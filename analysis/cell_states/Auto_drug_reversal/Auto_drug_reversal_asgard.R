####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_asgard.R
#   Methodology: not required (superseded exploratory drug-reversal branch)
#   Map: analysis/ANALYSIS_MAP.md
####################

####################
# Auto_drug_reversal_asgard.R
#
# Run ASGARD mono-drug reversal from prepared scRef state-vs-rest DEG lists.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

####################
# setup
####################

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

out_dir <- file.path("Auto_drug_reversal", "asgard")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_dir <- file.path("Auto_drug_reversal", "asgard_inputs")
gene_list_path <- file.path(input_dir, "Auto_asgard_gene_list.rds")

state_order <- c(
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)

####################
# helpers
####################

write_status <- function(status, detail) {
  fwrite(
    data.frame(
      step = "asgard",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_asgard_status.csv")
  )
}

first_matching_col <- function(df, candidates) {
  nms <- colnames(df)
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit) > 0) return(hit[1])
  pattern <- paste(candidates, collapse = "|")
  hit <- grep(pattern, nms, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) hit[1] else NA_character_
}

standardize_asgard_table <- function(x, state_name) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.data.frame(x)) {
    x <- tryCatch(as.data.frame(x), error = function(e) NULL)
  }
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- rownames_to_column(x, "row_id")

  drug_col <- first_matching_col(
    x,
    c("drug", "Drug", "drug_name", "DrugName", "Drug.Name", "compound", "pert_iname", "name")
  )
  if (is.na(drug_col)) drug_col <- colnames(x)[2]

  fdr_col <- first_matching_col(x, c("FDR", "fdr", "adj.P.Val", "adj_p", "qvalue", "q.value"))
  p_col <- first_matching_col(x, c("P.Value", "p.value", "pvalue", "P", "p_val"))
  score_col <- first_matching_col(x, c("drug_score", "DrugScore", "score", "OS", "connectivity", "tau"))
  target_col <- first_matching_col(x, c("target", "targets", "target_genes", "Drug_targets", "gene_target"))
  moa_col <- first_matching_col(x, c("moa", "MOA", "mechanism", "mechanism_of_action"))

  rank_metric <- seq_len(nrow(x))
  rank_direction <- "row_order"
  if (!is.na(fdr_col)) {
    rank_metric <- suppressWarnings(as.numeric(x[[fdr_col]]))
    rank_direction <- "fdr_ascending"
  } else if (!is.na(p_col)) {
    rank_metric <- suppressWarnings(as.numeric(x[[p_col]]))
    rank_direction <- "p_ascending"
  } else if (!is.na(score_col)) {
    rank_metric <- -suppressWarnings(as.numeric(x[[score_col]]))
    rank_direction <- "score_descending"
  }
  if (any(!is.finite(rank_metric))) {
    finite_max <- suppressWarnings(max(rank_metric[is.finite(rank_metric)], na.rm = TRUE))
    if (!is.finite(finite_max)) finite_max <- nrow(x)
    rank_metric[!is.finite(rank_metric)] <- finite_max + 1
  }

  tibble(
    state = state_name,
    pipeline = "ASGARD",
    drug = as.character(x[[drug_col]]),
    rank_metric = rank_metric,
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(x[[score_col]])) else NA_real_,
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(x[[p_col]])) else NA_real_,
    fdr = if (!is.na(fdr_col)) suppressWarnings(as.numeric(x[[fdr_col]])) else NA_real_,
    moa = if (!is.na(moa_col)) as.character(x[[moa_col]]) else NA_character_,
    target_genes = if (!is.na(target_col)) as.character(x[[target_col]]) else NA_character_,
    rank_direction = rank_direction
  ) %>%
    filter(!is.na(drug), nzchar(drug)) %>%
    arrange(rank_metric, drug) %>%
    mutate(rank = row_number()) %>%
    select(state, pipeline, drug, rank, score, p_value, fdr, moa, target_genes, rank_direction)
}

####################
# validation
####################

if (!file.exists(gene_list_path)) {
  write_status("missing_inputs", paste("Missing", gene_list_path, "Run Auto_drug_reversal_inputs.R first."))
  quit(save = "no", status = 0)
}

if (!requireNamespace("Asgard", quietly = TRUE)) {
  write_status(
    "missing_package",
    "R package Asgard is not installed in this environment. Run Auto_setup_drug_reversal_env.sh or install lanagarmire/Asgard."
  )
  quit(save = "no", status = 0)
}

drug_ref_rds <- Sys.getenv("AUTO_ASGARD_DRUG_REF_RDS", "")
drug_response_path <- Sys.getenv("AUTO_ASGARD_DRUG_RESPONSE", "")
gene_info_path <- Sys.getenv("AUTO_ASGARD_GENE_INFO", "")
drug_info_path <- Sys.getenv("AUTO_ASGARD_DRUG_INFO", "")

reference_config <- file.path("Auto_drug_reversal", "asgard_reference", "Auto_asgard_reference_paths.csv")
default_ref_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference"
if (!nzchar(drug_ref_rds) &&
    (!file.exists(drug_response_path) || !file.exists(gene_info_path) || !file.exists(drug_info_path)) &&
    file.exists(reference_config)) {
  cfg <- fread(reference_config)
  drug_response_path <- cfg$AUTO_ASGARD_DRUG_RESPONSE[1]
  gene_info_path <- cfg$AUTO_ASGARD_GENE_INFO[1]
  drug_info_path <- cfg$AUTO_ASGARD_DRUG_INFO[1]
}
if (!nzchar(drug_ref_rds) &&
    (!file.exists(drug_response_path) || !file.exists(gene_info_path) || !file.exists(drug_info_path))) {
  drug_response_path <- file.path(default_ref_dir, "stomach_rankMatrix.txt")
  gene_info_path <- file.path(default_ref_dir, "stomach_gene_info.txt")
  drug_info_path <- file.path(default_ref_dir, "stomach_drug_info.txt")
}

if (!nzchar(drug_ref_rds) &&
    (!file.exists(drug_response_path) || !file.exists(gene_info_path) || !file.exists(drug_info_path))) {
  write_status(
    "missing_reference",
    paste(
      "ASGARD requires either AUTO_ASGARD_DRUG_REF_RDS or",
      "AUTO_ASGARD_DRUG_RESPONSE + AUTO_ASGARD_GENE_INFO + AUTO_ASGARD_DRUG_INFO.",
      "These are tissue-specific L1000 rankMatrix/gene/drug reference files produced by ASGARD PrepareReference()."
    )
  )
  quit(save = "no", status = 0)
}

####################
# ASGARD execution
####################

message("Loading ASGARD gene list and drug reference.")
gene_list <- readRDS(gene_list_path)
gene_list <- gene_list[intersect(state_order, names(gene_list))]

if (nzchar(drug_ref_rds)) {
  drug_ref_profiles <- readRDS(drug_ref_rds)
} else {
  gene_info <- read.table(gene_info_path, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
  drug_info <- read.table(drug_info_path, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
  drug_ref_profiles <- Asgard::GetDrugRef(
    drug.response.path = drug_response_path,
    probe.to.genes = gene_info,
    drug.info = drug_info
  )
}

message("Running ASGARD GetDrug() with negative connectivity.")
drug_type <- Sys.getenv("AUTO_ASGARD_DRUG_TYPE", "all")
drug_ident_res <- Asgard::GetDrug(
  gene.data = gene_list,
  drug.ref.profiles = drug_ref_profiles,
  repurposing.unit = "drug",
  connectivity = "negative",
  drug.type = drug_type
)

saveRDS(drug_ident_res, file.path(out_dir, "Auto_asgard_raw_drug_identification.rds"))

ranked <- bind_rows(lapply(names(drug_ident_res), function(state_name) {
  standardize_asgard_table(drug_ident_res[[state_name]], state_name)
}))

if (nrow(ranked) == 0) {
  write_status("no_results", "ASGARD completed but no rankable drug rows were returned.")
  quit(save = "no", status = 0)
}

fwrite(ranked, file.path(out_dir, "Auto_asgard_ranked_drugs.csv"))
write_status("complete", paste("Ranked", nrow(ranked), "ASGARD drug-state rows."))

message("ASGARD drug reversal complete.")
