####################
# Auto_drug_reversal_local_cmap.R
#
# Local CMap-style transcriptomic reversal ranking from ASGARD tissue-specific
# rank matrices, used as a fallback secondary pipeline when scDrugPrio or
# CLUE batch submission is not operational.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
})

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

base_dir <- "Auto_drug_reversal"
out_dir <- file.path(base_dir, "clue_fallback")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

signature_path <- file.path(base_dir, "Auto_drug_reversal_signature_top150.csv")
ref_cfg_path <- file.path(base_dir, "asgard_reference", "Auto_asgard_reference_paths.csv")

write_status <- function(status, detail) {
  fwrite(
    data.frame(
      step = "clue_local_rankmatrix",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_clue_local_rankmatrix_status.csv")
  )
}

if (!file.exists(signature_path)) {
  write_status("missing_inputs", "Missing signature file.")
  quit(save = "no", status = 0)
}

signature_dt <- fread(signature_path)
default_ref_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference"
if (file.exists(ref_cfg_path)) {
  ref_cfg <- fread(ref_cfg_path)
  rank_path <- ref_cfg$AUTO_ASGARD_DRUG_RESPONSE[1]
  gene_path <- ref_cfg$AUTO_ASGARD_GENE_INFO[1]
  drug_path <- ref_cfg$AUTO_ASGARD_DRUG_INFO[1]
} else {
  rank_path <- file.path(default_ref_dir, "stomach_rankMatrix.txt")
  gene_path <- file.path(default_ref_dir, "stomach_gene_info.txt")
  drug_path <- file.path(default_ref_dir, "stomach_drug_info.txt")
}

if (!all(file.exists(c(rank_path, gene_path, drug_path)))) {
  write_status("missing_reference", "Rank matrix, gene info, or drug info path does not exist.")
  quit(save = "no", status = 0)
}

message("Loading stomach rank matrix for local reversal scoring.")
rank_dt <- fread(rank_path)
gene_dt <- fread(gene_path)
drug_dt <- fread(drug_path)

setnames(gene_dt, old = names(gene_dt)[1:2], new = c("probe_id", "gene_symbol"))
instance_cols <- setdiff(colnames(rank_dt), "probe_id")
rank_dt <- merge(rank_dt, gene_dt, by = "probe_id", all.x = FALSE, all.y = FALSE)
rank_dt <- rank_dt[!is.na(gene_symbol) & nzchar(gene_symbol)]

rank_gene <- rank_dt[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = gene_symbol,
  .SDcols = instance_cols
]

rank_mat <- as.matrix(rank_gene[, ..instance_cols])
rownames(rank_mat) <- rank_gene$gene_symbol
norm_mat <- (rank_mat - 1) / pmax(nrow(rank_mat) - 1, 1)

instance_map <- data.table(instance_id = as.character(drug_dt$instance_id))
if ("cmap_name" %in% names(drug_dt) && "catalog_name" %in% names(drug_dt)) {
  instance_map[, drug := mapply(
    FUN = function(cmap_name, catalog_name) {
      sub(paste0("_", catalog_name, "$"), "", cmap_name)
    },
    as.character(drug_dt$cmap_name),
    as.character(drug_dt$catalog_name),
    USE.NAMES = FALSE
  )]
} else if ("cmap_name" %in% names(drug_dt)) {
  instance_map[, drug := sub("_[^_]+$", "", as.character(drug_dt$cmap_name))]
} else {
  instance_map[, drug := as.character(instance_id)]
}

annotate_targets <- function(drugs) {
  api_key <- Sys.getenv("CLUE_API_KEY", Sys.getenv("CLUE_KEY", ""))
  if (!nzchar(api_key) || !requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    return(data.table(drug = character(), target_genes = character(), moa = character()))
  }

  drugs <- unique(drugs[!is.na(drugs) & nzchar(drugs)])
  if (length(drugs) == 0) return(data.table(drug = character(), target_genes = character(), moa = character()))

  chunks <- split(drugs, ceiling(seq_along(drugs) / 50))
  out <- rbindlist(lapply(chunks, function(chunk) {
    filter_txt <- jsonlite::toJSON(
      list(
        where = list(pert_iname = list(inq = as.list(chunk))),
        fields = list(pert_iname = TRUE, target = TRUE, moa = TRUE),
        limit = 1000
      ),
      auto_unbox = TRUE
    )
    res <- httr::GET(
      "https://api.clue.io/api/perts",
      query = list(filter = filter_txt),
      httr::add_headers(user_key = api_key, Accept = "application/json")
    )
    if (httr::status_code(res) >= 300) {
      return(data.table(drug = character(), target_gene = character(), moa = character()))
    }
    txt <- httr::content(res, as = "text", encoding = "UTF-8")
    x <- jsonlite::fromJSON(txt)
    if (length(x) == 0) {
      return(data.table(drug = character(), target_gene = character(), moa = character()))
    }

    per_pert <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
      targets_i <- x$target[[i]]
      moas_i <- x$moa[[i]]
      data.table(
        drug = as.character(x$pert_iname[[i]]),
        target_gene = if (length(targets_i) == 0) NA_character_ else as.character(targets_i),
        moa = if (length(moas_i) == 0) NA_character_ else paste(sort(unique(as.character(moas_i))), collapse = ";")
      )
    }), fill = TRUE)

    per_pert
  }), fill = TRUE)

  out[
    !is.na(drug) & nzchar(drug),
    .(
      target_genes = paste(sort(unique(target_gene[!is.na(target_gene) & nzchar(target_gene)])), collapse = ";"),
      moa = paste(sort(unique(moa[!is.na(moa) & nzchar(moa)])), collapse = ";")
    ),
    by = drug
  ]
}

state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)
state_order <- intersect(state_order, unique(signature_dt$state))
ranked_all <- rbindlist(lapply(state_order, function(state_name) {
  up_genes <- unique(signature_dt[state == state_name & direction == "up", gene])
  down_genes <- unique(signature_dt[state == state_name & direction == "down", gene])
  up_genes <- intersect(up_genes, rownames(norm_mat))
  down_genes <- intersect(down_genes, rownames(norm_mat))

  if (length(up_genes) < 10 || length(down_genes) < 10) return(NULL)

  up_score <- colMeans(norm_mat[up_genes, , drop = FALSE], na.rm = TRUE)
  down_score <- colMeans(norm_mat[down_genes, , drop = FALSE], na.rm = TRUE)
  instance_score <- up_score - down_score

  dt <- data.table(
    instance_id = names(instance_score),
    instance_score = as.numeric(instance_score)
  )
  dt <- merge(dt, instance_map, by = "instance_id", all.x = TRUE)
  dt[is.na(drug) | !nzchar(drug), drug := instance_id]

  dt <- dt[
    ,
    .(
      score = mean(instance_score, na.rm = TRUE),
      best_instance_score = max(instance_score, na.rm = TRUE),
      n_instances = .N
    ),
    by = drug
  ][order(-score, -best_instance_score, drug)]

  dt[, `:=`(
    state = state_name,
    pipeline = "CLUE_FALLBACK_LOCAL",
    rank = seq_len(.N),
    p_value = NA_real_,
    fdr = NA_real_,
    moa = NA_character_
  )]

  dt[, .(state, pipeline, drug, rank, score, p_value, fdr, moa, n_instances, best_instance_score)]
}), fill = TRUE)

if (nrow(ranked_all) == 0) {
  write_status("no_results", "No local reversal scores could be computed.")
  quit(save = "no", status = 0)
}

top_targets <- annotate_targets(unique(ranked_all[rank <= 100, drug]))
if (nrow(top_targets) > 0) {
  ranked_all[, drug_key := tolower(trimws(drug))]
  top_targets[, drug_key := tolower(trimws(drug))]
  ranked_all <- merge(
    ranked_all,
    top_targets[, .(drug_key, target_genes, moa_annotation = moa)],
    by = "drug_key",
    all.x = TRUE
  )
  ranked_all[
    is.na(moa) | !nzchar(moa),
    moa := moa_annotation
  ]
  ranked_all[, moa_annotation := NULL]
  ranked_all[, drug_key := NULL]
} else {
  ranked_all[, `:=`(target_genes = NA_character_, moa = NA_character_)]
}

setcolorder(ranked_all, c("state", "pipeline", "drug", "rank", "score", "p_value", "fdr", "moa", "target_genes", "n_instances", "best_instance_score"))
fwrite(ranked_all, file.path(out_dir, "Auto_clue_ranked_drugs.csv"))
write_status("complete", paste("Ranked", nrow(ranked_all), "local fallback drug-state rows."))

message("Local CMap fallback ranking complete.")
