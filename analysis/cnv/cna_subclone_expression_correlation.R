####################
# Analysis registry:
#   Status: active, terminal
#   Script: analysis/cnv/cna_subclone_expression_correlation.R
#   Methodology: analysis/methodology/cnv/cna_subclone_expression_correlation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# cna_subclone_expression_correlation.R
#
# Self-contained CNA subclone expression analysis (merged from v1+v2 split):
#   Phase 1: computation — arm-level CNA, dominant clone, pairwise distance,
#            consensus clustering, OAC/OCCAMS/TCGA annotation, recurrent events
#   Phase 2: presentation-quality plots — corrected event recounting, standardized
#            feature boxplots, per-sample heatmaps, association dot plots
#
# Inputs:
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_summary.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_mp_subclone_tests.csv
#   ref_outs/by_samples/<sample>/<sample>_outs.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/OAC_CNV_summary.xlsx
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/OAC_CNV_all.xlsx
#
# Outputs:
#   ref_outs/Auto_cna_subclone_expression/
#
# Cache: saves intermediate RDS after computation.
#   Set SCREF_REPLOT_ONLY=TRUE to skip computation and replot from cache.
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(scales)
  library(readxl)
  library(jsonlite)
  library(httr)
  library(grid)
  library(gridExtra)
  library(ggrepel)
  library(cluster)
  library(ggpubr)
})

options(stringsAsFactors = FALSE)
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
robust_effect_size_threshold <- if (length(args) >= 1 && nzchar(args[1])) as.numeric(args[1]) else 1.0

replot_only <- identical(Sys.getenv("SCREF_REPLOT_ONLY"), "TRUE")

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

out_dir <- file.path("Auto_cna_subclone_expression", paste0("threshold_", robust_effect_size_threshold))
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
live_rds_dir <- file.path(out_dir, "rds")
ephemeral_out_dir <- file.path("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_cna_subclone_expression", paste0("threshold_", robust_effect_size_threshold))
ephemeral_rds_dir <- file.path(ephemeral_out_dir, "rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(live_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ephemeral_rds_dir, recursive = TRUE, showWarnings = FALSE)

subclone_dir <- file.path("Auto_malignant_subclone_mp", paste0("threshold_", robust_effect_size_threshold))
cell_path <- file.path(subclone_dir, "Auto_malignant_subclone_cells.csv")
summary_path <- file.path(subclone_dir, "Auto_malignant_subclone_summary.csv")
mp_subclone_path <- file.path(subclone_dir, "Auto_malignant_subclone_mp_subclone_tests.csv")
oac_cnv_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/OAC_CNV_summary.xlsx"
occams_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/OAC_CNV_all.xlsx"

required_paths <- c(cell_path, summary_path, mp_subclone_path, oac_cnv_path, occams_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = ", "))
}

gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)

arm_call_threshold <- 0.10
dominant_min_fraction <- 0.50
dominant_min_gap <- 0.15
recurrent_min_sample_fraction <- 0.15
recurrent_min_samples <- 3L

chrom_levels <- c(paste0("chr", 1:22), "chrX")
arm_levels <- as.vector(rbind(paste0(chrom_levels, "p"), paste0(chrom_levels, "q")))
chr_from_arm <- sub("[pq]$", "", arm_levels)

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

mp_descriptions <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)

mp_order <- c("MP1", "MP5", "MP13+", "MP2+", "MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+",
              "MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP12", "MP15")

state_order <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia",
  "Stress adaptive",
  "Cancer-cell immune mimicry",
  "Unresolved",
  "Hybrid"
)

state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

mp_cols <- c(
  "MP1" = "#B0B0B0",
  "MP5" = "#C0C0C0",
  "MP13+" = "#999999",
  "MP2+" = "#E41A1C",
  "MP14" = "#8DA0CB",
  "MP3+" = "#4DAF4A",
  "MP6+" = "#66C2A5",
  "MP11+" = "#A6D854",
  "MP9+" = "#FC8D62",
  "MP10+" = "#FF7F00",
  "MP8+" = "#FFD92F",
  "MP8b" = "#E78AC3",
  "MP16" = "#FFD92F",
  "MP18b" = "#FF7F00",
  "MP17" = "#4DAF4A",
  "MP12" = "#984EA3",
  "MP15" = "#377EB8"
)

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(!is.na(x)) < 2) return(NA_real_)
  sd(x, na.rm = TRUE)
}

clean_feature_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  x
}

infer_study <- function(sample_id) {
  sub("^([^_]+_[0-9]{4}).*$", "\\1", as.character(sample_id))
}

write_table <- function(x, filename) {
  write.csv(x, file.path(table_dir, filename), row.names = FALSE)
}

save_rds <- function(x, filename) {
  saveRDS(x, file.path(ephemeral_rds_dir, filename))
  saveRDS(x, file.path(live_rds_dir, filename))
}

label_mp <- function(mp) {
  desc <- mp_descriptions[mp]
  desc[is.na(desc)] <- mp[is.na(desc)]
  paste0(mp, " ", desc)
}

feature_label <- function(feature) {
  out <- feature
  out <- ifelse(grepl("^mp__", out), label_mp(sub("^mp__", "", out)), out)
  out <- ifelse(grepl("^state__", out), paste0("State: ", gsub("_", " ", sub("^state__", "", out))), out)
  out <- gsub("_", " ", out)
  out <- gsub("^nCount RNA$", "nCount_RNA", out)
  out <- gsub("^nFeature RNA$", "nFeature_RNA", out)
  out <- gsub("^percent mt$", "percent.mt", out)
  out <- gsub("^cc score$", "cc_score", out)
  out <- gsub("^cs score$", "cs_score", out)
  out <- gsub("^cell cycle mp mean$", "Cell-cycle MP mean", out)
  out <- gsub("^cna burden mean abs$", "CNA burden", out)
  out
}

scale_rows <- function(mat) {
  mat <- as.matrix(mat)
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_weighted_mean <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

p_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values)]
  missing <- setdiff(values, names(cols))
  if (length(missing) > 0) {
    extra <- setNames(colorRampPalette(brewer.pal(8, palette))(length(missing)), missing)
    cols <- c(cols, extra)
  }
  cols[values]
}

subclone_palette <- c(
  "Subclone A" = "#D73027", "Subclone B" = "#4575B4", "Subclone C" = "#1A9850",
  "Subclone D" = "#984EA3", "Subclone E" = "#FF7F00", "Subclone F" = "#A65628"
)

subclone_colours <- function(values) {
  cols <- complete_palette(subclone_palette, values, "Set2")
  cols[!is.na(names(cols)) & nzchar(names(cols))]
}

state_display <- c(
  "Classic proliferation" = "state__Classic_proliferation",
  "Basal to intestinal metaplasia" = "state__Basal_to_intestinal_metaplasia",
  "SMG to intestinal metaplasia" = "state__SMG_to_intestinal_metaplasia",
  "Stress adaptive" = "state__Stress_adaptive",
  "Cancer-cell immune mimicry" = "state__Cancer_cell_immune_mimicry"
)

event_label <- function(event_id, known_genes = NULL) {
  direction <- sub("_.*$", "", event_id)
  arm <- sub("^[^_]+_", "", event_id)
  base <- paste0(ifelse(direction == "gain", "Gain ", "Loss "), arm)
  if (!is.null(known_genes)) {
    kg <- as.character(known_genes)
    kg[is.na(kg)] <- ""
    kg <- ifelse(nzchar(kg), paste0(" (", substr(kg, 1, 36), ")"), "")
    return(paste0(base, kg))
  }
  base
}

cytoband_to_arm <- function(cytoband) {
  cytoband <- as.character(cytoband)
  cytoband <- gsub("^chr", "", cytoband, ignore.case = TRUE)
  m <- regexec("^([0-9]+|X|Y)([pq])", cytoband)
  parts <- regmatches(cytoband, m)
  vapply(parts, function(x) {
    if (length(x) < 3) return(NA_character_)
    chr <- paste0("chr", x[2])
    if (!chr %in% chrom_levels) return(NA_character_)
    paste0(chr, x[3])
  }, character(1))
}

parse_peak_coord <- function(x) {
  x <- as.character(x)
  m <- regexec("(chr[0-9XY]+):([0-9]+)-([0-9]+)", x)
  parts <- regmatches(x, m)
  if (length(parts[[1]]) < 4) {
    return(data.frame(chr = NA_character_, start = NA_real_, end = NA_real_, arm = NA_character_))
  }
  chr <- parts[[1]][2]
  start <- as.numeric(parts[[1]][3])
  end <- as.numeric(parts[[1]][4])
  mid <- mean(c(start, end), na.rm = TRUE)
  arm <- if (chr %in% names(centromere_pos)) {
    paste0(chr, ifelse(mid <= centromere_pos[[chr]], "p", "q"))
  } else {
    NA_character_
  }
  data.frame(chr = chr, start = start, end = end, arm = arm)
}

split_gene_string <- function(x) {
  x <- as.character(x)
  x <- gsub("\\([^)]*\\)", "", x)
  x <- gsub("\\[|\\]", "", x)
  pieces <- unlist(strsplit(x, "[,;/ ]+"))
  pieces <- toupper(trimws(pieces))
  pieces <- pieces[nzchar(pieces) & !is.na(pieces)]
  pieces <- pieces[!grepl("^HSA-", pieces)]
  unique(pieces)
}

gene_order <- read.table(
  gene_order_path,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "end"),
  stringsAsFactors = FALSE
) %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start) %>%
  distinct(.data$gene, .keep_all = TRUE) %>%
  mutate(
    chromosome = as.character(.data$chromosome),
    arm = ifelse(.data$start <= centromere_pos[.data$chromosome], "p", "q"),
    chr_arm = paste0(.data$chromosome, .data$arm)
  )

gene_to_arm <- gene_order %>%
  select(gene = .data$gene, chr_arm = .data$chr_arm, chromosome = .data$chromosome, start = .data$start, end = .data$end)

cache_rds_path <- file.path(live_rds_dir, "Auto_cna_subclone_expression_results.rds")
if (!file.exists(cache_rds_path)) {
  ephemeral_cache_path <- file.path(ephemeral_rds_dir, "Auto_cna_subclone_expression_results.rds")
  if (file.exists(ephemeral_cache_path)) cache_rds_path <- ephemeral_cache_path
}

if (replot_only && file.exists(cache_rds_path)) {
  message("SCREF_REPLOT_ONLY=TRUE: loading cached computation results")
  res <- readRDS(cache_rds_path)
  features <- as.data.frame(res$features)
  arm_long <- as.data.frame(res$arm_long)
  arm_matrix <- as.matrix(res$arm_matrix)
  arm_call_matrix <- as.matrix(res$arm_call_matrix)
  dominance_df <- as.data.frame(res$dominance)
  pairwise_df <- as.data.frame(res$pairwise)
  event_summary_annotated <- as.data.frame(res$event_summary)
  cna_cluster <- res$cna_cluster
  cell_path <- file.path(subclone_dir, "Auto_malignant_subclone_cells.csv")
  cell_df <- read.csv(cell_path, check.names = FALSE)
} else {

message("Loading completed subclone outputs")
cell_df <- read.csv(cell_path, check.names = FALSE)
sample_summary <- read.csv(summary_path, check.names = FALSE)
mp_subclone_df <- read.csv(mp_subclone_path, check.names = FALSE)

cell_df <- cell_df %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = paste(.data$sample, .data$subclone, sep = "::"),
    study = infer_study(.data$sample),
    state_label = as.character(.data$state_label),
    top_mp = as.character(.data$top_mp)
  )

sample_summary <- sample_summary %>%
  mutate(sample = as.character(.data$sample), study = infer_study(.data$sample))

mp_subclone_df <- mp_subclone_df %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = paste(.data$sample, .data$subclone, sep = "::"),
    mp = as.character(.data$mp)
  )

analysed_samples <- sample_summary %>%
  filter(.data$status == "analysed", .data$n_subclones >= 1) %>%
  pull(.data$sample)

cell_df <- cell_df %>% filter(.data$sample %in% analysed_samples)
if (nrow(cell_df) == 0) stop("No analysed malignant subclone cells found.")

message("Computing chromosome-arm CNA means per subclone")

compute_sample_arm_cna <- function(sample_id) {
  cells_one <- cell_df %>% filter(.data$sample == sample_id)
  outs_path <- file.path("by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  if (!file.exists(outs_path)) {
    warning("Missing CNA matrix for ", sample_id)
    return(NULL)
  }
  outs <- readRDS(outs_path)
  cells <- intersect(cells_one$cell, colnames(outs))
  if (length(cells) == 0) return(NULL)

  common_genes <- intersect(rownames(outs), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE] %>%
    arrange(factor(.data$chromosome, levels = chrom_levels), .data$start)
  mat <- as.matrix(outs[go$gene, cells, drop = FALSE])
  mat[!is.finite(mat)] <- NA_real_
  arm_by_gene <- factor(go$chr_arm, levels = arm_levels)
  subclone_vec <- cells_one$subclone[match(cells, cells_one$cell)]
  names(subclone_vec) <- cells
  subclone_levels <- sort(unique(subclone_vec))

  rows <- list()
  row_i <- 1L
  for (cl in subclone_levels) {
    cl_cells <- names(subclone_vec)[subclone_vec == cl]
    arm_vals <- vapply(arm_levels, function(arm) {
      idx <- which(arm_by_gene == arm)
      if (length(idx) == 0 || length(cl_cells) == 0) return(NA_real_)
      safe_mean(mat[idx, cl_cells, drop = FALSE])
    }, numeric(1))
    rows[[row_i]] <- data.frame(
      sample = sample_id,
      study = infer_study(sample_id),
      subclone = cl,
      subclone_id = paste(sample_id, cl, sep = "::"),
      arm = arm_levels,
      chr = chr_from_arm,
      arm_mean = as.numeric(arm_vals),
      n_cells = length(cl_cells),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
  bind_rows(rows)
}

arm_long <- bind_rows(lapply(analysed_samples, function(s) {
  message("  ", s)
  compute_sample_arm_cna(s)
}))

if (nrow(arm_long) == 0) stop("No arm-level CNA profiles could be computed.")

arm_long <- arm_long %>%
  mutate(
    arm = factor(.data$arm, levels = arm_levels),
    cna_call = case_when(
      is.na(.data$arm_mean) ~ NA_integer_,
      .data$arm_mean >= arm_call_threshold ~ 1L,
      .data$arm_mean <= -arm_call_threshold ~ -1L,
      TRUE ~ 0L
    ),
    call_label = case_when(
      .data$cna_call == 1L ~ "gain",
      .data$cna_call == -1L ~ "loss",
      .data$cna_call == 0L ~ "neutral",
      TRUE ~ NA_character_
    )
  )

write_table(arm_long, "Auto_subclone_arm_cna_long.csv")
save_rds(arm_long, "Auto_subclone_arm_cna_long.rds")

cna_summary <- arm_long %>%
  group_by(.data$sample, .data$study, .data$subclone, .data$subclone_id) %>%
  summarise(
    cna_burden_mean_abs = mean(abs(.data$arm_mean), na.rm = TRUE),
    cna_burden_sd = safe_sd(.data$arm_mean),
    n_gain_arms = sum(.data$cna_call == 1L, na.rm = TRUE),
    n_loss_arms = sum(.data$cna_call == -1L, na.rm = TRUE),
    frac_gain_arms = mean(.data$cna_call == 1L, na.rm = TRUE),
    frac_loss_arms = mean(.data$cna_call == -1L, na.rm = TRUE),
    .groups = "drop"
  )

arm_wide <- arm_long %>%
  mutate(arm = as.character(.data$arm)) %>%
  select(.data$subclone_id, .data$arm, .data$arm_mean) %>%
  pivot_wider(names_from = .data$arm, values_from = .data$arm_mean, values_fill = 0)

arm_call_wide <- arm_long %>%
  mutate(arm = as.character(.data$arm)) %>%
  select(.data$subclone_id, .data$arm, .data$cna_call) %>%
  pivot_wider(names_from = .data$arm, values_from = .data$cna_call, values_fill = 0)

arm_matrix <- as.matrix(arm_wide[, arm_levels, drop = FALSE])
rownames(arm_matrix) <- arm_wide$subclone_id
arm_call_matrix <- as.matrix(arm_call_wide[, arm_levels, drop = FALSE])
rownames(arm_call_matrix) <- arm_call_wide$subclone_id
arm_matrix[!is.finite(arm_matrix)] <- 0
arm_call_matrix[!is.finite(arm_call_matrix)] <- 0

write_table(arm_wide, "Auto_subclone_arm_cna_matrix.csv")
write_table(arm_call_wide, "Auto_subclone_arm_cna_call_matrix.csv")
save_rds(arm_matrix, "Auto_subclone_arm_cna_matrix.rds")
save_rds(arm_call_matrix, "Auto_subclone_arm_cna_call_matrix.rds")

message("Building subclone feature table")

qc_cols <- intersect(
  c("cna_signal", "cna_cor", "nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score"),
  colnames(cell_df)
)

cell_summary <- cell_df %>%
  group_by(.data$sample, .data$study, .data$subclone, .data$subclone_id) %>%
  summarise(
    n_cells = n(),
    across(all_of(qc_cols), safe_mean),
    .groups = "drop"
  ) %>%
  group_by(.data$sample) %>%
  mutate(
    sample_malignant_cells = sum(.data$n_cells),
    subclone_fraction = .data$n_cells / .data$sample_malignant_cells
  ) %>%
  ungroup()

state_prop <- cell_df %>%
  mutate(state_label = ifelse(is.na(.data$state_label) | !nzchar(.data$state_label), "Unknown", .data$state_label)) %>%
  count(.data$sample, .data$subclone, .data$subclone_id, .data$state_label, name = "n") %>%
  group_by(.data$sample, .data$subclone, .data$subclone_id) %>%
  mutate(prop = .data$n / sum(.data$n)) %>%
  ungroup() %>%
  mutate(state_feature = paste0("state__", clean_feature_name(.data$state_label))) %>%
  select(.data$subclone_id, .data$state_feature, .data$prop) %>%
  pivot_wider(names_from = .data$state_feature, values_from = .data$prop, values_fill = 0)

top_mp_prop <- cell_df %>%
  mutate(top_mp = ifelse(is.na(.data$top_mp) | !nzchar(.data$top_mp), "Unknown", .data$top_mp)) %>%
  count(.data$sample, .data$subclone, .data$subclone_id, .data$top_mp, name = "n") %>%
  group_by(.data$sample, .data$subclone, .data$subclone_id) %>%
  mutate(prop = .data$n / sum(.data$n)) %>%
  ungroup() %>%
  mutate(top_mp_feature = paste0("topmp__", clean_feature_name(.data$top_mp))) %>%
  select(.data$subclone_id, .data$top_mp_feature, .data$prop) %>%
  pivot_wider(names_from = .data$top_mp_feature, values_from = .data$prop, values_fill = 0)

mp_wide <- mp_subclone_df %>%
  select(.data$sample, .data$subclone, .data$subclone_id, .data$mp, .data$mean_score) %>%
  distinct() %>%
  filter(.data$mp %in% mp_order) %>%
  mutate(mp_feature = paste0("mp__", .data$mp)) %>%
  select(.data$subclone_id, .data$mp_feature, .data$mean_score) %>%
  pivot_wider(names_from = .data$mp_feature, values_from = .data$mean_score)

features <- cell_summary %>%
  left_join(cna_summary, by = c("sample", "study", "subclone", "subclone_id")) %>%
  left_join(mp_wide, by = "subclone_id") %>%
  left_join(state_prop, by = "subclone_id") %>%
  left_join(top_mp_prop, by = "subclone_id")

cc_mp_cols <- intersect(paste0("mp__", c("MP1", "MP5", "MP13+")), colnames(features))
if (length(cc_mp_cols) > 0) {
  features$cell_cycle_mp_mean <- rowMeans(features[, cc_mp_cols, drop = FALSE], na.rm = TRUE)
}

features <- features %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

state_feature_cols <- grep("^state__", colnames(features), value = TRUE)
mp_feature_cols <- grep("^mp__", colnames(features), value = TRUE)
topmp_feature_cols <- grep("^topmp__", colnames(features), value = TRUE)

write_table(features, "Auto_subclone_feature_summary.csv")
save_rds(features, "Auto_subclone_feature_summary.rds")

message("Testing dominant clone associations")

dominance_rows <- lapply(split(features, features$sample), function(df) {
  df <- df %>% arrange(desc(.data$n_cells), .data$subclone)
  n_sub <- nrow(df)
  top_n <- df$n_cells[1]
  top_frac <- df$subclone_fraction[1]
  second_n <- if (n_sub >= 2) df$n_cells[2] else NA_real_
  second_frac <- if (n_sub >= 2) df$subclone_fraction[2] else NA_real_
  p_value <- if (n_sub >= 2) {
    tryCatch(binom.test(top_n, top_n + second_n, p = 0.5, alternative = "greater")$p.value,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }
  data.frame(
    sample = df$sample[1],
    study = df$study[1],
    n_subclones = n_sub,
    dominant_subclone = df$subclone[1],
    dominant_n_cells = top_n,
    dominant_fraction = top_frac,
    second_subclone = if (n_sub >= 2) df$subclone[2] else NA_character_,
    second_n_cells = second_n,
    second_fraction = second_frac,
    dominant_gap = if (n_sub >= 2) top_frac - second_frac else NA_real_,
    dominant_ratio = if (n_sub >= 2 && second_frac > 0) top_frac / second_frac else NA_real_,
    binom_p_value = p_value,
    stringsAsFactors = FALSE
  )
})

dominance_df <- bind_rows(dominance_rows) %>%
  mutate(
    binom_p_adj = ifelse(
      !is.na(.data$binom_p_value),
      p.adjust(.data$binom_p_value, method = "BH"),
      NA_real_
    ),
    dominance_class = case_when(
      .data$n_subclones == 1 ~ "single_subclone",
      !is.na(.data$binom_p_adj) &
        .data$binom_p_adj < 0.05 &
        .data$dominant_fraction >= dominant_min_fraction &
        .data$dominant_gap >= dominant_min_gap ~ "significant_dominant",
      TRUE ~ "largest_not_significant"
    )
  )

features <- features %>%
  left_join(dominance_df %>% select(.data$sample, .data$dominant_subclone, .data$dominance_class),
            by = "sample") %>%
  mutate(
    is_largest_subclone = .data$subclone == .data$dominant_subclone,
    is_significant_dominant = .data$is_largest_subclone & .data$dominance_class == "significant_dominant"
  )

write_table(dominance_df, "Auto_dominant_clone_summary.csv")
write_table(features, "Auto_subclone_feature_summary_with_dominance.csv")

classic_state_col <- paste0("state__", clean_feature_name("Classic proliferation"))
basal_state_col <- paste0("state__", clean_feature_name("Basal to intestinal metaplasia"))
smg_state_col <- paste0("state__", clean_feature_name("SMG to intestinal metaplasia"))
stress_state_col <- paste0("state__", clean_feature_name("Stress adaptive"))
immune_state_col <- paste0("state__", clean_feature_name("Cancer-cell immune mimicry"))
hybrid_state_col <- paste0("state__", clean_feature_name("Hybrid"))
unresolved_state_col <- paste0("state__", clean_feature_name("Unresolved"))

target_features <- unique(c(
  paste0("mp__", mp_order),
  "cell_cycle_mp_mean",
  classic_state_col,
  basal_state_col,
  smg_state_col,
  stress_state_col,
  immune_state_col,
  hybrid_state_col,
  unresolved_state_col,
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "cc_score",
  "cs_score",
  "cna_signal",
  "cna_cor",
  "cna_burden_mean_abs",
  "n_gain_arms",
  "n_loss_arms",
  "frac_gain_arms",
  "frac_loss_arms"
))
target_features <- intersect(target_features, colnames(features))

dominant_feature_tests <- lapply(target_features, function(feat) {
  deltas <- lapply(split(features, features$sample), function(df) {
    if (nrow(df) < 2) return(NULL)
    top <- df %>% filter(.data$is_largest_subclone)
    rest <- df %>% filter(!.data$is_largest_subclone)
    if (nrow(top) != 1 || nrow(rest) == 0) return(NULL)
    rest_value <- weighted.mean(rest[[feat]], w = rest$n_cells, na.rm = TRUE)
    data.frame(
      sample = top$sample,
      feature = feat,
      dominant_value = top[[feat]],
      rest_weighted_value = rest_value,
      delta = top[[feat]] - rest_value,
      dominance_class = top$dominance_class,
      stringsAsFactors = FALSE
    )
  })
  d <- bind_rows(deltas)
  p_value <- if (nrow(d) >= 3 && sum(!is.na(d$delta)) >= 3) {
    tryCatch(wilcox.test(d$delta, mu = 0)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  p_value_sigdom <- if (nrow(d) >= 3 && sum(d$dominance_class == "significant_dominant", na.rm = TRUE) >= 3) {
    tryCatch(wilcox.test(d$delta[d$dominance_class == "significant_dominant"], mu = 0)$p.value,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }
  centered_df <- features %>%
    group_by(.data$sample) %>%
    mutate(
      centered_fraction = .data$subclone_fraction - mean(.data$subclone_fraction, na.rm = TRUE),
      centered_feature = .data[[feat]] - mean(.data[[feat]], na.rm = TRUE)
    ) %>%
    ungroup()
  cor_res <- if (sum(complete.cases(centered_df[, c("centered_fraction", "centered_feature")])) >= 5) {
    tryCatch(cor.test(centered_df$centered_fraction, centered_df$centered_feature, method = "spearman"),
             error = function(e) NULL)
  } else {
    NULL
  }
  data.frame(
    feature = feat,
    feature_label = feature_label(feat),
    n_samples = nrow(d),
    median_delta = if (nrow(d) > 0) median(d$delta, na.rm = TRUE) else NA_real_,
    mean_delta = if (nrow(d) > 0) mean(d$delta, na.rm = TRUE) else NA_real_,
    pct_positive_delta = if (nrow(d) > 0) mean(d$delta > 0, na.rm = TRUE) else NA_real_,
    wilcox_p_value = p_value,
    significant_dominant_p_value = p_value_sigdom,
    sample_centered_spearman_rho = if (!is.null(cor_res)) unname(cor_res$estimate) else NA_real_,
    sample_centered_spearman_p = if (!is.null(cor_res)) cor_res$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(
    wilcox_p_adj = p.adjust(.data$wilcox_p_value, method = "BH"),
    significant_dominant_p_adj = p.adjust(.data$significant_dominant_p_value, method = "BH"),
    sample_centered_spearman_p_adj = p.adjust(.data$sample_centered_spearman_p, method = "BH")
  ) %>%
  arrange(.data$wilcox_p_adj, desc(abs(.data$median_delta)))

write_table(dominant_feature_tests, "Auto_dominant_clone_feature_tests.csv")

dominant_feature_deltas <- bind_rows(lapply(target_features, function(feat) {
  bind_rows(lapply(split(features, features$sample), function(df) {
    if (nrow(df) < 2) return(NULL)
    top <- df %>% filter(.data$is_largest_subclone)
    rest <- df %>% filter(!.data$is_largest_subclone)
    if (nrow(top) != 1 || nrow(rest) == 0) return(NULL)
    data.frame(
      sample = top$sample,
      feature = feat,
      feature_label = feature_label(feat),
      dominant_value = top[[feat]],
      rest_weighted_value = weighted.mean(rest[[feat]], w = rest$n_cells, na.rm = TRUE),
      delta = top[[feat]] - weighted.mean(rest[[feat]], w = rest$n_cells, na.rm = TRUE),
      dominance_class = top$dominance_class,
      stringsAsFactors = FALSE
    )
  }))
}))
write_table(dominant_feature_deltas, "Auto_dominant_clone_feature_deltas.csv")

message("Computing pairwise CNA-expression distances")

feature_by_id <- features
rownames(feature_by_id) <- feature_by_id$subclone_id
all_pair_rows <- list()
pair_i <- 1L

js_distance <- function(p, q) {
  p[is.na(p)] <- 0
  q[is.na(q)] <- 0
  if (sum(p) == 0 && sum(q) == 0) return(0)
  if (sum(p) > 0) p <- p / sum(p)
  if (sum(q) > 0) q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl <- function(a, b) {
    idx <- a > 0 & b > 0
    sum(a[idx] * log2(a[idx] / b[idx]))
  }
  sqrt(0.5 * kl(p, m) + 0.5 * kl(q, m))
}

for (sample_id in unique(features$sample)) {
  ids <- features$subclone_id[features$sample == sample_id]
  if (length(ids) < 2) next
  combos <- combn(ids, 2)
  for (j in seq_len(ncol(combos))) {
    id1 <- combos[1, j]
    id2 <- combos[2, j]
    a1 <- arm_matrix[id1, arm_levels]
    a2 <- arm_matrix[id2, arm_levels]
    c1 <- arm_call_matrix[id1, arm_levels]
    c2 <- arm_call_matrix[id2, arm_levels]
    f1 <- feature_by_id[id1, , drop = FALSE]
    f2 <- feature_by_id[id2, , drop = FALSE]
    mp1 <- as.numeric(f1[, mp_feature_cols, drop = TRUE])
    mp2 <- as.numeric(f2[, mp_feature_cols, drop = TRUE])
    st1 <- as.numeric(f1[, state_feature_cols, drop = TRUE])
    st2 <- as.numeric(f2[, state_feature_cols, drop = TRUE])
    row <- data.frame(
      sample = sample_id,
      study = f1$study,
      subclone_1 = f1$subclone,
      subclone_2 = f2$subclone,
      subclone_id_1 = id1,
      subclone_id_2 = id2,
      n_cells_1 = f1$n_cells,
      n_cells_2 = f2$n_cells,
      fraction_1 = f1$subclone_fraction,
      fraction_2 = f2$subclone_fraction,
      cna_euclidean = sqrt(mean((a1 - a2)^2, na.rm = TRUE)),
      cna_abs_mean = mean(abs(a1 - a2), na.rm = TRUE),
      cna_max_abs = max(abs(a1 - a2), na.rm = TRUE),
      cna_call_discordance = mean(c1 != c2, na.rm = TRUE),
      mp_euclidean = sqrt(mean((mp1 - mp2)^2, na.rm = TRUE)),
      mp_abs_mean = mean(abs(mp1 - mp2), na.rm = TRUE),
      state_l1 = 0.5 * sum(abs(st1 - st2), na.rm = TRUE),
      state_js = js_distance(st1, st2),
      fraction_abs_delta = abs(f1$subclone_fraction - f2$subclone_fraction),
      stringsAsFactors = FALSE
    )
    for (feat in target_features) {
      if (feat %in% colnames(f1)) {
        row[[paste0("abs_delta__", feat)]] <- abs(as.numeric(f1[[feat]]) - as.numeric(f2[[feat]]))
      }
    }
    all_pair_rows[[pair_i]] <- row
    pair_i <- pair_i + 1L
  }
}

pairwise_df <- bind_rows(all_pair_rows)
write_table(pairwise_df, "Auto_pairwise_cna_expression_distances.csv")

pair_endpoints <- c(
  "mp_euclidean",
  "mp_abs_mean",
  "state_l1",
  "state_js",
  "fraction_abs_delta",
  paste0("abs_delta__", intersect(c("mp__MP2+", "cell_cycle_mp_mean", classic_state_col, "cc_score", "cs_score", "nCount_RNA", "cna_burden_mean_abs"), target_features))
)
pair_endpoints <- intersect(pair_endpoints, colnames(pairwise_df))
cna_distance_metrics <- c("cna_abs_mean", "cna_euclidean", "cna_call_discordance")

pairwise_tests <- bind_rows(lapply(cna_distance_metrics, function(cna_metric) {
  bind_rows(lapply(pair_endpoints, function(endpoint) {
    df <- pairwise_df %>% filter(is.finite(.data[[cna_metric]]), is.finite(.data[[endpoint]]))
    cor_res <- if (nrow(df) >= 5) {
      tryCatch(cor.test(df[[cna_metric]], df[[endpoint]], method = "spearman"), error = function(e) NULL)
    } else {
      NULL
    }
    centered <- df %>%
      group_by(.data$sample) %>%
      mutate(
        cna_centered = .data[[cna_metric]] - mean(.data[[cna_metric]], na.rm = TRUE),
        endpoint_centered = .data[[endpoint]] - mean(.data[[endpoint]], na.rm = TRUE)
      ) %>%
      ungroup()
    centered_cor <- if (nrow(centered) >= 5 && sd(centered$cna_centered, na.rm = TRUE) > 0 && sd(centered$endpoint_centered, na.rm = TRUE) > 0) {
      tryCatch(cor.test(centered$cna_centered, centered$endpoint_centered, method = "spearman"),
               error = function(e) NULL)
    } else {
      NULL
    }
    lm_res <- if (nrow(df) >= length(unique(df$sample)) + 3 && sd(df[[cna_metric]], na.rm = TRUE) > 0) {
      tryCatch({
        fit <- lm(as.formula(paste0(endpoint, " ~ ", cna_metric, " + factor(sample)")), data = df)
        coef(summary(fit))[cna_metric, ]
      }, error = function(e) NULL)
    } else {
      NULL
    }
    data.frame(
      cna_metric = cna_metric,
      endpoint = endpoint,
      endpoint_label = feature_label(gsub("^abs_delta__", "", endpoint)),
      n_pairs = nrow(df),
      spearman_rho = if (!is.null(cor_res)) unname(cor_res$estimate) else NA_real_,
      spearman_p_value = if (!is.null(cor_res)) cor_res$p.value else NA_real_,
      sample_centered_rho = if (!is.null(centered_cor)) unname(centered_cor$estimate) else NA_real_,
      sample_centered_p_value = if (!is.null(centered_cor)) centered_cor$p.value else NA_real_,
      sample_fixed_beta = if (!is.null(lm_res)) unname(lm_res["Estimate"]) else NA_real_,
      sample_fixed_p_value = if (!is.null(lm_res)) unname(lm_res["Pr(>|t|)"]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  group_by(.data$cna_metric) %>%
  mutate(
    spearman_p_adj = p.adjust(.data$spearman_p_value, method = "BH"),
    sample_centered_p_adj = p.adjust(.data$sample_centered_p_value, method = "BH"),
    sample_fixed_p_adj = p.adjust(.data$sample_fixed_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(.data$cna_metric, .data$spearman_p_adj)

write_table(pairwise_tests, "Auto_pairwise_cna_expression_distance_tests.csv")

message("Consensus clustering subclone arm-level CNA profiles")

choose_consensus_k <- function(mat, max_k = 6L) {
  if (nrow(mat) < 4) return(1L)
  max_k <- min(max_k, nrow(mat) - 1L)
  d <- dist(mat)
  scores <- lapply(2:max_k, function(k) {
    cl <- cutree(hclust(d, method = "ward.D2"), k = k)
    if (any(table(cl) < 2)) return(data.frame(k = k, silhouette = NA_real_))
    sil <- tryCatch(mean(cluster::silhouette(cl, d)[, "sil_width"]), error = function(e) NA_real_)
    data.frame(k = k, silhouette = sil)
  }) %>% bind_rows()
  scores <- scores %>% filter(is.finite(.data$silhouette))
  if (nrow(scores) == 0) return(min(3L, nrow(mat)))
  scores$k[which.max(scores$silhouette)]
}

valid_rows <- rownames(arm_matrix)[rowSums(is.finite(arm_matrix)) == ncol(arm_matrix)]
arm_matrix_valid <- arm_matrix[valid_rows, , drop = FALSE]
hc <- hclust(dist(arm_matrix_valid), method = "ward.D2")
best_k <- choose_consensus_k(arm_matrix_valid, max_k = 6L)
cna_cluster <- if (best_k <= 1L) {
  setNames(rep("CNA cluster 1", nrow(arm_matrix_valid)), rownames(arm_matrix_valid))
} else {
  setNames(paste0("CNA cluster ", cutree(hc, k = best_k)), rownames(arm_matrix_valid))
}

cluster_df <- data.frame(
  subclone_id = names(cna_cluster),
  cna_cluster = unname(cna_cluster),
  stringsAsFactors = FALSE
) %>%
  left_join(features %>% select(.data$subclone_id, .data$sample, .data$study, .data$subclone, .data$n_cells,
                                .data$subclone_fraction, .data$is_largest_subclone, .data$is_significant_dominant),
            by = "subclone_id")

features <- features %>% left_join(cluster_df %>% select(.data$subclone_id, .data$cna_cluster), by = "subclone_id")
write_table(cluster_df, "Auto_cna_consensus_clusters.csv")
write_table(features, "Auto_subclone_feature_summary_with_clusters.csv")
save_rds(cna_cluster, "Auto_cna_consensus_clusters.rds")

cluster_feature_tests <- bind_rows(lapply(target_features, function(feat) {
  df <- features %>% filter(!is.na(.data$cna_cluster), is.finite(.data[[feat]]))
  p_value <- if (nrow(df) >= 5 && n_distinct(df$cna_cluster) >= 2) {
    tryCatch(kruskal.test(df[[feat]] ~ df$cna_cluster)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  data.frame(
    feature = feat,
    feature_label = feature_label(feat),
    n_subclones = nrow(df),
    n_clusters = n_distinct(df$cna_cluster),
    kruskal_p_value = p_value,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(kruskal_p_adj = p.adjust(.data$kruskal_p_value, method = "BH")) %>%
  arrange(.data$kruskal_p_adj)

write_table(cluster_feature_tests, "Auto_cna_cluster_feature_tests.csv")

####################
message("Plotting CNA consensus cluster MP/state associations")

cluster_k_diagnostics <- bind_rows(lapply(2:min(12L, nrow(arm_matrix_valid) - 1L), function(k) {
  cl <- cutree(hc, k = k)
  sil <- if (any(table(cl) < 2)) {
    NA_real_
  } else {
    tryCatch(mean(cluster::silhouette(cl, dist(arm_matrix_valid))[, "sil_width"]), error = function(e) NA_real_)
  }
  data.frame(
    k = k,
    min_cluster_size = min(table(cl)),
    max_cluster_size = max(table(cl)),
    mean_silhouette = sil,
    selected_current = k == length(unique(cna_cluster)),
    stringsAsFactors = FALSE
  )
}))
write_table(cluster_k_diagnostics, "Auto_cna_consensus_k_silhouette_diagnostics.csv")

cluster_size_summary <- features %>%
  count(.data$cna_cluster, name = "n_subclones") %>%
  arrange(.data$cna_cluster)
write_table(cluster_size_summary, "Auto_cna_consensus_cluster_sizes.csv")

cluster_plot_state_features <- setdiff(state_feature_cols, c(hybrid_state_col, unresolved_state_col))
cluster_plot_features <- intersect(c(mp_feature_cols, cluster_plot_state_features), colnames(features))
cluster_plot_group <- c(
  setNames(rep("Metaprogrammes", length(intersect(mp_feature_cols, cluster_plot_features))), intersect(mp_feature_cols, cluster_plot_features)),
  setNames(rep("Centred states", length(intersect(cluster_plot_state_features, cluster_plot_features))), intersect(cluster_plot_state_features, cluster_plot_features))
)
cluster_plot_features <- names(cluster_plot_group)

cluster_p_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

cluster_plot_long <- features %>%
  select(.data$sample, .data$subclone, .data$subclone_id, .data$n_cells,
         .data$subclone_fraction, .data$cna_cluster, all_of(cluster_plot_features)) %>%
  pivot_longer(all_of(cluster_plot_features), names_to = "feature", values_to = "value") %>%
  filter(!is.na(.data$cna_cluster), is.finite(.data$value)) %>%
  mutate(
    cna_cluster = factor(.data$cna_cluster, levels = sort(unique(as.character(features$cna_cluster)))),
    feature_group = unname(cluster_plot_group[.data$feature]),
    feature_label = feature_label(.data$feature)
  )

cluster_mp_state_tests <- cluster_plot_long %>%
  group_by(.data$feature_group, .data$feature, .data$feature_label) %>%
  summarise(
    n_subclones = n(),
    n_clusters = n_distinct(.data$cna_cluster),
    kruskal_p_value = if (n() >= 5 && n_distinct(.data$cna_cluster) >= 2) {
      tryCatch(kruskal.test(.data$value ~ .data$cna_cluster)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(.data$feature_group) %>%
  mutate(kruskal_p_adj_group = p.adjust(.data$kruskal_p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    kruskal_p_adj_global = p.adjust(.data$kruskal_p_value, method = "BH"),
    sig_label = cluster_p_to_stars(.data$kruskal_p_adj_group)
  ) %>%
  arrange(.data$feature_group, .data$kruskal_p_adj_group)
write_table(cluster_mp_state_tests, "Auto_cna_consensus_cluster_mp_state_tests.csv")

cluster_pairwise_tests <- bind_rows(lapply(split(cluster_plot_long, cluster_plot_long$feature), function(df) {
  if (n_distinct(df$cna_cluster) < 2) return(NULL)
  pw <- tryCatch(pairwise.wilcox.test(df$value, df$cna_cluster, p.adjust.method = "BH")$p.value,
                 error = function(e) NULL)
  if (is.null(pw)) return(NULL)
  as.data.frame(as.table(pw), stringsAsFactors = FALSE) %>%
    filter(!is.na(.data$Freq)) %>%
    transmute(
      feature = df$feature[1],
      feature_label = df$feature_label[1],
      feature_group = df$feature_group[1],
      cluster_1 = .data$Var1,
      cluster_2 = .data$Var2,
      pairwise_p_adj = .data$Freq
    )
})) %>%
  arrange(.data$feature_group, .data$pairwise_p_adj)
write_table(cluster_pairwise_tests, "Auto_cna_consensus_cluster_mp_state_pairwise_wilcox.csv")

cluster_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(cluster_plot_long$cna_cluster))),
  sort(unique(as.character(cluster_plot_long$cna_cluster)))
)

plot_cluster_feature_boxplots <- function(group_name, title_text) {
  df <- cluster_plot_long %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(.data$feature_label, levels = unique(feature_label(cluster_plot_features[cluster_plot_group == group_name]))))
  sig_df <- cluster_mp_state_tests %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(.data$feature_label, levels = unique(feature_label(cluster_plot_features[cluster_plot_group == group_name]))))
  
  pw_df <- cluster_pairwise_tests %>%
    filter(.data$feature_group == group_name, .data$pairwise_p_adj < 0.05) %>%
    mutate(feature_label = factor(.data$feature_label, levels = levels(df$feature_label))) %>%
    mutate(p.signif = cluster_p_to_stars(.data$pairwise_p_adj)) %>%
    rename(group1 = .data$cluster_1, group2 = .data$cluster_2)
  
  if (nrow(pw_df) > 0) {
    max_vals <- df %>% 
      group_by(.data$feature_label) %>% 
      summarise(max_val = max(.data$value, na.rm = TRUE), diff = max(.data$value, na.rm=TRUE) - min(.data$value, na.rm=TRUE), .groups="drop")
    pw_df <- pw_df %>%
      left_join(max_vals, by = "feature_label") %>%
      group_by(.data$feature_label) %>%
      mutate(y.position = .data$max_val + .data$diff * 0.1 * row_number()) %>%
      ungroup()
  }
  
  p <- ggplot(df, aes(.data$cna_cluster, .data$value, fill = .data$cna_cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.72, linewidth = 0.45) +
    geom_point(aes(size = .data$n_cells), position = position_jitter(width = 0.16, height = 0),
               alpha = 0.55, shape = 21, color = "black") +
    geom_text(data = sig_df, aes(x = Inf, y = Inf, label = .data$sig_label),
              inherit.aes = FALSE, hjust = 1.2, vjust = 1.2, size = 8, fontface = "bold")
              
  if (nrow(pw_df) > 0) {
    p <- p + stat_pvalue_manual(pw_df, label = "p.signif", size = 6, hide.ns = TRUE)
  }
  
  p + facet_wrap(~feature_label, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = cluster_cols, drop = FALSE) +
    scale_size_continuous(range = c(1.2, 4.5)) +
    labs(title = title_text, x = NULL, y = "Subclone feature value",
         fill = "CNA cluster", size = "Cells") +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 26),
      strip.text = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 16)
    )
}

cluster_k_plot <- ggplot(cluster_k_diagnostics, aes(.data$k, .data$mean_silhouette)) +
  geom_line(linewidth = 0.8, color = "grey35") +
  geom_point(aes(fill = .data$selected_current), shape = 21, size = 4, color = "black") +
  geom_text(aes(label = ifelse(.data$selected_current, "current", "")),
            vjust = -0.9, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#B2182B"), guide = "none") +
  scale_x_continuous(breaks = cluster_k_diagnostics$k) +
  labs(title = "Consensus CNA cluster k diagnostics",
       x = "Number of clusters", y = "Mean silhouette width") +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 22),
        axis.text = element_text(size = 15))
ggsave(file.path(figure_dir, "Auto_cna_consensus_k_silhouette_diagnostics.pdf"),
       cluster_k_plot, width = 9, height = 6, useDingbats = FALSE)

cluster_mp_plot <- plot_cluster_feature_boxplots(
  "Metaprogrammes",
  "CNA consensus clusters versus metaprogramme expression"
)
cluster_state_plot <- plot_cluster_feature_boxplots(
  "Centred states",
  "CNA consensus clusters versus state abundance"
)
ggsave(file.path(figure_dir, "Auto_cna_consensus_cluster_mp_boxplots.pdf"),
       cluster_mp_plot, width = 18, height = 12, useDingbats = FALSE)
ggsave(file.path(figure_dir, "Auto_cna_consensus_cluster_state_boxplots.pdf"),
       cluster_state_plot, width = 16, height = 9, useDingbats = FALSE)
pdf(file.path(figure_dir, "Auto_cna_consensus_cluster_mp_state_boxplots.pdf"),
    width = 18, height = 12, useDingbats = FALSE)
print(cluster_mp_plot)
print(cluster_state_plot)
dev.off()

dir.create("../updates/new_updates/summaries", recursive = TRUE, showWarnings = FALSE)
cluster_compact_summary <- bind_rows(
  cluster_k_diagnostics %>%
    transmute(summary_type = "k_diagnostic",
              item = paste0("k_", .data$k),
              metric = "mean_silhouette",
              value = signif(.data$mean_silhouette, 5),
              detail = paste0("min_cluster_size=", .data$min_cluster_size,
                              ";current=", .data$selected_current)),
  cluster_mp_state_tests %>%
    head(12) %>%
    transmute(summary_type = "top_cluster_feature",
              item = .data$feature,
              metric = "kruskal_group_fdr",
              value = signif(.data$kruskal_p_adj_group, 5),
              detail = paste(.data$feature_group, .data$feature_label, .data$sig_label, sep = " | "))
)
write.csv(cluster_compact_summary,
          "../updates/new_updates/summaries/Auto_cna_consensus_cluster_mp_state_summary.csv",
          row.names = FALSE)
####################

message("Summarising recurrent chromosome-arm CNA events")

event_summary <- bind_rows(
  arm_long %>%
    filter(.data$cna_call == 1L) %>%
    mutate(direction = "gain"),
  arm_long %>%
    filter(.data$cna_call == -1L) %>%
    mutate(direction = "loss")
) %>%
  mutate(event_id = paste(.data$direction, as.character(.data$arm), sep = "_")) %>%
  group_by(.data$direction, .data$arm, .data$event_id) %>%
  summarise(
    n_subclones_event = n_distinct(.data$subclone_id),
    n_samples_event = n_distinct(.data$sample),
    median_arm_mean = median(.data$arm_mean, na.rm = TRUE),
    mean_abs_arm_mean = mean(abs(.data$arm_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_subclones_total = nrow(features),
    n_samples_total = n_distinct(features$sample),
    frac_subclones_event = .data$n_subclones_event / .data$n_subclones_total,
    frac_samples_event = .data$n_samples_event / .data$n_samples_total,
    is_recurrent = .data$n_samples_event >= pmax(recurrent_min_samples, ceiling(recurrent_min_sample_fraction * .data$n_samples_total))
  ) %>%
  arrange(desc(.data$frac_samples_event), desc(.data$frac_subclones_event))

write_table(event_summary, "Auto_recurrent_cna_events.csv")

parse_oac_cnv <- function(path) {
  raw <- readxl::read_excel(path, sheet = 1, col_names = FALSE)
  rows <- list()
  current <- NA_character_
  row_i <- 1L
  for (i in seq_len(nrow(raw))) {
    first <- as.character(raw[[1]][i])
    if (is.na(first)) next
    if (grepl("^CNV gain", first, ignore.case = TRUE)) {
      current <- "gain"
      next
    }
    if (grepl("^CNV loss", first, ignore.case = TRUE)) {
      current <- "loss"
      next
    }
    rank_num <- suppressWarnings(as.integer(first))
    if (!is.na(rank_num) && !is.na(current)) {
      genes <- as.character(raw[[4]][i])
      rows[[row_i]] <- data.frame(
        source = "OAC_CNV_curated",
        direction = current,
        rank = rank_num,
        cytoband = as.character(raw[[2]][i]),
        region = as.character(raw[[3]][i]),
        genes = genes,
        frequency = as.character(raw[[5]][i]),
        pathway = as.character(raw[[6]][i]),
        clinical_relevance = as.character(raw[[7]][i]),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  bind_rows(rows) %>%
    mutate(
      arm = cytoband_to_arm(.data$cytoband),
      gene = lapply(.data$genes, split_gene_string)
    ) %>%
    tidyr::unnest(.data$gene, keep_empty = TRUE)
}

parse_occams_gistic_peaks <- function(path, sheet, direction) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE)
  first_col <- as.character(raw[[1]])
  cytoband_row <- which(tolower(first_col) == "cytoband")[1]
  q_row <- which(tolower(first_col) == "q value")[1]
  residual_q_row <- which(tolower(first_col) == "residual q value")[1]
  boundary_row <- which(tolower(first_col) == "wide peak boundaries")[1]
  gene_row <- which(tolower(first_col) == "genes in wide peak")[1]
  if (is.na(cytoband_row) || is.na(boundary_row) || is.na(gene_row)) return(data.frame())
  rows <- list()
  row_i <- 1L
  for (j in 2:ncol(raw)) {
    cytoband <- as.character(raw[[j]][cytoband_row])
    if (is.na(cytoband) || !nzchar(cytoband)) next
    genes <- as.character(raw[[j]][(gene_row + 1):nrow(raw)])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    coord <- parse_peak_coord(raw[[j]][boundary_row])
    rows[[row_i]] <- data.frame(
      source = "OCCAMS_GISTIC_peak",
      direction = direction,
      cytoband = cytoband,
      q_value = suppressWarnings(as.numeric(raw[[j]][q_row])),
      residual_q_value = suppressWarnings(as.numeric(raw[[j]][residual_q_row])),
      wide_peak_boundaries = as.character(raw[[j]][boundary_row]),
      chr = coord$chr,
      start = coord$start,
      end = coord$end,
      arm = ifelse(!is.na(coord$arm), coord$arm, cytoband_to_arm(cytoband)),
      genes = paste(genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
  bind_rows(rows)
}

parse_occams_driver_sheet <- function(path, sheet, direction, confidence) {
  x <- suppressMessages(readxl::read_excel(path, sheet = sheet, skip = 6))
  if (!"hgnc_symbol" %in% colnames(x)) return(data.frame())
  peak_col <- grep("GISTIC Peak", colnames(x), value = TRUE)[1]
  if (is.na(peak_col)) return(data.frame())
  out <- x %>%
    filter(!is.na(.data$hgnc_symbol), nzchar(as.character(.data$hgnc_symbol))) %>%
    mutate(
      source = paste0("OCCAMS_", confidence, "_driver"),
      direction = direction,
      gene = toupper(as.character(.data$hgnc_symbol)),
      peak = as.character(.data[[peak_col]]),
      cytoband = sub(" .*", "", .data$peak)
    )
  coords <- bind_rows(lapply(out$peak, parse_peak_coord))
  out <- bind_cols(out, coords) %>%
    mutate(arm = ifelse(!is.na(.data$arm), .data$arm, cytoband_to_arm(.data$cytoband))) %>%
    transmute(
      source = as.character(.data$source),
      direction = as.character(.data$direction),
      gene = as.character(.data$gene),
      cytoband = as.character(.data$cytoband),
      arm = as.character(.data$arm),
      peak = as.character(.data$peak)
    )
  as.data.frame(out)
}

message("Parsing local OAC and OCCAMS CNA annotation files")
oac_cnv <- parse_oac_cnv(oac_cnv_path)
occams_amp_peaks <- parse_occams_gistic_peaks(occams_path, "ST1 GISTIC amplification peaks", "gain")
occams_del_peaks <- parse_occams_gistic_peaks(occams_path, "ST2 GISTIC Deletion peaks", "loss")
occams_driver_rows <- bind_rows(
  parse_occams_driver_sheet(occams_path, "ST3 High Confidence Del Drivers", "loss", "high_confidence_deletion"),
  parse_occams_driver_sheet(occams_path, "ST4 Candidate Del Drivers", "loss", "candidate_deletion"),
  parse_occams_driver_sheet(occams_path, "ST5 High confidence Amp Drivers", "gain", "high_confidence_amplification"),
  parse_occams_driver_sheet(occams_path, "ST6 Candidate Amp Drivers", "gain", "candidate_amplification")
)

write_table(oac_cnv, "Auto_oac_cnv_curated.csv")
write_table(bind_rows(occams_amp_peaks, occams_del_peaks), "Auto_occams_gistic_peaks.csv")
write_table(occams_driver_rows, "Auto_occams_cna_driver_genes.csv")

fetch_cbioportal_tcga_cna <- function(symbols) {
  symbols <- unique(toupper(symbols[!is.na(symbols) & nzchar(symbols)]))
  if (length(symbols) == 0) return(data.frame())
  base_url <- "https://www.cbioportal.org/api"
  study_id <- "esca_tcga_pan_can_atlas_2018"
  profile_id <- "esca_tcga_pan_can_atlas_2018_gistic"
  sample_list_id <- "esca_tcga_pan_can_atlas_2018_cna"
  tryCatch({
    sample_resp <- httr::GET(file.path(base_url, "sample-lists", sample_list_id, "sample-ids"))
    httr::stop_for_status(sample_resp)
    sample_ids <- jsonlite::fromJSON(httr::content(sample_resp, as = "text", encoding = "UTF-8"))
    sample_count <- length(sample_ids)

    gene_resp <- httr::POST(
      paste0(base_url, "/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY"),
      body = jsonlite::toJSON(symbols, auto_unbox = TRUE),
      httr::content_type_json()
    )
    httr::stop_for_status(gene_resp)
    gene_df <- jsonlite::fromJSON(httr::content(gene_resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
    if (nrow(gene_df) == 0) return(data.frame())

    entrez_ids <- unique(gene_df$entrezGeneId)
    chunks <- split(entrez_ids, ceiling(seq_along(entrez_ids) / 150))
    cna_rows <- lapply(chunks, function(ids) {
      body <- list(entrezGeneIds = as.integer(ids), sampleListId = sample_list_id)
      resp <- httr::POST(
        paste0(base_url, "/molecular-profiles/", profile_id, "/discrete-copy-number/fetch?projection=SUMMARY"),
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        httr::content_type_json()
      )
      httr::stop_for_status(resp)
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      if (!nzchar(txt) || identical(txt, "[]")) return(data.frame())
      jsonlite::fromJSON(txt, flatten = TRUE)
    })
    cna_df <- bind_rows(cna_rows)
    if (nrow(cna_df) == 0) {
      return(gene_df %>%
               transmute(
                 source = "TCGA_cBioPortal",
                 study_id = study_id,
                 molecular_profile_id = profile_id,
                 sample_list_id = sample_list_id,
                 sample_count = sample_count,
                 gene = toupper(.data$hugoGeneSymbol),
                 entrezGeneId = .data$entrezGeneId,
                 n_gain = 0L,
                 n_amp = 0L,
                 n_any_gain = 0L,
                 n_loss = 0L,
                 n_deep_del = 0L,
                 n_any_loss = 0L,
                 pct_any_gain = 0,
                 pct_amp = 0,
                 pct_any_loss = 0,
                 pct_deep_del = 0
               ))
    }
    cna_df %>%
      left_join(gene_df %>% select(.data$entrezGeneId, gene = .data$hugoGeneSymbol), by = "entrezGeneId") %>%
      mutate(gene = toupper(.data$gene)) %>%
      group_by(.data$gene, .data$entrezGeneId) %>%
      summarise(
        n_gain = sum(.data$alteration == 1L, na.rm = TRUE),
        n_amp = sum(.data$alteration == 2L, na.rm = TRUE),
        n_any_gain = sum(.data$alteration %in% c(1L, 2L), na.rm = TRUE),
        n_loss = sum(.data$alteration == -1L, na.rm = TRUE),
        n_deep_del = sum(.data$alteration == -2L, na.rm = TRUE),
        n_any_loss = sum(.data$alteration %in% c(-1L, -2L), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        source = "TCGA_cBioPortal",
        study_id = study_id,
        molecular_profile_id = profile_id,
        sample_list_id = sample_list_id,
        sample_count = sample_count,
        pct_any_gain = .data$n_any_gain / sample_count,
        pct_amp = .data$n_amp / sample_count,
        pct_any_loss = .data$n_any_loss / sample_count,
        pct_deep_del = .data$n_deep_del / sample_count
      ) %>%
      select(.data$source, .data$study_id, .data$molecular_profile_id, .data$sample_list_id,
             .data$sample_count, .data$gene, .data$entrezGeneId, everything())
  }, error = function(e) {
    warning("cBioPortal fetch failed; continuing with local annotations only: ", conditionMessage(e))
    data.frame()
  })
}

annotation_symbols <- unique(c(
  oac_cnv$gene,
  occams_driver_rows$gene,
  "MYC", "ERBB2", "EGFR", "CCNE1", "CCND1", "CDK6", "MDM2", "MCL1",
  "GATA4", "GATA6", "KRAS", "MET", "PIK3CA", "VEGFA", "CDKN2A",
  "MTAP", "TP53", "SMAD4", "APC", "PTEN", "RB1", "ATM", "BRCA1",
  "BRCA2", "ARID1A", "FHIT", "CSMD1"
))
tcga_cna <- fetch_cbioportal_tcga_cna(annotation_symbols)
write_table(tcga_cna, "Auto_tcga_cbioportal_cna_gene_frequencies.csv")

tcga_annotation <- data.frame()
if (nrow(tcga_cna) > 0) {
  tcga_annotation <- tcga_cna %>%
    left_join(gene_to_arm, by = c("gene" = "gene")) %>%
    select(.data$source, .data$gene, .data$chr_arm, .data$pct_any_gain, .data$pct_amp,
           .data$pct_any_loss, .data$pct_deep_del) %>%
    pivot_longer(cols = c(.data$pct_any_gain, .data$pct_any_loss),
                 names_to = "tcga_event", values_to = "frequency") %>%
    mutate(
      direction = ifelse(.data$tcga_event == "pct_any_gain", "gain", "loss"),
      arm = .data$chr_arm,
      cytoband = NA_character_,
      genes = .data$gene,
      frequency_label = percent(.data$frequency, accuracy = 0.1),
      clinical_relevance = NA_character_,
      pathway = NA_character_
    ) %>%
    filter(!is.na(.data$arm), .data$frequency > 0) %>%
    select(.data$source, .data$direction, .data$arm, .data$cytoband, .data$gene,
           .data$genes, .data$frequency_label, .data$pathway, .data$clinical_relevance)
}

oac_annotation <- oac_cnv %>%
  transmute(
    source = .data$source,
    direction = .data$direction,
    arm = .data$arm,
    cytoband = .data$cytoband,
    gene = .data$gene,
    genes = .data$genes,
    frequency_label = .data$frequency,
    pathway = .data$pathway,
    clinical_relevance = .data$clinical_relevance
  )

occams_driver_annotation <- occams_driver_rows %>%
  transmute(
    source = .data$source,
    direction = .data$direction,
    arm = .data$arm,
    cytoband = .data$cytoband,
    gene = .data$gene,
    genes = .data$gene,
    frequency_label = NA_character_,
    pathway = NA_character_,
    clinical_relevance = NA_character_
  )

cnv_annotation_long <- bind_rows(oac_annotation, occams_driver_annotation, tcga_annotation) %>%
  filter(!is.na(.data$arm), .data$arm %in% arm_levels) %>%
  mutate(event_id = paste(.data$direction, .data$arm, sep = "_"))

event_annotation <- cnv_annotation_long %>%
  group_by(.data$event_id, .data$direction, .data$arm) %>%
  summarise(
    annotation_sources = paste(sort(unique(.data$source)), collapse = "; "),
    known_genes = paste(head(sort(unique(.data$gene[!is.na(.data$gene) & nzchar(.data$gene)])), 18), collapse = ", "),
    cytobands = paste(sort(unique(.data$cytoband[!is.na(.data$cytoband) & nzchar(.data$cytoband)])), collapse = ", "),
    pathways = paste(sort(unique(.data$pathway[!is.na(.data$pathway) & nzchar(.data$pathway)])), collapse = "; "),
    clinical_relevance = paste(sort(unique(.data$clinical_relevance[!is.na(.data$clinical_relevance) & nzchar(.data$clinical_relevance)])), collapse = "; "),
    .groups = "drop"
  )

event_summary_annotated <- event_summary %>%
  left_join(event_annotation, by = c("event_id", "direction", "arm")) %>%
  mutate(
    known_genes = ifelse(is.na(.data$known_genes), "", .data$known_genes),
    annotation_sources = ifelse(is.na(.data$annotation_sources), "", .data$annotation_sources)
  )

write_table(cnv_annotation_long, "Auto_cnv_annotation_long.csv")
write_table(event_annotation, "Auto_cnv_event_annotation_summary.csv")
write_table(event_summary_annotated, "Auto_recurrent_cna_events_annotated.csv")

message("Testing recurrent CNA event feature associations")

recurrent_events <- event_summary_annotated %>%
  filter(.data$is_recurrent) %>%
  arrange(desc(.data$frac_samples_event), desc(.data$frac_subclones_event)) %>%
  pull(.data$event_id)

if (length(recurrent_events) < 1) {
  recurrent_events <- event_summary_annotated %>%
    arrange(desc(.data$frac_samples_event), desc(.data$frac_subclones_event)) %>%
    head(12) %>%
    pull(.data$event_id)
}
recurrent_events <- head(recurrent_events, 30)

event_call_long <- arm_long %>%
  mutate(
    arm = as.character(.data$arm),
    gain_event = paste("gain", .data$arm, sep = "_"),
    loss_event = paste("loss", .data$arm, sep = "_")
  ) %>%
  select(.data$subclone_id, .data$sample, .data$subclone, .data$arm, .data$cna_call) %>%
  tidyr::crossing(direction = c("gain", "loss")) %>%
  mutate(
    event_id = paste(.data$direction, .data$arm, sep = "_"),
    event_present = case_when(
      .data$direction == "gain" ~ .data$cna_call == 1L,
      .data$direction == "loss" ~ .data$cna_call == -1L,
      TRUE ~ FALSE
    )
  ) %>%
  filter(.data$event_id %in% recurrent_events)

event_feature_tests <- bind_rows(lapply(recurrent_events, function(event_id) {
  event_df <- event_call_long %>%
    filter(.data$event_id == event_id) %>%
    select(.data$subclone_id, .data$event_present) %>%
    distinct() %>%
    left_join(features, by = "subclone_id")
  bind_rows(lapply(target_features, function(feat) {
    df <- event_df %>% filter(is.finite(.data[[feat]]), !is.na(.data$event_present))
    unpaired <- if (nrow(df) >= 5 && length(unique(df$event_present)) == 2) {
      tryCatch(wilcox.test(df[[feat]][df$event_present], df[[feat]][!df$event_present])$p.value,
               error = function(e) NA_real_)
    } else {
      NA_real_
    }
    paired_delta <- bind_rows(lapply(split(df, df$sample), function(sdf) {
      if (!all(c(TRUE, FALSE) %in% sdf$event_present)) return(NULL)
      data.frame(
        sample = sdf$sample[1],
        delta = weighted.mean(sdf[[feat]][sdf$event_present], sdf$n_cells[sdf$event_present], na.rm = TRUE) -
          weighted.mean(sdf[[feat]][!sdf$event_present], sdf$n_cells[!sdf$event_present], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    paired_p <- if (nrow(paired_delta) >= 3 && sum(is.finite(paired_delta$delta)) >= 3) {
      tryCatch(wilcox.test(paired_delta$delta, mu = 0)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    data.frame(
      event_id = event_id,
      feature = feat,
      feature_label = feature_label(feat),
      n_subclones_event = sum(df$event_present, na.rm = TRUE),
      n_subclones_no_event = sum(!df$event_present, na.rm = TRUE),
      mean_event = safe_mean(df[[feat]][df$event_present]),
      mean_no_event = safe_mean(df[[feat]][!df$event_present]),
      unpaired_delta = safe_mean(df[[feat]][df$event_present]) - safe_mean(df[[feat]][!df$event_present]),
      unpaired_p_value = unpaired,
      n_paired_samples = nrow(paired_delta),
      paired_median_delta = if (nrow(paired_delta) > 0) median(paired_delta$delta, na.rm = TRUE) else NA_real_,
      paired_p_value = paired_p,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  left_join(event_summary_annotated %>% select(.data$event_id, .data$direction, .data$arm,
                                               .data$frac_samples_event, .data$known_genes,
                                               .data$annotation_sources),
            by = "event_id") %>%
  group_by(.data$feature) %>%
  mutate(
    unpaired_p_adj_by_feature = p.adjust(.data$unpaired_p_value, method = "BH"),
    paired_p_adj_by_feature = p.adjust(.data$paired_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    unpaired_p_adj_global = p.adjust(.data$unpaired_p_value, method = "BH"),
    paired_p_adj_global = p.adjust(.data$paired_p_value, method = "BH")
  ) %>%
  arrange(.data$paired_p_adj_global, .data$unpaired_p_adj_global)

  write_table(event_feature_tests, "Auto_recurrent_cna_event_feature_tests.csv")
  save_rds(
    list(
      features = features,
      arm_long = arm_long,
      arm_matrix = arm_matrix,
      arm_call_matrix = arm_call_matrix,
      dominance = dominance_df,
      pairwise = pairwise_df,
      event_summary = event_summary_annotated,
      event_feature_tests = event_feature_tests,
      cna_cluster = cna_cluster
    ),
    "Auto_cna_subclone_expression_results.rds"
  )
} # end of else (!replot_only)

event_summary <- event_summary_annotated

message("Creating plots and publication-quality figures")

event_call_threshold <- 0.05
legacy_event_call_threshold <- 0.10
max_plot_cells <- 1200L
cna_colour_limit <- 0.15


mp_features <- paste0("mp__", mp_order)
mp_features <- mp_features[mp_features %in% colnames(features)]

state_features <- unname(state_display[state_display %in% colnames(features)])

qc_features <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score",
    "cna_signal", "cna_cor", "cna_burden_mean_abs", "n_gain_arms", "n_loss_arms"),
  colnames(features)
)

feature_group <- c(
  setNames(rep("Metaprogrammes", length(mp_features)), mp_features),
  setNames(rep("Centred states", length(state_features)), state_features),
  setNames(rep("QC / CNA metrics", length(qc_features)), qc_features)
)
plot_features <- names(feature_group)

feature_label <- function(feature) {
  out <- as.character(feature)
  is_mp <- grepl("^mp__", out)
  mp <- sub("^mp__", "", out[is_mp])
  out[is_mp] <- paste0(mp, " ", mp_descriptions[mp])
  state_lookup <- setNames(names(state_display), unname(state_display))
  is_state <- out %in% names(state_lookup)
  out[is_state] <- state_lookup[out[is_state]]
  out <- gsub("^nCount_RNA$", "nCount_RNA", out)
  out <- gsub("^nFeature_RNA$", "nFeature_RNA", out)
  out <- gsub("^percent.mt$", "percent.mt", out)
  out <- gsub("^cc_score$", "CC score", out)
  out <- gsub("^cs_score$", "CS score", out)
  out <- gsub("^cna_signal$", "CNA signal", out)
  out <- gsub("^cna_cor$", "CNA correlation", out)
  out <- gsub("^cna_burden_mean_abs$", "CNA burden", out)
  out <- gsub("^n_gain_arms$", "No. gained arms", out)
  out <- gsub("^n_loss_arms$", "No. lost arms", out)
  out
}

features <- features %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = as.character(.data$subclone_id)
  )

arm_long <- arm_long %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = as.character(.data$subclone_id),
    arm = as.character(.data$arm),
    cna_call_legacy = .data$cna_call,
    cna_call = case_when(
      .data$arm_mean >= event_call_threshold ~ 1L,
      .data$arm_mean <= -event_call_threshold ~ -1L,
      TRUE ~ 0L
    )
  )

event_known_genes <- event_summary %>%
  select(.data$event_id, .data$known_genes) %>%
  distinct()

summarise_events <- function(direction_label, call_value) {
  arm_long %>%
    filter(.data$cna_call == call_value) %>%
    group_by(.data$arm) %>%
    summarise(
      direction = direction_label,
      n_samples_event = n_distinct(.data$sample),
      n_subclones_event = n_distinct(.data$subclone_id),
      samples_event = paste(sort(unique(.data$sample)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(event_id = paste0(.data$direction, "_", .data$arm))
}

event_summary_recomputed <- bind_rows(
  summarise_events("gain", 1L),
  summarise_events("loss", -1L)
) %>%
  mutate(
    n_samples_total = n_distinct(features$sample),
    n_subclones_total = nrow(features),
    frac_samples_event = .data$n_samples_event / .data$n_samples_total,
    frac_subclones_event = .data$n_subclones_event / .data$n_subclones_total,
    is_recurrent = .data$n_samples_event >= recurrent_min_samples &
      .data$frac_samples_event >= recurrent_min_sample_fraction
  ) %>%
  left_join(event_known_genes, by = "event_id") %>%
  mutate(known_genes = ifelse(is.na(.data$known_genes), "", .data$known_genes)) %>%
  arrange(desc(.data$is_recurrent), desc(.data$frac_samples_event), desc(.data$frac_subclones_event))

threshold_sensitivity <- bind_rows(lapply(c(0, 0.02, 0.03, 0.05, 0.07, 0.08, 0.10, 0.12), function(thr) {
  bind_rows(
    arm_long %>%
      group_by(.data$arm) %>%
      summarise(
        threshold = thr,
        direction = "gain",
        n_subclones_event = sum(.data$arm_mean >= .env$thr, na.rm = TRUE),
        n_samples_event = n_distinct(.data$sample[.data$arm_mean >= .env$thr]),
        .groups = "drop"
      ),
    arm_long %>%
      group_by(.data$arm) %>%
      summarise(
        threshold = thr,
        direction = "loss",
        n_subclones_event = sum(.data$arm_mean <= -.env$thr, na.rm = TRUE),
        n_samples_event = n_distinct(.data$sample[.data$arm_mean <= -.env$thr]),
        .groups = "drop"
      )
  ) %>%
    mutate(event_id = paste0(.data$direction, "_", .data$arm))
})) %>%
  mutate(
    n_samples_total = n_distinct(features$sample),
    n_subclones_total = nrow(features),
    frac_samples_event = .data$n_samples_event / .data$n_samples_total,
    frac_subclones_event = .data$n_subclones_event / .data$n_subclones_total
  )

write_table(threshold_sensitivity, "Auto_v2_cna_event_threshold_sensitivity.csv")
write_table(event_summary_recomputed, "Auto_v2_recomputed_recurrent_cna_event_summary.csv")

ranked_events <- event_summary_recomputed %>%
  arrange(desc(.data$is_recurrent), desc(.data$frac_samples_event), desc(.data$frac_subclones_event))
n_events_to_plot <- min(nrow(ranked_events), 8L)
recurrent_events <- ranked_events %>%
  head(n_events_to_plot) %>%
  pull(.data$event_id)
boxplot_events <- head(recurrent_events, min(4L, length(recurrent_events)))

event_meta <- event_summary_recomputed %>%
  filter(.data$event_id %in% recurrent_events) %>%
  mutate(
    event_label = event_label(.data$event_id, .data$known_genes),
    event_label = factor(.data$event_label, levels = event_label(.data$event_id, .data$known_genes))
  )

event_presence_for <- function(event_id) {
  event_direction <- sub("_.*$", "", event_id)
  event_arm <- sub("^[^_]+_", "", event_id)
  arm_long %>%
    filter(.data$arm == .env$event_arm) %>%
    transmute(
      event_id = .env$event_id,
      sample = .data$sample,
      subclone = .data$subclone,
      subclone_id = .data$subclone_id,
      event_present = if (.env$event_direction == "gain") .data$cna_call == 1L else .data$cna_call == -1L
    ) %>%
    distinct(.data$event_id, .data$sample, .data$subclone, .data$subclone_id, .keep_all = TRUE)
}

event_presence <- bind_rows(lapply(recurrent_events, event_presence_for)) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

write_table(event_presence, "Auto_v2_recurrent_cna_event_subclone_presence.csv")

feature_sd <- vapply(plot_features, function(feat) safe_sd(features[[feat]]), numeric(1))
feature_sd[!is.finite(feature_sd) | feature_sd == 0] <- NA_real_
feature_mean <- vapply(plot_features, function(feat) safe_mean(features[[feat]]), numeric(1))

message("Computing corrected recurrent CNA event feature tests")

event_feature_values <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_df <- event_presence %>%
    filter(.data$event_id == ev) %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    value <- suppressWarnings(as.numeric(ev_df[[feat]]))
    data.frame(
      event_id = ev,
      sample = ev_df$sample,
      subclone = ev_df$subclone,
      subclone_id = ev_df$subclone_id,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      event_present = ev_df$event_present,
      feature_value = value,
      feature_z = (value - feature_mean[[feat]]) / feature_sd[[feat]],
      n_cells = ev_df$n_cells,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

event_sample_deltas <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_df <- event_presence %>%
    filter(.data$event_id == ev) %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    bind_rows(lapply(split(ev_df, ev_df$sample), function(sdf) {
      if (!all(c(TRUE, FALSE) %in% sdf$event_present)) return(NULL)
      event_value <- safe_weighted_mean(sdf[[feat]][sdf$event_present], sdf$n_cells[sdf$event_present])
      no_event_value <- safe_weighted_mean(sdf[[feat]][!sdf$event_present], sdf$n_cells[!sdf$event_present])
      data.frame(
        event_id = ev,
        sample = sdf$sample[1],
        feature = feat,
        feature_label = feature_label(feat),
        feature_group = feature_group[[feat]],
        event_value = event_value,
        no_event_value = no_event_value,
        delta = event_value - no_event_value,
        std_delta = (event_value - no_event_value) / feature_sd[[feat]],
        n_event_subclones = sum(sdf$event_present, na.rm = TRUE),
        n_no_event_subclones = sum(!sdf$event_present, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

event_feature_tests_v2 <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_presence <- event_presence %>% filter(.data$event_id == ev)
  ev_df <- ev_presence %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    d <- event_sample_deltas %>% filter(.data$event_id == ev, .data$feature == feat)
    event_values <- suppressWarnings(as.numeric(ev_df[[feat]][ev_df$event_present]))
    no_event_values <- suppressWarnings(as.numeric(ev_df[[feat]][!ev_df$event_present]))
    unpaired_p <- if (n_distinct(ev_df$event_present) == 2) {
      tryCatch(wilcox.test(event_values, no_event_values)$p.value,
               error = function(e) NA_real_)
    } else {
      NA_real_
    }
    paired_p <- if (nrow(d) >= 3 && sum(is.finite(d$delta)) >= 3) {
      tryCatch(wilcox.test(d$delta, mu = 0)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    data.frame(
      event_id = ev,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      n_subclones_event = sum(ev_df$event_present, na.rm = TRUE),
      n_subclones_no_event = sum(!ev_df$event_present, na.rm = TRUE),
      n_paired_samples = nrow(d),
      mean_event = safe_mean(event_values),
      mean_no_event = safe_mean(no_event_values),
      median_event = safe_median(event_values),
      median_no_event = safe_median(no_event_values),
      unpaired_delta = safe_mean(event_values) - safe_mean(no_event_values),
      unpaired_median_delta = safe_median(event_values) - safe_median(no_event_values),
      unpaired_median_std_delta = (safe_median(event_values) - safe_median(no_event_values)) / feature_sd[[feat]],
      paired_median_delta = if (nrow(d) > 0) median(d$delta, na.rm = TRUE) else NA_real_,
      paired_mean_delta = if (nrow(d) > 0) mean(d$delta, na.rm = TRUE) else NA_real_,
      paired_median_std_delta = if (nrow(d) > 0) median(d$std_delta, na.rm = TRUE) else NA_real_,
      unpaired_p_value = unpaired_p,
      paired_p_value = paired_p,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id") %>%
  group_by(.data$feature_group) %>%
  mutate(
    paired_p_adj_group = p.adjust(.data$paired_p_value, method = "BH"),
    unpaired_p_adj_group = p.adjust(.data$unpaired_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    paired_p_adj_global = p.adjust(.data$paired_p_value, method = "BH"),
    unpaired_p_adj_global = p.adjust(.data$unpaired_p_value, method = "BH"),
    sig_label = ifelse(.data$feature_group == "Metaprogrammes" & abs(.data$unpaired_median_delta) < 1.00, "", p_to_stars(.data$unpaired_p_adj_group)),
    neglog10_fdr = pmin(-log10(pmax(.data$unpaired_p_adj_group, 1e-12)), 12),
    primary_delta = .data$unpaired_median_std_delta,
    primary_p_adj_group = .data$unpaired_p_adj_group
  )

write_table(event_feature_values, "Auto_v2_recurrent_cna_event_feature_values.csv")
write_table(event_sample_deltas, "Auto_v2_recurrent_cna_event_per_sample_feature_deltas.csv")
write_table(event_feature_tests_v2, "Auto_v2_recurrent_cna_event_feature_tests.csv")

message("Computing largest-subclone effects with standardized x-axis")

dominant_deltas_v2 <- bind_rows(lapply(plot_features, function(feat) {
  bind_rows(lapply(split(features, features$sample), function(df) {
    if (nrow(df) < 2) return(NULL)
    top <- df %>% filter(.data$is_largest_subclone)
    rest <- df %>% filter(!.data$is_largest_subclone)
    if (nrow(top) != 1 || nrow(rest) == 0) return(NULL)
    rest_value <- safe_weighted_mean(rest[[feat]], rest$n_cells)
    data.frame(
      sample = top$sample,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      dominant_value = top[[feat]],
      rest_weighted_value = rest_value,
      delta = top[[feat]] - rest_value,
      std_delta = (top[[feat]] - rest_value) / feature_sd[[feat]],
      dominance_class = top$dominance_class,
      stringsAsFactors = FALSE
    )
  }))
}))

dominant_tests_v2 <- dominant_deltas_v2 %>%
  group_by(.data$feature, .data$feature_label, .data$feature_group) %>%
  summarise(
    n_samples = n(),
    median_delta = median(.data$delta, na.rm = TRUE),
    median_std_delta = median(.data$std_delta, na.rm = TRUE),
    pct_positive_delta = mean(.data$delta > 0, na.rm = TRUE),
    wilcox_p_value = if (n() >= 3 && sum(is.finite(.data$delta)) >= 3) {
      tryCatch(wilcox.test(.data$delta, mu = 0)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(.data$feature_group) %>%
  mutate(wilcox_p_adj_group = p.adjust(.data$wilcox_p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    wilcox_p_adj_global = p.adjust(.data$wilcox_p_value, method = "BH"),
    sig_label = ifelse(.data$feature_group == "Metaprogrammes" & abs(.data$median_delta) < 1.00, "", p_to_stars(.data$wilcox_p_adj_group)),
    neglog10_fdr = pmin(-log10(pmax(.data$wilcox_p_adj_group, 1e-12)), 12)
  )

write_table(dominant_deltas_v2, "Auto_v2_largest_subclone_per_sample_feature_deltas.csv")
write_table(dominant_tests_v2, "Auto_v2_largest_subclone_feature_tests.csv")

message("Computing pairwise CNA-distance tests for all feature groups")

cna_distance_metrics <- c("cna_abs_mean", "cna_euclidean", "cna_call_discordance")
cna_metric_labels <- c(
  cna_abs_mean = "Mean absolute CNA distance",
  cna_euclidean = "Euclidean CNA distance",
  cna_call_discordance = "Arm-call discordance"
)

pairwise_feature_tests_v2 <- bind_rows(lapply(cna_distance_metrics, function(cna_metric) {
  bind_rows(lapply(plot_features, function(feat) {
    endpoint <- paste0("abs_delta__", feat)
    if (!endpoint %in% colnames(pairwise_df)) return(NULL)
    df <- pairwise_df %>%
      filter(is.finite(.data[[cna_metric]]), is.finite(.data[[endpoint]]))
    cor_res <- if (nrow(df) >= 5 && safe_sd(df[[cna_metric]]) > 0 && safe_sd(df[[endpoint]]) > 0) {
      tryCatch(cor.test(df[[cna_metric]], df[[endpoint]], method = "spearman"), error = function(e) NULL)
    } else {
      NULL
    }
    centered <- df %>%
      group_by(.data$sample) %>%
      mutate(
        cna_centered = .data[[cna_metric]] - mean(.data[[cna_metric]], na.rm = TRUE),
        endpoint_centered = .data[[endpoint]] - mean(.data[[endpoint]], na.rm = TRUE)
      ) %>%
      ungroup()
    centered_cor <- if (nrow(centered) >= 5 &&
                        safe_sd(centered$cna_centered) > 0 &&
                        safe_sd(centered$endpoint_centered) > 0) {
      tryCatch(cor.test(centered$cna_centered, centered$endpoint_centered, method = "spearman"),
               error = function(e) NULL)
    } else {
      NULL
    }
    data.frame(
      cna_metric = cna_metric,
      cna_metric_label = cna_metric_labels[[cna_metric]],
      endpoint = endpoint,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      n_pairs = nrow(df),
      spearman_rho = if (!is.null(cor_res)) unname(cor_res$estimate) else NA_real_,
      spearman_p_value = if (!is.null(cor_res)) cor_res$p.value else NA_real_,
      sample_centered_rho = if (!is.null(centered_cor)) unname(centered_cor$estimate) else NA_real_,
      sample_centered_p_value = if (!is.null(centered_cor)) centered_cor$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  group_by(.data$feature_group, .data$cna_metric) %>%
  mutate(
    sample_centered_p_adj_group = p.adjust(.data$sample_centered_p_value, method = "BH"),
    spearman_p_adj_group = p.adjust(.data$spearman_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    sample_centered_p_adj_global = p.adjust(.data$sample_centered_p_value, method = "BH"),
    spearman_p_adj_global = p.adjust(.data$spearman_p_value, method = "BH"),
    sig_label = p_to_stars(.data$sample_centered_p_adj_group),
    neglog10_fdr = pmin(-log10(pmax(.data$sample_centered_p_adj_group, 1e-12)), 12)
  )

write_table(pairwise_feature_tests_v2, "Auto_v2_pairwise_cna_distance_all_feature_tests.csv")

message("Creating v2 consensus heatmap")

chrom_levels <- c(paste0("chr", 1:22), "chrX")
arm_levels <- colnames(arm_matrix)
chr_from_arm <- sub("[pq]$", "", arm_levels)
arm_matrix[!is.finite(arm_matrix)] <- 0
valid_rows <- rownames(arm_matrix)[rowSums(is.finite(arm_matrix)) == ncol(arm_matrix)]
arm_matrix_valid <- arm_matrix[valid_rows, , drop = FALSE]
arm_call_wide_v2 <- arm_long %>%
  select(.data$subclone_id, .data$arm, .data$cna_call) %>%
  pivot_wider(names_from = .data$arm, values_from = .data$cna_call, values_fill = 0) %>%
  as.data.frame()
rownames(arm_call_wide_v2) <- arm_call_wide_v2$subclone_id
arm_call_matrix_v2 <- as.matrix(arm_call_wide_v2[, setdiff(colnames(arm_call_wide_v2), "subclone_id"), drop = FALSE])
missing_call_arms <- setdiff(arm_levels, colnames(arm_call_matrix_v2))
if (length(missing_call_arms) > 0) {
  arm_call_matrix_v2 <- cbind(
    arm_call_matrix_v2,
    matrix(0L, nrow = nrow(arm_call_matrix_v2), ncol = length(missing_call_arms),
           dimnames = list(rownames(arm_call_matrix_v2), missing_call_arms))
  )
}
arm_call_matrix_v2 <- arm_call_matrix_v2[rownames(arm_matrix_valid), arm_levels, drop = FALSE]
hc <- hclust(dist(arm_matrix_valid), method = "ward.D2")
cluster_df <- data.frame(
  subclone_id = names(cna_cluster),
  cna_cluster = unname(cna_cluster),
  stringsAsFactors = FALSE
) %>%
  left_join(features %>% select(.data$subclone_id, .data$sample, .data$subclone,
                                .data$subclone_fraction, .data$dominance_class),
            by = "subclone_id")
cluster_df <- cluster_df[match(rownames(arm_matrix_valid), cluster_df$subclone_id), , drop = FALSE]

cluster_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(cluster_df$cna_cluster))),
  sort(unique(cluster_df$cna_cluster))
)
dominance_cols <- c(
  "single_subclone" = "grey70",
  "largest_not_significant" = "#FDB863",
  "significant_dominant" = "#B2182B"
)

plot_consensus_heatmap_v2 <- function() {
  row_meta <- cluster_df[hc$order, , drop = FALSE]
  mat <- arm_matrix_valid[hc$order, arm_levels, drop = FALSE]
  row_ha <- rowAnnotation(
    Cluster = row_meta$cna_cluster,
    Dominance = row_meta$dominance_class,
    `Clone fraction` = row_meta$subclone_fraction,
    col = list(
      Cluster = cluster_cols,
      Dominance = dominance_cols,
      `Clone fraction` = colorRamp2(c(0, 0.5, 1), c("white", "#FDB863", "#B2182B"))
    ),
    annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
    annotation_legend_param = list(
      Cluster = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      Dominance = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      `Clone fraction` = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12))
    ),
    simple_anno_size = unit(5, "mm")
  )
  Heatmap(
    mat,
    name = "Mean CNA",
    col = colorRamp2(c(-0.18, 0, 0.18), c("#2166AC", "white", "#B2182B")),
    left_annotation = row_ha,
    row_split = factor(row_meta$cna_cluster, levels = sort(unique(row_meta$cna_cluster))),
    row_title = NULL,
    row_title_gp = gpar(fontsize = 0, col = NA),
    row_gap = unit(1.2, "mm"),
    column_split = factor(chr_from_arm, levels = chrom_levels),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA)
  )
}

pdf(file.path(figure_dir, "Auto_cna_consensus_heatmap.pdf"),
    width = 17, height = 10, useDingbats = FALSE)
draw(plot_consensus_heatmap_v2(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

prepare_cna_matrix_v2 <- function(outs, cells) {
  cells <- intersect(cells, colnames(outs))
  common_genes <- intersect(rownames(outs), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(outs[go$gene, cells, drop = FALSE])
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  signal <- rowMeans(abs(mat), na.rm = TRUE)
  keep_signal <- signal >= as.numeric(quantile(signal, probs = 1 / 3, na.rm = TRUE))
  list(mat = mat[keep_signal, , drop = FALSE], gene_order = go[keep_signal, , drop = FALSE])
}

make_binned_cna_v2 <- function(cna_mat, go, bin_size = 100L) {
  go2 <- go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bin_levels <- unique(go2$bin)
  bins <- split(go2$.row, factor(go2$bin, levels = bin_levels))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  row_chr <- sub("_.*$", "", rownames(binned))
  list(mat = binned, chr = row_chr)
}

order_cells_by_subclone_v2 <- function(cna_mat, subclone) {
  split_cells <- split(names(subclone), factor(subclone, levels = unique(subclone)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    d <- dist(t(cna_mat[, cells, drop = FALSE]))
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

sample_plot_cells_v2 <- function(cells, subclone, max_cells = 1200L) {
  cells <- intersect(cells, names(subclone))
  if (length(cells) <= max_cells) return(cells)
  split_cells <- split(cells, factor(subclone[cells], levels = unique(subclone[cells])))
  target <- pmax(20L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

event_matrix_for_cells <- function(meta_plot, event_ids) {
  subclone_ids <- paste(meta_plot$sample, meta_plot$subclone, sep = "::")
  out <- do.call(rbind, lapply(event_ids, function(ev) {
    event_direction <- sub("_.*$", "", ev)
    event_arm <- sub("^[^_]+_", "", ev)
    calls <- arm_call_matrix_v2[subclone_ids, event_arm]
    present <- if (event_direction == "gain") calls == 1L else calls == -1L
    ifelse(present, ifelse(event_direction == "gain", 1L, -1L), 0L)
  }))
  rownames(out) <- as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])
  colnames(out) <- rownames(meta_plot)
  out
}

plot_sample_cell_heatmap <- function(sample_id) {
  sample_cells <- cell_df %>%
    filter(.data$sample == sample_id, .data$subclone != "Unresolved")
  if (nrow(sample_cells) == 0) return(NULL)
  outs_path <- file.path("by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  if (!file.exists(outs_path)) return(NULL)
  outs <- readRDS(outs_path)
  cells <- intersect(sample_cells$cell, colnames(outs))
  if (length(cells) == 0) return(NULL)
  sample_cells <- sample_cells[match(cells, sample_cells$cell), , drop = FALSE]
  subclone <- sample_cells$subclone
  names(subclone) <- sample_cells$cell
  cna_prepped <- prepare_cna_matrix_v2(outs, cells)
  cna_order <- order_cells_by_subclone_v2(cna_prepped$mat, subclone)
  plot_cells <- sample_plot_cells_v2(cna_order, subclone, max_plot_cells)
  binned <- make_binned_cna_v2(cna_prepped$mat[, cna_order, drop = FALSE], cna_prepped$gene_order)
  meta_plot <- sample_cells[match(plot_cells, sample_cells$cell), , drop = FALSE]
  rownames(meta_plot) <- meta_plot$cell
  meta_plot$subclone <- factor(meta_plot$subclone, levels = unique(subclone[cna_order]))
  meta_plot$top_mp_label <- factor(meta_plot$top_mp_label, levels = unique(meta_plot$top_mp_label))
  meta_plot$state_label <- factor(meta_plot$state_label, levels = unique(meta_plot$state_label))
  mat <- binned$mat[, rownames(meta_plot), drop = FALSE]
  row_chr <- factor(binned$chr, levels = unique(binned$chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr))), levels(row_chr))
  subclone_cols <- subclone_colours(meta_plot$subclone)
  topmp_cols <- mp_cols[names(mp_cols) %in% unique(as.character(meta_plot$top_mp_label))]
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(meta_plot$state_label))]
  missing_topmp <- setdiff(unique(as.character(meta_plot$top_mp_label)), names(topmp_cols))
  if (length(missing_topmp) > 0) topmp_cols <- c(topmp_cols, setNames(hue_pal()(length(missing_topmp)), missing_topmp))
  missing_states <- setdiff(unique(as.character(meta_plot$state_label)), names(local_state_cols))
  if (length(missing_states) > 0) local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  top_ha <- HeatmapAnnotation(
    Subclone = meta_plot$subclone,
    TopMP = meta_plot$top_mp_label,
    State = meta_plot$state_label,
    col = list(Subclone = subclone_cols, TopMP = topmp_cols, State = local_state_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(3, "mm"),
    na_col = "grey90"
  )
  left_ha <- rowAnnotation(
    Chr = row_chr,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )
  cna_ht <- Heatmap(
    mat,
    name = "CNA",
    col = colorRamp2(c(-cna_colour_limit, 0, cna_colour_limit), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_ha,
    left_annotation = left_ha,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    column_split = factor(meta_plot$subclone, levels = unique(meta_plot$subclone)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_title_gp = gpar(fontsize = 7),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA),
    column_title = sample_id
  )
  event_mat <- event_matrix_for_cells(meta_plot, recurrent_events)
  event_ht <- Heatmap(
    event_mat,
    name = "Event",
    col = c("-1" = "#2166AC", "0" = "white", "1" = "#B2182B"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      labels = c("Loss present", "Absent", "Gain present"),
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    height = unit(max(2.2, 0.32 * nrow(event_mat)), "cm"),
    border = TRUE,
    rect_gp = gpar(col = "grey70", lwd = 0.4)
  )
  cna_ht %v% event_ht
}

pdf(file.path(figure_dir, "Auto_per_sample_heatmap_recurrent_events.pdf"),
    width = 14, height = 9, useDingbats = FALSE)
for (sample_id in sort(unique(features$sample))) {
  ht <- plot_sample_cell_heatmap(sample_id)
  if (!is.null(ht)) draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()

message("Creating v2 recurrent CNA event plots")

event_order <- event_meta$event_label
event_feature_tests_v2 <- event_feature_tests_v2 %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features))))
  )
event_sample_deltas <- event_sample_deltas %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features))))
  )
event_feature_values <- event_feature_values %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features)))),
    event_group = factor(ifelse(.data$event_present, "Event-positive", "Event-negative"),
                         levels = c("Event-negative", "Event-positive"))
  )

plot_event_assoc <- function(group_name, title_suffix) {
  df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  ggplot(df, aes(.data$event_label, .data$feature_label)) +
    geom_point(aes(size = .data$neglog10_fdr, fill = .data$primary_delta),
               shape = 21, color = "black", stroke = 0.45, alpha = 0.95) +
    geom_text(aes(label = .data$sig_label), size = 6.2, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         na.value = "grey90") +
    scale_size_continuous(range = c(3.5, 11), limits = c(0, 12)) +
    labs(
      title = paste0("Recurrent CNA event associations: ", title_suffix),
      x = NULL,
      y = NULL,
      fill = "Median standardized\nevent delta",
      size = "-log10(FDR)"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15)
    )
}

event_page_chunks <- function(event_ids, page_size = 4L) {
  if (length(event_ids) == 0) return(list(character(0)))
  split(event_ids, ceiling(seq_along(event_ids) / page_size))
}

plot_event_boxplots <- function(group_name, title_suffix, event_ids) {
  df <- event_feature_values %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])))
  sig_df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)]))) %>%
    left_join(
      df %>%
        group_by(.data$event_id, .data$feature) %>%
        summarise(star_y = max(.data$feature_z, na.rm = TRUE) + 1.00, .groups = "drop"),
      by = c("event_id", "feature")
    )
  ggplot(df, aes(x = .data$feature_label, y = .data$feature_z, fill = .data$event_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.88,
                 position = position_dodge(width = 0.72), linewidth = 0.55) +
    geom_point(aes(color = .data$event_group),
               position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0.02,
                                               dodge.width = 0.72),
               alpha = 0.36, size = 1.3, show.legend = FALSE) +
    geom_text(data = sig_df, aes(x = .data$feature_label, y = .data$star_y, label = .data$sig_label),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~event_label, ncol = 2) +
    scale_fill_manual(values = c("Event-negative" = "grey72", "Event-positive" = "#B2182B")) +
    scale_color_manual(values = c("Event-negative" = "grey45", "Event-positive" = "#B2182B")) +
    labs(
      title = paste0("Recurrent CNA event feature distributions: ", title_suffix),
      x = NULL,
      y = "Standardized subclone feature value",
      fill = NULL
    ) +
    theme_classic(base_size = 20) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(face = "bold", size = 24),
      legend.position = "top",
      legend.text = element_text(size = 16)
    )
}

plot_event_sample_deltas <- function(group_name, title_suffix, event_ids) {
  df <- event_sample_deltas %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])))
  sig_df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)]))) %>%
    left_join(
      df %>%
        group_by(.data$event_id, .data$feature) %>%
        summarise(star_y = max(.data$std_delta, na.rm = TRUE) + 0.20, .groups = "drop"),
      by = c("event_id", "feature")
    )
  ggplot(df, aes(x = .data$feature_label, y = .data$std_delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.62, fill = "grey84", color = "black", linewidth = 0.55) +
    geom_point(position = position_jitter(width = 0.12, height = 0), alpha = 0.48,
               size = 1.8, color = "#B2182B") +
    geom_text(data = sig_df, aes(x = .data$feature_label, y = .data$star_y, label = .data$sig_label),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~event_label, ncol = 2) +
    labs(
      title = paste0("Per-sample recurrent CNA event deltas: ", title_suffix),
      x = NULL,
      y = "Standardized paired delta"
    ) +
    theme_classic(base_size = 20) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(face = "bold", size = 24)
    )
}

event_bar <- event_meta %>%
  mutate(event_label = factor(.data$event_label, levels = rev(as.character(.data$event_label)))) %>%
  ggplot(aes(.data$frac_samples_event, .data$event_label, fill = .data$direction)) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(aes(label = percent(.data$frac_samples_event, accuracy = 1)), hjust = -0.08, size = 5) +
  scale_x_continuous(labels = percent, limits = c(0, min(1, max(event_meta$frac_samples_event) * 1.18))) +
  scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
  labs(title = "Recurrent arm-level CNA events used for association testing",
       x = "Fraction of samples with event in at least one subclone", y = NULL, fill = NULL) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 24),
        axis.text = element_text(size = 16),
        legend.position = "top")

pdf(file.path(figure_dir, "Auto_recurrent_cna_event_associations_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(event_bar)
print(plot_event_assoc("Metaprogrammes", "all metaprogrammes"))
print(plot_event_assoc("Centred states", "five centred states excluding Hybrid and Unresolved"))
print(plot_event_assoc("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

pdf(file.path(figure_dir, "Auto_recurrent_cna_event_boxplots_all_features.pdf"),
    width = 22, height = 14, useDingbats = FALSE)
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Metaprogrammes", "all metaprogrammes", event_page))
}
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Centred states", "five centred states excluding Hybrid and Unresolved", event_page))
}
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("QC / CNA metrics", "QC and CNA metrics", event_page))
}
dev.off()

pdf(file.path(figure_dir, "Auto_recurrent_cna_event_per_sample_deltas.pdf"),
    width = 22, height = 14, useDingbats = FALSE)
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("Metaprogrammes", "all metaprogrammes", event_page))
}
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("Centred states", "five centred states excluding Hybrid and Unresolved", event_page))
}
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("QC / CNA metrics", "QC and CNA metrics", event_page))
}
dev.off()

chr8q_myc <- arm_long %>%
  filter(.data$arm == "chr8q") %>%
  select(.data$sample, .data$subclone, .data$subclone_id, .data$cna_call) %>%
  left_join(features %>% select(.data$sample, .data$subclone, .data$subclone_id, .data[["mp__MP2+"]], .data$n_cells),
            by = c("sample", "subclone", "subclone_id")) %>%
  mutate(
    chr8q_group = case_when(
      .data$cna_call == 1L ~ "8q gain",
      .data$cna_call == -1L ~ "8q loss",
      TRUE ~ "No 8q CNA"
    ),
    chr8q_group = factor(.data$chr8q_group, levels = c("8q loss", "No 8q CNA", "8q gain"))
  )
if (nrow(chr8q_myc) > 0) {
  group_n <- chr8q_myc %>%
    group_by(.data$chr8q_group) %>%
    summarise(n = n(), y = max(.data[["mp__MP2+"]], na.rm = TRUE) + 0.02, .groups = "drop")
  p_gain8q_myc <- ggplot(chr8q_myc, aes(.data$chr8q_group, .data[["mp__MP2+"]], fill = .data$chr8q_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.88, linewidth = 0.7, width = 0.62) +
    geom_point(aes(size = .data$n_cells), position = position_jitter(width = 0.14, height = 0),
               alpha = 0.42, color = "black") +
    geom_text(data = group_n, aes(x = .data$chr8q_group, y = .data$y, label = paste0("n=", .data$n)),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    scale_fill_manual(values = c("8q loss" = "#2166AC", "No 8q CNA" = "grey72", "8q gain" = "#B2182B")) +
    scale_size_continuous(range = c(1.8, 6)) +
    labs(
      title = "chr8q CNA status versus MYC-related proliferation MP",
      x = NULL,
      y = "Subclone mean MP2+ score",
      fill = "chr8q CNA",
      size = "Cells"
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 24),
      axis.text = element_text(size = 18),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15),
      legend.position = "right"
    )
  pdf(file.path(figure_dir, "Auto_gain_chr8q_myc_mp_per_sample.pdf"), width = 12, height = 8, useDingbats = FALSE)
  print(p_gain8q_myc)
  dev.off()
}

message("Creating v2 largest-subclone and pairwise-distance plots")

dominant_tests_v2 <- dominant_tests_v2 %>%
  mutate(
    feature_label = factor(.data$feature_label,
                           levels = rev(unique(feature_label(plot_features)))),
    feature_group = factor(.data$feature_group, levels = c("Metaprogrammes", "Centred states", "QC / CNA metrics"))
  )

plot_dominant_group <- function(group_name, title_suffix) {
  df <- dominant_deltas_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  sig_df <- dominant_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))),
           star_x = max(df$std_delta, na.rm = TRUE) + 0.15)
  ggplot(df, aes(.data$std_delta, .data$feature_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey86", color = "black", linewidth = 0.65) +
    geom_point(aes(color = .data$dominance_class),
               position = position_jitter(height = 0.12, width = 0), alpha = 0.58, size = 2.2) +
    geom_text(data = sig_df, aes(x = .data$star_x, y = .data$feature_label, label = .data$sig_label),
              inherit.aes = FALSE, size = 6.5, fontface = "bold") +
    scale_color_manual(
      values = dominance_cols,
      breaks = c("single_subclone", "largest_not_significant", "significant_dominant"),
      labels = c("Single subclone", "Largest not dominant", "Dominant largest")
    ) +
    labs(
      title = paste0("Largest subclone versus other subclones: ", title_suffix),
      x = "Standardized largest-minus-rest delta",
      y = NULL,
      color = "Largest-subclone class"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15),
      legend.position = "top"
    )
}

pdf(file.path(figure_dir, "Auto_largest_subclone_effects_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(plot_dominant_group("Metaprogrammes", "all metaprogrammes"))
print(plot_dominant_group("Centred states", "five centred states excluding Hybrid and Unresolved"))
print(plot_dominant_group("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

pairwise_feature_tests_v2 <- pairwise_feature_tests_v2 %>%
  mutate(
    cna_metric_label = factor(.data$cna_metric_label, levels = cna_metric_labels[cna_distance_metrics]),
    feature_group = factor(.data$feature_group, levels = c("Metaprogrammes", "Centred states", "QC / CNA metrics"))
  )

plot_pairwise_group <- function(group_name, title_suffix) {
  df <- pairwise_feature_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  ggplot(df, aes(.data$cna_metric_label, .data$feature_label)) +
    geom_point(aes(size = .data$neglog10_fdr, fill = .data$sample_centered_rho),
               shape = 21, color = "black", stroke = 0.45) +
    geom_text(aes(label = .data$sig_label), size = 6.1, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), na.value = "grey90") +
    scale_size_continuous(range = c(3.5, 11), limits = c(0, 12)) +
    labs(
      title = paste0("Pairwise CNA distance versus subclone divergence: ", title_suffix),
      x = NULL,
      y = NULL,
      fill = "Sample-centered rho",
      size = "-log10(FDR)"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15)
    )
}

pdf(file.path(figure_dir, "Auto_pairwise_cna_distance_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(plot_pairwise_group("Metaprogrammes", "all metaprogrammes"))
print(plot_pairwise_group("Centred states", "five centred states excluding Hybrid and Unresolved"))
print(plot_pairwise_group("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

run_summary_v2 <- data.frame(
  metric = c(
    "analysed_samples",
    "analysed_subclones",
    "recurrent_events_plotted",
    "features_plotted",
    "mp_features",
    "state_features_excluding_hybrid_unresolved",
    "qc_cna_features",
    "pairwise_subclone_comparisons"
  ),
  value = c(
    dplyr::n_distinct(features$sample),
    nrow(features),
    length(recurrent_events),
    length(plot_features),
    length(mp_features),
    length(state_features),
    length(qc_features),
    nrow(pairwise_df)
  )
)

write_table(run_summary_v2, "Auto_run_summary.csv")

summary_dir <- file.path("..", "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
event_count_summary <- event_feature_tests_v2 %>%
  distinct(.data$event_id, .data$event_label, .data$n_subclones_event,
           .data$n_subclones_no_event, .data$n_paired_samples) %>%
  mutate(
    summary_type = "recurrent_event_count",
    item = .data$event_id,
    metric = "event_positive_event_negative_paired_samples",
    value = paste(.data$n_subclones_event, .data$n_subclones_no_event,
                  .data$n_paired_samples, sep = "/"),
    detail = as.character(.data$event_label)
  ) %>%
  select(.data$summary_type, .data$item, .data$metric, .data$value, .data$detail)

top_assoc_summary <- event_feature_tests_v2 %>%
  filter(is.finite(.data$unpaired_p_adj_group)) %>%
  arrange(.data$unpaired_p_adj_group) %>%
  head(10) %>%
  mutate(
    summary_type = "top_recurrent_event_association",
    item = paste(.data$event_id, .data$feature, sep = "__"),
    metric = "unpaired_median_standardized_delta_group_fdr",
    value = paste(signif(.data$unpaired_median_std_delta, 4),
                  signif(.data$unpaired_p_adj_group, 4), sep = "/"),
    detail = paste(.data$event_label, .data$feature_label, sep = " | ")
  ) %>%
  select(.data$summary_type, .data$item, .data$metric, .data$value, .data$detail)

compact_summary <- bind_rows(
  run_summary_v2 %>%
    mutate(
      summary_type = "run",
      item = .data$metric,
      detail = ""
    ) %>%
    transmute(.data$summary_type, .data$item, metric = "value",
              value = as.character(.data$value), .data$detail),
  event_count_summary,
  top_assoc_summary
)

write.csv(compact_summary, file.path(summary_dir, "Auto_cna_subclone_expression_summary.csv"), row.names = FALSE)
write.csv(compact_summary, file.path(summary_dir, "Auto_cna_subclone_expression_v2_summary.csv"), row.names = FALSE)

save_rds(
  list(
    event_presence = event_presence,
    event_feature_values = event_feature_values,
    event_sample_deltas = event_sample_deltas,
    event_feature_tests = event_feature_tests_v2,
    dominant_deltas = dominant_deltas_v2,
    dominant_tests = dominant_tests_v2,
    pairwise_feature_tests = pairwise_feature_tests_v2
  ),
  "Auto_v2_visualisation_results.rds"
)

message("Done. Presentation-quality plots and figures successfully generated under: ", out_dir)
