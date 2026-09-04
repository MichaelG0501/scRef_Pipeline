####################
# Analysis registry:
#   Status: active heavy/intermediate workflow.
#   Script: analysis/trajectory/scatlas_velocity_metadata.R
#   Methodology: analysis/methodology/trajectory/scatlas_velocity_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     ref_outs/EAC_Ref_epi.rds
#     ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#     /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files/*_possorted_genome_bam.bam
#   Outputs:
#     ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_cell_metadata.csv
#     ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_sample_manifest.csv
#     ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_raw_bam_missing_samples.csv
#     ref_outs/Auto_velocity_scATLAS/barcodes/<sample>_qc_barcodes.tsv
#   Cache/replot behavior: rebuilds lightweight CSV/barcode inputs every run.
#   Run command:
#     Rscript analysis/trajectory/scatlas_velocity_metadata.R
#   Conda environment:
#     dmtcp
####################

####################
# Export scATLAS per-cell metadata and per-sample barcode lists for RNA
# velocity. Direction summaries use the five primary noreg Approach B states,
# while UMAP plots retain the finalized state labels.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(file.path(root_dir, "ref_outs"))

out_dir <- "Auto_velocity_scATLAS"
for (subdir in c("tables", "barcodes", "logs", "figures", "intermediate")) {
  dir.create(file.path(out_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

raw_bam_root <- Sys.getenv(
  "SCATLAS_RAW_BAM_ROOT",
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files"
)

primary_states <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia",
  "Stress adaptive",
  "Cancer-cell immune mimicry"
)

message("Loading scATLAS epithelial object and state vectors...")
epi <- readRDS("EAC_Ref_epi.rds")
state_final <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
state_noreg <- state_final

state_final_names <- names(state_final)
state_final <- as.character(state_final)
names(state_final) <- state_final_names

state_noreg_names <- names(state_noreg)
state_noreg <- as.character(state_noreg)
names(state_noreg) <- state_noreg_names

cells <- colnames(epi)
meta <- epi@meta.data[cells, , drop = FALSE]
meta$cell_id <- cells
meta$sample <- as.character(meta$orig.ident)

raw_bams <- list.files(
  raw_bam_root,
  pattern = "_possorted_genome_bam\\.bam$",
  full.names = TRUE
)
raw_manifest <- data.frame(
  sample = sub("_possorted_genome_bam\\.bam$", "", basename(raw_bams)),
  bam = normalizePath(raw_bams, mustWork = FALSE),
  stringsAsFactors = FALSE
) %>%
  arrange(.data$sample)

samples_in_object <- sort(unique(meta$sample))
samples_with_bam <- intersect(samples_in_object, raw_manifest$sample)
if (length(samples_with_bam) == 0) {
  stop("No epithelial samples have matching raw BAMs under: ", raw_bam_root)
}

message("Matched ", length(samples_with_bam), " epithelial samples with raw BAMs.")
meta <- meta %>% filter(.data$sample %in% samples_with_bam)
cells <- meta$cell_id

meta$raw_barcode <- mapply(
  function(cell_id, sample) {
    prefix <- paste0(sample, "_")
    if (startsWith(cell_id, prefix)) {
      substring(cell_id, nchar(prefix) + 1)
    } else {
      sub("^.*_", "", cell_id)
    }
  },
  meta$cell_id,
  meta$sample,
  USE.NAMES = FALSE
)
meta$state_final <- state_final[meta$cell_id]
meta$state_noreg <- state_noreg[meta$cell_id]

message("Computing dominant MPs for Basal and SMG states...")
mp_adj <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds")

basal_mps <- c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+")
smg_mps <- c("MP8+", "MP8b", "MP16", "MP18b", "MP17")

meta$dominant_basal_mp <- NA_character_
basal_cells <- meta$cell_id[meta$state_noreg == "Basal to intestinal metaplasia"]
basal_cells <- intersect(basal_cells, rownames(mp_adj))
if (length(basal_cells) > 0) {
  basal_scores <- mp_adj[basal_cells, intersect(basal_mps, colnames(mp_adj)), drop = FALSE]
  basal_best_idx <- max.col(basal_scores, ties.method = "first")
  basal_sorted <- t(apply(basal_scores, 1, sort, decreasing = TRUE))
  basal_gap <- basal_sorted[, 1] - basal_sorted[, 2]
  basal_dom <- colnames(basal_scores)[basal_best_idx]
  basal_dom[basal_gap < 0.1] <- "Hybrid"
  meta[basal_cells, "dominant_basal_mp"] <- basal_dom
}

meta$dominant_smg_mp <- NA_character_
smg_cells <- meta$cell_id[meta$state_noreg == "SMG to intestinal metaplasia"]
smg_cells <- intersect(smg_cells, rownames(mp_adj))
if (length(smg_cells) > 0) {
  smg_scores <- mp_adj[smg_cells, intersect(smg_mps, colnames(mp_adj)), drop = FALSE]
  smg_best_idx <- max.col(smg_scores, ties.method = "first")
  smg_sorted <- t(apply(smg_scores, 1, sort, decreasing = TRUE))
  smg_gap <- smg_sorted[, 1] - smg_sorted[, 2]
  smg_dom <- colnames(smg_scores)[smg_best_idx]
  smg_dom[smg_gap < 0.1] <- "Hybrid"
  meta[smg_cells, "dominant_smg_mp"] <- smg_dom
}

fallback_cols <- c(
  "cell_id", "sample", "raw_barcode", "study", "Author", "Year", "Patient",
  "Sample Name", "Sampling Time", "Treatment", "Treatment Strategy",
  "Diagnosis", "Tumor Type", "malignancy", "celltype_update",
  "state_noreg", "state_final", "dominant_basal_mp", "dominant_smg_mp",
  "nCount_RNA", "nFeature_RNA", "percent.mt"
)
keep_cols <- fallback_cols[fallback_cols %in% colnames(meta)]
meta_out <- meta[, keep_cols, drop = FALSE]
fwrite(meta_out, file.path(out_dir, "tables", "Auto_scatlas_velocity_cell_metadata.csv"))

for (sample_name in samples_with_bam) {
  barcodes <- sort(unique(meta_out$raw_barcode[meta_out$sample == sample_name]))
  writeLines(barcodes, file.path(out_dir, "barcodes", paste0(sample_name, "_qc_barcodes.tsv")))
}

sample_counts <- data.frame(
  sample = names(table(meta_out$sample)),
  n_cells = as.integer(table(meta_out$sample)),
  stringsAsFactors = FALSE
)

manifest <- raw_manifest %>%
  filter(.data$sample %in% samples_with_bam) %>%
  mutate(
    dataset = sub("^([^_]+_[0-9]{4}).*$", "\\1", .data$sample),
    study = .data$dataset,
    bai = paste0(.data$bam, ".bai"),
    barcodes_file = file.path(root_dir, "ref_outs", out_dir, "barcodes", paste0(.data$sample, "_qc_barcodes.tsv")),
    has_bam = file.exists(.data$bam),
    has_bai = file.exists(.data$bai)
  ) %>%
  left_join(sample_counts, by = "sample") %>%
  select(
    "sample", "dataset", "study", "bam", "bai",
    "barcodes_file", "n_cells", "has_bam", "has_bai"
  )
fwrite(manifest, file.path(out_dir, "tables", "Auto_scatlas_velocity_sample_manifest.csv"))

missing_tbl <- data.frame(sample = samples_in_object, stringsAsFactors = FALSE) %>%
  mutate(has_raw_bam = .data$sample %in% samples_with_bam) %>%
  filter(!.data$has_raw_bam)
fwrite(missing_tbl, file.path(out_dir, "tables", "Auto_scatlas_velocity_raw_bam_missing_samples.csv"))

summary_tbl <- meta_out %>%
  count(.data$sample, .data$study, .data$state_noreg, name = "cells") %>%
  arrange(.data$study, .data$sample, .data$state_noreg)
fwrite(summary_tbl, file.path(out_dir, "tables", "Auto_scatlas_velocity_state_summary.csv"))

message("Wrote velocity metadata for ", nrow(meta_out), " cells and ", nrow(manifest), " raw-BAM samples.")
