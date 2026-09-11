####################
# Auto_scatlas_numbat_export_inputs.R
#
# Analysis registry
# Status: active upstream.
# Script: analysis/cnv/Auto_scatlas_numbat_export_inputs.R
# Short description: export scATLAS Carroll/Alcindor count matrices, cell maps,
# and BAM/barcode paths for Numbat SNP pileup and subclone calling.
# Methodology: analysis/methodology/cnv/scatlas_numbat_methodology.md
# Inputs:
# - ref_outs/EAC_Ref_epi.rds
# - /rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/*/cellranger/<raw_sample>/outs/possorted_genome_bam.bam
# Outputs:
# - ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv
# - ref_outs/Auto_scatlas_numbat/by_samples/<sample>/input/Auto_<sample>_counts_raw_barcodes.rds
# - ref_outs/Auto_scatlas_numbat/by_samples/<sample>/input/Auto_<sample>_cell_map.csv
# Cache/replot behavior: rerun overwrites only the exported Auto_ inputs and
# manifest for samples with available BAM/barcode inputs.
# Run command: Rscript analysis/cnv/Auto_scatlas_numbat_export_inputs.R
# Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
out_root <- file.path(root_dir, "ref_outs")
raw_root <- Sys.getenv("SCATLAS_RAW_ROOT", "/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat")
setwd(root_dir)

out_dir <- file.path(out_root, "Auto_scatlas_numbat")
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "by_samples"), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("SCATLAS_NUMBAT_SAMPLES", "all")

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

raw_sample_from <- function(sample_name) {
  if (grepl("^Alcindor_2025_", sample_name)) return(sub("^Alcindor_2025_", "", sample_name))
  if (grepl("^Carroll_2023_", sample_name)) return(sub("^Carroll_2023_", "", sample_name))
  sample_name
}

dataset_from <- function(sample_name) {
  if (grepl("^Alcindor_2025_", sample_name)) return("Alcindor_2025")
  if (grepl("^Carroll_2023_", sample_name)) return("Carroll_2023")
  NA_character_
}

resolve_cellranger_out <- function(dataset, raw_sample) {
  env_key <- paste0("SCATLAS_", gsub("[^A-Z0-9]", "_", toupper(dataset)), "_BAM_CELLRANGER_ROOT")
  root <- Sys.getenv(env_key, file.path(raw_root, dataset, "cellranger"))
  file.path(root, raw_sample, "outs")
}

atlas_path <- file.path(out_root, "EAC_Ref_epi.rds")
if (!file.exists(atlas_path)) stop("Missing merged epithelial atlas: ", atlas_path)
message("Loading merged epithelial atlas: ", atlas_path)
atlas <- readRDS(atlas_path)
atlas_meta <- atlas@meta.data
if (!("orig.ident" %in% colnames(atlas_meta))) stop("Atlas metadata has no orig.ident column.")

sample_counts <- sort(table(atlas_meta$orig.ident), decreasing = TRUE)
samples <- names(sample_counts)[grepl("^Alcindor_2025_|^Carroll_2023_.*tumour", names(sample_counts))]
if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  samples <- intersect(samples, requested)
}
samples <- sort(samples)
if (length(samples) == 0) stop("No Carroll/Alcindor scATLAS samples selected.")

min_numbat_cells <- as.integer(Sys.getenv("SCATLAS_NUMBAT_MIN_CELLS", "50"))
skip_small <- data.frame(
  sample = names(sample_counts),
  n_cells = as.integer(sample_counts),
  stringsAsFactors = FALSE
) %>%
  filter(.data$sample %in% samples, .data$n_cells < min_numbat_cells) %>%
  arrange(.data$sample)
if (nrow(skip_small) > 0) {
  fwrite(skip_small, file.path(out_dir, "Auto_scatlas_numbat_manifest_skipped_low_cell_count.csv"))
  samples <- setdiff(samples, skip_small$sample)
}
if (length(samples) == 0) stop("No selected samples have at least ", min_numbat_cells, " cells.")

atlas_counts <- get_counts(atlas)
if (!inherits(atlas_counts, "dgCMatrix")) atlas_counts <- as(atlas_counts, "dgCMatrix")

manifest_rows <- list()
for (sample_name in samples) {
  dataset <- dataset_from(sample_name)
  raw_sample <- raw_sample_from(sample_name)
  cellranger_out <- resolve_cellranger_out(dataset, raw_sample)
  bam <- file.path(cellranger_out, "possorted_genome_bam.bam")
  if (!file.exists(bam)) {
    warning("Missing Cell Ranger BAM for ", sample_name, " under ", cellranger_out)
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      sample = sample_name,
      dataset = dataset,
      raw_sample = raw_sample,
      bam = normalizePath(bam, mustWork = FALSE),
      barcodes_file = NA_character_,
      count_rds = NA_character_,
      metadata_csv = NA_character_,
      sample_out_dir = file.path(out_dir, "by_samples", sample_name),
      numbat_dir = file.path(out_dir, "by_samples", sample_name, "numbat"),
      allele_file = file.path(out_dir, "by_samples", sample_name, paste0(sample_name, "_allele_counts.tsv.gz")),
      clone_post_file = file.path(out_dir, "by_samples", sample_name, "numbat", "clone_post_2.tsv"),
      joint_post_file = file.path(out_dir, "by_samples", sample_name, "numbat", "joint_post_2.tsv"),
      n_cells = NA_integer_,
      has_bam = file.exists(bam),
      has_barcode_file = FALSE,
      stringsAsFactors = FALSE
    )
    next
  }

  sample_out <- file.path(out_dir, "by_samples", sample_name)
  input_dir <- file.path(sample_out, "input")
  numbat_dir <- file.path(sample_out, "numbat")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(numbat_dir, recursive = TRUE, showWarnings = FALSE)

  message("Exporting Numbat inputs for ", sample_name)
  sample_cells <- rownames(atlas_meta)[atlas_meta$orig.ident == sample_name]
  counts <- atlas_counts[, sample_cells, drop = FALSE]
  raw_barcodes <- sub(paste0("^", sample_name, "_"), "", colnames(counts))
  colnames(counts) <- raw_barcodes
  cell_ids <- paste(sample_name, raw_barcodes, sep = "_")

  count_rds <- file.path(input_dir, paste0("Auto_", sample_name, "_counts_raw_barcodes.rds"))
  metadata_csv <- file.path(input_dir, paste0("Auto_", sample_name, "_cell_map.csv"))
  barcodes_file <- file.path(input_dir, paste0("Auto_", sample_name, "_barcodes_raw.tsv"))
  saveRDS(counts, count_rds)
  writeLines(raw_barcodes, barcodes_file)

  meta <- atlas_meta[sample_cells, , drop = FALSE]
  meta$raw_barcode <- raw_barcodes
  meta$cell_id <- cell_ids
  meta$sample <- sample_name
  meta$dataset <- dataset
  meta$raw_sample <- raw_sample
  keep_cols <- intersect(
    c("cell_id", "sample", "dataset", "raw_sample", "raw_barcode", "orig.ident",
      "nCount_RNA", "nFeature_RNA", "percent.mt", "celltype_initial", "malignancy"),
    colnames(meta)
  )
  fwrite(meta[, keep_cols, drop = FALSE], metadata_csv)

  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    sample = sample_name,
    dataset = dataset,
    raw_sample = raw_sample,
    bam = normalizePath(bam, mustWork = TRUE),
    barcodes_file = normalizePath(barcodes_file, mustWork = TRUE),
    count_rds = normalizePath(count_rds, mustWork = TRUE),
    metadata_csv = normalizePath(metadata_csv, mustWork = TRUE),
    sample_out_dir = normalizePath(sample_out, mustWork = TRUE),
    numbat_dir = normalizePath(numbat_dir, mustWork = TRUE),
    allele_file = normalizePath(file.path(sample_out, paste0(sample_name, "_allele_counts.tsv.gz")), mustWork = FALSE),
    clone_post_file = normalizePath(file.path(numbat_dir, "clone_post_2.tsv"), mustWork = FALSE),
    joint_post_file = normalizePath(file.path(numbat_dir, "joint_post_2.tsv"), mustWork = FALSE),
    n_cells = ncol(counts),
    has_bam = TRUE,
    has_barcode_file = TRUE,
    stringsAsFactors = FALSE
  )

  rm(counts, meta)
  gc()
}

manifest <- bind_rows(manifest_rows) %>% arrange(.data$dataset, .data$sample)
if (nrow(manifest) == 0) stop("No Numbat manifest rows were created.")

manifest_path <- file.path(out_dir, "Auto_scatlas_numbat_manifest.csv")
fwrite(manifest, manifest_path)

bad <- manifest %>% filter(!.data$has_bam | !.data$has_barcode_file | is.na(.data$count_rds))
if (nrow(bad) > 0) {
  fwrite(bad, file.path(out_dir, "Auto_scatlas_numbat_manifest_missing_inputs.csv"))
  stop("Some samples are missing BAM/barcode/count inputs. See Auto_scatlas_numbat_manifest_missing_inputs.csv")
}

summary_tbl <- manifest %>% count(.data$dataset, name = "n_samples") %>% arrange(.data$dataset)
fwrite(summary_tbl, file.path(out_dir, "Auto_scatlas_numbat_manifest_summary.csv"))

message("Wrote Numbat manifest: ", manifest_path)
message("Samples exported: ", nrow(manifest))
