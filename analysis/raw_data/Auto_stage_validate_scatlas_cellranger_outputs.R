####################
# Auto_stage_validate_scatlas_cellranger_outputs.R
#
# Analysis registry
# Status: active upstream validation.
# Short description: mirror the historical Cell Ranger post-processing for
# Carroll/Alcindor by staging filtered_feature_bc_matrix outputs as matrix_all,
# optionally exporting count CSVs with the original write.sh logic, and comparing
# new 10x matrices against the live reference matrices before Numbat.
# Methodology: analysis/methodology/raw_data/scatlas_raw_redownload_numbat_methodology.md
# Inputs:
# - /rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/<dataset>/cellranger*/<sample>/outs/filtered_feature_bc_matrix
# - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/<dataset>/matrix_all/<sample>_filtered
# Outputs:
# - /rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/<dataset>/matrix_all*/<sample>_filtered
# - /rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/00_counts_matrix_all/<dataset>_<sample>.csv when SCATLAS_EXPORT_COUNT_CSV=TRUE
# - /rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/validation/Auto_scatlas_cellranger_matrix_validation.csv
# Cache/replot behavior: existing staged files are left in place; set
# SCATLAS_EXPORT_COUNT_CSV=TRUE to write historical-style dense CSVs.
# Run command: Rscript analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.R
# Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

raw_root <- Sys.getenv("SCATLAS_RAW_ROOT", "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files")
live_root <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all"
datasets <- c("Alcindor_2025", "Carroll_2023")
export_count_csv <- identical(toupper(Sys.getenv("SCATLAS_EXPORT_COUNT_CSV", "FALSE")), "TRUE")
validation_prefix <- Sys.getenv("SCATLAS_VALIDATION_PREFIX", "Auto_scatlas_cellranger_matrix_validation")

dataset_cr_root <- function(dataset) {
  env_name <- paste0("SCATLAS_", toupper(gsub("[^A-Za-z0-9]", "_", dataset)), "_CELLRANGER_ROOT")
  Sys.getenv(env_name, file.path(raw_root, dataset, "cellranger"))
}

dataset_stage_root <- function(dataset) {
  env_name <- paste0("SCATLAS_", toupper(gsub("[^A-Za-z0-9]", "_", dataset)), "_MATRIX_STAGE_ROOT")
  Sys.getenv(env_name, file.path(raw_root, dataset, "matrix_all"))
}

validation_dir <- file.path(raw_root, "validation")
counts_out_dir <- file.path(raw_root, "00_counts_matrix_all")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
if (export_count_csv) dir.create(counts_out_dir, recursive = TRUE, showWarnings = FALSE)

copy_if_missing <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(to)) {
    ok <- file.copy(from, to, overwrite = FALSE, copy.date = TRUE)
    if (!ok) stop("Failed to copy ", from, " to ", to)
  }
}

stage_matrix <- function(dataset, sample_name) {
  src <- file.path(dataset_cr_root(dataset), sample_name, "outs", "filtered_feature_bc_matrix")
  dest <- file.path(dataset_stage_root(dataset), paste0(sample_name, "_filtered"))
  if (!dir.exists(src)) return(NA_character_)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  for (f in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")) {
    copy_if_missing(file.path(src, f), file.path(dest, f))
  }
  normalizePath(dest, mustWork = TRUE)
}

matrix_summary <- function(path) {
  counts <- Read10X(data.dir = path)
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  list(
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    nnzero = length(counts@x),
    total_umis = sum(counts@x),
    rownames = rownames(counts),
    colnames = colnames(counts),
    counts = counts
  )
}

same_sparse_matrix <- function(x, y) {
  if (!identical(dim(x), dim(y))) return(FALSE)
  if (!identical(rownames(x), rownames(y))) return(FALSE)
  if (!identical(colnames(x), colnames(y))) return(FALSE)
  diff <- x - y
  length(diff@x) == 0 || all(diff@x == 0)
}

export_counts_csv <- function(matrix_dir, dataset, sample_name) {
  output_file <- file.path(counts_out_dir, paste0(dataset, "_", sample_name, ".csv"))
  if (file.exists(output_file)) return(normalizePath(output_file, mustWork = TRUE))
  counts <- Read10X(data.dir = matrix_dir)
  fwrite(as.matrix(counts), output_file, row.names = TRUE)
  normalizePath(output_file, mustWork = TRUE)
}

rows <- list()
for (dataset in datasets) {
  cr_root <- dataset_cr_root(dataset)
  if (!dir.exists(cr_root)) next
  sample_dirs <- list.dirs(cr_root, recursive = FALSE, full.names = FALSE)
  sample_dirs <- sort(sample_dirs[nzchar(sample_dirs)])
  for (sample_name in sample_dirs) {
    message("Staging and validating ", dataset, " / ", sample_name)
    staged_dir <- stage_matrix(dataset, sample_name)
    live_dir <- file.path(live_root, dataset, "matrix_all", paste0(sample_name, "_filtered"))
    status <- "ok"
    note <- ""
    csv_path <- NA_character_

    if (is.na(staged_dir)) {
      status <- "missing_new_filtered_matrix"
      note <- "No new Cell Ranger filtered_feature_bc_matrix"
      rows[[length(rows) + 1L]] <- data.frame(dataset, sample = sample_name, status, note, stringsAsFactors = FALSE)
      next
    }
    if (!dir.exists(live_dir)) {
      status <- "missing_live_reference_matrix"
      note <- live_dir
      rows[[length(rows) + 1L]] <- data.frame(dataset, sample = sample_name, status, note, staged_dir, live_dir, stringsAsFactors = FALSE)
      next
    }

    new_sum <- matrix_summary(staged_dir)
    live_sum <- matrix_summary(live_dir)
    exact_match <- same_sparse_matrix(new_sum$counts, live_sum$counts)
    if (!exact_match) status <- "matrix_mismatch"

    if (export_count_csv && exact_match) {
      csv_path <- export_counts_csv(staged_dir, dataset, sample_name)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset,
      sample = sample_name,
      status = status,
      exact_sparse_match = exact_match,
      new_n_genes = new_sum$n_genes,
      live_n_genes = live_sum$n_genes,
      new_n_cells = new_sum$n_cells,
      live_n_cells = live_sum$n_cells,
      new_nnzero = new_sum$nnzero,
      live_nnzero = live_sum$nnzero,
      new_total_umis = new_sum$total_umis,
      live_total_umis = live_sum$total_umis,
      staged_dir = staged_dir,
      live_dir = live_dir,
      count_csv = csv_path,
      note = note,
      stringsAsFactors = FALSE
    )

    rm(new_sum, live_sum)
    gc()
  }
}

validation <- bind_rows(rows)
out_csv <- file.path(validation_dir, paste0(validation_prefix, ".csv"))
fwrite(validation, out_csv)

bad <- validation %>% filter(.data$status != "ok")
if (nrow(bad) > 0) {
  fwrite(bad, file.path(validation_dir, paste0(validation_prefix, "_failures.csv")))
  stop("Some staged matrices do not match historical live matrices. See ", out_csv)
}

message("All staged matrices match historical live matrices. Wrote: ", out_csv)
