#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_rctd_annotation.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: Visium HD square_016um filtered_feature_bc_matrix.h5 files,
#     spatial tissue positions, and ref_outs/EAC_Ref_merged.rds.
#   Outputs: ref_outs/visium_hd_outs/rctd/{intermediate,tables,logs}/.
#   Cache/replot: per-sample RCTD RDS objects cache the doublet weights; tables
#     are sufficient for downstream epithelial and malignancy filtering.
#   Run: Rscript analysis/spatial/visiumhd_rctd_annotation.R --inputs <...>
#     --sample-names <...> --output-dir ref_outs/visium_hd_outs/rctd
#   Environment: dmtcp.
####################

####################
# RCTD annotation of 16 um Visium HD bins. Doublet-mode results are retained
# in full, while only RCTD singlet epithelial bins pass to CNA/state mapping.
####################
suppressPackageStartupMessages({
  library(arrow)
  library(Matrix)
  library(Seurat)
  library(spacexr)
})

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

parse_cli <- function(args) {
  values <- list(inputs = character(), sample_names = character())
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% c("inputs", "sample-names", "output-dir", "reference-path", "min-umis", "max-cores")) {
      stop("Unknown argument: ", args[[i]])
    }
    i <- i + 1L
    start <- i
    while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1L
    if (start == i) stop("Missing value for --", key)
    value <- args[seq.int(start, i - 1L)]
    values[[gsub("-", "_", key)]] <- value
  }
  values
}

get_counts <- function(obj) {
  suppressWarnings(tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  ))
}

read_10x_counts <- function(path) {
  counts <- Read10X_h5(path)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  Matrix::Matrix(counts, sparse = TRUE)
}

write_csv_gz <- function(data, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(data, con, row.names = FALSE)
}

read_spatial_positions <- function(spatial_dir, barcodes) {
  parquet_path <- file.path(spatial_dir, "tissue_positions.parquet")
  csv_path <- file.path(spatial_dir, "tissue_positions.csv")
  positions <- if (file.exists(parquet_path)) {
    as.data.frame(arrow::read_parquet(parquet_path))
  } else if (file.exists(csv_path)) {
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("No tissue_positions.parquet or tissue_positions.csv in ", spatial_dir)
  }
  barcode_col <- intersect(c("barcode", "Barcode"), colnames(positions))
  if (length(barcode_col) != 1L) stop("Spatial positions need exactly one barcode column")
  rownames(positions) <- positions[[barcode_col]]
  x_col <- intersect(c("pxl_col_in_fullres", "array_col"), colnames(positions))
  y_col <- intersect(c("pxl_row_in_fullres", "array_row"), colnames(positions))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    stop("Spatial positions lack supported coordinate columns")
  }
  shared <- intersect(barcodes, rownames(positions))
  if (length(shared) == 0L) stop("No count barcodes matched spatial positions")
  coords <- data.frame(
    x = as.numeric(positions[shared, x_col[[1]]]),
    y = as.numeric(positions[shared, y_col[[1]]]),
    row.names = shared
  )
  list(coords = coords, barcodes = shared, positions = positions[shared, , drop = FALSE])
}

build_reference <- function(reference_path) {
  message("Loading RCTD reference: ", reference_path)
  ref_obj <- readRDS(reference_path)
  label_col <- intersect(c("celltype_update", "celltype_manual", "celltype_group"), colnames(ref_obj@meta.data))
  if (length(label_col) == 0L) stop("No recognised cell-type column in EAC_Ref_merged metadata")
  label_col <- label_col[[1]]
  labels <- as.character(ref_obj@meta.data[[label_col]])
  names(labels) <- rownames(ref_obj@meta.data)
  labels[labels %in% c("t.cell", "nk.cell")] <- "t_nk.cell"
  valid <- !is.na(labels) & nzchar(labels) & labels != "unresolved_inconsistent"
  labels <- labels[valid]
  set.seed(666)
  sampled_cells <- unlist(lapply(split(names(labels), labels), function(cells) {
    n_keep <- if (labels[[cells[[1]]]] == "epithelial") 6000L else 2000L
    sample(cells, min(length(cells), n_keep))
  }), use.names = FALSE)
  counts <- get_counts(ref_obj)[, sampled_cells, drop = FALSE]
  counts <- round(counts)
  ref_labels <- factor(labels[sampled_cells])
  names(ref_labels) <- sampled_cells
  list(
    reference = spacexr::Reference(counts[, names(ref_labels), drop = FALSE], ref_labels, Matrix::colSums(counts)),
    label_column = label_col,
    sampled_counts = table(ref_labels)
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (length(args$inputs) == 0L || length(args$sample_names) == 0L) {
  stop("--inputs and --sample-names are required")
}
if (length(args$inputs) != length(args$sample_names)) {
  stop("--inputs and --sample-names must have the same length")
}

output_dir <- if (length(args$output_dir)) args$output_dir[[1]] else file.path(WD, "ref_outs", "visium_hd_outs", "rctd")
reference_path <- if (length(args$reference_path)) args$reference_path[[1]] else file.path(WD, "ref_outs", "EAC_Ref_merged.rds")
min_umis <- if (length(args$min_umis)) as.integer(args$min_umis[[1]]) else 200L
max_cores <- if (length(args$max_cores)) as.integer(args$max_cores[[1]]) else 8L
if (!file.exists(reference_path)) stop("Missing RCTD reference: ", reference_path)
if (!is.finite(min_umis) || min_umis < 1L) stop("--min-umis must be a positive integer")
if (!is.finite(max_cores) || max_cores < 1L) stop("--max-cores must be a positive integer")

dirs <- file.path(output_dir, c("intermediate", "tables", "logs"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
run_log <- c(paste0("start=", format(Sys.time(), tz = "Europe/London")), paste0("reference_path=", reference_path), paste0("min_umis=", min_umis))

reference_info <- build_reference(reference_path)
write.csv(
  data.frame(celltype = names(reference_info$sampled_counts), n_cells = as.integer(reference_info$sampled_counts)),
  file.path(output_dir, "tables", "Auto_visiumhd_rctd_reference_composition.csv"),
  row.names = FALSE
)

summary_rows <- list()
for (i in seq_along(args$inputs)) {
  input_dir <- normalizePath(args$inputs[[i]], mustWork = TRUE)
  sample_name <- args$sample_names[[i]]
  message("RCTD: ", sample_name)
  counts_path <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(counts_path)) stop("Missing binned count matrix: ", counts_path)
  counts <- read_10x_counts(counts_path)
  n_umis <- Matrix::colSums(counts)
  keep <- n_umis >= min_umis
  counts <- counts[, keep, drop = FALSE]
  n_umis <- n_umis[keep]
  spatial <- read_spatial_positions(file.path(input_dir, "spatial"), colnames(counts))
  counts <- counts[, spatial$barcodes, drop = FALSE]
  n_umis <- n_umis[spatial$barcodes]
  puck <- spacexr::SpatialRNA(spatial$coords[colnames(counts), , drop = FALSE], counts, n_umis)
  rctd <- spacexr::create.RCTD(puck, reference_info$reference, max_cores = max_cores, test_mode = FALSE)
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = "doublet")
  result_df <- as.data.frame(rctd@results$results_df)
  result_df$barcode <- rownames(result_df)
  result_df$sample <- sample_name
  result_df$total_counts <- as.numeric(n_umis[result_df$barcode])
  result_df$pxl_col_in_fullres <- spatial$coords[result_df$barcode, 1]
  result_df$pxl_row_in_fullres <- spatial$coords[result_df$barcode, 2]
  result_df$Auto_annotation_method <- "RCTD_doublet"
  result_df$Auto_annotation_celltype <- as.character(result_df$first_type)
  result_df$Auto_annotation_pass_doublet_filter <- as.character(result_df$spot_class) == "singlet"
  result_df$Auto_annotation_keep_epithelial <- result_df$Auto_annotation_pass_doublet_filter & result_df$Auto_annotation_celltype == "epithelial"
  result_df <- result_df[, c("barcode", "sample", "total_counts", "pxl_row_in_fullres", "pxl_col_in_fullres", setdiff(colnames(result_df), c("barcode", "sample", "total_counts", "pxl_row_in_fullres", "pxl_col_in_fullres")))]
  saveRDS(rctd, file.path(output_dir, "intermediate", paste0("Auto_", sample_name, "_binned_rctd.rds")))
  write_csv_gz(result_df, file.path(output_dir, "tables", paste0("Auto_", sample_name, "_binned_rctd_annotations.csv.gz")))
  class_counts <- as.data.frame(table(result_df$spot_class, useNA = "ifany"), stringsAsFactors = FALSE)
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    n_input_barcodes = length(keep),
    n_qc_barcodes = sum(keep),
    n_rctd_barcodes = nrow(result_df),
    n_singlets = sum(result_df$Auto_annotation_pass_doublet_filter, na.rm = TRUE),
    n_epithelial_singlets = sum(result_df$Auto_annotation_keep_epithelial, na.rm = TRUE),
    rctd_classes = paste(paste0(class_counts$Var1, "=", class_counts$Freq), collapse = ";"),
    stringsAsFactors = FALSE
  )
  run_log <- c(run_log, paste0("sample=", sample_name, "; rctd_barcodes=", nrow(result_df)))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "tables", "Auto_visiumhd_binned_rctd_summary.csv"), row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(summary_dir, "visiumhd_binned_rctd_annotation_summary.csv"), row.names = FALSE)
run_log <- c(run_log, paste0("end=", format(Sys.time(), tz = "Europe/London")), paste0("reference_label_column=", reference_info$label_column))
writeLines(run_log, file.path(output_dir, "logs", "Auto_visiumhd_binned_rctd_run_summary.txt"))
####################
