#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_rctd_doublet_detection.R
#   Description: Detect singlet and non-singlet 16 um Visium HD bins with RCTD.
#   Methodology: analysis/methodology/spatial/visium_hd_final_annotation_methodology.md
#   Inputs: analysis/spatial/visium_hd_samples.tsv; square_016um count matrices
#     and positions; ref_outs/EAC_Ref_merged.rds.
#   Outputs:
#     intermediate/: per-sample RCTD RDS objects.
#     tables/: per-bin RCTD calls, reference composition, and run summary.
#     logs/: lightweight run summary.
#   Cache/replot: completed live outputs are reused. The archived live RCTD
#     cache is imported once for matching samples when available.
#   Run: Rscript analysis/spatial/visium_hd_rctd_doublet_detection.R
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

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
  values <- list()
  allowed <- c(
    "manifest", "output-dir", "reference-path", "legacy-cache-dir",
    "min-umis", "max-cores", "force-rebuild"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% allowed) stop("Unknown argument: ", args[[i]])
    i <- i + 1L
    if (i > length(args) || startsWith(args[[i]], "--")) {
      stop("Missing value for --", key)
    }
    values[[gsub("-", "_", key)]] <- args[[i]]
    i <- i + 1L
  }
  values
}

get_counts <- function(object) {
  suppressWarnings(tryCatch(
    GetAssayData(object, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(object, assay = "RNA", slot = "counts")
  ))
}

read_10x_counts <- function(path) {
  counts <- Read10X_h5(path)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) {
      counts[["Gene Expression"]]
    } else {
      counts[[1]]
    }
  }
  Matrix::Matrix(counts, sparse = TRUE)
}

write_csv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE)
}

read_spatial_positions <- function(spatial_dir, barcodes) {
  parquet_path <- file.path(spatial_dir, "tissue_positions.parquet")
  csv_path <- file.path(spatial_dir, "tissue_positions.csv")
  positions <- if (file.exists(parquet_path)) {
    as.data.frame(arrow::read_parquet(parquet_path))
  } else if (file.exists(csv_path)) {
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("No tissue positions file in ", spatial_dir)
  }
  barcode_col <- intersect(c("barcode", "Barcode"), colnames(positions))
  if (length(barcode_col) != 1L) {
    stop("Spatial positions need exactly one barcode column in ", spatial_dir)
  }
  rownames(positions) <- positions[[barcode_col]]
  x_col <- intersect(c("pxl_col_in_fullres", "array_col"), colnames(positions))
  y_col <- intersect(c("pxl_row_in_fullres", "array_row"), colnames(positions))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    stop("Spatial positions lack supported coordinate columns in ", spatial_dir)
  }
  shared <- intersect(barcodes, rownames(positions))
  if (length(shared) == 0L) stop("No count barcodes matched spatial positions")
  list(
    barcodes = shared,
    coords = data.frame(
      x = as.numeric(positions[shared, x_col[[1]]]),
      y = as.numeric(positions[shared, y_col[[1]]]),
      row.names = shared
    )
  )
}

build_reference <- function(reference_path) {
  message("Loading RCTD reference: ", reference_path)
  object <- readRDS(reference_path)
  label_col <- intersect(
    c("celltype_update", "celltype_manual", "celltype_group"),
    colnames(object@meta.data)
  )
  if (length(label_col) == 0L) {
    stop("No recognised cell-type column in the RCTD reference")
  }
  label_col <- label_col[[1]]
  labels <- as.character(object@meta.data[[label_col]])
  names(labels) <- rownames(object@meta.data)
  labels[labels %in% c("t.cell", "nk.cell")] <- "t_nk.cell"
  valid <- !is.na(labels) & nzchar(labels) &
    !labels %in% c("unresolved", "unresolved_inconsistent")
  labels <- labels[valid]
  set.seed(666)
  sampled_cells <- unlist(lapply(split(names(labels), labels), function(cells) {
    label <- labels[[cells[[1]]]]
    maximum <- if (label == "epithelial") 6000L else 2000L
    sample(cells, min(length(cells), maximum))
  }), use.names = FALSE)
  counts <- round(get_counts(object)[, sampled_cells, drop = FALSE])
  cell_types <- factor(labels[sampled_cells])
  names(cell_types) <- sampled_cells
  list(
    reference = spacexr::Reference(
      counts[, names(cell_types), drop = FALSE],
      cell_types,
      Matrix::colSums(counts)
    ),
    label_column = label_col,
    sampled_counts = table(cell_types)
  )
}

copy_legacy_cache <- function(sample_name, legacy_dir, output_dir) {
  if (is.null(legacy_dir) || !dir.exists(legacy_dir)) return(FALSE)
  legacy_table <- file.path(
    legacy_dir, "tables",
    paste0("Auto_", sample_name, "_binned_rctd_annotations.csv.gz")
  )
  legacy_rds <- file.path(
    legacy_dir, "intermediate",
    paste0("Auto_", sample_name, "_binned_rctd.rds")
  )
  output_table <- file.path(
    output_dir, "tables",
    paste0("Auto_", sample_name, "_binned_rctd_annotations.csv.gz")
  )
  output_rds <- file.path(
    output_dir, "intermediate",
    paste0("Auto_", sample_name, "_binned_rctd.rds")
  )
  if (!file.exists(legacy_table)) return(FALSE)
  table_ok <- file.copy(legacy_table, output_table, overwrite = FALSE)
  if (file.exists(legacy_rds)) file.copy(legacy_rds, output_rds, overwrite = FALSE)
  isTRUE(table_ok) || file.exists(output_table)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))

output_dir <- args$output_dir %||%
  file.path(WD, "ref_outs", "visium_hd_outs", "rctd")
reference_path <- args$reference_path %||%
  file.path(WD, "ref_outs", "EAC_Ref_merged.rds")
legacy_cache_dir <- args$legacy_cache_dir %||%
  file.path(WD, "ref_outs", "visium_hd_outs", "legacy_visiumhd", "rctd")
min_umis <- as.integer(args$min_umis %||% 100L)
max_cores <- as.integer(args$max_cores %||% 8L)
force_rebuild <- tolower(args$force_rebuild %||% "false") %in% c("true", "1", "yes")


if (!file.exists(reference_path)) stop("Missing RCTD reference: ", reference_path)
if (!is.finite(min_umis) || min_umis < 1L) {
  stop("--min-umis must be a positive integer")
}
if (!is.finite(max_cores) || max_cores < 1L) {
  stop("--max-cores must be a positive integer")
}

manifest <- data.frame(
  sample = c("SUR1122", "SUR1231", "FFPEA1", "FFPED1"),
  binned_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"),
  segmented_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs"),
  stringsAsFactors = FALSE
)
required_columns <- c("sample", "binned_input")
if (!all(required_columns %in% colnames(manifest))) {
  stop("Manifest must contain: ", paste(required_columns, collapse = ", "))
}
if (anyDuplicated(manifest$sample)) stop("Sample names must be unique")
missing_inputs <- !dir.exists(manifest$binned_input)
if (any(missing_inputs)) {
  stop(
    "Missing binned inputs: ",
    paste(manifest$binned_input[missing_inputs], collapse = "; ")
  )
}

for (directory in file.path(output_dir, c("intermediate", "tables", "logs"))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

run_log <- c(
  paste0("start=", format(Sys.time(), tz = "Europe/London")),

  paste0("reference=", normalizePath(reference_path)),
  paste0("min_umis=", min_umis),
  paste0("max_cores=", max_cores)
)

for (sample_name in manifest$sample) {
  output_table <- file.path(
    output_dir, "tables",
    paste0("Auto_", sample_name, "_binned_rctd_annotations.csv.gz")
  )
  if (file.exists(output_table) || force_rebuild) next
  imported <- copy_legacy_cache(sample_name, legacy_cache_dir, output_dir)
  if (imported) {
    message("Imported archived RCTD cache: ", sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=legacy_cache"))
  }
}

pending <- if (force_rebuild) manifest$sample else manifest$sample[!file.exists(file.path(
  output_dir, "tables",
  paste0("Auto_", manifest$sample, "_binned_rctd_annotations.csv.gz")
))]

reference_info <- NULL
if (length(pending) > 0L) {
  reference_info <- build_reference(reference_path)
  write.csv(
    data.frame(
      celltype = names(reference_info$sampled_counts),
      n_cells = as.integer(reference_info$sampled_counts)
    ),
    file.path(output_dir, "tables", "Auto_visiumhd_rctd_reference_composition.csv"),
    row.names = FALSE
  )
}

summary_rows <- list()
for (i in seq_len(nrow(manifest))) {
  sample_name <- manifest$sample[[i]]
  input_dir <- normalizePath(manifest$binned_input[[i]], mustWork = TRUE)
  output_table <- file.path(
    output_dir, "tables",
    paste0("Auto_", sample_name, "_binned_rctd_annotations.csv.gz")
  )

  if (!file.exists(output_table) || force_rebuild) {
    message("Running RCTD: ", sample_name)
    counts_path <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
    if (!file.exists(counts_path)) stop("Missing count matrix: ", counts_path)
    counts <- read_10x_counts(counts_path)
    n_input <- ncol(counts)
    n_umis <- Matrix::colSums(counts)
    keep <- n_umis >= min_umis
    counts <- counts[, keep, drop = FALSE]
    n_umis <- n_umis[keep]
    spatial <- read_spatial_positions(file.path(input_dir, "spatial"), colnames(counts))
    counts <- counts[, spatial$barcodes, drop = FALSE]
    n_umis <- n_umis[spatial$barcodes]
    puck <- spacexr::SpatialRNA(
      spatial$coords[colnames(counts), , drop = FALSE],
      counts,
      n_umis
    )
    rctd <- spacexr::create.RCTD(
      puck,
      reference_info$reference,
      max_cores = max_cores,
      test_mode = FALSE,
      UMI_min = min_umis
    )
    rctd <- spacexr::run.RCTD(rctd, doublet_mode = "doublet")
    result <- as.data.frame(rctd@results$results_df)
    result$barcode <- rownames(result)
    result$sample <- sample_name
    result$total_counts <- as.numeric(n_umis[result$barcode])
    result$pxl_col_in_fullres <- spatial$coords[result$barcode, "x"]
    result$pxl_row_in_fullres <- spatial$coords[result$barcode, "y"]
    result$Auto_rctd_is_singlet <- as.character(result$spot_class) == "singlet"
    result <- result[, c(
      "barcode", "sample", "total_counts",
      "pxl_row_in_fullres", "pxl_col_in_fullres",
      setdiff(
        colnames(result),
        c(
          "barcode", "sample", "total_counts",
          "pxl_row_in_fullres", "pxl_col_in_fullres"
        )
      )
    )]
    saveRDS(
      rctd,
      file.path(
        output_dir, "intermediate",
        paste0("Auto_", sample_name, "_binned_rctd.rds")
      )
    )
    write_csv_gz(result, output_table)
    run_log <- c(
      run_log,
      paste0(
        "sample=", sample_name,
        "; status=calculated; n_input=", n_input,
        "; n_qc=", ncol(counts)
      )
    )
  }

  result <- read.csv(output_table, stringsAsFactors = FALSE)
  singlet <- if ("Auto_rctd_is_singlet" %in% colnames(result)) {
    as.logical(result$Auto_rctd_is_singlet)
  } else {
    as.character(result$spot_class) == "singlet"
  }
  class_counts <- sort(table(result$spot_class, useNA = "ifany"), decreasing = TRUE)
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    n_rctd_bins = nrow(result),
    n_singlets = sum(singlet, na.rm = TRUE),
    n_non_singlets = sum(!singlet, na.rm = TRUE),
    pct_singlets = 100 * mean(singlet, na.rm = TRUE),
    rctd_classes = paste(
      paste0(names(class_counts), "=", as.integer(class_counts)),
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(
  summary_df,
  file.path(output_dir, "tables", "Auto_visiumhd_binned_rctd_summary.csv"),
  row.names = FALSE
)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  summary_df,
  file.path(summary_dir, "visiumhd_final_rctd_summary.csv"),
  row.names = FALSE
)
run_log <- c(
  run_log,
  paste0("end=", format(Sys.time(), tz = "Europe/London")),
  paste0(
    "reference_label_column=",
    if (is.null(reference_info)) "not_loaded_all_cached" else reference_info$label_column
  )
)
writeLines(
  run_log,
  file.path(output_dir, "logs", "Auto_visiumhd_rctd_run_summary.txt")
)
####################
