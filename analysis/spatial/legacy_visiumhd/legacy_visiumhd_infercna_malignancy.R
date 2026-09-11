#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_infercna_malignancy.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: per-sample Visium HD count matrices and canonical annotation CSVs
#     from process_visium_hd.py / visiumhd_rctd_annotation.R.
#   Outputs: ref_outs/visium_hd_outs/malignancy/{intermediate,tables,figures,logs}/.
#   Cache/replot: saved per-sample infercna matrices support replotting without
#     rerunning CNA inference.
#   Run: Rscript analysis/spatial/visiumhd_infercna_malignancy.R --mode <mode>
#     --inputs <...> --sample-names <...> --annotation-dir <...> --output-dir <...>
#   Environment: dmtcp.
####################

####################
# Per-sample spatial InferCNA. Epithelial observations are targets; exactly two
# normal reference compartments are selected from endothelial/macrophage/
# fibroblast according to abundance after the annotation/doublet filter.
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(infercna)
  library(ggplot2)
  library(scales)
})

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

parse_cli <- function(args) {
  out <- list(inputs = character(), sample_names = character())
  multi <- c("inputs", "sample-names")
  allowed <- c(multi, "mode", "annotation-dir", "output-dir", "min-reference-cells", "min-epithelial-cells", "cancer-signature-path", "cancer-signature-threshold", "cna-sd-k")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% allowed) stop("Unknown argument: ", args[[i]])
    i <- i + 1L
    if (i > length(args) || startsWith(args[[i]], "--")) stop("Missing value for --", key)
    if (key %in% multi) {
      start <- i
      while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1L
      out[[gsub("-", "_", key)]] <- args[start:(i - 1L)]
    } else {
      out[[gsub("-", "_", key)]] <- args[[i]]
      i <- i + 1L
    }
  }
  out
}

read_10x_counts <- function(path) {
  counts <- Read10X_h5(path)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }
  Matrix::Matrix(counts, sparse = TRUE)
}

to_cpm <- function(counts) {
  library_size <- Matrix::colSums(counts)
  library_size[!is.finite(library_size) | library_size <= 0] <- 1
  cpm <- Matrix::t(Matrix::t(counts) * (1e6 / library_size))
  dimnames(cpm) <- dimnames(counts)
  cpm
}

as_logical <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

nearest_binned_malignancy <- function(annotation, scatter, output_dir, sample_name) {
  result <- list(
    malignant = rep(NA, nrow(scatter)),
    barcode = rep(NA_character_, nrow(scatter)),
    distance = rep(NA_real_, nrow(scatter))
  )
  binned_path <- file.path(output_dir, "tables", paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz"))
  if (!file.exists(binned_path)) return(result)
  binned <- read.csv(binned_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres", "is_epithelial_target", "Auto_malignant")
  if (length(setdiff(required, colnames(binned)))) return(result)
  binned$Auto_malignant <- as_logical(binned$Auto_malignant)
  binned$is_epithelial_target <- as_logical(binned$is_epithelial_target)
  binned$pxl_col_in_fullres <- suppressWarnings(as.numeric(binned$pxl_col_in_fullres))
  binned$pxl_row_in_fullres <- suppressWarnings(as.numeric(binned$pxl_row_in_fullres))
  binned <- binned[
    binned$is_epithelial_target & is.finite(binned$pxl_col_in_fullres) & is.finite(binned$pxl_row_in_fullres),
    , drop = FALSE
  ]
  annotation <- annotation[match(scatter$barcode, annotation$barcode), , drop = FALSE]
  annotation$pxl_col_in_fullres <- suppressWarnings(as.numeric(annotation$pxl_col_in_fullres))
  annotation$pxl_row_in_fullres <- suppressWarnings(as.numeric(annotation$pxl_row_in_fullres))
  target <- scatter$is_epithelial_target & is.finite(annotation$pxl_col_in_fullres) & is.finite(annotation$pxl_row_in_fullres)
  if (nrow(binned) == 0L || !any(target)) return(result)
  nearest <- RANN::nn2(
    data = as.matrix(binned[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    query = as.matrix(annotation[target, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    k = 1L
  )
  target_idx <- which(target)
  result$malignant[target_idx] <- binned$Auto_malignant[nearest$nn.idx[, 1L]]
  result$barcode[target_idx] <- binned$barcode[nearest$nn.idx[, 1L]]
  result$distance[target_idx] <- nearest$nn.dists[, 1L]
  result
}

segmented_keratinocyte_projection <- function(annotation, scatter, annotation_dir, sample_name) {
  result <- data.frame(
    Auto_segmented_keratinocyte_n_cells = integer(nrow(scatter)),
    Auto_segmented_keratinocyte_fraction = numeric(nrow(scatter)),
    Auto_segmented_keratinocyte_nearest_distance = rep(NA_real_, nrow(scatter)),
    stringsAsFactors = FALSE
  )
  segmented_path <- file.path(annotation_dir, paste0("Auto_", sample_name, "_segmented_cell_annotations.csv.gz"))
  if (!file.exists(segmented_path)) return(result)
  segmented <- read.csv(segmented_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Auto_annotation_celltype", "pxl_col_in_fullres", "pxl_row_in_fullres")
  if (length(setdiff(required, colnames(segmented)))) return(result)
  annotation <- annotation[match(scatter$barcode, annotation$barcode), , drop = FALSE]
  annotation$pxl_col_in_fullres <- suppressWarnings(as.numeric(annotation$pxl_col_in_fullres))
  annotation$pxl_row_in_fullres <- suppressWarnings(as.numeric(annotation$pxl_row_in_fullres))
  segmented$pxl_col_in_fullres <- suppressWarnings(as.numeric(segmented$pxl_col_in_fullres))
  segmented$pxl_row_in_fullres <- suppressWarnings(as.numeric(segmented$pxl_row_in_fullres))
  bin_idx <- which(scatter$is_epithelial_target & is.finite(annotation$pxl_col_in_fullres) & is.finite(annotation$pxl_row_in_fullres))
  cell_idx <- which(is.finite(segmented$pxl_col_in_fullres) & is.finite(segmented$pxl_row_in_fullres))
  if (!length(bin_idx) || !length(cell_idx)) return(result)
  nearest <- RANN::nn2(
    data = as.matrix(annotation[bin_idx, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    query = as.matrix(segmented[cell_idx, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    k = 1L
  )
  matched <- data.frame(
    scatter_idx = bin_idx[nearest$nn.idx[, 1L]],
    distance = nearest$nn.dists[, 1L],
    keratinocyte = segmented$Auto_annotation_celltype[cell_idx] == "keratinocyte"
  )
  matched <- matched[is.finite(matched$distance) & matched$distance <= 25, , drop = FALSE]
  if (!nrow(matched)) return(result)
  composition <- aggregate(cbind(n_cells = rep(1L, nrow(matched)), n_keratinocyte = as.integer(matched$keratinocyte)), by = list(scatter_idx = matched$scatter_idx), FUN = sum)
  result$Auto_segmented_keratinocyte_n_cells[composition$scatter_idx] <- composition$n_cells
  result$Auto_segmented_keratinocyte_fraction[composition$scatter_idx] <- composition$n_keratinocyte / composition$n_cells
  nearest_distance <- aggregate(distance ~ scatter_idx, data = matched, FUN = median)
  result$Auto_segmented_keratinocyte_nearest_distance[nearest_distance$scatter_idx] <- nearest_distance$distance
  result
}

write_csv_gz <- function(data, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(data, con, row.names = FALSE)
}

plot_scatter <- function(scatter_df, threshold_signal, threshold_cor, cna_sd_k, sample_name, mode, output_dir) {
  scatter_df$plot_group <- ifelse(scatter_df$is_reference, "Reference", "Non-malignant epithelial")
  scatter_df$plot_group[scatter_df$Auto_malignancy == "malignant_level_1"] <- "Malignant level 1 (CNA)"
  scatter_df$plot_group[scatter_df$Auto_malignancy == "malignant_level_2"] <- "Malignant level 2 (CNA unresolved + signature)"
  scatter_df$plot_group <- factor(
    scatter_df$plot_group,
    levels = c("Reference", "Non-malignant epithelial", "Malignant level 1 (CNA)", "Malignant level 2 (CNA unresolved + signature)")
  )
  p <- ggplot(scatter_df, aes(x = cna.signal, y = cna.cor, colour = plot_group)) +
    geom_point(size = 0.75, alpha = 0.65) +
    geom_vline(xintercept = threshold_signal, linetype = "dashed", linewidth = 0.45, colour = "grey30") +
    geom_hline(yintercept = threshold_cor, linetype = "dashed", linewidth = 0.45, colour = "grey30") +
    scale_colour_manual(values = c(
      "Reference" = "#9E9E9E",
      "Non-malignant epithelial" = "#4DAF4A",
      "Malignant level 1 (CNA)" = "#E41A1C",
      "Malignant level 2 (CNA unresolved + signature)" = "#984EA3"
    ), drop = FALSE) +
    scale_x_continuous(labels = scales::scientific) +
    labs(
      title = paste0(sample_name, " ", mode, ": InferCNA malignancy classification"),
      subtitle = paste0("Level 1: both CNA metrics > reference mean + ", cna_sd_k, " SD; level 2: CNA unresolved plus malignant signature"),
      x = "CNA signal",
      y = "CNA correlation",
      colour = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  prefix <- file.path(output_dir, "figures", paste0("Auto_", sample_name, "_", mode, "_infercna_scatter"))
  ggsave(paste0(prefix, ".pdf"), p, width = 8.5, height = 6.5, useDingbats = FALSE)
  ggsave(paste0(prefix, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (length(args$inputs) == 0L || length(args$sample_names) == 0L || is.null(args$mode) || is.null(args$annotation_dir) || is.null(args$output_dir)) {
  stop("Required: --mode --inputs --sample-names --annotation-dir --output-dir")
}
####################
# Custom uses the same 16 um count matrix as binned mode, but its annotation
# table has RCTD-filtered doublets and manual marker calls for singlets.
####################
####################
# The spatial method uses segmented counts with a different annotation table.
if (!args$mode %in% c("binned", "segmented", "custom", "spatial")) {
  stop("--mode must be binned, segmented, custom, or spatial")
}
####################
####################
if (length(args$inputs) != length(args$sample_names)) stop("--inputs and --sample-names must have equal length")
####################
# Twenty cells per distinct normal type is sufficient to retain the RCTD
# singlet-gated SUR1231 macrophage compartment (n=21) while still rejecting
# one-compartment reference sets.
####################
min_reference_cells <- if (!is.null(args$min_reference_cells)) as.integer(args$min_reference_cells) else 20L
min_epithelial_cells <- if (!is.null(args$min_epithelial_cells)) as.integer(args$min_epithelial_cells) else 30L
cancer_signature_path <- if (!is.null(args$cancer_signature_path)) args$cancer_signature_path else file.path(WD, "ref_outs", "cancer_signatures.txt")
cancer_signature_threshold <- if (!is.null(args$cancer_signature_threshold)) as.numeric(args$cancer_signature_threshold) else 1
cna_sd_k <- if (!is.null(args$cna_sd_k)) as.numeric(args$cna_sd_k) else 1
if (!is.finite(min_reference_cells) || min_reference_cells < 1L) stop("--min-reference-cells must be positive")
if (!is.finite(min_epithelial_cells) || min_epithelial_cells < 1L) stop("--min-epithelial-cells must be positive")
if (!file.exists(cancer_signature_path)) stop("Missing cancer signature file: ", cancer_signature_path)
if (!is.finite(cancer_signature_threshold)) stop("--cancer-signature-threshold must be finite")
if (!is.finite(cna_sd_k) || cna_sd_k <= 0) stop("--cna-sd-k must be positive")
cancer_signature_genes <- unique(read.table(cancer_signature_path, header = FALSE, stringsAsFactors = FALSE)[[1]])

annotation_dir <- normalizePath(args$annotation_dir, mustWork = TRUE)
output_dir <- args$output_dir
for (subdir in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)
gene_order <- read.table(gene_order_path, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene", "chromosome", "start", "end"))
normal_types <- c("endothelial", "macrophage", "fibroblast")
summary_rows <- list()
run_log <- c(paste0("start=", format(Sys.time(), tz = "Europe/London")), paste0("mode=", args$mode), paste0("min_reference_cells=", min_reference_cells), paste0("min_epithelial_cells=", min_epithelial_cells), paste0("cancer_signature_path=", normalizePath(cancer_signature_path)), paste0("cancer_signature_threshold=", cancer_signature_threshold), paste0("cna_sd_k=", cna_sd_k))

for (i in seq_along(args$inputs)) {
  sample_name <- args$sample_names[[i]]
  input_dir <- normalizePath(args$inputs[[i]], mustWork = TRUE)
  annotation_path <- file.path(annotation_dir, paste0("Auto_", sample_name, "_", args$mode, "_cell_annotations.csv.gz"))
  if (!file.exists(annotation_path)) stop("Missing annotation table: ", annotation_path)
  annotation <- read.csv(annotation_path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("barcode", "Auto_annotation_celltype", "Auto_annotation_pass_doublet_filter", "Auto_annotation_keep_epithelial")
  missing <- setdiff(needed, colnames(annotation))
  if (length(missing)) stop("Annotation table missing columns: ", paste(missing, collapse = ", "))
  annotation$barcode <- as.character(annotation$barcode)
  annotation$Auto_annotation_pass_doublet_filter <- as_logical(annotation$Auto_annotation_pass_doublet_filter)
  annotation$Auto_annotation_keep_epithelial <- as_logical(annotation$Auto_annotation_keep_epithelial)
  annotation$Auto_annotation_celltype <- as.character(annotation$Auto_annotation_celltype)
  ####################
  # The two 16 um annotation methods must be compared on one CNA coordinate
  # system. Use the RCTD binned singlet universe and its normal compartments for
  # both binned modes; otherwise annotation-dependent references shift identical
  # expression profiles on the InferCNA scatter. Segmented cells retain their
  # independently measured and independently annotated normal reference.
  ####################
  reference_annotation <- annotation
  reference_annotation_mode <- args$mode
  if (args$mode == "custom") {
    reference_path <- file.path(annotation_dir, paste0("Auto_", sample_name, "_binned_cell_annotations.csv.gz"))
    if (!file.exists(reference_path)) stop("Missing common binned reference annotation: ", reference_path)
    reference_annotation <- read.csv(reference_path, stringsAsFactors = FALSE, check.names = FALSE)
    reference_missing <- setdiff(needed, colnames(reference_annotation))
    if (length(reference_missing)) stop("Binned reference annotation missing columns: ", paste(reference_missing, collapse = ", "))
    reference_annotation$barcode <- as.character(reference_annotation$barcode)
    reference_annotation$Auto_annotation_pass_doublet_filter <- as_logical(reference_annotation$Auto_annotation_pass_doublet_filter)
    reference_annotation$Auto_annotation_keep_epithelial <- as_logical(reference_annotation$Auto_annotation_keep_epithelial)
    reference_annotation$Auto_annotation_celltype <- as.character(reference_annotation$Auto_annotation_celltype)
    reference_annotation_mode <- "binned_rctd_common"
  }
  ####################
  usable <- annotation$Auto_annotation_pass_doublet_filter
  reference_usable <- reference_annotation$Auto_annotation_pass_doublet_filter
  annotation_type_counts <- table(annotation$Auto_annotation_celltype[usable])
  annotation_normal_counts <- setNames(as.integer(annotation_type_counts[normal_types]), normal_types)
  annotation_normal_counts[is.na(annotation_normal_counts)] <- 0L
  type_counts <- table(reference_annotation$Auto_annotation_celltype[reference_usable])
  normal_counts <- setNames(as.integer(type_counts[normal_types]), normal_types)
  normal_counts[is.na(normal_counts)] <- 0L
  epithelial_barcodes <- annotation$barcode[usable & annotation$Auto_annotation_keep_epithelial]
  eligible_types <- names(normal_counts)[normal_counts >= min_reference_cells]
  selected_types <- names(sort(normal_counts[eligible_types], decreasing = TRUE))[seq_len(min(2L, length(eligible_types)))]
  base_summary <- data.frame(
    sample = sample_name,
    mode = args$mode,
    n_annotated = nrow(annotation),
    n_pass_doublet_filter = sum(usable),
    n_epithelial = length(epithelial_barcodes),
    n_endothelial = annotation_normal_counts[["endothelial"]],
    n_macrophage = annotation_normal_counts[["macrophage"]],
    n_fibroblast = annotation_normal_counts[["fibroblast"]],
    n_reference_endothelial = normal_counts[["endothelial"]],
    n_reference_macrophage = normal_counts[["macrophage"]],
    n_reference_fibroblast = normal_counts[["fibroblast"]],
    selected_reference_types = paste(selected_types, collapse = ";"),
    reference_annotation_mode = reference_annotation_mode,
    stringsAsFactors = FALSE
  )
  if (length(selected_types) < 2L) {
    base_summary$status <- "insufficient_reference_compartments"
    base_summary$n_reference_cells <- sum(normal_counts[selected_types])
    base_summary$n_malignant <- NA_integer_
    base_summary$pct_malignant_epithelial <- NA_real_
    summary_rows[[sample_name]] <- base_summary
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=insufficient_reference_compartments"))
    next
  }
  if (length(epithelial_barcodes) < min_epithelial_cells) {
    base_summary$status <- "insufficient_epithelial_targets"
    base_summary$n_reference_cells <- sum(normal_counts[selected_types])
    base_summary$n_malignant <- NA_integer_
    base_summary$pct_malignant_epithelial <- NA_real_
    summary_rows[[sample_name]] <- base_summary
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=insufficient_epithelial_targets"))
    next
  }
  counts_path <- file.path(input_dir, if (args$mode %in% c("binned", "custom")) "filtered_feature_bc_matrix.h5" else "filtered_feature_cell_matrix.h5")
  if (!file.exists(counts_path)) stop("Missing count matrix: ", counts_path)
  counts <- read_10x_counts(counts_path)
  shared_barcodes <- intersect(colnames(counts), annotation$barcode)
  if (length(shared_barcodes) == 0L) stop("No count barcodes matched annotation for ", sample_name)
  annotation <- annotation[match(shared_barcodes, annotation$barcode), , drop = FALSE]
  usable <- annotation$Auto_annotation_pass_doublet_filter
  reference_shared <- intersect(shared_barcodes, reference_annotation$barcode)
  reference_annotation <- reference_annotation[match(reference_shared, reference_annotation$barcode), , drop = FALSE]
  reference_usable <- reference_annotation$Auto_annotation_pass_doublet_filter
  ref_barcodes <- lapply(selected_types, function(type) {
    reference_annotation$barcode[reference_usable & reference_annotation$Auto_annotation_celltype == type]
  })
  names(ref_barcodes) <- selected_types
  if (any(lengths(ref_barcodes) < min_reference_cells)) stop("Reference barcode matching failed for ", sample_name)
  reference_barcodes <- unlist(ref_barcodes, use.names = FALSE)
  epithelial_barcodes <- setdiff(intersect(epithelial_barcodes, shared_barcodes), reference_barcodes)
  if (args$mode %in% c("binned", "custom")) {
    use_barcodes <- reference_annotation$barcode[reference_usable]
  } else {
    use_barcodes <- unique(c(epithelial_barcodes, reference_barcodes))
  }
  base_summary$n_epithelial_assessed <- length(epithelial_barcodes)
  base_summary$n_common_binned_infercna_cells <- if (args$mode %in% c("binned", "custom")) length(use_barcodes) else NA_integer_
  ####################
  counts <- counts[, use_barcodes, drop = FALSE]
  keep_genes <- intersect(rownames(counts), gene_order$gene)
  keep_genes <- keep_genes[!duplicated(keep_genes)]
  if (length(keep_genes) < 5000L) {
    base_summary$status <- "too_few_genome_ordered_genes"
    base_summary$n_reference_cells <- sum(lengths(ref_barcodes))
    base_summary$n_malignant <- NA_integer_
    base_summary$pct_malignant_epithelial <- NA_real_
    summary_rows[[sample_name]] <- base_summary
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=too_few_genome_ordered_genes"))
    next
  }
  counts <- counts[keep_genes, , drop = FALSE]
  cpm <- to_cpm(counts)
  message("InferCNA: ", sample_name, " (targets=", length(epithelial_barcodes), ", references=", sum(lengths(ref_barcodes)), ")")
  ####################
  # Custom and RCTD binned modes deliberately share the complete InferCNA input
  # and reference. Reuse the finalized binned matrix for custom after strict
  # dimension-name validation, then recompute only annotation-dependent gates.
  ####################
  common_binned_cache <- file.path(output_dir, "intermediate", paste0("Auto_", sample_name, "_binned_infercna_outs.rds"))
  common_binned_cell_cache <- file.path(output_dir, "tables", paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz"))
  reused_common_binned_cache <- FALSE
  attempt_binned_cache <- args$mode == "custom" && file.exists(common_binned_cache)
  if (args$mode == "binned" && file.exists(common_binned_cache) && file.exists(common_binned_cell_cache)) {
    cached_cells <- read.csv(common_binned_cell_cache, stringsAsFactors = FALSE, check.names = FALSE)
    cached_finite_barcodes <- cached_cells$barcode[is.finite(cached_cells$cna.signal) & is.finite(cached_cells$cna.cor)]
    attempt_binned_cache <- setequal(cached_finite_barcodes, colnames(cpm))
    rm(cached_cells)
  }
  if (attempt_binned_cache) {
    cached_outs <- readRDS(common_binned_cache)
    cache_matches <- !is.list(cached_outs) &&
      nrow(cached_outs) >= 5000L &&
      setequal(colnames(cached_outs), colnames(cpm))
    if (args$mode == "custom" && !cache_matches) stop("Common binned InferCNA cache does not match custom input for ", sample_name)
    if (cache_matches) {
      outs <- cached_outs[, colnames(cpm), drop = FALSE]
      reused_common_binned_cache <- TRUE
    } else {
      rm(cached_outs)
      gc()
    }
  }
  if (!reused_common_binned_cache) {
    outs <- tryCatch(
      infercna::infercna(as.matrix(cpm), refCells = ref_barcodes, isLog = FALSE, verbose = TRUE),
      error = function(error) error
    )
  }
  base_summary$reused_common_binned_cna <- reused_common_binned_cache
  ####################
  if (inherits(outs, "error")) {
    base_summary$status <- "infercna_error"
    base_summary$n_reference_cells <- sum(lengths(ref_barcodes))
    base_summary$n_malignant <- NA_integer_
    base_summary$pct_malignant_epithelial <- NA_real_
    base_summary$error_message <- conditionMessage(outs)
    summary_rows[[sample_name]] <- base_summary
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=infercna_error; message=", conditionMessage(outs)))
    next
  }
  scatter <- as.data.frame(infercna::cnaScatterPlot(outs, excludeFromAvg = unlist(ref_barcodes, use.names = FALSE)))
  scatter$barcode <- rownames(scatter)
  scatter$is_reference <- scatter$barcode %in% unlist(ref_barcodes, use.names = FALSE)
  scatter$is_epithelial_target <- scatter$barcode %in% epithelial_barcodes
  ref_scatter <- scatter[scatter$is_reference, , drop = FALSE]
  ####################
  # Per-cell segmented data have broader reference scatter than 16 um bins.
  # A one-SD, two-metric gate preserves both independent CNA requirements and
  # avoids treating reference technical variance as a tumour-cell cutoff.
  ####################
  threshold_signal <- mean(ref_scatter$cna.signal, na.rm = TRUE) + cna_sd_k * sd(ref_scatter$cna.signal, na.rm = TRUE)
  threshold_cor <- mean(ref_scatter$cna.cor, na.rm = TRUE) + cna_sd_k * sd(ref_scatter$cna.cor, na.rm = TRUE)
  if (!is.finite(threshold_signal)) threshold_signal <- max(ref_scatter$cna.signal, na.rm = TRUE)
  if (!is.finite(threshold_cor)) threshold_cor <- max(ref_scatter$cna.cor, na.rm = TRUE)
  ####################
  # Keep the two-metric CNA tier, then rescue only CNA-unresolved epithelial
  # observations using the scATLAS cancer-signature procedure in Malignancy.R.
  ####################
  scatter$Auto_cna_class <- "not_assessed"
  scatter$Auto_cna_class[scatter$is_reference] <- "reference"
  target_idx <- scatter$is_epithelial_target
  both_cna <- scatter$cna.signal > threshold_signal & scatter$cna.cor > threshold_cor
  either_cna <- scatter$cna.signal > threshold_signal | scatter$cna.cor > threshold_cor
  scatter$Auto_cna_class[target_idx & both_cna] <- "cna_malignant"
  scatter$Auto_cna_class[target_idx & !both_cna & either_cna] <- "cna_unresolved"
  scatter$Auto_cna_class[target_idx & !either_cna] <- "cna_non_malignant"
  signature_genes_present <- intersect(cancer_signature_genes, rownames(cpm))
  target_in_cpm <- intersect(epithelial_barcodes, colnames(cpm))
  if (length(signature_genes_present) == 0L || length(target_in_cpm) == 0L) {
    stop("No cancer-signature genes or epithelial targets available after genome ordering for ", sample_name)
  }
  signature_gene_means <- Matrix::rowMeans(log1p(cpm[signature_genes_present, target_in_cpm, drop = FALSE] / 100), na.rm = TRUE)
  signature_genes_selected <- names(sort(signature_gene_means, decreasing = TRUE))[seq_len(min(50L, length(signature_gene_means)))]
  signature_scores <- Matrix::colMeans(log1p(cpm[signature_genes_selected, , drop = FALSE] / 100), na.rm = TRUE)
  scatter$Auto_cancer_signature_score <- unname(signature_scores[match(scatter$barcode, names(signature_scores))])
  scatter$Auto_cancer_signature_status <- "not_assessed"
  scatter$Auto_cancer_signature_status[target_idx] <- ifelse(
    scatter$Auto_cancer_signature_score[target_idx] >= cancer_signature_threshold,
    "cs_malignant", "cs_unresolved"
  )
  ####################
  # Segmented observations have lower CNA signal-to-noise than their matched
  # 16 um bins. A CNA-non-malignant segmented cell is rescued only when both
  # its cancer signature is positive and its nearest binned epithelial call
  # is malignant; signature alone remains insufficient.
  ####################
  scatter$Auto_binned_malignant <- NA
  scatter$Auto_nearest_binned_barcode <- NA_character_
  scatter$Auto_nearest_binned_distance <- NA_real_
  ####################
  # Spatial and custom segmented annotations use the same cell matrix and matched binned rescue.
  if (args$mode %in% c("segmented", "spatial")) {
    binned_projection <- nearest_binned_malignancy(annotation, scatter, output_dir, sample_name)
    scatter$Auto_binned_malignant <- binned_projection$malignant
    scatter$Auto_nearest_binned_barcode <- binned_projection$barcode
    scatter$Auto_nearest_binned_distance <- binned_projection$distance
  }
  ####################
  scatter$Auto_segmented_keratinocyte_n_cells <- 0L
  scatter$Auto_segmented_keratinocyte_fraction <- 0
  scatter$Auto_segmented_keratinocyte_nearest_distance <- NA_real_
  if (args$mode == "binned") {
    keratinocyte_projection <- segmented_keratinocyte_projection(annotation, scatter, annotation_dir, sample_name)
    scatter[, colnames(keratinocyte_projection)] <- keratinocyte_projection
  }
  scatter$Auto_malignancy <- "unresolved"
  scatter$Auto_malignancy_evidence <- "none"
  scatter$Auto_malignancy[target_idx & scatter$Auto_cna_class == "cna_malignant"] <- "malignant_level_1"
  scatter$Auto_malignancy_evidence[target_idx & scatter$Auto_cna_class == "cna_malignant"] <- "two_metric_cna"
  scatter$Auto_malignancy[target_idx & scatter$Auto_cna_class == "cna_unresolved" & scatter$Auto_cancer_signature_status == "cs_malignant"] <- "malignant_level_2"
  scatter$Auto_malignancy_evidence[target_idx & scatter$Auto_cna_class == "cna_unresolved" & scatter$Auto_cancer_signature_status == "cs_malignant"] <- "cna_unresolved_plus_signature"
  binned_signature_rescue <- target_idx & scatter$Auto_cna_class == "cna_non_malignant" & scatter$Auto_cancer_signature_status == "cs_malignant" & scatter$Auto_binned_malignant %in% TRUE
  scatter$Auto_malignancy[binned_signature_rescue] <- "malignant_level_2"
  scatter$Auto_malignancy_evidence[binned_signature_rescue] <- "binned_cna_plus_signature"
  ####################
  # RCTD has no keratinocyte reference and therefore retains these bins as
  # epithelial. Keep that annotation intact, but use the matched segmented
  # keratinocyte composition to classify a keratinocyte-dominant bin as normal
  # in the malignancy layer before state mapping.
  ####################
  normal_keratinocyte <- args$mode == "binned" & target_idx & scatter$Auto_segmented_keratinocyte_n_cells >= 1L & scatter$Auto_segmented_keratinocyte_fraction >= 0.5
  scatter$Auto_malignancy_before_keratinocyte_exclusion <- scatter$Auto_malignancy
  scatter$Auto_malignancy[normal_keratinocyte] <- "normal_keratinocyte"
  scatter$Auto_malignancy_evidence[normal_keratinocyte] <- "segmented_keratinocyte_dominant"
  scatter$Auto_malignant <- scatter$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")
  cell_table <- merge(annotation, scatter[, c("barcode", "cna.signal", "cna.cor", "is_reference", "is_epithelial_target", "Auto_cna_class", "Auto_cancer_signature_score", "Auto_cancer_signature_status", "Auto_binned_malignant", "Auto_nearest_binned_barcode", "Auto_nearest_binned_distance", "Auto_segmented_keratinocyte_n_cells", "Auto_segmented_keratinocyte_fraction", "Auto_segmented_keratinocyte_nearest_distance", "Auto_malignancy_before_keratinocyte_exclusion", "Auto_malignancy", "Auto_malignancy_evidence", "Auto_malignant")], by = "barcode", all.x = TRUE, sort = FALSE)
  cell_table$sample <- sample_name
  cell_table$mode <- args$mode
  cell_table$Auto_cna_signal_threshold <- threshold_signal
  cell_table$Auto_cna_cor_threshold <- threshold_cor
  cell_table$Auto_cna_sd_k <- cna_sd_k
  cell_table$Auto_cancer_signature_threshold <- cancer_signature_threshold
  cell_table$Auto_cancer_signature_n_genes <- length(signature_genes_selected)
  write_csv_gz(cell_table, file.path(output_dir, "tables", paste0("Auto_", sample_name, "_", args$mode, "_infercna_cells.csv.gz")))
  saveRDS(outs, file.path(output_dir, "intermediate", paste0("Auto_", sample_name, "_", args$mode, "_infercna_outs.rds")))
  plot_scatter(scatter, threshold_signal, threshold_cor, cna_sd_k, sample_name, args$mode, output_dir)
  target_scatter <- scatter[scatter$is_epithelial_target, , drop = FALSE]
  base_summary$status <- "complete"
  base_summary$n_reference_cells <- sum(lengths(ref_barcodes))
  base_summary$n_genome_ordered_genes <- length(keep_genes)
  base_summary$cna_signal_threshold <- threshold_signal
  base_summary$cna_cor_threshold <- threshold_cor
  base_summary$cna_sd_k <- cna_sd_k
  base_summary$cancer_signature_threshold <- cancer_signature_threshold
  base_summary$n_cancer_signature_genes <- length(signature_genes_selected)
  base_summary$n_signature_malignant_epithelial <- sum(target_scatter$Auto_cancer_signature_status == "cs_malignant", na.rm = TRUE)
  base_summary$n_signature_malignant_cna_malignant <- sum(target_scatter$Auto_cna_class == "cna_malignant" & target_scatter$Auto_cancer_signature_status == "cs_malignant", na.rm = TRUE)
  base_summary$n_signature_malignant_cna_unresolved <- sum(target_scatter$Auto_cna_class == "cna_unresolved" & target_scatter$Auto_cancer_signature_status == "cs_malignant", na.rm = TRUE)
  base_summary$n_signature_malignant_cna_non_malignant <- sum(target_scatter$Auto_cna_class == "cna_non_malignant" & target_scatter$Auto_cancer_signature_status == "cs_malignant", na.rm = TRUE)
  base_summary$n_segmented_binned_malignant_concordant <- sum(target_scatter$Auto_binned_malignant %in% TRUE, na.rm = TRUE)
  base_summary$n_malignant_level_2_from_binned_signature <- sum(target_scatter$Auto_malignancy_evidence == "binned_cna_plus_signature", na.rm = TRUE)
  base_summary$n_normal_keratinocyte <- sum(target_scatter$Auto_malignancy == "normal_keratinocyte", na.rm = TRUE)
  base_summary$n_keratinocyte_normal_from_level_1 <- sum(target_scatter$Auto_malignancy == "normal_keratinocyte" & target_scatter$Auto_malignancy_before_keratinocyte_exclusion == "malignant_level_1", na.rm = TRUE)
  base_summary$n_malignant_level_1 <- sum(target_scatter$Auto_malignancy == "malignant_level_1", na.rm = TRUE)
  base_summary$n_cna_unresolved <- sum(target_scatter$Auto_cna_class == "cna_unresolved", na.rm = TRUE)
  base_summary$n_malignant_level_2 <- sum(target_scatter$Auto_malignancy == "malignant_level_2", na.rm = TRUE)
  base_summary$n_malignant <- sum(target_scatter$Auto_malignant, na.rm = TRUE)
  base_summary$pct_malignant_epithelial <- 100 * mean(target_scatter$Auto_malignant, na.rm = TRUE)
  summary_rows[[sample_name]] <- base_summary
  run_log <- c(run_log, paste0("sample=", sample_name, "; status=complete; malignant=", base_summary$n_malignant))
}

summary_df <- as.data.frame(dplyr::bind_rows(summary_rows))
write.csv(summary_df, file.path(output_dir, "tables", paste0("Auto_visiumhd_", args$mode, "_infercna_malignancy_summary.csv")), row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(summary_dir, paste0("visiumhd_", args$mode, "_infercna_malignancy_summary.csv")), row.names = FALSE)
run_log <- c(run_log, paste0("end=", format(Sys.time(), tz = "Europe/London")))
writeLines(run_log, file.path(output_dir, "logs", paste0("Auto_visiumhd_", args$mode, "_infercna_run_summary.txt")))

failed_samples <- summary_df$sample[summary_df$status != "complete"]
if (length(failed_samples)) {
  stop("CNA classification did not complete for: ", paste(failed_samples, collapse = ", "), ". See the per-sample summary CSV before state mapping.")
}
####################
