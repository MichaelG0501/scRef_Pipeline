#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: terminal
#   Script: analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R
#   Description: sample-independent and pooled spatial-neighbourhood tests of
#     malignant centred-refined states/MPs against whole non-malignant cell
#     types and silhouette-retained non-malignant GeneNMF MPs in filtered
#     Visium HD 16 um bins from SUR1231, FFPEA1, and FFPED1 (D1).
#   Methodology:
#     analysis/methodology/spatial/visium_hd_spatial_cancer_tme_interactions_methodology.md
#   Inputs:
#     analysis/spatial/visium_hd_samples.tsv
#     analysis/shared/visium_hd_celltype_colours.tsv
#     ref_outs/visium_hd_outs/malignancy/tables/Auto_<sample>_binned_malignancy.csv.gz
#     ref_outs/visium_hd_outs/state_mapping/tables/Auto_visiumhd_binned_malignant_state_annotations.csv.gz
#     ref_outs/nmf_{fibroblast,endothelial,macrophage,nk,plasma,cd4,cd8}/MP_outs_default.rds
#     ref_outs/non_malignant_mp_correlations/{01,03,05,06}_*/Auto_celltype_correlations_all.csv
#     Space Ranger square_016um filtered count matrices from the manifest
#   Outputs:
#     ref_outs/visium_hd_outs/spatial_interactions/intermediate/: persistent
#       TME MP scores, gene sets, and replot-ready per-sample memberships
#     ref_outs/visium_hd_outs/spatial_interactions/tables/: sample-specific,
#       pooled, recurrence, sensitivity, scRNA-validation, and audit tables
#     ref_outs/visium_hd_outs/spatial_interactions/figures/: sample-side-by-side
#       and pooled dot maps, robustness plots, and spatial interaction maps
#     ref_outs/visium_hd_outs/spatial_interactions/logs/: run/session summary
#     updates/new_updates/summaries/visium_hd_spatial_cancer_tme_interactions_summary.csv
#   Cache/replot: SCREF_FORCE_REBUILD=TRUE rebuilds memberships and tests;
#     SCREF_REPLOT_ONLY=TRUE rebuilds figures from persistent live inputs.
#   Run: Rscript analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R
#   Environment: dmtcp
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(patchwork)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)

samples <- c("SUR1231", "FFPEA1", "FFPED1")
sample_display <- c("SUR1231" = "SUR1231", "FFPEA1" = "FFPEA1", "FFPED1" = "D1 (FFPED1)")
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "spatial_interactions")
output_tiers <- setNames(file.path(output_dir, c("intermediate", "tables", "figures", "logs", "reports")),
                         c("intermediate", "tables", "figures", "logs", "reports"))
for (path in output_tiers) dir.create(path, recursive = TRUE, showWarnings = FALSE)
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)


malignancy_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "malignancy", "tables")
state_path <- file.path(
  wd, "ref_outs", "visium_hd_outs", "state_mapping", "tables",
  "Auto_visiumhd_binned_malignant_state_annotations.csv.gz"
)

cache_version <- "2026-08-09_v1_spatial_cancer_tme"
force_rebuild <- tolower(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")) %in% c("true", "t", "1", "yes")
replot_only <- tolower(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")) %in% c("true", "t", "1", "yes")
n_permutations <- as.integer(Sys.getenv("SCREF_SPATIAL_PERMUTATIONS", "499"))
n_cores <- max(1L, as.integer(Sys.getenv(
  "SCREF_NCORES",
  if (nzchar(Sys.getenv("PBS_JOBID", ""))) Sys.getenv("NCPUS", "8") else "1"
)))
if (!is.finite(n_permutations) || n_permutations < 99L) stop("SCREF_SPATIAL_PERMUTATIONS must be >= 99")

neighbour_rings <- c(1L, 2L, 3L)
primary_ring <- 1L
mp_positive_threshold <- 0.5
top_genes_per_mp <- 100L
min_source_bins <- 20L
min_whole_target_bins <- 20L
min_parent_bins_for_mp <- 20L
min_mp_positive_bins <- 10L
min_observed_neighbour_slots <- 20L
sample_fdr_cutoff <- 0.05
pooled_fdr_cutoff <- 0.05
recurrence_fdr_cutoff <- 0.10
spatial_pages_per_direction <- 3L
maximum_plotted_contact_edges <- 2500L

state_crosswalk <- c(
  "Classic proliferation" = "Classic proliferation",
  "Basal to intestinal metaplasia" = "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia" = "SMG to intestinal metaplasia",
  "Stress adaptive" = "Stress adaptive",
  "Cancer-cell immune mimicry" = "Cancer-cell immune mimicry"
)
state_order <- c(
  "Classic proliferation", "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia", "Stress adaptive",
  "Cancer-cell immune mimicry"
)
state_colours <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

cancer_mp_descriptions <- c(
  "MP1" = "G2/M cell cycle", "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis", "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium", "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia", "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)
cancer_excluded_mps <- character(0)

celltype_colour_table <- data.frame(
  celltype = c("epithelial", "fibroblast", "endothelial", "macrophage", "mast", "t.cell", "b.cell", "nk.cell", "plasma", "dendritic", "lymph", "erythrocyte", "keratinocyte", "neutrophil", "unresolved", "combined"),
  colour = c("#D73027", "#8C564B", "#1F78B4", "#FF7F00", "#A65628", "#33A02C", "#377EB8", "#984EA3", "#E377C2", "#17BECF", "#6BAED6", "#7F7F7F", "#E6AB02", "#1B9E77", "#BDBDBD", "#555555"),
  stringsAsFactors = FALSE
)
celltype_colours <- setNames(celltype_colour_table$colour, celltype_colour_table$celltype)
celltype_colours <- celltype_colours[names(celltype_colours) != "combined"]
celltype_colours[c("nk", "cd4", "cd8")] <- celltype_colours[c("nk.cell", "t.cell", "t.cell")]

tme_mp_cfg <- data.table(
  compartment = c("fibroblast", "endothelial", "macrophage", "nk", "plasma", "cd4", "cd8"),
  annotation_celltype = c("fibroblast", "endothelial", "macrophage", "nk.cell", "plasma", "t.cell", "t.cell"),
  mp_path = file.path(
    wd, "ref_outs",
    c("nmf_fibroblast", "nmf_endothelial", "nmf_macrophage", "nmf_nk", "nmf_plasma", "nmf_cd4", "nmf_cd8"),
    "MP_outs_default.rds"
  )
)

mode_table <- data.table(
  mode_id = c(
    "01_cancer_mps_vs_tme_mps", "02_cancer_states_vs_tme_mps",
    "03_cancer_mps_vs_whole_celltypes", "04_cancer_states_vs_whole_celltypes"
  ),
  source_kind = c("cancer_mp", "cancer_state", "cancer_mp", "cancer_state"),
  target_kind = c("tme_mp", "tme_mp", "whole_celltype", "whole_celltype"),
  mode_label = c(
    "Malignant MPs vs non-malignant MPs", "Malignant states vs non-malignant MPs",
    "Malignant MPs vs whole cell types", "Malignant states vs whole cell types"
  ),
  scrna_subdir = c(
    "01_cancer_mps_cross_only", "03_cancer_states_cross_only",
    "05_cancer_mps_vs_whole_celltypes", "06_cancer_states_vs_whole_celltypes"
  )
)

run_start <- Sys.time()
run_log <- c(
  paste0("start=", format(run_start, tz = "Europe/London")),
  paste0("cache_version=", cache_version),
  paste0("samples=", paste(samples, collapse = ";")),
  "sample_alias=D1 is canonical manifest sample FFPED1",
  paste0("n_permutations=", n_permutations),
  paste0("n_cores=", n_cores),
  paste0("force_rebuild=", force_rebuild),
  paste0("replot_only=", replot_only)
)

as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

parse_grid_coordinates <- function(barcodes) {
  matched <- regexec("^s_016um_([0-9]+)_([0-9]+)-[0-9]+$", as.character(barcodes))
  pieces <- regmatches(as.character(barcodes), matched)
  valid <- lengths(pieces) == 3L
  row <- col <- rep(NA_integer_, length(barcodes))
  row[valid] <- as.integer(vapply(pieces[valid], `[[`, character(1), 2L))
  col[valid] <- as.integer(vapply(pieces[valid], `[[`, character(1), 3L))
  data.table(grid_row = row, grid_col = col)
}

normalise_log1p_cp10k <- function(counts) {
  libsize <- Matrix::colSums(counts)
  libsize[!is.finite(libsize) | libsize <= 0] <- 1
  normalised <- Matrix::t(Matrix::t(counts) * (1e4 / libsize))
  normalised@x <- log1p(normalised@x)
  normalised
}

read_10x_counts <- function(input_dir) {
  candidates <- file.path(input_dir, c("filtered_feature_bc_matrix.h5", "filtered_feature_cell_matrix.h5"))
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("No filtered 10x H5 matrix under ", input_dir)
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1L]]
  Matrix::Matrix(counts, sparse = TRUE)
}

filter_gene_sets_by_silhouette <- function(mp_outs) {
  mp_genes <- mp_outs$metaprograms.genes
  sil <- mp_outs$metaprograms.metrics$silhouette
  if (is.null(mp_genes) || is.null(sil)) stop("Malformed GeneNMF MP object")
  if (!is.null(names(sil)) && length(names(sil)) == length(mp_genes)) {
    sil_names <- names(sil)
  } else if (!is.null(names(mp_genes)) && length(names(mp_genes)) == length(sil)) {
    sil_names <- names(mp_genes)
  } else {
    sil_names <- paste0("MP", seq_along(sil))
  }
  names(sil) <- sil_names
  keep_names <- sil_names[!is.na(sil) & sil >= 0]
  mp_genes[intersect(names(mp_genes), keep_names)]
}

load_tme_gene_sets <- function() {
  rows <- list()
  gene_sets <- list()
  for (idx in seq_len(nrow(tme_mp_cfg))) {
    cfg <- tme_mp_cfg[idx]
    if (!file.exists(cfg$mp_path)) stop("Missing GeneNMF object: ", cfg$mp_path)
    message("Loading silhouette-filtered gene sets: ", cfg$compartment)
    mp_outs <- readRDS(cfg$mp_path)
    retained <- filter_gene_sets_by_silhouette(mp_outs)
    rm(mp_outs)
    invisible(gc())
    retained <- lapply(retained, function(genes) head(unique(as.character(genes)), top_genes_per_mp))
    for (mp_name in names(retained)) {
      feature_id <- paste(cfg$compartment, mp_name, sep = "__")
      gene_sets[[feature_id]] <- retained[[mp_name]]
      rows[[length(rows) + 1L]] <- data.table(
        feature_id = feature_id, compartment = cfg$compartment,
        annotation_celltype = cfg$annotation_celltype, mp = mp_name,
        feature_label = paste0(cfg$compartment, " ", mp_name),
        n_signature_genes = length(retained[[mp_name]]), mp_path = cfg$mp_path
      )
    }
  }
  list(gene_sets = gene_sets, feature_table = rbindlist(rows))
}

score_signature_matrix <- function(counts, gene_sets) {
  log_norm <- normalise_log1p_cp10k(counts)
  genes <- rownames(log_norm)
  means <- Matrix::rowMeans(log_norm)
  squared <- log_norm
  squared@x <- squared@x^2
  variances <- pmax(Matrix::rowMeans(squared) - means^2, 0)
  sds <- sqrt(variances)
  sds[!is.finite(sds) | sds == 0] <- 1
  score_list <- lapply(gene_sets, function(signature) {
    keep <- intersect(signature, genes)
    if (!length(keep)) return(rep(NA_real_, ncol(log_norm)))
    indices <- match(keep, genes)
    inv_sd <- 1 / sds[indices]
    scaled_sum <- Matrix::colSums(log_norm[indices, , drop = FALSE] * inv_sd)
    as.numeric(scaled_sum / length(indices) - mean(means[indices] * inv_sd))
  })
  raw <- do.call(cbind, score_list)
  colnames(raw) <- names(gene_sets)
  rownames(raw) <- colnames(counts)
  adjusted <- apply(raw, 2L, function(x) {
    centre <- mean(x, na.rm = TRUE)
    spread <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) spread <- 1
    (x - centre) / spread
  })
  if (is.null(dim(adjusted))) adjusted <- matrix(adjusted, ncol = 1L, dimnames = list(rownames(raw), colnames(raw)))
  adjusted
}

make_sparse_membership <- function(labels, levels) {
  column <- match(as.character(labels), levels)
  keep <- which(!is.na(column))
  Matrix::sparseMatrix(
    i = keep, j = column[keep], x = 1,
    dims = c(length(labels), length(levels)),
    dimnames = list(NULL, levels)
  )
}

make_adjacency <- function(source_meta, target_meta, ring) {
  target_keys <- paste(target_meta$grid_row, target_meta$grid_col, sep = ":")
  if (anyDuplicated(target_keys)) stop("Duplicate target 16 um grid positions")
  offsets <- CJ(dr = seq.int(-ring, ring), dc = seq.int(-ring, ring))
  offsets <- offsets[!(dr == 0L & dc == 0L)]
  edge_i <- edge_j <- vector("list", nrow(offsets))
  for (idx in seq_len(nrow(offsets))) {
    query <- paste(source_meta$grid_row + offsets$dr[idx], source_meta$grid_col + offsets$dc[idx], sep = ":")
    matched <- match(query, target_keys)
    keep <- which(!is.na(matched))
    edge_i[[idx]] <- keep
    edge_j[[idx]] <- matched[keep]
  }
  i <- unlist(edge_i, use.names = FALSE)
  j <- unlist(edge_j, use.names = FALSE)
  Matrix::sparseMatrix(
    i = i, j = j, x = 1,
    dims = c(nrow(source_meta), nrow(target_meta)),
    dimnames = list(source_meta$barcode, target_meta$barcode)
  )
}

safe_bh <- function(p) {
  result <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  result[keep] <- p.adjust(p[keep], method = "BH")
  result
}
####################

####################
# Persistent per-sample feature construction
####################
prepare_sample_membership <- function(sample_name, manifest, mapped_states, tme_reference) {
  message("Preparing spatial memberships for ", sample_name)
  malignancy_path <- file.path(malignancy_dir, paste0("Auto_", sample_name, "_binned_malignancy.csv.gz"))
  if (!file.exists(malignancy_path)) stop("Missing malignancy table: ", malignancy_path)
  cell_table <- fread(malignancy_path)
  required <- c(
    "barcode", "sample", "Auto_postfilter_celltype", "Auto_postfilter_keep",
    "pxl_col_in_fullres", "pxl_row_in_fullres", "Auto_malignancy"
  )
  missing <- setdiff(required, names(cell_table))
  if (length(missing)) stop(sample_name, " malignancy table lacks: ", paste(missing, collapse = ", "))
  cell_table <- cell_table[sample == sample_name]
  grid <- parse_grid_coordinates(cell_table$barcode)
  cell_table[, `:=`(grid_row = grid$grid_row, grid_col = grid$grid_col)]
  cell_table <- cell_table[
    as_bool(Auto_postfilter_keep) & is.finite(grid_row) & is.finite(grid_col) &
      is.finite(pxl_col_in_fullres) & is.finite(pxl_row_in_fullres)
  ]
  if (anyDuplicated(cell_table$barcode)) stop("Duplicate barcodes in ", sample_name)

  source_meta <- mapped_states[sample == sample_name]
  if (!nrow(source_meta)) stop("No mapped malignant bins for ", sample_name)
  source_grid <- parse_grid_coordinates(source_meta$barcode)
  source_meta[, `:=`(grid_row = source_grid$grid_row, grid_col = source_grid$grid_col)]
  source_meta <- source_meta[is.finite(grid_row) & is.finite(grid_col)]
  source_meta <- source_meta[barcode %in% cell_table$barcode]
  if (!nrow(source_meta)) stop("Mapped malignant bins do not overlap retained annotation for ", sample_name)
  setkey(cell_table, barcode)
  source_meta[, `:=`(
    pxl_col_in_fullres = cell_table[barcode, pxl_col_in_fullres],
    pxl_row_in_fullres = cell_table[barcode, pxl_row_in_fullres]
  )]

  raw_states <- as.character(source_meta$Auto_state_B)
  canonical_states <- unname(state_crosswalk[raw_states])
  state_membership <- make_sparse_membership(canonical_states, state_order)
  state_counts <- Matrix::colSums(state_membership)
  keep_states <- names(state_counts)[state_counts >= min_source_bins]
  state_membership <- state_membership[, keep_states, drop = FALSE]

  cancer_mp_columns <- grep("^Auto_adj_", names(source_meta), value = TRUE)
  cancer_mp_names <- sub("^Auto_adj_", "", cancer_mp_columns)
  keep_mp <- !cancer_mp_names %in% cancer_excluded_mps
  cancer_mp_columns <- cancer_mp_columns[keep_mp]
  cancer_mp_names <- cancer_mp_names[keep_mp]
  if (!length(cancer_mp_columns)) stop("No malignant adjusted MP columns for ", sample_name)
  cancer_mp_scores <- as.matrix(source_meta[, ..cancer_mp_columns])
  colnames(cancer_mp_scores) <- cancer_mp_names
  cancer_mp_membership <- Matrix::Matrix(cancer_mp_scores > mp_positive_threshold, sparse = TRUE)
  mp_counts <- Matrix::colSums(cancer_mp_membership)
  keep_cancer_mps <- names(mp_counts)[mp_counts >= min_source_bins]
  cancer_mp_membership <- cancer_mp_membership[, keep_cancer_mps, drop = FALSE]

  source_membership <- cbind(cancer_mp_membership, state_membership)
  source_features <- rbindlist(list(
    data.table(
      source_feature = keep_cancer_mps, source_kind = "cancer_mp",
      source_label = ifelse(
        keep_cancer_mps %in% names(cancer_mp_descriptions),
        paste0(keep_cancer_mps, " + ", cancer_mp_descriptions[keep_cancer_mps]), keep_cancer_mps
      ),
      source_validation_key = keep_cancer_mps,
      source_validation_key_whole = ifelse(
        keep_cancer_mps %in% names(cancer_mp_descriptions),
        unname(cancer_mp_descriptions[keep_cancer_mps]), keep_cancer_mps
      )
    ),
    data.table(
      source_feature = keep_states, source_kind = "cancer_state",
      source_label = keep_states, source_validation_key = keep_states,
      source_validation_key_whole = keep_states
    )
  ), fill = TRUE)
  colnames(source_membership) <- paste(source_features$source_kind, source_features$source_feature, sep = "__")
  source_features[, source_column := colnames(source_membership)]
  source_features[, n_positive_bins := as.integer(Matrix::colSums(source_membership))]

  target_meta <- cell_table[
    Auto_postfilter_celltype != "epithelial" &
      !is.na(Auto_postfilter_celltype) & Auto_postfilter_celltype != "unresolved"
  ]
  target_meta[, whole_celltype := as.character(Auto_postfilter_celltype)]
  whole_counts <- target_meta[, .N, by = whole_celltype]
  whole_keep <- whole_counts[N >= min_whole_target_bins, whole_celltype]
  whole_membership <- make_sparse_membership(target_meta$whole_celltype, whole_keep)
  whole_features <- data.table(
    target_feature = whole_keep, target_kind = "whole_celltype",
    target_compartment = whole_keep, target_annotation_celltype = whole_keep,
    target_label = whole_keep, permutation_block = "all_non_epithelial",
    n_positive_bins = as.integer(Matrix::colSums(whole_membership))
  )
  colnames(whole_membership) <- paste("whole_celltype", whole_keep, sep = "__")
  whole_features[, target_column := colnames(whole_membership)]

  input_dir <- manifest[sample == sample_name, binned_input][1L]
  if (is.na(input_dir) || !dir.exists(input_dir)) stop("Missing manifest input for ", sample_name)
  counts <- read_10x_counts(input_dir)
  missing_barcodes <- setdiff(target_meta$barcode, colnames(counts))
  if (length(missing_barcodes)) stop(sample_name, ": ", length(missing_barcodes), " retained target bins absent from counts")

  score_rows <- list()
  score_feature_rows <- list()
  target_mp_membership <- Matrix::Matrix(0, nrow = nrow(target_meta), ncol = 0, sparse = TRUE)
  target_index <- setNames(seq_len(nrow(target_meta)), target_meta$barcode)
  for (cfg_idx in seq_len(nrow(tme_mp_cfg))) {
    cfg <- tme_mp_cfg[cfg_idx]
    compartment_rows <- tme_reference$feature_table[compartment == cfg$compartment]
    bin_idx <- which(target_meta$whole_celltype == cfg$annotation_celltype)
    if (length(bin_idx) < min_parent_bins_for_mp) {
      score_feature_rows[[length(score_feature_rows) + 1L]] <- compartment_rows[, .(
        target_feature = feature_id, target_kind = "tme_mp", target_compartment = compartment,
        target_annotation_celltype = annotation_celltype, target_label = feature_label,
        permutation_block = annotation_celltype, n_parent_bins = length(bin_idx),
        n_positive_bins = 0L, eligible = FALSE, skip_reason = "insufficient_parent_celltype_bins"
      )]
      next
    }
    barcodes <- target_meta$barcode[bin_idx]
    feature_ids <- compartment_rows$feature_id
    sets <- tme_reference$gene_sets[feature_ids]
    signature_genes <- unique(unlist(sets, use.names = FALSE))
    keep_genes <- intersect(signature_genes, rownames(counts))
    if (!length(keep_genes)) stop("No TME signature genes detected for ", sample_name, " ", cfg$compartment)
    score_mat <- score_signature_matrix(counts[keep_genes, barcodes, drop = FALSE], sets)
    positive_mat <- Matrix::Matrix(score_mat > mp_positive_threshold, sparse = TRUE)
    positive_counts <- Matrix::colSums(positive_mat)
    keep_features <- names(positive_counts)[positive_counts >= min_mp_positive_bins]
    if (length(keep_features)) {
      full_block <- Matrix::sparseMatrix(
        i = rep(bin_idx, times = length(keep_features)),
        j = rep(seq_along(keep_features), each = length(bin_idx)),
        x = as.numeric(positive_mat[, keep_features, drop = FALSE]),
        dims = c(nrow(target_meta), length(keep_features)),
        dimnames = list(target_meta$barcode, paste0("tme_mp__", keep_features))
      )
      target_mp_membership <- cbind(target_mp_membership, full_block)
    }
    long_scores <- as.data.table(score_mat, keep.rownames = "barcode")
    long_scores <- melt(long_scores, id.vars = "barcode", variable.name = "feature_id", value.name = "adjusted_score")
    long_scores[, `:=`(
      sample = sample_name,
      compartment = tme_reference$feature_table$compartment[match(feature_id, tme_reference$feature_table$feature_id)],
      annotation_celltype = cfg$annotation_celltype,
      positive = adjusted_score > mp_positive_threshold
    )]
    score_rows[[length(score_rows) + 1L]] <- long_scores
    score_feature_rows[[length(score_feature_rows) + 1L]] <- compartment_rows[, .(
      target_feature = feature_id, target_kind = "tme_mp", target_compartment = compartment,
      target_annotation_celltype = annotation_celltype, target_label = feature_label,
      permutation_block = annotation_celltype, n_parent_bins = length(bin_idx),
      n_positive_bins = as.integer(positive_counts[feature_id]),
      eligible = feature_id %in% keep_features,
      skip_reason = ifelse(feature_id %in% keep_features, NA_character_, "insufficient_mp_positive_bins")
    )]
  }
  rm(counts)
  invisible(gc())

  mp_feature_audit <- rbindlist(score_feature_rows, fill = TRUE)
  eligible_mp_features <- mp_feature_audit[eligible %in% TRUE]
  if (ncol(target_mp_membership)) {
    expected_columns <- paste0("tme_mp__", eligible_mp_features$target_feature)
    target_mp_membership <- target_mp_membership[, expected_columns, drop = FALSE]
    eligible_mp_features[, target_column := expected_columns]
  } else {
    eligible_mp_features[, target_column := character()]
  }

  tme_scores <- rbindlist(score_rows, fill = TRUE)
  membership <- list(
    cache_version = cache_version, sample = sample_name,
    background = cell_table[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, Auto_postfilter_celltype, Auto_malignancy)],
    source_meta = source_meta[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, Auto_state_B)],
    target_meta = target_meta[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, whole_celltype)],
    source_membership = source_membership,
    source_features = source_features,
    target_membership = list(whole_celltype = whole_membership, tme_mp = target_mp_membership),
    target_features = list(whole_celltype = whole_features, tme_mp = eligible_mp_features),
    target_feature_audit = mp_feature_audit
  )
  list(membership = membership, scores = tme_scores)
}

membership_paths <- setNames(
  file.path(output_tiers[["intermediate"]], paste0("Auto_", samples, "_spatial_interaction_memberships.rds")),
  samples
)
tme_score_path <- file.path(output_tiers[["intermediate"]], "Auto_visiumhd_tme_mp_scores.rds")
tme_gene_set_path <- file.path(output_tiers[["intermediate"]], "Auto_visiumhd_tme_mp_gene_sets.rds")

if (!file.exists(state_path)) stop("Missing current mapped malignant state/MP table: ", state_path)
manifest <- data.table(
  sample = c("SUR1122", "SUR1231", "FFPEA1", "FFPED1"),
  binned_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"),
  segmented_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs")
)[sample %in% samples]
if (!setequal(manifest$sample, samples)) stop("Manifest lacks one or more requested samples")
mapped_states <- fread(state_path)
if (!setequal(unique(mapped_states$sample), samples)) stop("Mapped state table does not contain exactly the three requested samples")

if (replot_only) {
  required_replot <- c(unname(membership_paths), tme_score_path,
                       file.path(output_tiers[["tables"]], "Auto_spatial_interactions_by_sample.csv.gz"),
                       file.path(output_tiers[["tables"]], "Auto_spatial_interactions_pooled.csv.gz"))
  missing <- required_replot[!file.exists(required_replot)]
  if (length(missing)) stop("SCREF_REPLOT_ONLY missing: ", paste(missing, collapse = "; "))
  memberships <- lapply(membership_paths, readRDS)
  tme_scores_all <- readRDS(tme_score_path)
  run_log <- c(run_log, "membership_cache=replot_only")
} else {
  need_build <- force_rebuild || !all(file.exists(membership_paths)) || !file.exists(tme_score_path) || !file.exists(tme_gene_set_path)
  if (!need_build) {
    memberships <- lapply(membership_paths, readRDS)
    cache_valid <- all(vapply(memberships, function(x) identical(x$cache_version, cache_version), logical(1)))
    need_build <- !cache_valid
  }
  if (need_build) {
    tme_reference <- load_tme_gene_sets()
    saveRDS(tme_reference, tme_gene_set_path)
    prepared <- lapply(samples, prepare_sample_membership, manifest = manifest, mapped_states = mapped_states, tme_reference = tme_reference)
    memberships <- lapply(prepared, `[[`, "membership")
    names(memberships) <- samples
    for (sample_name in samples) saveRDS(memberships[[sample_name]], membership_paths[[sample_name]])
    tme_scores_all <- rbindlist(lapply(prepared, `[[`, "scores"), fill = TRUE)
    saveRDS(tme_scores_all, tme_score_path)
    fwrite(tme_scores_all, file.path(output_tiers[["intermediate"]], "Auto_visiumhd_tme_mp_scores.csv.gz"))
    run_log <- c(run_log, "membership_cache=rebuilt")
  } else {
    names(memberships) <- samples
    tme_scores_all <- readRDS(tme_score_path)
    run_log <- c(run_log, "membership_cache=reused")
  }
}

source_audit <- rbindlist(lapply(memberships, function(x) copy(x$source_features)[, sample := x$sample]), fill = TRUE)
target_audit <- rbindlist(lapply(memberships, function(x) {
  rbindlist(list(
    copy(x$target_features$whole_celltype),
    copy(x$target_feature_audit)
  ), fill = TRUE)[, sample := x$sample]
}), fill = TRUE)
fwrite(source_audit, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_source_feature_audit.csv"))
fwrite(target_audit, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_target_feature_audit.csv"))
####################

####################
# Degree-aware neighbourhood permutation tests
####################
permutation_summary_chunk <- function(chunk_n, source_membership, adjacency, target_membership,
                                      target_meta, target_features, target_kind, seed) {
  set.seed(seed)
  dims <- c(ncol(source_membership), ncol(target_membership))
  sum_null <- matrix(0, nrow = dims[1], ncol = dims[2])
  sumsq_null <- matrix(0, nrow = dims[1], ncol = dims[2])
  greater_equal <- matrix(0L, nrow = dims[1], ncol = dims[2])
  less_equal <- matrix(0L, nrow = dims[1], ncol = dims[2])
  observed <- as.matrix(Matrix::crossprod(source_membership, adjacency %*% target_membership))
  if (target_kind == "whole_celltype") {
    blocks <- list(all_non_epithelial = seq_len(nrow(target_meta)))
    feature_blocks <- list(all_non_epithelial = seq_len(ncol(target_membership)))
  } else {
    block_names <- unique(target_features$permutation_block)
    blocks <- setNames(lapply(block_names, function(block) which(target_meta$whole_celltype == block)), block_names)
    feature_blocks <- setNames(lapply(block_names, function(block) which(target_features$permutation_block == block)), block_names)
  }
  for (iteration in seq_len(chunk_n)) {
    permuted_counts <- matrix(0, nrow = dims[1], ncol = dims[2])
    for (block in names(blocks)) {
      row_idx <- blocks[[block]]
      feature_idx <- feature_blocks[[block]]
      if (!length(row_idx) || !length(feature_idx)) next
      permuted_rows <- sample(row_idx, length(row_idx), replace = FALSE)
      block_counts <- Matrix::crossprod(
        source_membership,
        adjacency[, row_idx, drop = FALSE] %*% target_membership[permuted_rows, feature_idx, drop = FALSE]
      )
      permuted_counts[, feature_idx] <- as.matrix(block_counts)
    }
    sum_null <- sum_null + permuted_counts
    sumsq_null <- sumsq_null + permuted_counts^2
    greater_equal <- greater_equal + (permuted_counts >= observed)
    less_equal <- less_equal + (permuted_counts <= observed)
  }
  list(n = chunk_n, sum = sum_null, sumsq = sumsq_null, ge = greater_equal, le = less_equal)
}

run_permutation_test <- function(membership, ring, target_kind) {
  adjacency <- make_adjacency(membership$source_meta, membership$target_meta, ring)
  source_membership <- membership$source_membership
  target_membership <- membership$target_membership[[target_kind]]
  source_features <- membership$source_features
  target_features <- membership$target_features[[target_kind]]
  if (!ncol(source_membership) || !ncol(target_membership) || !length(adjacency@x)) return(data.table())
  observed <- as.matrix(Matrix::crossprod(source_membership, adjacency %*% target_membership))
  source_slots <- as.numeric(Matrix::crossprod(source_membership, Matrix::rowSums(adjacency)))
  source_counts <- as.numeric(Matrix::colSums(source_membership))

  chunks <- rep(n_permutations %/% n_cores, n_cores)
  if (n_permutations %% n_cores) chunks[seq_len(n_permutations %% n_cores)] <- chunks[seq_len(n_permutations %% n_cores)] + 1L
  chunks <- chunks[chunks > 0]
  seeds <- 100000L + match(membership$sample, samples) * 10000L + ring * 100L +
    ifelse(target_kind == "whole_celltype", 1L, 2L) + seq_along(chunks)
  chunk_results <- parallel::mclapply(
    seq_along(chunks),
    function(idx) permutation_summary_chunk(
      chunks[idx], source_membership, adjacency, target_membership,
      membership$target_meta, target_features, target_kind, seeds[idx]
    ),
    mc.cores = min(n_cores, length(chunks)), mc.preschedule = TRUE
  )
  null_n <- sum(vapply(chunk_results, `[[`, integer(1), "n"))
  null_sum <- Reduce(`+`, lapply(chunk_results, `[[`, "sum"))
  null_sumsq <- Reduce(`+`, lapply(chunk_results, `[[`, "sumsq"))
  ge <- Reduce(`+`, lapply(chunk_results, `[[`, "ge"))
  le <- Reduce(`+`, lapply(chunk_results, `[[`, "le"))
  null_mean <- null_sum / null_n
  null_variance <- pmax((null_sumsq - null_sum^2 / null_n) / max(null_n - 1L, 1L), 0)
  null_sd <- sqrt(null_variance)
  z_score <- (observed - null_mean) / null_sd
  z_score[null_sd == 0 & observed == null_mean] <- 0
  z_score[null_sd == 0 & observed > null_mean] <- Inf
  z_score[null_sd == 0 & observed < null_mean] <- -Inf
  p_z_matrix <- matrix(
    pmax(2 * pnorm(-abs(z_score)), .Machine$double.xmin),
    nrow = nrow(observed), ncol = ncol(observed)
  )
  p_empirical_matrix <- matrix(
    pmin(1, 2 * pmin((ge + 1) / (null_n + 1), (le + 1) / (null_n + 1))),
    nrow = nrow(observed), ncol = ncol(observed)
  )

  grid <- CJ(source_index = seq_len(nrow(observed)), target_index = seq_len(ncol(observed)))
  matrix_lookup <- cbind(grid$source_index, grid$target_index)
  grid[, `:=`(
    sample = membership$sample, ring = ring, target_kind = target_kind,
    observed_contacts = as.numeric(observed[matrix_lookup]),
    null_mean_contacts = as.numeric(null_mean[matrix_lookup]),
    null_sd_contacts = as.numeric(null_sd[matrix_lookup]),
    z_score = as.numeric(z_score[matrix_lookup]),
    p_z = as.numeric(p_z_matrix[matrix_lookup]),
    p_empirical = as.numeric(p_empirical_matrix[matrix_lookup]),
    source_neighbour_slots = source_slots[source_index],
    n_source_positive_bins = source_counts[source_index],
    n_target_positive_bins = as.numeric(Matrix::colSums(target_membership))[target_index],
    n_permutations = null_n
  )]
  grid[, log2_enrichment := log2((observed_contacts + 0.5) / (null_mean_contacts + 0.5))]
  grid <- cbind(
    grid,
    source_features[grid$source_index, .(
      source_feature, source_kind, source_label,
      source_validation_key, source_validation_key_whole
    )],
    target_features[grid$target_index, .(
      target_feature, target_compartment, target_annotation_celltype,
      target_label, permutation_block
    )]
  )
  grid[, mode_id := mode_table$mode_id[match(paste(source_kind, target_kind), paste(mode_table$source_kind, mode_table$target_kind))]]
  grid[, mode_label := mode_table$mode_label[match(mode_id, mode_table$mode_id)]]
  grid[, eligible_pair :=
    n_source_positive_bins >= min_source_bins &
      n_target_positive_bins >= ifelse(target_kind == "whole_celltype", min_whole_target_bins, min_mp_positive_bins) &
      source_neighbour_slots >= min_observed_neighbour_slots]
  grid[eligible_pair == FALSE, `:=`(p_z = NA_real_, p_empirical = NA_real_)]
  grid[, fdr_sample := safe_bh(p_z), by = .(sample, mode_id, ring)]
  grid[, direction := fifelse(log2_enrichment > 0, "enriched", fifelse(log2_enrichment < 0, "depleted", "neutral"))]
  grid[]
}

sample_result_path <- file.path(output_tiers[["tables"]], "Auto_spatial_interactions_by_sample.csv.gz")
if (replot_only) {
  sample_results <- fread(sample_result_path)
} else {
  test_rows <- list()
  for (sample_name in samples) {
    membership <- memberships[[sample_name]]
    for (ring in neighbour_rings) {
      for (target_kind in c("whole_celltype", "tme_mp")) {
        message("Permutation test: ", sample_name, " ring ", ring, " ", target_kind)
        test_rows[[length(test_rows) + 1L]] <- run_permutation_test(membership, ring, target_kind)
      }
    }
  }
  sample_results <- rbindlist(test_rows, fill = TRUE)
  setcolorder(sample_results, c(
    "sample", "mode_id", "mode_label", "ring", "source_kind", "source_feature", "source_label",
    "target_kind", "target_compartment", "target_feature", "target_label"
  ))
  fwrite(sample_results, sample_result_path)
}
sample_results[is.finite(p_z), p_z := pmax(p_z, .Machine$double.xmin)]
sample_results[, fdr_sample := safe_bh(p_z), by = .(sample, mode_id, ring)]
fwrite(sample_results, sample_result_path)
####################

####################
# Cross-sample meta-analysis, recurrence, and scRNA validation
####################
pooled_results <- sample_results[eligible_pair %in% TRUE & is.finite(z_score), {
  weights <- sqrt(pmax(n_source_positive_bins, 1))
  pooled_z <- sum(weights * z_score) / sqrt(sum(weights^2))
  .(
    n_samples_tested = .N,
    samples_tested = paste(sample, collapse = ";"),
    observed_contacts = sum(observed_contacts),
    null_mean_contacts = sum(null_mean_contacts),
    pooled_z = pooled_z,
    pooled_p = pmax(2 * pnorm(-abs(pooled_z)), .Machine$double.xmin),
    pooled_log2_enrichment = log2((sum(observed_contacts) + 0.5) / (sum(null_mean_contacts) + 0.5)),
    min_sample_fdr = min(fdr_sample, na.rm = TRUE),
    n_samples_fdr_0_05 = sum(fdr_sample < sample_fdr_cutoff, na.rm = TRUE),
    n_samples_fdr_0_10 = sum(fdr_sample < recurrence_fdr_cutoff, na.rm = TRUE),
    n_samples_enriched = sum(log2_enrichment > 0),
    n_samples_depleted = sum(log2_enrichment < 0),
    sample_log2_enrichment = paste(paste0(sample, "=", sprintf("%.4f", log2_enrichment)), collapse = ";"),
    sample_fdr = paste(paste0(sample, "=", format(fdr_sample, scientific = TRUE, digits = 3)), collapse = ";")
  )
}, by = .(
  mode_id, mode_label, ring, source_kind, source_feature, source_label,
  source_validation_key, source_validation_key_whole,
  target_kind, target_compartment, target_annotation_celltype,
  target_feature, target_label
)]
pooled_results[n_samples_tested < 2L, `:=`(pooled_p = NA_real_, pooled_z = NA_real_)]
pooled_results[, pooled_fdr := safe_bh(pooled_p), by = .(mode_id, ring)]
pooled_results[, pooled_direction := fifelse(
  pooled_log2_enrichment > 0, "enriched",
  fifelse(pooled_log2_enrichment < 0, "depleted", "neutral")
)]
pooled_results[, direction_concordant :=
  (pooled_log2_enrichment > 0 & n_samples_enriched == n_samples_tested) |
    (pooled_log2_enrichment < 0 & n_samples_depleted == n_samples_tested)]

recurrence_results <- pooled_results[, .(
  mode_id, mode_label, ring, source_feature, source_label,
  target_compartment, target_feature, target_label,
  n_samples_tested, n_samples_fdr_0_05, n_samples_fdr_0_10,
  n_samples_enriched, n_samples_depleted, direction_concordant,
  pooled_log2_enrichment, pooled_z, pooled_p, pooled_fdr,
  recurrent_same_direction = direction_concordant & n_samples_fdr_0_10 >= 2L,
  pooled_significant = pooled_fdr < pooled_fdr_cutoff
)]

scrna_inputs <- file.path(
  wd, "ref_outs", "non_malignant_mp_correlations", mode_table$scrna_subdir,
  "Auto_celltype_correlations_all.csv"
)
names(scrna_inputs) <- mode_table$mode_id
load_scrna_validation <- function(mode_id) {
  path <- scrna_inputs[[mode_id]]
  if (!file.exists(path)) return(data.table())
  input <- fread(path)
  cancer_first <- input$compartment1 == "cancer"
  result <- data.table(
    mode_id = mode_id,
    scrna_source_feature = ifelse(cancer_first, input$mp1_name, input$mp2_name),
    scrna_source_display = ifelse(cancer_first, input$mp1_display, input$mp2_display),
    target_compartment = ifelse(cancer_first, input$celltype2_display, input$celltype1_display),
    scrna_target_feature = ifelse(cancer_first, input$mp2_name, input$mp1_name),
    scrna_target_display = ifelse(cancer_first, input$mp2_display, input$mp1_display),
    scrna_spearman_r = input$spearman_r,
    scrna_spearman_p = input$spearman_p,
    scrna_spearman_sig = as_bool(input$spearman_sig),
    scrna_shared_samples = input$shared_sample_n,
    scrna_edge_id = input$edge_id
  )
  if (mode_id %in% c("03_cancer_mps_vs_whole_celltypes", "04_cancer_states_vs_whole_celltypes")) {
    result[, scrna_target_feature := target_compartment]
  }
  result
}
scrna_validation <- rbindlist(lapply(mode_table$mode_id, load_scrna_validation), fill = TRUE)

validation_base <- copy(pooled_results)
validation_base[, spatial_source_key := ifelse(
  mode_id == "03_cancer_mps_vs_whole_celltypes",
  source_validation_key_whole, source_validation_key
)]
validation_base[, spatial_target_key := ifelse(
  target_kind == "whole_celltype", target_compartment,
  sub("^[^_]+__", "", target_feature)
)]
validation <- merge(
  validation_base,
  scrna_validation,
  by.x = c("mode_id", "spatial_source_key", "target_compartment", "spatial_target_key"),
  by.y = c("mode_id", "scrna_source_feature", "target_compartment", "scrna_target_feature"),
  all.x = TRUE
)
validation[, `:=`(
  scrna_available = is.finite(scrna_spearman_r),
  spatial_scrna_direction_concordant = is.finite(scrna_spearman_r) &
    sign(pooled_log2_enrichment) == sign(scrna_spearman_r),
  validation_class = fifelse(
    !is.finite(scrna_spearman_r), "No matching scRNA test",
    fifelse(
      pooled_fdr < pooled_fdr_cutoff & scrna_spearman_sig & sign(pooled_log2_enrichment) == sign(scrna_spearman_r),
      "Significant and direction-concordant in both",
      fifelse(
        pooled_fdr < pooled_fdr_cutoff & sign(pooled_log2_enrichment) == sign(scrna_spearman_r),
        "Spatial hit with concordant scRNA direction",
        fifelse(pooled_fdr < pooled_fdr_cutoff, "Spatial hit with discordant scRNA direction", "Not a pooled spatial hit")
      )
    )
  )
)]

fwrite(pooled_results, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_pooled.csv.gz"))
fwrite(recurrence_results, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_recurrence.csv"))
fwrite(validation, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_scrna_validation.csv.gz"))
fwrite(
  sample_results[eligible_pair %in% TRUE & fdr_sample < sample_fdr_cutoff],
  file.path(output_tiers[["tables"]], "Auto_spatial_interactions_significant_by_sample.csv")
)
fwrite(
  pooled_results[pooled_fdr < pooled_fdr_cutoff],
  file.path(output_tiers[["tables"]], "Auto_spatial_interactions_significant_pooled.csv")
)
####################

####################
# Presentation-readable interaction summaries
####################
interaction_scale_fill <- scale_fill_gradient2(
  low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
  limits = c(-max(0.5, max(abs(sample_results$log2_enrichment), na.rm = TRUE)),
              max(0.5, max(abs(sample_results$log2_enrichment), na.rm = TRUE))),
  oob = scales::squish, name = "log2 spatial\nenrichment"
)

empty_plot <- function(title, subtitle = NULL) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No interactions passed the stated threshold", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void(base_size = 14) + labs(title = title, subtitle = subtitle)
}

make_sample_dotmap <- function(plot_data, mode_row, target_compartment = NULL) {
  data <- plot_data[mode_id == mode_row$mode_id & ring == primary_ring & fdr_sample < sample_fdr_cutoff]
  if (!is.null(target_compartment)) {
    compartment_value <- target_compartment
    data <- data[target_compartment == compartment_value]
  }
  title <- mode_row$mode_label
  if (!is.null(target_compartment)) title <- paste0(title, " | ", target_compartment)
  if (!nrow(data)) return(empty_plot(title, "Immediate 16 um-bin contacts; sample BH FDR < 0.05"))
  data[, sample_plot := factor(sample_display[sample], levels = unname(sample_display[samples]))]
  data[, source_label := factor(source_label, levels = rev(unique(source_label[order(source_feature)])))]
  data[, target_label := factor(target_label, levels = unique(target_label[order(target_feature)]))]
  ggplot(data, aes(x = target_label, y = source_label)) +
    geom_point(aes(size = pmin(-log10(fdr_sample), 12), fill = log2_enrichment), shape = 21, colour = "grey20", stroke = 0.35) +
    facet_wrap(~sample_plot, nrow = 1, drop = FALSE) +
    interaction_scale_fill +
    scale_size_continuous(range = c(2.5, 10), name = "-log10 FDR\n(capped at 12)") +
    labs(title = title, subtitle = "Each sample tested independently | immediate 16 um-bin contacts", x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 13, face = "bold"), legend.text = element_text(size = 11),
      legend.title = element_text(size = 12), plot.title = element_text(size = 17, face = "bold")
    )
}

sample_dotmap_path <- file.path(output_tiers[["figures"]], "Auto_spatial_interaction_dotmap_with_neg.pdf")
pdf(sample_dotmap_path, width = 17, height = 9.5, onefile = TRUE, useDingbats = FALSE)
for (mode_idx in seq_len(nrow(mode_table))) {
  mode_row <- mode_table[mode_idx]
  if (mode_row$target_kind == "tme_mp") {
    compartments <- unique(sample_results[mode_id == mode_row$mode_id, target_compartment])
    for (compartment in compartments) print(make_sample_dotmap(sample_results, mode_row, compartment))
  } else {
    print(make_sample_dotmap(sample_results, mode_row))
  }
}
dev.off()

make_pooled_dotmap <- function(plot_data, mode_row, target_compartment = NULL) {
  data <- plot_data[mode_id == mode_row$mode_id & ring == primary_ring & pooled_fdr < pooled_fdr_cutoff]
  if (!is.null(target_compartment)) {
    compartment_value <- target_compartment
    data <- data[target_compartment == compartment_value]
  }
  title <- paste0("Three-sample pooled: ", mode_row$mode_label)
  if (!is.null(target_compartment)) title <- paste0(title, " | ", target_compartment)
  if (!nrow(data)) return(empty_plot(title, "Weighted Stouffer meta-analysis; BH FDR < 0.05"))
  data[, source_label := factor(source_label, levels = rev(unique(source_label[order(source_feature)])))]
  data[, target_label := factor(target_label, levels = unique(target_label[order(target_feature)]))]
  ggplot(data, aes(x = target_label, y = source_label)) +
    geom_point(aes(size = pmin(-log10(pooled_fdr), 12), fill = pooled_log2_enrichment), shape = 21, colour = "grey20", stroke = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "pooled log2\nenrichment") +
    scale_size_continuous(range = c(3, 11), name = "-log10 pooled FDR\n(capped at 12)") +
    labs(title = title, subtitle = "No cross-sample neighbours; sample z-scores combined after independent tests", x = NULL, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.text.y = element_text(size = 11),
      legend.text = element_text(size = 11), legend.title = element_text(size = 12),
      plot.title = element_text(size = 17, face = "bold")
    )
}

pooled_dotmap_path <- file.path(output_tiers[["figures"]], "Auto_spatial_interaction_dotmap_with_neg_pooled.pdf")
pdf(pooled_dotmap_path, width = 14, height = 9.5, onefile = TRUE, useDingbats = FALSE)
for (mode_idx in seq_len(nrow(mode_table))) {
  mode_row <- mode_table[mode_idx]
  if (mode_row$target_kind == "tme_mp") {
    compartments <- unique(pooled_results[mode_id == mode_row$mode_id, target_compartment])
    for (compartment in compartments) print(make_pooled_dotmap(pooled_results, mode_row, compartment))
  } else {
    print(make_pooled_dotmap(pooled_results, mode_row))
  }
}
dev.off()

recurrence_plot_data <- pooled_results[
  ring == primary_ring & (pooled_fdr < pooled_fdr_cutoff | n_samples_fdr_0_10 >= 2L)
]
if (nrow(recurrence_plot_data)) {
  recurrence_plot_data[, pair_label := paste(source_label, target_label, sep = "  <->  ")]
  recurrence_plot_data[, pair_label := factor(pair_label, levels = rev(unique(pair_label[order(mode_id, pooled_fdr)])))]
  recurrence_plot <- ggplot(recurrence_plot_data, aes(x = factor(n_samples_fdr_0_10), y = pair_label)) +
    geom_point(aes(size = pmin(-log10(pooled_fdr), 12), fill = pooled_log2_enrichment), shape = 21, colour = "grey20") +
    facet_wrap(~mode_label, scales = "free_y", ncol = 2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "pooled log2\nenrichment") +
    scale_size_continuous(range = c(3, 10), name = "-log10 pooled FDR") +
    labs(x = "Samples with same-pair FDR < 0.10", y = NULL, title = "Recurrent and pooled Visium HD interactions") +
    theme_classic(base_size = 13) +
    theme(axis.text.y = element_text(size = 8), strip.text = element_text(size = 12, face = "bold"), plot.title = element_text(size = 17, face = "bold"))
  recurrence_height <- max(9, min(24, 5 + 0.22 * nrow(recurrence_plot_data)))
  ggsave(file.path(output_tiers[["figures"]], "Auto_spatial_interaction_recurrence.pdf"), recurrence_plot, width = 16, height = recurrence_height, limitsize = FALSE)
} else {
  ggsave(file.path(output_tiers[["figures"]], "Auto_spatial_interaction_recurrence.pdf"), empty_plot("Recurrent and pooled Visium HD interactions"), width = 14, height = 8)
}

sensitivity_plot_data <- pooled_results[is.finite(pooled_log2_enrichment)]
sensitivity_plot_data[, primary_value := pooled_log2_enrichment[ring == primary_ring][match(
  paste(mode_id, source_feature, target_feature),
  paste(mode_id[ring == primary_ring], source_feature[ring == primary_ring], target_feature[ring == primary_ring])
)]]
sensitivity_plot <- ggplot(sensitivity_plot_data[ring != primary_ring & is.finite(primary_value)],
                           aes(x = primary_value, y = pooled_log2_enrichment, colour = factor(ring))) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
  geom_point(alpha = 0.45, size = 1.8) +
  facet_wrap(~mode_label, scales = "free", ncol = 2) +
  scale_colour_manual(values = c("2" = "#E69F00", "3" = "#0072B2"), labels = c("32 um", "48 um"), name = "Neighbour ring") +
  labs(x = "Immediate-contact pooled log2 enrichment (16 um)", y = "Broader-neighbourhood pooled log2 enrichment", title = "Neighbourhood-radius sensitivity") +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(size = 12, face = "bold"), plot.title = element_text(size = 17, face = "bold"))
ggsave(file.path(output_tiers[["figures"]], "Auto_spatial_interaction_neighbourhood_sensitivity.pdf"), sensitivity_plot, width = 14, height = 10)

validation_plot_data <- validation[
  ring == primary_ring & is.finite(scrna_spearman_r) & is.finite(pooled_fdr) &
    (pooled_fdr < pooled_fdr_cutoff | scrna_spearman_sig)
]
if (nrow(validation_plot_data)) {
  validation_plot <- ggplot(validation_plot_data, aes(x = scrna_spearman_r, y = pooled_log2_enrichment)) +
    geom_hline(yintercept = 0, colour = "grey65") + geom_vline(xintercept = 0, colour = "grey65") +
    geom_point(aes(fill = validation_class, size = pmin(-log10(pooled_fdr), 10)), shape = 21, colour = "grey20", alpha = 0.85) +
    facet_wrap(~mode_label, scales = "free", ncol = 2) +
    scale_fill_manual(values = c(
      "Significant and direction-concordant in both" = "#1B9E77",
      "Spatial hit with concordant scRNA direction" = "#66A61E",
      "Spatial hit with discordant scRNA direction" = "#D95F02",
      "Not a pooled spatial hit" = "#BDBDBD"
    ), drop = FALSE) +
    scale_size_continuous(range = c(2.5, 8), name = "-log10 spatial FDR") +
    labs(x = "scRNA sample-level Spearman rho", y = "Visium HD pooled log2 neighbourhood enrichment",
         fill = "Validation class", title = "Spatial validation of scRNA-predicted interactions") +
    theme_classic(base_size = 14) +
    theme(strip.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10), plot.title = element_text(size = 17, face = "bold"))
  ggsave(file.path(output_tiers[["figures"]], "Auto_spatial_interaction_scrna_validation.pdf"), validation_plot, width = 14, height = 10)
} else {
  ggsave(file.path(output_tiers[["figures"]], "Auto_spatial_interaction_scrna_validation.pdf"), empty_plot("Spatial validation of scRNA-predicted interactions"), width = 14, height = 8)
}

####################
# Spatial maps for significant or best-ranked audit interactions
####################
select_spatial_map_pairs <- function(results, sample_name, mode_id) {
  mode_value <- mode_id
  candidates <- results[
    sample == sample_name & mode_id == mode_value & ring == primary_ring & eligible_pair %in% TRUE & is.finite(fdr_sample)
  ]
  if (!nrow(candidates)) return(candidates)
  significant <- candidates[fdr_sample < sample_fdr_cutoff]
  if (nrow(significant)) {
    selected <- rbindlist(list(
      significant[log2_enrichment > 0][order(fdr_sample, -log2_enrichment)][seq_len(min(.N, spatial_pages_per_direction))],
      significant[log2_enrichment < 0][order(fdr_sample, log2_enrichment)][seq_len(min(.N, spatial_pages_per_direction))]
    ), fill = TRUE)
    selected[, selection_basis := "sample BH FDR < 0.05"]
  } else {
    selected <- rbindlist(list(
      candidates[log2_enrichment > 0][order(fdr_sample, -log2_enrichment)][1L],
      candidates[log2_enrichment < 0][order(fdr_sample, log2_enrichment)][1L]
    ), fill = TRUE)
    selected[, selection_basis := "exploratory best-ranked; no sample FDR hit"]
  }
  unique(selected, by = c("source_feature", "target_feature"))
}

make_spatial_interaction_plot <- function(membership, pair_row) {
  source_column <- paste(pair_row$source_kind, pair_row$source_feature, sep = "__")
  target_column <- if (pair_row$target_kind == "whole_celltype") {
    paste("whole_celltype", pair_row$target_feature, sep = "__")
  } else {
    paste("tme_mp", pair_row$target_feature, sep = "__")
  }
  source_idx <- which(as.numeric(membership$source_membership[, source_column]) > 0)
  target_idx <- which(as.numeric(membership$target_membership[[pair_row$target_kind]][, target_column]) > 0)
  adjacency <- make_adjacency(membership$source_meta, membership$target_meta, primary_ring)
  contact_summary <- summary(adjacency[source_idx, target_idx, drop = FALSE])
  if (nrow(contact_summary)) {
    contacts <- data.table(
      source_index = source_idx[contact_summary$i],
      target_index = target_idx[contact_summary$j]
    )
    if (nrow(contacts) > maximum_plotted_contact_edges) {
      set.seed(20260809)
      contacts <- contacts[sample(.N, maximum_plotted_contact_edges)]
    }
    segments <- data.table(
      x = membership$source_meta$pxl_col_in_fullres[contacts$source_index],
      y = membership$source_meta$pxl_row_in_fullres[contacts$source_index],
      xend = membership$target_meta$pxl_col_in_fullres[contacts$target_index],
      yend = membership$target_meta$pxl_row_in_fullres[contacts$target_index]
    )
  } else {
    segments <- data.table(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  }
  source_points <- membership$source_meta[source_idx]
  target_points <- membership$target_meta[target_idx]
  point_size <- max(0.22, min(1.1, 50000 / max(nrow(membership$background), 1)))
  ggplot() +
    geom_point(
      data = membership$background,
      aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres),
      colour = "grey87", size = point_size * 0.65, alpha = 0.55
    ) +
    geom_segment(
      data = segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#6A3D9A", linewidth = 0.18, alpha = 0.14
    ) +
    geom_point(
      data = target_points,
      aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres),
      colour = "#2166AC", size = point_size, alpha = 0.78
    ) +
    geom_point(
      data = source_points,
      aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres),
      colour = "#B2182B", size = point_size, alpha = 0.82
    ) +
    scale_y_reverse() + coord_fixed() +
    labs(
      title = paste0(sample_display[[membership$sample]], " | ", pair_row$source_label, " <-> ", pair_row$target_label),
      subtitle = paste0(
        pair_row$selection_basis, " | log2 enrichment=", sprintf("%.2f", pair_row$log2_enrichment),
        " | FDR=", format(pair_row$fdr_sample, digits = 2, scientific = TRUE),
        " | observed/null contacts=", pair_row$observed_contacts, "/", sprintf("%.1f", pair_row$null_mean_contacts)
      ),
      caption = "Red: malignant source-positive bins | blue: TME target-positive bins | purple: immediate contacts (subsampled when dense)",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11), plot.caption = element_text(size = 10, hjust = 0)
    )
}

spatial_selection_rows <- list()
for (mode_idx in seq_len(nrow(mode_table))) {
  mode_row <- mode_table[mode_idx]
  spatial_path <- file.path(output_tiers[["figures"]], paste0("Auto_", mode_row$mode_id, "_spatial_maps.pdf"))
  pdf(spatial_path, width = 14, height = 9.5, onefile = TRUE, useDingbats = FALSE)
  for (sample_name in samples) {
    selected <- select_spatial_map_pairs(sample_results, sample_name, mode_row$mode_id)
    if (!nrow(selected)) {
      print(empty_plot(paste0(sample_display[[sample_name]], " | ", mode_row$mode_label), "No eligible pair for a spatial audit map"))
      next
    }
    spatial_selection_rows[[length(spatial_selection_rows) + 1L]] <- selected
    for (row_idx in seq_len(nrow(selected))) {
      print(make_spatial_interaction_plot(memberships[[sample_name]], selected[row_idx]))
    }
  }
  dev.off()
}
spatial_selections <- rbindlist(spatial_selection_rows, fill = TRUE)
fwrite(spatial_selections, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_map_selections.csv"))

parameters <- data.table(
  parameter = c(
    "cache_version", "samples", "sample_alias", "bin_size_um", "primary_ring",
    "sensitivity_rings", "n_permutations", "mp_positive_threshold",
    "top_genes_per_mp", "min_source_bins", "min_whole_target_bins",
    "min_parent_bins_for_mp", "min_mp_positive_bins", "min_observed_neighbour_slots",
    "sample_fdr_cutoff", "pooled_fdr_cutoff", "pooled_method",
    "whole_celltype_null", "tme_mp_null", "state_crosswalk"
  ),
  value = c(
    cache_version, paste(samples, collapse = ";"), "D1=FFPED1", "16", as.character(primary_ring),
    paste(neighbour_rings, collapse = ";"), as.character(n_permutations), as.character(mp_positive_threshold),
    as.character(top_genes_per_mp), as.character(min_source_bins), as.character(min_whole_target_bins),
    as.character(min_parent_bins_for_mp), as.character(min_mp_positive_bins), as.character(min_observed_neighbour_slots),
    as.character(sample_fdr_cutoff), as.character(pooled_fdr_cutoff), "weighted Stouffer; weight=sqrt(source-positive bins)",
    "permute whole-celltype labels across retained non-epithelial bins within sample",
    "permute joint MP-positive matrix within its annotated parent celltype and sample",
    paste(paste(names(state_crosswalk), state_crosswalk, sep = "->"), collapse = ";")
  )
)
fwrite(parameters, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_parameters.csv"))

compact_summary <- rbindlist(lapply(seq_len(nrow(mode_table)), function(idx) {
  mode_row <- mode_table[idx]
  sample_primary <- sample_results[mode_id == mode_row$mode_id & ring == primary_ring & eligible_pair %in% TRUE]
  pooled_primary <- pooled_results[mode_id == mode_row$mode_id & ring == primary_ring]
  data.table(
    mode_id = mode_row$mode_id, mode_label = mode_row$mode_label,
    n_sample_pair_tests = nrow(sample_primary),
    n_sample_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff, na.rm = TRUE),
    n_sample_enriched_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff & sample_primary$log2_enrichment > 0, na.rm = TRUE),
    n_sample_depleted_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff & sample_primary$log2_enrichment < 0, na.rm = TRUE),
    n_pooled_pair_tests = nrow(pooled_primary),
    n_pooled_fdr_hits = sum(pooled_primary$pooled_fdr < pooled_fdr_cutoff, na.rm = TRUE),
    n_recurrent_same_direction = sum(pooled_primary$direction_concordant & pooled_primary$n_samples_fdr_0_10 >= 2L, na.rm = TRUE),
    n_spatial_scrna_concordant_both_significant = validation[
      mode_id == mode_row$mode_id & ring == primary_ring & validation_class == "Significant and direction-concordant in both", .N
    ]
  )
}))
fwrite(compact_summary, file.path(summary_dir, "visium_hd_spatial_cancer_tme_interactions_summary.csv"))
fwrite(compact_summary, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_compact_summary.csv"))

run_log <- c(
  run_log,
  paste0("n_sample_results=", nrow(sample_results)),
  paste0("n_pooled_results=", nrow(pooled_results)),
  paste0("n_sample_primary_fdr_hits=", sum(sample_results$ring == primary_ring & sample_results$fdr_sample < sample_fdr_cutoff, na.rm = TRUE)),
  paste0("n_pooled_primary_fdr_hits=", sum(pooled_results$ring == primary_ring & pooled_results$pooled_fdr < pooled_fdr_cutoff, na.rm = TRUE)),
  paste0("end=", format(Sys.time(), tz = "Europe/London")),
  "status=complete"
)
writeLines(run_log, file.path(output_tiers[["logs"]], "Auto_visium_hd_spatial_cancer_tme_interactions_run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_tiers[["logs"]], "Auto_visium_hd_spatial_cancer_tme_interactions_session_info.txt"))
####################
