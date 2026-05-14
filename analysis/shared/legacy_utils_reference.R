####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/shared/legacy_utils_reference.R
#   Methodology: analysis/methodology/README.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

# analysis/utils.R
# ==============================================================================
# REFERENCE FILE — NOT SOURCED BY ANY SCRIPT
#
# This file documents the canonical patterns used repeatedly across analysis/
# scripts. All scripts are self-contained (zero inter-script source() calls).
# This file exists for human reference and as a starting point if future
# refactoring adds shared sourcing.
#
# Patterns documented here:
#   1. filter_silhouette()     — used in 7+ scripts
#   2. score_ucell()           — used in 11+ scripts
#   3. load_epi_data()         — used in 8+ scripts
#   4. plot_correlation_heatmap() — used in 3+ scripts
# ==============================================================================

# ==============================================================================
# Utility Functions
# ==============================================================================

#' Filter metaprograms by silhouette score
#' @description Remove MPs with silhouette < threshold from a gene list.
#'   Used in: metaprograms/mp_ucell_scoring.R, metaprograms/mp_database_correlation.R,
#'            enrichment/enrichment_annotation.R, and 4 other scripts.
#' @param mp_genes Named list of MP gene vectors (e.g. geneNMF_obj$metaprograms.genes)
#' @param sil_scores Named numeric vector of silhouette scores (e.g. geneNMF_obj$metaprograms.metrics$silhouette)
#' @param threshold Minimum silhouette score to KEEP (default 0)
#' @return Filtered mp_genes list with low-quality MPs removed
#' @example
#'   bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
#'   bad_mp_names <- paste0("MP", bad_mps)
#'   mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
filter_silhouette <- function(mp_genes, sil_scores, threshold = 0) {
  bad_idx   <- which(sil_scores < threshold)
  bad_names <- names(sil_scores)[bad_idx]
  if (length(bad_names) == 0) bad_names <- paste0("MP", bad_idx)
  cat("Removing low-quality MPs (silhouette <", threshold, "):", bad_names, "\n")
  mp_genes[!names(mp_genes) %in% bad_names]
}

#' Plot MP correlation heatmap
#' @description Build a ComplexHeatmap for an MP x MP (or MP x term) correlation matrix.
#'   Used in: metaprograms/mp_correlation_sc.R, metaprograms/mp_correlation_pdo.R,
#'            metaprograms/mp_database_correlation.R.
#' @param cor_matrix Numeric matrix (MPs x MPs or MPs x reference terms)
#' @param col_fun circlize::colorRamp2 colour function (default: blue-white-red)
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#' @return ComplexHeatmap Heatmap object
#' @example
#'   col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#'   ht <- ComplexHeatmap::Heatmap(cor_mat, col = col_fun, name = "Spearman r",
#'           cluster_rows = TRUE, cluster_columns = TRUE,
#'           show_row_names = TRUE, show_column_names = TRUE)
plot_correlation_heatmap <- function(cor_matrix,
                                     col_fun = circlize::colorRamp2(c(-1, 0, 1),
                                                c("blue", "white", "red")),
                                     ...) {
  ComplexHeatmap::Heatmap(cor_matrix, col = col_fun, name = "Spearman r",
                          cluster_rows = TRUE, cluster_columns = TRUE,
                          show_row_names = TRUE, show_column_names = TRUE, ...)
}

#' UCell scoring wrapper
#' @description Score cells in a Seurat object against named gene lists using UCell.
#'   Used in: metaprograms/mp_ucell_scoring.R, metaprograms/mp_database_correlation.R,
#'            metaprograms/mp_external_scoring.R, and 8 other scripts.
#' @param seurat_obj Seurat object
#' @param gene_lists Named list of gene vectors (MP gene sets)
#' @param ncores Number of cores (default 2)
#' @return Matrix of UCell scores (cells x gene_lists)
#' @example
#'   tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp.genes, ncores = 2, name = "")
#'   ucell_scores <- tmdata_all@meta.data[, names(mp.genes), drop = FALSE]
score_ucell <- function(seurat_obj, gene_lists, ncores = 2) {
  seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, features = gene_lists,
                                            ncores = ncores, name = "")
  seurat_obj@meta.data[, names(gene_lists), drop = FALSE]
}

#' Load EAC epithelial Seurat object
#' @description Load the main epithelial Seurat object (75,348 OAC cells).
#'   Used in: metaprograms/mp_ucell_scoring.R, metaprograms/mp_correlation_sc.R,
#'            metaprograms/mp_database_correlation.R, and 5 other scripts.
#' @note All scripts call setwd("ref_outs/") before this, so the relative path works.
#'   Full absolute path: ref_outs/EAC_Ref_epi.rds
#' @param path Path to EAC_Ref_epi.rds (default: relative path after setwd to ref_outs/)
#' @return Seurat object (75,348 OAC epithelial cells, orig.ident = sample IDs)
#' @example
#'   setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
#'   tmdata_all <- readRDS("EAC_Ref_epi.rds")
load_epi_data <- function(path = "EAC_Ref_epi.rds") {
  readRDS(path)
}
