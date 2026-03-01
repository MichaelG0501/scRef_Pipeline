####################
# Moved from: analysis/Auto_MP19_analysis.R
# Reorganized as part of analysis/ restructuring
####################
####################
# Auto_MP19_analysis.R
# MP19 UCell scoring, silhouette filtering, correlation heatmap,
# Jaccard self-similarity heatmap, and ComplexHeatmap UCell score heatmap.
# Uses same patterns as geneNMF.R, MP_analysis_sc.R, compare_pdos_sc.R
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ============================================================================
# 1. Load data
# ============================================================================

geneNMF.metaprograms <- readRDS(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)
mp.genes <- geneNMF.metaprograms$metaprograms.genes

tmdata_all <- readRDS("EAC_Ref_epi.rds")

# ============================================================================
# 2. Filter MPs by silhouette score < 0
# ============================================================================

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  cat("Removing low-quality MPs (silhouette < 0):", bad_mp_names, "\n")
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
cat("Retained MPs:", names(mp.genes), "\n")

# ============================================================================
# 3. UCell scoring (similar to geneNMF.R line 64)
# ============================================================================

tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp.genes, ncores = 2, name = "")
ucell_scores <- tmdata_all@meta.data[, names(mp.genes), drop = FALSE]

saveRDS(ucell_scores, file = "Metaprogrammes_Results/UCell_nMP19_filtered.rds")
cat("Saved UCell scores: Metaprogrammes_Results/UCell_nMP19_filtered.rds\n")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")

# ============================================================================
# 4. Define mp_tree_order (after bad MP removal)
# ============================================================================

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)

# Keep only valid (non-bad, non-NA) cluster IDs
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp.genes)))
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

cat("MP tree order (after filtering):", mp_tree_order_names, "\n")

# ============================================================================
# 5. Name mapping (placeholder — fill in when annotations are done)
# ============================================================================

mp_descriptions <- setNames(names(mp.genes), names(mp.genes))
# e.g. mp_descriptions <- c(
#   "MP1"  = "MP1",
#   "MP2"  = "MP2",
#   "MP3"  = "MP3",
#   ...
# )
# Will be updated with final annotations

# ============================================================================
# 6. Scale UCell scores and prepare matrix
# ============================================================================

module_scores <- scale(as.matrix(ucell_scores))

# Order rows by mp_tree_order
mod_mat <- t(module_scores[, mp_tree_order_names, drop = FALSE])  # MPs x Cells

# Apply name mapping for display
rownames(mod_mat) <- mp_descriptions[rownames(mod_mat)]

# ============================================================================
# 7. Correlation Heatmap — per-sample Fisher Z meta-analysis
#    (similar to MP_analysis_sc.R lines 190-272, grouped by orig.ident)
# ============================================================================

# Extract sample name by removing last "_" + barcode
samples_vec <- sub("_[^_]+$", "", colnames(mod_mat))
samples <- unique(samples_vec)
mps <- rownames(mod_mat)
n_mps <- length(mps)

cat("Computing per-sample correlations across", length(samples), "samples...\n")

# Create 3D array: [MP x MP x Sample]
cor_array <- array(NA, dim = c(n_mps, n_mps, length(samples)),
                   dimnames = list(mps, mps, samples))

for (samp in samples) {
  cells_in_sample <- colnames(mod_mat)[samples_vec == samp]
  if (length(cells_in_sample) < 3) next  # skip tiny samples
  sub_mat <- mod_mat[, cells_in_sample, drop = FALSE]
  cor_array[, , samp] <- cor(t(sub_mat), method = "spearman")
}
# Fisher Z-transformation for averaging and t-testing
z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

mean_rho <- matrix(0, n_mps, n_mps, dimnames = list(mps, mps))
p_vals   <- matrix(1, n_mps, n_mps, dimnames = list(mps, mps))

for (i in 1:n_mps) {
  for (j in 1:n_mps) {
    if (i == j) {
      mean_rho[i, j] <- 1
      p_vals[i, j] <- 0
    } else {
      z_scores <- z_array[i, j, ]
      # Remove NaN from tiny/missing samples
      z_scores <- z_scores[is.finite(z_scores)]
      if (length(z_scores) < 3) {
        mean_rho[i, j] <- NA
        p_vals[i, j] <- NA
      } else {
        test_res <- t.test(z_scores)
        mean_rho[i, j] <- tanh(mean(z_scores))
        p_vals[i, j] <- test_res$p.value
      }
    }
  }
}

# Plot correlation heatmap
col_cor <- colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))

Heatmap(mean_rho,
        name = paste0("Mean Rho\n(", length(samples), " Samples)"),
        col = col_cor,
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,

        # Overlay significance stars
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- p_vals[i, j]
          rho <- mean_rho[i, j]
          if (is.na(p) || is.na(rho)) {
            grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "grey50"))
          } else if (p < 0.001) {
            grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 10))
          } else if (p < 0.01) {
            grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 10))
          } else if (p < 0.05) {
            grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 10))
          } else {
            grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 10))
          }
        },

        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold")))
cat("Saved: Metaprogrammes_Results/nMP19_correlation_heatmap_persample.pdf\n")

# ============================================================================
# 8. Jaccard Self-Similarity Heatmap
#    (similar to compare_pdos_sc.R lines 111-257, but self-comparison)
# ============================================================================

# Order gene lists by tree order
mp_list <- mp.genes[mp_tree_order_names]
mp_list <- lapply(mp_list, unique)

# Rename using mp_descriptions
names(mp_list) <- mp_descriptions[names(mp_list)]

mp_names <- names(mp_list)
universe <- unique(unlist(mp_list))

# Initialize matrices
jaccard_mat   <- matrix(NA_real_, length(mp_list), length(mp_list),
                        dimnames = list(mp_names, mp_names))
overlap_n_mat <- jaccard_mat
pval_mat      <- jaccard_mat

# Compute Jaccard, overlap counts, Fisher p-values
for (i in seq_along(mp_list)) {
  A <- mp_list[[i]]
  for (j in seq_along(mp_list)) {
    B <- mp_list[[j]]

    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))

    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j]   <- if (uni == 0) NA_real_ else inter / uni

    a <- inter
    b <- length(setdiff(A, B))
    cc <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))

    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_
    else fisher.test(matrix(c(a, b, cc, d), nrow = 2),
                     alternative = "greater")$p.value
  }
}

# Adjust p-values
padj_mat <- matrix(
  p.adjust(as.vector(pval_mat), method = "BH"),
  nrow = nrow(pval_mat), ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

# Build stars from adjusted p-values
stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat),
                    dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05]  <- "*"
stars_mat[padj_mat < 0.01]  <- "**"
stars_mat[padj_mat < 0.001] <- "***"

# Combine overlap count + stars for display
display_mat <- matrix(
  paste0(overlap_n_mat, "\n", stars_mat),
  nrow = nrow(overlap_n_mat),
  ncol = ncol(overlap_n_mat),
  dimnames = dimnames(overlap_n_mat)
)

pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "MP Gene Set Overlap (Jaccard Index) - nMP19",
  labels_row = rownames(jaccard_mat),
  labels_col = colnames(jaccard_mat),
  angle_col = "90",
  display_numbers = display_mat,
  fontsize_number = 8,
  number_color = "black",
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
)
cat("Saved: Metaprogrammes_Results/nMP19_jaccard_self_similarity_heatmap.pdf\n")

# ============================================================================
# 9. ComplexHeatmap — UCell Score Heatmap (subsampled, 20k cells)
#    (similar to MP_analysis_sc.R lines 130-188, cluster_rows = FALSE,
#     cluster_columns = TRUE)
# ============================================================================

set.seed(42)
n_sub <- min(20000, ncol(tmdata_all))
sub_idx <- sample(seq_len(ncol(tmdata_all)), size = n_sub)
tmdata_sub <- tmdata_all[, sub_idx]
cat("Subsampled to", ncol(tmdata_sub), "cells for heatmap\n")

# Subset module scores to match
module_scores_sub <- module_scores[Cells(tmdata_sub), mp_tree_order_names, drop = FALSE]

# --- kmeans clustering for MP_state (on subset) ---
K <- 300  # reduced from 500 proportional to 20k vs 75k
set.seed(1)
km <- kmeans(module_scores_sub, centers = K, nstart = 10, iter.max = 200,
             algorithm = "MacQueen")

clus <- factor(km$cluster)
centroids <- km$centers

# Hierarchical clustering of centroids
d <- as.dist(1 - cor(t(centroids), method = "pearson"))
hc <- hclust(d, method = "average")

# Determine optimal k for MP states using silhouette
library(factoextra)
library(cluster)

p1 <- fviz_nbclust(centroids, hcut, hc_method = "average", hc_metric = "pearson",
                    method = "silhouette", k.max = 20) +
  labs(title = "Silhouette Analysis (Centroids) - 20k subset")
p2 <- fviz_nbclust(centroids, hcut, hc_method = "average", hc_metric = "pearson",
                    method = "wss", k.max = 20) +
  labs(title = "Elbow Method (Centroids) - 20k subset")
library(patchwork)
print(p1 + p2)
cat("Saved: Metaprogrammes_Results/nMP19_subset_elbow_silhouette.pdf\n")
cat("CHECK the silhouette/elbow plot and set k below accordingly.\n")

# ---- SET k HERE after inspecting the elbow/silhouette plot ----
k <- 10  # <-- ADJUST based on silhouette/elbow result

mc_state <- cutree(hc, k = k)
names(mc_state) <- rownames(centroids)

tmdata_sub$km_cluster <- as.character(km$cluster)
tmdata_sub$MP_state <- as.character(mc_state[tmdata_sub$km_cluster])

# --- Cell ordering within states ---
mc_order <- rownames(centroids)[hc$order]
cell_scores_sub <- module_scores_sub
centroid_mat <- centroids[tmdata_sub$km_cluster, , drop = FALSE]
cell_to_centroid_cor <- vapply(seq_len(nrow(cell_scores_sub)), function(i) {
  cor(cell_scores_sub[i, ], centroid_mat[i, ], method = "pearson")
}, numeric(1))

order_df <- data.frame(
  cell_id = Cells(tmdata_sub),
  state = tmdata_sub$MP_state,
  state_order = match(tmdata_sub$MP_state, unique(mc_state[hc$order])),
  mc_order_val = match(tmdata_sub$km_cluster, mc_order),
  cor_val = cell_to_centroid_cor
)

order_df <- order_df[order(order_df$state_order, order_df$mc_order_val, -order_df$cor_val), ]
cell_ord <- match(order_df$cell_id, Cells(tmdata_sub))
tmdata_sub <- tmdata_sub[, cell_ord]

# --- Prepare heatmap matrix ---
mod_mat_sub <- t(as.matrix(tmdata_sub@meta.data[, mp_tree_order_names, drop = FALSE]))
rownames(mod_mat_sub) <- mp_descriptions[rownames(mod_mat_sub)]

# Scale per-row (z-score across cells for each MP)
mod_mat_sub <- t(scale(t(mod_mat_sub)))

# Color scale
max_val <- quantile(mod_mat_sub, 0.98, na.rm = TRUE)
col_fun <- colorRamp2(c(-max_val, 0, max_val), c("blue", "white", "red"))

# Annotations
state_names <- as.character(unique(tmdata_sub$MP_state))
state_cols  <- setNames(DiscretePalette(length(state_names), palette = "alphabet"), state_names)
study_cols  <- setNames(DiscretePalette(length(unique(tmdata_sub$study)), palette = "polychrome"),
                        unique(tmdata_sub$study))

col_ann <- HeatmapAnnotation(
  State = tmdata_sub$MP_state,
  Study = tmdata_sub$study,
  col = list(State = state_cols, Study = study_cols),
  annotation_name_side = "left",
  show_legend = TRUE
)

# Plot
ht <- Heatmap(
  mod_mat_sub,
  name = "UCell Score\n(z-scaled)",
  col = col_fun,

  # Grouping
  top_annotation = col_ann,
  column_split = tmdata_sub$MP_state,
  column_gap = unit(0, "mm"),

  # Row/Column Ordering
  cluster_rows = FALSE,             # Keep GeneNMF tree order
  row_order = rownames(mod_mat_sub),
  cluster_columns = TRUE,           # Cluster cells within each MP_state split
  clustering_method_columns = "ward.D2",

  # Aesthetics
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = paste0("MPs UCell scores (20k subset, k = ", k, ")"),
  show_column_names = FALSE,

  # Performance
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf("nMP19_UCellHeatmap_20k_subset.pdf",
    width = 18, height = 10, useDingbats = FALSE)
ht2 <- draw(ht, merge_legend = TRUE)
dev.off()
