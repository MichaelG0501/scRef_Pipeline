####################
# Auto_states_cluster.R
# Approach 1: Cluster-based cell state assignment via Louvain community detection
# on residualised, Z-normalised MP UCell scores.
#
# Input:  ref_outs/EAC_Ref_epi.rds, geneNMF_metaprograms_nMP_19.rds, UCell_nMP19_filtered.rds
# Output: ref_outs/Auto_cluster_states.rds          (named vector: barcode -> state)
#         ref_outs/Auto_cluster_states_umap.pdf
#         ref_outs/Auto_cluster_states_mean_heatmap.pdf
#         ref_outs/Auto_cluster_states_singlecell_heatmap.pdf
#         ref_outs/Auto_cluster_states_composition.pdf
#         ref_outs/Auto_cluster_states_silhouette.pdf
####################
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cluster)
library(scales)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

cat("=== Approach 1: Cluster-Based State Assignment ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

# ============================================================================
# 1. Load data
# ============================================================================
cat("Loading data...\n")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")

# ============================================================================
# 2. Silhouette filtering
# ============================================================================
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  cat("Removing low-quality MPs (silhouette < 0):", bad_mp_names, "\n")
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)
cat("Retained MPs:", retained_mps, "\n")

# Derive canonical MP ordering from geneNMF program tree
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)
cat("Canonical MP tree order:", mp_tree_order_names, "\n")

# ============================================================================
# 3. MP annotations and cell-cycle identification
# ============================================================================
mp_descriptions <- c(
  "MP1"  = "G2M_cycle",
  "MP2"  = "MYC_prolif",
  "MP5"  = "IFN_response",
  "MP7"  = "S_cycle",
  "MP8"  = "Intestinal_diff",
  "MP9"  = "G1S_cycle",
  "MP10" = "Columnar_diff",
  "MP12" = "Neuro_epithelial",
  "MP13" = "Partial_EMT",
  "MP14" = "Hypoxia_epithelial",
  "MP15" = "T_NK_infiltration",
  "MP16" = "Secretory_diff",
  "MP17" = "Squamous_transition",
  "MP18" = "Adaptive_secretory"
)

cc_mps     <- c("MP1", "MP7", "MP9")
non_cc_mps <- setdiff(retained_mps, cc_mps)
cat("Cell-cycle MPs:", cc_mps, "\n")
cat("Non-CC MPs for state definition:", non_cc_mps, "\n")

# ============================================================================
# 4. Align UCell scores to Seurat object
# ============================================================================
common_cells <- intersect(rownames(ucell_scores), Cells(tmdata_all))
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
tmdata_all   <- tmdata_all[, common_cells]
cat(sprintf("Aligned %d cells\n", length(common_cells)))

# Ensure all retained MPs are present in UCell scores
retained_in_ucell <- intersect(retained_mps, colnames(ucell_scores))
cc_in_ucell       <- intersect(cc_mps, colnames(ucell_scores))
non_cc_in_ucell   <- intersect(non_cc_mps, colnames(ucell_scores))
cat("Non-CC MPs available in UCell:", non_cc_in_ucell, "\n")

# ============================================================================
# 5. Regress cell-cycle signal from non-CC MPs
# ============================================================================
cat("Regressing cell-cycle signal...\n")
ucell_mat <- as.matrix(ucell_scores[, retained_in_ucell, drop = FALSE])

X_cc    <- ucell_mat[, cc_in_ucell, drop = FALSE]
Y_other <- ucell_mat[, non_cc_in_ucell, drop = FALSE]

# Efficient OLS residualization: Y_resid = Y - X (X'X)^-1 X'Y
X <- cbind(Intercept = 1, X_cc)
XtX_inv <- solve(crossprod(X))
B       <- XtX_inv %*% crossprod(X, Y_other)
Y_hat   <- X %*% B
Y_resid <- Y_other - Y_hat

cat("Residualization complete.\n")

# ============================================================================
# 6. Z-score normalisation: centre per sample, scale by within-study SD
# ============================================================================
cat("Z-normalising: centre per sample, scale by within-study SD...\n")
sample_var <- tmdata_all$orig.ident
study_var  <- tmdata_all$study

clust_df <- as.data.frame(Y_resid)
clust_df$.cell   <- rownames(Y_resid)
clust_df$.sample <- sample_var
clust_df$.study  <- study_var

# Compute per-study SD for each MP
study_sd <- clust_df %>%
  group_by(.study) %>%
  summarise(across(all_of(non_cc_in_ucell), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames(".study") %>%
  as.matrix()

# Replace 0 / NA SDs with 1 to avoid Inf
study_sd[is.na(study_sd) | study_sd == 0] <- 1

# Centre within sample
clust_centered <- clust_df %>%
  group_by(.sample) %>%
  mutate(across(all_of(non_cc_in_ucell), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Scale by within-study SD
mp_adj <- as.matrix(clust_centered[, non_cc_in_ucell])
rownames(mp_adj) <- clust_centered$.cell

for (mp in non_cc_in_ucell) {
  cell_studies <- clust_centered$.study
  mp_adj[, mp] <- mp_adj[, mp] / study_sd[cell_studies, mp]
}
mp_adj[!is.finite(mp_adj)] <- 0

cat("Z-normalisation complete.\n")

# ============================================================================
# 7. Create Seurat MP assay and cluster
# ============================================================================
cat("Creating MP assay and running PCA + Louvain...\n")

# MPs x cells for Seurat assay
mp_matrix <- t(mp_adj)

tmdata_all[["MPs"]] <- CreateAssayObject(data = mp_matrix)
DefaultAssay(tmdata_all) <- "MPs"

n_pcs <- min(30, nrow(mp_matrix) - 1)
tmdata_all <- ScaleData(tmdata_all, features = rownames(tmdata_all), verbose = FALSE)
tmdata_all <- RunPCA(tmdata_all, features = rownames(tmdata_all),
                     npcs = n_pcs, verbose = FALSE)
tmdata_all <- FindNeighbors(tmdata_all, reduction = "pca", dims = 1:n_pcs,
                            graph.name = "MPs_snn", verbose = FALSE)

# Store multiple resolutions for comparison
for (res in c(0.5, 0.8, 1.0)) {
  tmdata_all <- FindClusters(tmdata_all, graph.name = "MPs_snn",
                             resolution = res, verbose = FALSE)
  col_name <- paste0("MP_state_res", gsub("\\.", "", as.character(res)))
  tmdata_all@meta.data[[col_name]] <-
    paste0("State_", as.numeric(as.character(tmdata_all$seurat_clusters)) + 1)
}

# Default: resolution 0.8
tmdata_all$cluster_state <- tmdata_all@meta.data[["MP_state_res08"]]
cat(sprintf("\nResolution 0.8 → %d states\n", length(unique(tmdata_all$cluster_state))))
print(table(tmdata_all$cluster_state))

# ============================================================================
# 8. UMAP from MP PCA embedding
# ============================================================================
cat("Running UMAP...\n")
tmdata_all <- RunUMAP(tmdata_all, reduction = "pca", dims = 1:n_pcs,
                      reduction.name = "umap_mp", verbose = FALSE)

# ============================================================================
# 9. Visualisation: UMAP (state + study side-by-side)
# ============================================================================
cat("Generating UMAP plot...\n")
state_names <- sort(unique(tmdata_all$cluster_state))
state_cols  <- setNames(hue_pal()(length(state_names)), state_names)
study_cols  <- setNames(
  DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"),
  unique(tmdata_all$study)
)

p1 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "cluster_state",
              cols = state_cols, label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Cluster-based states (res=0.8)")
p2 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "study",
              cols = study_cols, pt.size = 0.1) +
  ggtitle("Coloured by study")

pdf("Auto_cluster_states_umap.pdf", width = 20, height = 8, useDingbats = FALSE)
print(p1 | p2)
dev.off()
cat("Saved: Auto_cluster_states_umap.pdf\n")

# ============================================================================
# 10. Mean Z-normalised MP score per state (ComplexHeatmap)
# ============================================================================
cat("Generating mean heatmap...\n")

meta <- tmdata_all@meta.data
mean_scores <- as.data.frame(mp_adj[rownames(meta), , drop = FALSE])
mean_scores$state <- meta$cluster_state

mean_mat <- mean_scores %>%
  group_by(state) %>%
  summarise(across(all_of(non_cc_in_ucell), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("state") %>%
  as.matrix()

# Reorder columns by canonical MP tree order
non_cc_tree_order <- mp_tree_order_names[mp_tree_order_names %in% non_cc_in_ucell]
mean_mat <- mean_mat[, non_cc_tree_order, drop = FALSE]
colnames(mean_mat) <- mp_descriptions[colnames(mean_mat)]

# Cell counts per state for annotation
cell_counts <- table(meta$cluster_state)
row_anno <- rowAnnotation(
  Cells = anno_barplot(as.numeric(cell_counts[rownames(mean_mat)]),
                       gp = gpar(fill = "steelblue"),
                       width = unit(2, "cm"))
)

col_fun_z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mean <- Heatmap(
  mean_mat,
  name = "Mean Z-score",
  col = col_fun_z,
  cluster_rows = TRUE, cluster_columns = FALSE,
  clustering_method_rows = "ward.D2",
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_rot = 45,
  column_title = "Mean Z-normalised MP score per cluster state",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mean_mat[i, j]),
              x, y, gp = gpar(fontsize = 7))
  }
)

pdf("Auto_cluster_states_mean_heatmap.pdf", width = 14, height = 8, useDingbats = FALSE)
draw(ht_mean, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
cat("Saved: Auto_cluster_states_mean_heatmap.pdf\n")

# ============================================================================
# 11. Single-cell heatmap (downsampled, split by state)
# ============================================================================
cat("Building single-cell heatmap...\n")
set.seed(42)
MAX_CELLS_TOTAL <- 8000

# Proportional downsampling: preserve relative state sizes
state_counts <- table(meta$cluster_state)
state_fracs  <- state_counts / sum(state_counts)
cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)  # min 20 per state

cells_to_plot <- unlist(
  mapply(function(cells, n) sample(cells, min(length(cells), n)),
         split(rownames(meta), meta$cluster_state),
         cells_per_state[names(split(rownames(meta), meta$cluster_state))],
         SIMPLIFY = FALSE),
  use.names = FALSE
)

# Reorder MP rows by canonical tree order
sub_scores <- t(mp_adj[cells_to_plot, , drop = FALSE])
mp_row_order <- mp_tree_order_names[mp_tree_order_names %in% rownames(sub_scores)]
sub_scores <- sub_scores[mp_row_order, , drop = FALSE]
rownames(sub_scores) <- mp_descriptions[rownames(sub_scores)]

sub_meta  <- meta[cells_to_plot, ]
split_vec <- factor(sub_meta$cluster_state, levels = sort(unique(sub_meta$cluster_state)))
study_vals <- sub_meta$study

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

col_ann <- HeatmapAnnotation(
  State = split_vec,
  Study = study_vals,
  col = list(
    State = state_cols,
    Study = study_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht_sc <- Heatmap(
  sub_scores, name = "Adj score", col = col_fun_sc,
  top_annotation = col_ann,
  column_split = split_vec, column_gap = unit(1.5, "mm"),
  cluster_rows = FALSE, cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  show_row_dend = FALSE, row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  show_column_names = FALSE,
  use_raster = TRUE, raster_quality = 5,
  border = FALSE, rect_gp = gpar(col = NA),
  column_title = "Single-cell MP heatmap (cluster states, proportional sampling)"
)

pdf("Auto_cluster_states_singlecell_heatmap.pdf", width = 18,
    height = max(7, length(non_cc_in_ucell) * 0.5), useDingbats = FALSE)
draw(ht_sc, merge_legend = TRUE)
dev.off()
cat("Saved: Auto_cluster_states_singlecell_heatmap.pdf\n")

# ============================================================================
# 12. Study composition barplot per state
# ============================================================================
cat("Generating composition barplot...\n")
comp_df <- meta %>%
  count(study, cluster_state) %>%
  group_by(study) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(x = study, y = pct, fill = cluster_state)) +
  geom_col(colour = "black", linewidth = 0.2) +
  labs(title = "Cluster state proportions per study",
       x = "Study", y = "% of cells", fill = "State") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Auto_cluster_states_proportion.pdf", p_comp, width = 12, height = 7)
cat("Saved: Auto_cluster_states_composition.pdf\n")

# ============================================================================
# 13. Silhouette plot
# ============================================================================
cat("Computing silhouette...\n")
# Use PCA embeddings for distance
pca_embeddings <- Embeddings(tmdata_all, reduction = "pca")[, 1:n_pcs]
cluster_ids    <- as.integer(factor(tmdata_all$cluster_state))

# Subsample for silhouette (full 75k is expensive)
set.seed(42)
n_sil <- min(10000, ncol(tmdata_all))
sil_idx <- sample(seq_len(ncol(tmdata_all)), n_sil)

dist_mat <- dist(pca_embeddings[sil_idx, ])
sil_res  <- silhouette(cluster_ids[sil_idx], dist_mat)

pdf("Auto_cluster_states_silhouette.pdf", width = 10, height = 8, useDingbats = FALSE)
plot(sil_res, col = state_cols[sort(unique(tmdata_all$cluster_state))],
     main = "Silhouette plot (cluster states, 10k subsample)",
     border = NA)
dev.off()
cat("Saved: Auto_cluster_states_silhouette.pdf\n")

# ============================================================================
# 14. Save state assignments
# ============================================================================
state_vec <- setNames(tmdata_all$cluster_state, Cells(tmdata_all))
saveRDS(state_vec, "Auto_cluster_states.rds")
cat("Saved: Auto_cluster_states.rds\n")

# Also save the UMAP embeddings for comparison script
umap_embeddings <- Embeddings(tmdata_all, reduction = "umap_mp")
saveRDS(umap_embeddings, "Auto_cluster_umap_embeddings.rds")
cat("Saved: Auto_cluster_umap_embeddings.rds\n")

# Save the normalised MP matrix for comparison script
saveRDS(mp_adj, "Auto_cluster_mp_adj.rds")
cat("Saved: Auto_cluster_mp_adj.rds\n")

cat("=== Approach 1 complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
