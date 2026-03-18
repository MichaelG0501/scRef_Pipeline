####################
# Auto_states_topmp.R
# Approach 2: Top-MP residual state assignment — each cell assigned to its
# dominant (highest Z-score) non-cell-cycle metaprogramme.
#
# Input:  ref_outs/EAC_Ref_epi.rds, geneNMF_metaprograms_nMP_19.rds, UCell_nMP19_filtered.rds
# Output: ref_outs/Auto_topmp_states.rds                (named vector: barcode -> state)
#         ref_outs/Auto_topmp_states_proportion.pdf
#         ref_outs/Auto_topmp_states_mean_heatmap.pdf
#         ref_outs/Auto_topmp_states_singlecell_heatmap.pdf
#         ref_outs/Auto_topmp_states_violin.pdf
####################
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

cat("=== Approach 2: Top-MP State Assignment ===\n")
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
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC Proliferation",
  "MP17" = "Basal-like Trans.",
  "MP14" = "Hypoxia Adapted",
  "MP5"  = "Epithelial IFN",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Intest.",
  "MP16" = "Secretory Gastric",
  "MP15" = "Immune Infilt.",
  "MP12" = "Neuro-responsive"
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

# Efficient OLS residualization
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
# 7. Assign states based on dominant MP
# ============================================================================
cat("Assigning states based on dominant MP...\n")

# For each cell, identify the MP with the highest Z-normalised score
max_mp_idx   <- max.col(mp_adj, ties.method = "first")
max_mp_names <- non_cc_in_ucell[max_mp_idx]
max_mp_vals  <- apply(mp_adj, 1, max)

# Map MP names to readable descriptions
state_labels <- mp_descriptions[max_mp_names]

# Apply minimum-signal threshold
THRESHOLD <- 0.5
state_labels[max_mp_vals < THRESHOLD] <- "Unresolved"

names(state_labels) <- rownames(mp_adj)
tmdata_all$topmp_state <- state_labels

cat("State distribution:\n")
print(sort(table(state_labels), decreasing = TRUE))
cat(sprintf("Unresolved cells: %d (%.1f%%)\n",
            sum(state_labels == "Unresolved"),
            100 * mean(state_labels == "Unresolved")))

# ============================================================================
# 8. Colour palettes
# ============================================================================
all_states <- sort(unique(state_labels))
non_unresolved <- setdiff(all_states, "Unresolved")

# Assign colours: use hue_pal for defined states, grey for Unresolved
defined_cols <- setNames(hue_pal()(length(non_unresolved)), non_unresolved)
state_cols <- c(defined_cols, "Unresolved" = "grey80")

study_cols <- setNames(
  DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"),
  unique(tmdata_all$study)
)

# ============================================================================
# 9. Proportion barplot (overall + per-study)
# ============================================================================
cat("Generating proportion barplot...\n")
meta <- tmdata_all@meta.data

# Overall proportions
overall_df <- data.frame(state = state_labels) %>%
  count(state) %>%
  mutate(pct = 100 * n / sum(n),
         study = "Overall")

# Per-study proportions
per_study_df <- data.frame(state = state_labels, study = meta$study) %>%
  count(state, study) %>%
  group_by(study) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

plot_df <- bind_rows(overall_df, per_study_df)
plot_df$study <- factor(plot_df$study, levels = c("Overall", sort(unique(meta$study))))

p_prop <- ggplot(plot_df, aes(x = study, y = pct, fill = state)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = state_cols) +
  labs(title = "Cell state proportions (Top-MP assignment)",
       x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Auto_topmp_states_proportion.pdf", p_prop, width = 14, height = 7)
cat("Saved: Auto_topmp_states_proportion.pdf\n")

# ============================================================================
# 10. Mean MP activity per state (ComplexHeatmap)
# ============================================================================
cat("Generating mean heatmap...\n")

mean_scores_df <- as.data.frame(mp_adj)
mean_scores_df$state <- state_labels

mean_mat <- mean_scores_df %>%
  group_by(state) %>%
  summarise(across(all_of(non_cc_in_ucell), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("state") %>%
  as.matrix()
# Reorder columns by canonical MP tree order
non_cc_tree_order <- mp_tree_order_names[mp_tree_order_names %in% non_cc_in_ucell]
mean_mat <- mean_mat[, non_cc_tree_order, drop = FALSE]
colnames(mean_mat) <- mp_descriptions[colnames(mean_mat)]

# Cell counts for annotation
cell_counts <- table(state_labels)
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
  column_title = "Mean Z-normalised MP score per top-MP state",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mean_mat[i, j]),
              x, y, gp = gpar(fontsize = 7))
  }
)

pdf("Auto_topmp_states_mean_heatmap.pdf", width = 14, height = 8, useDingbats = FALSE)
draw(ht_mean, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
cat("Saved: Auto_topmp_states_mean_heatmap.pdf\n")

# ============================================================================
# 11. Single-cell heatmap (downsampled, split by state)
# ============================================================================
cat("Building single-cell heatmap...\n")
set.seed(42)
MAX_CELLS_TOTAL <- 8000

# Proportional downsampling: preserve relative state sizes
state_counts <- table(state_labels)
state_fracs  <- state_counts / sum(state_counts)
cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)  # min 20 per state

cells_to_plot <- unlist(
  mapply(function(cells, n) sample(cells, min(length(cells), n)),
         split(rownames(meta), state_labels),
         cells_per_state[names(split(rownames(meta), state_labels))],
         SIMPLIFY = FALSE),
  use.names = FALSE
)

# Reorder MP rows by canonical tree order
sub_scores <- t(mp_adj[cells_to_plot, , drop = FALSE])
mp_row_order <- mp_tree_order_names[mp_tree_order_names %in% rownames(sub_scores)]
sub_scores <- sub_scores[mp_row_order, , drop = FALSE]
rownames(sub_scores) <- mp_descriptions[rownames(sub_scores)]

sub_meta  <- meta[cells_to_plot, ]
split_vec <- factor(state_labels[cells_to_plot],
                    levels = sort(unique(state_labels)))
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
  column_title = "Single-cell MP heatmap (top-MP states, proportional sampling)"
)

pdf("Auto_topmp_states_singlecell_heatmap.pdf", width = 18,
    height = max(7, length(non_cc_in_ucell) * 0.5), useDingbats = FALSE)
draw(ht_sc, merge_legend = TRUE)
dev.off()
cat("Saved: Auto_topmp_states_singlecell_heatmap.pdf\n")

# ============================================================================
# 12. Violin plot of dominant MP Z-score per state
# ============================================================================
cat("Generating violin plot...\n")

violin_df <- data.frame(
  state = state_labels,
  max_zscore = max_mp_vals,
  stringsAsFactors = FALSE
)

p_violin <- ggplot(violin_df, aes(x = state, y = max_zscore, fill = state)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  scale_fill_manual(values = state_cols) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed", colour = "red") +
  annotate("text", x = 0.5, y = THRESHOLD + 0.15, label = paste0("Threshold = ", THRESHOLD),
           hjust = 0, colour = "red", size = 3) +
  labs(title = "Dominant MP Z-score distribution per state",
       subtitle = "Confirms state separation quality",
       x = "State", y = "Max Z-score (dominant MP)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("Auto_topmp_states_violin.pdf", p_violin, width = 12, height = 7)
cat("Saved: Auto_topmp_states_violin.pdf\n")

# ============================================================================
# 13. Save state assignments
# ============================================================================
saveRDS(state_labels, "Auto_topmp_states.rds")
cat("Saved: Auto_topmp_states.rds\n")

# Save the normalised MP matrix for comparison script
saveRDS(mp_adj, "Auto_topmp_mp_adj.rds")
cat("Saved: Auto_topmp_mp_adj.rds\n")

cat("=== Approach 2 complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
