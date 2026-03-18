####################
# Auto_states_comparison.R
# Comparison script: evaluate cluster-based vs top-MP state assignment approaches.
# Produces confusion matrix, ARI/NMI, bootstrap stability, study-bias diagnostic,
# and biological coherence (DEGs + fgsea Hallmark enrichment).
#
# Input:  ref_outs/Auto_cluster_states.rds, Auto_topmp_states.rds,
#         Auto_cluster_umap_embeddings.rds, Auto_cluster_mp_adj.rds,
#         EAC_Ref_epi.rds, geneNMF_metaprograms_nMP_19.rds, UCell_nMP19_filtered.rds
# Output: ref_outs/Auto_states_confusion.pdf
#         ref_outs/Auto_states_comparison_umap.pdf
#         ref_outs/Auto_states_comparison_summary.csv
#         ref_outs/Auto_states_coherence_cluster.pdf
#         ref_outs/Auto_states_coherence_topmp.pdf
#         ref_outs/Auto_states_studybias.pdf
####################
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(mclust)
library(cluster)
library(scales)
library(grid)
library(msigdbr)
library(fgsea)
library(parallel)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

cat("=== State Assignment Comparison ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

# ============================================================================
# 1. Load both sets of state assignments and shared data
# ============================================================================
cat("Loading state assignments and shared data...\n")

cluster_states <- readRDS("Auto_cluster_states.rds")
topmp_states   <- readRDS("Auto_topmp_states.rds")
umap_emb       <- readRDS("Auto_cluster_umap_embeddings.rds")
mp_adj         <- readRDS("Auto_cluster_mp_adj.rds")

tmdata_all <- readRDS("EAC_Ref_epi.rds")

# Align cells across all objects
common_cells <- Reduce(intersect, list(
  names(cluster_states),
  names(topmp_states),
  rownames(umap_emb),
  Cells(tmdata_all)
))

cluster_states <- cluster_states[common_cells]
topmp_states   <- topmp_states[common_cells]
umap_emb       <- umap_emb[common_cells, ]
mp_adj         <- mp_adj[common_cells, ]
tmdata_all     <- tmdata_all[, common_cells]

cat(sprintf("Aligned %d cells across all objects\n", length(common_cells)))

# Add states to Seurat metadata
tmdata_all$cluster_state <- cluster_states
tmdata_all$topmp_state   <- topmp_states

# ============================================================================
# 2. Load MP annotations (for context)
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
non_cc_mps <- colnames(mp_adj)

# ============================================================================
# 3. Confusion matrix heatmap
# ============================================================================
cat("Generating confusion matrix heatmap...\n")

conf_table <- table(Cluster = cluster_states, TopMP = topmp_states)

# Column-normalised proportions
conf_prop <- sweep(conf_table, 2, colSums(conf_table), "/")
conf_prop[is.nan(conf_prop)] <- 0

# Color function
col_fun_conf <- colorRamp2(c(0, 0.5, 1), c("white", "orange", "darkred"))

ht_conf <- Heatmap(
  as.matrix(conf_prop),
  name = "Proportion",
  col = col_fun_conf,
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Cluster states vs Top-MP states (column-normalised)",
  row_title = "Cluster approach",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  cell_fun = function(j, i, x, y, width, height, fill) {
    count_val <- as.matrix(conf_table)[i, j]
    prop_val  <- as.matrix(conf_prop)[i, j]
    grid.text(sprintf("%d\n(%.0f%%)", count_val, prop_val * 100),
              x, y, gp = gpar(fontsize = 7))
  }
)

pdf("Auto_states_confusion.pdf", width = 14, height = 10, useDingbats = FALSE)
draw(ht_conf, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()
cat("Saved: Auto_states_confusion.pdf\n")

# ============================================================================
# 4. Concordance metrics: ARI and NMI
# ============================================================================
cat("Computing concordance metrics...\n")

ari_val <- adjustedRandIndex(cluster_states, topmp_states)
cat(sprintf("Adjusted Rand Index: %.4f\n", ari_val))

# NMI implementation (since aricode may not be installed)
compute_nmi <- function(x, y) {
  # Normalised Mutual Information via entropy
  tab <- table(x, y)
  n <- sum(tab)
  p_xy <- tab / n
  p_x  <- rowSums(p_xy)
  p_y  <- colSums(p_xy)

  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))

  mi <- 0
  for (i in seq_along(p_x)) {
    for (j in seq_along(p_y)) {
      if (p_xy[i, j] > 0) {
        mi <- mi + p_xy[i, j] * log(p_xy[i, j] / (p_x[i] * p_y[j]))
      }
    }
  }
  nmi <- mi / sqrt(h_x * h_y)
  return(nmi)
}

nmi_val <- compute_nmi(cluster_states, topmp_states)
cat(sprintf("Normalised Mutual Information: %.4f\n", nmi_val))

# ============================================================================
# 5. Side-by-side UMAP
# ============================================================================
cat("Generating comparison UMAPs...\n")

# Add UMAP embeddings to Seurat object
tmdata_all[["umap_mp"]] <- CreateDimReducObject(
  embeddings = umap_emb, key = "UMAPMP_", assay = "RNA"
)

# Colour palettes
cluster_names <- sort(unique(cluster_states))
cluster_cols  <- setNames(hue_pal()(length(cluster_names)), cluster_names)

topmp_names <- sort(unique(topmp_states))
non_unresolved <- setdiff(topmp_names, "Unresolved")
topmp_cols <- c(
  setNames(hue_pal()(length(non_unresolved)), non_unresolved),
  "Unresolved" = "grey80"
)

p1 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "cluster_state",
              cols = cluster_cols, label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Approach 1: Cluster-based states")
p2 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "topmp_state",
              cols = topmp_cols, label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Approach 2: Top-MP states")

pdf("Auto_states_comparison_umap.pdf", width = 22, height = 8, useDingbats = FALSE)
print(p1 | p2)
dev.off()
cat("Saved: Auto_states_comparison_umap.pdf\n")

# ============================================================================
# 6. Bootstrap stability assessment
# ============================================================================
cat("Running bootstrap stability (100 iterations, 80% subsample)...\n")

set.seed(42)
n_boot   <- 100
frac     <- 0.8
n_cells  <- length(common_cells)
n_sample <- round(frac * n_cells)

# Pre-compute full-data assignments
full_topmp   <- topmp_states
full_cluster <- cluster_states

ari_topmp_boot   <- numeric(n_boot)
ari_cluster_boot <- numeric(n_boot)

for (b in seq_len(n_boot)) {
  if (b %% 10 == 0) cat(sprintf("  Bootstrap %d/%d\n", b, n_boot))

  idx <- sample(seq_len(n_cells), n_sample)
  sub_cells <- common_cells[idx]
  sub_adj   <- mp_adj[sub_cells, , drop = FALSE]

  # Top-MP re-assignment on subsample
  max_idx  <- max.col(sub_adj, ties.method = "first")
  max_vals <- apply(sub_adj, 1, max)
  sub_topmp <- mp_descriptions[non_cc_mps[max_idx]]
  sub_topmp[max_vals < 0.5] <- "Unresolved"

  ari_topmp_boot[b] <- adjustedRandIndex(full_topmp[sub_cells], sub_topmp)

  # Cluster re-assignment on subsample (use PCA + Louvain)
  sub_mat <- t(sub_adj)
  sub_seurat <- CreateSeuratObject(counts = sub_mat, assay = "MPs")
  sub_seurat[["MPs"]] <- CreateAssayObject(data = sub_mat)
  DefaultAssay(sub_seurat) <- "MPs"
  n_pcs <- min(30, nrow(sub_mat) - 1)
  sub_seurat <- ScaleData(sub_seurat, features = rownames(sub_seurat), verbose = FALSE)
  sub_seurat <- RunPCA(sub_seurat, features = rownames(sub_seurat),
                       npcs = n_pcs, verbose = FALSE)
  sub_seurat <- FindNeighbors(sub_seurat, reduction = "pca", dims = 1:n_pcs,
                              graph.name = "MPs_snn", verbose = FALSE)
  sub_seurat <- FindClusters(sub_seurat, graph.name = "MPs_snn",
                             resolution = 0.8, verbose = FALSE)
  sub_cluster <- paste0("State_", as.numeric(as.character(sub_seurat$seurat_clusters)) + 1)
  names(sub_cluster) <- Cells(sub_seurat)

  ari_cluster_boot[b] <- adjustedRandIndex(full_cluster[sub_cells], sub_cluster)
}

cat(sprintf("\nBootstrap ARI (Top-MP):  median=%.3f, IQR=[%.3f, %.3f]\n",
            median(ari_topmp_boot),
            quantile(ari_topmp_boot, 0.25),
            quantile(ari_topmp_boot, 0.75)))

cat(sprintf("Bootstrap ARI (Cluster): median=%.3f, IQR=[%.3f, %.3f]\n",
            median(ari_cluster_boot),
            quantile(ari_cluster_boot, 0.25),
            quantile(ari_cluster_boot, 0.75)))

# ============================================================================
# 7. Study-bias diagnostic: Cramér's V
# ============================================================================
cat("Computing study-bias diagnostic (Cramér's V)...\n")

cramers_v <- function(x, y) {
  tab <- table(x, y)
  chi2 <- chisq.test(tab, simulate.p.value = TRUE, B = 2000)$statistic
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  sqrt(chi2 / (n * (k - 1)))
}

study_var <- tmdata_all$study

cv_cluster <- cramers_v(cluster_states, study_var)
cv_topmp   <- cramers_v(topmp_states, study_var)

cat(sprintf("Cramér's V (Cluster vs study): %.4f\n", cv_cluster))
cat(sprintf("Cramér's V (TopMP vs study):   %.4f\n", cv_topmp))

cv_df <- data.frame(
  Approach = c("Cluster-based", "Top-MP"),
  Cramers_V = c(cv_cluster, cv_topmp)
)

p_cv <- ggplot(cv_df, aes(x = Approach, y = Cramers_V, fill = Approach)) +
  geom_col(colour = "black", width = 0.5) +
  scale_fill_manual(values = c("Cluster-based" = "steelblue", "Top-MP" = "coral")) +
  labs(title = "Study-bias diagnostic: Cramér's V",
       subtitle = "Lower = less study-confounded",
       y = "Cramér's V") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("Auto_states_studybias.pdf", p_cv, width = 6, height = 5)
cat("Saved: Auto_states_studybias.pdf\n")

# ============================================================================
# 8. Biological coherence: DEGs + fgsea Hallmark enrichment
# ============================================================================
cat("Running biological coherence analysis...\n")

# Load Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Function to run DEGs + fgsea for one approach
run_coherence <- function(seurat_obj, state_col, approach_name) {
  cat(sprintf("  %s: FindAllMarkers...\n", approach_name))
  DefaultAssay(seurat_obj) <- "RNA"
  Idents(seurat_obj) <- state_col

  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )

  # Top 50 DEGs per state
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 50, order_by = avg_log2FC) %>%
    ungroup()

  # fgsea per state: use all DEGs ranked by avg_log2FC
  cat(sprintf("  %s: fgsea enrichment...\n", approach_name))
  states <- unique(markers$cluster)
  nes_list <- lapply(states, function(s) {
    state_markers <- markers %>% filter(cluster == s)
    ranked_genes  <- setNames(state_markers$avg_log2FC, state_markers$gene)
    ranked_genes  <- sort(ranked_genes, decreasing = TRUE)

    if (length(ranked_genes) < 5) {
      return(setNames(rep(NA, length(hallmark_list)), names(hallmark_list)))
    }

    ####################
# Fixed: Use scoreType = "pos" to handle all-positive stats
    fgsea_res <- fgsea(
      pathways  = hallmark_list,
      stats     = ranked_genes,
      minSize   = 5,
      maxSize   = 500,
      scoreType = "pos"
    )
####################

    nes <- setNames(fgsea_res$NES, fgsea_res$pathway)
    return(nes)
  })
  names(nes_list) <- states

  # Build NES matrix (pathways x states)
  all_pathways <- unique(unlist(lapply(nes_list, names)))
  nes_mat <- matrix(NA, nrow = length(all_pathways), ncol = length(states),
                    dimnames = list(all_pathways, states))
  for (s in states) {
    matched <- intersect(names(nes_list[[s]]), all_pathways)
    nes_mat[matched, s] <- nes_list[[s]][matched]
  }

  # Filter to pathways significant in at least one state
  nes_mat[is.na(nes_mat)] <- 0

  # Keep top variable pathways
  pathway_var <- apply(nes_mat, 1, var, na.rm = TRUE)
  ####################
  # Fixed: remove NA/NaN variances BEFORE sorting, then cap n_top to the
  # surviving length so seq_len() never overshoots → no subscript error
  pathway_var <- pathway_var[!is.na(pathway_var) & !is.nan(pathway_var)]
  n_top <- min(30L, length(pathway_var))
  if (n_top == 0L) {
    top_pathways <- character(0)
  } else {
    top_pathways <- names(sort(pathway_var, decreasing = TRUE))[seq_len(n_top)]
  }
  ####################
  nes_mat <- nes_mat[top_pathways, , drop = FALSE]

  # Clean pathway names for display
  rownames(nes_mat) <- gsub("^HALLMARK_", "", rownames(nes_mat))
  rownames(nes_mat) <- gsub("_", " ", rownames(nes_mat))

  ####################
  # Guard against empty NES matrix
  if (length(top_pathways) == 0 || ncol(nes_mat) == 0) {
    warning(sprintf("%s: fgsea returned empty NES matrix, skipping coherence heatmap.", approach_name))
    return(list(markers = top_markers, nes_mat = NULL))
  }
####################
  return(list(markers = top_markers, nes_mat = nes_mat))
}

# Run for both approaches
coherence_cluster <- run_coherence(tmdata_all, "cluster_state", "Cluster")
coherence_topmp   <- run_coherence(tmdata_all, "topmp_state", "TopMP")

# Plot NES heatmaps
plot_nes_heatmap <- function(nes_mat, title, filename) {
  max_abs <- max(abs(nes_mat), na.rm = TRUE)
  lim <- min(max_abs, 3)
  col_fun <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

  ht <- Heatmap(
    nes_mat,
    name = "NES",
    col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    column_title = title,
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    use_raster = FALSE
  )

  pdf(filename, width = 12,
      height = max(8, nrow(nes_mat) * 0.3), useDingbats = FALSE)
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
  dev.off()
  cat(sprintf("Saved: %s\n", filename))
}

####################
# Fixed: Guard heatmap plotting against NULL nes_mat
if (!is.null(coherence_cluster$nes_mat)) {
  plot_nes_heatmap(coherence_cluster$nes_mat,
                   "Hallmark NES: Cluster-based states",
                   "Auto_states_coherence_cluster.pdf")
}

if (!is.null(coherence_topmp$nes_mat)) {
  plot_nes_heatmap(coherence_topmp$nes_mat,
                   "Hallmark NES: Top-MP states",
                   "Auto_states_coherence_topmp.pdf")
}
####################

# ============================================================================
# 9. Summary table
# ============================================================================
cat("Writing summary table...\n")

summary_df <- data.frame(
  Metric = c(
    "Number of states",
    "ARI between approaches",
    "NMI between approaches",
    "Bootstrap ARI median (80% subsample)",
    "Bootstrap ARI IQR",
    "Cramer's V (study confounding)",
    "Number of DEGs (top 50/state)",
    "Hallmark pathways enriched (NES>1.5 in any state)"
  ),
  Cluster_approach = c(
    length(unique(cluster_states)),
    sprintf("%.4f", ari_val),
    sprintf("%.4f", nmi_val),
    sprintf("%.3f", median(ari_cluster_boot)),
    sprintf("[%.3f, %.3f]", quantile(ari_cluster_boot, 0.25), quantile(ari_cluster_boot, 0.75)),
    sprintf("%.4f", cv_cluster),
    nrow(coherence_cluster$markers),
    if (!is.null(coherence_cluster$nes_mat)) sum(apply(coherence_cluster$nes_mat, 1, function(x) any(x > 1.5, na.rm=TRUE))) else NA_integer_
  ),
  TopMP_approach = c(
    length(unique(topmp_states)),
    sprintf("%.4f", ari_val),
    sprintf("%.4f", nmi_val),
    sprintf("%.3f", median(ari_topmp_boot)),
    sprintf("[%.3f, %.3f]", quantile(ari_topmp_boot, 0.25), quantile(ari_topmp_boot, 0.75)),
    sprintf("%.4f", cv_topmp),
    nrow(coherence_topmp$markers),
    if (!is.null(coherence_topmp$nes_mat)) sum(apply(coherence_topmp$nes_mat, 1, function(x) any(x > 1.5, na.rm=TRUE))) else NA_integer_
  ),
  Interpretation = c(
    "Cluster: data-driven count; TopMP: one per non-CC MP + Unresolved",
    "Agreement between the two partitions (0=random, 1=identical)",
    "Mutual information normalised by entropy (0=independent, 1=identical)",
    "Stability under 80% subsampling; higher=more stable",
    "Interquartile range of bootstrap ARI",
    "Lower = less study-confounded; 0=independent, 1=perfect association",
    "Total DEGs across all states (top 50 per state)",
    "Pathways with strong enrichment signal"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "Auto_states_comparison_summary.csv", row.names = FALSE)
cat("Saved: Auto_states_comparison_summary.csv\n")

cat("=== Comparison complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
