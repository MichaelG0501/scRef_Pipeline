library(patchwork)

plot_df <- cbind(
  as.data.frame(Embeddings(tmdata, "umap")), # UMAP_1 and UMAP_2
  tmdata@meta.data                          # All metadata, including scores and clusters
)
plot_df$seurat_clusters <- as.factor(plot_df$seurat_clusters)

plot_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), # Centered, bold
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11)
  )

mod_plots <- imap(mod_cols, function(mod_col_name, plot_title) {
  
  ggplot(plot_df, aes(x = umap_1, y = umap_2, color = .data[[mod_col_name]])) +
    geom_point(size = 0.5) +
    scale_color_gradientn(colors = c("grey90", "blue"), 
                          limits = c(-0.2, 1.0), 
                          oob = scales::squish) +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Score") +
    ggtitle(plot_title) + # Use the nice name as the title
    plot_theme # Apply the shared theme
})

umap_cluster <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = seurat_clusters)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
  ggtitle("UMAP by cluster") +
  plot_theme +
  # Make legend dots bigger so they are visible
  guides(color = guide_legend(override.aes = list(size = 3))) 

umap_cluster2 <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = celltype_initial)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
  ggtitle("UMAP by cluster") +
  plot_theme +
  # Make legend dots bigger so they are visible
  guides(color = guide_legend(override.aes = list(size = 3)))

umap_cluster3 <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = marker_expression)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Expression") +
  ggtitle("UMAP by cluster") +
  plot_theme +
  # Make legend dots bigger so they are visible
  guides(color = guide_legend(override.aes = list(size = 3))) 

combined_plot <- wrap_plots(c(mod_plots), ncol = 4)
combined_plot <- wrap_plots(c(mod_plots, list(umap_cluster), list(umap_cluster2), list(umap_cluster3)), ncol = 4)

combined_plot





enrich_one <- function(cl) {
  de_genes <- de %>% filter(cluster == cl) %>% pull(gene) %>% unique()
  
  imap_dfr(setsN, ~ {
    a <- length(intersect(.x, de_genes))
    b <- length(setdiff(de_genes, .x))
    c <- length(setdiff(.x, de_genes))
    d <- length(universe) - a - b - c
    p <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")$p.value
    tibble(cell_type = .y, overlap = a, pval = p)
  }) %>%
    mutate(
      padj = p.adjust(pval, "BH"),
      cluster = cl
    ) %>%
    arrange(padj, desc(overlap))
}
enrichment_scores <- unique(de$cluster) %>% map_dfr(enrich_one)






expr <- GetAssayData(tmdata, slot = "data")
rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
cellsAUC <- AUCell_calcAUC(setsN, rankings, aucMaxRank = floor(0.05 * nrow(rankings)))

auc_df <- as.data.frame(t(getAUC(cellsAUC)))
colnames(auc_df) <- paste0("auc_", ct_names)

# add AUCs to metadata
tmdata@meta.data[, colnames(auc_df)] <- auc_df
mod_cols <- colnames(auc_df)  # e.g., "auc_b.cell", …

# --- pick module-score label per cluster (median AUC) ---
cluster_mod <- tmdata@meta.data %>%
  mutate(cluster = tmdata$leiden_clusters) %>%
  group_by(cluster) %>%
  summarize(across(all_of(mod_cols), median, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
  mutate(cell_type = gsub("^auc_", "", mod)) %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(cluster, step2 = cell_type)

tmdata$celltype_step2 <- cluster_mod$step2[match(tmdata$leiden_clusters, cluster_mod$cluster)]

# --- visualization: UMAP of each AUCell signature in a grid + clusters UMAP ---
mod_plots <- imap(setNames(mod_cols, gsub("^auc_", "", mod_cols)), function(mod_name, celltype) {
  FeaturePlot(tmdata, features = mod_name) + ggtitle(celltype)
})

umap_cluster <- DimPlot(tmdata, group.by = "leiden_clusters") + ggtitle("UMAP by cluster")

combined_plot <- wrap_plots(c(mod_plots, list(umap_cluster)), ncol = 4)
combined_plot

