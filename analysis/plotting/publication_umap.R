####################
# Moved from: analysis/beaut_umap.R
# Reorganized as part of analysis/ restructuring
####################
cosmx_list <- list()
for (sample in sample_dirs) {
  rds_path <- file.path("by_samples/", sample, paste0(sample, "_anno.rds"))
  cosmx_list[[sample]] <- readRDS(rds_path)
}


cosmx$celltype_final2 <- combined_metadata[colnames(cosmx), "celltype_final"]



Idents(cosmx) <- cosmx$seurat_clusters
clusters <- levels(Idents(cosmx))
available_genes <- rownames(cosmx)

markers_in_data <- lapply(markers_list, function(gene_set) intersect(gene_set, available_genes))
markers_in_data <- Filter(function(v) length(v) > 0, markers_in_data)
ct_names <- names(markers_in_data)

score_mat <- matrix(NA_real_,
                    nrow = ncol(cosmx),
                    ncol = length(ct_names),
                    dimnames = list(colnames(cosmx), ct_names))

for (cl in clusters) {
  cells_cl <- WhichCells(cosmx, idents = cl)
  if (length(cells_cl) == 0) next
  tm_sub <- subset(cosmx, cells = cells_cl)
  mtx <- GetAssayData(tm_sub, slot = "data")
  
  # for each cell type, pick top 8 expressed markers *within this cluster*
  cl_features <- lapply(markers_in_data, function(genes) {
    g <- intersect(genes, rownames(mtx))
    if (length(g) == 0) {
      character(0)
    } else {
      # mean expression across cells in this cluster
      m <- Matrix::rowMeans(mtx[g, , drop = FALSE])
      g[order(m, decreasing = TRUE)][seq_len(min(4, length(g)))]
    }
  })
  
  # keep only cell types with at least one expressed gene in this cluster
  keep <- vapply(cl_features, function(v) length(v) > 0, logical(1))
  if (!any(keep)) next
  cl_features_kept <- cl_features[keep]
  kept_ct <- names(cl_features_kept)
  
  tm_sub <- AddModuleScore(
    object   = tm_sub,
    features = cl_features_kept,
    name     = "mod_tmp_",
    assay    = DefaultAssay(tm_sub)
  )
  
  score_cols <- paste0("mod_tmp_", seq_along(cl_features_kept))
  scdf <- tm_sub@meta.data[, score_cols, drop = FALSE]
  colnames(scdf) <- kept_ct
  score_mat[cells_cl, kept_ct] <- as.matrix(scdf[rownames(scdf), kept_ct, drop = FALSE])
}

for (ct in ct_names) {
  cosmx@meta.data[[paste0("mod_", ct)]] <- score_mat[colnames(cosmx), ct]
}

mod_cols <- setNames(paste0("mod_", ct_names), ct_names)

cosmx@meta.data[mod_cols] <- lapply(cosmx@meta.data[mod_cols], function(x) {
  ifelse(x >= 1, 1,
         ifelse(x <= 0, 0, x))  # NA marks intermediate scores; replace with 0 if you want strict binary
})

# 1) per-cluster median module scores and z
scores_long <- cosmx@meta.data %>%
  mutate(cluster = cosmx$seurat_clusters) %>%
  group_by(cluster) %>%
  summarize(across(all_of(mod_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
  mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
  group_by(cluster) %>%
  mutate(z = as.numeric(scale(score))) %>%   # z across celltypes for THIS cluster
  ungroup() %>%
  select(cluster, cell_type, score, z)

# 2) step-2 call per cluster (new logic + pair exception)
# final rule: step2 = { ... }
step2_calls <- scores_long %>%
  group_by(cluster) %>%
  summarize(
    cell_types = list(cell_type),
    zs = list(z),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    all_z = list(unlist(zs)),
    all_ct = list(unlist(cell_types)),
    max_idx = which.max(all_z),
    top_ct  = all_ct[max_idx],
    top_z   = all_z[max_idx],
    # NEW: margin rule (top must be >= 0.8 higher than all others to call top_ct)
    step2 = {
      az  <- all_z
      act <- all_ct
      other_idx <- setdiff(seq_along(az), max_idx)
      
      if (length(other_idx) == 0) {
        # only one cell type scored in this cluster
        top_ct
      } else {
        margins <- top_z - az[other_idx]
        # "close" others are those within 0.8 of the top (i.e., margin < 0.8)
        close_idx <- other_idx[margins < gap_cut]
        
        if (length(close_idx) == 0) {
          # top is >= 0.8 higher than ALL others -> call top_ct
          top_ct
        } else if (length(close_idx) == 1) {
          # exactly one close competitor -> allow pair if in allowed_pairs, else unresolved
          other_ct <- act[close_idx]
          pair_vec <- c(top_ct, other_ct)
          is_allowed <- any(vapply(allowed_pairs, function(p) identical(p, pair_vec), logical(1)))
          if (is_allowed) paste0(top_ct, "|", other_ct) else "unresolved"
        } else {
          # multiple close competitors -> unresolved
          "unresolved"
        }
      }
    }
  ) %>%
  ungroup() %>%
  select(cluster, step2)

# 3) map back to cells
cl_map <- step2_calls$step2
names(cl_map) <- step2_calls$cluster

cosmx$celltype_final <- as.character(cl_map[as.character(cosmx$seurat_clusters)])




library(ggplot2)
library(dplyr)
library(scattermore) # Much faster than geom_point for millions of points

# Extract UMAP coordinates and metadata
plot_data <- data.frame(
  UMAP_1 = cosmx@reductions$umap@cell.embeddings[,1],
  UMAP_2 = cosmx@reductions$umap@cell.embeddings[,2],
  celltype = cosmx$celltype_final, 
  core = cosmx$coreid, 
  patient = as.character(cosmx$ID)
)
set.seed(0) # Set seed for reproducibility
plot_data <- plot_data[sample(nrow(plot_data)), ]

# Filter out unresolved cells
plot_data <- plot_data %>%
  filter(celltype != "unresolved") # Change "unresolved" to the exact string used in your data

cell_type_names <- c(
  "b.cell", "dendritic", "endothelial", "epithelial", 
  "fibroblast", "macrophage", "mast", "plasma", "t.cell"
)

my_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#A65628", "#F781BF", "#999999", "#00CED1"
)

# Create the named mapping
names(my_colors) <- cell_type_names
p1 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1, alpha = 0.8) + 
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    title = "CosMx 6k EAC",
    subtitle = paste0("Total annotated cells: ", format(nrow(plot_data), big.mark=",")),
    color = "Celltype", 
    x = "UMAP 1", 
    y = "UMAP 2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),
    
    # ADD THESE THREE LINES
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    aspect.ratio = 1
  )

library(pals)
unique_patients <- unique(plot_data$patient)
my_colors <- as.vector(pals::alphabet(n = 26))
my_colors <- c(my_colors, "#000000", "#555555", "#AAAAAA") 
names(my_colors) <- unique_patients
p2<-ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = patient)) +
  # scattermore is optimized for large datasets
  geom_point(size = 0.03, alpha = 0.8) + 
  scale_color_manual(values = my_colors) +
  theme_void() + # Cleanest look for UMAPs
  labs(
    subtitle = "By patient",
    color = "Patient ID"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Make legend dots readable
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )
p1/p2



library(ggplot2)
library(dplyr)
library(scattermore) # Much faster than geom_point for millions of points

# Extract UMAP coordinates and metadata
plot_data <- data.frame(
  UMAP_1 = tmdata@reductions$umap@cell.embeddings[,1],
  UMAP_2 = tmdata@reductions$umap@cell.embeddings[,2],
  celltype = tmdata$celltype_update
)

# Filter out unresolved cells
plot_data <- plot_data %>%
  filter(celltype != "unresolved_inconsistent") # Change "unresolved" to the exact string used in your data

# Updated cell type list (10 types now)
cell_type_names <- c(
  "b.cell", "dendritic", "endothelial", "epithelial", 
  "fibroblast", "macrophage", "mast", "nk.cell", "plasma", "t.cell"
)

# Added #FFD700 (Gold) for nk.cell - distinct from the existing oranges and yellows
my_colors <- c(
  "#E41A1C", # b.cell
  "#377EB8", # dendritic
  "#4DAF4A", # endothelial
  "#984EA3", # epithelial
  "#FF7F00", # fibroblast
  "#A65628", # macrophage
  "#F781BF", # mast
  "#FFD700", # nk.cell (NEW)
  "#999999", # plasma
  "#00CED1"  # t.cell
)

# Apply names to ensure consistency
names(my_colors) <- cell_type_names

ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  # scattermore is optimized for large datasets
  geom_scattermore(pointsize = 0.8, alpha = 0.8) + 
  scale_color_manual(values = my_colors) +
  theme_void() + # Cleanest look for UMAPs
  labs(
    title = "Single cell atals EAC",
    subtitle = paste0("Total annotated cells: ", format(nrow(plot_data), big.mark=",")),
    color = "Cell Type"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Make legend dots readable
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

#######################################

library(dplyr)
library(ggplot2)
library(scales)

library(dplyr)
library(ggplot2)
library(scales)

# 1) counts df
celltype_df <- as.data.frame(table(EAC_Ref_merged$celltype_update), stringsAsFactors = FALSE) %>%
  rename(celltype_update = Var1, n = Freq) %>%
  filter(celltype_update != "unresolved_inconsistent") %>%
  mutate(
    percent = 100 * n / sum(n),
    label = paste0("N = ", comma(n), "\n(", sprintf("%.1f", percent), "%)")
  )

# 2) colors (named -> always match regardless of order)
cell_type_names <- c(
  "b.cell", "dendritic", "endothelial", "epithelial",
  "fibroblast", "macrophage", "mast", "nk.cell", "plasma", "t.cell"
)
my_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#FFD700", "#999999", "#00CED1"
)
my_colors_named <- setNames(my_colors, cell_type_names)

# 3) order so MOST abundant is on TOP (after coord_flip)
celltype_df <- celltype_df %>%
  filter(celltype_update %in% names(my_colors_named)) %>%
  arrange(n) %>%  # smallest bottom, largest top (after flip)
  mutate(celltype_update = factor(celltype_update, levels = celltype_update))

# 4) place labels just ABOVE each bar (outside, not obscuring)
pad <- max(celltype_df$n) * 0.02  # spacing above bars

ggplot(celltype_df, aes(x = celltype_update, y = -n, fill = celltype_update)) +
  geom_col(width = 0.8, color = "black") +
  geom_text(
    aes(label = label),
    hjust = 1,                 # right-align text
    nudge_y = -pad * 1.8,      # shift LEFT so it clears the bar
    lineheight = 0.95,
    size = 4
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = my_colors_named, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) comma(abs(x)),
    expand = expansion(mult = c(0.20, 0.02))
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )






library(dplyr)
library(ggplot2)
library(cowplot)

# 1) Name the colors so mapping is stable regardless of order
my_colors_named <- setNames(my_colors, cell_type_names)

# 2) Choose legend order (example: by counts descending from celltype_df)
legend_order <- celltype_df %>%
  arrange(desc(n)) %>%
  pull(celltype_update)

# 3) Make a tiny "dummy" dataset just to generate the legend
legend_df <- celltype_df %>%
  mutate(celltype_update = factor(celltype_update, levels = legend_order))

p_leg <- ggplot(legend_df, aes(x = celltype_update, y = n, color = celltype_update)) +
  geom_point(size = 6) +
  scale_color_manual(values = my_colors_named, drop = FALSE) +
  guides(color = guide_legend(title = "Cell Type")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 18, face = "bold"),
    
    # 🔽 Increase distance between legend entries
    legend.key.height = unit(0.9, "cm"),
    legend.spacing.y = unit(0.6, "cm")
  )

leg <- cowplot::get_legend(p_leg)
cowplot::ggdraw(leg)



#########################



library(ggplot2)
library(dplyr)
library(scattermore) # Much faster than geom_point for millions of points

# Extract UMAP coordinates and metadata
plot_data <- data.frame(
  UMAP_1 = tmdata@reductions$umap@cell.embeddings[,1],
  UMAP_2 = tmdata@reductions$umap@cell.embeddings[,2],
  celltype = tmdata$SUR
)

my_colors <- colorRampPalette(brewer.pal(12, "Paired"))(16)

ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.03, alpha = 0.8) + 
  # Apply your 16+ colors
  scale_color_manual(values = my_colors) +
  # theme_classic provides the axis lines you need
  theme_classic() + 
  labs(
    title = "OAC PDOs scRNA-seq",
    subtitle = paste0("Total annotated cells: ", format(nrow(plot_data), big.mark=",")),
    color = "Patient",
    x = "UMAP 1", 
    y = "UMAP 2"
  ) +
  # Maintain your legend dots size
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(
    # Keeping your specific title and subtitle formatting
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    
    # Keeping your legend formatting
    legend.position = "right",
    legend.text = element_text(size = 10),
    
    # Hide axis numbers/ticks for a cleaner "UMAP" look while keeping the lines
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    # Maintain square proportions
    aspect.ratio = 1 
  )

