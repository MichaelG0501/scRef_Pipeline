####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_manual_state_umap.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/states_umap.R
# Reorganized as part of analysis/ restructuring
####################
library(dplyr)
library(ggplot2)
library(scales)

# --- 1. Data Prep ---
# Ensure manual_state is a factor so 'table' creates 0s for missing states
tmdata_all$manual_state <- as.factor(tmdata_all$manual_state)

md <- tmdata_all@meta.data[, c("orig.ident", "manual_state")]
# 'table' ensures we get rows with n=0 if a state is missing in a sample
counts_long <- as.data.frame(table(md$orig.ident, md$manual_state))
colnames(counts_long) <- c("orig.ident", "manual_state", "n")

props_long <- counts_long %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Calculate total N for plotting overlay
sample_totals <- props_long %>%
  group_by(orig.ident) %>%
  summarise(total_sample_n = sum(n), .groups = "drop")

target_states <- c("Classic Proliferative", "Basal to Intest. Meta", "Intestinal Metaplasia", "Stress-adaptive")

rank_df <- props_long %>%
  filter(manual_state %in% target_states) %>%
  group_by(orig.ident) %>%
  summarise(
    # 1. Total cells in just these target states
    target_n = sum(n),
    
    # 2. Geometric Mean of Counts: The "Strictest" Balance Score.
    # Formula: exp(mean(log(n + 1)))
    # If any state has n=0 or n=5, this score drops drastically.
    # It only becomes high if ALL target states are abundant.
    geo_mean_score = exp(mean(log(n + 1))),
    
    # Optional: Keep min_n just to see which state is the bottleneck
    min_n = min(n),
    .groups = "drop"
  ) %>%
  left_join(sample_totals, by = "orig.ident") %>%
  # Filter: Must have basic abundance to be considered
  filter(target_n > 20) %>%
  # Sort by the Geometric Mean Score
  arrange(desc(geo_mean_score)) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 40)

# --- 3. Visualization ---
plot_df <- props_long %>%
  filter(orig.ident %in% rank_df$orig.ident) %>%
  left_join(rank_df %>% select(orig.ident, rank, total_sample_n, geo_mean_score), by = "orig.ident")

scale_factor <- max(rank_df$total_sample_n) / 1.1

ggplot(plot_df) +
  # A. Stacked Bar Chart
  geom_col(aes(x = reorder(orig.ident, rank), y = prop, fill = manual_state), width = 0.75) +
  
  # B. Total Cell Count Overlay
  geom_point(data = rank_df, aes(x = reorder(orig.ident, rank), y = total_sample_n / scale_factor), 
             color = "black", size = 2.5, shape = 18) +
  geom_line(data = rank_df, aes(x = reorder(orig.ident, rank), y = total_sample_n / scale_factor, group = 1), 
            color = "black", alpha = 0.4, linetype = "dashed") +
  
  # C. Styling
  scale_fill_manual(values = manual_cols) +
  scale_y_continuous(
    name = "Proportion of Cells",
    expand = c(0,0),
    sec.axis = sec_axis(~ . * scale_factor, name = "Total Cell Count (N)", labels = comma)
  ) +
  theme_minimal() +
  labs(
    title = "Samples Ranked by Strict Balance (Geometric Mean)",
    subtitle = "High Rank = Substantial counts in ALL target states (No state is left behind)",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

rankplot = 3
#plot <- subset(tmdata_all, subset = orig.ident == rank_df$orig.ident[rankplot])
plot <- readRDS(paste0("by_samples/", rank_df$orig.ident[rankplot], "/", rank_df$orig.ident[rankplot], "_epi_f.rds"))
plot <- subset(plot, coexpression == "singlet" & marker_expression == "good" & (malignancy == "malignant_level_1" | malignancy == "malignant_level_2"))
plot <- subset(plot, malignancy == "malignant_level_1" | malignancy == "malignant_level_2")
plot$manual_state <- tmdata_all$manual_state[match(colnames(plot), colnames(tmdata_all))]


DimPlot(plot, 
        reduction = "umap", 
        group.by = "manual_state", 
        cols = manual_cols, 
        pt.size = 0.8) + 
  labs(title = "Cell States: Alcindor_2025_SRR27335929",
       subtitle = "UMAP of individual sample composition") +
  theme_minimal()


library(Seurat)
library(ggplot2)
library(cowplot) # or library(patchwork)

# 1. Setup the PDF Output
pdf("columnar_states_UMAP.pdf", width = 14, height = 18)

# 2. Loop through Top 12 Ranked Samples
plot_list <- list()

# We loop 1 to 12 to generate all plots first
for (i in 1:12) {
  
  # Get Sample ID
  sample_id <- rank_df$orig.ident[i]
  cat(sprintf("Processing Rank %d: %s ...\n", i, sample_id))
  
  # A. Subset
  # Note: ensure rank_df is character, not factor, to avoid errors
  sub_obj <- subset(tmdata_all, subset = orig.ident == as.character(sample_id))
  
  # B. Run Seurat Pipeline (Fast Version)
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
  # Using fewer dims (1:15) is usually sufficient for single samples and faster
  sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)
  
  # C. Generate Plot
  p <- DimPlot(sub_obj, 
               reduction = "umap", 
               group.by = "manual_state", 
               cols = manual_cols, 
               pt.size = 0.7) + # Slightly larger dots for clarity
    labs(title = paste0(sample_id),
         subtitle = paste0(rank_df$total_sample_n[i])) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8)
    )
  
  # Store in list
  plot_list[[i]] <- p
}

# 3. Print to PDF (6 per page)
# Page 1: Ranks 1-6
print(plot_grid(plotlist = plot_list[1:6], ncol = 2, nrow = 3))

# Page 2: Ranks 7-12
print(plot_grid(plotlist = plot_list[7:12], ncol = 2, nrow = 3))

# 4. Close PDF
dev.off()

cat("\nPDF Generated: Top12_Balanced_Samples_UMAP.pdf\n")
