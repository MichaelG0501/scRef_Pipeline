####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/top_diverse_sample_state_umap.R
#   Methodology: not required (sample ranking and per-sample UMAP plotting)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_top_sample_umap_all_states.R
# Plot a freshly fitted UMAP for the most state-diverse sample, including
# Unresolved and Hybrid cells. Diversity is the geometric mean of the five
# current centred-state counts plus one.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#
# Output:
#   ref_outs/top_diverse_sample_state_umap/figures/Auto_top1_sample_umap_all_states.pdf
#   ref_outs/top_diverse_sample_state_umap/tables/Auto_top1_sample_summary.csv
#   updates/new_updates/summaries/top_diverse_sample_state_umap_summary.csv
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

out_dir <- "top_diverse_sample_state_umap"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
summary_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/updates/new_updates/summaries"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
message("Loading data ...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")

# Align cells
common_cells <- intersect(names(state_B), Cells(tmdata_all))
tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
tmdata_all$state_B <- state_B[common_cells]

# Diversity metric constants
target_states <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "Stress adaptive",
  "SMG to intestinal metaplasia",
  "Cancer-cell immune mimicry"
)

# Colors provided by USER
group_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "Stress adaptive" = "#984EA3",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

# Identify Top 1 Sample by Diversity
message("Calculating diversity ranking ...")
state_df <- data.frame(
  cell = names(state_B),
  state = as.character(state_B),
  orig.ident = as.character(tmdata_all$orig.ident),
  stringsAsFactors = FALSE
)

sample_totals <- state_df %>% count(orig.ident, name = "total_n")

counts_long <- state_df %>%
  filter(state %in% target_states) %>%
  count(orig.ident, state, .drop = FALSE) %>%
  complete(orig.ident, state = target_states, fill = list(n = 0))

rank_df <- counts_long %>%
  group_by(orig.ident) %>%
  summarise(
    target_n = sum(n),
    geo_mean_score = exp(mean(log(n + 1))),
    .groups = "drop"
  ) %>%
  left_join(sample_totals, by = "orig.ident") %>%
  filter(target_n > 20) %>%
  arrange(desc(geo_mean_score)) %>%
  mutate(rank = row_number())

if (nrow(rank_df) == 0) stop("No sample has more than 20 current primary-state cells.")
top1_sample <- rank_df$orig.ident[1]
message(sprintf("Top 1 sample identified: %s (Diversity Rank 1, n_defined = %d, total_n = %d)", 
                top1_sample, rank_df$target_n[1], rank_df$total_n[1]))

# Subset to top 1 sample
sub_obj <- tmdata_all[, tmdata_all$orig.ident == top1_sample]
sub_obj$state_label <- factor(state_B[Cells(sub_obj)], levels = names(group_cols))

# Quality-of-life: add counts to legend labels
state_counts_sample <- table(sub_obj$state_label)
legend_labels <- setNames(
  paste0(names(group_cols), " (", as.integer(state_counts_sample[names(group_cols)]), ")"),
  names(group_cols)
)

# Per-sample UMAP processing
message(sprintf("Processing UMAP for %s ...", top1_sample))
sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sub_obj <- ScaleData(sub_obj, verbose = FALSE)
n_pcs <- min(30, ncol(sub_obj) - 1)
sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = n_pcs, verbose = FALSE)
dims_use <- min(15, n_pcs)
sub_obj <- RunUMAP(sub_obj, dims = 1:dims_use, verbose = FALSE)

# Plotting
library(ggplot2)

message("Generating plot ...")

# 1. Extract UMAP coordinates and the metadata label from the Seurat object
plot_data <- data.frame(
  UMAP_1 = sub_obj@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = sub_obj@reductions$umap@cell.embeddings[, 2],
  state_label = sub_obj$state_label
)

p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = state_label)) +
  geom_point(shape = 21, color = "black", stroke = 0.2, size = 1.2) + 
  scale_fill_manual(values = group_cols, labels = legend_labels, name = "State") +
  labs(title = paste0("Cell States UMAP - ", top1_sample),
       subtitle = paste0("Total cells: ", ncol(sub_obj)),
       x = "UMAP 1", 
       y = "UMAP 2") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +  # 👈 bigger legend dots
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 11)
  )

# Save output
pdf(file.path(out_dir, "figures", "Auto_top1_sample_umap_all_states.pdf"), width = 10, height = 8)
print(p)
dev.off()

####################
# Persist the ranking and a compact login-node-readable result summary.
####################
write.csv(rank_df, file.path(out_dir, "tables", "Auto_top1_sample_summary.csv"), row.names = FALSE)
summary_row <- rank_df[1, c("orig.ident", "target_n", "geo_mean_score", "total_n", "rank"), drop = FALSE]
summary_row$output_pdf <- file.path(out_dir, "figures", "Auto_top1_sample_umap_all_states.pdf")
write.csv(summary_row, file.path(summary_dir, "top_diverse_sample_state_umap_summary.csv"), row.names = FALSE)
message("Finished. Outputs saved under ref_outs/", out_dir)
