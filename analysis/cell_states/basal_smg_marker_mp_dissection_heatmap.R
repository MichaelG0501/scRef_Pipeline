####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/basal_smg_marker_mp_dissection_heatmap.R
#   Methodology: not required (fixed marker/MP aggregation and plotting)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Relates fixed basal/SMG marker expression to the dominant MP
#     within the two current centred refined metaplasia states.
#   Inputs: EAC_Ref_epi.rds, centred merged refined UCell scores, and centred states.
#   Outputs: ref_outs/Auto_basal_smg_marker_mp_dissection/ figures and summary tables;
#     updates/new_updates/summaries/basal_smg_marker_mp_dissection_summary.csv.
#   Cache/replot: deterministic plot-only rebuild from live inputs.
#   Run: qsub analysis/cell_states/basal_smg_marker_mp_dissection_heatmap.sh
#   Environment: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
})

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

out_dir <- "Auto_basal_smg_marker_mp_dissection"
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
summary_dir <- file.path(project_dir, "updates", "new_updates", "summaries")
for (path in c(fig_dir, table_dir, summary_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Inputs
epi_file <- "EAC_Ref_epi.rds"
ucell_file <- "Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds"
state_file <- "Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds"

message("Loading Seurat object...")
tmdata_epi <- readRDS(epi_file)

message("Loading UCell scores and state assignments...")
ucell_scores <- readRDS(ucell_file)
states <- readRDS(state_file)

# The target markers
basal_markers <- c("ISG20", "KRT17", "CEACAM6", "ADGRF1", "DUOX2")
smg_markers <- c("ADH1C", "OLFM4", "MT1G", "HMGCS2", "WFDC2", "PPP1R1B")
markers_to_plot <- c(basal_markers, smg_markers)
markers_to_plot <- rev(markers_to_plot) # Reverse so basal is on top in ggplot if needed

# Sub-group MPs
basal_mps <- c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+")
smg_mps <- c("MP8+", "MP8b", "MP16", "MP18b", "MP17")

target_states <- c("Squamous-to-intestinal", "Glandular-to-intestinal")

# Align cells
common_cells <- intersect(colnames(tmdata_epi), names(states))
common_cells <- intersect(common_cells, rownames(ucell_scores))

# Filter to only the states of interest
target_cells <- common_cells[states[common_cells] %in% target_states]

message("Extracting expression data...")
DefaultAssay(tmdata_epi) <- "RNA"
genes_use <- intersect(markers_to_plot, rownames(tmdata_epi))
if(length(genes_use) == 0) stop("None of the requested markers are in the Seurat object!")

expr_mat <- GetAssayData(tmdata_epi, assay = "RNA", slot = "data")[genes_use, target_cells, drop = FALSE]
expr_mat_dense <- as.matrix(expr_mat)

# Subgrouping logic
message("Subgrouping cells by Top MP...")
cell_meta <- data.frame(
  cell = target_cells,
  state = states[target_cells],
  stringsAsFactors = FALSE
)

# Function to get top MP from a specific list
get_top_mp <- function(cell_id, mp_list, scores_mat) {
  # Get scores for this cell for the given MPs
  s <- scores_mat[cell_id, mp_list, drop = FALSE]
  if(ncol(s) == 0) return(NA)
  top_idx <- which.max(s[1, ])
  if(length(top_idx) > 0) return(colnames(s)[top_idx])
  return(NA)
}

top_mps <- character(nrow(cell_meta))
for (i in seq_len(nrow(cell_meta))) {
  c_id <- cell_meta$cell[i]
  st <- cell_meta$state[i]
  
  if (st == "Squamous-to-intestinal") {
    top_mps[i] <- get_top_mp(c_id, basal_mps, ucell_scores)
  } else if (st == "Glandular-to-intestinal") {
    top_mps[i] <- get_top_mp(c_id, smg_mps, ucell_scores)
  } else {
    top_mps[i] <- NA
  }
}

cell_meta$Top_MP <- top_mps
cell_meta <- cell_meta[!is.na(cell_meta$Top_MP), ]

# Order levels
cell_meta$Top_MP <- factor(cell_meta$Top_MP, levels = c(basal_mps, smg_mps))
cell_meta$state <- factor(cell_meta$state, levels = target_states)

# Prepare for DotPlot
message("Preparing DotPlot data...")
plot_df <- cell_meta %>%
  group_by(state, Top_MP) %>%
  summarise(n_cells = n(), .groups = "drop")

# We need mean expression and pct detected per group for each gene
agg_list <- lapply(genes_use, function(g) {
  g_expr <- expr_mat_dense[g, cell_meta$cell]
  df <- data.frame(
    cell = cell_meta$cell,
    state = cell_meta$state,
    Top_MP = cell_meta$Top_MP,
    expr = g_expr
  )
  
  df %>%
    group_by(state, Top_MP) %>%
    summarise(
      mean_expr = mean(expr),
      pct_detected = mean(expr > 0),
      .groups = "drop"
    ) %>%
    mutate(gene = g)
})
agg_df <- bind_rows(agg_list)

# Add full MP descriptions
mp_desc_map <- c(
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor"
)

# Create marker state grouping with newlines to prevent obscuring
agg_df$marker_state <- ifelse(agg_df$gene %in% basal_markers, "Squamous-to-\nintestinal", 
                              ifelse(agg_df$gene %in% smg_markers, "Glandular-to-\nintestinal", "Other"))

# Keep all marker-state pairs so we can see cross-state expression for comparison
# (Filter removed)

# Apply MP mapping
mp_levels <- c(basal_mps, smg_mps)
mp_levels_desc <- paste0(mp_levels, "\n", mp_desc_map[mp_levels])
agg_df$Top_MP_desc <- paste0(agg_df$Top_MP, "\n", mp_desc_map[as.character(agg_df$Top_MP)])
agg_df$Top_MP_desc <- factor(agg_df$Top_MP_desc, levels = mp_levels_desc)

# Ordering
agg_df$gene <- factor(agg_df$gene, levels = rev(genes_use))
agg_df$state <- factor(agg_df$state, levels = target_states)
target_marker_states <- c("Squamous-to-\nintestinal", "Glandular-to-\nintestinal")
agg_df$marker_state <- factor(agg_df$marker_state, levels = target_marker_states)

message("Plotting...")
dotplot <- ggplot(agg_df, aes(x = Top_MP_desc, y = gene)) +
  geom_point(aes(size = pct_detected, fill = mean_expr), shape = 21, color = "grey35", stroke = 0.35) +
  facet_grid(marker_state ~ state, scales = "free", space = "free") +
  scale_size_continuous(range = c(0.5, 6), labels = scales::percent_format(accuracy = 1), name = "% Detected") +
  scale_fill_gradientn(
    colors = c("gray95", "khaki1", "orange", "red", "darkred"),
    values = scales::rescale(c(0, 0.5, 1.2, 2.0, 3.0)),
    limits = c(0, 3.0),
    oob = scales::squish,
    name = "Mean\nExpression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold.italic"),
    strip.text = element_text(face = "bold", size = 11),
    strip.text.y = element_text(angle = 270, size = 9, lineheight = 1.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    plot.margin = margin(t = 5, r = 15, b = 5, l = 5)
  ) +
  labs(x = "Top State-Defining MP", y = "Marker Genes", title = "Marker Expression by Top Metaprogram Subset")

out_pdf <- file.path(fig_dir, "marker_mp_dissection_dotplot.pdf")
ggsave(out_pdf, dotplot, width = 10.5, height = 5.5, useDingbats = FALSE)
message("Saved PDF to: ", out_pdf)

####################
# Persistent source table supports plot-only regeneration and compact review.
####################
write.csv(agg_df, file.path(table_dir, "marker_mp_dissection_plot_data.csv"), row.names = FALSE)
summary_df <- plot_df %>%
  mutate(output_pdf = out_pdf)
write.csv(summary_df, file.path(summary_dir, "basal_smg_marker_mp_dissection_summary.csv"), row.names = FALSE)

message("Done.")
