####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/centred/3ca_vs_refined_mp_correlation.R
#   Description: Correlation of 3CA MPs vs Centred Refined MPs for all cells and unresolved cells.
#   Map: analysis/ANALYSIS_MAP.md
####################

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(Seurat)
  library(ggplot2)
})

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- "Metaprogrammes_Results/centred/correlation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs
ucell_refined_file <- "Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds"
ucell_3ca_file <- "UCell_3CA_MPs.rds"
states_file <- "Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds"

ucell_refined <- readRDS(ucell_refined_file)
ucell_3ca <- readRDS(ucell_3ca_file)
states <- readRDS(states_file)

common_cells <- intersect(rownames(ucell_refined), rownames(ucell_3ca))
common_cells <- intersect(common_cells, names(states))

ucell_refined <- ucell_refined[common_cells, ]
ucell_3ca <- ucell_3ca[common_cells, ]
states <- states[common_cells]

unresolved_cells <- names(states)[states == "Unresolved"]

# Refined MP groups
cc_mps <- c("MP1", "MP5", "MP13+")
state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP2x"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)
excluded_mps <- c("MP11c", "MP18a")

plot_mp_order <- c(cc_mps, unlist(state_groups, use.names = FALSE), excluded_mps)
plot_mp_order <- intersect(plot_mp_order, colnames(ucell_refined))

# Clean 3CA names
clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}
colnames(ucell_3ca) <- clean_3ca_name(colnames(ucell_3ca))

CC_FIXED <- clean_3ca_name(c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
))

# Cor for all cells
cor_all <- cor(ucell_refined[, plot_mp_order], ucell_3ca, method = "pearson")

# Cor for unresolved
cor_unres <- cor(ucell_refined[unresolved_cells, plot_mp_order], ucell_3ca[unresolved_cells, ], method = "pearson")

# Row annotation (State)
state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)
mp_group_cols <- c(
  "Cell cycle" = "#6B7280",
  state_cols,
  "Excluded" = "grey80",
  "Other" = "grey70"
)

mp_to_group <- rep("Other", length(plot_mp_order))
names(mp_to_group) <- plot_mp_order
mp_to_group[intersect(cc_mps, names(mp_to_group))] <- "Cell cycle"
for (grp in names(state_groups)) {
  mp_to_group[intersect(state_groups[[grp]], names(mp_to_group))] <- grp
}
mp_to_group[intersect(excluded_mps, names(mp_to_group))] <- "Excluded"

row_ann <- rowAnnotation(
  Group = factor(mp_to_group, levels = names(mp_group_cols)),
  col = list(Group = mp_group_cols)
)

# MP descriptions for row names
mp_desc_map <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "Fanconi/HR repair-active glandular progenitor",
  "MP2+" = "MYC driven biosynthesis",
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
  "MP17" = "Immune-interactive glandular progenitor",
  "MP2x" = "Wnt-active glandular stem/progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry",
  "MP11c" = "Excluded",
  "MP18a" = "Excluded"
)

get_full_name <- function(mp_names) {
  desc <- mp_desc_map[mp_names]
  desc[is.na(desc)] <- "Unknown"
  paste0(mp_names, ": ", desc)
}
rownames(cor_all) <- get_full_name(rownames(cor_all))
rownames(cor_unres) <- get_full_name(rownames(cor_unres))

# Column grouping for 3CA MPs
col_group <- ifelse(colnames(cor_all) %in% CC_FIXED, "Cell cycle 3CA", "Other 3CA")
col_sums_all <- colSums(cor_all, na.rm = TRUE)
col_order <- order(factor(col_group, levels = c("Cell cycle 3CA", "Other 3CA")), -col_sums_all)

cor_all <- cor_all[, col_order]
cor_unres <- cor_unres[, col_order]
col_group_sorted <- factor(col_group[col_order], levels = c("Cell cycle 3CA", "Other 3CA"))

col_fun <- colorRamp2(c(-0.4, 0, 0.4), c("navy", "white", "firebrick3"))

ht_all <- Heatmap(
  cor_all,
  name = "Pearson R",
  col = col_fun,
  left_annotation = row_ann,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_group_sorted,
  row_split = factor(mp_to_group, levels = names(mp_group_cols)),
  row_title = NULL,
  column_title = paste0("All Cells (n = ", length(common_cells), ")"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 10),
  row_gap = unit(2, "mm"),
  column_gap = unit(4, "mm")
)

ht_unres <- Heatmap(
  cor_unres,
  name = "Pearson R",
  col = col_fun,
  left_annotation = row_ann,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_group_sorted,
  row_split = factor(mp_to_group, levels = names(mp_group_cols)),
  row_title = NULL,
  column_title = paste0("Unresolved Cells (n = ", length(unresolved_cells), ")"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 10),
  row_gap = unit(2, "mm"),
  column_gap = unit(4, "mm")
)

pdf(file.path(out_dir, "3ca_vs_refined_mp_correlation.pdf"), width = 16, height = 12)
draw(ht_all, column_title = "Correlation of Centred Refined MPs vs 3CA MPs (All Cells)", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
draw(ht_unres, column_title = "Correlation of Centred Refined MPs vs 3CA MPs (Unresolved Cells)", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

cat("Saved correlation heatmap to", file.path(out_dir, "3ca_vs_refined_mp_correlation.pdf"), "\n")

# --- Per-cell heatmap and barplot for unresolved cells ---

tmdata_all <- readRDS("EAC_Ref_epi.rds")
sample_var <- tmdata_all$orig.ident[common_cells]
study_var <- tmdata_all$study[common_cells]

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  
  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  
  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

top_10_least_3ca <- tail(colnames(cor_all)[col_group_sorted == "Other 3CA"], 10)
selected_refined <- c("MP13+", "MP11c", "MP18a")

# For heatmap 1 (All non-CC 3CA MPs)
all_non_cc_3ca <- setdiff(colnames(ucell_3ca), CC_FIXED)
combined_scores_all <- cbind(ucell_refined[, selected_refined, drop = FALSE], ucell_3ca[, all_non_cc_3ca, drop = FALSE])
mp_adj_all <- z_normalise(combined_scores_all, sample_var, study_var)
sub_scores <- t(mp_adj_all[unresolved_cells, ])

# For assignment (Top 10 + 3)
combined_scores_top10 <- cbind(ucell_refined[, selected_refined, drop = FALSE], ucell_3ca[, top_10_least_3ca, drop = FALSE])
mp_adj_top10 <- z_normalise(combined_scores_top10, sample_var, study_var)

# TopMP assignment calculation function
get_assignments <- function(mat, include_hybrid=TRUE) {
  group_max <- as.matrix(mat)
  best_group_idx <- max.col(group_max, ties.method = "first")
  best_group_val <- apply(group_max, 1, max, na.rm = TRUE)
  state_vec <- colnames(group_max)[best_group_idx]
  state_vec[!is.finite(best_group_val) | best_group_val < 0.5] <- "Unresolved"
  
  if (include_hybrid) {
    sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
    gap <- sorted_groups[, 1] - sorted_groups[, 2]
    state_vec[(gap < 0.3) & (state_vec != "Unresolved")] <- "Hybrid"
  }
  return(state_vec)
}

state_vec_hybrid <- get_assignments(mp_adj_top10[unresolved_cells, ], include_hybrid=TRUE)
state_vec_nohybrid <- get_assignments(mp_adj_top10[unresolved_cells, ], include_hybrid=FALSE)

split_levels <- c(colnames(combined_scores_top10), "Hybrid", "Unresolved")
split_vec <- factor(state_vec_hybrid, levels = split_levels)

# Color map
topmp_cols <- c(
  setNames(scales::hue_pal()(ncol(combined_scores_top10)), colnames(combined_scores_top10)),
  "Hybrid" = "black",
  "Unresolved" = "grey80"
)

# Column clustering for the heatmap
col_order_list <- lapply(levels(split_vec), function(lvl) {
  idx <- which(as.character(split_vec) == lvl)
  if (length(idx) <= 1) return(idx)
  hc <- hclust(dist(t(sub_scores[, idx, drop=FALSE])), method="ward.D2")
  idx[hc$order]
})
col_order <- unlist(col_order_list, use.names = FALSE)

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm=TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

top_ann <- HeatmapAnnotation(
  Assignment = split_vec,
  col = list(Assignment = topmp_cols)
)

ht_cells <- Heatmap(
  sub_scores,
  name = "Adj score",
  col = col_fun_sc,
  top_annotation = top_ann,
  column_split = split_vec,
  column_order = col_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_title = "Unresolved Cells (All non-CC 3CA MPs)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  column_title_rot = 0,
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE
)

# Page 2 heatmap (Top 10 least sum 3CA MPs)
page2_mps <- c(selected_refined, top_10_least_3ca)
sub_scores_page2 <- sub_scores[page2_mps, , drop = FALSE]

lim2 <- as.numeric(quantile(abs(sub_scores_page2), 0.98, na.rm=TRUE))
col_fun_sc2 <- colorRamp2(c(-lim2, 0, lim2), c("navy", "white", "firebrick3"))

ht_cells_page2 <- Heatmap(
  sub_scores_page2,
  name = "Adj score",
  col = col_fun_sc2,
  top_annotation = top_ann,
  column_split = split_vec,
  column_order = col_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_title = "Unresolved Cells (Top 10 least sum 3CA MPs)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  column_title_rot = 0,
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE
)

pdf(file.path(out_dir, "unresolved_cells_heatmap.pdf"), width = 20, height = 12, onefile = TRUE)
draw(ht_cells)
draw(ht_cells_page2)
dev.off()

# Barplots (using all non-CC 3CA MPs for assignment)
state_vec_hybrid_all <- get_assignments(mp_adj_all[unresolved_cells, ], include_hybrid=TRUE)
state_vec_nohybrid_all <- get_assignments(mp_adj_all[unresolved_cells, ], include_hybrid=FALSE)

topmp_cols_all <- c(
  setNames(scales::hue_pal()(ncol(combined_scores_all)), colnames(combined_scores_all)),
  "Hybrid" = "black",
  "Unresolved" = "grey80"
)

bar_df_hybrid <- data.frame(Assignment = state_vec_hybrid_all, stringsAsFactors = FALSE) %>%
  count(Assignment) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(desc(pct))
bar_df_hybrid$Assignment <- factor(bar_df_hybrid$Assignment, levels = bar_df_hybrid$Assignment)

p_bar_hybrid <- ggplot(bar_df_hybrid, aes(x = Assignment, y = pct, fill = Assignment)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = topmp_cols_all, drop = FALSE) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(title = "TopMP Assignment (With Hybrid)", y = "Proportion (%)", x = NULL)

bar_df_nohybrid <- data.frame(Assignment = state_vec_nohybrid_all, stringsAsFactors = FALSE) %>%
  count(Assignment) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(desc(pct))
bar_df_nohybrid$Assignment <- factor(bar_df_nohybrid$Assignment, levels = bar_df_nohybrid$Assignment)

p_bar_nohybrid <- ggplot(bar_df_nohybrid, aes(x = Assignment, y = pct, fill = Assignment)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = topmp_cols_all, drop = FALSE) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(title = "TopMP Assignment (No Hybrid)", y = "Proportion (%)", x = NULL)

pdf(file.path(out_dir, "unresolved_cells_barplot.pdf"), width = 12, height = 8, onefile = TRUE)
print(p_bar_hybrid)
print(p_bar_nohybrid)
dev.off()

cat("Saved per-cell heatmap and barplot to", out_dir, "\n")
