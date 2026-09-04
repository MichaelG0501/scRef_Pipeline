####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_hybrid_pairwise_percell_heatmap.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_task6_hybrid_pairwise_percell_heatmap.R
# Hybrid-only subdivision and visualisation based on Approach B style
# Matches states_topmp_hybrid.R logic (line 581+)
#
# Input:  ref_outs/EAC_Ref_epi.rds
#         ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#         ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#         ref_outs/Auto_topmp_v2_states_B.rds
#         ref_outs/Auto_topmp_v2_mp_adj.rds
#         /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#         ref_outs/cancer_signatures.txt
# Output: ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_subtypes.rds
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_states_expanded.rds
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_heatmap.pdf
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_heatmap.pdf
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_mean_heatmap.pdf
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_proportion.pdf
#         ref_outs/task6_hybrid_pairwise/Auto_task6_realstates_plushybrid_umap_top12.pdf
#         ref_outs/task6_hybrid_pairwise/Auto_task6_hybrid_umap_top12.pdf
#         updates/new_updates/summaries/Auto_task6_hybrid_pairwise_summary.csv
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "task6"
out_dir <- paste0(task_prefix, "_hybrid_pairwise")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Auto task6 hybrid pairwise subtyping ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load inputs
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
# Note: loading noreg states/mp_adj specifically as per original task context
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")

cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]

if (file.exists("cancer_signatures.txt")) {
  cs_genes_raw <- read.table("cancer_signatures.txt", header = FALSE)$V1
} else {
  warning("cancer_signatures.txt not found; CS score will be set to 0")
  cs_genes_raw <- character(0)
}

####################
# 2) Canonical MP setup and ordering
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

cc_mps <- c("MP1", "MP7", "MP9")

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

pair_levels <- combn(ordered_group_names, 2, simplify = FALSE)
pair_labels <- vapply(pair_levels, function(x) paste(x, collapse = "__"), character(1))

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

####################
# 3) Align matrices and compute helper scores
####################
common_cells <- intersect(intersect(Cells(tmdata_all), rownames(ucell_scores)), names(state_B))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))

tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))

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

sample_var <- tmdata_all$orig.ident
study_var <- tmdata_all$study
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)

mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc[common_cells, , drop = FALSE])

cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))

cs_genes <- intersect(cs_genes_raw, rownames(tmdata_all))
if (length(cs_genes) > 0) {
  cs_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cs_genes, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:min(50, length(cs_genes))]
  cs_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cs_top50, , drop = FALSE]))
} else {
  cs_score <- rep(0, ncol(tmdata_all))
  names(cs_score) <- colnames(tmdata_all)
}

####################
# 4) Hybrid-B subdivision
####################
HYBRID_GAP_B <- 0.3

group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_noncc))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[common_cells, mps_avail]))
  apply(mp_adj_noncc[common_cells, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)
rownames(group_max) <- common_cells

hybrid_cells <- names(state_B)[state_B == "Hybrid"]
cat(sprintf("Hybrid cells found: %d\n", length(hybrid_cells)))
if (length(hybrid_cells) == 0) stop("No Hybrid cells found in state_B")

assign_hybrid_subtype <- function(score_vec, order_vec, gap_thr = 0.3) {
  score_vec <- score_vec[order_vec]
  ord <- names(sort(score_vec, decreasing = TRUE))[1:2]
  ord <- order_vec[order_vec %in% ord]
  paste(ord, collapse = "__")
}

hybrid_subtype <- vapply(
  hybrid_cells,
  function(cl) assign_hybrid_subtype(group_max[cl, ], ordered_group_names, HYBRID_GAP_B),
  character(1)
)

hybrid_levels <- pair_labels
hybrid_subtype <- factor(hybrid_subtype, levels = hybrid_levels)
names(hybrid_subtype) <- hybrid_cells

expanded_state <- as.character(state_B)
names(expanded_state) <- names(state_B)
expanded_state[names(hybrid_subtype)] <- as.character(hybrid_subtype)
expanded_state <- factor(expanded_state)

saveRDS(hybrid_subtype, file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_subtypes.rds")))
saveRDS(expanded_state, file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_states_expanded.rds")))

####################
# 4b) 5 real states + Hybrid view (for requested UMAP display)
####################
state5_plus_hybrid <- as.character(state_B)
names(state5_plus_hybrid) <- names(state_B)
state5_plus_hybrid[state5_plus_hybrid %in% c("Unresolved")] <- NA_character_
state5_plus_hybrid <- factor(
  state5_plus_hybrid,
  levels = c(ordered_group_names, "Hybrid")
)

####################
# 5) Hybrid-only per-cell heatmap
####################
set.seed(42)
MAX_HYBRID_CELLS <- 6000
subtype_cells <- split(names(hybrid_subtype), hybrid_subtype)
subtype_counts <- table(hybrid_subtype)
subtype_fracs <- subtype_counts / sum(subtype_counts)
cells_per_subtype <- pmax(round(subtype_fracs * MAX_HYBRID_CELLS), 30)

cells_to_plot <- unlist(
  mapply(
    function(cells, n) sample(cells, min(length(cells), n)),
    subtype_cells,
    cells_per_subtype[names(subtype_cells)],
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)

sub_scores_orig <- t(mp_adj_all[cells_to_plot, , drop = FALSE])
cc_block_order <- cc_mps[cc_mps %in% rownames(sub_scores_orig)]
non_cc_block_order <- mp_tree_order_names[
  mp_tree_order_names %in% rownames(sub_scores_orig) & !(mp_tree_order_names %in% cc_mps)
]
mp_row_order <- c(cc_block_order, non_cc_block_order)
sub_scores <- sub_scores_orig[mp_row_order, , drop = FALSE]

mp_label_map <- mp_descriptions
missing_mps <- setdiff(rownames(sub_scores), names(mp_label_map))
if (length(missing_mps) > 0) mp_label_map[missing_mps] <- missing_mps
rownames(sub_scores) <- mp_label_map[rownames(sub_scores)]

split_vec <- factor(as.character(hybrid_subtype[cells_to_plot]), levels = hybrid_levels)
split_vec <- droplevels(split_vec)

meta <- tmdata_all@meta.data
study_vals <- as.character(meta[cells_to_plot, "study"])

state_df <- data.frame(
  cell = names(hybrid_subtype),
  subtype = as.character(hybrid_subtype),
  sample = as.character(meta[names(hybrid_subtype), "orig.ident"]),
  study = as.character(meta[names(hybrid_subtype), "study"]),
  stringsAsFactors = FALSE
)

total_samples <- length(unique(state_df$sample))
total_studies <- length(unique(state_df$study))
subtype_div_df <- state_df %>%
  group_by(subtype) %>%
  summarise(
    sample_cov = n_distinct(sample) / max(total_samples, 1),
    study_cov = n_distinct(study) / max(total_studies, 1),
    diversity_score = 0.5 * sample_cov + 0.5 * study_cov,
    .groups = "drop"
  )
subtype_div_map <- setNames(subtype_div_df$diversity_score, subtype_div_df$subtype)
div_vals <- subtype_div_map[as.character(split_vec)]
div_vals[is.na(div_vals)] <- 0
names(div_vals) <- cells_to_plot

hybrid_cols <- setNames(hue_pal()(length(pair_labels)), pair_labels)
local_hybrid_cols <- hybrid_cols[levels(split_vec)]

study_cols <- setNames(
  DiscretePalette(length(unique(meta$study)), palette = "polychrome"),
  unique(meta$study)
)

max_cc_val <- max(cc_score[cells_to_plot], na.rm = TRUE)
if (is.na(max_cc_val) || max_cc_val == 0) max_cc_val <- 1

max_cs_val <- max(cs_score[cells_to_plot], na.rm = TRUE)
if (is.na(max_cs_val) || max_cs_val == 0) max_cs_val <- 1

col_ann <- HeatmapAnnotation(
  HybridSubtype = split_vec,
  CC_score = cc_score[cells_to_plot],
  CS_score = cs_score[cells_to_plot],
  Diversity = div_vals,
  Study = study_vals,
  col = list(
    HybridSubtype = local_hybrid_cols,
    CC_score = colorRamp2(c(0, max_cc_val), c("white", "darkgreen")),
    CS_score = colorRamp2(c(0, max_cs_val), c("white", "darkorange")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Study = study_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE,
  na_col = "white"
)

mp_to_group <- rep("Other", length(mp_row_order))
names(mp_to_group) <- mp_row_order
mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
  mp_to_group[grp_mps] <- grp
}

group_cols <- setNames(
  c(hue_pal()(length(ordered_group_names)), "gold", "grey70"),
  c(ordered_group_names, "Cell_cycle", "Other")
)

mp_group_label <- mp_to_group
names(mp_group_label) <- rownames(sub_scores)

row_ann <- rowAnnotation(
  MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
  col = list(MP_group = group_cols),
  show_annotation_name = FALSE
)

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

ht_sc <- Heatmap(
  sub_scores,
  name = "Adj score",
  col = col_fun_sc,
  top_annotation = col_ann,
  left_annotation = row_ann,
  column_split = split_vec,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  cluster_rows = FALSE,
  row_split = factor(ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
                     levels = c("Cell_cycle_MPs", "Other_MPs")),
  row_gap = unit(2.5, "mm"),
  column_gap = unit(1.5, "mm"),
  show_row_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
  show_column_names = FALSE,
  column_title_rot = 30,
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_heatmap.pdf")), width = 19, height = max(7, length(mp_row_order) * 0.5), useDingbats = FALSE)
draw(ht_sc, merge_legend = TRUE)
dev.off()

####################
# 6) Hybrid-only proportion plot
####################
overall_df <- data.frame(subtype = hybrid_subtype) %>%
  count(subtype) %>%
  mutate(pct = 100 * n / sum(n), study = "Overall")

per_study_df <- data.frame(
  subtype = hybrid_subtype,
  study = meta[names(hybrid_subtype), "study"],
  stringsAsFactors = FALSE
) %>%
  count(subtype, study) %>%
  group_by(study) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

plot_df <- bind_rows(overall_df, per_study_df)
plot_df$subtype <- factor(plot_df$subtype, levels = hybrid_levels)
plot_df$study <- factor(plot_df$study, levels = c("Overall", sort(unique(as.character(meta$study)))))

p_prop <- ggplot(plot_df, aes(x = study, y = pct, fill = subtype)) +
  geom_col(colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = hybrid_cols) +
  labs(
    title = "Hybrid subtype proportions",
    subtitle = "Pairwise-only double-hybrid classes",
    x = NULL,
    y = "% of hybrid cells",
    fill = "Hybrid subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_proportion.pdf")), p_prop, width = 14, height = 7)

####################
# 6b) Pairwise-only hybrid matrix heatmap
####################
pair_mat <- matrix(0, nrow = length(ordered_group_names), ncol = length(ordered_group_names),
                   dimnames = list(ordered_group_names, ordered_group_names))
pair_tab <- table(as.character(hybrid_subtype))
for (lbl in names(pair_tab)) {
  if (!grepl("__", lbl)) next
  sp <- strsplit(lbl, "__")[[1]]
  if (length(sp) != 2) next
  a <- sp[1]
  b <- sp[2]
  if (a %in% ordered_group_names && b %in% ordered_group_names) {
    pair_mat[a, b] <- pair_tab[[lbl]]
    pair_mat[b, a] <- pair_tab[[lbl]]
  }
}
pair_pct <- 100 * pair_mat / max(length(state_B), 1)
pair_df <- as.data.frame(as.table(pair_pct), stringsAsFactors = FALSE)
colnames(pair_df) <- c("StateA", "StateB", "Pct")

p_pair <- ggplot(pair_df, aes(StateB, StateA, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Pct)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick3") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Pairwise-only Hybrid heatmap (% of all cells)", x = NULL, y = NULL, fill = "%")

ggsave(file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_heatmap.pdf")), p_pair, width = 8, height = 6)

####################
# 7) Mean heatmap per hybrid subtype (non-CC MPs)
####################
mean_scores_df <- as.data.frame(mp_adj_noncc[names(hybrid_subtype), , drop = FALSE])
mean_scores_df$subtype <- hybrid_subtype
mean_mat <- mean_scores_df %>%
  group_by(subtype) %>%
  summarise(across(all_of(colnames(mp_adj_noncc)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("subtype") %>%
  as.matrix()

ordered_subtypes_present <- hybrid_levels[hybrid_levels %in% rownames(mean_mat)]
mean_mat <- mean_mat[ordered_subtypes_present, , drop = FALSE]

non_cc_tree_order <- mp_tree_order_names[
  mp_tree_order_names %in% colnames(mp_adj_noncc) & !(mp_tree_order_names %in% cc_mps)
]
mean_mat <- mean_mat[, non_cc_tree_order, drop = FALSE]
colnames(mean_mat) <- mp_descriptions[colnames(mean_mat)]

cell_counts <- table(hybrid_subtype)
row_anno <- rowAnnotation(
  Cells = anno_barplot(as.numeric(cell_counts[rownames(mean_mat)]), gp = gpar(fill = "steelblue"), width = unit(2, "cm"))
)

col_fun_z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mean <- Heatmap(
  mean_mat,
  name = "Mean Z-score",
  col = col_fun_z,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9.5, fontface = "italic"),
  column_names_rot = 45,
  column_title = "Mean non-CC MP activity per hybrid subtype",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mean_mat[i, j]), x, y, gp = gpar(fontsize = 6.5))
  }
)

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_mean_heatmap.pdf")), width = 15, height = 8, useDingbats = FALSE)
draw(ht_mean, merge_legend = TRUE)
dev.off()

####################
# 8) UMAPs in per-sample style (separate real-state and hybrid views)
####################
if ("umap" %in% names(tmdata_all@reductions)) {
  # Real states + Hybrid (single Hybrid class), per-sample facets
  state_real_hybrid <- as.character(state_B)
  names(state_real_hybrid) <- names(state_B)
  state_real_hybrid[state_real_hybrid == "Unresolved"] <- NA_character_
  state_real_hybrid <- factor(state_real_hybrid, levels = c(ordered_group_names, "Hybrid"))
  names(state_real_hybrid) <- names(state_B)
  tmdata_all$Auto_state_real_hybrid <- state_real_hybrid[Cells(tmdata_all)]

  state5_cols <- setNames(
    c(hue_pal()(length(ordered_group_names)), "orchid"),
    c(ordered_group_names, "Hybrid")
  )

  # Rank samples by abundance of defined states (excluding unresolved)
  md_rank <- tmdata_all@meta.data
  rank_df <- md_rank %>%
    dplyr::mutate(state_plot = Auto_state_real_hybrid) %>%
    dplyr::filter(!is.na(state_plot)) %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(total_defined = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_defined))

  # safe slicing
  n_top <- min(12, nrow(rank_df))
  if (n_top == 0) n_top <- 0
  rank_df <- rank_df %>% dplyr::slice_head(n = n_top)

  top_ids <- as.character(rank_df$orig.ident)

  # Build one-page-per-sample list for real-state+hybrid view
  real_plots <- vector("list", length(top_ids))
  for (i in seq_along(top_ids)) {
    sid <- top_ids[i]
    sub_obj <- subset(tmdata_all, subset = orig.ident == sid)
    sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub_obj <- ScaleData(sub_obj, verbose = FALSE)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
    sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

    p <- DimPlot(
      sub_obj,
      reduction = "umap",
      group.by = "Auto_state_real_hybrid",
      cols = state5_cols,
      pt.size = 0.7
    ) +
      labs(
        title = sid,
        subtitle = sprintf("Defined cells: %d", rank_df$total_defined[i])
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    real_plots[[i]] <- p
  }

  pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_realstates_plushybrid_umap_top12.pdf")), width = 12, height = 9, useDingbats = FALSE)
  for (p in real_plots) print(p)
  dev.off()

  # Hybrid-only subtype view, per-sample facets
  hybrid_plot_label <- rep(NA_character_, length(common_cells))
  names(hybrid_plot_label) <- common_cells
  hybrid_plot_label[names(hybrid_subtype)] <- as.character(hybrid_subtype)
  hybrid_plot_label <- factor(hybrid_plot_label, levels = hybrid_levels)
  tmdata_all$Auto_hybrid_subtype <- hybrid_plot_label[Cells(tmdata_all)]

  # Build one-page-per-sample list for hybrid-only subtype view
  hybrid_plots <- vector("list", length(top_ids))
  for (i in seq_along(top_ids)) {
    sid <- top_ids[i]
    sub_obj <- subset(tmdata_all, subset = orig.ident == sid)
    sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub_obj <- ScaleData(sub_obj, verbose = FALSE)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
    sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

    p <- DimPlot(
      sub_obj,
      reduction = "umap",
      group.by = "Auto_hybrid_subtype",
      cols = hybrid_cols,
      pt.size = 0.7
    ) +
      labs(
        title = sid,
        subtitle = sprintf("Hybrid-labelled cells in sample"),
        colour = "Hybrid subtype"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    hybrid_plots[[i]] <- p
  }

  pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_umap_top12.pdf")), width = 12, height = 9, useDingbats = FALSE)
  for (p in hybrid_plots) print(p)
  dev.off()
}

####################
# 9) Summary output
####################
summary_df <- data.frame(
  subtype = names(table(hybrid_subtype)),
  cells = as.integer(table(hybrid_subtype)),
  stringsAsFactors = FALSE
)
summary_df$pct <- 100 * summary_df$cells / sum(summary_df$cells)

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(summary_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_summary.csv")), row.names = FALSE)

cat("Saved task6 pairwise-only per-cell hybrid heatmap and plots in:", out_dir, "\n")
