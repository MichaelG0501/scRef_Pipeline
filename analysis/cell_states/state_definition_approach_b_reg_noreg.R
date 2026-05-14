####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/state_definition_approach_b_reg_noreg.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# state_definition_approach_b_reg_noreg.R
# Unified Approach B workflow with parameterised mode:
#   - mode = reg: regress out CC MPs from non-CC MPs
#   - mode = noreg: do not regress out CC MPs
#
# This script runs BOTH modes in one pass and writes paired comparison figures.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/meta_full_epi.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#
# Output (RDS):
#   ref_outs/Auto_topmp_v2_reg_states_B.rds
#   ref_outs/Auto_topmp_v2_reg_mp_adj.rds
#   ref_outs/Auto_topmp_v2_reg_group_max.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#   ref_outs/Auto_topmp_v2_noreg_group_max.rds
#
# Output (figures):
#   ref_outs/Auto_topmp_v2_reg_noreg_heatmap_B_cconly.pdf
#   ref_outs/Auto_topmp_v2_reg_noreg_proportion_B_withpie.pdf
#   ref_outs/Auto_topmp_v2_reg_noreg_ccscore_boxplot_B.pdf
#
# Output (summary):
#   updates/new_updates/summaries/Auto_topmp_v2_reg_noreg_summary.csv
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

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")

tmdata_all <- readRDS("EAC_Ref_epi.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
meta_full_epi <- readRDS("meta_full_epi.rds")
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)

mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

cc_mps <- c("MP1", "MP7", "MP9")
non_cc_mps <- setdiff(retained_mps, cc_mps)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "Immune Infiltrating"   = c("MP15")
)

# Reorder MPs: CC MPs first (original relative order), then state MPs (state_groups order)
state_ordered_mps <- unlist(state_groups, use.names = FALSE)
# Identify original order from tree if available
orig_tree_order <- geneNMF.metaprograms$programs.tree$order
orig_clusters <- geneNMF.metaprograms$programs.clusters[orig_tree_order]
orig_order <- paste0("MP", unique(orig_clusters))
orig_order <- orig_order[orig_order %in% retained_mps]

reordered_mps <- c(
  orig_order[orig_order %in% cc_mps],
  state_ordered_mps[state_ordered_mps %in% retained_mps]
)
# Add any remaining ones that might not be in state_groups or cc_mps
reordered_mps <- unique(c(reordered_mps, orig_order))
mp_tree_order_names <- reordered_mps

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))
state_level_order <- c(ordered_group_names, "Unresolved", "Hybrid")

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

common_cells <- intersect(rownames(ucell_scores), Cells(tmdata_all))
tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]

sample_var <- tmdata_all$orig.ident
study_var <- tmdata_all$study

cna_cells <- intersect(rownames(meta_full_epi), common_cells)
cna_status <- rep(NA_character_, length(common_cells))
names(cna_status) <- common_cells
cna_status[cna_cells] <- as.character(meta_full_epi[cna_cells, "classification"])

cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))
names(cc_score) <- colnames(tmdata_all)

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

assign_states_mode <- function(mode = c("reg", "noreg")) {
  mode <- match.arg(mode)

  retained_in_ucell <- intersect(retained_mps, colnames(ucell_scores))
  cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
  non_cc_in_ucell <- intersect(non_cc_mps, colnames(ucell_scores))

  ucell_mat <- as.matrix(ucell_scores[, retained_in_ucell, drop = FALSE])
  X_cc <- ucell_mat[, cc_in_ucell, drop = FALSE]
  Y_other <- ucell_mat[, non_cc_in_ucell, drop = FALSE]

  if (mode == "reg") {
    X <- cbind(Intercept = 1, X_cc)
    XtX_inv <- solve(crossprod(X))
    B <- XtX_inv %*% crossprod(X, Y_other)
    Y_hat <- X %*% B
    Y_use <- Y_other - Y_hat
  } else {
    Y_use <- Y_other
  }

  mp_adj_noncc <- z_normalise(Y_use, sample_var, study_var)
  cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
  mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
  mp_adj_all <- cbind(mp_adj_noncc, mp_adj_cc)

  group_max <- sapply(state_groups, function(mps) {
    mps_avail <- intersect(mps, colnames(mp_adj_noncc))
    if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[, mps_avail]))
    apply(mp_adj_noncc[, mps_avail, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(mp_adj_noncc)

  THRESHOLD <- 0.5
  HYBRID_GAP_B <- 0.3

  best_group_idx <- max.col(group_max, ties.method = "first")
  best_group_val <- apply(group_max, 1, max)
  base_state <- names(state_groups)[best_group_idx]
  base_state[best_group_val < THRESHOLD] <- "Unresolved"

  sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
  gap <- sorted_groups[, 1] - sorted_groups[, 2]
  state_B <- base_state
  state_B[(gap < HYBRID_GAP_B) & (base_state != "Unresolved")] <- "Hybrid"
  names(state_B) <- rownames(group_max)

  list(state = state_B, mp_adj_noncc = mp_adj_noncc, mp_adj_all = mp_adj_all, group_max = group_max)
}

make_state_heatmap <- function(state_vec, mp_adj_all, mode_label) {
  set.seed(42)
  MAX_CELLS_TOTAL <- 8000
  state_counts <- table(state_vec)
  state_fracs <- state_counts / sum(state_counts)
  cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)
  state_cells <- split(names(state_vec), state_vec)
  cells_to_plot <- unlist(
    mapply(function(cells, n) sample(cells, min(length(cells), n)),
           state_cells,
           cells_per_state[names(state_cells)],
           SIMPLIFY = FALSE),
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

  present_states <- intersect(state_level_order, unique(as.character(state_vec[cells_to_plot])))
  if (length(present_states) == 0) present_states <- unique(as.character(state_vec[cells_to_plot]))
  split_vec <- factor(as.character(state_vec[cells_to_plot]), levels = present_states)

  state_df_full <- data.frame(
    cell = names(state_vec),
    state = as.character(state_vec),
    sample = as.character(tmdata_all@meta.data[names(state_vec), "orig.ident"]),
    study = as.character(tmdata_all@meta.data[names(state_vec), "study"]),
    stringsAsFactors = FALSE
  )
  total_samples <- length(unique(state_df_full$sample))
  total_studies <- length(unique(state_df_full$study))
  state_div_df <- state_df_full %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(
      sample_cov = dplyr::n_distinct(sample) / max(total_samples, 1),
      study_cov = dplyr::n_distinct(study) / max(total_studies, 1),
      diversity_score = 0.5 * sample_cov + 0.5 * study_cov,
      .groups = "drop"
    )
  state_div_map <- setNames(state_div_df$diversity_score, state_div_df$state)
  state_div_vals <- state_div_map[as.character(split_vec)]
  state_div_vals[is.na(state_div_vals)] <- 0
  names(state_div_vals) <- cells_to_plot

  local_group_cols <- group_cols[names(group_cols) %in% present_states]
  for (st in present_states) if (!st %in% names(local_group_cols)) local_group_cols[st] <- "grey80"
  local_group_cols <- local_group_cols[present_states]

  study_vals <- tmdata_all@meta.data[cells_to_plot, "study"]
  study_cols <- setNames(
    DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"),
    unique(tmdata_all$study)
  )

  max_cc <- max(cc_score[cells_to_plot], na.rm = TRUE)

  col_ann <- HeatmapAnnotation(
    State = split_vec,
    CNA = cna_status[cells_to_plot],
    CC_score = cc_score[cells_to_plot],
    Diversity = state_div_vals,
    Study = study_vals,
    col = list(
      State = local_group_cols,
      CNA = c(cna_malignant = "black", cna_unresolved = "grey70"),
      CC_score = colorRamp2(c(0, max_cc), c("white", "darkgreen")),
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
  group_colors_row <- c(group_cols[ordered_group_names], Cell_cycle = "gold", Other = "grey70")
  mp_group_label <- mp_to_group
  names(mp_group_label) <- rownames(sub_scores)

  row_ann <- rowAnnotation(
    MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
    col = list(MP_group = group_colors_row),
    show_annotation_name = FALSE
  )

  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

  Heatmap(
    sub_scores,
    name = "Adj score",
    col = col_fun_sc,
    top_annotation = col_ann,
    left_annotation = row_ann,
    column_split = split_vec,
    column_order = (function() {
      col_order_list <- lapply(levels(split_vec), function(lvl) {
        idx <- which(as.character(split_vec) == lvl)
        if (length(idx) <= 1) return(idx)
        mat_lvl <- sub_scores[, idx, drop = FALSE]
        dcols <- dist(t(mat_lvl))
        hc <- hclust(dcols, method = "ward.D2")
        idx[hc$order]
      })
      full_ord <- unlist(col_order_list, use.names = FALSE)
      if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
        return(seq_len(ncol(sub_scores)))
      }
      full_ord
    })(),
    column_gap = unit(1.5, "mm"),
    row_split = factor(ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
                       levels = c("Cell_cycle_MPs", "Other_MPs")),
    row_gap = unit(2.5, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 5,
    border = FALSE,
    rect_gp = gpar(col = NA)
  )
}

make_prop_and_pie <- function(state_vec, mode_label) {
  prop_df <- data.frame(
    state = as.character(state_vec),
    study = as.character(tmdata_all@meta.data[names(state_vec), "study"]),
    stringsAsFactors = FALSE
  )
  overall <- prop_df %>% count(state) %>% mutate(study = "Overall", pct = 100 * n / sum(n))
  per_study <- prop_df %>%
    count(study, state) %>%
    group_by(study) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  plot_df <- bind_rows(overall, per_study)
  plot_df$state <- factor(plot_df$state, levels = state_level_order)

  p_bar <- ggplot(plot_df, aes(study, pct, fill = state)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.6) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    labs(title = paste0("Approach B ", mode_label, ": state proportions"), x = NULL, y = "% of cells", fill = "State") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pie_df <- overall %>% mutate(label = paste0(state, "\n", sprintf("%.1f%%", pct)))
  p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    labs(title = paste0("Approach B ", mode_label, ": overall pie"), fill = "State") +
    theme_void(base_size = 11)

  list(bar = p_bar, pie = p_pie, overall = overall)
}

make_cc_box <- function(state_vec, mode_label) {
  cc_df <- data.frame(
    state = factor(as.character(state_vec[names(cc_score)]), levels = state_level_order),
    cc_score = as.numeric(cc_score[names(state_vec)]),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(state))

  ggplot(cc_df, aes(state, cc_score, fill = state)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.15, alpha = 0.2) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    labs(title = paste0("Approach B ", mode_label, ": cell-cycle score by state"), x = NULL, y = "Cell-cycle score") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
}

results <- list()
if ("reg" %in% requested_modes) results$reg <- assign_states_mode("reg")
if ("noreg" %in% requested_modes) results$noreg <- assign_states_mode("noreg")

if (!is.null(results$reg)) {
  saveRDS(results$reg$state, "Auto_topmp_v2_reg_states_B.rds")
  saveRDS(results$reg$mp_adj_noncc, "Auto_topmp_v2_reg_mp_adj.rds")
  saveRDS(results$reg$group_max, "Auto_topmp_v2_reg_group_max.rds")
}
if (!is.null(results$noreg)) {
  saveRDS(results$noreg$state, "Auto_topmp_v2_noreg_states_B.rds")
  saveRDS(results$noreg$mp_adj_noncc, "Auto_topmp_v2_noreg_mp_adj.rds")
  saveRDS(results$noreg$group_max, "Auto_topmp_v2_noreg_group_max.rds")
}

pdf("Auto_topmp_v2_reg_noreg_heatmap_B_cconly.pdf", width = 18, height = 8, useDingbats = FALSE)
if (!is.null(results$reg)) {
  draw(make_state_heatmap(results$reg$state, results$reg$mp_adj_all, "reg"), merge_legend = TRUE)
  grid.text("Approach B - Regressed", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"), 
            just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
}
if (!is.null(results$noreg)) {
  draw(make_state_heatmap(results$noreg$state, results$noreg$mp_adj_all, "noreg"), merge_legend = TRUE)
  grid.text("Approach B - No regress", x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"), 
            just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
}
dev.off()

prop_panels <- list()
pie_panels <- list()
cc_panels <- list()
if (!is.null(results$reg)) {
  p <- make_prop_and_pie(results$reg$state, "reg")
  prop_panels[["reg"]] <- p$bar
  pie_panels[["reg"]] <- p$pie
  cc_panels[["reg"]] <- make_cc_box(results$reg$state, "reg")
}
if (!is.null(results$noreg)) {
  p <- make_prop_and_pie(results$noreg$state, "noreg")
  prop_panels[["noreg"]] <- p$bar
  pie_panels[["noreg"]] <- p$pie
  cc_panels[["noreg"]] <- make_cc_box(results$noreg$state, "noreg")
}
if (length(prop_panels) == 2) {
  ggsave("Auto_topmp_v2_reg_noreg_proportion_B_withpie.pdf", (prop_panels[[1]] + prop_panels[[2]]) / (pie_panels[[1]] + pie_panels[[2]]), width = 18, height = 12)
  ggsave("Auto_topmp_v2_reg_noreg_ccscore_boxplot_B.pdf", cc_panels[[1]] + cc_panels[[2]], width = 16, height = 6)
} else if (length(prop_panels) == 1) {
  ggsave("Auto_topmp_v2_reg_noreg_proportion_B_withpie.pdf", prop_panels[[1]] + pie_panels[[1]], width = 16, height = 7)
  ggsave("Auto_topmp_v2_reg_noreg_ccscore_boxplot_B.pdf", cc_panels[[1]], width = 10, height = 6)
}

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

sum_rows <- list()
if (!is.null(results$reg)) {
  pr <- make_prop_and_pie(results$reg$state, "reg")
  sum_rows[["reg"]] <- pr$overall %>% transmute(mode = "reg", state, cells = n, pct)
}
if (!is.null(results$noreg)) {
  pn <- make_prop_and_pie(results$noreg$state, "noreg")
  sum_rows[["noreg"]] <- pn$overall %>% transmute(mode = "noreg", state, cells = n, pct)
}
write.csv(bind_rows(sum_rows), file.path(summary_dir, "Auto_topmp_v2_reg_noreg_summary.csv"), row.names = FALSE)

message("Saved unified reg+noreg Approach B outputs and comparison figures.")
