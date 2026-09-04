####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/state_mp_sample_abundance.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Summarises current centred refined state and top-MP prevalence by orig.ident.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#
# Output:
#   ref_outs/sample_abundance/Auto_sample_abundance.pdf
#   updates/new_updates/summaries/Auto_sample_abundance_summary.csv
####################

####################
# libraries
####################
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

####################
# paths and setup
####################
setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")
task_prefix <- "task3"
out_dir <- paste0(task_prefix, "_sample_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/updates/new_updates/summaries"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# load data
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
mp_adj_noncc <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds")
grouping <- read.csv("Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv", check.names = FALSE)

####################
# constants
####################
mp_descriptions <- setNames(grouping$description, grouping$mp)

state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

group_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

# Guard against unexpected extra labels while preserving technical states.
extra_states <- setdiff(unique(as.character(state_B)), names(group_cols))
if (length(extra_states) > 0) {
  new_cols <- setNames(scales::hue_pal()(length(extra_states)), extra_states)
  group_cols <- c(group_cols, new_cols)
}

mp_cols <- setNames(scales::hue_pal()(nrow(grouping)), paste0(grouping$mp, "_", grouping$description))

####################
# metaprogram filtering
####################
retained_mps <- grouping$mp
mp_tree_order_names <- retained_mps
cc_mps <- c("MP1", "MP5", "MP13+")
non_cc_mps <- setdiff(retained_mps, cc_mps)

####################
# helpers
####################
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

label_mp <- function(mp_vec) {
  desc <- mp_descriptions[mp_vec]
  desc[is.na(desc)] <- mp_vec[is.na(desc)]
  out <- paste0(mp_vec, "_", desc)
  names(out) <- names(mp_vec)
  out
}

make_prop_data <- function(label_vec, sample_vec, all_labels) {
  df <- data.frame(
    orig.ident = as.character(sample_vec),
    label = as.character(label_vec),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$label), , drop = FALSE]
  df <- df[df$label %in% all_labels, , drop = FALSE]
  all_samples <- sort(unique(as.character(sample_vec)))
  counts <- df %>%
    count(orig.ident, label, name = "n")

  out <- tidyr::expand_grid(
    orig.ident = all_samples,
    label = as.character(all_labels)
  ) %>%
    left_join(counts, by = c("orig.ident", "label")) %>%
    mutate(n = ifelse(is.na(n), 0L, as.integer(n))) %>%
    group_by(orig.ident) %>%
    mutate(pct = {
      total_n <- sum(n, na.rm = TRUE)
      if (total_n > 0) 100 * n / total_n else rep(0, dplyr::n())
    }) %>%
    ungroup()

  out
}

plot_abundance <- function(prop_data, sample_order, col_map, title_text, totals_df) {
  plot_df <- prop_data %>%
    filter(orig.ident %in% sample_order) %>%
    filter(!is.na(label), !is.na(pct), is.finite(pct), pct >= 0) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  totals_plot <- totals_df %>%
    filter(orig.ident %in% sample_order) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  scale_factor <- max(totals_plot$total_n, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

  # enforce top-to-bottom bar order to follow provided label order
  stack_levels <- rev(names(col_map))
  plot_df$label <- factor(plot_df$label, levels = stack_levels)

  ggplot(plot_df, aes(x = orig.ident, y = pct, fill = label)) +
    geom_col(width = 0.75) +
    geom_point(
      data = totals_plot,
      aes(x = orig.ident, y = total_n / scale_factor, fill = NULL),
      color = "black",
      size = 1.5,
      shape = 18,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = totals_plot,
      aes(x = orig.ident, y = total_n / scale_factor, group = 1, fill = NULL),
      color = "black",
      alpha = 0.4,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = col_map, breaks = names(col_map), drop = FALSE) +
    scale_y_continuous(
      name = "Proportion (%)",
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cell Count (N)", labels = comma)
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(title = title_text, x = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      legend.text = element_text(size = 7)
    ) +
    guides(fill = guide_legend(ncol = 4, reverse = FALSE))
}

####################
# align cells and reconstruct mp_adj_all
####################
common_cells <- intersect(rownames(ucell_scores), Cells(tmdata_all))
common_cells <- intersect(common_cells, names(state_B))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))

tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
state_B <- state_B[common_cells]
mp_adj_noncc <- as.matrix(mp_adj_noncc[common_cells, , drop = FALSE])

sample_var <- tmdata_all$orig.ident
study_var <- tmdata_all$study

cc_mps <- c("MP1", "MP5", "MP13+")
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc)

####################
# per-cell labels for 3 plot types
####################
sample_by_cell <- tmdata_all$orig.ident
names(sample_by_cell) <- Cells(tmdata_all)

noncc_avail <- intersect(non_cc_mps, colnames(mp_adj_noncc))
mp_noncc_use <- mp_adj_noncc[, noncc_avail, drop = FALSE]
topmp_noncc <- colnames(mp_noncc_use)[max.col(mp_noncc_use, ties.method = "first")]
names(topmp_noncc) <- rownames(mp_noncc_use)
topmp_noncc_label <- label_mp(topmp_noncc)

all_mps_avail <- intersect(retained_mps, colnames(mp_adj_all))
mp_all_use <- mp_adj_all[, all_mps_avail, drop = FALSE]
topmp_all <- colnames(mp_all_use)[max.col(mp_all_use, ties.method = "first")]
names(topmp_all) <- rownames(mp_all_use)
topmp_all_label <- label_mp(topmp_all)

state_label <- as.character(state_B)
names(state_label) <- names(state_B)

####################
# sample totals and sort orders
####################
totals_df <- data.frame(
  orig.ident = as.character(sample_by_cell),
  stringsAsFactors = FALSE
) %>%
  count(orig.ident, name = "total_n")

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))
target_states <- ordered_group_names

state_df <- data.frame(
  cell = names(state_B),
  state = as.character(state_B),
  stringsAsFactors = FALSE
)
state_df$orig.ident <- sample_by_cell[state_df$cell]

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
  left_join(totals_df, by = "orig.ident") %>%
  arrange(desc(geo_mean_score), orig.ident) %>%
  mutate(rank = row_number())

diversity_order <- rank_df$orig.ident
diversity_order <- c(diversity_order, setdiff(sort(unique(sample_by_cell)), diversity_order))

study_map <- tmdata_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(study = first(study), .groups = "drop")

study_order <- study_map %>%
  arrange(study, orig.ident) %>%
  pull(orig.ident)

####################
# proportions for each type
####################
# Define preferred MP order: CC MPs first, then defined state MPs, then any others.
state_ordered_mps <- unname(unlist(state_groups))
cc_mps_order <- c("MP1", "MP5", "MP13+")
preferred_mp_order <- c(cc_mps_order, state_ordered_mps)
# Ensure we include any retained MPs that might not be in the explicit lists
preferred_mp_order <- c(preferred_mp_order, setdiff(retained_mps, preferred_mp_order))

noncc_tree_order <- preferred_mp_order[preferred_mp_order %in% non_cc_mps & preferred_mp_order %in% noncc_avail]
noncc_levels <- label_mp(noncc_tree_order)
all_levels <- label_mp(preferred_mp_order[preferred_mp_order %in% all_mps_avail])
state_levels <- c(ordered_group_names, "Unresolved", "Hybrid")

# ensure bars are exactly 100% over labels-of-interest only
force_100 <- function(df) {
  df %>%
    group_by(orig.ident) %>%
    mutate(pct = {
      total_n <- sum(n, na.rm = TRUE)
      if (total_n > 0) 100 * n / total_n else rep(0, dplyr::n())
    }) %>%
    ungroup()
}

prop_type1 <- make_prop_data(topmp_noncc_label, sample_by_cell[names(topmp_noncc_label)], noncc_levels)
prop_type2 <- make_prop_data(topmp_all_label, sample_by_cell[names(topmp_all_label)], all_levels)
prop_type3 <- make_prop_data(state_label, sample_by_cell[names(state_label)], state_levels)
prop_type1 <- force_100(prop_type1)
prop_type2 <- force_100(prop_type2)
prop_type3 <- force_100(prop_type3)

col_type1 <- mp_cols[noncc_levels]
col_type2 <- mp_cols[all_levels]
col_type3 <- group_cols[state_levels]

if (any(is.na(col_type1))) {
  miss1 <- names(col_type1)[is.na(col_type1)]
  col_type1[miss1] <- setNames(scales::hue_pal()(length(miss1)), miss1)
}
if (any(is.na(col_type2))) {
  miss2 <- names(col_type2)[is.na(col_type2)]
  col_type2[miss2] <- setNames(scales::hue_pal()(length(miss2)), miss2)
}

####################
# generate 6 plots and save one multi-page pdf
####################
n_samples <- length(unique(sample_by_cell))
pdf_w <- max(16, 0.15 * n_samples)
pdf_h <- 8

p1 <- plot_abundance(
  prop_data = prop_type1,
  sample_order = diversity_order,
  col_map = col_type1,
  title_text = "Type 1: Top MP (non-CC MPs) | Sort: Diversity",
  totals_df = totals_df
)

p2 <- plot_abundance(
  prop_data = prop_type1,
  sample_order = study_order,
  col_map = col_type1,
  title_text = "Type 1: Top MP (non-CC MPs) | Sort: Study",
  totals_df = totals_df
)

p3 <- plot_abundance(
  prop_data = prop_type2,
  sample_order = diversity_order,
  col_map = col_type2,
  title_text = "Type 2: Top MP (all MPs including CC) | Sort: Diversity",
  totals_df = totals_df
)

p4 <- plot_abundance(
  prop_data = prop_type2,
  sample_order = study_order,
  col_map = col_type2,
  title_text = "Type 2: Top MP (all MPs including CC) | Sort: Study",
  totals_df = totals_df
)

p5 <- plot_abundance(
  prop_data = prop_type3,
  sample_order = diversity_order,
  col_map = col_type3,
  title_text = "Type 3: Approach B States (noreg) | Sort: Diversity",
  totals_df = totals_df
)

p6 <- plot_abundance(
  prop_data = prop_type3,
  sample_order = study_order,
  col_map = col_type3,
  title_text = "Type 3: Approach B States (noreg) | Sort: Study",
  totals_df = totals_df
)

####################
# p7: Basal Metaplasia Breakdown
####################
# Population: cells in the current basal-to-intestinal-metaplasia state
# Categories: the six current defining MPs
basal_cells <- names(state_B)[state_B == "Basal to intestinal metaplasia"]
basal_topmp_label <- topmp_all_label[basal_cells]
basal_sample_by_cell <- sample_by_cell[basal_cells]

basal_mp_ids <- c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+")
basal_levels <- label_mp(basal_mp_ids)

# Proportions scaled to 100% of these six MPs within the basal population
prop_basal <- make_prop_data(basal_topmp_label, basal_sample_by_cell, basal_levels)
prop_basal <- force_100(prop_basal)

# Secondary axis: total count of basal-to-intestinal-metaplasia cells per sample
totals_basal <- data.frame(
  orig.ident = as.character(basal_sample_by_cell),
  stringsAsFactors = FALSE
) %>%
  count(orig.ident, name = "total_n") %>%
  right_join(data.frame(orig.ident = unique(sample_by_cell), stringsAsFactors = FALSE), by = "orig.ident") %>%
  mutate(total_n = ifelse(is.na(total_n), 0, total_n))

# Sorting order: proportion of MP14 within this state
mp14_label <- label_mp("MP14")
basal_sort_order <- prop_basal %>%
  filter(label == mp14_label) %>%
  arrange(desc(pct), orig.ident) %>%
  pull(orig.ident)

# Specialized palette for basal MPs to ensure high separability and NO overlap with state colors
# (States use: Red, Green, Purple, Orange, Blue)
col_basal <- setNames(
  c("#FF69B4", "#FFD700", "#00CED1", "#A0522D", "#708090", "#7B68EE"),
  basal_levels
)

p7 <- plot_abundance(
  prop_data = prop_basal,
  sample_order = basal_sort_order,
  col_map = col_basal,
  title_text = "Basal Metaplasia MP Breakdown | Sort: MP14 Proportion",
  totals_df = totals_basal
)

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_sample_abundance.pdf")), width = pdf_w, height = pdf_h, onefile = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

####################
# summary csv
####################
make_summary <- function(type_name, label_vec, label_levels) {
  tmp <- data.frame(label = as.character(label_vec), stringsAsFactors = FALSE) %>%
    count(label, name = "total_cells") %>%
    right_join(data.frame(label = label_levels, stringsAsFactors = FALSE), by = "label") %>%
    mutate(total_cells = ifelse(is.na(total_cells), 0L, as.integer(total_cells)))

  total_all <- sum(tmp$total_cells)
  tmp %>%
    transmute(
      type = type_name,
      label = label,
      total_cells = total_cells,
      pct_of_total = ifelse(total_all > 0, 100 * total_cells / total_all, 0)
    )
}

summary_df <- bind_rows(
  make_summary("type1_mp_noncc", topmp_noncc_label, noncc_levels),
  make_summary("type2_mp_all", topmp_all_label, all_levels),
  make_summary("type3_state", state_label, state_levels)
)

write.csv(
  summary_df,
  file.path(summary_dir, "Auto_sample_abundance_summary.csv"),
  row.names = FALSE
)

message(sprintf("Saved: %s", file.path(out_dir, paste0("Auto_", task_prefix, "_sample_abundance.pdf"))))
message("Saved: updates/new_updates/summaries/Auto_sample_abundance_summary.csv")
