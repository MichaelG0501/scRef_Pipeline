####################
# Auto_sample_abundance.R
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
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
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
dir.create("sample_abundance/", recursive = TRUE, showWarnings = FALSE)

summary_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/updates/new_updates/summaries"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# load data
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")

####################
# constants
####################
mp_descriptions <- c(
  "MP1" = "G2M_cycle", "MP2" = "MYC_prolif", "MP5" = "IFN_response",
  "MP7" = "S_cycle", "MP8" = "Intestinal_diff", "MP9" = "G1S_cycle",
  "MP10" = "Columnar_diff", "MP12" = "Neuro_epithelial", "MP13" = "Partial_EMT",
  "MP14" = "Hypoxia_epithelial", "MP15" = "T_NK_infiltration",
  "MP16" = "Secretory_diff", "MP17" = "Squamous_transition", "MP18" = "Adaptive_secretory"
)

state_groups <- list(
  Classic_Proliferative = c("MP2"),
  Barretts_Metaplasia = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  EMT_related = c("MP13", "MP12"),
  Intestinal_Metaplasia = c("MP18", "MP16"),
  Immune_Infiltrated = c("MP15")
)

group_cols <- c(
  Classic_Proliferative = "#E41A1C", Barretts_Metaplasia = "#4DAF4A",
  EMT_related = "#984EA3", Intestinal_Metaplasia = "#FF7F00",
  Immune_Infiltrated = "#377EB8", Unresolved = "grey80", Hybrid = "black"
)

mp_cols <- c(
  "MP1_G2M_cycle" = "#B0B0B0",
  "MP2_MYC_prolif" = "#E41A1C",
  "MP5_IFN_response" = "#66C2A5",
  "MP7_S_cycle" = "#999999",
  "MP8_Intestinal_diff" = "#FC8D62",
  "MP9_G1S_cycle" = "#C0C0C0",
  "MP10_Columnar_diff" = "#A6D854",
  "MP12_Neuro_epithelial" = "#E78AC3",
  "MP13_Partial_EMT" = "#984EA3",
  "MP14_Hypoxia_epithelial" = "#8DA0CB",
  "MP15_T_NK_infiltration" = "#377EB8",
  "MP16_Secretory_diff" = "#FFD92F",
  "MP17_Squamous_transition" = "#4DAF4A",
  "MP18_Adaptive_secretory" = "#FF7F00"
)

####################
# metaprogram filtering
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
retained_mps <- names(mp.genes)

cc_mps <- c("MP1", "MP7", "MP9")
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
  df$label <- factor(df$label, levels = all_labels)

  out <- df %>%
    count(orig.ident, label, .drop = FALSE) %>%
    complete(orig.ident, label = factor(all_labels, levels = all_labels), fill = list(n = 0)) %>%
    group_by(orig.ident) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(label = as.character(label))

  out
}

plot_abundance <- function(prop_data, sample_order, col_map, title_text, totals_df) {
  plot_df <- prop_data %>%
    filter(orig.ident %in% sample_order) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  totals_plot <- totals_df %>%
    filter(orig.ident %in% sample_order) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  scale_factor <- max(totals_plot$total_n, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

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
    scale_fill_manual(values = col_map, drop = FALSE) +
    scale_y_continuous(
      name = "Proportion (%)",
      limits = c(0, 100),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cell Count (N)", labels = comma)
    ) +
    labs(title = title_text, x = NULL, fill = "Label") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      legend.text = element_text(size = 7)
    ) +
    guides(fill = guide_legend(ncol = 4))
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

cc_mps <- c("MP1", "MP7", "MP9")
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

target_states <- c("Classic_Proliferative", "Barretts_Metaplasia", "EMT_related", "Intestinal_Metaplasia", "Immune_Infiltrated")

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
noncc_levels <- label_mp(noncc_avail)
all_levels <- label_mp(all_mps_avail)
state_levels <- c(names(state_groups), "Unresolved", "Hybrid")

prop_type1 <- make_prop_data(topmp_noncc_label, sample_by_cell[names(topmp_noncc_label)], noncc_levels)
prop_type2 <- make_prop_data(topmp_all_label, sample_by_cell[names(topmp_all_label)], all_levels)
prop_type3 <- make_prop_data(state_label, sample_by_cell[names(state_label)], state_levels)

col_type1 <- mp_cols[noncc_levels]
col_type2 <- mp_cols[all_levels]
col_type3 <- group_cols[state_levels]

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
  title_text = "Type 1: Top MP (non-CC MPs) | Sort: Diversity (Geometric Mean)",
  totals_df = totals_df
)

p2 <- plot_abundance(
  prop_data = prop_type1,
  sample_order = study_order,
  col_map = col_type1,
  title_text = "Type 1: Top MP (non-CC MPs) | Sort: Study + Alphabetical",
  totals_df = totals_df
)

p3 <- plot_abundance(
  prop_data = prop_type2,
  sample_order = diversity_order,
  col_map = col_type2,
  title_text = "Type 2: Top MP (all MPs including CC) | Sort: Diversity (Geometric Mean)",
  totals_df = totals_df
)

p4 <- plot_abundance(
  prop_data = prop_type2,
  sample_order = study_order,
  col_map = col_type2,
  title_text = "Type 2: Top MP (all MPs including CC) | Sort: Study + Alphabetical",
  totals_df = totals_df
)

p5 <- plot_abundance(
  prop_data = prop_type3,
  sample_order = diversity_order,
  col_map = col_type3,
  title_text = "Type 3: Approach B States (noreg) | Sort: Diversity (Geometric Mean)",
  totals_df = totals_df
)

p6 <- plot_abundance(
  prop_data = prop_type3,
  sample_order = study_order,
  col_map = col_type3,
  title_text = "Type 3: Approach B States (noreg) | Sort: Study + Alphabetical",
  totals_df = totals_df
)

pdf("sample_abundance/Auto_sample_abundance.pdf", width = pdf_w, height = pdf_h, onefile = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
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

message("Saved: sample_abundance/Auto_sample_abundance.pdf")
message("Saved: updates/new_updates/summaries/Auto_sample_abundance_summary.csv")
