####################
# Auto_unresolved_relabel.R
#
# Relabel unresolved Approach-B (noreg) cells by top pan-cancer (3CA) MP,
# retain 3-5 pan-cancer MPs by sample/study coverage thresholds, and regenerate:
#   - updated proportions (bar + pie)
#   - per-cell MP heatmap (all cells; split by updated states)
#   - TCGA survival volcano plots (EAC, ESCC)
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#   ref_outs/UCell_3CA_MPs.rds
#   ref_outs/meta_full_epi.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#
# Output:
#   ref_outs/unresolved_states/Auto_unresolved_relabel_states.rds
#   ref_outs/unresolved_states/Auto_unresolved_relabel_proportion.pdf
#   ref_outs/unresolved_states/Auto_unresolved_relabel_heatmap.pdf
#   ref_outs/unresolved_states/Auto_unresolved_relabel_volcano.pdf
#   ref_outs/unresolved_states/Auto_unresolved_relabel_cox_results.csv
#   ref_outs/unresolved_states/Auto_unresolved_relabel_mp_coverage.csv
#   updates/new_updates/summaries/Auto_unresolved_relabel_summary.csv
####################

####################
# libraries
####################
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(survival)
library(survminer)
library(GSVA)
library(patchwork)
library(data.table)

####################
# setup
####################
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
dir.create("unresolved_states/", recursive = TRUE, showWarnings = FALSE)

####################
# constants
####################
state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intest. Meta" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrated"    = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrated"    = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

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
  "MP15" = "Immune Attracting",
  "MP12" = "Stressed-basal"
)
cc_mps <- c("MP1", "MP7", "MP9")
non_cc_mps <- c("MP2", "MP5", "MP8", "MP10", "MP12", "MP13", "MP14", "MP15", "MP16", "MP17", "MP18")

####################
# helper functions
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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}

run_cox_for_group <- function(df, features, cohort_name) {
  out <- list()
  for (feat in features) {
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`")), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      cohort = cohort_name,
      feature = feat,
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, title_text) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(sig = P_value < 0.05, neglog10 = -log10(P_value), log2HR = log2(HR))
  ggplot(df, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.85) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(aes(label = feature), size = 3.1, max.overlaps = 100) +
    theme_minimal(base_size = 13) +
    labs(title = title_text, x = "log2(HR)", y = "-log10(p)")
}

####################
# data loading
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")
meta_full_epi <- readRDS("meta_full_epi.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")

meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

####################
# STEP 1: relabel unresolved cells by top 3CA MP
####################
common_cells <- intersect(names(state_B), Cells(tmdata_all))
common_cells <- intersect(common_cells, rownames(ucell_3ca))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))
common_cells <- intersect(common_cells, rownames(ucell_scores))

tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]

unresolved_cells <- names(state_B)[state_B == "Unresolved"]
unresolved_3ca <- ucell_3ca[unresolved_cells, , drop = FALSE]

sample_var <- tmdata_all$orig.ident
study_var <- tmdata_all$study
names(sample_var) <- Cells(tmdata_all)
names(study_var) <- Cells(tmdata_all)

# Kept normalization here in case you still need z_3ca_unresolved for heatmap plotting
z_3ca <- z_normalise(ucell_3ca, sample_var, study_var)
z_3ca_unresolved <- z_3ca[unresolved_cells, , drop = FALSE]

CC_FIXED <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)

# remove CC_FIXED MPs from top-label calculation
unresolved_3ca_nocc <- unresolved_3ca[, !colnames(unresolved_3ca) %in% CC_FIXED, drop = FALSE]

# Calculate max column using raw UCell scores, ignoring CC_FIXED MPs
top_3ca_mp <- colnames(unresolved_3ca_nocc)[max.col(unresolved_3ca_nocc, ties.method = "first")]
names(top_3ca_mp) <- unresolved_cells

####################
# STEP 2: coverage-threshold selection of 3CA MPs (3-5)
####################
unresolved_meta <- data.frame(
  cell = unresolved_cells,
  mp_label = top_3ca_mp,
  orig.ident = as.character(tmdata_all$orig.ident[unresolved_cells]),
  study = as.character(tmdata_all$study[unresolved_cells]),
  stringsAsFactors = FALSE
)

total_samples <- length(unique(tmdata_all$orig.ident[common_cells]))
total_studies <- length(unique(tmdata_all$study[common_cells]))

mp_coverage <- unresolved_meta %>%
  group_by(mp_label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(orig.ident),
    n_studies = n_distinct(study),
    pct_samples = 100 * n_distinct(orig.ident) / total_samples,
    pct_studies = 100 * n_distinct(study) / total_studies,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

message("Pan-cancer MP coverage in unresolved cells:")
print(mp_coverage)

candidates <- mp_coverage %>%
  filter(n_samples >= 50 & n_studies >= 6) %>%
  arrange(desc(n_cells))

retained_3ca <- candidates %>% pull(mp_label)

message(paste("Retained pan-cancer MPs:", paste(retained_3ca, collapse = ", ")))
write.csv(mp_coverage, "unresolved_states/Auto_unresolved_relabel_mp_coverage.csv", row.names = FALSE)

####################
# STEP 3: update state labels
####################
state_updated <- state_B

for (cell in unresolved_cells) {
  mp <- top_3ca_mp[cell]
  if (mp %in% retained_3ca) {
    state_updated[cell] <- clean_3ca_name(mp)
  }
}

new_state_names <- unique(clean_3ca_name(retained_3ca))
state_level_order_updated <- c(
  "Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive",
  "SMG-like Metaplasia", "Immune Infiltrated",
  sort(new_state_names),
  "Unresolved", "Hybrid"
)

new_state_cols <- setNames(
  scales::hue_pal()(length(new_state_names)),
  new_state_names
)
group_cols_updated <- c(group_cols, new_state_cols)

saveRDS(state_updated, "unresolved_states/Auto_unresolved_relabel_states.rds")

####################
# STEP 4: updated proportion plot
####################
prop_df <- data.frame(
  state = as.character(state_updated[common_cells]),
  study = as.character(tmdata_all$study[common_cells]),
  stringsAsFactors = FALSE
)
overall <- prop_df %>% count(state) %>% mutate(study = "Overall", pct = 100 * n / sum(n))
per_study <- prop_df %>%
  count(study, state) %>%
  group_by(study) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
plot_df <- bind_rows(overall, per_study)
plot_df$state <- factor(plot_df$state, levels = state_level_order_updated)

p_bar <- ggplot(plot_df, aes(study, pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions (5 original + pan-cancer relabeled)", x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pie_df <- overall %>% 
  mutate(label = ifelse(pct < 5, "", paste0(state, "\n", sprintf("%.1f%%", pct))))

p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions: overall pie", fill = "State") +
  theme_void(base_size = 11)

ggsave(
  "unresolved_states/Auto_unresolved_relabel_proportion.pdf",
  p_bar + p_pie + plot_layout(widths = c(2, 1)),
  width = 18,
  height = 8
)

####################
# STEP 5: per-cell heatmap (all cells)
####################
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc)

cna_cells <- intersect(rownames(meta_full_epi), common_cells)
cna_status <- rep(NA_character_, length(common_cells))
names(cna_status) <- common_cells
cna_status[cna_cells] <- as.character(meta_full_epi[cna_cells, "classification"])

cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]
cc_consensus <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
retained_mps <- names(mp.genes)
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

set.seed(42)
MAX_CELLS_TOTAL <- 8000
state_counts <- table(state_updated)
state_fracs <- state_counts / sum(state_counts)
cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)
state_cells <- split(names(state_updated), state_updated)
cells_to_plot <- unlist(
  mapply(
    function(cells, n) sample(cells, min(length(cells), n)),
    state_cells,
    cells_per_state[names(state_cells)],
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)
if (length(cells_to_plot) > MAX_CELLS_TOTAL) {
  cells_to_plot <- sample(cells_to_plot, MAX_CELLS_TOTAL)
}

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

present_states <- intersect(state_level_order_updated, unique(as.character(state_updated[cells_to_plot])))
if (length(present_states) == 0) present_states <- unique(as.character(state_updated[cells_to_plot]))
split_vec <- factor(as.character(state_updated[cells_to_plot]), levels = present_states)

state_df_full <- data.frame(
  cell = names(state_updated),
  state = as.character(state_updated),
  sample = as.character(tmdata_all@meta.data[names(state_updated), "orig.ident"]),
  study = as.character(tmdata_all@meta.data[names(state_updated), "study"]),
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

local_group_cols <- group_cols_updated[names(group_cols_updated) %in% present_states]
for (st in present_states) if (!st %in% names(local_group_cols)) local_group_cols[st] <- "grey80"
local_group_cols <- local_group_cols[present_states]

study_vals <- tmdata_all@meta.data[cells_to_plot, "study"]
study_levels <- unique(as.character(tmdata_all$study))
study_cols <- setNames(
  DiscretePalette(length(study_levels), palette = "polychrome"),
  study_levels
)

max_cc <- max(cc_score[cells_to_plot], na.rm = TRUE)

cna_levels <- sort(unique(cna_status[!is.na(cna_status)]))
if (length(cna_levels) == 0) cna_levels <- "unknown"
cna_palette <- setNames(rep("grey70", length(cna_levels)), cna_levels)
if ("cna_malignant" %in% names(cna_palette)) cna_palette["cna_malignant"] <- "black"
if ("cna_unresolved" %in% names(cna_palette)) cna_palette["cna_unresolved"] <- "grey70"

col_ann <- HeatmapAnnotation(
  State = split_vec,
  CNA = cna_status[cells_to_plot],
  CC_score = cc_score[cells_to_plot],
  Diversity = state_div_vals,
  Study = study_vals,
  col = list(
    State = local_group_cols,
    CNA = cna_palette,
    CC_score = colorRamp2(c(0, max_cc), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Study = study_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE,
  na_col = "white"
)

ordered_group_names <- names(state_groups)
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

ht <- Heatmap(
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
  column_gap = grid::unit(1.5, "mm"),
  row_split = factor(
    ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
    levels = c("Cell_cycle_MPs", "Other_MPs")
  ),
  row_gap = grid::unit(2.5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 9, fontface = "italic"),
  show_column_names = FALSE,
  column_title_rot = 30,
  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = grid::gpar(col = NA)
)

pdf("unresolved_states/Auto_unresolved_relabel_heatmap.pdf", width = 18, height = 8, useDingbats = FALSE)
draw(ht, merge_legend = TRUE)
dev.off()

####################
# STEP 6: survival volcano using GSVA
####################
set.seed(42)

# 1. Extract original MP genes from geneNMF
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
gsva_sets <- lapply(mp.genes, unique)

# 2. Load and format 3CA pan-cancer genes
MP_df <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", check.names = FALSE)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])

# Replicate the exact naming used when creating the UCell scores
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))

# Filter to only the retained 3CA MPs
new_state_sigs <- MP_list[retained_3ca]

# Rename them using the clean function so they exactly match new_state_names
names(new_state_sigs) <- clean_3ca_name(names(new_state_sigs))

# Append the 3CA gene sets to the original GSVA sets
gsva_sets <- c(gsva_sets, new_state_sigs)

# 3. Filter genes present in TCGA and run GSVA
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]

gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

# 4. Merge with TCGA metadata
surv_data <- meta_tcga %>%
  inner_join(gsva_df, by = "sample_barcode") %>%
  filter(sample_type_code == "01")

surv_data$HistologyGroup <- infer_histology(surv_data$type)

# 5. Aggregate original MP scores into State scores (taking max MP score for each state group)
for (nm in names(state_groups)) {
  mps <- intersect(state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}

# 6. Combine original states and the newly added 3CA states for Cox regression
state_cols <- intersect(c(names(state_groups), new_state_names), colnames(surv_data))

all_cox <- list()

# ONLY plot EAC
pdf("unresolved_states/Auto_unresolved_relabel_volcano.pdf", width = 9, height = 7)
for (coh in c("EAC")) {
  cox_df <- run_cox_for_group(
    surv_data %>% filter(HistologyGroup == coh),
    state_cols,
    cohort_name = coh
  )
  all_cox[[coh]] <- cox_df
  p <- plot_volcano(cox_df, paste0("Updated states: TCGA survival volcano (", coh, ")"))
  if (!is.null(p)) print(p)
}
dev.off()

cox_res <- bind_rows(all_cox)
write.csv(cox_res, "unresolved_states/Auto_unresolved_relabel_cox_results.csv", row.names = FALSE)