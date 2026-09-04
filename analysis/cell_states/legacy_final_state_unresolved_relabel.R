####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_final_state_unresolved_relabel.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

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
#   ref_outs/unresolved_states/Auto_unresolved_relabel_cc_boxplot.pdf
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
task_prefix <- "task4"
out_dir <- paste0(task_prefix, "_unresolved_states")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

####################
# constants
####################
state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "Immune Infiltrating"   = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Stress-adaptive"       = "#984EA3",
  "Immune Infiltrating"   = "#377EB8",
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
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
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

# MP tree and group order logic (moved up for earlier use)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
retained_mps <- names(mp.genes)
# Reorder MPs: CC MPs first (original relative order), then state MPs (state_groups order)
state_ordered_mps <- unlist(state_groups, use.names = FALSE)
# To keep original relative order of CC, we identify their current sequence
orig_tree_order <- geneNMF.metaprograms$programs.tree$order
orig_clusters <- geneNMF.metaprograms$programs.clusters[orig_tree_order]
orig_order <- paste0("MP", unique(orig_clusters))
orig_order <- orig_order[orig_order %in% retained_mps]

reordered_mps <- c(
  orig_order[orig_order %in% cc_mps],
  state_ordered_mps[state_ordered_mps %in% retained_mps]
)
# Add any remaining ones (like 3CA ones if they were added to retained_mps, 
# although they are typically handled separately in the rows below)
reordered_mps <- unique(c(reordered_mps, orig_order))
mp_tree_order_names <- reordered_mps

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

tcga_meta_path <- "tcga_esca_meta.rds"
tcga_tpm_path <- "cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
has_tcga_inputs <- file.exists(tcga_meta_path) && file.exists(tcga_tpm_path)
if (has_tcga_inputs) {
  meta_tcga <- readRDS(tcga_meta_path)
  tpm_df <- data.table::fread(tcga_tpm_path)
  tpm_mat <- as.matrix(tpm_df[, -1])
  rownames(tpm_mat) <- tpm_df$GeneSymbol
} else {
  meta_tcga <- NULL
  tpm_mat <- NULL
  message("TCGA survival inputs not found; unresolved relabeling will skip STEP 7 survival volcano.")
}

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
total_cells <- length(common_cells)   # or ncol(tmdata_all)

mp_coverage <- unresolved_meta %>%
  group_by(mp_label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(orig.ident),
    n_studies = n_distinct(study),
    pct_cells = 100 * n() / total_cells,
    pct_samples = 100 * n_distinct(orig.ident) / total_samples,
    pct_studies = 100 * n_distinct(study) / total_studies,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

message("Pan-cancer MP coverage in unresolved cells:")
print(mp_coverage)

candidates <- mp_coverage %>%
  filter(
    n_samples >= 50,
    n_studies >= 6,
    pct_cells >= 1
  ) %>%
  arrange(desc(n_cells))

retained_3ca <- candidates %>% pull(mp_label)

message(paste("Retained pan-cancer MPs:", paste(retained_3ca, collapse = ", ")))
write.csv(mp_coverage, file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_mp_coverage.csv")), row.names = FALSE)

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

# Apply Merges
# 1. Merge 3CA_mp_30 Respiration 1 (from unresolved) with Classic Proliferative
state_updated[state_updated == "3CA_mp_30 Respiration 1"] <- "Classic Proliferative"

# 2. Merge 3CA_mp_12 Protein maturation and 3CA_mp_17 EMT III into 3CA_EMT_and_Protein_maturation
state_updated[state_updated %in% c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")] <- "3CA_EMT_and_Protein_maturation"

# Define final levels and colors
new_retained_3ca_names <- clean_3ca_name(retained_3ca)
# Remove the ones that were merged into others
new_retained_3ca_names <- setdiff(new_retained_3ca_names, c("3CA_mp_30 Respiration 1", "3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III"))
# Add the new merged one if its inputs were present in candidate 3CA MPs
if (any(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III") %in% clean_3ca_name(retained_3ca))) {
  new_retained_3ca_names <- unique(c(new_retained_3ca_names, "3CA_EMT_and_Protein_maturation"))
}

state_level_order_updated <- c(
  ordered_group_names,
  new_retained_3ca_names,
  "Unresolved", "Hybrid"
)

# Custom colors for new states to look premium
custom_3ca_cols <- c(
  "3CA_EMT_and_Protein_maturation" = "#666666", # Dark grey/brown for EMT-maturation
  "3CA_mp_1 Epithelial-1"          = "#F781BF", # Pink
  "3CA_mp_5 Epithelial-5"          = "#A65628", # Brown
  "3CA_mp_21 Epithelial-21"         = "#FFFF33"  # Yellow
)

# For any others, use a palette
remaining_3ca <- setdiff(new_retained_3ca_names, names(custom_3ca_cols))
if (length(remaining_3ca) > 0) {
  new_state_cols <- c(
    custom_3ca_cols[intersect(names(custom_3ca_cols), new_retained_3ca_names)],
    setNames(scales::hue_pal()(length(remaining_3ca)), remaining_3ca)
  )
} else {
  new_state_cols <- custom_3ca_cols[intersect(names(custom_3ca_cols), new_retained_3ca_names)]
}

group_cols_updated <- c(group_cols, new_state_cols)

saveRDS(state_updated, file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_states.rds")))
saveRDS(state_updated, "Auto_final_states.rds")

####################
# STEP 4: updated proportion plot
####################
make_prop_and_pie <- function(state_vec, mode_label, lvls, cols) {
  prop_df <- data.frame(
    state = as.character(state_vec),
    study = as.character(tmdata_all$study[common_cells]),
    stringsAsFactors = FALSE
  )
  overall <- prop_df %>% count(state) %>% mutate(study = "Total", pct = 100 * n / sum(n))
  per_study <- prop_df %>%
    count(study, state) %>%
    group_by(study) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  
  plot_df <- bind_rows(per_study, overall)
  plot_df$state <- factor(plot_df$state, levels = lvls)
  
  study_levels <- c(sort(unique(per_study$study)), "Total")
  plot_df$study <- factor(plot_df$study, levels = study_levels)
  plot_df$is_total <- factor(ifelse(plot_df$study == "Total", "Total", "Studies"), levels = c("Studies", "Total"))

  p_bar <- ggplot(plot_df, aes(study, pct, fill = state)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct > 3, sprintf("%.1f%%", pct), "")), position = position_stack(vjust = 0.5), size = 4.5) +
    scale_fill_manual(values = cols, drop = FALSE) +
    facet_grid(~ is_total, scales = "free_x", space = "free_x") +
    labs(title = paste0(mode_label, ": state proportions"), x = NULL, y = "% of cells", fill = "State") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.title.y = element_text(size = 16, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      panel.spacing = unit(1, "lines")
    )

  pie_df <- overall %>% mutate(label = paste0(state, "\n", sprintf("%.1f%%", pct)))
  p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = ifelse(pct > 3, label, "")), position = position_stack(vjust = 0.5), size = 5, fontface = "bold") +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(title = paste0(mode_label, ": overall pie"), fill = "State") +
    theme_void(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  list(bar = p_bar, pie = p_pie)
}

p_after <- make_prop_and_pie(state_updated[common_cells], "Updated state proportions (5 original + pan-cancer relabeled)", state_level_order_updated, group_cols_updated)

ggsave(
  file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_proportion.pdf")),
  plot = p_after$bar / p_after$pie,
  width = 10.55,
  height = 18
)

####################
# STEP 5: cell-cycle boxplot
####################
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]
cc_consensus <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))

cc_df <- data.frame(
  state = factor(as.character(state_updated[names(cc_score)]), levels = state_level_order_updated),
  cc_score = as.numeric(cc_score),
  stringsAsFactors = FALSE
) %>% filter(!is.na(state))

p_cc <- ggplot(cc_df, aes(state, cc_score, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.2) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Finalised states: Cell-cycle score by state", x = NULL, y = "Cell-cycle score") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

ggsave(
  file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_cc_boxplot.pdf")),
  p_cc,
  width = 12,
  height = 7
)

####################
# STEP 6: per-cell heatmap (all cells)
####################
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc)

# normalized pan-cancer MP expression (only retained threshold MPs)
retained_3ca <- retained_3ca[retained_3ca %in% colnames(ucell_3ca)]
z_3ca_all <- z_normalise(ucell_3ca[common_cells, retained_3ca, drop = FALSE], sample_var, study_var)
# raw pan-cancer MP UCell scores (same retained MPs)
raw_3ca_all <- as.matrix(ucell_3ca[common_cells, retained_3ca, drop = FALSE])

# CNA status logic unchanged
cna_cells <- intersect(rownames(meta_full_epi), common_cells)
cna_status <- rep(NA_character_, length(common_cells))
names(cna_status) <- common_cells
cna_status[cna_cells] <- as.character(meta_full_epi[cna_cells, "classification"])

# Logic moved up to constants/loading section

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

# append normalized pan-cancer rows (retained order)
pan_mat <- t(z_3ca_all[cells_to_plot, , drop = FALSE])
# keep retained order, no extra sorting
pan_mat <- pan_mat[retained_3ca[retained_3ca %in% rownames(pan_mat)], , drop = FALSE]
rownames(pan_mat) <- clean_3ca_name(rownames(pan_mat))

# append raw pan-cancer rows (retained order)
raw_pan_mat <- t(raw_3ca_all[cells_to_plot, , drop = FALSE])
raw_pan_mat <- raw_pan_mat[retained_3ca[retained_3ca %in% rownames(raw_pan_mat)], , drop = FALSE]
rownames(raw_pan_mat) <- paste0(clean_3ca_name(rownames(raw_pan_mat)), " (raw)")

sub_scores <- rbind(sub_scores, pan_mat, raw_pan_mat)

# keep raw row IDs for robust MP-group annotation before display renaming
row_ids_raw <- rownames(sub_scores)

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
  na_col = "white",
  annotation_legend_param = list(
    State = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CNA = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CC_score = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Diversity = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Study = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
  )
)

mp_to_group <- rep("Other", length(row_ids_raw))
names(mp_to_group) <- row_ids_raw
mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
  mp_to_group[grp_mps] <- grp
}
pan_rows <- clean_3ca_name(retained_3ca)
pan_rows_raw <- paste0(pan_rows, " (raw)")
mp_to_group[intersect(names(mp_to_group), pan_rows)] <- "PanCancer"
mp_to_group[intersect(names(mp_to_group), pan_rows_raw)] <- "PanCancerRaw"

# Define labels and colors mapping
mp_group_names <- c(
  "Cell_cycle" = "Cell cycle",
  "Classic Proliferative" = "Classic\nProlif",
  "Basal to Intestinal Metaplasia" = "Basal-IM",
  "Stress-adaptive" = "Stress\nadaptive",
  "SMG-like Metaplasia" = "SMG-like\nMeta",
  "Immune Infiltrating" = "Immune\nInfiltrating",
  "PanCancer" = "Pan-Cancer\n(norm)",
  "PanCancerRaw" = "Pan-Cancer\n(raw)",
  "Other" = "Other"
)

# Colors according to labels
group_colors_row <- c(
  setNames(group_cols[names(state_groups)], mp_group_names[names(state_groups)]),
  "Cell cycle" = "gold",
  "Pan-Cancer\n(norm)" = "#4B0082",
  "Pan-Cancer\n(raw)" = "#7F0000",
  "Other" = "grey70"
)

# Map labels using original IDs, then rename the vector to match current sub_scores rows
mp_group_labels <- mp_group_names[mp_to_group[row_ids_raw]]
names(mp_group_labels) <- mp_label_map[row_ids_raw]

# Split heatmap into Adjusted and Raw sections
adj_rows <- rownames(sub_scores)[!grepl("\\(raw\\)$", rownames(sub_scores))]
raw_rows <- rownames(sub_scores)[grepl("\\(raw\\)$", rownames(sub_scores))]

# Separate row annotations for each segment
row_ann_adj <- rowAnnotation(
  MP_group = factor(mp_group_labels[adj_rows], levels = mp_group_names),
  col = list(MP_group = group_colors_row),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(MP_group = list(title = "MP Group", 
                                                 title_gp = gpar(fontsize = 9, fontface = "bold"), 
                                                 labels_gp = gpar(fontsize = 8)))
)

row_ann_raw <- rowAnnotation(
  MP_group = factor(mp_group_labels[raw_rows], levels = mp_group_names),
  col = list(MP_group = group_colors_row),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

# Row split for adjusted part
# use the display names for matching
cc_labels <- mp_descriptions[cc_mps]
pan_labels <- clean_3ca_name(retained_3ca)

row_split_adj <- factor(
  ifelse(adj_rows %in% pan_labels,
         "Pan-Cancer\n(norm)",
         ifelse(adj_rows %in% cc_labels, "Cell cycle\nMPs", "State\nMPs")),
  levels = c("Cell cycle\nMPs", "State\nMPs", "Pan-Cancer\n(norm)")
)
  
  # Column order logic from state_definition_approach_b_reg_noreg.R
  col_order <- (function() {
    col_order_list <- lapply(levels(split_vec), function(lvl) {
      idx <- which(as.character(split_vec) == lvl)
      if (length(idx) <= 1) return(idx)
      mat_lvl <- sub_scores[adj_rows, idx, drop = FALSE]
      dcols <- dist(t(mat_lvl))
      hc <- hclust(dcols, method = "ward.D2")
      idx[hc$order]
    })
    full_ord <- unlist(col_order_list, use.names = FALSE)
    if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
      return(seq_len(ncol(sub_scores)))
    }
    full_ord
  })()

  ht_adj <- Heatmap(
    sub_scores[adj_rows, , drop = FALSE],
    name = "Adj score",
    col = col_fun_sc,
    top_annotation = col_ann,
    left_annotation = row_ann_adj,
    column_split = split_vec,
    row_split = row_split_adj,
    column_order = col_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_gap = unit(2, "mm"),
    column_gap = unit(1.5, "mm"),
    use_raster = TRUE,
    raster_quality = 3,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), 
                                labels_gp = gpar(fontsize = 8))
  )
  
  # Raw scores color scale
  raw_lim <- as.numeric(quantile(sub_scores[raw_rows, ], 0.98, na.rm = TRUE))
  if (is.na(raw_lim) || raw_lim == 0) raw_lim <- 0.15
  col_fun_raw <- colorRamp2(c(0, raw_lim), c("white", "firebrick3"))
  
  ht_raw <- Heatmap(
    sub_scores[raw_rows, , drop = FALSE],
    name = "Raw score",
    col = col_fun_raw,
    left_annotation = row_ann_raw,
    column_split = split_vec,
    column_order = col_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_gap = unit(1.5, "mm"),
    use_raster = TRUE,
    raster_quality = 3,
    show_heatmap_legend = TRUE,
    column_title = NULL,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), 
                                labels_gp = gpar(fontsize = 8))
  )

  ht <- ht_adj %v% ht_raw

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_heatmap.pdf")), width = 18, height = 10, useDingbats = FALSE)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
grid.text("Unresolved cells subclass: 3CA-based relabeling (noreg)", x = unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"), 
          just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

####################
# STEP 7: survival volcano using GSVA
####################
if (has_tcga_inputs) {
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

# Rename them using the clean function
names(new_state_sigs) <- clean_3ca_name(names(new_state_sigs))

# Handle merges for GSVA sets (same as state labels)
# Keep originals in new_state_sigs so they can be aggregated into groups, 
# but they will be excluded from the final standalone Cox columns.

# EMT + Protein Maturation merge into a new set
if (any(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III") %in% names(new_state_sigs))) {
  merged_genes <- unique(unlist(new_state_sigs[intersect(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III"), names(new_state_sigs))]))
  new_state_sigs[["3CA_EMT_and_Protein_maturation"]] <- merged_genes
}

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
# Update state groups to include Respiration for Classic Proliferative
local_state_groups <- state_groups
local_state_groups[["Classic Proliferative"]] <- c(local_state_groups[["Classic Proliferative"]], "3CA_mp_30 Respiration 1")

for (nm in names(local_state_groups)) {
  mps <- intersect(local_state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}

# 6. Combine original states and the newly added 3CA states for Cox regression
# new_state_names in Step 3 was updated to new_retained_3ca_names
final_3ca_names <- setdiff(clean_3ca_name(retained_3ca), c("3CA_mp_30 Respiration 1", "3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III"))
if (any(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III") %in% clean_3ca_name(retained_3ca))) {
  final_3ca_names <- unique(c(final_3ca_names, "3CA_EMT_and_Protein_maturation"))
}

state_cols <- intersect(c(names(local_state_groups), final_3ca_names), colnames(surv_data))

all_cox <- list()

# Plot both EAC and ESCC
p_list <- list()
for (coh in c("EAC", "ESCC")) {
  cox_df <- run_cox_for_group(
    surv_data %>% filter(HistologyGroup == coh),
    state_cols,
    cohort_name = coh
  )
  if (nrow(cox_df) > 0) {
    all_cox[[coh]] <- cox_df
    p_list[[coh]] <- plot_volcano(cox_df, paste0("Finalised states: TCGA survival (", coh, ")"))
  }
}

if (length(p_list) > 0) {
  pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_volcano.pdf")), width = 16, height = 7)
  print(wrap_plots(p_list, ncol = 2))
  dev.off()
}

cox_res <- bind_rows(all_cox)
write.csv(cox_res, file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_cox_results.csv")), row.names = FALSE)
} else {
  write.csv(
    data.frame(),
    file.path(out_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_cox_results.csv")),
    row.names = FALSE
  )
}

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_df <- data.frame(
  task = task_prefix,
  n_cells_total = length(common_cells),
  n_unresolved_input = length(unresolved_cells),
  n_retained_3ca = length(retained_3ca),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(summary_dir, paste0("Auto_", task_prefix, "_unresolved_relabel_summary.csv")), row.names = FALSE)

message(sprintf("Saved unresolved relabel outputs in %s", out_dir))
