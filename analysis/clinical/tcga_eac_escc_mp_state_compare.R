library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(GSVA)
library(patchwork)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "task8"
out_dir <- paste0(task_prefix, "_tcga_eac_escc_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}

run_gsva <- function(expr_mat, gene_sets) {
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

# Load retained 3CA MPs
coverage_path <- "task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv"
if (!file.exists(coverage_path)) coverage_path <- "unresolved_states/Auto_unresolved_relabel_mp_coverage.csv"

retained_3ca <- character(0)
if (file.exists(coverage_path)) {
  cov_df <- read.csv(coverage_path)
  retained_3ca <- cov_df %>%
    filter(n_samples >= 50, n_studies >= 6, pct_cells >= 1) %>%
    pull(mp_label)
}

# 3CA Gene sets
nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
pan_sets <- MP_list[retained_3ca]

# Combine all sets for GSVA
all_gsva_sets <- c(mp.genes, pan_sets)

tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
meta_tcga <- meta_tcga %>% 
  filter(sample_type_code == "01", HistologyGroup %in% c("EAC", "ESCC")) %>%
  mutate(HistologyGroup = factor(HistologyGroup, levels = c("EAC", "ESCC")))

mp_gs <- run_gsva(tpm_mat, all_gsva_sets)
if (is.null(mp_gs)) stop("No valid MP GSVA results")
mp_df <- as.data.frame(t(mp_gs))
mp_df$sample_barcode <- rownames(mp_df)

plot_df <- meta_tcga %>% inner_join(mp_df, by = "sample_barcode")

# State scores as max MP-in-group expression (with finalized 3CA-merge: Respiration -> Classic Prolif)
local_state_groups <- state_groups
local_state_groups[["Classic Proliferative"]] <- c(local_state_groups[["Classic Proliferative"]], "3CA_mp_30 Respiration 1")

for (nm in names(local_state_groups)) {
  mps <- intersect(local_state_groups[[nm]], colnames(plot_df))
  if (length(mps) == 0) next
  plot_df[[nm]] <- apply(as.matrix(plot_df[, mps, drop = FALSE]), 1, max)
}

# New states: just the 3CA scores (with finalized 3CA-merge: Protein Maturation + EMT)
retained_3ca_clean <- clean_3ca_name(retained_3ca)
kept_3ca_names <- setdiff(retained_3ca_clean, c("3CA_mp_30 Respiration 1", "3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III"))

# 1. As-is states
for (orig in retained_3ca) {
  cln <- clean_3ca_name(orig)
  if (cln %in% kept_3ca_names && orig %in% colnames(plot_df)) {
    plot_df[[cln]] <- plot_df[[orig]]
  }
}

# 2. Merged 3CA: EMT + Protein Maturation
emt_prot_sources <- intersect(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III"), retained_3ca_clean)
if (length(emt_prot_sources) > 0) {
  orig_map <- setNames(retained_3ca_clean, retained_3ca)
  orig_sources <- names(orig_map)[orig_map %in% emt_prot_sources]
  mps_sources <- intersect(orig_sources, colnames(plot_df))
  if (length(mps_sources) > 0) {
    plot_df[["3CA_EMT_and_Protein_maturation"]] <- apply(as.matrix(plot_df[, mps_sources, drop = FALSE]), 1, max)
  }
}

final_new_state_names <- unique(c(kept_3ca_names, intersect("3CA_EMT_and_Protein_maturation", colnames(plot_df))))

mp_features <- intersect(names(all_gsva_sets), colnames(plot_df))
state_features <- intersect(c(names(local_state_groups), final_new_state_names), colnames(plot_df))
all_features <- c(mp_features, state_features)

# Define proper ordering
state_level_order <- c(names(local_state_groups), final_new_state_names)
state_levels <- intersect(state_level_order, colnames(plot_df))

mp_ideal_order <- c(unlist(state_groups), setdiff(names(mp.genes), unlist(state_groups)), retained_3ca)
mp_levels_raw <- intersect(mp_ideal_order, colnames(plot_df))

# Map labels for MPs: Description
label_mp <- function(id) {
  if (id %in% names(mp_descriptions)) {
    paste0(id, ": ", mp_descriptions[id])
  } else if (grepl("3CA_mp", id)) {
    clean_3ca_name(id)
  } else {
    id
  }
}
mp_labels <- setNames(sapply(mp_levels_raw, label_mp), mp_levels_raw)

# Correlation of cohort means per feature (EAC mean vs ESCC mean)
mean_by_cohort <- plot_df %>%
  group_by(HistologyGroup) %>%
  summarise(across(all_of(all_features), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(all_features), names_to = "feature", values_to = "mean_score") %>%
  pivot_wider(names_from = HistologyGroup, values_from = mean_score)

# Split into MP and State versions
mean_mp <- mean_by_cohort %>% filter(feature %in% mp_features) %>% mutate(label = mp_labels[feature])
mean_state <- mean_by_cohort %>% filter(feature %in% state_features) %>% mutate(label = feature)

cor_mp <- cor.test(mean_mp$EAC, mean_mp$ESCC, method = "spearman")
cor_state <- cor.test(mean_state$EAC, mean_state$ESCC, method = "spearman")

# Define common limits for correlation plots (Page 1)
common_min <- min(c(mean_mp$EAC, mean_mp$ESCC, mean_state$EAC, mean_state$ESCC), na.rm=TRUE)
common_max <- max(c(mean_mp$EAC, mean_mp$ESCC, mean_state$EAC, mean_state$ESCC), na.rm=TRUE)
padding <- (common_max - common_min) * 0.05
common_lims <- c(common_min - padding, common_max + padding)

p_mean_cor_mp <- ggplot(mean_mp, aes(EAC, ESCC, label = label)) +
  geom_point(color = "darkgrey", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_text_repel(size = 4, max.overlaps = 30, box.padding = 0.5) +
  scale_x_continuous(limits = common_lims) +
  scale_y_continuous(limits = common_lims) +
  labs(
    title = "Clinical MP Correlation: EAC vs ESCC",
    subtitle = paste0("Spearman rho = ", round(unname(cor_mp$estimate), 3), ", p = ", signif(cor_mp$p.value, 3)),
    x = "EAC mean", y = "ESCC mean"
  ) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 14))

p_mean_cor_state <- ggplot(mean_state, aes(EAC, ESCC, label = label)) +
  geom_point(color = "navy", size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  geom_text_repel(size = 5.5, max.overlaps = 30, box.padding = 0.7) +
  scale_x_continuous(limits = common_lims) +
  scale_y_continuous(limits = common_lims) +
  labs(
    title = "Clinical State Correlation: EAC vs ESCC",
    subtitle = paste0("Spearman rho = ", round(unname(cor_state$estimate), 3), ", p = ", signif(cor_state$p.value, 3)),
    x = "EAC mean", y = "ESCC mean"
  ) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 14))

p_page1 <- (p_mean_cor_mp | p_mean_cor_state)

# Removed presence plots as requested

# Comparison Plot (Side-by-side MP and State)

# Long format for MPs
mp_long <- plot_df %>%
  select(sample_barcode, HistologyGroup, any_of(mp_levels_raw)) %>%
  pivot_longer(cols = any_of(mp_levels_raw), names_to = "MP", values_to = "Score") %>%
  mutate(MP_label = factor(mp_labels[MP], levels = mp_labels))

# Compute stats for MPs
mp_stats <- mp_long %>%
  group_by(MP_label) %>%
  do({
    test <- wilcox.test(Score ~ HistologyGroup, data = .)
    data.frame(p_val = test$p.value, max_val = max(.$Score, na.rm = TRUE))
  }) %>%
  mutate(stars = cut(p_val, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                     labels = c("***", "**", "*", "ns"))) %>%
  mutate(y_pos = max_val + (max(max_val) * 0.05))

# Long format for States
state_long <- plot_df %>%
  select(sample_barcode, HistologyGroup, any_of(state_levels)) %>%
  pivot_longer(cols = any_of(state_levels), names_to = "State", values_to = "Score") %>%
  mutate(State_label = factor(paste0("State: ", State), levels = paste0("State: ", state_levels)))

# Compute stats for States
state_stats <- state_long %>%
  group_by(State_label) %>%
  do({
    test <- wilcox.test(Score ~ HistologyGroup, data = .)
    data.frame(p_val = test$p.value, max_val = max(.$Score, na.rm = TRUE))
  }) %>%
  mutate(stars = cut(p_val, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                     labels = c("***", "**", "*", "ns"))) %>%
  mutate(y_pos = max_val + (max(max_val) * 0.05))

# Define common Y limits for boxplots (Page 2)
common_y_min <- min(c(mp_long$Score, state_long$Score), na.rm=TRUE)
common_y_max <- max(c(mp_stats$y_pos, state_stats$y_pos), na.rm=TRUE)
common_y_lims <- c(common_y_min, common_y_max * 1.05)

p_mp <- ggplot(mp_long, aes(x = MP_label, y = Score, fill = HistologyGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7, color = "black", linewidth = 0.7) +
  geom_text(data = mp_stats, aes(x = MP_label, y = y_pos, label = stars), 
            inherit.aes = FALSE, size = 6, vjust = -0.2, fontface = "bold") +
  scale_fill_manual(values = c("EAC" = "#E41A1C", "ESCC" = "#377EB8")) +
  scale_y_continuous(limits = common_y_lims) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        plot.title = element_text(face = "bold")) +
  labs(title = "Clinical Metaprogram Scores", x = NULL, y = "GSVA Score")

p_state <- ggplot(state_long, aes(x = State_label, y = Score, fill = HistologyGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7, color = "black", linewidth = 0.7) +
  geom_text(data = state_stats, aes(x = State_label, y = y_pos, label = stars), 
            inherit.aes = FALSE, size = 6, vjust = -0.2, fontface = "bold") +
  scale_fill_manual(values = c("EAC" = "#E41A1C", "ESCC" = "#377EB8")) +
  scale_y_continuous(limits = common_y_lims) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        plot.title = element_text(face = "bold")) +
  labs(title = "Clinical State Scores", x = NULL, y = "Score")

p_combined <- (p_mp | p_state) + 
  plot_layout(guides = "collect", widths = c(1.8, 1)) & 
  theme(legend.position = "bottom", legend.text = element_text(size = 16))

####################
# Define Page 3: Delta Plots (EAC - ESCC) with Sample-level Stats
# 1. Prepare statistics
temp_mp_stats <- mp_stats %>% ungroup() %>% select(MP_label, stars, p_val)
temp_state_stats <- state_stats %>% ungroup() %>% 
  mutate(feature_name = gsub("^State: ", "", State_label)) %>%
  select(feature_name, stars, p_val)

mean_mp <- mean_mp %>%
  mutate(delta = EAC - ESCC) %>%
  left_join(temp_mp_stats, by = c("label" = "MP_label")) %>%
  mutate(Direction = ifelse(delta > 0, "High in EAC", "High in ESCC"),
         Direction = factor(Direction, levels = c("High in EAC", "High in ESCC")))

mean_state <- mean_state %>%
  mutate(delta = EAC - ESCC) %>%
  left_join(temp_state_stats, by = c("feature" = "feature_name")) %>%
  mutate(Direction = ifelse(delta > 0, "High in EAC", "High in ESCC"),
         Direction = factor(Direction, levels = c("High in EAC", "High in ESCC")))

# Determine limits for delta axes with more padding for labels/stars
common_delta_lims <- c(min(c(mean_mp$delta, mean_state$delta), na.rm=TRUE),
                       max(c(mean_mp$delta, mean_state$delta), na.rm=TRUE))
padding_delta <- (common_delta_lims[2] - common_delta_lims[1]) * 0.25
common_delta_lims <- c(common_delta_lims[1] - padding_delta, common_delta_lims[2] + padding_delta)

p_delta_mp <- ggplot(mean_mp, aes(x = reorder(label, delta), y = delta, fill = Direction)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = stars, hjust = ifelse(delta >= 0, -0.2, 1.2)), 
            size = 5.5, fontface = "bold", vjust = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("High in EAC" = "#E41A1C", "High in ESCC" = "#377EB8")) +
  scale_y_continuous(limits = common_delta_lims) +
  theme_classic(base_size = 18) +
  labs(title = "MP Delta: EAC - ESCC", 
       subtitle = "Stars: Sample-level Wilcoxon test", x = NULL, y = "Delta Mean Score") +
  theme(axis.text = element_text(color = "black", size = 13),
        plot.title = element_text(face = "bold"))

p_delta_state <- ggplot(mean_state, aes(x = reorder(label, delta), y = delta, fill = Direction)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = stars, hjust = ifelse(delta >= 0, -0.2, 1.2)), 
            size = 5.5, fontface = "bold", vjust = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("High in EAC" = "#E41A1C", "High in ESCC" = "#377EB8"), guide = "none") +
  scale_y_continuous(limits = common_delta_lims) +
  theme_classic(base_size = 18) +
  labs(title = "State Delta: EAC - ESCC", 
       subtitle = "Stars: Sample-level Wilcoxon test", x = NULL, y = "Delta Mean Score") +
  theme(axis.text = element_text(color = "black", size = 13),
        plot.title = element_text(face = "bold"))

p_page3 <- (p_delta_mp | p_delta_state) + 
  plot_layout(guides = "collect", widths = c(1.8, 1)) & 
  theme(legend.position = "bottom", legend.title = element_blank())
####################

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_tcga_eac_escc_compare_plots.pdf")), width = 16, height = 9, onefile = TRUE)
print(p_page1)
print(p_combined)
print(p_page3)
dev.off()

write.csv(mean_by_cohort, file.path(out_dir, paste0("Auto_", task_prefix, "_feature_means_eac_escc.csv")), row.names = FALSE)
# presence_summary removed

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  task = task_prefix,
  n_eac = sum(meta_tcga$HistologyGroup == "EAC"),
  n_escc = sum(meta_tcga$HistologyGroup == "ESCC"),
  presence_threshold = 0,
  mean_corr_mp_rho = unname(cor_mp$estimate),
  mean_corr_mp_p = cor_mp$p.value,
  mean_corr_state_rho = unname(cor_state$estimate),
  mean_corr_state_p = cor_state$p.value,
  stringsAsFactors = FALSE
), file.path(summary_dir, paste0("Auto_", task_prefix, "_tcga_eac_escc_compare_summary.csv")), row.names = FALSE)

cat("Saved task8 outputs in:", out_dir, "\n")
