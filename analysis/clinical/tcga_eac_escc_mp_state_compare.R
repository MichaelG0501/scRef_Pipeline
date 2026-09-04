####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/tcga_eac_escc_mp_state_compare.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Scores the final 17 centred refined MPs and five state gene
#     sets in reconstructed TCGA-ESCA bulk expression and compares EAC vs ESCC.
#   Inputs: centred merged refined MP genes and reconstructed TCGA ESCA TPM/metadata.
#   Outputs: ref_outs/task8_tcga_eac_escc_compare/ figures/tables and compact summary.
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(GSVA)
library(patchwork)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

task_prefix <- "task8"
out_dir <- paste0(task_prefix, "_tcga_eac_escc_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Updated TCGA reconstruction paths and histology helper
####################
tcga_recon_dir <- "TCGA/esca_gdc_reconstruction"
tcga_meta_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_meta.rds")
tcga_matrix_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_tpm_matrix.rds")
tcga_mixture_path <- file.path(tcga_recon_dir, "tables", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tcga_compat_meta_path <- "tcga_esca_meta.rds"
tcga_compat_mixture_path <- file.path("cibersortx", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
expression_transform <- "log2(TPM + 1)"

infer_histology <- function(type_vec, detailed_vec = NA_character_) {
  t <- tolower(paste(as.character(type_vec), as.character(detailed_vec)))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}
####################

run_gsva <- function(expr_mat, gene_sets) {
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

mp.genes <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")

state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

grouping_current <- read.csv("Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv", check.names = FALSE)
mp_descriptions <- setNames(grouping_current$description, grouping_current$mp)

# Final centred panel only; historical 3CA relabel features are excluded.
all_gsva_sets <- mp.genes

####################
# Load updated TCGA reconstruction and prepare GSVA matrix
####################
if (file.exists(tcga_meta_path)) {
  meta_tcga <- readRDS(tcga_meta_path)
} else if (file.exists(tcga_compat_meta_path)) {
  meta_tcga <- readRDS(tcga_compat_meta_path)
} else {
  stop("Missing TCGA metadata. Run analysis/TCGA/tcga_esca_reconstruct_data.R first.")
}

if (file.exists(tcga_matrix_path)) {
  tpm_mat <- readRDS(tcga_matrix_path)
} else if (file.exists(tcga_mixture_path)) {
  tpm_df <- data.table::fread(tcga_mixture_path)
  tpm_mat <- as.matrix(tpm_df[, -1, with = FALSE])
  rownames(tpm_mat) <- tpm_df[[1]]
} else if (file.exists(tcga_compat_mixture_path)) {
  tpm_df <- data.table::fread(tcga_compat_mixture_path)
  tpm_mat <- as.matrix(tpm_df[, -1, with = FALSE])
  rownames(tpm_mat) <- tpm_df[[1]]
} else {
  stop("Missing TCGA TPM matrix. Run analysis/TCGA/tcga_esca_reconstruct_data.R first.")
}

if (!"sample_barcode" %in% colnames(meta_tcga)) {
  stop("TCGA metadata lacks sample_barcode.")
}
if (is.null(colnames(tpm_mat))) {
  stop("TCGA TPM matrix lacks sample barcode column names.")
}

detailed_vec <- if ("Cancer_Type_Detailed" %in% colnames(meta_tcga)) meta_tcga$Cancer_Type_Detailed else NA_character_
type_vec <- if ("type" %in% colnames(meta_tcga)) meta_tcga$type else detailed_vec
inferred_histology <- infer_histology(type_vec, detailed_vec)
if ("HistologyGroup" %in% colnames(meta_tcga)) {
  meta_tcga$HistologyGroup <- ifelse(
    is.na(meta_tcga$HistologyGroup) | meta_tcga$HistologyGroup == "Other",
    inferred_histology,
    as.character(meta_tcga$HistologyGroup)
  )
} else {
  meta_tcga$HistologyGroup <- inferred_histology
}

common_samples <- intersect(colnames(tpm_mat), meta_tcga$sample_barcode)
if (length(common_samples) < 20) {
  stop("Too few TCGA samples overlap between metadata and TPM matrix: ", length(common_samples))
}
tpm_mat <- tpm_mat[, common_samples, drop = FALSE]
expr_mat <- log2(tpm_mat + 1)
expr_mat[!is.finite(expr_mat)] <- 0

meta_tcga <- meta_tcga %>% 
  filter(sample_barcode %in% common_samples) %>%
  filter(sample_type_code == "01", HistologyGroup %in% c("EAC", "ESCC")) %>%
  mutate(HistologyGroup = factor(HistologyGroup, levels = c("EAC", "ESCC")))

if (nrow(meta_tcga) == 0) {
  stop("No primary EAC/ESCC TCGA samples remain after filtering.")
}

mp_gs <- run_gsva(expr_mat, all_gsva_sets)
if (is.null(mp_gs)) stop("No valid MP GSVA results")
mp_df <- as.data.frame(t(mp_gs))
mp_df$sample_barcode <- rownames(mp_df)

plot_df <- meta_tcga %>% inner_join(mp_df, by = "sample_barcode")
if (nrow(plot_df) == 0) {
  stop("No TCGA metadata rows joined to GSVA scores.")
}
####################

# State scores are the maximum GSVA score among each current state's MPs.
local_state_groups <- state_groups

for (nm in names(local_state_groups)) {
  mps <- intersect(local_state_groups[[nm]], colnames(plot_df))
  if (length(mps) == 0) next
  plot_df[[nm]] <- apply(as.matrix(plot_df[, mps, drop = FALSE]), 1, max)
}

mp_features <- intersect(names(all_gsva_sets), colnames(plot_df))
state_features <- intersect(names(local_state_groups), colnames(plot_df))
all_features <- c(mp_features, state_features)

# Define proper ordering
state_level_order <- names(local_state_groups)
state_levels <- intersect(state_level_order, colnames(plot_df))

mp_ideal_order <- c(unlist(state_groups), setdiff(names(mp.genes), unlist(state_groups)))
mp_levels_raw <- intersect(mp_ideal_order, colnames(plot_df))

# Map labels for MPs: Description
label_mp <- function(id) {
  if (id %in% names(mp_descriptions)) {
    paste0(id, ": ", mp_descriptions[id])
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

####################
# Classify mean-score quadrants for EAC/ESCC comparison
####################
classify_mean_region <- function(eac, escc) {
  dplyr::case_when(
    eac > 0 & escc < 0 ~ "EAC-specific",
    eac < 0 & escc > 0 ~ "ESCC-high_EAC-low",
    eac > 0 & escc > 0 ~ "High-in-both",
    eac < 0 & escc < 0 ~ "Low-in-both",
    TRUE ~ "Boundary"
  )
}

mean_mp <- mean_mp %>%
  mutate(region = classify_mean_region(EAC, ESCC))
mean_state <- mean_state %>%
  mutate(region = classify_mean_region(EAC, ESCC))

region_summary <- bind_rows(
  mean_mp %>% count(region, name = "n_features") %>% mutate(feature_type = "MP"),
  mean_state %>% count(region, name = "n_features") %>% mutate(feature_type = "State")
) %>%
  select(feature_type, region, n_features) %>%
  arrange(feature_type, region)
####################

cor_mp <- cor.test(mean_mp$EAC, mean_mp$ESCC, method = "spearman")
cor_state <- cor.test(mean_state$EAC, mean_state$ESCC, method = "spearman")

# Define common limits for correlation plots (Page 1)
# Allow X and Y to have different ranges so we don't have empty space
x_min <- min(c(mean_mp$EAC, mean_state$EAC, 0), na.rm=TRUE)
x_max <- max(c(mean_mp$EAC, mean_state$EAC, 0), na.rm=TRUE)
y_min <- min(c(mean_mp$ESCC, mean_state$ESCC, 0), na.rm=TRUE)
y_max <- max(c(mean_mp$ESCC, mean_state$ESCC, 0), na.rm=TRUE)

x_pad <- (x_max - x_min) * 0.10
y_pad <- (y_max - y_min) * 0.10

x_lims <- c(x_min - x_pad, x_max + x_pad)
y_lims <- c(y_min - y_pad, y_max + y_pad)

n_eac <- sum(meta_tcga$HistologyGroup == "EAC")
n_escc <- sum(meta_tcga$HistologyGroup == "ESCC")
x_lab <- paste0("EAC mean (n = ", n_eac, ")")
y_lab <- paste0("ESCC mean (n = ", n_escc, ")")

p_mean_cor_mp <- ggplot(mean_mp, aes(EAC, ESCC, label = label)) +
  annotate("rect", xmin = 0, xmax = x_lims[2], ymin = y_lims[1], ymax = 0,
           fill = NA, color = "#E41A1C", linetype = "dashed", linewidth = 1.5) +
  annotate("rect", xmin = x_lims[1], xmax = 0, ymin = 0, ymax = y_lims[2],
           fill = NA, color = "#377EB8", linetype = "dashed", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 1) +
  geom_point(color = "darkgrey", size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  geom_text_repel(size = 5, fontface = "bold", max.overlaps = Inf, box.padding = 1.5, point.padding = 0.8, force = 5, segment.color = "grey50") +
  scale_x_continuous(limits = x_lims) +
  scale_y_continuous(limits = y_lims) +
  labs(
    title = "Clinical MP Correlation: EAC vs ESCC",
    subtitle = paste0("Spearman rho = ", round(unname(cor_mp$estimate), 3), ", p = ", signif(cor_mp$p.value, 3),
                      "; dotted boxes: EAC-specific and ESCC-high/EAC-low regions"),
    x = x_lab, y = y_lab
  ) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 22),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(color = "black", size = 18))

p_mean_cor_state <- ggplot(mean_state, aes(EAC, ESCC, label = label)) +
  annotate("rect", xmin = 0, xmax = x_lims[2], ymin = y_lims[1], ymax = 0,
           fill = NA, color = "#E41A1C", linetype = "dashed", linewidth = 1.5) +
  annotate("rect", xmin = x_lims[1], xmax = 0, ymin = 0, ymax = y_lims[2],
           fill = NA, color = "#377EB8", linetype = "dashed", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 1) +
  geom_point(color = "navy", size = 4.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  geom_text_repel(size = 5.5, fontface = "bold", max.overlaps = Inf, box.padding = 1.5, point.padding = 0.8, force = 5, segment.color = "grey50") +
  scale_x_continuous(limits = x_lims) +
  scale_y_continuous(limits = y_lims) +
  labs(
    title = "Clinical State Correlation: EAC vs ESCC",
    subtitle = paste0("Spearman rho = ", round(unname(cor_state$estimate), 3), ", p = ", signif(cor_state$p.value, 3),
                      "; dotted boxes: EAC-specific and ESCC-high/EAC-low regions"),
    x = x_lab, y = y_lab
  ) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 22),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(color = "black", size = 18))

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

####################
# Avoid patchwork '&' operator for current patchwork/S7 method dispatch
####################
p_mp <- p_mp + theme(legend.position = "bottom", legend.text = element_text(size = 16))
p_state <- p_state + theme(legend.position = "bottom", legend.text = element_text(size = 16))
p_combined <- (p_mp | p_state) + 
  plot_layout(guides = "collect", widths = c(1.8, 1))
####################

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

p_delta_mp <- p_delta_mp + theme(legend.position = "bottom", legend.title = element_blank())
p_delta_state <- p_delta_state + theme(legend.position = "bottom", legend.title = element_blank())
p_page3 <- (p_delta_mp | p_delta_state) + 
  plot_layout(guides = "collect", widths = c(1.8, 1))
####################

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_tcga_eac_escc_compare_plots.pdf")), width = 16, height = 9, onefile = TRUE)
print(p_page1)
print(p_combined)
print(p_page3)
dev.off()

write.csv(mean_by_cohort, file.path(out_dir, paste0("Auto_", task_prefix, "_feature_means_eac_escc.csv")), row.names = FALSE)
write.csv(
  bind_rows(
    mean_mp %>% mutate(feature_type = "MP"),
    mean_state %>% mutate(feature_type = "State")
  ) %>%
    select(feature_type, feature, label, EAC, ESCC, region, everything()),
  file.path(out_dir, paste0("Auto_", task_prefix, "_feature_mean_regions_eac_escc.csv")),
  row.names = FALSE
)
write.csv(region_summary, file.path(out_dir, paste0("Auto_", task_prefix, "_feature_mean_region_summary.csv")), row.names = FALSE)
# presence_summary removed

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
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
  tcga_expression_transform = expression_transform,
  n_high_in_both = sum(region_summary$region == "High-in-both"),
  n_low_in_both = sum(region_summary$region == "Low-in-both"),
  stringsAsFactors = FALSE
), file.path(summary_dir, paste0("Auto_", task_prefix, "_tcga_eac_escc_compare_summary.csv")), row.names = FALSE)

cat("Saved task8 outputs in:", out_dir, "\n")
