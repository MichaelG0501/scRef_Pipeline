####################
# Auto_mp18_cna_investigation.R
# 
# Deep dive into the significant association between MP18 (Secretory Diff. Intest.)
# and pairwise CNA distance across subclones.
#
# Key finding from Auto_v2_pairwise_cna_distance_all_features.pdf:
#   MP18 is the most significant MP associated with CNA distance
#   (sample-centered rho=0.42, p=7.9e-06 for cna_abs_mean)
#
# Questions addressed:
#   1. What is the direction of the correlation?
#   2. Which chromosome arms/CNA events drive this?
#   3. Which samples contribute most?
#   4. Is there a specific subclone CNA profile linked to MP18?
#   5. Visualisations of results
#
# Inputs:
#   ref_outs/Auto_cna_subclone_expression/rds/Auto_cna_subclone_expression_results.rds
#   ref_outs/Auto_cna_subclone_expression/rds/Auto_v2_visualisation_results.rds
#   ref_outs/Auto_cna_subclone_expression/rds/Auto_subclone_arm_cna_long.rds
#   ref_outs/Auto_cna_subclone_expression/rds/Auto_subclone_feature_summary.rds
#
# Outputs:
#   ref_outs/Auto_mp18_cna_investigation/
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
})

options(stringsAsFactors = FALSE)
set.seed(42)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- "Auto_mp18_cna_investigation"
fig_dir <- file.path(out_dir, "figures")
tbl_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ─── Load cached data ────────────────────────────────────────────────
message("Loading cached data")
res <- readRDS("Auto_cna_subclone_expression/rds/Auto_cna_subclone_expression_results.rds")
vis <- readRDS("Auto_cna_subclone_expression/rds/Auto_v2_visualisation_results.rds")
arm_long <- readRDS("Auto_cna_subclone_expression/rds/Auto_subclone_arm_cna_long.rds")
features <- readRDS("Auto_cna_subclone_expression/rds/Auto_subclone_feature_summary.rds")

pairwise_df <- as.data.frame(res$pairwise)
arm_matrix <- as.matrix(res$arm_matrix)

mp_descriptions <- c(
  "MP1" = "G2M Cell Cycle", "MP7" = "DNA Damage Repair", "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation", "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.", "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.", "MP8" = "Intestinal Diff.",
  "MP18" = "Secretory Diff. (Intest.)", "MP16" = "Secretory Diff. (Gastric)",
  "MP13" = "Hypoxic Inflam. Epi.", "MP12" = "Neuro-responsive Epi",
  "MP15" = "Immune Infiltration"
)

mp_order <- c("MP1","MP7","MP9","MP2","MP17","MP14","MP5","MP10",
              "MP8","MP18","MP16","MP13","MP12","MP15")

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Unresolved" = "grey80", "Hybrid" = "black"
)

infer_study <- function(x) sub("^([^_]+_[0-9]{4}).*$", "\\1", as.character(x))

# ─── 1. Overall scatter: CNA distance vs MP18 delta ──────────────────
message("1. Creating CNA distance vs MP18 delta scatter plot")

pw <- pairwise_df %>%
  mutate(study = infer_study(sample))

# Confirm correlation stats
mp18_col <- "abs_delta__mp__MP18"
if (!mp18_col %in% colnames(pw)) stop("MP18 pairwise column missing")

# Global correlation
global_cor <- cor.test(pw$cna_abs_mean, pw[[mp18_col]], method = "spearman")
message(sprintf("  Global Spearman rho=%.3f, p=%.2e", global_cor$estimate, global_cor$p.value))

# Sample-centered correlation
pw_centered <- pw %>%
  group_by(sample) %>%
  mutate(
    cna_c = cna_abs_mean - mean(cna_abs_mean, na.rm = TRUE),
    mp18_c = .data[[mp18_col]] - mean(.data[[mp18_col]], na.rm = TRUE)
  ) %>%
  ungroup()

centered_cor <- cor.test(pw_centered$cna_c, pw_centered$mp18_c, method = "spearman")
message(sprintf("  Sample-centered Spearman rho=%.3f, p=%.2e", centered_cor$estimate, centered_cor$p.value))

# Scatter: raw
p1_raw <- ggplot(pw, aes(cna_abs_mean, .data[[mp18_col]])) +
  geom_point(aes(color = study), alpha = 0.65, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(
    title = "Pairwise CNA distance vs MP18 expression difference",
    subtitle = sprintf("Spearman ρ = %.3f, p = %.2e (n = %d subclone pairs)",
                       global_cor$estimate, global_cor$p.value, nrow(pw)),
    x = "Mean absolute CNA distance\n(arm-level, between subclone pair)",
    y = "Absolute delta MP18 score\n(Secretory Diff. Intestinal)",
    color = "Study"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# Scatter: sample-centered
p1_centered <- ggplot(pw_centered, aes(cna_c, mp18_c)) +
  geom_point(aes(color = study), alpha = 0.65, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(
    title = "Sample-centered: CNA distance vs MP18 divergence",
    subtitle = sprintf("Sample-centered Spearman ρ = %.3f, p = %.2e",
                       centered_cor$estimate, centered_cor$p.value),
    x = "Sample-centered CNA distance",
    y = "Sample-centered MP18 delta",
    color = "Study"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# ─── 2. Compare MP18 to all other MPs ────────────────────────────────
message("2. Comparing all MP correlations with CNA distance")

mp_delta_cols <- grep("^abs_delta__mp__MP", colnames(pw), value = TRUE)

all_mp_cors <- bind_rows(lapply(mp_delta_cols, function(col) {
  mp_name <- sub("^abs_delta__mp__", "", col)
  cr <- tryCatch(cor.test(pw$cna_abs_mean, pw[[col]], method = "spearman"),
                 error = function(e) NULL)
  # also sample-centered
  pw_c <- pw %>%
    group_by(sample) %>%
    mutate(
      cna_c = cna_abs_mean - mean(cna_abs_mean, na.rm = TRUE),
      feat_c = .data[[col]] - mean(.data[[col]], na.rm = TRUE)
    ) %>%
    ungroup()
  cr_c <- tryCatch(cor.test(pw_c$cna_c, pw_c$feat_c, method = "spearman"),
                   error = function(e) NULL)
  data.frame(
    mp = mp_name,
    mp_label = paste0(mp_name, " ", mp_descriptions[mp_name]),
    global_rho = if (!is.null(cr)) unname(cr$estimate) else NA_real_,
    global_p = if (!is.null(cr)) cr$p.value else NA_real_,
    centered_rho = if (!is.null(cr_c)) unname(cr_c$estimate) else NA_real_,
    centered_p = if (!is.null(cr_c)) cr_c$p.value else NA_real_,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    global_p_adj = p.adjust(global_p, method = "BH"),
    centered_p_adj = p.adjust(centered_p, method = "BH"),
    sig_global = ifelse(global_p_adj < 0.05, "*", ""),
    sig_centered = ifelse(centered_p_adj < 0.05, "*", ""),
    is_mp18 = mp == "MP18"
  )

write.csv(all_mp_cors, file.path(tbl_dir, "Auto_all_mp_cna_distance_correlations.csv"), row.names = FALSE)

# Bar chart of sample-centered rho by MP
all_mp_cors <- all_mp_cors %>%
  mutate(mp_label = factor(mp_label, levels = mp_label[order(centered_rho)]))

p2_bar <- ggplot(all_mp_cors, aes(centered_rho, mp_label, fill = is_mp18)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.72) +
  geom_text(aes(label = ifelse(centered_p_adj < 0.001, "***",
                                ifelse(centered_p_adj < 0.01, "**",
                                       ifelse(centered_p_adj < 0.05, "*", "")))),
            hjust = ifelse(all_mp_cors$centered_rho >= 0, -0.3, 1.3),
            size = 6, fontface = "bold") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("TRUE" = "#FF7F00", "FALSE" = "grey70"), guide = "none") +
  labs(
    title = "Sample-centered Spearman ρ: CNA distance vs MP expression divergence",
    subtitle = "Stars indicate FDR < 0.05 after BH correction. MP18 highlighted in orange.",
    x = "Sample-centered Spearman ρ",
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 13, color = "grey30"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  )

# ─── 3. Per-arm CNA vs MP18: which arms drive the association? ──────
message("3. Testing which chromosome arms are associated with MP18 expression")

arm_levels <- colnames(arm_matrix)

# For each arm, test correlation between arm CNA value and MP18 score across subclones
features <- features %>%
  mutate(subclone_id = paste(sample, subclone, sep = "::"))

arm_mp18_tests <- bind_rows(lapply(arm_levels, function(arm_name) {
  if (!arm_name %in% colnames(arm_matrix)) return(NULL)
  ids <- intersect(rownames(arm_matrix), features$subclone_id)
  if (length(ids) < 5) return(NULL)
  
  arm_val <- arm_matrix[ids, arm_name]
  mp18_val <- features$mp__MP18[match(ids, features$subclone_id)]
  
  ok <- is.finite(arm_val) & is.finite(mp18_val)
  if (sum(ok) < 5) return(NULL)
  
  cr <- tryCatch(cor.test(arm_val[ok], mp18_val[ok], method = "spearman"),
                 error = function(e) NULL)
  if (is.null(cr)) return(NULL)
  
  data.frame(
    arm = arm_name,
    chr = sub("[pq]$", "", arm_name),
    rho = unname(cr$estimate),
    p_value = cr$p.value,
    n = sum(ok),
    mean_arm_cna = mean(arm_val[ok], na.rm = TRUE),
    sd_arm_cna = sd(arm_val[ok], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = ifelse(p_adj < 0.001, "***",
                 ifelse(p_adj < 0.01, "**",
                        ifelse(p_adj < 0.05, "*", ""))),
    direction_label = ifelse(rho > 0, "positive", "negative")
  ) %>%
  arrange(p_adj)

write.csv(arm_mp18_tests, file.path(tbl_dir, "Auto_arm_vs_mp18_correlations.csv"), row.names = FALSE)

message("  Top arms associated with MP18:")
print(head(arm_mp18_tests %>% select(arm, rho, p_value, p_adj, sig), 15))

# Genome-wide Manhattan-style plot
arm_mp18_tests <- arm_mp18_tests %>%
  mutate(arm = factor(arm, levels = arm_levels))

p3_manhattan <- ggplot(arm_mp18_tests, aes(arm, -log10(p_value))) +
  geom_col(aes(fill = rho), color = "black", linewidth = 0.25, width = 0.72) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text(aes(label = sig), vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman ρ\n(arm CNA vs MP18)") +
  labs(
    title = "Chromosome arm CNA vs MP18 score (across all subclones)",
    subtitle = "Red dashed line = p = 0.05 (uncorrected). Stars = BH-adjusted FDR < 0.05.",
    x = "Chromosome arm",
    y = "-log10(p-value)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# ─── 4. Per-sample MP18 vs CNA distance ─────────────────────────────
message("4. Per-sample breakdown of MP18 vs CNA distance")

# Which samples have multiple subclones with informative MP18 variation?
sample_mp18_stats <- features %>%
  filter(!is.na(mp__MP18)) %>%
  group_by(sample) %>%
  summarise(
    n_subclones = n(),
    mp18_mean = mean(mp__MP18, na.rm = TRUE),
    mp18_sd = sd(mp__MP18, na.rm = TRUE),
    mp18_range = max(mp__MP18, na.rm = TRUE) - min(mp__MP18, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_subclones >= 2) %>%
  arrange(desc(mp18_range))

write.csv(sample_mp18_stats, file.path(tbl_dir, "Auto_sample_mp18_subclone_variation.csv"), row.names = FALSE)

# Per-sample pairwise plot
pw_per_sample <- pw %>%
  filter(sample %in% sample_mp18_stats$sample) %>%
  mutate(sample_short = sub("^[^_]+_[0-9]{4}_", "", sample))

# Keep top samples with most variation
top_samples <- head(sample_mp18_stats$sample, 12)

pw_top <- pw_per_sample %>%
  filter(sample %in% top_samples) %>%
  mutate(sample = factor(sample, levels = top_samples))

if (nrow(pw_top) > 0) {
  p4_persample <- ggplot(pw_top, aes(cna_abs_mean, .data[[mp18_col]])) +
    geom_point(color = "#FF7F00", alpha = 0.7, size = 3.5) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
    facet_wrap(~sample, scales = "free", ncol = 4) +
    labs(
      title = "Per-sample: CNA distance vs MP18 divergence (top 12 by MP18 range)",
      x = "Mean absolute CNA distance",
      y = "Absolute MP18 score delta"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 9)
    )
}

# ─── 5. Subclone-level: high vs low MP18 subclones CNA profile ──────
message("5. Comparing CNA profiles of high vs low MP18 subclones")

# Classify subclones as high/low MP18 within their sample
features_classified <- features %>%
  filter(!is.na(mp__MP18)) %>%
  group_by(sample) %>%
  mutate(
    mp18_sample_median = median(mp__MP18, na.rm = TRUE),
    mp18_group = ifelse(mp__MP18 > mp18_sample_median, "High MP18", "Low MP18"),
    mp18_rank = rank(mp__MP18, ties.method = "first"),
    n_sub = n()
  ) %>%
  ungroup() %>%
  filter(n_sub >= 2)  # need at least 2 subclones per sample

ids_classified <- intersect(features_classified$subclone_id, rownames(arm_matrix))
features_classified <- features_classified %>% filter(subclone_id %in% ids_classified)

# For each arm, compare CNA values between high and low MP18 subclones
arm_high_low_tests <- bind_rows(lapply(arm_levels, function(arm_name) {
  d <- features_classified %>%
    mutate(arm_cna = arm_matrix[subclone_id, arm_name])
  
  high_vals <- d$arm_cna[d$mp18_group == "High MP18"]
  low_vals <- d$arm_cna[d$mp18_group == "Low MP18"]
  
  if (length(high_vals) < 3 || length(low_vals) < 3) return(NULL)
  
  wt <- tryCatch(wilcox.test(high_vals, low_vals), error = function(e) NULL)
  if (is.null(wt)) return(NULL)
  
  data.frame(
    arm = arm_name,
    mean_high = mean(high_vals, na.rm = TRUE),
    mean_low = mean(low_vals, na.rm = TRUE),
    delta = mean(high_vals, na.rm = TRUE) - mean(low_vals, na.rm = TRUE),
    wilcox_p = wt$p.value,
    n_high = length(high_vals),
    n_low = length(low_vals),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    p_adj = p.adjust(wilcox_p, method = "BH"),
    sig = ifelse(p_adj < 0.001, "***",
                 ifelse(p_adj < 0.01, "**",
                        ifelse(p_adj < 0.05, "*", "")))
  ) %>%
  arrange(p_adj)

write.csv(arm_high_low_tests, file.path(tbl_dir, "Auto_high_vs_low_mp18_arm_cna.csv"), row.names = FALSE)

message("  Top arms with CNA difference between high vs low MP18 subclones:")
print(head(arm_high_low_tests %>% select(arm, delta, wilcox_p, p_adj, sig), 15))

# Genome-wide delta plot
arm_high_low_tests <- arm_high_low_tests %>%
  mutate(arm = factor(arm, levels = arm_levels))

p5_delta <- ggplot(arm_high_low_tests, aes(arm, delta)) +
  geom_col(aes(fill = delta), color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(aes(label = sig), vjust = ifelse(arm_high_low_tests$delta >= 0, -0.3, 1.3),
            size = 4.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Mean CNA delta\n(High − Low MP18)") +
  labs(
    title = "CNA difference: High MP18 vs Low MP18 subclones",
    x = "Chromosome arm",
    y = "Mean CNA delta (High − Low MP18)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8)
  )

p5b_significance <- ggplot(arm_high_low_tests, aes(arm, ifelse(delta > 0, -log10(wilcox_p), log10(wilcox_p)))) +
  geom_col(aes(fill = delta > 0), color = "black", linewidth = 0.25, width = 0.72) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_text(aes(label = sig), vjust = ifelse(arm_high_low_tests$delta >= 0, -0.3, 1.3),
            size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
  labs(
    title = "Significance of CNA differences: High vs Low MP18",
    x = "Chromosome arm",
    y = "Signed -log10(p-value)\n(Up = Gain in High MP18, Down = Loss)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, face = "bold"),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold")
  )

# ─── 5b. TCGA Validation of Top Arms ────────────────────────────────
message("5b. TCGA Validation of Top Arms")

tcga_arm <- readRDS("TCGA/cna_recurrent_event_validation/intermediate/Auto_tcga_weighted_arm_cna_calls.rds")
tcga_gsva <- readRDS("TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores.rds")

tcga_mp18 <- data.frame(
  sample_barcode = rownames(tcga_gsva$mp_scores),
  sample_key = substr(rownames(tcga_gsva$mp_scores), 1, 15),
  mp18_score = tcga_gsva$mp_scores[, "MP18"],
  stringsAsFactors = FALSE
)

tcga_data <- tcga_arm %>%
  inner_join(tcga_mp18, by = "sample_key") %>%
  mutate(
    cna_call = case_when(
      arm_mean >= 0.1 ~ "Gain",
      arm_mean <= -0.1 ~ "Loss",
      TRUE ~ "Neutral"
    )
  )

top4_arms <- head(arm_high_low_tests$arm, 4)

sc_validation_data <- arm_long %>%
  filter(arm %in% top4_arms) %>%
  left_join(features %>% select(subclone_id, mp18_score = mp__MP18), by = "subclone_id") %>%
  filter(!is.na(mp18_score)) %>%
  mutate(
    cna_call = case_when(
      cna_call == 1 ~ "Gain",
      cna_call == -1 ~ "Loss",
      TRUE ~ "Neutral"
    ),
    Dataset = "scATLAS subclones"
  )

tcga_validation_data <- tcga_data %>%
  filter(arm %in% top4_arms) %>%
  mutate(Dataset = "TCGA")

combined_validation <- bind_rows(
  sc_validation_data %>% select(arm, cna_call, mp18_score, Dataset),
  tcga_validation_data %>% select(arm, cna_call, mp18_score, Dataset)
) %>%
  mutate(
    cna_call = factor(cna_call, levels = c("Loss", "Neutral", "Gain")),
    Dataset = factor(Dataset, levels = c("scATLAS subclones", "TCGA")),
    arm = factor(arm, levels = top4_arms)
  )

p_vals <- combined_validation %>%
  group_by(Dataset, arm) %>%
  summarise(
    p_val = tryCatch(kruskal.test(mp18_score ~ cna_call)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = ifelse(is.na(p_val), "NA", sprintf("p = %.2e", p_val))
  ) %>%
  left_join(
    combined_validation %>%
      group_by(Dataset) %>%
      summarise(max_y = max(mp18_score, na.rm=TRUE) + 0.15 * diff(range(mp18_score, na.rm=TRUE)), .groups="drop"),
    by = "Dataset"
  )

p10_tcga_validation <- ggplot(combined_validation, aes(x = cna_call, y = mp18_score, fill = cna_call)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
  geom_text(data = p_vals, aes(x = 2, y = max_y, label = p_label), 
            inherit.aes = FALSE, size = 6, fontface = "italic") +
  facet_grid(Dataset ~ arm, scales = "free_y") +
  scale_fill_manual(values = c("Loss" = "#2166AC", "Neutral" = "grey80", "Gain" = "#B2182B")) +
  labs(
    title = "Validation of Top 4 CNA Arms Associated with MP18",
    x = "Chromosome Arm CNA Status",
    y = "MP18 Expression Score",
    fill = "CNA Call"
  ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    strip.text = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

# ─── 6. Pairwise arm-specific delta vs MP18 delta ────────────────────
message("6. Testing per-arm CNA delta vs MP18 delta at pairwise level")

# For each subclone pair, compute per-arm CNA delta and see which arm deltas
# correlate most with MP18 delta
pair_arm_cors <- bind_rows(lapply(arm_levels, function(arm_name) {
  id1s <- pw$subclone_id_1
  id2s <- pw$subclone_id_2
  
  # arm delta for each pair
  ok1 <- id1s %in% rownames(arm_matrix)
  ok2 <- id2s %in% rownames(arm_matrix)
  ok <- ok1 & ok2
  
  if (sum(ok) < 5) return(NULL)
  
  arm_delta <- abs(arm_matrix[id1s[ok], arm_name] - arm_matrix[id2s[ok], arm_name])
  mp18_delta <- pw[[mp18_col]][ok]
  
  cr <- tryCatch(cor.test(arm_delta, mp18_delta, method = "spearman"),
                 error = function(e) NULL)
  if (is.null(cr)) return(NULL)
  
  data.frame(
    arm = arm_name,
    rho = unname(cr$estimate),
    p_value = cr$p.value,
    n = sum(ok),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = ifelse(p_adj < 0.001, "***",
                 ifelse(p_adj < 0.01, "**",
                        ifelse(p_adj < 0.05, "*", ""))),
    arm = factor(arm, levels = arm_levels)
  ) %>%
  arrange(p_adj)

write.csv(pair_arm_cors, file.path(tbl_dir, "Auto_pairwise_arm_delta_vs_mp18_delta.csv"), row.names = FALSE)

message("  Top arms whose pairwise CNA delta correlates with MP18 delta:")
print(head(pair_arm_cors %>% select(arm, rho, p_value, p_adj, sig), 15))

p6_pairarm <- ggplot(pair_arm_cors, aes(arm, rho)) +
  geom_col(aes(fill = rho), color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(aes(label = sig), vjust = ifelse(pair_arm_cors$rho >= 0, -0.3, 1.3),
            size = 4.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman ρ") +
  labs(
    title = "Pairwise arm CNA delta vs MP18 delta (all subclone pairs)",
    subtitle = "Which arm-level CNA differences drive MP18 expression differences? Stars = BH FDR < 0.05.",
    x = "Chromosome arm",
    y = "Spearman ρ (arm CNA Δ vs MP18 Δ)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8)
  )

# ─── 7. Top driving arms: focused scatter plots ─────────────────────
message("7. Creating focused scatter plots for top driving arms")

top_driving_arms <- pair_arm_cors %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  head(6)

if (nrow(top_driving_arms) > 0) {
  scatter_data_list <- lapply(seq_len(nrow(top_driving_arms)), function(i) {
    arm_name <- as.character(top_driving_arms$arm[i])
    id1s <- pw$subclone_id_1
    id2s <- pw$subclone_id_2
    ok <- id1s %in% rownames(arm_matrix) & id2s %in% rownames(arm_matrix)
    
    data.frame(
      arm = arm_name,
      arm_delta = abs(arm_matrix[id1s[ok], arm_name] - arm_matrix[id2s[ok], arm_name]),
      mp18_delta = pw[[mp18_col]][ok],
      sample = pw$sample[ok],
      rho = top_driving_arms$rho[i],
      p_adj = top_driving_arms$p_adj[i],
      stringsAsFactors = FALSE
    )
  })
  
  scatter_data <- bind_rows(scatter_data_list) %>%
    mutate(arm_label = paste0(arm, " (ρ=", round(rho, 3), ", FDR=", signif(p_adj, 2), ")"))
  
  p7_topscatter <- ggplot(scatter_data, aes(arm_delta, mp18_delta)) +
    geom_point(alpha = 0.5, size = 2.5, color = "#FF7F00") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
    facet_wrap(~arm_label, scales = "free_x", ncol = 3) +
    labs(
      title = "Top chromosome arms driving MP18-CNA distance association",
      subtitle = "Arm-specific CNA delta vs MP18 expression delta between subclone pairs",
      x = "Absolute arm-level CNA delta",
      y = "Absolute MP18 score delta"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 12, color = "grey30"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

# ─── 8. Subclone CNA burden vs MP18 ─────────────────────────────────
message("8. CNA burden vs MP18 at subclone level")

features_plot <- features %>%
  filter(!is.na(mp__MP18), !is.na(cna_burden_mean_abs))

burden_cor <- cor.test(features_plot$cna_burden_mean_abs, features_plot$mp__MP18, method = "spearman")
message(sprintf("  CNA burden vs MP18: rho=%.3f, p=%.2e", burden_cor$estimate, burden_cor$p.value))

# Determine dominant state for each subclone
state_cols_in_data <- grep("^state__", colnames(features_plot), value = TRUE)
if (length(state_cols_in_data) > 0) {
  features_plot$dominant_state <- apply(features_plot[, state_cols_in_data], 1, function(x) {
    state_names <- sub("^state__", "", state_cols_in_data)
    state_names <- gsub("_", " ", state_names)
    state_names[which.max(x)]
  })
} else {
  features_plot$dominant_state <- "Unknown"
}

p8_burden <- ggplot(features_plot, aes(cna_burden_mean_abs, mp__MP18)) +
  geom_point(aes(color = dominant_state, size = n_cells), alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  scale_color_manual(values = c(
    "Classic Proliferative" = "#E41A1C",
    "Basal to Intestinal Metaplasia" = "#4DAF4A",
    "SMG like Metaplasia" = "#FF7F00",
    "Stress adaptive" = "#984EA3",
    "Immune Infiltrating" = "#377EB8",
    "3CA EMT and Protein maturation" = "#666666",
    "Unresolved" = "grey80", "Hybrid" = "black"
  )) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    title = "CNA burden vs MP18 score per subclone",
    subtitle = sprintf("Spearman ρ = %.3f, p = %.2e", burden_cor$estimate, burden_cor$p.value),
    x = "CNA burden (mean absolute arm-level CNA)",
    y = "Mean MP18 score (Secretory Diff. Intestinal)",
    color = "Dominant state",
    size = "Cells"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 11)
  )

# ─── 9. CNA cluster-level MP18 comparison ───────────────────────────
message("9. CNA cluster-level MP18 comparison")

cna_cluster <- res$cna_cluster
cluster_mp18 <- features %>%
  filter(!is.na(mp__MP18)) %>%
  mutate(cna_cluster = cna_cluster[subclone_id]) %>%
  filter(!is.na(cna_cluster))

if (n_distinct(cluster_mp18$cna_cluster) >= 2) {
  kw <- kruskal.test(mp__MP18 ~ cna_cluster, data = cluster_mp18)
  message(sprintf("  Kruskal-Wallis MP18 by CNA cluster: p=%.4e", kw$p.value))
  
  # Pairwise Wilcoxon
  pw_test <- pairwise.wilcox.test(cluster_mp18$mp__MP18, cluster_mp18$cna_cluster, 
                                   p.adjust.method = "BH")
  
  cluster_n <- cluster_mp18 %>%
    group_by(cna_cluster) %>%
    summarise(n = n(), mp18_median = median(mp__MP18, na.rm = TRUE), .groups = "drop")
  
  p9_cluster <- ggplot(cluster_mp18, aes(cna_cluster, mp__MP18, fill = cna_cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6, linewidth = 0.7) +
    geom_jitter(aes(size = n_cells), alpha = 0.45, width = 0.15, color = "black") +
    geom_text(data = cluster_n, aes(x = cna_cluster, y = max(cluster_mp18$mp__MP18, na.rm = TRUE) + 0.02,
                                     label = paste0("n=", n)),
              inherit.aes = FALSE, size = 5.5, fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    scale_size_continuous(range = c(1.5, 5.5)) +
    labs(
      title = "MP18 expression by CNA consensus cluster",
      subtitle = sprintf("Kruskal-Wallis p = %.2e (k = %d clusters)", kw$p.value, n_distinct(cluster_mp18$cna_cluster)),
      x = NULL,
      y = "Mean MP18 score (Secretory Diff. Intestinal)",
      fill = "CNA cluster",
      size = "Cells"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 14, color = "grey30"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )
}

# ─── 10. Combined summary figure ─────────────────────────────────────
message("10. Saving all figures")

pdf(file.path(fig_dir, "Auto_mp18_cna_investigation.pdf"),
    width = 18, height = 12, useDingbats = FALSE)

# Page 1: Pairwise scatter (raw + centered)
print(p1_raw)
print(p1_centered)

# Page 2: All MPs comparison bar chart
print(p2_bar)

# Page 3: Genome-wide arm CNA vs MP18 correlation (Manhattan-style)
print(p3_manhattan)

# Page 4: High vs Low MP18 CNA profile
print(p5_delta)
print(p5b_significance)
print(p10_tcga_validation)

# Page 5: Pairwise arm CNA delta vs MP18 delta
print(p6_pairarm)

# Page 6: Top driving arm scatters
if (exists("p7_topscatter")) print(p7_topscatter)

# Page 7: CNA burden vs MP18
print(p8_burden)

# Page 8: Per-sample breakdown
if (exists("p4_persample") && nrow(pw_top) > 0) print(p4_persample)

# Page 9: CNA cluster comparison
if (exists("p9_cluster")) print(p9_cluster)

dev.off()

# ─── Summary statistics ──────────────────────────────────────────────
message("\n=== SUMMARY ===")
message(sprintf("MP18 is POSITIVELY correlated with CNA distance (rho=%.3f)", global_cor$estimate))
message("Interpretation: Subclone pairs with MORE divergent CNA profiles tend to")
message("have LARGER differences in MP18 (Secretory Diff. Intestinal) expression.")
message("")
message("Top 5 arms driving the association (pairwise arm CNA Δ vs MP18 Δ):")
top5 <- pair_arm_cors %>% arrange(p_adj) %>% head(5)
for (i in seq_len(nrow(top5))) {
  message(sprintf("  %s: rho=%.3f, FDR=%.2e %s", top5$arm[i], top5$rho[i], top5$p_adj[i], top5$sig[i]))
}
message("")
message("Top 5 arms with CNA difference between High vs Low MP18 subclones:")
top5b <- arm_high_low_tests %>% arrange(p_adj) %>% head(5)
for (i in seq_len(nrow(top5b))) {
  message(sprintf("  %s: CNA delta=%.4f (high-low), FDR=%.2e %s", 
                  top5b$arm[i], top5b$delta[i], top5b$p_adj[i], top5b$sig[i]))
}
message("")
message(sprintf("CNA burden vs MP18 subclone-level: rho=%.3f, p=%.2e", 
                burden_cor$estimate, burden_cor$p.value))

message("\nDone. All outputs in: ref_outs/Auto_mp18_cna_investigation/")
