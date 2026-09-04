####################
# Analysis registry:
#   Status: active
#   Script: analysis/TCGA/tcga_gender_state_mp_validation.R
#   Methodology: analysis/methodology/TCGA/tcga_reconstruction_and_gender_validation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#     - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds
#     - ref_outs/meta_full_epi.rds
#     - canonical centred refined states, UCell scores, MP genes, and grouping table
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx
#   Outputs:
#     - ref_outs/TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores_centred17.rds
#     - ref_outs/TCGA/gender_validation/tables/Auto_tcga_gender_*csv
#     - ref_outs/TCGA/gender_validation/figures/Auto_tcga_gender_scRef_concordance.pdf/png
#     - updates/new_updates/summaries/Auto_tcga_gender_scRef_concordance_summary.csv
#   Cache/replot behavior:
#     - GSVA scores are cached and reused unless SCREF_FORCE_REBUILD=TRUE.
#     - Plot/stat regeneration is always performed from cached scores and current scRef inputs.
#   Run:
#     Rscript analysis/TCGA/tcga_gender_state_mp_validation.R
#   Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GSVA)
library(patchwork)
library(readxl)
library(scales)
library(stringr)
library(tidyr)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_helpers.R")

setwd(SCREF_PROJECT_DIR)

####################
# 1) Paths and constants
####################
script_path <- file.path(SCREF_ANALYSIS_DIR, "TCGA", "tcga_gender_state_mp_validation.R")
tcga_recon_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "esca_gdc_reconstruction")
out_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "gender_validation")
tiers <- ensure_output_tiers(out_dir)
summary_dir <- SCREF_SUMMARY_DIR
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

tcga_meta_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_meta.rds")
tcga_matrix_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_tpm_matrix.rds")
tcga_mixture_path <- file.path(tcga_recon_dir, "tables", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
gsva_cache_path <- file.path(tiers[["intermediate"]], "Auto_tcga_gender_gsva_scores_centred17.rds")

state_groups <- SCREF_STATE_GROUPS
core_state_order <- SCREF_PRIMARY_STATE_ORDER
sex_palette <- c("Female" = "#B24745", "Male" = "#2F7FB8")

run_summary <- start_run_summary(
  script = script_path,
  input_files = c(
    tcga_meta_path,
    tcga_matrix_path,
    SCREF_META_FULL_EPI_RDS,
    SCREF_FINAL_STATE_RDS,
    SCREF_MP_UCELL_RDS,
    SCREF_MP_GENES_RDS,
    SCREF_MP_GROUPING_CSV,
    SCREF_CLINICAL_XLSX
  ),
  output_files = c(
    file.path(tiers[["figures"]], "Auto_tcga_gender_scRef_concordance.pdf"),
    file.path(tiers[["tables"]], "Auto_tcga_gender_scRef_concordance.csv"),
    file.path(summary_dir, "Auto_tcga_gender_scRef_concordance_summary.csv")
  ),
  parameters = list(
    tcga_validation_cohort = "EAC primary tumors",
    expression_transform = "log2(TPM + 1)",
    state_method = SCREF_PREFERRED_STATE_METHOD
  )
)

####################
# 2) Helpers
####################
normalise_sex <- function(x) {
  x <- tolower(as.character(x))
  dplyr::case_when(
    x %in% c("female", "f") ~ "Female",
    x %in% c("male", "m") ~ "Male",
    TRUE ~ NA_character_
  )
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }
  diffs <- outer(x, y, "-")
  (sum(diffs > 0) - sum(diffs < 0)) / (length(x) * length(y))
}

feature_test_table <- function(long_df, cohort_label) {
  long_df |>
    filter(!is.na(sex), sex %in% c("Female", "Male"), is.finite(value)) |>
    group_by(feature_type, feature) |>
    group_modify(function(.x, .y) {
      female_vals <- .x$value[.x$sex == "Female"]
      male_vals <- .x$value[.x$sex == "Male"]
      if (length(female_vals) < 2 || length(male_vals) < 2) {
        return(tibble(
          cohort = cohort_label,
          n_female = length(female_vals),
          n_male = length(male_vals),
          median_female = median(female_vals, na.rm = TRUE),
          median_male = median(male_vals, na.rm = TRUE),
          delta_female_minus_male = NA_real_,
          cliffs_delta = NA_real_,
          p_value = NA_real_
        ))
      }
      test <- suppressWarnings(wilcox.test(female_vals, male_vals, exact = FALSE))
      tibble(
        cohort = cohort_label,
        n_female = length(female_vals),
        n_male = length(male_vals),
        median_female = median(female_vals, na.rm = TRUE),
        median_male = median(male_vals, na.rm = TRUE),
        delta_female_minus_male = median_female - median_male,
        cliffs_delta = cliffs_delta(female_vals, male_vals),
        p_value = test$p.value
      )
    }) |>
    ungroup() |>
    group_by(feature_type) |>
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>
    ungroup()
}

make_feature_labels <- function(mp_features, state_features) {
  mp_labels <- label_mps(mp_features)
  names(mp_labels) <- mp_features
  state_labels <- setNames(state_features, state_features)
  c(mp_labels, state_labels)
}

run_gsva_scores <- function(expr_mat, gene_sets) {
  gene_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(expr_mat)))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]
  if (length(gene_sets) == 0) {
    stop("No gene sets retained at >=5 genes for GSVA.")
  }
  suppressMessages(gsva(expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian"))
}

plot_tcga_boxes <- function(tcga_long_eac, tcga_stats_eac, features, feature_labels, title_text) {
  features <- features[features %in% tcga_long_eac$feature]
  if (length(features) == 0) {
    return(NULL)
  }
  plot_df <- tcga_long_eac |>
    filter(feature %in% features) |>
    mutate(
      feature = factor(feature, levels = features),
      sex = factor(sex, levels = c("Female", "Male"))
    )

  legend_counts <- plot_df |>
    distinct(sample_barcode, sex) |>
    count(sex)
  legend_labels <- setNames(
    paste0(legend_counts$sex, " (n=", legend_counts$n, ")"),
    legend_counts$sex
  )

  # Calculate y position for text annotations above boxplots per feature
  y_span <- diff(range(plot_df$value, na.rm = TRUE))
  y_offset <- 0.08 * y_span
  annot_df <- plot_df |>
    group_by(feature) |>
    summarise(y_pos = max(value, na.rm = TRUE) + y_offset, .groups = "drop") |>
    left_join(tcga_stats_eac |> select(feature, p_value) |> distinct(), by = "feature") |>
    mutate(
      label = ifelse(
        is.na(p_value),
        "",
        paste0("p=", scales::pvalue(p_value, accuracy = 0.001))
      )
    )

  ggplot(plot_df, aes(x = feature, y = value, fill = sex, color = sex)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.8,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      alpha = 0.7,
      size = 1.0,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = annot_df,
      aes(x = feature, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = sex_palette, labels = legend_labels, drop = FALSE) +
    scale_color_manual(values = sex_palette, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = feature_labels) +
    labs(
      title = title_text,
      x = NULL,
      y = "TCGA EAC GSVA score",
      fill = "Sex"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(12, 16, 12, 12)
    )
}

####################
# 3) Gene sets for TCGA scoring
####################
mp_genes <- readRDS(SCREF_MP_GENES_RDS)

preferred_mp_order <- c(
  "MP1", "MP5", "MP13+", "MP2+", "MP14", "MP3+", "MP6+", "MP11+",
  "MP9+", "MP10+", "MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP12", "MP15"
)
mp_genes <- mp_genes[intersect(c(preferred_mp_order, setdiff(names(mp_genes), preferred_mp_order)), names(mp_genes))]
mp_sets <- mp_genes

state_sets <- lapply(state_groups, function(mps) {
  unique(c(
    unlist(mp_genes[intersect(mps, names(mp_genes))], use.names = FALSE)
  ))
})
state_sets <- state_sets[sapply(state_sets, length) >= 5]
state_sets <- state_sets[intersect(core_state_order, names(state_sets))]
feature_labels <- make_feature_labels(names(mp_sets), names(state_sets))

####################
# 4) TCGA bulk score table
####################
if (!file.exists(tcga_meta_path)) {
  tcga_meta_path <- file.path(SCREF_REF_OUTS_DIR, "tcga_esca_meta.rds")
}
if (!file.exists(tcga_meta_path)) {
  stop("Missing TCGA metadata. Run analysis/TCGA/tcga_esca_reconstruct_data.R first.")
}
meta_tcga <- readRDS(tcga_meta_path) |>
  mutate(
    Gender = normalise_sex(Gender),
    HistologyGroup = ifelse(is.na(HistologyGroup), "Other", HistologyGroup)
  )

if (file.exists(tcga_matrix_path)) {
  tpm_mat <- readRDS(tcga_matrix_path)
} else if (file.exists(tcga_mixture_path)) {
  tpm_df <- fread(tcga_mixture_path)
  tpm_mat <- as.matrix(tpm_df[, -1, with = FALSE])
  rownames(tpm_mat) <- tpm_df[[1]]
} else {
  stop("Missing TCGA TPM matrix. Run analysis/TCGA/tcga_esca_reconstruct_data.R first.")
}

common_samples <- intersect(colnames(tpm_mat), meta_tcga$sample_barcode)
if (length(common_samples) < 20) {
  stop("Too few TCGA samples overlap between metadata and TPM matrix: ", length(common_samples))
}
meta_tcga <- meta_tcga |> filter(sample_barcode %in% common_samples)
tpm_mat <- tpm_mat[, meta_tcga$sample_barcode, drop = FALSE]
expr_log <- log2(tpm_mat + 1)
expr_log[!is.finite(expr_log)] <- 0

gsva_cache <- load_or_build_cache(
  gsva_cache_path,
  build_fun = function() {
    list(
      mp_scores = t(run_gsva_scores(expr_log, mp_sets)),
      state_scores = t(run_gsva_scores(expr_log, state_sets)),
      mp_sets = mp_sets,
      state_sets = state_sets,
      transform = "log2(TPM + 1)"
    )
  },
  label = "tcga_gender_gsva_scores"
)
run_summary <- add_cache_status(run_summary, gsva_cache$label, gsva_cache$path, gsva_cache$reused)
gsva_obj <- gsva_cache$value

tcga_score_df <- meta_tcga |>
  select(sample_barcode, case_barcode, HistologyGroup, sample_type_code, sample_type, sex = Gender) |>
  bind_cols(
    as.data.frame(gsva_obj$mp_scores[meta_tcga$sample_barcode, names(mp_sets), drop = FALSE]),
    as.data.frame(gsva_obj$state_scores[meta_tcga$sample_barcode, names(state_sets), drop = FALSE])
  )

tcga_long_all <- bind_rows(
  tcga_score_df |>
    select(sample_barcode, case_barcode, HistologyGroup, sample_type_code, sample_type, sex, all_of(names(mp_sets))) |>
    pivot_longer(cols = all_of(names(mp_sets)), names_to = "feature", values_to = "value") |>
    mutate(feature_type = "MP"),
  tcga_score_df |>
    select(sample_barcode, case_barcode, HistologyGroup, sample_type_code, sample_type, sex, all_of(names(state_sets))) |>
    pivot_longer(cols = all_of(names(state_sets)), names_to = "feature", values_to = "value") |>
    mutate(feature_type = "State")
) |>
  mutate(
    feature_label = feature_labels[feature],
    sex = factor(sex, levels = c("Female", "Male"))
  )

tcga_long_eac <- tcga_long_all |>
  filter(sample_type_code == "01", HistologyGroup == "EAC", !is.na(sex))
tcga_long_escc <- tcga_long_all |>
  filter(sample_type_code == "01", HistologyGroup == "ESCC", !is.na(sex))
tcga_long_esca <- tcga_long_all |>
  filter(sample_type_code == "01", HistologyGroup %in% c("EAC", "ESCC"), !is.na(sex))

tcga_stats <- bind_rows(
  feature_test_table(tcga_long_eac, "TCGA_EAC_primary"),
  feature_test_table(tcga_long_escc, "TCGA_ESCC_primary"),
  feature_test_table(tcga_long_esca, "TCGA_ESCA_primary")
) |>
  mutate(feature_label = feature_labels[feature])

####################
# 5) scRef sample-level gender association table
####################
meta_full_epi <- readRDS(SCREF_META_FULL_EPI_RDS)
final_states <- readRDS(SCREF_FINAL_STATE_RDS)
ucell_scores <- readRDS(SCREF_MP_UCELL_RDS)

clinical_sheet <- read_excel(SCREF_CLINICAL_XLSX, sheet = 3, skip = 1) |>
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

cell_ids <- Reduce(intersect, list(rownames(meta_full_epi), rownames(ucell_scores), names(final_states)))
cell_meta <- data.frame(
  cell = cell_ids,
  orig.ident = as.character(meta_full_epi[cell_ids, "orig.ident"]),
  stringsAsFactors = FALSE
) |>
  left_join(clinical_sheet |> select(orig.ident, Gender), by = "orig.ident") |>
  mutate(sex = normalise_sex(Gender))

scref_mp_cols <- intersect(names(mp_sets), colnames(ucell_scores))
scref_mp_df <- as.data.frame(ucell_scores[cell_ids, scref_mp_cols, drop = FALSE])
scref_mp_df$cell <- rownames(scref_mp_df)
scref_mp_long <- scref_mp_df |>
  pivot_longer(cols = all_of(scref_mp_cols), names_to = "feature", values_to = "value") |>
  left_join(cell_meta |> select(cell, orig.ident, sex), by = "cell") |>
  filter(!is.na(sex)) |>
  group_by(orig.ident, sex, feature) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(feature_type = "MP")

scref_state_cell <- cell_meta |>
  mutate(state = as.character(final_states[cell])) |>
  filter(!is.na(sex), !is.na(state), !state %in% c("Unresolved", "Hybrid"))

scref_state_levels <- intersect(names(state_sets), unique(scref_state_cell$state))
scref_state_totals <- scref_state_cell |>
  count(orig.ident, sex, name = "total_cells")
scref_state_long <- scref_state_cell |>
  count(orig.ident, sex, state, name = "n_cells") |>
  right_join(
    scref_state_cell |> distinct(orig.ident, sex) |> tidyr::crossing(state = scref_state_levels),
    by = c("orig.ident", "sex", "state")
  ) |>
  mutate(n_cells = replace_na(n_cells, 0L)) |>
  left_join(scref_state_totals, by = c("orig.ident", "sex")) |>
  transmute(
    orig.ident = orig.ident,
    sex = sex,
    feature = state,
    value = 100 * n_cells / pmax(total_cells, 1),
    feature_type = "State"
  )

scref_long <- bind_rows(scref_mp_long, scref_state_long) |>
  mutate(feature_label = feature_labels[feature])
scref_stats <- feature_test_table(scref_long, "scRef_OAC") |>
  mutate(feature_label = feature_labels[feature])

####################
# 6) Concordance tables and summaries
####################
tcga_stats_eac <- tcga_stats |> filter(cohort == "TCGA_EAC_primary")
concordance <- scref_stats |>
  select(
    feature_type,
    feature,
    feature_label,
    scref_n_female = n_female,
    scref_n_male = n_male,
    scref_median_female = median_female,
    scref_median_male = median_male,
    scref_delta = delta_female_minus_male,
    scref_cliffs_delta = cliffs_delta,
    scref_p_value = p_value,
    scref_p_adj = p_adj
  ) |>
  inner_join(
    tcga_stats_eac |>
      select(
        feature_type,
        feature,
        tcga_n_female = n_female,
        tcga_n_male = n_male,
        tcga_median_female = median_female,
        tcga_median_male = median_male,
        tcga_delta = delta_female_minus_male,
        tcga_cliffs_delta = cliffs_delta,
        tcga_p_value = p_value,
        tcga_p_adj = p_adj
      ),
    by = c("feature_type", "feature")
  ) |>
  mutate(
    direction = case_when(
      is.na(scref_cliffs_delta) | is.na(tcga_cliffs_delta) ~ "Not tested",
      sign(scref_cliffs_delta) == sign(tcga_cliffs_delta) ~ "Concordant",
      TRUE ~ "Discordant"
    ),
    scref_nominal = scref_p_value < 0.05,
    tcga_nominal = tcga_p_value < 0.05,
    scRef_direction = ifelse(scref_cliffs_delta > 0, "Female higher", "Male higher"),
    TCGA_direction = ifelse(tcga_cliffs_delta > 0, "Female higher", "Male higher")
  ) |>
  arrange(feature_type, scref_p_value)

summary_df <- concordance |>
  group_by(feature_type) |>
  summarise(
    n_features = n(),
    n_scref_nominal = sum(scref_nominal, na.rm = TRUE),
    n_tcga_nominal = sum(tcga_nominal, na.rm = TRUE),
    n_concordant_all = sum(direction == "Concordant", na.rm = TRUE),
    n_discordant_all = sum(direction == "Discordant", na.rm = TRUE),
    n_concordant_scref_nominal = sum(direction == "Concordant" & scref_nominal, na.rm = TRUE),
    n_discordant_scref_nominal = sum(direction == "Discordant" & scref_nominal, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(tcga_stats, file.path(tiers[["tables"]], "Auto_tcga_gender_tcga_feature_stats.csv"), row.names = FALSE)
write.csv(scref_stats, file.path(tiers[["tables"]], "Auto_tcga_gender_scref_feature_stats.csv"), row.names = FALSE)
write.csv(concordance, file.path(tiers[["tables"]], "Auto_tcga_gender_scRef_concordance.csv"), row.names = FALSE)
write.csv(summary_df, file.path(tiers[["tables"]], "Auto_tcga_gender_scRef_concordance_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(summary_dir, "Auto_tcga_gender_scRef_concordance_summary.csv"), row.names = FALSE)

####################
# 7) Figures
####################
plot_conc <- concordance |>
  filter(is.finite(scref_cliffs_delta), is.finite(tcga_cliffs_delta)) |>
  filter(direction != "Not tested") |>
  mutate(
    label = ifelse(scref_p_value < 0.10 | tcga_p_value < 0.10 | feature_type == "State", feature_labels[feature], ""),
    direction = factor(direction, levels = c("Concordant", "Discordant"))
  )

n_conc_candidates <- sum(plot_conc$direction == "Concordant" & plot_conc$scref_nominal, na.rm = TRUE)
n_candidates <- sum(plot_conc$scref_nominal, na.rm = TRUE)
plot_lim <- max(abs(c(plot_conc$scref_cliffs_delta, plot_conc$tcga_cliffs_delta)), na.rm = TRUE)
plot_lim <- min(1, max(0.45, ceiling((plot_lim + 0.08) * 10) / 10))

p_concordance <- ggplot(plot_conc, aes(scref_cliffs_delta, tcga_cliffs_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey45") +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey45") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey65") +
  geom_point(aes(fill = direction, shape = feature_type, size = -log10(pmax(scref_p_value, 1e-6))), colour = "black", alpha = 0.9) +
  geom_text_repel(aes(label = label), size = 3.0, max.overlaps = 80, min.segment.length = 0, box.padding = 0.35) +
  scale_fill_manual(values = c("Concordant" = "#1B9E77", "Discordant" = "#D95F02")) +
  scale_shape_manual(values = c("MP" = 21, "State" = 24)) +
  scale_size_continuous(range = c(2.4, 5.2), name = "-log10(scRef p)") +
  coord_cartesian(xlim = c(-plot_lim, plot_lim), ylim = c(-plot_lim, plot_lim)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4))) +
  labs(
    title = "Sex-associated scRef MP/state signals in TCGA EAC",
    x = "scRef OAC sample-level effect size",
    y = "TCGA EAC bulk effect size",
    fill = "Direction",
    shape = "Feature type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(colour = "black"),
    legend.position = "right"
  )

candidate_mps <- concordance |>
  filter(feature_type == "MP") |>
  arrange(match(feature, preferred_mp_order)) |>
  pull(feature)
candidate_states <- concordance |>
  filter(feature_type == "State") |>
  arrange(match(feature, core_state_order)) |>
  pull(feature)

p_mp_boxes <- plot_tcga_boxes(
  tcga_long_eac,
  tcga_stats_eac,
  candidate_mps,
  feature_labels,
  "TCGA EAC MP scores by sex"
)
p_state_boxes <- plot_tcga_boxes(
  tcga_long_eac,
  tcga_stats_eac,
  candidate_states,
  feature_labels,
  "TCGA EAC state scores by sex"
)

effect_bar_df <- concordance |>
  mutate(
    feature_label = factor(feature_labels[feature], levels = rev(feature_labels[c(intersect(preferred_mp_order, feature), intersect(core_state_order, feature))])),
    TCGA_direction = factor(TCGA_direction, levels = c("Female higher", "Male higher"))
  )

p_effect_bar <- ggplot(effect_bar_df, aes(feature_label, tcga_cliffs_delta, fill = TCGA_direction)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25, alpha = 0.9) +
  geom_point(aes(y = scref_cliffs_delta, shape = "scRef OAC"), fill = "white", colour = "black", size = 2.5, stroke = 0.7) +
  coord_flip() +
  facet_grid(feature_type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Female higher" = sex_palette[["Female"]], "Male higher" = sex_palette[["Male"]]), drop = FALSE) +
  scale_shape_manual(name = "", values = c("scRef OAC" = 21)) +
  labs(
    title = "TCGA effect direction with scRef effect overlay",
    x = NULL,
    y = "Cliff's delta: Female vs Male",
    fill = "TCGA direction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    axis.text = element_text(colour = "black"),
    strip.text.y = element_text(face = "bold")
  )

pdf_path <- file.path(tiers[["figures"]], "Auto_tcga_gender_scRef_concordance.pdf")
png_path <- file.path(tiers[["figures"]], "Auto_tcga_gender_scRef_concordance.png")
grDevices::cairo_pdf(pdf_path, width = 13, height = 8.5, onefile = TRUE)
print(p_concordance)
if (!is.null(p_mp_boxes)) print(p_mp_boxes)
if (!is.null(p_state_boxes)) print(p_state_boxes)
print(p_effect_bar)
dev.off()

ggsave(png_path, p_concordance, width = 13, height = 8.5, dpi = 300)

write_run_summary(
  run_summary,
  file.path(tiers[["logs"]], "Auto_tcga_gender_scRef_concordance_run_summary.rds")
)

message("TCGA sex concordance validation complete.")
message("Concordance table: ", file.path(tiers[["tables"]], "Auto_tcga_gender_scRef_concordance.csv"))
message("Figure: ", pdf_path)
