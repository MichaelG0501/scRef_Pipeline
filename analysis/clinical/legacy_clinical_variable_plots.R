####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_clinical_variable_plots.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_clinical_variable_plots_topmp_v2B.R
# Clinical association plots using saved Top-MP v2 Approach B states
# (no tmdata object required).
#
# Input:
#   ref_outs/meta_full_epi.rds
#   ref_outs/Auto_topmp_v2_states_B.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx (sheet 3)
#
# Output:
#   ref_outs/Auto_clinical_assoc_topmp_v2B_combined.pdf
#   ref_outs/Auto_clinical_assoc_topmp_v2B_per_study.pdf
#   new_update/summaries/Auto_clinical_assoc_topmp_v2B_summary.csv
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(readxl)

####################
# 1) Setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# 2) Load metadata + Approach B states
####################
meta_full_epi <- readRDS("meta_full_epi.rds")
state_B <- readRDS("Auto_topmp_v2_states_B.rds")

####################
# 3) Load sample-level clinical sheet and build sample key
####################
clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
  sheet = 3,
  skip = 1
) %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

####################
# 4) Build cell-level table from saved states and metadata
####################
cell_ids <- intersect(names(state_B), rownames(meta_full_epi))

####################
# `meta_full_epi.rds` does not always contain a dedicated `study` column.
# Derive study from `orig.ident` as Author_Year when missing.
####################
study_vec <- if ("study" %in% colnames(meta_full_epi)) {
  as.character(meta_full_epi[cell_ids, "study"])
} else {
  orig_tmp <- as.character(meta_full_epi[cell_ids, "orig.ident"])
  sub("^(([^_]+)_[^_]+).*$", "\\1", orig_tmp)
}

cell_df <- data.frame(
  cell = cell_ids,
  state = as.character(state_B[cell_ids]),
  orig.ident = as.character(meta_full_epi[cell_ids, "orig.ident"]),
  study = study_vec,
  stringsAsFactors = FALSE
)

####################
# 5) Attach sample-level clinical covariates to each cell via orig.ident
####################
cell_df <- cell_df %>%
  left_join(clinical_sheet, by = "orig.ident")

####################
# 6) Clinical cleaning
####################
cell_df <- cell_df %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age)),
    Age_Group = ifelse(!is.na(Age) & Age > 60, ">60", "<=60"),
    Gender = ifelse(Gender %in% c("F", "female", "Female"), "Female", "Male"),
    Treatment = factor(Treatment, levels = c("Tx-naïve", "Post")),
    `Clinical response` = factor(
      ifelse(`Clinical response` == "R", "Responder", "Nonresponder"),
      levels = c("Responder", "Nonresponder")
    )
  )

####################
# 7) State order + colours (Approach B)
####################
state_levels <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "f, Immune Infiltrating",
  "Unresolved",
  "Hybrid"
)

state_colors <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

####################
# 8) Per-sample proportions -> group mean -> rescale to 100%
####################
compute_plot_data <- function(data, group_var, filter_expr = NULL, studyname = NULL) {
  if (!is.null(studyname)) {
    data <- data %>% filter(study == studyname)
  }
  
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }
  
  data <- data %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(state), !is.na(study))
  
  if (nrow(data) == 0) {
    return(list(plot_data = NULL, label_data = NULL, total_samples = 0L))
  }
  
  sample_prop <- data %>%
    group_by(orig.ident, .data[[group_var]], state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident, .data[[group_var]]) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()
  
  mean_prop <- sample_prop %>%
    group_by(.data[[group_var]], state) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")
  
  group_sum <- mean_prop %>%
    group_by(.data[[group_var]]) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")
  
  plot_data <- mean_prop %>%
    left_join(group_sum, by = group_var) %>%
    mutate(
      pct = ifelse(group_total > 0, 100 * mean_prop / group_total, 0),
      state = factor(state, levels = state_levels)
    )
  
  total_samples <- data %>% distinct(orig.ident) %>% nrow()
  total_studies <- data %>% distinct(study) %>% nrow()
  
  n_sample_group <- data %>%
    distinct(orig.ident, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_samples = n_distinct(orig.ident), .groups = "drop")
  
  n_cell_group <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  study_group <- data %>%
    distinct(orig.ident, study, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n_studies = n_distinct(study),
      studies_list = paste(sort(unique(study)), collapse = "\n"),
      .groups = "drop"
    )
  
  label_data <- n_sample_group %>%
    left_join(n_cell_group, by = group_var) %>%
    left_join(study_group, by = group_var) %>%
    mutate(
      label = if (!is.null(studyname)) {
        paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples)
      } else {
        paste0(
          "N=", comma(n_cells),
          "\n", n_samples, "/", total_samples,
          "\n(", n_studies, "/", total_studies, ")",
          "\n\n", studies_list
        )
      }
    )
  
  list(plot_data = plot_data, label_data = label_data, total_samples = total_samples)
}

####################
# 9) Combined-style plot helper
####################
plot_clinical_assoc <- function(data, group_var, title, filter_expr = NULL, studyname = NULL) {
  prep <- compute_plot_data(data, group_var, filter_expr, studyname)
  
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) {
    return(NULL)
  }
  
  show_details <- is.null(studyname)
  max_lines <- max(stringr::str_count(prep$label_data$label, "\\n") + 1)
  expansion_mult <- if (show_details) {
    0.04 + (max_lines * 0.04)
  } else {
    0.15
  }
  
  ggplot(prep$plot_data, aes(x = .data[[group_var]], y = pct / 100, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    geom_text(
      data = prep$label_data,
      aes(x = .data[[group_var]], y = 1.02, label = label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3,
      lineheight = 0.9,
      color = "black"
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, expansion_mult))
    ) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")) +
    theme_minimal(base_size = 16) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 18)
    )
}

####################
# 10) Per-study faceted helper (style like previous)
####################
plot_variable_facet <- function(data, group_var, title, filter_expr = NULL) {
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }
  
  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(study), !is.na(state))
  if (nrow(data) == 0) return(NULL)
  
  sample_prop <- data %>%
    group_by(study, orig.ident, .data[[group_var]], state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(study, orig.ident, .data[[group_var]]) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()
  
  mean_prop <- sample_prop %>%
    group_by(study, .data[[group_var]], state) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")
  
  group_sum <- mean_prop %>%
    group_by(study, .data[[group_var]]) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")
  
  plot_data <- mean_prop %>%
    left_join(group_sum, by = c("study", group_var)) %>%
    mutate(
      freq = ifelse(group_total > 0, mean_prop / group_total, 0),
      state = factor(state, levels = state_levels)
    )
  
  study_stats <- data %>%
    group_by(study) %>%
    summarise(
      n_cells = n(),
      n_samples = n_distinct(orig.ident),
      .groups = "drop"
    ) %>%
    mutate(study_label = paste0(study, "\n(N=", comma(n_cells), "; samples=", n_samples, ")"))
  
  label_data <- data %>%
    distinct(study, orig.ident, .data[[group_var]]) %>%
    group_by(study, .data[[group_var]]) %>%
    summarise(n_samples_level = n_distinct(orig.ident), .groups = "drop") %>%
    left_join(study_stats, by = "study") %>%
    mutate(
      label = paste0(n_samples_level, "/", n_samples),
      study_label = paste0(study, "\n(N=", comma(n_cells), "; samples=", n_samples, ")")
    )
  
  plot_data <- plot_data %>%
    left_join(study_stats[, c("study", "study_label")], by = "study")
  
  ggplot(plot_data, aes(x = .data[[group_var]], y = freq, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    geom_text(
      data = label_data,
      aes(x = .data[[group_var]], y = 1.02, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    facet_wrap(~ study_label, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")) +
    theme_bw(base_size = 14) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

####################
# 11) Plot configuration
####################
plot_configs <- list(
  list(var = "Gender", title = "Gender Distribution by Study", filter = NULL),
  list(var = "Age_Group", title = "Age (>60) Distribution by Study", filter = NULL),
  list(var = "Race/Ethnicity", title = "Race Distribution by Study", filter = NULL),
  list(var = "Tumor Location", title = "Tumor Location by Study", filter = NULL),
  list(var = "Tumor Type", title = "Tumor Type by Study", filter = NULL),
  list(var = "Grade (Differentiation)", title = "Grade by Study", filter = NULL),
  list(var = "Stage AJCC", title = "Stage (All Samples) by Study", filter = NULL),
  list(var = "Stage AJCC", title = "Stage (Tx-Naive Only) by Study", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Stage AJCC", title = "Stage (Post-Tx Only) by Study", filter = "Treatment == 'Post'"),
  list(var = "Technology", title = "Technology Used by Study", filter = NULL),
  list(var = "Treatment", title = "Treatment Timing by Study", filter = NULL),
  list(var = "Clinical response", title = "Response (Predictive / Tx-Naive) by Study", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Clinical response", title = "Response (Assessment / Post-Tx) by Study", filter = "Treatment == 'Post'")
)

####################
# 12) Combined cohort report (old multi-panel layout)
####################
p1 <- plot_clinical_assoc(cell_df, "Gender", "Gender")
p2 <- plot_clinical_assoc(cell_df, "Age_Group", "Age (>60)")
p3 <- plot_clinical_assoc(cell_df, "Race/Ethnicity", "Race")

p4 <- plot_clinical_assoc(cell_df, "Tumor Location", "Location")
p5 <- plot_clinical_assoc(cell_df, "Tumor Type", "Type")
p6 <- plot_clinical_assoc(cell_df, "Grade (Differentiation)", "Grade")

p7 <- plot_clinical_assoc(cell_df, "Stage AJCC", "Stage (All)")
p8 <- plot_clinical_assoc(cell_df, "Stage AJCC", "Stage (Tx-Naive)", "Treatment == 'Tx-naïve'")
p9 <- plot_clinical_assoc(cell_df, "Stage AJCC", "Stage (Post-Tx)", "Treatment == 'Post'")

p10 <- plot_clinical_assoc(cell_df, "Technology", "Technology")
p11 <- plot_clinical_assoc(cell_df, "Treatment", "Tx Timing")

p12 <- plot_clinical_assoc(cell_df, "Clinical response", "Response (Predictive)", "Treatment == 'Tx-naïve'")
p13 <- plot_clinical_assoc(cell_df, "Clinical response", "Response (Post-Tx)", "Treatment == 'Post'")

pdf("Auto_clinical_assoc_topmp_v2B_combined.pdf", width = 14, height = 16)
print(((p1 + p2) / p3) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = "Demographics Overview - All"))
print(((p4 + p5) / p6) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = "Tumor Characteristics - All"))
print((p7 / p8 / p9) + plot_annotation(title = "Staging Analysis - All"))
print((p10 / p11) + plot_layout(heights = c(1.5, 1)) + plot_annotation(title = "Technology & Treatment - All"))
print((p12 + p13) + plot_annotation(title = "Clinical Response - All"))
dev.off()

####################
# 13) Per-study report (old faceted style)
####################
pdf("Auto_clinical_assoc_topmp_v2B_per_study.pdf", width = 16, height = 10)
for (cfg in plot_configs) {
  p <- plot_variable_facet(
    data = cell_df,
    group_var = cfg$var,
    filter_expr = cfg$filter,
    title = cfg$title
  )
  if (!is.null(p)) print(p)
}
dev.off()

####################
# 14) Save machine-readable summary
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "new_update",
  "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- lapply(plot_configs, function(cfg) {
  prep <- compute_plot_data(
    data = cell_df,
    group_var = cfg$var,
    filter_expr = cfg$filter,
    studyname = NULL
  )
  
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) {
    return(NULL)
  }
  
  prep$plot_data %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      clinical_variable = cfg$var,
      level = .data[[cfg$var]],
      state = as.character(state),
      mean_sample_prop_pct = pct,
      samples_in_level = n_samples,
      total_samples = prep$total_samples
    )
})

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(summary_dir, "Auto_clinical_assoc_topmp_v2B_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

message("Saved: Auto_clinical_assoc_topmp_v2B_combined.pdf")
message("Saved: Auto_clinical_assoc_topmp_v2B_per_study.pdf")
message(sprintf("Saved: %s", summary_path))