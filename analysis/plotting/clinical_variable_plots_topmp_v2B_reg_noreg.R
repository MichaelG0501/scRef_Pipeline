####################
# Auto_clinical_variable_plots_topmp_v2B_reg_noreg.R
# Unified clinical association plotting for Approach B reg and noreg states.
# Keeps the full variable coverage from clinical_variable_plots.R and writes
# paired pages (reg then noreg) into the same PDF output.
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(readxl)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")

meta_full_epi <- readRDS("meta_full_epi.rds")
state_map <- list(
  reg = readRDS("Auto_topmp_v2_reg_states_B.rds"),
  noreg = readRDS("Auto_topmp_v2_noreg_states_B.rds")
)

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
  sheet = 3,
  skip = 1
) %>% mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

build_cell_df <- function(state_vec, mode_name) {
  cell_ids <- intersect(names(state_vec), rownames(meta_full_epi))
  study_vec <- if ("study" %in% colnames(meta_full_epi)) {
    as.character(meta_full_epi[cell_ids, "study"])
  } else {
    sub("^(([^_]+)_[^_]+).*$", "\\1", as.character(meta_full_epi[cell_ids, "orig.ident"]))
  }

  data.frame(
    cell = cell_ids,
    mode = mode_name,
    state = as.character(state_vec[cell_ids]),
    orig.ident = as.character(meta_full_epi[cell_ids, "orig.ident"]),
    study = study_vec,
    stringsAsFactors = FALSE
  ) %>%
    left_join(clinical_sheet, by = "orig.ident") %>%
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
}

cell_df <- bind_rows(
  build_cell_df(state_map$reg, "reg"),
  build_cell_df(state_map$noreg, "noreg")
)

state_levels <- c(
  "Classic_Proliferative", "Barretts_Metaplasia", "EMT_related",
  "Intestinal_Metaplasia", "Immune_Infiltrated", "Unresolved", "Hybrid"
)
state_colors <- c(
  Classic_Proliferative = "#E41A1C",
  Barretts_Metaplasia = "#4DAF4A",
  EMT_related = "#984EA3",
  Intestinal_Metaplasia = "#FF7F00",
  Immune_Infiltrated = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

compute_plot_data <- function(data, group_var, filter_expr = NULL, studyname = NULL) {
  if (!is.null(studyname)) {
    data <- data %>% filter(study == studyname)
  }
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }
  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(state), !is.na(study))
  if (nrow(data) == 0) return(list(plot_data = NULL, label_data = NULL, total_samples = 0L))

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

  n_sample_group <- data %>%
    distinct(orig.ident, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_samples = n_distinct(orig.ident), .groups = "drop")

  n_cell_group <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_cells = n(), .groups = "drop")

  label_data <- n_sample_group %>%
    left_join(n_cell_group, by = group_var) %>%
    mutate(label = paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples))

  list(plot_data = plot_data, label_data = label_data, total_samples = total_samples)
}

plot_clinical_assoc <- function(data, group_var, title, filter_expr = NULL) {
  prep <- compute_plot_data(data, group_var, filter_expr)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

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
    geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.2, color = "black") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
}

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

pdf("Auto_clinical_assoc_topmp_v2B_reg_noreg_combined.pdf", width = 14, height = 16)
for (mode_name in requested_modes) {
  mode_df <- cell_df %>% filter(mode == mode_name)

  p1 <- plot_clinical_assoc(mode_df, "Gender", paste0("[", mode_name, "] Gender"))
  p2 <- plot_clinical_assoc(mode_df, "Age_Group", paste0("[", mode_name, "] Age (>60)"))
  p3 <- plot_clinical_assoc(mode_df, "Race/Ethnicity", paste0("[", mode_name, "] Race"))
  p4 <- plot_clinical_assoc(mode_df, "Tumor Location", paste0("[", mode_name, "] Location"))
  p5 <- plot_clinical_assoc(mode_df, "Tumor Type", paste0("[", mode_name, "] Type"))
  p6 <- plot_clinical_assoc(mode_df, "Grade (Differentiation)", paste0("[", mode_name, "] Grade"))
  p7 <- plot_clinical_assoc(mode_df, "Stage AJCC", paste0("[", mode_name, "] Stage (All)"))
  p8 <- plot_clinical_assoc(mode_df, "Stage AJCC", paste0("[", mode_name, "] Stage (Tx-Naive)"), "Treatment == 'Tx-naïve'")
  p9 <- plot_clinical_assoc(mode_df, "Stage AJCC", paste0("[", mode_name, "] Stage (Post-Tx)"), "Treatment == 'Post'")
  p10 <- plot_clinical_assoc(mode_df, "Technology", paste0("[", mode_name, "] Technology"))
  p11 <- plot_clinical_assoc(mode_df, "Treatment", paste0("[", mode_name, "] Tx Timing"))
  p12 <- plot_clinical_assoc(mode_df, "Clinical response", paste0("[", mode_name, "] Response (Predictive)"), "Treatment == 'Tx-naïve'")
  p13 <- plot_clinical_assoc(mode_df, "Clinical response", paste0("[", mode_name, "] Response (Post-Tx)"), "Treatment == 'Post'")

  print(((p1 + p2) / p3) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = paste0("[", mode_name, "] Demographics Overview")))
  print(((p4 + p5) / p6) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = paste0("[", mode_name, "] Tumor Characteristics")))
  print((p7 / p8 / p9) + plot_annotation(title = paste0("[", mode_name, "] Staging Analysis")))
  print((p10 / p11) + plot_layout(heights = c(1.5, 1)) + plot_annotation(title = paste0("[", mode_name, "] Technology & Treatment")))
  print((p12 + p13) + plot_annotation(title = paste0("[", mode_name, "] Clinical Response")))
}
dev.off()

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- lapply(requested_modes, function(mode_name) {
  mode_df <- cell_df %>% filter(mode == mode_name)
  bind_rows(lapply(plot_configs, function(cfg) {
    prep <- compute_plot_data(mode_df, cfg$var, cfg$filter)
    if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
    prep$plot_data %>%
      transmute(
        mode = mode_name,
        clinical_variable = cfg$var,
        level = .data[[cfg$var]],
        state = as.character(state),
        mean_sample_prop_pct = pct
      )
  }))
})

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(summary_dir, "Auto_clinical_assoc_topmp_v2B_reg_noreg_summary.csv"), row.names = FALSE)

message("Saved unified reg+noreg clinical association outputs.")
