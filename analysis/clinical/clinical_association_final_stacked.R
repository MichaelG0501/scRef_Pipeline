####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/clinical_association_final_stacked.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_clinical_assoc_stacked_final.R
# Clean final-state stacked clinical association workflow for scRef.
#
# Consolidates the useful parts of the older stacked-bar scripts:
#   - Uses Auto_final_states.rds when available.
#   - Keeps full clinical-variable coverage and per-study faceting.
#   - Writes outputs to the current updates/new_updates/summaries path.
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(readxl)
library(stringr)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

meta_full_epi <- readRDS("meta_full_epi.rds")
state_path <- "Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds"
state_vec <- readRDS(state_path)

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
  sheet = 3,
  skip = 1
) %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

cell_ids <- intersect(names(state_vec), rownames(meta_full_epi))
study_vec <- if ("study" %in% colnames(meta_full_epi)) {
  as.character(meta_full_epi[cell_ids, "study"])
} else {
  sub("^(([^_]+)_[^_]+).*$", "\\1", as.character(meta_full_epi[cell_ids, "orig.ident"]))
}

cell_df <- data.frame(
  cell = cell_ids,
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

core_state_levels <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia",
  "Stress adaptive",
  "Cancer-cell immune mimicry"
)
preferred_extra_states <- c()
trailing_state_levels <- c("Unresolved", "Hybrid")
present_states <- unique(as.character(cell_df$state))
other_extra_states <- setdiff(present_states, c(core_state_levels, preferred_extra_states, trailing_state_levels))
state_levels <- c(
  core_state_levels[core_state_levels %in% present_states],
  preferred_extra_states[preferred_extra_states %in% present_states],
  sort(other_extra_states),
  trailing_state_levels[trailing_state_levels %in% present_states]
)

state_colors <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

missing_color_states <- setdiff(state_levels, names(state_colors))
if (length(missing_color_states) > 0) {
  state_colors <- c(state_colors, setNames(scales::hue_pal()(length(missing_color_states)), missing_color_states))
}

plot_configs <- list(
  list(var = "Gender", title = "Gender", filter = NULL),
  list(var = "Age_Group", title = "Age (>60)", filter = NULL),
  list(var = "Race/Ethnicity", title = "Race", filter = NULL),
  list(var = "Tumor Location", title = "Tumor Location", filter = NULL),
  list(var = "Tumor Type", title = "Tumor Type", filter = NULL),
  list(var = "Grade (Differentiation)", title = "Grade", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (All)", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (Tx-Naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Stage AJCC", title = "Stage AJCC (Post-Tx)", filter = "Treatment == 'Post'"),
  list(var = "Technology", title = "Sequencing Technology", filter = NULL),
  list(var = "Treatment", title = "Treatment Timing", filter = NULL),
  list(var = "Clinical response", title = "Clinical Response (Predictive, Tx-Naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Clinical response", title = "Clinical Response (Post-Tx)", filter = "Treatment == 'Post'")
)

apply_optional_filter <- function(data, filter_expr = NULL) {
  if (is.null(filter_expr)) return(data)
  data %>% filter(!!rlang::parse_expr(filter_expr))
}

compute_plot_data <- function(data, group_var, filter_expr = NULL, facet_var = NULL) {
  data <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(state), !is.na(study))
  if (nrow(data) == 0) return(list(plot_data = NULL, label_data = NULL))

  sample_group_cols <- unique(c(facet_var, "orig.ident", group_var, "state"))
  sample_total_cols <- unique(c(facet_var, "orig.ident", group_var))
  mean_group_cols <- unique(c(facet_var, group_var, "state"))
  total_group_cols <- unique(c(facet_var, group_var))

  sample_prop <- data %>%
    group_by(across(all_of(sample_group_cols))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(sample_total_cols))) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  plot_data <- sample_prop %>%
    group_by(across(all_of(mean_group_cols))) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(total_group_cols))) %>%
    mutate(group_total = sum(mean_prop), pct = ifelse(group_total > 0, 100 * mean_prop / group_total, 0)) %>%
    ungroup() %>%
    mutate(state = factor(state, levels = state_levels))

  label_data <- data %>%
    distinct(across(all_of(unique(c(facet_var, "orig.ident", group_var, "study"))))) %>%
    group_by(across(all_of(total_group_cols))) %>%
    summarise(
      n_samples = n_distinct(orig.ident),
      n_studies = n_distinct(study),
      studies_list = paste(sort(unique(study)), collapse = "\n"),
      .groups = "drop"
    ) %>%
    left_join(
      data %>%
        group_by(across(all_of(total_group_cols))) %>%
        summarise(n_cells = n(), .groups = "drop"),
      by = total_group_cols
    )

  if (is.null(facet_var)) {
    total_samples <- n_distinct(data$orig.ident)
    total_studies <- n_distinct(data$study)
    label_data <- label_data %>%
      mutate(label = paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples, "\n(", n_studies, "/", total_studies, ")\n\n", studies_list))
  } else {
    facet_stats <- data %>%
      group_by(across(all_of(facet_var))) %>%
      summarise(facet_cells = n(), facet_samples = n_distinct(orig.ident), .groups = "drop") %>%
      mutate(facet_label = paste0(.data[[facet_var]], "\n(N=", comma(facet_cells), "; samples=", facet_samples, ")"))
    plot_data <- plot_data %>% left_join(facet_stats[, c(facet_var, "facet_label")], by = facet_var)
    label_data <- label_data %>%
      left_join(facet_stats, by = facet_var) %>%
      mutate(label = paste0(n_samples, "/", facet_samples), facet_label = paste0(.data[[facet_var]], "\n(N=", comma(facet_cells), "; samples=", facet_samples, ")"))
  }

  list(plot_data = plot_data, label_data = label_data)
}

plot_stacked_assoc <- function(data, group_var, title, filter_expr = NULL, facet_var = NULL) {
  prep <- compute_plot_data(data, group_var, filter_expr, facet_var)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

  x_levels <- unique(as.character(prep$plot_data[[group_var]]))
  detail_lines <- if (is.null(facet_var)) max(stringr::str_count(prep$label_data$label, "\\n") + 1) else 2
  expansion_mult <- if (is.null(facet_var)) 0.05 + 0.035 * detail_lines else 0.16

  p <- ggplot(prep$plot_data, aes(x = .data[[group_var]], y = pct / 100, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.72) +
    geom_text(
      data = prep$plot_data %>% mutate(label_text = ifelse(pct >= 8, sprintf("%.0f%%", pct), "")),
      aes(label = label_text, group = state),
      position = position_stack(vjust = 0.5),
      size = 3.8,
      fontface = "bold",
      colour = "black",
      show.legend = FALSE
    ) +
    geom_text(
      data = prep$label_data,
      aes(x = .data[[group_var]], y = 1.02, label = label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4,
      lineheight = 1.0,
      fontface = "bold"
    ) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, expansion_mult))) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"), width = 16)) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = ifelse(length(x_levels) > 3, 45, 25), hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 15),
      plot.margin = margin(10, 20, 10, 10)
    )

  if (!is.null(facet_var)) {
    p <- p +
      facet_wrap(~ facet_label, scales = "free_x", ncol = 4) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 10)
      )
  }

  p
}

safe_plot <- function(cfg, facet_var = NULL) {
  tryCatch(
    plot_stacked_assoc(cell_df, cfg$var, cfg$title, cfg$filter, facet_var),
    error = function(e) {
      message("Skipping ", cfg$title, ": ", e$message)
      NULL
    }
  )
}

combined_plots <- lapply(plot_configs, safe_plot)
names(combined_plots) <- vapply(plot_configs, `[[`, character(1), "title")
combined_plots <- combined_plots[!vapply(combined_plots, is.null, logical(1))]

pdf("Auto_clinical_assoc_stacked_final_combined.pdf", width = 17, height = 12, useDingbats = FALSE)
for (plot_name in names(combined_plots)) {
  print(combined_plots[[plot_name]])
}
dev.off()

pdf("Auto_clinical_assoc_stacked_final_per_study.pdf", width = 17, height = 11, useDingbats = FALSE)
for (cfg in plot_configs) {
  p <- safe_plot(cfg, facet_var = "study")
  if (!is.null(p)) print(p)
}
dev.off()

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- bind_rows(lapply(plot_configs, function(cfg) {
  prep <- compute_plot_data(cell_df, cfg$var, cfg$filter)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  prep$plot_data %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      source_state_file = state_path,
      clinical_variable = cfg$var,
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
      level = .data[[cfg$var]],
      state = as.character(state),
      mean_sample_prop_pct = pct,
      samples_in_level = n_samples,
      cells_in_level = n_cells
    )
}))

write.csv(
  summary_df,
  file.path(summary_dir, "Auto_clinical_assoc_stacked_final_summary.csv"),
  row.names = FALSE
)

message("Saved scRef final stacked clinical association outputs.")
