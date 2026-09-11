####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/Auto_occams_bulk_clinical_associations.R
#   Description: Tests current centred refined OCCAMS bulk MP/state GSVA scores
#     across source-column-specific clinical strata and plots subject-level boxplots.
#   Methodology: analysis/methodology/clinical/Auto_occams_bulk_clinical_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
# Inputs:
#   - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_mp_scores.rds
#   - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_state_scores.rds
#   - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
# Outputs:
#   intermediate/:
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_clinical_association_plot_data.rds
#   tables/:
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_variable_inventory.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_stats.csv
#   figures/:
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_clinical_score_boxplots_MP.pdf
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_clinical_score_boxplots_State.pdf
#   logs/:
#     - ref_outs/OCCAMS/clinical/logs/Auto_occams_bulk_clinical_associations.stdout.log
#     - ref_outs/OCCAMS/clinical/logs/Auto_occams_bulk_clinical_associations.stderr.log
#   reports/ and summaries/:
#     - ref_outs/OCCAMS/clinical/reports/Auto_OCCAMS_clinical_association_run_summary.txt
#     - updates/new_updates/summaries/Auto_OCCAMS_clinical_association_summary.csv
# Cache/replot behavior:
#   SCREF_FORCE_REBUILD=TRUE rebuilds persistent plot/statistics data.
#   SCREF_REPLOT_ONLY=TRUE redraws both PDFs from the persistent plot-data RDS.
# Run command:
#   /opt/pbs/bin/qsub analysis/clinical/Auto_occams_bulk_clinical_associations.sh
# Conda env: dmtcp
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

out_root <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "clinical")
out_dirs <- setNames(file.path(out_root, c("intermediate", "tables", "figures", "logs", "reports")),
                     c("intermediate", "tables", "figures", "logs", "reports"))
invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(SCREF_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

mp_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_mp_scores.rds")
state_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_state_scores.rds")
metadata_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_qc_pass_subject_metadata.csv")
grouping_path <- SCREF_MP_GROUPING_CSV
plot_data_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_clinical_association_plot_data.rds")
inventory_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_variable_inventory.csv")
stats_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_stats.csv")
boxplot_path_mp <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_clinical_score_boxplots_MP.pdf")
boxplot_path_state <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_clinical_score_boxplots_State.pdf")
report_path <- file.path(out_dirs[["reports"]], "Auto_OCCAMS_clinical_association_run_summary.txt")
summary_path <- file.path(SCREF_SUMMARY_DIR, "Auto_OCCAMS_clinical_association_summary.csv")

force_rebuild <- identical(toupper(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")), "TRUE")
replot_only <- identical(toupper(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")), "TRUE")
minimum_level_n <- 10L

variable_config <- data.frame(
  variable = c(
    "gender", "age_group", "diagnostic_combined_tumour_site", "pretreatment_combined_tumour_site",
    "pathology_location", "diagnostic_histology", "diagnostic_siewert_classification",
    "pretreatment_siewert_classification", "pathology_siewert_classification",
    "pretreatment_t_stage", "pretreatment_n_stage_tnm7", "pretreatment_m_stage",
    "treatment_plan_performance_status_group", "treatment_plan_intent",
    "treatment_plan_curative_modality", "surgery_main", "treatment_plan_comorbidity_group",
    "pathology_t_stage", "pathology_n_stage_tnm7", "pathology_m_stage",
    "pathology_positive_node_group", "pathology_grade_group", "pathology_signet_ring_cells",
    "pathology_barrett_adjacent_dysplasia", "pathology_neoadjuvant_therapy",
    "surgery_procedure_group", "pathology_mandard_score", "pathology_mandard_response_group",
    "treatment_response_group", "recurrence_ever_recorded"
  ),
  label = c(
    "Gender", "Age at diagnosis (>60)", "Diagnostic combined tumour site",
    "Final pretreatment combined tumour site", "Resection pathology location", "Diagnostic histology",
    "Diagnostic Siewert classification", "Final pretreatment Siewert classification",
    "Resection pathology Siewert classification", "Final pretreatment T stage",
    "Final pretreatment N stage (TNM7)", "Final pretreatment M stage", "Treatment-plan performance status",
    "Treatment-plan intent", "Treatment-plan curative modality", "Main surgery",
    "Treatment-plan comorbidity burden", "Resection pathology T stage",
    "Resection pathology N stage (TNM7)", "Resection pathology M stage",
    "Resection positive-node burden", "Resection pathology grade", "Resection signet-ring cells",
    "Resection adjacent Barrett-associated dysplasia", "Recorded neoadjuvant therapy at resection",
    "Operation type", "Mandard tumour regression grade", "Mandard response grouping",
    "Recorded treatment response grouping", "Disease recurrence ever recorded"
  ),
  category = c(
    rep("Demographic", 2), "Diagnostic", "Pretreatment", "Resection pathology",
    "Diagnostic", "Diagnostic", "Pretreatment", "Resection pathology",
    rep("Pretreatment", 3), rep("Treatment plan", 3), "Surgery", "Treatment plan",
    rep("Resection pathology", 8), "Surgery", rep("Resection pathology", 2),
    "Treatment response", "Longitudinal follow-up"
  ),
  rationale = c(
    "Demographic variable", "Demographic/prognostic variable", "RD anatomical site only",
    "PS anatomical site only", "RP anatomical location only", "RD disease class only",
    "RD GOJ anatomical class only", "PS GOJ anatomical class only", "RP GOJ anatomical class only",
    "PS pretreatment local extent", "PS pretreatment nodal extent using TNM7 only",
    "PS pretreatment metastatic status", "TP recorded fitness", "TP recorded treatment intent",
    "TP recorded curative modality", "ST recorded surgery", "TP recorded comorbidity context",
    "RP resection local extent", "RP resection nodal extent using TNM7 only",
    "RP resection metastatic status", "RP positive-node count grouping", "RP differentiation",
    "RP histological phenotype", "RP dysplasia field context", "RP neoadjuvant exposure record",
    "ST surgical approach", "RP pathological response", "RP grouped pathological response",
    "TR recorded response", "Any Yes recorded in the single recurrence source column"
  ),
  source_column = c(
    "DI Patient Gender", "DI Age At Diagnosis", "RD Combined Tumour Site",
    "PS Pre Treatment Combined Tumour Site", "RP Location", "RD Histology",
    "RD Siewert Classification", "PS Siewert Classification", "RP Siewert Classification",
    "PS TStage Primary Tumour Final Pretreatment Staging",
    "PS NStage Primary Tumour Final Pretreatment Staging TNM7",
    "PS MStage Primary Tumour Final Pretreatment Staging", "TP Performance Status",
    "TP Treatment Intent", "TP Curative Treatment Modality", "ST Main Surgery", "TP Comorbidities",
    "RP TStage Primary Tumour", "RP Nstage RP TNM7", "RP MStage Distant Metastasis",
    "RP NStage Number Of Positive Nodes", "RP Tumour Grading Differentiation Status",
    "RP Signet Ring Cells Present", "RP Baretts Adjacent To Tumour Microscopic Dysplasia",
    "RP History Of Neo Adjuvant Therapy", "ST Procedure", "RP Mandard Score For Response",
    "RP Mandard Score For Response", "TR Response", "Has Original Disease Reoccurred"
  ),
  timepoint = c(
    "DI", "DI", "RD", "PS", "RP", "RD", "RD", "PS", "RP", "PS", "PS", "PS",
    "TP", "TP", "TP", "ST", "TP", "RP", "RP", "RP", "RP", "RP", "RP", "RP",
    "RP", "ST", "RP", "RP", "TR", "longitudinal source field"
  ),
  stringsAsFactors = FALSE
)

feature_long <- function(score_matrix, feature_type) {
  x <- as.data.frame(score_matrix, check.names = FALSE)
  x$subject_id <- rownames(x)
  x %>% pivot_longer(-subject_id, names_to = "feature", values_to = "score") %>%
    mutate(feature_type = feature_type)
}

prepare_variable_data <- function(metadata, scores_long, variable_name) {
  x <- metadata %>%
    select(subject_id, group = all_of(variable_name)) %>%
    mutate(group = trimws(as.character(group)), group = na_if(group, "")) %>%
    filter(!is.na(group))
  level_counts <- x %>% count(group, name = "level_n")
  retained <- level_counts %>% filter(level_n >= minimum_level_n) %>% pull(group)
  x <- x %>% filter(group %in% retained)
  if (n_distinct(x$group) < 2) return(NULL)
  x %>% inner_join(scores_long, by = "subject_id") %>%
    mutate(clinical_variable = variable_name)
}

compute_stats <- function(plot_data) {
  plot_data %>%
    group_by(clinical_variable, feature_type, feature) %>%
    group_modify(~ {
      d <- .x %>% filter(!is.na(group), !is.na(score))
      groups <- sort(unique(d$group))
      group_summary <- d %>% group_by(group) %>%
        summarise(n = n_distinct(subject_id), median = median(score), mean = mean(score), .groups = "drop")
      test_name <- if (length(groups) == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"
      test <- tryCatch(
        if (length(groups) == 2) wilcox.test(score ~ group, data = d, exact = FALSE) else kruskal.test(score ~ group, data = d),
        error = function(e) NULL
      )
      data.frame(
        test = test_name, n_groups = length(groups), n_subjects = n_distinct(d$subject_id),
        p_value = if (is.null(test)) NA_real_ else test$p.value,
        group_summary = paste0(group_summary$group, " (n=", group_summary$n,
                               ", median=", sprintf("%.3f", group_summary$median), ")", collapse = " | "),
        stringsAsFactors = FALSE
      )
    }) %>% ungroup() %>%
    group_by(clinical_variable, feature_type) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
}

build_plot_objects <- function() {
  required <- c(mp_score_path, state_score_path, metadata_path, grouping_path)
  if (any(!file.exists(required))) stop("Missing upstream OCCAMS survival output(s): ", paste(required[!file.exists(required)], collapse = ", "))
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_scores <- readRDS(mp_score_path)
  state_scores <- readRDS(state_score_path)
  if (nrow(metadata) != 281 || nrow(mp_scores) != 281 || nrow(state_scores) != 281) {
    stop("Expected exactly 281 QC-pass subjects in metadata, MP scores, and state scores")
  }
  if (!setequal(metadata$subject_id, rownames(mp_scores)) || !setequal(metadata$subject_id, rownames(state_scores))) {
    stop("Subject identities differ between metadata and score matrices")
  }

  grouping <- read.csv(grouping_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_labels <- setNames(grouping$plot_label, grouping$mp)
  score_long <- bind_rows(feature_long(mp_scores, "MP"), feature_long(state_scores, "State"))
  plot_data <- bind_rows(lapply(variable_config$variable, function(variable_name) {
    prepare_variable_data(metadata, score_long, variable_name)
  })) %>%
    left_join(variable_config, by = c("clinical_variable" = "variable")) %>%
    mutate(feature_label = ifelse(feature_type == "MP", unname(mp_labels[feature]), feature))

  tested_variables <- unique(plot_data$clinical_variable)
  inventory <- bind_rows(lapply(variable_config$variable, function(variable_name) {
    values <- metadata[[variable_name]]
    values <- trimws(as.character(values))
    values[values == ""] <- NA_character_
    counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
    data.frame(
      variable = variable_name,
      n_available = sum(!is.na(values)), n_missing = sum(is.na(values)),
      n_levels_observed = length(counts), n_levels_at_least_10 = sum(counts >= minimum_level_n),
      levels_and_counts = paste0(names(counts), "=", as.integer(counts), collapse = " | "),
      included_in_tests = variable_name %in% tested_variables,
      exclusion_rule = ifelse(variable_name %in% tested_variables, "", "fewer than two levels with n>=10"),
      stringsAsFactors = FALSE
    )
  })) %>% left_join(variable_config, by = "variable")

  stats <- compute_stats(plot_data)

  list(plot_data = plot_data, stats = stats,
       inventory = inventory, mp_order = grouping$mp,
       state_order = intersect(SCREF_PRIMARY_STATE_ORDER, colnames(state_scores)))
}

plot_boxplots <- function(objects) {
  for (type_name in c("MP", "State")) {
    out_path <- if (type_name == "MP") boxplot_path_mp else boxplot_path_state
    grDevices::cairo_pdf(out_path, width = 18, height = 9, onefile = TRUE)
    for (variable_name in unique(objects$plot_data$clinical_variable)) {
      variable_data <- objects$plot_data %>% filter(clinical_variable == variable_name)
      d <- variable_data %>% filter(feature_type == type_name)
      if (nrow(d) == 0) next
      variable_label <- unique(variable_data$label)
      
      feature_order <- if (type_name == "MP") objects$mp_order else objects$state_order
      d$feature_label <- factor(d$feature_label,
                                levels = unique(d$feature_label[match(feature_order, d$feature)]))
      d$group <- as.factor(d$group)
      
      legend_counts <- d %>% distinct(subject_id, group) %>% count(group, name = "n_samples")
      legend_labels <- setNames(
        paste0(legend_counts$group, " (n=", legend_counts$n_samples, ")"),
        as.character(legend_counts$group)
      )
      
      group_levels <- levels(d$group)
      palette <- scales::hue_pal()(length(group_levels))
      names(palette) <- group_levels
      
      stats_subset <- objects$stats %>% filter(clinical_variable == variable_name, feature_type == type_name)
      annot_df <- d %>%
        group_by(feature, feature_label) %>%
        summarise(y_pos = max(score, na.rm = TRUE) + max(0.08 * diff(range(score, na.rm=TRUE)), 0.02), .groups = "drop") %>%
        left_join(stats_subset, by = c("feature" = "feature")) %>%
        mutate(sig_label = case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )) %>%
        mutate(sig_label = ifelse(p_adj_bh < 0.05, sig_label, ""))
      
      p <- ggplot(d, aes(x = feature_label, y = score, fill = group, color = group)) +
        geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, linewidth = 0.4, color = "black") +
        geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), alpha = 0.7, size = 1.0, stroke = 0, show.legend = FALSE) +
        geom_text(data = annot_df %>% filter(sig_label != ""), aes(x = feature_label, y = y_pos, label = sig_label), inherit.aes = FALSE, size = 4, fontface = "bold") +
        scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
        scale_color_manual(values = palette, guide = "none", drop = FALSE) +
        labs(title = paste0(variable_label, " - ", type_name, " scores"),
             subtitle = "Sample-level scores; stars mark BH-adjusted p < 0.05 across clinical groups.",
             x = NULL, y = "GSVA score", fill = "Clinical group") +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 10, colour = "grey35"),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.line.x = element_blank(),
              legend.position = "top",
              legend.title = element_text(face = "bold"),
              plot.margin = margin(12, 16, 12, 12))
      print(p)
    }
    dev.off()
  }
}



if (replot_only) {
  if (!file.exists(plot_data_path)) stop("SCREF_REPLOT_ONLY requires ", plot_data_path)
  objects <- readRDS(plot_data_path)
} else {
  if (!force_rebuild && file.exists(plot_data_path)) {
    objects <- readRDS(plot_data_path)
  } else {
    objects <- build_plot_objects()
    saveRDS(objects, plot_data_path, compress = "xz")
  }
  write.csv(objects$inventory, inventory_path, row.names = FALSE, na = "")
  write.csv(objects$stats, stats_path, row.names = FALSE, na = "")
}

plot_boxplots(objects)

summary_df <- objects$stats %>% group_by(feature_type) %>%
  summarise(qc_pass_subjects = 281L, n_clinical_variables = n_distinct(clinical_variable),
            n_tests = n(), n_nominal_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
            n_bh_fdr_lt_0_05 = sum(p_adj_bh < 0.05, na.rm = TRUE), .groups = "drop")
write.csv(summary_df, summary_path, row.names = FALSE)
writeLines(c(
  "OCCAMS current-centred bulk clinical association run",
  paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "qc_pass_subjects=281",
  paste0("minimum_subjects_per_displayed_level=", minimum_level_n),
  paste0("clinical_variables_tested=", n_distinct(objects$plot_data$clinical_variable)),
  paste0("association_tests=", nrow(objects$stats)),
  "result=PASS"
), report_path)
message("Completed OCCAMS QC-pass clinical association workflow.")
