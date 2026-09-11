#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active diagnostic
#   Script: analysis/spatial/visium_hd_binned_reference_selection_audit.R
#   Description: compare legacy RCTD, current random-balanced, and current
#     manual/RCTD-concordant normal-reference definitions for binned InferCNA.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_filter_malignancy_methodology.md
#   Inputs:
#     current filtered annotation and malignancy tables; legacy binned
#     malignancy cell tables.
#   Outputs:
#     ref_outs/visium_hd_outs/malignancy/reference_audit/{tables,figures,logs}/
#     and updates/new_updates/summaries/visium_hd_binned_reference_audit.csv
#   Cache/replot: inexpensive table-only diagnostic; always recomputed.
#   Run: Rscript analysis/spatial/visium_hd_binned_reference_selection_audit.R
#   Environment: dmtcp
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)
samples <- c("SUR1231", "FFPEA1", "FFPED1")
normal_types <- c("endothelial", "macrophage", "fibroblast")
filter_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "post_annotation_filter", "tables")
malignancy_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "malignancy", "tables")
legacy_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "legacy_visiumhd", "malignancy", "tables")
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "malignancy", "reference_audit")
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
for (subdir in c("tables", "figures", "logs")) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

metric_output_path <- file.path(output_dir, "tables", "Auto_visium_hd_binned_reference_metric_comparison.csv")
count_output_path <- file.path(output_dir, "tables", "Auto_visium_hd_binned_reference_definition_counts.csv")
previous_metric_rows <- data.table()
previous_count_rows <- data.table()
if (file.exists(metric_output_path)) {
  previous_metric_rows <- fread(metric_output_path)[
    reference_definition %in% c("Current random-balanced", "Previous current random-balanced")
  ]
  previous_metric_rows[, reference_definition := "Previous current random-balanced"]
}
if (file.exists(count_output_path)) {
  previous_count_rows <- fread(count_output_path)
}

count_rows <- list()
metric_rows <- list()
classification_rows <- list()
rctd_cross_rows <- list()
run_log <- c(paste0("start=", format(Sys.time(), tz = "Europe/London")))

pdf_path <- file.path(output_dir, "figures", "Auto_visium_hd_binned_reference_selection_audit.pdf")
pdf(pdf_path, width = 12, height = 8.5, onefile = TRUE, useDingbats = FALSE)

for (sample_name in samples) {
  annotation <- fread(file.path(filter_dir, paste0("Auto_", sample_name, "_binned_filtered_annotations.csv.gz")))
  current <- fread(file.path(malignancy_dir, paste0("Auto_", sample_name, "_binned_malignancy.csv.gz")))
  legacy <- fread(file.path(legacy_dir, paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz")))

  annotation[, Auto_reference_exact_manual_label :=
    !grepl("\\|", Auto_annotation_celltype) & Auto_annotation_celltype == Auto_postfilter_celltype]
  annotation[, Auto_reference_rctd_concordant :=
    Auto_postfilter_keep %in% TRUE &
      Auto_reference_exact_manual_label %in% TRUE &
      Auto_rctd_is_singlet %in% TRUE &
      Auto_rctd_first_type == Auto_postfilter_celltype &
      Auto_postfilter_celltype %in% normal_types]

  for (reference_type in normal_types) {
    sample_value <- sample_name
    reference_value <- reference_type
    count_rows[[paste(sample_name, reference_type, sep = "_")]] <- data.table(
      sample = sample_name,
      reference_type = reference_type,
      n_current_postfilter = annotation[
        Auto_postfilter_keep %in% TRUE & Auto_postfilter_celltype == reference_type, .N
      ],
      n_current_exact_label = annotation[
        Auto_postfilter_keep %in% TRUE & Auto_postfilter_celltype == reference_type &
          Auto_reference_exact_manual_label %in% TRUE, .N
      ],
      n_current_rctd_concordant = annotation[
        Auto_reference_rctd_concordant %in% TRUE & Auto_postfilter_celltype == reference_type, .N
      ],
      n_current_rctd_singlet_reference = annotation[
        Auto_rctd_is_singlet %in% TRUE & Auto_rctd_first_type == reference_type, .N
      ],
      n_previous_random_balanced = if (nrow(previous_count_rows)) {
        value <- previous_count_rows[
          sample == sample_value & reference_type == reference_value,
          n_current_selected
        ]
        if (length(value)) value[[1L]] else NA_integer_
      } else {
        NA_integer_
      },
      n_current_selected = current[
        is_reference %in% TRUE & Auto_postfilter_celltype == reference_type, .N
      ],
      n_legacy_rctd_reference = legacy[
        is_reference %in% TRUE & Auto_annotation_celltype == reference_type, .N
      ]
    )
  }

  rctd_cross_rows[[sample_name]] <- annotation[
    Auto_rctd_is_singlet %in% TRUE & Auto_rctd_first_type %in% normal_types,
    .(n_bins = .N),
    by = .(
      rctd_reference_type = Auto_rctd_first_type,
      current_postfilter_celltype = Auto_postfilter_celltype,
      current_postfilter_keep = Auto_postfilter_keep
    )
  ]
  rctd_cross_rows[[sample_name]][, sample := sample_name]
  setcolorder(rctd_cross_rows[[sample_name]], c(
    "sample", "rctd_reference_type", "current_postfilter_celltype",
    "current_postfilter_keep", "n_bins"
  ))

  current_reference <- current[is_reference %in% TRUE]
  current_reference[, reference_definition := "Final current RCTD/manual concordant"]
  legacy_reference <- legacy[is_reference %in% TRUE]
  legacy_reference[, reference_definition := "Legacy RCTD"]

  metric_sets <- list(
    "Final current RCTD/manual concordant" = current_reference,
    "Legacy RCTD" = legacy_reference
  )
  for (definition in names(metric_sets)) {
    values <- metric_sets[[definition]]
    if (!nrow(values)) next
    metric_rows[[paste(sample_name, definition, sep = "_")]] <- data.table(
      sample = sample_name,
      reference_definition = definition,
      n_reference = nrow(values),
      signal_mean = mean(values$cna.signal, na.rm = TRUE),
      signal_sd = sd(values$cna.signal, na.rm = TRUE),
      signal_threshold_mean_plus_1sd = mean(values$cna.signal, na.rm = TRUE) + sd(values$cna.signal, na.rm = TRUE),
      correlation_mean = mean(values$cna.cor, na.rm = TRUE),
      correlation_sd = sd(values$cna.cor, na.rm = TRUE),
      correlation_threshold_mean_plus_1sd = mean(values$cna.cor, na.rm = TRUE) + sd(values$cna.cor, na.rm = TRUE)
    )
  }

  if (nrow(current_reference) >= 20L) {
    signal_threshold <- mean(current_reference$cna.signal, na.rm = TRUE) +
      sd(current_reference$cna.signal, na.rm = TRUE)
    correlation_threshold <- mean(current_reference$cna.cor, na.rm = TRUE) +
      sd(current_reference$cna.cor, na.rm = TRUE)
    targets <- current[is_epithelial_target %in% TRUE]
    classification_rows[[sample_name]] <- data.table(
      sample = sample_name,
      threshold_basis = "final_current_RCTD_manual_concordant_reference",
      signal_threshold = signal_threshold,
      correlation_threshold = correlation_threshold,
      n_epithelial = nrow(targets),
      n_both_above = sum(targets$cna.signal > signal_threshold & targets$cna.cor > correlation_threshold),
      pct_both_above = 100 * mean(targets$cna.signal > signal_threshold & targets$cna.cor > correlation_threshold)
    )
  }

  plot_data <- copy(current_reference)
  p <- ggplot(plot_data, aes(x = cna.signal, y = cna.cor, colour = reference_definition)) +
    geom_point(size = 1.25, alpha = 0.6) +
    scale_colour_manual(values = c("Final current RCTD/manual concordant" = "#0072B2"), drop = FALSE) +
    labs(title = sample_name, x = "CNA signal", y = "CNA correlation", colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  print(p)
  run_log <- c(run_log, paste0(
    "sample=", sample_name,
    "; current_reference=", nrow(current_reference),
    "; legacy_reference=", nrow(legacy_reference)
  ))
}
dev.off()

count_table <- rbindlist(count_rows, use.names = TRUE, fill = TRUE)
metric_table <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
if (nrow(previous_metric_rows)) {
  metric_table <- rbindlist(list(previous_metric_rows, metric_table), use.names = TRUE, fill = TRUE)
}
classification_table <- rbindlist(classification_rows, use.names = TRUE, fill = TRUE)
rctd_cross_table <- rbindlist(rctd_cross_rows, use.names = TRUE, fill = TRUE)
fwrite(count_table, count_output_path)
fwrite(metric_table, metric_output_path)
fwrite(classification_table, file.path(output_dir, "tables", "Auto_visium_hd_binned_reference_threshold_sensitivity.csv"))
fwrite(rctd_cross_table, file.path(output_dir, "tables", "Auto_visium_hd_binned_rctd_reference_current_annotation_crosstab.csv"))
fwrite(count_table, file.path(summary_dir, "visium_hd_binned_reference_audit.csv"))
run_log <- c(run_log, paste0("end=", format(Sys.time(), tz = "Europe/London")))
writeLines(run_log, file.path(output_dir, "logs", "Auto_visium_hd_binned_reference_audit_run_summary.txt"))
####################
