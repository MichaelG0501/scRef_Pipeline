#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_compare_annotation_infercna.R
#   Description: Compare InferCNA scatter distributions across binned RCTD,
#     binned custom, segmented custom, and spatial-corrected annotation on shared per-sample axes.
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: ref_outs/visium_hd_outs/malignancy/tables/Auto_<sample>_<mode>_infercna_cells.csv.gz
#   Outputs: ref_outs/visium_hd_outs/malignancy/figures/Auto_visiumhd_four_method_infercna_comparison.pdf;
#     ref_outs/visium_hd_outs/malignancy/tables/Auto_visiumhd_four_method_infercna_comparison_summary.csv;
#     updates/new_updates/summaries/visiumhd_four_method_infercna_comparison_summary.csv
#   Cache/replot: Plot-only; always reuses completed malignancy cell tables.
#   Run: qsub analysis/spatial/visiumhd_compare_annotation_infercna.sh
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

####################
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)
base_dir <- file.path(wd, "ref_outs/visium_hd_outs/malignancy")
table_dir <- file.path(base_dir, "tables")
figure_dir <- file.path(base_dir, "figures")
summary_dir <- file.path(wd, "updates/new_updates/summaries")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("SUR1231", "FFPEA1", "FFPED1")
mode_labels <- c(
  binned = "Binned RCTD",
  custom = "Binned custom",
  segmented = "Segmented custom",
  spatial = "Segmented spatial"
)

read_mode <- function(sample_name, mode) {
  path <- file.path(table_dir, paste0("Auto_", sample_name, "_", mode, "_infercna_cells.csv.gz"))
  if (!file.exists(path)) stop("Missing malignancy table: ", path)
  x <- read.csv(path, check.names = FALSE)
  required <- c("barcode", "cna.signal", "cna.cor", "is_reference", "is_epithelial_target", "Auto_malignancy",
                "Auto_cna_signal_threshold", "Auto_cna_cor_threshold")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x$is_reference <- as.logical(x$is_reference)
  x$is_epithelial_target <- as.logical(x$is_epithelial_target)
  x <- x[is.finite(x$cna.signal) & is.finite(x$cna.cor) &
           (x$is_reference | x$is_epithelial_target), , drop = FALSE]
  x$sample <- sample_name
  x$mode <- mode
  x$method <- unname(mode_labels[mode])
  x$plot_group <- ifelse(x$is_reference, "Reference", "Non-malignant epithelial")
  x$plot_group[x$Auto_malignancy == "malignant_level_1"] <- "Malignant level 1 (CNA)"
  x$plot_group[x$Auto_malignancy == "malignant_level_2"] <- "Malignant level 2 (signature)"
  x[, c(
    "barcode", "cna.signal", "cna.cor", "is_reference", "is_epithelial_target",
    "Auto_malignancy", "Auto_cna_signal_threshold", "Auto_cna_cor_threshold",
    "sample", "mode", "method", "plot_group"
  ), drop = FALSE]
}

all_data <- do.call(rbind, lapply(samples, function(sample_name) {
  do.call(rbind, lapply(names(mode_labels), function(mode) read_mode(sample_name, mode)))
}))
all_data$method <- factor(all_data$method, levels = unname(mode_labels))
all_data$plot_group <- factor(all_data$plot_group, levels = c(
  "Reference", "Non-malignant epithelial", "Malignant level 1 (CNA)", "Malignant level 2 (signature)"
))

summary_rows <- list()
for (sample_name in samples) {
  sample_data <- all_data[all_data$sample == sample_name, , drop = FALSE]
  binned <- sample_data[sample_data$mode == "binned", c("barcode", "cna.signal", "cna.cor")]
  custom <- sample_data[sample_data$mode == "custom", c("barcode", "cna.signal", "cna.cor")]
  paired <- merge(binned, custom, by = "barcode", suffixes = c("_rctd", "_custom"))
  for (mode in names(mode_labels)) {
    x <- sample_data[sample_data$mode == mode, , drop = FALSE]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      sample = sample_name,
      mode = mode,
      method = unname(mode_labels[mode]),
      n_points = nrow(x),
      n_reference = sum(as.logical(x$is_reference), na.rm = TRUE),
      cna_signal_threshold = unique(x$Auto_cna_signal_threshold)[1],
      cna_cor_threshold = unique(x$Auto_cna_cor_threshold)[1],
      paired_binned_n = if (mode == "custom") nrow(paired) else NA_integer_,
      paired_binned_signal_cor = if (mode == "custom") cor(paired$cna.signal_rctd, paired$cna.signal_custom) else NA_real_,
      paired_binned_cna_cor_cor = if (mode == "custom") cor(paired$cna.cor_rctd, paired$cna.cor_custom) else NA_real_,
      paired_binned_signal_max_abs_diff = if (mode == "custom") max(abs(paired$cna.signal_rctd - paired$cna.signal_custom)) else NA_real_,
      paired_binned_cna_cor_max_abs_diff = if (mode == "custom") max(abs(paired$cna.cor_rctd - paired$cna.cor_custom)) else NA_real_
    )
  }
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(table_dir, "Auto_visiumhd_four_method_infercna_comparison_summary.csv"), row.names = FALSE)
write.csv(summary_table, file.path(summary_dir, "visiumhd_four_method_infercna_comparison_summary.csv"), row.names = FALSE)

pdf_path <- file.path(figure_dir, "Auto_visiumhd_four_method_infercna_comparison.pdf")
pdf(pdf_path, width = 20, height = 8, useDingbats = FALSE)
for (sample_name in samples) {
  x <- all_data[all_data$sample == sample_name, , drop = FALSE]
  thresholds <- unique(x[, c("method", "Auto_cna_signal_threshold", "Auto_cna_cor_threshold")])
  x_upper <- unname(quantile(x$cna.signal, 0.995, na.rm = TRUE)) * 1.05
  y_limits <- unname(quantile(x$cna.cor, c(0.005, 0.995), na.rm = TRUE))
  p <- ggplot(x, aes(cna.signal, cna.cor, colour = plot_group)) +
    geom_point(size = 0.55, alpha = 0.55) +
    geom_vline(data = thresholds, aes(xintercept = Auto_cna_signal_threshold),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 0.45, colour = "grey25") +
    geom_hline(data = thresholds, aes(yintercept = Auto_cna_cor_threshold),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 0.45, colour = "grey25") +
    facet_wrap(~method, nrow = 1) +
    coord_cartesian(xlim = c(0, x_upper), ylim = y_limits) +
    scale_colour_manual(values = c(
      "Reference" = "#8A8A8A",
      "Non-malignant epithelial" = "#4DAF4A",
      "Malignant level 1 (CNA)" = "#E41A1C",
      "Malignant level 2 (signature)" = "#984EA3"
    ), drop = FALSE) +
    scale_x_continuous(labels = scientific) +
    labs(title = sample_name, x = "CNA signal", y = "CNA correlation", colour = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 11)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1)))
  print(p)
}
dev.off()
####################
