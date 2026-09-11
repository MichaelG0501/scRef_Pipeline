#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_keratinocyte_evidence_audit.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: final binned InferCNA cell tables with segmented keratinocyte and
#     optional complete-CNA-profile fields.
#   Outputs: ref_outs/visium_hd_outs/malignancy/figures/
#     Auto_visiumhd_keratinocyte_vs_epithelial_evidence_audit.pdf and summary.
#   Cache/replot: plot-only; no InferCNA calculation is performed.
#   Run: Rscript analysis/spatial/visiumhd_keratinocyte_evidence_audit.R
#     --samples SUR1231 FFPEA1 FFPED1 --malignancy-dir <dir>
#   Environment: dmtcp.
####################

####################
# Compact evidence audit for keratinocyte-dominant RCTD epithelial bins versus
# non-keratinocyte malignant epithelial bins. It reports overlapping scalar
# distributions and profile-correlation evidence without relabelling cells.
####################
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

parse_cli <- function(args) {
  out <- list(samples = character())
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% c("samples", "malignancy-dir")) stop("Unknown argument: ", args[[i]])
    i <- i + 1L
    if (i > length(args) || startsWith(args[[i]], "--")) stop("Missing value for --", key)
    if (key == "samples") {
      start <- i
      while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1L
      out$samples <- args[start:(i - 1L)]
    } else {
      out[[gsub("-", "_", key)]] <- args[[i]]
      i <- i + 1L
    }
  }
  out
}

as_logical <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")

auc_keratinocyte_high <- function(keratinocyte, epithelial) {
  keratinocyte <- keratinocyte[is.finite(keratinocyte)]
  epithelial <- epithelial[is.finite(epithelial)]
  if (!length(keratinocyte) || !length(epithelial)) return(NA_real_)
  ranks <- rank(c(keratinocyte, epithelial), ties.method = "average")
  (sum(ranks[seq_along(keratinocyte)]) - length(keratinocyte) * (length(keratinocyte) + 1) / 2) / (length(keratinocyte) * length(epithelial))
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!length(args$samples) || is.null(args$malignancy_dir)) stop("Required: --samples --malignancy-dir")
malignancy_dir <- normalizePath(args$malignancy_dir, mustWork = TRUE)
figure_dir <- file.path(malignancy_dir, "figures")
####################
# Create output tiers for standalone PBS execution.
####################
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(malignancy_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
####################
summary_rows <- list()
pdf_path <- file.path(figure_dir, "Auto_visiumhd_keratinocyte_vs_epithelial_evidence_audit.pdf")
pdf(pdf_path, width = 15, height = 10, useDingbats = FALSE)

for (sample_name in args$samples) {
  cell_path <- file.path(malignancy_dir, "tables", paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz"))
  if (!file.exists(cell_path)) stop("Missing binned malignancy table: ", cell_path)
  ####################
  # Use full-height pages and the established InferCNA scatter presentation.
  # Scalar evidence is intentionally restricted to signal, correlation, and
  # the malignancy signature; full-profile fields are not displayed here.
  ####################
  full_data <- read.csv(cell_path, stringsAsFactors = FALSE, check.names = FALSE)
  full_data$is_epithelial_target <- as_logical(full_data$is_epithelial_target)
  full_data$is_reference <- as_logical(full_data$is_reference)
  keratinocyte <- full_data$is_epithelial_target & full_data$Auto_segmented_keratinocyte_n_cells >= 1L & full_data$Auto_segmented_keratinocyte_fraction >= 0.5
  anchor_label <- if ("Auto_malignancy_before_keratinocyte_profile" %in% colnames(full_data)) full_data$Auto_malignancy_before_keratinocyte_profile else full_data$Auto_malignancy
  epithelial <- full_data$is_epithelial_target & !keratinocyte & anchor_label == "malignant_level_1"
  full_data$plot_group <- "Other epithelial"
  full_data$plot_group[full_data$is_reference] <- "Reference"
  full_data$plot_group[epithelial] <- "Non-keratinocyte malignant epithelial"
  full_data$plot_group[keratinocyte] <- "Keratinocyte-dominant RCTD epithelial"
  full_data$plot_group <- factor(full_data$plot_group, levels = c("Reference", "Other epithelial", "Non-keratinocyte malignant epithelial", "Keratinocyte-dominant RCTD epithelial"))
  threshold_signal <- full_data$Auto_cna_signal_threshold[which(is.finite(full_data$Auto_cna_signal_threshold))[1]]
  threshold_cor <- full_data$Auto_cna_cor_threshold[which(is.finite(full_data$Auto_cna_cor_threshold))[1]]
  threshold_signature <- full_data$Auto_cancer_signature_threshold[which(is.finite(full_data$Auto_cancer_signature_threshold))[1]]
  cna_sd_k <- full_data$Auto_cna_sd_k[which(is.finite(full_data$Auto_cna_sd_k))[1]]
  if (!is.finite(cna_sd_k)) cna_sd_k <- 1
  scatter_data <- do.call(rbind, list(
    transform(full_data, facet_group = "All"),
    transform(full_data[full_data$is_reference, , drop = FALSE], facet_group = "Reference"),
    transform(full_data[epithelial, , drop = FALSE], facet_group = "Epithelial level 1"),
    transform(full_data[keratinocyte, , drop = FALSE], facet_group = "Keratinocyte dominant")
  ))
  scatter_data$facet_group <- factor(scatter_data$facet_group, levels = c("All", "Reference", "Epithelial level 1", "Keratinocyte dominant"))
  colours <- c("Reference" = "#9E9E9E", "Other epithelial" = "#4DAF4A", "Non-keratinocyte malignant epithelial" = "#E41A1C", "Keratinocyte-dominant RCTD epithelial" = "#FFD92F")
  p_scatter <- ggplot(scatter_data, aes(cna.signal, cna.cor, colour = plot_group)) +
    geom_point(size = 0.35, alpha = 0.65) +
    geom_vline(xintercept = threshold_signal, linetype = "dotted", linewidth = 0.6, colour = "grey20") +
    geom_hline(yintercept = threshold_cor, linetype = "dotted", linewidth = 0.6, colour = "grey20") +
    facet_wrap(~ facet_group, ncol = 2) +
    scale_colour_manual(values = colours, drop = FALSE) +
    scale_x_continuous(labels = scales::scientific) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = paste(sample_name, "binned InferCNA"), x = "CNA signal", y = "CNA correlation", colour = NULL) +
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 12))
  print(p_scatter)
  data <- full_data[keratinocyte | epithelial, , drop = FALSE]
  data$comparison_group <- factor(ifelse(keratinocyte[keratinocyte | epithelial], "Keratinocyte dominant", "Epithelial level 1"), levels = c("Epithelial level 1", "Keratinocyte dominant"))
  metrics <- c("cna.signal" = "CNA signal", "cna.cor" = "CNA correlation", "Auto_cancer_signature_score" = "Cancer signature")
  long <- do.call(rbind, lapply(names(metrics), function(metric) data.frame(comparison_group = data$comparison_group, metric = unname(metrics[[metric]]), value = data[[metric]], stringsAsFactors = FALSE)))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$metric <- factor(long$metric, levels = unname(metrics))
  thresholds <- data.frame(
    metric = factor(unname(metrics), levels = unname(metrics)),
    threshold = c(threshold_signal, threshold_cor, threshold_signature),
    label = c(paste0("threshold = ", format(threshold_signal, scientific = TRUE, digits = 2)), paste0("threshold = ", sprintf("%.3f", threshold_cor)), paste0("threshold = ", sprintf("%.2f", threshold_signature))),
    stringsAsFactors = FALSE
  )
  comparison_colours <- c("Epithelial level 1" = "#E41A1C", "Keratinocyte dominant" = "#FFD92F")
  p_density <- ggplot(long, aes(value, fill = comparison_group, colour = comparison_group)) +
    geom_density(alpha = 0.25, linewidth = 0.8, adjust = 1.1) +
    geom_vline(data = thresholds, aes(xintercept = threshold), inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7, colour = "grey20") +
    geom_text(data = thresholds, aes(x = threshold, y = Inf, label = label), inherit.aes = FALSE, angle = 90, hjust = 1.05, vjust = -0.35, size = 3.6) +
    scale_fill_manual(values = comparison_colours) + scale_colour_manual(values = comparison_colours) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    labs(title = paste(sample_name, "scalar evidence"), x = NULL, y = "Density", fill = NULL, colour = NULL) +
    theme_classic(base_size = 15) + theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12))
  p_box <- ggplot(long, aes(comparison_group, value, fill = comparison_group)) +
    geom_boxplot(outlier.size = 0.3, width = 0.58) +
    geom_hline(data = thresholds, aes(yintercept = threshold), inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7, colour = "grey20") +
    geom_text(data = thresholds, aes(x = 1.5, y = threshold, label = label), inherit.aes = FALSE, vjust = -0.45, size = 3.6) +
    scale_fill_manual(values = comparison_colours) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10), strip.text = element_text(face = "bold", size = 12))
  print((p_density / p_box) + plot_layout(heights = c(1.05, 1)))
  summary_rows[[sample_name]] <- data.frame(sample = sample_name, n_keratinocyte_dominant = sum(keratinocyte), n_nonkeratinocyte_malignant = sum(epithelial), stringsAsFactors = FALSE)
  ####################
}
dev.off()
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(malignancy_dir, "tables", "Auto_visiumhd_keratinocyte_vs_epithelial_evidence_audit_summary.csv"), row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries"); dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(summary_dir, "visiumhd_keratinocyte_vs_epithelial_evidence_audit_summary.csv"), row.names = FALSE)
