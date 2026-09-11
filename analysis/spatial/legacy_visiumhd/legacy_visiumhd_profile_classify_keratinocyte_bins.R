#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_profile_classify_keratinocyte_bins.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: cached binned InferCNA matrices/cell tables with segmented
#     keratinocyte composition fields.
#   Outputs: corrected binned malignancy tables/scatters, profile diagnostics,
#     and a keratinocyte profile-classification summary.
#   Cache/replot: reuses saved InferCNA matrices; no CNA inference is rerun.
#   Run: Rscript analysis/spatial/visiumhd_profile_classify_keratinocyte_bins.R
#     --samples SUR1231 FFPEA1 FFPED1 --malignancy-dir <dir>
#   Environment: dmtcp.
####################

####################
# Classify keratinocyte-dominant RCTD epithelial bins from their complete CNA
# profile. Scalar CNA signal and cancer-signature score are retained as audits
# but do not define normal versus malignant keratinocyte calls.
####################
suppressPackageStartupMessages({
  library(ggplot2)
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

correlation_to_profile <- function(cna, profile) {
  profile <- as.numeric(profile)
  profile <- profile - mean(profile, na.rm = TRUE)
  column_mean <- colMeans(cna, na.rm = TRUE)
  column_sum_sq <- colSums(cna^2, na.rm = TRUE)
  numerator <- as.numeric(crossprod(profile, cna))
  denominator <- sqrt(sum(profile^2, na.rm = TRUE) * pmax(column_sum_sq - nrow(cna) * column_mean^2, .Machine$double.eps))
  numerator / denominator
}

plot_profile_diagnostics <- function(data, sample_name, figure_dir) {
  plot_data <- data[data$Auto_keratinocyte_cna_profile_class %in% c("keratinocyte_cna_like", "keratinocyte_normal_like", "keratinocyte_indeterminate") | data$Auto_malignancy == "malignant_level_1", , drop = FALSE]
  plot_data$plot_group <- ifelse(plot_data$Auto_malignancy == "malignant_level_1" & plot_data$Auto_keratinocyte_cna_profile_class == "not_keratinocyte", "non_keratinocyte_malignant", plot_data$Auto_keratinocyte_cna_profile_class)
  plot_data$plot_group <- factor(plot_data$plot_group, levels = c("non_keratinocyte_malignant", "keratinocyte_cna_like", "keratinocyte_normal_like", "keratinocyte_indeterminate"))
  colours <- c(non_keratinocyte_malignant = "#E41A1C", keratinocyte_cna_like = "#984EA3", keratinocyte_normal_like = "#4DAF4A", keratinocyte_indeterminate = "#FFD92F")
  p1 <- ggplot(plot_data, aes(plot_group, Auto_keratinocyte_tumour_profile_correlation, fill = plot_group)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.75) + geom_boxplot(width = 0.16, outlier.size = 0.25) +
    scale_fill_manual(values = colours, drop = FALSE) + labs(title = paste(sample_name, "CNA-profile classification"), x = NULL, y = "Correlation to non-keratinocyte malignant CNA centroid") +
    theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  p2 <- ggplot(plot_data, aes(cna.signal, cna.cor, colour = plot_group)) +
    geom_point(size = 0.65, alpha = 0.55) + scale_colour_manual(values = colours, drop = FALSE) + scale_x_continuous(labels = scales::scientific) +
    labs(x = "CNA signal", y = "CNA correlation", colour = NULL) + theme_classic(base_size = 13) + theme(legend.position = "bottom")
  pdf(file.path(figure_dir, paste0("Auto_", sample_name, "_binned_keratinocyte_cna_profile_diagnostics.pdf")), width = 15, height = 6)
  print(p1); print(p2); dev.off()
  ggsave(file.path(figure_dir, paste0("Auto_", sample_name, "_binned_keratinocyte_cna_profile_scatter.png")), p2, width = 8, height = 6, dpi = 300)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!length(args$samples) || is.null(args$malignancy_dir)) stop("Required: --samples --malignancy-dir")
malignancy_dir <- normalizePath(args$malignancy_dir, mustWork = TRUE)
####################
# Create output tiers for standalone PBS execution.
####################
dir.create(file.path(malignancy_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(malignancy_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
####################
summary_rows <- list()

for (sample_name in args$samples) {
  cell_path <- file.path(malignancy_dir, "tables", paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz"))
  cna_path <- file.path(malignancy_dir, "intermediate", paste0("Auto_", sample_name, "_binned_infercna_outs.rds"))
  if (!file.exists(cell_path) || !file.exists(cna_path)) stop("Missing cached binned CNA inputs for ", sample_name)
  cells <- read.csv(cell_path, stringsAsFactors = FALSE, check.names = FALSE)
  cna <- readRDS(cna_path)
  if (is.list(cna)) stop("Expected infercna matrix, found list for ", sample_name)
  shared <- intersect(colnames(cna), cells$barcode)
  if (length(shared) < 100L) stop("Too few matching cached CNA profiles for ", sample_name)
  cna <- cna[, shared, drop = FALSE]
  cells <- cells[match(shared, cells$barcode), , drop = FALSE]
  cells$is_reference <- as_logical(cells$is_reference)
  cells$is_epithelial_target <- as_logical(cells$is_epithelial_target)
  if (!"Auto_segmented_keratinocyte_fraction" %in% colnames(cells)) stop("Missing segmented keratinocyte fields for ", sample_name)
  keratinocyte <- cells$is_epithelial_target & cells$Auto_segmented_keratinocyte_n_cells >= 1L & cells$Auto_segmented_keratinocyte_fraction >= 0.5
  prior_malignancy <- if ("Auto_malignancy_before_keratinocyte_exclusion" %in% colnames(cells)) cells$Auto_malignancy_before_keratinocyte_exclusion else cells$Auto_malignancy
  anchors <- cells$is_epithelial_target & !keratinocyte & prior_malignancy == "malignant_level_1"
  references <- cells$is_reference
  if (sum(anchors) < 30L || sum(references) < 30L) stop("Insufficient malignant anchors or references for ", sample_name)
  ####################
  # A sample with no segmented keratinocyte-dominant bins requires no profile
  # correction. Retain its malignancy calls and emit an explicit valid status
  # instead of stopping the complete three-sample workflow.
  ####################
  if (sum(keratinocyte) < 1L) {
    cells$Auto_keratinocyte_cna_profile_class <- "not_keratinocyte"
    cells$Auto_keratinocyte_tumour_profile_correlation <- NA_real_
    cells$Auto_keratinocyte_normal_profile_correlation <- NA_real_
    cells$Auto_keratinocyte_tumour_normal_margin <- NA_real_
    cells$Auto_keratinocyte_tumour_profile_threshold <- NA_real_
    cells$Auto_keratinocyte_margin_threshold <- NA_real_
    cells$Auto_malignancy_before_keratinocyte_profile <- cells$Auto_malignancy
    cells$Auto_malignant <- cells$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")
    con <- gzfile(cell_path, open = "wt"); write.csv(cells, con, row.names = FALSE); close(con)
    target <- cells[cells$is_epithelial_target, , drop = FALSE]
    summary_rows[[sample_name]] <- data.frame(
      sample = sample_name, mode = "binned", profile_status = "no_keratinocyte_targets",
      n_keratinocyte_cna_like = 0L, n_keratinocyte_normal_like = 0L,
      n_keratinocyte_indeterminate = 0L, tumour_profile_threshold = NA_real_,
      margin_threshold = NA_real_, n_malignant_level_1 = sum(target$Auto_malignancy == "malignant_level_1"),
      n_malignant_level_2 = sum(target$Auto_malignancy == "malignant_level_2"),
      n_malignant = sum(target$Auto_malignant), pct_malignant_epithelial = 100 * mean(target$Auto_malignant),
      stringsAsFactors = FALSE
    )
    next
  }
  ####################
  tumour_profile <- rowMeans(cna[, anchors, drop = FALSE], na.rm = TRUE)
  normal_profile <- rowMeans(cna[, references, drop = FALSE], na.rm = TRUE)
  tumour_cor <- correlation_to_profile(cna, tumour_profile)
  normal_cor <- correlation_to_profile(cna, normal_profile)
  margin <- tumour_cor - normal_cor
  tumour_threshold <- as.numeric(quantile(tumour_cor[references], 0.99, na.rm = TRUE))
  margin_threshold <- as.numeric(quantile(margin[references], 0.95, na.rm = TRUE))
  profile_class <- rep("not_keratinocyte", nrow(cells))
  profile_class[keratinocyte & tumour_cor <= tumour_threshold & margin <= margin_threshold] <- "keratinocyte_normal_like"
  profile_class[keratinocyte & tumour_cor > tumour_threshold & margin > margin_threshold] <- "keratinocyte_cna_like"
  profile_class[keratinocyte & profile_class == "not_keratinocyte"] <- "keratinocyte_indeterminate"
  cells$Auto_keratinocyte_cna_profile_class <- profile_class
  cells$Auto_keratinocyte_tumour_profile_correlation <- tumour_cor
  cells$Auto_keratinocyte_normal_profile_correlation <- normal_cor
  cells$Auto_keratinocyte_tumour_normal_margin <- margin
  cells$Auto_keratinocyte_tumour_profile_threshold <- tumour_threshold
  cells$Auto_keratinocyte_margin_threshold <- margin_threshold
  cells$Auto_malignancy_before_keratinocyte_profile <- cells$Auto_malignancy
  cells$Auto_malignancy[profile_class == "keratinocyte_normal_like"] <- "normal_keratinocyte"
  cells$Auto_malignancy[profile_class == "keratinocyte_indeterminate"] <- "keratinocyte_indeterminate"
  cells$Auto_malignancy[profile_class == "keratinocyte_cna_like"] <- "malignant_level_1"
  cells$Auto_malignancy_evidence[profile_class == "keratinocyte_normal_like"] <- "keratinocyte_normal_cna_profile"
  cells$Auto_malignancy_evidence[profile_class == "keratinocyte_indeterminate"] <- "keratinocyte_indeterminate_cna_profile"
  cells$Auto_malignancy_evidence[profile_class == "keratinocyte_cna_like"] <- "keratinocyte_malignant_cna_profile"
  cells$Auto_malignant <- cells$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")
  con <- gzfile(cell_path, open = "wt"); write.csv(cells, con, row.names = FALSE); close(con)
  plot_profile_diagnostics(cells, sample_name, file.path(malignancy_dir, "figures"))
  target <- cells[cells$is_epithelial_target, , drop = FALSE]
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name, mode = "binned", profile_status = "complete",
    n_keratinocyte_cna_like = sum(target$Auto_keratinocyte_cna_profile_class == "keratinocyte_cna_like"),
    n_keratinocyte_normal_like = sum(target$Auto_keratinocyte_cna_profile_class == "keratinocyte_normal_like"),
    n_keratinocyte_indeterminate = sum(target$Auto_keratinocyte_cna_profile_class == "keratinocyte_indeterminate"),
    tumour_profile_threshold = tumour_threshold, margin_threshold = margin_threshold,
    n_malignant_level_1 = sum(target$Auto_malignancy == "malignant_level_1"),
    n_malignant_level_2 = sum(target$Auto_malignancy == "malignant_level_2"),
    n_malignant = sum(target$Auto_malignant), pct_malignant_epithelial = 100 * mean(target$Auto_malignant),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(malignancy_dir, "tables", "Auto_visiumhd_binned_infercna_malignancy_summary.csv")
existing <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
for (column in setdiff(colnames(summary_df), c("sample", "mode"))) existing[[column]] <- summary_df[[column]][match(existing$sample, summary_df$sample)]
write.csv(existing, summary_path, row.names = FALSE)
write.csv(summary_df, file.path(malignancy_dir, "tables", "Auto_visiumhd_binned_keratinocyte_cna_profile_summary.csv"), row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries"); dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(existing, file.path(summary_dir, "visiumhd_binned_infercna_malignancy_summary.csv"), row.names = FALSE)
writeLines(c(paste0("samples=", paste(args$samples, collapse = ",")), paste0("end=", format(Sys.time(), tz = "Europe/London"))), file.path(malignancy_dir, "logs", "Auto_visiumhd_binned_keratinocyte_cna_profile_run_summary.txt"))
