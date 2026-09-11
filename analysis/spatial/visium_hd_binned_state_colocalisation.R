#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: terminal
#   Script: analysis/spatial/visium_hd_binned_state_colocalisation.R
#   Description: State abundance and six-nearest-neighbour colocalisation for
#     malignant level-1/level-2 epithelial Visium HD 16 um bins.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_state_mapping_methodology.md
#   Inputs:
#     ref_outs/visium_hd_outs/state_mapping/tables/
#       Auto_visiumhd_binned_malignant_state_annotations.csv.gz
#   Outputs:
#     tables/: per-bin colocalisation, sample/state and pooled summaries
#     figures/: abundance, raw colocalisation, and multi-page audit report
#     updates/new_updates/summaries/: compact colocalisation summary
#   Cache/replot: inexpensive relative to state scoring; always recomputed from
#     the live mapped state table.
#   Run: Rscript analysis/spatial/visium_hd_binned_state_colocalisation.R
#   Environment: dmtcp
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
  library(ggplot2)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "state_mapping")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
log_dir <- file.path(output_dir, "logs")
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
for (path in c(table_dir, figure_dir, log_dir, summary_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

input_path <- file.path(table_dir, "Auto_visiumhd_binned_malignant_state_annotations.csv.gz")
if (!file.exists(input_path)) stop("Missing mapped state table: ", input_path)
mapped <- fread(input_path)
required <- c("barcode", "sample", "Auto_state_B", "pxl_row_in_fullres", "pxl_col_in_fullres", "Auto_malignancy")
missing <- setdiff(required, names(mapped))
if (length(missing)) stop("Mapped state table lacks: ", paste(missing, collapse = ", "))

samples <- c("SUR1231", "FFPEA1", "FFPED1")
primary_states <- c(
  "Classic proliferation", "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia", "Stress adaptive",
  "Cancer-cell immune mimicry"
)
state_order <- c(primary_states, "Unresolved", "Hybrid")
state_colours <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved" = "#9E9E9E", "Hybrid" = "#111111"
)
knn_k <- 6L
jitter_height <- 0.4 / knn_k

####################
# Match the legacy Visium HD colocalisation presentation without importing its
# outdated state constants or its publication helper side effects.
legacy_colocalisation_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "#111827", family = ""),
      axis.text = element_text(size = rel(0.85), colour = "#111827"),
      axis.title = element_text(size = rel(1), face = "bold"),
      panel.border = element_rect(colour = "#CBD5E1", fill = NA, linewidth = 0.4),
      axis.line = element_line(colour = "#111827", linewidth = 0.4),
      plot.margin = margin(8, 8, 8, 8)
    )
}
####################

if (!all(mapped$sample %in% samples)) stop("Unexpected sample in mapped state table")
if (!all(mapped$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2"))) {
  stop("State table contains bins outside malignant levels 1 and 2")
}

abundance <- mapped[, .(n_bins = .N), by = .(sample, state = Auto_state_B)]
abundance[, pct_malignant_bins := 100 * n_bins / sum(n_bins), by = sample]
abundance[, state := factor(state, levels = state_order)]
setorder(abundance, sample, state)
fwrite(abundance, file.path(table_dir, "Auto_visiumhd_binned_malignant_state_abundance_for_colocalisation.csv"))

colocalisation_rows <- list()
for (sample_name in samples) {
  sample_data <- mapped[
    sample == sample_name & Auto_state_B %in% primary_states &
      is.finite(pxl_row_in_fullres) & is.finite(pxl_col_in_fullres)
  ]
  if (nrow(sample_data) < 8L) stop("Too few biological-state bins for ", sample_name)
  effective_k <- min(knn_k, nrow(sample_data) - 1L)
  coordinates <- as.matrix(sample_data[, .(pxl_row_in_fullres, pxl_col_in_fullres)])
  neighbours <- FNN::get.knn(coordinates, k = effective_k)
  state_values <- as.character(sample_data$Auto_state_B)
  same_score <- vapply(seq_len(nrow(sample_data)), function(index) {
    mean(state_values[neighbours$nn.index[index, ]] == state_values[[index]])
  }, numeric(1))
  baseline <- table(state_values) / length(state_values)
  baseline_score <- as.numeric(baseline[state_values])
  colocalisation_rows[[sample_name]] <- data.table(
    barcode = sample_data$barcode,
    sample = sample_name,
    state = state_values,
    k = effective_k,
    same_neighbor_score = same_score,
    sample_state_fraction = baseline_score,
    same_neighbor_excess = same_score - baseline_score,
    same_neighbor_ratio = fifelse(baseline_score > 0, same_score / baseline_score, NA_real_),
    mean_neighbor_distance = rowMeans(neighbours$nn.dist),
    maximum_neighbor_distance = apply(neighbours$nn.dist, 1L, max)
  )
}
colocalisation <- rbindlist(colocalisation_rows, use.names = TRUE)
colocalisation[, state := factor(state, levels = primary_states)]
fwrite(
  colocalisation,
  file.path(table_dir, "Auto_visiumhd_binned_malignant_state_colocalisation_per_bin.csv.gz")
)

sample_summary <- colocalisation[, .(
  n_bins = .N,
  mean_same_neighbor_score = mean(same_neighbor_score),
  median_same_neighbor_score = median(same_neighbor_score),
  q25_same_neighbor_score = quantile(same_neighbor_score, 0.25),
  q75_same_neighbor_score = quantile(same_neighbor_score, 0.75),
  mean_sample_state_fraction = mean(sample_state_fraction),
  mean_same_neighbor_excess = mean(same_neighbor_excess),
  median_neighbor_distance = median(mean_neighbor_distance),
  q95_maximum_neighbor_distance = quantile(maximum_neighbor_distance, 0.95)
), by = .(sample, state)]
setorder(sample_summary, sample, state)
fwrite(sample_summary, file.path(table_dir, "Auto_visiumhd_binned_malignant_state_colocalisation_by_sample.csv"))
fwrite(sample_summary, file.path(summary_dir, "visium_hd_binned_malignant_state_colocalisation_summary.csv"))

pooled_summary <- colocalisation[, .(
  n_bins = .N,
  mean_same_neighbor_score = mean(same_neighbor_score),
  median_same_neighbor_score = median(same_neighbor_score),
  mean_same_neighbor_excess = mean(same_neighbor_excess)
), by = state]
setorder(pooled_summary, state)
fwrite(pooled_summary, file.path(table_dir, "Auto_visiumhd_binned_malignant_state_colocalisation_pooled.csv"))

abundance_plot <- ggplot(abundance, aes(x = sample, y = pct_malignant_bins, fill = state)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = state_colours, drop = FALSE) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Malignant epithelial bins", fill = "State") +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 11),
    legend.position = "right", legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

raw_colocalisation_plot <- ggplot(
  colocalisation, aes(x = state, y = same_neighbor_score, fill = state, colour = state)
) +
  geom_boxplot(
    width = 0.5, outlier.shape = NA, alpha = 0.8,
    linewidth = 0.6, colour = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = jitter_height),
    size = 0.5, alpha = 0.05, colour = "black"
  ) +
  scale_fill_manual(values = state_colours, guide = "none", drop = FALSE) +
  scale_colour_manual(values = state_colours, guide = "none", drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1), limits = c(0, 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = NULL, y = "Same-state neighbours") +
  legacy_colocalisation_theme(14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )

sample_colocalisation_plot <- raw_colocalisation_plot +
  facet_wrap(~sample, nrow = 1) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 9))

####################
excess_colocalisation_plot <- ggplot(
  colocalisation,
  aes(x = state, y = same_neighbor_excess, fill = state, colour = state)
) +
  geom_hline(
    yintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey30"
  ) +
  geom_boxplot(
    width = 0.5, outlier.shape = NA, alpha = 0.8,
    linewidth = 0.6, colour = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = jitter_height),
    size = 0.5, alpha = 0.05, colour = "black"
  ) +
  scale_fill_manual(values = state_colours, guide = "none", drop = FALSE) +
  scale_colour_manual(values = state_colours, guide = "none", drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Excess same-state neighbours") +
  legacy_colocalisation_theme(14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )
####################

ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_abundance.pdf"), abundance_plot, width = 10, height = 7)
ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_abundance.png"), abundance_plot, width = 10, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_colocalisation.pdf"), raw_colocalisation_plot, width = 7.5, height = 6)
ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_colocalisation.png"), raw_colocalisation_plot, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_colocalisation_abundance_adjusted.pdf"), excess_colocalisation_plot, width = 7.5, height = 6)
ggsave(file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_colocalisation_abundance_adjusted.png"), excess_colocalisation_plot, width = 7.5, height = 6, dpi = 300)

report_path <- file.path(figure_dir, "Auto_visiumhd_binned_malignant_state_spatial_summary.pdf")
pdf(report_path, width = 16, height = 9, onefile = TRUE, useDingbats = FALSE)
print(abundance_plot)
print(raw_colocalisation_plot)
print(sample_colocalisation_plot)
print(excess_colocalisation_plot)
dev.off()

writeLines(
  c(
    paste0("status=complete"),
    paste0("samples=", paste(samples, collapse = ";")),
    paste0("input=", input_path),
    paste0("primary_states=", paste(primary_states, collapse = ";")),
    paste0("knn_k=", knn_k),
    "neighbour_pool=primary biological states only; Hybrid and Unresolved excluded before kNN",
    "distance_limit=none; exact legacy six-nearest-neighbour definition",
    "same_neighbor_excess=same_neighbor_score minus within-sample state fraction; included as report page and standalone figure",
    paste0("n_colocalisation_bins=", nrow(colocalisation))
  ),
  file.path(log_dir, "Auto_visiumhd_binned_state_colocalisation_run_summary.txt")
)
####################
