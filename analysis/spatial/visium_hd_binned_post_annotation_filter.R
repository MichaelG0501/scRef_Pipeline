#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_binned_post_annotation_filter.R
#   Description: permissive marker-expression and unrelated-lineage
#     coexpression filtering of final 16 um Visium HD annotations.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_filter_malignancy_methodology.md
#   Inputs:
#     analysis/shared/visium_hd_celltype_colours.tsv
#     ref_outs/visium_hd_outs/tables/Auto_<sample>_binned_cell_annotations.csv.gz
#   Outputs:
#     intermediate/: none
#     tables/: filtered per-bin annotations, before/after counts, filter summary,
#       and parameter table under ref_outs/visium_hd_outs/post_annotation_filter/
#     figures/: landscape before/after count PDF, segmented-versus-binned
#       comparison PDFs before and after filtering, per-sample count PNG files,
#       and post-filter spatial/UMAP diagnostics under
#       figures/annotation_diagnostics/after_filtering/
#     logs/: run summary
#     updates/new_updates/summaries/: compact filter summary
#   Cache/replot: inexpensive; always recomputed from the annotation tables.
#   Run: Rscript analysis/spatial/visium_hd_binned_post_annotation_filter.R
#   Environment: dmtcp
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)

samples <- c("SUR1231", "FFPEA1", "FFPED1")
annotation_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "tables")
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "post_annotation_filter")
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
after_filter_figure_dir <- file.path(
  wd, "ref_outs", "visium_hd_outs", "figures",
  "annotation_diagnostics", "after_filtering"
)
for (subdir in c("tables", "figures", "logs")) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(after_filter_figure_dir, recursive = TRUE, showWarnings = FALSE)

# The annotation score is the mean log1p(CP10K) expression of the fixed marker
# panel. A positive home score therefore means at least one home-panel marker
# was detected. Structural epithelial/fibroblast signal is ignored as an
# off-lineage competitor because it is the expected dominant spillover source.
home_score_min <- 0
competitor_score_min <- 1
competitor_to_home_min <- 1
structural_types <- c("epithelial", "fibroblast")
related_types <- list(
  "t.cell" = c("nk.cell"),
  "nk.cell" = c("t.cell"),
  "b.cell" = c("plasma"),
  "plasma" = c("b.cell"),
  "macrophage" = c("dendritic"),
  "dendritic" = c("macrophage"),
  "keratinocyte" = c("epithelial"),
  "epithelial" = c("keratinocyte")
)

celltype_order <- c(
  "epithelial", "keratinocyte", "fibroblast", "endothelial", "macrophage",
  "dendritic", "mast", "neutrophil", "t.cell", "nk.cell", "b.cell",
  "plasma", "lymph", "erythrocyte", "unresolved"
)
celltype_colour_table <- data.frame(
  celltype = c("epithelial", "fibroblast", "endothelial", "macrophage", "mast", "t.cell", "b.cell", "nk.cell", "plasma", "dendritic", "lymph", "erythrocyte", "keratinocyte", "neutrophil", "unresolved", "combined"),
  colour = c("#D73027", "#8C564B", "#1F78B4", "#FF7F00", "#A65628", "#33A02C", "#377EB8", "#984EA3", "#E377C2", "#17BECF", "#6BAED6", "#7F7F7F", "#E6AB02", "#1B9E77", "#BDBDBD", "#555555"),
  stringsAsFactors = FALSE
)
celltype_colours <- setNames(celltype_colour_table$colour, celltype_colour_table$celltype)
combined_celltype_colour <- unname(celltype_colours[["combined"]])
celltype_colours <- celltype_colours[names(celltype_colours) != "combined"]
celltype_plot_colours <- function(labels) {
  labels <- unique(as.character(labels))
  colours <- unname(celltype_colours[labels])
  colours[is.na(colours) | grepl("\\|", labels)] <- combined_celltype_colour
  setNames(colours, labels)
}

####################
# Apply the exact production filter to segmented annotations for the count
# comparison only. Segmented results are not written as malignancy inputs.
calculate_comparison_counts <- function(annotation, sample_name, method_name) {
  annotation <- copy(annotation)
  annotation[, Auto_prefilter_celltype := as.character(Auto_annotation_celltype)]
  annotation[is.na(Auto_prefilter_celltype) | Auto_prefilter_celltype == "", Auto_prefilter_celltype := "unresolved"]
  annotation[, Auto_postfilter_celltype := sub("\\|.*$", "", as.character(Auto_annotation_celltype))]
  annotation[is.na(Auto_postfilter_celltype) | Auto_postfilter_celltype == "", Auto_postfilter_celltype := "unresolved"]

  score_columns <- intersect(paste0(setdiff(celltype_order, "unresolved"), "_score"), names(annotation))
  score_types <- sub("_score$", "", score_columns)
  names(score_columns) <- score_types
  missing_home <- setdiff(unique(annotation$Auto_postfilter_celltype), c(score_types, "unresolved"))
  if (length(missing_home)) stop("No marker-score column for: ", paste(missing_home, collapse = ", "))

  annotation[, `:=`(
    Auto_postfilter_home_score = NA_real_,
    Auto_postfilter_competing_score = NA_real_
  )]
  for (home_type in setdiff(unique(annotation$Auto_postfilter_celltype), "unresolved")) {
    idx <- which(annotation$Auto_postfilter_celltype == home_type)
    annotation$Auto_postfilter_home_score[idx] <- as.numeric(annotation[[score_columns[[home_type]]]][idx])
    excluded <- unique(c(home_type, structural_types, related_types[[home_type]]))
    candidate_types <- setdiff(score_types, excluded)
    if (!length(candidate_types)) next
    candidate_columns <- unname(score_columns[candidate_types])
    candidate_matrix <- as.matrix(annotation[idx, ..candidate_columns])
    candidate_matrix[!is.finite(candidate_matrix)] <- -Inf
    top_score <- candidate_matrix[cbind(seq_len(nrow(candidate_matrix)), max.col(candidate_matrix, ties.method = "first"))]
    top_score[!is.finite(top_score)] <- NA_real_
    annotation$Auto_postfilter_competing_score[idx] <- top_score
  }

  annotation[, Auto_postfilter_expression_pass :=
    Auto_postfilter_celltype != "unresolved" &
      is.finite(Auto_postfilter_home_score) &
      Auto_postfilter_home_score > home_score_min]
  annotation[, Auto_postfilter_coexpression_conflict :=
    Auto_postfilter_expression_pass &
      is.finite(Auto_postfilter_competing_score) &
      Auto_postfilter_competing_score >= competitor_score_min &
      Auto_postfilter_competing_score >= competitor_to_home_min * Auto_postfilter_home_score]
  annotation[, Auto_postfilter_keep :=
    Auto_postfilter_expression_pass & !Auto_postfilter_coexpression_conflict]

  before_counts <- annotation[, .(count = .N), by = .(celltype = Auto_prefilter_celltype)]
  before_counts[, stage := "Before filtering"]
  after_counts <- annotation[Auto_postfilter_keep %in% TRUE, .(count = .N), by = .(celltype = Auto_postfilter_celltype)]
  after_counts[, stage := "After filtering"]
  result <- rbind(before_counts, after_counts, fill = TRUE)
  result[, `:=`(sample = sample_name, method = method_name)]
  result
}
####################

count_rows <- list()
filter_rows <- list()
run_log <- c(
  paste0("start=", format(Sys.time(), tz = "Europe/London")),
  paste0("samples=", paste(samples, collapse = ";")),
  paste0("home_score_min=", home_score_min),
  paste0("competitor_score_min=", competitor_score_min),
  paste0("competitor_to_home_min=", competitor_to_home_min),
  paste0("ignored_structural_competitors=", paste(structural_types, collapse = ";"))
)

pdf_path <- file.path(output_dir, "figures", "Auto_visium_hd_binned_celltype_counts_before_after_filtering.pdf")
pdf(pdf_path, width = 15, height = 8.5, onefile = TRUE, useDingbats = FALSE)

for (sample_name in samples) {
  annotation_path <- file.path(annotation_dir, paste0("Auto_", sample_name, "_binned_cell_annotations.csv.gz"))
  if (!file.exists(annotation_path)) stop("Missing annotation table: ", annotation_path)
  annotation <- data.table::fread(annotation_path)
  needed <- c("barcode", "Auto_annotation_celltype")
  missing <- setdiff(needed, names(annotation))
  if (length(missing)) stop("Annotation table missing columns: ", paste(missing, collapse = ", "))

  annotation[, Auto_prefilter_celltype := as.character(Auto_annotation_celltype)]
  annotation[is.na(Auto_prefilter_celltype) | Auto_prefilter_celltype == "", Auto_prefilter_celltype := "unresolved"]
  annotation[, Auto_postfilter_celltype := sub("\\|.*$", "", as.character(Auto_annotation_celltype))]
  annotation[is.na(Auto_postfilter_celltype) | Auto_postfilter_celltype == "", Auto_postfilter_celltype := "unresolved"]

  score_columns <- intersect(paste0(setdiff(celltype_order, "unresolved"), "_score"), names(annotation))
  score_types <- sub("_score$", "", score_columns)
  names(score_columns) <- score_types
  missing_home <- setdiff(unique(annotation$Auto_postfilter_celltype), c(score_types, "unresolved"))
  if (length(missing_home)) stop("No marker-score column for: ", paste(missing_home, collapse = ", "))

  annotation[, `:=`(
    Auto_postfilter_home_score = NA_real_,
    Auto_postfilter_competing_celltype = NA_character_,
    Auto_postfilter_competing_score = NA_real_
  )]

  for (home_type in setdiff(unique(annotation$Auto_postfilter_celltype), "unresolved")) {
    idx <- which(annotation$Auto_postfilter_celltype == home_type)
    annotation$Auto_postfilter_home_score[idx] <- as.numeric(annotation[[score_columns[[home_type]]]][idx])
    excluded <- unique(c(home_type, structural_types, related_types[[home_type]]))
    candidate_types <- setdiff(score_types, excluded)
    if (!length(candidate_types)) next
    candidate_columns <- unname(score_columns[candidate_types])
    candidate_matrix <- as.matrix(annotation[idx, ..candidate_columns])
    candidate_matrix[!is.finite(candidate_matrix)] <- -Inf
    top_index <- max.col(candidate_matrix, ties.method = "first")
    top_score <- candidate_matrix[cbind(seq_len(nrow(candidate_matrix)), top_index)]
    top_score[!is.finite(top_score)] <- NA_real_
    annotation$Auto_postfilter_competing_celltype[idx] <- candidate_types[top_index]
    annotation$Auto_postfilter_competing_score[idx] <- top_score
  }

  annotation[, Auto_postfilter_expression_pass :=
    Auto_postfilter_celltype != "unresolved" &
      is.finite(Auto_postfilter_home_score) &
      Auto_postfilter_home_score > home_score_min]
  annotation[, Auto_postfilter_coexpression_conflict :=
    Auto_postfilter_expression_pass &
      is.finite(Auto_postfilter_competing_score) &
      Auto_postfilter_competing_score >= competitor_score_min &
      Auto_postfilter_competing_score >= competitor_to_home_min * Auto_postfilter_home_score]
  annotation[, Auto_postfilter_keep :=
    Auto_postfilter_expression_pass & !Auto_postfilter_coexpression_conflict]
  annotation[, Auto_postfilter_reason := fifelse(
    Auto_postfilter_celltype == "unresolved", "unresolved_annotation",
    fifelse(!Auto_postfilter_expression_pass, "no_home_marker_expression",
      fifelse(Auto_postfilter_coexpression_conflict,
        "strong_unrelated_nonstructural_coexpression", "pass")
    )
  )]
  annotation[, `:=`(
    Auto_postfilter_home_score_min = home_score_min,
    Auto_postfilter_competitor_score_min = competitor_score_min,
    Auto_postfilter_competitor_to_home_min = competitor_to_home_min,
    Auto_postfilter_structural_competitors_ignored = paste(structural_types, collapse = ";")
  )]

  output_path <- file.path(output_dir, "tables", paste0("Auto_", sample_name, "_binned_filtered_annotations.csv.gz"))
  data.table::fwrite(annotation, output_path)

  before_counts <- annotation[, .(count = .N), by = Auto_prefilter_celltype]
  before_counts[, stage := "Before filtering"]
  after_counts <- annotation[Auto_postfilter_keep == TRUE, .(count = .N), by = Auto_postfilter_celltype]
  after_counts[, stage := "After filtering"]
  setnames(before_counts, "Auto_prefilter_celltype", "celltype")
  setnames(after_counts, "Auto_postfilter_celltype", "celltype")
  sample_counts <- rbind(before_counts, after_counts, fill = TRUE)
  all_types <- unique(sample_counts$celltype)
  count_grid <- CJ(celltype = all_types, stage = c("Before filtering", "After filtering"), unique = TRUE)
  sample_counts <- merge(count_grid, sample_counts, by = c("celltype", "stage"), all.x = TRUE)
  sample_counts[is.na(count), count := 0L]
  sample_counts[, sample := sample_name]
  count_rows[[sample_name]] <- copy(sample_counts)

  reason_summary <- annotation[, .(n_bins = .N), by = .(Auto_postfilter_celltype, Auto_postfilter_reason)]
  reason_summary[, sample := sample_name]
  filter_rows[[sample_name]] <- reason_summary

  plot_order <- before_counts[order(-count), celltype]
  sample_counts[, celltype_plot := factor(celltype, levels = plot_order)]
  sample_counts[, stage := factor(stage, levels = c("Before filtering", "After filtering"))]
  max_count <- max(sample_counts$count)
  p <- ggplot(sample_counts, aes(x = celltype_plot, y = count, fill = celltype)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = scales::comma(count)), vjust = -0.35, size = 4.0) +
    facet_grid(. ~ stage, scales = "fixed", space = "fixed") +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = celltype_plot_colours(sample_counts$celltype), drop = FALSE) +
    scale_y_continuous(
      limits = c(0, max_count * 1.18),
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(title = sample_name, x = NULL, y = "Bin count") +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(12, 34, 12, 12)
    )
  print(p)
  ggsave(
    file.path(output_dir, "figures", paste0("Auto_", sample_name, "_binned_celltype_counts_before_after_filtering.png")),
    p, width = 15, height = 8.5, dpi = 300
  )

  ####################
  # Preserve the canonical pre-filter diagnostics and create a separate,
  # directly comparable view after expression/coexpression filtering. Grey
  # observations show where bins were removed; retained bins use final labels.
  diagnostic_columns <- c(
    "pxl_col_in_fullres", "pxl_row_in_fullres", "UMAP_1", "UMAP_2"
  )
  missing_diagnostic_columns <- setdiff(diagnostic_columns, names(annotation))
  if (length(missing_diagnostic_columns)) {
    stop(
      "Cannot create post-filter diagnostics; missing columns: ",
      paste(missing_diagnostic_columns, collapse = ", ")
    )
  }
  retained <- annotation[
    Auto_postfilter_keep %in% TRUE &
      is.finite(pxl_col_in_fullres) & is.finite(pxl_row_in_fullres) &
      is.finite(UMAP_1) & is.finite(UMAP_2)
  ]
  removed <- annotation[
    !Auto_postfilter_keep %in% TRUE &
      is.finite(pxl_col_in_fullres) & is.finite(pxl_row_in_fullres) &
      is.finite(UMAP_1) & is.finite(UMAP_2)
  ]
  retained[, display_celltype := factor(
    Auto_postfilter_celltype,
    levels = intersect(celltype_order, unique(Auto_postfilter_celltype))
  )]
  retained_counts <- retained[, .N, by = Auto_postfilter_celltype]
  retained_count_labels <- setNames(
    paste0(retained_counts$Auto_postfilter_celltype, " (", comma(retained_counts$N), ")"),
    retained_counts$Auto_postfilter_celltype
  )
  retained_colours <- celltype_colours[levels(retained$display_celltype)]

  p_spatial <- ggplot() +
    geom_point(
      data = removed,
      aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres),
      colour = "#D9D9D9", size = 0.7, alpha = 0.45
    ) +
    geom_point(
      data = retained,
      aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, colour = display_celltype),
      size = 1.35, alpha = 0.88
    ) +
    scale_colour_manual(
      values = retained_colours,
      labels = retained_count_labels,
      drop = FALSE
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = paste(sample_name, "spatial"), x = NULL, y = NULL, colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    theme_void(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(0.7, "cm"),
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5)
    )

  p_umap <- ggplot() +
    geom_point(
      data = removed, aes(x = UMAP_1, y = UMAP_2),
      colour = "#D9D9D9", size = 0.7, alpha = 0.45
    ) +
    geom_point(
      data = retained, aes(x = UMAP_1, y = UMAP_2, colour = display_celltype),
      size = 1.35, alpha = 0.88
    ) +
    scale_colour_manual(
      values = retained_colours,
      labels = retained_count_labels,
      drop = FALSE
    ) +
    coord_equal() +
    labs(title = paste(sample_name, "UMAP"), x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(0.7, "cm"),
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5)
    )

  diagnostic_plot <- gridExtra::arrangeGrob(
    p_spatial,
    p_umap,
    ncol = 2,
    widths = c(1, 1)
  )
  diagnostic_stem <- file.path(
    after_filter_figure_dir,
    paste0("Auto_", sample_name, "_binned_annotation_diagnostics_after_filtering")
  )
  ggsave(
    paste0(diagnostic_stem, ".pdf"), diagnostic_plot,
    width = 18, height = 10.5, device = cairo_pdf
  )
  ggsave(
    paste0(diagnostic_stem, ".png"), diagnostic_plot,
    width = 18, height = 10.5, dpi = 300
  )
  ####################

  run_log <- c(run_log, paste0(
    "sample=", sample_name,
    "; n_before=", nrow(annotation),
    "; n_after=", sum(annotation$Auto_postfilter_keep),
    "; pct_after=", round(100 * mean(annotation$Auto_postfilter_keep), 2),
    "; diagnostic_pdf=", paste0(diagnostic_stem, ".pdf")
  ))
  rm(annotation)
  gc()
}
dev.off()

all_counts <- rbindlist(count_rows, use.names = TRUE, fill = TRUE)
setcolorder(all_counts, c("sample", "stage", "celltype", "count"))
fwrite(all_counts, file.path(output_dir, "tables", "Auto_visium_hd_binned_celltype_counts_before_after_filtering.csv"))

####################
# Compare segmented cells (left) with 16 um bins (right) using a shared count
# scale within each sample. Separate PDFs show the unfiltered and identically
# filtered annotations.
comparison_rows <- list()
binned_comparison_counts <- copy(all_counts)
binned_comparison_counts[, method := "16 um bins"]
comparison_rows[["binned"]] <- binned_comparison_counts

for (sample_name in samples) {
  segmented_path <- file.path(
    annotation_dir,
    paste0("Auto_", sample_name, "_segmented_cell_annotations.csv.gz")
  )
  if (!file.exists(segmented_path)) stop("Missing segmented annotation table: ", segmented_path)
  segmented_annotation <- fread(segmented_path)
  needed_segmented <- c("barcode", "Auto_annotation_celltype")
  missing_segmented <- setdiff(needed_segmented, names(segmented_annotation))
  if (length(missing_segmented)) {
    stop("Segmented annotation table missing columns: ", paste(missing_segmented, collapse = ", "))
  }
  comparison_rows[[paste0(sample_name, "_segmented")]] <- calculate_comparison_counts(
    segmented_annotation, sample_name, "Cell segmentation"
  )
  rm(segmented_annotation)
  gc()
}

comparison_counts <- rbindlist(comparison_rows, use.names = TRUE, fill = TRUE)
comparison_counts[, method := factor(method, levels = c("Cell segmentation", "16 um bins"))]
comparison_counts[, stage := factor(stage, levels = c("Before filtering", "After filtering"))]
comparison_grid <- CJ(
  sample = samples,
  stage = levels(comparison_counts$stage),
  method = levels(comparison_counts$method),
  celltype = unique(comparison_counts$celltype),
  unique = TRUE
)
comparison_grid[, `:=`(
  stage = factor(stage, levels = levels(comparison_counts$stage)),
  method = factor(method, levels = levels(comparison_counts$method))
)]
comparison_counts <- merge(
  comparison_grid,
  comparison_counts,
  by = c("sample", "stage", "method", "celltype"),
  all.x = TRUE
)
comparison_counts[is.na(count), count := 0L]
fwrite(
  comparison_counts,
  file.path(output_dir, "tables", "Auto_visium_hd_celltype_counts_segmented_vs_binned_before_after_filtering.csv")
)

for (stage_name in levels(comparison_counts$stage)) {
  stage_slug <- if (stage_name == "Before filtering") "before_filtering" else "after_filtering"
  comparison_pdf <- file.path(
    output_dir, "figures",
    paste0("Auto_visium_hd_celltype_counts_segmented_vs_binned_", stage_slug, ".pdf")
  )
  pdf(comparison_pdf, width = 16, height = 9, onefile = TRUE, useDingbats = FALSE)
  for (sample_name in samples) {
    plot_counts <- comparison_counts[
      sample == sample_name & as.character(stage) == stage_name
    ]
    present_celltypes <- plot_counts[, .(maximum_count = max(count)), by = celltype][
      maximum_count > 0, celltype
    ]
    plot_counts <- plot_counts[celltype %in% present_celltypes]
    plot_order <- plot_counts[, .(maximum_count = max(count)), by = celltype][
      order(-maximum_count), celltype
    ]
    plot_counts[, celltype_plot := factor(celltype, levels = plot_order)]
    shared_maximum <- max(plot_counts$count)
    comparison_plot <- ggplot(
      plot_counts,
      aes(x = celltype_plot, y = count, fill = celltype)
    ) +
      geom_col(width = 0.72) +
      geom_text(aes(label = comma(count)), vjust = -0.35, size = 3.9) +
      facet_grid(. ~ method, scales = "fixed", space = "fixed") +
      coord_cartesian(clip = "off") +
      scale_fill_manual(values = celltype_plot_colours(plot_counts$celltype), drop = FALSE) +
      scale_y_continuous(
        limits = c(0, max(1, shared_maximum * 1.18)),
        labels = comma,
        expand = expansion(mult = c(0, 0.01))
      ) +
      labs(title = sample_name, x = NULL, y = "Observation count") +
      theme_classic(base_size = 15) +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.margin = margin(12, 36, 12, 12)
      )
    print(comparison_plot)
  }
  dev.off()
  run_log <- c(run_log, paste0("comparison_pdf_", stage_slug, "=", comparison_pdf))
}
####################

all_reasons <- rbindlist(filter_rows, use.names = TRUE, fill = TRUE)
setnames(all_reasons, c("Auto_postfilter_celltype", "Auto_postfilter_reason"), c("celltype", "filter_reason"))
setcolorder(all_reasons, c("sample", "celltype", "filter_reason", "n_bins"))
fwrite(all_reasons, file.path(output_dir, "tables", "Auto_visium_hd_binned_post_annotation_filter_summary.csv"))
fwrite(all_reasons, file.path(summary_dir, "visium_hd_binned_post_annotation_filter_summary.csv"))

parameters <- data.table(
  parameter = c(
    "samples", "home_score_min_strictly_greater_than", "competitor_score_min",
    "competitor_to_home_min", "ignored_structural_competitors", "related_pairs_ignored"
  ),
  value = c(
    paste(samples, collapse = ";"), as.character(home_score_min),
    as.character(competitor_score_min), as.character(competitor_to_home_min),
    paste(structural_types, collapse = ";"),
    paste(vapply(names(related_types), function(x) paste0(x, ":", paste(related_types[[x]], collapse = ",")), character(1)), collapse = ";")
  )
)
fwrite(parameters, file.path(output_dir, "tables", "Auto_visium_hd_binned_post_annotation_filter_parameters.csv"))

run_log <- c(run_log, paste0("end=", format(Sys.time(), tz = "Europe/London")))
writeLines(run_log, file.path(output_dir, "logs", "Auto_visium_hd_binned_post_annotation_filter_run_summary.txt"))
####################
