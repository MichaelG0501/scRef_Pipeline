#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_binned_malignancy.R
#   Description: binned-only Visium HD epithelial malignancy classification
#     using InferCNA followed by scATLAS cancer-signature rescue.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_filter_malignancy_methodology.md
#   Inputs:
#     analysis/spatial/visium_hd_samples.tsv
#     ref_outs/visium_hd_outs/post_annotation_filter/tables/
#       Auto_<sample>_binned_filtered_annotations.csv.gz
#     ref_outs/cancer_signatures.txt
#     Space Ranger 16 um filtered_feature_bc_matrix.h5 files in the manifest
#     /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/
#       hg38_gencode_v27.txt
#   Outputs:
#     intermediate/: InferCNA matrix caches under the corresponding ephemeral
#       scRef_Pipeline path only
#     tables/: full per-bin malignancy tables, selected references/signature
#       genes, parameters, and sample summaries under
#       ref_outs/visium_hd_outs/malignancy/
#     figures/: InferCNA scatter and spatial malignancy diagnostics, including
#       a standalone spatial PDF and PNG for each sample, plus a three-sample
#       horizontal stacked malignancy-count plot
#     logs/: run summary and session information
#     updates/new_updates/summaries/: compact malignancy summary
#   Cache/replot: reuses an InferCNA cache only after exact target/reference
#     barcode validation. Set SCREF_FORCE_REBUILD=TRUE to ignore caches or
#     SCREF_REPLOT_ONLY=TRUE to rebuild plots from live per-bin tables.
#   Run: Rscript analysis/spatial/visium_hd_binned_malignancy.R
#   Environment: dmtcp
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
  library(infercna)
  library(ggplot2)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)

samples <- c("SUR1231", "FFPEA1", "FFPED1")

filter_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "post_annotation_filter")
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "malignancy")
ephemeral_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/visium_hd_outs/malignancy/intermediate"
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
for (subdir in c("tables", "figures", "logs")) {
  dir.create(file.path(output_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}
dir.create(ephemeral_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

normal_types <- c("endothelial", "macrophage", "fibroblast")
min_reference_cells <- 15L
min_epithelial_cells <- 30L
cna_sd_k <- 1
cancer_signature_threshold <- 1
cancer_signature_top_n <- 50L
cancer_signature_path <- file.path(wd, "ref_outs", "cancer_signatures.txt")
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
force_rebuild <- tolower(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")) %in% c("true", "t", "1", "yes")
replot_only <- tolower(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")) %in% c("true", "t", "1", "yes")

if (!file.exists(cancer_signature_path)) stop("Missing cancer signature: ", cancer_signature_path)
if (!file.exists(gene_order_path)) stop("Missing gene order: ", gene_order_path)
manifest <- data.table(
  sample = c("SUR1122", "SUR1231", "FFPEA1", "FFPED1"),
  binned_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"),
  segmented_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs")
)
manifest <- manifest[sample %in% samples]
if (!setequal(manifest$sample, samples)) stop("Manifest does not contain exactly the requested samples")
manifest[, sample_order := match(sample, samples)]
setorder(manifest, sample_order)
manifest[, sample_order := NULL]
cancer_signature_genes <- unique(scan(cancer_signature_path, what = character(), quiet = TRUE))
gene_order <- fread(gene_order_path, header = FALSE, fill = TRUE, select = 1L)
genome_genes <- unique(gene_order[[1L]])

read_10x_counts <- function(path) {
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1L]]
  }
  Matrix::Matrix(counts, sparse = TRUE)
}

to_cpm <- function(counts) {
  library_size <- Matrix::colSums(counts)
  library_size[!is.finite(library_size) | library_size <= 0] <- 1
  Matrix::t(Matrix::t(counts) * (1e6 / library_size))
}

make_scatter_plot <- function(cell_table, sample_name) {
  plot_data <- cell_table[is_reference %in% TRUE | is_epithelial_target %in% TRUE]
  plot_data[, plot_group := fifelse(
    is_reference %in% TRUE, "Reference",
    fifelse(Auto_malignancy == "malignant_level_1", "Malignant level 1",
      fifelse(Auto_malignancy == "malignant_level_2", "Malignant level 2",
        fifelse(Auto_cna_class == "cna_unresolved", "CNA unresolved", "CNA non-malignant")
      )
    )
  )]
  plot_data[, plot_group := factor(
    plot_group,
    levels = c("Reference", "CNA non-malignant", "CNA unresolved", "Malignant level 2", "Malignant level 1")
  )]
  signal_threshold <- unique(na.omit(plot_data$Auto_cna_signal_threshold))[1L]
  cor_threshold <- unique(na.omit(plot_data$Auto_cna_cor_threshold))[1L]
  ggplot(plot_data, aes(x = cna.signal, y = cna.cor, colour = plot_group)) +
    geom_point(size = 0.9, alpha = 0.68) +
    geom_vline(xintercept = signal_threshold, linetype = "dashed", linewidth = 0.6, colour = "grey25") +
    geom_hline(yintercept = cor_threshold, linetype = "dashed", linewidth = 0.6, colour = "grey25") +
    scale_colour_manual(values = c(
      "Reference" = "#969696", "CNA non-malignant" = "#4DAF4A",
      "CNA unresolved" = "#FDB863", "Malignant level 2" = "#984EA3",
      "Malignant level 1" = "#D73027"
    ), drop = FALSE) +
    scale_x_continuous(labels = scales::scientific) +
    labs(title = sample_name, x = "CNA signal", y = "CNA correlation", colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1), nrow = 2, byrow = TRUE)) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
}

make_spatial_plot <- function(cell_table, sample_name) {
  plot_data <- cell_table[
    Auto_postfilter_keep %in% TRUE &
      is.finite(pxl_col_in_fullres) & is.finite(pxl_row_in_fullres)
  ]
  plot_data[, plot_group := "Other retained cell type"]
  plot_data[is_epithelial_target %in% TRUE, plot_group := "Epithelial non-malignant/unresolved"]
  plot_data[Auto_malignancy == "malignant_level_2", plot_group := "Malignant level 2"]
  plot_data[Auto_malignancy == "malignant_level_1", plot_group := "Malignant level 1"]
  plot_data[, plot_group := factor(
    plot_group,
    levels = c("Other retained cell type", "Epithelial non-malignant/unresolved", "Malignant level 2", "Malignant level 1")
  )]
  setorder(plot_data, plot_group)
  ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, colour = plot_group)) +
    geom_point(size = 0.45, alpha = 0.72) +
    scale_colour_manual(values = c(
      "Other retained cell type" = "#D0D0D0",
      "Epithelial non-malignant/unresolved" = "#4DAF4A",
      "Malignant level 2" = "#984EA3", "Malignant level 1" = "#D73027"
    ), drop = FALSE) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = sample_name, x = NULL, y = NULL, colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1), nrow = 2, byrow = TRUE)) +
    theme_void(base_size = 15) +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
}

save_diagnostic_plots <- function(cell_table, sample_name) {
  scatter_plot <- make_scatter_plot(cell_table, sample_name)
  spatial_plot <- make_spatial_plot(cell_table, sample_name)
  pdf_path <- file.path(output_dir, "figures", paste0("Auto_", sample_name, "_binned_malignancy_diagnostics.pdf"))
  pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE, useDingbats = FALSE)
  print(scatter_plot)
  print(spatial_plot)
  dev.off()
  ggsave(
    file.path(output_dir, "figures", paste0("Auto_", sample_name, "_binned_infercna_scatter.png")),
    scatter_plot, width = 10, height = 7.5, dpi = 300
  )
  ggsave(
    file.path(output_dir, "figures", paste0("Auto_", sample_name, "_binned_malignancy_spatial.png")),
    spatial_plot, width = 11, height = 8.5, dpi = 300
  )
  ggsave(
    file.path(output_dir, "figures", paste0("Auto_", sample_name, "_binned_malignancy_spatial.pdf")),
    spatial_plot, width = 11, height = 8.5, device = cairo_pdf
  )
}

####################
# Summarise every post-filter epithelial target so the stacked counts reconcile
# exactly to the per-bin malignancy tables. Keep unresolved separate rather than
# treating uncertain CNA evidence as non-malignant.
save_malignancy_count_plot <- function() {
  malignancy_levels <- c(
    "non_malignant", "unresolved", "malignant_level_1", "malignant_level_2"
  )
  malignancy_labels <- c(
    "non_malignant" = "Non-malignant",
    "unresolved" = "Unresolved",
    "malignant_level_1" = "Malignant level 1",
    "malignant_level_2" = "Malignant level 2"
  )
  malignancy_colours <- c(
    "Non-malignant" = "#4DAF4A",
    "Unresolved" = "#FDB863",
    "Malignant level 1" = "#D73027",
    "Malignant level 2" = "#984EA3"
  )

  count_rows <- lapply(samples, function(sample_name) {
    table_path <- file.path(
      output_dir, "tables", paste0("Auto_", sample_name, "_binned_malignancy.csv.gz")
    )
    if (!file.exists(table_path)) stop("Missing malignancy table for count plot: ", table_path)
    cell_table <- fread(table_path, select = c("is_epithelial_target", "Auto_malignancy"))
    observed <- cell_table[is_epithelial_target %in% TRUE, .N, by = Auto_malignancy]
    unknown <- setdiff(observed$Auto_malignancy, malignancy_levels)
    if (length(unknown)) {
      stop(sample_name, " has unexpected epithelial malignancy classes: ", paste(unknown, collapse = ", "))
    }
    completed <- observed[
      data.table(Auto_malignancy = malignancy_levels), on = "Auto_malignancy"
    ]
    completed[is.na(N), N := 0L]
    completed[, sample := sample_name]
    completed
  })
  count_table <- rbindlist(count_rows, use.names = TRUE)
  count_table[, `:=`(
    sample = factor(sample, levels = samples),
    malignancy = factor(
      unname(malignancy_labels[Auto_malignancy]),
      levels = unname(malignancy_labels[malignancy_levels])
    )
  )]
  count_table[, total_epithelial_bins := sum(N), by = sample]
  count_table[, percent_epithelial_bins := 100 * N / total_epithelial_bins]
  count_table[, class_order := match(Auto_malignancy, malignancy_levels)]
  setorder(count_table, sample, class_order)
  count_table[, `:=`(
    segment_end = cumsum(N),
    segment_start = shift(cumsum(N), fill = 0)
  ), by = sample]
  count_table[, segment_mid := (segment_start + segment_end) / 2]
  count_table[, label_inside := percent_epithelial_bins >= 7]
  count_table[, label_y := fifelse(
    label_inside, 1,
    fifelse(class_order %% 2L == 0L, 0.48, 1.52)
  )]

  plot_data <- count_table[N > 0]
  count_plot <- ggplot(plot_data, aes(fill = malignancy)) +
    geom_rect(
      aes(xmin = segment_start, xmax = segment_end, ymin = 0.71, ymax = 1.29),
      linewidth = 0.3, colour = "white"
    ) +
    geom_segment(
      data = plot_data[label_inside == FALSE],
      aes(x = segment_mid, xend = segment_mid, y = 0.72, yend = label_y),
      inherit.aes = FALSE, linewidth = 0.45, colour = "grey35"
    ) +
    geom_text(
      aes(x = segment_mid, y = label_y, label = scales::comma(N)),
      inherit.aes = FALSE, size = 4.2, fontface = "bold", colour = "black"
    ) +
    facet_wrap(vars(sample), ncol = 1, scales = "free_x") +
    scale_fill_manual(values = malignancy_colours, drop = FALSE) +
    scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.06))) +
    scale_y_continuous(limits = c(0.25, 1.75), breaks = NULL) +
    labs(x = "Post-filter epithelial bins", y = NULL, fill = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_classic(base_size = 15) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.7, "cm"),
      legend.text = element_text(size = 12),
      panel.spacing.y = grid::unit(0.45, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(size = 16, face = "bold", hjust = 0)
    )

  count_table[, sample := as.character(sample)]
  count_table[, malignancy := as.character(malignancy)]
  setcolorder(
    count_table,
    c(
      "sample", "Auto_malignancy", "malignancy", "N", "percent_epithelial_bins",
      "total_epithelial_bins", "class_order", "segment_start", "segment_end",
      "segment_mid", "label_inside", "label_y"
    )
  )
  fwrite(
    count_table,
    file.path(output_dir, "tables", "Auto_visium_hd_binned_malignancy_class_counts.csv")
  )
  fwrite(
    count_table[, .(sample, Auto_malignancy, malignancy, N, percent_epithelial_bins, total_epithelial_bins)],
    file.path(summary_dir, "visium_hd_binned_malignancy_class_counts.csv")
  )
  ggsave(
    file.path(output_dir, "figures", "Auto_visium_hd_binned_malignancy_class_counts.pdf"),
    count_plot, width = 10, height = 9, device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, "figures", "Auto_visium_hd_binned_malignancy_class_counts.png"),
    count_plot, width = 10, height = 9, dpi = 300
  )
  invisible(count_table)
}
####################

summary_rows <- list()
selected_reference_rows <- list()
signature_gene_rows <- list()
failed_samples <- character()
run_log <- c(
  paste0("start=", format(Sys.time(), tz = "Europe/London")),
  paste0("samples=", paste(samples, collapse = ";")),
  paste0("normal_types=", paste(normal_types, collapse = ";")),
  paste0("min_reference_cells=", min_reference_cells),
  paste0("min_epithelial_cells=", min_epithelial_cells),
  paste0("cna_sd_k=", cna_sd_k),
  paste0("cancer_signature_threshold=", cancer_signature_threshold),
  paste0("cancer_signature_top_n=", cancer_signature_top_n),
  paste0("force_rebuild=", force_rebuild),
  paste0("replot_only=", replot_only)
)

for (sample_name in samples) {
  output_table_path <- file.path(output_dir, "tables", paste0("Auto_", sample_name, "_binned_malignancy.csv.gz"))
  if (replot_only) {
    if (!file.exists(output_table_path)) stop("Missing live malignancy table for replot: ", output_table_path)
    cell_table <- fread(output_table_path)
    save_diagnostic_plots(cell_table, sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=replotted"))
    next
  }

  annotation_path <- file.path(
    filter_dir, "tables", paste0("Auto_", sample_name, "_binned_filtered_annotations.csv.gz")
  )
  if (!file.exists(annotation_path)) stop("Missing filtered annotation: ", annotation_path)
  annotation <- fread(annotation_path)
  needed <- c(
    "barcode", "Auto_postfilter_celltype", "Auto_postfilter_keep",
    "pxl_col_in_fullres", "pxl_row_in_fullres"
  )
  missing <- setdiff(needed, names(annotation))
  if (length(missing)) stop("Filtered annotation missing columns: ", paste(missing, collapse = ", "))

  input_dir <- manifest[sample == sample_name, binned_input][1L]
  counts_path <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(counts_path)) stop("Missing count matrix: ", counts_path)
  counts <- read_10x_counts(counts_path)
  annotation <- annotation[barcode %in% colnames(counts)]
  annotation <- annotation[order(match(barcode, colnames(counts)))]

  normal_counts <- setNames(
    vapply(normal_types, function(type) {
      annotation[Auto_postfilter_keep %in% TRUE & Auto_postfilter_celltype == type, .N]
    }, integer(1)),
    normal_types
  )
  annotation[, Auto_reference_rctd_manual_concordant :=
    Auto_postfilter_keep %in% TRUE &
      !grepl("\\|", Auto_annotation_celltype) &
      Auto_annotation_celltype == Auto_postfilter_celltype &
      Auto_rctd_is_singlet %in% TRUE &
      Auto_rctd_first_type == Auto_postfilter_celltype &
      Auto_postfilter_celltype %in% normal_types]
  reference_counts <- setNames(
    vapply(normal_types, function(type) {
      annotation[
        Auto_reference_rctd_manual_concordant %in% TRUE &
          Auto_postfilter_celltype == type, .N
      ]
    }, integer(1)),
    normal_types
  )
  eligible_types <- names(reference_counts)[reference_counts >= min_reference_cells]
  selected_types <- names(sort(reference_counts[eligible_types], decreasing = TRUE))[seq_len(min(2L, length(eligible_types)))]
  epithelial_barcodes <- annotation[
    Auto_postfilter_keep %in% TRUE & Auto_postfilter_celltype == "epithelial", barcode
  ]

  base_summary <- data.table(
    sample = sample_name,
    status = "pending",
    n_annotated = nrow(annotation),
    n_postfilter_keep = sum(annotation$Auto_postfilter_keep %in% TRUE),
    n_epithelial_targets = length(epithelial_barcodes),
    n_endothelial_available = normal_counts[["endothelial"]],
    n_macrophage_available = normal_counts[["macrophage"]],
    n_fibroblast_available = normal_counts[["fibroblast"]],
    n_endothelial_high_confidence = reference_counts[["endothelial"]],
    n_macrophage_high_confidence = reference_counts[["macrophage"]],
    n_fibroblast_high_confidence = reference_counts[["fibroblast"]],
    selected_reference_types = paste(selected_types, collapse = ";")
  )
  if (length(selected_types) < 2L) {
    base_summary[, status := "insufficient_reference_compartments"]
    summary_rows[[sample_name]] <- base_summary
    failed_samples <- c(failed_samples, sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=insufficient_reference_compartments"))
    rm(annotation, counts)
    gc()
    next
  }
  if (length(epithelial_barcodes) < min_epithelial_cells) {
    base_summary[, status := "insufficient_epithelial_targets"]
    summary_rows[[sample_name]] <- base_summary
    failed_samples <- c(failed_samples, sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=insufficient_epithelial_targets"))
    rm(annotation, counts)
    gc()
    next
  }

  # InferCNA constructs a separate mean profile for each reference group and
  # uses the range between those group means, so unequal group sizes do not
  # numerically weight the correction. Retain every independently concordant
  # reference bin to estimate each group mean as precisely as possible.
  ref_barcodes <- lapply(selected_types, function(type) {
    sort(annotation[
      Auto_reference_rctd_manual_concordant %in% TRUE &
        Auto_postfilter_celltype == type, barcode
    ])
  })
  names(ref_barcodes) <- selected_types
  reference_barcodes <- unlist(ref_barcodes, use.names = FALSE)
  use_barcodes <- unique(c(epithelial_barcodes, reference_barcodes))

  selected_reference_rows[[sample_name]] <- rbindlist(lapply(selected_types, function(type) {
    data.table(
      sample = sample_name,
      reference_type = type,
      reference_definition = "postfilter_exact_manual_label_and_RCTD_singlet_concordant",
      barcode = ref_barcodes[[type]]
    )
  }))

  keep_genes <- intersect(rownames(counts), genome_genes)
  if (length(keep_genes) < 5000L) {
    base_summary[, status := "too_few_genome_ordered_genes"]
    base_summary[, n_genome_ordered_genes := length(keep_genes)]
    summary_rows[[sample_name]] <- base_summary
    failed_samples <- c(failed_samples, sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=too_few_genome_ordered_genes"))
    rm(annotation, counts)
    gc()
    next
  }
  counts <- counts[keep_genes, use_barcodes, drop = FALSE]
  cpm <- to_cpm(counts)
  dimnames(cpm) <- dimnames(counts)
  rm(counts)
  gc()

  cache_path <- file.path(ephemeral_dir, paste0("Auto_", sample_name, "_binned_infercna_outs.rds"))
  cache_metadata_path <- file.path(ephemeral_dir, paste0("Auto_", sample_name, "_binned_infercna_cache_metadata.rds"))
  reused_cache <- FALSE
  if (!force_rebuild && file.exists(cache_path) && file.exists(cache_metadata_path)) {
    cache_metadata <- readRDS(cache_metadata_path)
    cache_matches <- identical(sort(cache_metadata$use_barcodes), sort(use_barcodes)) &&
      identical(cache_metadata$ref_barcodes, ref_barcodes) &&
      identical(cache_metadata$genes, rownames(cpm))
    if (cache_matches) {
      outs <- readRDS(cache_path)
      cache_matches <- !is.list(outs) && setequal(colnames(outs), colnames(cpm))
      if (cache_matches) {
        outs <- outs[, colnames(cpm), drop = FALSE]
        reused_cache <- TRUE
      } else {
        rm(outs)
      }
    }
  }

  if (!reused_cache) {
    message(
      "InferCNA: ", sample_name, " (epithelial targets=", length(epithelial_barcodes),
      ", references=", length(reference_barcodes), ")"
    )
    outs <- tryCatch(
      infercna::infercna(as.matrix(cpm), refCells = ref_barcodes, isLog = FALSE, verbose = TRUE),
      error = function(error) error
    )
    if (!inherits(outs, "error")) {
      saveRDS(outs, cache_path)
      saveRDS(
        list(
          sample = sample_name, use_barcodes = use_barcodes,
          ref_barcodes = ref_barcodes, genes = rownames(cpm),
          created = format(Sys.time(), tz = "Europe/London")
        ),
        cache_metadata_path
      )
    }
  }
  if (inherits(outs, "error")) {
    base_summary[, status := "infercna_error"]
    base_summary[, error_message := conditionMessage(outs)]
    summary_rows[[sample_name]] <- base_summary
    failed_samples <- c(failed_samples, sample_name)
    run_log <- c(run_log, paste0("sample=", sample_name, "; status=infercna_error; message=", conditionMessage(outs)))
    rm(annotation, cpm, outs)
    gc()
    next
  }

  scatter <- as.data.table(infercna::cnaScatterPlot(outs, refCells = ref_barcodes), keep.rownames = "barcode")
  scatter[, is_reference := barcode %in% reference_barcodes]
  scatter[, is_epithelial_target := barcode %in% epithelial_barcodes]
  ref_scatter <- scatter[is_reference %in% TRUE]
  threshold_signal <- mean(ref_scatter$cna.signal, na.rm = TRUE) + cna_sd_k * sd(ref_scatter$cna.signal, na.rm = TRUE)
  threshold_cor <- mean(ref_scatter$cna.cor, na.rm = TRUE) + cna_sd_k * sd(ref_scatter$cna.cor, na.rm = TRUE)
  if (!is.finite(threshold_signal)) threshold_signal <- max(ref_scatter$cna.signal, na.rm = TRUE)
  if (!is.finite(threshold_cor)) threshold_cor <- max(ref_scatter$cna.cor, na.rm = TRUE)

  scatter[, Auto_cna_class := "not_assessed"]
  scatter[is_reference %in% TRUE, Auto_cna_class := "reference"]
  scatter[
    is_epithelial_target %in% TRUE & cna.signal > threshold_signal & cna.cor > threshold_cor,
    Auto_cna_class := "cna_malignant"
  ]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "not_assessed" &
      (cna.signal > threshold_signal | cna.cor > threshold_cor),
    Auto_cna_class := "cna_unresolved"
  ]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "not_assessed",
    Auto_cna_class := "cna_non_malignant"
  ]

  signature_genes_present <- intersect(cancer_signature_genes, rownames(cpm))
  if (!length(signature_genes_present)) stop("No cancer signature genes present for ", sample_name)
  signature_gene_means <- Matrix::rowMeans(
    log1p(cpm[signature_genes_present, epithelial_barcodes, drop = FALSE] / 100),
    na.rm = TRUE
  )
  signature_genes_selected <- names(sort(signature_gene_means, decreasing = TRUE))[
    seq_len(min(cancer_signature_top_n, length(signature_gene_means)))
  ]
  signature_scores <- Matrix::colMeans(
    log1p(cpm[signature_genes_selected, , drop = FALSE] / 100),
    na.rm = TRUE
  )
  signature_gene_rows[[sample_name]] <- data.table(
    sample = sample_name,
    rank = seq_along(signature_genes_selected),
    gene = signature_genes_selected,
    mean_log1p_cp10k_in_epithelial = unname(signature_gene_means[signature_genes_selected])
  )
  scatter[, Auto_cancer_signature_score := unname(signature_scores[match(barcode, names(signature_scores))])]
  scatter[, Auto_cancer_signature_status := "not_assessed"]
  scatter[
    is_epithelial_target %in% TRUE,
    Auto_cancer_signature_status := fifelse(
      Auto_cancer_signature_score >= cancer_signature_threshold,
      "cs_malignant", "cs_unresolved"
    )
  ]

  scatter[, `:=`(Auto_malignancy = "not_assessed", Auto_malignancy_evidence = "not_assessed")]
  scatter[is_reference %in% TRUE, `:=`(Auto_malignancy = "reference", Auto_malignancy_evidence = "normal_reference")]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "cna_non_malignant",
    `:=`(Auto_malignancy = "non_malignant", Auto_malignancy_evidence = "both_cna_metrics_below_threshold")
  ]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "cna_unresolved",
    `:=`(Auto_malignancy = "unresolved", Auto_malignancy_evidence = "one_cna_metric_above_threshold")
  ]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "cna_malignant",
    `:=`(Auto_malignancy = "malignant_level_1", Auto_malignancy_evidence = "both_cna_metrics_above_threshold")
  ]
  scatter[
    is_epithelial_target %in% TRUE & Auto_cna_class == "cna_unresolved" &
      Auto_cancer_signature_status == "cs_malignant",
    `:=`(Auto_malignancy = "malignant_level_2", Auto_malignancy_evidence = "one_cna_metric_plus_cancer_signature")
  ]
  scatter[, Auto_malignant := is_epithelial_target %in% TRUE & Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")]
  scatter[, `:=`(
    Auto_cna_signal_threshold = threshold_signal,
    Auto_cna_cor_threshold = threshold_cor,
    Auto_cna_sd_k = cna_sd_k,
    Auto_cancer_signature_threshold = cancer_signature_threshold,
    Auto_cancer_signature_n_genes = length(signature_genes_selected)
  )]

  annotation[, `:=`(
    cna.signal = NA_real_, cna.cor = NA_real_, is_reference = FALSE,
    is_epithelial_target = FALSE, Auto_cna_class = NA_character_,
    Auto_cancer_signature_score = NA_real_, Auto_cancer_signature_status = NA_character_,
    Auto_malignancy = fifelse(Auto_postfilter_keep %in% TRUE, "not_assessed", "filtered_before_malignancy"),
    Auto_malignancy_evidence = fifelse(Auto_postfilter_keep %in% TRUE, "not_assessed", Auto_postfilter_reason),
    Auto_malignant = NA,
    Auto_cna_signal_threshold = threshold_signal,
    Auto_cna_cor_threshold = threshold_cor,
    Auto_cna_sd_k = cna_sd_k,
    Auto_cancer_signature_threshold = cancer_signature_threshold,
    Auto_cancer_signature_n_genes = length(signature_genes_selected)
  )]
  scatter_match <- match(annotation$barcode, scatter$barcode)
  scatter_columns <- c(
    "cna.signal", "cna.cor", "is_reference", "is_epithelial_target", "Auto_cna_class",
    "Auto_cancer_signature_score", "Auto_cancer_signature_status", "Auto_malignancy",
    "Auto_malignancy_evidence", "Auto_malignant"
  )
  for (column_name in scatter_columns) {
    matched <- !is.na(scatter_match)
    set(annotation, i = which(matched), j = column_name, value = scatter[[column_name]][scatter_match[matched]])
  }
  fwrite(annotation, output_table_path)

  target_scatter <- scatter[is_epithelial_target %in% TRUE]
  base_summary[, `:=`(
    status = "complete",
    selected_reference_cell_counts = paste(
      paste0(selected_types, "=", lengths(ref_barcodes)), collapse = ";"
    ),
    n_reference_cells_total = length(reference_barcodes),
    n_genome_ordered_genes = nrow(outs),
    reused_infercna_cache = reused_cache,
    cna_signal_threshold = threshold_signal,
    cna_cor_threshold = threshold_cor,
    cna_sd_k = cna_sd_k,
    cancer_signature_threshold = cancer_signature_threshold,
    n_cancer_signature_genes = length(signature_genes_selected),
    n_cna_malignant = sum(target_scatter$Auto_cna_class == "cna_malignant"),
    n_cna_unresolved = sum(target_scatter$Auto_cna_class == "cna_unresolved"),
    n_cna_non_malignant = sum(target_scatter$Auto_cna_class == "cna_non_malignant"),
    n_malignant_level_1 = sum(target_scatter$Auto_malignancy == "malignant_level_1"),
    n_malignant_level_2 = sum(target_scatter$Auto_malignancy == "malignant_level_2"),
    n_malignant_total = sum(target_scatter$Auto_malignant),
    pct_malignant_epithelial = 100 * mean(target_scatter$Auto_malignant),
    n_signature_positive = sum(target_scatter$Auto_cancer_signature_status == "cs_malignant")
  )]
  summary_rows[[sample_name]] <- base_summary
  save_diagnostic_plots(annotation, sample_name)
  run_log <- c(run_log, paste0(
    "sample=", sample_name,
    "; status=complete",
    "; targets=", nrow(target_scatter),
    "; references=", length(reference_barcodes),
    "; malignant=", sum(target_scatter$Auto_malignant),
    "; pct_malignant=", round(100 * mean(target_scatter$Auto_malignant), 2),
    "; reused_cache=", reused_cache
  ))

  rm(annotation, cpm, outs, scatter, target_scatter, ref_scatter, signature_scores)
  gc()
}

####################
if (!length(failed_samples)) {
  save_malignancy_count_plot()
}
####################

if (!replot_only) {
  summary_table <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  fwrite(summary_table, file.path(output_dir, "tables", "Auto_visium_hd_binned_malignancy_summary.csv"))
  fwrite(summary_table, file.path(summary_dir, "visium_hd_binned_malignancy_summary.csv"))
  if (length(selected_reference_rows)) {
    fwrite(
      rbindlist(selected_reference_rows, use.names = TRUE, fill = TRUE),
      file.path(output_dir, "tables", "Auto_visium_hd_binned_selected_reference_bins.csv.gz")
    )
  }
  if (length(signature_gene_rows)) {
    fwrite(
      rbindlist(signature_gene_rows, use.names = TRUE, fill = TRUE),
      file.path(output_dir, "tables", "Auto_visium_hd_binned_selected_cancer_signature_genes.csv")
    )
  }
  parameters <- data.table(
    parameter = c(
      "samples", "normal_reference_candidates", "normal_reference_types_required",
      "min_reference_cells_available", "reference_confidence_definition",
      "reference_group_size_handling", "min_epithelial_cells", "cna_threshold",
      "cancer_signature_source", "cancer_signature_top_n", "cancer_signature_threshold",
      "signature_rescue_scope", "infercna_intermediate_storage"
    ),
    value = c(
      paste(samples, collapse = ";"), paste(normal_types, collapse = ";"), "2",
      as.character(min_reference_cells),
      "postfilter exact manual label concordant with an RCTD singlet first type",
      "all qualifying bins retained; InferCNA corrects against separate group means",
      as.character(min_epithelial_cells),
      paste0("both metrics > pooled reference mean + ", cna_sd_k, " SD"),
      normalizePath(cancer_signature_path), as.character(cancer_signature_top_n),
      as.character(cancer_signature_threshold), "CNA-unresolved epithelial bins only",
      ephemeral_dir
    )
  )
  fwrite(parameters, file.path(output_dir, "tables", "Auto_visium_hd_binned_malignancy_parameters.csv"))
}

run_log <- c(
  run_log,
  paste0("end=", format(Sys.time(), tz = "Europe/London")),
  "session_info:", capture.output(sessionInfo())
)
writeLines(run_log, file.path(output_dir, "logs", "Auto_visium_hd_binned_malignancy_run_summary.txt"))

if (length(failed_samples)) {
  stop("Malignancy classification failed for: ", paste(unique(failed_samples), collapse = ", "),
       ". See the sample summary for reference and epithelial counts.")
}
####################
