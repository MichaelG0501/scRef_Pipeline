####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/bulk_tcga_geo_qc.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_bulk_tcga_geo_qc.R
#
# Cross-platform bulk QC for TCGA-ESCA whole-bulk RNA-seq and GEO GSE19417
# microarray. Harmonizes genes, applies dataset-wise transforms, visualizes PCA
# structure, adds sample-level expression-strength metrics, and flags samples for
# review/removal using a deliberately conservative rule set.
#
# Inputs:
#   ref_outs/tcga_esca_meta.rds
#   ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#
# Outputs:
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_review.pdf
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_pca_before.pdf
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_pca_after.pdf
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_sample_table.csv
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_removed_samples.csv
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_retained_samples.csv
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_preprocessing_summary.csv
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_objects.rds
#   updates/new_updates/summaries/Auto_bulk_tcga_geo_qc_summary.csv
####################

library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- "bulk_crossplatform"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("../updates/new_updates/summaries", recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

write_grob_pdf <- function(path, grob_list, width, height) {
  grob_list <- grob_list[!vapply(grob_list, is.null, logical(1))]
  if (length(grob_list) == 0) stop("No grobs available to write: ", path)
  open_pdf_device(path, width = width, height = height)
  for (g in grob_list) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  dev.off()
}

infer_tcga_histology <- function(type_vec) {
  type_vec <- tolower(as.character(type_vec))
  out <- rep("Other", length(type_vec))
  out[grepl("adeno", type_vec)] <- "EAC"
  out[grepl("squamous", type_vec)] <- "ESCC"
  out
}

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1
  scaled <- sweep(mat, 1, row_means, "-")
  scaled <- sweep(scaled, 1, row_sds, "/")
  scaled[is.na(scaled)] <- 0
  scaled
}

select_variable_genes <- function(mat, top_n = 5000) {
  gene_mad <- apply(mat, 1, mad, na.rm = TRUE)
  gene_mad[is.na(gene_mad)] <- 0
  ord <- order(gene_mad, decreasing = TRUE)
  keep_n <- min(top_n, sum(gene_mad > 0))
  if (keep_n < 500) keep_n <- min(nrow(mat), top_n)
  rownames(mat)[ord[seq_len(keep_n)]]
}

collapse_reasons <- function(...) {
  vals <- c(...)
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) return("")
  paste(vals, collapse = "; ")
}

calc_expression_metrics <- function(expr_mat, raw_mat = NULL) {
  expr_threshold <- as.numeric(quantile(expr_mat, probs = 0.35, na.rm = TRUE))
  out <- data.frame(
    sample_id = colnames(expr_mat),
    sample_median_expression = apply(expr_mat, 2, median, na.rm = TRUE),
    sample_iqr_expression = apply(expr_mat, 2, IQR, na.rm = TRUE),
    sample_sd_expression = apply(expr_mat, 2, sd, na.rm = TRUE),
    sample_expression_breadth = colMeans(expr_mat > expr_threshold, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if (!is.null(raw_mat)) {
    raw_mat <- raw_mat[rownames(expr_mat), colnames(expr_mat), drop = FALSE]
    out$sample_detected_genes_tpm1 <- colSums(raw_mat > 1, na.rm = TRUE)
  } else {
    out$sample_detected_genes_tpm1 <- NA_real_
  }
  out
}

add_lower_tail_flags <- function(df, metric_name, n_mad = 4) {
  threshold_name <- paste0(metric_name, "_lower_threshold")
  flag_name <- paste0("flag_low_", metric_name)
  df %>%
    group_by(dataset) %>%
    mutate(
      !!threshold_name := {
        med_x <- median(.data[[metric_name]], na.rm = TRUE)
        mad_x <- mad(.data[[metric_name]], na.rm = TRUE)
        if (is.na(mad_x) || mad_x == 0) med_x - n_mad else med_x - n_mad * mad_x
      },
      !!flag_name := .data[[metric_name]] < .data[[threshold_name]]
    ) %>%
    ungroup()
}

histology_colors <- c(
  "EAC" = "#B2182B",
  "ESCC" = "#2166AC",
  "Gastric" = "#E69F00",
  "Other" = "grey70"
)

qc_status_colors <- c(
  "Keep" = "#333333",
  "Review" = "#E69F00",
  "Remove" = "#D73027"
)

tcga_meta <- readRDS("tcga_esca_meta.rds") %>%
  transmute(
    sample_id = sample_barcode,
    dataset = "TCGA",
    histology_group = infer_tcga_histology(type),
    source_histology = type,
    analysis_ready_for_survival = histology_group == "EAC" & !is.na(OS_time) & !is.na(OS_event),
    OS_time = OS_time,
    OS_event = OS_event
  )

tcga_df <- fread("cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tcga_expr_raw <- as.matrix(tcga_df[, -1, with = FALSE])
rownames(tcga_expr_raw) <- tcga_df[[1]]
storage.mode(tcga_expr_raw) <- "numeric"
tcga_expr_raw <- tcga_expr_raw[, tcga_meta$sample_id[tcga_meta$sample_id %in% colnames(tcga_expr_raw)], drop = FALSE]
tcga_meta <- tcga_meta %>% filter(sample_id %in% colnames(tcga_expr_raw))
tcga_expr <- log2(tcga_expr_raw + 1)

geo_meta <- readRDS(file.path("geo_survival", "Auto_GSE19417_meta.rds")) %>%
  transmute(
    sample_id = sample_geo_accession,
    dataset = "GEO_GSE19417",
    histology_group = HistologyGroup,
    source_histology = HistologyGroup,
    analysis_ready_for_survival = analysis_ready_for_survival,
    OS_time = OS_time,
    OS_event = OS_event
  )

geo_expr <- readRDS(file.path("geo_survival", "Auto_GSE19417_expr_gene.rds"))
geo_expr <- geo_expr[, geo_meta$sample_id[geo_meta$sample_id %in% colnames(geo_expr)], drop = FALSE]
geo_meta <- geo_meta %>% filter(sample_id %in% colnames(geo_expr))
geo_expr <- as.matrix(geo_expr)
storage.mode(geo_expr) <- "numeric"

shared_genes <- intersect(rownames(tcga_expr), rownames(geo_expr))
shared_genes <- sort(shared_genes)
tcga_shared <- tcga_expr[shared_genes, , drop = FALSE]
geo_shared <- geo_expr[shared_genes, , drop = FALSE]
tcga_raw_shared <- tcga_expr_raw[shared_genes, , drop = FALSE]

tcga_scaled <- row_zscore(tcga_shared)
geo_scaled <- row_zscore(geo_shared)
combined_scaled <- cbind(tcga_scaled, geo_scaled)

meta_all <- bind_rows(tcga_meta, geo_meta) %>%
  mutate(
    dataset = factor(dataset, levels = c("TCGA", "GEO_GSE19417")),
    histology_group = factor(histology_group, levels = c("EAC", "ESCC", "Gastric", "Other"))
  )

meta_all <- meta_all %>% filter(sample_id %in% colnames(combined_scaled))
combined_scaled <- combined_scaled[, meta_all$sample_id, drop = FALSE]

expression_metrics <- bind_rows(
  calc_expression_metrics(tcga_shared, tcga_raw_shared) %>% mutate(dataset = "TCGA"),
  calc_expression_metrics(geo_shared, NULL) %>% mutate(dataset = "GEO_GSE19417")
)

variable_genes <- select_variable_genes(combined_scaled, top_n = 5000)
pca <- prcomp(t(combined_scaled[variable_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
pc_keep <- min(10, ncol(pca$x))
pc_df <- as.data.frame(pca$x[, seq_len(pc_keep), drop = FALSE])
pc_df$sample_id <- rownames(pc_df)
pc_df <- pc_df %>%
  left_join(meta_all, by = "sample_id") %>%
  left_join(expression_metrics, by = c("sample_id", "dataset"))

qc_df <- pc_df %>% filter(histology_group %in% c("EAC", "ESCC"))
pc_cols <- paste0("PC", seq_len(pc_keep))
pc_mat <- as.matrix(qc_df[, pc_cols, drop = FALSE])

centroids <- lapply(split(seq_len(nrow(qc_df)), qc_df$histology_group), function(idx) {
  colMeans(pc_mat[idx, , drop = FALSE])
})

predicted_histology <- vapply(seq_len(nrow(qc_df)), function(i) {
  dists <- vapply(names(centroids), function(cl) {
    sqrt(sum((pc_mat[i, ] - centroids[[cl]])^2))
  }, numeric(1))
  names(which.min(dists))
}, character(1))

own_centroid_distance <- vapply(seq_len(nrow(qc_df)), function(i) {
  cls <- as.character(qc_df$histology_group[i])
  sqrt(sum((pc_mat[i, ] - centroids[[cls]])^2))
}, numeric(1))

dist_mat <- as.matrix(dist(pc_mat))
k_neighbors <- min(15, nrow(qc_df) - 1)
nn_same_histology_fraction <- rep(NA_real_, nrow(qc_df))
if (k_neighbors >= 1) {
  for (i in seq_len(nrow(qc_df))) {
    ord <- order(dist_mat[i, ])
    ord <- ord[ord != i]
    near_idx <- ord[seq_len(min(k_neighbors, length(ord)))]
    nn_same_histology_fraction[i] <- mean(
      as.character(qc_df$histology_group[near_idx]) == as.character(qc_df$histology_group[i])
    )
  }
}

distance_thresholds <- tapply(own_centroid_distance, qc_df$histology_group, function(x) {
  med_x <- median(x, na.rm = TRUE)
  mad_x <- mad(x, na.rm = TRUE)
  if (is.na(mad_x) || mad_x == 0) med_x + 2.5 else med_x + 2.5 * mad_x
})

qc_df <- qc_df %>%
  mutate(
    predicted_histology = predicted_histology,
    own_centroid_distance = own_centroid_distance,
    nn_same_histology_fraction = nn_same_histology_fraction,
    flag_histology_mismatch = predicted_histology != as.character(histology_group),
    flag_knn_mismatch = !is.na(nn_same_histology_fraction) & nn_same_histology_fraction < 0.6,
    flag_strong_knn_mismatch = !is.na(nn_same_histology_fraction) & nn_same_histology_fraction < 0.5,
    flag_distance_outlier = own_centroid_distance > unname(distance_thresholds[as.character(histology_group)])
  )

qc_df <- qc_df %>%
  mutate(
    qc_remove = flag_knn_mismatch & flag_distance_outlier
  )

pc_df <- pc_df %>%
  left_join(
    qc_df %>%
      select(
        sample_id,
        predicted_histology,
        own_centroid_distance,
        nn_same_histology_fraction,
        starts_with("flag_"),
        qc_remove
      ),
    by = "sample_id"
  ) %>%
  mutate(
    predicted_histology = ifelse(is.na(predicted_histology), as.character(histology_group), predicted_histology),
    own_centroid_distance = ifelse(is.na(own_centroid_distance), NA_real_, own_centroid_distance),
    nn_same_histology_fraction = ifelse(is.na(nn_same_histology_fraction), NA_real_, nn_same_histology_fraction),
    distance_threshold = unname(distance_thresholds[as.character(histology_group)]),
    flag_histology_mismatch = ifelse(is.na(flag_histology_mismatch), FALSE, flag_histology_mismatch),
    flag_knn_mismatch = ifelse(is.na(flag_knn_mismatch), FALSE, flag_knn_mismatch),
    flag_strong_knn_mismatch = ifelse(is.na(flag_strong_knn_mismatch), FALSE, flag_strong_knn_mismatch),
    flag_distance_outlier = ifelse(is.na(flag_distance_outlier), FALSE, flag_distance_outlier),
    qc_remove = ifelse(is.na(qc_remove), FALSE, qc_remove)
  )

pc_df$qc_status <- ifelse(pc_df$qc_remove, "Remove", "Keep")
pc_df$qc_status <- factor(pc_df$qc_status, levels = c("Keep", "Remove"))

pc_df$qc_reason <- vapply(seq_len(nrow(pc_df)), function(i) {
  collapse_reasons(
    if (pc_df$flag_knn_mismatch[i]) "knn_mismatch" else "",
    if (pc_df$flag_distance_outlier[i]) "distance_outlier" else ""
  )
}, character(1))

pc_df$integration_keep <- !pc_df$qc_remove

plot_pca_combined <- function(df, title_text, keep_only = FALSE) {
  plot_df <- if (keep_only) df %>% filter(!qc_remove) else df
  p <- ggplot(plot_df, aes(PC1, PC2, color = histology_group, shape = dataset)) +
    geom_point(aes(alpha = qc_status), size = 3) +
    scale_color_manual(values = histology_colors, drop = FALSE) +
    scale_alpha_manual(values = c("Keep" = 1, "Remove" = 1), guide = "none") +
    theme_classic(base_size = 12) +
    labs(title = title_text, color = "Histology", shape = "Dataset")
  p
}

plot_dataset_pca <- function(df, dataset_name) {
  plot_df <- df %>% filter(dataset == dataset_name)
  ggplot(plot_df, aes(PC1, PC2, color = histology_group, shape = qc_status)) +
    geom_point(size = 3) +
    scale_color_manual(values = histology_colors, drop = FALSE) +
    scale_shape_manual(values = c("Keep" = 16, "Remove" = 8), drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      title = paste0(dataset_name, ": PCA review"),
      subtitle = "Joint PCA coordinates, faceted by dataset for manual QC review",
      color = "Histology",
      shape = "QC status"
    )
}

plot_dataset_histology_metrics <- function(df, dataset_name) {
  plot_df <- df %>% filter(dataset == dataset_name, histology_group %in% c("EAC", "ESCC"))
  ggplot(plot_df, aes(own_centroid_distance, nn_same_histology_fraction, color = histology_group, shape = qc_status)) +
    geom_vline(
      data = plot_df %>% distinct(histology_group, distance_threshold),
      aes(xintercept = distance_threshold, color = histology_group),
      linetype = "dotted",
      linewidth = 0.5,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = 0.6,
      linetype = "dotted",
      linewidth = 0.5,
      color = "grey40"
    ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = histology_colors, drop = FALSE) +
    scale_shape_manual(values = c("Keep" = 16, "Remove" = 8), drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      title = paste0(dataset_name, ": histology-consistency metrics"),
      subtitle = "Dotted lines show the QC thresholds used for removal",
      x = "Distance to own histology centroid",
      y = "Fraction of 15 nearest neighbors with same histology",
      color = "Histology",
      shape = "QC status"
    )
}

plot_dataset_expression_metrics <- function(df, dataset_name) {
  plot_df <- df %>% filter(dataset == dataset_name)
  ggplot(plot_df, aes(sample_median_expression, sample_expression_breadth, color = histology_group, shape = qc_status)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = histology_colors, drop = FALSE) +
    scale_shape_manual(values = c("Keep" = 16, "Remove" = 8), drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      title = paste0(dataset_name, ": overall expression strength"),
      x = "Median expression",
      y = "Expression breadth",
      color = "Histology",
      shape = "QC status"
    )
}

plot_dataset_spread_metrics <- function(df, dataset_name) {
  plot_df <- df %>% filter(dataset == dataset_name)
  ggplot(plot_df, aes(sample_median_expression, sample_iqr_expression, color = histology_group, shape = qc_status)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = histology_colors, drop = FALSE) +
    scale_shape_manual(values = c("Keep" = 16, "Remove" = 8), drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(
      title = paste0(dataset_name, ": expression spread"),
      x = "Median expression",
      y = "IQR expression",
      color = "Histology",
      shape = "QC status"
    )
}

review_pages <- lapply(c("TCGA", "GEO_GSE19417"), function(dataset_name) {
  plot_df <- pc_df %>% filter(dataset == dataset_name)
  gridExtra::arrangeGrob(
    plot_dataset_pca(plot_df, dataset_name),
    plot_dataset_histology_metrics(plot_df, dataset_name),
    plot_dataset_expression_metrics(plot_df, dataset_name),
    plot_dataset_spread_metrics(plot_df, dataset_name),
    ncol = 2,
    top = grid::textGrob(
      paste0(dataset_name, " QC review page"),
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
})

unlink(file.path(out_dir, "Auto_bulk_crossplatform_qc_metrics.pdf"), force = TRUE)
write_grob_pdf(
  file.path(out_dir, "Auto_bulk_crossplatform_qc_review.pdf"),
  grob_list = review_pages,
  width = 13,
  height = 10
)

open_pdf_device(file.path(out_dir, "Auto_bulk_crossplatform_qc_pca_before.pdf"), width = 9.5, height = 7.5)
print(plot_pca_combined(pc_df, "Cross-platform PCA Before QC", keep_only = FALSE))
dev.off()

open_pdf_device(file.path(out_dir, "Auto_bulk_crossplatform_qc_pca_after.pdf"), width = 9.5, height = 7.5)
print(plot_pca_combined(pc_df, "Cross-platform PCA After QC", keep_only = TRUE))
dev.off()

preprocessing_summary <- bind_rows(
  data.frame(
    dataset = "TCGA",
    input_scale = "TPM",
    transform_applied = "log2(TPM + 1)",
    standardization = "row-wise z-score within dataset",
    qc_expression_metrics = "median expression, IQR expression, expression breadth, detected genes > TPM1",
    sample_count = ncol(tcga_shared),
    eac_count = sum(tcga_meta$histology_group == "EAC"),
    escc_count = sum(tcga_meta$histology_group == "ESCC"),
    gastric_count = 0,
    survival_ready_count = sum(tcga_meta$analysis_ready_for_survival),
    shared_gene_count = length(shared_genes),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "GEO_GSE19417",
    input_scale = "processed microarray log-scale",
    transform_applied = "kept as supplied by GEO",
    standardization = "row-wise z-score within dataset",
    qc_expression_metrics = "median expression, IQR expression, expression breadth",
    sample_count = ncol(geo_shared),
    eac_count = sum(geo_meta$histology_group == "EAC"),
    escc_count = sum(geo_meta$histology_group == "ESCC"),
    gastric_count = sum(geo_meta$histology_group == "Gastric"),
    survival_ready_count = sum(geo_meta$analysis_ready_for_survival),
    shared_gene_count = length(shared_genes),
    stringsAsFactors = FALSE
  )
)

fwrite(pc_df, file.path(out_dir, "Auto_bulk_crossplatform_qc_sample_table.csv"))
fwrite(pc_df %>% filter(qc_remove), file.path(out_dir, "Auto_bulk_crossplatform_qc_removed_samples.csv"))
fwrite(pc_df %>% filter(!qc_remove), file.path(out_dir, "Auto_bulk_crossplatform_qc_retained_samples.csv"))
fwrite(preprocessing_summary, file.path(out_dir, "Auto_bulk_crossplatform_preprocessing_summary.csv"))

saveRDS(
  list(
    sample_table = pc_df,
    preprocessing_summary = preprocessing_summary,
    shared_genes = shared_genes,
    variable_genes = variable_genes,
    pca_rotation = pca$rotation,
    pca_variance = (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  ),
  file.path(out_dir, "Auto_bulk_crossplatform_qc_objects.rds")
)

summary_df <- data.frame(
  metric = c(
    "shared_genes",
    "tcga_samples",
    "geo_samples",
    "samples_removed_total",
    "samples_removed_tcga",
    "samples_removed_geo",
    "eac_survival_retained_tcga",
    "eac_survival_retained_geo"
  ),
  value = c(
    length(shared_genes),
    ncol(tcga_shared),
    ncol(geo_shared),
    sum(pc_df$qc_remove),
    sum(pc_df$qc_remove & pc_df$dataset == "TCGA"),
    sum(pc_df$qc_remove & pc_df$dataset == "GEO_GSE19417"),
    sum(!pc_df$qc_remove & pc_df$dataset == "TCGA" & pc_df$histology_group == "EAC" & pc_df$analysis_ready_for_survival),
    sum(!pc_df$qc_remove & pc_df$dataset == "GEO_GSE19417" & pc_df$histology_group == "EAC" & pc_df$analysis_ready_for_survival)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_df,
  "../updates/new_updates/summaries/Auto_bulk_tcga_geo_qc_summary.csv"
)
