####################
# Analysis registry:
#   Status: active terminal/validation workflow.
#   Script: analysis/cnv/scatlas_numbat_raw_expression_concordance.R
#   Methodology: analysis/methodology/cnv/scatlas_numbat_raw_expression_concordance_methodology.md
#   Inputs:
#     ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv
#     ref_outs/Auto_scatlas_numbat/by_samples/<sample>/numbat/gexp_roll_wide.tsv.gz
#     ref_outs/Auto_scatlas_numbat/by_samples/<sample>/numbat/Auto_<sample>_numbat_cell_map.csv
#     ref_outs/by_samples/<sample>/<sample>_outs.rds
#     optional ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv
#     /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt
#   Outputs:
#     ref_outs/Auto_scatlas_numbat/raw_expression_concordance/figures/Auto_scatlas_numbat_raw_expression_infercna_matched_heatmaps.pdf
#     ref_outs/Auto_scatlas_numbat/raw_expression_concordance/figures/per_sample/Auto_<sample>_raw_expression_infercna_concordance.{pdf,png}
#     ref_outs/Auto_scatlas_numbat/raw_expression_concordance/tables/Auto_scatlas_numbat_raw_expression_infercna_summary.csv
#     ref_outs/Auto_scatlas_numbat/raw_expression_concordance/tables/Auto_scatlas_numbat_raw_expression_cell_clusters.csv
#     ref_outs/Auto_scatlas_numbat/raw_expression_concordance/logs/Auto_scatlas_numbat_raw_expression_concordance_run_summary.txt
#   Cache/replot behavior:
#     Rebuilds plots/tables from the raw Numbat and InferCNA matrices each run.
#     Use first CLI argument to restrict to comma-separated sample IDs.
#   Run command:
#     Rscript analysis/cnv/scatlas_numbat_raw_expression_concordance.R all 1200
#   Conda environment:
#     dmtcp
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
  library(cluster)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
out_root <- file.path(root_dir, "ref_outs")
setwd(out_root)

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "all"
max_plot_cells <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 1200L
gene_bin_size <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 100L

manifest_path <- "Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
infer_cell_path <- "Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv"

if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)

out_dir <- "Auto_scatlas_numbat/raw_expression_concordance"
fig_dir <- file.path(out_dir, "figures")
per_sample_dir <- file.path(fig_dir, "per_sample")
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

combined_pdf <- file.path(fig_dir, "Auto_scatlas_numbat_raw_expression_infercna_matched_heatmaps.pdf")
summary_csv <- file.path(table_dir, "Auto_scatlas_numbat_raw_expression_infercna_summary.csv")
cell_cluster_csv <- file.path(table_dir, "Auto_scatlas_numbat_raw_expression_cell_clusters.csv")
run_summary_txt <- file.path(log_dir, "Auto_scatlas_numbat_raw_expression_concordance_run_summary.txt")

run_start <- Sys.time()

manifest <- fread(manifest_path)
if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  manifest <- manifest[sample %in% requested]
}
if (nrow(manifest) == 0) stop("No Numbat samples selected by argument: ", sample_arg)

chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- fread(gene_order_path, header = FALSE, fill = TRUE)
gene_order <- gene_order[, seq_len(4), with = FALSE]
setnames(gene_order, c("gene", "chromosome", "start", "end"))
gene_order <- gene_order[!is.na(gene) & nzchar(gene)]
gene_order <- gene_order %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(
    start = as.numeric(.data$start),
    end = as.numeric(.data$end),
    chromosome = factor(.data$chromosome, levels = chrom_levels)
  ) %>%
  filter(is.finite(.data$start), is.finite(.data$end)) %>%
  arrange(.data$chromosome, .data$start) %>%
  distinct(.data$gene, .keep_all = TRUE)

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

infer_cells <- NULL
if (file.exists(infer_cell_path)) {
  infer_cells <- fread(infer_cell_path) %>%
    mutate(
      cell_id = as.character(.data$cell),
      sample = as.character(.data$sample),
      infercna_subclone = as.character(.data$subclone),
      state_label = as.character(.data$state_label),
      top_mp_label = as.character(.data$top_mp_label)
    ) %>%
    select(.data$cell_id, .data$sample, .data$infercna_subclone, .data$state_label, .data$top_mp_label)
}

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  setNames(colorRampPalette(base)(length(values)), values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values[is.na(values) | !nzchar(values)] <- "NA"
  values <- values[!is.na(values) & nzchar(values)]
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) cols <- c(cols, make_palette(missing_values, palette))
  cols[values]
}

comb2 <- function(x) x * (x - 1) / 2

adjusted_rand <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.character(x[ok])
  y <- as.character(y[ok])
  n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y)
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_i * sum_j / total
  max_index <- (sum_i + sum_j) / 2
  denom <- max_index - expected
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (sum_ij - expected) / denom
}

parse_done_status <- function(sample_id, numbat_dir) {
  done_file <- file.path(numbat_dir, paste0("Auto_", sample_id, "_numbat_done.txt"))
  if (!file.exists(done_file)) {
    return(list(status = "missing_done_file", terminal_no_subclone = NA, iteration = NA_character_))
  }
  done_lines <- readLines(done_file, warn = FALSE)
  get_field <- function(key) {
    hit <- done_lines[grepl(paste0("^", key, "\\t"), done_lines)]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", key, "\\t"), "", hit[1])
  }
  list(
    status = get_field("status"),
    terminal_no_subclone = identical(get_field("terminal_no_subclone"), "TRUE"),
    iteration = get_field("numbat_iteration")
  )
}

resolve_infer_columns <- function(sample_id, cell_ids, raw_barcodes, infer_cols) {
  resolved <- vapply(seq_along(cell_ids), function(i) {
    cell_id <- cell_ids[i]
    raw_barcode <- raw_barcodes[i]
    candidates <- unique(c(
      cell_id,
      raw_barcode,
      paste(sample_id, raw_barcode, sep = "__"),
      paste(sample_id, raw_barcode, sep = "_")
    ))
    hit <- candidates[candidates %in% infer_cols]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }, character(1))
  resolved
}

read_clone_labels <- function(sample_id, numbat_dir, cell_map) {
  clone_file <- file.path(numbat_dir, paste0("Auto_", sample_id, "_numbat_clone_post.csv"))
  if (!file.exists(clone_file)) {
    return(setNames(rep("No accepted clone", nrow(cell_map)), cell_map$cell_id))
  }
  clone_post <- fread(clone_file)
  if (nrow(clone_post) == 0 || !all(c("cell", "clone_opt") %in% colnames(clone_post))) {
    return(setNames(rep("No accepted clone", nrow(cell_map)), cell_map$cell_id))
  }
  map_use <- cell_map %>%
    select(.data$cell_id, .data$numbat_cell, .data$raw_barcode, .data$numbat_cell_prefixed)
  clone_post <- clone_post %>%
    left_join(map_use, by = c("cell" = "numbat_cell"))
  missing_cell <- is.na(clone_post$cell_id)
  if (any(missing_cell)) {
    fallback <- map_use$cell_id[match(clone_post$cell[missing_cell], map_use$raw_barcode)]
    clone_post$cell_id[missing_cell] <- fallback
  }
  missing_cell <- is.na(clone_post$cell_id)
  if (any(missing_cell)) {
    fallback <- map_use$cell_id[match(clone_post$cell[missing_cell], map_use$numbat_cell_prefixed)]
    clone_post$cell_id[missing_cell] <- fallback
  }
  out <- setNames(rep("No accepted clone", nrow(cell_map)), cell_map$cell_id)
  ok <- !is.na(clone_post$cell_id)
  out[clone_post$cell_id[ok]] <- paste0("Numbat clone ", clone_post$clone_opt[ok])
  out
}

make_gene_bins <- function(go, bin_size = gene_bin_size) {
  go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(
      bin_index = ((row_number() - 1L) %/% bin_size) + 1L,
      bin = paste0(.data$chromosome, "_", .data$bin_index)
    ) %>%
    ungroup()
}

bin_gene_matrix <- function(mat, go, bin_size = gene_bin_size) {
  go2 <- make_gene_bins(go, bin_size)
  bins <- split(seq_len(nrow(go2)), factor(go2$bin, levels = unique(go2$bin)))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  bins_df <- go2 %>%
    group_by(.data$bin, .data$chromosome, .data$bin_index) %>%
    summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      midpoint = (.data$start + .data$end) / 2,
      arm = ifelse(.data$midpoint <= centromere_pos[as.character(.data$chromosome)], "p", "q"),
      arm_label = paste0(as.character(.data$chromosome), .data$arm)
    )
  bins_df <- bins_df[match(rownames(binned), bins_df$bin), , drop = FALSE]
  list(mat = binned, bins = bins_df)
}

arm_mean_matrix <- function(mat, bins_df) {
  arm_levels <- unique(bins_df$arm_label)
  out <- do.call(rbind, lapply(arm_levels, function(arm) {
    rows <- which(bins_df$arm_label == arm)
    colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(out) <- arm_levels
  out[!is.finite(out)] <- 0
  out
}

cluster_raw_expression <- function(numbat_binned) {
  n_cells <- ncol(numbat_binned)
  if (n_cells < 3 || nrow(numbat_binned) < 2) {
    out <- setNames(rep("RawExpr cluster 1", n_cells), colnames(numbat_binned))
    attr(out, "silhouette") <- NA_real_
    return(out)
  }
  mat <- numbat_binned
  mat[!is.finite(mat)] <- 0
  d <- dist(t(mat))
  hc <- hclust(d, method = "ward.D2")
  k <- min(2L, n_cells)
  raw <- cutree(hc, k = k)
  counts <- sort(table(raw), decreasing = TRUE)
  relabel <- setNames(paste0("RawExpr cluster ", seq_along(counts)), names(counts))
  out <- relabel[as.character(raw)]
  names(out) <- names(raw)
  if (length(unique(raw)) > 1) {
    sil <- tryCatch(mean(cluster::silhouette(raw, d)[, "sil_width"], na.rm = TRUE), error = function(e) NA_real_)
  } else {
    sil <- NA_real_
  }
  attr(out, "silhouette") <- sil
  out
}

sample_cells_for_plot <- function(cells, group, max_cells = max_plot_cells) {
  if (length(cells) <= max_cells) return(cells)
  group <- group[cells]
  group[is.na(group) | !nzchar(group)] <- "NA"
  split_cells <- split(cells, factor(group, levels = unique(group)))
  target <- pmax(10L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

order_cells_by_numbat <- function(numbat_binned, cells) {
  cells <- intersect(cells, colnames(numbat_binned))
  if (length(cells) <= 2) return(cells)
  mat <- numbat_binned[, cells, drop = FALSE]
  mat[!is.finite(mat)] <- 0
  d <- dist(t(mat))
  hc <- hclust(d, method = "ward.D2")
  cells[hc$order]
}

make_heatmap_grob <- function(mat, go, meta_plot, value_name, title, value_limit, show_annotation_legend = TRUE) {
  cells <- rownames(meta_plot)
  mat <- mat[, cells, drop = FALSE]
  row_chr <- factor(as.character(go$chromosome), levels = chrom_levels)
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(chrom_levels)), chrom_levels)

  raw_cols <- complete_palette(c(
    "RawExpr cluster 1" = "#2B8CBE",
    "RawExpr cluster 2" = "#D95F0E"
  ), meta_plot$raw_expr_cluster, "Set2")
  infer_cols <- complete_palette(c(
    "Subclone A" = "#D73027",
    "Subclone B" = "#4575B4",
    "Subclone C" = "#1A9850",
    "Subclone D" = "#984EA3",
    "Subclone E" = "#FF7F00",
    "Subclone F" = "#A65628"
  ), meta_plot$infercna_subclone, "Set2")
  clone_cols <- complete_palette(make_palette(meta_plot$numbat_clone, "Dark2"), meta_plot$numbat_clone, "Dark2")
  status_cols <- complete_palette(c(
    "Success" = "#1A9850",
    "No clones remain after filtering by size. Consider reducing min_cells." = "#B2182B",
    "No CNV remains after filtering by LLR in pseudobulks." = "#B2182B",
    "terminal_no_subclone" = "#B2182B"
  ), meta_plot$numbat_status, "Set3")

  top_ha <- HeatmapAnnotation(
    RawExpr = meta_plot$raw_expr_cluster,
    NumbatClone = meta_plot$numbat_clone,
    NumbatStatus = meta_plot$numbat_status_short,
    InferCNA = meta_plot$infercna_subclone,
    State = meta_plot$state_label,
    TopMP = meta_plot$top_mp_label,
    col = list(
      RawExpr = raw_cols,
      NumbatClone = clone_cols,
      NumbatStatus = complete_palette(c(Success = "#1A9850", NoFinalClone = "#B2182B"), meta_plot$numbat_status_short, "Set3"),
      InferCNA = infer_cols,
      State = complete_palette(make_palette(meta_plot$state_label, "Set3"), meta_plot$state_label, "Set3"),
      TopMP = complete_palette(make_palette(meta_plot$top_mp_label, "Paired"), meta_plot$top_mp_label, "Paired")
    ),
    simple_anno_size = unit(3.5, "mm"),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    show_legend = TRUE,
    na_col = "grey90"
  )

  left_ha <- rowAnnotation(
    Chr = row_chr,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )

  ht <- Heatmap(
    mat,
    name = value_name,
    col = colorRamp2(c(-value_limit, 0, value_limit), c("#2166AC", "white", "#B2182B")),
    left_annotation = left_ha,
    top_annotation = top_ha,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 6),
    use_raster = TRUE,
    raster_quality = 2,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title = value_name, legend_height = unit(25, "mm"))
  )

  grid.grabExpr(draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    show_annotation_legend = show_annotation_legend
  ))
}

make_arm_scatter_plot <- function(infer_binned, numbat_binned, bins_df) {
  common_cells <- intersect(colnames(infer_binned), colnames(numbat_binned))
  infer_arm <- arm_mean_matrix(infer_binned[, common_cells, drop = FALSE], bins_df)
  numbat_arm <- arm_mean_matrix(numbat_binned[, common_cells, drop = FALSE], bins_df)
  common_arms <- intersect(rownames(infer_arm), rownames(numbat_arm))
  scatter_df <- data.frame(
    cell_id = rep(common_cells, each = length(common_arms)),
    arm = rep(common_arms, times = length(common_cells)),
    infercna = as.vector(infer_arm[common_arms, common_cells, drop = FALSE]),
    numbat_raw_expression = as.vector(numbat_arm[common_arms, common_cells, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  scatter_df <- scatter_df[is.finite(scatter_df$infercna) & is.finite(scatter_df$numbat_raw_expression), , drop = FALSE]
  if (nrow(scatter_df) < 3) {
    return(ggplot() + theme_void() + ggtitle("Too few arm values"))
  }
  rho <- suppressWarnings(cor(scatter_df$infercna, scatter_df$numbat_raw_expression, method = "spearman", use = "complete.obs"))
  ggplot(scatter_df, aes(.data$infercna, .data$numbat_raw_expression)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
    geom_point(color = "#2B8CBE", alpha = 0.12, size = 0.18) +
    geom_smooth(method = "lm", se = FALSE, color = "#B2182B", linewidth = 0.35) +
    labs(
      title = "Arm-level agreement",
      subtitle = paste0("cells x arms=", nrow(scatter_df), ", Spearman rho=", signif(rho, 3)),
      x = "InferCNA arm mean",
      y = "Numbat raw expression-roll arm mean"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_overlap_plot <- function(meta_plot) {
  plot_df <- meta_plot %>%
    count(.data$raw_expr_cluster, .data$infercna_subclone, name = "n") %>%
    group_by(.data$raw_expr_cluster) %>%
    mutate(frac = .data$n / sum(.data$n)) %>%
    ungroup()
  ggplot(plot_df, aes(.data$raw_expr_cluster, .data$frac, fill = .data$infercna_subclone)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.78) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = complete_palette(make_palette(plot_df$infercna_subclone, "Set2"), plot_df$infercna_subclone, "Set2")) +
    labs(title = "Raw expression cluster overlap", x = "Numbat raw cluster", y = "InferCNA fraction", fill = "InferCNA") +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      legend.key.size = unit(2.5, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

render_page <- function(page_title, left, right, overlap_plot, scatter_plot) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 3,
    ncol = 3,
    heights = unit.c(unit(0.35, "in"), unit(1, "null"), unit(1, "null")),
    widths = unit(c(4.2, 4.2, 2.3), "null")
  )))
  grid.text(page_title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3),
            gp = gpar(fontsize = 13, fontface = "bold"))
  pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 1))
  grid.draw(left)
  popViewport()
  pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 2))
  grid.draw(right)
  popViewport()
  print(overlap_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(scatter_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
  popViewport()
}

summary_rows <- list()
cell_cluster_rows <- list()

pdf(combined_pdf, width = 21, height = 11.5, useDingbats = FALSE)
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample[i]
  numbat_dir <- manifest$numbat_dir[i]
  message("Processing ", sample_id)

  gexp_file <- file.path(numbat_dir, "gexp_roll_wide.tsv.gz")
  map_file <- file.path(numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  infer_file <- file.path("by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  missing <- c(gexp_file, map_file, infer_file)[!file.exists(c(gexp_file, map_file, infer_file))]
  if (length(missing) > 0) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nMissing inputs:\n", paste(missing, collapse = "\n")), gp = gpar(fontsize = 13))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_input", missing = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
    next
  }

  done_status <- parse_done_status(sample_id, numbat_dir)
  infer_mat_full <- readRDS(infer_file)
  if (!is.matrix(infer_mat_full)) infer_mat_full <- as.matrix(infer_mat_full)

  gexp_header <- names(fread(gexp_file, nrows = 0))
  numbat_genes <- setdiff(gexp_header, "cell")
  common_genes <- intersect(intersect(numbat_genes, rownames(infer_mat_full)), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE] %>%
    arrange(.data$chromosome, .data$start)
  common_genes <- go$gene
  if (length(common_genes) < 100) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nToo few common ordered genes: ", length(common_genes)), gp = gpar(fontsize = 13))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "too_few_common_genes", n_common_genes = length(common_genes), stringsAsFactors = FALSE)
    next
  }

  gexp_dt <- fread(gexp_file, select = c("cell", common_genes))
  cell_map <- fread(map_file)
  gexp_dt <- gexp_dt %>%
    left_join(cell_map %>% select(.data$raw_barcode, .data$numbat_cell, .data$cell_id, .data$malignancy), by = c("cell" = "numbat_cell"))
  missing_cell <- is.na(gexp_dt$cell_id)
  if (any(missing_cell)) {
    fallback <- cell_map$cell_id[match(gexp_dt$cell[missing_cell], cell_map$raw_barcode)]
    gexp_dt$cell_id[missing_cell] <- fallback
  }
  gexp_dt <- gexp_dt[!is.na(gexp_dt$cell_id), ]
  if (nrow(gexp_dt) < 20) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nToo few mapped Numbat cells: ", nrow(gexp_dt)), gp = gpar(fontsize = 13))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "too_few_mapped_numbat_cells", n_mapped_cells = nrow(gexp_dt), stringsAsFactors = FALSE)
    next
  }

  infer_lookup <- resolve_infer_columns(sample_id, gexp_dt$cell_id, gexp_dt$cell, colnames(infer_mat_full))
  keep_cells <- !is.na(infer_lookup)
  gexp_dt <- gexp_dt[keep_cells, ]
  infer_lookup <- infer_lookup[keep_cells]
  if (nrow(gexp_dt) < 20) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nToo few matched Numbat/InferCNA cells: ", nrow(gexp_dt)), gp = gpar(fontsize = 13))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "too_few_matched_cells", n_matched_cells = nrow(gexp_dt), stringsAsFactors = FALSE)
    next
  }

  gexp_df <- as.data.frame(gexp_dt)
  numbat_cell_mat <- as.matrix(gexp_df[, common_genes, drop = FALSE])
  storage.mode(numbat_cell_mat) <- "numeric"
  rownames(numbat_cell_mat) <- gexp_dt$cell_id
  numbat_gene_mat <- t(numbat_cell_mat)
  infer_gene_mat <- as.matrix(infer_mat_full[common_genes, infer_lookup, drop = FALSE])
  colnames(infer_gene_mat) <- gexp_dt$cell_id
  rownames(infer_gene_mat) <- common_genes
  rm(infer_mat_full, numbat_cell_mat)
  gc()

  finite_rows <- rowSums(is.finite(numbat_gene_mat)) == ncol(numbat_gene_mat) &
    rowSums(is.finite(infer_gene_mat)) == ncol(infer_gene_mat)
  numbat_gene_mat <- numbat_gene_mat[finite_rows, , drop = FALSE]
  infer_gene_mat <- infer_gene_mat[finite_rows, , drop = FALSE]
  go <- go[finite_rows, , drop = FALSE]
  if (nrow(numbat_gene_mat) < 100 || ncol(numbat_gene_mat) < 20) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nToo few finite matrix values."), gp = gpar(fontsize = 13))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "too_few_finite_values", stringsAsFactors = FALSE)
    next
  }

  numbat_binned <- bin_gene_matrix(numbat_gene_mat, go)
  infer_binned <- bin_gene_matrix(infer_gene_mat, go)
  raw_cluster <- cluster_raw_expression(numbat_binned$mat)
  raw_silhouette <- attr(raw_cluster, "silhouette")
  raw_cluster_sizes <- paste(paste0(names(table(raw_cluster)), "=", as.integer(table(raw_cluster))), collapse = ";")

  meta_plot_all <- data.frame(
    cell_id = colnames(numbat_gene_mat),
    sample = sample_id,
    raw_barcode = gexp_dt$cell,
    raw_expr_cluster = raw_cluster[colnames(numbat_gene_mat)],
    raw_expr_k2_silhouette = raw_silhouette,
    numbat_status = done_status$status,
    numbat_status_short = ifelse(isTRUE(done_status$terminal_no_subclone), "NoFinalClone", "Success"),
    stringsAsFactors = FALSE
  )
  clone_labels <- read_clone_labels(sample_id, numbat_dir, cell_map)
  meta_plot_all$numbat_clone <- clone_labels[meta_plot_all$cell_id]
  meta_plot_all$numbat_clone[is.na(meta_plot_all$numbat_clone) | !nzchar(meta_plot_all$numbat_clone)] <- "No accepted clone"
  if (!is.null(infer_cells)) {
    meta_plot_all <- meta_plot_all %>%
      left_join(infer_cells, by = c("cell_id", "sample"))
  } else {
    meta_plot_all$infercna_subclone <- NA_character_
    meta_plot_all$state_label <- NA_character_
    meta_plot_all$top_mp_label <- NA_character_
  }
  meta_plot_all <- meta_plot_all %>%
    mutate(
      infercna_subclone = ifelse(is.na(.data$infercna_subclone) | !nzchar(.data$infercna_subclone), "InferCNA unassigned", .data$infercna_subclone),
      state_label = ifelse(is.na(.data$state_label) | !nzchar(.data$state_label), "NA", .data$state_label),
      top_mp_label = ifelse(is.na(.data$top_mp_label) | !nzchar(.data$top_mp_label), "NA", .data$top_mp_label)
    )

  set.seed(42)
  plot_cells <- sample_cells_for_plot(meta_plot_all$cell_id, setNames(meta_plot_all$raw_expr_cluster, meta_plot_all$cell_id))
  plot_cells <- order_cells_by_numbat(numbat_binned$mat, plot_cells)
  meta_plot <- meta_plot_all[match(plot_cells, meta_plot_all$cell_id), , drop = FALSE]
  rownames(meta_plot) <- meta_plot$cell_id

  infer_arm <- arm_mean_matrix(infer_binned$mat, infer_binned$bins)
  numbat_arm <- arm_mean_matrix(numbat_binned$mat, numbat_binned$bins)
  common_arms <- intersect(rownames(infer_arm), rownames(numbat_arm))
  arm_cor <- suppressWarnings(cor(
    as.vector(infer_arm[common_arms, , drop = FALSE]),
    as.vector(numbat_arm[common_arms, , drop = FALSE]),
    method = "spearman",
    use = "complete.obs"
  ))
  ari_raw_vs_infer <- adjusted_rand(meta_plot_all$raw_expr_cluster, meta_plot_all$infercna_subclone)

  page_title <- paste0(
    sample_id,
    " | matched cells=", ncol(numbat_gene_mat),
    " | genes=", nrow(numbat_gene_mat),
    " | raw k2: ", raw_cluster_sizes,
    " | arm rho=", signif(arm_cor, 3)
  )
  if (isTRUE(done_status$terminal_no_subclone)) {
    page_title <- paste0(page_title, " | Numbat final: no accepted clone")
  }

  left <- make_heatmap_grob(
    numbat_gene_mat[, plot_cells, drop = FALSE],
    go,
    meta_plot,
    "Expression magnitude",
    "Numbat raw expression-roll matrix (gexp_roll_wide)",
    0.8,
    TRUE
  )
  right <- make_heatmap_grob(
    infer_gene_mat[, plot_cells, drop = FALSE],
    go,
    meta_plot,
    "InferCNA",
    "InferCNA unfiltered expression-CNA matrix (_outs.rds)",
    0.15,
    FALSE
  )
  overlap_plot <- make_overlap_plot(meta_plot_all)
  scatter_plot <- make_arm_scatter_plot(infer_binned$mat, numbat_binned$mat, infer_binned$bins)

  render_page(page_title, left, right, overlap_plot, scatter_plot)

  sample_pdf <- file.path(per_sample_dir, paste0("Auto_", sample_id, "_raw_expression_infercna_concordance.pdf"))
  pdf(sample_pdf, width = 21, height = 11.5, useDingbats = FALSE)
  render_page(page_title, left, right, overlap_plot, scatter_plot)
  dev.off()

  sample_png <- file.path(per_sample_dir, paste0("Auto_", sample_id, "_raw_expression_infercna_concordance.png"))
  png(sample_png, width = 4200, height = 2300, res = 200)
  render_page(page_title, left, right, overlap_plot, scatter_plot)
  dev.off()

  summary_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "ok",
    numbat_status = done_status$status,
    terminal_no_subclone = isTRUE(done_status$terminal_no_subclone),
    numbat_iteration = done_status$iteration,
    n_matched_cells = ncol(numbat_gene_mat),
    n_plot_cells = length(plot_cells),
    n_common_genes = nrow(numbat_gene_mat),
    raw_expr_k2_sizes = raw_cluster_sizes,
    raw_expr_k2_silhouette = raw_silhouette,
    infercna_subclone_count = length(unique(meta_plot_all$infercna_subclone)),
    numbat_final_clone_count = length(setdiff(unique(meta_plot_all$numbat_clone), "No accepted clone")),
    raw_expr_vs_infercna_ari = ari_raw_vs_infer,
    arm_spearman_rho = arm_cor,
    per_sample_pdf = normalizePath(sample_pdf, mustWork = TRUE),
    per_sample_png = normalizePath(sample_png, mustWork = TRUE),
    stringsAsFactors = FALSE
  )

  cell_cluster_rows[[sample_id]] <- meta_plot_all %>%
    select(
      .data$sample,
      .data$cell_id,
      .data$raw_barcode,
      .data$raw_expr_cluster,
      .data$numbat_clone,
      .data$numbat_status,
      .data$numbat_status_short,
      .data$infercna_subclone,
      .data$state_label,
      .data$top_mp_label
    )

  rm(numbat_gene_mat, infer_gene_mat, numbat_binned, infer_binned, left, right)
  gc()
}
dev.off()

summary_tbl <- bind_rows(summary_rows)
cell_cluster_tbl <- bind_rows(cell_cluster_rows)
fwrite(summary_tbl, summary_csv)
fwrite(cell_cluster_tbl, cell_cluster_csv)

run_end <- Sys.time()
writeLines(
  c(
    paste0("script\tanalysis/cnv/scatlas_numbat_raw_expression_concordance.R"),
    paste0("started\t", run_start),
    paste0("finished\t", run_end),
    paste0("elapsed_minutes\t", round(as.numeric(difftime(run_end, run_start, units = "mins")), 3)),
    paste0("sample_arg\t", sample_arg),
    paste0("max_plot_cells\t", max_plot_cells),
    paste0("gene_bin_size\t", gene_bin_size),
    paste0("combined_pdf\t", normalizePath(combined_pdf, mustWork = TRUE)),
    paste0("summary_csv\t", normalizePath(summary_csv, mustWork = TRUE)),
    paste0("cell_cluster_csv\t", normalizePath(cell_cluster_csv, mustWork = TRUE)),
    paste0("n_samples_ok\t", sum(summary_tbl$status == "ok", na.rm = TRUE)),
    paste0("n_samples_total\t", nrow(summary_tbl))
  ),
  run_summary_txt
)

message("Wrote combined PDF: ", combined_pdf)
message("Wrote summary table: ", summary_csv)
message("Wrote cell cluster table: ", cell_cluster_csv)
