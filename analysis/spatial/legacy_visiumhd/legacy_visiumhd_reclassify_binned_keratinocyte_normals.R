#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_reclassify_binned_keratinocyte_normals.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: cached binned InferCNA cell tables and canonical segmented
#     annotation tables under ref_outs/visium_hd_outs/.
#   Outputs: corrected binned InferCNA cell tables/scatters and summary tables.
#   Cache/replot: does not rerun InferCNA; uses cached CNA scatter coordinates.
#   Run: Rscript analysis/spatial/visiumhd_reclassify_binned_keratinocyte_normals.R
#     --samples SUR1231 FFPEA1 FFPED1 --annotation-dir <tables> --malignancy-dir <dir>
#   Environment: dmtcp.
####################

####################
# Reclassify RCTD epithelial bins as normal in the malignancy layer only when
# their matched segmented cells are keratinocyte-dominant. The RCTD annotation
# remains epithelial because keratinocytes are absent from its reference.
####################
suppressPackageStartupMessages({
  library(RANN)
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
    if (!key %in% c("samples", "annotation-dir", "malignancy-dir")) stop("Unknown argument: ", args[[i]])
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

keratinocyte_projection <- function(binned, segmented) {
  result <- data.frame(
    Auto_segmented_keratinocyte_n_cells = integer(nrow(binned)),
    Auto_segmented_keratinocyte_fraction = numeric(nrow(binned)),
    Auto_segmented_keratinocyte_nearest_distance = rep(NA_real_, nrow(binned))
  )
  required <- c("pxl_col_in_fullres", "pxl_row_in_fullres", "is_epithelial_target")
  if (length(setdiff(required, colnames(binned)))) return(result)
  binned$pxl_col_in_fullres <- suppressWarnings(as.numeric(binned$pxl_col_in_fullres))
  binned$pxl_row_in_fullres <- suppressWarnings(as.numeric(binned$pxl_row_in_fullres))
  segmented$pxl_col_in_fullres <- suppressWarnings(as.numeric(segmented$pxl_col_in_fullres))
  segmented$pxl_row_in_fullres <- suppressWarnings(as.numeric(segmented$pxl_row_in_fullres))
  binned$is_epithelial_target <- as_logical(binned$is_epithelial_target)
  bin_idx <- which(binned$is_epithelial_target & is.finite(binned$pxl_col_in_fullres) & is.finite(binned$pxl_row_in_fullres))
  cell_idx <- which(is.finite(segmented$pxl_col_in_fullres) & is.finite(segmented$pxl_row_in_fullres))
  if (!length(bin_idx) || !length(cell_idx)) return(result)
  nearest <- RANN::nn2(
    data = as.matrix(binned[bin_idx, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    query = as.matrix(segmented[cell_idx, c("pxl_col_in_fullres", "pxl_row_in_fullres")]),
    k = 1L
  )
  matched <- data.frame(
    binned_idx = bin_idx[nearest$nn.idx[, 1L]],
    distance = nearest$nn.dists[, 1L],
    keratinocyte = segmented$Auto_annotation_celltype[cell_idx] == "keratinocyte"
  )
  matched <- matched[is.finite(matched$distance) & matched$distance <= 25, , drop = FALSE]
  if (!nrow(matched)) return(result)
  composition <- aggregate(cbind(n_cells = rep(1L, nrow(matched)), n_keratinocyte = as.integer(matched$keratinocyte)), by = list(binned_idx = matched$binned_idx), FUN = sum)
  result$Auto_segmented_keratinocyte_n_cells[composition$binned_idx] <- composition$n_cells
  result$Auto_segmented_keratinocyte_fraction[composition$binned_idx] <- composition$n_keratinocyte / composition$n_cells
  distance <- aggregate(distance ~ binned_idx, data = matched, FUN = median)
  result$Auto_segmented_keratinocyte_nearest_distance[distance$binned_idx] <- distance$distance
  result
}

plot_scatter <- function(data, sample_name, figure_dir) {
  data <- data[is.finite(data$cna.signal) & is.finite(data$cna.cor), , drop = FALSE]
  data$plot_group <- ifelse(data$is_reference, "Reference", as.character(data$Auto_malignancy))
  data$plot_group[data$plot_group == "unresolved"] <- "Unresolved epithelial"
  data$plot_group <- factor(data$plot_group, levels = c("Reference", "normal_keratinocyte", "Unresolved epithelial", "malignant_level_1", "malignant_level_2"))
  p <- ggplot(data, aes(cna.signal, cna.cor, colour = plot_group)) +
    geom_point(size = 0.75, alpha = 0.65) +
    scale_colour_manual(values = c(Reference = "#9E9E9E", normal_keratinocyte = "#FFD92F", `Unresolved epithelial` = "#80B1D3", malignant_level_1 = "#E41A1C", malignant_level_2 = "#984EA3"), drop = FALSE) +
    scale_x_continuous(labels = scales::scientific) +
    labs(title = paste0(sample_name, " binned: malignancy after keratinocyte-normal correction"), x = "CNA signal", y = "CNA correlation", colour = NULL) +
    theme_classic(base_size = 14) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  prefix <- file.path(figure_dir, paste0("Auto_", sample_name, "_binned_infercna_scatter"))
  ggsave(paste0(prefix, ".pdf"), p, width = 8.5, height = 6.5, useDingbats = FALSE)
  ggsave(paste0(prefix, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!length(args$samples) || is.null(args$annotation_dir) || is.null(args$malignancy_dir)) stop("Required: --samples --annotation-dir --malignancy-dir")
annotation_dir <- normalizePath(args$annotation_dir, mustWork = TRUE)
malignancy_dir <- normalizePath(args$malignancy_dir, mustWork = TRUE)
summary_rows <- list()

for (sample_name in args$samples) {
  binned_path <- file.path(malignancy_dir, "tables", paste0("Auto_", sample_name, "_binned_infercna_cells.csv.gz"))
  segmented_path <- file.path(annotation_dir, paste0("Auto_", sample_name, "_segmented_cell_annotations.csv.gz"))
  if (!file.exists(binned_path) || !file.exists(segmented_path)) stop("Missing binned CNA or segmented annotation input for ", sample_name)
  binned <- read.csv(binned_path, stringsAsFactors = FALSE, check.names = FALSE)
  segmented <- read.csv(segmented_path, stringsAsFactors = FALSE, check.names = FALSE)
  projection <- keratinocyte_projection(binned, segmented)
  binned[, colnames(projection)] <- projection
  binned$Auto_malignancy_before_keratinocyte_exclusion <- binned$Auto_malignancy
  normal_keratinocyte <- as_logical(binned$is_epithelial_target) & binned$Auto_segmented_keratinocyte_n_cells >= 1L & binned$Auto_segmented_keratinocyte_fraction >= 0.5
  binned$Auto_malignancy[normal_keratinocyte] <- "normal_keratinocyte"
  binned$Auto_malignancy_evidence[normal_keratinocyte] <- "segmented_keratinocyte_dominant"
  binned$Auto_malignant <- binned$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")
  con <- gzfile(binned_path, open = "wt")
  write.csv(binned, con, row.names = FALSE)
  close(con)
  plot_scatter(binned, sample_name, file.path(malignancy_dir, "figures"))
  targets <- binned[as_logical(binned$is_epithelial_target), , drop = FALSE]
  summary_rows[[sample_name]] <- data.frame(
    sample = sample_name,
    mode = "binned",
    n_normal_keratinocyte = sum(targets$Auto_malignancy == "normal_keratinocyte", na.rm = TRUE),
    n_keratinocyte_normal_from_level_1 = sum(targets$Auto_malignancy == "normal_keratinocyte" & targets$Auto_malignancy_before_keratinocyte_exclusion == "malignant_level_1", na.rm = TRUE),
    n_malignant_level_1 = sum(targets$Auto_malignancy == "malignant_level_1", na.rm = TRUE),
    n_malignant_level_2 = sum(targets$Auto_malignancy == "malignant_level_2", na.rm = TRUE),
    n_malignant = sum(targets$Auto_malignant, na.rm = TRUE),
    pct_malignant_epithelial = 100 * mean(targets$Auto_malignant, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(malignancy_dir, "tables", "Auto_visiumhd_binned_infercna_malignancy_summary.csv")
existing <- if (file.exists(summary_path)) read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame(sample = character())
for (column in setdiff(colnames(summary_df), c("sample", "mode"))) {
  existing[[column]] <- summary_df[[column]][match(existing$sample, summary_df$sample)]
}
write.csv(existing, summary_path, row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(existing, file.path(summary_dir, "visiumhd_binned_infercna_malignancy_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(malignancy_dir, "tables", "Auto_visiumhd_binned_keratinocyte_normal_reclassification_summary.csv"), row.names = FALSE)
writeLines(c(paste0("samples=", paste(args$samples, collapse = ",")), paste0("end=", format(Sys.time(), tz = "Europe/London"))), file.path(malignancy_dir, "logs", "Auto_visiumhd_binned_keratinocyte_normal_reclassification_run_summary.txt"))
