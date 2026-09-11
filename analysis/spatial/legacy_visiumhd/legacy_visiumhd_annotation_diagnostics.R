#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_annotation_diagnostics.R
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: Visium HD count H5 files, canonical annotation CSVs, and optional
#     InferCNA/malignancy CSVs from ref_outs/visium_hd_outs/.
#   Outputs: ref_outs/visium_hd_outs/{intermediate,figures,tables,logs}/.
#   Cache/replot: cached UMAP coordinates are reused unless
#     SCREF_FORCE_REBUILD=TRUE; plots can therefore be regenerated cheaply.
#   Run: Rscript analysis/spatial/visiumhd_annotation_diagnostics.R --mode <mode>
#     --inputs <...> --sample-names <...> --annotation-dir <...> --output-dir <...>
#   Environment: dmtcp.
####################

####################
# Per-sample spatial and UMAP diagnostics for the canonical annotation and
# malignancy tiers. Space Ranger UMAPs are used when present; otherwise a
# Seurat UMAP is computed once and cached under the live output directory.
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(patchwork)
})

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

parse_cli <- function(args) {
  out <- list(inputs = character(), sample_names = character())
  multi <- c("inputs", "sample-names")
  allowed <- c(multi, "mode", "annotation-dir", "malignancy-dir", "output-dir")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% allowed) stop("Unknown argument: ", args[[i]])
    i <- i + 1L
    if (i > length(args) || startsWith(args[[i]], "--")) stop("Missing value for --", key)
    if (key %in% multi) {
      start <- i
      while (i <= length(args) && !startsWith(args[[i]], "--")) i <- i + 1L
      out[[gsub("-", "_", key)]] <- args[start:(i - 1L)]
    } else {
      out[[gsub("-", "_", key)]] <- args[[i]]
      i <- i + 1L
    }
  }
  out
}

celltype_colours <- c(
  epithelial = "#E41A1C", fibroblast = "#8C564B", endothelial = "#1F78B4",
  macrophage = "#FF7F00", mast = "#A65628", `t.cell` = "#4DAF4A",
  `b.cell` = "#377EB8", `nk.cell` = "#984EA3", plasma = "#F781BF",
  dendritic = "#17BECF", lymph = "#6BAED6", erythrocyte = "#BDBDBD",
  keratinocyte = "#FFD92F", neutrophil = "#66C2A5", Unknown = "#BDBDBD",
  unresolved = "#BDBDBD"
)
malignancy_colours <- c(
  "Background / non-epithelial" = "#D9D9D9", unresolved = "#80B1D3",
  normal_keratinocyte = "#4DAF4A", keratinocyte_indeterminate = "#FFD92F", malignant_level_1 = "#E41A1C", malignant_level_2 = "#984EA3"
)

read_counts <- function(input_dir, mode) {
  h5 <- file.path(input_dir, if (mode %in% c("binned", "custom")) "filtered_feature_bc_matrix.h5" else "filtered_feature_cell_matrix.h5")
  if (!file.exists(h5)) stop("Missing count matrix: ", h5)
  counts <- Read10X_h5(h5)
  if (is.list(counts)) counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  counts
}

get_umap_coordinates <- function(annotation, input_dir, sample_name, mode, intermediate_dir) {
  cache_path <- file.path(intermediate_dir, paste0("Auto_", sample_name, "_", mode, "_annotation_umap_coords.csv.gz"))
  force_rebuild <- identical(tolower(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")), "true")
  if (file.exists(cache_path) && !force_rebuild) return(fread(cache_path, data.table = FALSE))
  spaceranger_umap <- file.path(input_dir, "analysis", "umap", "gene_expression_2_components", "projection.csv")
  if (file.exists(spaceranger_umap)) {
    coords <- fread(spaceranger_umap, data.table = FALSE)
    colnames(coords)[match("Barcode", colnames(coords))] <- "barcode"
    colnames(coords)[match("UMAP-1", colnames(coords))] <- "UMAP_1"
    colnames(coords)[match("UMAP-2", colnames(coords))] <- "UMAP_2"
    coords <- coords[, c("barcode", "UMAP_1", "UMAP_2"), drop = FALSE]
  } else {
    counts <- read_counts(input_dir, mode)
    barcodes <- intersect(colnames(counts), annotation$barcode)
    if (length(barcodes) < 30L) stop("Too few annotated observations to calculate UMAP for ", sample_name)
    obj <- CreateSeuratObject(counts = counts[, barcodes, drop = FALSE], min.cells = 0, min.features = 0)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    npcs <- min(30L, length(VariableFeatures(obj)), ncol(obj) - 1L)
    if (npcs < 2L) stop("Too few principal components to calculate UMAP for ", sample_name)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)
    obj <- RunUMAP(obj, dims = seq_len(npcs), seed.use = 0, verbose = FALSE)
    coords <- as.data.frame(Embeddings(obj, "umap"))
    coords$barcode <- rownames(coords)
    colnames(coords)[1:2] <- c("UMAP_1", "UMAP_2")
    coords <- coords[, c("barcode", "UMAP_1", "UMAP_2"), drop = FALSE]
  }
  fwrite(coords, cache_path, compress = "gzip")
  coords
}

annotation_plot <- function(data, x, y, sample_name, mode, coordinate_title) {
  data$Auto_annotation_celltype <- as.character(data$Auto_annotation_celltype)
  ####################
  # Retain explicit legend entries for observed multi-label calls rather than
  # allowing an unmatched named-vector lookup to turn their names into NA.
  observed <- unique(c(names(celltype_colours), data$Auto_annotation_celltype))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  palette <- setNames(rep("#636363", length(observed)), observed)
  matched <- intersect(observed, names(celltype_colours))
  palette[matched] <- celltype_colours[matched]
  ####################
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], colour = Auto_annotation_celltype)) +
    geom_point(size = if (mode %in% c("binned", "custom")) 0.08 else 0.18, alpha = 0.75) +
    scale_colour_manual(values = palette, na.value = "#BDBDBD") +
    guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    coord_fixed() +
    labs(title = paste(sample_name, mode, "annotation"), subtitle = coordinate_title, colour = "Cell type") +
    theme_classic(base_size = 13) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right")
}

malignancy_plot <- function(data, x, y, sample_name, mode, coordinate_title) {
  data$plot_malignancy <- "Background / non-epithelial"
  epithelial <- data$Auto_annotation_keep_epithelial %in% TRUE
  data$plot_malignancy[epithelial] <- "unresolved"
  data$plot_malignancy[epithelial & data$Auto_malignancy == "normal_keratinocyte"] <- "normal_keratinocyte"
  data$plot_malignancy[epithelial & data$Auto_malignancy == "keratinocyte_indeterminate"] <- "keratinocyte_indeterminate"
  data$plot_malignancy[epithelial & data$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")] <- data$Auto_malignancy[epithelial & data$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2")]
  data$plot_malignancy <- factor(data$plot_malignancy, levels = names(malignancy_colours))
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], colour = plot_malignancy)) +
    geom_point(size = if (mode %in% c("binned", "custom")) 0.08 else 0.18, alpha = 0.75) +
    scale_colour_manual(values = malignancy_colours, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    coord_fixed() +
    labs(title = paste(sample_name, mode, "malignancy"), subtitle = coordinate_title, colour = NULL) +
    theme_classic(base_size = 13) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right")
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (length(args$inputs) == 0L || length(args$inputs) != length(args$sample_names) || is.null(args$mode) || is.null(args$annotation_dir) || is.null(args$output_dir)) {
  stop("Required: --mode --inputs --sample-names --annotation-dir --output-dir")
}
####################
# Custom diagnostics reuse binned geometry and cached binned expression UMAPs.
####################
####################
# Spatial annotations use the segmented count/coordinate branch.
if (!args$mode %in% c("binned", "segmented", "custom", "spatial")) {
  stop("--mode must be binned, segmented, custom, or spatial")
}
####################
####################
intermediate_dir <- file.path(args$output_dir, "intermediate")
figure_dir <- file.path(args$output_dir, "figures", "annotation_diagnostics")
table_dir <- file.path(args$output_dir, "tables")
log_dir <- file.path(args$output_dir, "logs")
for (dir_path in c(intermediate_dir, figure_dir, table_dir, log_dir)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()

for (i in seq_along(args$inputs)) {
  sample_name <- args$sample_names[[i]]
  annotation_path <- file.path(args$annotation_dir, paste0("Auto_", sample_name, "_", args$mode, "_cell_annotations.csv.gz"))
  if (!file.exists(annotation_path)) stop("Missing annotation table: ", annotation_path)
  annotation <- fread(annotation_path, data.table = FALSE)
  annotation$barcode <- as.character(annotation$barcode)
  annotation$Auto_annotation_keep_epithelial <- as.logical(annotation$Auto_annotation_keep_epithelial)
  annotation$Auto_malignancy <- NA_character_
  if (!is.null(args$malignancy_dir)) {
    malignancy_path <- file.path(args$malignancy_dir, "tables", paste0("Auto_", sample_name, "_", args$mode, "_infercna_cells.csv.gz"))
    if (file.exists(malignancy_path)) {
      malignancy <- fread(malignancy_path, data.table = FALSE)
      malignancy$barcode <- as.character(malignancy$barcode)
      annotation <- annotation %>% left_join(malignancy %>% select(barcode, Auto_malignancy), by = "barcode", suffix = c("", ".from_cna"))
      if ("Auto_malignancy.from_cna" %in% colnames(annotation)) annotation$Auto_malignancy <- annotation$Auto_malignancy.from_cna
    }
  }
  coords <- get_umap_coordinates(annotation, args$inputs[[i]], sample_name, args$mode, intermediate_dir)
  plot_data <- annotation %>% left_join(coords, by = "barcode")
  spatial_data <- plot_data[is.finite(plot_data$pxl_col_in_fullres) & is.finite(plot_data$pxl_row_in_fullres), , drop = FALSE]
  umap_data <- plot_data[is.finite(plot_data$UMAP_1) & is.finite(plot_data$UMAP_2), , drop = FALSE]
  spatial <- annotation_plot(spatial_data, "pxl_col_in_fullres", "pxl_row_in_fullres", sample_name, args$mode, "Spatial coordinates") + scale_y_reverse()
  spatial_malignancy <- malignancy_plot(spatial_data, "pxl_col_in_fullres", "pxl_row_in_fullres", sample_name, args$mode, "Spatial coordinates") + scale_y_reverse()
  umap <- annotation_plot(umap_data, "UMAP_1", "UMAP_2", sample_name, args$mode, "Expression UMAP")
  umap_malignancy <- malignancy_plot(umap_data, "UMAP_1", "UMAP_2", sample_name, args$mode, "Expression UMAP")
  spatial_combined <- spatial + spatial_malignancy + plot_layout(guides = "collect")
  umap_combined <- umap + umap_malignancy + plot_layout(guides = "collect")
  ggsave(file.path(figure_dir, paste0("Auto_", sample_name, "_", args$mode, "_annotation_spatial.pdf")), spatial_combined, width = 16, height = 8, useDingbats = FALSE)
  ggsave(file.path(figure_dir, paste0("Auto_", sample_name, "_", args$mode, "_annotation_spatial.png")), spatial_combined, width = 16, height = 8, dpi = 300)
  ggsave(file.path(figure_dir, paste0("Auto_", sample_name, "_", args$mode, "_annotation_umap.pdf")), umap_combined, width = 16, height = 8, useDingbats = FALSE)
  ggsave(file.path(figure_dir, paste0("Auto_", sample_name, "_", args$mode, "_annotation_umap.png")), umap_combined, width = 16, height = 8, dpi = 300)
  summary_rows[[sample_name]] <- data.frame(sample = sample_name, mode = args$mode, n_spatial = nrow(spatial_data), n_umap = nrow(umap_data), stringsAsFactors = FALSE)
}
summary_df <- dplyr::bind_rows(summary_rows)
write.csv(summary_df, file.path(table_dir, paste0("Auto_visiumhd_", args$mode, "_annotation_diagnostics_summary.csv")), row.names = FALSE)
summary_dir <- file.path(WD, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(summary_dir, paste0("visiumhd_", args$mode, "_annotation_diagnostics_summary.csv")), row.names = FALSE)
writeLines(c(paste0("mode=", args$mode), paste0("samples=", paste(args$sample_names, collapse = ",")), paste0("end=", format(Sys.time(), tz = "Europe/London"))), file.path(log_dir, paste0("Auto_visiumhd_", args$mode, "_annotation_diagnostics_run_summary.txt")))
