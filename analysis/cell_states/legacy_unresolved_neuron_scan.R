####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_unresolved_neuron_scan.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_unresolved_neuron_scan.R
# Scan per-sample annotation outputs for neuron-like unresolved clusters.
#
# Input:
#   ref_outs/by_samples/<sample>/<sample>_anno.rds
#
# Output:
#   ref_outs/Auto_unresolved_neuron_scan/Auto_unresolved_neuron_cluster_summary.csv
#   ref_outs/Auto_unresolved_neuron_scan/Auto_unresolved_neuron_sample_flags.csv
#   ref_outs/Auto_unresolved_neuron_scan/Auto_unresolved_neuron_samples_plots.pdf
#   updates/new_updates/summaries/Auto_unresolved_neuron_scan_summary.csv
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(Matrix)
####################
library(readxl)
####################

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

pan_neuronal_markers <- c(
  "TUBB3",
  "STMN2",
  "UCHL1",
  "SNAP25",
  "SYT1",
  "MAP2",
  "RBFOX3",
  "NEFM"
)

# Increased threshold a little bit
min_pct_expr <- 20
min_avg_expr <- 0.4
min_markers_positive <- 3

celltype_map <- c(
  "Fibroblast" = "fibroblast", "Macrophage" = "macrophage", "Mast" = "mast",
  "B cell" = "b.cell", "T cell" = "t.cell", "Dendritic" = "dendritic",
  "Endothelial" = "endothelial", "Epithelial" = "epithelial",
  "NK cell" = "nk.cell", "Plasma" = "plasma"
)

celltype_order <- c(
  "epithelial", "endothelial", "fibroblast", "macrophage",
  "b.cell", "plasma", "t.cell", "nk.cell", "mast", "dendritic"
)

celltype_display <- c(
  epithelial = "Epithelial", endothelial = "Endothelial",
  fibroblast = "Fibroblast", macrophage = "Macrophage",
  "b.cell" = "B cell", plasma = "Plasma",
  "t.cell" = "T cell", "nk.cell" = "NK cell",
  mast = "Mast", dendritic = "Dendritic"
)

celltype_colors <- c(
  epithelial = "#D73027", endothelial = "#2C7BB6",
  fibroblast = "#A6761D", macrophage = "#1B9E77",
  "b.cell" = "#7570B3", plasma = "#E7298A",
  "t.cell" = "#66A61E", "nk.cell" = "#1F78B4",
  mast = "#E6AB02", dendritic = "#6A3D9A"
)

canonical_markers <- list(
  fibroblast = c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"),
  macrophage = c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"),
  mast = c("MS4A2", "CPA3", "TPSB2", "TPSAB1"),
  epithelial = c("KRT7", "MUC1", "KRT19", "EPCAM"),
  t.cell = c("CD3E", "CD3D", "CD2", "CD3G"),
  b.cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  nk.cell = c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"),
  plasma = c("MZB1", "JCHAIN", "DERL3"),
  dendritic = c("CLEC10A", "CCR7", "CD86"),
  endothelial = c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5")
)

get_expr_layer <- function(obj, layer = "data") {
  layer_data <- tryCatch(obj@assays$RNA@layers[[layer]], error = function(e) NULL)
  if (!is.null(layer_data)) {
    if (length(rownames(layer_data)) == 0) rownames(layer_data) <- rownames(obj)
    if (length(colnames(layer_data)) == 0) colnames(layer_data) <- colnames(obj)
    return(layer_data)
  }
  GetAssayData(obj, slot = layer)
}

add_score_matrix <- function(obj, score_mat, prefix) {
  for (ct in celltype_order) {
    obj[[paste0(prefix, ct)]] <- score_mat[colnames(obj), ct]
  }
  obj
}

build_canonical_score_matrix <- function(obj) {
  expr <- get_expr_layer(obj, "data")
  score_mat <- matrix(0, nrow = ncol(obj), ncol = length(celltype_order),
                      dimnames = list(colnames(obj), celltype_order))
  for (ct in celltype_order) {
    genes <- intersect(canonical_markers[[ct]], rownames(expr))
    if (length(genes) == 0) next
    score_mat[, ct] <- as.numeric(Matrix::colMeans(expr[genes, , drop = FALSE]))
  }
  score_mat
}

get_point_size <- function(n_cells) {
  if (n_cells > 40000) 0.35
  else if (n_cells > 20000) 0.8
  else if (n_cells > 10000) 0.9
  else if (n_cells > 5000) 1.0
  else if (n_cells > 2000) 1.2
  else 1.4
}

base_dimplot <- function(obj, group_by, title_text, label = FALSE,
                         palette = NULL, show_legend = FALSE) {
  point_size <- get_point_size(ncol(obj))
  p <- DimPlot(object = obj, group.by = group_by, label = label, repel = label,
               raster = FALSE, pt.size = point_size,
               label.size = 4, combine = FALSE)[[1]]
  if (!is.null(palette)) {
    present_levels <- unique(as.character(obj@meta.data[[group_by]]))
    present_levels <- present_levels[!is.na(present_levels)]
    p <- p + scale_color_manual(
      values = palette,
      breaks = intersect(names(palette), present_levels),
      labels = if (!is.null(names(palette))) {
        vals <- intersect(names(palette), present_levels)
        if (all(vals %in% names(celltype_display))) celltype_display[vals] else vals
      },
      drop = TRUE
    )
  }
  p + ggtitle(title_text) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm"),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(), axis.title = element_blank(),
      panel.grid = element_blank()
    )
}

score_featureplot_sequential <- function(obj, feature_name, title_text, limits) {
  point_size <- get_point_size(ncol(obj))
  coords <- as.data.frame(Embeddings(obj, reduction = "umap")[, 1:2, drop = FALSE])
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  values <- as.numeric(obj@meta.data[[feature_name]])
  plot_df <- cbind(coords, score = values)
  plot_df$score_plot <- pmin(pmax(plot_df$score, limits[1]), limits[2])
  plot_df <- plot_df[order(plot_df$score_plot), , drop = FALSE]

  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(color = "grey88", size = point_size + 0.1, stroke = 0) +
    geom_point(aes(color = score_plot), size = point_size, stroke = 0, alpha = 0.95) +
    scale_color_gradient(
      low = "grey88", high = "#D73027",
      limits = limits, oob = scales::squish, name = NULL,
      guide = guide_colorbar(barwidth = unit(0.35, "cm"), barheight = unit(1.8, "cm"))
    ) +
    ggtitle(title_text) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(), axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.margin = margin(0, 0, 0, -4),
      panel.grid = element_blank()
    )
}

build_expression_page <- function(obj, score_prefix, page_title,
                                  label_column, label_title,
                                  label_palette = celltype_colors) {
  all_vals <- unlist(lapply(celltype_order, function(ct) {
    obj@meta.data[[paste0(score_prefix, ct)]]
  }))
  shared_max <- max(all_vals, na.rm = TRUE)
  shared_max <- max(shared_max, 0.1)
  limits <- c(0, ceiling(shared_max * 10) / 10)

  plots <- list(base_dimplot(obj, "seurat_clusters", "Clusters", label = TRUE))
  for (ct in celltype_order) {
    plots[[length(plots) + 1]] <- score_featureplot_sequential(
      obj, paste0(score_prefix, ct), celltype_display[[ct]], limits = limits
    )
  }
  plots[[length(plots) + 1]] <- base_dimplot(
    obj, label_column, label_title,
    label = FALSE, palette = label_palette, show_legend = TRUE
  )
  wrap_plots(plots, ncol = 4, nrow = 3) +
    plot_annotation(title = page_title)
}

get_annotation_column <- function(obj) {
  if ("celltype_update" %in% colnames(obj@meta.data)) {
    return("celltype_update")
  }
  if ("celltype_final" %in% colnames(obj@meta.data)) {
    return("celltype_final")
  }
  return(NA_character_)
}

is_unresolved_cell <- function(x) {
  x <- tolower(as.character(x))
  grepl("^unresolved", x)
}

sample_dirs <- list.dirs(path = "by_samples", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

out_dir <- "Auto_unresolved_neuron_scan"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("../updates/new_updates/summaries", recursive = TRUE, showWarnings = FALSE)

cluster_summary_all <- list()
sample_flags <- list()

# Open PDF early and print directly to it to save memory
pdf(file.path(out_dir, "Auto_unresolved_neuron_samples_plots.pdf"), width = 18, height = 13)

for (sample in sample_dirs) {
  anno_path <- file.path("by_samples", sample, paste0(sample, "_anno.rds"))
  if (!file.exists(anno_path)) {
    next
  }

  obj <- readRDS(anno_path)
  anno_col <- get_annotation_column(obj)

  if (is.na(anno_col) || !"seurat_clusters" %in% colnames(obj@meta.data)) {
    rm(obj); next
  }

  unresolved_idx <- is_unresolved_cell(obj@meta.data[[anno_col]])
  unresolved_cells <- colnames(obj)[unresolved_idx]

  # Scan ALL clusters in the sample
  cluster_ids <- unique(as.character(obj$seurat_clusters))
  sample_cluster_summary <- list()

  # Define markers present in this sample
  genes_in_obj <- intersect(pan_neuronal_markers, rownames(obj))
  if (length(genes_in_obj) == 0) {
    sample_flags[[sample]] <- data.frame(
      sample = sample,
      has_neuron_like_unresolved_cluster = FALSE,
      n_unresolved_cells = length(unresolved_cells),
      n_neuron_like_clusters = 0,
      neuron_like_clusters = "",
      stringsAsFactors = FALSE
    )
    rm(obj); next
  }

  # Pre-calculate expression for all cells to be faster
  expr_all <- GetAssayData(obj, slot = "data")[genes_in_obj, , drop = FALSE]

  for (cl in cluster_ids) {
    cl_cells <- colnames(obj)[obj$seurat_clusters == cl]
    if (length(cl_cells) < 5) {
      next
    }

    cl_expr <- expr_all[, cl_cells, drop = FALSE]
    pct_expr <- Matrix::rowMeans(cl_expr > 0) * 100
    avg_expr <- Matrix::rowMeans(cl_expr)

    marker_tbl <- data.frame(
      sample = sample,
      cluster = cl,
      marker = names(pct_expr),
      pct_expr = as.numeric(pct_expr),
      avg_expr = as.numeric(avg_expr),
      stringsAsFactors = FALSE
    )

    marker_tbl <- marker_tbl %>%
      mutate(marker_positive = pct_expr >= min_pct_expr & avg_expr >= min_avg_expr)

    n_pos <- sum(marker_tbl$marker_positive)
    marker_tbl$cluster_is_neuron_like <- n_pos >= min_markers_positive

    sample_cluster_summary[[cl]] <- marker_tbl
  }

  if (length(sample_cluster_summary) == 0) {
    post_annotation_sample_flags[[sample]] <- data.frame(
      sample = sample,
      has_neuron_like_unresolved_cluster = FALSE,
      n_post_annotation_unresolved_cells = sum(unresolved_idx),
      n_post_annotation_unresolved_clusters = length(unresolved_clusters),
      n_neuron_like_clusters = 0,
      neuron_like_clusters = "",
      stringsAsFactors = FALSE
    )
    rm(obj); next
  }

  sample_cluster_df <- bind_rows(sample_cluster_summary)
  cluster_summary_all[[sample]] <- sample_cluster_df

  neuron_like_clusters <- sample_cluster_df %>%
    group_by(cluster) %>%
    summarize(cluster_is_neuron_like = any(cluster_is_neuron_like), .groups = "drop") %>%
    filter(cluster_is_neuron_like) %>%
    pull(cluster)

  has_neuron <- length(neuron_like_clusters) > 0

  sample_flags[[sample]] <- data.frame(
    sample = sample,
    has_neuron_like_unresolved_cluster = has_neuron,
    n_unresolved_cells = length(unresolved_cells),
    n_neuron_like_clusters = length(neuron_like_clusters),
    neuron_like_clusters = paste(neuron_like_clusters, collapse = ";"),
    stringsAsFactors = FALSE
  )

  if (has_neuron) {
    message(paste("Flagged sample:", sample, "| neuron-like clusters:", paste(neuron_like_clusters, collapse = ",")))
    
    # Label entire clusters as neuron-like
    obj$auto_neuron_scan_label <- "other_cluster"
    neuron_like_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% neuron_like_clusters]
    obj$auto_neuron_scan_label[neuron_like_cells] <- "neuron_like_cluster"
    
    obj$auto_neuron_scan_label <- factor(
      obj$auto_neuron_scan_label,
      levels = c("other_cluster", "neuron_like_cluster")
    )

    if (!"umap" %in% names(obj@reductions)) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      obj <- ScaleData(obj, verbose = FALSE)
      n_pcs <- min(30, max(5, ncol(obj) - 1))
      obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
      dims_use <- 1:min(20, n_pcs)
      obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
    }

    # --- PAGE 1: Neuronal markers (3x4 panel) ---
    
    # 1. Compute shared scale limits for the 8 markers
    # Use the full list of markers to maintain position
    marker_expr_list <- lapply(pan_neuronal_markers, function(g) {
      if (g %in% rownames(obj)) FetchData(obj, vars = g)[, 1] else rep(0, ncol(obj))
    })
    shared_max_neuron <- max(unlist(marker_expr_list), na.rm = TRUE)
    shared_max_neuron <- max(shared_max_neuron, 0.1)
    limits_neuron <- c(0, ceiling(shared_max_neuron * 10) / 10)

    # 2. Build marker plots for the 8 markers
    marker_plots <- lapply(pan_neuronal_markers, function(g) {
      if (g %in% rownames(obj)) {
        obj@meta.data[[paste0("neuronal_scan_", g)]] <- FetchData(obj, vars = g)[, 1]
        score_featureplot_sequential(obj, paste0("neuronal_scan_", g), g, limits_neuron)
      } else {
        ggplot() + theme_void() + ggtitle(paste(g, "(not in sample)"))
      }
    })

    # 3. Third row, Column 1: auto_neuron_scan_label UMAP
    p_scan_umap <- base_dimplot(
      obj, "auto_neuron_scan_label", "Neuron-like Cluster Scan", 
      label = FALSE, 
      palette = c("other_cluster" = "grey85", "neuron_like_cluster" = "#D55E00"), 
      show_legend = TRUE
    )

    # 4. Third row, Column 2: celltype annotation UMAP
    p_anno_umap <- base_dimplot(
      obj, anno_col, "Celltype Annotation", 
      label = FALSE, 
      palette = celltype_colors, 
      show_legend = TRUE
    )

    # Assemble Page 1
    plots_page1 <- marker_plots # 8 plots
    plots_page1[[9]] <- p_scan_umap
    plots_page1[[10]] <- p_anno_umap
    plots_page1[[11]] <- plot_spacer()
    plots_page1[[12]] <- plot_spacer()

    combined_page1 <- wrap_plots(plots_page1, ncol = 4, nrow = 3) +
      plot_annotation(title = paste0(sample, " | Neuronal Markers"))

    # --- PAGE 2: Canonical marker expression panel ---
    obj <- add_score_matrix(obj, build_canonical_score_matrix(obj), prefix = "canonical_score_")
    combined_page2 <- build_expression_page(
      obj = obj, score_prefix = "canonical_score_",
      page_title = paste0(sample, " | Canonical Marker Expression"),
      label_column = anno_col, label_title = "Celltype Annotation"
    )

    print(combined_page1)
    print(combined_page2)
  }
  
  rm(obj)
  gc()
}

dev.off()

cluster_summary_df <- if (length(cluster_summary_all) > 0) bind_rows(cluster_summary_all) else data.frame()
sample_flags_df <- if (length(sample_flags) > 0) bind_rows(sample_flags) else data.frame()

if (nrow(cluster_summary_df) > 0) {
  write.csv(
    cluster_summary_df,
    file = file.path(out_dir, "Auto_unresolved_neuron_cluster_summary.csv"),
    row.names = FALSE
  )
}

if (nrow(sample_flags_df) > 0) {
  write.csv(
    sample_flags_df,
    file = file.path(out_dir, "Auto_unresolved_neuron_sample_flags.csv"),
    row.names = FALSE
  )
}

summary_df <- sample_flags_df %>%
  summarize(
    n_samples_scanned = n(),
    n_samples_flagged = sum(has_neuron_like_unresolved_cluster, na.rm = TRUE),
    median_unresolved_cells = median(n_unresolved_cells, na.rm = TRUE)
  )

write.csv(
  summary_df,
  file = "../updates/new_updates/summaries/Auto_unresolved_neuron_scan_summary.csv",
  row.names = FALSE
)

####################
# Additional scan: reconstruct the immediate post-Annotation.R calls from
# by-sample <sample>.rds inputs, then scan only clusters called unresolved
# before Expr_filtering.R overwrites <sample>_anno.rds.
####################

combine_marker_scores_post_annotation <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }
  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)
  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

get_excluded_off_cts_post_annotation <- function(home_ct) {
  excluded <- character(0)
  if (grepl("^b.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "plasma", "dendritic")
  }
  if (grepl("^plasma", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "b.cell", "dendritic")
  }
  if (grepl("^t.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "nk.cell", "dendritic")
  }
  if (grepl("^nk.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "t.cell", "dendritic")
  }
  if (grepl("^macrophage", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "dendritic")
  }
  excluded
}

build_post_annotation_inputs <- function(sample_dirs) {
  tmdata_list <- list()
  for (sample in sample_dirs) {
    rds_path <- file.path("by_samples", sample, paste0(sample, ".rds"))
    if (file.exists(rds_path)) {
      tmdata_list[[sample]] <- readRDS(rds_path)
      if (ncol(tmdata_list[[sample]]) < 50) {
        tmdata_list[[sample]]$celltype_update <- rep("unresolved", ncol(tmdata_list[[sample]]))
      }
    }
  }

  tmdata_list <- tmdata_list[vapply(tmdata_list, function(obj) {
    "celltype_initial" %in% colnames(obj@meta.data) &&
      "seurat_clusters" %in% colnames(obj@meta.data)
  }, logical(1))]

  if (length(tmdata_list) == 0) {
    return(list(tmdata_list = list(), markers_list = list()))
  }

  markers <- read_excel(
    "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx",
    sheet = 1
  )
  markers <- markers[markers$specificity > 0.2 & markers$cell_type != "Malignant", ]
  markers_list_initial <- markers %>%
    mutate(cell_type = recode(cell_type, !!!celltype_map)) %>%
    split(.$cell_type)
  markers_ranked <- lapply(markers_list_initial, function(df) {
    combine_marker_scores_post_annotation(df, w_specificity = 0.2, w_sensitivity = 0.8)
  })

  setsN <- lapply(markers_ranked, function(df) {
    df %>% arrange(desc(Combined)) %>% slice_head(n = 100) %>% pull(gene)
  })

  celltypes <- names(setsN)
  all_genes <- unique(unlist(setsN))

  pct_tbl <- bind_rows(lapply(names(tmdata_list), function(sid) {
    obj <- tmdata_list[[sid]]
    meta_calls <- obj@meta.data$celltype_initial
    genes <- intersect(all_genes, rownames(obj))
    if (length(genes) == 0) return(NULL)
    expr <- GetAssayData(obj, slot = "data")[genes, , drop = FALSE]
    split_cells <- split(colnames(obj), meta_calls)
    bind_rows(lapply(names(split_cells), function(ct) {
      cells <- split_cells[[ct]]
      if (length(cells) == 0) return(NULL)
      pct <- Matrix::rowMeans(expr[, cells, drop = FALSE] > 1) * 100
      data.frame(sample = sid, celltype = ct, gene = genes, pct_expr = pct, stringsAsFactors = FALSE)
    }))
  }))

  pct_tbl <- pct_tbl %>% filter(celltype %in% celltypes)
  markers_long <- tibble(celltype = names(setsN), gene = setsN) %>% unnest(gene)

  pct_with_homeflag <- pct_tbl %>%
    left_join(markers_long %>% mutate(is_home = TRUE), by = c("gene", "celltype")) %>%
    mutate(
      is_home = coalesce(is_home, FALSE),
      pass = ifelse(is_home, pct_expr >= 30, pct_expr >= 15)
    )

  gene_ct_long <- pct_with_homeflag %>%
    group_by(gene, celltype, is_home) %>%
    summarise(
      n_total = sum(!is.na(pct_expr)),
      n_pass = sum(pass, na.rm = TRUE),
      pct_samples = ifelse(n_total > 0, 100 * n_pass / n_total, 0),
      .groups = "drop"
    ) %>%
    select(gene, celltype, pct_samples, is_home)

  home_tbl <- markers_long %>%
    left_join(gene_ct_long, by = c("gene", "celltype")) %>%
    mutate(home_pct = coalesce(pct_samples, 0)) %>%
    select(gene, celltype, home_pct)

  off_tbl_raw <- gene_ct_long %>%
    rename(off_celltype = celltype, off_pct = pct_samples) %>%
    inner_join(markers_long, by = "gene")

  off_tbl_filtered <- off_tbl_raw %>%
    rowwise() %>%
    mutate(
      drop_these = list(get_excluded_off_cts_post_annotation(celltype)),
      keep_row = !(off_celltype %in% drop_these)
    ) %>%
    ungroup() %>%
    filter(off_celltype != celltype, keep_row) %>%
    select(-drop_these, -keep_row)

  off_tbl <- off_tbl_filtered %>%
    group_by(gene, celltype) %>%
    summarise(off_pct_max = max(off_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(off_pct_max = ifelse(is.finite(off_pct_max), off_pct_max, 0))

  exclusive_tbl <- home_tbl %>%
    left_join(off_tbl, by = c("gene", "celltype")) %>%
    mutate(
      off_pct_max = coalesce(off_pct_max, 0),
      is_exclusive = (home_pct >= 15) & (off_pct_max <= 30)
    ) %>%
    arrange(desc(is_exclusive), desc(home_pct), off_pct_max)

  save <- exclusive_tbl %>%
    filter(is_exclusive) %>%
    select(gene, home_celltype = celltype, home_pct, off_pct_max)

  list(tmdata_list = tmdata_list, markers_list = split(save$gene, save$home_celltype))
}

annotate_post_annotation_sample <- function(tmdata, markers_list) {
  allowed_pairs <- list(
    c("t.cell", "nk.cell"),
    c("nk.cell", "t.cell"),
    c("b.cell", "plasma"),
    c("plasma", "b.cell")
  )

  Idents(tmdata) <- tmdata$seurat_clusters
  clusters <- levels(Idents(tmdata))
  available_genes <- rownames(tmdata)
  markers_in_data <- lapply(markers_list, function(gene_set) intersect(gene_set, available_genes))
  markers_in_data <- Filter(function(v) length(v) > 0, markers_in_data)
  ct_names <- names(markers_in_data)
  if (length(ct_names) == 0) return(tmdata)

  score_mat <- matrix(
    NA_real_,
    nrow = ncol(tmdata),
    ncol = length(ct_names),
    dimnames = list(colnames(tmdata), ct_names)
  )

  for (cl in clusters) {
    cells_cl <- WhichCells(tmdata, idents = cl)
    if (length(cells_cl) == 0) next
    tm_sub <- subset(tmdata, cells = cells_cl)
    mtx <- GetAssayData(tm_sub, slot = "data")
    cl_features <- lapply(markers_in_data, function(genes) {
      g <- intersect(genes, rownames(mtx))
      if (length(g) == 0) {
        character(0)
      } else {
        m <- Matrix::rowMeans(mtx[g, , drop = FALSE])
        g[order(m, decreasing = TRUE)][seq_len(min(4, length(g)))]
      }
    })
    keep <- vapply(cl_features, function(v) length(v) > 0, logical(1))
    if (!any(keep)) next
    cl_features_kept <- cl_features[keep]
    kept_ct <- names(cl_features_kept)
    tm_sub <- AddModuleScore(
      object = tm_sub,
      features = cl_features_kept,
      name = "mod_tmp_",
      assay = DefaultAssay(tm_sub)
    )
    score_cols <- paste0("mod_tmp_", seq_along(cl_features_kept))
    scdf <- tm_sub@meta.data[, score_cols, drop = FALSE]
    colnames(scdf) <- kept_ct
    score_mat[cells_cl, kept_ct] <- as.matrix(scdf[rownames(scdf), kept_ct, drop = FALSE])
  }

  for (ct in ct_names) {
    tmdata@meta.data[[paste0("mod_", ct)]] <- score_mat[colnames(tmdata), ct]
  }

  mod_cols <- setNames(paste0("mod_", ct_names), ct_names)
  scores_long <- tmdata@meta.data %>%
    mutate(cluster = tmdata$seurat_clusters) %>%
    group_by(cluster) %>%
    summarize(across(all_of(mod_cols), median, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
    mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
    group_by(cluster) %>%
    mutate(z = as.numeric(scale(score))) %>%
    ungroup() %>%
    select(cluster, cell_type, score, z)

  step2_calls <- scores_long %>%
    group_by(cluster) %>%
    summarize(cell_types = list(cell_type), zs = list(z), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      all_z = list(unlist(zs)),
      all_ct = list(unlist(cell_types)),
      max_idx = which.max(all_z),
      top_ct = all_ct[max_idx],
      top_z = all_z[max_idx],
      step2 = {
        az <- all_z
        act <- all_ct
        other_idx <- setdiff(seq_along(az), max_idx)
        if (length(other_idx) == 0) {
          top_ct
        } else {
          margins <- top_z - az[other_idx]
          close_idx <- other_idx[margins < 0.8]
          if (length(close_idx) == 0) {
            top_ct
          } else if (length(close_idx) == 1) {
            other_ct <- act[close_idx]
            pair_vec <- c(top_ct, other_ct)
            is_allowed <- any(vapply(allowed_pairs, function(p) identical(p, pair_vec), logical(1)))
            if (is_allowed) paste0(top_ct, "|", other_ct) else "unresolved"
          } else {
            "unresolved"
          }
        }
      }
    ) %>%
    ungroup() %>%
    select(cluster, step2)

  cl_map <- step2_calls$step2
  names(cl_map) <- step2_calls$cluster
  tmdata$celltype_update <- as.character(cl_map[as.character(tmdata$seurat_clusters)])
  tmdata
}

post_annotation_inputs <- build_post_annotation_inputs(sample_dirs)
post_annotation_cluster_summary_all <- list()
post_annotation_sample_flags <- list()

pdf(file.path(out_dir, "Auto_post_annotation_unresolved_neuron_samples_plots.pdf"), width = 18, height = 13)

for (sample in names(post_annotation_inputs$tmdata_list)) {
  obj <- annotate_post_annotation_sample(
    post_annotation_inputs$tmdata_list[[sample]],
    post_annotation_inputs$markers_list
  )
  anno_col <- get_annotation_column(obj)
  if (is.na(anno_col) || !"seurat_clusters" %in% colnames(obj@meta.data)) {
    rm(obj); next
  }

  unresolved_idx <- is_unresolved_cell(obj@meta.data[[anno_col]])
  unresolved_clusters <- unique(as.character(obj$seurat_clusters[unresolved_idx]))
  genes_in_obj <- intersect(pan_neuronal_markers, rownames(obj))

  if (length(genes_in_obj) == 0 || length(unresolved_clusters) == 0) {
    post_annotation_sample_flags[[sample]] <- data.frame(
      sample = sample,
      has_neuron_like_unresolved_cluster = FALSE,
      n_post_annotation_unresolved_cells = sum(unresolved_idx),
      n_post_annotation_unresolved_clusters = length(unresolved_clusters),
      n_neuron_like_clusters = 0,
      neuron_like_clusters = "",
      stringsAsFactors = FALSE
    )
    rm(obj); next
  }

  expr_all <- GetAssayData(obj, slot = "data")[genes_in_obj, , drop = FALSE]
  sample_cluster_summary <- list()
  for (cl in unresolved_clusters) {
    cl_cells <- colnames(obj)[as.character(obj$seurat_clusters) == cl]
    if (length(cl_cells) < 5) next
    cl_expr <- expr_all[, cl_cells, drop = FALSE]
    pct_expr <- Matrix::rowMeans(cl_expr > 0) * 100
    avg_expr <- Matrix::rowMeans(cl_expr)
    marker_tbl <- data.frame(
      sample = sample,
      cluster = cl,
      marker = names(pct_expr),
      pct_expr = as.numeric(pct_expr),
      avg_expr = as.numeric(avg_expr),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        scan_scope = "post_annotation_unresolved_clusters",
        marker_positive = pct_expr >= min_pct_expr & avg_expr >= min_avg_expr
      )
    marker_tbl$cluster_is_neuron_like <- sum(marker_tbl$marker_positive) >= min_markers_positive
    sample_cluster_summary[[cl]] <- marker_tbl
  }

  if (length(sample_cluster_summary) == 0) {
    rm(obj); next
  }

  sample_cluster_df <- bind_rows(sample_cluster_summary)
  post_annotation_cluster_summary_all[[sample]] <- sample_cluster_df
  neuron_like_clusters <- sample_cluster_df %>%
    group_by(cluster) %>%
    summarize(cluster_is_neuron_like = any(cluster_is_neuron_like), .groups = "drop") %>%
    filter(cluster_is_neuron_like) %>%
    pull(cluster)

  has_neuron <- length(neuron_like_clusters) > 0
  post_annotation_sample_flags[[sample]] <- data.frame(
    sample = sample,
    has_neuron_like_unresolved_cluster = has_neuron,
    n_post_annotation_unresolved_cells = sum(unresolved_idx),
    n_post_annotation_unresolved_clusters = length(unresolved_clusters),
    n_neuron_like_clusters = length(neuron_like_clusters),
    neuron_like_clusters = paste(neuron_like_clusters, collapse = ";"),
    stringsAsFactors = FALSE
  )

  if (has_neuron) {
    message(paste(
      "Flagged post-Annotation unresolved sample:",
      sample,
      "| neuron-like clusters:",
      paste(neuron_like_clusters, collapse = ",")
    ))

    obj$auto_post_annotation_neuron_scan_label <- "other_cluster"
    neuron_like_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% neuron_like_clusters]
    obj$auto_post_annotation_neuron_scan_label[neuron_like_cells] <- "neuron_like_unresolved_cluster"
    obj$auto_post_annotation_neuron_scan_label <- factor(
      obj$auto_post_annotation_neuron_scan_label,
      levels = c("other_cluster", "neuron_like_unresolved_cluster")
    )

    if (!"umap" %in% names(obj@reductions)) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      obj <- ScaleData(obj, verbose = FALSE)
      n_pcs <- min(30, max(5, ncol(obj) - 1))
      obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
      dims_use <- 1:min(20, n_pcs)
      obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
    }

    marker_expr_list <- lapply(pan_neuronal_markers, function(g) {
      if (g %in% rownames(obj)) FetchData(obj, vars = g)[, 1] else rep(0, ncol(obj))
    })
    shared_max_neuron <- max(unlist(marker_expr_list), na.rm = TRUE)
    shared_max_neuron <- max(shared_max_neuron, 0.1)
    limits_neuron <- c(0, ceiling(shared_max_neuron * 10) / 10)

    marker_plots <- lapply(pan_neuronal_markers, function(g) {
      if (g %in% rownames(obj)) {
        obj@meta.data[[paste0("post_annotation_neuronal_scan_", g)]] <- FetchData(obj, vars = g)[, 1]
        score_featureplot_sequential(obj, paste0("post_annotation_neuronal_scan_", g), g, limits_neuron)
      } else {
        ggplot() + theme_void() + ggtitle(paste(g, "(not in sample)"))
      }
    })

    p_scan_umap <- base_dimplot(
      obj,
      "auto_post_annotation_neuron_scan_label",
      "Post-Annotation Unresolved Scan",
      label = FALSE,
      palette = c("other_cluster" = "grey85", "neuron_like_unresolved_cluster" = "#D55E00"),
      show_legend = TRUE
    )

    p_anno_umap <- base_dimplot(
      obj,
      anno_col,
      "Post-Annotation Celltype",
      label = FALSE,
      palette = celltype_colors,
      show_legend = TRUE
    )

    plots_page1 <- marker_plots
    plots_page1[[9]] <- p_scan_umap
    plots_page1[[10]] <- p_anno_umap
    plots_page1[[11]] <- plot_spacer()
    plots_page1[[12]] <- plot_spacer()

    combined_page1 <- wrap_plots(plots_page1, ncol = 4, nrow = 3) +
      plot_annotation(title = paste0(sample, " | Post-Annotation Unresolved Neuronal Markers"))

    obj <- add_score_matrix(obj, build_canonical_score_matrix(obj), prefix = "post_annotation_canonical_score_")
    combined_page2 <- build_expression_page(
      obj = obj,
      score_prefix = "post_annotation_canonical_score_",
      page_title = paste0(sample, " | Post-Annotation Canonical Marker Expression"),
      label_column = anno_col,
      label_title = "Post-Annotation Celltype"
    )

    print(combined_page1)
    print(combined_page2)
  }

  rm(obj)
  gc()
}

dev.off()

post_annotation_cluster_summary_df <- if (length(post_annotation_cluster_summary_all) > 0) {
  bind_rows(post_annotation_cluster_summary_all)
} else {
  data.frame()
}
post_annotation_sample_flags_df <- if (length(post_annotation_sample_flags) > 0) {
  bind_rows(post_annotation_sample_flags)
} else {
  data.frame()
}

if (nrow(post_annotation_cluster_summary_df) > 0) {
  write.csv(
    post_annotation_cluster_summary_df,
    file = file.path(out_dir, "Auto_post_annotation_unresolved_neuron_cluster_summary.csv"),
    row.names = FALSE
  )
}

if (nrow(post_annotation_sample_flags_df) > 0) {
  write.csv(
    post_annotation_sample_flags_df,
    file = file.path(out_dir, "Auto_post_annotation_unresolved_neuron_sample_flags.csv"),
    row.names = FALSE
  )
}

post_annotation_summary_df <- if (nrow(post_annotation_sample_flags_df) > 0) {
  post_annotation_sample_flags_df %>%
    summarize(
      n_samples_scanned = n(),
      n_samples_flagged = sum(has_neuron_like_unresolved_cluster, na.rm = TRUE),
      median_post_annotation_unresolved_cells = median(n_post_annotation_unresolved_cells, na.rm = TRUE)
    )
} else {
  data.frame(
    n_samples_scanned = 0,
    n_samples_flagged = 0,
    median_post_annotation_unresolved_cells = NA_real_
  )
}

write.csv(
  post_annotation_summary_df,
  file = "../updates/new_updates/summaries/Auto_post_annotation_unresolved_neuron_scan_summary.csv",
  row.names = FALSE
)

####################
message("Finished Auto_unresolved_neuron_scan.R")
