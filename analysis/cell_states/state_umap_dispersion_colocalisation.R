####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/state_umap_dispersion_colocalisation.R
#   Methodology: analysis/methodology/cell_states/state_umap_dispersion_colocalisation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# state_umap_dispersion_colocalisation.R
#
# Quantifies how epithelial noreg Approach B states occupy per-sample UMAP
# space and how locally mixed or self-colocalised their cells are.
#
# Metrics:
#   1. Dispersion: distance from each cell to the centroid of its assigned
#      state or top-MP label in the relevant per-sample UMAP. The primary
#      comparable value is sample-normalised by the median all-cell distance
#      from the sample centroid.
#   2. Co-localisation: fraction of the k nearest UMAP neighbours assigned to
#      the same state or top-MP label.
#
# State-level scope:
#   Uses the saved cell sets from:
#   ref_outs/state_distance_pseudotime/sample_state_trajectories/
#   and rebuilds per-sample UMAPs on the primary state cells only.
#
# Basal-to-intestinal-metaplasia MP scope:
#   Rebuilds a basal-only per-sample UMAP and labels each basal cell by the
#   top adjusted current MP among MP14, MP3+, MP6+, MP11+, MP9+, and MP10+.
#
# Inputs:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds
#   ref_outs/state_distance_pseudotime/sample_state_trajectories/*.rds
#
# Regeneration:
#   If sample_state_trajectories are missing, this script first runs
#   analysis/cell_states/pseudotime_state_distance_matrix.R in the active R
#   environment, then rechecks for the required trajectory RDS files.
#
# Output:
#   ref_outs/state_umap_dispersion_colocalisation/intermediate/
#     Auto_state_umap_coordinates.rds
#     Auto_basal_topmp_umap_coordinates.rds
#     Auto_state_umap_cell_metrics.rds
#     Auto_basal_topmp_umap_cell_metrics.rds
#   ref_outs/state_umap_dispersion_colocalisation/tables/
#     Auto_state_umap_cell_metrics.csv.gz
#     Auto_state_umap_label_metrics_per_sample.csv
#     Auto_state_umap_label_metrics_overall.csv
#     Auto_basal_topmp_umap_cell_metrics.csv.gz
#     Auto_basal_topmp_label_metrics_per_sample.csv
#     Auto_basal_topmp_label_metrics_overall.csv
#     Auto_state_top_mp_composition.csv
#   ref_outs/state_umap_dispersion_colocalisation/figures/
#     Auto_state_umap_colocalisation_boxplot.pdf
#     Auto_state_umap_dispersion_boxplot.pdf
#     Auto_state_umap_dispersion_vs_colocalisation.pdf
#     Auto_state_umap_colocalisation_umap_pages.pdf
#     Auto_basal_topmp_colocalisation_boxplot.pdf
#     Auto_basal_topmp_dispersion_boxplot.pdf
#     Auto_basal_topmp_dispersion_vs_colocalisation.pdf
#     Auto_basal_topmp_colocalisation_umap_pages.pdf
#   ref_outs/state_umap_dispersion_colocalisation/logs/
#     Auto_state_umap_dispersion_colocalisation_run_summary.rds
#     Auto_state_umap_dispersion_colocalisation_run_summary.txt
#   updates/new_updates/summaries/
#     Auto_state_umap_dispersion_colocalisation_summary.csv
#
# Cache/replot:
#   SCREF_FORCE_REBUILD=TRUE ignores cached UMAP coordinates and metrics.
#   SCREF_REPLOT_ONLY=TRUE reuses cached metrics and regenerates plots/tables.
#   SCREF_COLOCAL_K controls nearest-neighbour k; default 15.
#
# Run:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/cell_states/state_umap_dispersion_colocalisation.R
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(FNN)
  library(patchwork)
  library(scales)
})

set.seed(12345)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_helpers.R")

setwd(SCREF_PROJECT_DIR)

####################
# Configuration
####################
out_root <- file.path(SCREF_REF_OUTS_DIR, "state_umap_dispersion_colocalisation")
out_paths <- ensure_output_tiers(out_root)
dir.create(SCREF_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

env_true <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = ifelse(default, "TRUE", "FALSE"))
  tolower(value) %in% c("true", "1", "yes", "y")
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) == 0 || is.na(value)) default else value
}

COLOCAL_K <- env_int("SCREF_COLOCAL_K", 15)
STATE_MIN_CELLS <- env_int("SCREF_STATE_UMAP_MIN_CELLS", 50)
BASAL_MIN_CELLS <- env_int("SCREF_BASAL_UMAP_MIN_CELLS", 50)
FORCE_REBUILD <- env_true("SCREF_FORCE_REBUILD", FALSE)
REPLOT_ONLY <- env_true("SCREF_REPLOT_ONLY", FALSE)

primary_states <- SCREF_PRIMARY_STATE_ORDER
basal_state <- "Basal to intestinal metaplasia"
basal_mps <- SCREF_STATE_GROUPS[[basal_state]]
basal_mp_levels <- label_mps(basal_mps)

state_colours <- state_colours_current(primary_states)
state_colours <- state_colours[primary_states]

basal_mp_colours <- setNames(
  c("#D95F02", "#7570B3", "#E7298A", "#66A61E", "#1B9E77", "#6A3D9A"),
  basal_mp_levels
)

trajectory_dir <- file.path(SCREF_REF_OUTS_DIR, "state_distance_pseudotime", "sample_state_trajectories")

input_files <- c(
  SCREF_EPI_RDS,
  SCREF_STATE_NOREG_RDS,
  SCREF_STATE_NOREG_MP_ADJ_RDS,
  trajectory_dir
)

output_files <- c(
  file.path(out_paths[["intermediate"]], "Auto_state_umap_coordinates.rds"),
  file.path(out_paths[["intermediate"]], "Auto_basal_topmp_umap_coordinates.rds"),
  file.path(out_paths[["intermediate"]], "Auto_state_umap_cell_metrics.rds"),
  file.path(out_paths[["intermediate"]], "Auto_basal_topmp_umap_cell_metrics.rds"),
  file.path(out_paths[["tables"]], "Auto_state_umap_label_metrics_overall.csv"),
  file.path(out_paths[["tables"]], "Auto_basal_topmp_label_metrics_overall.csv"),
  file.path(out_paths[["figures"]], "Auto_state_umap_colocalisation_boxplot.pdf"),
  file.path(out_paths[["figures"]], "Auto_basal_topmp_colocalisation_boxplot.pdf"),
  file.path(SCREF_SUMMARY_DIR, "Auto_state_umap_dispersion_colocalisation_summary.csv")
)

run_summary <- start_run_summary(
  script = "analysis/cell_states/state_umap_dispersion_colocalisation.R",
  input_files = input_files,
  output_files = output_files,
  parameters = list(
    colocal_k = COLOCAL_K,
    state_min_cells = STATE_MIN_CELLS,
    basal_min_cells = BASAL_MIN_CELLS,
    force_rebuild = FORCE_REBUILD,
    replot_only = REPLOT_ONLY,
    state_source = "centred_refined_noreg_states.rds"
  )
)

####################
# Utility functions
####################
cache_or_build <- function(path, build_fun, label) {
  if (file.exists(path) && !FORCE_REBUILD) {
    message("Reusing cache: ", path)
    run_summary <<- add_cache_status(run_summary, label, path, reused = TRUE)
    return(readRDS(path))
  }

  if (REPLOT_ONLY) {
    stop("SCREF_REPLOT_ONLY=TRUE but cache is missing: ", path)
  }

  message("Building cache: ", path)
  value <- build_fun()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path)
  run_summary <<- add_cache_status(run_summary, label, path, reused = FALSE)
  value
}

sample_from_trajectory_file <- function(path) {
  sub("_pseudotime_states\\.rds$", "", basename(path))
}

get_trajectory_files <- function() {
  files <- list.files(trajectory_dir, pattern = "\\.rds$", full.names = TRUE)
  files <- files[grepl("_pseudotime_states\\.rds$", basename(files))]
  sort(files)
}

ensure_trajectory_files <- function() {
  files <- get_trajectory_files()
  if (length(files) > 0) {
    return(files)
  }

  if (REPLOT_ONLY) {
    stop("No trajectory RDS files found and SCREF_REPLOT_ONLY=TRUE: ", trajectory_dir)
  }

  message("Missing sample_state_trajectories; regenerating pseudotime_state_distance_matrix.R ...")
  status <- system2(
    "Rscript",
    args = file.path(SCREF_ANALYSIS_DIR, "cell_states", "pseudotime_state_distance_matrix.R")
  )
  if (!identical(status, 0L)) {
    stop("Regeneration failed for pseudotime_state_distance_matrix.R with status: ", status)
  }

  files <- get_trajectory_files()
  if (length(files) == 0) {
    stop("Regeneration completed but no trajectory RDS files were found in: ", trajectory_dir)
  }

  run_summary$parameters$trajectory_regenerated <<- TRUE
  files
}

read_trajectory_cells <- function(files) {
  message("Reading trajectory cell sets: ", length(files), " files")
  map_dfr(files, function(path) {
    x <- readRDS(path)
    tibble(
      sample = sample_from_trajectory_file(path),
      cell = names(x),
      basal_root_pseudotime = as.numeric(x)
    )
  }) |>
    filter(!is.na(cell), nzchar(cell)) |>
    distinct(sample, cell, .keep_all = TRUE)
}

add_top_state_mp <- function(meta_df, mp_adj, state_groups) {
  top_mp <- rep(NA_character_, nrow(meta_df))

  for (state_name in names(state_groups)) {
    state_cells <- meta_df$cell[meta_df$state == state_name]
    state_mps <- intersect(state_groups[[state_name]], colnames(mp_adj))

    if (length(state_cells) == 0 || length(state_mps) == 0) {
      next
    }

    state_scores <- mp_adj[state_cells, state_mps, drop = FALSE]
    top_idx <- max.col(state_scores, ties.method = "first")
    top_mp[match(state_cells, meta_df$cell)] <- colnames(state_scores)[top_idx]
  }

  meta_df |>
    mutate(
      top_state_mp = top_mp,
      top_state_mp_label = label_mps(top_state_mp)
    )
}

prepare_umap_for_cells <- function(seurat_obj, cells, sample_id, scope_name) {
  if (length(cells) < 3) {
    return(tibble())
  }

  sample_obj <- seurat_obj[, cells]
  if ("RNA" %in% Assays(sample_obj)) {
    DefaultAssay(sample_obj) <- "RNA"
  }

  sample_obj <- NormalizeData(sample_obj, verbose = FALSE)
  sample_obj <- FindVariableFeatures(sample_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sample_obj <- ScaleData(sample_obj, verbose = FALSE)

  n_pcs <- min(30, ncol(sample_obj) - 1)
  if (n_pcs < 2) {
    return(tibble())
  }

  sample_obj <- RunPCA(
    sample_obj,
    features = VariableFeatures(sample_obj),
    npcs = n_pcs,
    verbose = FALSE
  )
  dims_use <- min(15, n_pcs)
  sample_obj <- RunUMAP(sample_obj, dims = 1:dims_use, seed.use = 12345, verbose = FALSE)

  coords <- Embeddings(sample_obj, "umap")
  out <- tibble(
    sample = sample_id,
    cell = rownames(coords),
    UMAP_1 = coords[, 1],
    UMAP_2 = coords[, 2],
    umap_scope = scope_name,
    n_umap_cells = nrow(coords)
  )

  rm(sample_obj)
  gc(verbose = FALSE)
  out
}

build_umap_coordinates <- function(seurat_obj, meta_df, scope_name, min_cells) {
  sample_counts <- meta_df |>
    count(sample, name = "n_cells") |>
    filter(n_cells >= min_cells) |>
    arrange(desc(n_cells), sample)

  message("Building ", scope_name, " UMAPs for ", nrow(sample_counts), " samples")

  rows <- vector("list", nrow(sample_counts))
  for (i in seq_len(nrow(sample_counts))) {
    sample_id <- sample_counts$sample[i]
    sample_cells <- meta_df$cell[meta_df$sample == sample_id]
    message("  [", i, "/", nrow(sample_counts), "] ", sample_id, " | cells=", length(sample_cells))

    rows[[i]] <- tryCatch(
      prepare_umap_for_cells(seurat_obj, sample_cells, sample_id, scope_name),
      error = function(e) {
        message("    Skipping UMAP after error: ", e$message)
        tibble()
      }
    )
  }

  bind_rows(rows)
}

polygon_area <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) {
    return(NA_real_)
  }

  hull_idx <- chull(x, y)
  hx <- x[hull_idx]
  hy <- y[hull_idx]
  hx_next <- c(hx[-1], hx[1])
  hy_next <- c(hy[-1], hy[1])
  abs(sum(hx * hy_next - hy * hx_next)) / 2
}

safe_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok])
}

compute_sample_cell_metrics <- function(df, k) {
  if (nrow(df) == 0) {
    return(tibble())
  }

  coords <- as.matrix(df[, c("UMAP_1", "UMAP_2")])
  labels <- as.character(df$label)
  n_cells <- nrow(df)

  sample_centroid <- colMeans(coords)
  sample_centroid_dist <- sqrt(rowSums((coords - matrix(sample_centroid, n_cells, 2, byrow = TRUE)) ^ 2))
  sample_scale <- median(sample_centroid_dist[sample_centroid_dist > 0], na.rm = TRUE)
  if (!is.finite(sample_scale) || sample_scale <= 0) {
    sample_scale <- max(diff(range(coords[, 1], na.rm = TRUE)), diff(range(coords[, 2], na.rm = TRUE)), na.rm = TRUE)
  }
  if (!is.finite(sample_scale) || sample_scale <= 0) {
    sample_scale <- 1
  }

  sample_hull_area <- polygon_area(coords[, 1], coords[, 2])

  k_use <- min(k, n_cells - 1)
  if (k_use >= 1) {
    nn <- FNN::get.knn(coords, k = k_use)
    same_mat <- matrix(FALSE, nrow = n_cells, ncol = k_use)
    for (j in seq_len(k_use)) {
      same_mat[, j] <- labels[nn$nn.index[, j]] == labels
    }
    same_neighbor_n <- rowSums(same_mat, na.rm = TRUE)
    same_neighbor_score <- same_neighbor_n / k_use
  } else {
    same_neighbor_n <- rep(NA_real_, n_cells)
    same_neighbor_score <- rep(NA_real_, n_cells)
  }

  centroids <- df |>
    group_by(label) |>
    summarise(
      label_centroid_x = mean(UMAP_1, na.rm = TRUE),
      label_centroid_y = mean(UMAP_2, na.rm = TRUE),
      label_cells_in_sample = n(),
      .groups = "drop"
    )

  df |>
    left_join(centroids, by = "label") |>
    mutate(
      neighbor_k = k_use,
      same_neighbor_n = same_neighbor_n,
      same_neighbor_score = same_neighbor_score,
      distance_to_label_centroid = sqrt((UMAP_1 - label_centroid_x) ^ 2 + (UMAP_2 - label_centroid_y) ^ 2),
      dispersion_norm = distance_to_label_centroid / sample_scale,
      sample_centroid_x = sample_centroid[1],
      sample_centroid_y = sample_centroid[2],
      sample_scale = sample_scale,
      sample_hull_area = sample_hull_area
    )
}

compute_label_metrics <- function(coord_df, meta_df, label_col, label_type, label_levels, k) {
  metric_df <- coord_df |>
    inner_join(meta_df, by = c("sample", "cell")) |>
    mutate(
      label = as.character(.data[[label_col]]),
      label_type = label_type
    ) |>
    filter(!is.na(label), nzchar(label))

  sample_rows <- split(metric_df, metric_df$sample)
  cell_metrics <- map_dfr(sample_rows, compute_sample_cell_metrics, k = k) |>
    mutate(
      label = factor(label, levels = label_levels),
      label = as.character(label)
    ) |>
    filter(!is.na(label))

  per_sample <- cell_metrics |>
    group_by(label_type, sample, label) |>
    summarise(
      n_cells = n(),
      neighbor_k = max(neighbor_k, na.rm = TRUE),
      mean_same_neighbor_score = mean(same_neighbor_score, na.rm = TRUE),
      median_same_neighbor_score = median(same_neighbor_score, na.rm = TRUE),
      mean_dispersion_norm = mean(dispersion_norm, na.rm = TRUE),
      median_dispersion_norm = median(dispersion_norm, na.rm = TRUE),
      p90_dispersion_norm = quantile(dispersion_norm, 0.9, na.rm = TRUE, names = FALSE),
      mean_distance_to_centroid = mean(distance_to_label_centroid, na.rm = TRUE),
      median_distance_to_centroid = median(distance_to_label_centroid, na.rm = TRUE),
      label_hull_area = polygon_area(UMAP_1, UMAP_2),
      sample_hull_area = first(sample_hull_area),
      label_hull_area_norm = label_hull_area / sample_hull_area,
      .groups = "drop"
    ) |>
    mutate(
      label_hull_area_norm = ifelse(is.finite(label_hull_area_norm), label_hull_area_norm, NA_real_)
    )

  overall <- per_sample |>
    group_by(label_type, label) |>
    summarise(
      n_samples = n_distinct(sample),
      n_cells_total = sum(n_cells),
      mean_same_neighbor_score = safe_weighted_mean(mean_same_neighbor_score, n_cells),
      median_sample_same_neighbor_score = median(median_same_neighbor_score, na.rm = TRUE),
      mean_dispersion_norm = safe_weighted_mean(mean_dispersion_norm, n_cells),
      median_sample_dispersion_norm = median(median_dispersion_norm, na.rm = TRUE),
      mean_p90_dispersion_norm = safe_weighted_mean(p90_dispersion_norm, n_cells),
      mean_distance_to_centroid = safe_weighted_mean(mean_distance_to_centroid, n_cells),
      median_sample_distance_to_centroid = median(median_distance_to_centroid, na.rm = TRUE),
      mean_label_hull_area_norm = safe_weighted_mean(label_hull_area_norm, n_cells),
      .groups = "drop"
    ) |>
    arrange(match(label, label_levels))

  list(cell_metrics = cell_metrics, per_sample = per_sample, overall = overall)
}

plot_colocalisation_boxplot <- function(cell_metrics, per_sample, label_levels, colours, output_path, x_title = NULL) {
  cell_metrics <- cell_metrics |>
    mutate(label = factor(label, levels = label_levels))
  per_sample <- per_sample |>
    mutate(label = factor(label, levels = label_levels))

  p <- ggplot(cell_metrics, aes(x = label, y = same_neighbor_score, fill = label, color = label)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.82, linewidth = 0.55, color = "black") +
    geom_point(
      data = per_sample,
      aes(y = mean_same_neighbor_score),
      position = position_jitter(width = 0.15, height = 0),
      size = 1.7,
      alpha = 0.75,
      stroke = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
    scale_color_manual(values = colours, guide = "none", drop = FALSE) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(x = x_title, y = "Same-label neighbours") +
    ppt_theme(14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25)
    )

  ggsave(output_path, p, width = 7.4, height = 5.8)
}

plot_dispersion_boxplot <- function(cell_metrics, per_sample, label_levels, colours, output_path, x_title = NULL) {
  cell_metrics <- cell_metrics |>
    mutate(label = factor(label, levels = label_levels))
  per_sample <- per_sample |>
    mutate(label = factor(label, levels = label_levels))

  p <- ggplot(cell_metrics, aes(x = label, y = dispersion_norm, fill = label, color = label)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.82, linewidth = 0.55, color = "black") +
    geom_point(
      data = per_sample,
      aes(y = mean_dispersion_norm),
      position = position_jitter(width = 0.15, height = 0),
      size = 1.7,
      alpha = 0.75,
      stroke = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
    scale_color_manual(values = colours, guide = "none", drop = FALSE) +
    labs(x = x_title, y = "Dispersion from label centroid\n(sample-normalised UMAP distance)") +
    ppt_theme(14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25)
    )

  ggsave(output_path, p, width = 7.4, height = 5.8)
}

plot_dispersion_colocalisation_scatter <- function(overall, label_levels, colours, output_path) {
  plot_df <- overall |>
    mutate(label = factor(label, levels = label_levels))

  p <- ggplot(
    plot_df,
    aes(
      x = mean_dispersion_norm,
      y = mean_same_neighbor_score,
      fill = label,
      size = n_cells_total
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.9, stroke = 0.45) +
    geom_text(aes(label = label), size = 3.4, nudge_y = 0.025, check_overlap = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_size_continuous(range = c(4, 10), name = "Cells") +
    labs(
      x = "Mean dispersion from label centroid\n(sample-normalised UMAP distance)",
      y = "Mean same-label neighbours"
    ) +
    ppt_theme(14) +
    ####################
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold")
    )
    ####################

  ggsave(output_path, p, width = 7.2, height = 5.8)
}

plot_umap_pages <- function(cell_metrics, label_levels, colours, output_path, title_prefix) {
  samples <- cell_metrics |>
    count(sample, name = "n_cells") |>
    arrange(desc(n_cells), sample) |>
    pull(sample)

  if (length(samples) == 0) {
    return(invisible(NULL))
  }

  grDevices::cairo_pdf(output_path, width = 12.5, height = 5.8, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (sample_id in samples) {
    df <- cell_metrics |>
      filter(sample == sample_id) |>
      mutate(label = factor(label, levels = label_levels)) |>
      arrange(label)

    ####################
    point_size <- if (nrow(df) > 4000) {
      0.35
    } else if (nrow(df) > 1500) {
      0.55
    } else if (nrow(df) > 500) {
      1.0
    } else if (nrow(df) > 150) {
      1.8
    } else {
      3.0
    }
    ####################

    p_label <- ggplot(df, aes(UMAP_1, UMAP_2, color = label)) +
      geom_point(size = point_size, alpha = 0.82, stroke = 0) +
      scale_color_manual(values = colours, drop = FALSE, name = NULL) +
      coord_equal() +
      theme_void(base_size = 12) +
      ####################
      theme(
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_blank()
      ) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      ####################
      labs(title = "Label")

    p_coloc <- ggplot(df, aes(UMAP_1, UMAP_2, color = same_neighbor_score)) +
      geom_point(size = point_size, alpha = 0.9, stroke = 0) +
      scale_color_viridis_c(
        option = "C",
        limits = c(0, 1),
        labels = scales::percent_format(accuracy = 1),
        name = "Same-label\nneighbours"
      ) +
      coord_equal() +
      theme_void(base_size = 12) +
      ####################
      theme(
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold")
      ) +
      guides(color = guide_colorbar(barwidth = 1.2, barheight = 8)) +
      ####################
      labs(title = "Co-localisation")

    print(
      (p_label + p_coloc) +
        plot_annotation(
          title = paste0(title_prefix, ": ", sample_id),
          subtitle = paste0("Cells: ", nrow(df), " | k = ", max(df$neighbor_k, na.rm = TRUE)),
          theme = theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 11)
          )
        )
    )
  }

  invisible(output_path)
}

####################
# Load and align data
####################
message("Checking upstream trajectory files ...")
trajectory_files <- ensure_trajectory_files()
trajectory_cells <- read_trajectory_cells(trajectory_files)

message("Loading epithelial object and state/MP assignments ...")
tmdata_all <- readRDS(SCREF_EPI_RDS)
state_vec <- readRDS(SCREF_STATE_NOREG_RDS)
mp_adj <- readRDS(SCREF_STATE_NOREG_MP_ADJ_RDS)

common_cells <- Reduce(
  intersect,
  list(Cells(tmdata_all), names(state_vec), rownames(mp_adj), trajectory_cells$cell)
)
if (length(common_cells) == 0) {
  stop("No overlapping cells across epithelial object, state vector, MP matrix, and trajectory outputs.")
}

trajectory_cells <- trajectory_cells |>
  filter(cell %in% common_cells) |>
  distinct(cell, .keep_all = TRUE)

common_cells <- trajectory_cells$cell
tmdata_all <- tmdata_all[, common_cells]
state_vec <- state_vec[common_cells]
mp_adj <- mp_adj[common_cells, , drop = FALSE]

meta_base <- tibble(
  cell = Cells(tmdata_all),
  sample = as.character(tmdata_all$orig.ident),
  study = as.character(tmdata_all$study),
  state = as.character(state_vec[Cells(tmdata_all)])
) |>
  inner_join(trajectory_cells |> select(cell, basal_root_pseudotime), by = "cell") |>
  filter(state %in% primary_states)

meta_primary <- add_top_state_mp(meta_base, mp_adj, SCREF_STATE_GROUPS)

state_top_mp_composition <- meta_primary |>
  count(state, top_state_mp, top_state_mp_label, name = "n_cells") |>
  group_by(state) |>
  mutate(pct_cells = 100 * n_cells / sum(n_cells)) |>
  ungroup() |>
  arrange(match(state, primary_states), desc(n_cells))

write_csv(
  state_top_mp_composition,
  file.path(out_paths[["tables"]], "Auto_state_top_mp_composition.csv")
)

basal_meta <- meta_primary |>
  filter(state == basal_state, top_state_mp %in% basal_mps) |>
  mutate(top_state_mp_label = factor(top_state_mp_label, levels = basal_mp_levels)) |>
  filter(!is.na(top_state_mp_label)) |>
  mutate(top_state_mp_label = as.character(top_state_mp_label))

if (nrow(basal_meta) == 0) {
  stop("No basal-to-intestinal-metaplasia cells with basal MP labels were available.")
}

####################
# Build or reuse UMAP coordinate caches
####################
state_coord_path <- file.path(out_paths[["intermediate"]], "Auto_state_umap_coordinates.rds")
basal_coord_path <- file.path(out_paths[["intermediate"]], "Auto_basal_topmp_umap_coordinates.rds")

state_coords <- cache_or_build(
  state_coord_path,
  function() build_umap_coordinates(tmdata_all, meta_primary, "primary_states", STATE_MIN_CELLS),
  "state_umap_coordinates"
)

basal_coords <- cache_or_build(
  basal_coord_path,
  function() build_umap_coordinates(tmdata_all, basal_meta, "basal_topmp", BASAL_MIN_CELLS),
  "basal_topmp_umap_coordinates"
)

if (nrow(state_coords) == 0) {
  stop("No state-level UMAP coordinates were generated.")
}
if (nrow(basal_coords) == 0) {
  stop("No basal-only UMAP coordinates were generated.")
}

####################
# Build or reuse metric caches
####################
state_metric_path <- file.path(out_paths[["intermediate"]], "Auto_state_umap_cell_metrics.rds")
basal_metric_path <- file.path(out_paths[["intermediate"]], "Auto_basal_topmp_umap_cell_metrics.rds")

state_metrics <- cache_or_build(
  state_metric_path,
  function() compute_label_metrics(
    coord_df = state_coords,
    meta_df = meta_primary,
    label_col = "state",
    label_type = "state",
    label_levels = primary_states,
    k = COLOCAL_K
  ),
  "state_metrics"
)

basal_metrics <- cache_or_build(
  basal_metric_path,
  function() compute_label_metrics(
    coord_df = basal_coords,
    meta_df = basal_meta,
    label_col = "top_state_mp_label",
    label_type = "basal_top_mp",
    label_levels = basal_mp_levels,
    k = COLOCAL_K
  ),
  "basal_topmp_metrics"
)

####################
# Write tables
####################
message("Writing tables ...")
write_csv(
  state_metrics$cell_metrics,
  file.path(out_paths[["tables"]], "Auto_state_umap_cell_metrics.csv.gz")
)
write_csv(
  state_metrics$per_sample,
  file.path(out_paths[["tables"]], "Auto_state_umap_label_metrics_per_sample.csv")
)
write_csv(
  state_metrics$overall,
  file.path(out_paths[["tables"]], "Auto_state_umap_label_metrics_overall.csv")
)

write_csv(
  basal_metrics$cell_metrics,
  file.path(out_paths[["tables"]], "Auto_basal_topmp_umap_cell_metrics.csv.gz")
)
write_csv(
  basal_metrics$per_sample,
  file.path(out_paths[["tables"]], "Auto_basal_topmp_label_metrics_per_sample.csv")
)
write_csv(
  basal_metrics$overall,
  file.path(out_paths[["tables"]], "Auto_basal_topmp_label_metrics_overall.csv")
)

compact_summary <- bind_rows(
  state_metrics$overall |>
    mutate(analysis_scope = "state", display_label = label),
  basal_metrics$overall |>
    mutate(analysis_scope = "basal_top_mp", display_label = label)
) |>
  transmute(
    analysis_scope,
    display_label,
    n_samples,
    n_cells_total,
    mean_same_neighbor_score,
    median_sample_same_neighbor_score,
    mean_dispersion_norm,
    median_sample_dispersion_norm,
    mean_label_hull_area_norm
  )

write_csv(
  compact_summary,
  file.path(SCREF_SUMMARY_DIR, "Auto_state_umap_dispersion_colocalisation_summary.csv")
)

####################
# Write figures
####################
message("Writing figures ...")
plot_colocalisation_boxplot(
  state_metrics$cell_metrics,
  state_metrics$per_sample,
  primary_states,
  state_colours,
  file.path(out_paths[["figures"]], "Auto_state_umap_colocalisation_boxplot.pdf")
)

plot_dispersion_boxplot(
  state_metrics$cell_metrics,
  state_metrics$per_sample,
  primary_states,
  state_colours,
  file.path(out_paths[["figures"]], "Auto_state_umap_dispersion_boxplot.pdf")
)

plot_dispersion_colocalisation_scatter(
  state_metrics$overall,
  primary_states,
  state_colours,
  file.path(out_paths[["figures"]], "Auto_state_umap_dispersion_vs_colocalisation.pdf")
)

plot_umap_pages(
  state_metrics$cell_metrics,
  primary_states,
  state_colours,
  file.path(out_paths[["figures"]], "Auto_state_umap_colocalisation_umap_pages.pdf"),
  "Primary-state UMAP"
)

plot_colocalisation_boxplot(
  basal_metrics$cell_metrics,
  basal_metrics$per_sample,
  basal_mp_levels,
  basal_mp_colours,
  file.path(out_paths[["figures"]], "Auto_basal_topmp_colocalisation_boxplot.pdf")
)

plot_dispersion_boxplot(
  basal_metrics$cell_metrics,
  basal_metrics$per_sample,
  basal_mp_levels,
  basal_mp_colours,
  file.path(out_paths[["figures"]], "Auto_basal_topmp_dispersion_boxplot.pdf")
)

plot_dispersion_colocalisation_scatter(
  basal_metrics$overall,
  basal_mp_levels,
  basal_mp_colours,
  file.path(out_paths[["figures"]], "Auto_basal_topmp_dispersion_vs_colocalisation.pdf")
)

plot_umap_pages(
  basal_metrics$cell_metrics,
  basal_mp_levels,
  basal_mp_colours,
  file.path(out_paths[["figures"]], "Auto_basal_topmp_colocalisation_umap_pages.pdf"),
  "Basal-to-intestinal-metaplasia top-MP UMAP"
)

####################
# Finish
####################
write_run_summary(
  run_summary,
  file.path(out_paths[["logs"]], "Auto_state_umap_dispersion_colocalisation_run_summary.rds")
)

message("State UMAP dispersion/co-localisation complete: ", out_root)
