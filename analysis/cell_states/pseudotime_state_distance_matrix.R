####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/pseudotime_state_distance_matrix.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# pseudotime_state_distance_matrix.R
# Monocle3-based state-to-state distance comparison for the 5 primary
# current centred refined noreg states across all valid samples.
#
# Methods:
#   1. Directed pseudotime gaps (median)
#   2. Directed pseudotime gaps (mean)
#   3. Principal-graph geodesic distance between state medoids
#   4. Principal-graph geodesic distance between state centroids
#   5. UMAP centroid Euclidean distance (root-free baseline)
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#
# Output:
#   ref_outs/state_distance_pseudotime/Auto_state_distance_sample_summary.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_directed_pseudotime.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_geodesic_medoid.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_geodesic_centroid.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_umap_centroid.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_summary.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_matrices.rds
#   ref_outs/state_distance_pseudotime/Auto_*_matrix.csv
#   ref_outs/state_distance_pseudotime/Auto_state_distance_method_comparison_heatmap.pdf
#   ref_outs/state_distance_pseudotime/Auto_state_pseudotime_all_valid_samples.pdf
#   ref_outs/state_distance_pseudotime/sample_state_trajectories/<sample>_pseudotime_states.rds
#   updates/new_updates/summaries/Auto_state_distance_summary.csv
####################

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(igraph)
library(patchwork)

live_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
setwd(live_dir)

####################
# Output directories and constants
####################
out_dir <- "state_distance_pseudotime"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "sample_state_trajectories"), recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(live_dir, "../updates/new_updates/summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

primary_states <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "Stress adaptive",
  "SMG to intestinal metaplasia",
  "Cancer-cell immune mimicry"
)

group_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "Stress adaptive" = "#984EA3",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Cancer-cell immune mimicry" = "#377EB8"
)

####################
# Sample inclusion thresholds
# - 20 root cells: enough root support for stable root assignment
# - 20 cells in at least one second state: avoid one-state samples
# - 80 total primary-state cells: enough cells to learn a useful graph
####################
ROOT_MIN_CELLS <- 20
OTHER_MIN_CELLS <- 20
TOTAL_MIN_CELLS <- 80
MIN_STATES_OVER_THRESHOLD <- 2

####################
# Helper functions
####################
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

make_empty_matrix <- function(states) {
  mat <- matrix(NA_real_, nrow = length(states), ncol = length(states))
  rownames(mat) <- states
  colnames(mat) <- states
  diag(mat) <- 0
  mat
}

matrix_to_long <- function(mat, method_name) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(state_a = Var1, state_b = Var2, distance = Freq) %>%
    mutate(method = method_name)
}

make_symmetric_matrix <- function(summary_df, method_name, states) {
  out <- make_empty_matrix(states)

  df_use <- summary_df %>%
    filter(method == method_name) %>%
    mutate(
      state_min = pmin(state_a, state_b),
      state_max = pmax(state_a, state_b)
    ) %>%
    group_by(state_min, state_max) %>%
    summarise(distance = mean(mean_distance, na.rm = TRUE), .groups = "drop")

  if (nrow(df_use) == 0) return(out)

  for (i in seq_len(nrow(df_use))) {
    a <- df_use$state_min[i]
    b <- df_use$state_max[i]
    out[a, b] <- df_use$distance[i]
    out[b, a] <- df_use$distance[i]
  }

  diag(out) <- 0
  out
}

build_heatmap_pdf <- function(long_df, file_path) {
  if (nrow(long_df) == 0) return(invisible(NULL))

  plot_df <- long_df %>%
    mutate(
      state_a = factor(state_a, levels = primary_states),
      state_b = factor(state_b, levels = primary_states)
    )

  p <- ggplot(plot_df, aes(x = state_a, y = state_b, fill = distance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(distance), "NA", sprintf("%.2f", distance))), size = 2.7) +
    scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "grey90") +
    facet_wrap(~method) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave(file_path, p, width = 13, height = 8)
}

build_state_pseudotime_page <- function(cds, sample_id, sample_meta) {
  state_counts <- table(factor(sample_meta$state, levels = primary_states))
  legend_labels <- setNames(
    paste0(primary_states, " (", as.integer(state_counts[primary_states]), ")"),
    primary_states
  )

  p_states <- plot_cells(
    cds,
    color_cells_by = "state_B",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.8
  ) +
    scale_color_manual(
      values = group_cols[primary_states],
      breaks = primary_states,
      labels = legend_labels[primary_states],
      name = "State",
      na.value = "grey80",
      drop = FALSE,
      guide = guide_legend(override.aes = list(size = 4))
    ) +
    labs(
      title = paste0("Cell states - ", sample_id, " (n = ", nrow(sample_meta), ")"),
      color = NULL
    ) +
    theme_minimal(base_size = 11)

  p_pseudotime <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.8
  ) +
    scale_color_viridis_c() +
    labs(
      title = paste0("Pseudotime - ", sample_id, " | Basal root"),
      color = "Pseudotime"
    ) +
    theme_minimal(base_size = 11)

  p_states + p_pseudotime + plot_layout(guides = "collect")
}

get_state_medoid_cell <- function(umap_df) {
  if (nrow(umap_df) == 1) return(umap_df$cell[1])
  dist_mat <- as.matrix(dist(umap_df[, c("UMAP_1", "UMAP_2"), drop = FALSE]))
  umap_df$cell[which.min(rowSums(dist_mat))]
}

compute_graph_weights <- function(graph_obj, graph_coords) {
  edge_df <- igraph::as_data_frame(graph_obj, what = "edges")
  if (nrow(edge_df) == 0) return(graph_obj)

  edge_weights <- purrr::map2_dbl(edge_df$from, edge_df$to, function(from_node, to_node) {
    from_xy <- graph_coords[, from_node]
    to_xy <- graph_coords[, to_node]
    sqrt(sum((from_xy - to_xy) ^ 2))
  })

  igraph::E(graph_obj)$weight <- edge_weights
  graph_obj
}

nearest_graph_vertex <- function(point_xy, graph_coords) {
  vertex_names <- colnames(graph_coords)
  dists <- apply(graph_coords, 2, function(node_xy) {
    sqrt(sum((point_xy - node_xy) ^ 2))
  })
  vertex_names[which.min(dists)]
}

coerce_graph_vertex_name <- function(raw_vertex, graph_obj) {
  graph_nodes <- igraph::V(graph_obj)$name
  raw_vertex <- as.character(raw_vertex)[1]

  if (length(raw_vertex) == 0 || is.na(raw_vertex) || raw_vertex == "") {
    return(NA_character_)
  }

  if (raw_vertex %in% graph_nodes) {
    return(raw_vertex)
  }

  raw_num <- suppressWarnings(as.numeric(raw_vertex))
  if (!is.na(raw_num)) {
    raw_idx <- as.integer(raw_num)
    if (raw_idx >= 1 && raw_idx <= length(graph_nodes)) {
      return(graph_nodes[raw_idx])
    }
  }

  NA_character_
}

get_root_pr_node <- function(cds, root_cells) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(closest_vertex)) return(NULL)

  if (is.matrix(closest_vertex)) {
    closest_vertex <- closest_vertex[, 1, drop = TRUE]
  }

  root_cells <- intersect(root_cells, names(closest_vertex))
  if (length(root_cells) == 0) return(NULL)

  root_vertex_raw <- names(sort(table(as.character(closest_vertex[root_cells])), decreasing = TRUE))[1]
  graph_obj <- principal_graph(cds)[["UMAP"]]
  root_vertex <- coerce_graph_vertex_name(root_vertex_raw, graph_obj)

  if (!is.na(root_vertex)) {
    return(root_vertex)
  }

  NULL
}

prepare_sample_for_monocle <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  n_pcs <- min(30, max(2, ncol(seurat_obj) - 1))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, verbose = FALSE)

  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  cds
}

extract_graph_bits <- function(cds) {
  graph_obj <- principal_graph(cds)[["UMAP"]]
  graph_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

  if (is.null(graph_obj) || is.null(graph_coords) || is.null(closest_vertex)) {
    return(NULL)
  }

  if (is.matrix(closest_vertex)) {
    closest_vertex <- closest_vertex[, 1, drop = TRUE]
  }

  graph_obj <- compute_graph_weights(graph_obj, graph_coords)

  list(
    graph = graph_obj,
    graph_coords = graph_coords,
    closest_vertex = closest_vertex
  )
}

####################
# Load data and align cells
####################
message("Loading data ...")
tmdata_all <- readRDS(file.path(live_dir, "EAC_Ref_epi.rds"))
state_vec <- readRDS(file.path(live_dir, "Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds"))

common_cells <- intersect(Cells(tmdata_all), names(state_vec))
tmdata_all <- tmdata_all[, common_cells]
state_vec <- state_vec[common_cells]
tmdata_all$state_B <- state_vec[Cells(tmdata_all)]

meta_df <- data.frame(
  cell = Cells(tmdata_all),
  orig.ident = as.character(tmdata_all$orig.ident),
  state = as.character(tmdata_all$state_B),
  stringsAsFactors = FALSE
) %>%
  filter(state %in% primary_states)

####################
# Sample inclusion summary
####################
sample_state_counts <- meta_df %>%
  count(orig.ident, state, name = "n_cells") %>%
  complete(orig.ident, state = primary_states, fill = list(n_cells = 0))

sample_summary <- sample_state_counts %>%
  group_by(orig.ident) %>%
  summarise(
    total_primary_cells = sum(n_cells),
    n_states_over_threshold = sum(n_cells >= OTHER_MIN_CELLS),
    qualifies_total = total_primary_cells >= TOTAL_MIN_CELLS,
    qualifies_multistate = n_states_over_threshold >= MIN_STATES_OVER_THRESHOLD,
    valid_sample = qualifies_total & qualifies_multistate,
    .groups = "drop"
  )

root_summary <- sample_state_counts %>%
  left_join(sample_summary %>% select(orig.ident, valid_sample), by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(other_state_max = vapply(state, function(this_state) {
    max(n_cells[state != this_state], na.rm = TRUE)
  }, numeric(1))) %>%
  ungroup() %>%
  mutate(
    valid_root = valid_sample &
      n_cells >= ROOT_MIN_CELLS &
      other_state_max >= OTHER_MIN_CELLS
  ) %>%
  rename(root_state = state, root_state_cells = n_cells)

write_csv(root_summary, file.path(out_dir, "Auto_state_distance_sample_summary.csv"))

valid_samples <- sample_summary %>%
  filter(valid_sample) %>%
  pull(orig.ident)

message("Valid samples detected: ", length(valid_samples))

####################
# Per-sample processing
####################
directed_rows <- list()
geodesic_medoid_rows <- list()
geodesic_centroid_rows <- list()
umap_centroid_rows <- list()

pdf(
  file.path(out_dir, "Auto_state_pseudotime_all_valid_samples.pdf"),
  width = 14,
  height = 6,
  onefile = TRUE
)

for (sample_id in valid_samples) {
  message("Processing sample: ", sample_id)

  sample_cells <- meta_df %>%
    filter(orig.ident == sample_id) %>%
    pull(cell)

  if (length(sample_cells) < TOTAL_MIN_CELLS) next

  sample_obj <- tmdata_all[, sample_cells]
  sample_meta <- meta_df %>%
    filter(cell %in% sample_cells)

  cds <- tryCatch(
    prepare_sample_for_monocle(sample_obj),
    error = function(e) {
      message("Skipping sample after monocle3 prep failure: ", sample_id, " | ", e$message)
      NULL
    }
  )

  if (is.null(cds)) next

  graph_bits <- tryCatch(
    extract_graph_bits(cds),
    error = function(e) {
      message("Skipping sample after graph extraction failure: ", sample_id, " | ", e$message)
      NULL
    }
  )

  if (is.null(graph_bits)) next

  umap_mat <- reducedDims(cds)$UMAP
  umap_df <- data.frame(
    cell = rownames(umap_mat),
    UMAP_1 = umap_mat[, 1],
    UMAP_2 = umap_mat[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_meta, by = "cell")

  rep_df <- sample_meta %>%
    count(state, name = "n_cells") %>%
    mutate(
      medoid_cell = vapply(state, function(state_name) {
        get_state_medoid_cell(umap_df %>% filter(state == state_name))
      }, character(1))
    ) %>%
    left_join(umap_df %>% select(cell, UMAP_1, UMAP_2), by = c("medoid_cell" = "cell")) %>%
    rowwise() %>%
    mutate(
      centroid_x = mean(umap_df$UMAP_1[umap_df$state == state], na.rm = TRUE),
      centroid_y = mean(umap_df$UMAP_2[umap_df$state == state], na.rm = TRUE),
      medoid_vertex = coerce_graph_vertex_name(graph_bits$closest_vertex[medoid_cell], graph_bits$graph),
      centroid_vertex = nearest_graph_vertex(c(centroid_x, centroid_y), graph_bits$graph_coords)
    ) %>%
    ungroup()

  ####################
  # Root-free distances for this sample
  ####################
  present_states <- intersect(primary_states, rep_df$state)

  for (state_a in present_states) {
    for (state_b in present_states) {
      row_a <- rep_df %>% filter(state == state_a)
      row_b <- rep_df %>% filter(state == state_b)
      if (nrow(row_a) == 0 || nrow(row_b) == 0) next

      medoid_distance <- suppressWarnings(
        igraph::distances(
          graph_bits$graph,
          v = row_a$medoid_vertex[1],
          to = row_b$medoid_vertex[1],
          weights = igraph::E(graph_bits$graph)$weight
        )[1, 1]
      )
      medoid_distance <- if (is.finite(medoid_distance)) as.numeric(medoid_distance) else NA_real_

      centroid_distance <- suppressWarnings(
        igraph::distances(
          graph_bits$graph,
          v = row_a$centroid_vertex[1],
          to = row_b$centroid_vertex[1],
          weights = igraph::E(graph_bits$graph)$weight
        )[1, 1]
      )
      centroid_distance <- if (is.finite(centroid_distance)) as.numeric(centroid_distance) else NA_real_

      euclidean_distance <- sqrt(
        (row_a$centroid_x[1] - row_b$centroid_x[1]) ^ 2 +
          (row_a$centroid_y[1] - row_b$centroid_y[1]) ^ 2
      )

      geodesic_medoid_rows[[length(geodesic_medoid_rows) + 1]] <- data.frame(
        sample = sample_id,
        state_a = state_a,
        state_b = state_b,
        distance = medoid_distance,
        stringsAsFactors = FALSE
      )

      geodesic_centroid_rows[[length(geodesic_centroid_rows) + 1]] <- data.frame(
        sample = sample_id,
        state_a = state_a,
        state_b = state_b,
        distance = centroid_distance,
        stringsAsFactors = FALSE
      )

      umap_centroid_rows[[length(umap_centroid_rows) + 1]] <- data.frame(
        sample = sample_id,
        state_a = state_a,
        state_b = state_b,
        distance = euclidean_distance,
        stringsAsFactors = FALSE
      )
    }
  }

  ####################
  # Directed pseudotime distances for each valid root
  ####################
  sample_root_rows <- root_summary %>%
    filter(orig.ident == sample_id, valid_root)

  if (nrow(sample_root_rows) == 0) next

  basal_root_cells <- sample_meta$cell[sample_meta$state == "Basal to intestinal metaplasia"]
  if (length(basal_root_cells) >= ROOT_MIN_CELLS) {
    basal_ordered_cds <- tryCatch(
      order_cells(cds, root_cells = basal_root_cells),
      error = function(e) {
        message("Basal-root order_cells failed for sample plot: ", sample_id, " | ", e$message)
        NULL
      }
    )

    if (!is.null(basal_ordered_cds)) {
      basal_pseudotime <- pseudotime(basal_ordered_cds)
      basal_pseudotime[is.infinite(basal_pseudotime)] <- NA_real_

      print(build_state_pseudotime_page(
        cds = basal_ordered_cds,
        sample_id = sample_id,
        sample_meta = sample_meta
      ))

      saveRDS(
        basal_pseudotime,
        file.path(out_dir, "sample_state_trajectories", paste0(sample_id, "_pseudotime_states.rds"))
      )
    }
  }

  for (i in seq_len(nrow(sample_root_rows))) {
    root_state <- sample_root_rows$root_state[i]
    root_cells <- sample_meta$cell[sample_meta$state == root_state]
    root_node <- get_root_pr_node(cds, root_cells)

    if (is.null(root_node)) {
      message("Skipping root for sample due to missing root node: ", sample_id, " | ", root_state)
      next
    }

    ordered_cds <- tryCatch(
      order_cells(cds, root_pr_nodes = root_node),
      error = function(e) {
        message("order_cells failed for sample/root: ", sample_id, " | ", root_state, " | ", e$message)
        NULL
      }
    )

    if (is.null(ordered_cds)) next

    pt_df <- data.frame(
      cell = names(pseudotime(ordered_cds)),
      pseudotime = as.numeric(pseudotime(ordered_cds)),
      stringsAsFactors = FALSE
    ) %>%
      left_join(sample_meta, by = "cell")

    state_pt <- pt_df %>%
      group_by(state) %>%
      summarise(
        median_pseudotime = safe_median(pseudotime),
        mean_pseudotime = safe_mean(pseudotime),
        n_cells_with_pt = sum(is.finite(pseudotime)),
        .groups = "drop"
      )

    root_median <- state_pt$median_pseudotime[state_pt$state == root_state]
    root_mean <- state_pt$mean_pseudotime[state_pt$state == root_state]
    if (length(root_median) == 0 || length(root_mean) == 0) next

    for (target_state in primary_states) {
      target_row <- state_pt %>% filter(state == target_state)
      if (nrow(target_row) == 0) next

      directed_rows[[length(directed_rows) + 1]] <- data.frame(
        sample = sample_id,
        root_state = root_state,
        target_state = target_state,
        root_state_cells = sample_root_rows$root_state_cells[i],
        target_state_cells = target_row$n_cells_with_pt[1],
        median_distance = target_row$median_pseudotime[1] - root_median[1],
        mean_distance = target_row$mean_pseudotime[1] - root_mean[1],
        stringsAsFactors = FALSE
      )
    }
  }
}

dev.off()

directed_df <- bind_rows(directed_rows)
geodesic_medoid_df <- bind_rows(geodesic_medoid_rows)
geodesic_centroid_df <- bind_rows(geodesic_centroid_rows)
umap_centroid_df <- bind_rows(umap_centroid_rows)

write_csv(directed_df, file.path(out_dir, "Auto_state_distance_directed_pseudotime.csv"))
write_csv(geodesic_medoid_df, file.path(out_dir, "Auto_state_distance_geodesic_medoid.csv"))
write_csv(geodesic_centroid_df, file.path(out_dir, "Auto_state_distance_geodesic_centroid.csv"))
write_csv(umap_centroid_df, file.path(out_dir, "Auto_state_distance_umap_centroid.csv"))

####################
# Summaries and final matrices
####################
directed_pairwise_median <- directed_df %>%
  filter(root_state != target_state) %>%
  mutate(
    state_a = pmin(root_state, target_state),
    state_b = pmax(root_state, target_state),
    distance = abs(median_distance)
  ) %>%
  group_by(sample, state_a, state_b) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop")

directed_pairwise_mean <- directed_df %>%
  filter(root_state != target_state) %>%
  mutate(
    state_a = pmin(root_state, target_state),
    state_b = pmax(root_state, target_state),
    distance = abs(mean_distance)
  ) %>%
  group_by(sample, state_a, state_b) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop")

geodesic_medoid_sym <- geodesic_medoid_df %>%
  mutate(
    state_min = pmin(state_a, state_b),
    state_max = pmax(state_a, state_b)
  ) %>%
  group_by(sample, state_min, state_max) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  rename(state_a = state_min, state_b = state_max)

geodesic_centroid_sym <- geodesic_centroid_df %>%
  mutate(
    state_min = pmin(state_a, state_b),
    state_max = pmax(state_a, state_b)
  ) %>%
  group_by(sample, state_min, state_max) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  rename(state_a = state_min, state_b = state_max)

umap_centroid_sym <- umap_centroid_df %>%
  mutate(
    state_min = pmin(state_a, state_b),
    state_max = pmax(state_a, state_b)
  ) %>%
  group_by(sample, state_min, state_max) %>%
  summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  rename(state_a = state_min, state_b = state_max)

summary_df <- bind_rows(
  directed_pairwise_median %>%
    group_by(state_a, state_b) %>%
    summarise(
      method = "directed_pseudotime_median",
      mean_distance = safe_mean(distance),
      median_distance = safe_median(distance),
      n_samples_used = sum(is.finite(distance)),
      .groups = "drop"
    ),
  directed_pairwise_mean %>%
    group_by(state_a, state_b) %>%
    summarise(
      method = "directed_pseudotime_mean",
      mean_distance = safe_mean(distance),
      median_distance = safe_median(distance),
      n_samples_used = sum(is.finite(distance)),
      .groups = "drop"
    ),
  geodesic_medoid_sym %>%
    group_by(state_a, state_b) %>%
    summarise(
      method = "principal_graph_geodesic_medoid",
      mean_distance = safe_mean(distance),
      median_distance = safe_median(distance),
      n_samples_used = sum(is.finite(distance)),
      .groups = "drop"
    ),
  geodesic_centroid_sym %>%
    group_by(state_a, state_b) %>%
    summarise(
      method = "principal_graph_geodesic_centroid",
      mean_distance = safe_mean(distance),
      median_distance = safe_median(distance),
      n_samples_used = sum(is.finite(distance)),
      .groups = "drop"
    ),
  umap_centroid_sym %>%
    group_by(state_a, state_b) %>%
    summarise(
      method = "umap_centroid_euclidean",
      mean_distance = safe_mean(distance),
      median_distance = safe_median(distance),
      n_samples_used = sum(is.finite(distance)),
      .groups = "drop"
    )
) %>%
  arrange(method, state_a, state_b)

write_csv(summary_df, file.path(out_dir, "Auto_state_distance_summary.csv"))
write_csv(summary_df, file.path(summary_dir, "Auto_state_distance_summary.csv"))

matrix_list <- list(
  directed_pseudotime_median = make_symmetric_matrix(summary_df, "directed_pseudotime_median", primary_states),
  directed_pseudotime_mean = make_symmetric_matrix(summary_df, "directed_pseudotime_mean", primary_states),
  principal_graph_geodesic_medoid = make_symmetric_matrix(summary_df, "principal_graph_geodesic_medoid", primary_states),
  principal_graph_geodesic_centroid = make_symmetric_matrix(summary_df, "principal_graph_geodesic_centroid", primary_states),
  umap_centroid_euclidean = make_symmetric_matrix(summary_df, "umap_centroid_euclidean", primary_states)
)

saveRDS(matrix_list, file.path(out_dir, "Auto_state_distance_matrices.rds"))

for (method_name in names(matrix_list)) {
  write.csv(
    matrix_list[[method_name]],
    file.path(out_dir, paste0("Auto_", method_name, "_matrix.csv")),
    quote = FALSE
  )
}

heatmap_df <- bind_rows(
  matrix_to_long(matrix_list$directed_pseudotime_median, "directed_pseudotime_median"),
  matrix_to_long(matrix_list$directed_pseudotime_mean, "directed_pseudotime_mean"),
  matrix_to_long(matrix_list$principal_graph_geodesic_medoid, "principal_graph_geodesic_medoid"),
  matrix_to_long(matrix_list$principal_graph_geodesic_centroid, "principal_graph_geodesic_centroid"),
  matrix_to_long(matrix_list$umap_centroid_euclidean, "umap_centroid_euclidean")
)

build_heatmap_pdf(
  heatmap_df,
  file.path(out_dir, "Auto_state_distance_method_comparison_heatmap.pdf")
)

message("Saved state distance outputs to: ", file.path(getwd(), out_dir))
