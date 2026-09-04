####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/hybrid_pairwise_distance_nodeplot.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# hybrid_pairwise_distance_nodeplot.R
# Pairwise hybrid node plots using biologically derived state-distance matrices.
#
# Behaviour:
#   - loops through all requested distance methods in one run
#   - writes per-method CSV summaries
#   - writes one combined nodeplot PDF page (3 columns x 2 rows)
#   - writes one combined heatmap PDF page (3 columns x 2 rows)
#
# Layout:
#   - unconstrained 2D MDS layout for better distance fidelity
#   - orientation aligned consistently across methods
#   - outward label placement to avoid overlap with nodes
#
# Input:
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_group_max.rds
#   ref_outs/state_distance_pseudotime/Auto_state_distance_matrices.rds
#
# Output:
#   ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_nodeplot_all_methods.pdf
#   ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_distance_heatmap_all_methods.pdf
#   ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_edges_<method>.csv
#   ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_layout_<method>.csv
#   updates/new_updates/summaries/Auto_task6_hybrid_pairwise_distance_summary.csv
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(patchwork)
library(MASS)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

####################
# Output directories and constants
####################
task_prefix <- "task6"
out_dir <- paste0(task_prefix, "_hybrid_pairwise_distance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/updates/new_updates/summaries"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

real_states <- c(
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
  "Cancer-cell immune mimicry" = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

method_labels <- c(
  "directed_pseudotime_median" = "Directed pseudotime median",
  "directed_pseudotime_mean" = "Directed pseudotime mean",
  "principal_graph_geodesic_medoid" = "Graph geodesic medoid",
  "principal_graph_geodesic_centroid" = "Graph geodesic centroid",
  "umap_centroid_euclidean" = "UMAP centroid Euclidean"
)

####################
# Helper functions
####################
make_safe_method <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

blank_panel <- function() {
  ggplot() + theme_void()
}

scale_target_matrix <- function(distance_mat) {
  mat <- distance_mat
  diag(mat) <- NA_real_
  vals <- mat[upper.tri(mat)]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(mat)
  }
  mat / stats::median(vals)
}

layout_distance_matrix <- function(layout_df, state_order) {
  coords <- as.matrix(layout_df[match(state_order, layout_df$state), c("x", "y")])
  out <- as.matrix(stats::dist(coords))
  rownames(out) <- state_order
  colnames(out) <- state_order
  out
}

make_mds_layout <- function(distance_mat, state_order) {
  mat_use <- scale_target_matrix(distance_mat[state_order, state_order, drop = FALSE])
  mds <- cmdscale(as.dist(mat_use), k = 2, eig = TRUE)

  coords <- as.data.frame(mds$points, stringsAsFactors = FALSE)
  colnames(coords) <- c("x", "y")
  coords$state <- state_order

  coords$x <- coords$x - mean(coords$x)
  coords$y <- coords$y - mean(coords$y)
  coords$radius <- sqrt(coords$x ^ 2 + coords$y ^ 2)
  coords
}

####################
# Layout diagnostics and candidate selection
####################
layout_fit_summary <- function(distance_mat, layout_df, state_order) {
  true_mat <- distance_mat[state_order, state_order, drop = FALSE]
  plotted_mat <- layout_distance_matrix(layout_df, state_order)
  pair_idx <- which(upper.tri(true_mat), arr.ind = TRUE)
  true_vals <- true_mat[upper.tri(true_mat)]
  plotted_vals <- plotted_mat[upper.tri(plotted_mat)]
  scale_factor <- sum(true_vals * plotted_vals) / sum(plotted_vals ^ 2)
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    scale_factor <- 1
  }
  plotted_scaled_mat <- plotted_mat * scale_factor

  pair_df <- data.frame(
    state_a = rownames(true_mat)[pair_idx[, 1]],
    state_b = colnames(true_mat)[pair_idx[, 2]],
    true_distance = true_vals,
    plotted_distance = plotted_vals,
    plotted_distance_scaled = plotted_scaled_mat[upper.tri(plotted_scaled_mat)],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      abs_error = abs(plotted_distance_scaled - true_distance),
      rel_error = ifelse(true_distance > 0, abs_error / true_distance, NA_real_)
    )

  nearest_true <- vapply(seq_len(nrow(true_mat)), function(i) {
    ord <- order(replace(true_mat[i, ], i, Inf))
    colnames(true_mat)[ord[1]]
  }, character(1))

  nearest_plot <- vapply(seq_len(nrow(plotted_mat)), function(i) {
    ord <- order(replace(plotted_mat[i, ], i, Inf))
    colnames(plotted_mat)[ord[1]]
  }, character(1))

  list(
    pair_df = pair_df,
    metrics = data.frame(
      global_scale_factor = scale_factor,
      mean_abs_error = mean(pair_df$abs_error),
      max_abs_error = max(pair_df$abs_error),
      mean_rel_error = mean(pair_df$rel_error, na.rm = TRUE),
      nearest_neighbor_mismatches = sum(nearest_true != nearest_plot),
      stringsAsFactors = FALSE
    ),
    nearest_df = data.frame(
      state = state_order,
      nearest_true = nearest_true,
      nearest_plot = nearest_plot,
      match = nearest_true == nearest_plot,
      stringsAsFactors = FALSE
    )
  )
}

make_best_layout <- function(distance_mat, state_order, reference_layout = NULL) {
  base_layout <- make_mds_layout(distance_mat, state_order)
  candidates <- list(cmdscale = base_layout)

  iso_layout <- tryCatch({
    mat_use <- scale_target_matrix(distance_mat[state_order, state_order, drop = FALSE])
    iso_fit <- MASS::isoMDS(
      as.dist(mat_use),
      y = as.matrix(base_layout[, c("x", "y")]),
      k = 2,
      trace = FALSE
    )

    coords <- as.data.frame(iso_fit$points, stringsAsFactors = FALSE)
    colnames(coords) <- c("x", "y")
    coords$state <- state_order
    coords$x <- coords$x - mean(coords$x)
    coords$y <- coords$y - mean(coords$y)
    coords$radius <- sqrt(coords$x ^ 2 + coords$y ^ 2)
    coords
  }, error = function(e) {
    NULL
  })

  if (!is.null(iso_layout)) {
    candidates$isoMDS <- iso_layout
  }

  candidate_results <- lapply(names(candidates), function(candidate_name) {
    cand <- candidates[[candidate_name]]
    if (!is.null(reference_layout)) {
      cand <- align_layout_to_reference(cand, reference_layout)
    }
    fit <- layout_fit_summary(distance_mat, cand, state_order)
    list(name = candidate_name, layout = cand, fit = fit)
  })

  score_df <- bind_rows(lapply(candidate_results, function(x) {
    x$fit$metrics %>%
      mutate(candidate = x$name)
  })) %>%
    arrange(nearest_neighbor_mismatches, mean_abs_error, max_abs_error)

  best_name <- score_df$candidate[1]
  best_result <- candidate_results[[match(best_name, vapply(candidate_results, `[[`, character(1), "name"))]]

  best_result$layout$layout_method <- best_name
  best_result
}

align_layout_to_reference <- function(layout_df, ref_df) {
  state_order <- ref_df$state
  coords <- as.matrix(layout_df[match(state_order, layout_df$state), c("x", "y")])
  ref_coords <- as.matrix(ref_df[, c("x", "y")])

  transforms <- list(
    diag(c(1, 1)),
    diag(c(-1, 1)),
    diag(c(1, -1)),
    diag(c(-1, -1)),
    matrix(c(0, 1, 1, 0), nrow = 2),
    matrix(c(0, -1, 1, 0), nrow = 2),
    matrix(c(0, 1, -1, 0), nrow = 2),
    matrix(c(0, -1, -1, 0), nrow = 2)
  )

  best_loss <- Inf
  best_coords <- coords

  for (tr in transforms) {
    cand <- coords %*% tr
    loss <- sum((cand - ref_coords) ^ 2)
    if (loss < best_loss) {
      best_loss <- loss
      best_coords <- cand
    }
  }

  layout_df$x <- best_coords[, 1]
  layout_df$y <- best_coords[, 2]
  layout_df$radius <- sqrt(layout_df$x ^ 2 + layout_df$y ^ 2)
  layout_df$align_loss <- best_loss
  layout_df
}

expand_layout <- function(layout_df, factor = 1.18) {
  layout_df$x <- layout_df$x * factor
  layout_df$y <- layout_df$y * factor
  layout_df$radius <- sqrt(layout_df$x ^ 2 + layout_df$y ^ 2)
  layout_df
}

add_label_positions <- function(layout_df) {
  center_x <- mean(layout_df$x)
  center_y <- mean(layout_df$y)

  dx <- layout_df$x - center_x
  dy <- layout_df$y - center_y
  norms <- sqrt(dx ^ 2 + dy ^ 2)
  norms[norms == 0] <- 1

  ux <- dx / norms
  uy <- dy / norms

  span <- max(diff(range(layout_df$x)), diff(range(layout_df$y)))
  if (!is.finite(span) || span <= 0) {
    span <- 1
  }
  label_offset <- 0.10 * span

  layout_df$label_x <- layout_df$x + label_offset * ux
  layout_df$label_y <- layout_df$y + label_offset * uy
  layout_df$hjust <- ifelse(ux >= 0.15, 0, ifelse(ux <= -0.15, 1, 0.5))
  layout_df$vjust <- ifelse(uy >= 0.15, 0, ifelse(uy <= -0.15, 1, 0.5))
  layout_df
}

wrap_state_label <- function(state_name, width = 18) {
  label_map <- c(
    "Classic proliferation" = "Classic\nProliferation",
    "Basal to intestinal metaplasia" = "Basal to\nIntestinal\nMetaplasia",
    "Stress adaptive" = "Stress\nAdaptive",
    "SMG to intestinal metaplasia" = "SMG to\nIntestinal\nMetaplasia",
    "Cancer-cell immune mimicry" = "Cancer-cell\nImmune Mimicry"
  )

  vapply(state_name, function(x) {
    if (x %in% names(label_map)) {
      return(unname(label_map[[x]]))
    }
    paste(strwrap(x, width = width), collapse = "\n")
  }, character(1))
}

build_pair_df <- function(group_max_df, hybrid_cells, real_states) {
  hybrid_pair_calls <- group_max_df %>%
    filter(cell %in% hybrid_cells)

  if (nrow(hybrid_pair_calls) == 0) {
    return(data.frame(
      pair = character(0),
      from = character(0),
      to = character(0),
      hybrid_cells = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  hybrid_score_mat <- as.matrix(hybrid_pair_calls[, real_states, drop = FALSE])
  top_pair_mat <- t(apply(hybrid_score_mat, 1, function(x) {
    ord <- sort(names(sort(x, decreasing = TRUE))[1:2])
    c(from = ord[1], to = ord[2])
  }))

  data.frame(
    cell = hybrid_pair_calls$cell,
    from = top_pair_mat[, "from"],
    to = top_pair_mat[, "to"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(pair = paste(from, to, sep = "__")) %>%
    count(pair, from, to, name = "hybrid_cells")
}

make_nodeplot <- function(edge_df, node_df, method_name, fit_metrics) {
  node_df <- node_df %>%
    mutate(
      label_text = paste0(
        wrap_state_label(state),
        "\n",
        sprintf("%.1f%%", pct)
      )
    )

  ggplot() +
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct),
      color = "grey35",
      alpha = 0.8
    ) +
    geom_point(
      data = node_df,
      aes(x = x, y = y, size = pct, color = state)
    ) +
    geom_text(
      data = node_df,
      aes(
        x = label_x,
        y = label_y,
        label = label_text,
        hjust = hjust,
        vjust = vjust
      ),
      size = 3.3,
      fontface = "bold",
      lineheight = 0.95
    ) +
    geom_label(
      data = edge_df,
      aes(
        x = (x + xend) / 2,
        y = (y + yend) / 2,
        label = paste0(sprintf("%.1f", pct), "%")
      ),
      size = 2.0,
      fill = "white",
      label.size = 0,
      fontface = "bold"
    ) +
    scale_color_manual(values = group_cols) +
    scale_size(range = c(6, 18), guide = "none") +
    scale_linewidth(range = c(0.5, 5.5), guide = "none") +
    coord_equal(clip = "off") +
    expand_limits(
      x = {
        pad <- 0.10 * max(diff(range(node_df$label_x)), diff(range(node_df$x)))
        c(min(node_df$label_x) - pad, max(node_df$label_x) + pad)
      },
      y = {
        pad <- 0.10 * max(diff(range(node_df$label_y)), diff(range(node_df$y)))
        c(min(node_df$label_y) - pad, max(node_df$label_y) + pad)
      }
    ) +
    theme_void(base_size = 12) +
    labs(
      title = method_labels[[method_name]]#,
      # subtitle = paste0(
      #   "2D fit: ", fit_metrics$layout_method,
      #   " | NN mismatches ", fit_metrics$nearest_neighbor_mismatches, "/5",
      #   " | mean abs error ", sprintf("%.2f", fit_metrics$mean_abs_error),
      #   " | scale ", sprintf("%.2f", fit_metrics$global_scale_factor)
      # )
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(t = 8, b = 4)),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, margin = margin(b = 8)),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
}

make_heatmap <- function(distance_mat, method_name, state_order) {
  heatmap_df <- as.data.frame(as.table(distance_mat), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("StateA", "StateB", "Distance")

  heatmap_df$StateA <- factor(heatmap_df$StateA, levels = state_order)
  heatmap_df$StateB <- factor(heatmap_df$StateB, levels = state_order)

  ggplot(heatmap_df, aes(StateB, StateA, fill = Distance)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(Distance), "NA", sprintf("%.2f", Distance))), size = 2.7) +
    scale_fill_gradient(low = "white", high = "firebrick3", na.value = "grey90") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    ) +
    labs(
      title = method_labels[[method_name]],
      x = NULL,
      y = NULL,
      fill = "Distance"
    )
}

####################
# Load inputs
####################
state_B <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
group_max <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_group_max.rds")
distance_list <- readRDS(file.path("state_distance_pseudotime", "Auto_state_distance_matrices.rds"))

requested_methods <- if (length(args) >= 1 && nzchar(args[1])) {
  intersect(names(distance_list), unlist(strsplit(args[1], ",")))
} else {
  names(distance_list)
}

if (length(requested_methods) == 0) {
  stop("No valid distance methods requested. Available methods: ", paste(names(distance_list), collapse = ", "))
}

####################
# Reference orientation from the average scaled distance matrix
####################
average_scaled_distance <- Reduce(
  `+`,
  lapply(requested_methods, function(method_name) {
    scale_target_matrix(distance_list[[method_name]][real_states, real_states, drop = FALSE])
  })
) / length(requested_methods)

reference_layout <- make_best_layout(average_scaled_distance, real_states)$layout %>%
  add_label_positions()

state_order_global <- reference_layout$state

group_max_df <- as.data.frame(group_max, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("cell")

if (!all(real_states %in% colnames(group_max_df))) {
  stop("Expected current centred-state columns were not found in centred_refined_noreg_group_max.rds")
}

hybrid_cells <- names(state_B)[state_B == "Hybrid"]
pair_df_base <- build_pair_df(group_max_df, hybrid_cells, real_states)

tot_cells <- length(state_B)
state_df <- data.frame(
  state = as.character(state_B),
  stringsAsFactors = FALSE
) %>%
  filter(state %in% real_states) %>%
  count(state, name = "cells") %>%
  mutate(pct = 100 * cells / tot_cells)

pair_df_base <- pair_df_base %>%
  mutate(pct = 100 * hybrid_cells / tot_cells)

####################
# Loop through methods
####################
node_plots <- list()
heatmap_plots <- list()
summary_rows <- list()

for (method_name in requested_methods) {
  distance_mat <- distance_list[[method_name]]
  distance_mat <- distance_mat[state_order_global, state_order_global, drop = FALSE]
  diag(distance_mat) <- 0

  layout_result <- make_best_layout(distance_mat, state_order_global, reference_layout = reference_layout)
  layout_df <- layout_result$layout %>%
    expand_layout(factor = 1.18) %>%
    add_label_positions()
  fit_metrics <- layout_result$fit$metrics %>%
    mutate(
      method = method_name,
      layout_method = layout_df$layout_method[1]
    )

  node_df <- layout_df %>%
    left_join(state_df, by = "state")
  node_df$cells[is.na(node_df$cells)] <- 0
  node_df$pct[is.na(node_df$pct)] <- 0

  edge_df <- pair_df_base %>%
    left_join(layout_df, by = c("from" = "state")) %>%
    rename(x = x, y = y) %>%
    left_join(layout_df, by = c("to" = "state"), suffix = c("", "_to")) %>%
    rename(xend = x_to, yend = y_to) %>%
    left_join(
      as.data.frame(as.table(distance_mat), stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        rename(from = Var1, to = Var2, distance = Freq),
      by = c("from", "to")
    )

  safe_method <- make_safe_method(method_name)

  write_csv(
    edge_df,
    file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_edges_", safe_method, ".csv"))
  )

  write_csv(
    node_df,
    file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_layout_", safe_method, ".csv"))
  )

  write_csv(
    layout_result$fit$pair_df,
    file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_fit_pairs_", safe_method, ".csv"))
  )

  write_csv(
    layout_result$fit$nearest_df,
    file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_fit_nearest_", safe_method, ".csv"))
  )

  summary_rows[[method_name]] <- bind_rows(
    fit_metrics %>%
      transmute(
        method = method,
        type = "layout_fit",
        label = layout_method,
        cells = nearest_neighbor_mismatches,
        pct = mean_abs_error,
        distance = max_abs_error
      ),
    node_df %>%
      transmute(
        method = method_name,
        type = "state",
        label = state,
        cells = cells,
        pct = pct,
        distance = NA_real_
      ),
    edge_df %>%
      transmute(
        method = method_name,
        type = "pairwise_hybrid",
        label = pair,
        cells = hybrid_cells,
        pct = pct,
        distance = distance
      )
  )

  node_plots[[method_name]] <- make_nodeplot(edge_df, node_df, method_name, fit_metrics)
  heatmap_plots[[method_name]] <- make_heatmap(distance_mat, method_name, state_order_global)
}

write_csv(
  bind_rows(summary_rows),
  file.path(summary_dir, "Auto_task6_hybrid_pairwise_distance_summary.csv")
)

####################
# Combined panels: 3 columns x 2 rows
####################
while (length(node_plots) < 6) {
  node_plots[[paste0("blank_node_", length(node_plots) + 1)]] <- blank_panel()
}

while (length(heatmap_plots) < 6) {
  heatmap_plots[[paste0("blank_heat_", length(heatmap_plots) + 1)]] <- blank_panel()
}

combined_nodeplot <- wrap_plots(node_plots[1:6], ncol = 3, nrow = 2) +
  plot_annotation(
    title = "Hybrid Pairwise Nodeplots Across Distance Methods",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

combined_heatmap <- wrap_plots(heatmap_plots[1:6], ncol = 3, nrow = 2) +
  plot_annotation(
    title = "State Distance Heatmaps Across Methods",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

ggsave(
  file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_nodeplot_all_methods.pdf")),
  combined_nodeplot,
  width = 20,
  height = 13,
  limitsize = FALSE
)

ggsave(
  file.path(out_dir, paste0("Auto_", task_prefix, "_hybrid_pairwise_distance_heatmap_all_methods.pdf")),
  combined_heatmap,
  width = 20,
  height = 13,
  limitsize = FALSE
)

message("Saved combined hybrid nodeplot and heatmap panels for methods: ", paste(requested_methods, collapse = ", "))
