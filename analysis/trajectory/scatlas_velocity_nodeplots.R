####################
# Analysis registry:
#   Status: terminal plotting workflow.
#   Script: analysis/trajectory/scatlas_velocity_nodeplots.R
#   Methodology: analysis/methodology/trajectory/scatlas_velocity_methodology.md
#   Inputs:
#     ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_state_nodes.csv
#     ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_state_direction_edges.csv
#   Outputs:
#     ref_outs/Auto_velocity_scATLAS/figures/Auto_scatlas_velocity_nodeplot_by_dataset.pdf
#     updates/new_updates/summaries/Auto_scatlas_velocity_direction_summary.csv
#   Cache/replot behavior: plot-only; reruns from scVelo CSV tables.
#   Run command:
#     Rscript analysis/trajectory/scatlas_velocity_nodeplots.R
#   Conda environment:
#     dmtcp
####################

####################
# Directed node plots for scATLAS RNA velocity state transitions. Nodes are
# five primary state abundances; arrows are positive velocity alignment from
# source-state mean velocity toward target-state centroids.
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(data.table)
})

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

out_dir <- "Auto_velocity_scATLAS"
node_path <- file.path(out_dir, "tables", "Auto_scatlas_velocity_state_nodes.csv")
edge_path <- file.path(out_dir, "tables", "Auto_scatlas_velocity_state_direction_edges.csv")
if (!file.exists(node_path) || !file.exists(edge_path)) {
  stop("Missing scVelo state transition tables.")
}

nodes <- fread(node_path)
edges <- fread(edge_path)
if (!"pct_total_direction_states" %in% names(nodes)) {
  nodes$pct_total_direction_states <- nodes$pct_major
}

state_levels <- c(
  "Classic proliferation",
  "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia",
  "Stress adaptive",
  "Cancer-cell immune mimicry"
)
state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A",
  "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)
layout_df <- data.frame(
  state = state_levels,
  x = c(-1.18, 1.18, 1.18, -1.18, 0),
  y = c(0.76, 0.76, -0.76, -0.76, -1.34),
  stringsAsFactors = FALSE
)

nodes <- nodes %>%
  mutate(
    state = as.character(.data$state),
    dataset = as.character(.data$dataset),
    pct_plot = .data$pct_total_direction_states
  ) %>%
  left_join(layout_df, by = c("state" = "state")) %>%
  mutate(state = factor(.data$state, levels = state_levels))

edges <- edges %>%
  mutate(
    source = as.character(.data$source),
    target = as.character(.data$target),
    dataset = as.character(.data$dataset),
    plotted = .data$velocity_alignment > 0.10
  ) %>%
  left_join(layout_df, by = c("source" = "state")) %>%
  left_join(layout_df, by = c("target" = "state"), suffix = c("", "_to")) %>%
  rename(xend = .data$x_to, yend = .data$y_to) %>%
  mutate(
    source = factor(.data$source, levels = state_levels),
    target = factor(.data$target, levels = state_levels)
  )

combine_nodes <- function(df) {
  df %>%
    group_by(.data$state, .data$x, .data$y) %>%
    summarise(
      cells = sum(.data$cells, na.rm = TRUE),
      total_cells = sum(unique(.data$total_cells), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pct_major = ifelse(.data$total_cells > 0, 100 * .data$cells / .data$total_cells, 0),
      pct_plot = .data$pct_major
    )
}

combine_edges <- function(df) {
  df %>%
    group_by(.data$source, .data$target, .data$x, .data$y, .data$xend, .data$yend) %>%
    summarise(
      velocity_alignment = mean(.data$velocity_alignment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(.data$velocity_alignment > 0.10)
}

draw_nodeplot <- function(ndf, edf, title_text, max_node, max_edge, node_levels, node_cols, labels_map = NULL) {
  ndf <- ndf %>%
    mutate(
      label_raw = as.character(.data$state),
      label_txt = if (!is.null(labels_map)) dplyr::recode(label_raw, !!!labels_map) else label_raw,
      label = paste0(.data$label_txt, "\n", sprintf("%.1f%%", .data$pct_plot))
    )
  edf <- edf %>%
    filter(.data$velocity_alignment > 0.10) %>%
    mutate(
      edge_dx = .data$xend - .data$x,
      edge_dy = .data$yend - .data$y,
      edge_len = sqrt(.data$edge_dx^2 + .data$edge_dy^2),
      edge_ux = ifelse(.data$edge_len > 0, .data$edge_dx / .data$edge_len, 0),
      edge_uy = ifelse(.data$edge_len > 0, .data$edge_dy / .data$edge_len, 0),
      source_idx = match(as.character(.data$source), node_levels),
      target_idx = match(as.character(.data$target), node_levels),
      curve_side = ifelse(.data$source_idx < .data$target_idx, 1, -1),
      perp_x = -.data$edge_uy,
      perp_y = .data$edge_ux,
      x_plot = .data$x + 0.30 * .data$edge_ux,
      y_plot = .data$y + 0.30 * .data$edge_uy,
      xend_plot = .data$xend - 0.36 * .data$edge_ux,
      yend_plot = .data$yend - 0.36 * .data$edge_uy,
      label_x = (.data$x + .data$xend) / 2 + 0.20 * .data$curve_side * .data$perp_x,
      label_y = (.data$y + .data$yend) / 2 + 0.20 * .data$curve_side * .data$perp_y
    )

  ggplot() +
    {
      if (nrow(edf) > 0) {
        geom_curve(
          data = edf,
          aes(x = .data$x_plot, y = .data$y_plot, xend = .data$xend_plot, yend = .data$yend_plot, linewidth = .data$velocity_alignment),
          curvature = 0.22,
          color = "grey10",
          alpha = 0.95,
          arrow = arrow(length = unit(0.24, "inches"), type = "closed"),
          lineend = "round"
        )
      }
    } +
    {
      if (nrow(edf) > 0) {
        geom_label(
          data = edf,
          aes(x = .data$label_x, y = .data$label_y, label = sprintf("%.2f", .data$velocity_alignment)),
          size = 3.2,
          fill = "white",
          linewidth = 0.2,
          fontface = "bold"
        )
      }
    } +
    geom_point(data = ndf, aes(x = .data$x, y = .data$y, size = .data$pct_plot, color = .data$state)) +
    geom_text(data = ndf, aes(x = .data$x, y = .data$y, label = .data$label), size = 3.2, fontface = "bold", lineheight = 0.86, color = "black") +
    scale_color_manual(values = node_cols, drop = FALSE) +
    scale_size(limits = c(0, max_node), range = c(10, 30), guide = "none") +
    scale_linewidth(limits = c(0, max_edge), range = c(0.8, 5.4), guide = "none") +
    coord_equal() +
    expand_limits(x = c(-1.85, 1.85), y = c(-1.70, 1.22)) +
    theme_void(base_size = 13) +
    labs(title = title_text) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
}

safe_max <- function(x, fallback = 1) {
  out <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(out) || out <= 0) fallback else out
}

pdf_path <- file.path(out_dir, "figures", "Auto_scatlas_velocity_nodeplot_by_dataset.pdf")
summary_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/updates/new_updates/summaries"
dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

dataset_levels <- sort(unique(nodes$dataset))
panel_nodes <- list(All = combine_nodes(nodes))
panel_edges <- list(All = combine_edges(edges))
for (dataset_id in dataset_levels) {
  panel_nodes[[dataset_id]] <- combine_nodes(nodes %>% filter(.data$dataset == dataset_id))
  panel_edges[[dataset_id]] <- combine_edges(edges %>% filter(.data$dataset == dataset_id))
}

max_node <- safe_max(unlist(lapply(panel_nodes, function(x) x$pct_plot)))
max_edge <- safe_max(unlist(lapply(panel_edges, function(x) x$velocity_alignment)), 0.2)

state_labels_map <- c(
  "Classic proliferation" = "Classic\nproliferation",
  "Basal to intestinal metaplasia" = "Basal to\nintestinal meta",
  "SMG to intestinal metaplasia" = "SMG to\nintestinal meta",
  "Stress adaptive" = "Stress\nadaptive",
  "Cancer-cell immune mimicry" = "Cancer-cell\nimmune mimicry"
)

message("Writing: ", pdf_path)
pdf(pdf_path, width = 14, height = 8)
print(draw_nodeplot(panel_nodes$All, panel_edges$All, "All raw-BAM scATLAS samples", max_node, max_edge, state_levels, state_cols, state_labels_map))
if (length(dataset_levels) > 0) {
  plots <- lapply(dataset_levels, function(dataset_id) {
    draw_nodeplot(panel_nodes[[dataset_id]], panel_edges[[dataset_id]], dataset_id, max_node, max_edge, state_levels, state_cols, state_labels_map)
  })
  while (length(plots) > 0) {
    print(wrap_plots(plots[seq_len(min(2, length(plots)))], nrow = 1))
    plots <- plots[-seq_len(min(2, length(plots)))]
  }
}
dev.off()

draw_mp_nodeplot_chain <- function(node_path, edge_path, mp_levels, mp_cols, layout_df, pdf_path, title_prefix, labels_map = NULL) {
  if (!file.exists(node_path) || !file.exists(edge_path)) return(NULL)

  nodes <- fread(node_path)
  edges <- fread(edge_path)
  if (!"pct_total_direction_states" %in% names(nodes)) {
    nodes$pct_total_direction_states <- nodes$pct_major
  }

  nodes <- nodes %>%
    mutate(
      state = as.character(.data$state),
      dataset = as.character(.data$dataset),
      pct_plot = .data$pct_total_direction_states
    ) %>%
    left_join(layout_df, by = c("state" = "state")) %>%
    mutate(state = factor(.data$state, levels = mp_levels))

  edges <- edges %>%
    mutate(
      source = as.character(.data$source),
      target = as.character(.data$target),
      dataset = as.character(.data$dataset),
      plotted = .data$velocity_alignment > 0.10
    ) %>%
    left_join(layout_df, by = c("source" = "state")) %>%
    left_join(layout_df, by = c("target" = "state"), suffix = c("", "_to")) %>%
    rename(xend = .data$x_to, yend = .data$y_to) %>%
    mutate(
      source = factor(.data$source, levels = mp_levels),
      target = factor(.data$target, levels = mp_levels)
    )

  dataset_levels <- sort(unique(nodes$dataset))
  panel_nodes <- list(All = combine_nodes(nodes))
  panel_edges <- list(All = combine_edges(edges))
  for (dataset_id in dataset_levels) {
    panel_nodes[[dataset_id]] <- combine_nodes(nodes %>% filter(.data$dataset == dataset_id))
    panel_edges[[dataset_id]] <- combine_edges(edges %>% filter(.data$dataset == dataset_id))
  }

  max_node <- safe_max(unlist(lapply(panel_nodes, function(x) x$pct_plot)))
  max_edge <- safe_max(unlist(lapply(panel_edges, function(x) x$velocity_alignment)), 0.2)

  pdf(pdf_path, width = 14, height = 8)
  print(draw_nodeplot(panel_nodes$All, panel_edges$All, paste(title_prefix, "- All samples"), max_node, max_edge, mp_levels, mp_cols, labels_map))
  if (length(dataset_levels) > 0) {
    plots <- lapply(dataset_levels, function(dataset_id) {
      draw_nodeplot(panel_nodes[[dataset_id]], panel_edges[[dataset_id]], paste(title_prefix, "-", dataset_id), max_node, max_edge, mp_levels, mp_cols, labels_map)
    })
    while (length(plots) > 0) {
      print(wrap_plots(plots[seq_len(min(2, length(plots)))], nrow = 1))
      plots <- plots[-seq_len(min(2, length(plots)))]
    }
  }
  dev.off()
}

basal_mps <- c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+")
smg_mps <- c("MP8+", "MP8b", "MP16", "MP18b", "MP17")

mp_cols <- c(
  "MP14" = "#4DAF4A", "MP3+" = "#8DA0CB", "MP6+" = "#66C2A5",
  "MP11+" = "#FC8D62", "MP9+" = "#A6D854", "MP10+" = "#E78AC3",
  "MP8+" = "#FF7F00", "MP8b" = "#FDBF6F", "MP16" = "#FF9896",
  "MP18b" = "#C5B0D5", "MP17" = "#C49C94"
)

mp_labels_map <- c(
  "MP14" = "MP14\nSquamoid/basal transition",
  "MP3+" = "MP3+\nBasal-columnar invasive epi",
  "MP6+" = "MP6+\nStress-reactive columnar epi",
  "MP11+" = "MP11+\nEpithelial antiviral IFN",
  "MP9+" = "MP9+\nMetabolic columnar epi",
  "MP10+" = "MP10+\nIntestinal metaplasia",
  "MP8+" = "MP8+\nGlandular intestinal meta",
  "MP8b" = "MP8b\nMetabolic intestinal meta",
  "MP16" = "MP16\nMucous-secretory glandular",
  "MP18b" = "MP18b\nMucous-secretory diff",
  "MP17" = "MP17\nImmune-interactive progenitor"
)

basal_layout_df <- data.frame(
  state = basal_mps,
  x = c(-1.18, -0.59, 0.59, 1.18, 0.59, -0.59),
  y = c(0, 0.98, 0.98, 0, -0.98, -0.98),
  stringsAsFactors = FALSE
)

smg_theta <- seq(pi, pi + 2 * pi, length.out = length(smg_mps) + 1L)[seq_along(smg_mps)]
smg_layout_df <- data.frame(
  state = smg_mps,
  x = 1.18 * cos(smg_theta),
  y = 1.18 * sin(smg_theta),
  stringsAsFactors = FALSE
)

draw_mp_nodeplot_chain(
  file.path(out_dir, "tables", "Auto_scatlas_velocity_basal_mp_nodes.csv"),
  file.path(out_dir, "tables", "Auto_scatlas_velocity_basal_mp_direction_edges.csv"),
  basal_mps,
  mp_cols,
  basal_layout_df,
  file.path(out_dir, "figures", "Auto_scatlas_velocity_nodeplot_basal_mps.pdf"),
  "Basal to intestinal metaplasia MPs",
  labels_map = mp_labels_map
)

draw_mp_nodeplot_chain(
  file.path(out_dir, "tables", "Auto_scatlas_velocity_smg_mp_nodes.csv"),
  file.path(out_dir, "tables", "Auto_scatlas_velocity_smg_mp_direction_edges.csv"),
  smg_mps,
  mp_cols,
  smg_layout_df,
  file.path(out_dir, "figures", "Auto_scatlas_velocity_nodeplot_smg_mps.pdf"),
  "SMG to intestinal metaplasia MPs",
  labels_map = mp_labels_map
)

summary_tbl <- edges %>%
  group_by(.data$dataset, .data$source, .data$target) %>%
  summarise(
    sample_n = n_distinct(.data$sample),
    mean_alignment = mean(.data$velocity_alignment, na.rm = TRUE),
    median_alignment = median(.data$velocity_alignment, na.rm = TRUE),
    positive_fraction = mean(.data$velocity_alignment > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(.data$dataset, desc(.data$mean_alignment))
fwrite(summary_tbl, file.path(summary_dir, "Auto_scatlas_velocity_direction_summary.csv"))

message("Node plot PDF done.")
