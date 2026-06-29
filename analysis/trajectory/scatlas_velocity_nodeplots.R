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

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

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
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating"
)
state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8"
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

draw_nodeplot <- function(ndf, edf, title_text, max_node, max_edge) {
  ndf <- ndf %>%
    mutate(
      label = dplyr::recode(
        as.character(.data$state),
        "Classic Proliferative" = "Classic\nProlif.",
        "Basal to Intestinal Metaplasia" = "Basal to\nIntest. Meta",
        "SMG-like Metaplasia" = "SMG-like\nMetaplasia",
        "Stress-adaptive" = "Stress-\nadaptive",
        "Immune Infiltrating" = "Immune\nInfiltrating"
      ),
      label = paste0(.data$label, "\n", sprintf("%.1f%%", .data$pct_plot))
    )
  edf <- edf %>%
    filter(.data$velocity_alignment > 0.10) %>%
    mutate(
      edge_dx = .data$xend - .data$x,
      edge_dy = .data$yend - .data$y,
      edge_len = sqrt(.data$edge_dx^2 + .data$edge_dy^2),
      edge_ux = ifelse(.data$edge_len > 0, .data$edge_dx / .data$edge_len, 0),
      edge_uy = ifelse(.data$edge_len > 0, .data$edge_dy / .data$edge_len, 0),
      source_idx = match(as.character(.data$source), state_levels),
      target_idx = match(as.character(.data$target), state_levels),
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
    scale_color_manual(values = state_cols, drop = FALSE) +
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
summary_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/updates/new_updates/summaries"
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

message("Writing: ", pdf_path)
pdf(pdf_path, width = 14, height = 8)
print(draw_nodeplot(panel_nodes$All, panel_edges$All, "All raw-BAM scATLAS samples", max_node, max_edge))
if (length(dataset_levels) > 0) {
  plots <- lapply(dataset_levels, function(dataset_id) {
    draw_nodeplot(panel_nodes[[dataset_id]], panel_edges[[dataset_id]], dataset_id, max_node, max_edge)
  })
  while (length(plots) > 0) {
    print(wrap_plots(plots[seq_len(min(2, length(plots)))], nrow = 1))
    plots <- plots[-seq_len(min(2, length(plots)))]
  }
}
dev.off()

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
