####################
# Auto_states_hybrid_pairwise_nodeplot.R
# Node plot for real states and pairwise hybrid proportions.
# Multi-class hybrids (>2 groups) are excluded from edge construction.
####################

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")

real_states <- c(
  "Classic_Proliferative",
  "Barretts_Metaplasia",
  "EMT_related",
  "Intestinal_Metaplasia",
  "Immune_Infiltrated"
)

state_groups <- list(
  Classic_Proliferative = c("MP2"),
  Barretts_Metaplasia = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  EMT_related = c("MP13", "MP12"),
  Intestinal_Metaplasia = c("MP18", "MP16"),
  Immune_Infiltrated = c("MP15")
)

group_cols <- c(
  Classic_Proliferative = "#E41A1C",
  Barretts_Metaplasia = "#4DAF4A",
  EMT_related = "#984EA3",
  Intestinal_Metaplasia = "#FF7F00",
  Immune_Infiltrated = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

make_nodeplot_for_mode <- function(mode_name) {
  state_file <- paste0("Auto_topmp_v2_", mode_name, "_states_B.rds")
  mp_adj_file <- paste0("Auto_topmp_v2_", mode_name, "_mp_adj.rds")
  
  if (!file.exists(state_file) || !file.exists(mp_adj_file)) {
    message("Files missing for mode: ", mode_name)
    return(NULL)
  }
  
  state_B <- readRDS(state_file)
  mp_adj <- readRDS(mp_adj_file)

  group_max <- sapply(state_groups, function(mps) {
    mps_avail <- intersect(mps, colnames(mp_adj))
    if (length(mps_avail) == 0) return(rep(NA_real_, nrow(mp_adj)))
    if (length(mps_avail) == 1) return(as.numeric(mp_adj[, mps_avail]))
    apply(mp_adj[, mps_avail, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(mp_adj)

  hybrid_cells <- names(state_B)[state_B == "Hybrid"]
  
  assign_pair <- function(x, names_vec) {
    ord <- names(sort(x, decreasing = TRUE))[1:2]
    ord <- names_vec[names_vec %in% ord]
    paste(ord, collapse = "__")
  }

  pair_labels <- vapply(hybrid_cells, function(cl) assign_pair(group_max[cl, real_states], real_states), character(1))

  pair_df <- data.frame(pair = pair_labels, stringsAsFactors = FALSE) %>%
    count(pair, name = "hybrid_cells") %>%
    separate(pair, into = c("from", "to"), sep = "__", remove = FALSE)

  state_df <- data.frame(state = state_B, stringsAsFactors = FALSE) %>%
    filter(state %in% real_states) %>%
    count(state, name = "cells")

  tot_cells <- length(state_B)
  state_df <- state_df %>% mutate(pct = 100 * cells / tot_cells)
  pair_df <- pair_df %>% mutate(pct = 100 * hybrid_cells / tot_cells)

  n <- length(real_states)
  theta <- seq(0, 2 * pi, length.out = n + 1)[1:n]
  layout_df <- data.frame(
    state = real_states,
    x = cos(theta),
    y = sin(theta),
    stringsAsFactors = FALSE
  )

  node_df <- left_join(layout_df, state_df, by = c("state" = "state"))
  node_df$cells[is.na(node_df$cells)] <- 0
  node_df$pct[is.na(node_df$pct)] <- 0
  # Radial label placement to avoid overlap and obscuring edges
  node_df$label_x <- node_df$x * 1.25
  node_df$label_y <- node_df$y * 1.25

  edge_df <- pair_df %>%
    left_join(layout_df, by = c("from" = "state")) %>%
    rename(x = x, y = y) %>%
    left_join(layout_df, by = c("to" = "state"), suffix = c("", "_to")) %>%
    rename(xend = x_to, yend = y_to)

  p <- ggplot() +
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
      aes(x = label_x, y = label_y, label = paste0(state, "\n", sprintf("%.1f%%", pct))),
      size = 3.5,
      fontface = "bold"
    ) +
    geom_label(
      data = edge_df,
      aes(
        x = (x + xend) / 2,
        y = (y + yend) / 2,
        label = sprintf("%.1f%%", pct)
      ),
      size = 2.6,
      fill = "white",
      label.size = 0,
      fontface = "bold"
    ) +
    scale_color_manual(values = group_cols) +
    scale_size(range = c(8, 22), guide = "none") +
    scale_linewidth(range = c(0.6, 6), guide = "none") +
    coord_equal() +
    expand_limits(x = c(-1.4, 1.4), y = c(-1.4, 1.4)) +
    theme_void(base_size = 14) +
    labs(title = paste0("Pairwise hybrid network - ", mode_name)) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t = 10, b = 10)))
    
  list(plot = p, state_df = state_df, pair_df = pair_df)
}

results <- lapply(requested_modes, make_nodeplot_for_mode)
names(results) <- requested_modes
results <- results[!sapply(results, is.null)]

if (length(results) > 0) {
  pdf("Auto_topmp_v2_hybrid_pairwise_nodeplot_reg_noreg.pdf", width = 6, height = 6)
  for (mode_res in results) {
    print(mode_res$plot)
  }
  dev.off()
}

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- lapply(names(results), function(mode_name) {
  res <- results[[mode_name]]
  bind_rows(
    res$state_df %>% transmute(mode = mode_name, type = "state", label = state, cells = cells, pct = pct),
    res$pair_df %>% transmute(mode = mode_name, type = "pairwise_hybrid", label = pair, cells = hybrid_cells, pct = pct)
  )
})

write.csv(
  bind_rows(summary_rows),
  file.path(summary_dir, "Auto_topmp_v2_hybrid_pairwise_nodeplot_summary.csv"),
  row.names = FALSE
)

message("Saved unified pairwise hybrid node plot and summary.")
