####################
# Auto_pseudotime_trajectory_reports_improved.R
#
# Unified trajectory analysis and visualization for Part A samples.
# Rebuilds Monocle3 trajectories, extracts exact principal graphs,
# and generates multi-panel per-sample reports.
#
# Improvements:
#   1. Removed the graph-skeleton-only panel
#   2. Removed the 1D linear progression panel
#   3. Root is explicitly labeled on principal-graph panels
#   4. Ridge heights are weighted by relative state abundance
#   5. Layout adjusted to a clean 2 x 2 report
#
# Panels included:
#   1. Principal Graph + projected cells colored by state + ROOT label
#   2. Principal Graph + projected cells colored by pseudotime + ROOT label
#   3. UMAP colored by pseudotime (cell coordinates)
#   4. Count-weighted pseudotime density ridges
#
# Output:
#   ref_outs/Auto_pseudotime_trajectory_summary_reports_partA.pdf
#   ref_outs/pseudotime_trajectory_assets/<sample>_cds.rds
#   ref_outs/pseudotime_trajectory_assets/<sample>_projections.csv
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(igraph)
  library(readr)
  library(ggridges)
  library(scales)
})

set.seed(12345)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# Configuration
####################
group_cols <- c(
  "Classic Proliferative"          = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"                = "#984EA3",
  "SMG-like Metaplasia"            = "#FF7F00",
  "Immune Infiltrating"            = "#377EB8"
)

state_order <- c(
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "Classic Proliferative"
)

target_states <- names(group_cols)
root_state <- "Basal to Intestinal Metaplasia"

partA_dir <- "task1_pseudotime_updated_naming/partA"
asset_dir <- "pseudotime_trajectory_assets"
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- "Auto_pseudotime_trajectory_summary_reports_partA.pdf"

####################
# Helper Functions
####################

# Project a point (px, py) onto a line segment (ax, ay) -> (bx, by)
project_point_to_segment <- function(px, py, ax, ay, bx, by) {
  abx <- bx - ax
  aby <- by - ay
  ab2 <- abx * abx + aby * aby

  if (!is.finite(ab2) || ab2 == 0) {
    return(list(x = ax, y = ay, t = 0, dist2 = (px - ax)^2 + (py - ay)^2))
  }

  t <- ((px - ax) * abx + (py - ay) * aby) / ab2
  t <- max(0, min(1, t))

  proj_x <- ax + t * abx
  proj_y <- ay + t * aby
  dist2 <- (px - proj_x)^2 + (py - proj_y)^2

  list(x = proj_x, y = proj_y, t = t, dist2 = dist2)
}

# Prepare trajectory for an individual sample
prepare_sample_trajectory <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  n_pcs <- min(30, ncol(seurat_obj) - 1)
  if (n_pcs < 2) stop("Too few cells/features for PCA.")

  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, verbose = FALSE)

  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE)

  cds
}

# Extract nodes and edges from the CDS graph
extract_graph_structure <- function(cds) {
  graph_obj <- principal_graph(cds)[["UMAP"]]
  graph_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst

  edge_df <- igraph::as_data_frame(graph_obj, what = "edges") %>%
    mutate(
      x    = graph_coords[1, from],
      y    = graph_coords[2, from],
      xend = graph_coords[1, to],
      yend = graph_coords[2, to]
    )

  node_df <- data.frame(
    node = colnames(graph_coords),
    x = graph_coords[1, ],
    y = graph_coords[2, ],
    stringsAsFactors = FALSE
  )

  list(edges = edge_df, nodes = node_df)
}

# Identify the principal graph root node used by Monocle3 ordering
get_root_label_position <- function(cds, graph_bits, root_cells) {
  aux <- cds@principal_graph_aux[["UMAP"]]

  # Try to identify the principal node most associated with the root cells
  root_node <- NA_character_

  if (!is.null(aux$pr_graph_cell_proj_closest_vertex)) {
    closest_vertex <- aux$pr_graph_cell_proj_closest_vertex

    # Handle both matrix/data.frame structures safely
    if (is.matrix(closest_vertex) || is.data.frame(closest_vertex)) {
      closest_vertex <- as.data.frame(closest_vertex, stringsAsFactors = FALSE)

      # Standard monocle3 format usually has one column
      if (ncol(closest_vertex) >= 1) {
        common_root <- intersect(root_cells, rownames(closest_vertex))
        if (length(common_root) > 0) {
          vals <- closest_vertex[common_root, 1, drop = TRUE]
          vals <- as.character(vals)
          vals <- vals[!is.na(vals)]
          if (length(vals) > 0) {
            root_node <- names(sort(table(vals), decreasing = TRUE))[1]
          }
        }
      }
    }
  }

  # If a graph node was found, use that location
  if (!is.na(root_node) && root_node %in% graph_bits$nodes$node) {
    out <- graph_bits$nodes %>%
      filter(node == root_node) %>%
      slice(1) %>%
      mutate(label = "ROOT")
    return(out)
  }

  # Fallback: use centroid of root cells in UMAP space
  umap_mat <- reducedDims(cds)$UMAP
  common_root <- intersect(root_cells, rownames(umap_mat))
  if (length(common_root) > 0) {
    xy <- umap_mat[common_root, , drop = FALSE]
    out <- data.frame(
      node = "root_centroid",
      x = median(xy[, 1], na.rm = TRUE),
      y = median(xy[, 2], na.rm = TRUE),
      label = "ROOT",
      stringsAsFactors = FALSE
    )
    return(out)
  }

  NULL
}

# Project all cells in a CDS onto the principal graph
project_cells_to_graph <- function(cds, edges_df, state_map) {
  umap_mat <- reducedDims(cds)$UMAP
  cell_names <- rownames(umap_mat)
  pt_vec <- pseudotime(cds)

  proj_res <- lapply(seq_along(cell_names), function(i) {
    px <- umap_mat[i, 1]
    py <- umap_mat[i, 2]

    best_dist2 <- Inf
    best_p <- NULL
    best_e <- 1L

    for (e in seq_len(nrow(edges_df))) {
      pr <- project_point_to_segment(
        px, py,
        edges_df$x[e], edges_df$y[e],
        edges_df$xend[e], edges_df$yend[e]
      )
      if (pr$dist2 < best_dist2) {
        best_dist2 <- pr$dist2
        best_p <- pr
        best_e <- e
      }
    }

    data.frame(
      cell = cell_names[i],
      umap_x = px,
      umap_y = py,
      graph_x = best_p$x,
      graph_y = best_p$y,
      projection_distance = sqrt(best_dist2),
      edge_index = best_e,
      state = factor(state_map[cell_names[i]], levels = state_order),
      pseudotime = as.numeric(pt_vec[cell_names[i]]),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(proj_res)
}

# Build a count-weighted ridge dataframe
# Ridge height = density * (n_state / max_state_n)
# This makes the tallest ridges correspond to the most abundant states.
build_weighted_ridge_df <- function(proj_df, state_order, n_points = 512) {
  df <- proj_df %>%
    filter(is.finite(pseudotime), !is.na(state)) %>%
    mutate(state = factor(as.character(state), levels = state_order))

  if (nrow(df) == 0) return(NULL)

  counts <- df %>%
    count(state, name = "n") %>%
    complete(state = factor(state_order, levels = state_order), fill = list(n = 0)) %>%
    mutate(weight = ifelse(max(n) > 0, n / max(n), 0),
           y = rev(seq_along(state_order)))

  pt_range <- range(df$pseudotime, na.rm = TRUE)
  if (!all(is.finite(pt_range))) return(NULL)

  # Expand tiny ranges so density() does not break
  if (diff(pt_range) == 0) {
    pt_range <- pt_range + c(-1e-6, 1e-6)
  }

  ridge_list <- lapply(seq_len(nrow(counts)), function(i) {
    st <- counts$state[i]
    n_st <- counts$n[i]
    weight_st <- counts$weight[i]
    y_st <- counts$y[i]

    sub <- df %>% filter(state == st)

    if (n_st < 2) {
      x_grid <- seq(pt_range[1], pt_range[2], length.out = n_points)
      dens_y <- rep(0, length(x_grid))
    } else {
      dens <- density(
        sub$pseudotime,
        from = pt_range[1],
        to = pt_range[2],
        n = n_points,
        na.rm = TRUE,
        bw = "nrd0"
      )
      x_grid <- dens$x
      dens_y <- dens$y
    }

    data.frame(
      state = factor(as.character(st), levels = state_order),
      x = x_grid,
      height = dens_y * weight_st,
      y = y_st,
      n_cells = n_st,
      rel_weight = weight_st,
      stringsAsFactors = FALSE
    )
  })

  ridge_df <- bind_rows(ridge_list)

  # Scale ridge heights to a nice visual range
  max_h <- max(ridge_df$height, na.rm = TRUE)
  if (is.finite(max_h) && max_h > 0) {
    ridge_df <- ridge_df %>%
      mutate(height = height / max_h * 0.9)
  }

  state_labels <- counts %>%
    mutate(
      label = paste0(as.character(state), " (n=", n, ")")
    ) %>%
    select(state, y, label)

  list(ridge_df = ridge_df, state_labels = state_labels, counts = counts)
}

####################
# Main Execution Loop
####################

message("Loading datasets...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states_path <- "Auto_final_states.rds"
state_vec <- readRDS(ifelse(file.exists(final_states_path), final_states_path, "Auto_topmp_v2_noreg_states_B.rds"))

# Align indices
common_cells <- intersect(Cells(tmdata_all), names(state_vec))
tmdata_all <- tmdata_all[, common_cells]
state_vec <- state_vec[common_cells]

files <- list.files(partA_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No Part A pseudotime RDS files found.")

pdf(pdf_path, width = 16, height = 12, onefile = TRUE)

for (f in files) {
  sample_id <- basename(f) %>% sub("_pseudotime.*", "", .)
  message("Processing sample: ", sample_id)

  # Load saved pseudotime to get cells belonging to this sample's model
  pt_saved <- readRDS(f)
  sample_cells <- intersect(names(pt_saved), Cells(tmdata_all))
  if (length(sample_cells) < 50) {
    message("  Skipping: < 50 overlapping cells.")
    next
  }

  # Prepare sample subset and trajectory
  sample_obj <- tmdata_all[, sample_cells]
  sample_obj$state_label <- factor(as.character(state_vec[sample_cells]), levels = target_states)

  root_cells <- colnames(sample_obj)[as.character(sample_obj$state_label) == root_state]
  if (length(root_cells) == 0) {
    message("  Skipping: No root cells for root-state.")
    next
  }

  cds <- tryCatch(
    prepare_sample_trajectory(sample_obj),
    error = function(e) {
      message("  Prep error: ", e$message)
      NULL
    }
  )
  if (is.null(cds)) next

  cds <- tryCatch(
    order_cells(cds, root_cells = root_cells),
    error = function(e) {
      message("  Order error: ", e$message)
      NULL
    }
  )
  if (is.null(cds)) next

  # Graph components
  graph_bits <- extract_graph_structure(cds)
  proj_df <- project_cells_to_graph(cds, graph_bits$edges, state_vec)
  root_label_df <- get_root_label_position(cds, graph_bits, root_cells)
  ridge_bits <- build_weighted_ridge_df(proj_df, state_order)

  # State legend labels with counts
  state_counts <- table(factor(proj_df$state, levels = target_states))
  legend_labels <- setNames(
    paste0(target_states, " (", as.integer(state_counts), ")"),
    target_states
  )

  # Dynamic text offset for root label
  x_span <- diff(range(c(graph_bits$nodes$x, proj_df$graph_x), na.rm = TRUE))
  y_span <- diff(range(c(graph_bits$nodes$y, proj_df$graph_y), na.rm = TRUE))
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1

  root_n <- sum(as.character(proj_df$state) == root_state, na.rm = TRUE)

  # Panel 1: Principal graph + projected cells by state + ROOT label
  p_proj_state <- ggplot() +
    geom_segment(
      data = graph_bits$edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.9, color = "grey55", alpha = 0.65
    ) +
    geom_point(
      data = proj_df,
      aes(x = graph_x, y = graph_y, color = state),
      size = 1.1, alpha = 0.85
    ) +
    {
      if (!is.null(root_label_df))
        geom_point(
          data = root_label_df,
          aes(x = x, y = y),
          shape = 8, size = 4.2, stroke = 1.0, color = "black"
        )
    } +
    {
      if (!is.null(root_label_df))
        geom_label(
          data = root_label_df,
          aes(x = x + 0.03 * x_span, y = y + 0.04 * y_span,
              label = paste0("ROOT\n", root_state, "\n(n=", root_n, ")")),
          size = 3.2, label.size = 0.25, hjust = 0, vjust = 0,
          fill = "white", color = "black", alpha = 0.95
        )
    } +
    scale_color_manual(values = group_cols, labels = legend_labels, drop = FALSE) +
    coord_equal() +
    theme_classic() +
    labs(
      title = "Principal Graph: State-projected cells",
      x = "Dim 1", y = "Dim 2", color = "State"
    )

  # Panel 2: Principal graph + projected cells by pseudotime + ROOT label
  p_proj_pt <- ggplot() +
    geom_segment(
      data = graph_bits$edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.9, color = "grey55", alpha = 0.65
    ) +
    geom_point(
      data = proj_df,
      aes(x = graph_x, y = graph_y, color = pseudotime),
      size = 1.1, alpha = 0.85
    ) +
    {
      if (!is.null(root_label_df))
        geom_point(
          data = root_label_df,
          aes(x = x, y = y),
          shape = 8, size = 4.2, stroke = 1.0, color = "black"
        )
    } +
    {
      if (!is.null(root_label_df))
        geom_label(
          data = root_label_df,
          aes(x = x + 0.03 * x_span, y = y + 0.04 * y_span, label = "ROOT"),
          size = 3.2, label.size = 0.25, hjust = 0, vjust = 0,
          fill = "white", color = "black", alpha = 0.95
        )
    } +
    scale_color_viridis_c(option = "D", na.value = "grey85") +
    coord_equal() +
    theme_classic() +
    labs(
      title = "Principal Graph: Pseudotime-projected cells",
      x = "Dim 1", y = "Dim 2", color = "Pseudotime"
    )

  # Panel 3: Raw UMAP colored by pseudotime
  p_umap_pt <- ggplot(
    proj_df %>% filter(is.finite(umap_x), is.finite(umap_y)),
    aes(x = umap_x, y = umap_y, color = pseudotime)
  ) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_color_viridis_c(option = "D", na.value = "grey85") +
    coord_equal() +
    theme_classic() +
    labs(
      title = "UMAP: Pseudotime",
      x = "UMAP 1", y = "UMAP 2", color = "Pseudotime"
    )

  # Panel 4: Count-weighted pseudotime density ridges
  if (!is.null(ridge_bits) && nrow(ridge_bits$ridge_df) > 0) {
    ridge_axis_breaks <- ridge_bits$state_labels$y
    ridge_axis_labels <- ridge_bits$state_labels$label

    p_ridges <- ggplot(
      ridge_bits$ridge_df,
      aes(x = x, y = y, height = height, group = state, fill = state)
    ) +
      geom_ridgeline(
        stat = "identity",
        scale = 1,
        alpha = 0.85,
        color = "white",
        linewidth = 0.35
      ) +
      scale_fill_manual(values = group_cols, drop = FALSE) +
      scale_y_continuous(
        breaks = ridge_axis_breaks,
        labels = ridge_axis_labels,
        expand = expansion(mult = c(0.03, 0.1))
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.y = element_blank()
      ) +
      labs(
        title = "Pseudotime density ridges (height weighted by cell count)",
        x = "Pseudotime"
      )
  } else {
    p_ridges <- ggplot() +
      theme_void() +
      labs(title = "Pseudotime density ridges") +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient finite pseudotime values", size = 5)
  }

  # Composite 2 x 2 report
  report_page <- (p_proj_state + p_proj_pt) /
    (p_umap_pt + p_ridges) +
    plot_annotation(
      title = paste0("Trajectory Report: ", sample_id),
      subtitle = paste0(
        "Total cells: ", nrow(proj_df),
        " | Root state: ", root_state,
        " | Root cells: ", root_n
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 11)
      )
    )

  print(report_page)

  # Save assets
  saveRDS(cds, file.path(asset_dir, paste0(sample_id, "_cds.rds")))
  write_csv(proj_df, file.path(asset_dir, paste0(sample_id, "_projections.csv")))
}

dev.off()
message("Reports complete. Output: ref_outs/", pdf_path)