####################
# Analysis registry:
#   Status: active
#   Script: analysis/publication/poster_section5_lineage_validation.R
#   Inputs:
#     ref_outs/Auto_six_state_markers/Auto_six_state_markers_final.csv
#     ref_outs/Auto_basal_smg_mp_signature_heatmap/Auto_basal_smg_mp_signature_heatmap_data.rds
#     ref_outs/state_distance_pseudotime/sample_state_trajectories/Alcindor_2025_SRR27335937_pseudotime_states.rds
#     ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_layout_principal_graph_geodesic_medoid.csv
#     ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_edges_principal_graph_geodesic_medoid.csv
#     ref_outs/visium_scatlas_states/Auto_visium_state_summary_by_sample.csv
#     ref_outs/visium_scatlas_states/Auto_P8_A_state_map.png
#     ref_outs/visium_scatlas_states/Auto_P8_B_state_map.png
#     updates/new_updates/summaries/Auto_clinical_assoc_state_boxplots_final_stats.csv
#   Outputs:
#     ref_outs/publication/section5/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section5_lineage_validation.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(scales)
  library(Seurat)
  library(patchwork)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section5"
out_dir <- pub_section_dir(section)
lineage_states <- c("Basal to Intestinal Metaplasia", "SMG-like Metaplasia")
mp_group_df <- tibble(mp = PUB_MP_ORDER, group = pub_mp_state(PUB_MP_ORDER)) |>
  mutate(group = factor(group, levels = PUB_MP_STATE_ORDER)) |>
  arrange(group, match(mp, PUB_MP_ORDER)) |>
  mutate(group_id = as.integer(group),
         x = row_number() + cumsum(c(0, diff(group_id) != 0)) * 0.7)

# ==============================================================================
# FIGURE 1: Signature bubble plot — lineage MPs show corresponding signatures
# ==============================================================================
cat("Generating signature bubble plot...\n")
sig_rds <- file.path(SCREF_REF_OUTS_DIR, "Auto_basal_smg_mp_signature_heatmap",
                      "Auto_basal_smg_mp_signature_heatmap_data.rds")
if (file.exists(sig_rds)) {
  sig <- readRDS(sig_rds)
  agg <- sig$agg_scores |>
    mutate(across(-top_mp_label, ~ as.numeric(scale(.)))) |>
    pivot_longer(-top_mp_label, names_to = "signature", values_to = "z")
    
  # Calculate percent of cells expressing signature (using score_norm > 0 as threshold)
  scores <- t(sig$score_norm)
  pct_expr <- cbind(top_mp_label = sig$cell_metadata$top_mp_label, as.data.frame(scores, check.names = FALSE)) |>
    group_by(top_mp_label) |>
    summarise(across(everything(), ~ mean(. > 0) * 100)) |>
    pivot_longer(-top_mp_label, names_to = "signature", values_to = "pct_expr")

  agg <- agg |> left_join(pct_expr, by = c("top_mp_label", "signature"))

  sig_order <- c("Squamous", "Gastric Columnar", "Intestinal Metaplasia",
                  "SMG-like Secretory", "Interferon Alpha Response",
                  "Interferon Gamma Response", "Buffa hypoxia A mean", "Cell-cycle score")
                  
  agg <- agg |>
    mutate(signature = factor(signature, levels = rev(sig_order[sig_order %in% unique(signature)])),
           mp = sub(" .*", "", as.character(top_mp_label)),
           state_group = factor(pub_mp_state(mp), levels = PUB_MP_STATE_ORDER),
           top_mp_label = factor(top_mp_label, levels = levels(sig$agg_scores$top_mp_label)))

  mp_group_df_sub <- agg |> select(mp, state_group, top_mp_label) |> distinct() |>
    arrange(state_group, match(mp, PUB_MP_ORDER)) |>
    mutate(group_id = as.integer(state_group),
           x = row_number() + cumsum(c(0, diff(group_id) != 0)) * 0.7)
           
  agg <- agg |> left_join(mp_group_df_sub |> select(mp, x), by = "mp")

  signature_levels <- levels(agg$signature)
  plot_sig_levels <- c(signature_levels, " ")

  agg <- agg |> mutate(signature = factor(as.character(signature), levels = plot_sig_levels))

  ann <- mp_group_df_sub |>
    mutate(signature = factor(" ", levels = plot_sig_levels))
           
  p <- ggplot(agg, aes(x, signature)) +
    geom_tile(data = ann, aes(x = x, y = signature),
              inherit.aes = FALSE, width = 0.92, height = 0.18,
              fill = PUB_MP_STATE_COLOURS[as.character(ann$state_group)], colour = NA) +
    geom_point(aes(size = pct_expr, colour = z), alpha = 0.95) +
    scale_colour_gradient(low = "white", high = "#B2182B", name = "Scaled\nscore") +
    scale_size_continuous(range = c(2.2, 9), name = "% Expressing\n(z > 0)") +
    scale_x_continuous(breaks = mp_group_df_sub$x,
                       labels = mp_group_df_sub$top_mp_label,
                       expand = expansion(mult = c(0.01, 0.01)),
                       position = "top") +
    scale_y_discrete(labels = function(x) ifelse(x == " ", "", x)) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    pub_theme(12) +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, size = 10))
    
  save_pub_gg(p, section, "s5_signature_bubble", width = 8, height = 5)
} else {
  abort_missing_figure(section, "s5_signature_bubble", "Signature bubble",
                   "Missing Auto_basal_smg_mp_signature_heatmap_data.rds.")
}

# ==============================================================================
# FIGURE 3: Pseudotime single-sample UMAP (Alcindor SRR27335937)
# ==============================================================================
cat("Generating pseudotime trajectory UMAP...\n")
sample_id <- "Alcindor_2025_SRR27335937"
sample_seurat_path <- file.path(SCREF_REF_OUTS_DIR, "by_samples", sample_id, paste0(sample_id, "_epi_f.rds"))

if (file.exists(sample_seurat_path)) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(monocle3)
    library(SeuratWrappers)
  })
  
  sub_obj <- readRDS(sample_seurat_path)
  states <- readRDS(SCREF_STATE_NOREG_RDS)
  
  common_cells <- intersect(colnames(sub_obj), names(states))
  sub_obj <- sub_obj[, common_cells]
  sub_obj$state_clean <- clean_state(states[common_cells])
  
  target_states <- setdiff(names(PUB_STATE_COLOURS), c("Unresolved", "Hybrid"))
  cells_keep <- colnames(sub_obj)[sub_obj$state_clean %in% target_states]
  sub_obj <- sub_obj[, cells_keep]
  sub_obj$state_label <- factor(sub_obj$state_clean, levels = target_states)
  
  # Exact pipeline from pseudotime_top_diverse_samples.R
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  n_pcs <- min(30, ncol(sub_obj) - 1)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  sub_obj <- RunUMAP(sub_obj, dims = 1:dims_use, verbose = FALSE)
  
  cds <- as.cell_data_set(sub_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  
  root_cells <- colnames(sub_obj)[sub_obj$state_label == "Basal to Intestinal Metaplasia"]
  cds <- order_cells(cds, root_cells = root_cells)
  
  p <- plot_cells(cds, color_cells_by = "state_label", show_trajectory_graph = TRUE,
                  label_cell_groups = FALSE, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
       scale_color_manual(values = PUB_STATE_COLOURS, name = "State") +
       labs(title = "Pseudotime trajectory (Alcindor SRR27335937)", x = NULL, y = NULL) +
       pub_theme(13) +
       theme(axis.text = element_blank(), axis.ticks = element_blank(),
             axis.title = element_blank(), legend.position = "right") +
       coord_fixed()
       
  save_pub_gg(p, section, "s5_pseudotime_umap", width = 7, height = 5.5)
} else {
  abort_missing_figure(section, "s5_pseudotime_umap", "Pseudotime UMAP",
                   "Missing Alcindor sample seurat object for monocle3 rebuilding.")
}

# ==============================================================================
# FIGURE 4: Geodesic medoid node plot
# ==============================================================================
cat("Generating geodesic nodeplot...\n")
layout_file <- file.path(SCREF_REF_OUTS_DIR, "task6_hybrid_pairwise_distance",
                          "Auto_task6_hybrid_pairwise_layout_principal_graph_geodesic_medoid.csv")
edge_file <- file.path(SCREF_REF_OUTS_DIR, "task6_hybrid_pairwise_distance",
                        "Auto_task6_hybrid_pairwise_edges_principal_graph_geodesic_medoid.csv")
if (file.exists(layout_file) && file.exists(edge_file)) {
  layout_df <- read_csv(layout_file, show_col_types = FALSE) |>
    mutate(state = clean_state(state)) |>
    filter(state %in% PUB_STATE_ORDER) |>
    mutate(
      label_x = case_when(
        state == "Classic Proliferative" ~ label_x + 0.05,
        state == "Basal to Intestinal Metaplasia" ~ label_x + 0.15,
        state == "SMG-like Metaplasia" ~ label_x + 0.15,
        TRUE ~ label_x
      ),
      label_y = case_when(
        state == "Classic Proliferative" ~ label_y - 0.10,
        state == "Basal to Intestinal Metaplasia" ~ label_y + 0.02,
        TRUE ~ label_y
      ),
      hjust = case_when(
        state == "Classic Proliferative" ~ 1,
        state == "Basal to Intestinal Metaplasia" ~ 0,
        state == "SMG-like Metaplasia" ~ 0,
        TRUE ~ hjust
      ),
      vjust = case_when(
        state == "Classic Proliferative" ~ 1,
        state == "Basal to Intestinal Metaplasia" ~ 0,
        TRUE ~ vjust
      )
    )
  edge_df <- read_csv(edge_file, show_col_types = FALSE) |>
    mutate(from = clean_state(from), to = clean_state(to)) |>
    filter(from %in% PUB_STATE_ORDER, to %in% PUB_STATE_ORDER)

  p <- ggplot() +
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct),
                 colour = "#94A3B8", alpha = 0.75, lineend = "round") +
    geom_text(data = edge_df,
              aes(x = (x + xend) / 2, y = (y + yend) / 2, label = paste0(round(pct, 1), "%")),
              size = 3.2, fontface = "bold", colour = "#111827") +
    geom_point(data = layout_df, aes(x, y, size = pct, fill = state),
               shape = 21, colour = "#111827", stroke = 0.6) +
    geom_text(data = layout_df,
              aes(label_x, label_y, label = state, hjust = hjust, vjust = vjust),
              size = 4, fontface = "bold", colour = "#111827") +
    scale_fill_manual(values = PUB_STATE_COLOURS, guide = "none") +
    scale_size_continuous(range = c(9, 18), guide = "none") +
    scale_linewidth_continuous(range = c(0.8, 4.2), guide = "none") +
    coord_equal(clip = "off") +
    labs(title = NULL) +
    theme_void(base_size = 13) +
    theme(plot.margin = margin(20, 20, 20, 20))
  save_pub_gg(p, section, "s5_geodesic_nodeplot", width = 7, height = 5.5)
} else {
  abort_missing_figure(section, "s5_geodesic_nodeplot", "Geodesic nodeplot",
                   "Missing task6_hybrid_pairwise_distance outputs.")
}

# ==============================================================================
# FIGURE 5: Visium stacked barplot, P8A/P8B state maps, and colocalisation score
# ==============================================================================
cat("Generating Visium plots...\n")
visium_summary <- file.path(SCREF_REF_OUTS_DIR, "visium_scatlas_states",
                            "Auto_visium_state_summary_by_sample.csv")
spot_file <- file.path(SCREF_REF_OUTS_DIR, "visium_scatlas_states", "Auto_visium_spot_annotations.csv.gz")

visium_colours <- c(
  PUB_STATE_COLOURS,
  "Hybrid" = "#000000",
  "Unresolved" = "#999999",
  "Normal/Mixed" = "#E6E6E6"
)

if (file.exists(visium_summary)) {
  vis <- read_csv(visium_summary, show_col_types = FALSE) |>
    mutate(state = clean_state(Auto_state_B)) |>
    mutate(state = ifelse(!state %in% c(PUB_STATE_ORDER, "Unresolved", "Hybrid", "Normal/Mixed"), "Normal/Mixed", state)) |>
    group_by(state) |>
    summarise(spots = sum(spots, na.rm = TRUE), .groups = "drop") |>
    mutate(pct = 100 * spots / sum(spots)) |>
    mutate(state = factor(state, levels = rev(c(PUB_STATE_ORDER, "Hybrid", "Unresolved", "Normal/Mixed"))))

  p_pie <- ggplot(vis, aes(x = 1, y = pct, fill = state)) +
    geom_bar(stat = "identity", width = 1, colour = "white", linewidth = 0.5) +
    geom_text(aes(x = 1.15, label = ifelse(pct > 1, paste0(round(pct, 1), "%"), ""),
                  color = ifelse(state == "Hybrid", "white", "black"),
                  group = state), 
              position = position_stack(vjust = 0.5), 
              size = 2.7, # Changed from 5.0 to make the text smaller
              fontface = "bold", 
              show.legend = FALSE) +
    scale_color_identity() +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = visium_colours, breaks = c(PUB_STATE_ORDER, "Hybrid", "Unresolved"), name = "State", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 14) +
    theme(legend.position = "right",
          plot.margin = margin(15, 15, 15, 15))
          
  save_pub_gg(p_pie, section, "s5_visium_pie_chart", width = 6.0, height = 4.5)

  if (file.exists(spot_file)) {

    spots <- read_csv(spot_file, show_col_types = FALSE) |>
      mutate(state = clean_state(Auto_state_B)) |>
      mutate(state = ifelse(!state %in% c(PUB_STATE_ORDER, "Unresolved", "Hybrid", "Normal/Mixed"), "Normal/Mixed", state)) |>
      mutate(state = factor(state, levels = c("Normal/Mixed", "Unresolved", "Hybrid", PUB_STATE_ORDER))) |>
      filter(sample %in% c("P8_A", "P8_B", "P8A", "P8B", "P8 A", "P8 B")) |>
      arrange(state)

    spots_a <- spots |> filter(sample %in% c("P8_A", "P8A", "P8 A"))
    p_a <- ggplot(spots_a, aes(pxl_col_in_fullres, -pxl_row_in_fullres, color=state)) +
      geom_point(size = 1.8, shape = 16) +
      scale_color_manual(values = visium_colours, name = "State", drop = FALSE) +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none")

    spots_b <- spots |> filter(sample %in% c("P8_B", "P8B", "P8 B"))
    p_b <- ggplot(spots_b, aes(pxl_col_in_fullres, -pxl_row_in_fullres, color=state)) +
      geom_point(size = 1.8, shape = 16) +
      scale_color_manual(values = visium_colours, name = "State", drop = FALSE) +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none")

    save_pub_gg(p_a, section, "s5_P8_A_state_map", width = 5, height = 5)
    save_pub_gg(p_b, section, "s5_P8_B_state_map", width = 5, height = 5)
  }
} else {
  abort_missing_figure(section, "s5_visium_stacked_bar", "Visium stacked bar",
                   "Missing Visium state summary CSV.")
}

if (file.exists(spot_file)) {
  spots <- read_csv(spot_file, show_col_types = FALSE) |>
    mutate(state = clean_state(Auto_state_B)) |>
    filter(state %in% PUB_STATE_ORDER, sample %in% c("P8_A", "P8_B", "P8A", "P8B", "P8 A", "P8 B"))
  if (nrow(spots) == 0) {
    spots <- read_csv(spot_file, show_col_types = FALSE) |>
      mutate(state = clean_state(Auto_state_B)) |>
      filter(state %in% PUB_STATE_ORDER)
  }
  coloc <- spots |>
    group_by(sample) |>
    group_modify(function(.x, .y) {
      coords <- as.matrix(.x[, c("array_row", "array_col")])
      if (nrow(coords) < 8) return(tibble())
      d <- as.matrix(dist(coords))
      diag(d) <- Inf
      neigh <- apply(d, 1, function(z) order(z)[1:min(6, length(z) - 1)])
      same <- vapply(seq_len(nrow(.x)), function(i) mean(.x$state[neigh[, i]] == .x$state[i], na.rm = TRUE), numeric(1))
      tibble(state = .x$state, same_neighbor_score = same)
    }) |>
    ungroup() |>
    mutate(state = factor(state, levels = PUB_STATE_ORDER))
    
  p <- ggplot(coloc, aes(x = state, y = same_neighbor_score, fill = state, color = state)) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      alpha = 0.8,
      linewidth = 0.6,
      color = "black"
    ) +
    geom_point(
      position = position_jitter(width = 0.15),
      alpha = 0.7,
      size = 1.5,
      stroke = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = PUB_STATE_COLOURS, guide = "none") +
    scale_color_manual(values = PUB_STATE_COLOURS, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Same-state neighbours") +
    pub_theme(14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 13, face = "bold"),
          plot.margin = margin(15, 15, 15, 15))
  save_pub_gg(p, section, "s5_visium_colocalisation", width = 6.0, height = 5.5)
}

cat("Section 5 complete.\n")

