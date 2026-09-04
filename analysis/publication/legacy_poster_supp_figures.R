####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/publication/legacy_poster_supp_figures.R
#   Methodology: not required (legacy poster assembly)
#   Inputs:
#     various
#   Outputs:
#     ref_outs/publication/supp_figures/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_supp_figures.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(igraph)
  library(ggraph)
  library(Seurat)
  library(patchwork)
  library(scales)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "supp_figures"
out_dir <- pub_section_dir(section)


# Extracted from poster_section1_atlas_metaprograms.R
# FIGURE 3: MP Correlation Heatmap (from final_mp_correlation.R)
# ==============================================================================
cat("Generating MP correlation heatmap...\n")
ucell_file <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "UCell_nMP19_filtered.rds")
if (file.exists(ucell_file)) {
  ucell <- readRDS(ucell_file)

  use_mps <- intersect(PUB_MP_ORDER, colnames(ucell))
  cor_mat <- cor(ucell[, use_mps, drop = FALSE], method = "spearman", use = "pairwise.complete.obs")
  mp_names <- setNames(gsub("\n", " ", short_mp_label(use_mps)), use_mps)
  rownames(cor_mat) <- mp_names[rownames(cor_mat)]
  colnames(cor_mat) <- mp_names[colnames(cor_mat)]

  state_vec <- factor(pub_mp_state(use_mps), levels = PUB_MP_STATE_ORDER)
  box_list <- split(seq_along(use_mps), state_vec)
  box_list <- box_list[lengths(box_list) > 0]

  ha_left <- rowAnnotation(
    State = state_vec,
    col = list(State = PUB_MP_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  ha_top <- HeatmapAnnotation(
    State = state_vec,
    col = list(State = PUB_MP_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  col_fun_cor <- colorRamp2(c(-0.4, 0, 0.4), c("#2166AC", "white", "#B2182B"))

  ht_cor <- Heatmap(
    cor_mat,
    name = "Spearman ρ",
    col = col_fun_cor,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 1),
    left_annotation = ha_left,
    top_annotation = ha_top,
    row_split = state_vec,
    column_split = state_vec,
    row_names_side = "right",
    column_names_side = "bottom",
    row_names_gp = gpar(fontsize = 8.5, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8.5, fontface = "bold"),
    column_names_rot = -30,
    row_title = NULL,
    column_title = NULL,
    width = unit(7.2, "inch"),
    height = unit(7.2, "inch"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 7.5, fontface = "bold"))
      for (idx in box_list) {
        if (i %in% idx && j %in% idx) {
          if (i == min(idx)) grid.lines(c(x - width / 2, x + width / 2), c(y + height / 2, y + height / 2), gp = gpar(col = "#111827", lwd = 1.5, lty = "dashed"))
          if (i == max(idx)) grid.lines(c(x - width / 2, x + width / 2), c(y - height / 2, y - height / 2), gp = gpar(col = "#111827", lwd = 1.5, lty = "dashed"))
          if (j == min(idx)) grid.lines(c(x - width / 2, x - width / 2), c(y - height / 2, y + height / 2), gp = gpar(col = "#111827", lwd = 1.5, lty = "dashed"))
          if (j == max(idx)) grid.lines(c(x + width / 2, x + width / 2), c(y - height / 2, y + height / 2), gp = gpar(col = "#111827", lwd = 1.5, lty = "dashed"))
        }
      }
    },
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                                labels_gp = gpar(fontsize = 9))
  )

  draw_expr <- quote(draw(ht_cor, merge_legend = TRUE, padding = unit(c(4, 18, 2, 2), "mm")))
  save_pub_grid(draw_expr, section, "s1_mp_correlation_heatmap", width = 8, height = 8)
} else {
  make_placeholder(section, "s1_mp_correlation_heatmap", "MP Correlation Heatmap", "Filtered MP scores missing")
}

# ==============================================================================

# Extracted from poster_section1_atlas_metaprograms.R
# FIGURE 5: Per cell heatmap (from state_definition_approach_b_reg_noreg.R)
# ==============================================================================
cat("Generating per-cell state heatmap...\n")
noreg_file <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_mp_adj.rds")
state_file <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_states_B.rds")

if (file.exists(noreg_file) && file.exists(state_file)) {
  mat_noncc <- readRDS(noreg_file)
  ucell_all <- readRDS(ucell_file)
  states <- readRDS(state_file)

  common <- Reduce(intersect, list(rownames(mat_noncc), rownames(ucell_all), names(states)))
  mat <- as.matrix(mat_noncc[common, , drop = FALSE])
  cc_mps <- intersect(c("MP1", "MP7", "MP9"), colnames(ucell_all))
  if (length(cc_mps) > 0) {
    cc_mat <- scale(as.matrix(ucell_all[common, cc_mps, drop = FALSE]))
    mat <- cbind(cc_mat, mat)
  }
  mat <- mat[, intersect(PUB_MP_ORDER, colnames(mat)), drop = FALSE]
  states <- states[common]

  valid_cells <- states %in% PUB_STATE_ORDER
  mat <- mat[valid_cells, ]
  states <- factor(states[valid_cells], levels = PUB_STATE_ORDER)

  set.seed(42)
  if (nrow(mat) > 6000) {
    samp <- unlist(lapply(split(seq_len(nrow(mat)), states), function(idx) {
      sample(idx, min(length(idx), ceiling(6000 * length(idx) / nrow(mat))))
    }), use.names = FALSE)
    mat <- mat[samp, ]
    states <- states[samp]
  }

  col_order <- unlist(lapply(levels(states), function(st) {
    idx <- which(states == st)
    if (length(idx) <= 2) return(idx)
    idx[hclust(dist(mat[idx, , drop = FALSE]), method = "ward.D2")$order]
  }), use.names = FALSE)
  mat <- mat[col_order, , drop = FALSE]
  states <- states[col_order]
  sub_scores <- t(mat)
  mp_group <- factor(pub_mp_state(rownames(sub_scores)), levels = PUB_MP_STATE_ORDER)
  col_ann <- HeatmapAnnotation(
    State = states,
    col = list(State = PUB_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = TRUE
  )
  row_ann <- rowAnnotation(
    MP_group = mp_group,
    col = list(MP_group = PUB_MP_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  ht_cell <- Heatmap(
    sub_scores,
    name = "Adj score",
    top_annotation = col_ann,
    left_annotation = row_ann,
    column_split = states,
    row_split = factor(ifelse(rownames(sub_scores) %in% cc_mps, "Cell-cycle MPs", "State MPs"),
                       levels = c("Cell-cycle MPs", "State MPs")),
    col = colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_labels = short_mp_label(rownames(sub_scores)),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
    show_column_names = FALSE,
    column_gap = unit(1.5, "mm"),
    row_gap = unit(2, "mm"),
    column_title = NULL,
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE
  )
  prop_df <- tibble(state = states) |>
    count(state) |>
    mutate(pct = n / sum(n), ypos = cumsum(pct) - pct / 2)
  pie <- ggplot(prop_df, aes(x = 1, y = pct, fill = state)) +
    geom_col(width = 1, colour = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = PUB_STATE_COLOURS, guide = "none") +
    theme_void() +
    geom_text(aes(y = ypos, label = paste0(round(100 * pct), "%")),
              x = 1.35, size = 3.2, fontface = "bold")
  draw_expr <- quote({
    pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2,
                                               widths = unit(c(0.82, 0.18), "npc"))))
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    draw(ht_cell, merge_legend = TRUE, newpage = FALSE)
    upViewport()
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    grid.draw(ggplotGrob(pie))
    upViewport(2)
  })
  save_pub_grid(draw_expr, section, "s1_percell_state_heatmap", width = 9, height = 5.6)
} else {
  make_placeholder(section, "s1_percell_state_heatmap", "Per-cell Heatmap", "State data missing")
}

# ==============================================================================

# Extracted from poster_section2_genetics_regulons.R
# FIGURE 1: CNA recurrent event consensus heatmap (clean)
# ==============================================================================
cna_rds <- file.path(SCREF_REF_OUTS_DIR, "Auto_cna_subclone_expression_v2", "rds", "Auto_v2_visualisation_results.rds")
if (!file.exists(cna_rds)) {
  make_placeholder(section, "s2_cna_consensus_heatmap", "CNA consensus heatmap",
                   "Missing Auto_v2_visualisation_results.rds.\nRun legacy_cna_subclone_expression_visuals_v2.R first.")
  make_placeholder(section, "s2_cna_mp_association", "CNA-MP association",
                   "Missing Auto_v2_visualisation_results.rds.")
} else {
  cna <- readRDS(cna_rds)

  # --- Consensus heatmap ---
  presence <- cna$event_presence |>
    mutate(event_short = gsub(" \\(.*", "", as.character(event_label)),
           event_short = gsub("Gain ", "+", event_short),
           event_short = gsub("Loss ", "-", event_short),
           event_status = case_when(
             !event_present ~ "Absent",
             grepl("^\\+", event_short) ~ "Gain",
             grepl("^-", event_short) ~ "Loss",
             TRUE ~ "Present"
           )) |>
    group_by(event_id) |>
    mutate(frac = mean(event_present, na.rm = TRUE)) |>
    ungroup() |>
    mutate(event_short = fct_reorder(event_short, frac, .desc = TRUE))

  # Simple subclone ordering by total event count
  sc_order <- presence |>
    group_by(subclone_id) |>
    summarise(n_events = sum(event_present, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(n_events)) |>
    pull(subclone_id)
  presence <- presence |> mutate(subclone_id = factor(subclone_id, levels = sc_order))

  cna_heat <- ggplot(presence, aes(subclone_id, event_short, fill = event_status)) +
    geom_tile(colour = "white", linewidth = 0.05) +
    scale_fill_manual(values = c("Absent" = "#F1F5F9", "Gain" = "#B2182B",
                                 "Loss" = "#2166AC", "Present" = "#64748B"),
                      breaks = c("Absent", "Gain", "Loss"), name = "CNA event") +
    labs(x = "Subclones", y = NULL) +
    pub_theme(13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "right")
  save_pub_gg(cna_heat, section, "s2_cna_consensus_heatmap", width = 9, height = 4.5)
  pub_copy_to_assets(file.path(out_dir, "figures", "s2_cna_consensus_heatmap.png"),
                     "s2_cna_state_consensus_heatmap.png", section)

  # ==============================================================================
  # FIGURE 2: CNA-MP association — showing no significant drivers
  # ==============================================================================
  tests <- cna$event_feature_tests |>
    filter(feature_group == "Metaprogrammes",
           feature %in% paste0("mp__", PUB_MP_ORDER)) |>
    mutate(MP = sub("^mp__", "", feature),
           MP = factor(MP, levels = PUB_MP_ORDER),
           event_short = gsub(" \\(.*", "", as.character(event_label)),
           event_short = gsub("Gain ", "+", event_short),
           event_short = gsub("Loss ", "-", event_short),
           event_short = factor(event_short, levels = levels(presence$event_short)),
           state_group = factor(PUB_MP_TO_STATE[as.character(MP)], levels = PUB_MP_STATE_ORDER),
           fdr = primary_p_adj_group,
           sig = ifelse(!is.na(fdr) & fdr < 0.05, "FDR < 0.05", "n.s."))

  delta_lim <- max(1.5, ceiling(max(abs(tests$primary_delta), na.rm = TRUE) * 3))
  cna_assoc <- ggplot(tests, aes(MP, fct_rev(event_short), fill = primary_delta)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    geom_point(aes(shape = sig), size = 2.2, colour = "#111827") +
    scale_x_discrete(labels = function(x) gsub("\n", " ", short_mp_label(x))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-delta_lim, delta_lim),
                         name = "Paired\nmedian Δ") +
    scale_shape_manual(values = c("FDR < 0.05" = 8, "n.s." = 1), name = "Significance") +
    labs(x = NULL, y = "Recurrent CNA event") +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), legend.position = "right")
  save_pub_gg(cna_assoc, section, "s2_cna_mp_association", width = 9, height = 5)
  pub_copy_to_assets(file.path(out_dir, "figures", "s2_cna_mp_association.png"),
                     "s2_cna_recurrent_mp_association.png", section)
  write_csv(tests, file.path(out_dir, "tables", "s2_cna_mp_association.csv"))
}

# ==============================================================================

# Extracted from poster_section2_genetics_regulons.R
# FIGURE 3: SCENIC Top 100 regulon network (placeholder if not available)
# ==============================================================================
scenic_dir <- file.path(SCREF_REF_OUTS_DIR, "final_mp_scenic")
network_csv <- list.files(scenic_dir, pattern = "regulon|network|edge", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
network_csv <- network_csv[grepl("\\.csv$", network_csv, ignore.case = TRUE)]

# Also check for the output/regulon tables
output_dir <- file.path(scenic_dir, "output")
if (dir.exists(output_dir)) {
  extra <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  network_csv <- c(network_csv, extra)
}

edge_file <- network_csv[1]

if (!is.na(edge_file) && file.exists(edge_file)) {
  edges <- read_csv(edge_file, show_col_types = FALSE)
  if (all(c("from", "to") %in% colnames(edges))) {
    # Exclude 3CA_EMT state
    edges <- edges |>
      filter(!grepl("3CA_EMT", from, ignore.case = TRUE),
             !grepl("3CA_EMT", to, ignore.case = TRUE))

    weight_col <- intersect(c("weight", "score", "importance", "n_regulons"), colnames(edges))[1]
    if (!is.na(weight_col)) {
      edges <- edges |> mutate(weight = .data[[weight_col]])
    } else {
      edges$weight <- 1
    }
    edges <- edges |>
      slice_max(abs(weight), n = 100, with_ties = FALSE)

    graph <- graph_from_data_frame(edges, directed = FALSE)
    node_names <- V(graph)$name
    is_state <- node_names %in% PUB_STATE_ORDER
    V(graph)$node_type <- ifelse(is_state, "State", "Regulon")
    V(graph)$state <- ifelse(is_state, node_names, NA_character_)
    V(graph)$degree_val <- degree(graph)

    # Label: states always, TFs only top degree
    top_tf_n <- 15
    tf_nodes <- which(!is_state)
    if (length(tf_nodes) > 0) {
      tf_deg <- degree(graph, v = tf_nodes)
      top_tfs <- names(sort(tf_deg, decreasing = TRUE))[seq_len(min(top_tf_n, length(tf_deg)))]
    } else {
      top_tfs <- character(0)
    }
    V(graph)$label_show <- ifelse(is_state | node_names %in% top_tfs, node_names, "")

    regulon_plot <- ggraph(graph, layout = "fr") +
      geom_edge_link(aes(alpha = abs(weight)), colour = "#64748B",
                     linewidth = 0.3, show.legend = FALSE) +
      geom_node_point(aes(size = ifelse(node_type == "State", 14, 3 + degree_val * 0.5),
                          fill = state),
                      shape = 21, colour = "#111827", stroke = 0.4, alpha = 0.92) +
      geom_node_text(aes(label = label_show), repel = TRUE,
                     size = ifelse(is_state[order(V(graph))], 5, 3.5),
                     fontface = ifelse(is_state[order(V(graph))], "bold", "plain"),
                     colour = "#111827", max.overlaps = 20) +
      scale_size_identity() +
      scale_fill_manual(values = c(PUB_STATE_COLOURS, "NA" = "#94A3B8"),
                        na.value = "#94A3B8", guide = "none") +
      labs(title = NULL) +
      theme_void(base_size = 13) +
      theme(plot.title = element_blank())
    save_pub_gg(regulon_plot, section, "s2_scenic_regulon_network", width = 8, height = 6.5)
  } else {
    make_placeholder(section, "s2_scenic_regulon_network", "SCENIC regulon network",
                     "Edge CSV found but missing from/to columns.", width = 8, height = 6.5)
  }
} else {
  make_placeholder(section, "s2_scenic_regulon_network", "SCENIC regulon network",
                   "SCENIC still running.\nPlaceholder — rerun after completion.", width = 8, height = 6.5)
  write_pub_status(section, "s2_scenic_regulon_network", "scenic_pending",
                   "No network CSV found in final_mp_scenic/")
}

cat("Section 2 complete.\n")


# Extracted from poster_section3_tme_interactions.R
# FIGURE 1: TME-cancer MP interaction dotmap with relaxed threshold
# ==============================================================================
cat("Generating TME-cancer MP dotmap...\n")

dot <- read_csv(dot_file, show_col_types = FALSE) |>
  mutate(mp1_name = as.character(mp1_name),
         cancer_mp = ifelse(focal_celltype == "cancer", mp1_name, mp2_name),
         tme_celltype = ifelse(focal_celltype == "cancer", partner_celltype, focal_celltype),
         tme_mp = ifelse(focal_celltype == "cancer", mp2_name, mp1_name),
         cancer_label = SCREF_MP_DESCRIPTIONS[cancer_mp],
         cancer_label = ifelse(is.na(cancer_label), cancer_mp,
                                paste0(cancer_mp, " ", cancer_label))) |>
  filter(cancer_mp %in% PUB_MP_ORDER)

# Relaxed threshold: p < 0.30 or |rho| >= 0.15 (lowered from original)
dot <- dot |>
  mutate(significant = spearman_p < 0.05,
         relaxed = spearman_p < 0.30 | abs(spearman_r) >= 0.15) |>
  filter(relaxed) |>
  group_by(tme_celltype, cancer_mp) |>
  slice_max(abs(spearman_r), n = 1, with_ties = FALSE) |>
  ungroup()

tme_order <- c("fibroblast", "macrophage", "endothelial", "cd8", "cd4", "nk", "plasma")
dot <- dot |>
  mutate(tme_celltype = factor(tme_celltype, levels = tme_order[tme_order %in% unique(tme_celltype)]),
         cancer_mp = factor(cancer_mp, levels = PUB_MP_ORDER),
         state_group = factor(pub_mp_state(as.character(cancer_mp)), levels = PUB_MP_STATE_ORDER),
         sig_class = ifelse(significant, "Significant (P < 0.05)", "Trend"))

dot_plot <- ggplot(dot, aes(cancer_mp, tme_celltype)) +
  geom_point(aes(size = co_positive_sample_pct, fill = spearman_r, alpha = sig_class),
             shape = 21, colour = "#111827", stroke = 0.35) +
  scale_x_discrete(labels = function(x) gsub("\n", " ", short_mp_label(x))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Spearman ρ") +
  scale_size_continuous(range = c(2, 8), name = "Co-positive\nsamples (%)") +
  scale_alpha_manual(values = c("Significant (P < 0.05)" = 1.0, "Trend" = 0.3),
                     name = "Evidence") +
  labs(x = "Cancer metaprogramme", y = "TME compartment") +
  pub_theme(13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right")
save_pub_gg(dot_plot, section, "s3_tme_dotmap", width = 9.5, height = 5)
write_csv(dot, file.path(out_dir, "tables", "s3_tme_dotmap_data.csv"))

# ==============================================================================

# Extracted from poster_section3_tme_interactions.R
# FIGURE 2: Simplified LR-annotated positive network
# ==============================================================================
if (!file.exists(lr_file)) {
  make_placeholder(section, "s3_lr_network", "LR network",
                   "Missing Auto_celltype_positive_edges_lr_annotated.csv.")
  quit(save = "no")
}

cat("Generating LR network...\n")

lr <- read_csv(lr_file, show_col_types = FALSE) |>
  filter(compartment1 == "cancer" | compartment2 == "cancer") |>
  mutate(cancer_mp = ifelse(compartment1 == "cancer", mp1_name, mp2_name),
         cancer_label = ifelse(cancer_mp %in% names(SCREF_MP_DESCRIPTIONS),
                                paste0(cancer_mp, " ", SCREF_MP_DESCRIPTIONS[cancer_mp]),
                                cancer_mp),
         tme_label = ifelse(compartment1 == "cancer", node2_label, node1_label),
         lr_label = ifelse(is.na(top_lr_label), "", top_lr_label)) |>
  filter(cancer_mp %in% PUB_MP_ORDER)

# Score edges and take top ones for clean network
lr <- lr |>
  mutate(score = spearman_significance + pmin(n_lr_pairs, 20) / 4 + abs(spearman_r) * 2) |>
  arrange(desc(score)) |>
  slice_head(n = 20)  # Fewer edges for clarity

if (nrow(lr) == 0) {
  make_placeholder(section, "s3_lr_network", "LR network",
                   "No cancer-TME LR edges passed filter.")
} else {
  edges <- lr |>
    transmute(from = cancer_label, to = tme_label,
              weight = spearman_r, lr_label = lr_label, n_lr = n_lr_pairs)
  graph <- graph_from_data_frame(edges, directed = FALSE)
  node_names <- V(graph)$name
  V(graph)$type <- ifelse(node_names %in% edges$from, "Cancer MP", "TME MP")

  # Assign colours by state for cancer MPs
  mp_in_name <- sub(" .*", "", node_names)
  V(graph)$state_col <- ifelse(
    V(graph)$type == "Cancer MP" & mp_in_name %in% names(PUB_MP_TO_STATE),
    PUB_MP_STATE_COLOURS[PUB_MP_TO_STATE[mp_in_name]],
    "#F59E0B"  # TME colour
  )
  V(graph)$state_col[is.na(V(graph)$state_col)] <- "#94A3B8"

  # Show only a few key LR annotations to avoid overcrowding
  key_lr <- unique(edges$lr_label[nzchar(edges$lr_label)])[seq_len(min(2, sum(nzchar(unique(edges$lr_label)))))]
  edges$edge_label <- ifelse(edges$lr_label %in% key_lr, edges$lr_label, "")

  set.seed(42)
  coords <- as.data.frame(igraph::layout_with_fr(graph))
  colnames(coords) <- c("x", "y")
  coords$name <- V(graph)$name
  coords$type <- V(graph)$type
  coords$state_col <- V(graph)$state_col
  edge_coords <- edges |>
    left_join(coords |> select(from = name, x, y), by = "from") |>
    left_join(coords |> select(to = name, xend = x, yend = y), by = "to") |>
    mutate(label_x = (x + xend) / 2,
           label_y = (y + yend) / 2,
           edge_label = ifelse(lr_label %in% key_lr, lr_label, ""))

  network_plot <- ggplot() +
    geom_segment(data = edge_coords,
                 aes(x = x, y = y, xend = xend, yend = yend, linewidth = pmax(weight, 0.05)),
                 colour = "#64748B", alpha = 0.55, lineend = "round") +
    geom_text(data = edge_coords |> filter(nzchar(edge_label)),
              aes(label_x, label_y, label = edge_label),
              size = 3, colour = "#111827", fontface = "bold", check_overlap = TRUE) +
    geom_point(data = coords,
               aes(x, y, size = ifelse(type == "Cancer MP", 10, 5)),
               fill = coords$state_col, shape = 21,
               colour = "#111827", stroke = 0.4) +
    geom_text(data = coords,
              aes(x, y, label = name, fontface = ifelse(type == "Cancer MP", "bold", "plain")),
              size = 3.2, colour = "#111827", check_overlap = TRUE, vjust = -0.8) +
    scale_size_identity() +
    scale_linewidth_continuous(range = c(0.35, 2.2), guide = "none") +
    labs(title = NULL, caption = NULL) +
    theme_void(base_size = 13) +
    theme(plot.title = element_blank(), plot.caption = element_blank())
  save_pub_gg(network_plot, section, "s3_lr_network", width = 8, height = 6)
  write_csv(lr, file.path(out_dir, "tables", "s3_lr_network_edges.csv"))
}

cat("Section 3 complete.\n")


# Extracted from poster_section4_pdo_concordance.R
# FIGURE 1: PDO NMF composition heatmap with dotted selection squares
# ==============================================================================
cat("Generating PDO NMF heatmap...\n")
pdo_nmf_file <- file.path(pdo_dir, "Metaprogrammes_Results", "geneNMF_metaprograms_nMP_13.rds")
if (file.exists(pdo_nmf_file)) {
  pdo_nmf <- readRDS(pdo_nmf_file)
  sil <- pdo_nmf$metaprograms.metrics$silhouette
  names(sil) <- names(pdo_nmf$metaprograms.genes)[seq_along(sil)]
  retained <- names(sil)[sil >= 0]
  programme_order <- pdo_nmf$programs.tree$order
  sim <- pdo_nmf$programs.similarity[programme_order, programme_order]
  clusters <- paste0("MP", pdo_nmf$programs.clusters[programme_order])
  runs <- rle(clusters)
  run_end <- cumsum(runs$lengths)
  retained_runs <- data.frame(mp = runs$values, start = run_end - runs$lengths + 1,
                              end = run_end, stringsAsFactors = FALSE) |>
    filter(mp %in% retained)
  n_prog <- nrow(sim)
  ht <- Heatmap(
    sim,
    name = "Jaccard\nsimilarity",
    col = colorRamp2(c(0, 0.5, 1), c("#F8FAFC", "#60A5FA", "#111827")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,
    column_title = NULL,
    use_raster = TRUE,
    raster_quality = 3,
    border = FALSE
  )
  draw_expr <- quote({
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
    decorate_heatmap_body("Jaccard\nsimilarity", {
      for (rr in seq_len(nrow(retained_runs))) {
        s <- retained_runs$start[rr]; e <- retained_runs$end[rr]
        centre <- (s + e - 1) / (2 * n_prog)
        extent <- (e - s + 1) / n_prog
        grid.rect(x = unit(centre, "npc"), y = unit(1 - centre, "npc"),
                  width = unit(extent, "npc"), height = unit(extent, "npc"),
                  gp = gpar(fill = NA, col = "#111827", lwd = 1.4, lty = "dashed"))
      }
    })
  })
  save_pub_grid(draw_expr, section, "s4_pdo_nmf_heatmap", width = 7.2, height = 7.0)
} else {
  make_placeholder(section, "s4_pdo_nmf_heatmap", "PDO NMF heatmap",
                   "Missing geneNMF_metaprograms_nMP_13.rds.")
}

# ==============================================================================

# Extracted from poster_section4_pdo_concordance.R
# FIGURE 2: PDO-to-scRef MP relationship heatmap (ordered by state)
# ==============================================================================
cat("Generating PDO MP correlation heatmap...\n")
if (file.exists(pdo_nmf_file) && file.exists(SCREF_LEGACY_MP_OBJECT_RDS)) {
  pdo_nmf <- readRDS(pdo_nmf_file)
  sc_nmf <- readRDS(SCREF_LEGACY_MP_OBJECT_RDS)
  pdo_sil <- pdo_nmf$metaprograms.metrics$silhouette
  names(pdo_sil) <- paste0("MP", seq_along(pdo_sil))
  pdo_use <- pdo_mp_order[pdo_mp_order %in% names(pdo_sil) & pdo_sil[pdo_mp_order] >= 0]
  sc_use <- PUB_MP_ORDER[PUB_MP_ORDER %in% names(sc_nmf$metaprograms.genes)]
  overlap <- outer(pdo_use, sc_use, Vectorize(function(a, b) {
    length(intersect(pdo_nmf$metaprograms.genes[[a]], sc_nmf$metaprograms.genes[[b]])) /
      length(union(pdo_nmf$metaprograms.genes[[a]], sc_nmf$metaprograms.genes[[b]]))
  }))
  dimnames(overlap) <- list(pdo_use, sc_use)
  ha_top <- HeatmapAnnotation(
    State = factor(pub_mp_state(sc_use), levels = PUB_MP_STATE_ORDER),
    col = list(State = PUB_MP_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  ha_left <- rowAnnotation(
    State = factor(pdo_mp_to_state[pdo_use], levels = PUB_MP_STATE_ORDER),
    col = list(State = PUB_MP_STATE_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  ht <- Heatmap(
    overlap,
    name = "Gene-set\nJaccard",
    col = colorRamp2(c(0, max(overlap, na.rm = TRUE) / 2, max(overlap, na.rm = TRUE)),
                     c("white", "#FCA5A5", "#B2182B")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    row_labels = paste0(pdo_use, " ", pdo_mp_desc[pdo_use]),
    column_labels = short_mp_label(sc_use),
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = -30,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", overlap[i, j]), x, y, gp = gpar(fontsize = 7))
    }
  )
  save_pub_grid(quote(draw(ht, padding = unit(c(3, 16, 2, 2), "mm"))),
                section, "s4_pdo_mp_correlation", width = 9, height = 5.8)
} else {
  make_placeholder(section, "s4_pdo_mp_correlation", "PDO MP correlation",
                   "Missing UCell_scores_filtered.rds or NMF object.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 1: Marker specificity bubble plot for two lineage states
# ==============================================================================
cat("Generating marker specificity bubble plot...\n")
marker_file <- file.path(SCREF_REF_OUTS_DIR, "Auto_six_state_markers", "Auto_six_state_markers_final.csv")
if (file.exists(marker_file)) {
  markers <- read_csv(marker_file, show_col_types = FALSE) |>
    mutate(state = clean_state(state)) |>
    filter(state %in% lineage_states, best_state_match) |>
    group_by(state) |>
    slice_max(ranking_score, n = 8, with_ties = FALSE) |>
    ungroup()

  plot_df <- markers |>
    mutate(state = factor(state, levels = lineage_states),
           gene = factor(gene, levels = rev(unlist(lapply(lineage_states, function(st) {
             markers |> filter(state == st) |> arrange(desc(ranking_score)) |> pull(gene)
           })))))

  p <- ggplot(plot_df, aes(state, gene)) +
    geom_point(aes(size = global_pct_state, colour = global_mean_diff), alpha = 0.92) +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           name = "Mean\nexpression\ndifference") +
    scale_size_continuous(range = c(3, 9), labels = percent_format(accuracy = 1),
                          name = "Cells\nexpressing") +
    labs(x = NULL, y = NULL) +
    pub_theme(13) +
    facet_grid(state ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 11),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          panel.spacing.y = unit(3, "mm"))
  save_pub_gg(p, section, "s5_marker_specificity_bubble", width = 6.5, height = 5.5)
  write_csv(markers, file.path(out_dir, "tables", "s5_marker_specificity.csv"))
} else {
  make_placeholder(section, "s5_marker_specificity_bubble", "Marker bubble",
                   "Missing Auto_six_state_markers_final.csv.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 2: Signature bubble plot — lineage MPs show corresponding signatures
# ==============================================================================
cat("Generating signature bubble plot...\n")
sig_rds <- file.path(SCREF_REF_OUTS_DIR, "Auto_basal_smg_mp_signature_heatmap",
                      "Auto_basal_smg_mp_signature_heatmap_data.rds")
if (file.exists(sig_rds)) {
  sig <- readRDS(sig_rds)
  agg <- sig$agg_scores |>
    pivot_longer(-top_mp_label, names_to = "signature", values_to = "z")

  # Clean signature order
  sig_order <- c("Squamous", "Gastric Columnar", "Intestinal Metaplasia",
                  "SMG-like Secretory", "Interferon Alpha Response",
                  "Interferon Gamma Response", "Buffa hypoxia A mean", "Cell-cycle score")
  agg <- agg |>
    mutate(signature = factor(signature, levels = sig_order[sig_order %in% unique(signature)]),
           mp = sub(" .*", "", as.character(top_mp_label)),
           state_group = factor(pub_mp_state(mp), levels = PUB_MP_STATE_ORDER),
           top_mp_label = factor(top_mp_label, levels = rev(levels(sig$agg_scores$top_mp_label))))

  p <- ggplot(agg, aes(signature, top_mp_label)) +
    geom_point(aes(size = pmin(abs(z), 2.5), colour = z), alpha = 0.95) +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           name = "Scaled\nscore") +
    scale_size_continuous(range = c(2.2, 9), name = "|Score|") +
    labs(x = NULL, y = NULL) +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10))
  save_pub_gg(p, section, "s5_signature_bubble", width = 8, height = 5)
} else {
  make_placeholder(section, "s5_signature_bubble", "Signature bubble",
                   "Missing Auto_basal_smg_mp_signature_heatmap_data.rds.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 3: Pseudotime single-sample UMAP (Alcindor SRR27335937)
# ==============================================================================
cat("Generating pseudotime trajectory UMAP...\n")

sample_id <- "Alcindor_2025_SRR27335937"
epi <- readRDS(file.path(SCREF_REF_OUTS_DIR, "EAC_Ref_epi.rds"))
states <- readRDS(SCREF_STATE_NOREG_RDS)
epi$state_B <- states[colnames(epi)]

sample_cells <- rownames(epi@meta.data)[epi$orig.ident == sample_id]
sample_cells <- intersect(sample_cells, names(states))

primary_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia", 
                    "Stress-adaptive", "SMG-like Metaplasia", "Immune Infiltrating")
sample_cells <- sample_cells[states[sample_cells] %in% primary_states]

if (length(sample_cells) > 0) {
  sub_obj <- epi[, sample_cells]
  sample_meta <- data.frame(
    cell = sample_cells,
    state = states[sample_cells],
    stringsAsFactors = FALSE
  )
  
  # Ensure all necessary packages for Monocle3 are loaded
  suppressPackageStartupMessages({
    library(monocle3)
    library(patchwork)
    library(SeuratWrappers)
  })
  
  # Seurat prep
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  n_pcs <- min(30, max(2, ncol(sub_obj) - 1))
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  sub_obj <- RunUMAP(sub_obj, dims = 1:dims_use, verbose = FALSE)
  
  cds <- as.cell_data_set(sub_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  
  root_cells <- sample_meta$cell[sample_meta$state == "Basal to Intestinal Metaplasia"]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.matrix(closest_vertex)) closest_vertex <- closest_vertex[, 1]
  root_cells <- intersect(root_cells, names(closest_vertex))
  
  if (length(root_cells) > 0) {
    root_vertex_raw <- names(sort(table(as.character(closest_vertex[root_cells])), decreasing = TRUE))[1]
    
    graph_obj <- principal_graph(cds)[["UMAP"]]
    graph_nodes <- igraph::V(graph_obj)$name
    root_vertex <- as.character(root_vertex_raw)[1]
    if (!is.na(suppressWarnings(as.numeric(root_vertex))) && !(root_vertex %in% graph_nodes)) {
      root_vertex <- graph_nodes[as.integer(root_vertex)]
    }
    
    cds <- order_cells(cds, root_pr_nodes = root_vertex)
  }
  
  # Group colours for primary states exactly as in state distance pseudotime script
  group_cols <- c(
    "Classic Proliferative" = "#E41A1C",
    "Basal to Intestinal Metaplasia" = "#4DAF4A",
    "Stress-adaptive" = "#984EA3",
    "SMG-like Metaplasia" = "#FF7F00",
    "Immune Infiltrating" = "#377EB8"
  )
  
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
  
  p <- p_states + p_pseudotime + plot_layout(guides = "collect")
  
  # Use 14x6 to match the two-plot PDF dimensions used in the source script
  save_pub_gg(p, section, "s5_pseudotime_umap", width = 14, height = 6)
  
  # Cleanup memory
  rm(epi, sub_obj, cds)
  gc()
} else {
  make_placeholder(section, "s5_pseudotime_umap", "Pseudotime UMAP",
                   "No cells found for Alcindor_2025_SRR27335937.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
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
    filter(state %in% PUB_STATE_ORDER)
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
    coord_equal() +
    labs(title = NULL) +
    theme_void(base_size = 13) +
    theme(plot.margin = margin(12, 30, 12, 30))
  save_pub_gg(p, section, "s5_geodesic_nodeplot", width = 7, height = 5.5)
} else {
  make_placeholder(section, "s5_geodesic_nodeplot", "Geodesic nodeplot",
                   "Missing task6_hybrid_pairwise_distance outputs.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 5: Tumour location — Distal vs GOJ for lineage states
# ==============================================================================
cat("Generating tumour location plot...\n")
clin_file <- file.path(SCREF_SUMMARY_DIR, "Auto_clinical_assoc_state_boxplots_final_stats.csv")
if (file.exists(clin_file)) {
  clin <- read_csv(clin_file, show_col_types = FALSE) |>
    filter(clinical_variable == "Tumor Location",
           feature %in% lineage_states)

  # Parse medians from group_summary
  parse_median <- function(summary_str, group) {
    m <- regmatches(summary_str, regexpr(paste0(group, " \\(n=[0-9]+, median=([0-9.]+)\\)"),
                                          summary_str))
    if (length(m) == 0) return(NA_real_)
    as.numeric(sub(".*median=([0-9.]+)\\).*", "\\1", m))
  }

  clin_plot <- clin |>
    rowwise() |>
    mutate(Distal = parse_median(group_summary, "Distal Oesophagus"),
           GOJ = parse_median(group_summary, "GOJ")) |>
    ungroup() |>
    select(feature, p_value, p_adj, significance, Distal, GOJ) |>
    pivot_longer(c(Distal, GOJ), names_to = "location", values_to = "median_pct") |>
    mutate(location = recode(location, Distal = "Distal\noesophagus"),
           sig_label = ifelse(location == "Distal\noesophagus", pub_sig_label(p_adj), ""),
           feature = factor(feature, levels = lineage_states))

  p <- ggplot(clin_plot, aes(feature, median_pct, fill = location)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.55,
             colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sig_label), position = position_dodge(width = 0.65),
              vjust = -0.5, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("Distal\noesophagus" = "#B2182B", "GOJ" = "#2166AC"),
                      name = "Location") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Median state abundance (%)") +
    pub_theme(13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_pub_gg(p, section, "s5_tumour_location", width = 6.2, height = 4.5)
  write_csv(clin_plot, file.path(out_dir, "tables", "s5_tumour_location.csv"))
} else {
  make_placeholder(section, "s5_tumour_location", "Tumour location",
                   "Missing clinical stats CSV.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 6: Visium stacked barplot + copy P8A/P8B spatial maps
# ==============================================================================
cat("Generating Visium plots...\n")
visium_summary <- file.path(SCREF_REF_OUTS_DIR, "visium_scatlas_states",
                             "Auto_visium_state_summary_by_sample.csv")
if (file.exists(visium_summary)) {
  vis <- read_csv(visium_summary, show_col_types = FALSE) |>
    mutate(state = clean_state(Auto_state_B)) |>
    group_by(state) |>
    summarise(spots = sum(spots, na.rm = TRUE), .groups = "drop") |>
    mutate(pct = 100 * spots / sum(spots),
           sample = "Combined Visium",
           state = factor(state, levels = rev(c(PUB_STATE_ORDER, "Unresolved", "Hybrid", "Normal/Mixed"))))

  p <- ggplot(vis, aes(sample, pct, fill = state)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(PUB_STATE_COLOURS, "Unresolved" = "grey80",
                                 "Hybrid" = "#111827", "Normal/Mixed" = "#CBD5E1"),
                      name = "State", drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Assigned spots (%)") +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10))
  save_pub_gg(p, section, "s5_visium_stacked_bar", width = 7, height = 4.5)

  # Copy P8A and P8B spatial maps
  visium_dir <- file.path(SCREF_REF_OUTS_DIR, "visium_scatlas_states")
  for (f in c("Auto_P8_A_state_map.png", "Auto_P8_B_state_map.png")) {
    src <- file.path(visium_dir, f)
    if (file.exists(src)) {
      dest_name <- sub("Auto_", "s5_", f)
      file.copy(src, file.path(PUB_ASSET_DIR, dest_name), overwrite = TRUE)
      file.copy(src, file.path(POSTER_ASSET_DIR, dest_name), overwrite = TRUE)
      file.copy(src, file.path(out_dir, "figures", dest_name), overwrite = TRUE)
    }
  }
} else {
  make_placeholder(section, "s5_visium_stacked_bar", "Visium stacked bar",
                   "Missing Visium state summary CSV.")
}

# ==============================================================================

# Extracted from poster_section5_lineage_validation.R
# FIGURE 7: Lineage state survival (placeholder — needs TCGA results)
# ==============================================================================
# This will be filled by the TCGA survival publication replotting
survival_csv <- file.path(SCREF_REF_OUTS_DIR, "task2_survival", "Auto_task2_survival_mp_state_cox_methods_splits.csv")

if (file.exists(survival_csv)) {
  cox <- read_csv(survival_csv, show_col_types = FALSE) |>
    filter(feature %in% lineage_states,
           cohort == "EAC",
           mode == "noreg",
           split_method == "continuous")

  # Use consistent method across all survival plotting
  best_method <- "whole_tcga"

  cox_best <- cox |> filter(method == best_method)

  if (nrow(cox_best) > 0) {
    cox_best <- cox_best |>
      mutate(feature = factor(feature, levels = lineage_states),
             log2HR = log2(HR),
             sig = pub_sig_label(P_value),
             direction = ifelse(HR < 1, "Protective", "Risk"))

    p <- ggplot(cox_best, aes(feature, log2HR, fill = feature)) +
      geom_col(width = 0.5, colour = "#111827", linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_text(aes(label = paste0(sig, "\nHR=", round(HR, 2))),
                vjust = ifelse(cox_best$log2HR > 0, -0.3, 1.3), size = 4) +
      scale_fill_manual(values = PUB_STATE_COLOURS[lineage_states], guide = "none") +
      labs(x = NULL, y = "log2(Hazard Ratio)",
           title = "Lineage states and better survival",
           subtitle = paste0("Method: ", best_method)) +
      pub_theme(13) +
      theme(axis.text.x = element_text(angle = 15, hjust = 1))
    save_pub_gg(p, section, "s5_lineage_survival", width = 5.5, height = 4.5)
  }
} else {
  make_placeholder(section, "s5_lineage_survival", "Lineage survival",
                   "Run TCGA survival analysis first,\nthen rerun this script.", width = 5.5, height = 4.5)
}

cat("Section 5 complete.\n")


# Extracted from poster_section6_flot_remodelling.R
# FIGURE 3: MP expression paired boxplots for MP5, MP4, MP8, MP10 only
# ==============================================================================
cat("Generating MP expression paired plot...\n")
mp_file <- file.path(flot_dir, "Auto_pdo_flot_matched_mp_expression_changes.csv")
if (file.exists(mp_file)) {
  pdo_mp_desc_sel <- c("MP5" = "MYC proliferation", "MP4" = "Intestinal metaplasia",
                        "MP8" = "Columnar progenitor", "MP10" = "Inflammatory stress")
  mp <- read_csv(mp_file, show_col_types = FALSE) |>
    filter(MP %in% names(pdo_mp_desc_sel)) |>
    mutate(MP_label = factor(
      paste0(MP, "\n", pdo_mp_desc_sel[MP]),
      levels = paste0(c("MP5", "MP4", "MP8", "MP10"), "\n", pdo_mp_desc_sel[c("MP5", "MP4", "MP8", "MP10")])
    )) |>
    dplyr::select(Patient, Patient_label, MP, MP_label, Untreated, Treated) |>
    pivot_longer(c(Untreated, Treated), names_to = "condition", values_to = "score") |>
    mutate(condition = factor(condition, levels = c("Untreated", "Treated")))

  p <- ggplot(mp, aes(condition, score, group = Patient)) +
    geom_line(colour = "#64748B", linewidth = 0.55, alpha = 0.6) +
    geom_point(aes(fill = condition), shape = 21, size = 3, colour = "#111827", stroke = 0.3) +
    facet_wrap(~ MP_label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("Untreated" = "#F8FAFC", "Treated" = "#B2182B"), guide = "none") +
    labs(x = NULL, y = "Mean UCell score") +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
          strip.text = element_text(size = 9.5))
  save_pub_gg(p, section, "s6_mp_paired_expression", width = 9, height = 4)
}

cat("Section 6 complete.\n")


# Extracted from poster_section7_survival_targeting.R
# FIGURE 1: Survival association for stress-adaptive and classic proliferative
# ==============================================================================
cat("Generating survival figure...\n")
survival_csv <- file.path(SCREF_REF_OUTS_DIR, "task2_survival", "Auto_task2_survival_mp_state_cox_methods_splits.csv")

if (file.exists(survival_csv)) {
  cox <- read_csv(survival_csv, show_col_types = FALSE) |>
    filter(cohort == "EAC",
           mode == "noreg",
           split_method == "continuous")

  # Get state-level results for poor-prognosis states
  cox_states <- cox |>
    filter(feature_type == "State", feature %in% poor_states)

  # Also get MP-level results for the MPs belonging to these states
  poor_mps <- c("MP2", "MP13", "MP12")  # MYC, Hypoxic Inflam, Neuro-responsive
  cox_mps <- cox |>
    filter(feature_type == "MP", feature %in% poor_mps)

  # Pick best method
  combined <- bind_rows(cox_states, cox_mps)
  if (nrow(combined) > 0) {
    best_method <- "whole_tcga"

    # Bar plot for the selected features
    barplot_data <- combined |>
      filter(method == best_method) |>
      mutate(log2HR = log2(HR),
             neglog10 = -log10(P_value),
             sig = P_value < 0.05,
             label = case_when(
               feature %in% names(SCREF_MP_DESCRIPTIONS) ~
                 paste0(feature, " ", SCREF_MP_DESCRIPTIONS[feature]),
               TRUE ~ feature
             ),
             label = fct_reorder(label, log2HR))

    p <- ggplot(barplot_data, aes(x = log2HR, y = label, fill = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
      geom_col(width = 0.7, colour = "#111827", linewidth = 0.3) +
      geom_text(aes(label = sprintf("P=%.3f", P_value)), 
                hjust = ifelse(barplot_data$log2HR > 0, -0.1, 1.1),
                size = 3.5, colour = "#111827") +
      scale_fill_manual(values = c("FALSE" = "#CBD5E1", "TRUE" = "#DC2626"),
                          labels = c("n.s.", "P < 0.05"), name = "Significance") +
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = "log2(Hazard Ratio)", y = NULL) +
      pub_theme(13) +
      theme(legend.position = "bottom")
    save_pub_gg(p, section, "s7_survival_volcano", width = 7.5, height = 4.5)
    write_csv(barplot_data, file.path(out_dir, "tables", "s7_survival_volcano.csv"))
  }
} else {
  make_placeholder(section, "s7_survival_volcano", "Survival volcano",
                   "Run TCGA survival analysis first,\nthen rerun this script.", width = 7.5, height = 5.5)
}

# ==============================================================================

# Extracted from poster_section7_survival_targeting.R
# FIGURE 2: scAtlas and PDO L1000 signature reversal profiles
# ==============================================================================
cat("Generating L1000 reversal profiles...\n")
profile_files <- c(
  "scAtlas" = file.path(SCREF_REF_OUTS_DIR, "Auto_drug_reversal", "method_visuals",
                        "Auto_drug_reversal_l1000_signature_reversal_profiles.csv"),
  "PDOs" = "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_drug_reversal/method_visuals/Auto_drug_reversal_l1000_signature_reversal_profiles.csv"
)
profiles <- bind_rows(lapply(names(profile_files), function(dataset) {
  if (!file.exists(profile_files[[dataset]])) return(NULL)
  read_csv(profile_files[[dataset]], show_col_types = FALSE) |>
    mutate(dataset = dataset)
}))
if (nrow(profiles) > 0) {
  prof_plot <- profiles |>
    filter(state %in% poor_states) |>
    group_by(dataset, state, drug_label) |>
    mutate(best_rank = min(mean_l1000_rank, na.rm = TRUE)) |>
    ungroup() |>
    group_by(dataset, state) |>
    slice_min(best_rank, n = 7, with_ties = FALSE) |>
    ungroup() |>
    mutate(state = factor(state, levels = poor_states),
           dataset = factor(dataset, levels = c("scAtlas", "PDOs")),
           direction = factor(direction, levels = c("State-up genes", "State-down genes")),
           drug_label = fct_reorder(drug_label, mean_l1000_rank, .fun = min))
  p <- ggplot(prof_plot, aes(mean_l1000_rank, drug_label, colour = direction)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#CBD5E1") +
    geom_point(size = 2.6, alpha = 0.92) +
    geom_path(aes(group = interaction(dataset, state, drug_label)), colour = "#CBD5E1", linewidth = 0.35) +
    facet_grid(dataset ~ state, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = c("State-up genes" = "#B2182B", "State-down genes" = "#2166AC"),
                        name = "Signature") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
    labs(x = "L1000 reversal rank (lower is stronger)", y = NULL) +
    pub_theme(11) +
    theme(strip.text = element_text(size = 10.5, face = "bold"),
          axis.text.y = element_text(size = 8.5),
          legend.position = "bottom")
  save_pub_gg(p, section, "s7_l1000_reversal_profiles", width = 9, height = 6.5)
  pub_copy_to_assets(file.path(out_dir, "figures", "s7_l1000_reversal_profiles.png"),
                     "s7_inhibitor_overlap.png", section)
  write_csv(prof_plot, file.path(out_dir, "tables", "s7_l1000_reversal_profiles.csv"))
}

cat("Section 7 complete.\n")
