####################
# delete_refined_mp_abundance.R
#
# Temporary script: compare per-sample refined MP abundance between
# centred and uncentred GeneNMF workflows.
#
# Input:
#   ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   ref_outs/by_samples/
#
# Output:
#   ref_outs/Metaprogrammes_Results/centred_comparison/figures/Auto_refined_mp_sample_abundance.pdf
#
# Run:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript /rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG/delete_refined_mp_abundance.R
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(grid)

project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
ref_dir <- file.path(project_dir, "ref_outs")
setwd(ref_dir)

outdir <- file.path(ref_dir, "Metaprogrammes_Results", "centred_comparison", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load UCell scores ────────────────────────────────────────────────────────
unc_ucell <- as.data.frame(readRDS("Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_ucell_scores.rds"), check.names = FALSE)
cent_ucell <- as.data.frame(readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds"), check.names = FALSE)

common_cells <- intersect(rownames(unc_ucell), rownames(cent_ucell))
unc_ucell <- unc_ucell[common_cells, , drop = FALSE]
cent_ucell <- cent_ucell[common_cells, , drop = FALSE]
cat("Common cells:", length(common_cells), "\n")

# ── Infer sample identity ────────────────────────────────────────────────────
infer_samples <- function(cells) {
  sample_dirs <- list.dirs("by_samples", full.names = FALSE, recursive = FALSE)
  sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]
  sample_dirs <- sample_dirs[order(nchar(sample_dirs), decreasing = TRUE)]
  out <- rep(NA_character_, length(cells))
  for (sample in sample_dirs) {
    idx <- is.na(out) & startsWith(cells, paste0(sample, "_"))
    out[idx] <- sample
  }
  if (anyNA(out)) {
    warning("Could not infer sample for ", sum(is.na(out)), " cells")
    out[is.na(out)] <- cells[is.na(out)]
  }
  out
}

sample_vec <- infer_samples(common_cells)
names(sample_vec) <- common_cells
study_vec <- sub("^([A-Za-z]+_[0-9]{4}).*", "\\1", sample_vec)
cat("Unique samples:", length(unique(sample_vec)), "\n")

# ── Description labels ───────────────────────────────────────────────────────

# Uncentred refined MP labels (from centred_method_comparison_figures.R)
submp_desc_map <- c(
  "MP7+"  = "Fanconi/HR repair progenitor",
  "MP7h"  = "Replication-stress dormant epi.",
  "MP7j"  = "DNA damage response",
  "MP7r"  = "Stem-like glandular duct progenitor",
  "MP7v"  = "Mucous secretory progenitor",
  "MP10+" = "Metabolic columnar epi.",
  "MP10e" = "Inflam. mucous-secr. columnar epi.",
  "MP8+"  = "Intestinal metaplasia",
  "MP8b"  = "Quiescent glandular-metabolic progenitor",
  "MP8c"  = "NF-kB inflam. cycling glandular progenitor",
  "MP8e"  = "Cycling intestinal-columnar progenitor",
  "MP12a" = "Enteroendocrine-primed progenitor",
  "MP12b" = "Enteroendocrine differentiation",
  "MP12c" = "Cycling glandular-intestinal progenitor",
  "MP2+"  = "MYC proliferation",
  "MP2v"  = "EMT-V cycling invasive progenitor",
  "MP15a" = "T/NK-like epi. immune mimicry",
  "MP15b" = "Type I IFN-activated EMT-primed epi.",
  "MP15c" = "Type II IFN / NF-kB peak inflam. epi.",
  "MP5+"  = "Epithelial IFN Resp.",
  "MP16+" = "Secretory Diff. (Gastric)",
  "MP1"   = "G2M Cell Cycle",
  "MP9"   = "G1S Cell Cycle",
  "MP13"  = "Hypoxic Inflam. Epi.",
  "MP14"  = "Hypoxia Adapted Epi.",
  "MP17"  = "Basal-like Transition",
  "MP18"  = "Secretory Diff. (Intest.)"
)

# Centred refined MP labels
cent_desc_map <- c(
  "MP10+" = "Intestinal metaplasia",
  "MP9+" = "Metabolic columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP14" = "Squamoid/basal transition",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP18b" = "Mucous-secretory differentiation",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP2x" = "Wnt-active glandular stem/progenitor",
  "MP18a" = "MP18a",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP11c" = "MP11c",
  "MP15" = "T/NK-like cancer-cell immune mimicry",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP2+" = "MYC driven biosynthesis"
)

make_display_label <- function(mp_name, desc_map) {
  desc <- desc_map[mp_name]
  if (is.na(desc) || desc == mp_name) return(mp_name)
  # If adapted_label already includes the MP prefix, use it directly
  if (startsWith(desc, mp_name)) return(desc)
  paste0(mp_name, ": ", desc)
}

unc_labels <- setNames(
  vapply(colnames(unc_ucell), make_display_label, character(1), desc_map = submp_desc_map),
  colnames(unc_ucell)
)
cent_labels <- setNames(
  vapply(colnames(cent_ucell), make_display_label, character(1), desc_map = cent_desc_map),
  colnames(cent_ucell)
)

# ── Top-MP assignment and abundance computation ─────────────────────────────

compute_topmp_abundance <- function(ucell_mat, sample_vec, mp_labels) {
  # Assign each cell to its top-scoring MP
  topmp <- colnames(ucell_mat)[max.col(as.matrix(ucell_mat), ties.method = "first")]
  names(topmp) <- rownames(ucell_mat)
  
  # Build abundance per sample
  df <- data.frame(
    cell = names(topmp),
    mp = topmp,
    sample = sample_vec[names(topmp)],
    stringsAsFactors = FALSE
  )
  
  all_mps <- colnames(ucell_mat)
  all_samples <- sort(unique(df$sample))
  
  counts <- df %>%
    count(sample, mp, name = "n")
  
  full <- expand_grid(sample = all_samples, mp = all_mps) %>%
    left_join(counts, by = c("sample", "mp")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    group_by(sample) %>%
    mutate(
      total = sum(n),
      pct = if_else(total > 0, 100 * n / total, 0)
    ) %>%
    ungroup()
  
  # Add display labels
  full$mp_label <- mp_labels[full$mp]
  full$mp_label[is.na(full$mp_label)] <- full$mp[is.na(full$mp_label)]
  
  # Totals
  totals <- full %>% distinct(sample, total)
  
  list(abundance = full, totals = totals, topmp = topmp, all_mps = all_mps)
}

unc_abund <- compute_topmp_abundance(unc_ucell, sample_vec, unc_labels)
cent_abund <- compute_topmp_abundance(cent_ucell, sample_vec, cent_labels)

# ── Normalised versions ──────────────────────────────────────────────────────
z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

unc_ucell_norm <- z_normalise(unc_ucell, sample_vec, study_vec)
cent_ucell_norm <- z_normalise(cent_ucell, sample_vec, study_vec)

unc_abund_norm <- compute_topmp_abundance(unc_ucell_norm, sample_vec, unc_labels)
cent_abund_norm <- compute_topmp_abundance(cent_ucell_norm, sample_vec, cent_labels)

# ── Expression (Mean UCell) ──────────────────────────────────────────────────
compute_mean_expression <- function(ucell_mat, sample_vec, mp_labels) {
  df <- ucell_mat
  df$sample <- sample_vec[rownames(df)]
  
  mean_expr <- df %>%
    group_by(sample) %>%
    summarise(across(all_of(colnames(ucell_mat)), mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(cols = -sample, names_to = "mp", values_to = "pct")
    
  totals <- data.frame(
    sample = unique(sample_vec),
    total = as.integer(table(sample_vec)[unique(sample_vec)]),
    stringsAsFactors = FALSE
  )
    
  list(abundance = mean_expr, totals = totals)
}

unc_expr <- compute_mean_expression(unc_ucell, sample_vec, unc_labels)
cent_expr <- compute_mean_expression(cent_ucell, sample_vec, cent_labels)

# ── Sort order: by diversity (geometric mean across MPs) ─────────────────────
compute_diversity_order <- function(abund_data) {
  abund_data$abundance %>%
    group_by(sample) %>%
    summarise(
      geo_mean = exp(mean(log(n + 1))),
      total = first(total),
      .groups = "drop"
    ) %>%
    arrange(desc(geo_mean), sample) %>%
    pull(sample)
}

diversity_order <- compute_diversity_order(unc_abund)

# ── Colour palettes ─────────────────────────────────────────────────────────
make_mp_palette <- function(mp_names) {
  n <- length(mp_names)
  if (n <= 12) {
    cols <- brewer.pal(max(n, 3), "Paired")
  } else if (n <= 20) {
    cols <- c(brewer.pal(12, "Paired"), brewer.pal(max(n - 12, 3), "Set3"))
  } else {
    cols <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set3"), viridis(n - 20, option = "turbo"))
  }
  setNames(cols[seq_len(n)], mp_names)
}

unc_mp_cols <- make_mp_palette(unc_abund$all_mps)
cent_mp_cols <- make_mp_palette(cent_abund$all_mps)

# ── Parent MP extraction (for grouping) ──────────────────────────────────────
parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]+$", "", x))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Tile heatmap — MP proportion per sample (both methods side-by-side)
#
# A heatmap is far better than stacked bars for 20-27 MPs × 117 samples.
# Rows = MPs, columns = samples (sorted by diversity), fill = proportion.
# ═══════════════════════════════════════════════════════════════════════════════

make_abundance_heatmap <- function(abund_data, sample_order, mp_labels, title, is_expression = FALSE, row_order_param = NULL) {
  # Build proportion matrix: rows = MPs, cols = samples
  mat <- abund_data$abundance %>%
    select(mp, sample, pct) %>%
    pivot_wider(names_from = sample, values_from = pct, values_fill = 0) %>%
    tibble::column_to_rownames("mp") %>%
    as.matrix()
  
  # Reorder
  mat <- mat[, intersect(sample_order, colnames(mat)), drop = FALSE]
  
  # Row labels
  row_display <- mp_labels[rownames(mat)]
  row_display[is.na(row_display)] <- rownames(mat)[is.na(row_display)]
  
  # Sort rows by total percentage sum across samples or use provided order
  if (!is.null(row_order_param)) {
    ordered_rows <- intersect(row_order_param, rownames(mat))
    mat <- mat[ordered_rows, , drop = FALSE]
  } else {
    row_sums <- rowSums(mat)
    mat <- mat[order(row_sums, decreasing = FALSE), , drop = FALSE]
  }
  row_display <- row_display[rownames(mat)]
  
  # Parent MP annotation for row split
  parent_vec <- parent_id(rownames(mat))
  
  # Study annotation for columns
  study_names <- sub("^([A-Za-z]+_[0-9]{4}).*", "\\1", colnames(mat))
  n_studies <- length(unique(study_names))
  qual_cols <- c(brewer.pal(max(3, min(n_studies, 9)), "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
  
  study_cols <- setNames(
    qual_cols[1:n_studies],
    sort(unique(study_names))
  )
  
  # Sample size annotation
  totals <- abund_data$totals
  sample_sizes <- setNames(totals$total, totals$sample)
  
  top_ha <- HeatmapAnnotation(
    `Cell Count` = anno_barplot(
      sample_sizes[colnames(mat)],
      gp = gpar(fill = "grey40", col = NA),
      height = unit(15, "mm"),
      border = FALSE
    ),
    Study = study_names,
    col = list(Study = study_cols),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = FALSE,
    simple_anno_size = unit(3, "mm")
  )
  
  if (is_expression) {
    col_fun <- colorRamp2(
      c(0, 0.05, 0.1, 0.2, 0.4),
      c("#F7F7F7", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D")
    )
    legend_name <- "Mean\nUCell"
  } else {
    col_fun <- colorRamp2(
      c(0, 5, 15, 40, 80),
      c("#F7F7F7", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D")
    )
    legend_name <- "% of cells"
  }
  
  rownames(mat) <- row_display
  
  Heatmap(
    mat,
    name = legend_name,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = study_names,
    column_gap = unit(0.2, "mm"),
    cluster_column_slices = FALSE,
    show_column_names = TRUE,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 7),
    row_names_max_width = unit(6, "cm"),
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    top_annotation = top_ha,
    rect_gp = gpar(col = "white", lwd = 0.3),
    width = unit(18, "inch"),
    height = unit(0.35 * nrow(mat), "inch"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(35, "mm")
    ),
    use_raster = TRUE,
    raster_quality = 3
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Faceted dot-plot — Mean UCell score per sample per MP
#
# Dot size = cell count, colour = mean UCell score.
# Faceted by parent MP group for readability.
# ═══════════════════════════════════════════════════════════════════════════════

make_dotplot <- function(ucell_mat, sample_vec, mp_labels, sample_order, title) {
  # Compute per-sample per-MP mean UCell and cell count
  long <- as.data.frame(ucell_mat, check.names = FALSE) %>%
    mutate(sample = sample_vec[rownames(ucell_mat)]) %>%
    pivot_longer(cols = -sample, names_to = "mp", values_to = "score") %>%
    group_by(sample, mp) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    )
  
  # Add display labels and parent groups
  long$mp_label <- mp_labels[long$mp]
  long$mp_label[is.na(long$mp_label)] <- long$mp[is.na(long$mp_label)]
  long$parent <- parent_id(long$mp)
  
  # Factor for ordering
  long$sample <- factor(long$sample, levels = sample_order)
  
  # Order MPs: by parent, then by suffix
  mp_order <- sort(unique(long$mp))
  mp_num <- as.integer(sub("^MP", "", sub("[a-z\\+]*$", "", mp_order)))
  mp_order <- mp_order[order(mp_num, mp_order)]
  label_order <- mp_labels[mp_order]
  label_order[is.na(label_order)] <- mp_order[is.na(label_order)]
  long$mp_label <- factor(long$mp_label, levels = label_order)
  
  ggplot(long, aes(x = sample, y = mp_label, size = mean_score, colour = mean_score)) +
    geom_point() +
    scale_colour_viridis_c(option = "inferno", name = "Mean\nUCell") +
    scale_size_continuous(range = c(0.2, 3.5), name = "Mean\nUCell") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_line(colour = "grey95"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "right"
    )
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Per-parent-MP stacked bar chart
#
# Split the 20-27 MPs by their parent MP (e.g. MP7, MP8, MP12, etc.) and
# show stacked bars of only the sub-MPs within each parent. This keeps each
# panel readable (2-5 colours max per panel).
# ═══════════════════════════════════════════════════════════════════════════════

make_parent_faceted_bars <- function(abund_data, sample_order, mp_labels, mp_cols, title) {
  df <- abund_data$abundance %>%
    mutate(
      parent = parent_id(mp),
      mp_label = mp_labels[mp]
    )
  df$mp_label[is.na(df$mp_label)] <- df$mp[is.na(df$mp_label)]
  df$sample <- factor(df$sample, levels = sample_order)
  
  # Only keep parents with >1 sub-MP for the faceted view
  parent_counts <- df %>% distinct(parent, mp) %>% count(parent, name = "n_sub")
  multi_parents <- parent_counts %>% filter(n_sub > 1) %>% pull(parent)
  
  if (length(multi_parents) == 0) {
    message("No parent MPs with multiple sub-MPs for faceted bars.")
    return(NULL)
  }
  
  df_multi <- df %>% filter(parent %in% multi_parents)
  
  # Re-normalise within each parent per sample
  df_multi <- df_multi %>%
    group_by(sample, parent) %>%
    mutate(
      parent_total = sum(n),
      pct_within_parent = if_else(parent_total > 0, 100 * n / parent_total, 0)
    ) %>%
    ungroup()
  
  # Order parents by their numeric ID
  parent_num <- as.integer(sub("^MP", "", multi_parents))
  multi_parents <- multi_parents[order(parent_num)]
  df_multi$parent <- factor(df_multi$parent, levels = multi_parents)
  
  # Re-key colour palette by display labels (mp_label) instead of raw MP names
  label_cols <- setNames(mp_cols, mp_labels[names(mp_cols)])
  # Ensure all labels in data have a colour
  all_labels <- unique(df_multi$mp_label)
  missing <- setdiff(all_labels, names(label_cols))
  if (length(missing) > 0) {
    extra <- setNames(scales::hue_pal()(length(missing)), missing)
    label_cols <- c(label_cols, extra)
  }
  
  ggplot(df_multi, aes(x = sample, y = pct_within_parent, fill = mp_label)) +
    geom_col(width = 0.85) +
    facet_wrap(~ parent, ncol = 1, scales = "free_y",
               strip.position = "left",
               labeller = labeller(parent = function(x) paste("Parent", x))) +
    scale_fill_manual(values = label_cols, name = "Refined MP") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = title, x = NULL, y = "% within parent MP") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      panel.spacing = unit(0.3, "lines"),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 13, face = "bold")
    ) +
    guides(fill = guide_legend(ncol = 4))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Side-by-side stacked bars (top-MP assignment)
#
# Like state_mp_sample_abundance.R but with one page per method.
# For 20-27 MPs the legend is large, so we keep it at the bottom in multiple cols.
# ═══════════════════════════════════════════════════════════════════════════════

make_stacked_bar <- function(abund_data, sample_order, mp_labels, mp_cols, title) {
  df <- abund_data$abundance
  df$mp_label <- mp_labels[df$mp]
  df$mp_label[is.na(df$mp_label)] <- df$mp[is.na(df$mp_label)]
  df$sample <- factor(df$sample, levels = sample_order)
  
  # Stack order: reverse so top of legend matches top of bar
  label_order <- mp_labels[abund_data$all_mps]
  label_order[is.na(label_order)] <- abund_data$all_mps[is.na(label_order)]
  df$mp_label <- factor(df$mp_label, levels = rev(label_order))
  
  # named colour vector by label
  col_vec <- mp_cols
  names(col_vec) <- mp_labels[names(col_vec)]
  
  totals_df <- abund_data$totals %>%
    filter(sample %in% sample_order) %>%
    mutate(sample = factor(sample, levels = sample_order))
  
  scale_factor <- max(totals_df$total, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1
  
  ggplot(df, aes(x = sample, y = pct, fill = mp_label)) +
    geom_col(width = 0.75) +
    geom_point(
      data = totals_df,
      aes(x = sample, y = total / scale_factor),
      color = "black", size = 1.2, shape = 18, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = col_vec, breaks = rev(label_order), drop = FALSE) +
    scale_y_continuous(
      name = "Proportion (%)",
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cell Count (N)", labels = comma)
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(title = title, x = NULL, fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
      legend.position = "bottom",
      legend.text = element_text(size = 5.5),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 13, face = "bold")
    ) +
    guides(fill = guide_legend(ncol = 4, reverse = FALSE))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Shannon diversity per sample comparison (centred vs uncentred)
# ═══════════════════════════════════════════════════════════════════════════════

compute_shannon <- function(abund_data) {
  abund_data$abundance %>%
    filter(total > 0) %>%
    mutate(p = n / total) %>%
    filter(p > 0) %>%
    group_by(sample) %>%
    summarise(
      shannon = -sum(p * log(p)),
      total = first(total),
      .groups = "drop"
    )
}

unc_shannon <- compute_shannon(unc_abund) %>% mutate(method = "Uncentred")
cent_shannon <- compute_shannon(cent_abund) %>% mutate(method = "Centred")

shannon_df <- bind_rows(unc_shannon, cent_shannon)

p_shannon <- ggplot(shannon_df, aes(x = method, y = shannon, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = c("Uncentred" = "#4DAF4A", "Centred" = "#377EB8")) +
  labs(
    title = "Shannon Diversity of Refined MP Proportions per Sample",
    y = "Shannon Diversity (H')",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

shannon_scatter_df <- inner_join(
  unc_shannon %>% select(sample, shannon_unc = shannon),
  cent_shannon %>% select(sample, shannon_cent = shannon),
  by = "sample"
)

p_shannon_scatter <- ggplot(shannon_scatter_df, aes(x = shannon_unc, y = shannon_cent)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(
    title = "Shannon Diversity: Uncentred vs Centred",
    x = "Uncentred Shannon H'",
    y = "Centred Shannon H'"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Cross-method agreement — per-cell top-MP concordance
# ═══════════════════════════════════════════════════════════════════════════════

# For matched MPs (via Jaccard from label transfer), check if cells get
# assigned to equivalent MPs in both methods
cross_jaccard_csv <- file.path(ref_dir, "Metaprogrammes_Results", "centred_comparison", "tables", "Auto_refined_centred_label_transfer.csv")
if (file.exists(cross_jaccard_csv)) {
  lt <- read.csv(cross_jaccard_csv, stringsAsFactors = FALSE)
  
  # Build mapping: centred -> uncentred
  cent_to_unc <- setNames(lt$matched_uncentred_mp, lt$centred_mp)
  cent_to_unc_jac <- setNames(lt$best_jaccard, lt$centred_mp)
  
  # Per sample: fraction of cells where top centred MP maps to same parent as top uncentred MP
  unc_top <- colnames(unc_ucell)[max.col(as.matrix(unc_ucell), ties.method = "first")]
  names(unc_top) <- rownames(unc_ucell)
  cent_top <- colnames(cent_ucell)[max.col(as.matrix(cent_ucell), ties.method = "first")]
  names(cent_top) <- rownames(cent_ucell)
  
  # Match by parent: uncentred parent vs centred's best uncentred match's parent
  unc_parent <- parent_id(unc_top)
  cent_match_parent <- parent_id(cent_to_unc[cent_top])
  
  agree_df <- data.frame(
    sample = sample_vec,
    agree = unc_parent == cent_match_parent,
    stringsAsFactors = FALSE
  ) %>%
    group_by(sample) %>%
    summarise(
      pct_agree = 100 * mean(agree, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(pct_agree))
  
  agree_df$sample <- factor(agree_df$sample, levels = agree_df$sample)
  
  p_agree <- ggplot(agree_df, aes(x = sample, y = pct_agree)) +
    geom_col(aes(fill = pct_agree), width = 0.8) +
    scale_fill_viridis_c(option = "viridis", name = "% Agree") +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "red", linewidth = 0.5) +
    labs(
      title = "Per-Sample Cross-Method Agreement (Parent-MP Level)",
      subtitle = "Fraction of cells where top centred refined MP maps to same parent as top uncentred refined MP",
      x = NULL,
      y = "% Cells Agreeing"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9)
    )
} else {
  p_agree <- NULL
}

# ═══════════════════════════════════════════════════════════════════════════════
# Render PDF
# ═══════════════════════════════════════════════════════════════════════════════

n_samples <- length(unique(sample_vec))
pdf_w <- max(22, 0.2 * n_samples)
pdf_h_bar <- 10

pdf_file <- file.path(outdir, "Auto_refined_mp_sample_abundance.pdf")
pdf(pdf_file, width = pdf_w, height = 14, useDingbats = FALSE, onefile = TRUE)

# Compute global row orders based on raw abundance sums (to keep exactly the same across the 3 plots)
compute_global_order <- function(abund_data) {
  mat <- abund_data$abundance %>%
    select(mp, sample, pct) %>%
    pivot_wider(names_from = sample, values_from = pct, values_fill = 0) %>%
    tibble::column_to_rownames("mp") %>%
    as.matrix()
  
  row_sums <- rowSums(mat)
  rownames(mat)[order(row_sums, decreasing = TRUE)]
}

unc_global_order <- compute_global_order(unc_abund)
cent_global_order <- compute_global_order(cent_abund)

# Page 1: Tile heatmap — uncentred (raw abundance)
cat("Rendering page 1: Uncentred tile heatmap (raw abundance)\n")
draw(
  make_abundance_heatmap(unc_abund, diversity_order, unc_labels,
                         "Uncentred Refined MP Abundance (Raw UCell)", row_order_param = unc_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Page 2: Tile heatmap — centred (raw abundance)
cat("Rendering page 2: Centred tile heatmap (raw abundance)\n")
draw(
  make_abundance_heatmap(cent_abund, diversity_order, cent_labels,
                         "Centred Refined MP Abundance (Raw UCell)", row_order_param = cent_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Page 3: Tile heatmap — uncentred (normalised abundance)
cat("Rendering page 3: Uncentred tile heatmap (normalised abundance)\n")
draw(
  make_abundance_heatmap(unc_abund_norm, diversity_order, unc_labels,
                         "Uncentred Refined MP Abundance (Normalised UCell)", row_order_param = unc_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Page 4: Tile heatmap — centred (normalised abundance)
cat("Rendering page 4: Centred tile heatmap (normalised abundance)\n")
draw(
  make_abundance_heatmap(cent_abund_norm, diversity_order, cent_labels,
                         "Centred Refined MP Abundance (Normalised UCell)", row_order_param = cent_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Page 5: Tile heatmap — uncentred (expression)
cat("Rendering page 5: Uncentred tile heatmap (expression)\n")
draw(
  make_abundance_heatmap(unc_expr, diversity_order, unc_labels,
                         "Uncentred Refined MP Expression (Mean UCell)", is_expression = TRUE, row_order_param = unc_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# Page 6: Tile heatmap — centred (expression)
cat("Rendering page 6: Centred tile heatmap (expression)\n")
draw(
  make_abundance_heatmap(cent_expr, diversity_order, cent_labels,
                         "Centred Refined MP Expression (Mean UCell)", is_expression = TRUE, row_order_param = cent_global_order),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()
cat("Saved:", pdf_file, "\n")
cat("Done.\n")
