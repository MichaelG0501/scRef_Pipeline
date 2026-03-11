####################
# clinical_mp_ucell_plots.R
# Clinical association plots using raw per-cell MP UCell scores.
#
# Visualisation approach:
#   A) Dot-plot heatmap:  MPs (y) × clinical levels (x)
#      - Dot colour = z-scored MP activity (blue → red)
#      - Dot size  = absolute z-score (emphasises strong signals)
#      - MPs grouped by biological function (state groups)
#   B) Lollipop difference plot (for 2-level variables only):
#      - Shows directional change between groups
#
# Aggregation:
#   1) per-sample mean MP score
#   2) mean of sample-means within each clinical group
#   3) z-score across clinical levels per MP (row-wise) for colour mapping
#
# Input:
#   ref_outs/meta_full_epi.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx (sheet 3)
#
# Output:
#   ref_outs/Auto_clinical_assoc_mp_ucell_combined.pdf
#   ref_outs/Auto_clinical_assoc_mp_ucell_per_study.pdf
#   new_update/summaries/Auto_clinical_assoc_mp_ucell_summary.csv
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(scales)
library(patchwork)
library(ggnewscale)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# 1) Load data
####################
meta_full_epi <- readRDS("meta_full_epi.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
  sheet = 3,
  skip = 1
) %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

####################
# 2) MP order and annotations (matched to Auto_states_topmp_v2.R)
####################
mp_descriptions <- c(
  "MP1" = "G2M_cycle",
  "MP2" = "MYC_prolif",
  "MP5" = "IFN_response",
  "MP7" = "S_cycle",
  "MP8" = "Intestinal_diff",
  "MP9" = "G1S_cycle",
  "MP10" = "Columnar_diff",
  "MP12" = "Neuro_epithelial",
  "MP13" = "Partial_EMT",
  "MP14" = "Hypoxia_epithelial",
  "MP15" = "T_NK_infiltration",
  "MP16" = "Secretory_diff",
  "MP17" = "Squamous_transition",
  "MP18" = "Adaptive_secretory"
)

# State groups from Auto_states_topmp_v2.R
state_groups <- list(
  "Cell Cycle"              = c("MP1", "MP7", "MP9"),
  "Classic Proliferative"   = c("MP2"),
  "Barretts Metaplasia"   = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "EMT-related"             = c("MP13", "MP12"),
  "Intestinal Metaplasia"  = c("MP18", "MP16"),
  "Immune Infiltrated"      = c("MP15")
)

# Silhouette filtering
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)

# Derive tree-order
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

mp_cols <- intersect(mp_tree_order_names, colnames(ucell_scores))
mp_labels <- mp_descriptions
mp_labels[setdiff(mp_cols, names(mp_labels))] <- setdiff(mp_cols, names(mp_labels))

# Build display label: "MP# Description", ordered by state group
mp_display <- setNames(
  paste0(names(mp_labels[mp_cols]), " ", mp_labels[mp_cols]),
  mp_cols
)

# Build MP → group mapping for y-axis grouping
mp_to_group <- setNames(rep("Other", length(mp_cols)), mp_cols)
for (grp in names(state_groups)) {
  for (mp in intersect(state_groups[[grp]], mp_cols)) {
    mp_to_group[mp] <- grp
  }
}

# Order MPs: within each group, keep tree order; groups ordered as defined
group_order <- names(state_groups)
mp_ordered <- c()
for (grp in group_order) {
  grp_mps <- mp_cols[mp_cols %in% state_groups[[grp]]]
  mp_ordered <- c(mp_ordered, grp_mps)
}
remaining <- setdiff(mp_cols, mp_ordered)
mp_ordered <- c(mp_ordered, remaining)

# Group colour palette (for strip/annotation)
group_palette <- c(
  "Cell Cycle"            = "#FFD700",
  "Classic Proliferative" = "#E41A1C",
  "Intestinal Metaplasia" = "#4DAF4A",
  "EMT-related"           = "#984EA3",
  "Secretory"             = "#FF7F00",
  "Immune Infiltrated"    = "#377EB8",
  "Other"                 = "grey60"
)

####################
# 3) Cell alignment and metadata
####################
cell_ids <- Reduce(intersect, list(rownames(ucell_scores), rownames(meta_full_epi)))

study_vec <- if ("study" %in% colnames(meta_full_epi)) {
  as.character(meta_full_epi[cell_ids, "study"])
} else {
  orig_tmp <- as.character(meta_full_epi[cell_ids, "orig.ident"])
  sub("^(([^_]+)_[^_]+).*$", "\\1", orig_tmp)
}

cell_meta <- data.frame(
  cell = cell_ids,
  orig.ident = as.character(meta_full_epi[cell_ids, "orig.ident"]),
  study = study_vec,
  stringsAsFactors = FALSE
) %>%
  left_join(clinical_sheet, by = "orig.ident") %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age)),
    Age_Group = ifelse(!is.na(Age) & Age > 60, ">60", "<=60"),
    Gender = ifelse(Gender %in% c("F", "female", "Female"), "Female", "Male"),
    Treatment = factor(Treatment, levels = c("Tx-naïve", "Post")),
    `Clinical response` = factor(
      ifelse(`Clinical response` == "R", "Responder", "Nonresponder"),
      levels = c("Responder", "Nonresponder")
    )
  )

cell_long <- as.data.frame(ucell_scores[cell_ids, mp_cols, drop = FALSE]) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(mp_cols), names_to = "MP", values_to = "score") %>%
  left_join(cell_meta, by = "cell")

####################
# 4) Core aggregation helper
####################
compute_mp_group_means <- function(data, group_var, filter_expr = NULL, studyname = NULL) {
  if (!is.null(studyname)) data <- data %>% filter(study == studyname)
  if (!is.null(filter_expr)) data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  
  data <- data %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(MP), !is.na(score), !is.na(study))
  
  if (nrow(data) == 0) return(list(plot_data = NULL, label_data = NULL, total_samples = 0L))
  
  sample_mean <- data %>%
    group_by(orig.ident, .data[[group_var]], MP) %>%
    summarise(sample_mean = mean(score, na.rm = TRUE), .groups = "drop")
  
  group_mean <- sample_mean %>%
    group_by(.data[[group_var]], MP) %>%
    summarise(group_mean = mean(sample_mean, na.rm = TRUE), .groups = "drop")
  
  total_samples <- data %>% distinct(orig.ident) %>% nrow()
  n_sample_group <- data %>%
    distinct(orig.ident, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_samples = n_distinct(orig.ident), .groups = "drop")
  
  n_cell_group <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_cells = n_distinct(cell), .groups = "drop")
  
  label_data <- n_sample_group %>%
    left_join(n_cell_group, by = group_var) %>%
    mutate(
      label = paste0(.data[[group_var]], "\n(n=", n_samples, " samples; ", comma(n_cells), " cells)")
    )
  
  list(plot_data = group_mean, label_data = label_data, total_samples = total_samples)
}

####################
# 5) Dot-plot heatmap (main visualisation)
#    Uses facet_grid on MP biological groups for clean y-axis annotation.
####################
plot_dot_heatmap <- function(data, group_var, title, filter_expr = NULL, studyname = NULL) {
  prep <- compute_mp_group_means(data, group_var, filter_expr, studyname)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  
  plot_df <- prep$plot_data
  
  # Z-score across clinical levels per MP (row-wise standardisation)
  plot_df <- plot_df %>%
    group_by(MP) %>%
    mutate(
      mp_mean = mean(group_mean, na.rm = TRUE),
      mp_sd   = sd(group_mean, na.rm = TRUE),
      z_score = ifelse(mp_sd > 0, (group_mean - mp_mean) / mp_sd, 0)
    ) %>%
    ungroup()
  
  # Clamp z-scores for colour scale
  z_lim <- max(abs(plot_df$z_score), na.rm = TRUE)
  z_lim <- max(z_lim, 0.5)
  
  # Display labels and grouping
  plot_df$mp_display <- factor(mp_display[plot_df$MP], levels = rev(mp_display[mp_ordered]))
  plot_df$mp_group <- factor(mp_to_group[plot_df$MP], levels = names(state_groups))
  
  # X-axis: clinical levels with sample info
  level_labels <- prep$label_data$label
  names(level_labels) <- as.character(prep$label_data[[group_var]])
  plot_df$x_label <- level_labels[as.character(plot_df[[group_var]])]
  plot_df$x_label <- factor(plot_df$x_label, levels = unique(level_labels))
  
  n_x <- length(unique(plot_df$x_label))
  
  # MP name within group (short label for y-axis inside each facet)
  plot_df$mp_short <- sub("^MP\\d+\\s+", "", as.character(plot_df$mp_display))
  # Preserve ordering within facets
  mp_short_levels <- sub("^MP\\d+\\s+", "", rev(mp_display[mp_ordered]))
  plot_df$mp_short <- factor(plot_df$mp_short, levels = unique(mp_short_levels))
  
  p <- ggplot(plot_df, aes(x = x_label, y = mp_short)) +
    # Dots (fixed size for clarity)
    geom_point(
      aes(fill = z_score),
      size = 7,
      shape = 21, colour = "grey30", stroke = 0.4
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-z_lim, z_lim),
      name = "Relative\nactivity\n(Z-score)"
    ) +
    # Raw score labels next to dots for small N
    {if (n_x <= 5) {
      geom_text(
        aes(label = sprintf("%.3f", group_mean)),
        size = 2.0, colour = "grey30", fontface = "plain",
        nudge_y = -0.3
      )
    }} +
    # Facet by MP group — creates labelled strips on the right
    facet_grid(
      mp_group ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = paste0("n=", prep$total_samples, " samples")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, margin = margin(b = 0)),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey45", margin = margin(b = 8)),
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 10, lineheight = 0.9),
      axis.text.y = element_text(size = 9.5),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0.3, "lines"),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 0, face = "bold", size = 8, hjust = 1,
        margin = margin(r = 0)
      ),
      strip.background = element_rect(fill = "grey97", colour = NA),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 15, 10, 10)
    )
  
  p
}


####################
# 6) Lollipop difference plot (for binary variables)
####################
plot_lollipop_diff <- function(data, group_var, title, filter_expr = NULL) {
  prep <- compute_mp_group_means(data, group_var, filter_expr)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  
  n_levels <- length(unique(prep$plot_data[[group_var]]))
  if (n_levels != 2) return(NULL)
  
  levels_vec <- sort(unique(as.character(prep$plot_data[[group_var]])))
  ref_level <- levels_vec[1]
  alt_level <- levels_vec[2]
  
  wide_df <- prep$plot_data %>%
    mutate(level = as.character(.data[[group_var]])) %>%
    dplyr::select(MP, level, group_mean) %>%
    pivot_wider(names_from = level, values_from = group_mean)
  
  wide_df$diff <- wide_df[[alt_level]] - wide_df[[ref_level]]
  wide_df$direction <- ifelse(wide_df$diff > 0, paste0("Higher in ", alt_level), paste0("Higher in ", ref_level))
  wide_df$mp_display <- factor(mp_display[wide_df$MP], levels = rev(mp_display[mp_ordered]))
  
  # Get sample counts for subtitle
  lbl_df <- prep$label_data
  
  ggplot(wide_df, aes(x = diff, y = mp_display, colour = direction)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = diff, y = mp_display, yend = mp_display), linewidth = 0.8) +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("%+.4f", diff)), 
              hjust = ifelse(wide_df$diff >= 0, -0.3, 1.3), size = 2.8, show.legend = FALSE) +
    scale_colour_manual(
      values = setNames(c("#B2182B", "#2166AC"), 
                        c(paste0("Higher in ", alt_level), paste0("Higher in ", ref_level))),
      name = "Direction"
    ) +
    labs(
      x = paste0("Δ Mean UCell score (", alt_level, " − ", ref_level, ")"),
      y = NULL,
      title = paste0(title, " — Difference"),
      subtitle = paste0("n=", prep$total_samples, " samples")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      axis.text.y = element_text(size = 9, face = "italic"),
      panel.grid.major.y = element_line(colour = "grey95"),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}


####################
# 7) Faceted per-study dot plot (smaller version)
####################
plot_dot_faceted <- function(data, group_var, title, filter_expr = NULL) {
  if (!is.null(filter_expr)) data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(study), !is.na(MP), !is.na(score))
  if (nrow(data) == 0) return(NULL)
  
  total_samples <- n_distinct(data$orig.ident)
  
  sample_mean <- data %>%
    group_by(study, orig.ident, .data[[group_var]], MP) %>%
    summarise(sample_mean = mean(score, na.rm = TRUE), .groups = "drop")
  
  plot_data <- sample_mean %>%
    group_by(study, .data[[group_var]], MP) %>%
    summarise(group_mean = mean(sample_mean, na.rm = TRUE), .groups = "drop")
  
  # Z-score per study per MP
  plot_data <- plot_data %>%
    group_by(study, MP) %>%
    mutate(
      mp_mean = mean(group_mean, na.rm = TRUE),
      mp_sd   = sd(group_mean, na.rm = TRUE),
      z_score = ifelse(mp_sd > 0, (group_mean - mp_mean) / mp_sd, 0)
    ) %>%
    ungroup()
  
  z_lim <- max(abs(plot_data$z_score), na.rm = TRUE)
  z_lim <- max(z_lim, 0.5)
  
  # Study stats for facet labels
  study_stats <- data %>%
    group_by(study) %>%
    summarise(n_cells = n_distinct(cell), n_samples = n_distinct(orig.ident), .groups = "drop") %>%
    mutate(study_label = paste0(study, " (n=", n_samples, "; ", comma(n_cells), " cells)"))
  
  plot_data <- plot_data %>% left_join(study_stats[, c("study", "study_label")], by = "study")
  
  # Display labels
  plot_data$mp_display <- factor(mp_display[plot_data$MP], levels = rev(mp_display[mp_ordered]))
  
  ggplot(plot_data, aes(x = .data[[group_var]], y = mp_display)) +
    geom_point(
      aes(fill = z_score),
      size = 5.5,
      shape = 21, colour = "grey30", stroke = 0.3
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-z_lim, z_lim),
      name = "Z"
    ) +
    # Show raw value inside bubbles for clarity (especially when no comparison group)
    geom_text(
      aes(label = sprintf("%.3f", group_mean)),
      size = 1.3, colour = "grey20"
    ) +
    facet_wrap(~ study_label, scales = "free_x", ncol = 3) +
    labs(x = NULL, y = NULL, title = title, subtitle = paste0("n=", total_samples, " samples")) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey45", margin = margin(b = 8)),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7, face = "italic"),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.key.size = unit(0.35, "cm")
    )
}


####################
# 8) Variables
####################
plot_configs <- list(
  list(var = "Gender", title = "Gender", filter = NULL),
  list(var = "Age_Group", title = "Age (>60)", filter = NULL),
  list(var = "Race/Ethnicity", title = "Race / Ethnicity", filter = NULL),
  list(var = "Tumor Location", title = "Tumor Location", filter = NULL),
  list(var = "Tumor Type", title = "Tumor Type", filter = NULL),
  list(var = "Grade (Differentiation)", title = "Grade (Differentiation)", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (All)", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (Tx-Naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Stage AJCC", title = "Stage AJCC (Post-Tx)", filter = "Treatment == 'Post'"),
  list(var = "Technology", title = "Sequencing Technology", filter = NULL),
  list(var = "Treatment", title = "Treatment Timing", filter = NULL),
  list(var = "Clinical response", title = "Clinical Response (Predictive, Tx-naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Clinical response", title = "Clinical Response (Post-Tx)", filter = "Treatment == 'Post'")
)

####################
# 9) Combined report — dot-plot heatmap + lollipop for binary vars
####################
cat("Generating combined report...\n")
pdf("Auto_clinical_assoc_mp_ucell_combined.pdf", width = 14, height = 10)
for (cfg in plot_configs) {
  # Main dot-plot heatmap
  p_dot <- tryCatch(
    plot_dot_heatmap(cell_long, cfg$var, cfg$title, cfg$filter),
    error = function(e) { message("  Skipping dot plot for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p_dot)) print(p_dot)
  
  # Lollipop difference plot (only for 2-level variables)
  p_lol <- tryCatch(
    plot_lollipop_diff(cell_long, cfg$var, cfg$title, cfg$filter),
    error = function(e) { message("  Skipping lollipop for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p_lol)) print(p_lol)
}
dev.off()

####################
# 10) Per-study report — faceted dot plots
####################
cat("Generating per-study report...\n")
pdf("Auto_clinical_assoc_mp_ucell_per_study.pdf", width = 18, height = 12)
for (cfg in plot_configs) {
  p <- tryCatch(
    plot_dot_faceted(cell_long, cfg$var, paste0(cfg$title, " — per Study"), cfg$filter),
    error = function(e) { message("  Skipping faceted plot for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p)) print(p)
}
dev.off()

####################
# 11) Summary CSV
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "new_update",
  "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- lapply(plot_configs, function(cfg) {
  prep <- compute_mp_group_means(cell_long, cfg$var, cfg$filter)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  
  prep$plot_data %>%
    mutate(MP_label = mp_labels[MP]) %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      clinical_variable = cfg$var,
      level = .data[[cfg$var]],
      MP = MP,
      MP_label = MP_label,
      MP_group = mp_to_group[MP],
      mean_of_sample_means = group_mean,
      samples_in_level = n_samples,
      cells_in_level = n_cells,
      total_samples = prep$total_samples
    )
})

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(summary_dir, "Auto_clinical_assoc_mp_ucell_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

message("Saved: Auto_clinical_assoc_mp_ucell_combined.pdf")
message("Saved: Auto_clinical_assoc_mp_ucell_per_study.pdf")
message(sprintf("Saved: %s", summary_path))