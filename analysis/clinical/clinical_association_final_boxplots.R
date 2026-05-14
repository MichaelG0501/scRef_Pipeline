####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/clinical_association_final_boxplots.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_clinical_assoc_boxplots_final.R
# Alternative clinical association plotting for MPs and finalized states.
#
# Visualisation:
#   - One page per clinical variable
#   - MPs or states on the x-axis
#   - Thin, side-by-side sample-level boxplots for the clinical groups
#   - Separate PDFs for MP activity and finalized state proportions
#
# Statistics:
#   - Wilcoxon rank-sum for 2-level variables
#   - Kruskal-Wallis for >2-level variables
#   - Benjamini-Hochberg adjustment within each clinical variable / feature type
#
# Input:
#   ref_outs/meta_full_epi.rds
#   ref_outs/Auto_final_states.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx (sheet 3)
#
# Output:
#   ref_outs/Auto_clinical_assoc_mp_boxplots_final.pdf           (combined)
#   ref_outs/Auto_clinical_assoc_mp_boxplots_final_per_study.pdf (per-study)
#   ref_outs/Auto_clinical_assoc_state_boxplots_final.pdf           (combined)
#   ref_outs/Auto_clinical_assoc_state_boxplots_final_per_study.pdf (per-study)
#   updates/new_updates/summaries/Auto_clinical_assoc_mp_boxplots_final_stats.csv
#   updates/new_updates/summaries/Auto_clinical_assoc_mp_boxplots_final_per_study_stats.csv
#   updates/new_updates/summaries/Auto_clinical_assoc_state_boxplots_final_stats.csv
#   updates/new_updates/summaries/Auto_clinical_assoc_state_boxplots_final_per_study_stats.csv
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(scales)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# 1) Load data
####################
meta_full_epi <- readRDS("meta_full_epi.rds")
final_states <- readRDS("Auto_final_states.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")

common_cells <- intersect(rownames(ucell_scores), rownames(ucell_3ca))
ucell_scores <- cbind(
  ucell_scores[common_cells, , drop = FALSE], 
  ucell_3ca[common_cells, c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1"), drop = FALSE]
)

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
  sheet = 3,
  skip = 1
) %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

####################
# 2) Constants
####################
####################
# Clinical variable plot configurations
# Matching clinical_variable_plots_topmp_v2B_reg_noreg.R (except Race)
####################
plot_configs <- list(
  list(var = "Gender", title = "Gender", filter = NULL),
  list(var = "Age_Group", title = "Age (>60)", filter = NULL),
  list(var = "Tumor Location", title = "Tumor Location", filter = NULL),
  list(var = "Tumor Type", title = "Tumor Type", filter = NULL),
  list(var = "Grade (Differentiation)", title = "Grade (Differentiation)", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (All)", filter = NULL),
  list(var = "Stage AJCC", title = "Stage AJCC (Tx-Naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Stage AJCC", title = "Stage AJCC (Post-Tx)", filter = "Treatment == 'Post'"),
  list(var = "Technology", title = "Sequencing Technology", filter = NULL),
  list(var = "Treatment", title = "Treatment Timing", filter = NULL),
  list(var = "Clinical response", title = "Clinical Response (Predictive, Tx-Naive)", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Clinical response", title = "Clinical Response (Post-Tx)", filter = "Treatment == 'Post'")
)

clinical_level_orders <- list(
  Gender = c("Female", "Male"),
  Age_Group = c("<=60", ">60"),
  Treatment = c("Tx-naïve", "Post"),
  `Clinical response` = c("Responder", "Nonresponder")
)

mp_descriptions <- c(
  "MP1" = "G2M cycle",
  "MP2" = "MYC prolif",
  "MP5" = "IFN response",
  "MP7" = "S cycle",
  "MP8" = "Intestinal diff",
  "MP9" = "G1S cycle",
  "MP10" = "Columnar diff",
  "MP12" = "Neuro-epithelial",
  "MP13" = "Partial EMT",
  "MP14" = "Hypoxia epithelial",
  "MP15" = "Immune infiltration",
  "MP16" = "Secretory diff",
  "MP17" = "Squamous transition",
  "MP18" = "Adaptive secretory"
)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15"),
  "Cell Cycle" = c("MP1", "MP7", "MP9")
)

state_axis_labels <- c(
  "Classic Proliferative" = "Classic\nProlif",
  "Basal to Intestinal Metaplasia" = "Basal to\nIntestinal\nMetaplasia",
  "SMG-like Metaplasia" = "SMG-like\nMetaplasia",
  "Stress-adaptive" = "Stress\nadaptive",
  "Immune Infiltrating" = "Immune\nInfiltrating",
  "3CA_EMT_and_Protein_maturation" = "3CA EMT +\nProtein\nmaturation",
  "Unresolved" = "Unresolved",
  "Hybrid" = "Hybrid"
)

####################
# 3) Helpers
####################
make_group_order <- function(values, group_var) {
  present_levels <- unique(as.character(values))
  preferred_levels <- clinical_level_orders[[group_var]]
  ordered_levels <- character(0)
  if (!is.null(preferred_levels)) {
    ordered_levels <- preferred_levels[preferred_levels %in% present_levels]
  }
  remaining_levels <- setdiff(present_levels, ordered_levels)
  c(ordered_levels, sort(remaining_levels))
}

make_group_palette <- function(levels_use) {
  base_cols <- c(
    "Female" = "#E64B35",
    "Male" = "#4DBBD5",
    "<=60" = "#00A087",
    ">60" = "#3C5488",
    "Tx-naïve" = "#F39B7F",
    "Post" = "#8491B4",
    "Responder" = "#91D1C2",
    "Nonresponder" = "#DC0000"
  )
  palette <- base_cols[levels_use]
  missing_levels <- setdiff(levels_use, names(base_cols))
  if (length(missing_levels) > 0) {
    palette[missing_levels] <- scales::hue_pal(l = 65, c = 100)(length(missing_levels))
  }
  palette[levels_use]
}

significance_label <- function(p_val) {
  dplyr::case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

make_state_levels <- function(state_vec) {
  core_states <- c(
    "Classic Proliferative",
    "Basal to Intestinal Metaplasia",
    "SMG-like Metaplasia",
    "Stress-adaptive",
    "Immune Infiltrating"
  )
  preferred_extra_states <- c("3CA_EMT_and_Protein_maturation")
  present_states <- unique(as.character(state_vec))
  other_extra_states <- setdiff(
    present_states,
    c(core_states, preferred_extra_states, "Unresolved", "Hybrid")
  )
  c(
    core_states[core_states %in% present_states],
    preferred_extra_states[preferred_extra_states %in% present_states],
    sort(other_extra_states)
  )
}

apply_optional_filter <- function(data, filter_expr = NULL) {
  if (is.null(filter_expr)) {
    return(data)
  }
  data %>% filter(!!rlang::parse_expr(filter_expr))
}

compute_feature_stats <- function(sample_df, clinical_var, feature_type, feature_label_map, digits = 3) {
  feature_splits <- split(sample_df, as.character(sample_df$feature))
  stats_rows <- lapply(feature_splits, function(df_feature) {
    df_feature <- df_feature %>%
      filter(!is.na(group), !is.na(value))

    present_groups <- unique(as.character(df_feature$group))
    if (length(present_groups) < 2 || n_distinct(df_feature$orig.ident) < 3) {
      return(data.frame(
        clinical_variable = clinical_var,
        feature_type = feature_type,
        feature = as.character(df_feature$feature[1]),
        feature_label = unname(feature_label_map[as.character(df_feature$feature[1])]),
        test = NA_character_,
        n_groups = length(present_groups),
        n_samples = n_distinct(df_feature$orig.ident),
        p_value = NA_real_,
        group_summary = NA_character_,
        median_range = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    test_name <- if (length(present_groups) == 2) "wilcoxon" else "kruskal"
    test_res <- tryCatch(
      suppressWarnings(
        if (test_name == "wilcoxon") {
          wilcox.test(value ~ group, data = df_feature, exact = FALSE)
        } else {
          kruskal.test(value ~ group, data = df_feature)
        }
      ),
      error = function(e) NULL
    )

    median_df <- df_feature %>%
      group_by(group) %>%
      summarise(
        n_samples = n_distinct(orig.ident),
        median_value = median(value, na.rm = TRUE),
        mean_value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

    data.frame(
      clinical_variable = clinical_var,
      feature_type = feature_type,
      feature = as.character(df_feature$feature[1]),
      feature_label = unname(feature_label_map[as.character(df_feature$feature[1])]),
      test = test_name,
      n_groups = length(present_groups),
      n_samples = n_distinct(df_feature$orig.ident),
      p_value = if (is.null(test_res)) NA_real_ else test_res$p.value,
      group_summary = paste0(
        as.character(median_df$group),
        " (n=",
        median_df$n_samples,
        ", median=",
        formatC(median_df$median_value, format = "f", digits = digits),
        ")",
        collapse = " | "
      ),
      median_range = diff(range(median_df$median_value, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })

  stats_df <- bind_rows(stats_rows)
  if (nrow(stats_df) == 0) {
    return(stats_df)
  }

  stats_df %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significance = significance_label(p_value)
    )
}

plot_feature_boxplot <- function(sample_df, stats_df, title_text, subtitle_text, y_label, feature_labels, feature_type) {
  sample_df <- sample_df %>%
    filter(!is.na(group), !is.na(feature), !is.na(value))

  if (nrow(sample_df) == 0 || n_distinct(sample_df$group) < 2) {
    return(NULL)
  }

  group_levels <- levels(sample_df$group)
  if (is.null(group_levels)) {
    group_levels <- unique(as.character(sample_df$group))
  }

  legend_counts <- sample_df %>%
    distinct(orig.ident, group) %>%
    count(group, name = "n_samples")

  legend_labels <- setNames(
    paste0(
      iconv(as.character(legend_counts$group), from = "UTF-8", to = "ASCII//TRANSLIT"),
      " (n=",
      legend_counts$n_samples,
      ")"
    ),
    as.character(legend_counts$group)
  )

  palette <- make_group_palette(group_levels)
  box_width <- min(0.22, 0.8 / max(length(group_levels), 3))
  y_range <- range(sample_df$value, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- if (feature_type == "state") 5 else 0.02
  }
  y_offset <- if (feature_type == "state") max(2, 0.08 * y_span) else max(0.01, 0.08 * y_span)

  annot_df <- sample_df %>%
    group_by(feature) %>%
    summarise(y_pos = max(value, na.rm = TRUE) + y_offset, .groups = "drop") %>%
    left_join(stats_df %>% select(feature, significance, p_value), by = "feature") %>%
    mutate(label = ifelse(!is.na(p_value) & p_value < 0.05, significance, ""))

  p <- ggplot(sample_df, aes(x = feature, y = value, fill = group, color = group)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.8,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      alpha = 0.7,
      size = 1.0,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = annot_df %>% filter(label != ""),
      aes(x = feature, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
    scale_color_manual(values = palette, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = feature_labels) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = y_label,
      fill = "Clinical group"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "grey35"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(12, 16, 12, 12)
    )

  if (feature_type == "state") {
    p <- p +
      scale_y_continuous(
        labels = function(x) paste0(formatC(x, format = "f", digits = 0), "%"),
        expand = expansion(mult = c(0.02, 0.18))
      )
  } else {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }

  p
}

####################
# 4) Cell-level metadata tables
####################
cell_ids <- Reduce(intersect, list(
  rownames(meta_full_epi),
  rownames(ucell_scores),
  names(final_states)
))

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

####################
# 5) MP order and labels
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

retained_mps <- names(mp.genes)
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

mp_cols <- intersect(mp_tree_order_names, colnames(ucell_scores))
missing_mps <- setdiff(mp_tree_order_names, colnames(ucell_scores))
extra_mps <- intersect(c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1"), colnames(ucell_scores))
mp_cols <- c(mp_cols, extra_mps)

mp_ordered <- c()
for (grp in names(state_groups)) {
  grp_mps <- mp_cols[mp_cols %in% state_groups[[grp]]]
  mp_ordered <- c(mp_ordered, grp_mps)
}
mp_ordered <- unique(c(mp_ordered, setdiff(mp_cols, mp_ordered)))

extended_mp_descriptions <- c(
  mp_descriptions, 
  "X3CA_mp_12.Protein.maturation" = "Protein maturation", 
  "X3CA_mp_17.EMT.III" = "EMT III", 
  "X3CA_mp_30.Respiration.1" = "Respiration 1"
)

mp_feature_labels <- setNames(
  ifelse(!is.na(extended_mp_descriptions[mp_ordered]), paste0(mp_ordered, " ", extended_mp_descriptions[mp_ordered]), mp_ordered),
  mp_ordered
)

mp_axis_labels <- setNames(
  ifelse(!is.na(extended_mp_descriptions[mp_ordered]), paste0(mp_ordered, "\n", extended_mp_descriptions[mp_ordered]), mp_ordered),
  mp_ordered
)

mp_to_group <- setNames(rep("Other", length(mp_ordered)), mp_ordered)
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], mp_ordered)
  mp_to_group[grp_mps] <- grp
}

####################
# 6) MP and state sample-level tables
####################
cell_long <- as.data.frame(ucell_scores[cell_ids, mp_cols, drop = FALSE]) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(mp_cols), names_to = "MP", values_to = "score") %>%
  left_join(cell_meta, by = "cell")

state_cell_df <- cell_meta %>%
  mutate(state = as.character(final_states[cell]))

build_mp_sample_df <- function(data, group_var, filter_expr = NULL) {
  data_use <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(MP), !is.na(score), !is.na(study))

  if (nrow(data_use) == 0) {
    return(NULL)
  }

  group_levels <- make_group_order(data_use[[group_var]], group_var)

  data_use %>%
    mutate(group = factor(as.character(.data[[group_var]]), levels = group_levels)) %>%
    group_by(orig.ident, study, group, MP) %>%
    summarise(value = mean(score, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      orig.ident = orig.ident,
      study = study,
      group = group,
      feature = factor(MP, levels = mp_ordered),
      value = value
    )
}

build_state_sample_df <- function(data, group_var, filter_expr = NULL) {
  data_use <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(state), !is.na(study), !(state %in% c("Unresolved", "Hybrid")))

  if (nrow(data_use) == 0) {
    return(NULL)
  }

  group_levels <- make_group_order(data_use[[group_var]], group_var)
  state_levels <- make_state_levels(data_use$state)

  sample_meta <- data_use %>%
    transmute(
      orig.ident = orig.ident,
      study = study,
      group = as.character(.data[[group_var]])
    ) %>%
    distinct()

  sample_totals <- data_use %>%
    transmute(
      orig.ident = orig.ident,
      group = as.character(.data[[group_var]])
    ) %>%
    count(orig.ident, group, name = "total_cells")

  state_counts <- data_use %>%
    transmute(
      orig.ident = orig.ident,
      group = as.character(.data[[group_var]]),
      state = state
    ) %>%
    count(orig.ident, group, state, name = "n_cells")

  sample_meta %>%
    tidyr::crossing(state = state_levels) %>%
    left_join(state_counts, by = c("orig.ident", "group", "state")) %>%
    mutate(n_cells = replace_na(n_cells, 0L)) %>%
    left_join(sample_totals, by = c("orig.ident", "group")) %>%
    mutate(
      group = factor(group, levels = group_levels),
      feature = factor(state, levels = state_levels),
      value = 100 * n_cells / pmax(total_cells, 1)
    ) %>%
    transmute(
      orig.ident = orig.ident,
      study = study,
      group = group,
      feature = feature,
      value = value
    )
}

####################
# 7) Output paths
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates",
  "new_updates",
  "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# 8) MP clinical boxplots
####################
mp_stats_all <- list()
pdf("Auto_clinical_assoc_mp_boxplots_final.pdf", width = 18, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  mp_sample_df <- build_mp_sample_df(cell_long, cfg$var, cfg$filter)
  if (is.null(mp_sample_df) || nrow(mp_sample_df) == 0 || n_distinct(mp_sample_df$group) < 2) {
    next
  }

  stats_df <- compute_feature_stats(
    mp_sample_df,
    clinical_var = cfg$var,
    feature_type = "MP",
    feature_label_map = mp_feature_labels,
    digits = 3
  ) %>%
    mutate(
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
      feature_group = mp_to_group[as.character(feature)]
    )

  mp_stats_all[[cfg$title]] <- stats_df

  p <- plot_feature_boxplot(
    mp_sample_df,
    stats_df,
    title_text = paste0(cfg$title, " - MP activity"),
    subtitle_text = "Sample-level mean UCell scores; stars mark BH-adjusted p < 0.05 across clinical groups.",
    y_label = "Mean sample UCell score",
    feature_labels = mp_axis_labels,
    feature_type = "mp"
  )
  if (!is.null(p)) {
    print(p)
  }
}
dev.off()

mp_stats_df <- bind_rows(mp_stats_all)
write.csv(
  mp_stats_df,
  file.path(summary_dir, "Auto_clinical_assoc_mp_boxplots_final_stats.csv"),
  row.names = FALSE
)

####################
# 9) Final-state clinical boxplots
####################
state_levels_all <- make_state_levels(state_cell_df$state)
state_feature_labels <- setNames(state_levels_all, state_levels_all)
state_axis_labels[names(state_axis_labels) %in% state_levels_all] <- state_axis_labels[names(state_axis_labels) %in% state_levels_all]
missing_state_labels <- setdiff(state_levels_all, names(state_axis_labels))
if (length(missing_state_labels) > 0) {
  state_axis_labels[missing_state_labels] <- missing_state_labels
}

state_stats_all <- list()
pdf("Auto_clinical_assoc_state_boxplots_final.pdf", width = 16, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  state_sample_df <- build_state_sample_df(state_cell_df, cfg$var, cfg$filter)
  if (is.null(state_sample_df) || nrow(state_sample_df) == 0 || n_distinct(state_sample_df$group) < 2) {
    next
  }

  stats_df <- compute_feature_stats(
    state_sample_df,
    clinical_var = cfg$var,
    feature_type = "state",
    feature_label_map = state_feature_labels,
    digits = 1
  ) %>%
    mutate(
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter)
    )

  state_stats_all[[cfg$title]] <- stats_df

  p <- plot_feature_boxplot(
    state_sample_df,
    stats_df,
    title_text = paste0(cfg$title, " - final state proportions"),
    subtitle_text = "Sample-level state proportions; stars mark BH-adjusted p < 0.05 across clinical groups.",
    y_label = "Sample state proportion",
    feature_labels = state_axis_labels,
    feature_type = "state"
  )
  if (!is.null(p)) {
    print(p)
  }
}
dev.off()

state_stats_df <- bind_rows(state_stats_all)
write.csv(
  state_stats_df,
  file.path(summary_dir, "Auto_clinical_assoc_state_boxplots_final_stats.csv"),
  row.names = FALSE
)

message("Saved combined MP/state clinical boxplot PDFs and statistics.")

####################
# 10) Per-study MP clinical boxplots
####################
all_studies <- unique(cell_long$study)
all_studies <- all_studies[!is.na(all_studies) & nzchar(all_studies)]
all_studies <- sort(all_studies)

mp_stats_per_study_all <- list()
pdf("Auto_clinical_assoc_mp_boxplots_final_per_study.pdf", width = 18, height = 9, useDingbats = FALSE)
# Iterate by clinical variable first, then by study within each variable
for (cfg in plot_configs) {
  for (study_name in all_studies) {
    study_data <- cell_long %>% filter(study == study_name)
    
    mp_sample_df <- build_mp_sample_df(study_data, cfg$var, cfg$filter)
    if (is.null(mp_sample_df) || nrow(mp_sample_df) == 0 || n_distinct(mp_sample_df$group) < 2) {
      next
    }
    
    stats_df <- compute_feature_stats(
      mp_sample_df,
      clinical_var = cfg$var,
      feature_type = "MP",
      feature_label_map = mp_feature_labels,
      digits = 3
    ) %>%
      mutate(
        study = study_name,
        plot_title = cfg$title,
        filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
        feature_group = mp_to_group[as.character(feature)]
      )
    
    mp_stats_per_study_all[[paste0(cfg$title, "_", study_name)]] <- stats_df
    
    p <- plot_feature_boxplot(
      mp_sample_df,
      stats_df,
      title_text = paste0("[", study_name, "] ", cfg$title, " - MP activity"),
      subtitle_text = "Sample-level mean UCell scores; stars mark BH-adjusted p < 0.05 across clinical groups.",
      y_label = "Mean sample UCell score",
      feature_labels = mp_axis_labels,
      feature_type = "mp"
    )
    if (!is.null(p)) {
      print(p)
    }
  }
}
dev.off()

mp_stats_per_study_df <- bind_rows(mp_stats_per_study_all)
write.csv(
  mp_stats_per_study_df,
  file.path(summary_dir, "Auto_clinical_assoc_mp_boxplots_final_per_study_stats.csv"),
  row.names = FALSE
)

message("Saved per-study MP clinical boxplot PDF and statistics.")

####################
# 11) Per-study final-state clinical boxplots
####################
state_stats_per_study_all <- list()
pdf("Auto_clinical_assoc_state_boxplots_final_per_study.pdf", width = 16, height = 9, useDingbats = FALSE)
# Iterate by clinical variable first, then by study within each variable
for (cfg in plot_configs) {
  for (study_name in all_studies) {
    study_state_data <- state_cell_df %>% filter(study == study_name)
    
    state_sample_df <- build_state_sample_df(study_state_data, cfg$var, cfg$filter)
    if (is.null(state_sample_df) || nrow(state_sample_df) == 0 || n_distinct(state_sample_df$group) < 2) {
      next
    }
    
    stats_df <- compute_feature_stats(
      state_sample_df,
      clinical_var = cfg$var,
      feature_type = "state",
      feature_label_map = state_feature_labels,
      digits = 1
    ) %>%
      mutate(
        study = study_name,
        plot_title = cfg$title,
        filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter)
      )
    
    state_stats_per_study_all[[paste0(cfg$title, "_", study_name)]] <- stats_df
    
    p <- plot_feature_boxplot(
      state_sample_df,
      stats_df,
      title_text = paste0("[", study_name, "] ", cfg$title, " - final state proportions"),
      subtitle_text = "Sample-level state proportions; stars mark BH-adjusted p < 0.05 across clinical groups.",
      y_label = "Sample state proportion",
      feature_labels = state_axis_labels,
      feature_type = "state"
    )
    if (!is.null(p)) {
      print(p)
    }
  }
}
dev.off()

state_stats_per_study_df <- bind_rows(state_stats_per_study_all)
write.csv(
  state_stats_per_study_df,
  file.path(summary_dir, "Auto_clinical_assoc_state_boxplots_final_per_study_stats.csv"),
  row.names = FALSE
)

message("Saved per-study state clinical boxplot PDF and statistics.")
message("All outputs complete.")

####################
# 12) Optional stacked clinical association companion
####################
stacked_script <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "analysis",
  "clinical",
  "clinical_association_final_stacked.R"
)
skip_stacked <- tolower(Sys.getenv("SCREF_SKIP_STACKED_ASSOC", unset = "FALSE")) %in% c("true", "1", "yes", "y")
if (skip_stacked) {
  message("Skipping stacked clinical association companion because SCREF_SKIP_STACKED_ASSOC is TRUE.")
} else if (file.exists(stacked_script)) {
  message("Running stacked clinical association companion: ", stacked_script)
  source(stacked_script, local = TRUE)
} else {
  message("Stacked clinical association companion not found; boxplot outputs are complete.")
}
