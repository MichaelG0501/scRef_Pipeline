#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visium_hd_summary.R
#   Methodology: not required (legacy summary aggregation)
#   Map: analysis/ANALYSIS_MAP.md
####################
library(readr)
library(dplyr)
library(ggplot2)
library(FNN)

source("analysis/shared/scRef_config.R")
source("analysis/shared/scRef_helpers.R")
source("analysis/publication/publication_helpers.R")

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

OUT_DIR <- file.path("ref_outs", "visium_hd_outs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

binned_file_all <- file.path(OUT_DIR, "Auto_visiumhd_binned_spot_annotations.csv.gz")
segmented_file_all <- file.path(OUT_DIR, "Auto_visiumhd_segmented_spot_annotations.csv.gz")

####################
# State summaries are also generated for the malignant-only subset
####################
binned_file_malignant <- file.path(OUT_DIR, "tables", "Auto_visiumhd_binned_malignant_epithelial_state_annotations.csv.gz")
segmented_file_malignant <- file.path(OUT_DIR, "tables", "Auto_visiumhd_segmented_malignant_epithelial_state_annotations.csv.gz")
####################

visium_colours <- c(
  PUB_STATE_COLOURS,
  "Hybrid" = "#000000",
  "Unresolved" = "#999999",
  "Normal/Mixed" = "#E6E6E6"
)

clean_state_visium <- function(x) {
  dplyr::case_when(
    x == "Classic proliferation" ~ "Classic Proliferative",
    x == "Basal to intestinal metaplasia" ~ "Basal to Intestinal Metaplasia",
    x == "SMG to intestinal metaplasia" ~ "SMG-like Metaplasia",
    x == "Stress adaptive" ~ "Stress-adaptive",
    x == "Cancer-cell immune mimicry" ~ "Immune Infiltrating",
    x == "Unresolved" ~ "Unresolved",
    x == "Hybrid" ~ "Hybrid",
    TRUE ~ "Normal/Mixed"
  )
}

process_summary <- function(spot_file, mode_name, output_suffix = "") {
  if (!file.exists(spot_file)) {
    cat(sprintf("File missing: %s\n", spot_file))
    return(NULL)
  }
  
  cat(sprintf("Processing summary for %s%s...\n", mode_name, output_suffix))
  spots <- read_csv(spot_file, show_col_types = FALSE)
  
  # Ensure state column exists depending on what's available
  if ("Auto_state_B" %in% colnames(spots)) {
      spots <- spots |> mutate(state = clean_state_visium(Auto_state_B))
  } else if ("Auto_annotation_celltype" %in% colnames(spots)) {
      spots <- spots |> mutate(state = clean_state_visium(Auto_annotation_celltype))
  } else {
      cat(sprintf("No valid state column found in %s\n", spot_file))
      return(NULL)
  }
  
  # 1. State Abundance (Stacked Barplot)
  vis <- spots |>
    group_by(sample, state) |>
    summarise(spots_count = n(), .groups = "drop") |>
    group_by(sample) |>
    mutate(pct = 100 * spots_count / sum(spots_count)) |>
    ungroup() |>
    mutate(state = factor(state, levels = rev(c(PUB_STATE_ORDER, "Hybrid", "Unresolved", "Normal/Mixed"))))
  
  p_bar <- ggplot(vis, aes(x = sample, y = pct, fill = state)) +
    geom_bar(stat = "identity", width = 0.8, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = visium_colours) +
    labs(x = "Sample", y = "Proportion of Spots/Cells (%)", fill = "State") +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  ggsave(file.path(OUT_DIR, sprintf("Auto_visiumhd_%s%s_abundance_barplot.pdf", mode_name, output_suffix)), p_bar, width = 7, height = 6)
  ggsave(file.path(OUT_DIR, sprintf("Auto_visiumhd_%s%s_abundance_barplot.png", mode_name, output_suffix)), p_bar, width = 7, height = 6, dpi = 300)
  
  # 2. Colocalisation Score (kNN = 6)
  cat(sprintf("Calculating colocalisation scores for %s...\n", mode_name))
  coloc <- spots |>
    filter(state %in% PUB_STATE_ORDER) |>
    group_by(sample) |>
    group_modify(function(.x, .y) {
      coords <- as.matrix(.x[, c("pxl_row_in_fullres", "pxl_col_in_fullres")])
      if (nrow(coords) < 8) return(tibble())
      
      # FNN for fast k-nearest neighbours (k=6)
      knn_res <- FNN::get.knn(coords, k = min(6, nrow(coords) - 1))
      neigh <- knn_res$nn.index
      
      # Calculate proportion of neighbours sharing the same state
      same <- vapply(seq_len(nrow(.x)), function(i) {
         mean(.x$state[neigh[i, ]] == .x$state[i], na.rm = TRUE)
      }, numeric(1))
      
      tibble(state = .x$state, same_neighbor_score = same)
    }) |>
    ungroup() |>
    mutate(state = factor(state, levels = PUB_STATE_ORDER))
    
  p_coloc <- ggplot(coloc, aes(x = state, y = same_neighbor_score, fill = state, color = state)) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      alpha = 0.8,
      linewidth = 0.6,
      color = "black"
    ) +
    geom_point(
      position = position_jitter(width = 0.15),
      size = 0.5,
      alpha = 0.05,
      color = "black"
    ) +
    scale_fill_manual(values = visium_colours, guide = "none") +
    scale_color_manual(values = visium_colours, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = "Same-state neighbours") +
    pub_theme(14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 13, face = "bold"),
          plot.margin = margin(15, 15, 15, 15))
          
  ggsave(file.path(OUT_DIR, sprintf("Auto_visiumhd_%s%s_colocalisation.pdf", mode_name, output_suffix)), p_coloc, width = 6.0, height = 5.5)
  ggsave(file.path(OUT_DIR, sprintf("Auto_visiumhd_%s%s_colocalisation.png", mode_name, output_suffix)), p_coloc, width = 6.0, height = 5.5, dpi = 300)
}

process_summary(binned_file_all, "binned", "")
process_summary(segmented_file_all, "segmented", "")
process_summary(binned_file_malignant, "binned", "_malignant_only")
process_summary(segmented_file_malignant, "segmented", "_malignant_only")

cat("Visium HD summary generation complete.\n")
