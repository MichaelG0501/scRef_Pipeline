####################
# Analysis registry:
#   Status: active
#   Script: analysis/publication/poster_section6_flot_remodelling.R
#   Inputs:
#     PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_pathway_deltas.csv
#     PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_mp_expression_changes.csv
#     PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_matched_response/pseudobulk_deg/Auto_pdo_flot_matched_deg_*.csv
#   Outputs:
#     ref_outs/publication/section6/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section6_flot_remodelling.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(clusterProfiler)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(tibble)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section6"
out_dir <- pub_section_dir(section)
flot_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_matched_response"

# States to include (excluding 3CA_EMT)
flot_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia",
                  "SMG-like Metaplasia", "Stress-adaptive")
state_file_map <- c(
  "Classic Proliferative" = "Classic_Proliferative",
  "Basal to Intestinal Metaplasia" = "Basal_to_Intest_Meta",
  "SMG-like Metaplasia" = "SMG_like_Metaplasia",
  "Stress-adaptive" = "Stress_adaptive"
)

# ==============================================================================
# FIGURE 1: DEG-based pathway enrichment per state after FLOT
# ==============================================================================
cat("Performing DEG pathway enrichment per state...\n")

deg_dir <- file.path(flot_dir, "pseudobulk_deg")
all_deg_enrichments <- list()

# Hallmark enrichment keeps the poster at pathway-level resolution.
hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") |>
  dplyr::transmute(gs_name, gene_symbol)

for (state_name in flot_states) {
  file_suffix <- state_file_map[state_name]
  deg_file <- file.path(deg_dir, paste0("Auto_pdo_flot_matched_deg_", file_suffix, ".csv"))
  if (!file.exists(deg_file)) next

  degs <- read_csv(deg_file, show_col_types = FALSE)

  # Get upregulated genes (FLOT vs untreated)
  up_genes <- degs |>
    filter(logFC > 0.25, FDR < 0.05) |>
    arrange(desc(logFC)) |>
    pull(gene) |>
    head(200)
    
  # Get downregulated genes (FLOT vs untreated)
  down_genes <- degs |>
    filter(logFC < -0.25, FDR < 0.05) |>
    arrange(logFC) |>
    pull(gene) |>
    head(200)

  if (length(up_genes) >= 10) {
    ego_up <- tryCatch(
      enricher(gene = up_genes, TERM2GENE = hallmark, pvalueCutoff = 0.20, qvalueCutoff = 0.50),
      error = function(e) NULL
    )
    if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
      all_deg_enrichments[[paste0(state_name, "_Up")]] <- as.data.frame(ego_up) |>
        slice_min(p.adjust, n = 5, with_ties = FALSE) |>
        mutate(state = state_name, direction = "Up")
    }
  }

  if (length(down_genes) >= 10) {
    ego_down <- tryCatch(
      enricher(gene = down_genes, TERM2GENE = hallmark, pvalueCutoff = 0.20, qvalueCutoff = 0.50),
      error = function(e) NULL
    )
    if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
      all_deg_enrichments[[paste0(state_name, "_Down")]] <- as.data.frame(ego_down) |>
        slice_min(p.adjust, n = 5, with_ties = FALSE) |>
        mutate(state = state_name, direction = "Down")
    }
  }
}

if (length(all_deg_enrichments) > 0) {
  enrich_combined <- bind_rows(all_deg_enrichments) |>
    mutate(signed_log10 = ifelse(direction == "Up", -log10(p.adjust), log10(p.adjust)),
           abs_log10 = -log10(p.adjust),
           state = factor(state, levels = flot_states),
           # Clean up hallmark names
           Description = gsub("^HALLMARK_", "", ID),
           Description = gsub("_", " ", Description),
           Description = tolower(Description),
           Description = paste0(toupper(substr(Description, 1, 1)), substring(Description, 2)),
           term_short = substr(Description, 1, 45))

  # Determine order: group each pathway by the state/direction where it is most significant
  pathway_ordering <- enrich_combined |>
    group_by(term_short) |>
    slice_max(abs_log10, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      state_num = as.numeric(factor(state, levels = flot_states)),
      dir_num = ifelse(direction == "Up", 1, 2)
    ) |>
    arrange(state_num, dir_num, desc(abs_log10)) |>
    pull(term_short)

  enrich_combined <- enrich_combined |>
    mutate(term_short = factor(term_short, levels = rev(pathway_ordering)))

  plot_term_levels <- c(levels(enrich_combined$term_short), " ")
  enrich_combined <- enrich_combined |> mutate(term_short = factor(as.character(term_short), levels = plot_term_levels))

  ann <- expand_grid(state = factor(flot_states, levels = flot_states)) |>
         mutate(term_short = factor(" ", levels = plot_term_levels))

  pw_clip <- max(1, quantile(abs(enrich_combined$signed_log10), 0.95, na.rm = TRUE))

  p <- ggplot(enrich_combined, aes(x = state, y = term_short)) +
    geom_tile(data = ann, aes(x = state, y = term_short), inherit.aes = FALSE,
              width = 0.92, height = 0.6, fill = PUB_STATE_COLOURS[as.character(ann$state)], colour = NA) +
    geom_point(aes(size = Count, colour = signed_log10), alpha = 0.9) +
    scale_colour_gradient2(low = "#245F7B", mid = "grey90", high = "#B63E2F", 
                           limits = c(-pw_clip, pw_clip), oob = squish, midpoint = 0,
                           name = "Signed\n-log10(FDR)") +
    scale_size_continuous(range = c(2.5, 7), name = "Genes") +
    scale_x_discrete(position = "top", expand = expansion(mult = c(0.01, 0.01)), drop = FALSE) +
    scale_y_discrete(drop = FALSE, labels = function(x) ifelse(x == " ", "", as.character(x))) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    pub_theme(14) +
    theme(
      axis.text.x.top = element_text(angle = 30, hjust = 0, vjust = 0, size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.position = "right",
      plot.margin = margin(10, 40, 10, 10)
    )
  save_pub_gg(p, section, "s6_deg_pathway_enrichment", width = 9, height = 7)
  write_csv(enrich_combined, file.path(out_dir, "tables", "s6_deg_pathway_enrichment.csv"))
} else {
  make_placeholder(section, "s6_deg_pathway_enrichment", "DEG pathway enrichment",
                   "No significant enrichments found from FLOT DEGs.")
}

# ==============================================================================
# FIGURE 2: Mean pathway response matrix
# ==============================================================================
cat("Generating pathway response matrix...\n")
path_file <- file.path(flot_dir, "Auto_pdo_flot_matched_pathway_deltas.csv")
if (file.exists(path_file)) {
  pathway_delta <- read_csv(file.path(flot_dir, "Auto_pdo_flot_matched_comb_delta.csv"), show_col_types = FALSE) |>
    mutate(state = clean_state(state),
           pathway = ifelse(pathway == "CCSIG", "cell cycle signatures", pathway)) |>
    filter(state %in% flot_states)

  # Compute mean delta
  comb_mean <- pathway_delta |>
    group_by(state, pathway) |>
    summarise(delta_score = mean(delta_score, na.rm = TRUE), .groups = "drop")

  sig_label <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    ""
  }

  comb_sig <- pathway_delta |>
    group_by(state, pathway) |>
    summarise(
      p = tryCatch(t.test(delta_score, mu = 0)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) |>
    mutate(label = sapply(p, sig_label))

  pathway_groups <- c("Proliferation / \nDNA repair", "Stress / \nimmune response", "Stress adaptation", "Metabolism", "Lineage")
  pw_order <- c("DNA repair", "E2F targets", "G2M checkpoint", "cell cycle signatures",
                "Apoptosis", "Interferon response", "p53 pathway", "TNFa/NFkB",
                "EMT", "Hypoxia", "Unfolded protein response",
                "Oxidative phosphorylation", "Xenobiotic metabolism",
                "Intestinal Metaplasia", "Columnar Progenitor")
  
  group_map <- c(
    rep("Proliferation / \nDNA repair", 4), 
    rep("Stress / \nimmune response", 4), 
    rep("Stress adaptation", 3), 
    rep("Metabolism", 2), 
    rep("Lineage", 2)
  )
  names(group_map) <- pw_order

  df_plot <- comb_mean |>
    filter(pathway %in% pw_order) |>
    mutate(
      pathway = factor(pathway, levels = rev(pw_order)),
      state = factor(state, levels = flot_states),
      group = factor(group_map[as.character(pathway)], levels = pathway_groups)
    )

  plot_term_levels <- c(levels(df_plot$pathway), " ")
  df_plot <- df_plot |> mutate(pathway = factor(as.character(pathway), levels = plot_term_levels))

  ann <- expand_grid(state = factor(flot_states, levels = flot_states),
                     group = factor(pathway_groups[1], levels = pathway_groups)) |>
         mutate(pathway = factor(" ", levels = plot_term_levels))
         
  pw_clip <- max(0.2, quantile(abs(df_plot$delta_score), 0.95, na.rm = TRUE))

  p <- ggplot(df_plot, aes(x = state, y = pathway, fill = delta_score)) +
    geom_tile(width = 0.92, height = 0.92, colour = "white", linewidth = 0.18) +
    geom_tile(data = ann, aes(x = state, y = pathway), inherit.aes = FALSE,
              width = 0.92, height = 0.6, fill = PUB_STATE_COLOURS[as.character(ann$state)], colour = NA) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#245F7B", mid = "white", high = "#B63E2F", 
                         limits = c(-pw_clip, pw_clip), oob = squish,
                         name = "Treated - Untreated\n\u0394 score") +
    scale_x_discrete(position = "top", expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_discrete(labels = function(x) ifelse(x == " ", "", as.character(x))) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    pub_theme(10.5) +
    theme(
      axis.text.x.top = element_text(angle = 30, hjust = 0, vjust = 0, size = 9),
      axis.text.y = element_text(size = 9.5),
      strip.placement = "outside",
      strip.background = element_rect(fill = "#F1F5F9", colour = "#CBD5E1"),
      strip.text.y.left = element_text(angle = 0, size = 10, face = "bold"),
      panel.spacing.y = unit(3, "mm"),
      legend.position = "right",
      plot.margin = margin(10, 40, 10, 10)
    )

  sig_df <- comb_sig |> 
    filter(pathway %in% pw_order) |>
    mutate(pathway = factor(pathway, levels = plot_term_levels),
           state = factor(state, levels = flot_states),
           group = factor(group_map[as.character(pathway)], levels = pathway_groups))

  p <- p + geom_text(data = sig_df, aes(label = label, x = state, y = pathway), 
                     inherit.aes = FALSE, vjust = 0.75, size = 4.5, fontface = "bold")

  save_pub_gg(p, section, "s6_pathway_response_matrix", width = 8.5, height = 5.5)

  write_csv(comb_mean, file.path(out_dir, "tables", "s6_pathway_response_matrix.csv"))
}

cat("Section 6 complete.\n")
