####################
# Analysis registry:
#   Status: active
#   Script: analysis/publication/poster_section2_genetics_regulons.R
#   Inputs:
#     ref_outs/Auto_cna_subclone_expression_v2/rds/Auto_v2_visualisation_results.rds
#     ref_outs/final_mp_scenic/ when available
#   Outputs:
#     ref_outs/publication/section2/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section2_genetics_regulons.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(AUCell)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section2"
out_dir <- pub_section_dir(section)

# ==============================================================================
# FIGURE 1: SCENIC state-linked regulon network
# ==============================================================================
scenic_dir <- file.path(SCREF_REF_OUTS_DIR, "final_mp_scenic")
auc_file <- file.path(scenic_dir, "Auto_final_mp_scenic_regulon_auc.rds")
cells_file <- file.path(scenic_dir, "Auto_final_mp_scenic_selected_cells.csv")

if (file.exists(auc_file) && file.exists(cells_file)) {
  regulon_auc <- readRDS(auc_file)
  auc_mat <- AUCell::getAUC(regulon_auc)
  
  selected_df <- read_csv(cells_file, show_col_types = FALSE)
  
  state_df <- selected_df |> 
    mutate(final_state = clean_state(final_state)) |>
    filter(final_state %in% PUB_STATE_ORDER) |>
    mutate(final_state = factor(final_state, levels = PUB_STATE_ORDER))
    
  state_cells <- as.character(state_df$cell)
  state_label_map <- setNames(as.character(state_df$final_state), state_cells)
  
  common_cells <- intersect(colnames(auc_mat), state_cells)
  state_auc_mat <- auc_mat[, common_cells, drop = FALSE]
  state_label_map <- state_label_map[colnames(state_auc_mat)]
  
  state_mean_auc_mat <- sapply(split(names(state_label_map), state_label_map), function(cells) {
    rowMeans(state_auc_mat[, cells, drop = FALSE], na.rm = TRUE)
  })
  state_mean_auc_mat <- as.matrix(state_mean_auc_mat)
  state_mean_auc_mat <- state_mean_auc_mat[, PUB_STATE_ORDER, drop = FALSE]
  
  n_top_global <- 100
  min_per_category <- 5
  
  state_per_cat_regs <- lapply(PUB_STATE_ORDER, function(st) {
    vals <- sort(state_mean_auc_mat[, st], decreasing = TRUE)
    names(vals)[seq_len(min(min_per_category, length(vals)))]
  })
  state_guaranteed_regs <- unique(unlist(state_per_cat_regs))
  state_guaranteed_regs <- state_guaranteed_regs[!is.na(state_guaranteed_regs)]
  
  state_global_top <- names(sort(apply(state_mean_auc_mat[, PUB_STATE_ORDER, drop = FALSE], 1, max, na.rm = TRUE), decreasing = TRUE))
  state_remaining <- setdiff(state_global_top, state_guaranteed_regs)
  state_n_fill <- max(0, n_top_global - length(state_guaranteed_regs))
  state_auc_top_regulons <- unique(c(state_guaranteed_regs, head(state_remaining, state_n_fill)))
  
  format_regulon_name <- function(x) {
    x <- gsub("_extended$", "", x)
    x <- gsub(" \\([0-9]+g\\)$", "", x)
    x <- gsub(" \\([0-9]+ genes\\)$", "", x)
    x
  }
  
  state_auc_edge_df <- bind_rows(lapply(PUB_STATE_ORDER, function(st) {
    vals <- state_mean_auc_mat[state_auc_top_regulons, st]
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) == 0) return(NULL)
    threshold <- quantile(vals, 0.5, na.rm = TRUE)
    keep <- names(vals)[vals >= threshold]
    if (length(keep) == 0) return(NULL)
    data.frame(
      regulon_label = format_regulon_name(keep),
      state_label = st,
      weight = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) %>% distinct(regulon_label, state_label, .keep_all = TRUE)
  
  state_auc_node_df <- data.frame(
    name = unique(c(state_auc_edge_df$regulon_label, state_auc_edge_df$state_label)),
    stringsAsFactors = FALSE
  ) %>% mutate(
    node_type = ifelse(name %in% PUB_STATE_ORDER, "State", "Regulon"),
    state_group = ifelse(node_type == "State", name, "Regulon")
  )
  
  commonly_known_tfs <- c(
    "MYC", "TP53", "HIF1A", "STAT1", "STAT3", "STAT5A", "STAT5B", "NFKB1", "RELA", "SOX2", "SOX4", "SOX6", "SOX9", 
    "CDX1", "CDX2", "EZH2", "BCL6", "E2F1", "E2F2", "E2F3", "E2F4", "E2F7", "E2F8", "FOXA1", "FOXA2", "FOXQ1",
    "GATA4", "GATA6", "KLF4", "KLF5", "ASCL2", "CEBPB", "CEBPG", "JUN", "JUND", "FOS", "FOSL1", "FOSL2", 
    "SMAD3", "SMAD4", "HNF4A", "IRF7", "HMGB1", "MTA3", "NFIA", "NFYB", "NFYC", "PBX1", "POU3F1", "PRDM1", 
    "SP1", "SP2", "TEAD1", "TEAD4", "ZBTB7B", "ZFHX3", "ZNF143", "YBX1", "BRCA1", "AR", "ESR1", "HOXA9", "HOXB6", 
    "ATF3", "ATF4", "NRF2", "NFE2L2", "REST", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "YAP1", "TAZ", "XBP1", "FOXM1"
  )
  
  top_regulons_fallback <- state_auc_edge_df |> 
    group_by(state_label) |> 
    slice_max(weight, n = 2, with_ties = FALSE) |> 
    pull(regulon_label) |> 
    unique()
  
  state_reg_degree <- state_auc_edge_df %>% count(regulon_label, name = "n_states")
  state_auc_node_df <- state_auc_node_df %>%
    left_join(state_reg_degree, by = c("name" = "regulon_label")) %>%
    mutate(is_shared = !is.na(n_states) & n_states > 1) %>%
    mutate(label_show = ifelse(node_type == "State" | name %in% commonly_known_tfs | name %in% top_regulons_fallback, name, ""))
  
  state_auc_network_graph <- tbl_graph(
    nodes = state_auc_node_df,
    edges = state_auc_edge_df %>% transmute(from = regulon_label, to = state_label, weight = weight),
    directed = FALSE
  )
  
  set.seed(42)
  regulon_plot <- ggraph(state_auc_network_graph, layout = "stress") +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "#64748B",
                   show.legend = FALSE, lineend = "round") +
    scale_edge_width(range = c(0.2, 1.5)) +
    scale_edge_alpha(range = c(0.15, 0.6)) +
    geom_node_point(
      aes(fill = state_group, shape = node_type, 
          size = ifelse(node_type == "State", 9, ifelse(is_shared, 5, 3.5))),
      colour = "#111827", stroke = 0.4
    ) +
    scale_shape_manual(values = c(State = 21, Regulon = 22), guide = "none") +
    scale_size_identity() +
    scale_fill_manual(values = c(PUB_STATE_COLOURS, Regulon = "grey35"), na.value = "#94A3B8", guide = "none") +
    geom_node_text(
      aes(label = ifelse(node_type == "Regulon", label_show, "")),
      fontface = "plain", size = 7.5, repel = TRUE, 
      point.padding = 0.3, force = 2,
      colour = "#111827", max.overlaps = 100, bg.color = "white", bg.r = 0.15
    ) +
    geom_node_text(
      aes(label = ifelse(node_type == "State", label_show, "")),
      fontface = "bold", size = 9.0, repel = TRUE, 
      point.padding = 0.8, force = 0.5, min.segment.length = Inf,
      colour = "#111827", max.overlaps = 100, bg.color = "white", bg.r = 0.15
    ) +
    theme_void(base_size = 16) +
    theme(plot.margin = margin(10, 10, 10, 10))
  
  save_pub_gg(regulon_plot, section, "s2_scenic_regulon_network", width = 22, height = 14)
  write_csv(state_auc_edge_df, file.path(out_dir, "tables", "s2_scenic_regulon_network_edges.csv"))
} else {
  abort_missing_figure(section, "s2_scenic_regulon_network", "SCENIC regulon network",
                   "Missing AUC or selected cells outputs.", width = 10, height = 8)
}

cat("Section 2 complete.\n")
