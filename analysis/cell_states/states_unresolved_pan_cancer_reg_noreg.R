####################
# Auto_states_unresolved_pan_cancer_reg_noreg.R
# Unified unresolved-cell pan-cancer subclassification for reg and noreg modes.
#
# For each mode:
#   - classify unresolved cells by top 3CA UCell signature
#   - plot per-cell heatmap in states_topmp_hybrid style (keep CNA + CC + Study)
#
# Output:
#   ref_outs/Auto_topmp_v2_reg_unresolved_pan_cancer.{csv,rds}
#   ref_outs/Auto_topmp_v2_noreg_unresolved_pan_cancer.{csv,rds}
#   ref_outs/Auto_topmp_v2_reg_noreg_unresolved_pan_cancer_heatmap.pdf
#   ref_outs/Auto_topmp_v2_reg_noreg_unresolved_pan_cancer_barplot.pdf
#   updates/new_updates/summaries/Auto_topmp_v2_reg_noreg_unresolved_pan_cancer_summary.csv
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")

tmdata_all <- readRDS("EAC_Ref_epi.rds")
meta_full_epi <- readRDS("meta_full_epi.rds")
state_reg <- readRDS("Auto_topmp_v2_reg_states_B.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")

cc_genes_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
cell_cycle_genes <- read.csv(cc_genes_path, header = TRUE, stringsAsFactors = FALSE)[, 1:3]

common_cells <- intersect(Cells(tmdata_all), rownames(ucell_3ca))
tmdata_all <- tmdata_all[, common_cells]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]
state_reg <- state_reg[common_cells]
state_noreg <- state_noreg[common_cells]

cna_cells <- intersect(rownames(meta_full_epi), common_cells)
cna_status <- rep(NA_character_, length(common_cells))
names(cna_status) <- common_cells
cna_status[cna_cells] <- as.character(meta_full_epi[cna_cells, "classification"])

cc_consensus <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))

THRESHOLD <- 0.10

classify_mode <- function(state_vec, mode_name) {
  unresolved_cells <- names(state_vec)[state_vec == "Unresolved"]
  unresolved_cells <- intersect(unresolved_cells, rownames(ucell_3ca))
  if (length(unresolved_cells) == 0) return(NULL)

  sub_scores <- ucell_3ca[unresolved_cells, , drop = FALSE]
  
  sample_var   <- tmdata_all$orig.ident[unresolved_cells]
  study_var    <- tmdata_all$study[unresolved_cells]
  
  clust_df <- as.data.frame(sub_scores)
  clust_df$.cell   <- rownames(sub_scores)
  clust_df$.sample <- sample_var
  clust_df$.study  <- study_var

  mps <- colnames(sub_scores)

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(mps), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()

  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(mps), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, mps])
  rownames(mp_adj) <- clust_centered$.cell

  for (mp in mps) {
    cell_studies <- clust_centered$.study
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[cell_studies, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0

  mp_matrix <- t(mp_adj)
  temp_obj <- CreateSeuratObject(counts = mp_matrix, assay = "MPs")
  temp_obj <- ScaleData(temp_obj, assay = "MPs", layer = "counts", features = rownames(temp_obj), verbose = FALSE)
  n_pcs <- min(30, nrow(mp_matrix) - 1)
  temp_obj <- RunPCA(temp_obj, features = rownames(temp_obj), npcs = n_pcs, verbose = FALSE)
  temp_obj <- FindNeighbors(temp_obj, reduction = "pca", dims = 1:n_pcs, graph.name = "MPs_snn", verbose = FALSE)
  temp_obj <- FindClusters(temp_obj, graph.name = "MPs_snn", resolution = 1, verbose = FALSE)
  
  subclass <- paste0("State_", as.numeric(as.character(temp_obj$seurat_clusters)) + 1)

  out <- data.frame(
    cell = unresolved_cells,
    mode = mode_name,
    unresolved_subclass = subclass,
    stringsAsFactors = FALSE
  )

  write.csv(out, paste0("Auto_topmp_v2_", mode_name, "_unresolved_subclass.csv"), row.names = FALSE)
  saveRDS(out, paste0("Auto_topmp_v2_", mode_name, "_unresolved_subclass.rds"))
  saveRDS(mp_adj, paste0("Auto_topmp_v2_", mode_name, "_unresolved_mp_adj.rds"))
  list(df = out, mp_adj = mp_adj)
}

make_unresolved_heatmap <- function(res_list, mode_name) {
  if (is.null(res_list) || is.null(res_list$df)) return(NULL)

  df_mode <- res_list$df
  mp_adj <- res_list$mp_adj
  cells <- df_mode$cell
  scores <- mp_adj[cells, , drop = FALSE]
  split_vec <- factor(df_mode$unresolved_subclass)

  set.seed(42)
  max_cells <- 6000
  subtype_cells <- split(cells, split_vec)
  subtype_counts <- table(split_vec)
  subtype_fracs <- subtype_counts / sum(subtype_counts)
  cells_per_subtype <- pmax(round(subtype_fracs * max_cells), 20)
  cells_to_plot <- unlist(mapply(
    function(cset, n) sample(cset, min(length(cset), n)),
    subtype_cells,
    cells_per_subtype[names(subtype_cells)],
    SIMPLIFY = FALSE
  ), use.names = FALSE)

  sub_scores <- t(scores[cells_to_plot, , drop = FALSE])
  rownames(sub_scores) <- gsub("\\.", " ", gsub("^X3CA_", "3CA_", rownames(sub_scores)))
  split_plot <- factor(df_mode$unresolved_subclass[match(cells_to_plot, df_mode$cell)])
  split_plot <- droplevels(split_plot)

  subtype_cols <- setNames(hue_pal()(length(levels(split_plot))), levels(split_plot))
  study_vals <- tmdata_all@meta.data[cells_to_plot, "study"]
  study_cols <- setNames(
    DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"),
    unique(tmdata_all$study)
  )

  col_ann <- HeatmapAnnotation(
    Subclass = split_plot,
    CNA = cna_status[cells_to_plot],
    CC_score = cc_score[cells_to_plot],
    Study = study_vals,
    col = list(
      Subclass = subtype_cols,
      CNA = c(cna_malignant = "black", cna_unresolved = "grey70"),
      CC_score = colorRamp2(c(0, max(cc_score[cells_to_plot], na.rm = TRUE)), c("white", "darkgreen")),
      Study = study_cols
    ),
    annotation_name_side = "left",
    na_col = "white"
  )

  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  Heatmap(
    sub_scores,
    name = "Adj score",
    col = colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3")),
    top_annotation = col_ann,
    column_split = split_plot,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    use_raster = TRUE,
    raster_quality = 5,
    column_title = paste0("Unresolved subclassification (", mode_name, ")")
  )
}

res_reg <- if ("reg" %in% requested_modes) classify_mode(state_reg, "reg") else NULL
res_noreg <- if ("noreg" %in% requested_modes) classify_mode(state_noreg, "noreg") else NULL

pdf("Auto_topmp_v2_reg_noreg_unresolved_heatmap.pdf", width = 18, height = 10, useDingbats = FALSE)
if (!is.null(res_reg)) {
  ht_reg <- make_unresolved_heatmap(res_reg, "reg")
  if (!is.null(ht_reg)) draw(ht_reg, merge_legend = TRUE)
}
if (!is.null(res_noreg)) {
  ht_noreg <- make_unresolved_heatmap(res_noreg, "noreg")
  if (!is.null(ht_noreg)) draw(ht_noreg, merge_legend = TRUE)
}
dev.off()



summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- bind_rows(res_reg$df, res_noreg$df) %>%
  count(mode, unresolved_subclass, name = "cells") %>%
  group_by(mode) %>%
  mutate(pct = 100 * cells / sum(cells)) %>%
  ungroup()
write.csv(summary_df, file.path(summary_dir, "Auto_topmp_v2_reg_noreg_unresolved_summary.csv"), row.names = FALSE)

message("Saved unified unresolved outputs (reg+noreg).")
