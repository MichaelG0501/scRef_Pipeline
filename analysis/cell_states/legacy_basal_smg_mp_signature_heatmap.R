####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_basal_smg_mp_signature_heatmap.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_basal_smg_mp_signature_heatmap.R
# Basal vs SMG-like state heatmap grouped by the top state-defining MP.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#
# Output:
#   ref_outs/Auto_basal_smg_mp_signature_heatmap/
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(Matrix)
  library(grid)
  library(msigdbr)
})

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- "Auto_basal_smg_mp_signature_heatmap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MAX_CELLS_PER_MP <- Inf
NCORES_UCELL <- 4

target_states <- c("Basal to Intestinal Metaplasia", "SMG-like Metaplasia")
basal_mps <- c("MP17", "MP14", "MP5", "MP10", "MP8")
smg_mps <- c("MP18", "MP16")
target_mp_order <- c(basal_mps, smg_mps)

mp_descriptions <- c(
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)"
)

state_cols <- c(
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00"
)

mp_cols <- c(
  "MP17" = "#4DAF4A",
  "MP14" = "#F781BF",
  "MP5"  = "#377EB8",
  "MP10" = "#A65628",
  "MP8"  = "#984EA3",
  "MP18" = "#FF7F00",
  "MP16" = "#FFD92F"
)

label_mp <- function(mp) {
  desc <- mp_descriptions[mp]
  desc[is.na(desc)] <- mp[is.na(desc)]
  paste0(mp, " ", desc)
}

match_genes <- function(genes, universe) {
  genes <- unique(toupper(trimws(genes)))
  universe_upper <- toupper(universe)
  hits <- universe[match(genes, universe_upper)]
  unique(hits[!is.na(hits)])
}

get_data_layer <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
}

z_scale_rows <- function(mat, zlim = 2.5) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z[z > zlim] <- zlim
  z[z < -zlim] <- -zlim
  z
}

normalize_within_sample <- function(score_mat, sample_vector, zlim = 2.5) {
  samples <- unique(sample_vector)
  out_mat <- matrix(NA_real_, nrow = nrow(score_mat), ncol = ncol(score_mat), 
                    dimnames = dimnames(score_mat))
  
  for (s in samples) {
    idx <- which(sample_vector == s)
    if (length(idx) == 0) next
    if (length(idx) < 2) {
      out_mat[, idx] <- 0 
      next
    }
    sub_mat <- score_mat[, idx, drop = FALSE]
    sub_z <- t(scale(t(sub_mat)))
    sub_z[!is.finite(sub_z)] <- 0
    sub_z[sub_z > zlim] <- zlim
    sub_z[sub_z < -zlim] <- -zlim
    out_mat[, idx] <- sub_z
  }
  out_mat
}

score_mean_signature <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) {
    out <- rep(NA_real_, ncol(expr_mat))
    names(out) <- colnames(expr_mat)
    return(out)
  }
  out <- Matrix::colMeans(expr_mat[genes, , drop = FALSE])
  names(out) <- colnames(expr_mat)
  out
}

# Function to build the multi-row heatmap list
generate_heatmap_list <- function(score_mat, plot_meta, top_ann, column_split, column_order, is_norm = TRUE) {
  ht_list <- NULL
  for (i in 1:nrow(score_mat)) {
    row_name <- rownames(score_mat)[i]
    row_data <- score_mat[i, , drop = FALSE]
    
    if (is_norm) {
      # Z-score color scale
      lim <- max(abs(quantile(row_data, c(0.01, 0.99), na.rm = TRUE)))
      if (is.na(lim) || lim <= 0) lim <- 2.0
      col_fun <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))
    } else {
      # Raw color scale
      max_val <- quantile(row_data, 0.99, na.rm = TRUE)
      if (is.na(max_val) || max_val <= 0) max_val <- max(row_data, na.rm = TRUE)
      if (max_val <= 0) max_val <- 1
      col_fun <- colorRamp2(c(0, max_val), c("white", "firebrick3"))
    }
    
    ht_row <- Heatmap(
      row_data,
      name = row_name,
      col = col_fun,
      top_annotation = if (i == 1) top_ann else NULL,
      column_split = column_split,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_order = column_order,
      cluster_column_slices = FALSE,
      show_column_names = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 12),
      column_title = NULL,
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = row_name,
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        direction = "horizontal"
      ),
      use_raster = TRUE,
      raster_quality = 5,
      border = FALSE,
      height = unit(12, "mm")
    )
    
    if (is.null(ht_list)) {
      ht_list <- ht_row
    } else {
      ht_list <- ht_list %v% ht_row
    }
  }
  ht_list
}

generate_agg_heatmap_list <- function(agg_mat, target_mp_order, is_norm = TRUE) {
  ht_agg_list <- NULL
  for (i in 1:nrow(agg_mat)) {
    row_name <- rownames(agg_mat)[i]
    row_data <- agg_mat[, label_mp(target_mp_order), drop = FALSE][i, , drop = FALSE]
    
    if (is_norm) {
      max_val <- max(abs(row_data), na.rm = TRUE)
      if (max_val == 0) max_val <- 1
      col_fun <- colorRamp2(c(-max_val, 0, max_val), c("navy", "white", "firebrick3"))
    } else {
      max_val <- max(row_data, na.rm = TRUE)
      if (max_val <= 0) max_val <- 1
      col_fun <- colorRamp2(c(0, max_val), c("white", "firebrick3"))
    }
    
    ht_row <- Heatmap(
      row_data,
      name = paste0("Agg_", row_name),
      col = col_fun,
      top_annotation = NULL,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_side = "top",
      show_column_names = (i == 1),
      column_names_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_rot = 60,
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = row_name, 
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        direction = "horizontal"
      ),
      border = TRUE,
      height = unit(18, "mm")
    )
    
    if (is.null(ht_agg_list)) {
      ht_agg_list <- ht_row
    } else {
      ht_agg_list <- ht_agg_list %v% ht_row
    }
  }
  ht_agg_list
}

message("Loading state assignments and MP scores...")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")

common_cells <- intersect(names(state_B), rownames(mp_adj_noncc))
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

target_cells <- names(state_B)[state_B %in% target_states]
if (length(target_cells) == 0) stop("No Basal or SMG-like cells found in Auto_topmp_v2_noreg_states_B.rds")

assign_top_mp <- function(cells, state_name, mp_set) {
  if (length(cells) == 0) return(setNames(character(0), character(0)))
  mp_set <- intersect(mp_set, colnames(mp_adj_noncc))
  if (length(mp_set) == 0) stop("No MPs from ", state_name, " are present in mp_adj_noncc")
  sub_mat <- mp_adj_noncc[cells, mp_set, drop = FALSE]
  top_mp <- colnames(sub_mat)[max.col(sub_mat, ties.method = "first")]
  names(top_mp) <- cells
  top_mp
}

basal_cells <- target_cells[state_B[target_cells] == "Basal to Intestinal Metaplasia"]
smg_cells <- target_cells[state_B[target_cells] == "SMG-like Metaplasia"]

top_mp <- c(
  assign_top_mp(basal_cells, "Basal to Intestinal Metaplasia", basal_mps),
  assign_top_mp(smg_cells, "SMG-like Metaplasia", smg_mps)
)
target_cells <- names(top_mp)

top_mp_score <- vapply(
  target_cells,
  function(cell) mp_adj_noncc[cell, top_mp[cell]],
  numeric(1)
)

message("Loading epithelial Seurat object and subsetting to target cells...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]

expr_all <- get_data_layer(tmdata_all, assay = "RNA")
cc_consensus <- match_genes(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(expr_all))
if (length(cc_consensus) == 0) stop("No consensus cell-cycle genes found in the RNA data")
cc_gene_means <- Matrix::rowMeans(expr_all[cc_consensus, , drop = FALSE])
cc_top50 <- names(sort(cc_gene_means, decreasing = TRUE))[seq_len(min(50, length(cc_gene_means)))]
rm(expr_all, cc_gene_means)
gc()

target_cells <- intersect(target_cells, Cells(tmdata_all))
tmdata_all <- tmdata_all[, target_cells]
top_mp <- top_mp[target_cells]
top_mp_score <- top_mp_score[target_cells]
state_B <- state_B[target_cells]

plot_meta <- data.frame(
  cell = target_cells,
  state = as.character(state_B[target_cells]),
  top_mp = as.character(top_mp[target_cells]),
  top_mp_score = as.numeric(top_mp_score[target_cells]),
  sample = as.character(tmdata_all$orig.ident[target_cells]),
  study = as.character(tmdata_all$study[target_cells]),
  stringsAsFactors = FALSE
)

plot_meta$top_mp <- factor(plot_meta$top_mp, levels = target_mp_order)
plot_meta <- plot_meta %>%
  filter(!is.na(top_mp)) %>%
  group_by(top_mp) %>%
  arrange(desc(top_mp_score), study, sample, .by_group = TRUE) %>%
  { if (is.finite(MAX_CELLS_PER_MP)) slice_head(., n = MAX_CELLS_PER_MP) else . } %>%
  ungroup()

plot_meta$top_mp <- factor(as.character(plot_meta$top_mp), levels = target_mp_order)
plot_meta$top_mp_label <- factor(label_mp(as.character(plot_meta$top_mp)), levels = label_mp(target_mp_order))
plot_meta$state <- factor(plot_meta$state, levels = target_states)
cells_to_plot <- plot_meta$cell

tmdata_all <- tmdata_all[, cells_to_plot]
expr_mat <- get_data_layer(tmdata_all, assay = "RNA")

message("Preparing gene signatures...")
hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_alpha <- hallmark_tbl$gene_symbol[hallmark_tbl$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]
hallmark_gamma <- hallmark_tbl$gene_symbol[hallmark_tbl$gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE"]

ucell_gene_sets_input <- list(
  Squamous = c("TP63", "KRT5", "KRT14", "SOX2", "KRT4", "KRT13", "DSG3", "IVL", "SPRR1A", "KRT15"),
  Gastric_Columnar = c("FOXA1", "HNF4A", "GATA4", "KRT8", "KRT18", "VSIG1", "TFF1", "MUC5AC", "GKN1", "GKN2"),
  Intestinal_Metaplasia = c("CDX2", "MUC2", "TFF3", "REG4", "KRT20", "ATOH1", "SPINK4", "FCGBP", "VIL1", "CDH17"),
  SMG_like_Secretory = c("MUC6", "MUC5B", "TFF2", "BPIFB1", "LTF", "AQP5", "SPINK1", "ERO1B", "SOX9", "PIGR"),
  Interferon_Alpha_Response = hallmark_alpha,
  Interferon_Gamma_Response = hallmark_gamma
)

ucell_gene_sets <- lapply(ucell_gene_sets_input, match_genes, universe = rownames(expr_mat))
gene_presence <- bind_rows(lapply(names(ucell_gene_sets_input), function(sig) {
  input <- unique(toupper(ucell_gene_sets_input[[sig]]))
  matched <- ucell_gene_sets[[sig]]
  data.frame(
    signature = sig,
    n_input = length(input),
    n_present = length(matched),
    present_genes = paste(matched, collapse = ";"),
    missing_genes = paste(setdiff(input, toupper(matched)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(gene_presence, file.path(out_dir, "Auto_basal_smg_signature_gene_presence.csv"), row.names = FALSE)

ucell_gene_sets <- ucell_gene_sets[lengths(ucell_gene_sets) > 0]
if (length(ucell_gene_sets) == 0) stop("None of the requested UCell gene sets have genes present in the object")

message("Scoring marker signatures with UCell...")
tmdata_all <- AddModuleScore_UCell(tmdata_all, features = ucell_gene_sets, ncores = NCORES_UCELL, name = "")
ucell_scores <- tmdata_all@meta.data[cells_to_plot, names(ucell_gene_sets), drop = FALSE]

message("Scoring cell-cycle and Buffa hypoxia signature A from normalized RNA data...")
cc_score <- score_mean_signature(expr_mat, cc_top50)

buffa_hypoxia_A <- c(
  "ACOT7", "ADM", "AK4", "ALDOA", "ANKRD37", "ANLN", "BNIP3", "MRGBP", "CA9", "CDKN3",
  "CHCHD2", "CORO1C", "CTSV", "DDIT4", "ENO1", "ESRP1", "GAPDH", "GPI", "HILPDA", "HK2",
  "KIF20A", "KIF4A", "LDHA", "LRRC42", "MAD2L2", "MAP7D1", "MCTS1", "MIF", "MRPL13", "MRPL15",
  "MRPS17", "NDRG1", "ZNF384", "P4HA1", "PFKP", "PGAM1", "PGK1", "PSMA7", "PSRC1", "SEC61G",
  "SHCBP1", "SLC16A1", "SLC25A32", "SLC2A1", "TPI1", "TUBA1B", "TUBA1C", "TUBB6", "UTP11",
  "VEGFA", "YKT6"
)
hypoxia_genes <- match_genes(buffa_hypoxia_A, rownames(expr_mat))
hypoxia_score <- score_mean_signature(expr_mat, hypoxia_genes)

score_raw <- cbind(
  as.matrix(ucell_scores),
  Cell_Cycle_Score = cc_score[cells_to_plot],
  Buffa_Hypoxia_A_Mean = hypoxia_score[cells_to_plot]
)
score_raw <- t(score_raw)

row_display <- c(
  Squamous = "Squamous",
  Gastric_Columnar = "Gastric Columnar",
  Intestinal_Metaplasia = "Intestinal Metaplasia",
  SMG_like_Secretory = "SMG-like Secretory",
  Interferon_Alpha_Response = "Interferon Alpha Response",
  Interferon_Gamma_Response = "Interferon Gamma Response",
  Cell_Cycle_Score = "Cell-cycle score",
  Buffa_Hypoxia_A_Mean = "Buffa hypoxia A mean"
)
score_raw <- score_raw[names(row_display), , drop = FALSE]
rownames(score_raw) <- row_display[rownames(score_raw)]

message("Normalizing rows within each sample...")
score_norm <- normalize_within_sample(score_raw, plot_meta$sample)

plot_meta <- plot_meta[match(colnames(score_norm), plot_meta$cell), , drop = FALSE]
if (!identical(plot_meta$cell, colnames(score_norm))) stop("Column metadata and heatmap scores are out of sync")

study_vals <- unique(as.character(plot_meta$study))
study_cols <- setNames(Seurat::DiscretePalette(length(study_vals), palette = "polychrome"), study_vals)

top_ann <- HeatmapAnnotation(
  State = plot_meta$state,
  Top_MP = plot_meta$top_mp_label,
  Study = plot_meta$study,
  col = list(
    State = state_cols,
    Top_MP = setNames(mp_cols[target_mp_order], label_mp(target_mp_order)),
    Study = study_cols
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(
    State = list(title_gp = gpar(fontsize = 11, fontface = "bold"), labels_gp = gpar(fontsize = 10)),
    Top_MP = list(title_gp = gpar(fontsize = 11, fontface = "bold"), labels_gp = gpar(fontsize = 10)),
    Study = list(title_gp = gpar(fontsize = 11, fontface = "bold"), labels_gp = gpar(fontsize = 10))
  ),
  show_legend = TRUE,
  na_col = "white"
)

column_split <- factor(plot_meta$top_mp_label, levels = label_mp(target_mp_order))

## Manual column order calculation for unsupervised clustering within each MP group
message("Computing hierarchical clustering within each MP group using all signatures...")
column_split_vec <- as.character(column_split)
splits <- levels(column_split)
final_column_order <- integer(0)

# We use score_norm for clustering on all pages to keep cell order consistent
for (s in splits) {
  idx_in_slice <- which(column_split_vec == s)
  if (length(idx_in_slice) > 1) {
    # Cluster within slice using all 8 signatures from score_norm
    # dist() calculates distance between columns of score_norm (rows of transposed)
    d <- dist(t(score_norm[, idx_in_slice, drop = FALSE]))
    hc_sub <- hclust(d, method = "ward.D2")
    final_column_order <- c(final_column_order, idx_in_slice[hc_sub$order])
  } else {
    final_column_order <- c(final_column_order, idx_in_slice)
  }
}

####################
# Bubble plot helper function — nature-figure style (UPGRADED)
# Shows mean signature score (color) and detection rate (size) per MP group
# Each signature now has its own relative color scale for maximum clarity
####################
generate_bubble_plot <- function(score_mat, plot_meta, target_mp_order, mp_cols, mp_descriptions, is_norm = TRUE) {
  library(ggplot2)
  library(dplyr)
  
  score_df <- as.data.frame(t(score_mat))
  score_df$top_mp <- plot_meta$top_mp
  
  sig_names <- rownames(score_mat)
  
  # Compute mean score and detection rate per MP group
  bubble_data <- do.call(rbind, lapply(sig_names, function(sig) {
    do.call(rbind, lapply(target_mp_order, function(mp) {
      mp_lab <- paste0(mp, " ", mp_descriptions[mp])
      idx <- which(as.character(score_df$top_mp) == mp)
      if (length(idx) == 0) return(NULL)
      vals <- score_df[idx, sig]
      
      if (is_norm) {
        # For z-scored data, detection = fraction of cells with |z| > 0.5
        det_rate <- mean(abs(vals) > 0.5, na.rm = TRUE)
      } else {
        # For raw data, detection = fraction of cells with score > 0
        det_rate <- mean(vals > 0, na.rm = TRUE)
      }
      
      data.frame(
        signature = sig,
        mp_group = mp_lab,
        mp_id = mp,
        mean_score = mean(vals, na.rm = TRUE),
        detection_rate = det_rate,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Perform per-signature normalization of mean_score for the color mapping
  # This makes each row have its own relative scale, highlighting differences within the signature
  bubble_data <- bubble_data %>%
    group_by(signature) %>%
    mutate(
      mean_score_scaled = if(is_norm) {
        # Scale to [-1, 1] based on max absolute value within signature
        mx <- max(abs(mean_score), na.rm = TRUE)
        if(mx > 1e-6) mean_score / mx else 0
      } else {
        # Scale to [0, 1]
        mn <- min(mean_score, na.rm = TRUE)
        mx <- max(mean_score, na.rm = TRUE)
        if(mx > mn + 1e-6) (mean_score - mn) / (mx - mn) else 0
      }
    ) %>%
    ungroup()
  
  # Factor ordering
  mp_order_labels <- paste0(target_mp_order, " ", mp_descriptions[target_mp_order])
  bubble_data$mp_group <- factor(bubble_data$mp_group, levels = mp_order_labels)
  bubble_data$signature <- factor(bubble_data$signature, levels = rev(sig_names))
  
  # Color scale — using the SCALED score
  if (is_norm) {
    fill_scale <- scale_color_gradient2(
      low = "navy", mid = "white", high = "firebrick3",
      midpoint = 0, limits = c(-1, 1),
      name = "Relative\nZ-score"
    )
  } else {
    fill_scale <- scale_color_gradient(
      low = "grey95", high = "firebrick3",
      limits = c(0, 1),
      name = "Relative\nScore"
    )
  }
  
  p <- ggplot(bubble_data, aes(x = mp_group, y = signature)) +
    geom_point(aes(size = detection_rate, color = mean_score_scaled), stroke = 0.5) +
    fill_scale +
    scale_size_continuous(
      range = c(2, 16), # Larger dots
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      labels = c("25%", "50%", "75%", "100%"),
      name = "Detection\nRate"
    ) +
    # Basal vs SMG separator
    geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    annotate("text", x = 3, y = length(sig_names) + 0.8, label = "Basal to Intestinal Metaplasia",
             hjust = 0.5, size = 5.5, fontface = "bold", color = "#4DAF4A") +
    annotate("text", x = 6.5, y = length(sig_names) + 0.8, label = "SMG-like Metaplasia",
             hjust = 0.5, size = 5.5, fontface = "bold", color = "#FF7F00") +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 14) + # Larger base font
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right",
      legend.margin = margin(l = 10),
      panel.grid = element_blank(),
      plot.margin = margin(t = 40, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2, barwidth = 1.2, barheight = 6)
    )
  
  p
}

pdf(file.path(out_dir, "Auto_basal_smg_mp_signature_heatmap.pdf"), width = 18, height = 12, useDingbats = FALSE)

# Page 1: Per-cell Heatmap (Normalized within sample)
message("Generating Page 1: Per-cell Normalized Heatmap...")
ht_list_norm <- generate_heatmap_list(score_norm, plot_meta, top_ann, column_split, final_column_order, is_norm = TRUE)
draw(ht_list_norm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "right", ht_gap = unit(1.5, "mm"))

# Page 2: Combined Heatmap (Normalized)
message("Generating Page 2: Aggregated Normalized Heatmap...")
score_df_norm <- as.data.frame(t(score_norm))
score_df_norm$top_mp_label <- plot_meta$top_mp_label
agg_scores_norm <- score_df_norm %>%
  group_by(top_mp_label) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()
rownames(agg_scores_norm) <- agg_scores_norm$top_mp_label
agg_mat_norm <- t(as.matrix(agg_scores_norm[, -1]))

top_ann_agg <- HeatmapAnnotation(
  Top_MP = label_mp(target_mp_order),
  col = list(
    Top_MP = setNames(mp_cols[target_mp_order], label_mp(target_mp_order))
  ),
  show_legend = FALSE
)

ht_agg_norm <- generate_agg_heatmap_list(agg_mat_norm, target_mp_order, is_norm = TRUE)
draw(ht_agg_norm, merge_legend = TRUE, heatmap_legend_side = "bottom", ht_gap = unit(3, "mm"),
     padding = unit(c(20, 2, 2, 2), "mm"))

####################
# Page 2b: Bubble plot (Normalized) — combined signature expression per MP group
####################
message("Generating Page 2b: Bubble Plot (Normalized)...")
p_bubble_norm <- generate_bubble_plot(score_norm, plot_meta, target_mp_order, mp_cols, mp_descriptions, is_norm = TRUE)
print(p_bubble_norm)

# Page 3: Per-cell Heatmap (Raw scores)
message("Generating Page 3: Per-cell Raw Heatmap...")
ht_list_raw <- generate_heatmap_list(score_raw, plot_meta, top_ann, column_split, final_column_order, is_norm = FALSE)
draw(ht_list_raw, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "right", ht_gap = unit(2, "mm"))

# Page 4: Combined Heatmap (Raw)
message("Generating Page 4: Aggregated Raw Heatmap...")
score_df_raw <- as.data.frame(t(score_raw))
score_df_raw$top_mp_label <- plot_meta$top_mp_label
agg_scores_raw <- score_df_raw %>%
  group_by(top_mp_label) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()
rownames(agg_scores_raw) <- agg_scores_raw$top_mp_label
agg_mat_raw <- t(as.matrix(agg_scores_raw[, -1]))

ht_agg_raw <- generate_agg_heatmap_list(agg_mat_raw, target_mp_order, is_norm = FALSE)
draw(ht_agg_raw, merge_legend = TRUE, heatmap_legend_side = "bottom", ht_gap = unit(3, "mm"),
     padding = unit(c(20, 2, 2, 2), "mm"))

####################
# Page 4b: Bubble plot (Raw) — combined signature expression per MP group
####################
message("Generating Page 4b: Bubble Plot (Raw)...")
p_bubble_raw <- generate_bubble_plot(score_raw, plot_meta, target_mp_order, mp_cols, mp_descriptions, is_norm = FALSE)
print(p_bubble_raw)

dev.off()

####################
# SVG export of bubble plots (nature-figure primary vector format)
####################
message("Exporting bubble plot SVGs...")
tryCatch({
  if (requireNamespace("svglite", quietly = TRUE)) {
    # Using more compact aspect ratio (8x6) for SVGs to make fonts relatively larger
    svglite::svglite(file.path(out_dir, "Auto_basal_smg_bubble_normalized.svg"), width = 8, height = 6)
    print(p_bubble_norm)
    dev.off()
    svglite::svglite(file.path(out_dir, "Auto_basal_smg_bubble_raw.svg"), width = 8, height = 6)
    print(p_bubble_raw)
    dev.off()
    message("SVG exports saved.")
  } else {
    message("svglite not available, skipping SVG export.")
  }
}, error = function(e) message("SVG export skipped: ", e$message))

score_long <- as.data.frame(t(score_norm), check.names = FALSE)
score_long$cell <- rownames(score_long)
summary_df <- plot_meta %>%
  select(cell, state, top_mp, top_mp_label, study, sample) %>%
  left_join(score_long, by = "cell") %>%
  group_by(state, top_mp, top_mp_label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    n_studies = n_distinct(study),
    across(all_of(rownames(score_norm)), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write.csv(summary_df, file.path(out_dir, "Auto_basal_smg_mp_signature_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    score_raw = score_raw,
    score_norm = score_norm,
    cell_metadata = plot_meta,
    gene_presence = gene_presence,
    ucell_gene_sets = ucell_gene_sets,
    cell_cycle_top50 = cc_top50,
    hypoxia_genes = hypoxia_genes,
    target_mp_order = target_mp_order,
    agg_scores = agg_scores_norm,
    agg_scores_raw = agg_scores_raw
  ),
  file.path(out_dir, "Auto_basal_smg_mp_signature_heatmap_data.rds")
)

message("Saved Basal/SMG MP signature heatmap outputs to: ", file.path(getwd(), out_dir))
