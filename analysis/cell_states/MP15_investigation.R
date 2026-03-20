####################
# Auto_MP15_investigation.R
#
# Investigate whether epithelial cells contributing to MP15 (T_NK_infiltration)
# express non-epithelial cell type markers — indicating contamination or doublets.
#
# Two approaches (both use z-normalised UCell scores):
#   A) Top-3 samples contributing to MP15 → expression heatmap, column split
#      by MP15 binary activity (high vs low z-score)
#   B) All cells where MP15 z-score is "high" → expression heatmap
#
# Quick check: QC heatmap for Alcindor_2025_SRR27335939 (all cell types)
#
# Z-normalisation: center within sample, scale by within-study SD
# (as in mp_external_scoring.R lines 177-223)
#
# Plotting style: qc_heatmap.R (per-celltype column grobs, gridExtra layout)
# Run interactively in dmtcp environment.
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
####################
# Added dependencies for MP15 CNV/UMAP plots
library(ggplot2)
####################

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ============================================================================
# 1. Load data
# ============================================================================
cat("Loading geneNMF metaprograms (nMP=19)...\n")
geneNMF.metaprograms <- readRDS(
  "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)

cat("Loading EAC_Ref_epi.rds...\n")
tmdata_all <- readRDS("EAC_Ref_epi.rds")

cat("Loading UCell scores (nMP19 filtered)...\n")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")

cat(sprintf("Seurat object: %d cells | UCell scores: %d cells x %d MPs\n",
            ncol(tmdata_all), nrow(ucell_scores), ncol(ucell_scores)))

# ============================================================================
# 2. Updated MP annotations
# ============================================================================
mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

# ============================================================================
# 3. Silhouette filtering — retain only good MPs
# ============================================================================
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
na_mps <- which(is.na(geneNMF.metaprograms$metaprograms.metrics$silhouette))
na_mp_names <- paste0("MP", na_mps)
remove_mps <- unique(c(bad_mp_names, na_mp_names))

all_mp_names <- names(geneNMF.metaprograms$metaprograms.genes)
retained_mp_names <- all_mp_names[!all_mp_names %in% remove_mps]
cat(sprintf("Removed %d MPs (silhouette < 0 or NA): %s\n",
            length(remove_mps), paste(remove_mps, collapse = ", ")))
cat(sprintf("Retained %d MPs: %s\n", length(retained_mp_names),
            paste(retained_mp_names, collapse = ", ")))

# Filter UCell scores to retained MPs only
ucell_retained <- ucell_scores[, intersect(retained_mp_names, colnames(ucell_scores)),
                               drop = FALSE]

# Align cells
common_cells <- intersect(rownames(ucell_retained), colnames(tmdata_all))
ucell_retained <- ucell_retained[common_cells, , drop = FALSE]
tmdata_all <- tmdata_all[, common_cells]
cat(sprintf("Aligned: %d common cells\n", length(common_cells)))

# ============================================================================
# 4. Z-normalise UCell scores: center within sample, scale by within-study SD
#    (following mp_external_scoring.R lines 177-223)
# ============================================================================
cat("Z-normalising UCell scores (center per sample, scale by within-study SD)...\n")

final_mps <- colnames(ucell_retained)
sample_var <- as.character(tmdata_all$orig.ident)
names(sample_var) <- colnames(tmdata_all)
study_var <- as.character(tmdata_all$study)
names(study_var) <- colnames(tmdata_all)

clust_df <- as.data.frame(ucell_retained)
clust_df$.cell   <- rownames(ucell_retained)
clust_df$.sample <- sample_var[clust_df$.cell]
clust_df$.study  <- study_var[clust_df$.cell]

# Per-study SD for each MP
study_sd <- clust_df %>%
  group_by(.study) %>%
  summarise(across(all_of(final_mps), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames(".study") %>%
  as.matrix()
study_sd[is.na(study_sd) | study_sd == 0] <- 1

# Center within sample
clust_centered <- clust_df %>%
  group_by(.sample) %>%
  mutate(across(all_of(final_mps), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Scale by within-study SD
mp_adj <- as.matrix(clust_centered[, final_mps])
rownames(mp_adj) <- clust_centered$.cell

for (mp in final_mps) {
  cell_studies <- clust_centered$.study
  mp_adj[, mp] <- mp_adj[, mp] / study_sd[cell_studies, mp]
}
mp_adj[!is.finite(mp_adj)] <- 0

cat(sprintf("Z-normalised matrix: %d cells x %d MPs\n",
            nrow(mp_adj), ncol(mp_adj)))

# MP15 z-scores
mp15_z <- mp_adj[, "MP15"]
cat(sprintf("MP15 z-score range: [%.2f, %.2f], mean: %.4f, sd: %.4f\n",
            min(mp15_z), max(mp15_z), mean(mp15_z), sd(mp15_z)))

# Define MP15 "high" threshold: z > 1 (1 SD above mean, conservative)
MP15_Z_THRESHOLD <- 1
mp15_high <- mp15_z > MP15_Z_THRESHOLD
cat(sprintf("MP15 high (z > %.1f): %d cells (%.1f%%)\n",
            MP15_Z_THRESHOLD, sum(mp15_high),
            100 * sum(mp15_high) / length(mp15_high)))

# ============================================================================
# 5. Identify which samples contributed NMF programs to MP15
# ============================================================================
programs_clusters <- geneNMF.metaprograms$programs.clusters
programs_clusters <- programs_clusters[!is.na(programs_clusters)]

mp15_programs <- names(programs_clusters)[programs_clusters == 15]
cat(sprintf("\nMP15 programs (%d total):\n", length(mp15_programs)))

mp15_samples <- mp15_programs
mp15_samples <- sub("\\.[0-9]+\\.k[0-9]+$", "", mp15_samples)
mp15_samples <- sub("\\.k[0-9]+\\.[0-9]+$", "", mp15_samples)
mp15_samples_unique <- unique(mp15_samples)

cat(sprintf("MP15 derives from %d unique samples:\n", length(mp15_samples_unique)))
for (s in mp15_samples_unique) cat(sprintf("  - %s\n", s))

# Rank by cell count
cell_counts <- table(tmdata_all$orig.ident)
mp15_cell_counts <- cell_counts[intersect(mp15_samples_unique, names(cell_counts))]
mp15_cell_counts <- sort(mp15_cell_counts, decreasing = TRUE)

cat("\nMP15 samples ranked by cell count:\n")
for (i in seq_along(mp15_cell_counts)) {
  cat(sprintf("  %d. %s : %d cells\n", i, names(mp15_cell_counts)[i],
              mp15_cell_counts[i]))
}

top3_samples <- names(mp15_cell_counts)[1:min(3, length(mp15_cell_counts))]
cat(sprintf("\nTop 3 samples selected: %s\n", paste(top3_samples, collapse = ", ")))

# ============================================================================
# 6. Define cell type markers (from Expr_filtering.R)
# ============================================================================
fibroblast  <- c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN")
macrophage  <- c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68")
mast        <- c("MS4A2", "CPA3", "TPSB2", "TPSAB1")
epithelial  <- c("KRT7", "MUC1", "KRT19", "EPCAM")
t.cell      <- c("CD3E", "CD3D", "CD2", "CD3G")
b.cell      <- c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1")
nk.cell     <- c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1")
plasma      <- c("MZB1", "JCHAIN", "DERL3")
dendritic   <- c("CLEC10A", "CCR7", "CD86")
endothelial <- c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5")

housekeeping <- c("ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16",
                   "RPS18", "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18")

markers <- c(fibroblast, macrophage, mast, epithelial, t.cell, nk.cell,
             b.cell, plasma, dendritic, endothelial, housekeeping)

markers_list <- list(
  b.cell      = b.cell,
  dendritic   = dendritic,
  endothelial = endothelial,
  epithelial  = epithelial,
  fibroblast  = fibroblast,
  macrophage  = macrophage,
  mast        = mast,
  nk.cell     = nk.cell,
  plasma      = plasma,
  t.cell      = t.cell,
  housekeeping = housekeeping
)

# ============================================================================
# 7. QC-style heatmap function (from qc_heatmap.R)
# ============================================================================
plot_qc_heatmap <- function(tmdata, sampleid, identity, reorder = FALSE,
                            ct_reorder = NULL) {

  expr_data <- FetchData(tmdata, vars = markers, layer = "data")
  missing_markers <- setdiff(markers, colnames(expr_data))
  for (gene in missing_markers) {
    expr_data[[gene]] <- 0
  }
  expr_data <- expr_data[, markers]
  local_markers_list <- markers_list
  for (name in names(local_markers_list)) {
    local_markers_list[[name]] <- local_markers_list[[name]][
      local_markers_list[[name]] %in% colnames(expr_data)
    ]
  }
  expr_data$celltype <- tmdata@meta.data[[identity]]
  expr_data$celltype <- sapply(expr_data$celltype, function(x) gsub(" ", "\n", x))

  expr_data <- expr_data[order(expr_data$celltype), ]
  hk_avg <- rowMeans(expr_data[, colnames(expr_data) %in%
    local_markers_list$housekeeping, drop = FALSE])
  hk_avg <- matrix(hk_avg, nrow = 1, dimnames = list("avg_hk", names(hk_avg)))
  nCounts <- tmdata$nFeature_RNA[rownames(expr_data)]
  nCounts <- matrix(nCounts, nrow = 1, dimnames = list("nGenes", names(nCounts)))

  heatplot <- t(as.matrix(expr_data[, 1:(ncol(expr_data) - 1)]))

  custom_colors <- colorRamp2(
    c(0, round(0.6 * max(heatplot), 1), ceiling(max(heatplot))),
    c("#D0D0D0", "red4", "red4")
  )

  heatmap_grobs <- list()
  len <- length(unique(expr_data$celltype))

  if (!reorder || is.null(ct_reorder)) {
    desired_order <- unique(expr_data$celltype)
  } else {
    desired_order <- ct_reorder
  }
  present_celltypes <- intersect(desired_order, unique(expr_data$celltype))

  min_ngenes_val <- 200
  min_hk_expr_val <- 0
  max_mt_val <- 15

  for (i in seq_along(local_markers_list)) {
    marker <- local_markers_list[[i]]
    marker_names <- rownames(heatplot[marker, , drop = FALSE])
    temp <- list()
    for (j in seq_along(present_celltypes)) {
      celltype <- present_celltypes[j]
      cells <- rownames(expr_data)[which(expr_data$celltype == celltype)]

      ht <- Heatmap(
        heatplot[marker, cells, drop = FALSE],
        col = custom_colors,
        show_column_names = FALSE, show_row_names = FALSE,
        row_names_side = "left", row_names_gp = gpar(fontsize = 40),
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE, use_raster = FALSE
      )
      ht_grob <- grid.grabExpr(
        draw(ht, newpage = FALSE, padding = unit(c(2, 1, 2, 1), "mm"))
      )
      temp[[j]] <- ht_grob
    }

    gene_label_grobs <- lapply(marker_names, function(name) {
      textGrob(label = name, x = unit(1, "npc"), just = "right",
               gp = gpar(fontsize = 40))
    })
    gene_label_col <- arrangeGrob(grobs = gene_label_grobs, ncol = 1)
    temp <- c(list(gene_label_col), temp)

    if (j == length(present_celltypes)) {
      text_grob <- textGrob(
        names(local_markers_list)[i],
        gp = gpar(fontsize = 50, fontface = "bold")
      )
      rect_grob <- rectGrob(gp = gpar(fill = "grey", col = NA))
      merged_grob <- gTree(children = gList(rect_grob, text_grob))
      temp[[j + 2]] <- merged_grob
      combined_grob <- do.call(
        arrangeGrob,
        c(temp, list(ncol = len + 2,
                     widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
      )
      heatmap_grobs[[length(heatmap_grobs) + 1]] <- combined_grob
    }
  }

  stats_grobs <- list()
  for (i in seq_along(list(nCounts, hk_avg))) {
    data <- list(nCounts, hk_avg)[[i]]
    data_names <- rownames(data)
    temp <- list()
    for (j in seq_along(present_celltypes)) {
      celltype <- present_celltypes[j]
      cells <- rownames(expr_data)[which(expr_data$celltype == celltype)]

      if (i == 1) {
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name = "Number of\ngenes detected",
          col = colorRamp2(
            c(min_ngenes_val, round(0.7 * max(nCounts), -3),
              ceiling(max(nCounts) / 1000) * 1000),
            c("#D0D0D0", "blue3", "blue3")
          ),
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = FALSE, row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE, show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      } else {
        rounded_value <- ceiling(max(hk_avg) * 10) / 10
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name = "Average\nhousekeeping\nexpression",
          col = colorRamp2(
            c(0, max(min_hk_expr_val, 0.01), rounded_value),
            c("#D0D0D0", "#A0B8E6", "blue3")
          ),
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = FALSE, row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE, show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      }
      ht_grob <- grid.grabExpr(
        draw(ht_add, newpage = FALSE, padding = unit(c(6, 1.5, 6, 1.5), "mm"))
      )
      temp[[j]] <- ht_grob
    }

    gene_label_grobs <- lapply(data_names, function(name) {
      textGrob(label = name, just = "center", gp = gpar(fontsize = 40))
    })
    gene_label_col <- arrangeGrob(grobs = gene_label_grobs, ncol = 1)
    temp <- c(list(gene_label_col), temp)

    if (j == length(present_celltypes)) {
      text_grob <- textGrob(
        c("Number\nof genes", "Average\nhousekeeping\nexpression")[i],
        gp = gpar(fontsize = 40, fontface = "bold")
      )
      rect_grob <- rectGrob(gp = gpar(fill = "white", col = NA))
      merged_grob <- gTree(children = gList(rect_grob, text_grob))
      temp[[j + 2]] <- merged_grob
      combined_grob <- do.call(
        arrangeGrob,
        c(temp, list(ncol = len + 2,
                     widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
      )
      stats_grobs[[length(stats_grobs) + 1]] <- combined_grob
    }
  }

  celltype_labels <- c("", as.character(present_celltypes), "", "")
  text_grobs <- lapply(celltype_labels, function(ct) {
    textGrob(ct, gp = gpar(fontsize = 50, fontface = "bold"),
             just = "center", rot = 0)
  })
  text_grob_row <- do.call(
    arrangeGrob,
    c(text_grobs, list(ncol = len + 2,
                       widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  celltype_num <- c("",
                    as.numeric(table(expr_data$celltype)[present_celltypes]),
                    "Markers", "")
  num_grobs <- lapply(celltype_num, function(x) {
    textGrob(x, gp = gpar(fontsize = 40, fontface = "bold"), just = "center")
  })
  num_grob_row <- do.call(
    arrangeGrob,
    c(num_grobs, list(ncol = len + 2,
                      widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  title_grob <- textGrob(
    paste0("Mito DNA < ", max_mt_val,
           "  &&  Min genes > ", min_ngenes_val,
           "  &&  Min HK expr > ", min_hk_expr_val),
    gp = gpar(fontsize = 40), just = "center"
  )

  legend_obj1 <- Legend(
    title = "\nExpression\nlevels (E)\n",
    at = pretty(c(0, ceiling(max(heatplot))), n = 5),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"), legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"), grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(
      c(0, round(0.6 * max(heatplot), 1), ceiling(max(heatplot))),
      c("#D0D0D0", "red4", "red4")
    )
  )

  legend_obj2 <- Legend(
    title = "\nNumber\nof genes\n",
    at = round(seq(min_ngenes_val,
                   ceiling(max(nCounts) / 1000) * 1000,
                   length.out = 4), -3),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"), legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"), grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(
      c(min_ngenes_val, round(0.7 * max(nCounts), -3),
        ceiling(max(nCounts) / 1000) * 1000),
      c("#D0D0D0", "blue3", "blue3")
    )
  )

  rounded_value <- ceiling(max(hk_avg) * 10) / 10
  legend_obj3 <- Legend(
    title = "\nAverage\nhousekeeping\nexpression\n",
    at = seq(floor(min_hk_expr_val), rounded_value, length.out = 4) %>% round(0),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"), legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"), grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(
      c(0, max(min_hk_expr_val, 0.01), rounded_value),
      c("#D0D0D0", "#A0B8E6", "blue3")
    )
  )

  main_content <- arrangeGrob(
    grobs = c(list(title_grob), list(text_grob_row), list(num_grob_row),
              heatmap_grobs, stats_grobs),
    ncol = 1,
    heights = c(2.5, 2, 1.5, rep(4, length(heatmap_grobs)), 4, 4)
  )

  legend_grob <- textGrob(
    label = sampleid, gp = gpar(fontsize = 35, fontface = "bold")
  )
  legend_grob0 <- textGrob(
    label = paste0("Total cell count: ",
                   sum(expr_data$celltype %in% present_celltypes)),
    gp = gpar(fontsize = 35)
  )
  legend_grob1 <- grid.grabExpr(draw(legend_obj1))
  legend_grob2 <- grid.grabExpr(draw(legend_obj2))
  legend_grob3 <- grid.grabExpr(draw(legend_obj3))
  legend_column <- arrangeGrob(
    legend_grob, legend_grob0, legend_grob1, legend_grob2, legend_grob3,
    ncol = 1, heights = c(0.3, 0.2, 2, 2, 2)
  )

  final_layout <- grid.arrange(
    main_content, legend_column,
    ncol = 2, widths = c(9.3, 1)
  )

  return(final_layout)
}

# ============================================================================
# 8. APPROACH A: Top-3 MP15 samples — column split by MP15 z-score (high/low)
# ============================================================================
cat("\n========================================\n")
cat("APPROACH A: Top-3 MP15 samples, split by MP15 activity\n")
cat("========================================\n")

for (samp in top3_samples) {
  cat(sprintf("\nProcessing sample: %s\n", samp))

  samp_cells <- common_cells[tmdata_all$orig.ident[common_cells] == samp]
  cat(sprintf("  Cells in this sample: %d\n", length(samp_cells)))

  if (length(samp_cells) < 10) {
    cat(sprintf("  Skipping: only %d cells\n", length(samp_cells)))
    next
  }

  # Subsample if too many
  set.seed(42)
  if (length(samp_cells) > 10000) {
    samp_cells <- sample(samp_cells, 10000)
    cat(sprintf("  Subsampled to %d cells\n", length(samp_cells)))
  }

  tmdata_samp <- subset(tmdata_all, cells = samp_cells)

  # Assign MP15 binary activity as the "celltype" for column split
  samp_mp15_z <- mp_adj[samp_cells, "MP15"]
  tmdata_samp$MP15_activity <- ifelse(samp_mp15_z > MP15_Z_THRESHOLD,
                                       "MP15_high", "MP15_low")

  cat(sprintf("  MP15 high: %d | MP15 low: %d\n",
              sum(tmdata_samp$MP15_activity == "MP15_high"),
              sum(tmdata_samp$MP15_activity == "MP15_low")))

  safe_name <- gsub("[^A-Za-z0-9_]", "_", samp)

  png(paste0("Auto_MP15_top3_", safe_name, ".png"),
      width = 80, height = 50, units = "in", res = 400)
  plot_qc_heatmap(tmdata_samp, sampleid = samp,
                  identity = "MP15_activity",
                  reorder = TRUE,
                  ct_reorder = c("MP15_high", "MP15_low"))
  dev.off()
  cat(sprintf("  Saved: Auto_MP15_top3_%s.png\n", safe_name))
}

####################
# 8B. Top-3 MP15 samples — CNV + UMAP comparisons (high vs low)
####################
cat("\n========================================\n")
cat("APPROACH A2: Top-3 MP15 samples — CNV & UMAP comparisons\n")
cat("========================================\n")

####################
# Helper: compute CNA signal from infercna outs
compute_cna_signal <- function(outs_mat, use_quantile = 0.9) {
  infercna::cnaSignal(outs_mat, gene.quantile = use_quantile)
}
####################

####################
# Helper: build binned CNV matrix (genes x cells) using 100-gene bins
bin_cnv_matrix <- function(outs_mat, gene_order, bin_size = 100L) {
  chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
  go <- gene_order %>%
    filter(gene_id %in% rownames(outs_mat), chromosome %in% chrom_levels) %>%
    mutate(chromosome = factor(chromosome, levels = chrom_levels)) %>%
    arrange(chromosome, start)
  outs_mat <- outs_mat[go$gene_id, , drop = FALSE]
  go <- go %>%
    group_by(chromosome) %>%
    mutate(
      g_rank = row_number(),
      bin_in_chr = ((g_rank - 1L) %/% bin_size) + 1L,
      bin_key = paste(chromosome, bin_in_chr, sep = "_")
    ) %>%
    ungroup()
  ordered_bin_keys <- unique(go$bin_key)
  bins_idx <- split(seq_len(nrow(go)), factor(go$bin_key, levels = ordered_bin_keys))
  binned_mat <- do.call(rbind, lapply(bins_idx, function(ix) colMeans(outs_mat[ix, , drop = FALSE])))
  rownames(binned_mat) <- names(bins_idx)
  list(binned_mat = binned_mat, go = go)
}
####################

####################
# Load gene order (hg38) for CNV binning
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
if (!file.exists(gene_order_path)) {
  stop("Gene order file not found: ", gene_order_path)
}
gene_order <- read.table(
  gene_order_path,
  header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")
)
####################

####################
# Prepare summary collector
summary_rows <- list()
####################

for (samp in top3_samples) {
  cat(sprintf("\n[CNV/UMAP] Processing sample: %s\n", samp))

  samp_cells <- common_cells[tmdata_all$orig.ident[common_cells] == samp]
  if (length(samp_cells) < 20) {
    cat(sprintf("  Skipping CNV/UMAP: too few cells (%d)\n", length(samp_cells)))
    next
  }

  ####################
  # MP15 high/low within sample
  samp_mp15_z <- mp_adj[samp_cells, "MP15"]
  mp15_group <- ifelse(samp_mp15_z > MP15_Z_THRESHOLD, "MP15_high", "MP15_low")
  names(mp15_group) <- samp_cells
  ####################

  ####################
  # CNV: load infercna outs + CNA signal distribution
  ####################
  outs_path <- file.path("by_samples", samp, paste0(samp, "_outs.rds"))
  if (!file.exists(outs_path)) {
    cat(sprintf("  No infercna outs found: %s (skipping CNV plots)\n", outs_path))
  } else {
    outs <- readRDS(outs_path)
    common_cnv_cells <- intersect(colnames(outs), samp_cells)
    if (length(common_cnv_cells) < 20) {
      cat("  Too few cells with CNV outs for this sample.\n")
    } else {
      outs <- outs[, common_cnv_cells, drop = FALSE]

      ####################
      # CNA signal
      cna_signal <- compute_cna_signal(outs, use_quantile = 0.9)
      cna_df <- data.frame(
        cell = names(cna_signal),
        cna_signal = as.numeric(cna_signal),
        mp15_group = mp15_group[names(cna_signal)],
        stringsAsFactors = FALSE
      )
      cna_df <- cna_df[!is.na(cna_df$mp15_group), ]

      ####################
      # Save CNA signal distribution plot
      safe_name <- gsub("[^A-Za-z0-9_]", "_", samp)
      p_cna <- ggplot(cna_df, aes(x = mp15_group, y = cna_signal, fill = mp15_group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
        scale_fill_manual(values = c(MP15_high = "#D73027", MP15_low = "#4575B4")) +
        labs(title = paste0("CNA signal: ", samp), x = NULL, y = "CNA signal (mean sq.)") +
        theme_classic() +
        theme(legend.position = "none")

      ggsave(filename = paste0("Auto_MP15_CNA_signal_", safe_name, ".pdf"),
             plot = p_cna, width = 6, height = 4)

      ####################
      # Binned CNV heatmap (subset to 2000 cells max)
      set.seed(42)
      if (ncol(outs) > 2000) {
        keep_cells <- sample(colnames(outs), 2000)
        outs_heat <- outs[, keep_cells, drop = FALSE]
      } else {
        outs_heat <- outs
      }

      bin_res <- bin_cnv_matrix(outs_heat, gene_order, bin_size = 100L)
      binned_mat <- bin_res$binned_mat
      binned_mat <- binned_mat[, order(mp15_group[colnames(binned_mat)])]

      ####################
      # Heatmap annotation by MP15 group
      row_ann <- rowAnnotation(
        MP15 = mp15_group[rownames(t(binned_mat))],
        col = list(MP15 = c(MP15_high = "#D73027", MP15_low = "#4575B4")),
        show_annotation_name = TRUE
      )

      col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

      pdf(paste0("Auto_MP15_CNV_heatmap_", safe_name, ".pdf"), width = 8, height = 6)
      draw(
        Heatmap(
          t(binned_mat),
          name = "CNV",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          left_annotation = row_ann,
          use_raster = TRUE
        )
      )
      dev.off()

      cat(sprintf("  Saved CNV plots for %s\n", samp))

      ####################
      # Summary row
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        sample = samp,
        n_cells_cnv = ncol(outs),
        n_cells_heat = ncol(outs_heat),
        cna_signal_median_high = median(cna_df$cna_signal[cna_df$mp15_group == "MP15_high"], na.rm = TRUE),
        cna_signal_median_low  = median(cna_df$cna_signal[cna_df$mp15_group == "MP15_low"], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  ####################
  # UMAP: re-run per-sample and color by MP15 high/low
  ####################
  sub_obj <- subset(tmdata_all, cells = samp_cells)
  sub_obj$MP15_activity <- mp15_group[colnames(sub_obj)]

  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

  p_umap <- DimPlot(
    sub_obj,
    reduction = "umap",
    group.by = "MP15_activity",
    cols = c(MP15_high = "#D73027", MP15_low = "#4575B4"),
    pt.size = 0.4,
    shuffle = TRUE
  ) +
    labs(title = paste0("UMAP: ", samp), subtitle = "MP15 high vs low") +
    theme_minimal()

  safe_name <- gsub("[^A-Za-z0-9_]", "_", samp)
  ggsave(filename = paste0("Auto_MP15_UMAP_", safe_name, ".pdf"),
         plot = p_umap, width = 6, height = 5)
  cat(sprintf("  Saved UMAP plot for %s\n", samp))
}

####################
# Save summary CSV if any
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df, file = "Auto_MP15_CNV_summary.csv", row.names = FALSE)
  updates_dir <- "../updates/02mar/summaries"
  if (dir.exists(updates_dir)) {
    write.csv(summary_df,
              file = file.path(updates_dir, "Auto_MP15_CNV_summary.csv"),
              row.names = FALSE)
  } else {
    cat(sprintf("WARNING: updates summaries dir not found: %s\n", updates_dir))
  }
  cat("Saved: Auto_MP15_CNV_summary.csv\n")
}
####################

####################
# 8C. Combined PDF: UMAP → CNV → CNA signal (Top-3 samples)
####################
cat("\n========================================\n")
cat("APPROACH A3: Combined PDF for Top-3 samples\n")
cat("========================================\n")

combined_pdf <- "Auto_MP15_top3_combined.pdf"
pdf(combined_pdf, width = 9, height = 7)
####################
#################### 20#################### 20####################
# Combined PDF - ordered: all UMAPs, then all CNV heatmaps, then all CNA signal plots
#################### 20#################### 20####################

# 1) UMAPs for all top3 samples
for (samp in top3_samples) {
  cat(sprintf("\n[Combined-UMAP] Processing sample: %s\n", samp))
  samp_cells <- common_cells[tmdata_all$orig.ident[common_cells] == samp]
  if (length(samp_cells) < 20) next
  samp_mp15_z <- mp_adj[samp_cells, "MP15"]
  mp15_group <- ifelse(samp_mp15_z > MP15_Z_THRESHOLD, "MP15_high", "MP15_low")
  names(mp15_group) <- samp_cells
  n_high <- sum(mp15_group == "MP15_high")
  n_low  <- sum(mp15_group == "MP15_low")

  sub_obj <- subset(tmdata_all, cells = samp_cells)
  sub_obj$MP15_activity <- mp15_group[colnames(sub_obj)]
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

  p_umap <- DimPlot(
    sub_obj,
    reduction = "umap",
    group.by = "MP15_activity",
    cols = c(MP15_high = "#D73027", MP15_low = "#4575B4"),
    pt.size = 0.4,
    shuffle = TRUE
  ) + labs(title = paste0("UMAP: ", samp), subtitle = paste0("MP15_high: ", n_high, " | MP15_low: ", n_low)) + theme_minimal()
  print(p_umap)
}

# 2) CNV heatmaps for all top3 samples (fallback style from cnv_subsetting.R)
for (samp in top3_samples) {
  cat(sprintf("\n[Combined-CNV] Processing sample: %s\n", samp))
  samp_cells <- common_cells[tmdata_all$orig.ident[common_cells] == samp]
  if (length(samp_cells) < 20) next

  samp_mp15_z <- mp_adj[samp_cells, "MP15"]
  mp15_group <- ifelse(samp_mp15_z > MP15_Z_THRESHOLD, "MP15_high", "MP15_low")
  names(mp15_group) <- samp_cells
  n_high <- sum(mp15_group == "MP15_high")
  n_low  <- sum(mp15_group == "MP15_low")

  outs_path <- file.path("by_samples", samp, paste0(samp, "_outs.rds"))
  epi_file  <- file.path("by_samples", samp, paste0(samp, "_epi_f.rds"))
  if (!file.exists(outs_path)) { cat("  outs not found, skipping\n"); next }
  outs <- readRDS(outs_path)
  common_cnv_cells <- intersect(colnames(outs), samp_cells)
  if (length(common_cnv_cells) < 20) { cat("  too few CNV cells, skipping\n"); next }
  outs <- outs[, common_cnv_cells, drop = FALSE]

  ####################
  # Balance MP15_high vs MP15_low cell counts for comparable heatmap width
  mp15_group_cnv <- mp15_group[colnames(outs)]
  keep_high <- names(mp15_group_cnv)[mp15_group_cnv == "MP15_high"]
  keep_low  <- names(mp15_group_cnv)[mp15_group_cnv == "MP15_low"]
  n_each <- min(length(keep_high), length(keep_low))
  if (n_each > 0) {
    set.seed(42)
    keep_high <- sample(keep_high, n_each)
    keep_low  <- sample(keep_low, n_each)
    keep_cells <- c(keep_high, keep_low)
    outs <- outs[, keep_cells, drop = FALSE]
  }
  ####################

  # Binning (reuse helper)
  bin_res <- bin_cnv_matrix(outs, gene_order, bin_size = 100L)
  binned_mat <- bin_res$binned_mat

  # chromosome per bin
  row_chr_labels <- sub("_.*$", "", rownames(binned_mat))
  row_chr <- factor(row_chr_labels, levels = unique(row_chr_labels))

  # load epi meta for top annotation if available
  meta <- NULL
  if (file.exists(epi_file)) {
    epi <- readRDS(epi_file)
    meta <- epi@meta.data
    meta <- meta[colnames(outs), , drop = FALSE]
  } else {
    meta <- data.frame(row.names = colnames(outs))
  }

  # build top annotation fields if present
  ann_list <- list()
  ann_cols <- list()
  if ("cs_score" %in% colnames(meta)) { ann_list$cancer_signature <- as.numeric(meta$cs_score); ann_cols$cancer_signature <- colorRamp2(c(min(meta$cs_score, na.rm=TRUE), median(meta$cs_score, na.rm=TRUE), max(meta$cs_score, na.rm=TRUE)), c("white","grey80","black")) }
  if ("cc_score" %in% colnames(meta)) { ann_list$cell_cycling <- as.numeric(meta$cc_score); ann_cols$cell_cycling <- colorRamp2(c(min(meta$cc_score, na.rm=TRUE), median(meta$cc_score, na.rm=TRUE), max(meta$cc_score, na.rm=TRUE)), c("white","grey80","black")) }
  if ("study" %in% colnames(meta)) { ann_list$study <- factor(meta$study); study_levels <- levels(droplevels(ann_list$study)); ann_cols$study <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(study_levels)), study_levels) }

  if (length(ann_list) > 0) {
    top_ha <- do.call(HeatmapAnnotation, c(ann_list, list(col = ann_cols, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9), annotation_name_offset = unit(2, "mm"), annotation_height = unit(rep(4, length(ann_list)), "mm"), show_legend = TRUE)))
  } else {
    top_ha <- NULL
  }

  # left chr bar
  chr_used <- levels(droplevels(row_chr))
  base_cols <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))
  chr_cols <- setNames(base_cols[seq_along(chr_used)], chr_used)
  left_chr_bar <- rowAnnotation(chr = row_chr, col = list(chr = chr_cols), show_annotation_name = FALSE, show_legend = FALSE, gp = gpar(col = NA), width = unit(4, "mm"))

  # column split by mp15_group
  col_split <- factor(mp15_group[colnames(binned_mat)], levels = c("MP15_high", "MP15_low"))

  # chromosome boundaries
  chr_bounds <- which(head(row_chr_labels, -1L) != tail(row_chr_labels, -1L))
  line_gp <- gpar(col = "black", lwd = 1, lineend = "square")

  # stronger contrast for CNV
  # compute robust color limits using quantiles to avoid faint colours
  flat_vals <- as.numeric(binned_mat)
  if (all(is.na(flat_vals))) {
    qlims <- c(-2, 2)
  } else {
    qlims <- quantile(flat_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    if (!is.finite(qlims[1]) || !is.finite(qlims[2]) || qlims[1] == qlims[2]) {
      qlims <- c(min(flat_vals, na.rm = TRUE), max(flat_vals, na.rm = TRUE))
    }
    # expand towards symmetric around zero if required
    qmax <- max(abs(qlims))
    qlims <- c(-qmax, qmax)
  }
  col_fun <- colorRamp2(c(qlims[1], 0, qlims[2]), c("#08306B", "white", "#99000D"))

  ht <- Heatmap(
    binned_mat,
    name = "CNV",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    column_split = col_split,
    column_title_rot = 0,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    column_gap = unit(2, "mm"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ha,
    left_annotation = left_chr_bar,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    row_title_rot = 0,
    rect_gp = gpar(col = NA),
    border = NA,
    layer_fun = function(j, i, x, y, w, h, fill) {
      hits <- intersect(i, chr_bounds)
      if (length(hits)) {
        id <- match(hits, i)
        yy <- y[id] - h[id] / 2
        grid.segments(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = yy, y1 = yy, gp = line_gp)
      }
    }
  )

  # Remove explicit bold annotation on the heatmap canvas (sample name & counts)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
}

# 3) CNA signal plots for all top3 samples
for (samp in top3_samples) {
  cat(sprintf("\n[Combined-CNA] Processing sample: %s\n", samp))
  samp_cells <- common_cells[tmdata_all$orig.ident[common_cells] == samp]
  if (length(samp_cells) < 20) next
  samp_mp15_z <- mp_adj[samp_cells, "MP15"]
  mp15_group <- ifelse(samp_mp15_z > MP15_Z_THRESHOLD, "MP15_high", "MP15_low")
  names(mp15_group) <- samp_cells
  outs_path <- file.path("by_samples", samp, paste0(samp, "_outs.rds"))
  if (!file.exists(outs_path)) { cat("  outs not found, skipping CNA plot\n"); next }
  outs <- readRDS(outs_path)
  common_cnv_cells <- intersect(colnames(outs), samp_cells)
  if (length(common_cnv_cells) < 20) { cat("  too few CNV cells, skipping\n"); next }
  outs <- outs[, common_cnv_cells, drop = FALSE]
  cna_signal <- compute_cna_signal(outs, use_quantile = 0.9)
  cna_df <- data.frame(cell = names(cna_signal), cna_signal = as.numeric(cna_signal), mp15_group = mp15_group[names(cna_signal)], stringsAsFactors = FALSE)
  cna_df <- cna_df[!is.na(cna_df$mp15_group), ]
  n_high <- sum(cna_df$mp15_group == "MP15_high", na.rm = TRUE)
  n_low <- sum(cna_df$mp15_group == "MP15_low", na.rm = TRUE)
  p_cna <- ggplot(cna_df, aes(x = mp15_group, y = cna_signal, fill = mp15_group)) + geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) + scale_fill_manual(values = c(MP15_high = "#D73027", MP15_low = "#4575B4")) + labs(title = paste0("CNA signal: ", samp), subtitle = paste0("MP15_high: ", n_high, " | MP15_low: ", n_low), x = NULL, y = "CNA signal (mean sq.)") + theme_classic() + theme(legend.position = "none")
  print(p_cna)
}

dev.off()
cat(sprintf("Saved combined PDF: %s\n", combined_pdf))

#################### 20#################### 20####################
####################

# ============================================================================
# 9. APPROACH B: Compare all MP15_high cells vs an equal number of MP15_low cells
# ============================================================================
cat("\n========================================\n")
cat("APPROACH B: Compare all MP15_high cells vs MP15_low cells\n")
cat("========================================\n")

cells_mp15_high <- common_cells[mp15_high[common_cells]]
cells_mp15_low <- common_cells[!mp15_high[common_cells]]

cat(sprintf("Total MP15 high cells: %d\n", length(cells_mp15_high)))
cat(sprintf("Total MP15 low cells: %d\n", length(cells_mp15_low)))

set.seed(42)
# Sample MP15 low cells to match the count of MP15 high cells
n_high <- length(cells_mp15_high)
if (length(cells_mp15_low) > n_high) {
  cells_mp15_low_sub <- sample(cells_mp15_low, n_high)
} else {
  cells_mp15_low_sub <- cells_mp15_low
}

cells_plot_b <- c(cells_mp15_high, cells_mp15_low_sub)

# Subsample combined if > 20000 to avoid memory issues with huge heatmaps
if (length(cells_plot_b) > 20000) {
  # Subsample high and low equally
  n_each <- 10000
  cells_mp15_high_sub <- sample(cells_mp15_high, n_each)
  cells_mp15_low_sub <- sample(cells_mp15_low_sub, n_each)
  cells_plot_b <- c(cells_mp15_high_sub, cells_mp15_low_sub)
  cat(sprintf("Subsampled combined to %d cells for plotting\n", length(cells_plot_b)))
}

tmdata_mp15_b <- subset(tmdata_all, cells = cells_plot_b)

# Assign MP15 activity
b_mp15_z <- mp_adj[cells_plot_b, "MP15"]
tmdata_mp15_b$MP15_activity <- ifelse(b_mp15_z > MP15_Z_THRESHOLD,
                                      "MP15_high", "MP15_low")

cat(sprintf("Cells in plot - MP15 high: %d | MP15 low: %d\n",
            sum(tmdata_mp15_b$MP15_activity == "MP15_high"),
            sum(tmdata_mp15_b$MP15_activity == "MP15_low")))

# Use MP15_activity as identity for column split
png("Auto_MP15_allSamples_high_vs_low.png",
    width = 80, height = 50, units = "in", res = 400)
plot_qc_heatmap(tmdata_mp15_b,
                sampleid = "All_Samples_MP15_Activity",
                identity = "MP15_activity",
                reorder = TRUE,
                ct_reorder = c("MP15_high", "MP15_low"))
dev.off()
cat("Saved: Auto_MP15_allSamples_high_vs_low.png\n")

# ============================================================================
# 10. QUICK CHECK: QC heatmap for Alcindor_2025_SRR27335939 (all cell types)
# ============================================================================
cat("\n========================================\n")
cat("QUICK CHECK: Alcindor_2025_SRR27335939 QC heatmap (all cell types)\n")
cat("========================================\n")

qc_sample <- "Alcindor_2025_SRR27335939"
qc_path <- file.path("by_samples", qc_sample, paste0(qc_sample, "_anno.rds"))

if (!file.exists(qc_path)) {
  cat(sprintf("WARNING: %s not found. Skipping.\n", qc_path))
} else {
  cat(sprintf("Loading %s...\n", qc_path))
  tmdata_qc <- readRDS(qc_path)
  cat(sprintf("  Raw cells: %d\n", ncol(tmdata_qc)))

  # Standard filtering: coexpression_loose == "singlet" & marker_expression == "good"
  has_coexp <- "coexpression_loose" %in% colnames(tmdata_qc@meta.data)
  has_marker <- "marker_expression" %in% colnames(tmdata_qc@meta.data)

  if (has_coexp && has_marker) {
    tmdata_qc_filt <- tryCatch({
      subset(tmdata_qc,
             subset = coexpression_loose == "singlet" &
                      marker_expression == "good")
    }, error = function(e) {
      cat(sprintf("  Filter error: %s\n", e$message))
      tmdata_qc
    })
  } else if (has_coexp) {
    tmdata_qc_filt <- subset(tmdata_qc, subset = coexpression_loose == "singlet")
    cat("  Note: marker_expression column not found, filtering by coexpression_loose only\n")
  } else {
    tmdata_qc_filt <- tmdata_qc
    cat("  Note: filtering columns not found, using all cells\n")
  }

  cat(sprintf("  After filtering: %d cells\n", ncol(tmdata_qc_filt)))

  if (ncol(tmdata_qc_filt) > 10) {
    cat("\n  Cell type composition:\n")
    print(table(tmdata_qc_filt$celltype_update))

    ct_order_qc <- c("b.cell", "dendritic", "endothelial", "epithelial",
                     "fibroblast", "macrophage", "mast", "nk.cell",
                     "nk.cell|t.cell", "plasma", "t.cell", "t.cell|nk.cell",
                     "unresolved", "unresolved_inconsistent")

    png(paste0("Auto_QC_", qc_sample, ".png"),
        width = 80, height = 50, units = "in", res = 400)
    plot_qc_heatmap(tmdata_qc_filt,
                    sampleid = qc_sample,
                    identity = "celltype_update",
                    reorder = TRUE,
                    ct_reorder = ct_order_qc)
    dev.off()
    cat(sprintf("  Saved: Auto_QC_%s.png\n", qc_sample))
  } else {
    cat("  Too few cells after filtering, skipping plot.\n")
  }

  rm(tmdata_qc, tmdata_qc_filt)
  gc()
}

cat("\n====== DONE ======\n")
