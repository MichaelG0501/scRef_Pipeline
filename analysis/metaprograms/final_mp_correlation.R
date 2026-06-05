####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/final_mp_correlation.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_final_mp_correlation.R
# Finalized MP Correlation (Fisher Z) and Jaccard similarity plots.
# Order follows state definitions as requested by USER.
# Excludes CC MPs and keeps only the curated state-associated MPs.
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ============================================================================
# 1. Load data
# ============================================================================

geneNMF.metaprograms <- readRDS(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)
mp.genes <- geneNMF.metaprograms$metaprograms.genes

tmdata_all <- readRDS("EAC_Ref_epi.rds")
relabeled_state_path <- "unresolved_states/Auto_unresolved_relabel_states.rds"
relabeled_states <- NULL
if (file.exists(relabeled_state_path)) {
  relabeled_states <- readRDS(relabeled_state_path)
}

# ============================================================================
# 2. Define Final Order and Subset (Excluding CC MPs)
# ============================================================================

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)

final_mps <- unlist(state_groups, use.names = FALSE)
cat("Requested MPs for final plots:", paste(final_mps, collapse = ", "), "\n")

# Filter mp.genes to only include these (and optionally check silhouette if desired, 
# but user requested exactly these in this order)
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)

# Log if any requested MPs are low quality
low_qual_requested <- intersect(final_mps, bad_mp_names)
if (length(low_qual_requested) > 0) {
  cat("Warning: Some requested MPs have silhouette < 0:", paste(low_qual_requested, collapse = ", "), "\n")
}

# Subset and enforce order
mp.genes <- mp.genes[final_mps]
# Ensure no NAs in case some MP wasn't in mp.genes (though for nMP=19 they should be)
mp.genes <- mp.genes[!is.na(names(mp.genes))]
final_mps_ordered <- names(mp.genes)

cat("Retained and ordered MPs:", paste(final_mps_ordered, collapse = ", "), "\n")

# ============================================================================
# 3. UCell scoring (retrieve or recalculate)
# ============================================================================

# We can reuse the filtered UCell scores if they exist, but we might need 
# to ensure we have all 11. The original script saved to Metaprogrammes_Results/UCell_nMP19_filtered.rds
ucell_path <- "Metaprogrammes_Results/UCell_nMP19_filtered.rds"
if (file.exists(ucell_path)) {
  cat("Loading UCell scores from:", ucell_path, "\n")
  all_ucell_scores <- readRDS(ucell_path)
} else {
  cat("Recalculating UCell scores...\n")
  tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp.genes, ncores = 1, name = "")
  all_ucell_scores <- tmdata_all@meta.data[, names(mp.genes), drop = FALSE]
}

# Only keep the final_mps_ordered in the scores matrix
ucell_scores <- all_ucell_scores[, final_mps_ordered, drop = FALSE]

# ============================================================================
# 4. UCell scores and matrix preparation
# ============================================================================

module_scores <- scale(as.matrix(ucell_scores))

# Order rows
mod_mat <- t(module_scores) # MPs x Cells (already ordered by final_mps_ordered)

# Descriptions (just use names for now as in original script)
mp_descriptions <- setNames(rownames(mod_mat), rownames(mod_mat))

# ============================================================================
# 5. Name mapping
# ============================================================================

mp_desc_map <- c(
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
  "MP12" = "Neuro-responsive Epi."
)

# Apply name mapping: "MPX Description"
mp_names_with_desc <- setNames(
  paste(final_mps_ordered, mp_desc_map[final_mps_ordered]),
  final_mps_ordered
)
# Fallback for any missing descriptions
missing <- is.na(mp_desc_map[final_mps_ordered])
if (any(missing)) {
  mp_names_with_desc[final_mps_ordered[missing]] <- final_mps_ordered[missing]
}

# Mapping from MP to State Group for annotations
mp_to_state <- unlist(lapply(names(state_groups), function(gn) {
  setNames(rep(gn, length(state_groups[[gn]])), state_groups[[gn]])
}))

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8"
)

# ============================================================================
# 6. Correlation Heatmap — per-sample Fisher Z meta-analysis
# ============================================================================

samples_vec <- tmdata_all$orig.ident[match(colnames(mod_mat), Cells(tmdata_all))]
samples_vec[is.na(samples_vec)] <- colnames(mod_mat)[is.na(samples_vec)]
samples <- unique(samples_vec)
mps <- rownames(mod_mat)
n_mps <- length(mps)

cat("Computing per-sample correlations across", length(samples), "samples...\n")

cor_array <- array(NA, dim = c(n_mps, n_mps, length(samples)),
                   dimnames = list(mps, mps, samples))

for (samp in samples) {
  cells_in_sample <- colnames(mod_mat)[samples_vec == samp]
  if (length(cells_in_sample) < 10) next  # skip tiny samples for stability
  sub_mat <- mod_mat[, cells_in_sample, drop = FALSE]
  cor_array[, , samp] <- cor(t(sub_mat), method = "spearman")
}

z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

mean_rho <- matrix(0, n_mps, n_mps, dimnames = list(mps, mps))
p_vals   <- matrix(1, n_mps, n_mps, dimnames = list(mps, mps))

for (i in 1:n_mps) {
  for (j in 1:n_mps) {
    if (i == j) {
      mean_rho[i, j] <- 1
      p_vals[i, j] <- 0
    } else {
      z_scores <- z_array[i, j, ]
      z_scores <- z_scores[is.finite(z_scores)]
      if (length(z_scores) < 3) {
        mean_rho[i, j] <- NA
        p_vals[i, j] <- NA
      } else {
        test_res <- t.test(z_scores)
        mean_rho[i, j] <- tanh(mean(z_scores))
        p_vals[i, j] <- test_res$p.value
      }
    }
  }
}

# Update names for visualization
rownames(mean_rho) <- mp_names_with_desc[rownames(mean_rho)]
colnames(mean_rho) <- mp_names_with_desc[colnames(mean_rho)]

# Build annotations
state_vec_for_mps <- factor(mp_to_state[final_mps_ordered], levels = names(state_groups))

# State colors only, no legend
ha_left <- rowAnnotation(
  State = state_vec_for_mps,
  col = list(State = group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)
ha_top <- HeatmapAnnotation(
  State = state_vec_for_mps,
  col = list(State = group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Common heatmap settings for square shape and labels
hm_width <- unit(9, "inch")
hm_height <- unit(9, "inch")

col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

pdf("Metaprogrammes_Results/Auto_final_mp_correlation_heatmap.pdf", width = 16, height = 16, useDingbats = FALSE)
ht_cor <- Heatmap(mean_rho,
        name = paste0("Mean Rho\n(", length(samples), " Samples)"),
        col = col_cor,
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        left_annotation = ha_left,
        top_annotation = ha_top,
        row_split = state_vec_for_mps,
        column_split = state_vec_for_mps,
        
        # Labels and rotation
        column_title_rot = 20,
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        row_title = NULL, # Remove state labels on the left
        
        # MP names on bottom and right
        row_names_side = "right",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        
        # Square shape
        width = hm_width,
        height = hm_height,
        
        # Overlay significance stars
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- p_vals[i, j]
          rho <- mean_rho[i, j]
          if (is.na(p) || is.na(rho)) {
            grid.text("NA", x, y, gp = gpar(fontsize = 12, col = "grey50"))
          } else if (p < 0.001) {
            grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 14))
          } else if (p < 0.01) {
            grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 14))
          } else if (p < 0.05) {
            grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 14))
          } else {
            grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 14))
          }
        },
        heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14)))
draw(ht_cor, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()
cat("Saved: Metaprogrammes_Results/Auto_final_mp_correlation_heatmap.pdf\n")

# ============================================================================
# 7. Jaccard Self-Similarity Heatmap
# ============================================================================

mp_list <- mp.genes # Already ordered
mp_list <- lapply(mp_list, unique)
mp_names <- names(mp_list)
universe <- unique(unlist(mp_list))

jaccard_mat   <- matrix(NA_real_, length(mp_list), length(mp_list),
                        dimnames = list(mp_names, mp_names))
overlap_n_mat <- jaccard_mat
pval_mat      <- jaccard_mat

for (i in seq_along(mp_list)) {
  A <- mp_list[[i]]
  for (j in seq_along(mp_list)) {
    B <- mp_list[[j]]
    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j]   <- if (uni == 0) NA_real_ else inter / uni
    a <- inter
    b <- length(setdiff(A, B))
    cc <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_
    else fisher.test(matrix(c(a, b, cc, d), nrow = 2),
                     alternative = "greater")$p.value
  }
}

padj_mat <- matrix(
  p.adjust(as.vector(pval_mat), method = "BH"),
  nrow = nrow(pval_mat), ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat),
                    dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05]  <- "*"
stars_mat[padj_mat < 0.01]  <- "**"
stars_mat[padj_mat < 0.001] <- "***"

display_mat <- matrix(
  paste0(overlap_n_mat, "\n", stars_mat),
  nrow = nrow(overlap_n_mat),
  ncol = ncol(overlap_n_mat),
  dimnames = dimnames(overlap_n_mat)
)

# Update names for visualization
rownames(jaccard_mat) <- mp_names_with_desc[rownames(jaccard_mat)]
colnames(jaccard_mat) <- mp_names_with_desc[colnames(jaccard_mat)]

pdf("Metaprogrammes_Results/Auto_final_mp_jaccard_heatmap.pdf", width = 16, height = 16, useDingbats = FALSE)
ht_jaccard <- Heatmap(
  jaccard_mat,
  name = "Jaccard Index",
  col = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100),
  rect_gp = gpar(col = "grey85", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = ha_left,
  top_annotation = ha_top,
  row_split = state_vec_for_mps,
  column_split = state_vec_for_mps,
  
  # Labels and rotation
  column_title_rot = 20,
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_title = NULL, # Remove state labels on the left
  
  # MP names on bottom and right
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 16, fontface = "bold"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  
  # Square shape
  width = hm_width,
  height = hm_height,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 12))
  },
  
  heatmap_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14))
)
draw(ht_jaccard, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()
cat("Saved: Metaprogrammes_Results/Auto_final_mp_jaccard_heatmap.pdf\n")