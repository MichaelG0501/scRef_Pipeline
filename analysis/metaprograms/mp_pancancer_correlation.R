####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/mp_pancancer_correlation.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_task5_mp_ucell_scoring_unresolved_pancancer.R
# Separate unresolved-aware MP UCell visualization script:
# - Includes canonical MPs (mp_tree_order) + retained pan-cancer MPs used in unresolved relabel
# - Names are always "MPx description" (for canonical MPs) and clean pan-cancer names
# - Produces per-sample correlation heatmap and Jaccard overlap heatmap
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(dplyr)
library(grid)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "task5"
out_dir <- file.path("Metaprogrammes_Results", paste0(task_prefix, "_mp_ucell"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
u3ca <- readRDS("UCell_3CA_MPs.rds")

coverage_candidates <- c(
  "task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv",
  "unresolved_states/Auto_unresolved_relabel_mp_coverage.csv"
)
coverage_path <- coverage_candidates[file.exists(coverage_candidates)][1]
if (is.na(coverage_path)) {
  stop("Missing unresolved MP coverage csv in task4 or unresolved_states folder")
}
cov_df <- read.csv(coverage_path, stringsAsFactors = FALSE)
retained_3ca <- cov_df %>%
  dplyr::filter(n_samples >= 50, n_studies >= 6, pct_cells >= 1) %>%
  dplyr::arrange(desc(n_cells)) %>%
  dplyr::pull(mp_label)

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp.genes)))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP5"  = "Epithelial IFN Resp.",
  "MP7"  = "DNA Damage Repair",
  "MP8"  = "Intestinal Diff.",
  "MP9"  = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)"
)

mp_display <- paste0(mp_tree_order_names, " ", ifelse(is.na(mp_descriptions[mp_tree_order_names]), mp_tree_order_names, mp_descriptions[mp_tree_order_names]))

common_cells <- intersect(Cells(tmdata_all), rownames(ucell_scores))
common_cells <- intersect(common_cells, rownames(u3ca))
tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
u3ca <- u3ca[common_cells, , drop = FALSE]

module_scores <- scale(as.matrix(ucell_scores))
canon_mps_avail <- mp_tree_order_names[mp_tree_order_names %in% colnames(module_scores)]

retained_3ca_corr <- retained_3ca[retained_3ca %in% colnames(u3ca)]
u3ca_scaled <- scale(as.matrix(u3ca[, retained_3ca_corr, drop = FALSE]))

mod_mat_canon <- t(module_scores[, canon_mps_avail, drop = FALSE])
rownames(mod_mat_canon) <- paste0(canon_mps_avail, ": ", ifelse(is.na(mp_descriptions[canon_mps_avail]), canon_mps_avail, mp_descriptions[canon_mps_avail]))

mod_mat_pan <- t(u3ca_scaled[, retained_3ca_corr, drop = FALSE])
rownames(mod_mat_pan) <- clean_3ca_name(rownames(mod_mat_pan))

mod_mat <- rbind(mod_mat_canon, mod_mat_pan)

samples_vec <- tmdata_all$orig.ident[match(colnames(mod_mat), Cells(tmdata_all))]
samples_vec[is.na(samples_vec)] <- colnames(mod_mat)[is.na(samples_vec)]
samples <- unique(samples_vec)

mp_names <- rownames(mod_mat)
n_mps <- length(mp_names)
cor_array <- array(NA_real_, dim = c(n_mps, n_mps, length(samples)), dimnames = list(mp_names, mp_names, samples))

for (samp in samples) {
  cells_in_sample <- colnames(mod_mat)[samples_vec == samp]
  if (length(cells_in_sample) < 3) next
  sub_mat <- mod_mat[, cells_in_sample, drop = FALSE]
  cor_array[, , samp] <- cor(t(sub_mat), method = "spearman")
}

z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

mean_rho <- matrix(0, n_mps, n_mps, dimnames = list(mp_names, mp_names))
p_vals <- matrix(1, n_mps, n_mps, dimnames = list(mp_names, mp_names))
for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    if (i == j) {
      mean_rho[i, j] <- 1
      p_vals[i, j] <- 0
    } else {
      z_scores <- z_array[i, j, ]
      z_scores <- z_scores[is.finite(z_scores)]
      if (length(z_scores) < 3) {
        mean_rho[i, j] <- NA_real_
        p_vals[i, j] <- NA_real_
      } else {
        test_res <- t.test(z_scores)
        mean_rho[i, j] <- tanh(mean(z_scores))
        p_vals[i, j] <- test_res$p.value
      }
    }
  }
}

col_cor <- colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_nMP19_correlation_heatmap_persample.pdf")), width = 10, height = 9, useDingbats = FALSE)
Heatmap(
  mean_rho,
        name = paste0("Mean Rho\n(", length(samples), " Samples)"),
        col = col_cor,
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- p_vals[i, j]
          rho <- mean_rho[i, j]
          if (is.na(p) || is.na(rho)) {
            grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "grey50"))
          } else if (p < 0.001) {
            grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 10))
          } else if (p < 0.01) {
            grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 10))
          } else if (p < 0.05) {
            grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 10))
          } else {
            grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 10))
          }
        },
        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"))
)
dev.off()

retained_3ca <- retained_3ca[retained_3ca %in% colnames(u3ca)]

canon_gene_sets <- mp.genes[mp_tree_order_names[mp_tree_order_names %in% names(mp.genes)]]
names(canon_gene_sets) <- paste0(names(canon_gene_sets), ": ", ifelse(is.na(mp_descriptions[names(canon_gene_sets)]), names(canon_gene_sets), mp_descriptions[names(canon_gene_sets)]))

nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
pan_sets <- MP_list[retained_3ca]
names(pan_sets) <- clean_3ca_name(names(pan_sets))

all_sets <- c(canon_gene_sets, pan_sets)
all_sets <- all_sets[sapply(all_sets, length) > 0]

all_names <- names(all_sets)
universe <- unique(unlist(all_sets))
jaccard_mat <- matrix(NA_real_, length(all_names), length(all_names), dimnames = list(all_names, all_names))
overlap_n_mat <- matrix(NA_real_, length(all_names), length(all_names), dimnames = list(all_names, all_names))
pval_mat <- matrix(NA_real_, length(all_names), length(all_names), dimnames = list(all_names, all_names))

for (i in seq_along(all_sets)) {
  A <- unique(all_sets[[i]])
  for (j in seq_along(all_sets)) {
    B <- unique(all_sets[[j]])
    inter <- length(intersect(A, B))
    uni <- length(union(A, B))
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j] <- if (uni == 0) NA_real_ else inter / uni

    a <- inter
    b <- length(setdiff(A, B))
    cc <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_ else fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")$p.value
  }
}

padj_mat <- matrix(
  p.adjust(as.vector(pval_mat), method = "BH"),
  nrow = nrow(pval_mat),
  ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)
stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05] <- "*"
stars_mat[padj_mat < 0.01] <- "**"
stars_mat[padj_mat < 0.001] <- "***"
display_mat <- matrix(
  paste0(overlap_n_mat, "\n", stars_mat),
  nrow = nrow(overlap_n_mat),
  ncol = ncol(overlap_n_mat),
  dimnames = dimnames(overlap_n_mat)
)

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_nMP19_jaccard_with_pancancer_heatmap.pdf")), width = 12, height = 10, useDingbats = FALSE)
pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "MP Gene Set Overlap (Jaccard Index) - nMP19",
  labels_row = rownames(jaccard_mat),
  labels_col = colnames(jaccard_mat),
  angle_col = "90",
  display_numbers = display_mat,
  fontsize_number = 8,
  number_color = "black",
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
)
dev.off()

####################
# NEW STYLED PDF PLOTS (matching final_mp_correlation.R)
####################
# 1. State definitions and colors
styled_group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Repeiration" = "#FF9999",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Other" = "grey80"
)

styled_mp_to_state <- c(
  "MP2" = "Classic Proliferative",
  "MP17" = "Basal to Intestinal Metaplasia",
  "MP14" = "Basal to Intestinal Metaplasia",
  "MP5" = "Basal to Intestinal Metaplasia",
  "MP10" = "Basal to Intestinal Metaplasia",
  "MP8" = "Basal to Intestinal Metaplasia",
  "MP13" = "Stress-adaptive",
  "MP12" = "Stress-adaptive",
  "MP18" = "SMG-like Metaplasia",
  "MP16" = "SMG-like Metaplasia",
  "MP15" = "Immune Infiltrating",
  "MP1" = "Cell cycle",
  "MP9" = "Cell cycle",
  "MP7" = "Cell cycle",
  "3CA_mp_30 Respiration 1" = "Repeiration",
  "3CA_mp_12 Protein maturation" = "3CA_EMT_and_Protein_maturation",
  "3CA_mp_17 EMT III" = "3CA_EMT_and_Protein_maturation"
)

raw_mp_keys <- c(canon_mps_avail, clean_3ca_name(retained_3ca_corr))
state_assign <- styled_mp_to_state[raw_mp_keys]
state_assign[is.na(state_assign)] <- "Other"

# Remove Cell cycle MPs
non_cc_idx <- state_assign != "Cell cycle"
state_assign <- state_assign[non_cc_idx]

styled_mean_rho <- mean_rho[non_cc_idx, non_cc_idx, drop = FALSE]
styled_p_vals <- p_vals[non_cc_idx, non_cc_idx, drop = FALSE]
styled_jaccard_mat <- jaccard_mat[non_cc_idx, non_cc_idx, drop = FALSE]
styled_display_mat <- display_mat[non_cc_idx, non_cc_idx, drop = FALSE]

state_vec_for_mps <- factor(state_assign, levels = c(
  "Classic Proliferative",
  "Repeiration",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation",
  "Other"
))

# Row/Col annotations
ha_left <- rowAnnotation(
  State = state_vec_for_mps,
  col = list(State = styled_group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)
ha_top <- HeatmapAnnotation(
  State = state_vec_for_mps,
  col = list(State = styled_group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm_width <- unit(9, "inch")
hm_height <- unit(9, "inch")
styled_col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_nMP19_correlation_heatmap_persample_styled.pdf")), width = 16, height = 16, useDingbats = FALSE)
ht_cor_styled <- Heatmap(styled_mean_rho,
        name = paste0("Mean Rho\n(", length(samples), " Samples)"),
        col = styled_col_cor,
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        left_annotation = ha_left,
        top_annotation = ha_top,
        row_split = state_vec_for_mps,
        column_split = state_vec_for_mps,
        
        column_title_rot = 20,
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        row_title = NULL,
        
        row_names_side = "right",
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        
        width = hm_width,
        height = hm_height,
        show_heatmap_legend = FALSE,
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- styled_p_vals[i, j]
          rho <- styled_mean_rho[i, j]
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
        })
draw(ht_cor_styled, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_nMP19_jaccard_with_pancancer_heatmap_styled.pdf")), width = 16, height = 16, useDingbats = FALSE)
ht_jaccard_styled <- Heatmap(
  styled_jaccard_mat,
  name = "Jaccard Index",
  col = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100),
  rect_gp = gpar(col = "grey85", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = ha_left,
  top_annotation = ha_top,
  row_split = state_vec_for_mps,
  column_split = state_vec_for_mps,
  
  column_title_rot = 20,
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_title = NULL,
  
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 16, fontface = "bold"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  
  width = hm_width,
  height = hm_height,
  show_heatmap_legend = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(styled_display_mat[i, j], x, y, gp = gpar(fontsize = 12))
  }
)
draw(ht_jaccard_styled, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  task = task_prefix,
  n_samples = length(samples),
  n_canonical_mps = length(canon_gene_sets),
  n_retained_pancancer_mps = length(pan_sets),
  stringsAsFactors = FALSE
), file.path(summary_dir, paste0("Auto_", task_prefix, "_mp_ucell_unresolved_pancancer_summary.csv")), row.names = FALSE)

cat("Saved task5 MP-UCell unresolved pan-cancer plots in:", out_dir, "\n")
