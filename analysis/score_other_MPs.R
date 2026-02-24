library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ══════════════════════════════════════════════════════════════════════════════
# USER OPTIONS — toggle these before running
# ══════════════════════════════════════════════════════════════════════════════
REMOVE_CELL_CYCLE    <- TRUE    # ← set TRUE to drop CC MPs and regress CC from others
ETA_THRESHOLD        <- 1    # ← η² cutoff for study-confounding filter
SIG_THRESHOLD        <- 0.1     # ← 99th-percentile cutoff for signal filter
CHOSEN_RES           <- 1     # ← Louvain resolution
MAX_CELLS_PER_STATE  <- 500     # ← downsample for single-cell heatmap
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Load data ─────────────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
#tmdata_all <- readRDS("EAC_Ref_epi.rds")

# cat("Reading New_NMFs.csv...\n")
# MP_df <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", check.names = FALSE)
# MP_list <- as.list(MP_df)
# 
# # Clean up empty strings and change names to be valid column names in Seurat
# MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
# names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
#
# tmdata_all <- AddModuleScore_UCell(tmdata_all, features = MP_list, ncores = 4, name = "")
# 
# # Extract scores
# score_cols <- grep("^X3CA_mp|^3CA_mp", colnames(tmdata_all@meta.data), value = TRUE)
# ucell_scores <- tmdata_all@meta.data[, score_cols, drop = FALSE]
# 
# # Save raw scores
# saveRDS(ucell_scores, file = "UCell_3CA_MPs.rds")

ucell_scores <- readRDS("UCell_3CA_MPs.rds")
score_cols   <- colnames(ucell_scores)
study_var    <- tmdata_all$study
sample_var   <- tmdata_all$orig.ident

# ══════════════════════════════════════════════════════════════════════════════
# CHANGE 1: Z-normalise within study BEFORE any filtering
#
#   Why: Studies with globally higher/lower UCell scores inflate the 99th
#   percentile and bias the ANOVA η². By z-scoring within study first, every
#   study contributes equally to the signal and confounding filters.
#   We keep the RAW scores for final interpretability heatmaps.
# ══════════════════════════════════════════════════════════════════════════════
cat("Using RAW UCell scores for all filtering (signal + eta^2)...\n")

# Raw UCell matrix (cells x MPs)
ucell_raw <- as.data.frame(ucell_scores)
ucell_raw$.cell <- rownames(ucell_raw)

# Ensure study_var is aligned to cells
study_var <- tmdata_all$study
names(study_var) <- rownames(tmdata_all@meta.data)
ucell_raw$.study <- study_var[ucell_raw$.cell]

# Raw matrix (cells x MPs)
ucell_mat <- as.matrix(ucell_raw[, score_cols, drop = FALSE])
rownames(ucell_mat) <- ucell_raw$.cell

# ── 1. Signal filter (RAW scores) ───────────────────────────────────────────
quantiles_99 <- apply(ucell_mat, 2, function(x) quantile(x, 0.99, na.rm = TRUE))

sig_threshold <- SIG_THRESHOLD
sig_mps <- names(quantiles_99)[quantiles_99 > sig_threshold]

cat(sprintf("After signal filter: %d / %d MPs (99th pctl of study-Z > %.2f)\n",
            length(sig_mps), length(score_cols), sig_threshold))

# ── 2. Study-confounding filter (η² on study-z-normalised scores) ────────────
cat("Computing eta-squared on study-Z scores...\n")
set.seed(42)
n_sub   <- min(15000, nrow(ucell_mat))
idx_sub <- sample(seq_len(nrow(ucell_mat)), n_sub)

eta_sq <- sapply(sig_mps, function(mp) {
  fit <- aov(ucell_mat[idx_sub, mp] ~ factor(study_var[idx_sub]))
  ss  <- summary(fit)[[1]]
  ss["Sum Sq"][1, 1] / sum(ss["Sum Sq"][, 1])
})

eta_df <- data.frame(MP = sig_mps, eta_sq = round(eta_sq, 4),
                     stringsAsFactors = FALSE) %>%
  arrange(desc(eta_sq))
eta_df$MP_clean <- gsub("\\.", " ", gsub("^X3CA_", "3CA_", eta_df$MP))
eta_df$keep     <- ifelse(eta_df$eta_sq <= ETA_THRESHOLD, "Keep", "Remove")

p_eta <- ggplot(eta_df, aes(x = reorder(MP_clean, eta_sq), y = eta_sq, fill = keep)) +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = ETA_THRESHOLD, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("Keep" = "steelblue", "Remove" = "tomato")) +
  coord_flip() +
  labs(title = "Study-confounding filter (computed on within-study Z-scores)",
       subtitle = paste0("Remove MPs with η² > ", ETA_THRESHOLD),
       x = NULL, y = expression(eta^2), fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
p_eta

robust_mps <- eta_df$MP[eta_df$eta_sq <= ETA_THRESHOLD]
cat(sprintf("After η² filter: %d / %d MPs (η² ≤ %.2f)\n",
            length(robust_mps), length(sig_mps), ETA_THRESHOLD))
cat("Removed (study-confounded):\n")
cat(paste(" ", eta_df$MP_clean[eta_df$eta_sq > ETA_THRESHOLD]), sep = "\n")

# ══════════════════════════════════════════════════════════════════════════════
# CHANGE 3 (OPTIONAL): Remove cell-cycle MPs & regress CC from remaining MPs
#
#   Step A: Identify cell-cycle MPs by grepping for "cell.cycle" (case-insensitive)
#   Step B: Remove them from the clustering feature set
#   Step C: For every remaining MP, fit   MP ~ CC_MP1 + CC_MP2 + ...
#           and replace the score with the residual, so that cell-cycle-
#           correlated variance is removed from all other MPs.
#   This runs on the WITHIN-STUDY Z-scored matrix so regressions aren't
#   biased by study-level shifts.
# ═════════════════════════════════════════════════════════════════════════���════
# ── Fixed CC MP list (use these five only) ──────────────────────────────────
CC_FIXED <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)

if (REMOVE_CELL_CYCLE) {
  cc_mps <- intersect(CC_FIXED, robust_mps)
  non_cc_mps <- setdiff(robust_mps, cc_mps)
  
  cat(sprintf("\nCell-cycle MPs (fixed list) found in robust_mps: %d\n", length(cc_mps)))
  if (length(cc_mps)) cat(paste(" ", cc_mps), sep = "\n")
  
  if (length(cc_mps) == 0) {
    cat("None of the fixed CC MPs present — skipping regression.\n")
    final_mps <- robust_mps
    ucell_for_clustering <- ucell_mat[, final_mps, drop = FALSE]
    
  } else {
    cat(sprintf("Regressing out CC signal using %d CC MPs; keeping %d MPs.\n",
                length(cc_mps), length(non_cc_mps)))
    
    X_cc    <- as.matrix(ucell_mat[, cc_mps, drop = FALSE])
    Y_other <- as.matrix(ucell_mat[, non_cc_mps, drop = FALSE])
    
    # Efficient OLS residualization: Y_resid = Y - X (X'X)^-1 X'Y
    X <- cbind(Intercept = 1, X_cc)
    XtX_inv <- solve(crossprod(X))
    B <- XtX_inv %*% crossprod(X, Y_other)
    Y_hat <- X %*% B
    Y_resid <- Y_other - Y_hat
    
    final_mps <- non_cc_mps
    ucell_for_clustering <- Y_resid
  }
  
} else {
  cat("\nSkipping cell-cycle removal (REMOVE_CELL_CYCLE = FALSE).\n")
  final_mps <- robust_mps
  ucell_for_clustering <- ucell_mat[, final_mps, drop = FALSE]
}

cat(sprintf("\nFinal MP set for clustering: %d MPs\n", length(final_mps)))


# ══════════════════════════════════════════════════════════════════════════════
# CHANGE 2: Center within sample, scale by WITHIN-STUDY global SD
#
#   Previously: each cell's score was centered by sample mean, then divided
#   by a single global SD computed across ALL cells/studies.
#   Problem: studies with inherently higher variance dominate the SD.
#
#   Now: for each MP, we compute the pooled SD separately for each study
#   (across all samples within that study). Each cell is centered by its
#   sample mean, then divided by its study's SD. This way the per-study
#   scaling matches each study's own dynamic range.
# ══════════════════════════════════════════════════════════════════════════════
cat("Preparing clustering matrix (center per sample, scale by within-study SD)...\n")

clust_df <- as.data.frame(ucell_for_clustering)
clust_df$.cell   <- rownames(ucell_for_clustering)
clust_df$.sample <- sample_var
clust_df$.study  <- study_var

# Compute per-study SD for each MP
study_sd <- clust_df %>%
  group_by(.study) %>%
  summarise(across(all_of(final_mps), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames(".study") %>%
  as.matrix()

# Replace 0 / NA SDs with 1 to avoid Inf
study_sd[is.na(study_sd) | study_sd == 0] <- 1

# Center within sample
clust_centered <- clust_df %>%
  group_by(.sample) %>%
  mutate(across(all_of(final_mps), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Scale by within-study SD: look up each cell's study → get that study's SD
mp_adj <- as.matrix(clust_centered[, final_mps])
rownames(mp_adj) <- clust_centered$.cell

for (mp in final_mps) {
  cell_studies <- clust_centered$.study
  mp_adj[, mp] <- mp_adj[, mp] / study_sd[cell_studies, mp]
}
mp_adj[!is.finite(mp_adj)] <- 0

# MPs × cells for Seurat assay
mp_matrix <- t(mp_adj)
# selected <- intersect(selected, rownames(mp_matrix))
# 
# mp_matrix <- mp_matrix[selected, , drop = FALSE]
# mp_adj    <- mp_adj[, selected, drop = FALSE]
# final_mps <- selected

# ── 4. Cluster ───────────────────────────────────────────────────────────────
cat("Clustering...\n")
orig_assay <- DefaultAssay(tmdata_all)

tmdata_all[["MPs"]] <- CreateAssayObject(data = mp_matrix)
DefaultAssay(tmdata_all) <- "MPs"

n_pcs <- min(30, nrow(mp_matrix) - 1)
tmdata_all <- ScaleData(tmdata_all, features = rownames(tmdata_all), verbose = FALSE)
tmdata_all <- RunPCA(tmdata_all, features = rownames(tmdata_all),
                     npcs = n_pcs, verbose = FALSE)
tmdata_all <- FindNeighbors(tmdata_all, reduction = "pca", dims = 1:n_pcs,
                            graph.name = "MPs_snn", verbose = FALSE)

# Store multiple resolutions for comparison
for (res in c(1)) {
  tmdata_all <- FindClusters(tmdata_all, graph.name = "MPs_snn",
                             resolution = res, verbose = FALSE)
  tmdata_all@meta.data[[paste0("MP_state_res", res)]] <-
    paste0("State_", as.numeric(as.character(tmdata_all$seurat_clusters)) + 1)
}

tmdata_all$MP_state <- tmdata_all@meta.data[[paste0("MP_state_res", CHOSEN_RES)]]
cat(sprintf("\nResolution %.1f → %d states\n", CHOSEN_RES,
            length(unique(tmdata_all$MP_state))))
print(table(tmdata_all$MP_state))

# ── 5. Study composition per state ──────────────────────────────────────────
comp_df <- tmdata_all@meta.data %>%
  count(MP_state, study) %>%
  group_by(MP_state) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(x = MP_state, y = pct, fill = study)) +
  geom_col(colour = "black", linewidth = 0.2) +
  labs(title = "Study composition per state (v3: study-Z → η² filter → CC regress → sample-center / study-SD)",
       x = "State", y = "% of cells", fill = "Study") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_comp

# ── 6. UMAP ─────────────────────────────────────────────────────────────────
tmdata_all <- RunUMAP(tmdata_all, reduction = "pca", dims = 1:n_pcs,
                      reduction.name = "umap_mp", verbose = FALSE)
p1 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "MP_state",
              label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("MP-based states (v3)")
p2 <- DimPlot(tmdata_all, reduction = "umap_mp", group.by = "study",
              pt.size = 0.1) +
  ggtitle("Coloured by study")
p1 | p2

# ── 7. Mean-score heatmap (RAW UCell scores for interpretability) ────────────
tmdata_all@meta.data <- cbind(tmdata_all@meta.data,
                              ucell_scores[rownames(tmdata_all@meta.data), final_mps, drop = FALSE])

# De-duplicate column names if cbind added duplicates
dup_cols <- duplicated(colnames(tmdata_all@meta.data))
if (any(dup_cols)) tmdata_all@meta.data <- tmdata_all@meta.data[, !dup_cols]

mean_scores <- tmdata_all@meta.data %>%
  select(MP_state, all_of(final_mps)) %>%
  group_by(MP_state) %>%
  summarise(across(all_of(final_mps), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("MP_state") %>%
  as.matrix()

colnames(mean_scores) <- gsub("\\.", " ", gsub("^X3CA_", "3CA_", colnames(mean_scores)))

mean_z    <- scale(mean_scores)
col_fun_z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mean <- Heatmap(
  t(mean_z),
  name = "Z-score",
  col  = col_fun_z,
  cluster_rows    = TRUE, cluster_columns = TRUE,
  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
  row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title    = "Mean MP enrichment per state (RAW scores, Z-scored per MP)",
  rect_gp         = gpar(col = "grey80", lwd = 0.3),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", t(mean_scores)[i, j]),
              x, y, gp = gpar(fontsize = 6))
  }
)

draw(ht_mean, merge_legend = TRUE)

# ── 8. Single-cell heatmap (adjusted scores, split by state) ────────────────
cat("Building single-cell heatmap...\n")
set.seed(42)
meta <- tmdata_all@meta.data
cells_to_plot <- unlist(
  lapply(split(rownames(meta), meta$MP_state),
         function(x) sample(x, min(length(x), MAX_CELLS_PER_STATE))),
  use.names = FALSE
)

sub_scores <- t(mp_adj[cells_to_plot, , drop = FALSE])
rownames(sub_scores) <- gsub("\\.", " ", gsub("^X3CA_", "3CA_", rownames(sub_scores)))

sub_meta  <- meta[cells_to_plot, ]
split_vec <- factor(sub_meta$MP_state, levels = sort(unique(sub_meta$MP_state)))
study_vals <- sub_meta$study

cell_order <- order(split_vec, apply(sub_scores, 2, which.max))
sub_scores <- sub_scores[, cell_order, drop = FALSE]
split_vec  <- split_vec[cell_order]
study_vals <- study_vals[cell_order]

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

col_ann <- HeatmapAnnotation(
  State = split_vec,
  Study = study_vals,
  col = list(
    State = setNames(scales::hue_pal()(length(levels(split_vec))), levels(split_vec)),
    Study = setNames(scales::hue_pal()(length(unique(study_vals))), unique(study_vals))
  ),
  annotation_name_side = "left"
)

# --- Add external annotation (state_temp) aligned to plotted cells ---
# state_temp: named character vector, names are cell IDs
state_temp <- setNames(as.character(state_temp), names(state_temp))
state_temp[state_temp == "Intest_Diff_Columnar"] <- "Intest_Metaplasia"
anno_vals <- state_temp[cells_to_plot]                 # align to sampled cells
anno_vals[is.na(anno_vals)] <- "Unassigned"            # fill missing

# reorder to match heatmap columns
anno_vals <- anno_vals[cell_order]

# make factor for stable legend ordering
anno_vals <- factor(anno_vals, levels = sort(unique(anno_vals)))

anno_cols <- c(
  "Classic_Prolif"        = "#E41A1C",  # Red
  "Squamous_Transition"   = "#4DAF4A",  # Green
  "Intest_Metaplasia"  = "#984EA3",  # Purple
  "Plastic_Tolerant"      = "#FF7F00",  # Orange
  #  "Intest_Diff"           = "#377EB8",  # Blue (Added)
  #  "IFN_columnar"          = "#A65628",  # Brown (Added)
  "Unresolved" = "grey80",   # Light Grey
  "Unassigned"            = "grey50"    # Darker Grey (Distinguishes noise from quienscent)
)

col_ann <- HeatmapAnnotation(
  State    = split_vec,
  Study    = study_vals,
  our_states = anno_vals,                                 # <-- NEW BAR
  col = list(
    State    = setNames(scales::hue_pal()(length(levels(split_vec))), levels(split_vec)),
    Study    = setNames(scales::hue_pal()(length(unique(study_vals))), unique(study_vals)),
    our_states = anno_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht_sc <- Heatmap(
  sub_scores, name = "Adj score", col = col_fun_sc,
  top_annotation = col_ann,
  column_split = split_vec, column_gap = unit(1.5, "mm"),
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE, row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  show_column_names = FALSE,
  use_raster = TRUE, raster_quality = 5,
  border = FALSE, rect_gp = gpar(col = NA),
  column_title = paste0("Single-cell MP heatmap (v3",
                        ifelse(REMOVE_CELL_CYCLE, "; CC removed & regressed", ""),
                        ")")
)

pdf("MP_state_selected.pdf", width = 18,
    height = max(7, length(final_mps) * 0.25), useDingbats = FALSE)
draw(ht_sc, merge_legend = TRUE)
dev.off()