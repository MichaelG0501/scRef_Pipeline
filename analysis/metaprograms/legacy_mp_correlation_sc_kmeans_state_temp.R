####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_mp_correlation_sc_kmeans_state_temp.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/MP_analysis_sc.R
# Reorganized as part of analysis/ restructuring
####################
library(Seurat)
library(ggplot2)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
module_scores <- readRDS("UCell_default.rds")
module_scores <- scale(as.matrix(module_scores))
geneNMF.metaprograms <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/MP_outs_default.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes

mp.genes$MP2 <- NULL
mp.genes$MP4 <- NULL
module_scores <- module_scores[, -c(2,4)]

mp.genes$MP1 <- NULL
mp.genes$MP7 <- NULL
module_scores <- module_scores[, c("MP3", "MP5", "MP6", "MP8", "MP9"), drop = FALSE]

tmdata_all@meta.data[, grepl("^MP", colnames(tmdata_all@meta.data))] <- NULL
tmdata_all@meta.data <- cbind(tmdata_all@meta.data, module_scores)


matrix <- tmdata_all@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
tmdata_all@reductions[["MPsignatures"]] <- new("DimReduc",
                                               cell.embeddings = dimred,
                                               assay.used = "RNA",
                                               key = "MP_",
                                               global = FALSE)
tmdata_all <- RunUMAP(tmdata_all, reduction="MPsignatures", dims=1:length(tmdata_all@reductions[["MPsignatures"]]),
                      metric = "euclidean", reduction.name = "umap_MP")
DimPlot(tmdata_all, reduction = "umap_MP", group.by = "Treatment") + theme(aspect.ratio = 1)

mp.emb <- Embeddings(tmdata_all, "MPsignatures")   # cells x 10 (or whatever)

set.seed(1)
K <- 500   # pick e.g. 200–2000 depending on how fine you want “micro-clusters”
km <- kmeans(mp.emb, centers = K, nstart = 10, iter.max = 200, algorithm = "MacQueen")
saveRDS(km, "km_500_subset_z.rds")

clus <- factor(km$cluster)
centroids <- km$centers   # K x 10 matrix (already the centroids)

d <- as.dist(1 - cor(t(centroids), method = "pearson"))
hc <- hclust(d, method = "average")

# after kmeans:
tmdata_all$km_cluster <- as.character(km$cluster)
library(factoextra)
library(cluster)
p1 <- fviz_nbclust(centroids, hcut, hc_method = "average", hc_metric = "pearson", 
                   method = "silhouette", k.max = 20) +
  labs(title = "Silhouette Analysis (Centroids)")
p2 <- fviz_nbclust(centroids, hcut, hc_method = "average", hc_metric = "pearson", 
                   method = "wss", k.max = 20) +
  labs(title = "Elbow Method (Centroids)")
library(patchwork)
p1 + p2

k <- 4
# name mc_state by kmeans cluster IDs (1..K)
mc_state <- cutree(hc, k = k)
names(mc_state) <- rownames(centroids)  # rownames(centroids) are "1","2",... by default

tmdata_all$MP_state <- as.character(mc_state[tmdata_all$km_cluster])

# --- 2. Define Row Order from GeneNMF ---
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Ensure naming matches module_scores (adjust "_" if your names use "MP1" instead of "MP_1")
mp_tree_order <- paste0("MP", mp_tree_order) 
# Remove the 10th element as per your requirement
mp_tree_order <- mp_tree_order[-8:-10]
mp_tree_order <- mp_tree_order[-1:-2]

# --- 3. Cell Ordering Logic (Within States) ---
rownames(centroids) <- as.character(seq_len(nrow(centroids)))  # if not already set

cell_clusters <- tmdata_all$km_cluster  # <-- use kmeans ids, not seurat_clusters
mc_order <- rownames(centroids)[hc$order]

cell_scores <- Embeddings(tmdata_all, "MPsignatures")

centroid_mat <- centroids[cell_clusters, , drop = FALSE]  # align per-cell centroid
cell_to_centroid_cor <- vapply(seq_len(nrow(cell_scores)), function(i) {
  cor(cell_scores[i, ], centroid_mat[i, ], method = "pearson")
}, numeric(1))

order_df <- data.frame(
  cell_id = colnames(tmdata_all),
  state = tmdata_all$MP_state,
  state_order = match(tmdata_all$MP_state, unique(mc_state[hc$order])),
  mc_order_val = match(cell_clusters, mc_order),
  cor_val = cell_to_centroid_cor
)

# Sort cells: State > Cluster Tree > Similarity to Centroid
order_df <- order_df[order(order_df$state_order, order_df$mc_order_val, -order_df$cor_val), ]
cell_ord <- match(order_df$cell_id, colnames(tmdata_all))

# Reorder the Seurat object
tmdata_all <- tmdata_all[, cell_ord]

# --- 4. Prepare Heatmap Matrix ---
# Subset and order the matrix rows by the MP tree
mod_mat <- t(as.matrix(tmdata_all@meta.data[, mp_tree_order, drop = FALSE]))
# Create a mapping vector based on your table
mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC Proliferation",
  "MP17" = "Basal-like Trans.",
  "MP14" = "Hypoxia Adapted",
  "MP5"  = "Epithelial IFN",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Intest.",
  "MP16" = "Secretory Gastric",
  "MP15" = "Immune Infilt.",
  "MP12" = "Neuro-responsive"
)
rownames(mod_mat) <- mp_descriptions[rownames(mod_mat)]

# --- 5. Visualization with ComplexHeatmap ---
library(ComplexHeatmap)
library(circlize)
library(Seurat)

# Color scale (White to Red)
max_val <- quantile(mod_mat, 0.98, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_val), c("white", "red"))

# Discrete Annotation Colors
state_names <- as.character(unique(tmdata_all$MP_state))
state_cols  <- setNames(DiscretePalette(length(state_names), palette = "alphabet"), state_names)
study_cols  <- setNames(DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"), unique(tmdata_all$study))

# Define Annotations
col_ann <- HeatmapAnnotation(
  State = tmdata_all$MP_state,
  Study = tmdata_all$study,
  col = list(State = state_cols, Study = study_cols),
  annotation_name_side = "left",
  show_legend = TRUE
)

# Plot
ht <- Heatmap(
  mod_mat,
  name = "UCell Score",
  col = col_fun,
  
  # Grouping
  top_annotation = col_ann,
  column_split = tmdata_all$MP_state,
  column_gap = unit(0, "mm"),
  
  # Row/Column Ordering
  cluster_rows = FALSE,         # Keep GeneNMF tree order
  row_order = rownames(mod_mat),
  cluster_columns = TRUE,       # Clusters cells LOCALLY within each MP_state split
  clustering_method_columns = "ward.D2",
  
  # Aesthetics
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = "MPs UCell scores clustering k = 8",
  show_column_names = FALSE,
  
  # Performance
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf("MP_heatmap_subset_z.pdf", width = 18, height = 10, useDingbats = FALSE)
ht2 <- draw(ht, merge_legend = TRUE)
ord_list <- column_order(ht2)     # usually a list (one per slice)
ord <- unlist(ord_list)
true_final_state_order <- unique(as.character(tmdata_all$MP_state[ord]))
dev.off()

library(ComplexHeatmap)
library(circlize)

# --- 1. Setup Data & Variables ---
studies <- unique(tmdata_all$study)
mps <- rownames(mod_mat)
n_mps <- length(mps)

# Create 3D array to store correlations: [MP x MP x Study]
cor_array <- array(0, dim = c(n_mps, n_mps, length(studies)), 
                   dimnames = list(mps, mps, studies))

# --- 2. Loop per Study (Compute 8 Matrices) ---
for (st in studies) {
  # Subset matrix for this study
  # Note: mod_mat is [MPs x Cells], verify column alignment with study vector
  cells_in_study <- colnames(mod_mat)[tmdata_all$study == st]
  sub_mat <- mod_mat[, cells_in_study, drop=FALSE]
  
  # Compute Spearman correlation
  cor_array[,,st] <- cor(t(sub_mat), method = "spearman")
}

# --- 3. Meta-Analysis (Compute Mean Rho & Significance) ---
# We use Fisher Z-transformation for accurate averaging and t-testing
z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999)) # Clip to avoid Inf

# Initialize result matrices
mean_rho <- matrix(0, n_mps, n_mps, dimnames = list(mps, mps))
p_vals   <- matrix(1, n_mps, n_mps, dimnames = list(mps, mps))

for (i in 1:n_mps) {
  for (j in 1:n_mps) {
    if (i == j) {
      mean_rho[i,j] <- 1
      p_vals[i,j] <- 0
    } else {
      # Extract the 8 Z-scores for this MP pair
      z_scores <- z_array[i, j, ]
      
      # T-test: Is the mean Z-score significantly different from 0?
      test_res <- t.test(z_scores)
      
      # Convert mean Z back to Rho for plotting
      mean_rho[i,j] <- tanh(mean(z_scores))
      p_vals[i,j]   <- test_res$p.value
    }
  }
}

# --- 4. Plot Beautiful Heatmap ---
library(ComplexHeatmap)
library(circlize)

col_cor <- colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")) # Adjusted scale for visibility

Heatmap(mean_rho,
        name = "Mean Rho\n(8 Studies)",
        col = col_cor,
        rect_gp = gpar(col = "white", lwd = 1), # Clean white borders
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        
        # Overlay Significance Stars
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- p_vals[i, j]
          cor <- mean_rho[i, j]
          # If significant, add text
          if(p < 0.001) {
            grid.text(paste0(round(cor, 2), "\n***"), x, y, gp = gpar(fontsize = 10))
          } else if(p < 0.01) {
            grid.text(paste0(round(cor, 2), "\n**"), x, y, gp = gpar(fontsize = 10))
          } else if(p < 0.05) {
            grid.text(paste0(round(cor, 2), "\n*"), x, y, gp = gpar(fontsize = 10))
          } else {
            grid.text(round(cor, 2), x, y, gp = gpar(fontsize = 10))
          }
        },
        
        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"))
)

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- mean MP scores per state ----
state_mp_means <- tmdata_all@meta.data %>%
  filter(!is.na(MP_state)) %>%
  group_by(MP_state) %>%
  summarise(
    across(all_of(mp_tree_order), mean, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# counts (named by state)
cell_counts <- state_mp_means$n_cells
names(cell_counts) <- as.character(state_mp_means$MP_state)

# matrix: MPs as rows, states as columns
state_mp_means_df <- as.data.frame(state_mp_means)
rownames(state_mp_means_df) <- as.character(state_mp_means_df$MP_state)
state_mp_means_df$MP_state <- NULL
state_mp_means_df$n_cells <- NULL

summary_mat <- t(as.matrix(state_mp_means_df))  # rows=MPs, cols=states

# ---- order rows + columns ----
summary_mat <- summary_mat[mp_tree_order, , drop = FALSE]
rownames(summary_mat) <- mp_descriptions[rownames(summary_mat)]

# true_final_state_order must match column names (states)
true_final_state_order <- intersect(true_final_state_order, colnames(summary_mat))
summary_mat <- summary_mat[, true_final_state_order, drop = FALSE]

# ---- colors ----
max_mean_val <- max(summary_mat, na.rm = TRUE)
col_fun_mean <- colorRamp2(c(0, max_mean_val), c("white", "red"))

# ---- top annotation: n_cells ----
state_cell_counts <- cell_counts[true_final_state_order]

col_ann_summary <- HeatmapAnnotation(
  n_cells = anno_barplot(state_cell_counts, height = unit(2, "cm")),
  annotation_name_side = "left"
)

# ---- plot ----
ht_summary <- Heatmap(
  summary_mat,
  name = "Mean UCell\nScore",
  col = col_fun_mean,
  top_annotation = col_ann_summary,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", summary_mat[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  
  column_title = "Average Meta-Program Expression per State",
  border = TRUE,
  rect_gp = gpar(col = "grey30", lwd = 1)
)

draw(ht_summary, heatmap_legend_side = "right")

#########################

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
# Define your 4 known states (NA = either 0 or 1)
template_list <- rbind(
  "Intestinal Metaplasia" = c(0, NA, 1, 1, 0),
  "Stress-adaptive"       = c(0, 0, 0, 0, 1),
  "Classic Proliferative" = c(1, 0, 0, 0, 0),
  "Basal to Intest. Meta" = c(0, 1, NA, 0, 0)
)

colnames(template_list) <- rownames(mod_mat)
plot_mat <- as.matrix(template_list)
col_logic <- c("1" = "red", "0" = "white")
ht_template <- Heatmap(
  t(plot_mat), # Transpose to keep MPs as rows and States as columns
  name = "Logic",
  col = col_logic,
  
  # Handle NAs explicitly with a grey color
  na_col = "grey80",
  
  # Grid and Style
  rect_gp = gpar(col = "grey30", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Text Labels
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  # Add text overlay to confirm values
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- t(plot_mat)[i, j]
    txt <- ifelse(is.na(v), "NA", as.character(v))
    grid.text(txt, x, y, gp = gpar(fontsize = 10, 
                                   col = ifelse(!is.na(v) && v == 1, "white", "black")))
  },
  
  column_title = "Manual State Assignment Templates",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Legend Customization
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Logic",
    at = c(0, 1),
    labels = c("Low (0)", "High (1)"),
    color_bar = "discrete"
  )
)
draw(ht_template)

library(proxy)

# Expand templates with NA to all explicit combinations they represent
expand_template <- function(template) {
  na_idx <- which(is.na(template))
  if (length(na_idx) == 0) {
    return(matrix(template, nrow = 1))
  }
  # Generate all 2^n combinations for NA positions
  n_na <- length(na_idx)
  combos <- as.matrix(expand.grid(rep(list(0:1), n_na)))
  expanded <- matrix(rep(template, nrow(combos)), nrow = nrow(combos), byrow = TRUE)
  for (i in seq_along(na_idx)) {
    expanded[, na_idx[i]] <- combos[, i]
  }
  expanded
}

# Get all explicit patterns covered by defined states
defined_patterns <- list()
for (state in rownames(template_list)) {
  expanded <- expand_template(template_list[state, ])
  for (i in 1:nrow(expanded)) {
    key <- paste(expanded[i, ], collapse = ",")
    defined_patterns[[key]] <- state
  }
}

# Generate ALL 32 possible binary combinations (2^5)
all_combos <- as.matrix(expand.grid(rep(list(0:1), 5)))
colnames(all_combos) <- rownames(mod_mat)

# Build full template list including "Unassigned" for undefined patterns
full_template_list <- list()
for (i in 1:nrow(all_combos)) {
  pattern <- all_combos[i, ]
  key <- paste(pattern, collapse = ",")
  if (key %in% names(defined_patterns)) {
    state_name <- defined_patterns[[key]]
  } else {
    state_name <- "Unassigned"
  }
  full_template_list[[key]] <- list(pattern = pattern, state = state_name)
}

# Convert to matrix for matching
all_templates <- t(sapply(full_template_list, function(x) x$pattern))
template_states <- sapply(full_template_list, function(x) x$state)

template_states["0,0,1,0,0"] <- template_states["0,0,0,1,0"] <- "Intest_Diff_Columnar" 

# Prepare cell matrix
cell_mat <- t(mod_mat)

# Standard cosine similarity (no masking needed - all templates are explicit now)
cosine_similarity <- function(X, tvec) {
  if (all(tvec == 0)) return(rep(0, nrow(X)))
  
  num <- as.vector(X %*% tvec)
  den <- sqrt(rowSums(X^2)) * sqrt(sum(tvec^2))
  out <- num / den
  out[is.na(out)] <- 0
  out
}

# Compute similarity to all 32 explicit templates
sim_matrix <- sapply(1:nrow(all_templates), function(i) {
  cosine_similarity(cell_mat, all_templates[i, ])
})

# Assign best matching template per cell
best_match_idx <- max.col(sim_matrix, ties.method = "first")
tmdata_all$manual_state <- as.vector(template_states[best_match_idx])

# Override: poor directional signal = Unassigned/Quiescent
max_score <- apply(sim_matrix, 1, max)
tmdata_all$manual_state[max_score < 0.1] <- "Unresolved"

max_score <- apply(cell_mat, 1, max)
tmdata_all$manual_state[max_score < 0.1] <- "Unresolved"


state_order <- c("Classic Proliferative", "Basal to Intest. Meta", "Intestinal Metaplasia", "Stress-adaptive", "Unresolved", "Unassigned")

tmdata_all$manual_state <- factor(tmdata_all$manual_state, levels = state_order)

manual_names <- levels(tmdata_all$manual_state)
manual_cols <- c(
  "Classic Proliferative" = "#E41A1C",  # Red
  "Basal to Intest. Meta" = "#4DAF4A",  # Green
  "Intestinal Metaplasia" = "#984EA3",  # Purple
  "Stress-adaptive"       = "#FF7F00",  # Orange
  "Unresolved"            = "grey80",   # Light Grey
  "Unassigned"            = "grey50"    # Darker Grey
)

sample_col <- "orig.ident"
study_col  <- "study" # Ensure this matches your metadata column name

# 1. Basic Metadata Setup
meta <- tmdata_all@meta.data
total_samples_n <- length(unique(meta[[sample_col]]))
total_studies_n <- length(unique(meta[[study_col]]))

# Create a mapping of Sample to Study (needed for counting studies later)
sample_to_study <- unique(meta[, c(sample_col, study_col)])
rownames(sample_to_study) <- sample_to_study[[sample_col]]

# 2. Create State x Sample contingency table
counts_matrix <- table(tmdata_all$manual_state, meta[[sample_col]])

# 3. Calculate Robust Sample Coverage (> 10 cells)
# Logical matrix: TRUE if state has > 10 cells in that sample
robust_mask <- counts_matrix > 10
robust_sample_coverage <- rowSums(robust_mask)

# 4. Calculate Robust Study Coverage
# A study is "robust" for a state if at least ONE sample in that study has > 10 cells
robust_study_coverage <- apply(robust_mask, 1, function(row) {
  # Get the names of samples that passed the > 10 threshold
  robust_samples <- colnames(robust_mask)[row]
  # Count unique studies associated with those samples
  length(unique(sample_to_study[robust_samples, study_col]))
})

# 5. Create the Final Summary Table
coverage_summary <- data.frame(
  State = rownames(counts_matrix),
  Total_Cell_Count = as.numeric(rowSums(counts_matrix)),
  Sample_N = as.numeric(robust_sample_coverage),
  Sample_Pct = round((as.numeric(robust_sample_coverage) / total_samples_n) * 100, 1),
  Study_N = as.numeric(robust_study_coverage),
  Study_Pct = round((as.numeric(robust_study_coverage) / total_studies_n) * 100, 1)
)

coverage_summary <- coverage_summary[order(-coverage_summary$Total_Cell_Count), ]

library(circlize)
max_val <- quantile(mod_mat, 0.98, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_val), c("white", "red"))
study_cols  <- setNames(DiscretePalette(length(unique(tmdata_all$study)), palette = "polychrome"), unique(tmdata_all$study))

col_ann_manual <- HeatmapAnnotation(
  State = tmdata_all$manual_state,
  Study = tmdata_all$study,
  col = list(State = manual_cols, Study = study_cols),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht_manual <- Heatmap(
  mod_mat,
  name = "MPs UCell scores clustering by states",
  col = col_fun,
  
  # split in fixed order (factor levels) + add gaps between states
  top_annotation = col_ann_manual,
  column_split = tmdata_all$manual_state,
  column_gap = unit(2, "mm"),        # <-- gap between state blocks (adjust)
  
  # keep rows fixed, cluster columns WITHIN each state
  cluster_rows = FALSE,
  row_order = rownames(mod_mat),
  
  cluster_columns = TRUE,            # cluster within slices
  clustering_method_columns = "ward.D2",
  cluster_column_slices = FALSE,     # <-- IMPORTANT: do NOT reorder state slices
  
  # aesthetics like your other heatmap
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = "MPs UCell scores clustering by states",
  show_column_names = FALSE,
  
  # performance + clean look
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf("MP_heatmap_states_subset_z_tweak3.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(
  ht_manual,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- mean MP scores per state ----
state_mp_means <- tmdata_all@meta.data %>%
  filter(!is.na(manual_state)) %>%
  group_by(manual_state) %>%
  summarise(
    across(all_of(mp_tree_order), mean, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# counts (named by state)
cell_counts <- state_mp_means$n_cells
names(cell_counts) <- as.character(state_mp_means$manual_state)

# matrix: MPs as rows, states as columns
state_mp_means_df <- as.data.frame(state_mp_means)
rownames(state_mp_means_df) <- as.character(state_mp_means_df$manual_state)
state_mp_means_df$manual_state <- NULL
state_mp_means_df$n_cells <- NULL

summary_mat <- t(as.matrix(state_mp_means_df))  # rows=MPs, cols=states

# ---- order rows + columns ----
summary_mat <- summary_mat[mp_tree_order, , drop = FALSE]

# ---- colors ----
max_mean_val <- max(summary_mat, na.rm = TRUE)
col_fun_mean <- colorRamp2(c(0, max_mean_val), c("white", "red"))

# ---- top annotation: n_cells ----
state_cell_counts <- cell_counts

col_ann_summary <- HeatmapAnnotation(
  n_cells = anno_barplot(state_cell_counts, height = unit(2, "cm")),
  annotation_name_side = "left"
)

# ---- plot ----
ht_summary <- Heatmap(
  summary_mat,
  name = "Mean UCell\nScore",
  col = col_fun_mean,
  top_annotation = col_ann_summary,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", summary_mat[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  
  column_title = "Average Meta-Program Expression per State",
  border = TRUE,
  rect_gp = gpar(col = "grey30", lwd = 1)
)

draw(ht_summary, heatmap_legend_side = "right")



library(dplyr)

tmdata_all$Treatment_B <- ifelse(tmdata_all$Treatment == "Post", "Post", "Pre")
tmdata_all$Treatment_B <- factor(tmdata_all$Treatment_B, levels = c("Pre", "Post"))

library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)

# 1. Prepare and Filter Metadata
comp_df <- tmdata_all@meta.data %>%
  # Basic cleaning
  filter(!is.na(manual_state), !is.na(Treatment_B), !is.na(orig.ident), !is.na(study)) %>%
  # CRITICAL FILTER: Keep only studies that have both "Pre" and "Post"
  # group_by(study) %>%
  # filter(all(c("Pre", "Post") %in% Treatment_B)) %>%
  # ungroup() %>%
  
  # 2. Calculate Sample-level proportions (Avoid cell abundance bias)
  count(orig.ident, Treatment_B, manual_state, study) %>%
  group_by(orig.ident) %>%
  mutate(sample_prop = n / sum(n)) %>%
  ungroup() %>%
  
  # 3. Average proportions across Treatment groups
  group_by(Treatment_B, manual_state) %>%
  summarize(mean_prop = mean(sample_prop), .groups = 'drop') %>%
  
  # 4. Re-normalize to 100% for visualization
  group_by(Treatment_B) %>%
  mutate(prop = mean_prop / sum(mean_prop)) %>%
  ungroup()

# 5. Plot
ggplot(comp_df, aes(x = Treatment_B, y = prop, fill = manual_state)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = manual_cols, drop = FALSE) +
  labs(
    y = "Cell Fraction (Mean of Samples)", 
    x = "", 
    fill = "Malignant States",
  ) +
  theme_classic(base_size = 18) +   # 👈 increases ALL text
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.height = unit(0.9, "cm")
  )


ggplot(comp_df, aes(x = Treatment_B, y = prop, color = Treatment_B)) +
  geom_point(size = 3) +
  geom_line(aes(group = manual_state), color = "grey70") +
  facet_wrap(~ manual_state, scales = "free_y") +
  labs(
    y = "Cell fraction",
    title = "MP-state proportions per treatment"
  ) +
  theme_classic(base_size = 13)


library(ggpubr) # Required for stat_compare_means
comp_sample <- tmdata_all@meta.data %>%
  filter(!is.na(manual_state), !is.na(Treatment_B)) %>%
  count(orig.ident, Treatment_B, manual_state, study) %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>% filter(manual_state == "Intest_Metaplasia")

# comp_sample <- comp_sample %>%
#   group_by(study) %>%
#   filter(all(c("Pre", "Post") %in% Treatment_B)) %>%
#   ungroup()

strip_states <- unique(comp_sample$manual_state)
strip_bg_colors <- manual_cols[strip_states]

ggplot(comp_sample, aes(x = Treatment_B, y = prop, fill = Treatment_B)) +
  stat_summary(fun = "mean", geom = "bar", color = "black", alpha = 0.7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  facet_wrap(~ manual_state, scales = "free_y") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    label.x = 1.45,
    size = 8,
    fontface = "bold",
    vjust = 0.5,
    hide.ns = FALSE
  ) +
  scale_fill_manual(values = c("Pre" = "#1B9E77", "Post" = "#7F0000")) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(
    y = "Cell fraction (%)",
    x = NULL,
    title = "State Proportions Pre vs Post"
  ) +
  theme_classic(base_size = 19) +
  theme(
    strip.background = element_rect(fill = strip_bg_colors, color = "black"),
    strip.text = element_text(face = "bold", color = "white", size = 16),
    legend.position = "none"
  )
  

# Example of how to structure the data if not already done:
comp_df <- tmdata_all@meta.data %>%
  group_by(orig.ident, Treatment_B, manual_state, study) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n))


library(ggplot2)
library(dplyr)

# 1. Filter for the specific state and ensure we are at the sample level
plot_data <- comp_df %>% 
  filter(manual_state == "Intest_Metaplasia")

# 2. Plot
ggplot(plot_data, aes(x = Treatment_B, y = prop, fill = Treatment_B)) +
  # Violin shows the density/shape of the distribution
  geom_violin(alpha = 0.3, color = NA) +
  # Boxplot shows median and quartiles
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  # Jitter shows every individual orig.ident
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.8, aes(color = Treatment_B)) +
  
  # Split into panels by Study
  facet_wrap(~ study, scales = "free_y") + 
  
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = "Cell fraction (%)",
    x = NULL,
    title = "Proportion of MP-State 1",
    subtitle = "Broken down by Study and Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
