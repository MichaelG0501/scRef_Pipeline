####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_manual_state_assignment.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/residual.R
# Purpose: Interactive state assignment — MP z-score residuals (regressing out cell-cycle MPs),
#          cosine-similarity-based manual state assignment for SC and PDO cells,
#          state coverage summaries, and ComplexHeatmap visualisations.
# Note: This is an interactive scratch script; objects (module_scores, tmdata_all,
#       tmdata_pdos, mod_mat, geneNMF.metaprograms) must be loaded in session.
####################

module_scores_raw <- as.matrix(module_scores)  # if module_scores currently holds raw
# (If module_scores already scaled in your session, reload the raw RDS first)

cycle_predictors <- intersect(c("MP1", "MP7"), colnames(module_scores_raw))
target_mps       <- intersect(c("MP3", "MP5", "MP6", "MP8", "MP9"), colnames(module_scores_raw))

X <- module_scores_raw[, cycle_predictors, drop = FALSE]

Y <- module_scores_raw[, target_mps, drop = FALSE]
Y_resid <- Y

for (mp in colnames(Y)) {
  fit <- lm(Y[, mp] ~ X)
  Y_resid[, mp] <- resid(fit)
}

module_scores[, colnames(Y_resid)] <- Y_resid
module_scores <- scale(as.matrix(module_scores))

############################

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
    state_name <- key#"Unassigned"
  }
  full_template_list[[key]] <- list(pattern = pattern, state = state_name)
}

# Convert to matrix for matching
all_templates <- t(sapply(full_template_list, function(x) x$pattern))
template_states <- sapply(full_template_list, function(x) x$state)

template_states["0,0,1,0,0"] <- template_states["0,0,0,1,0"] <- "Intestinal Metaplasia"
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
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Intestinal Metaplasia" = "#984EA3",
  "Stress-adaptive"       = "#FF7F00",
  "Unresolved"            = "grey80",
  "Unassigned"            = "grey50"
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

pdf("MP_heatmap_states_subset_z_residual.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(
  ht_manual,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

######################################

MP6_original <- scale(module_scores)[, "MP6"]
cycle_predictors <- c("MP1", "MP3", "MP4") 
fmla <- as.formula(paste("MP6 ~", paste(cycle_predictors, collapse = " + ")))
fit <- lm(fmla, data = as.data.frame(module_scores))
MP6_resid <- residuals(fit)

module_scores[, "MP6"] <- MP6_resid
module_scores <- scale(as.matrix(module_scores))

cell_cycle_signal <- rowMeans(module_scores[, c("MP1", "MP7")])
MP6_residual <- module_scores[, "MP6"]

# 3. Create Comparison Dataframe
# We combine Raw and Residual scores
score_df <- data.frame(
  Cell = rownames(module_scores),
  MP6_Raw = module_scores[, "MP6"],
  MP6_Residual = MP6_resid,
  Cycle_Signal = if("MP1" %in% colnames(module_scores)) module_scores[, "MP1"] else 0
)

# Merge with Metadata for Aggregation
score_df$Sample <- tmdata_pdos@meta.data[rownames(score_df), "orig.ident"]
score_df$State  <- tmdata_pdos@meta.data[rownames(score_df), "manual_state"]


# 1. Aggregate Scores by Sample
plot_data <- score_df %>%
  group_by(Sample) %>%
  summarise(
    Mean_MP6_Raw = mean(MP6_Raw, na.rm = TRUE),
    Mean_MP6_Resid = mean(MP6_Residual, na.rm = TRUE),
    Mean_Cycle = mean(Cycle_Signal, na.rm = TRUE)
  ) %>%
  as.data.frame()

rownames(plot_data) <- plot_data$Sample
plot_mat <- as.matrix(plot_data[, -1]) # Remove Sample Name column

# 2. Scale for Heatmap (Z-score)
# This makes the patterns comparable even if raw values differ in range
scaled_mat <- scale(plot_mat)

# 3. Plot Heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

Heatmap(scaled_mat,
        name = "Z-Score",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = FALSE, # Keep Raw -> Resid order
        
        # Split columns to separate Raw/Cycle from the Result (Residual)
        column_split = factor(c("Raw", "Result", "Raw"), levels = c("Raw", "Result")),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 8))
        },
        
        column_title = "Impact of Regressing Out Cycle (MP1/7) from MP6",
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        row_names_gp = gpar(fontsize = 8),
        border = TRUE
)


###

library(ComplexHeatmap)
library(circlize)

# 1. Prepare Data & Subset (60k limit for clustering)
# ---------------------------------------------------
# Create a matrix with human-readable column names
plot_mat <- cbind(
  "Original Signal" = scale(MP6_original),
  "Cell Cycle (Confounder)" = scale(cell_cycle_signal),
  "Refined Signal (Residual)" = scale(MP6_residual)
)

set.seed(123) # Ensure reproducible subset
if (nrow(plot_mat) > 60000) {
  subset_indices <- sample(seq_len(nrow(plot_mat)), 60000)
  plot_mat_subset <- plot_mat[subset_indices, ]
} else {
  plot_mat_subset <- plot_mat
}

# 2. Styling and Colors
# ---------------------
# Use a balanced Red-Blue divergence for Z-scores
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))

# Create a top annotation to categorize the columns
column_groups <- c("Raw Input", "Confounder", "Final Output")
ha_top <- HeatmapAnnotation(
  Type = column_groups,
  col = list(Type = c("Raw Input" = "grey80", 
                      "Confounder" = "#ffeda0", 
                      "Final Output" = "#a1d99b")),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.3, "cm")
)

# 3. Generate the Heatmap
# -----------------------
ht_resid <- Heatmap(
  plot_mat_subset,
  name = "Z-score",
  col = col_fun,
  
  # Clustering: Cluster rows to reveal patterns, but hide the tree
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  
  # Columns: Strict order (Input -> Confounder -> Output)
  cluster_columns = FALSE,
  top_annotation = ha_top,
  
  # Split the heatmap into three distinct blocks
  column_split = factor(column_groups, levels = column_groups),
  column_gap = unit(2, "mm"),
  
  # Labelling
  column_title = "Decomposition of MP6: Removing Cell Cycle Effects",
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Optimization for large data
  use_raster = TRUE,
  raster_quality = 5
)

# 4. Save to PDF
# --------------
pdf("Heatmap_temp2.pdf", width = 10, height = 8, useDingbats = FALSE)
draw(
  ht_resid,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


###################################

genes_mp1 <- geneNMF.metaprograms$metaprograms.genes[["MP1"]]
genes_mp3 <- geneNMF.metaprograms$metaprograms.genes[["MP3"]]
genes_mp4 <- geneNMF.metaprograms$metaprograms.genes[["MP4"]]
genes_mp6 <- geneNMF.metaprograms$metaprograms.genes[["MP6"]]
genes_to_exclude <- unique(c(genes_mp1, genes_mp3, genes_mp4))
selected_genes <- setdiff(genes_mp6, genes_to_exclude)
discarded_genes <- intersect(genes_mp6, genes_to_exclude)
gene_sets <- list(
  "Cell cycle genes MP6" = discarded_genes, 
  "Other genes MP6" = selected_genes
)
gene2module <- setNames(
  rep(names(gene_sets), times = lengths(gene_sets)),
  unlist(gene_sets, use.names = FALSE)
)
features <- names(gene2module)
DefaultAssay(tmdata_pdos) <- "RNA"
features_present <- features[features %in% rownames(tmdata_pdos)]

library(ComplexHeatmap)
library(circlize)
library(Seurat)
# Define the states you want to compare and their order
state_order <- c("Prolif. Columnar", "Classic Proliferative", "Stress Columnar", "Columnar", "Stress Plastic","Intestinal Diff.", "Unassigned/Quienscent", "Unassigned")

manual_names <- levels(tmdata_pdos$manual_state)
manual_cols <- c(
  "Prolif. Columnar"    = "#FB8072",
  "Classic Proliferative" = "#E41A1C",
  "Stress Columnar"     = "#BEBADA",
  "Columnar"            = "#984EA3",
  "Stress Plastic"      = "#FF7F00",
  "Intestinal Diff."    = "#377EB8",
  "Unassigned/Quienscent" = "grey80",
  "Unassigned"          = "grey50"
)


# --- 3. Prepare Matrix & Rows ---
# Extract expression matrix

features <- unlist(gene_sets)
features_present <- features[features %in% rownames(tmdata_pdos)]
mat <- as.matrix(GetAssayData(tmdata_pdos, slot = "data")[features_present, ])

# Create Row Split (Gene Modules)
# Assuming 'gene2module' exists and maps genes to module names
row_split <- factor(gene2module[features_present], levels = names(gene_sets))

# Define Color Function
col_fun <- colorRamp2(c(0, 1.5, 3, 6), c("#FCFDBF", "#FEB078", "#B73779", "#000004"))

# --- 4. Annotation (Top Bars) ---
top_anno <- HeatmapAnnotation(
  # Bar 1: The State (Main Grouping)
  State = tmdata_pdos$manual_state,
  
  # Bar 2: The Sample (To see patient composition within each state)
  Sample = tmdata_pdos$orig.ident,
  
  col = list(
    State = manual_cols
  ),
  
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm"),
  gap = unit(1, "mm")
)

# --- 5. Generate Heatmap ---
ht <- Heatmap(
  mat,
  name = "Log-Expr",
  col = col_fun,
  
  # --- KEY CHANGE: SPLIT BY STATE ---
  column_split = tmdata_pdos$manual_state, 
  column_gap = unit(2, "mm"), # Clean gap between states
  
  column_title_rot = 30, # Keep state names horizontal
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_side = "top",
  
  # Row Settings
  row_split = row_split,
  row_gap = unit(2.5, "mm"),
  cluster_rows = FALSE,      # Keep genes in your set order
  cluster_columns = TRUE,    # Cluster cells within each state
  show_column_dend = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  # Style
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_names_gp = gpar(fontsize = 9, fontitalic = TRUE),
  top_annotation = top_anno,
  border = TRUE, 
  use_raster = TRUE, 
  raster_quality = 5
)

# --- 6. Save ---
pdf("Heatmap_temp.pdf", width = 16, height = 8)
draw(ht, 
     merge_legends = TRUE, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm")) # Adjusted padding
dev.off()

###################################

# Start from raw UCell matrix
module_scores_raw <- as.matrix(module_scores)  # if module_scores currently holds raw
# (If module_scores already scaled in your session, reload the raw RDS first)

cycle_predictors <- intersect(c("MP1", "MP3", "MP4"), colnames(module_scores_raw))
target_mps       <- intersect(c("MP5", "MP6", "MP7", "MP8"), colnames(module_scores_raw))

X <- module_scores_raw[, cycle_predictors, drop = FALSE]

Y <- module_scores_raw[, target_mps, drop = FALSE]
Y_resid <- Y

for (mp in colnames(Y)) {
  fit <- lm(Y[, mp] ~ X)
  Y_resid[, mp] <- resid(fit)
}

module_scores[, colnames(Y_resid)] <- Y_resid
module_scores <- scale(as.matrix(module_scores))


#################

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Define your 4 known states (NA = either 0 or 1)
template_list <- rbind(
  "Stress Plastic"      = c(0, 0, 1, 0),
  "Stress Columnar"     = c(0, 1, 1, 0),
  "Stress Proliferative"= c(1, 0, 1, 0),
  "Classic Proliferative"= c(1, 0, 0, 0),
  "Columnar"            = c(0, 1, 0, 0),
  "Intestinal Diff."    = c(0, 0, 0, 1)
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
all_combos <- as.matrix(expand.grid(rep(list(0:1), 4)))
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
tmdata_pdos$manual_state <- as.vector(template_states[best_match_idx])

# Override: poor directional signal = Unassigned/Quiescent
max_score <- apply(sim_matrix, 1, max)
tmdata_pdos$manual_state[max_score < 0.1] <- "Unresolved"

max_score <- apply(cell_mat, 1, max)
tmdata_pdos$manual_state[max_score < 0.1] <- "Unresolved"


state_order <- c("Stress Columnar", "Classic Proliferative", "Stress Proliferative", "Columnar", "Stress Plastic","Intestinal Diff.", "Unresolved", "Unassigned")

tmdata_pdos$manual_state <- factor(tmdata_pdos$manual_state, levels = state_order)

manual_names <- levels(tmdata_pdos$manual_state)
manual_cols <- c(
  "Stress Columnar"     = "green",
  "Classic Proliferative" = "#E41A1C",
  "Stress Proliferative" = "#FF4D4D",
  "Columnar"            = "#984EA3",
  "Stress Plastic"      = "#FF7F00",
  "Intestinal Diff."    = "#377EB8",
  "Unresolved"          = "grey80",
  "Unassigned"          = "grey50"
)

sample_col <- "orig.ident"
study_col  <- "Batch" # Ensure this matches your metadata column name

# 1. Basic Metadata Setup
meta <- tmdata_pdos@meta.data
total_samples_n <- length(unique(meta[[sample_col]]))
total_studies_n <- length(unique(meta[[study_col]]))

# Create a mapping of Sample to Study (needed for counting studies later)
sample_to_study <- unique(meta[, c(sample_col, study_col)])
rownames(sample_to_study) <- sample_to_study[[sample_col]]

# 2. Create State x Sample contingency table
counts_matrix <- table(tmdata_pdos$manual_state, meta[[sample_col]])

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

max_val <- quantile(mod_mat, 0.98, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_val), c("white", "red"))

patient_levels <- unique(tmdata_pdos$SUR)
patient_cols <- setNames(
  DiscretePalette(length(patient_levels), palette = "alphabet"), 
  patient_levels
)
tmdata_pdos$batch <- ifelse(tmdata_pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia")
batch_names  <- as.character(unique(tmdata_pdos$batch))
batch_cols <- setNames(
  DiscretePalette(length(batch_names), palette = "polychrome"),
  batch_names
)

col_ann_manual <- HeatmapAnnotation(
  State = tmdata_pdos$manual_state,
  batch = tmdata_pdos$batch,
  Patient = tmdata_pdos$SUR, 
  col = list(State = manual_cols, batch = batch_cols, Patient = patient_cols),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht_manual <- Heatmap(
  mod_mat,
  name = "MPs UCell scores clustering by states",
  col = col_fun,
  
  # split in fixed order (factor levels) + add gaps between states
  top_annotation = col_ann_manual,
  column_split = tmdata_pdos$manual_state,
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

pdf("MP_heatmap_states_subset_z_residual.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(
  ht_manual,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()