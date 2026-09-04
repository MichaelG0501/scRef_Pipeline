####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_mp_pdo_sc_crossdataset_correlation.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/compare_pdos_sc.R
# Reorganized as part of analysis/ restructuring
####################
signature_A <- c(
  "ACOT7", "ADM", "AK4", "ALDOA", "ANKRD37", "ANLN", "BNIP3", "MRGBP", "CA9", "CDKN3",
  "CHCHD2", "CORO1C", "CTSV", "DDIT4", "ENO1", "ESRP1", "GAPDH", "GPI", "HILPDA", "HK2",
  "KIF20A", "KIF4A", "LDHA", "LRRC42", "MAD2L2", "MAP7D1", "MCTS1", "MIF", "MRPL13", "MRPL15",
  "MRPS17", "NDRG1", "ZNF384", "P4HA1", "PFKP", "PGAM1", "PGK1", "PSMA7", "PSRC1", "SEC61G",
  "SHCBP1", "SLC16A1", "SLC25A32", "SLC2A1", "TPI1", "TUBA1B", "TUBA1C", "TUBB6", "UTP11",
  "VEGFA", "YKT6"
)

signature_A <- intersect(signature_A, rownames(tmdata_pdos))

signature_B <- c(
  "ASF1B", "ADM", "ASPM", "BIRC5", "BUB3", "CENPE", "CENPU", "CMTM3", "DDIT4", "DONSON",
  "DSP", "DTL", "FER1L4", "FOXM1", "G6PD", "HILPDA", "HJURP", "MCM2", "MEP1A", "MTMR2",
  "P4HA1", "PGAM4", "PKM", "RIMKLA", "RNASE4", "SCD", "SPAG4", "TDG", "TRIP13", "UNG",
  "XRCC6", "ZWINT"
)

signature_B <- intersect(signature_B, rownames(tmdata_pdos))

# Interquartile mean function
iqm <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  mean(x[x >= q[1] & x <= q[2]], na.rm = TRUE)
}

score_A <- colMeans(tmdata_pdos@assays$RNA$data[signature_A, , drop = FALSE])             # mean
score_B <- apply(tmdata_pdos@assays$RNA$data[signature_B, , drop = FALSE], 2, iqm)        # iqm

# Identify MP columns
mp_cols <- grep("^MP", colnames(tmdata_pdos@meta.data), value = TRUE)

# Pull MP columns and force them numeric
mp_meta <- tmdata_pdos@meta.data[, mp_cols, drop = FALSE]
mp_meta[] <- lapply(mp_meta, function(x) as.numeric(as.character(x)))

# Combine
meta <- cbind(mp_meta, score_A = as.numeric(score_A), score_B = as.numeric(score_B))
cor_matrix <- cor(meta[, mp_cols], meta[, "score_A"], method = "spearman")
ordered_cor_matrix <- cor_matrix[pdo_mp_tree_order, , drop = FALSE]

# Define symmetric breaks around 0
max_abs <- max(abs(ordered_cor_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(
  ordered_cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = breaks,          # <- ensures 0 is white
  display_numbers = TRUE,
  fontsize = 14,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

cor_matrix <- cor(meta[, mp_cols], meta[, "score_B"], method = "spearman")
ordered_cor_matrix <- cor_matrix[pdo_mp_tree_order, , drop = FALSE]
# Define symmetric breaks around 0
max_abs <- max(abs(ordered_cor_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
pheatmap(
  ordered_cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = breaks,          # <- ensures 0 is white
  display_numbers = TRUE,
  fontsize = 14,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

#######################################

cell_cycle_genes <- read.csv("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
                             header = TRUE, stringsAsFactors = FALSE)[, 1:3]
cc_genes <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_genes <- intersect(cc_genes, rownames(tmdata_all))
gene_means <- rowMeans(tmdata_all@assays$RNA$data[cc_genes, , drop = FALSE], na.rm = TRUE)
cc_genes <- names(sort(gene_means, decreasing = TRUE))[1:50]

cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_genes, , drop = FALSE]))
cc_score <- cc_score[rownames(tmdata_all@meta.data)]   # align to meta.data rows

# Identify MP columns
mp_cols <- grep("^MP", colnames(tmdata_all@meta.data), value = TRUE)

# Pull MP columns and force them numeric
mp_meta <- tmdata_all@meta.data[, mp_cols, drop = FALSE]
mp_meta[] <- lapply(mp_meta, function(x) as.numeric(as.character(x)))

# Combine
meta <- cbind(mp_meta, cc_score = as.numeric(cc_score))
cor_matrix <- cor(meta[, mp_cols], meta[, "cc_score"], method = "spearman")
ordered_cor_matrix <- cor_matrix[sc_mp_tree_order, , drop = FALSE]

# Define symmetric breaks around 0
max_abs <- max(abs(ordered_cor_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(
  ordered_cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = breaks,          # <- ensures 0 is white
  display_numbers = TRUE,
  fontsize = 14,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

#####################

MP_pdo <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/MP_outs_default.rds")
MP_sc  <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/MP_outs_default.rds")

### 1) Extract lists
pdo_list <- MP_pdo$metaprograms.genes
sc_list  <- MP_sc$metaprograms.genes

pdo_list <- lapply(pdo_list, unique)
sc_list  <- lapply(sc_list, unique)

### 2) Derive PDO order from its own tree
pdo_tree_order <- MP_pdo$programs.tree$order
pdo_ordered_clusters <- MP_pdo$programs.clusters[pdo_tree_order]
pdo_mp_tree_order <- rev(unique(pdo_ordered_clusters))[c(-5,-8,-9)]
pdo_mp_tree_order <- paste0("MP", pdo_mp_tree_order) 

### 3) Derive sc order from its own tree
sc_tree_order <- MP_sc$programs.tree$order
sc_ordered_clusters <- MP_sc$programs.clusters[sc_tree_order]
sc_mp_tree_order <- unique(sc_ordered_clusters)[-8:-10]
sc_mp_tree_order <- paste0("MP", sc_mp_tree_order) 

### 4) Reorder lists
pdo_list <- pdo_list[pdo_mp_tree_order]
sc_list  <- sc_list[sc_mp_tree_order]

pdo_mp_descriptions <- c(
  "MP3"  = "G2M Cell Cycle",
  "MP6"  = "MYC Proliferation",
  "MP1"  = "G2M Checkpoint",
  "MP4"  = "G1S Cell Cycle",
  "MP7"  = "Columnar Progen.",
  "MP8"  = "Hypoxic Inflam.",
  "MP5"  = "Intestinal Diff."
)

sc_mp_descriptions <- c(
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

# Rename PDO list
names(pdo_list) <- paste0("PDOs_", pdo_mp_descriptions[names(pdo_list)])
names(sc_list) <- paste0("scATLAS_", sc_mp_descriptions[names(sc_list)])

### 5) Universe of genes
universe <- unique(c(unlist(pdo_list), unlist(sc_list)))

pdo_names <- names(pdo_list)
sc_names  <- names(sc_list)

### 6) Initialize matrices
jaccard_mat   <- matrix(NA_real_, length(pdo_list), length(sc_list),
                        dimnames = list(pdo_names, sc_names))
overlap_n_mat <- jaccard_mat
pval_mat      <- jaccard_mat

### 7) Compute Jaccard, overlap counts, Fisher p-values
for (i in seq_along(pdo_list)) {
  A <- pdo_list[[i]]
  for (j in seq_along(sc_list)) {
    B <- sc_list[[j]]
    
    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))
    
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j]   <- if (uni == 0) NA_real_ else inter / uni
    
    a <- inter
    b <- length(setdiff(A, B))
    c <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))
    
    pval_mat[i, j] <- if (any(c(a,b,c,d) < 0)) NA_real_
    else fisher.test(matrix(c(a, b, c, d), nrow = 2),
                     alternative = "greater")$p.value
  }
}

### 8) Adjust p-values
padj_mat <- matrix(
  p.adjust(as.vector(pval_mat), method = "BH"),
  nrow = nrow(pval_mat), ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

neglog10_fdr <- -log10(padj_mat)
neglog10_fdr[is.infinite(neglog10_fdr)] <-
  max(neglog10_fdr[is.finite(neglog10_fdr)], na.rm = TRUE) + 1

### 9) Heatmap (NO clustering, ordered by trees)
library(pheatmap)
# Transpose the Jaccard and Number matrices for sc (rows) vs PDO (cols)
jaccard_mat_t <- t(jaccard_mat)
overlap_n_mat_t <- t(overlap_n_mat)

# --- build stars from adjusted p-values (padj_mat) ---
stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat),
                    dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05]  <- "*"
stars_mat[padj_mat < 0.01]  <- "**"
stars_mat[padj_mat < 0.001] <- "***"

# transpose stars to match jaccard_mat_t / overlap_n_mat_t
stars_mat_t <- t(stars_mat)

# combine: number on first line, stars just below
display_mat <- matrix(
  paste0(overlap_n_mat_t, "\n", stars_mat_t),
  nrow = nrow(overlap_n_mat_t),
  ncol = ncol(overlap_n_mat_t),
  dimnames = dimnames(overlap_n_mat_t)
)

pheatmap(
  jaccard_mat_t,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "Gene sets overlap",
  # Names are already prefixed in your names(sc_list) / names(pdo_list)
  labels_row = rownames(jaccard_mat_t), 
  labels_col = colnames(jaccard_mat_t),
  angle_col = "90",           # Rotate column names 45 degrees
  display_numbers = display_mat,
  fontsize_number = 8, 
  number_color = "black",
  fontsize_row = 12, 
  fontsize_col = 12, 
  color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
)

#######################

mod_mat <- t(readRDS("UCell_default.rds"))
#mod_mat <- mod_mat[-c(2,9,10), ]
mod_mat <- mod_mat[-c(2,9,10), ]
save  <- t(readRDS("UCell_ref.rds"))
#save <- save[-c(2,4), ]
save <- save[-c(2,4), ]

#mod_mat <- mod_mat[, colnames(save)]

mod_mat <- mod_mat[pdo_mp_tree_order, ]
save <- save[sc_mp_tree_order, ]

rownames(mod_mat) <- paste0("PDOs_", pdo_mp_descriptions[rownames(mod_mat)])
rownames(save) <- paste0("scATLAS_", sc_mp_descriptions[rownames(save)])

cor_cross <- cor(t(mod_mat), t(save), method = "spearman")

# Transpose the correlation matrix: Rows = scATLAS, Cols = PDOs
# cor_cross was [PDO x sc], so t(cor_cross) is [sc x PDO]
cor_cross_t <- t(cor_cross)

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

col_cor <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht_cross <- Heatmap(
  cor_cross_t,
  name = "Spearman\nCorrelation",
  col = col_cor,
  
  # Clustering F
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  
  # Aesthetics
  rect_gp = gpar(col = "white", lwd = 1),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  # Column Name Rotation
  column_names_rot = 45,      # Rotate 45 degrees anticlockwise
  
  # Add correlation values (indexing cor_cross_t)
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- cor_cross_t[i, j]
    grid.text(sprintf("%.2f", val), x, y, 
              gp = gpar(fontsize = 8, 
                        col = ifelse(abs(val) > 0.5, "white", "black")))
  },
  
  column_title = "Patient-Derived Organoid (PDO) MPs (n=10)",
  row_title = "Reference scATLAS MPs (n=9)",
  column_title_side = "top"
)

draw(ht_cross)


sample_ids <- sub("_[^_]+$", "", colnames(mod_mat))
samples <- unique(sample_ids)

# sample_ids <- sapply(strsplit(colnames(mod_mat), "_"), function(x) {
#   paste(x[1:2], collapse = "_")
# })
# samples <- unique(sample_ids)

n_pdo <- nrow(mod_mat)
n_sc  <- nrow(save)

# Create 3D array: [scATLAS_MPs x PDO_MPs x Samples]
cor_array <- array(NA, dim = c(n_sc, n_pdo, length(samples)), 
                   dimnames = list(rownames(save), rownames(mod_mat), samples))

# --- 2. Loop per Sample (Compute Cross-Correlations) ---
for (smp in samples) {
  idx <- which(sample_ids == smp)
  
  if(length(idx) > 5) { # Minimum cell threshold
    # Compute cross-correlation for this sample
    # Result is [PDO x sc], so we transpose to get [sc x PDO]
    cor_sample <- cor(t(mod_mat[, idx, drop=FALSE]), 
                      t(save[, idx, drop=FALSE]), 
                      method = "spearman")
    cor_array[,,smp] <- t(cor_sample)
  }
}

# --- 3. Meta-Analysis (Averaging & T-test) ---
z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

mean_rho <- matrix(0, n_sc, n_pdo, dimnames = list(rownames(save), rownames(mod_mat)))
p_vals   <- matrix(1, n_sc, n_pdo, dimnames = list(rownames(save), rownames(mod_mat)))

for (i in 1:n_sc) {
  for (j in 1:n_pdo) {
    z_scores <- z_array[i, j, ]
    z_scores <- z_scores[!is.na(z_scores)]
    
    if(length(z_scores) > 1) {
      test_res <- t.test(z_scores)
      mean_rho[i,j] <- tanh(mean(z_scores))
      p_vals[i,j]   <- test_res$p.value
    }
  }
}

# --- 4. Beautiful Cross-Correlation Heatmap ---
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(mean_rho,
        name = "Mean Rho\n(Cross-Cor)",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        rect_gp = gpar(col = "white", lwd = 1),
        row_title = "scATLAS Meta-Programs",
        column_title = "PDO Meta-Programs (Sample-Averaged)",
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          p <- p_vals[i, j]
          r_val <- mean_rho[i, j]
          
          # Star logic
          lvl <- if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else ""
          
          grid.text(sprintf("%.2f\n%s", r_val, lvl), x, y, 
                    gp = gpar(fontsize = 9, fontface = "bold"))
        },
        
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))

################################################################
##################Plot expression of combined pdos and atlas MPs ###############

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

# --- 2. Define Row Order from GeneNMF ---
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Ensure naming matches module_scores (adjust "_" if your names use "MP1" instead of "MP_1")
mp_tree_order <- paste0("MP", mp_tree_order) 
# Remove the 10th element as per your requirement
mp_tree_order <- mp_tree_order[-8:-10]
mp_tree_order <- mp_tree_order[-1:-2]

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
rownames(mod_mat) <- paste0("scATLAS_", rownames(mod_mat))



# Define your 4 known states (NA = either 0 or 1)
template_list <- rbind(
  "Intestinal Metaplasia" = c(0, NA, 1, 1, 0),
  "Stress-adaptive"       = c(0, 0, 0, 0, 1),
  "Classic Proliferative" = c(1, 0, 0, 0, 0),
  "Basal to Intest. Meta" = c(0, 1, NA, 0, 0)
)

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
tmdata_all$manual_state[max_score < 0.1] <- "Unassigned/Quienscent"

max_score <- apply(cell_mat, 1, max)
tmdata_all$manual_state[max_score < 0.1] <- "Unassigned/Quienscent"


state_order <- c("Classic Proliferative", "Basal to Intest. Meta", "Intestinal Metaplasia", "Stress-adaptive", "Unassigned/Quienscent", "Unassigned")

tmdata_all$manual_state <- factor(tmdata_all$manual_state, levels = state_order)

manual_names <- levels(tmdata_all$manual_state)
manual_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Intestinal Metaplasia" = "#984EA3",
  "Stress-adaptive"       = "#FF7F00",
  "Unassigned/Quienscent" = "grey80",
  "Unassigned"            = "grey50"
)

add <- readRDS("UCell_pdos.rds")
add <- scale(as.matrix(add))
add <- t(add)
add <- add[-c(1,2,3,4,9,10), ]

MP_pdo <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/MP_outs_default.rds")
tree_order <- MP_pdo$programs.tree$order
ordered_clusters <- MP_pdo$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Ensure naming matches module_scores (adjust "_" if your names use "MP1" instead of "MP_1")
mp_tree_order <- paste0("MP", mp_tree_order) 
# Remove the 10th element as per your requirement
mp_tree_order <- rev(mp_tree_order)[c(-5,-8,-9)]
mp_tree_order <- mp_tree_order[c(-1,-3,-4)]

pdo_mp_descriptions <- c(
  "MP3"  = "G2M Cell Cycle",
  "MP6"  = "MYC Proliferation",
  "MP1"  = "G2M Checkpoint",
  "MP4"  = "G1S Cell Cycle",
  "MP7"  = "Columnar Progen.",
  "MP8"  = "Hypoxic Inflam.",
  "MP5"  = "Intestinal Diff."
)

add <- add[mp_tree_order, ]
rownames(add) <- paste0("PDOs_", pdo_mp_descriptions[rownames(add)])

add <- add[, colnames(mod_mat)]
mod_mat <- rbind(mod_mat, add)


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

pdf("MP_heatmap_states_subset_z_pdos.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(
  ht_manual,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

##################################################

library(Seurat)
library(ggplot2)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
tmdata_pdos <- readRDS("PDOs_final.rds")
module_scores <- readRDS("UCell_default.rds")
module_scores <- scale(as.matrix(module_scores))
geneNMF.metaprograms <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/MP_outs_default.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes

mp.genes$MP2 <- NULL
mp.genes$MP9 <- NULL
mp.genes$MP10 <- NULL
module_scores <- module_scores[, -c(2,9,10)]

mp.genes$MP1 <- NULL
mp.genes$MP3 <- NULL
mp.genes$MP4 <- NULL
module_scores <- module_scores[, c("MP5", "MP6", "MP7", "MP8"), drop = FALSE]

tmdata_pdos@meta.data[, grepl("^MP", colnames(tmdata_pdos@meta.data))] <- NULL
tmdata_pdos@meta.data <- cbind(tmdata_pdos@meta.data, module_scores)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Ensure naming matches module_scores (adjust "_" if your names use "MP1" instead of "MP_1")
mp_tree_order <- paste0("MP", mp_tree_order) 
# Remove the 10th element as per your requirement
mp_tree_order <- rev(mp_tree_order)[c(-5,-8,-9)]
mp_tree_order <- mp_tree_order[c(-1,-3,-4)]

mod_mat <- t(as.matrix(tmdata_pdos@meta.data[, mp_tree_order, drop = FALSE]))
# Create a mapping vector based on your table
mp_descriptions <- c(
  "MP3"  = "MP3_G2M_mitotic",
  "MP6"  = "MP6_MYC Biosynth",
  "MP1"  = "MP1_G2M_checkpoint",
  "MP4"  = "MP4_G1S Cycle",
  "MP7"  = "MP7_Columnar progenitor",
  "MP8"  = "MP8_Stress-induced plasticity",
  "MP5"  = "MP5_Intest diff"
)
rownames(mod_mat) <- mp_descriptions[rownames(mod_mat)]
rownames(mod_mat) <- paste0("PDOs_", rownames(mod_mat))


library(proxy)

# Define your 4 known states (NA = either 0 or 1)
template_list <- rbind(
  "Stress Plastic"      = c(0, 0, 1, 0),
  "Stress Columnar"     = c(0, 1, 1, 0),
  "Stress Proliferative"= c(1, 0, 1, 0),
  "Prolif. Columnar"    = c(1, 1, 1, 0),
  "Columnar"            = c(0, 1, 0, 0),
  "Intestinal Diff."    = c(0, 0, 0, 1)
)

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
    state_name <- key#"Unassigned"
  }
  full_template_list[[key]] <- list(pattern = pattern, state = state_name)
}

# Convert to matrix for matching
all_templates <- t(sapply(full_template_list, function(x) x$pattern))
template_states <- sapply(full_template_list, function(x) x$state)

#template_states["0,0,1,0,0"] <- template_states["0,0,0,1,0"] <- "Intest_Diff_Columnar" 

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
tmdata_pdos$manual_state[max_score < 0.1] <- "Unassigned/Quienscent"

max_score <- apply(cell_mat, 1, max)
tmdata_pdos$manual_state[max_score < 0.1] <- "Unassigned/Quienscent"


state_order <- c("Prolif. Columnar", "Classic Proliferative", "Stress Columnar", "Columnar", "Stress Plastic","Intestinal Diff.", "Unassigned/Quienscent", "Unassigned")

tmdata_pdos$manual_state <- factor(tmdata_pdos$manual_state, levels = state_order)

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

add <- readRDS("UCell_ref.rds")
add <- scale(as.matrix(add))
add <- t(add)
add <- add[-c(1,2,4,7), ]

MP_sc  <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/MP_outs_default.rds")
tree_order <- MP_sc$programs.tree$order
ordered_clusters <- MP_sc$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Ensure naming matches module_scores (adjust "_" if your names use "MP1" instead of "MP_1")
mp_tree_order <- paste0("MP", mp_tree_order) 
# Remove the 10th element as per your requirement
mp_tree_order <- mp_tree_order[-8:-10]
mp_tree_order <- mp_tree_order[-1:-2]

sc_mp_descriptions <- c(
  "MP1" = "MP1_G2M Cycle",
  "MP7" = "MP7_G1S Cycle",
  "MP5" = "MP5_MYC Biosynth",
  "MP9" = "MP9_Squamous transition",
  "MP3" = "MP3_IFN_act columnar",
  "MP6" = "MP6_Intest diff",
  "MP8" = "MP8_Stress-induced plasticity"
)

add <- add[mp_tree_order, ]
rownames(add) <- paste0("scATLAS_", sc_mp_descriptions[rownames(add)])

add <- add[, colnames(mod_mat)]
mod_mat <- rbind(mod_mat, add)


library(circlize)

max_val <- quantile(mod_mat, 0.98, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_val), c("white", "red"))

patient_levels <- unique(tmdata_pdos$SUR)
patient_cols <- setNames(
  DiscretePalette(length(patient_levels), palette = "alphabet"), 
  patient_levels
)

tmdata_pdos$batch <- ifelse(tmdata_pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia")
batch_cols  <- setNames(DiscretePalette(length(unique(tmdata_pdos$batch)), palette = "polychrome"), unique(tmdata_pdos$batch))

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

pdf("MP_heatmap_states_subset_z_ref.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(
  ht_manual,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
