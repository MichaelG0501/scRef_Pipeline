####################
# Moved from: analysis/robust_NMF.R
# Reorganized as part of analysis/ restructuring
####################
library(Seurat)
library(NMF)
library(reshape2)
library(ggplot2)
library(scales)
library(viridis)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

Genes_nmf_w_basis <- list()

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("_PDO$", sample_dirs)]

for (sample in sample_dirs[-21]) {
  Genes_nmf_w_basis[[sample]] <- readRDS(file.path("by_samples", sample, paste0(sample, "_rank4_9_nrun10.RDS")))
}

saveRDS(Genes_nmf_w_basis, "Genes_nmf_w_basis_all.rds")

robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  a <- inter_intersect
  b <- sort(apply(a, 2, function(x) sort(x, decreasing = TRUE)[2]), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
  if(length(b) > 1) {
    c <- names(b[1]) 
    for(y in 2:length(b)) {
      if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
    }
    final_filter <- c(final_filter, c)
  } else {
    final_filter <- c(final_filter, names(b))
  }
  return(final_filter)                                                      
}

intra_min_parameter <- 30
intra_max_parameter <- 20

names(Genes_nmf_w_basis) <- gsub("nruns", "nrun", names(Genes_nmf_w_basis))

nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs          <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
inter_intersect         <- apply(do.call(cbind, nmf_programs) , 2, function(x) apply(do.call(cbind, nmf_programs) , 2, function(y) length(intersect(x,y)))) 

# not filtering inter sample
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 

nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

library(pheatmap)

highlight <- colnames(inter_intersect) %in% nmf_filter_ccle
annotation_col <- data.frame(
  Highlight = ifelse(highlight, "Yes", "No"),
  Sample = sub("_rank.*", "", colnames(inter_intersect))
)

annotation_row <- data.frame(
  Highlight = ifelse(rownames(inter_intersect) %in% nmf_filter_ccle, "Yes", "No"),
  Sample = sub("_rank.*", "", rownames(inter_intersect))
)

# Set row names for annotations
rownames(annotation_col) <- colnames(inter_intersect)
rownames(annotation_row) <- rownames(inter_intersect)

# Plot with annotations
pheatmap(
  inter_intersect,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(Highlight = c(Yes = "red", No = "white")),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Highlighted NMF Programs"
)

heatmap(
  nmf_intersect,
  Rowv = nmf_intersect_hc,
  Colv = nmf_intersect_hc,
  scale = "none",
  col = colorRampPalette(c("white", "blue"))(100),
  margins = c(10, 10),
  main = "NMF Program Similarity"
)




overlap_list <- lapply(nmf_filter_ccle, function(prog) {
  overlaps <- inter_intersect[prog, ]
  others <- names(overlaps[overlaps > 20])
  return(others)
})

names(overlap_list) <- nmf_filter_ccle


get_prefix <- function(x) sub("_rank.*", "", x)

unique_sample_counts <- sapply(overlap_list, function(progs) {
  length(unique(get_prefix(progs)))
})
final_nmf <- unique_sample_counts[unique_sample_counts >= 3]

library(clusterProfiler)
library(org.Hs.eg.db)

enrich_res <- lapply(as.data.frame(nmf_programs[, names(final_nmf)]), function(genes) {
  enrichGO(gene          = genes,
           OrgDb         = org.Hs.eg.db,
           keyType       = "SYMBOL",
           ont           = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff  = 0.05)
})

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

# 1. Get Hallmark gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

# 2. Build TERM2GENE and TERM2NAME objects
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]  # pathway names

# 3. Run enrichment
enrich_res <- lapply(as.data.frame(nmf_programs[, names(final_nmf)]), function(genes) {
  list(
    GO = enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    ),
    Hallmark = enricher(
      gene          = genes,
      TERM2GENE     = hallmark_term2gene,
      TERM2NAME     = hallmark_term2name,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    )
  )
})


saveRDS(enrich_res, "enrichGOtempshort.rds")

enrich_res_sub <- enrich_res[names(final_nmf)]

top_terms <- lapply(enrich_res_sub, function(df) {
  df <- df[order(df$p.adjust), ]
  head(df$Description, 5)
})

all_terms <- unique(unlist(top_terms))
mat <- sapply(top_terms, function(terms) all_terms %in% terms)
mat <- t(mat) * 1  # binary matrix
rownames(mat) <- names(top_terms)
colnames(mat) <- all_terms

dist_mat <- dist(mat, method = "binary")  # Jaccard-like
hc <- hclust(dist_mat, method = "average")


#################

selected <- unique(unlist(overlap_list[names(final_nmf)], use.names = FALSE))
inter_intersect <- inter_intersect[selected, selected]

representative <- colnames(inter_intersect) %in% names(final_nmf)
annotation_col <- data.frame(
  Representative = ifelse(representative, "Yes", "No"),
  Sample    = sub("_rank.*", "", colnames(inter_intersect)))
annotation_col$Batch <- ifelse(
  grepl("_PDO", annotation_col$Sample) & 
    !grepl("Untreated_PDO|Treated_PDO", annotation_col$Sample),
  "batch_cynthia",
  ifelse(grepl("Untreated_PDO|Treated_PDO", annotation_col$Sample), "batch2", NA))
annotation_row <- data.frame(
  Representative = ifelse(rownames(inter_intersect) %in% nmf_filter_ccle, "Yes", "No"),
  Sample    = sub("_rank.*", "", rownames(inter_intersect)))
annotation_row$Batch <- ifelse(
  grepl("_PDO", annotation_row$Sample) & 
    !grepl("Untreated_PDO|Treated_PDO", annotation_row$Sample),
  "batch_cynthia",
  ifelse(grepl("Untreated_PDO|Treated_PDO", annotation_row$Sample), "batch2", NA))

# Set row names for annotations
rownames(annotation_col) <- colnames(inter_intersect)
rownames(annotation_row) <- rownames(inter_intersect)

sample_levels <- unique(c(annotation_col$Sample, annotation_row$Sample))
sample_palette <- setNames(
  viridis(length(sample_levels), option = "D"),
  sample_levels
)
pheatmap(
  inter_intersect,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(
    Representative = c(Yes = "red", No = "white"),
    Batch = c(batch_cynthia = "brown", batch2 = "darkgreen"), Sample = sample_palette
  ),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Highlighted NMF Programs"
)


nmf_intersect <- nmf_intersect[names(final_nmf), names(final_nmf)]
nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# Compute distance and clustering separately for rows and columns
row_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
col_hc <- hclust(as.dist(t(50 - nmf_intersect)), method = "average")
row_dend <- reorder(as.dendrogram(row_hc), rowMeans(nmf_intersect))
col_dend <- reorder(as.dendrogram(col_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(row_dend), order.dendrogram(col_dend)]
heatmap(
  nmf_intersect,
  Rowv = row_dend,
  Colv = col_dend,
  scale = "none",
  col = colorRampPalette(c("white", "blue"))(100),
  margins = c(10, 10),
  main = "NMF Program Similarity"
)



nmf_programs <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
#nmf_programs <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
nmf_programs <- do.call(cbind, nmf_programs)
all_genes <- nmf_programs[, selected]
all_genes <- as.vector(all_genes)
all_genes <- unique(all_genes)






MP_list <- readRDS("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/MP_list.RDS")
MP_list <- MP_list$Cancer

prog_names <- names(final_nmf)
nmf_sub    <- nmf_programs[, prog_names, drop = FALSE]

prog_list <- setNames(lapply(prog_names, function(pn) {
  unique(toupper(na.omit(nmf_sub[, pn])))
}), prog_names)

mp_list_u <- lapply(MP_list, function(x) unique(toupper(x)))
mp_names  <- names(mp_list_u)
mp_labels <- sub("^MP\\d+\\s+", "", mp_names)       # human-readable labels
names(mp_labels) <- mp_names

# Gene universe (for optional enrichment)
U <- unique(c(unlist(prog_list), unlist(mp_list_u)))

# ---------------------------
# 2) Similarity (Jaccard) and overlaps
# ---------------------------
jaccard <- function(A, B) {
  inter <- length(intersect(A, B)); union <- length(union(A, B))
  if (union == 0) 0 else inter / union
}
overlap_n <- function(A, B) length(intersect(A, B))

J  <- matrix(0, nrow = length(mp_names), ncol = length(prog_names),
             dimnames = list(mp_labels, prog_names))
OV <- matrix(0, nrow = length(mp_names), ncol = length(prog_names),
             dimnames = list(mp_labels, prog_names))

for (i in seq_along(mp_names)) {
  A <- mp_list_u[[mp_names[i]]]
  for (j in seq_along(prog_names)) {
    B <- prog_list[[prog_names[j]]]
    J[i, j]  <- jaccard(A, B)
    OV[i, j] <- overlap_n(A, B)
  }
}

# Best MP label per program
best_idx     <- apply(J, 2, which.max)
best_label   <- rownames(J)[best_idx]             # plain MP label
best_overlap <- apply(OV, 2, max)
best_jaccard <- apply(J, 2, max)

annotation_col <- data.frame(
  BestLabel = best_label,
  Overlap   = best_overlap,
  Jaccard   = round(best_jaccard, 3),
  row.names = prog_names,
  check.names = FALSE
)

# ---------------------------
# 3) Heatmap: MP vs NMF (Jaccard), with overlaps as numbers
# ---------------------------
cols <- colorRampPalette(c("#ffffff", "#cce5ff", "#80bfff", "#3385ff", "#0047b3"))(100)

pheatmap(
  J,
  color = cols,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = OV,           # show |A∩B|
  number_format = "%d",
  number_color = "black",
  fontsize_row = 10,
  fontsize_col = 9,
  main = "Jaccard similarity: NMF programs vs MP signatures",
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  border_color = NA
)

# ---------------------------
# 4) Relabel nmf_intersect by best MP labels and plot
# ---------------------------
# Map program → best label (plain). If duplicates, make them unique for display.
prog_to_label <- setNames(best_label, prog_names)

# Keep the original clustering order you prefer, then relabel names
# Recompute dendrograms concisely (average linkage on a distance derived from nmf_intersect)
# If nmf_intersect already aligned to prog_names both ways:
nmf_intersect <- nmf_intersect[prog_names, prog_names]

# Distance: larger shared → smaller distance; here using (max - value)
Dval <- max(nmf_intersect, na.rm = TRUE)
row_hc <- hclust(as.dist(Dval - nmf_intersect), method = "average")
col_hc <- hclust(as.dist(Dval - t(nmf_intersect)), method = "average")
row_dend <- reorder(as.dendrogram(row_hc), rowMeans(nmf_intersect))
col_dend <- reorder(as.dendrogram(col_hc), colMeans(nmf_intersect))

nmf_intersect <- nmf_intersect[order.dendrogram(row_dend), order.dendrogram(col_dend)]

# Replace program IDs with best labels (ensure uniqueness)
lab_rows <- make.unique(prog_to_label[rownames(nmf_intersect)], sep = " · ")
lab_cols <- make.unique(prog_to_label[colnames(nmf_intersect)], sep = " · ")
rownames(nmf_intersect) <- lab_rows
colnames(nmf_intersect) <- lab_cols

# Plot (pheatmap for cleaner labeling)
pheatmap(
  nmf_intersect,
  color = colorRampPalette(c("white", "blue"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,  # already ordered by dendrograms above
  fontsize_row = 10, fontsize_col = 9,
  main = "NMF Program Similarity (best MP labels)",
  border_color = NA
)

