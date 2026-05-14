####################
# Analysis registry:
#   Status: active
#   Script: analysis/enrichment/nmf_enrichment.R
#   Methodology: analysis/methodology/enrichment/enrichment_external_reference_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/nmf_plot.R
# Reorganized as part of analysis/ restructuring
####################
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)  # for dotplot
library(ggplot2)

##-----------------------------
## 0. Hallmark setup (you already had this)
##-----------------------------
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

##-----------------------------
## 1. Function to get merged genes for a cluster
##    (representative + all programs with >20 gene overlap)
##-----------------------------
get_cluster_genes <- function(rep_prog, nmf_programs, overlap_list) {
  # all programs in this cluster: representative + its overlaps
  cluster_progs <- unique(c(rep_prog, overlap_list[[rep_prog]]))
  
  # assuming nmf_programs has one column per program and entries are gene symbols
  # (if it's something else, adapt this extraction step)
  genes_cluster <- as.vector(unique(na.omit(unlist(
    nmf_programs[, cluster_progs, drop = FALSE]
  ))))
  
  list(
    rep_prog = rep_prog,
    members  = cluster_progs,
    genes    = genes_cluster
  )
}

nmf_programs_all          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs_all          <- lapply(nmf_programs_all,toupper) ## convert all genes to uppercase 
nmf_filter_ccle       <- unique(unlist(lapply(nmf_programs_all, colnames)))
nmf_programs_all          <- lapply(nmf_programs_all, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs_all          <- do.call(cbind, nmf_programs_all)

cluster_gene_sets <- lapply(names(final_nmf), get_cluster_genes,
                            nmf_programs = nmf_programs_all,
                            overlap_list = overlap_list)
names(cluster_gene_sets) <- names(final_nmf)

##-----------------------------
## 3. Run GO BP + Hallmark enrichment on each *cluster* gene set
##-----------------------------
cluster_enrich <- lapply(cluster_gene_sets, function(cl) {
  genes <- cl$genes
  
  res_GO <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  res_H <- enricher(
    gene          = genes,
    TERM2GENE     = hallmark_term2gene,
    TERM2NAME     = hallmark_term2name,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  list(
    rep_prog = cl$rep_prog,
    members  = cl$members,
    GO       = res_GO,
    Hallmark = res_H
  )
})

##-----------------------------
## 4. Plot enrichment for each representative program
##    - One GO dotplot + one Hallmark dotplot per cluster
##    - showCategory controls how many top terms you see
##-----------------------------
for (rp in names(cluster_enrich)) {
  cat("Plotting enrichment for cluster represented by:", rp, "\n")
  enr <- cluster_enrich[[rp]]
  
  ## GO BP
  if (!is.null(enr$GO) && nrow(enr$GO@result) > 0) {
    p_go <- dotplot(enr$GO, showCategory = 10) +
      ggtitle(paste0(rp, " cluster – GO BP"))
    print(p_go)
  } else {
    message("No significant GO terms for ", rp)
  }
  
  ## Hallmark
  if (!is.null(enr$Hallmark) && nrow(enr$Hallmark@result) > 0) {
    p_h <- dotplot(enr$Hallmark, showCategory = 10) +
      ggtitle(paste0(rp, " cluster – Hallmark"))
    print(p_h)
  } else {
    message("No significant Hallmark terms for ", rp)
  }
}


for (rp in names(enrich_res)) {
  cat("Plotting enrichment for cluster represented by:", rp, "\n")
  enr <- enrich_res[[rp]]
  
  ## GO BP
  if (!is.null(enr$GO) && nrow(enr$GO@result[enr$GO@result$p.adjust < 0.05, ]) > 0) {
    p_go <- barplot(enr$GO, showCategory = 10) +
      ggtitle(rp) +
      theme(
        plot.title = element_text(size = 16),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17, face = "bold"),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14)
      )
    print(p_go)
  } else {
    message("No significant GO terms for ", rp)
  }
  
  ## Hallmark
  if (!is.null(enr$Hallmark) && nrow(enr$Hallmark@result[enr$Hallmark@result$p.adjust < 0.05, ]) > 0) {
    p_h <- barplot(enr$Hallmark, showCategory = 10) +
      ggtitle(rp) +
      theme(
        plot.title = element_text(size = 16),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17, face = "bold"),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14)
      )
    print(p_h)
  } else {
    message("No significant Hallmark terms for ", rp)
  }
}



library(patchwork)

for (rp in names(enrich_res)) {
  cat("Plotting enrichment for cluster represented by:", rp, "\n")
  enr <- enrich_res[[rp]]
  
  ## GO BP
  p_go <- NULL
  if (!is.null(enr$GO) && nrow(enr$GO@result[enr$GO@result$p.adjust < 0.05, ]) > 0) {
    p_go <- barplot(enr$GO, showCategory = 6) +
      theme(
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17, face = "bold"),
        axis.title  = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title= element_text(size = 14)
      )
  } else {
    message("No significant GO terms for ", rp)
  }
  
  ## Hallmark
  p_h <- NULL
  if (!is.null(enr$Hallmark) && nrow(enr$Hallmark@result[enr$Hallmark@result$p.adjust < 0.05, ]) > 0) {
    p_h <- barplot(enr$Hallmark, showCategory = 6) +
      theme(
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17, face = "bold"),
        axis.title  = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title= element_text(size = 14)
      )
  } else {
    message("No significant Hallmark terms for ", rp)
  }
  
  ## Combine into one plot with a single title
  if (!is.null(p_go) && !is.null(p_h)) {
    combined <- (p_go / p_h) + plot_annotation(title = rp)
    print(combined)
  } else if (!is.null(p_go)) {
    print(p_go + ggtitle(rp))
  } else if (!is.null(p_h)) {
    print(p_h + ggtitle(rp))
  }
}








library(pheatmap)

# 1) Cluster columns to obtain the order (ignore previous dendrograms)
cols <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

p_tmp <- pheatmap(
  J,
  color = cols,
  cluster_rows = FALSE,
  cluster_cols = TRUE,    # let pheatmap compute clustering
  display_numbers = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  silent = TRUE
)

# 2) Extract clustered order and original program IDs IN THAT ORDER
col_order <- p_tmp$tree_col$order
orig_prog_order <- colnames(J)[col_order]   # store originals before renaming

# 3) Reorder matrices to the clustered order
J_ord  <- J[,  col_order, drop = FALSE]
OV_ord <- OV[, col_order, drop = FALSE]

# 4) Create sequential names in that order and apply them
new_names <- paste0("NMF_", seq_along(col_order))
colnames(J_ord)  <- new_names
colnames(OV_ord) <- new_names

# 5) Annotation (BestLabel aligned to the clustered order, keyed by new names)
annotation_col <- data.frame(
  BestLabel = best_label[col_order],
  row.names = new_names,
  check.names = FALSE
)

# 6) Mapping table (correct: NMF_* ↔ original program ↔ best label)
matching_tbl <- data.frame(
  NMF       = new_names,
  Program   = orig_prog_order,                 # original program IDs
  BestLabel = best_label[col_order],
  Overlap   = best_overlap[col_order],
  Jaccard   = round(best_jaccard[col_order], 3),
  check.names = FALSE
)
print(matching_tbl)

# 7) Final heatmap: show only NMF_* labels, keep the clustered order
pheatmap(
  J_ord,
  color = cols,
  cluster_rows = FALSE,
  cluster_cols = FALSE,     # do NOT recluster; keep the order you set
  display_numbers = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 9,
  main = "Jaccard similarity: NMF programs vs MP signatures",
  border_color = NA
)


library(writexl)

# Page 1 = nmf_programs subset
page1 <- nmf_programs[, names(final_nmf), drop = FALSE]

# Build mapping: original IDs -> NMF_x
nmf_map <- setNames(matching_tbl[,1], rownames(matching_tbl))

# Rename columns
colnames(page1) <- nmf_map[colnames(page1)]

# Explicit order NMF_1 ... NMF_12
nmf_order <- paste0("NMF_", 1:12)
page1 <- page1[, nmf_order, drop = FALSE]

# Page 2 = cluster_gene_sets
page2 <- lapply(cluster_gene_sets, function(x) x$genes)
max_len <- max(sapply(page2, length))
page2_mat <- sapply(page2, function(g) c(g, rep(NA, max_len - length(g))))

# Rename and reorder columns
colnames(page2_mat) <- nmf_map[colnames(page2_mat)]
page2_df <- as.data.frame(page2_mat, stringsAsFactors = FALSE)
page2_df <- page2_df[, nmf_order, drop = FALSE]

# Write both sheets in correct order
write_xlsx(
  list(
    Page1 = as.data.frame(page1, stringsAsFactors = FALSE),
    Page2 = page2_df
  ),
  path = "nmf_gene_programs.xlsx"
)
