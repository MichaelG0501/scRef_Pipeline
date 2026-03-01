####################
# Moved from: analysis/terms_overlap.R
# Reorganized as part of analysis/ restructuring
####################
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)

#-----------------------------
# 1) Build gene-set lists from TERM2GENE
#-----------------------------
t2g <- save_list$TERM2GENE %>%
  mutate(
    term = as.character(term),
    gene = as.character(gene)
  )

# Normal development: stomach subset
dev_stomach_terms <- t2g %>%
  filter(grepl("Normal_Development", term),
         grepl("Stomach", term)) %>%
  group_by(term) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

# Adult epithelium: all
adult_epi_terms <- t2g %>%
  filter(grepl("Adult_Epithelium", term)) %>%
  group_by(term) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

# Named lists (term -> vector of genes)
dev_list   <- setNames(dev_stomach_terms$genes, dev_stomach_terms$term)
adult_list <- setNames(adult_epi_terms$genes,  adult_epi_terms$term)

# Optional: shorten labels (keeps heatmap readable)
names(dev_list) <- gsub("\\.\\.", "_", names(dev_list))
names(adult_list) <- gsub("\\.\\.", "_", names(adult_list))

#-----------------------------
# 2) Universe of genes
#-----------------------------
universe <- unique(c(unlist(dev_list), unlist(adult_list)))

dev_names   <- names(dev_list)
adult_names <- names(adult_list)

#-----------------------------
# 3) Initialize matrices
#    rows = dev_stomach, cols = adult_epi
#-----------------------------
jaccard_mat <- matrix(NA_real_, length(dev_list), length(adult_list),
                      dimnames = list(dev_names, adult_names))
overlap_n_mat <- jaccard_mat
pval_mat <- jaccard_mat

#-----------------------------
# 4) Compute Jaccard, overlap counts, Fisher p-values
#-----------------------------
for (i in seq_along(dev_list)) {
  A <- dev_list[[i]]
  for (j in seq_along(adult_list)) {
    B <- adult_list[[j]]
    
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

#-----------------------------
# 5) Adjust p-values (BH/FDR) + stars
#-----------------------------
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

#-----------------------------
# 6) Heatmap
#-----------------------------
pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "Normal Development (Stomach) vs Adult Epithelium",
  angle_col = 90,
  display_numbers = display_mat,
  fontsize_number = 7,
  number_color = "black",
  fontsize_row = 10,
  fontsize_col = 9,
  color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
)
