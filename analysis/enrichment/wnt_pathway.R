####################
# Moved from: analysis/wnt_enrich.R
# Reorganized as part of analysis/ restructuring
# Note: WNT-specific pathway analysis with custom WNT gene sets
####################
library(pheatmap)
library(dplyr)


####################################

wnt_sets <- list(
  WNT_CM = c(
    "LGR5","ASCL2","OLFM4","PROM1","SMOC2","AXIN2","SOX9","TROY",
    "CD44","MSI1","MYC","CCND1","CCND2","CCND3","MKI67","FOXM1",
    "PCNA","RSPO2","RSPO3","DKK1","SFRP1","SFRP2","GREM1",
    "ITGA6","ITGB1","LRP6","FZD7","BCL9","PYGO2"), 
  WNT_Canonical = c("LEF1","TCF7","TCF7L2","NKD1","SNAI1","SNAI2","VIM","MMP7",
                    "MMP9","LGR5","ASCL2","PROM1","SOX9","FOSL1","BIRC5",
                    "REG1A","REG3A")
)

TERM2GENE <- data.frame(
  term = c(
    rep("WNT conditioned medium related genes", 29),
    rep("Canonical WNT genes (intrinsic WNT related genes)", 17)
  ),
  gene = c(
    "LGR5","ASCL2","OLFM4","PROM1","SMOC2","AXIN2","SOX9","TROY",
    "CD44","MSI1","MYC","CCND1","CCND2","CCND3","MKI67","FOXM1",
    "PCNA","RSPO2","RSPO3","DKK1","SFRP1","SFRP2","GREM1",
    "ITGA6","ITGB1","LRP6","FZD7","BCL9","PYGO2", 
    ###
    "LEF1","TCF7","TCF7L2","NKD1","SNAI1","SNAI2","VIM","MMP7",
    "MMP9","LGR5","ASCL2","PROM1","SOX9","FOSL1","BIRC5",
    "REG1A","REG3A"
  )
  ,
  stringsAsFactors = FALSE
)

TERM2NAME <- data.frame(
  term = c(
    "WNT conditioned medium related genes",
    "Canonical WNT genes (intrinsic WNT related genes)"
  ),
  name = c(
    "WNT conditioned medium related genes",
    "Canonical WNT genes (intrinsic WNT related genes)"
  ),
  stringsAsFactors = FALSE
)

enricher(
  gene          = mp_gene_lists$MP2,
  TERM2GENE     = TERM2GENE,
  TERM2NAME     = TERM2NAME,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

###################################

# 1. Define WNT sets
wnt_sets <- list(
  Wnt_conditioned_medium = c(
    "LGR5","ASCL2","OLFM4","PROM1","SMOC2","AXIN2","SOX9","TROY",
    "CD44","MSI1","MYC","CCND1","CCND2","CCND3","MKI67","FOXM1",
    "PCNA","RSPO2","RSPO3","DKK1","SFRP1","SFRP2","GREM1",
    "ITGA6","ITGB1","LRP6","FZD7","BCL9","PYGO2"
  ),
  Wnt_canonical = c(
    "LEF1","TCF7","TCF7L2","NKD1","SNAI1","SNAI2","VIM","MMP7",
    "MMP9","LGR5","ASCL2","PROM1","SOX9","FOSL1","BIRC5",
    "REG1A","REG3A"
  )
)

# 2. Setup names and dimensions
mp_names <- names(mp_gene_lists)
wnt_names <- names(wnt_sets)
if(!exists("universe")) universe <- unique(c(unlist(mp_gene_lists), unlist(wnt_sets)))

# 3. Initialize Matrices (MPs as columns, WNT as rows)
jaccard_mat   <- matrix(NA_real_, length(wnt_names), length(mp_names), dimnames = list(wnt_names, mp_names))
label_mat     <- matrix("", length(wnt_names), length(mp_names), dimnames = list(wnt_names, mp_names))

# 4. Compute Stats
for (j in seq_along(wnt_sets)) {
  B <- intersect(unique(wnt_sets[[j]]), universe)
  wnt_size <- length(B)
  
  for (i in seq_along(mp_gene_lists)) {
    A <- intersect(unique(mp_gene_lists[[i]]), universe)
    
    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))
    
    # Calculate Jaccard and cap it at 0.3 for the color scale if desired, 
    # but better to do it via 'breaks' in pheatmap
    jaccard_mat[j, i]   <- if (uni == 0) 0 else inter / uni
    
    # Create the label: overlap/reference size
    label_mat[j, i] <- paste0(inter, "/", wnt_size)
  }
}

# 5. Apply Tree Ordering
if(exists("mp_tree_order")){
  col_order <- rev(mp_tree_order)
  col_order <- col_order[col_order %in% colnames(jaccard_mat)]
  jaccard_mat <- jaccard_mat[, col_order, drop=FALSE]
  label_mat <- label_mat[, col_order, drop=FALSE]
}

# Create column labels: "MP1\nn=30"
mp_sizes <- sapply(mp_gene_lists[colnames(jaccard_mat)], length)
new_col_labels <- paste0(colnames(jaccard_mat), "\nn=", mp_sizes)

# 6. Define the Custom Scale (0 to 0.3)
# Any value above 0.3 will be colored with the highest color in the palette
my_breaks <- seq(0, 0.3, length.out = 101) 
my_colors <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

# 7. Render Heatmap

pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "MP overlap with Wnt genes",
  labels_row = rownames(jaccard_mat),
  labels_col = new_col_labels,
  display_numbers = label_mat, 
  fontsize_number = 8, 
  number_color = "black",
  fontsize_row = 12, 
  fontsize_col = 10, 
  angle_col = 0, 
  color = my_colors,
  breaks = my_breaks,
  legend_breaks = seq(0, 0.3, by = 0.1) # Specific ticks for the legend
)
