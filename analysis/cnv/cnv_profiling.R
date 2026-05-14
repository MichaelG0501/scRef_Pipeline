####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/cnv_profiling.R
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/cnv_profile_sc.R
# Reorganized as part of analysis/ restructuring
####################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(Seurat)
})


merged_obj <- readRDS("EAC_Ref_merged.rds")

library(infercna)
library(dplyr)
library(ggplot2)

# data <- subset(
#   merged_obj,
#   (orig.ident == "sampleC" & celltype_group %in% c("epithelial")) |
#     celltype_group %in% c("endothelial", "macrophage")
# )
# 
# data <- subset(
#   merged_obj,
#   celltype_manual %in% c("epithelial")
# )
# 
# data <- FindVariableFeatures(data)
# data <- ScaleData(data)
# data <- RunPCA(data)
# data <- FindNeighbors(data, dims = 1:20)
# 
# data <- FindClusters(data, resolution = 0.8, algorithm = 1)  # Leiden
# data$leiden_clusters <- Idents(data)
# 
# data <- RunUMAP(data, dims = 1:20)
# saveRDS(data, "all_epithelial.rds")
# 
p1 <- DimPlot(data, group.by = "leiden_clusters", label = TRUE) + ggtitle("Louvain Clustering")
p3 <- DimPlot(data, group.by = "study", label = FALSE) + ggtitle("Study Origin")
p3 <- DimPlot(data, group.by = "orig.ident", label = FALSE) + ggtitle("Identity") + NoLegend()

combined_plot <- p1 + p3

ggsave("clustering_comparison_epithelial.png", plot = combined_plot, width = 14, height = 6, dpi = 300)

data <- subset(
  merged_obj,
  celltype_update %in% c("epithelial") |
    celltype_update %in% c("macrophage", "b.cell")
)

meta <- data@meta.data[, c("orig.ident", "celltype_update")]

matrix <- as.matrix(data@assays$RNA$CPM)

ref_ct <- c("b.cell", "macrophage")
ct_col <- "celltype_manual"

ref <- lapply(ref_ct, function(ct) {
  colnames(data)[data@meta.data[[ct_col]] == ct]  # Extract cell IDs for each type
})
names(ref) <- ref_ct

outs <- infercna(matrix, refCells = ref)
#outs <- readRDS("cynthia_H_T1_outs.rds")

saveRDS(outs, "all_outs.rds")
