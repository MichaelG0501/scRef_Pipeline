library(data.table)
library(dplyr)
library(circlize)
library(grid)
library(Seurat)
library(ggplot2)
library(purrr)
library(tidyr)
library(readxl)
library(stringr)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
print(sample)
epi <- readRDS(paste0("by_samples/", sample, "/", sample, "_epi.rds"))

###################################

cell_cycle_genes <- read.csv("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv", 
                             header = TRUE, stringsAsFactors = FALSE)[, 1:3]
cc_genes <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_genes <- intersect(cc_genes, rownames(epi))
gene_means <- rowMeans(epi@assays$RNA$data[cc_genes, , drop = FALSE], na.rm = TRUE)
cc_genes <- names(sort(gene_means, decreasing = TRUE))[1:10]

cc_score <- colMeans(as.matrix(epi@assays$RNA$data[cc_genes, , drop = FALSE]))
epi$cc_score <- cc_score
cc_thr <- 1
cc_status <- rep("cc_unresolved", ncol(epi))
names(cc_status) <- colnames(epi)
cc_status[cc_score >= cc_thr] <- "cc_malignant"
epi$cc_status <- cc_status

################################

cs_genes <- read.table("cancer_signatures.txt", header = FALSE)$V1
cs_genes <- intersect(cs_genes, rownames(epi))
gene_means <- rowMeans(epi@assays$RNA$data[cs_genes, , drop = FALSE], na.rm = TRUE)
cs_genes <- names(sort(gene_means, decreasing = TRUE))[1:50]
cs_score <- colMeans(as.matrix(epi@assays$RNA$data[cs_genes, , drop = FALSE]), na.rm = TRUE)
epi$cs_score <- cs_score
cs_thr <- 1
cs_status <- rep("cs_unresolved", ncol(epi))
names(cs_status) <- colnames(epi)
cs_status[cs_score >= cs_thr] <- "cs_malignant"
epi$cs_status <- cs_status

epi$malignancy <- rep("unresolved", ncol(epi))
cluster_labels <- tapply(epi$malignant_clus, epi$seurat_clusters, function(x) {
  if (all(x == "malignant_clus")) {
    "malignant_clus"
  } else if (all(x == "non_malignant_clus")) {
    "non_malignant_clus"
  } else if (all(x == "unresolved_clus")) {
    "unresolved_clus"
  }
})
for (cl in names(cluster_labels)) {
  if (cluster_labels[cl] == "non_malignant_clus") {
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_non_malignant"] <- "non_malignant_level_1"
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_unresolved"] <- "non_malignant_level_2"
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_malignant"] <- "non_malignant_unresolved"
    
  } else {
    #if ("classification" %in% colnames(epi@meta.data)) {
    malignant <- epi$seurat_clusters == cl & (epi$classification == "cna_malignant" | epi$cc_status == "cc_malignant" | epi$cs_status == "cs_malignant")
    #} else {
    #  malignant <- epi$seurat_clusters == cl & (epi$cc_status == "cc_malignant" | epi$cs_status == "cs_malignant")
    #}
    if (sum(malignant) / sum(epi$seurat_clusters == cl) >= 0.5) {
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_malignant"] <- "malignant_level_1"
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_unresolved"] <- "malignant_level_2"
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_non_malignant"] <- "malignant_unresolved"
    } else {
      epi$malignancy[epi$seurat_clusters == cl] <- "unresolved"
    }
  }
}

saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi_f.rds"))