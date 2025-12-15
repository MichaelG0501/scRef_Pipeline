library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(Seurat)
library(ggplot2)
library(purrr)
library(tidyr)
library(readxl)
library(stringr)
library(infercna)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
print(sample)

tmdata <- readRDS(paste0("by_samples/", sample, "/", sample, "_anno.rds"))
if (sum(tmdata$coexpression_loose == "singlet" & tmdata$marker_expression == "good") > 0) {
  tmdata <- subset(tmdata, coexpression_loose == "singlet" & marker_expression == "good")
} else {
  saveRDS("No cells", paste0("by_samples/", sample, "/no_cell"))
  stop("No cells")
}

if (sum(tmdata$celltype_update == "epithelial") > 0) {
  data <- subset(
    tmdata,
    celltype_update %in% c("epithelial")
  )
} else {
  saveRDS("No epithelial cells", paste0("by_samples/", sample, "/no_epi"))
  stop("No epithelial cells")
}

reference <- readRDS(paste0(sub("^([^_]+_[^_]+).*", "\\1", tmdata@meta.data$orig.ident[1]), "_reference.rds"))
meta <- data@meta.data[, c("orig.ident", "celltype_update")]
meta <- rbind(meta, reference$meta)
matrix <- as.matrix(data@assays$RNA$CPM)
matrix <- cbind(matrix, as.matrix(reference$matrix[rownames(matrix), ]))
ref <- reference$ref

outs <- infercna(matrix, refCells = ref)
saveRDS(outs, paste0("by_samples/", sample, "/", sample, "_outs.rds"))

sd_k_cor <- 2
sd_k_sig <- 2

classify_one_ident <- function(df_ident) {
  ref <- df_ident %>% filter(celltype_update != "epithelial")
  epi <- df_ident %>% filter(celltype_update == "epithelial")
  
  # start with all labeled as reference
  df_ident$classification <- ifelse(df_ident$celltype_update != "epithelial",
                                    "non_epithelial_ref", NA)
  
  if (nrow(epi) == 0 | nrow(ref) == 0) return(df_ident)
  
  mean_cor <- mean(ref$cna_cor, na.rm = TRUE)
  sd_cor   <- sd(ref$cna_cor, na.rm = TRUE)
  mean_sig <- mean(ref$cna_signal, na.rm = TRUE)
  sd_sig   <- sd(ref$cna_signal, na.rm = TRUE)
  
  thr_cor <- mean_cor + sd_k_cor * sd_cor
  thr_sig <- mean_sig + sd_k_sig * sd_sig
  
  df_ident$classification[df_ident$celltype_update == "epithelial"] <-
    ifelse(epi$cna_cor > thr_cor & epi$cna_signal > thr_sig, "cna_malignant",
           ifelse(epi$cna_cor > thr_cor | epi$cna_signal > thr_sig, "cna_unresolved", "cna_non_malignant")
    )
  
  df_ident
}

coord <- cnaScatterPlot(outs)
meta <- meta %>%
  as.data.frame() %>%
  dplyr::mutate(
    cna_cor    = coord$cna.cor,
    cna_signal = coord$cna.signal
  ) %>%
  classify_one_ident()

tmdata@meta.data <- bind_cols(
  tmdata@meta.data,
  meta[rownames(tmdata@meta.data), c("cna_cor", "cna_signal", "classification"), drop = FALSE]
)

######################################################

epi <- subset(tmdata, celltype_update == "epithelial")
ncells <- ncol(epi)

for (do in 1) {
  if (ncells > 30) {
    
    epi <- FindVariableFeatures(epi, nfeatures = 3000)
    epi <- ScaleData(epi, verbose = FALSE)
    epi <- RunPCA(epi, npcs = min(50, ncells - 10), verbose = FALSE)
    k_use  <- min(30, max(10, round(sqrt(ncells) - 2)))
    pc_use <- 1:min(30, ncells - 10)
    epi <- FindNeighbors(epi, dims = pc_use, k.param = k_use, verbose = FALSE)
    epi <- FindClusters(epi, resolution = 0.3, algorithm = 1, verbose = FALSE)
    epi <- RunUMAP(epi, dims = pc_use, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
    
    tab <- table(epi$classification, epi$seurat_clusters)
    clusters <- colnames(tab)
    cluster_labels <- setNames(rep(NA, length(clusters)), clusters)
    for (cl in clusters) {
      malignant     <- if ("cna_malignant"     %in% rownames(tab)) tab["cna_malignant", cl]     else 0
      non_malignant <- if ("cna_non_malignant" %in% rownames(tab)) tab["cna_non_malignant", cl] else 0
      unresolved    <- if ("cna_unresolved"    %in% rownames(tab)) tab["cna_unresolved", cl]    else 0
      total         <- sum(tab[, cl])
      
      pct_non_malignant <- non_malignant / total
      pct_malignant <- malignant / total
      
      if (pct_non_malignant > 0.50) {
        cluster_labels[cl] <- "non_malignant_clus"
      } else {
        if (pct_malignant > 0.50) {
          cluster_labels[cl] <- "malignant_clus"
        } else {
          cluster_labels[cl] <- "unresolved_clus"
        }
      }
    }
    epi$malignant_clus <- as.character(cluster_labels[as.character(epi$seurat_clusters)]) 
    
    malignant_ref <- colnames(epi)[epi$malignant_clus == "malignant_clus" & epi$classification == "cna_malignant"]
    
    cell_cycle_genes <- read.csv("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv", 
                                 header = TRUE, stringsAsFactors = FALSE)[, 1:3]
    cc_candidates <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
    cc_candidates <- intersect(cc_candidates, rownames(epi))
    rna <- GetAssayData(epi, assay = "RNA", slot = "data")
    
    if (length(malignant_ref) > 0) {
      avg_expr_malignant <- rowMeans(rna[cc_candidates, malignant_ref, drop = FALSE])
      cc_genes <- names(sort(avg_expr_malignant, decreasing = TRUE))[1:min(20, length(avg_expr_malignant))]
    } else {
      epi@meta.data$malignancy <- rep("unresolved", ncol(epi))
      print("Not enough malignant cells, assigning all as unresolved")
      saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi.rds"))
      break
    }
    
    cc_score <- colMeans(rna[cc_genes, , drop = FALSE])
    epi$CCscore1 <- cc_score
    ref_cc_scores <- cc_score[malignant_ref]
    cc_mean <- mean(ref_cc_scores, na.rm = TRUE)
    cc_sd   <- sd(ref_cc_scores, na.rm = TRUE)
    cc_thr <- cc_mean - 1 * cc_sd
    cc_status <- rep("cc_unresolved", ncol(epi))
    names(cc_status) <- colnames(epi)
    cc_status[malignant_ref] <- "cc_reference"
    other_cells <- setdiff(colnames(epi), malignant_ref)
    cc_status[other_cells[cc_score[other_cells] >= cc_thr]] <- "cc_malignant"
    epi$cc_status <- cc_status
    
    epi$malignancy <- rep("unresolved", ncol(epi))
    for (cl in names(cluster_labels)) {
      if (cluster_labels[cl] == "malignant_clus") {
        epi$malignancy[epi$seurat_clusters == cl] <- "malignant"
      } else if (cluster_labels[cl] == "non_malignant_clus") {
        epi$malignancy[epi$seurat_clusters == cl] <- "non_malignant"
      } else if (cluster_labels[cl] == "unresolved_clus") {
        malignant <- epi$seurat_clusters == cl & (epi$classification == "cna_malignant" | epi$cc_status == "cc_malignant")
        if (sum(malignant) / sum(epi$seurat_clusters == cl) >= 0.8) {
          epi$malignancy[epi$seurat_clusters == cl] <- "malignant"
        } else {
          epi$malignancy[epi$seurat_clusters == cl] <- "unresolved"
        }
      }
    }
    
    saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi.rds"))
    
    #######################################
    
    nonmalignant_ref <- colnames(epi)[
      epi$malignant_clus == "non_malignant_clus" &
        epi$classification == "cna_non_malignant" &
        epi$cc_status != "cc_malignant"
    ]
    cs_genes <- character(0)
    if (length(malignant_ref) > 10 && length(nonmalignant_ref) > 10) {
      epi$signature <- rep("sig_unresolved", ncol(epi))
      epi$signature[malignant_ref]    <- "sig_malignant_ref"
      epi$signature[nonmalignant_ref] <- "sig_non_malignant_ref"
      Idents(epi) <- "signature"
      dge <- FindMarkers(
        epi,
        ident.1 = "sig_malignant_ref",
        ident.2 = "sig_non_malignant_ref",
        assay = "RNA",
        logfc.threshold = 1,
        min.pct = 0.8,
        only.pos = TRUE
      )
      cs_genes <- rownames(dge[dge$p_val_adj < 0.05, ])
      saveRDS(cs_genes, paste0("by_samples/", sample, "/", sample, "_signatures.rds"))
    } else {
      print("Not enough reference cells for cancer signature scoring")
      break
    }
    
  } else {
    print("Not enough cells, assigning all as unresolved")
    epi@meta.data$malignancy <- rep("unresolved", ncol(epi))
    saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi.rds"))
  }
}
