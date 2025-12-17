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
colnames(tmdata) <- paste0(tmdata$orig.ident, "_", colnames(tmdata))
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

if (ncells > 30) {
  nfeat <- min(3000, max(2, nrow(epi)))
  epi <- FindVariableFeatures(epi, nfeatures = nfeat, verbose = FALSE)
  epi <- ScaleData(epi, features = VariableFeatures(epi), verbose = FALSE)
  max_pcs <- max(1, min(length(VariableFeatures(epi)), ncells - 1))
  npcs <- min(50, max(1, round(sqrt(ncells))))
  npcs <- min(npcs, max_pcs)
  epi <- RunPCA(epi, npcs = npcs, verbose = FALSE)
  pc_use <- 1:npcs
  if (ncells <= 10) {
    k_use <- max(1, ncells - 1)
  } else if (ncells <= 500) {
    k_use <- max(5, min(15, ncells - 1))
  } else if (ncells <= 5000) {
    k_use <- 30
  } else if (ncells <= 50000) {
    k_use <- 30
  } else {
    k_use <- min(50, round(sqrt(ncells)))   # LOWER cap than before
  }
  nn_method <- if (ncells > 20000) "annoy" else "rann"
  resolution <- if (ncells < 1000) 0.8 else if (ncells < 5000) 1.0 else 1.2
  epi <- FindNeighbors(
    epi,
    dims = pc_use,
    k.param = max(1, k_use),
    nn.method = nn_method,
    annoy.metric = "euclidean",
    verbose = FALSE
  )
  if (ncells >= 3) {
    epi <- FindClusters(epi, resolution = resolution, algorithm = 1, verbose = FALSE)
  }
  
  epi <- RunUMAP(
    epi,
    dims = pc_use,
    n.neighbors = max(2, k_use),
    min.dist = if (ncells > 5000) 0.15 else 0.3,
    spread = if (ncells > 5000) 1.5 else 1.0,
    verbose = FALSE
  )
} else {
  print("Not enough cells, assigning same X clusters")
  epi$seurat_clusters <- rep("X", ncol(epi))
  Idents(epi) <- "seurat_clusters"
}

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
nonmalignant_ref <- colnames(epi)[epi$malignant_clus == "non_malignant_clus" & epi$classification == "cna_non_malignant" ]

if (length(malignant_ref) > 50 && length(nonmalignant_ref) > 50) {
  cs_genes <- character(0)
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
    min.pct = 0.5,
    only.pos = TRUE
  )
  cs_genes <- rownames(dge[dge$p_val_adj < 0.01 & dge$pct.2 < 0.1, ])
  if (length(cs_genes) > 0) {
    saveRDS(cs_genes, paste0("by_samples/", sample, "/", sample, "_signatures.rds"))
  } else {
    print("Not enough significant cancer signatures")
  }
  saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi.rds"))
} else {
  saveRDS(epi, paste0("by_samples/", sample, "/", sample, "_epi.rds"))
  print("Not enough reference cells for cancer signature scoring")
}