library(Seurat)
library(NMF)
library(reshape2)
library(ggplot2)
library(scales)
library(doParallel)
library(foreach)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
print(sample)
epi <- readRDS(paste0("by_samples/", sample, "/", sample, "_epi_f.rds"))

epi <- subset(epi, malignancy == "malignant_good" | malignancy == "malignant_ok")
if (ncol(epi) < 10) {
  saveRDS("Not enough malignant cells for NMF", paste0("by_samples/", sample, "/no_cancer"))
  stop("Not enough cells for NMF")
}

nmf.options(parallel = 6)

epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
var_genes <- VariableFeatures(epi)

expr_matrix <- GetAssayData(epi, layer = "CPM")
expr_var <- expr_matrix[var_genes, ] / 10
expr_norm <- t(scale(t(expr_var), center = TRUE, scale = FALSE))
expr_norm <- pmax(expr_norm, 0)
expr_norm <- expr_norm[rowSums(is.na(expr_norm)) == 0 & rowSums(expr_norm) > 0, ]

rank_all <- list()
for (rank in 4:9) {
  message(paste("Running NMF for", sample, "rank", rank))
  nmf_result <- nmf(expr_norm, rank = rank, nrun = 10, method = "brunet")
  basis_matrix <- basis(nmf_result)
  colnames(basis_matrix) <- paste0(sample, "_rank4_9_nrun10.RDS.", rank, ".", seq_len(rank))
  rank_all[[rank]] <- basis_matrix
  gc()
}
result <- do.call(cbind, rank_all)

saveRDS(result, paste0("by_samples/", sample, "/", sample, "_rank4_9_nrun10.RDS"))
