####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/centred/01_centred_geneNMF.R
#   Description:
#     Replicates geneNMF.R but uses center=TRUE in multiNMF. This natively
#     transforms the log normalised matrix by making it centered per gene 
#     (subtracting gene mean) and sets negative values to zero before running 
#     NMF factorization. Extracts metaprograms for k=8:30.
####################

library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

outdir <- "Metaprogrammes_Results/centred"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tmdata_annotated <- list()

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

for (sample in sample_dirs) {
  if (!file.exists(file.path("by_samples", sample, paste0(sample, "_epi_f.rds")))) {
    next
  }
  tmdata <- readRDS(file.path("by_samples", sample, paste0(sample, "_epi_f.rds")))
  cancer <- sum(tmdata$malignancy == "malignant_level_1" | tmdata$malignancy == "malignant_level_2")
  if (cancer < 10) {
    print(paste0("Not enough cells for NMF - ", sample, ": ", sample))
    next
  }
  tmdata <- subset(tmdata, malignancy == "malignant_level_1" | malignancy == "malignant_level_2")
  tmdata_annotated[[sample]] <- tmdata
}

# Run multiNMF with center=TRUE. 
# Internally, GeneNMF:::getDataMatrix applies: 
#   mat <- t(scale(Matrix::t(mat), center=TRUE, scale=FALSE))
#   mat[mat < 0] <- 0
rds_outs <- file.path(outdir, "geneNMF_outs.rds")
if (file.exists(rds_outs)) {
  message("Loading existing geneNMF.programs...")
  geneNMF.programs <- readRDS(rds_outs)
} else {
  geneNMF.programs <- multiNMF(tmdata_annotated, assay="RNA", k=4:9, min.exp = 0.05, center = TRUE)
  saveRDS(geneNMF.programs, file=rds_outs)
}

# Extract Metaprograms for k = 8 to 30
for (k in 8:30) {
  rds_path <- file.path(outdir, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))
  png_path <- file.path(outdir, paste0("metaprograms_heatmap_nMP_", k, ".png"))
  
  if (file.exists(rds_path)) {
    message(paste("Loading existing Metaprograms for nMP =", k))
    geneNMF.metaprograms <- readRDS(rds_path)
  } else {
    message(paste("Extracting Metaprograms for nMP =", k))
    geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                            metric = "cosine",
                                            specificity.weight = 5,
                                            weight.explained = 0.5,
                                            nMP=k, 
                                            min.confidence = 0.5)
    saveRDS(geneNMF.metaprograms, file=rds_path)
  }
  
  if (!file.exists(png_path)) {
    message(paste("Plotting heatmap for nMP =", k))
    png(png_path, width = 3000, height = 2500, res = 300)
    plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0,1))
    dev.off()
  }
}
