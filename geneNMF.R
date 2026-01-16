library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

tmdata_annotated <- list()

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

for (sample in sample_dirs) {
  if (!file.exists(file.path("by_samples", sample, paste0(sample, "_epi_f.rds")))) {
    next
  }
  tmdata <- readRDS(file.path("by_samples", sample, paste0(sample, "_epi_f.rds")))
  cancer <- sum(tmdata$malignancy == "malignant_level_1" | tmdata$malignancy == "malignant_level_2"    )
  if (cancer < 10) {
    print(paste0("Not enough cells for NMF - ", sample, ": ", sample))
    next
  }
  tmdata <- subset(tmdata, malignancy == "malignant_level_1" | malignancy == "malignant_level_2")
  tmdata_annotated[[sample]] <- tmdata
}

geneNMF.programs <- multiNMF(tmdata_annotated, assay="RNA", k=4:9, min.exp = 0.05)
saveRDS(geneNMF.programs, file="geneNMF_outs.rds")

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        specificity.weight = 5,
                                        weight.explained = 0.5,
                                        nMP=10, 
                                        min.confidence = 0.5)
saveRDS(geneNMF.metaprograms, file="MP_outs_default.rds")

anno_colors <- brewer.pal(n=10, name="Paired")
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

png("metaprograms_heatmap.png",
    width = 3000,       # pixels
    height = 2500,      # pixels
    res = 300)          # DPI (good for heatmaps)
plotMetaPrograms(
  geneNMF.metaprograms,
  annotation_colors = anno_colors,
  similarity.cutoff = c(0,1)
)
dev.off()

tmdata_all <- readRDS("EAC_Ref_merged.rds")
cells_to_extract <- as.vector(unlist(lapply(tmdata_annotated, Cells)))
tmdata_all <- subset(tmdata_all, cells = cells_to_extract)

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(tmdata_all), category = "C5", subcategory = "GO:BP")
})
saveRDS(top_p, "GO_outs.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp.genes, ncores=4, name = "")
ucell_scores <- tmdata_all@meta.data[, grep("^MP", colnames(tmdata_all@meta.data))]
saveRDS(ucell_scores, file = "UCell_default.rds")

png("vln_origident.png",
    width = 5000,   # wide
    height = 2500,
    res = 300)
VlnPlot(
  tmdata_all,
  features = names(mp.genes),
  group.by = "study",
  pt.size = 0,
  ncol = 5
)
dev.off()

png("vln_clinical_response.png",
    width = 3500,
    height = 2500,
    res = 300)
VlnPlot(
  tmdata_all,
  features = names(mp.genes),
  group.by = "Treatment",
  pt.size = 0,
  ncol = 5
)
dev.off()
