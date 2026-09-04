####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/nmf_celltype_geneNMF.R
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Per-cell-type GeneNMF, consensus MP construction, enrichment,
#     and UCell scoring for macrophage/fibroblast/endothelial/plasma/NK/CD4/CD8 cells.
#   Inputs: ref_outs/by_samples/<sample>/<sample>_anno.rds; cell type is the required CLI argument.
#   Outputs: ref_outs/nmf_<celltype>/{geneNMF_outs.rds,MP_outs_default.rds,tmdata_all.rds,GO_outs.rds,UCell_default.rds}.
#   Cache/replot: full rebuild; outputs are persistent inputs to cross-cell-type analysis.
#   Run: qsub analysis/non_malignant_nmf/submit_geneNMF_all.sh
#   Conda env: gnmf
####################

library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
celltype_arg <- args[1]

base_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
setwd(base_dir)

# ---- Cell type mapping ----
celltype_map <- list(
  macrophage  = list(ct = "macrophage",  folder = "nmf_macrophage"),
  fibroblast  = list(ct = "fibroblast",  folder = "nmf_fibroblast"),
  endothelial = list(ct = "endothelial", folder = "nmf_endothelial"),
  nk.cell     = list(ct = "nk.cell",     folder = "nmf_nk"),
  plasma      = list(ct = "plasma",      folder = "nmf_plasma"),
  cd4         = list(ct = "t.cell",      folder = "nmf_cd4",  subtype = "cd4"),
  cd8         = list(ct = "t.cell",      folder = "nmf_cd8",  subtype = "cd8")
)

ct_info <- celltype_map[[celltype_arg]]
if (is.null(ct_info)) stop("Unknown celltype: ", celltype_arg)

out_dir <- file.path(base_dir, ct_info$folder)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- CD4/CD8 classification (vectorised) ----
classify_tcell <- function(seurat_obj) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  cd4_e  <- if ("CD4"  %in% rownames(expr)) as.numeric(expr["CD4",  ]) else rep(0, ncol(expr))
  cd8a_e <- if ("CD8A" %in% rownames(expr)) as.numeric(expr["CD8A", ]) else rep(0, ncol(expr))
  cd8b_e <- if ("CD8B" %in% rownames(expr)) as.numeric(expr["CD8B", ]) else rep(0, ncol(expr))

  subtype <- rep("cd8", length(cd4_e))

  # All three > 0
  all_pos <- cd4_e > 0 & cd8a_e > 0 & cd8b_e > 0
  subtype[all_pos] <- ifelse(cd4_e[all_pos] > (cd8a_e[all_pos] + cd8b_e[all_pos]) / 2, "cd4", "cd8")

  # CD8A == 0 (covers both-zero case too)
  c2 <- !all_pos & cd8a_e == 0
  subtype[c2] <- ifelse(cd4_e[c2] > cd8b_e[c2], "cd4", "cd8")

  # CD8B == 0 and CD8A != 0
  c3 <- !all_pos & !c2 & cd8b_e == 0
  subtype[c3] <- ifelse(cd4_e[c3] > cd8a_e[c3], "cd4", "cd8")

  # Remaining (CD4 == 0, both CD8s > 0) stays "cd8"
  return(subtype)
}

# ---- Load per-sample data ----
tmdata_annotated <- list()
sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

for (sample in sample_dirs) {
  rds_path <- file.path("by_samples", sample, paste0(sample, "_anno.rds"))
  if (!file.exists(rds_path)) next

  tmdata <- readRDS(rds_path)

  # Quality filter (same as Expr_filtering.R)
  tmdata <- tryCatch({
    tmp <- subset(tmdata,
                  subset = (marker_expression == "good" | is.na(marker_expression)) &
                           coexpression_loose == "singlet")
    if (ncol(tmp) <= 1) NULL else tmp
  }, error = function(e) NULL)
  if (is.null(tmdata)) { rm(tmdata); gc(); next }

  # Subset to target cell type
  ct_mask <- tmdata$celltype_update == ct_info$ct
  if (sum(ct_mask) < 1) { rm(tmdata); gc(); next }
  tmdata <- subset(tmdata, cells = colnames(tmdata)[ct_mask])

  # T cell CD4/CD8 split
  if (celltype_arg %in% c("cd4", "cd8")) {
    ####################
    # Skip single-cell T-cell subsets that can produce malformed
    # Assay5 layer dimensions after subsetting in Seurat v5.
    ####################
    if (ncol(tmdata) < 2) {
      print(paste0("Skipping ", sample, ": only ", ncol(tmdata), " t.cell cells before CD4/CD8 split"))
      rm(tmdata); gc(); next
    }
    ####################
    # Robust fallback: if layer extraction fails for this sample,
    # skip it instead of aborting the full celltype run.
    ####################
    subtype_try <- tryCatch(
      classify_tcell(tmdata),
      error = function(e) {
        print(paste0("Skipping ", sample, ": classify_tcell failed - ", e$message))
        NULL
      }
    )
    if (is.null(subtype_try)) {
      rm(tmdata); gc(); next
    }
    tmdata$tcell_subtype <- subtype_try
    keep <- tmdata$tcell_subtype == ct_info$subtype
    if (sum(keep) < 1) { rm(tmdata); gc(); next }
    tmdata <- subset(tmdata, cells = colnames(tmdata)[keep])
  }

  # Minimum cell check
  if (ncol(tmdata) < 10) {
    print(paste0("Skipping ", sample, ": only ", ncol(tmdata), " ", celltype_arg, " cells"))
    rm(tmdata); gc(); next
  }

  # Strip existing sample prefix if present to avoid double-prefixing by multiNMF
  tmdata <- RenameCells(tmdata, new.names = gsub(paste0("^", sample, "_"), "", colnames(tmdata)))
  tmdata_annotated[[sample]] <- tmdata
  print(paste0("Loaded ", sample, ": ", ncol(tmdata), " ", celltype_arg, " cells"))
  rm(tmdata); gc()
}

cat("\n=== Total samples with sufficient ", celltype_arg, " cells: ",
    length(tmdata_annotated), " ===\n")

if (length(tmdata_annotated) < 2) {
  stop("Not enough samples (need >= 2) with sufficient cells for multiNMF")
}

# ---- Run multiNMF ----
geneNMF.programs <- multiNMF(tmdata_annotated, assay = "RNA", k = 4:9, min.exp = 0.05, center = TRUE)
saveRDS(geneNMF.programs, file = file.path(out_dir, "geneNMF_outs.rds"))

# ---- Get MetaPrograms ----
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                         metric = "cosine",
                                         specificity.weight = 5,
                                         weight.explained = 0.5,
                                         nMP = 10,
                                         min.confidence = 0.5)
saveRDS(geneNMF.metaprograms, file = file.path(out_dir, "MP_outs_default.rds"))

# ---- Heatmap ----
n_mp <- length(geneNMF.metaprograms$metaprograms.genes)
if (n_mp > 0) {
  pal_n <- max(min(n_mp, 12), 3)
  anno_colors <- brewer.pal(n = pal_n, name = "Paired")[1:min(n_mp, 12)]
  if (n_mp > 12) anno_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_mp)
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png(file.path(out_dir, "metaprograms_heatmap.png"),
      width = 3000, height = 2500, res = 300)
  plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0, 1))
  dev.off()
}

# ---- Merge per-sample objects for downstream ----
for (nm in names(tmdata_annotated)) {
  tmdata_annotated[[nm]] <- RenameCells(tmdata_annotated[[nm]], add.cell.id = nm)
}
if (length(tmdata_annotated) > 1) {
  tmdata_all <- merge(tmdata_annotated[[1]],
                      y = tmdata_annotated[2:length(tmdata_annotated)])
} else {
  tmdata_all <- tmdata_annotated[[1]]
}
rm(tmdata_annotated); gc()

# Add study from orig.ident
tmdata_all@meta.data$study <- sapply(
  strsplit(tmdata_all@meta.data$orig.ident, "_"),
  function(x) paste(x[1:2], collapse = "_"))

# Add Treatment from Excel
excel_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Summary_EAC_Ref.xlsx"
if (file.exists(excel_path) && !("Treatment" %in% colnames(tmdata_all@meta.data))) {
  tryCatch({
    library(readxl)
    xl <- read_excel(excel_path, sheet = 2, skip = 1)
    xl <- xl %>%
      mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_")) %>%
      select(orig.ident, Treatment) %>%
      distinct()
    cell_names <- colnames(tmdata_all)
    tmdata_all@meta.data <- tmdata_all@meta.data %>% left_join(xl, by = "orig.ident")
    rownames(tmdata_all@meta.data) <- cell_names
  }, error = function(e) message("Could not load Excel metadata: ", e$message))
}
saveRDS(tmdata_all, file = file.path(out_dir, "tmdata_all.rds"))

# ---- GSEA ----
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe = rownames(tmdata_all), category = "C5", subcategory = "GO:BP")
})
saveRDS(top_p, file.path(out_dir, "GO_outs.rds"))

# ---- UCell scoring ----
mp.genes <- geneNMF.metaprograms$metaprograms.genes
tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp.genes, ncores = 2, name = "")
ucell_scores <- tmdata_all@meta.data[, grep("^MP", colnames(tmdata_all@meta.data))]
saveRDS(ucell_scores, file = file.path(out_dir, "UCell_default.rds"))

# ---- Violin plots ----
png(file.path(out_dir, "vln_origident.png"),
    width = 5000, height = 2500, res = 300)
VlnPlot(tmdata_all, features = names(mp.genes),
        group.by = "study", pt.size = 0, ncol = 5)
dev.off()

if ("Treatment" %in% colnames(tmdata_all@meta.data)) {
  png(file.path(out_dir, "vln_clinical_response.png"),
      width = 3500, height = 2500, res = 300)
  VlnPlot(tmdata_all, features = names(mp.genes),
          group.by = "Treatment", pt.size = 0, ncol = 5)
  dev.off()
}

cat("\n=== NMF analysis complete for ", celltype_arg, " ===\n")
