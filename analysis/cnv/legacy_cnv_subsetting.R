####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cnv/legacy_cnv_subsetting.R
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Historical one-off Ju_2025 reference-spiking CNV heatmap; retained without execution because it contains fixed sampling and incomplete scratch code.
#   Inputs: historical Ju_2025 reference objects and cancer_signatures.txt.
#   Outputs: a historical inferCNV-style PDF in the current working directory.
#   Run: not part of the active workflow; use cnv_profiling.R and cnv_malignant_subclone_mp_heatmap.R instead.
####################

####################
# Moved from: analysis/cnv_subset.R
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


library(infercna)
library(dplyr)
library(ggplot2)

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

all_meta <- list()
all_outs <- list()

for (sample in sample_dirs) {
  
  epi_file  <- file.path("by_samples", sample, paste0(sample, "_epi_f.rds"))
  outs_file <- file.path("by_samples", sample, paste0(sample, "_outs.rds"))
  
  # skip if outs file is missing
  if (!file.exists(outs_file)) {
    message("Skipping ", sample, " — no _outs.rds file")
    next
  }
  
  epi  <- readRDS(epi_file)
  outs <- readRDS(outs_file)
  meta <- epi@meta.data
  outs <- outs[, rownames(meta)]
  
  # rename malignancy levels
  meta$status <- plyr::revalue(
    meta$malignancy,
    c(
      "malignant_level_1"        = "malignant",
      "malignant_level_2"        = "malignant",
      "non_malignant_level_1"    = "unresolved",
      "non_malignant_level_2"    = "unresolved",
      "malignant_unresolved"     = "unresolved",
      "non_malignant_unresolved" = "unresolved"
    )
  )
  
  # store results
  all_meta[[sample]] <- meta
  all_outs[[sample]] <- outs
  print(paste0("Finished processing sample: ", sample))
}

# 1. Find intersection of rownames across all outs
common_rows <- Reduce(intersect, lapply(all_outs, rownames))
all_outs_intersect <- lapply(all_outs, function(m) m[common_rows, , drop = FALSE])
common_cols <- Reduce(intersect, lapply(all_meta, colnames))
all_meta_intersect <- lapply(all_meta, function(m) m[, common_cols, drop = FALSE])
meta <- do.call(rbind, all_meta_intersect)
outs <- do.call(cbind, all_outs_intersect)


## ---- load ----
ref <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/by_samples/Ju_2025_patient1/Ju_2025_patient1_outs.rds")
ref_all <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Ju_2025_reference.rds")
EAC_Ref_merged <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_merged.rds")

## ---- cells to keep ----
ref_cells <- unlist(ref_all$ref, use.names = FALSE)
ref <- ref[, ref_cells]

## ---- expression matrix (genes x cells) ----
expr_mat <- as.matrix(EAC_Ref_merged@assays$RNA$data[, ref_cells, drop = FALSE])

## ---- CC score (top 10 consensus cell cycle genes by mean expression) ----
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE, stringsAsFactors = FALSE
)[, 1:3]

cc_genes <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(expr_mat))
stopifnot(length(cc_genes) > 0)

cc_top10 <- c(
  "ARL6IP1", "UBE2S", "SLBP", "DUT", "ECT2",
  "KPNA2", "DNMT1", "GMNN", "PCNA", "HELLS"
)
cc_score <- colMeans(expr_mat[cc_top10, , drop = FALSE])

## ---- CS score (top 50 signature genes by mean expression) ----
cs_genes <- read.table("cancer_signatures.txt", header = FALSE, stringsAsFactors = FALSE)$V1
cs_genes <- intersect(cs_genes, rownames(expr_mat))
stopifnot(length(cs_genes) > 0)

cs_top50 <- c(
  "PTPRK", "AKAP13", "MAGI1", "MCU", "EXT1", "MECOM",
  "LMO7", "LPP", "ARHGAP26", "CYP3A5", "RBM47", "KLF5",
  "HSPA1B", "DNAJB1", "ANKRD36C", "HSPH1", "JMJD1C", "TCF7L2",
  "TMC5", "FNDC3B", "PPP1R15A", "FOSB", "EPS8", "SHROOM3",
  "ZBTB20", "KLF6", "FOXP1", "PRKCA", "ZFP36L2", "GAREM1",
  "CTNND1", "SYNE2", "CHD2", "ERRFI1", "ATF3", "CEMIP2",
  "TJP1", "SIK3", "PARD3", "AKAP9", "SLC38A2", "JUNB",
  "LRBA", "SIPA1L1", "SH3RF1", "RAB10", "SMURF1", "EGR1",
  "IGF2BP2", "CDC42BPA"
)
cs_score <- colMeans(expr_mat[cs_top50, , drop = FALSE])

## ---- CNA scatter scores (optional; keep as separate columns) ----
coord <- cnaScatterPlot(ref)
studies <- c(
  "Alcindor_2025", "Ju_2025", "Carroll_2023", "Strasser_2025",
  "Baek_2025", "Yates_2025", "Walker_2025", "Croft_2022"
)

meta_ref <- data.frame(
  cell       = ref_cells,
  cs_score   = cs_score[ref_cells],
  cc_score   = cc_score[ref_cells],
  cna_signal = coord$cna.signal, 
  cna.cor    = coord$cna.cor, 
  status     = rep("reference", length(ref_cells)), 
  orig.ident      = paste0(sample(studies, length(ref_cells), replace = TRUE), "_XXX"), 
  row.names  = ref_cells,
  stringsAsFactors = FALSE
)
meta_ref <- meta_ref[sample(rownames(meta_ref), 810), ]
ref <- ref[, rownames(meta_ref)]

gene_order <- read.tarefgene_order <- read.table(
  "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")
)

save <- outs
meta_save <- meta

meta_save$status <- ifelse(meta_save$status == "malignant", "malignant", "unresolved")
meta_save$status <- ifelse(meta_save$malignancy %in% c("malignant_level_1", "malignant_level_2"), meta_save$malignancy, "unresolved")

set.seed(1)
ratio <- c(malignant_level_1 = 4000, malignant_level_2 = 3000, unresolved = 810)
ratio <- ratio / sum(ratio)

n_total <- 7810
n_per_group <- round(n_total * ratio)
statuses <- names(n_per_group)

keep_cells <- unlist(
  lapply(statuses, function(s) {
    idx <- which(meta_save$status == s)
    n_target <- n_per_group[[s]]
    
    if (length(idx) >= n_target) sample(idx, n_target) else idx
  }),
  use.names = FALSE
)

meta <- meta_save[keep_cells, ]
outs <- save[, keep_cells]


common_cols <- intersect(colnames(meta_save), colnames(meta_ref))
meta <- rbind(
  meta[, common_cols, drop = FALSE],
  meta_ref[,  common_cols, drop = FALSE]
)

common_rows <- intersect(rownames(outs), rownames(ref))
outs <- cbind(
  outs[common_rows, , drop = FALSE],
  ref[common_rows, , drop = FALSE]
)

meta$study <- sapply(
  strsplit(meta$orig.ident, "_"),
  function(x) paste(x[1:2], collapse = "_")
)
study_levels <- unique(meta$study)
study_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(study_levels)),
  study_levels
)

chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
common_genes <- intersect(rownames(outs), gene_order$gene_id)

go <- gene_order %>%
  filter(gene_id %in% common_genes, chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(chromosome, levels = chrom_levels)) %>%
  arrange(chromosome, start)

outs <- outs[go$gene_id, , drop = FALSE]
stopifnot(identical(rownames(outs), go$gene_id))

## --------------------------- BINNING ----------------------------- ##
bin_size <- 100L

go <- go %>%
  group_by(chromosome) %>%
  mutate(
    g_rank = row_number(),
    bin_in_chr = ((g_rank - 1L) %/% bin_size) + 1L,
    bin_key = paste(chromosome, bin_in_chr, sep = "_")
  ) %>%
  ungroup()


ordered_bin_keys <- unique(go$bin_key)
bins_idx <- split(seq_len(nrow(go)), factor(go$bin_key, levels = ordered_bin_keys))
# ----------------------------- FIX 1 END ---------------------------- ##

binned_mat <- do.call(rbind, lapply(bins_idx, function(ix) colMeans(outs[ix, , drop = FALSE])))
rownames(binned_mat) <- names(bins_idx)

# chromosome per bin, in row order; keep order-of-appearance to match matrix
row_chr_labels <- sub("_.*$", "", rownames(binned_mat))
row_chr <- factor(row_chr_labels, levels = unique(row_chr_labels))

top_ha <- HeatmapAnnotation(
  cancer_signature = as.numeric(meta$cs_score),
  cell_cycling     = as.numeric(meta$cc_score),
  study            = factor(meta$study),
  
  col = list(
    cancer_signature = colorRamp2(c(0,3,4),   c("white","grey80","black")),
    cell_cycling     = colorRamp2(c(0,1,1.5), c("white","grey80","black")),
    study            = study_colors
  ),
  
  annotation_name_side   = "left",
  annotation_name_gp     = gpar(fontsize = 9),      # optional: slightly smaller
  annotation_name_offset = unit(2, "mm"),           # adds spacing from bars

  # makes each annotation row taller
  annotation_height = unit(c(4, 4, 4), "mm"),
  
  show_legend = c(cancer_signature = TRUE, cell_cycling = TRUE, study = TRUE),
  annotation_legend_param = list(
    cancer_signature = list(title = "CS score"),
    cell_cycling     = list(title = "CC score"),
    study            = list(title = "Study")
  )
)


## ---------------------- CHR COLOR BAR (LEFT) --------------------- ##
chr_used <- levels(droplevels(row_chr))
base_cols <- c(brewer.pal(12, "Paired"),
               brewer.pal(8,  "Dark2"),
               brewer.pal(9,  "Set1"),
               brewer.pal(12, "Set3"))
chr_cols <- setNames(base_cols[seq_along(chr_used)], chr_used)

left_chr_bar <- rowAnnotation(
  chr = row_chr,
  col = list(chr = chr_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  gp = gpar(col = NA),
  width = unit(4, "mm")
)

## ---------------------- BOUNDARIES --------------------------------##
# Horizontal boundaries from actual row order
chr_bounds <- which(head(row_chr_labels, -1L) != tail(row_chr_labels, -1L))

class_order <- c("malignant_level_1", "malignant_level_2", "unresolved", "reference")

# make sure classification has the right order
meta$classification <- factor(meta$status, levels = class_order)

# use it directly (no need to pre-order)
col_split <- meta$classification
line_gp <- gpar(col = "black", lwd = 2, lineend = "square")

ht <- Heatmap(
  binned_mat,
  name = "CNV",
  cluster_rows = FALSE,
  
  cluster_columns = TRUE,
  column_split = col_split, 
  column_title_rot = 30, 
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  column_gap = unit(2, "mm"),
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation  = top_ha,
  left_annotation = left_chr_bar,
  
  row_split = row_chr,
  row_gap = unit(0, "mm"),
  row_title_rot  = 0,
  rect_gp = gpar(col = NA),
  border = NA,
  
  layer_fun = function(j, i, x, y, w, h, fill) {
    hits <- intersect(i, chr_bounds)
    if (length(hits)) {
      id <- match(hits, i)
      yy <- y[id] - h[id]/2
      grid.segments(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = yy, y1 = yy,
        gp = line_gp
      )
    }
  }
)


## ----------------------------- SAVE -------------------------------##
pdf(
  paste0(unique(meta$ident)[1], "_T_CNV_profile_infercnv_style4.pdf"),
  width = 8, height = 8
)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
