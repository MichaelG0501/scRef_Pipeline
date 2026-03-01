####################
# Moved from: analysis/plot_CNV.R
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
sample_dirs <- sample_dirs[grepl("^[^/]+_PDO$", sample_dirs)][-21]

all_outs <- list()

for (sample in sample_dirs) {
  
  outs_file <- file.path("by_samples", sample, paste0(sample, "_outs.rds"))

  outs <- readRDS(outs_file)
  outs <- outs[, grepl("PDO", colnames(outs)), drop = FALSE]

  all_outs[[sample]] <- outs
  print(paste0("Finished processing sample: ", sample))
}

# 1. Find intersection of rownames across all outs
common_rows <- Reduce(intersect, lapply(all_outs, rownames))
all_outs_intersect <- lapply(all_outs, function(m) m[common_rows, , drop = FALSE])
outs <- do.call(cbind, all_outs_intersect)

ref <- readRDS(outs_file)
ref <- ref[, !grepl("PDO", colnames(ref)), drop = FALSE]

outs <- cbind(outs, ref)

gene_order <- read.tarefgene_order <- read.table(
  "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")
)

save <- outs
set.seed(123)  # for reproducibility
cell_idx <- sample(ncol(outs), 3000)
outs <- outs[, cell_idx]

sample_ids <- sub("_[^_]+$", "", colnames(outs))
sample_ids2 <- ifelse(grepl("^Strasser_2025_", sample_ids), "reference", sample_ids)
study_levels <- unique(sample_ids2)
study_colors <- setNames(
  hcl.colors(21, palette = "Dark 3"),
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
bin_size <- 200L

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
  sample            = factor(sample_ids2),
  
  col = list(
    sample            = study_colors
  ),
  
  annotation_name_side   = "left",
  annotation_name_gp     = gpar(fontsize = 16),      # optional: slightly smaller
  annotation_name_offset = unit(2, "mm"),           # adds spacing from bars
  
  # makes each annotation row taller
  annotation_height = unit(c(4, 4, 4), "mm"),
  annotation_legend_param = list(
    sample            = list(title = "Sample")
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


classification <- ifelse(grepl("^Strasser_2025_", sample_ids), "Reference", "PDOs")
class_order <- c("PDOs", "Reference")

classification <- factor(classification, levels = class_order)

# use it directly (no need to pre-order)
col_split <- classification
line_gp <- gpar(col = "black", lwd = 2, lineend = "square")

ht <- Heatmap(
  binned_mat,
  name = "CNV",
  cluster_rows = FALSE,
  
  cluster_columns = TRUE,
  column_split = col_split, 
  column_title_rot = 0, 
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
pdf("PDOs_cnv.pdf", width = 14, height = 8)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
