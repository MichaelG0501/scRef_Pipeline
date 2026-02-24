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


outs <- readRDS("all_outs.rds")

# coord <- cnaScatterPlot(outs)
# 
# saveRDS(coord, "coord.rds")

load("meta.RData")


cell_summary <- as.data.frame(cell_id = rownames(meta), ident = meta$orig.ident, 
                              celltypes = meta$celltype_manual, cna_signal = coord$cna.signal, 
                              cna_cor = coord$cna.cor)


library(readxl)
data <- read_excel("/rds/general/project/tumourheterogeneity1/ephemeral/EAC_Ref_all/00_merged/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
data <- data %>%
  mutate(ident = paste(Author, Year, `Sample Name`, sep = "_"))
cell_names <- rownames(meta)
meta <- meta %>%
  left_join(data, by = "ident")
rownames(meta) <- cell_names


gene_order <- read.table(
  "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")
)

meta <- readRDS("cell_summary.rds")
outs <- readRDS("all_filtered_outs.rds")
save <- outs

## -------------------------- PREPROCESS --------------------------- ##
outs <- outs[, rownames(meta)[meta$celltypes == "epithelial" & meta$T_Status == "Tumour" &
                                meta$classification == "malignant"], drop = FALSE]

meta <- meta[colnames(outs), , drop = FALSE]
filter_levels <- c("malignant", "nonmalignant", "unresolved")
filter_colors <- setNames(c("red", "green", "grey60"), filter_levels)
meta$classification <- factor(meta$classification, levels = filter_levels)

celltype_levels <- unique(meta$celltypes)
celltype_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(celltype_levels)),
  celltype_levels
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

# ---------------------------- FIX 1 START --------------------------- ##
# PROBLEM: Chromosome sorting. `split()` on a character vector sorts alphabetically
# (e.g., chr1, chr10, chr11, chr2...).
# SOLUTION: Convert `bin_key` to a factor with levels ordered correctly based on
# the already-sorted `go` data frame. This ensures `split` respects the
# desired numeric chromosome order.
ordered_bin_keys <- unique(go$bin_key)
bins_idx <- split(seq_len(nrow(go)), factor(go$bin_key, levels = ordered_bin_keys))
# ----------------------------- FIX 1 END ---------------------------- ##

binned_mat <- do.call(rbind, lapply(bins_idx, function(ix) colMeans(outs[ix, , drop = FALSE])))
rownames(binned_mat) <- names(bins_idx)

# chromosome per bin, in row order; keep order-of-appearance to match matrix
row_chr_labels <- sub("_.*$", "", rownames(binned_mat))
row_chr <- factor(row_chr_labels, levels = unique(row_chr_labels))

## -------------------- TOP (COLUMN) ANNOTATIONS ------------------- ##
cell_clust <- hclust(dist(t(binned_mat)), method = "ward.D2")
saveRDS(cell_clust, "malignant_T_clust.RDS")
#cell_clust <- readRDS("cynthia_H_T1_clust.RDS")
k <- 3L
meta$cluster <- factor(cutree(cell_clust, k = k))

cluster_mean_signal <- tapply(meta$cna_signal, meta$cluster, mean, na.rm = TRUE)
cluster_mean_cor    <- tapply(meta$cna_cor,    meta$cluster, mean, na.rm = TRUE)
meta$mean_signal_per_cluster <- cluster_mean_signal[meta$cluster]
meta$mean_cor_per_cluster    <- cluster_mean_cor[meta$cluster]

cluster_ids  <- sort(unique(meta$cluster))
cluster_cols <- setNames(
  brewer.pal(max(3, length(cluster_ids)), "Set2")[seq_along(cluster_ids)],
  cluster_ids
)

top_ha <- HeatmapAnnotation(
  cluster         = meta$cluster,
  mean_cna_signal = meta$mean_signal_per_cluster,
  mean_cna_cor    = meta$mean_cor_per_cluster,
  filter_status   = meta$classification,
  celltype_group  = meta$celltypes, 
  ident = meta$ident, 
  col = list(
    cluster         = cluster_cols,
    mean_cna_signal = colorRamp2(c(0, max(meta$mean_signal_per_cluster, na.rm = TRUE)), c("white", "red")),
    mean_cna_cor    = colorRamp2(c(0, max(meta$mean_cor_per_cluster,    na.rm = TRUE)), c("white", "blue")),
    filter_status   = filter_colors,
    celltype_group  = celltype_colors
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    cluster         = list(title = "Cluster"),
    mean_cna_signal = list(title = "Mean CNA Signal"),
    mean_cna_cor    = list(title = "Mean CNA Correlation"),
    filter_status   = list(title = "Filter Status"),
    celltype_group  = list(title = "Cell Type")
  ),
  simple_anno_size = unit(2.5, "mm")
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

# Column order from clustering
cell_order <- cell_clust$order
clusters_in_order <- as.character(meta$cluster)[cell_order]

# Vertical boundaries from actual cluster transitions in this order
cluster_bounds <- which(head(clusters_in_order, -1L) != tail(clusters_in_order, -1L))


# ---------------------------- FIX 2 START --------------------------- ##
# PROBLEM: Line visibility. The original lines were too thin.
# SOLUTION: Increased `lwd` (line width) from 1.2 to 2 for much clearer,
# more visible boundary lines.
line_gp <- gpar(col = "black", lwd = 2, lineend = "square")
# ----------------------------- FIX 2 END ---------------------------- ##


## ----------------------------- HEATMAP --------------------------- ##
ht <- Heatmap(
  binned_mat,
  name = "CNV",
  #  col = cnv_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = cell_order,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation  = top_ha,
  left_annotation = left_chr_bar,
  
  row_split = row_chr,
  row_title_side = "left",
  row_title_gp   = gpar(fontsize = 7),
  row_title_rot  = 0,
  
  row_gap = unit(0, "mm"),
  column_gap = unit(0, "mm"),
  rect_gp = gpar(col = NA),
  use_raster = FALSE,
  border = NA,
  
  # ---------------------------- FIX 3 START --------------------------- ##
  # PROBLEM: Incorrect cluster boundary alignment. The original logic incorrectly
  # mixed indices from the original matrix with indices from the ordered matrix.
  # SOLUTION: Redraw the vertical lines using "npc" (Normalized Parent Coordinates).
  # We calculate the fractional position of each cluster boundary across the
  # *entire* heatmap width. This is robust and ensures perfect alignment.
  layer_fun = function(j, i, x, y, w, h, fill) {
    ## Horizontal chromosome lines — bottom edge between row groups
    hits_r <- intersect(i, chr_bounds)
    if (length(hits_r) > 0) {
      id <- match(hits_r, i)
      yy <- y[id] - h[id] / 2
      grid.segments(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                    y0 = yy, y1 = yy, gp = line_gp)
    }
    
    ## Vertical cluster lines
    total_cols <- length(cell_order) # Total number of columns in the heatmap
    if (length(cluster_bounds) > 0) {
      # Convert boundary indices to fractional positions (0 to 1)
      boundary_positions_npc <- cluster_bounds / total_cols
      # Draw segments at these fractional positions
      grid.segments(
        y0 = unit(0, "npc"), y1 = unit(1, "npc"),
        x0 = unit(boundary_positions_npc, "npc"),
        x1 = unit(boundary_positions_npc, "npc"),
        gp = line_gp
      )
    }
  }
  # ----------------------------- FIX 3 END ---------------------------- ##
)

## ----------------------------- SAVE -------------------------------##
pdf(
  paste0(unique(meta$ident)[1], "_T_CNV_profile_infercnv_style.pdf"),
  width = 15, height = 8
)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


#############################################################

library(dplyr)
library(purrr)

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
    ifelse(epi$cna_cor > thr_cor & epi$cna_signal > thr_sig, "malignant",
           ifelse(epi$cna_cor > thr_cor | epi$cna_signal > thr_sig, "cna_unresolved", "non_malignant")
    )
  
  df_ident
}

meta <- meta %>%
  mutate(
    cna_cor    = as.numeric(cna_cor),
    cna_signal = as.numeric(cna_signal)
  ) %>%
  group_split(ident, .keep = TRUE) %>%
  map_df(classify_one_ident)

meta <- meta[meta$T_Status == "Tumour" & meta$celltype_update == "epithelial", ]

#meta$classification <- ifelse(meta$cna_cor > 0.3 & meta$cna_signal > 0.001, "malignant", "non_malignant")

summary_tbl <- meta %>%
  group_by(study, orig.ident) %>%
  summarise(
    n_total = n(),
    n_malignant = sum(classification == "malignant", na.rm = TRUE),
    pct_malignant = 100 * n_malignant / n_total,
    .groups = "drop"
  )

print(summary_tbl)

library(ggplot2)

p <- ggplot(summary_tbl, aes(x = orig.ident, y = pct_malignant, fill = study)) +
  geom_col() +
  geom_text(aes(label = n_total), vjust = -0.3, size = 3) +
  facet_wrap(~ study, scales = "free_x") +
  labs(
    x = NULL,
    y = "% Malignant cells",
    title = "Fraction of malignant cells per sample by study"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave("malignant_fraction_plot.png", plot = p, width = 10, height = 6, dpi = 300)

##########################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

## ---------- build the working df (your code, kept simple) ----------
df <- data.frame(
  x = meta$cna_cor, 
  y = meta$cna_signal,
  color = meta$celltype_manual,
  color2 = meta$orig.ident,
  study = meta$study
)

celltype_counts <- df %>% count(color) %>% arrange(desc(n))
df <- df %>%
  mutate(color = factor(color, levels = celltype_counts$color)) %>%
  arrange(color)

df <- df %>%
  mutate(color2 = factor(color2, levels = sort(unique(color2)))) %>%
  arrange(color2)

## ---------- per-sample summary label (same logic as yours) ----------
make_summary <- function(d) {
  d %>%
    mutate(pass_threshold = (y > 0.001) & (x > 0.3)) %>%
    summarise(
      n_total = sum(color %in% c("epithelial")),
      n_pass  = sum(pass_threshold & color %in% c("epithelial")),
      pct_pass = ifelse(n_total > 0, n_pass / n_total * 100, 0),
      mean_x = mean(x, na.rm = TRUE),
      mean_y = mean(y, na.rm = TRUE)
    ) %>% 
    mutate(label = paste0("Passed : ", n_pass, " (", sprintf("%.3f", pct_pass), "%)"))
}

## ---------- one sample plot (only epithelial) ----------
plot_sample <- function(d_samp, smp) {
  d_epi <- d_samp %>% filter(color %in% c("epithelial"))
  lab_df <- make_summary(d_samp)
  
  ggplot(d_epi, aes(x = y, y = x)) +
    geom_point(size = 0.7, alpha = 1) +
    labs(
      x = "CNA Signal", y = "CNA Correlation",
      title = smp
    ) +
    theme_minimal() +
    geom_vline(xintercept = 0.001, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0.3, color = "red", linetype = "dashed") +
    geom_text(
      data = lab_df,
      aes(x = 0.008, y = 0.6, label = label),
      color = "black", size = 3.8, hjust = 0, vjust = 1
    ) +
    xlim(c(-0.001, 0.005))
}

## ---------- samples per study, 2x5 per page ----------
# distinct sample–study map
samples_tbl <- df %>% distinct(color2, study)
by_study <- split(samples_tbl$color2, samples_tbl$study)

pdf("all.pdf", width = 16, height = 9, onefile = TRUE)

for (stud in names(by_study)) {
  smps <- sort(by_study[[stud]])
  if (length(smps) == 0) next
  
  # chunk into pages of 10
  idx <- seq(1, length(smps), by = 10)
  for (s in idx) {
    chunk <- smps[s:min(s + 9, length(smps))]
    
    # make plots for this page
    grobs <- lapply(chunk, function(smp) {
      d_samp <- df %>% filter(color2 == smp)
      plot_sample(d_samp, smp)
    })
    
    # pad to exactly 10 if needed
    if (length(grobs) < 10) {
      blanks <- replicate(10 - length(grobs), ggplot() + theme_void(), simplify = FALSE)
      grobs <- c(grobs, blanks)
    }
    
    # arrange 2 rows x 5 columns, with a study-level title
    grid.arrange(
      grobs = grobs, nrow = 2, ncol = 5,
      top = textGrob(paste0("Study: ", stud), gp = gpar(fontface = "bold", cex = 1.2))
    )
  }
}

dev.off()
