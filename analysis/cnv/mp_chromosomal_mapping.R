####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/mp_chromosomal_mapping.R
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_mp_chromosomal_mapping.R
#
# Visualise metaprogramme gene locations across chromosomes/arms.
# Produces 5 plot styles in a single PDF.
#
# Inputs:
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt
#
# Outputs:
#   ref_outs/Auto_mp_chromosomal_mapping.pdf
#   ref_outs/Auto_mp_chromosomal_mapping_summary.csv
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(scales)
  library(RColorBrewer)
})

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ── 1. Load data ──
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
cat("Retained MPs:", length(mp.genes), "\n")

gene_order <- read.table(
  "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  header = FALSE, col.names = c("gene", "chromosome", "start", "end"),
  stringsAsFactors = FALSE
)
chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- gene_order %>%
  filter(chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(chromosome, levels = chrom_levels)) %>%
  arrange(chromosome, start) %>%
  distinct(gene, .keep_all = TRUE)

# ── 2. MP descriptions & colours ──
mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",       "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC Proliferation",    "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted",      "MP5"  = "Epithelial IFN",
  "MP10" = "Columnar Diff.",       "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam.",      "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory (Intest.)",  "MP16" = "Secretory (Gastric)",
  "MP15" = "Immune Infiltration",  "MP12" = "Neuro-responsive"
)

# Fixed colour palette — 14 distinct colours
mp_pal <- c(
  "MP1"  = "#B0B0B0", "MP9"  = "#C0C0C0", "MP7"  = "#999999",
  "MP2"  = "#E41A1C", "MP17" = "#4DAF4A", "MP14" = "#8DA0CB",
  "MP5"  = "#66C2A5", "MP10" = "#A6D854", "MP8"  = "#FC8D62",
  "MP18" = "#FF7F00", "MP16" = "#FFD92F", "MP13" = "#984EA3",
  "MP12" = "#E78AC3", "MP15" = "#377EB8"
)

mp_order <- c("MP1","MP9","MP7","MP2","MP17","MP14","MP5","MP10",
              "MP8","MP18","MP16","MP13","MP12","MP15")
mp_order <- intersect(mp_order, names(mp.genes))

label_mp <- function(x) {
  d <- mp_descriptions[x]; d[is.na(d)] <- x[is.na(d)]
  paste0(x, " ", d)
}

# ── 3. Map genes to coordinates ──
centromere_pos <- c(
  chr1=121700000, chr2=91800000, chr3=87900000, chr4=50600000,
  chr5=48400000, chr6=61000000, chr7=59900000, chr8=45600000,
  chr9=49000000, chr10=40200000, chr11=53400000, chr12=35500000,
  chr13=17700000, chr14=17200000, chr15=19000000, chr16=36800000,
  chr17=25100000, chr18=18500000, chr19=26200000, chr20=28100000,
  chr21=12000000, chr22=15000000, chrX=61000000
)

mp_gene_df <- do.call(rbind, lapply(names(mp.genes), function(mp) {
  genes <- mp.genes[[mp]]
  matched <- gene_order[gene_order$gene %in% genes, ]
  if (nrow(matched) == 0) return(NULL)
  matched$mp <- mp
  matched$mp_label <- label_mp(mp)
  matched
}))
mp_gene_df$arm <- ifelse(mp_gene_df$start <= centromere_pos[as.character(mp_gene_df$chromosome)], "p", "q")
mp_gene_df$chr_arm <- paste0(mp_gene_df$chromosome, mp_gene_df$arm)
mp_gene_df$chr_num <- as.integer(gsub("chr|X", "", gsub("chrX", "23", mp_gene_df$chromosome)))
mp_gene_df$midpoint <- (mp_gene_df$start + mp_gene_df$end) / 2
mp_gene_df$mp <- factor(mp_gene_df$mp, levels = mp_order)

cat("Total mapped genes:", nrow(mp_gene_df), "\n")
cat("Unmapped genes per MP:\n")
for (mp in mp_order) {
  total <- length(mp.genes[[mp]])
  mapped <- sum(mp_gene_df$mp == mp, na.rm = TRUE)
  cat("  ", mp, ":", total - mapped, "of", total, "unmapped\n")
}

# ── 4. Chromosome arm summary ──
chr_arm_levels <- paste0(rep(chrom_levels, each = 2), c("p","q"))
arm_summary <- mp_gene_df %>%
  group_by(mp, chr_arm) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(
    mp_gene_df %>% group_by(mp) %>% summarise(total = n(), .groups = "drop"),
    by = "mp"
  ) %>%
  mutate(pct = 100 * n / total)

# Fisher's exact test for enrichment
total_genes_in_ref <- nrow(gene_order)
arm_gene_counts <- gene_order %>%
  mutate(arm = ifelse(start <= centromere_pos[as.character(chromosome)], "p", "q"),
         chr_arm = paste0(chromosome, arm)) %>%
  count(chr_arm, name = "arm_total")

enrich_df <- arm_summary %>%
  left_join(arm_gene_counts, by = "chr_arm") %>%
  left_join(
    mp_gene_df %>% group_by(mp) %>% summarise(mp_total = n(), .groups = "drop"),
    by = "mp"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch({
      mat <- matrix(c(n, arm_total - n, mp_total - n, total_genes_in_ref - arm_total - mp_total + n),
                    nrow = 2)
      mat[mat < 0] <- 0
      fisher.test(mat, alternative = "greater")$p.value
    }, error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         log10p = -log10(p_adj),
         sig = ifelse(!is.na(p_adj) & p_adj < 0.05, "*", ""))

write.csv(
  enrich_df %>% select(mp, chr_arm, n, total, pct, p_value, p_adj, sig) %>% arrange(mp, chr_arm),
  "Auto_mp_chromosomal_mapping_summary.csv", row.names = FALSE
)

# ── 5. Chromosome sizes (approx hg38) ──
chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468, chrX=156040895
)

# ══════════════════════════════════════════════════════════
# PLOT 1: Circos plot
# ══════════════════════════════════════════════════════════
plot_circos <- function() {
  # Prepare circos sectors
  sector_df <- data.frame(
    chr = chrom_levels,
    start = 0,
    end = chr_sizes[chrom_levels],
    stringsAsFactors = FALSE
  )

  circos.clear()
  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.005, 0.005),
    start.degree = 90,
    gap.degree = c(rep(1, 22), 4),
    track.height = 0.06
  )

  circos.genomicInitialize(
    sector_df,
    plotType = c("axis", "labels"),
    axis.labels.cex = 0.35,
    labels.cex = 0.55
  )

  # Centromere markers
  circos.track(
    ylim = c(0, 1), track.height = 0.03, bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      cen <- centromere_pos[chr]
      if (!is.na(cen)) {
        circos.rect(cen - 2e6, 0, cen + 2e6, 1, col = "#333333", border = NA)
      }
    }
  )

  # One track per MP (stacked)
  for (mp in mp_order) {
    sub <- mp_gene_df[mp_gene_df$mp == mp, ]
    col <- mp_pal[mp]

    circos.track(
      ylim = c(0, 1), track.height = 0.035, bg.border = NA, bg.col = "grey97",
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        genes_here <- sub[as.character(sub$chromosome) == chr, ]
        if (nrow(genes_here) > 0) {
          for (i in seq_len(nrow(genes_here))) {
            circos.rect(
              genes_here$start[i], 0.05,
              genes_here$end[i], 0.95,
              col = col, border = NA
            )
          }
        }
      }
    )
  }

  # Legend
  mp_labels_legend <- label_mp(mp_order)
  lgd <- Legend(
    labels = mp_labels_legend,
    legend_gp = gpar(fill = mp_pal[mp_order]),
    title = "Metaprogramme",
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 5),
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    ncol = 2
  )
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.13, "npc"), just = "center")

  title("Metaprogramme Gene Locations - Circos View", cex.main = 0.9, line = -1)
  circos.clear()
}

# ══════════════════════════════════════════════════════════
# PLOT 2: Chromosome arm enrichment heatmap
# ══════════════════════════════════════════════════════════
plot_arm_heatmap <- function() {
  # Build matrix: MP x chr_arm (% of MP genes)
  hmat <- enrich_df %>%
    select(mp, chr_arm, pct) %>%
    pivot_wider(names_from = chr_arm, values_from = pct, values_fill = 0) %>%
    tibble::column_to_rownames("mp") %>%
    as.matrix()

  # Filter to arms with at least one gene
  keep_arms <- colnames(hmat)[colSums(hmat) > 0]
  # Order arms properly
  keep_arms <- intersect(chr_arm_levels, keep_arms)
  hmat <- hmat[mp_order[mp_order %in% rownames(hmat)], keep_arms, drop = FALSE]
  rownames(hmat) <- label_mp(rownames(hmat))

  # Significance overlay
  sig_mat <- enrich_df %>%
    select(mp, chr_arm, sig) %>%
    pivot_wider(names_from = chr_arm, values_from = sig, values_fill = "") %>%
    tibble::column_to_rownames("mp") %>%
    as.matrix()
  sig_mat <- sig_mat[mp_order[mp_order %in% rownames(sig_mat)], keep_arms, drop = FALSE]
  rownames(sig_mat) <- label_mp(rownames(sig_mat))

  # Chr grouping
  arm_chr <- sub("[pq]$", "", keep_arms)
  arm_chr_short <- gsub("chr", "", arm_chr)

  col_split <- factor(arm_chr, levels = unique(arm_chr))

  ht <- Heatmap(
    hmat,
    name = "% of\nMP genes",
    col = colorRamp2(c(0, 5, 15, 30), c("white", "#FEE08B", "#F46D43", "#A50026")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = col_split,
    column_gap = unit(0.5, "mm"),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 60,
    row_names_side = "left",
    column_title_gp = gpar(fontsize = 7),
    column_title_rot = 0,
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3),
    cell_fun = function(j, i, x, y, w, h, fill) {
      v <- hmat[i, j]
      s <- sig_mat[i, j]
      if (!is.na(v) && v >= 1) {
        grid.text(sprintf("%.0f", v), x, y, gp = gpar(fontsize = 5, col = ifelse(v > 15, "white", "black")))
      }
      if (!is.na(s) && nzchar(s)) {
        grid.text(s, x, unit.c(y + unit(2.5, "mm")), gp = gpar(fontsize = 7, col = "black", fontface = "bold"))
      }
    },
    width = unit(min(250, 5.5 * length(keep_arms)), "mm"),
    height = unit(5.5 * nrow(hmat), "mm"),
    column_title = "Chromosome Arm"
  )

  draw(ht, column_title = "MP Gene Distribution Across Chromosome Arms (% + Fisher enrichment *)",
       column_title_gp = gpar(fontsize = 10, fontface = "bold"),
       padding = unit(c(5, 5, 5, 5), "mm"))
}

# ══════════════════════════════════════════════════════════
# PLOT 3: Karyotype-style ideogram
# ══════════════════════════════════════════════════════════
plot_karyotype <- function() {
  # Build chromosome backbone
  chr_backbone <- data.frame(
    chr = factor(chrom_levels, levels = chrom_levels),
    chr_num = seq_along(chrom_levels),
    size = chr_sizes[chrom_levels],
    centromere = centromere_pos[chrom_levels]
  )

  # Offset each MP within a chromosome band
  mp_gene_df2 <- mp_gene_df %>%
    mutate(mp_idx = as.integer(factor(mp, levels = mp_order)),
           y_offset = (mp_idx - 1) / length(mp_order))

  chr_labels <- gsub("chr", "", chrom_levels)

  p <- ggplot() +
    # Chromosome rectangles
    geom_rect(data = chr_backbone,
              aes(xmin = 0, xmax = size, ymin = chr_num - 0.4, ymax = chr_num + 0.4),
              fill = "grey95", color = "grey60", linewidth = 0.3) +
    # Centromere
    geom_segment(data = chr_backbone,
                 aes(x = centromere, xend = centromere,
                     y = chr_num - 0.4, yend = chr_num + 0.4),
                 color = "grey30", linewidth = 0.6, linetype = "dashed") +
    # Gene tick marks (colored by MP)
    geom_segment(data = mp_gene_df2,
                 aes(x = midpoint, xend = midpoint,
                     y = chr_num - 0.38, yend = chr_num + 0.38,
                     color = mp),
                 linewidth = 0.4, alpha = 0.8) +
    scale_color_manual(values = mp_pal[mp_order], labels = label_mp(mp_order),
                       name = "Metaprogramme") +
    scale_y_continuous(breaks = seq_along(chrom_levels), labels = chr_labels,
                       trans = "reverse") +
    scale_x_continuous(labels = function(x) paste0(x / 1e6, "Mb"),
                       breaks = seq(0, 250e6, 50e6)) +
    labs(title = "Metaprogramme Gene Positions - Karyotype View",
         x = "Genomic position", y = "Chromosome") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      legend.key.size = unit(3, "mm"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(color = guide_legend(ncol = 3, override.aes = list(linewidth = 2, alpha = 1)))

  print(p)
}

# ══════════════════════════════════════════════════════════
# PLOT 4: Per-MP stacked bar (chromosome arm composition)
# ══════════════════════════════════════════════════════════
plot_stacked_bar <- function() {
  # Top arms per MP
  bar_df <- arm_summary %>%
    mutate(mp = factor(mp, levels = mp_order),
           mp_label = label_mp(as.character(mp)))

  # Colour by chromosome (not arm)
  chr_cols <- c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Pastel2"),
    brewer.pal(3, "Dark2")
  )
  chr_short <- gsub("chr", "", chrom_levels)
  names(chr_cols) <- chrom_levels[seq_along(chr_cols)]

  bar_df$chr <- sub("[pq]$", "", bar_df$chr_arm)
  bar_df$arm_label <- gsub("chr", "", bar_df$chr_arm)

  # Keep only arms with >=2 genes for readability
  bar_df_filt <- bar_df %>% filter(n >= 2)

  p <- ggplot(bar_df_filt, aes(x = mp_label, y = n, fill = chr)) +
    geom_col(color = "white", linewidth = 0.2, position = "stack") +
    geom_text(aes(label = ifelse(n >= 3, arm_label, "")),
              position = position_stack(vjust = 0.5), size = 2.2, fontface = "bold") +
    scale_fill_manual(values = chr_cols, name = "Chromosome") +
    coord_flip() +
    labs(title = "Per-MP Gene Count by Chromosome Arm",
         subtitle = "Labels shown for arms with >=3 genes",
         x = NULL, y = "Number of genes") +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.key.size = unit(3, "mm"),
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 7)
    ) +
    guides(fill = guide_legend(nrow = 2))

  print(p)
}

# ══════════════════════════════════════════════════════════
# PLOT 5: Cluster hotspot — top enriched arms per MP
# ══════════════════════════════════════════════════════════
plot_hotspot_lollipop <- function() {
  # Top 3 arms per MP by gene count
  top_arms <- arm_summary %>%
    group_by(mp) %>%
    slice_max(n, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(mp = factor(mp, levels = rev(mp_order)),
           mp_label = label_mp(as.character(mp)),
           arm_label = gsub("chr", "", chr_arm))

  top_arms$mp_label <- factor(top_arms$mp_label, levels = label_mp(rev(mp_order)))

  p <- ggplot(top_arms, aes(x = pct, y = mp_label)) +
    geom_segment(aes(x = 0, xend = pct, yend = mp_label, color = mp),
                 linewidth = 1, show.legend = FALSE) +
    geom_point(aes(size = n, color = mp), show.legend = TRUE) +
    geom_text(aes(label = arm_label), hjust = -0.3, size = 2.5, fontface = "bold") +
    scale_color_manual(values = mp_pal, guide = "none") +
    scale_size_continuous(range = c(2, 7), name = "Gene count") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = "Top Chromosome Arm Hotspots per Metaprogramme",
         subtitle = "Top 3 arms by gene count; labels = chr arm",
         x = "% of MP genes in arm", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_line(colour = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom"
    )

  print(p)
}

# ══════════════════════════════════════════════════════════
# PLOT 6: Per-MP faceted chromosome density ridge
# ══════════════════════════════════════════════════════════
plot_density_facet <- function() {
  # Normalise position to 0-1 within each chromosome
  df <- mp_gene_df %>%
    mutate(chr_size = chr_sizes[as.character(chromosome)],
           norm_pos = midpoint / chr_size,
           chr_label = gsub("chr", "", as.character(chromosome)),
           mp_label = label_mp(as.character(mp)))

  df$mp_label <- factor(df$mp_label, levels = label_mp(mp_order))

  p <- ggplot(df, aes(x = chr_label, fill = mp)) +
    geom_bar(color = "white", linewidth = 0.15, width = 0.7) +
    facet_wrap(~ mp_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = mp_pal, guide = "none") +
    labs(title = "Gene Count per Chromosome - Faceted by Metaprogramme",
         x = "Chromosome", y = "Gene count") +
    theme_minimal(base_size = 8) +
    theme(
      strip.text = element_text(size = 6.5, face = "bold"),
      axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )

  print(p)
}

# ══════════════════════════════════════════════════════════
# Generate all plots
# ══════════════════════════════════════════════════════════
out_pdf <- "Auto_mp_chromosomal_mapping.pdf"
pdf(out_pdf, width = 14, height = 11)

# Page 1: Circos
plot_circos()

# Page 2: Arm heatmap
plot_arm_heatmap()

# Page 3: Karyotype ideogram
plot_karyotype()

# Page 4: Stacked bar
plot_stacked_bar()

# Page 5: Hotspot lollipop
plot_hotspot_lollipop()

# Page 6: Faceted bar
plot_density_facet()

dev.off()
cat("Done. Output:", out_pdf, "\n")
