####################
# Auto_malignant_subclone_mp_heatmap.R
#
# CNA subclone versus cancer MP/state heterogeneity.
#
# Inputs:
#   ref_outs/by_samples/<sample>/<sample>_outs.rds
#   ref_outs/by_samples/<sample>/<sample>_epi_f.rds, or <sample>_epi.rds fallback
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   optional ref_outs/Auto_final_states.rds
#   optional ref_outs/Auto_topmp_v2_noreg_states_B.rds
#
# Outputs:
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_mp_sample_pages.pdf
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_mp_cohort_summary.pdf
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_summary.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_mp_tests.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_mp_subclone_tests.csv
#   ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_state_tests.csv
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
})

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "all"
min_cells <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 40L
min_subclone_cells <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 20L
min_subclone_frac <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 0.15
max_plot_cells <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1200L
max_subclones <- 6L
min_distinct_arm_delta <- 0.03
min_subclone_silhouette <- 0.12
mp_score_limit <- 2

out_dir <- "Auto_malignant_subclone_mp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)

gene_order <- read.table(
  gene_order_path,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "end"),
  stringsAsFactors = FALSE
)
chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- gene_order %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start)

ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
meta_full_epi <- readRDS("meta_full_epi.rds")

mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15")
)

cc_mps <- c("MP1", "MP7", "MP9")
extra_state_order <- c("3CA_EMT_and_Protein_maturation")
state_level_order <- c(names(state_groups), extra_state_order, "Unresolved", "Hybrid")

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

mp_cols <- c(
  "MP1_G2M Cell Cycle" = "#B0B0B0",
  "MP7_DNA Damage Repair" = "#999999",
  "MP9_G1S Cell Cycle" = "#C0C0C0",
  "MP2_MYC-related Proliferation" = "#E41A1C",
  "MP17_Basal-like Transition" = "#4DAF4A",
  "MP14_Hypoxia Adapted Epi." = "#8DA0CB",
  "MP5_Epithelial IFN Resp." = "#66C2A5",
  "MP10_Columnar Diff." = "#A6D854",
  "MP8_Intestinal Diff." = "#FC8D62",
  "MP18_Secretory Diff. (Intest.)" = "#FF7F00",
  "MP16_Secretory Diff. (Gastric)" = "#FFD92F",
  "MP13_Hypoxic Inflam. Epi." = "#984EA3",
  "MP12_Neuro-responsive Epi" = "#E78AC3",
  "MP15_Immune Infiltration" = "#377EB8"
)

label_mp <- function(mps) {
  desc <- mp_descriptions[mps]
  desc[is.na(desc)] <- mps[is.na(desc)]
  paste0(mps, "_", desc)
}

if (file.exists("Auto_final_states.rds")) {
  state_vec <- readRDS("Auto_final_states.rds")
  state_source <- "Auto_final_states.rds"
} else if (file.exists("Auto_topmp_v2_noreg_states_B.rds")) {
  state_vec <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
  state_source <- "Auto_topmp_v2_noreg_states_B.rds"
} else {
  stop("Missing state labels. Run analysis/cell_states/states_topmpB_reg_noreg.R",
       " and, when UCell_3CA_MPs.rds is available, analysis/cell_states/states_unresolved_relabel.R first.")
}

if (!file.exists("Auto_topmp_v2_noreg_mp_adj.rds")) {
  stop("Missing Auto_topmp_v2_noreg_mp_adj.rds. Run analysis/cell_states/states_topmpB_reg_noreg.R first.")
}
mp_adj_noncc <- as.matrix(readRDS("Auto_topmp_v2_noreg_mp_adj.rds"))

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- paste0("MP", unique(ordered_clusters))
mp_tree_order <- mp_tree_order[mp_tree_order %in% retained_mps]
mp_names <- unique(c(
  mp_tree_order[mp_tree_order %in% cc_mps],
  unlist(lapply(state_groups, function(mps) mp_tree_order[mp_tree_order %in% mps]), use.names = FALSE),
  mp_tree_order
))
mp_names <- mp_names[mp_names %in% names(mp_descriptions)]
mp_labels <- setNames(label_mp(mp_names), mp_names)

infer_study <- function(sample_id) {
  sub("^([^_]+_[0-9]{4}).*$", "\\1", as.character(sample_id))
}

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df %>%
    group_by(.data$.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df %>%
    group_by(.data$.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

score_sample_var <- as.character(meta_full_epi$orig.ident)
names(score_sample_var) <- rownames(meta_full_epi)
if ("study" %in% colnames(meta_full_epi)) {
  score_study_var <- as.character(meta_full_epi$study)
} else {
  score_study_var <- infer_study(score_sample_var)
}
names(score_study_var) <- rownames(meta_full_epi)

cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_common_cells <- intersect(rownames(ucell_scores), names(score_sample_var))
mp_adj_cc <- matrix(nrow = length(cc_common_cells), ncol = 0, dimnames = list(cc_common_cells, character(0)))
if (length(cc_in_ucell) > 0) {
  mp_adj_cc <- z_normalise(
    as.matrix(ucell_scores[cc_common_cells, cc_in_ucell, drop = FALSE]),
    score_sample_var,
    score_study_var
  )
}
score_common_cells <- intersect(rownames(mp_adj_noncc), cc_common_cells)
mp_score_mat <- cbind(
  mp_adj_noncc[score_common_cells, , drop = FALSE],
  mp_adj_cc[score_common_cells, , drop = FALSE]
)
mp_score_mat <- mp_score_mat[, intersect(mp_names, colnames(mp_score_mat)), drop = FALSE]
mp_names <- mp_names[mp_names %in% colnames(mp_score_mat)]
mp_labels <- setNames(label_mp(mp_names), mp_names)
topmp_mps <- setdiff(mp_names, cc_mps)
topmp_mps <- topmp_mps[topmp_mps %in% colnames(mp_adj_noncc)]
if (length(topmp_mps) == 0) stop("No non-cell-cycle MP scores available for top-MP assignment.")

sample_dirs <- basename(list.dirs("by_samples", recursive = FALSE, full.names = TRUE))
sample_dirs <- sample_dirs[file.exists(file.path("by_samples", sample_dirs, paste0(sample_dirs, "_outs.rds")))]
if (!identical(sample_arg, "all")) {
  requested <- unlist(strsplit(sample_arg, ","))
  sample_dirs <- intersect(sample_dirs, requested)
}
if (length(sample_dirs) == 0) stop("No samples found for argument: ", sample_arg)

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  cols <- colorRampPalette(base)(length(values))
  stats::setNames(cols, values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  cols <- as.character(cols)
  names(cols) <- names(cols)
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) {
    cols <- c(cols, make_palette(missing_values, palette))
  }
  cols[values]
}

assert_named_palette <- function(cols, label) {
  if (length(cols) == 0 || is.null(names(cols)) || any(is.na(names(cols))) || any(!nzchar(names(cols)))) {
    stop("Invalid annotation palette for ", label, ": ",
         paste0("[", paste(names(cols), collapse = ","), "]"), call. = FALSE)
  }
  cols
}

get_epi_path <- function(sample_id) {
  epi_f <- file.path("by_samples", sample_id, paste0(sample_id, "_epi_f.rds"))
  epi <- file.path("by_samples", sample_id, paste0(sample_id, "_epi.rds"))
  if (file.exists(epi_f)) return(epi_f)
  if (file.exists(epi)) return(epi)
  NA_character_
}

select_malignant_level1_cells <- function(epi) {
  meta <- epi@meta.data
  if ("malignancy" %in% colnames(meta)) {
    cells <- rownames(meta)[as.character(meta$malignancy) == "malignant_level_1"]
  } else if ("classification" %in% colnames(meta)) {
    cells <- rownames(meta)[as.character(meta$classification) == "cna_malignant"]
  } else {
    cells <- character(0)
  }
  intersect(cells, colnames(epi))
}

compute_sample_mp_scores <- function(cells) {
  cells <- intersect(cells, rownames(mp_score_mat))
  z <- as.matrix(mp_score_mat[cells, mp_names, drop = FALSE])
  top_use <- as.matrix(mp_adj_noncc[cells, topmp_mps, drop = FALSE])
  top_mp <- colnames(top_use)[max.col(top_use, ties.method = "first")]
  names(top_mp) <- rownames(top_use)
  list(z = z, top_mp = top_mp)
}

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

make_arm_labels <- function(go) {
  chr <- as.character(go$chromosome)
  arm <- ifelse(go$start <= centromere_pos[chr], "p", "q")
  paste0(chr, arm)
}

compute_arm_means <- function(cna_mat, cluster, arm_labels) {
  cl_levels <- sort(unique(as.character(cluster)))
  arm_levels <- unique(arm_labels)
  out <- do.call(rbind, lapply(cl_levels, function(cl) {
    cells <- names(cluster)[cluster == cl]
    vals <- tapply(seq_len(nrow(cna_mat)), factor(arm_labels, levels = arm_levels), function(ix) {
      mean(cna_mat[ix, cells, drop = FALSE], na.rm = TRUE)
    })
    vals[arm_levels]
  }))
  rownames(out) <- cl_levels
  out[is.na(out)] <- 0
  out
}

merge_small_clusters <- function(cluster, pcs) {
  cluster_names <- names(cluster)
  cluster <- as.character(cluster)
  names(cluster) <- cluster_names
  repeat {
    counts <- table(cluster)
    valid <- names(counts)[counts >= min_subclone_cells & as.numeric(counts) / length(cluster) >= min_subclone_frac]
    small <- setdiff(names(counts), valid)
    if (length(small) == 0) break
    if (length(valid) == 0) {
      cluster[] <- "1"
      break
    }
    centroids <- t(vapply(names(counts), function(cl) {
      colMeans(pcs[names(cluster)[cluster == cl], , drop = FALSE], na.rm = TRUE)
    }, numeric(ncol(pcs))))
    for (cl in small) {
      d <- vapply(valid, function(valid_cl) {
        sum((centroids[valid_cl, ] - centroids[cl, ])^2, na.rm = TRUE)
      }, numeric(1))
      if (length(d) == 0) next
      cluster[cluster == cl] <- valid[which.min(d)]
    }
  }
  cluster
}

merge_indistinct_clusters <- function(cna_mat, cluster, arm_labels) {
  cluster_names <- names(cluster)
  cluster <- as.character(cluster)
  names(cluster) <- cluster_names
  repeat {
    cl_levels <- sort(unique(cluster))
    if (length(cl_levels) <= 1) break
    arm_mean <- compute_arm_means(cna_mat, cluster, arm_labels)
    arm_call <- matrix(0L, nrow = nrow(arm_mean), ncol = ncol(arm_mean), dimnames = dimnames(arm_mean))
    arm_call[arm_mean > 0.15] <- 1L
    arm_call[arm_mean < -0.15] <- -1L

    merge_pair <- NULL
    merge_score <- Inf
    for (i in seq_len(length(cl_levels) - 1L)) {
      for (j in seq.int(i + 1L, length(cl_levels))) {
        same_calls <- identical(unname(arm_call[cl_levels[i], ]), unname(arm_call[cl_levels[j], ]))
        max_delta <- max(abs(arm_mean[cl_levels[i], ] - arm_mean[cl_levels[j], ]), na.rm = TRUE)
        rms_delta <- sqrt(mean((arm_mean[cl_levels[i], ] - arm_mean[cl_levels[j], ])^2, na.rm = TRUE))
        if (same_calls && is.finite(max_delta) && max_delta < min_distinct_arm_delta && rms_delta < 0.03) {
          if (max_delta < merge_score) {
            merge_pair <- c(cl_levels[i], cl_levels[j])
            merge_score <- max_delta
          }
        }
      }
    }
    if (is.null(merge_pair)) break
    counts <- table(cluster[cluster %in% merge_pair])
    keep <- names(which.max(counts))
    cluster[cluster %in% merge_pair] <- keep
  }
  cluster
}

infer_cna_subclones <- function(cna_mat, go) {
  n_cells <- ncol(cna_mat)
  if (n_cells < min_cells) {
    return(setNames(rep("Subclone A", n_cells), colnames(cna_mat)))
  }
  binned <- make_binned_cna(cna_mat, go, bin_size = 100L)$mat
  pc_n <- min(20L, n_cells - 1L, nrow(binned))
  pca <- prcomp(t(binned), center = TRUE, scale. = FALSE, rank. = pc_n)
  pcs <- pca$x[, seq_len(pc_n), drop = FALSE]
  hc <- hclust(dist(pcs), method = "ward.D2")
  arm_labels <- make_arm_labels(go)

  best <- list(cluster = setNames(rep("1", n_cells), colnames(cna_mat)), k = 1L, silhouette = NA_real_)
  d_pcs <- dist(pcs)
  for (k in seq.int(2L, min(max_subclones, n_cells - 1L))) {
    cluster <- cutree(hc, k = k)
    names(cluster) <- rownames(pcs)
    cluster <- merge_small_clusters(cluster, pcs)
    cluster <- merge_indistinct_clusters(cna_mat, cluster, arm_labels)
    counts <- table(cluster)
    valid <- all(counts >= min_subclone_cells & as.numeric(counts) / length(cluster) >= min_subclone_frac)
    if (!valid || length(counts) < 2) next
    sil <- if (requireNamespace("cluster", quietly = TRUE)) {
      mean(cluster::silhouette(as.integer(factor(cluster)), d_pcs)[, 3], na.rm = TRUE)
    } else {
      length(counts)
    }
    if (is.na(sil) || sil < min_subclone_silhouette) next
    if (is.na(best$silhouette) || sil > best$silhouette + 0.02 ||
        (abs(sil - best$silhouette) <= 0.02 && length(counts) < best$k)) {
      best <- list(cluster = cluster, k = length(counts), silhouette = sil)
    }
  }

  counts <- sort(table(best$cluster), decreasing = TRUE)
  message("  inferred CNA clusters: ",
          paste(paste0(names(counts), "=", as.integer(counts)), collapse = "; "),
          "; silhouette=", ifelse(is.na(best$silhouette), "NA", round(best$silhouette, 3)))
  clone_letters <- LETTERS[seq_along(counts)]
  names(clone_letters) <- names(counts)
  out <- paste0("Subclone ", clone_letters[best$cluster])
  names(out) <- names(best$cluster)
  out <- out[colnames(cna_mat)]
  attr(out, "silhouette") <- best$silhouette
  out
}

prepare_cna_matrix <- function(outs, cells) {
  cells <- intersect(cells, colnames(outs))
  common_genes <- intersect(rownames(outs), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(outs[go$gene, cells, drop = FALSE])
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  signal <- rowMeans(abs(mat), na.rm = TRUE)
  keep_signal <- signal >= as.numeric(quantile(signal, probs = 1 / 3, na.rm = TRUE))
  list(mat = mat[keep_signal, , drop = FALSE], gene_order = go[keep_signal, , drop = FALSE])
}

make_binned_cna <- function(cna_mat, go, bin_size = 100L) {
  go2 <- go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bin_levels <- unique(go2$bin)
  bins <- split(go2$.row, factor(go2$bin, levels = bin_levels))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  row_chr <- sub("_.*$", "", rownames(binned))
  list(mat = binned, chr = row_chr)
}

order_cells_by_subclone <- function(cna_mat, subclone) {
  split_cells <- split(names(subclone), factor(subclone, levels = unique(subclone)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    d <- dist(t(cna_mat[, cells, drop = FALSE]))
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

sample_plot_cells <- function(cells, subclone, max_cells = 1200L) {
  cells <- intersect(cells, names(subclone))
  if (length(cells) <= max_cells) return(cells)
  split_cells <- split(cells, factor(subclone[cells], levels = unique(subclone[cells])))
  target <- pmax(20L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

make_cna_heatmap <- function(binned, meta_plot, sample_id) {
  ord <- rownames(meta_plot)
  mat <- binned$mat[, ord, drop = FALSE]
  row_chr <- factor(binned$chr, levels = unique(binned$chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr))), levels(row_chr))
  subclone_cols <- complete_palette(make_palette(meta_plot$subclone, "Set2"), meta_plot$subclone, "Set2")
  topmp_cols <- mp_cols[names(mp_cols) %in% unique(as.character(meta_plot$top_mp_label))]
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(meta_plot$state_label))]
  missing_states <- setdiff(unique(meta_plot$state_label), names(local_state_cols))
  if (length(missing_states) > 0) {
    local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  }
  missing_topmp <- setdiff(unique(as.character(meta_plot$top_mp_label)), names(topmp_cols))
  if (length(missing_topmp) > 0) {
    topmp_cols <- c(topmp_cols, setNames(hue_pal()(length(missing_topmp)), missing_topmp))
  }
  subclone_cols <- assert_named_palette(subclone_cols, "Subclone")
  topmp_cols <- assert_named_palette(topmp_cols, "TopMP")
  local_state_cols <- assert_named_palette(local_state_cols, "State")

  top_ha <- HeatmapAnnotation(
    Subclone = meta_plot$subclone,
    TopMP = meta_plot$top_mp_label,
    State = meta_plot$state_label,
    col = list(Subclone = subclone_cols, TopMP = topmp_cols, State = local_state_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(3, "mm"),
    na_col = "grey90"
  )

  left_ha <- rowAnnotation(
    Chr = row_chr,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )

  Heatmap(
    mat,
    name = "CNA",
    col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_ha,
    left_annotation = left_ha,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    column_split = factor(meta_plot$subclone, levels = unique(meta_plot$subclone)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_title_gp = gpar(fontsize = 7),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA),
    column_title = sample_id
  )
}

make_mean_mp_heatmap <- function(mp_z, subclone) {
  mean_mat <- sapply(unique(subclone), function(cl) colMeans(mp_z[names(subclone)[subclone == cl], , drop = FALSE], na.rm = TRUE))
  if (is.null(dim(mean_mat))) mean_mat <- matrix(mean_mat, ncol = 1, dimnames = list(colnames(mp_z), unique(subclone)))
  mean_mat <- mean_mat[intersect(mp_names, rownames(mean_mat)), , drop = FALSE]
  rownames(mean_mat) <- mp_labels[rownames(mean_mat)]
  Heatmap(
    mean_mat,
    name = "Mean MP z",
    col = colorRamp2(c(-mp_score_limit, 0, mp_score_limit), c("#2166AC", "white", "#B2182B")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3)
  )
}

make_corr_heatmap <- function(mp_z, subclone) {
  mean_mat <- sapply(unique(subclone), function(cl) colMeans(mp_z[names(subclone)[subclone == cl], , drop = FALSE], na.rm = TRUE))
  if (is.null(dim(mean_mat)) || ncol(mean_mat) < 2) {
    return(textGrob("Only one CNA subclone", gp = gpar(fontsize = 12)))
  }
  cm <- suppressWarnings(cor(mean_mat, method = "spearman", use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- NA
  Heatmap(
    cm,
    name = "rho",
    col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
    cluster_rows = nrow(cm) > 2,
    cluster_columns = ncol(cm) > 2,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3)
  )
}

make_boxplot <- function(score_df, mp_test_df, sample_id) {
  set.seed(42)
  point_df <- score_df %>%
    group_by(.data$mp_label, .data$subclone) %>%
    group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), 200L)), , drop = FALSE]) %>%
    ungroup()

  label_df <- mp_test_df %>%
    mutate(sig_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "NS"
    ))
  label_df$mp_label <- factor(label_df$mp_label, levels = unique(score_df$mp_label))

  y_pos <- score_df %>%
    group_by(.data$mp_label) %>%
    summarise(y = max(.data$score_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
    left_join(label_df, by = "mp_label")

  ggplot(score_df, aes(.data$subclone, .data$score_z, fill = .data$subclone)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
    geom_jitter(data = point_df, width = 0.12, size = 0.15, alpha = 0.15) +
    geom_text(data = y_pos, aes(x = 1, y = .data$y, label = .data$sig_label), inherit.aes = FALSE, size = 2.4) +
    facet_wrap(~mp_label, scales = "free_y", ncol = 4) +
    labs(title = paste0(sample_id, ": MP scores by CNA subclone"), x = NULL, y = "MP score z") +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6)
    )
}

make_qc_boxplot <- function(meta_plot, sample_id) {
  top_subs <- meta_plot %>% count(subclone) %>% arrange(desc(n)) %>% head(2) %>% pull(subclone)
  
  if (length(top_subs) < 2) {
    return(ggplot() + theme_void() + ggtitle("Only one subclone"))
  }
  
  plot_data <- meta_plot %>%
    filter(subclone %in% top_subs) %>%
    select(subclone, nCount_RNA, nFeature_RNA, percent.mt, cc_score, cs_score) %>%
    pivot_longer(cols = c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score"), names_to = "QC_Metric", values_to = "Value") %>%
    mutate(QC_Metric = factor(QC_Metric, levels = c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score")))
  
  stats_df <- plot_data %>%
    group_by(QC_Metric) %>%
    summarise(
      p_value = tryCatch({
        sub1_vals <- Value[subclone == top_subs[1]]
        sub2_vals <- Value[subclone == top_subs[2]]
        if (length(na.omit(sub1_vals)) > 2 && length(na.omit(sub2_vals)) > 2) {
          wilcox.test(sub1_vals, sub2_vals)$p.value
        } else { NA_real_ }
      }, error = function(e) NA_real_),
      max_val = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = ifelse(!is.na(p_value) & p_value < 0.05, sprintf("p=%.2e", p_value), "NS"))
  
  ggplot(plot_data, aes(x = subclone, y = Value, fill = subclone)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8, linewidth = 0.3) +
    facet_wrap(~ QC_Metric, scales = "free_y", nrow = 2) +
    geom_text(data = stats_df, aes(x = 1.5, y = max_val * 1.05, label = label), inherit.aes = FALSE, size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 8) +
    labs(title = "QC Metrics (Top 2 Subclones)", x = NULL, y = "Value") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", strip.text = element_text(size = 7))
}

make_state_distribution_plot <- function(meta_plot, state_test_df = NULL) {
  df <- meta_plot %>%
    count(.data$subclone, .data$state_label, name = "n") %>%
    group_by(.data$subclone) %>%
    mutate(pct = 100 * .data$n / sum(.data$n)) %>%
    ungroup()
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(df$state_label))]
  missing_states <- setdiff(unique(as.character(df$state_label)), names(local_state_cols))
  if (length(missing_states) > 0) local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  subtitle <- NULL
  if (!is.null(state_test_df) && nrow(state_test_df) > 0 && !is.na(state_test_df$p_value[1])) {
    subtitle <- paste0("State x subclone p=", signif(state_test_df$p_value[1], 3),
                       ", Cramer's V=", signif(state_test_df$cramers_v[1], 3))
  }
  ggplot(df, aes(.data$subclone, .data$pct, fill = .data$state_label)) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_fill_manual(values = local_state_cols, breaks = state_level_order, drop = FALSE) +
    labs(title = "State distribution", subtitle = subtitle, x = NULL, y = "% malignant level 1 cells", fill = "State") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
}

test_states_by_subclone <- function(meta_all, sample_id) {
  tab <- table(meta_all$subclone, meta_all$state_label)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.frame(sample = sample_id, test = NA_character_, p_value = NA_real_, cramers_v = NA_real_, stringsAsFactors = FALSE))
  }
  test_name <- "chisq"
  p_val <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5, na.rm = TRUE)) {
    test_name <- "chisq_simulated"
    p_val <- tryCatch(chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA_real_)
  }
  chi <- suppressWarnings(chisq.test(tab)$statistic)
  n <- sum(tab)
  denom <- n * (min(dim(tab)) - 1)
  v <- if (is.finite(chi) && denom > 0) as.numeric(sqrt(chi / denom)) else NA_real_
  data.frame(sample = sample_id, test = test_name, p_value = p_val, cramers_v = v, stringsAsFactors = FALSE)
}

test_mps_by_subclone <- function(mp_z, subclone, sample_id) {
  df <- as.data.frame(mp_z)
  df$cell <- rownames(df)
  df$subclone <- subclone[df$cell]
  long <- df %>%
    pivot_longer(all_of(colnames(mp_z)), names_to = "mp", values_to = "score_z") %>%
    mutate(mp_label = mp_labels[.data$mp])
  test_long <- long %>% filter(.data$subclone != "Unresolved")

  tests <- test_long %>%
    group_by(.data$mp, .data$mp_label) %>%
    summarise(
      p_value = if (n_distinct(.data$subclone) >= 2) {
        ms <- tapply(.data$score_z, .data$subclone, mean, na.rm = TRUE)
        hi_cl <- names(ms)[which.max(ms)]
        lo_cl <- names(ms)[which.min(ms)]
        tryCatch(wilcox.test(.data$score_z[.data$subclone == hi_cl], .data$score_z[.data$subclone == lo_cl])$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(sample = sample_id, p_adj = p.adjust(.data$p_value, method = "BH"))

  sub_rows <- list()
  row_i <- 1L
  for (mp_id in unique(test_long$mp)) {
    mp_df <- test_long[test_long$mp == mp_id, , drop = FALSE]
    mp_label <- unique(mp_df$mp_label)[1]
    for (cl in unique(mp_df$subclone)) {
      x <- mp_df$score_z[mp_df$subclone == cl]
      y <- mp_df$score_z[mp_df$subclone != cl]
      mean_score <- mean(x, na.rm = TRUE)
      rest_mean <- mean(y, na.rm = TRUE)
      delta_mean <- mean_score - rest_mean
      p_val <- if (length(x) >= 2 && length(y) >= 2) {
        tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      sub_rows[[row_i]] <- data.frame(
        sample = sample_id,
        mp = mp_id,
        mp_label = mp_label,
        subclone = cl,
        n_cells = length(x),
        mean_score = mean_score,
        rest_mean = rest_mean,
        delta_mean = delta_mean,
        p_value = p_val,
        frac_large_effect = ifelse(delta_mean >= 0,
                                   mean(x > 1, na.rm = TRUE),
                                   mean(x < -1, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  sub_tests <- bind_rows(sub_rows)
  sub_tests$sample <- sample_id
  sub_tests <- sub_tests %>%
    ungroup() %>%
    group_by(.data$sample) %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05 & abs(.data$delta_mean) >= 0.25) %>%
    ungroup()

  list(long = long, tests = tests, sub_tests = sub_tests)
}

blank_page <- function(sample_id, reason) {
  grid.newpage()
  grid.text(sample_id, x = 0.03, y = 0.96, just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
  grid.text(reason, x = 0.03, y = 0.88, just = c("left", "top"), gp = gpar(fontsize = 11))
}

cell_rows <- list()
sample_rows <- list()
mp_tests_all <- list()
sub_tests_all <- list()
state_tests_all <- list()
qc_tests_all <- list()

sample_pdf <- file.path(out_dir, "Auto_malignant_subclone_mp_sample_pages.pdf")
pdf(sample_pdf, width = 18, height = 12, useDingbats = FALSE)

for (sample_id in sample_dirs) {
  message("Processing ", sample_id)
  epi_path <- get_epi_path(sample_id)
  if (is.na(epi_path)) {
    sample_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      status = "skipped_missing_epi_rds",
      n_malignant_level1 = NA_integer_,
      n_subclones = 0,
      stringsAsFactors = FALSE
    )
    next
  }

  outs_path <- file.path("by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  epi <- readRDS(epi_path)
  malignant_cells <- select_malignant_level1_cells(epi)
  malignant_cells <- Reduce(intersect, list(malignant_cells, rownames(mp_score_mat), names(state_vec)))
  if (length(malignant_cells) < min_cells) {
    sample_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      status = "skipped_low_malignant_level1_cells",
      n_malignant_level1 = length(malignant_cells),
      n_subclones = 0,
      stringsAsFactors = FALSE
    )
    next
  }

  outs <- readRDS(outs_path)
  cna_prepped <- prepare_cna_matrix(outs, malignant_cells)
  if (nrow(cna_prepped$mat) < 100 || ncol(cna_prepped$mat) < min_cells) {
    sample_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      status = "skipped_low_cna_signal_matrix",
      n_malignant_level1 = length(malignant_cells),
      n_subclones = 0,
      stringsAsFactors = FALSE
    )
    next
  }

  subclone <- infer_cna_subclones(cna_prepped$mat, cna_prepped$gene_order)
  subclone_silhouette <- attr(subclone, "silhouette")
  if (is.null(subclone_silhouette)) subclone_silhouette <- NA_real_

  mp <- compute_sample_mp_scores(names(subclone))
  keep_cells <- Reduce(intersect, list(names(subclone), rownames(mp$z), names(state_vec)))
  if (length(keep_cells) < min_cells) {
    sample_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      status = "skipped_low_cells_with_state_and_mp_scores",
      n_malignant_level1 = length(malignant_cells),
      n_subclones = length(unique(subclone)),
      subclone_silhouette = subclone_silhouette,
      stringsAsFactors = FALSE
    )
    next
  }
  subclone <- subclone[keep_cells]
  meta_epi <- epi@meta.data[keep_cells, , drop = FALSE]
  state_label <- as.character(state_vec[keep_cells])
  names(state_label) <- keep_cells
  top_mp_label <- mp_labels[mp$top_mp[keep_cells]]
  names(top_mp_label) <- keep_cells
  subclone_order <- paste0("Subclone ", LETTERS[seq_along(unique(subclone))])
  subclone_order <- subclone_order[subclone_order %in% unique(subclone)]
  state_order <- c(state_level_order, sort(setdiff(unique(state_label), state_level_order)))
  state_order <- state_order[state_order %in% unique(state_label)]
  topmp_order <- mp_labels[topmp_mps]
  topmp_order <- topmp_order[topmp_order %in% unique(top_mp_label)]

  meta_all <- data.frame(
    cell = keep_cells,
    sample = sample_id,
    subclone = factor(subclone[keep_cells], levels = subclone_order),
    top_mp = mp$top_mp[keep_cells],
    top_mp_label = factor(top_mp_label[keep_cells], levels = topmp_order),
    state_label = factor(state_label[keep_cells], levels = state_order),
    cna_signal = if ("cna_signal" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "cna_signal"]) else NA_real_,
    cna_cor = if ("cna_cor" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "cna_cor"]) else NA_real_,
    nCount_RNA = if ("nCount_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nCount_RNA"]) else NA_real_,
    nFeature_RNA = if ("nFeature_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nFeature_RNA"]) else NA_real_,
    percent.mt = if ("percent.mt" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "percent.mt"]) else NA_real_,
    cc_score = if ("cc_score" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "cc_score"]) else NA_real_,
    cs_score = if ("cs_score" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "cs_score"]) else NA_real_,
    stringsAsFactors = FALSE,
    row.names = keep_cells
  )

  plot_cells <- sample_plot_cells(keep_cells, subclone, max_plot_cells)
  cna_order <- order_cells_by_subclone(cna_prepped$mat[, plot_cells, drop = FALSE], factor(subclone[plot_cells], levels = subclone_order))
  meta_plot <- meta_all[cna_order, , drop = FALSE]
  binned <- make_binned_cna(cna_prepped$mat[, cna_order, drop = FALSE], cna_prepped$gene_order)
  test_res <- test_mps_by_subclone(mp$z[keep_cells, , drop = FALSE], subclone, sample_id)
  state_test <- test_states_by_subclone(meta_all, sample_id)

  cell_rows[[sample_id]] <- meta_all %>% as.data.frame() %>% select(-cell) %>% mutate(cell = rownames(meta_all), .before = 1)
  sample_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "analysed",
    n_malignant_level1 = length(malignant_cells),
    n_cells_used = length(keep_cells),
    n_subclones = length(unique(subclone)),
    subclone_silhouette = subclone_silhouette,
    state_source = state_source,
    stringsAsFactors = FALSE
  )
  mp_tests_all[[sample_id]] <- test_res$tests
  sub_tests_all[[sample_id]] <- test_res$sub_tests
  state_tests_all[[sample_id]] <- state_test
  
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score")
  qc_test_rows <- list()
  for (q in qc_metrics) {
    if (q %in% colnames(meta_all) && n_distinct(meta_all$subclone) >= 2) {
      ms <- tapply(meta_all[[q]], meta_all$subclone, mean, na.rm = TRUE)
      hi_cl <- names(ms)[which.max(ms)]
      lo_cl <- names(ms)[which.min(ms)]
      p <- tryCatch(wilcox.test(meta_all[[q]][meta_all$subclone == hi_cl], 
                                meta_all[[q]][meta_all$subclone == lo_cl])$p.value, 
                    error = function(e) NA_real_)
      qc_test_rows[[q]] <- data.frame(sample = sample_id, metric = q, p_value = p, stringsAsFactors = FALSE)
    }
  }
  qc_tests_all[[sample_id]] <- bind_rows(qc_test_rows)

  score_df <- test_res$long %>%
    mutate(subclone = factor(.data$subclone, levels = subclone_order),
           mp_label = factor(.data$mp_label, levels = mp_labels[mp_names]))

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 4,
    widths = unit(c(5.9, 2.7, 3.6, 3.6), "null"),
    heights = unit(c(5.8, 5.2), "null")
  )))

  cna_ht <- make_cna_heatmap(binned, meta_plot, sample_id)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  draw(cna_ht, newpage = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right")
  popViewport()

  mean_ht <- make_mean_mp_heatmap(mp$z[keep_cells, , drop = FALSE], subclone)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  draw(mean_ht, newpage = FALSE, heatmap_legend_side = "right")
  popViewport()

  corr_obj <- make_corr_heatmap(mp$z[keep_cells, , drop = FALSE], subclone)
  if (inherits(corr_obj, "Heatmap")) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
    draw(corr_obj, newpage = FALSE, heatmap_legend_side = "right")
    popViewport()
  } else {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
    grid.draw(corr_obj)
    popViewport()
  }

  print(make_boxplot(score_df, test_res$tests, sample_id), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  print(make_state_distribution_plot(meta_all, state_test), vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(make_qc_boxplot(meta_plot, sample_id), vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
  popViewport()
}

dev.off()

cell_df <- bind_rows(cell_rows)
sample_df <- bind_rows(sample_rows)
mp_tests_df <- bind_rows(mp_tests_all)
sub_tests_df <- bind_rows(sub_tests_all)
state_tests_df <- bind_rows(state_tests_all)
qc_tests_df <- bind_rows(qc_tests_all)
if (nrow(qc_tests_df) > 0) {
  qc_tests_df <- qc_tests_df %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
}
if (nrow(state_tests_df) > 0 && "p_value" %in% colnames(state_tests_df)) {
  state_tests_df <- state_tests_df %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05)
}

write.csv(cell_df, file.path(out_dir, "Auto_malignant_subclone_cells.csv"), row.names = FALSE)
write.csv(sample_df, file.path(out_dir, "Auto_malignant_subclone_summary.csv"), row.names = FALSE)
write.csv(mp_tests_df, file.path(out_dir, "Auto_malignant_subclone_mp_tests.csv"), row.names = FALSE)
write.csv(sub_tests_df, file.path(out_dir, "Auto_malignant_subclone_mp_subclone_tests.csv"), row.names = FALSE)
write.csv(state_tests_df, file.path(out_dir, "Auto_malignant_subclone_state_tests.csv"), row.names = FALSE)

if (nrow(sub_tests_df) > 0) {
  multi_subclone_samples <- sample_df$sample[sample_df$n_subclones >= 2]
  
  sig_counts_sample <- mp_tests_df %>%
    filter(sample %in% multi_subclone_samples) %>%
    mutate(significant = !is.na(p_adj) & p_adj < 0.05) %>%
    group_by(sample) %>%
    summarise(n_sig_mps = sum(significant, na.rm = TRUE), .groups = "drop") %>%
    mutate(category = case_when(
      n_sig_mps == 0 ~ "None",
      n_sig_mps == 1 ~ "One significant",
      TRUE ~ "More than one"
    ))

  p_counts <- sig_counts_sample %>%
    count(category, name = "n") %>%
    mutate(category = factor(category, levels = c("None", "One significant", "More than one")),
           pct = 100 * n / sum(n)) %>%
    ggplot(aes(category, pct, fill = category)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c("None" = "grey70", "One significant" = "#FDB863", "More than one" = "#B2182B")) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = "Significant MP differences per sample", x = NULL, y = "Percentage of samples") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  target_mps <- c("MP1", "MP7", "MP9", "MP2", "MP17", "MP14", "MP5", "MP10", "MP8", "MP18", "MP16", "MP13", "MP12", "MP15")
  # Filter to those present in mp_labels
  target_mps <- target_mps[target_mps %in% names(mp_labels)]
  
  mp_plot_df <- mp_tests_df %>%
    filter(sample %in% multi_subclone_samples, mp %in% target_mps) %>%
    mutate(mp_label = factor(mp_label, levels = mp_labels[target_mps]),
           val = -log10(p_adj))
  
  mp_pcts <- mp_plot_df %>%
    group_by(mp_label) %>%
    summarise(pct = 100 * mean(p_adj < 0.05, na.rm = TRUE), .groups = "drop")
  
  p_mp <- ggplot(mp_plot_df, aes(mp_label, val, fill = mp_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = mp_pcts, aes(x = mp_label, y = max(mp_plot_df$val, na.rm = TRUE) * 1.15, label = sprintf("%.1f%%", pct)), inherit.aes = FALSE, size = 3) +
    scale_y_log10(expand = expansion(mult = c(0.1, 0.3))) +
    scale_fill_manual(values = mp_cols[mp_labels[target_mps]]) +
    labs(title = NULL, x = NULL, y = "-log10(p_adj)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  # Page 2: QC Summary across the cohort
  qc_summary_rows <- list()
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score")
  
  for (samp in multi_subclone_samples) {
    samp_df <- cell_df %>% filter(sample == samp, subclone != "Unresolved")
    top_subs <- samp_df %>% count(subclone) %>% arrange(desc(n)) %>% head(2) %>% pull(subclone)
    if (length(top_subs) < 2) next
    
    row_data <- data.frame(sample = samp, clone1 = top_subs[1], clone2 = top_subs[2])
    for (q in qc_metrics) {
      if (q %in% colnames(samp_df)) {
        cl_means <- samp_df %>% group_by(subclone) %>% summarise(m = mean(.data[[q]], na.rm = TRUE), .groups = "drop")
        v1 <- max(cl_means$m, na.rm = TRUE)
        v2 <- min(cl_means$m, na.rm = TRUE)
        row_data[[paste0("X_", q)]] <- v1
        row_data[[paste0("Y_", q)]] <- v2
      }
    }
    qc_summary_rows[[length(qc_summary_rows) + 1]] <- row_data
  }
  qc_summary_df <- bind_rows(qc_summary_rows)

  plot_list_qc <- list()
  for (q in qc_metrics) {
    if (paste0("X_", q) %in% colnames(qc_summary_df)) {
      x_col <- paste0("X_", q)
      y_col <- paste0("Y_", q)
      wt <- tryCatch(wilcox.test(qc_summary_df[[x_col]], qc_summary_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
      p_val <- wt$p.value
      diff_stat <- mean(qc_summary_df[[x_col]] - qc_summary_df[[y_col]], na.rm = TRUE)
      
      p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
      subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
      
      qc_all_vals <- unlist(qc_summary_df[, c(x_col, y_col)])
      qc_lims <- quantile(qc_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
      
      p <- ggplot(qc_summary_df, aes(.data[[x_col]], .data[[y_col]])) +
        geom_point(alpha = 0.5, size = 1.2, color = "black") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        coord_cartesian(xlim = qc_lims, ylim = qc_lims) +
        labs(title = q, subtitle = subtitle, x = "Highest subclone", y = "Lowest subclone") +
        theme_classic(base_size = 9) +
        theme(plot.title = element_text(size = 8, face = "bold"),
              plot.subtitle = element_text(size = 7.5))
      plot_list_qc[[q]] <- p
    }
  }
  
  if (nrow(qc_tests_df) > 0) {
    qc_plot_df <- qc_tests_df %>% 
      filter(sample %in% multi_subclone_samples) %>%
      mutate(val = -log10(p_adj))
    qc_pcts <- qc_plot_df %>% group_by(metric) %>% summarise(pct = 100 * mean(p_adj < 0.05, na.rm = TRUE), .groups = "drop")
    p_qc_sig <- ggplot(qc_plot_df, aes(metric, val, fill = metric)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text(data = qc_pcts, aes(x = metric, y = max(qc_plot_df$val, na.rm = TRUE) * 1.15, label = sprintf("%.1f%%", pct)), inherit.aes = FALSE, size = 2.5) +
      scale_y_log10(expand = expansion(mult = c(0.1, 0.3))) +
      labs(title = "QC Significance", x = NULL, y = "-log10(p_adj)") +
      theme_classic(base_size = 9) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(size = 8, face = "bold"))
    plot_list_qc[["Sig"]] <- p_qc_sig
  }

  valid_cells <- intersect(cell_df$cell, rownames(mp_score_mat))
  cell_df_valid <- cell_df[match(valid_cells, cell_df$cell), ]
  mp_scores_valid <- mp_score_mat[valid_cells, target_mps, drop = FALSE]
  
  subclone_means <- cell_df_valid %>%
    select(cell, sample, subclone) %>%
    bind_cols(as.data.frame(mp_scores_valid)) %>%
    filter(sample %in% multi_subclone_samples, subclone != "Unresolved") %>%
    group_by(sample, subclone) %>%
    summarise(across(all_of(target_mps), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    
  pair_rows <- list()
  for (samp in unique(subclone_means$sample)) {
    samp_df <- subclone_means %>% filter(sample == samp) %>% arrange(subclone)
    if (nrow(samp_df) >= 2) {
      combos <- combn(nrow(samp_df), 2)
      for (i in 1:ncol(combos)) {
        idx1 <- combos[1, i]
        idx2 <- combos[2, i]
        row_data <- data.frame(
          sample = samp,
          clone1 = samp_df$subclone[idx1],
          clone2 = samp_df$subclone[idx2]
        )
        for (mp in target_mps) {
          val1 <- samp_df[[mp]][idx1]
          val2 <- samp_df[[mp]][idx2]
          row_data[[paste0("X_", mp)]] <- max(val1, val2, na.rm = TRUE)
          row_data[[paste0("Y_", mp)]] <- min(val1, val2, na.rm = TRUE)
        }
        pair_rows[[length(pair_rows) + 1]] <- row_data
      }
    }
  }
  pairs_df <- bind_rows(pair_rows)
  
  mp_x_cols <- paste0("X_", target_mps)
  mp_y_cols <- paste0("Y_", target_mps)
  mp_all_vals <- unlist(pairs_df[, c(mp_x_cols, mp_y_cols)])
  mp_global_lims <- quantile(mp_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
  
  plot_list_mp <- list()
  for (mp in target_mps) {
    x_col <- paste0("X_", mp)
    y_col <- paste0("Y_", mp)
    wt <- tryCatch(wilcox.test(pairs_df[[x_col]], pairs_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
    p_val <- wt$p.value
    diff_stat <- mean(pairs_df[[x_col]] - pairs_df[[y_col]], na.rm = TRUE)
    
    title_label <- mp_labels[mp]
    p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
    subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
    
    p <- ggplot(pairs_df, aes(.data[[x_col]], .data[[y_col]])) +
      geom_point(alpha = 0.5, color = mp_cols[mp_labels[mp]], size = 1.2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      labs(title = title_label, subtitle = subtitle, x = "Higher-expressed clone", y = "Lower-expressed clone") +
      coord_cartesian(xlim = mp_global_lims, ylim = mp_global_lims) +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7.5))
    plot_list_mp[[mp]] <- p
  }

  state_counts <- cell_df %>%
    filter(sample %in% multi_subclone_samples, subclone != "Unresolved") %>%
    count(sample, subclone, state_label) %>%
    group_by(sample, subclone) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  target_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia", "SMG-like Metaplasia", "Stress-adaptive", "Immune Infiltrating", "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid")
  
  pair_rows_state <- list()
  for (samp in unique(state_counts$sample)) {
    samp_df_st <- state_counts %>% filter(sample == samp)
    subs <- sort(unique(samp_df_st$subclone))
    if (length(subs) >= 2) {
      combos <- combn(length(subs), 2)
      for (i in 1:ncol(combos)) {
        sub1 <- subs[combos[1, i]]
        sub2 <- subs[combos[2, i]]
        row_data <- data.frame(sample = samp, clone1 = sub1, clone2 = sub2)
        for (st in target_states) {
          p1 <- samp_df_st$prop[samp_df_st$subclone == sub1 & samp_df_st$state_label == st]
          p2 <- samp_df_st$prop[samp_df_st$subclone == sub2 & samp_df_st$state_label == st]
          val1 <- if (length(p1) > 0) p1[1] else 0
          val2 <- if (length(p2) > 0) p2[1] else 0
          row_data[[paste0("X_", st)]] <- max(val1, val2, na.rm = TRUE)
          row_data[[paste0("Y_", st)]] <- min(val1, val2, na.rm = TRUE)
        }
        pair_rows_state[[length(pair_rows_state) + 1]] <- row_data
      }
    }
  }
  pairs_state_df <- bind_rows(pair_rows_state)

  state_x_cols <- paste0("X_", target_states)
  state_y_cols <- paste0("Y_", target_states)
  state_all_vals <- unlist(pairs_state_df[, c(state_x_cols, state_y_cols)])
  state_global_lims <- quantile(state_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)

  plot_list_state <- list()
  for (st in target_states) {
    x_col <- paste0("X_", st)
    y_col <- paste0("Y_", st)
    wt <- tryCatch(wilcox.test(pairs_state_df[[x_col]], pairs_state_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
    p_val <- wt$p.value
    diff_stat <- mean(pairs_state_df[[x_col]] - pairs_state_df[[y_col]], na.rm = TRUE)
    
    st_col <- if (st %in% names(state_cols)) state_cols[[st]] else "grey50"
    p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
    subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
    
    p <- ggplot(pairs_state_df, aes(.data[[x_col]], .data[[y_col]])) +
      geom_point(alpha = 0.5, color = st_col, size = 1.2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      labs(title = st, subtitle = subtitle, x = "Higher-abundance clone", y = "Lower-abundance clone") +
      coord_cartesian(xlim = state_global_lims, ylim = state_global_lims) +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7.5))
    plot_list_state[[st]] <- p
  }

  write.csv(sig_counts_sample, file.path(out_dir, "Auto_malignant_subclone_sig_count_summary.csv"), row.names = FALSE)
  write.csv(mp_summary_sample, file.path(out_dir, "Auto_malignant_subclone_mp_cohort_summary.csv"), row.names = FALSE)

  pdf(file.path(out_dir, "Auto_malignant_subclone_mp_cohort_summary.pdf"), width = 15, height = 9, useDingbats = FALSE)
  # Page 1
  grid.arrange(p_counts, p_mp, ncol = 2, widths = c(1, 2))
  # Page 2: QC Summary
  grid.arrange(grobs = plot_list_qc, ncol = 3)
  # Page 3
  grid.arrange(grobs = plot_list_mp, ncol = 4)
  # Page 4
  grid.arrange(grobs = plot_list_state, ncol = 4)
  dev.off()
}

message("Saved sample pages to: ", sample_pdf)
message("Saved tables and cohort summary to: ", out_dir)
