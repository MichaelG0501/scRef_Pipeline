####################
# Analysis registry:
#   Status: terminal
#   Script: analysis/metaprograms/centred_method_comparison_figures.R
#   Description:
#     Generate side-by-side comparison PDFs for the historical uncentred GeneNMF
#     workflow and the centred GeneNMF workflow.
#   Methodology: analysis/methodology/metaprograms/centred_method_comparison_figures_methodology.md
#   Inputs:
#     - /rds/general/project/tumourheterogeneity1/live/temp_save/geneNMF_metaprograms_nMP_19.rds
#     - /rds/general/project/tumourheterogeneity1/live/temp_save/UCell_nMP19_filtered.rds
#     - ref_outs/geneNMF_outs.rds
#     - ref_outs/meta_full_epi.rds (authoritative orig.ident metadata)
#     - ref_outs/by_samples/*/*_epi_f.rds
#     - ref_outs/Metaprogrammes_Results/centred/geneNMF_metaprograms_nMP_*.rds
#     - ref_outs/Metaprogrammes_Results/{centred/,}mp_refinement/intermediate/*.rds
#   Outputs:
#     - ref_outs/Metaprogrammes_Results/centred_comparison/figures/Auto_centred_vs_uncentred_metaprogramme_comparison_all.pdf
#     - ref_outs/Metaprogrammes_Results/centred_comparison/intermediate/
#     - ref_outs/Metaprogrammes_Results/centred_comparison/tables/
#     - ref_outs/Metaprogrammes_Results/centred_comparison/logs/
#   Cache/replot behavior:
#     Run with first CLI argument prepare_uncentred in the gnmf env to create
#     missing uncentred nMP caches. Run with render in dmtcp to draw figures.
#     SCREF_FORCE_REBUILD=TRUE rebuilds comparison caches.
#   Run command:
#     eval "$(~/miniforge3/bin/conda shell.bash hook)"
#     source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#     Rscript analysis/metaprograms/centred_method_comparison_figures.R render
#   Conda env: dmtcp for render; gnmf for prepare_uncentred
####################

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[[1]] else "render"

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
ref_dir <- file.path(project_dir, "ref_outs")
outdir <- file.path(ref_dir, "Metaprogrammes_Results", "centred_comparison")
for (subdir in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
}
summary_dir <- file.path(project_dir, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

force_rebuild <- Sys.getenv("SCREF_FORCE_REBUILD", "FALSE") == "TRUE"

uncentred_nmp19_source <- "/rds/general/project/tumourheterogeneity1/live/temp_save/geneNMF_metaprograms_nMP_19.rds"
uncentred_ucell_source <- "/rds/general/project/tumourheterogeneity1/live/temp_save/UCell_nMP19_filtered.rds"

cache_uncentred_mp <- function(k) {
  file.path(outdir, "intermediate", paste0("Auto_uncentred_geneNMF_metaprograms_nMP_", k, ".rds"))
}

if (identical(mode, "prepare_uncentred")) {
  suppressPackageStartupMessages({
    library(GeneNMF)
  })

  if (!file.exists(file.path(ref_dir, "geneNMF_outs.rds"))) {
    stop("Missing ref_outs/geneNMF_outs.rds")
  }
  if (!file.exists(uncentred_nmp19_source)) {
    stop("Missing historical uncentred nMP19 source: ", uncentred_nmp19_source)
  }

  historical_uncentred_mp <- readRDS(uncentred_nmp19_source)
  for (k in 8:30) {
    out_file <- cache_uncentred_mp(k)
    if (file.exists(out_file) && !force_rebuild) {
      message("Using cached uncentred nMP ", k, ": ", out_file)
      next
    }
    message("Caching historical uncentred tree/similarity object for nMP = ", k)
    saveRDS(historical_uncentred_mp, out_file)
  }
  quit(save = "no")
}

if (!identical(mode, "render")) {
  stop("Unknown mode: ", mode)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ComplexHeatmap)
  library(circlize)
  library(cluster)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(patchwork)
  library(RColorBrewer)
  library(tidyr)
  library(viridis)
})

setwd(ref_dir)

required_files <- c(
  uncentred_nmp19_source,
  uncentred_ucell_source,
  file.path(ref_dir, "Metaprogrammes_Results", "centred", "geneNMF_metaprograms_nMP_19.rds"),
  file.path(ref_dir, "Metaprogrammes_Results", "centred", "optimal_nMP.rds"),
  file.path(ref_dir, "Metaprogrammes_Results", "mp_refinement", "intermediate", "merged_refined_mp_genes.rds"),
  file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_mp_genes.rds"),
  file.path(ref_dir, "Metaprogrammes_Results", "mp_refinement", "intermediate", "merged_refined_ucell_scores.rds"),
  file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_ucell_scores.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

missing_rank_caches <- cache_uncentred_mp(8:30)[!file.exists(cache_uncentred_mp(8:30))]
if (length(missing_rank_caches) > 0) {
  stop(
    "Missing uncentred rank cache(s). Run in gnmf first: Rscript ",
    "analysis/metaprograms/centred_method_comparison_figures.R prepare_uncentred"
  )
}

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

mp_desc_map <- c(
  "MP1" = "G2M Cell Cycle",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8" = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7" = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi."
)

submp_desc_map <- c(
  "MP7+" = "Fanconi/HR repair progenitor",
  "MP7h" = "Replication-stress dormant epithelial",
  "MP7j" = "DNA damage response",
  "MP7r" = "Stem-like glandular duct progenitor",
  "MP7v" = "Mucous secretory progenitor",
  "MP10+" = "Metabolic columnar epithelium",
  "MP10e" = "Inflammatory mucous-secretory columnar epithelium",
  "MP8+" = "Intestinal metaplasia",
  "MP8b" = "Quiescent glandular-metabolic progenitor",
  "MP8c" = "NF-kB inflammatory cycling glandular progenitor",
  "MP8e" = "Cycling intestinal-columnar progenitor",
  "MP12a" = "Enteroendocrine-primed progenitor",
  "MP12b" = "Enteroendocrine differentiation",
  "MP12c" = "Cycling glandular-intestinal progenitor",
  "MP2+" = "MYC proliferation",
  "MP2v" = "EMT-V cycling invasive progenitor",
  "MP15a" = "T/NK-like epithelial immune mimicry",
  "MP15b" = "Type I IFN-activated EMT-primed epithelial",
  "MP15c" = "Type II IFN / NF-kB peak inflammatory epithelial",
  "MP5+" = "Epithelial IFN Resp.",
  "MP16+" = "Secretory Diff. (Gastric)"
)

display_uncentred_label <- function(x) {
  vapply(x, function(xi) {
    if (xi %in% names(submp_desc_map)) return(paste0(xi, ": ", submp_desc_map[[xi]]))
    parent <- parent_id(xi)
    if (parent %in% names(mp_desc_map)) return(paste0(xi, ": ", mp_desc_map[[parent]]))
    xi
  }, character(1), USE.NAMES = FALSE)
}

jaccard_one <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

jaccard_matrix <- function(row_genes, col_genes) {
  mat <- outer(
    names(row_genes),
    names(col_genes),
    Vectorize(function(r, c) jaccard_one(row_genes[[r]], col_genes[[c]]))
  )
  rownames(mat) <- names(row_genes)
  colnames(mat) <- names(col_genes)
  mat
}

transfer_labels <- function(source_names, reference_names, jac_mat,
                            threshold = 0.25,
                            reference_label_fun = display_uncentred_label) {
  best_col <- max.col(jac_mat, ties.method = "first")
  best_reference <- colnames(jac_mat)[best_col]
  best_jaccard <- jac_mat[cbind(seq_len(nrow(jac_mat)), best_col)]
  reference_labels <- reference_label_fun(best_reference)
  assigned <- ifelse(
    is.finite(best_jaccard) & best_jaccard >= threshold,
    paste0(source_names, ": ", sub("^[^:]+: ?", "", reference_labels)),
    source_names
  )
  data.frame(
    centred_mp = source_names,
    matched_uncentred_mp = best_reference,
    best_jaccard = best_jaccard,
    label_threshold = threshold,
    adapted_label = assigned,
    stringsAsFactors = FALSE
  )
}

order_centered_by_uncentred <- function(jac_mat, uncentred_order) {
  best_col <- max.col(jac_mat, ties.method = "first")
  best_uncentred <- colnames(jac_mat)[best_col]
  ordered <- unlist(lapply(uncentred_order, function(ref) {
    candidates <- rownames(jac_mat)[best_uncentred == ref]
    if (length(candidates) == 0) return(character(0))
    candidates[order(jac_mat[candidates, ref], decreasing = TRUE)]
  }), use.names = FALSE)
  c(ordered, setdiff(rownames(jac_mat), ordered))
}

orient_cross_matrix <- function(mat, centered_order, uncentred_order,
                                centered_labels, uncentred_labels) {
  out <- t(mat[centered_order, uncentred_order, drop = FALSE])
  rownames(out) <- uncentred_labels[rownames(out)]
  colnames(out) <- centered_labels[colnames(out)]
  out
}

filter_initial_genes <- function(mp_obj) {
  genes <- mp_obj$metaprograms.genes
  metrics <- mp_obj$metaprograms.metrics
  keep <- rownames(metrics)[metrics$silhouette >= 0]
  genes[names(genes) %in% keep]
}

tree_ordered_initial_mps <- function(mp_obj, keep_names) {
  ordered_clusters <- mp_obj$programs.clusters[mp_obj$programs.tree$order]
  ordered_clusters <- ordered_clusters[!is.na(ordered_clusters)]
  ordered_mps <- paste0("MP", unique(ordered_clusters))
  ordered_mps[ordered_mps %in% keep_names]
}

infer_samples <- function(cells) {
  ####################
  # Sample identity must come from metadata rather than barcode parsing.
  meta <- readRDS(file.path(ref_dir, "meta_full_epi.rds"))
  if (inherits(meta, "Seurat")) meta <- meta@meta.data
  if (!is.data.frame(meta) || !"orig.ident" %in% colnames(meta)) {
    stop("meta_full_epi.rds must be a cell-indexed data frame containing orig.ident")
  }
  out <- unname(setNames(as.character(meta$orig.ident), rownames(meta))[cells])
  if (anyNA(out)) stop("Missing orig.ident for ", sum(is.na(out)), " comparison cells")
  ####################
  out
}

compute_fisher_cor_p <- function(score_mat, sample_vec, min_cells = 10) {
  score_mat <- scale(as.matrix(score_mat))
  feature_names <- colnames(score_mat)
  samples <- unique(sample_vec)
  n_features <- length(feature_names)
  cor_array <- array(
    NA_real_,
    dim = c(n_features, n_features, length(samples)),
    dimnames = list(feature_names, feature_names, samples)
  )
  for (samp in samples) {
    idx <- which(sample_vec == samp)
    if (length(idx) < min_cells) next
    cor_array[, , samp] <- cor(score_mat[idx, , drop = FALSE], method = "spearman")
  }
  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(NA_real_, n_features, n_features, dimnames = list(feature_names, feature_names))
  p_vals <- matrix(NA_real_, n_features, n_features, dimnames = list(feature_names, feature_names))
  for (i in seq_len(n_features)) {
    for (j in seq_len(n_features)) {
      if (i == j) {
        mean_rho[i, j] <- 1
        p_vals[i, j] <- 0
        next
      }
      zs <- z_array[i, j, ]
      zs <- zs[is.finite(zs)]
      if (length(zs) < 3) next
      mean_rho[i, j] <- tanh(mean(zs))
      tt <- tryCatch(t.test(zs), error = function(e) NULL)
      p_vals[i, j] <- if (!is.null(tt)) tt$p.value else NA_real_
    }
  }
  list(mean_rho = mean_rho, p_values = p_vals, samples = samples)
}

rank_metrics <- function(results_dir = NULL, cache_fun = NULL) {
  k_vals <- 8:30
  avg_sil_widths <- numeric(length(k_vals))
  wss_vals <- numeric(length(k_vals))
  for (i in seq_along(k_vals)) {
    k <- k_vals[[i]]
    rds_path <- if (!is.null(cache_fun)) {
      cache_fun(k)
    } else {
      file.path(results_dir, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))
    }
    if (!file.exists(rds_path)) {
      avg_sil_widths[[i]] <- NA_real_
      wss_vals[[i]] <- NA_real_
      next
    }
    mp_res <- readRDS(rds_path)
    dist_mat <- as.dist(1 - mp_res$programs.similarity)
    cluster_assignments <- cutree(mp_res$programs.tree, k = k)
    sil <- silhouette(cluster_assignments, dist = dist_mat)
    avg_sil_widths[[i]] <- summary(sil)$avg.width
    wss_k <- 0
    dist_m <- as.matrix(dist_mat)
    for (clust_id in unique(cluster_assignments)) {
      idx <- which(cluster_assignments == clust_id)
      if (length(idx) > 1) {
        cluster_dist <- dist_m[idx, idx]
        wss_k <- wss_k + sum(cluster_dist^2) / (2 * length(idx))
      }
    }
    wss_vals[[i]] <- wss_k
  }
  data.frame(nMP = k_vals, Silhouette = avg_sil_widths, WSS = wss_vals)
}

find_knee <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x1 <- x_norm[[1]]
  y1 <- y_norm[[1]]
  x2 <- x_norm[[length(x_norm)]]
  y2 <- y_norm[[length(y_norm)]]
  dists <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  x[which.max(dists)]
}

plot_rank_panel <- function(df, title_prefix) {
  sil_knee <- find_knee(df$nMP, df$Silhouette)
  wss_knee <- find_knee(df$nMP, df$WSS)
  p1 <- ggplot(df, aes(x = nMP, y = Silhouette)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 1.8) +
    geom_vline(xintercept = sil_knee, linetype = "dashed", color = "red3", linewidth = 0.6) +
    annotate("text", x = sil_knee + 0.5, y = max(df$Silhouette, na.rm = TRUE),
             label = paste0("Inflection: ", sil_knee), hjust = 0, color = "red3", size = 3) +
    theme_minimal(base_size = 10) +
    scale_x_continuous(breaks = seq(8, 30, 2)) +
    labs(title = paste0(title_prefix, " silhouette"), x = "nMP", y = "Average silhouette")
  p2 <- ggplot(df, aes(x = nMP, y = WSS)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.8) +
    geom_vline(xintercept = wss_knee, linetype = "dashed", color = "red3", linewidth = 0.6) +
    annotate("text", x = wss_knee + 0.5, y = max(df$WSS, na.rm = TRUE) * 0.96,
             label = paste0("Elbow: ", wss_knee), hjust = 0, color = "red3", size = 3) +
    theme_minimal(base_size = 10) +
    scale_x_continuous(breaks = seq(8, 30, 2)) +
    labs(title = paste0(title_prefix, " WSS elbow"), x = "nMP", y = "Total WSS")
  p1 / p2
}

make_program_heatmap <- function(mp_obj, title) {
  sim_matrix <- mp_obj$programs.similarity
  mp_clusters <- mp_obj$programs.clusters
  keep_names <- names(mp_clusters)[!is.na(mp_clusters)]
  ordered_names <- mp_obj$programs.tree$labels[mp_obj$programs.tree$order]
  final_ordered_names <- ordered_names[ordered_names %in% keep_names]
  sim_matrix <- sim_matrix[final_ordered_names, final_ordered_names, drop = FALSE]
  annotation_df <- data.frame(
    Metaprogram = paste0("MP", mp_clusters[final_ordered_names]),
    Study = vapply(strsplit(sub("\\..*$", "", final_ordered_names), "_"),
                   function(x) paste(head(x, 2), collapse = "_"), character(1)),
    row.names = final_ordered_names,
    stringsAsFactors = FALSE
  )
  annotation_df$Metaprogram <- factor(annotation_df$Metaprogram, levels = unique(annotation_df$Metaprogram))
  mp_cols <- setNames(
    colorRampPalette(brewer.pal(8, "Paired"))(length(levels(annotation_df$Metaprogram))),
    levels(annotation_df$Metaprogram)
  )
  study_cols <- setNames(viridis::viridis(length(unique(annotation_df$Study)), option = "turbo"),
                         unique(annotation_df$Study))
  top_ha <- HeatmapAnnotation(
    df = annotation_df[, c("Metaprogram", "Study"), drop = FALSE],
    col = list(Metaprogram = mp_cols, Study = study_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2, "mm")
  )
  left_ha <- rowAnnotation(
    df = annotation_df[, c("Metaprogram", "Study"), drop = FALSE],
    col = list(Metaprogram = mp_cols, Study = study_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2, "mm")
  )
  col_fun <- colorRamp2(
    c(0.00, 0.12, 0.22, 0.70, 1.00),
    c("#FFFFFF", "#F6E8A6", "#E76F51", "#5E2A84", "#000000")
  )
  Heatmap(
    sim_matrix,
    name = "Similarity",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = annotation_df$Metaprogram,
    column_split = annotation_df$Metaprogram,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    rect_gp = gpar(col = NA),
    border = FALSE,
    row_gap = unit(0.4, "mm"),
    column_gap = unit(0.4, "mm"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    top_annotation = top_ha,
    left_annotation = left_ha,
    use_raster = TRUE,
    raster_quality = 3,
    width = unit(5.6, "inch"),
    height = unit(5.6, "inch")
  )
}

run_3ca_enrichment <- function(mp_genes, ordered_mps, method_label) {
  cache_file <- file.path(outdir, "intermediate", paste0("Auto_", method_label, "_initial_nMP19_3CA_enrichment.rds"))
  if (file.exists(cache_file) && !force_rebuild) {
    return(readRDS(cache_file))
  }
  mp_ref <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv")
  mp_ref <- as.list(mp_ref)
  mp_term2gene <- data.frame(
    term = rep(names(mp_ref), lengths(mp_ref)),
    gene = unlist(mp_ref),
    row.names = NULL
  )
  mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
  mp_term2name <- data.frame(term = unique(mp_term2gene$term), name = unique(mp_term2gene$term))
  rows <- list()
  for (mp in names(mp_genes)) {
    er <- enricher(
      gene = unique(mp_genes[[mp]]),
      TERM2GENE = mp_term2gene,
      TERM2NAME = mp_term2name,
      qvalueCutoff = 0.05
    )
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) next
    r <- r[r$p.adjust < 0.05 & r$p.adjust > 0, , drop = FALSE]
    if (nrow(r) == 0) next
    rows[[mp]] <- data.frame(
      Program = mp,
      Term = r$Description,
      padj = r$p.adjust,
      Overlap = r$GeneRatio,
      stringsAsFactors = FALSE
    )
  }
  df <- dplyr::bind_rows(rows)
  result <- list(df = df, ordered_mps = ordered_mps)
  saveRDS(result, cache_file)
  result
}

make_3ca_heatmap <- function(enrich_res, title, col_labels = NULL) {
  df <- enrich_res$df
  ordered_mps <- enrich_res$ordered_mps
  if (nrow(df) == 0) {
    mat <- matrix(0, nrow = 1, ncol = length(ordered_mps), dimnames = list("No significant terms", ordered_mps))
    return(Heatmap(mat, name = "-log10 padj", column_title = title))
  }
  terms_use <- df |>
    dplyr::group_by(Program) |>
    dplyr::arrange(padj, .by_group = TRUE) |>
    dplyr::slice_head(n = 8) |>
    dplyr::ungroup() |>
    dplyr::distinct(Term) |>
    dplyr::pull(Term)
  if (length(terms_use) > 80) {
    terms_use <- df |>
      dplyr::filter(Term %in% terms_use) |>
      dplyr::group_by(Term) |>
      dplyr::summarize(min_p = min(padj), .groups = "drop") |>
      dplyr::arrange(min_p) |>
      dplyr::slice_head(n = 80) |>
      dplyr::pull(Term)
  }
  ordered_mps <- ordered_mps[ordered_mps %in% unique(df$Program)]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  final_df <- full_grid |>
    dplyr::left_join(df, by = c("Term", "Program")) |>
    dplyr::mutate(score = tidyr::replace_na(pmin(-log10(padj), 7), 0))
  mat <- final_df |>
    dplyr::select(Term, Program, score) |>
    tidyr::pivot_wider(names_from = Program, values_from = score) |>
    as.data.frame()
  rownames(mat) <- mat$Term
  mat <- as.matrix(mat[, setdiff(colnames(mat), "Term"), drop = FALSE])
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), dimnames = dimnames(mat))
  best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
  mat <- mat[order(match(best_mp, colnames(mat)), -rowSums(mat)), , drop = FALSE]
  if (!is.null(col_labels)) {
    colnames(mat) <- col_labels[colnames(mat)]
  }
  Heatmap(
    mat,
    name = "-log10 padj",
    col = colorRamp2(c(0, 3.5, 7), c("#FFFFFF", "#FB6A4A", "#67000D")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_rot = 45,
    column_title = title,
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    rect_gp = gpar(col = "grey90", lwd = 0.2),
    width = unit(5.2, "inch"),
    height = unit(7.6, "inch")
  )
}

make_matrix_heatmap <- function(mat, title, col_fun, cluster_rows = FALSE,
                                cluster_columns = FALSE, show_numbers = TRUE,
                                width = 5.8, height = 5.8,
                                row_font = 8, col_font = 8,
                                legend_title = title,
                                legend_height_mm = 40) {
  display_mat <- mat
  Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_names_gp = gpar(fontsize = row_font, fontface = "bold"),
    column_names_gp = gpar(fontsize = col_font, fontface = "bold"),
    column_names_rot = 45,
    row_names_side = "right",
    column_names_side = "bottom",
    column_title = title,
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    width = unit(width, "inch"),
    height = unit(height, "inch"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(legend_height_mm, "mm")
    ),
    cell_fun = if (show_numbers) {
      function(j, i, x, y, w, h, fill) {
        v <- display_mat[i, j]
        if (is.finite(v)) grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 6))
      }
    } else {
      NULL
    }
  )
}

draw_two_heatmaps <- function(left_ht, right_ht, page_title) {
  grid.newpage()
  grid.text(page_title, x = unit(0.5, "npc"), y = unit(0.985, "npc"),
            gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 1), "null"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(left_ht, newpage = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(right_ht, newpage = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right")
  popViewport(2)
}

draw_two_ggplots <- function(left_plot, right_plot, page_title) {
  grid.newpage()
  grid.text(page_title, x = unit(0.5, "npc"), y = unit(0.985, "npc"),
            gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 1), "null"))))
  print(left_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(right_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  popViewport()
}

score_centered_initial_ucell <- function(centered_genes, target_cells) {
  cache_file <- file.path(outdir, "intermediate", "Auto_centred_nMP19_UCell_filtered.rds")
  if (file.exists(cache_file) && !force_rebuild) {
    scores <- readRDS(cache_file)
    if (all(target_cells %in% rownames(scores)) && all(names(centered_genes) %in% colnames(scores))) {
      return(scores[target_cells, names(centered_genes), drop = FALSE])
    }
  }
  sample_dirs <- list.dirs("by_samples", full.names = FALSE, recursive = FALSE)
  sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]
  all_scores <- vector("list", length(sample_dirs))
  names(all_scores) <- sample_dirs
  for (sample in sample_dirs) {
    sample_cells <- target_cells[startsWith(target_cells, paste0(sample, "_"))]
    if (length(sample_cells) == 0) next
    rds_file <- file.path("by_samples", sample, paste0(sample, "_epi_f.rds"))
    if (!file.exists(rds_file)) next
    message("Scoring centred initial MPs for ", sample, " (", length(sample_cells), " cells)")
    obj <- readRDS(rds_file)
    keep_cells <- intersect(sample_cells, Cells(obj))
    if (length(keep_cells) == 0) next
    obj <- subset(obj, cells = keep_cells)
    obj <- AddModuleScore_UCell(obj, features = centered_genes, ncores = 1, name = "")
    all_scores[[sample]] <- obj@meta.data[keep_cells, names(centered_genes), drop = FALSE]
    rm(obj)
    invisible(gc())
  }
  all_scores <- all_scores[!vapply(all_scores, is.null, logical(1))]
  scores <- do.call(rbind, unname(all_scores))
  scores <- as.data.frame(scores, check.names = FALSE)
  missing_cells <- setdiff(target_cells, rownames(scores))
  if (length(missing_cells) > 0) {
    stop("Failed to compute centred UCell scores for ", length(missing_cells), " target cells.")
  }
  scores <- scores[target_cells, names(centered_genes), drop = FALSE]
  saveRDS(scores, cache_file)
  scores
}

message("Loading GeneNMF and refinement inputs.")
unc_mp <- readRDS(uncentred_nmp19_source)
cent_mp <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "centred", "geneNMF_metaprograms_nMP_19.rds"))
unc_genes <- filter_initial_genes(unc_mp)
cent_genes <- filter_initial_genes(cent_mp)
unc_order <- tree_ordered_initial_mps(unc_mp, names(unc_genes))
cent_order <- tree_ordered_initial_mps(cent_mp, names(cent_genes))
unc_genes <- unc_genes[unc_order]
cent_genes <- cent_genes[cent_order]

initial_jaccard <- jaccard_matrix(cent_genes, unc_genes)
initial_label_transfer <- transfer_labels(rownames(initial_jaccard), colnames(initial_jaccard), initial_jaccard)
write.csv(initial_label_transfer, file.path(outdir, "tables", "Auto_initial_centred_label_transfer.csv"), row.names = FALSE)
initial_cent_labels <- setNames(initial_label_transfer$adapted_label, initial_label_transfer$centred_mp)
initial_unc_labels <- setNames(display_uncentred_label(names(unc_genes)), names(unc_genes))
initial_cent_cross_order <- order_centered_by_uncentred(initial_jaccard, unc_order)
initial_jaccard_cross <- orient_cross_matrix(
  initial_jaccard,
  initial_cent_cross_order,
  unc_order,
  initial_cent_labels,
  initial_unc_labels
)
write.csv(initial_jaccard_cross, file.path(outdir, "tables", "Auto_initial_centred_vs_uncentred_jaccard.csv"))

unc_ucell <- readRDS(uncentred_ucell_source)
unc_ucell <- as.data.frame(unc_ucell, check.names = FALSE)
common_initial_cells <- rownames(unc_ucell)
cent_ucell <- score_centered_initial_ucell(cent_genes, common_initial_cells)
sample_vec_initial <- infer_samples(common_initial_cells)

initial_expr_cor_raw <- cor(
  as.matrix(cent_ucell[common_initial_cells, names(cent_genes), drop = FALSE]),
  as.matrix(unc_ucell[common_initial_cells, names(unc_genes), drop = FALSE]),
  method = "spearman",
  use = "pairwise.complete.obs"
)
initial_expr_cor_cross <- orient_cross_matrix(
  initial_expr_cor_raw,
  initial_cent_cross_order,
  unc_order,
  initial_cent_labels,
  initial_unc_labels
)
saveRDS(initial_expr_cor_cross, file.path(outdir, "intermediate", "Auto_initial_centred_vs_uncentred_ucell_spearman.rds"))
write.csv(initial_expr_cor_cross, file.path(outdir, "tables", "Auto_initial_centred_vs_uncentred_ucell_spearman.csv"))

state_groups <- list(
  "Classic proliferation" = c("MP2"),
  "Basal to intestinal metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress adaptive" = c("MP13", "MP12"),
  "SMG to intestinal metaplasia" = c("MP18", "MP16"),
  "Cancer-cell immune mimicry" = c("MP15")
)
final_unc_order <- unlist(state_groups, use.names = FALSE)
final_unc_order <- final_unc_order[final_unc_order %in% colnames(unc_ucell)]
unc_final_cor <- compute_fisher_cor_p(
  unc_ucell[common_initial_cells, final_unc_order, drop = FALSE],
  sample_vec_initial
)
unc_final_rho <- unc_final_cor$mean_rho
rownames(unc_final_rho) <- initial_unc_labels[rownames(unc_final_rho)]
colnames(unc_final_rho) <- initial_unc_labels[colnames(unc_final_rho)]

cent_order_for_cor <- colnames(cent_ucell)
cent_initial_cor <- compute_fisher_cor_p(
  cent_ucell[common_initial_cells, cent_order_for_cor, drop = FALSE],
  sample_vec_initial
)
cent_initial_rho <- cent_initial_cor$mean_rho
rownames(cent_initial_rho) <- initial_cent_labels[rownames(cent_initial_rho)]
colnames(cent_initial_rho) <- initial_cent_labels[colnames(cent_initial_rho)]

unc_refined_genes <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "mp_refinement", "intermediate", "merged_refined_mp_genes.rds"))
cent_refined_genes <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_mp_genes.rds"))
unc_refined_ucell <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "mp_refinement", "intermediate", "merged_refined_ucell_scores.rds"))
cent_refined_ucell <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_ucell_scores.rds"))
unc_refined_ucell <- as.data.frame(unc_refined_ucell, check.names = FALSE)
cent_refined_ucell <- as.data.frame(cent_refined_ucell, check.names = FALSE)

strict_refined_order <- c(
  "MP7j", "MP9", "MP1", "MP2+", "MP17", "MP8+", "MP10+", "MP14", "MP5+",
  "MP7r", "MP7v", "MP10e", "MP16+", "MP18",
  "MP8c", "MP15c", "MP12c", "MP2v", "MP8e", "MP12a", "MP13",
  "MP7+", "MP7h", "MP8b", "MP12b", "MP15a", "MP15b"
)
unc_refined_order <- strict_refined_order[strict_refined_order %in% names(unc_refined_genes)]
cent_refined_order <- names(cent_refined_genes)
mp_num <- as.integer(sub("^MP", "", sub("[a-z\\+]*$", "", cent_refined_order)))
cent_refined_order <- cent_refined_order[order(mp_num, cent_refined_order)]

refined_jaccard <- jaccard_matrix(cent_refined_genes[cent_refined_order], unc_refined_genes[unc_refined_order])
refined_label_transfer <- transfer_labels(rownames(refined_jaccard), colnames(refined_jaccard), refined_jaccard)
write.csv(refined_label_transfer, file.path(outdir, "tables", "Auto_refined_centred_label_transfer.csv"), row.names = FALSE)
refined_cent_labels <- setNames(refined_label_transfer$adapted_label, refined_label_transfer$centred_mp)
refined_unc_labels <- setNames(display_uncentred_label(unc_refined_order), unc_refined_order)
refined_cent_cross_order <- order_centered_by_uncentred(refined_jaccard, unc_refined_order)
refined_jaccard_cross <- orient_cross_matrix(
  refined_jaccard,
  refined_cent_cross_order,
  unc_refined_order,
  refined_cent_labels,
  refined_unc_labels
)
write.csv(refined_jaccard_cross, file.path(outdir, "tables", "Auto_refined_centred_vs_uncentred_jaccard.csv"))

common_refined_cells <- intersect(rownames(cent_refined_ucell), rownames(unc_refined_ucell))
refined_expr_cor_raw <- cor(
  as.matrix(cent_refined_ucell[common_refined_cells, cent_refined_order, drop = FALSE]),
  as.matrix(unc_refined_ucell[common_refined_cells, unc_refined_order, drop = FALSE]),
  method = "spearman",
  use = "pairwise.complete.obs"
)
refined_expr_cor_cross <- orient_cross_matrix(
  refined_expr_cor_raw,
  refined_cent_cross_order,
  unc_refined_order,
  refined_cent_labels,
  refined_unc_labels
)
saveRDS(refined_expr_cor_cross, file.path(outdir, "intermediate", "Auto_refined_centred_vs_uncentred_ucell_spearman.rds"))
write.csv(refined_expr_cor_cross, file.path(outdir, "tables", "Auto_refined_centred_vs_uncentred_ucell_spearman.csv"))

rank_unc <- rank_metrics(cache_fun = cache_uncentred_mp) |>
  mutate(method = "uncentred")
rank_cent <- rank_metrics(results_dir = file.path(ref_dir, "Metaprogrammes_Results", "centred")) |>
  mutate(method = "centred")
write.csv(bind_rows(rank_unc, rank_cent),
          file.path(outdir, "tables", "Auto_rank_selection_metrics_uncentred_vs_centred.csv"),
          row.names = FALSE)

unc_3ca <- run_3ca_enrichment(unc_genes, unc_order, "uncentred")
cent_3ca <- run_3ca_enrichment(cent_genes, cent_order, "centred")

unc_refined_cor <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "mp_refinement", "intermediate", "merged_refined_mp_correlation_matrices.rds"))
cent_refined_cor <- readRDS(file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_mp_correlation_matrices.rds"))
unc_refined_rho <- unc_refined_cor$mean_rho[unc_refined_order, unc_refined_order, drop = FALSE]
cent_refined_order <- cent_refined_order[cent_refined_order %in% rownames(cent_refined_cor$mean_rho)]
cent_refined_rho <- cent_refined_cor$mean_rho[cent_refined_order, cent_refined_order, drop = FALSE]
rownames(unc_refined_rho) <- refined_unc_labels[rownames(unc_refined_rho)]
colnames(unc_refined_rho) <- refined_unc_labels[colnames(unc_refined_rho)]
rownames(cent_refined_rho) <- refined_cent_labels[rownames(cent_refined_rho)]
colnames(cent_refined_rho) <- refined_cent_labels[colnames(cent_refined_rho)]

col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))
col_cross_cor <- colorRamp2(c(-0.2, 0.35, 0.9), c("blue", "white", "red"))
col_jacc <- colorRamp2(c(0, 0.45, 0.9), c("#FFFFFF", "#FB6A4A", "#67000D"))

pdf_file <- file.path(outdir, "figures", "Auto_centred_vs_uncentred_metaprogramme_comparison_all.pdf")
pdf(pdf_file, width = 28, height = 16, useDingbats = FALSE, onefile = TRUE)

message("Rendering page 1.")
draw_two_heatmaps(
  make_program_heatmap(unc_mp, "Uncentred nMP19 NMF programmes"),
  make_program_heatmap(cent_mp, "Centred nMP19 NMF programmes"),
  "1. Unsupervised clustered NMF programme similarity heatmaps"
)

message("Rendering page 2.")
draw_two_ggplots(
  plot_rank_panel(rank_unc, "Uncentred"),
  plot_rank_panel(rank_cent, "Centred"),
  "2. Optimal nMP silhouette and WSS elbow diagnostics"
)

message("Rendering page 3.")
draw_two_heatmaps(
  make_3ca_heatmap(unc_3ca, "Uncentred initial nMP19 3CA enrichment", initial_unc_labels),
  make_3ca_heatmap(cent_3ca, "Centred initial nMP19 3CA enrichment", initial_cent_labels),
  "3. Initial nMP19 3CA enrichment heatmaps"
)

message("Rendering page 4.")
draw_two_heatmaps(
  make_matrix_heatmap(initial_jaccard_cross, "Initial uncentred vs centred Jaccard", col_jacc,
                      width = 8.2, height = 7.4, row_font = 7, col_font = 7),
  make_matrix_heatmap(initial_expr_cor_cross, "Initial uncentred vs centred UCell Spearman", col_cross_cor,
                      width = 8.2, height = 7.4, row_font = 7, col_font = 7),
  "4. Initial MP gene overlap and expression-score correlation"
)

message("Rendering page 5.")
draw_two_heatmaps(
  make_matrix_heatmap(unc_final_rho, "Uncentred final-MP correlation, original order", col_cor,
                      cluster_rows = FALSE, cluster_columns = FALSE, width = 6.4, height = 6.4,
                      row_font = 8, col_font = 8),
  make_matrix_heatmap(cent_initial_rho, "Centred nMP19 MP correlation, unsupervised order", col_cor,
                      cluster_rows = TRUE, cluster_columns = TRUE, width = 7.4, height = 7.4,
                      row_font = 7, col_font = 7),
  "5. Initial MP expression correlation heatmaps"
)

message("Rendering page 6.")
draw_two_heatmaps(
  make_matrix_heatmap(unc_refined_rho, "Uncentred refined MP correlation, requested order", col_cor,
                      cluster_rows = FALSE, cluster_columns = FALSE, width = 8.6, height = 8.6,
                      row_font = 6.5, col_font = 6.5, legend_title = "rho", legend_height_mm = 28),
  make_matrix_heatmap(cent_refined_rho, "Centred refined MP correlation, unsupervised order", col_cor,
                      cluster_rows = TRUE, cluster_columns = TRUE, width = 7.8, height = 7.8,
                      row_font = 6.5, col_font = 6.5, legend_title = "rho", legend_height_mm = 28),
  "6. Refined MP correlation heatmaps after correlated sub-MP merge"
)

message("Rendering page 7.")
draw_two_heatmaps(
  make_matrix_heatmap(refined_jaccard_cross, "Refined uncentred vs centred Jaccard", col_jacc,
                      width = 8.4, height = 8.2, row_font = 6.5, col_font = 6.5),
  make_matrix_heatmap(refined_expr_cor_cross, "Refined uncentred vs centred UCell Spearman", col_cross_cor,
                      width = 8.4, height = 8.2, row_font = 6.5, col_font = 6.5),
  "7. Refined MP gene overlap and expression-score correlation"
)

dev.off()
message("Saved: ", pdf_file)

summary_table <- data.frame(
  output = c(
    "Auto_centred_vs_uncentred_metaprogramme_comparison_all.pdf",
    "Auto_initial_centred_vs_uncentred_jaccard.csv",
    "Auto_initial_centred_vs_uncentred_ucell_spearman.csv",
    "Auto_refined_centred_vs_uncentred_jaccard.csv",
    "Auto_refined_centred_vs_uncentred_ucell_spearman.csv"
  ),
  path = c(
    file.path(outdir, "figures", "Auto_centred_vs_uncentred_metaprogramme_comparison_all.pdf"),
    file.path(outdir, "tables", "Auto_initial_centred_vs_uncentred_jaccard.csv"),
    file.path(outdir, "tables", "Auto_initial_centred_vs_uncentred_ucell_spearman.csv"),
    file.path(outdir, "tables", "Auto_refined_centred_vs_uncentred_jaccard.csv"),
    file.path(outdir, "tables", "Auto_refined_centred_vs_uncentred_ucell_spearman.csv")
  ),
  note = c(
    "Seven-page side-by-side comparison PDF.",
    "Rows are uncentred filtered initial nMP19 MPs; columns are centred filtered initial nMP19 MPs ordered by best uncentred match.",
    "Rows are uncentred filtered initial nMP19 MPs; columns are centred filtered initial nMP19 MPs ordered by best uncentred match.",
    "Rows are uncentred merged refined MPs; columns are centred merged refined MPs ordered by best uncentred match.",
    "Rows are uncentred merged refined MPs; columns are centred merged refined MPs ordered by best uncentred match."
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table,
          file.path(summary_dir, "centred_method_comparison_figures_summary.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "logs", "Auto_centred_method_comparison_sessionInfo.txt"))

################################################################################
# Addition: Centred-method specific Enrichment, UCell reference scoring, and Excel
################################################################################

cat("\n=== Running Enrichment Analysis for Centred MPs ===\n")
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

MP_list <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv")
MP_list <- as.list(MP_list)
mp_term2gene <- data.frame(
  term = rep(names(MP_list), lengths(MP_list)),
  gene = unlist(MP_list),
  row.names = NULL
)
mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
mp_term2name <- data.frame(
  term = unique(mp_term2gene$term),
  name = unique(mp_term2gene$term)
)

individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
custom_refs <- lapply(custom_files, readRDS)
names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)

for (nm in c("Normal_Development_long", "Normal_Development_short")) {
  if (nm %in% names(custom_refs)) {
    custom_refs[[nm]]$TERM2GENE <- custom_refs[[nm]]$TERM2GENE %>% dplyr::filter(grepl("_Stomach\\.\\.", term))
    custom_refs[[nm]]$TERM2NAME <- custom_refs[[nm]]$TERM2NAME %>% dplyr::filter(grepl("_Stomach\\.\\.", term))
  }
}

heatmap_gaps_col <- 0
mp_gene_lists <- cent_refined_genes
final_col_order <- c("MP10+", "MP9+", "MP11+", "MP6+", "MP3+", "MP14", "MP17", "MP18b", "MP8+", "MP16", "MP2x", "MP18a", "MP12", "MP13+", "MP11c", "MP15", "MP8b", "MP1", "MP5", "MP2+")
merged_display_labels <- refined_cent_labels

cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {
  genes <- mp_gene_lists[[mp_name]]
  cat("Processing MP for enrichment: ", mp_name, "\n")
  res_GO <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                     ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
  res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene, 
                    TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)
  res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene, 
                    TERM2NAME = mp_term2name, qvalueCutoff = 0.05)
  res_custom_list <- lapply(names(custom_refs), function(ref_name) {
    enricher(gene = genes, TERM2GENE = custom_refs[[ref_name]]$TERM2GENE,
             TERM2NAME = custom_refs[[ref_name]]$TERM2NAME,
             pAdjustMethod = "BH", qvalueCutoff = 0.05)
  })
  names(res_custom_list) <- names(custom_refs)
  base_results <- list(rep_prog = mp_name, genes = genes, GO = res_GO, Hallmark = res_H, MPs_3CA = res_M)
  return(c(base_results, res_custom_list))
})
names(cluster_enrich) <- names(mp_gene_lists)

enrich_heatmap <- function(cluster_enrich, element, top_per_program = 8, top_n = 80, cap = 7, 
                           cols = viridis::magma(100, direction = -1), fontsize_row = 7, fontsize_col = 9) {
  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")
  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    data_source <- if(is_custom) r else r_sig
    if (nrow(data_source) == 0 && !is_custom) return(NULL)
    term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
    data.frame(Program = prog, Term = term, padj = data_source$p.adjust, Overlap = data_source$GeneRatio, stringsAsFactors = FALSE)
  })
  df <- dplyr::bind_rows(df_list)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  if (is_custom) {
    if (!element %in% names(custom_refs)) return(invisible(NULL))
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    terms_use <- df %>% dplyr::filter(padj < 0.05) %>% dplyr::arrange(Program, padj) %>% dplyr::group_by(Program) %>% dplyr::slice_head(n = top_per_program) %>% dplyr::ungroup() %>% dplyr::distinct(Term) %>% dplyr::pull(Term)
    if (length(terms_use) > top_n) {
      terms_use <- df %>% dplyr::filter(Term %in% terms_use) %>% dplyr::group_by(Term) %>% dplyr::summarise(min_p = min(padj), .groups = "drop") %>% dplyr::arrange(min_p) %>% dplyr::slice_head(n = top_n) %>% dplyr::pull(Term)
    }
  }
  
  ordered_mps <- final_col_order
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  final_df <- full_grid %>% dplyr::left_join(df, by = c("Term", "Program")) %>% dplyr::mutate(score = tidyr::replace_na(pmin(-log10(padj), cap), 0), display_text = if(element %in% c("Hallmark","GO","MPs_3CA") || is_custom) tidyr::replace_na(Overlap, "") else "")
  mat <- final_df %>% dplyr::select(Term, Program, score) %>% tidyr::pivot_wider(names_from = Program, values_from = score) %>% as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  text_mat <- final_df %>% dplyr::select(Term, Program, display_text) %>% tidyr::pivot_wider(names_from = Program, values_from = display_text) %>% as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  mp_sizes <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
  col_labels <- paste0(merged_display_labels[colnames(mat)], " (n=", mp_sizes, ")")
  
  cluster_rows_param <- FALSE; row_gaps <- NULL
  if (is_custom) {
    mat <- mat[terms_use, , drop = FALSE]
    text_mat <- text_mat[terms_use, , drop = FALSE]
  } else {
    best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_order <- order(match(best_mp, colnames(mat)), -rowSums(mat))
    mat <- mat[row_order, , drop = FALSE]
    text_mat <- text_mat[row_order, , drop = FALSE]
    groups <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_gaps <- which(groups[-length(groups)] != groups[-1])
  }
  
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  hm_name <- paste0("enrich_", element)
  ht <- ComplexHeatmap::pheatmap(mat, name = hm_name, display_numbers = text_mat, number_color = "black", fontsize_number = fontsize_row * 1.1, labels_col = col_labels, color = cols, breaks = breaks, cluster_rows = cluster_rows_param, cluster_cols = FALSE, gaps_row = row_gaps, border_color = NA, show_colnames = TRUE, angle_col = "45", fontsize_row = fontsize_row, fontsize_col = fontsize_col, main = paste0(element, " Enrichment (-log10 padj)"))
  ht@column_names_param$rot <- 35
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 35, 2, 2), "mm"))
  
  target_mps <- c("MP17", "MP18a", "MP11c", "MP8b")
  num_slices <- if (is.null(row_gaps)) 1 else length(row_gaps) + 1
  for (tmp in target_mps) {
    idx <- match(tmp, colnames(mat))
    if (!is.na(idx) && idx > 1) {
      for (i in seq_len(num_slices)) {
        ComplexHeatmap::decorate_heatmap_body(hm_name, row_slice = i, {
          x_pos <- (idx - 1) / ncol(mat)
          grid::grid.lines(c(x_pos, x_pos), c(0, 1), gp = grid::gpar(lty = 2, lwd = 1.5, col = "black"))
        })
      }
    }
  }
  return(invisible(mat))
}

cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

cat("Saving enrichment heatmaps to figures/...\n")
pdf(file.path(outdir, "figures", "Auto_centred_refined_mp_enrichment_anno.pdf"), width = 14, height = 12)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "GO",       top_per_program = 6, top_n = 60, cols = cols_palette)
enrich_heatmap(cluster_enrich, "MPs_3CA",        top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

cat("\n=== UCell Scoring in External References ===\n")
external_cache <- file.path(ref_dir, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_mp_external_ucell_scores.rds")
if (file.exists(external_cache)) {
  external_summary <- readRDS(external_cache)
} else {
  external_summary <- data.frame()
  warning("External UCell cache not found. UCell heatmaps will be blank.")
}

ucell_heatmap <- function(external_summary, element, cols = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100), fontsize_row = 7, fontsize_col = 9) {
  if (!element %in% names(custom_refs)) { return(invisible(NULL)) }
  terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  df <- if (!is.null(external_summary) && nrow(external_summary) > 0) external_summary %>% dplyr::filter(term %in% terms_use) else data.frame()
  mat <- matrix(NA_real_, nrow = length(terms_use), ncol = length(final_col_order), dimnames = list(terms_use, final_col_order))
  
  if (nrow(df) > 0) {
    mat_df <- df %>% dplyr::select(term, dplyr::all_of(intersect(final_col_order, colnames(df)))) %>% as.data.frame()
    rownames(mat_df) <- mat_df$term; mat_df$term <- NULL
    common_cols <- intersect(final_col_order, colnames(mat_df)); common_rows <- intersect(terms_use, rownames(mat_df))
    if (length(common_rows) > 0 && length(common_cols) > 0) mat[common_rows, common_cols] <- as.matrix(mat_df[common_rows, common_cols])
    text_mat <- matrix(ifelse(!is.na(mat) & mat > 0, sprintf("%.2f", mat), ""), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    main_title <- paste0(element, " Mean UCell in Reference Cells")
  } else {
    mat[] <- NA; text_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)); main_title <- paste0(element, " Mean UCell in Reference Cells\n(No expression cells available)")
  }
  
  mp_sizes <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
  col_labels <- paste0(merged_display_labels[colnames(mat)], " (n=", mp_sizes, ")")
  
  cap <- quantile(mat[mat>0], 0.98, na.rm=TRUE)
  if (is.na(cap) || cap <= 0) cap <- max(mat, na.rm=TRUE)
  if (is.na(cap) || cap <= 0) cap <- 1
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  
  hm_name <- paste0("ucell_", element)
  ht <- ComplexHeatmap::pheatmap(mat, name = hm_name, display_numbers = text_mat, number_color = "black", fontsize_number = fontsize_row * 1.1, labels_col = col_labels, color = cols, na_col = "#F0F0F0", breaks = breaks, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, show_colnames = TRUE, angle_col = "45", fontsize_row = fontsize_row, fontsize_col = fontsize_col, main = main_title)
  ht@column_names_param$rot <- 35
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 35, 2, 2), "mm"))
  
  target_mps <- c("MP17", "MP18a", "MP11c", "MP8b")
  for (tmp in target_mps) {
    idx <- match(tmp, colnames(mat))
    if (!is.na(idx) && idx > 1) {
      ComplexHeatmap::decorate_heatmap_body(hm_name, {
        x_pos <- (idx - 1) / ncol(mat)
        grid::grid.lines(c(x_pos, x_pos), c(0, 1), gp = grid::gpar(lty = 2, lwd = 1.5, col = "black"))
      })
    }
  }
  return(invisible(mat))
}

cat("Saving UCell heatmaps to figures/...\n")
pdf(file.path(outdir, "figures", "Auto_centred_refined_mp_ucell_reference_anno.pdf"), width = 14, height = 12)
elements_to_plot <- c("Early_Embryogenesis", "Normal_Development_long", "Normal_Development_short", "Organogenesis_major", "Organogenesis_sub", "Adult_Epithelium", "Barretts_Oesophagus")
for (element in elements_to_plot) { ucell_heatmap(external_summary, element, cols = cols_palette) }
dev.off()

cat("\n=== Generating Excel MP Gene Summary ===\n")
library(openxlsx)

get_desc_cent <- function(mp_name) {
  desc <- refined_cent_labels[mp_name]
  if (is.na(desc) || is.null(desc) || desc == "") {
    return("unknown")
  }
  return(desc)
}

all_ordered_excel <- final_col_order
page1_order <- all_ordered_excel

seen_parents <- character(0); page2_order <- character(0)
for (feat in all_ordered_excel) {
  p <- parent_id(feat)
  if (!p %in% seen_parents) {
    if (length(seen_parents) > 0) page2_order <- c(page2_order, "GAP")
    seen_parents <- c(seen_parents, p)
  }
  page2_order <- c(page2_order, feat)
}

build_mp_matrix_cent <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)
  get_genes <- function(mp) {
    if (mp == "GAP") return(character(0))
    res <- mp_gene_lists[[mp]]
    if (is.null(res)) { nm <- gsub("\\+", "", mp); if (!is.null(cent_mp$metaprograms.genes[[nm]])) res <- cent_mp$metaprograms.genes[[nm]] }
    return(res)
  }
  max_g <- max(sapply(mp_names_vec, function(x) length(get_genes(x))))
  n_mp <- length(mp_names_vec)
  mat <- matrix(NA_character_, nrow = max_g + 2, ncol = n_mp)
  for (i in seq_along(mp_names_vec)) {
    mp <- mp_names_vec[i]
    if (mp == "GAP") { mat[1:2, i] <- "" } else {
      mat[1, i] <- mp; mat[2, i] <- get_desc_cent(mp)
      genes <- get_genes(mp); if (length(genes) > 0) mat[3:(length(genes)+2), i] <- genes
    }
  }
  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

df_p1 <- build_mp_matrix_cent(page1_order)
df_p2 <- build_mp_matrix_cent(page2_order)

wb <- createWorkbook()
mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
desc_style <- createStyle(fgFill = "#F2F2F2")
add_sheet <- function(wb, sheet_name, df, order_vec) {
  addWorksheet(wb, sheet_name)
  sheet_idx <- length(names(wb))
  if (!is.null(df)) {
    writeData(wb, sheet = sheet_idx, x = df, startCol = 1, startRow = 1, colNames = FALSE)
    for (i in seq_along(order_vec)) {
      if (order_vec[i] == "GAP") { setColWidths(wb, sheet_idx, cols = i, widths = 3) } else {
        setColWidths(wb, sheet_idx, cols = i, widths = 25)
        addStyle(wb, sheet = sheet_idx, mp_name_style, rows = 1, cols = i, gridExpand = TRUE)
        addStyle(wb, sheet = sheet_idx, desc_style, rows = 2, cols = i, gridExpand = TRUE)
      }
    }
  }
}

add_sheet(wb, "refined_mps", df_p1, page1_order)
add_sheet(wb, "Grouped_By_Parent", df_p2, page2_order)

output_excel <- file.path(outdir, "tables", "Auto_centred_refined_MP_genes_summary.xlsx")
saveWorkbook(wb, output_excel, overwrite = TRUE)
message("Saved: ", output_excel)
