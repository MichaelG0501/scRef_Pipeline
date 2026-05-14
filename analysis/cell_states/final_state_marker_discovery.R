####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/final_state_marker_discovery.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_six_state_markers.R
#
# Rebuild a six-state epithelial embedding from the finalized states,
# then derive robust state markers using a sample-aware recurrent DGE workflow.
#
# Inputs:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_final_states.rds
#
# Outputs:
#   ref_outs/Auto_six_state_markers/Auto_six_state_umap.pdf
#   ref_outs/Auto_six_state_markers/Auto_six_state_umap_embeddings.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_cluster_state_table.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_global_marker_screen.csv.gz
#   ref_outs/Auto_six_state_markers/Auto_six_state_sample_state_eligibility.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_per_sample_dge.csv.gz
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_summary.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_ranked.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_final.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_top5_recurrence_summary.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_top5_sample_support.csv.gz
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_top5_study_support.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_heatmap_top.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_heatmap_matrix.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_heatmap.pdf
#   analysis/methodology/cell_states/final_state_marker_discovery_methodology.md
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(parallel)
  library(data.table)
  library(grid)
})

####################
# setup
####################
project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

out_dir <- "Auto_six_state_markers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1234)

####################
# constants
####################
state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

params <- list(
  n_variable_features = 3000,
  n_pcs = 30,
  cluster_resolution = 0.5,
  min_cells_feature = 20,
  min_cells_state = 20,
  min_cells_rest = 20,
  candidate_pool_per_state = 1500,
  top_markers_per_state = 5,
  mc_cores = 1L
)

####################
# helpers
####################
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

pick_logfc_col <- function(df) {
  candidates <- c("avg_log2FC", "avg_logFC")
  out <- candidates[candidates %in% colnames(df)][1]
  if (is.na(out)) stop("Could not find a logFC column in the marker table.")
  out
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

run_sample_state_markers <- function(sample_id, sample_cells, obj, state_levels, candidate_map, min_cells_state, min_cells_rest) {
  sample_obj <- obj[, sample_cells]
  sample_meta <- sample_obj@meta.data
  sample_study <- unique(as.character(sample_meta$study))
  if (length(sample_study) == 0) sample_study <- NA_character_
  sample_study <- sample_study[1]
  Idents(sample_obj) <- "final_state6"

  eligibility_rows <- vector("list", length(state_levels))
  marker_rows <- list()

  for (i in seq_along(state_levels)) {
    state_name <- state_levels[i]
    state_cells <- colnames(sample_obj)[sample_obj$final_state6 == state_name]
    other_cells <- colnames(sample_obj)[sample_obj$final_state6 != state_name]
    state_n <- length(state_cells)
    other_n <- length(other_cells)
    eligible <- state_n >= min_cells_state && other_n >= min_cells_rest

    eligibility_rows[[i]] <- data.frame(
      sample = sample_id,
      study = sample_study,
      state = state_name,
      state_cell_n = state_n,
      other_cell_n = other_n,
      eligible = eligible,
      stringsAsFactors = FALSE
    )

    if (!eligible) next

    features_use <- unique(candidate_map[[state_name]] %||% character())
    if (length(features_use) == 0) next
    features_use <- intersect(features_use, rownames(sample_obj))
    if (length(features_use) == 0) next
    other_states <- intersect(
      setdiff(state_levels, state_name),
      unique(as.character(sample_obj$final_state6))
    )
    if (length(other_states) == 0) next

    res <- tryCatch(
      FindMarkers(
        object = sample_obj,
        ident.1 = state_name,
        ident.2 = other_states,
        features = features_use,
        logfc.threshold = 0,
        min.pct = 0,
        test.use = "wilcox",
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0) next

    logfc_col <- pick_logfc_col(res)
    res$gene <- rownames(res)
    res$avg_log2FC <- res[[logfc_col]]
    res$sample <- sample_id
    res$study <- sample_study
    res$state <- state_name
    res$state_cell_n <- state_n
    res$other_cell_n <- other_n
    res$pct_state <- res$pct.1 %||% NA_real_
    res$pct_other <- res$pct.2 %||% NA_real_
    res$pct_delta <- res$pct_state - res$pct_other

    marker_rows[[state_name]] <- res %>%
      select(
        gene,
        sample,
        study,
        state,
        p_val,
        p_val_adj,
        avg_log2FC,
        pct_state,
        pct_other,
        pct_delta,
        state_cell_n,
        other_cell_n
      )
  }

  list(
    eligibility = bind_rows(eligibility_rows),
    markers = bind_rows(marker_rows)
  )
}

####################
# data loading
####################
message("Loading Seurat object and finalized state labels.")

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_labels <- readRDS("Auto_final_states.rds")

DefaultAssay(tmdata_all) <- "RNA"

common_cells <- intersect(colnames(tmdata_all), names(state_labels))
state_labels <- as.character(state_labels[common_cells])
keep_cells <- common_cells[state_labels %in% state_order]
state_labels <- state_labels[state_labels %in% state_order]

message("Building a lean six-state count matrix.")

meta_state6 <- tmdata_all@meta.data[keep_cells, c("orig.ident", "study"), drop = FALSE]
meta_state6$final_state_full <- state_labels

counts_mat <- GetAssayData(tmdata_all, assay = "RNA", slot = "counts")[, keep_cells, drop = FALSE]
gene_detect_n <- Matrix::rowSums(counts_mat > 0)
keep_features <- rownames(counts_mat)[gene_detect_n >= params$min_cells_feature]
counts_mat <- counts_mat[keep_features, , drop = FALSE]

rm(tmdata_all, state_labels, common_cells, gene_detect_n, keep_features)
invisible(gc())

tmdata_state6 <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_state6,
  assay = "RNA"
)

rm(counts_mat, meta_state6)
invisible(gc())

tmdata_state6$final_state6 <- factor(as.character(tmdata_state6$final_state_full), levels = state_order)
Idents(tmdata_state6) <- "final_state6"

state_counts <- table(tmdata_state6$final_state6)
sample_counts <- tmdata_state6@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  count(final_state6, orig.ident, study, name = "n_cells")

####################
# six-state embedding
####################
tmdata_state6_cache <- file.path(cache_dir, "tmdata_state6_embedded.rds")

if (file.exists(tmdata_state6_cache)) {
  message("Loading cached six-state embedded Seurat object.")
  tmdata_state6 <- readRDS(tmdata_state6_cache)
} else {
  message("Rebuilding six-state embedding and clustering.")

  tmdata_state6 <- NormalizeData(tmdata_state6, verbose = FALSE)
  tmdata_state6 <- FindVariableFeatures(
    tmdata_state6,
    selection.method = "vst",
    nfeatures = params$n_variable_features,
    verbose = FALSE
  )
  tmdata_state6 <- ScaleData(
    tmdata_state6,
    features = VariableFeatures(tmdata_state6),
    verbose = FALSE
  )
  tmdata_state6 <- RunPCA(
    tmdata_state6,
    features = VariableFeatures(tmdata_state6),
    npcs = params$n_pcs,
    verbose = FALSE
  )
  tmdata_state6 <- FindNeighbors(
    tmdata_state6,
    dims = seq_len(params$n_pcs),
    verbose = FALSE
  )
  tmdata_state6 <- FindClusters(
    tmdata_state6,
    resolution = params$cluster_resolution,
    verbose = FALSE
  )
  tmdata_state6 <- RunUMAP(
    tmdata_state6,
    dims = seq_len(params$n_pcs),
    verbose = FALSE
  )
  
  message("Caching six-state embedded Seurat object.")
  saveRDS(tmdata_state6, tmdata_state6_cache)
}

umap_df <- Embeddings(tmdata_state6, reduction = "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    state = as.character(tmdata_state6$final_state6[match(cell, colnames(tmdata_state6))]),
    sample = as.character(tmdata_state6$orig.ident[match(cell, colnames(tmdata_state6))]),
    study = as.character(tmdata_state6$study[match(cell, colnames(tmdata_state6))]),
    cluster = as.character(tmdata_state6$seurat_clusters[match(cell, colnames(tmdata_state6))])
  )

fwrite(
  umap_df,
  file.path(out_dir, "Auto_six_state_umap_embeddings.csv")
)

cluster_state_table <- tmdata_state6@meta.data %>%
  count(seurat_clusters, final_state6, name = "n_cells") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

fwrite(
  cluster_state_table,
  file.path(out_dir, "Auto_six_state_cluster_state_table.csv")
)

p_state <- DimPlot(
  tmdata_state6,
  reduction = "umap",
  group.by = "final_state6",
  cols = state_cols,
  pt.size = 0.20,
  raster = TRUE
) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "Finalized six-state epithelial UMAP",
    subtitle = paste0(
      ncol(tmdata_state6),
      " cells across ",
      length(unique(tmdata_state6$orig.ident)),
      " samples / ",
      length(unique(tmdata_state6$study)),
      " studies"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_cluster <- DimPlot(
  tmdata_state6,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.20,
  raster = TRUE
) +
  labs(
    title = "Unsupervised clusters on the six-state subset",
    subtitle = paste0("resolution = ", params$cluster_resolution)
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(out_dir, "Auto_six_state_umap.pdf"),
  plot = p_state + p_cluster + plot_layout(widths = c(1.05, 1)),
  width = 16,
  height = 7,
  useDingbats = FALSE
)

####################
# global marker screen
####################
global_markers_cache <- file.path(cache_dir, "global_marker_screen.rds")

if (file.exists(global_markers_cache)) {
  message("Loading cached global marker screen.")
  global_markers <- readRDS(global_markers_cache)
} else {
  message("Running pooled six-state descriptive enrichment screen.")

  expr_global <- GetAssayData(tmdata_state6, assay = "RNA", slot = "data")

  global_markers <- bind_rows(lapply(state_order, function(state_name) {
    state_cells <- colnames(tmdata_state6)[tmdata_state6$final_state6 == state_name]
    other_cells <- colnames(tmdata_state6)[tmdata_state6$final_state6 != state_name]

    mean_state <- Matrix::rowMeans(expr_global[, state_cells, drop = FALSE])
    mean_other <- Matrix::rowMeans(expr_global[, other_cells, drop = FALSE])
    pct_state <- Matrix::rowMeans(expr_global[, state_cells, drop = FALSE] > 0)
    pct_other <- Matrix::rowMeans(expr_global[, other_cells, drop = FALSE] > 0)

    data.frame(
      gene = rownames(expr_global),
      state = state_name,
      global_mean_state = mean_state,
      global_mean_other = mean_other,
      global_mean_diff = mean_state - mean_other,
      global_pct_state = pct_state,
      global_pct_other = pct_other,
      global_pct_delta = pct_state - pct_other,
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(state, desc(global_mean_diff), desc(global_pct_delta), desc(global_pct_state))
  
  message("Caching global marker screen.")
  saveRDS(global_markers, global_markers_cache)
}

fwrite(
  global_markers,
  file.path(out_dir, "Auto_six_state_global_marker_screen.csv.gz")
)

candidate_map <- global_markers %>%
  filter(global_mean_diff > 0) %>%
  group_by(state) %>%
  arrange(desc(global_mean_diff), desc(global_pct_delta), desc(global_pct_state), .by_group = TRUE) %>%
  slice_head(n = params$candidate_pool_per_state) %>%
  summarise(candidate_genes = list(unique(gene)), .groups = "drop")

candidate_map <- setNames(candidate_map$candidate_genes, candidate_map$state)
candidate_map <- candidate_map[state_order]
candidate_map[sapply(candidate_map, is.null)] <- list(character())

if (exists("expr_global")) rm(expr_global)
invisible(gc())

####################
# per-sample DGE
####################
sample_dge_cache <- file.path(cache_dir, "per_sample_dge_ranked_top5.rds")

if (file.exists(sample_dge_cache)) {
  message("Loading cached per-sample DGE results.")
  cached_sample_res <- readRDS(sample_dge_cache)
  eligibility_df <- cached_sample_res$eligibility
  per_sample_dge <- cached_sample_res$markers
} else {
  message("Running per-sample recurrent DGE validation.")

  sample_cells_map <- split(colnames(tmdata_state6), tmdata_state6$orig.ident)
  sample_ids <- names(sample_cells_map)

  sample_res <- mclapply(
    sample_ids,
    function(sample_id) {
      run_sample_state_markers(
        sample_id = sample_id,
        sample_cells = sample_cells_map[[sample_id]],
        obj = tmdata_state6,
        state_levels = state_order,
        candidate_map = candidate_map,
        min_cells_state = params$min_cells_state,
        min_cells_rest = params$min_cells_rest
      )
    },
    mc.cores = params$mc_cores
  )

  eligibility_df <- bind_rows(lapply(sample_res, `[[`, "eligibility")) %>%
    arrange(state, sample)
  
  per_sample_dge <- bind_rows(lapply(sample_res, `[[`, "markers"))

  if (nrow(per_sample_dge) == 0) {
    stop("Per-sample DGE produced no marker rows. Check sample/state eligibility and FindMarkers inputs.")
  }
  
  message("Caching per-sample DGE results.")
  saveRDS(list(eligibility = eligibility_df, markers = per_sample_dge), sample_dge_cache)
}

# Always write the eligibility CSV for downstream reference
fwrite(
  eligibility_df,
  file.path(out_dir, "Auto_six_state_sample_state_eligibility.csv")
)

per_sample_dge <- per_sample_dge %>%
  mutate(
    hit_flag = !is.na(p_val_adj) &
      p_val_adj < 0.05 &
      avg_log2FC > 0
  ) %>%
  arrange(state, gene, sample)

fwrite(
  per_sample_dge,
  file.path(out_dir, "Auto_six_state_per_sample_dge.csv.gz")
)

state_coverage <- eligibility_df %>%
  filter(eligible) %>%
  group_by(state) %>%
  summarise(
    eligible_sample_n = n_distinct(sample),
    eligible_study_n = n_distinct(study),
    .groups = "drop"
  )

marker_summary <- per_sample_dge %>%
  group_by(state, gene) %>%
  summarise(
    tested_sample_n = n_distinct(sample),
    tested_study_n = n_distinct(study),
    hit_sample_n = n_distinct(sample[hit_flag]),
    hit_study_n = n_distinct(study[hit_flag]),
    hit_sample_pct = hit_sample_n / tested_sample_n,
    median_log2FC_hit = safe_median(avg_log2FC[hit_flag]),
    median_pct_state_hit = safe_median(pct_state[hit_flag]),
    median_pct_other_hit = safe_median(pct_other[hit_flag]),
    median_pct_delta_hit = safe_median(pct_delta[hit_flag]),
    max_log2FC_hit = safe_max(avg_log2FC[hit_flag]),
    min_p_adj_hit = safe_min(p_val_adj[hit_flag]),
    .groups = "drop"
  ) %>%
  left_join(
    global_markers %>%
      select(state, gene, global_mean_state, global_mean_other, global_mean_diff,
             global_pct_state, global_pct_other, global_pct_delta),
    by = c("state", "gene")
  ) %>%
  left_join(state_coverage, by = "state") %>%
  mutate(
    sample_recurrence = ifelse(eligible_sample_n > 0, hit_sample_n / eligible_sample_n, 0),
    study_recurrence = ifelse(eligible_study_n > 0, hit_study_n / eligible_study_n, 0),
    reproducibility_score = 0.5 * sample_recurrence + 0.5 * study_recurrence
  )

####################
# state-level specificity check
####################
specificity_cache <- file.path(cache_dir, "state_specificity.rds")

if (file.exists(specificity_cache)) {
  message("Loading cached state specificity results.")
  specificity_long <- readRDS(specificity_cache)
} else {
  message("Computing sample-aware state specificity summaries for candidate genes.")

  candidate_union <- sort(unique(marker_summary$gene))
  expr_data <- GetAssayData(tmdata_state6, assay = "RNA", slot = "data")[candidate_union, , drop = FALSE]

  state_expr_list <- lapply(state_order, function(state_name) {
    eligible_samples <- eligibility_df %>%
      filter(state == state_name, eligible) %>%
      pull(sample) %>%
      unique()

    if (length(eligible_samples) == 0) {
      return(
        data.frame(
          gene = candidate_union,
          state = state_name,
          median_sample_state_expr = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    state_cells_by_sample <- split(
      colnames(tmdata_state6)[tmdata_state6$final_state6 == state_name & tmdata_state6$orig.ident %in% eligible_samples],
      tmdata_state6$orig.ident[tmdata_state6$final_state6 == state_name & tmdata_state6$orig.ident %in% eligible_samples]
    )

    state_cells_by_sample <- state_cells_by_sample[lengths(state_cells_by_sample) >= params$min_cells_state]

    if (length(state_cells_by_sample) == 0) {
      return(
        data.frame(
          gene = candidate_union,
          state = state_name,
          median_sample_state_expr = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    sample_means <- vapply(
      state_cells_by_sample,
      function(cells_use) Matrix::rowMeans(expr_data[, cells_use, drop = FALSE]),
      FUN.VALUE = numeric(length(candidate_union))
    )

    if (is.null(dim(sample_means))) {
      sample_means <- matrix(sample_means, ncol = 1)
    }

    data.frame(
      gene = candidate_union,
      state = state_name,
      median_sample_state_expr = apply(sample_means, 1, median, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  state_expr_wide <- bind_rows(state_expr_list) %>%
    tidyr::pivot_wider(
      names_from = state,
      values_from = median_sample_state_expr
    )

  state_expr_mat <- state_expr_wide %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  best_state_idx <- max.col(state_expr_mat, ties.method = "first")
  best_state <- colnames(state_expr_mat)[best_state_idx]
  # max_expr <- apply(state_expr_mat, 1, max, na.rm = TRUE) # Not strictly needed if not used below

  specificity_long <- bind_rows(lapply(state_order, function(state_name) {
    state_vals <- state_expr_mat[, state_name]
    other_vals <- state_expr_mat[, setdiff(state_order, state_name), drop = FALSE]
    other_max <- apply(other_vals, 1, max, na.rm = TRUE)

    data.frame(
      gene = rownames(state_expr_mat),
      state = state_name,
      state_median_expr = state_vals,
      off_state_max_median_expr = other_max,
      specificity_gap = state_vals - other_max,
      best_state = best_state,
      stringsAsFactors = FALSE
    )
  }))
  
  message("Caching state specificity results.")
  saveRDS(specificity_long, specificity_cache)
}

# Re-derive state_expr_mat from specificity_long for downstream use (heatmap_expr relies on it)
state_expr_mat <- specificity_long %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

marker_summary <- marker_summary %>%
  left_join(specificity_long, by = c("state", "gene")) %>%
  mutate(
    best_state_match = best_state == state
  )

fwrite(
  marker_summary,
  file.path(out_dir, "Auto_six_state_marker_summary.csv")
)

####################
# final marker definition
####################
message("Selecting top publication-facing markers by reproducibility, effect size, and specificity.")

ranked_markers <- marker_summary %>%
  filter(
    hit_sample_n > 0,
    best_state_match,
    specificity_gap > 0
  ) %>%
  group_by(state) %>%
  mutate(
    reproducibility_rank = dplyr::percent_rank(reproducibility_score),
    effect_rank = dplyr::percent_rank(median_log2FC_hit),
    specificity_rank = dplyr::percent_rank(specificity_gap),
    ranking_score = reproducibility_rank + effect_rank + specificity_rank
  ) %>%
  ungroup() %>%
  arrange(
    state,
    desc(ranking_score),
    desc(reproducibility_score),
    desc(median_log2FC_hit),
    desc(specificity_gap),
    gene
  )

fwrite(
  ranked_markers,
  file.path(out_dir, "Auto_six_state_markers_ranked.csv")
)

final_markers <- ranked_markers %>%
  group_by(state) %>%
  slice_head(n = params$top_markers_per_state) %>%
  ungroup()

fwrite(
  final_markers,
  file.path(out_dir, "Auto_six_state_markers_final.csv")
)

####################
# recurrence reporting
####################
message("Writing sample/study recurrence tables for the top markers.")

top_marker_recurrence_summary <- final_markers %>%
  ####################
  # support classification
  ####################
  mutate(
    legacy_required_sample_n = pmax(8, ceiling(0.20 * eligible_sample_n)),
    legacy_required_study_n = pmax(2, ceiling(0.35 * eligible_study_n)),
    is_multi_sample = hit_sample_n >= 2,
    is_multi_study = hit_study_n >= 2,
    support_class = case_when(
      is_multi_study ~ "multi-study",
      is_multi_sample ~ "multi-sample_single-study",
      hit_sample_n == 1 ~ "single-sample_single-study",
      TRUE ~ "no-positive-hit"
    ),
    passes_legacy_strict_recurrence =
      hit_sample_n >= legacy_required_sample_n &
      hit_study_n >= legacy_required_study_n
  ) %>%
  select(
    state,
    gene,
    hit_sample_n,
    eligible_sample_n,
    sample_recurrence,
    hit_study_n,
    eligible_study_n,
    study_recurrence,
    reproducibility_score,
    median_log2FC_hit,
    specificity_gap,
    ranking_score,
    is_multi_sample,
    is_multi_study,
    support_class,
    legacy_required_sample_n,
    legacy_required_study_n,
    passes_legacy_strict_recurrence
  )

fwrite(
  top_marker_recurrence_summary,
  file.path(out_dir, "Auto_six_state_markers_top5_recurrence_summary.csv")
)

top_marker_sample_support <- per_sample_dge %>%
  semi_join(
    final_markers %>% select(state, gene),
    by = c("state", "gene")
  ) %>%
  select(
    state,
    gene,
    sample,
    study,
    hit_flag,
    p_val_adj,
    avg_log2FC,
    pct_state,
    pct_other,
    pct_delta,
    state_cell_n,
    other_cell_n
  ) %>%
  arrange(state, gene, desc(hit_flag), study, sample)

fwrite(
  top_marker_sample_support,
  file.path(out_dir, "Auto_six_state_markers_top5_sample_support.csv.gz")
)

top_marker_study_support <- top_marker_sample_support %>%
  group_by(state, gene, study) %>%
  summarise(
    eligible_sample_n_in_study = n_distinct(sample),
    hit_sample_n_in_study = n_distinct(sample[hit_flag]),
    study_hit_flag = hit_sample_n_in_study > 0,
    best_log2FC_in_study = safe_max(avg_log2FC[hit_flag]),
    best_specificity_delta_in_study = safe_max(pct_delta[hit_flag]),
    best_p_adj_in_study = safe_min(p_val_adj[hit_flag]),
    .groups = "drop"
  ) %>%
  arrange(state, gene, desc(study_hit_flag), study)

fwrite(
  top_marker_study_support,
  file.path(out_dir, "Auto_six_state_markers_top5_study_support.csv")
)

####################
# state-level sensitivity summary
####################
top_marker_state_recurrence_summary <- top_marker_recurrence_summary %>%
  group_by(state) %>%
  summarise(
    top_marker_n = n(),
    multi_sample_marker_n = sum(is_multi_sample, na.rm = TRUE),
    multi_study_marker_n = sum(is_multi_study, na.rm = TRUE),
    median_sample_recurrence = median(sample_recurrence, na.rm = TRUE),
    median_study_recurrence = median(study_recurrence, na.rm = TRUE),
    n_passing_legacy_strict_recurrence = sum(passes_legacy_strict_recurrence, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  top_marker_state_recurrence_summary,
  file.path(out_dir, "Auto_six_state_markers_top5_state_recurrence_summary.csv")
)

heatmap_markers <- final_markers

fwrite(
  heatmap_markers,
  file.path(out_dir, "Auto_six_state_markers_heatmap_top.csv")
)

####################
# heatmap matrix
####################
message("Building marker heatmap.")

if (nrow(heatmap_markers) == 0) {
  stop("No final markers passed the current recurrent marker thresholds.")
}

heatmap_genes <- heatmap_markers$gene
heatmap_expr <- state_expr_mat[heatmap_genes, state_order, drop = FALSE]

dup_counter <- ave(seq_along(heatmap_genes), heatmap_genes, FUN = seq_along)
heatmap_row_names <- ifelse(
  dup_counter == 1,
  heatmap_genes,
  paste0(heatmap_genes, "__", dup_counter)
)

rownames(heatmap_expr) <- heatmap_row_names
heatmap_z <- row_zscore(heatmap_expr)

heatmap_matrix_out <- as.data.frame(heatmap_expr) %>%
  tibble::rownames_to_column("row_id") %>%
  left_join(
    heatmap_markers %>%
      mutate(row_id = heatmap_row_names) %>%
      select(row_id, gene, state, hit_sample_n, hit_sample_pct, hit_study_n, median_log2FC_hit, specificity_gap),
    by = "row_id"
  ) %>%
  relocate(gene, state, hit_sample_n, hit_sample_pct, hit_study_n, median_log2FC_hit, specificity_gap, .after = row_id)

fwrite(
  heatmap_matrix_out,
  file.path(out_dir, "Auto_six_state_marker_heatmap_matrix.csv")
)

heatmap_state_factor <- factor(heatmap_markers$state, levels = state_order)
names(heatmap_state_factor) <- heatmap_row_names

####################
# publication-facing heatmap
####################
# FIX: Use explicitly named vectors to prevent ComplexHeatmap from getting confused
row_ann <- rowAnnotation(
  State = heatmap_state_factor,
  Sample_support = anno_barplot(
    heatmap_markers$hit_sample_pct * 100,
    gp = gpar(fill = "grey35", col = NA),
    border = FALSE,
    axis_param = list(side = "bottom", at = c(0, 50), labels = c("0", "50")),
    width = unit(18, "mm")
  ),
  Study_support = anno_barplot(
    heatmap_markers$hit_study_n,
    gp = gpar(fill = "grey55", col = NA),
    border = FALSE,
    axis_param = list(side = "bottom"),
    width = unit(15, "mm")
  ),
  col = list(State = state_cols),
  # Use named vectors so the package knows exactly which rule applies to which column
  show_annotation_name = c(State = FALSE, Sample_support = TRUE, Study_support = TRUE), 
  annotation_label = c(State = "State", Sample_support = "Sample\nSupport", Study_support = "Study\nSupport"), 
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9, lineheight = 0.9),
  gap = unit(c(2, 6), "mm"), 
  simple_anno_size = unit(4, "mm")
)

# FIX: Renamed to 'Top_State' to completely prevent the legend-merging bug
top_ann <- HeatmapAnnotation(
  Top_State = factor(state_order, levels = state_order),
  col = list(Top_State = state_cols),
  show_annotation_name = FALSE, 
  show_legend = FALSE, # We only need the legend from row_ann
  simple_anno_size = unit(4, "mm")
)

col_fun <- colorRamp2(
  c(-2.0, 0, 2.5),
  c("#1D4E89", "#F8F4EC", "#B22222")
)

ht <- Heatmap(
  heatmap_z,
  name = "Row Z-score",
  top_annotation = top_ann,
  left_annotation = row_ann,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = heatmap_state_factor,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, 
  border = TRUE,
  heatmap_legend_param = list(
    title = "Relative\nexpression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

pdf(
  file.path(out_dir, "Auto_six_state_marker_heatmap.pdf"),
  width = 11.5,
  height = 12, 
  useDingbats = FALSE
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(2.5, 1), c("cm", "null")))))

# Title block
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text(
  "Recurrent state markers across six finalized epithelial states",
  x = unit(0.5, "npc"),
  y = unit(0.70, "npc"), 
  gp = gpar(fontsize = 16, fontface = "bold")
)
grid.text(
  "Rows are per-state recurrent markers ranked by recurrence, effect size, and state specificity. Columns are median sample-level expression per state, scaled by row.",
  x = unit(0.5, "npc"),
  y = unit(0.30, "npc"), 
  gp = gpar(fontsize = 10)
)
popViewport()

# Heatmap block
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht, newpage = FALSE)
popViewport(2)
dev.off()

####################
# methodology document
####################
message("Writing simplified methodology document.")

state_summary_tbl <- data.frame(
  state = state_order,
  cell_n = as.integer(state_counts[state_order]),
  stringsAsFactors = FALSE
) %>%
  left_join(state_coverage, by = "state") %>%
  left_join(final_markers %>% count(state, name = "top_marker_n"), by = "state") %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))

state_summary_print <- capture.output(print(state_summary_tbl, row.names = FALSE))

state_recurrence_print <- capture.output(print(top_marker_state_recurrence_summary, row.names = FALSE))

####################
# 3CA recurrence note
####################
threeca_recurrence_tbl <- top_marker_recurrence_summary %>%
  filter(state == "3CA_EMT_and_Protein_maturation") %>%
  select(
    gene,
    hit_sample_n,
    hit_study_n,
    sample_recurrence,
    study_recurrence,
    support_class,
    passes_legacy_strict_recurrence
  )

threeca_recurrence_print <- capture.output(print(threeca_recurrence_tbl, row.names = FALSE))

method_lines <- c(
  "# Auto Six-State Marker Methodology",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 1. Goal and Scope",
  "",
  "The workflow aims to identify the top 5 most robust and specific markers for each finalized state. It prioritizes genes that are not only highly expressed but also consistently reproducible across the heterogeneous cohort of samples and studies, while maintaining clear specificity for a single state.",
  "",
  "**Core Inputs:**",
  "- `ref_outs/EAC_Ref_epi.rds`: The main epithelial Seurat object (75,348 cells).",
  "- `ref_outs/Auto_final_states.rds`: Finalized six-state labels, with `Unresolved` and `Hybrid` excluded.",
  "",
  "**Finalized States Retained:**",
  paste0("- `", state_order, "`"),
  "",
  "---",
  "",
  "## 2. Six-State Subset and Re-embedding",
  "",
  "To ensure marker analysis is focused on the finalized transcriptional landscape, a clean subset and re-embedding are performed.",
  "",
  "### 2.1 Cell and Feature Selection",
  "1. The global epithelial count matrix is subset to cells present in both inputs.",
  "2. Cells with `Unresolved` or `Hybrid` labels are removed.",
  paste0("3. Genes detected in fewer than **", params$min_cells_feature, " cells** within this six-state subset are discarded."),
  "",
  "### 2.2 Re-embedding Pipeline",
  "The subsetted object is processed through a standard Seurat pipeline:",
  "- **Normalization:** `NormalizeData` (standard log-normalization).",
  paste0("- **HVGs:** `FindVariableFeatures` (vst method, **", params$n_variable_features, " features**)."),
  "- **Scaling:** `ScaleData` on the HVGs.",
  paste0("- **PCA:** `RunPCA` (**", params$n_pcs, " PCs**)."),
  paste0("- **Neighbors:** `FindNeighbors` on ", params$n_pcs, " PCs."),
  paste0("- **Clustering:** `FindClusters` (Leiden algorithm, **resolution ", params$cluster_resolution, "**)."),
  paste0("- **UMAP:** `RunUMAP` on ", params$n_pcs, " PCs."),
  "",
  "---",
  "",
  "## 3. Candidate Gene Screening (Pooled)",
  "",
  "Because sample-wise DGE is computationally expensive, a global descriptive screen is used to define a tractable candidate set for detailed recurrence testing.",
  "",
  "For each state:",
  "1. Total cells are divided into **State** vs. **Rest**.",
  "2. Mean expression and detection frequency (pct) are computed for both groups.",
  "3. Genes are ranked by `global_mean_difference` (State - Rest).",
  paste0("4. The **top ", params$candidate_pool_per_state, "** genes with positive mean difference are retained as candidates for that state."),
  "",
  "---",
  "",
  "## 4. Sample-Wise Recurrent DGE",
  "",
  "The core of the reproducibility analysis is running DGE within individual samples to see how often a gene is a \"hit.\"",
  "",
  "### 4.1 Sample Eligibility",
  "A sample is considered \"eligible\" to test a specific state only if:",
  paste0("- The sample contains **at least ", params$min_cells_state, " cells** belonging to the target state."),
  paste0("- The sample contains **at least ", params$min_cells_rest, " cells** belonging to the other five states combined (the \"Rest\")."),
  "",
  "### 4.2 Differential Expression Testing",
  "For each eligible sample and each of its qualified states:",
  "- **Test:** Seurat `FindMarkers` (Wilcoxon rank-sum test).",
  "- **Universe:** One state versus the other 5 states within that sample.",
  paste0("- **Genes:** Only the ", params$candidate_pool_per_state, " candidate genes identified in the pooled screen for that state."),
  "- **Thresholds:** No hard expression or logFC gates are applied at the testing stage (`logfc.threshold = 0`, `min.pct = 0`) to capture all valid statistics.",
  "",
  "### 4.3 Hit Definition",
  "A gene is defined as a \"hit\" within a specific sample if:",
  "- `p_val_adj < 0.05` (FDR corrected).",
  "- `avg_log2FC > 0` (statistically higher in the target state).",
  "",
  "---",
  "",
  "## 5. Sample-Aware State Specificity",
  "",
  "To ensure markers are globally specific across the atlas, a \"specificity gap\" is computed using sample-level medians.",
  "",
  "1. For every gene in the marker summary, the mean expression is calculated within the target state for every sample eligible for that state.",
  "2. The **median of these sample-level means** is computed, representing the \"typical\" expression of that gene in that state.",
  "3. This is repeated for all six states.",
  "4. **Specificity Gap:** The typical expression in the target state minus the maximum typical expression seen in any of the other five states.",
  "5. A gene is considered a \"best state match\" only if the target state has the highest median expression.",
  "",
  "---",
  "",
  "## 6. Ranking and Final Selection",
  "",
  "Genes are ranked within each state using a multi-component reproducibility and specificity score.",
  "",
  "### 6.1 Ranking Metrics",
  "Three metrics are computed per gene/state:",
  "1. **Reproducibility Score:** ",
  "   - `sample_recurrence`: Fraction of eligible samples that are DGE hits.",
  "   - `study_recurrence`: Fraction of eligible studies that have at least one DGE hit sample.",
  "   - `score = 0.5 * sample_recurrence + 0.5 * study_recurrence`.",
  "2. **Effect Size:** Median `avg_log2FC` across all samples that were DGE hits.",
  "3. **Specificity Gap:** The gap computed in Section 5.",
  "",
  "### 6.2 The Ranking Score",
  "Within each state, genes are assigned a `percent_rank` (0 to 1) for each of the three metrics above. The final **Ranking Score** is the sum of these three ranks.",
  "",
  "### 6.3 Hard Selection Rules",
  "To be considered for the final top 5, a gene MUST:",
  "- Have at least one significant positive DGE hit in at least one sample (`hit_sample_n > 0`).",
  "- Have the target state as its highest-expressing state globally (`best_state_match == TRUE`).",
  "- Have a positive specificity gap (`specificity_gap > 0`).",
  "",
  paste0("The top ", params$top_markers_per_state, " markers by **Ranking Score** per state are selected for the final panel."),
  "",
  "---",
  "",
  "## 7. Recurrence and Support Classification",
  "",
  "- **Support Class:**",
  "  - `multi-study`: Hit detected in samples from 2 or more studies.",
  "  - `multi-sample_single-study`: Hit detected in 2 or more samples, but all within one study.",
  "  - `single-sample_single-study`: Hit detected in only 1 sample.",
  "- **Legacy Strict Check:** The markers are also checked against a \"Legacy Strict\" rule (`hit_sample_n >= 20% of eligible` and `hit_study_n >= 35% of eligible`). This highlights which public markers are extremely robust vs. those that are specific but lower sensitivity.",
  "",
  "---",
  "",
  "## 8. Heatmap Construction",
  "",
  "- **Data:** Median of sample-level means per state (as computed in Section 5).",
  "- **Z-scoring:** Values are Z-scored per row across the six states to highlight state-specific enrichment.",
  "",
  "## 9. Current 3CA recurrence profile",
  "",
  "```text",
  threeca_recurrence_print,
  "```",
  "",
  "## 10. State-level summary",
  "",
  "```text",
  state_summary_print,
  "```",
  "",
  "## 11. Recurrence summary by state",
  "",
  "```text",
  state_recurrence_print,
  "```",
  "",
  "## 12. Output Files",
  "",
  "- `Auto_six_state_markers_final.csv`: The final top 5 markers per state with their ranking and recurrence stats.",
  "- `Auto_six_state_markers_ranked.csv`: The full table of candidate genes ranked by the workflow.",
  "- `Auto_six_state_markers_top5_recurrence_summary.csv`: Summary of hit counts and support classes.",
  "- `Auto_six_state_markers_top5_sample_support.csv.gz`: per-sample support table for the final top-5 markers.",
  "- `Auto_six_state_markers_top5_study_support.csv`: per-study support table for the final top-5 markers.",
  "- `Auto_six_state_markers_top5_state_recurrence_summary.csv`: state-level sensitivity summary for the final top-5 markers.",
  "- `Auto_six_state_marker_heatmap.pdf`: final publication-facing heatmap.",
  "- `Auto_six_state_umap.pdf`: UMAP visualizations of the six-state subset.",
  ""
)

writeLines(
  method_lines, 
  file.path("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/methodology/cell_states/final_state_marker_discovery_methodology.md")
)
