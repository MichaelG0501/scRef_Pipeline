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
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_final.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_markers_heatmap_top.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_heatmap_matrix.csv
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_heatmap.pdf
#   ref_outs/Auto_six_state_markers/Auto_six_state_marker_methodology.md
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
  min_cells_state = 20,
  min_cells_rest = 20,
  global_min_pct = 0.05,
  global_candidate_min_pct_state = 0.10,
  global_candidate_min_delta_pct = 0.05,
  global_candidate_min_logfc = 0.10,
  candidate_pool_per_state = 1500,
  per_sample_min_logfc = 0.25,
  per_sample_min_pct_state = 0.25,
  per_sample_min_delta_pct = 0.10,
  final_min_logfc = 0.25,
  final_min_pct_state = 0.25,
  final_min_delta_pct = 0.10,
  min_recurrence_fraction = 0.20,
  min_recurrence_floor = 8,
  min_recurrence_study_fraction = 0.35,
  min_recurrence_study_floor = 2,
  heatmap_genes_per_state = 15,
  min_cells_feature = 20,
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
  filter(
    global_mean_diff >= params$global_candidate_min_logfc,
    global_pct_state >= params$global_candidate_min_pct_state,
    global_pct_delta >= params$global_candidate_min_delta_pct
  ) %>%
  group_by(state) %>%
  arrange(desc(global_mean_diff), desc(global_pct_delta), desc(global_pct_state), .by_group = TRUE) %>%
  slice_head(n = params$candidate_pool_per_state) %>%
  summarise(candidate_genes = list(unique(gene)), .groups = "drop")

candidate_map <- setNames(candidate_map$candidate_genes, candidate_map$state)
candidate_map <- candidate_map[state_order]
candidate_map[sapply(candidate_map, is.null)] <- list(character())

rm(expr_global)
invisible(gc())

####################
# per-sample DGE
####################
sample_dge_cache <- file.path(cache_dir, "per_sample_dge_full.rds")

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
      avg_log2FC >= params$per_sample_min_logfc &
      pct_state >= params$per_sample_min_pct_state &
      pct_delta >= params$per_sample_min_delta_pct
  ) %>%
  arrange(state, gene, sample)

fwrite(
  per_sample_dge,
  file.path(out_dir, "Auto_six_state_per_sample_dge.csv.gz")
)

recurrence_requirements <- eligibility_df %>%
  filter(eligible) %>%
  group_by(state) %>%
  summarise(
    eligible_sample_n = n_distinct(sample),
    eligible_study_n = n_distinct(study),
    required_sample_n = pmax(
      params$min_recurrence_floor,
      ceiling(params$min_recurrence_fraction * eligible_sample_n)
    ),
    required_study_n = pmax(
      params$min_recurrence_study_floor,
      ceiling(params$min_recurrence_study_fraction * eligible_study_n)
    ),
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
    median_log2FC_hit = median(avg_log2FC[hit_flag], na.rm = TRUE),
    median_pct_state_hit = median(pct_state[hit_flag], na.rm = TRUE),
    median_pct_other_hit = median(pct_other[hit_flag], na.rm = TRUE),
    median_pct_delta_hit = median(pct_delta[hit_flag], na.rm = TRUE),
    max_log2FC_hit = max(avg_log2FC[hit_flag], na.rm = TRUE),
    min_p_adj_hit = min(p_val_adj[hit_flag], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        median_log2FC_hit,
        median_pct_state_hit,
        median_pct_other_hit,
        median_pct_delta_hit,
        max_log2FC_hit,
        min_p_adj_hit
      ),
      ~ ifelse(is.infinite(.x), NA_real_, .x)
    )
  ) %>%
  left_join(
    global_markers %>%
      select(state, gene, global_mean_state, global_mean_other, global_mean_diff,
             global_pct_state, global_pct_other, global_pct_delta),
    by = c("state", "gene")
  ) %>%
  left_join(recurrence_requirements, by = "state")

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

best_state_map <- specificity_long %>%
  select(gene, best_state) %>%
  distinct() %>%
  tibble::deframe()

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
message("Selecting final recurrent and state-specific markers.")

final_markers <- marker_summary %>%
  filter(
    global_mean_diff >= params$global_candidate_min_logfc,
    global_pct_state >= params$final_min_pct_state,
    global_pct_delta >= params$final_min_delta_pct,
    hit_sample_n >= required_sample_n,
    hit_study_n >= required_study_n,
    median_log2FC_hit >= params$final_min_logfc,
    median_pct_state_hit >= params$final_min_pct_state,
    median_pct_delta_hit >= params$final_min_delta_pct,
    best_state_match,
    specificity_gap > 0
  ) %>%
  mutate(
    ranking_score = (2 * hit_sample_pct) + median_log2FC_hit + global_mean_diff + pmax(specificity_gap, 0)
  ) %>%
  arrange(state, desc(ranking_score), desc(hit_sample_pct), desc(median_log2FC_hit), desc(specificity_gap), gene)

fwrite(
  final_markers,
  file.path(out_dir, "Auto_six_state_markers_final.csv")
)

heatmap_markers <- final_markers %>%
  group_by(state) %>%
  slice_head(n = params$heatmap_genes_per_state) %>%
  ungroup()

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