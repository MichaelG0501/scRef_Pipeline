####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_inspect_expression.R
#   Description:
#     Combined script for TF expression inspection:
#     1) Heatmap of candidate genes across states (scAtlas vs PDO).
#     2) Per-sample scatter plots correlating primary TFs with state abundance & DGE UCell.
#     3) Regulon-level scatter: for every SCENIC regulon, plots mean AUCell activity (x)
#        vs Spearman rho with state metrics (y). 4-page PDF, one page per overlapping
#        state pair, 8 panels per page (top=scAtlas, bottom=PDO). Top 10 regulons
#        labelled; shared top 10 across datasets highlighted in red.
#   Outputs:
#     figures/:
#       final_mp_scenic/inspect_expression/Auto_gene_expression_state_heatmap.pdf
#       final_mp_scenic/inspect_expression/Auto_E2F8_expression_correlation.pdf
#       final_mp_scenic/inspect_expression/Auto_SNAI1_expression_correlation.pdf
#       final_mp_scenic/inspect_expression/Auto_regulon_state_correlation.pdf
#     tables/:
#       final_mp_scenic/inspect_expression/Auto_gene_expression_state_heatmap_values.csv
#       final_mp_scenic/inspect_expression/Auto_expression_correlation_data.csv
#       final_mp_scenic/inspect_expression/Auto_regulon_state_correlation_data.csv
#   Run command: qsub analysis/cell_states/Auto_inspect_expression.sh
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(UCell)
  library(data.table)
  library(AUCell)
  library(ggrepel)
})

sc_root  <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
pdo_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(file.path(sc_root, "ref_outs"))
out_dir <- "final_mp_scenic/inspect_expression"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

####################
# Configuration (Heatmap)
####################
gene_groups <- list(
  "Cell cycle"            = c("PTTG1", "FEN1"),
  "Classic proliferative" = c("E2F8", "MAZ"),
  "Stress-adaptive"       = c("SNAI1", "RXRB")
)
all_genes <- unlist(gene_groups, use.names = FALSE)

sc_state_order <- c("Classic proliferation", "Squamous-to-intestinal",
                     "Glandular-to-intestinal", "Stress-adaptive",
                     "Cancer-cell immune mimicry")
sc_state_cols <- c("Classic proliferation" = "#E41A1C",
                   "Squamous-to-intestinal" = "#4DAF4A",
                   "Glandular-to-intestinal" = "#FF7F00",
                   "Stress-adaptive" = "#984EA3",
                   "Cancer-cell immune mimicry" = "#377EB8")

pdo_state_order <- c("Classic proliferation", "Columnar-to-intestinal",
                     "Glandular differentiation", "Stress-adaptive",
                     "ECM-remodelling", "Motile-cilia differentiation")
pdo_state_cols <- c("Classic proliferation" = "#E41A1C",
                    "Columnar-to-intestinal" = "#4DAF4A",
                    "Glandular differentiation" = "#FF7F00",
                    "Stress-adaptive" = "#984EA3",
                    "ECM-remodelling" = "#A65628",
                    "Motile-cilia differentiation" = "#F781BF")

group_cols <- c("Cell cycle" = "#6B7280",
                "Classic proliferative" = "#E41A1C",
                "Stress-adaptive" = "#984EA3")

####################
# Configuration (Gene Correlation)
####################
primary_genes <- list(E2F8 = "Classic proliferation", SNAI1 = "Stress-adaptive")
neg_ctrl <- "IRF6"
tasks_to_run <- data.frame(
  gene = c("E2F8", "IRF6", "SNAI1", "IRF6"),
  target_state = c("Classic proliferation", "Classic proliferation",
                    "Stress-adaptive", "Stress-adaptive"),
  stringsAsFactors = FALSE
) %>% distinct()
top_n_markers <- 5

####################
# Configuration (Regulon Correlation)
####################
state_pairs <- list(
  list(sc_state = "Classic proliferation",
       pdo_state = "Classic proliferation",
       label = "Classic proliferation"),
  list(sc_state = "Squamous-to-intestinal",
       pdo_state = "Columnar-to-intestinal",
       label = "Intestinal transition"),
  list(sc_state = "Glandular-to-intestinal",
       pdo_state = "Glandular differentiation",
       label = "Glandular"),
  list(sc_state = "Stress-adaptive",
       pdo_state = "Stress-adaptive",
       label = "Stress-adaptive")
)

####################
# Helper functions (shared)
####################
get_assay_matrix <- function(seurat_obj, slot_name = c("counts", "data")) {
  slot_name <- match.arg(slot_name)
  mat <- tryCatch(GetAssayData(seurat_obj, assay = "RNA", slot = slot_name),
                  error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  mat <- tryCatch(LayerData(seurat_obj, assay = "RNA", layer = slot_name),
                  error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  stop("Unable to retrieve RNA matrix.")
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

format_regulon_name <- function(x) {
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  x <- gsub(" \\([0-9]+ genes\\)$", "", x)
  x <- gsub("_extended$", "", x)
  x
}

safe_spearman <- function(x, y, min_n = 5) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < min_n) return(NA_real_)
  suppressWarnings(cor(x[valid], y[valid], method = "spearman"))
}

####################
# Helper: Heatmap expression computation
####################
compute_state_expr <- function(seurat_obj, state_vec, state_order,
                               sample_col = "orig.ident", min_cells = 10) {
  DefaultAssay(seurat_obj) <- "RNA"
  norm_data <- tryCatch(get_assay_matrix(seurat_obj, "data"), error = function(e) NULL)
  if (is.null(norm_data) || max(norm_data, na.rm = TRUE) == 0) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    norm_data <- get_assay_matrix(seurat_obj, "data")
  }
  common_cells <- intersect(colnames(norm_data), names(state_vec))
  state_vec <- as.character(state_vec[common_cells])
  keep <- state_vec %in% state_order
  common_cells <- common_cells[keep]
  state_vec <- state_vec[keep]

  genes_present <- intersect(all_genes, rownames(norm_data))
  norm_data <- norm_data[genes_present, common_cells, drop = FALSE]
  sample_ids <- as.character(seurat_obj@meta.data[common_cells, sample_col])

  expr_mat <- matrix(NA_real_, nrow = length(genes_present), ncol = length(state_order),
                     dimnames = list(genes_present, state_order))
  for (st in state_order) {
    st_cells <- common_cells[state_vec == st]
    st_samples <- sample_ids[state_vec == st]
    if (length(st_cells) < min_cells) next
    sample_split <- split(st_cells, st_samples)
    sample_split <- sample_split[lengths(sample_split) >= min_cells]
    if (length(sample_split) == 0) {
      expr_mat[, st] <- Matrix::rowMeans(norm_data[, st_cells, drop = FALSE])
      next
    }
    sample_means <- vapply(sample_split,
                           function(cells) Matrix::rowMeans(norm_data[, cells, drop = FALSE]),
                           FUN.VALUE = numeric(length(genes_present)))
    if (is.null(dim(sample_means))) sample_means <- matrix(sample_means, ncol = 1)
    expr_mat[, st] <- apply(sample_means, 1, median, na.rm = TRUE)
  }
  expr_mat
}

ensure_gene_order <- function(expr_mat, genes, state_order) {
  out <- matrix(NA_real_, nrow = length(genes), ncol = length(state_order),
                dimnames = list(genes, state_order))
  present <- intersect(genes, rownames(expr_mat))
  out[present, ] <- expr_mat[present, state_order, drop = FALSE]
  out
}

####################
# Helper: Gene-level per-sample correlation computation
####################
compute_sample_data <- function(seurat_obj, state_vec, state_order,
                                abundance_df, auc_mat, gene_name,
                                target_state, dataset_label) {
  dge_genes <- dge_gene_sets[[target_state]]
  norm_data <- tryCatch(get_assay_matrix(seurat_obj, "data"), error = function(e) NULL)
  if (is.null(norm_data) || max(norm_data, na.rm = TRUE) == 0) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    norm_data <- get_assay_matrix(seurat_obj, "data")
  }
  common_cells <- intersect(colnames(norm_data), names(state_vec))
  state_vec <- as.character(state_vec[common_cells])
  keep <- state_vec %in% state_order
  common_cells <- common_cells[keep]
  sample_ids <- as.character(seurat_obj@meta.data[common_cells, "orig.ident"])
  names(sample_ids) <- common_cells

  if (gene_name %in% rownames(norm_data)) {
    gene_expr <- as.numeric(norm_data[gene_name, common_cells])
    names(gene_expr) <- common_cells
  } else {
    gene_expr <- setNames(rep(NA_real_, length(common_cells)), common_cells)
  }
  sample_mean_expr <- tapply(gene_expr, sample_ids, mean, na.rm = TRUE)

  if (!is.null(auc_mat)) {
    regulon_names <- rownames(auc_mat)
    clean_names <- format_regulon_name(regulon_names)
    matching_idx <- which(clean_names == gene_name)
    if (length(matching_idx) > 0) {
      non_ext <- matching_idx[!grepl("extended", regulon_names[matching_idx])]
      regulon_use <- if (length(non_ext) > 0) regulon_names[non_ext[1]] else regulon_names[matching_idx[1]]
      auc_cells <- intersect(colnames(auc_mat), common_cells)
      if (length(auc_cells) > 0) {
        regulon_scores <- as.numeric(auc_mat[regulon_use, auc_cells])
        names(regulon_scores) <- auc_cells
        sample_mean_auc <- tapply(regulon_scores, sample_ids[auc_cells], mean, na.rm = TRUE)
      } else {
        sample_mean_auc <- setNames(rep(NA_real_, length(unique(sample_ids))), unique(sample_ids))
      }
    } else {
      sample_mean_auc <- setNames(rep(NA_real_, length(unique(sample_ids))), unique(sample_ids))
    }
  } else {
    sample_mean_auc <- setNames(rep(NA_real_, length(unique(sample_ids))), unique(sample_ids))
  }

  state_abundance <- abundance_df %>%
    filter(label == target_state) %>%
    select(orig.ident, pct) %>%
    tibble::deframe()

  dge_genes_present <- intersect(dge_genes, rownames(norm_data))
  if (length(dge_genes_present) >= 2) {
    ucell_scores <- tryCatch(
      ScoreSignatures_UCell(matrix = norm_data[, common_cells, drop = FALSE],
                            features = list(DGE_sig = dge_genes_present), name = ""),
      error = function(e) NULL)
    if (!is.null(ucell_scores)) {
      sample_mean_ucell <- tapply(ucell_scores[, 1], sample_ids, mean, na.rm = TRUE)
    } else {
      sample_mean_ucell <- setNames(rep(NA_real_, length(unique(sample_ids))), unique(sample_ids))
    }
  } else {
    sample_mean_ucell <- setNames(rep(NA_real_, length(unique(sample_ids))), unique(sample_ids))
  }

  all_samples <- unique(sample_ids)
  res <- data.frame(
    sample = all_samples, dataset = dataset_label, gene = gene_name,
    target_state = target_state,
    mean_log_expr = as.numeric(sample_mean_expr[all_samples]),
    mean_regulon_auc = as.numeric(sample_mean_auc[all_samples]),
    state_abundance = as.numeric(state_abundance[all_samples]),
    mean_dge_ucell = as.numeric(sample_mean_ucell[all_samples]),
    stringsAsFactors = FALSE
  )
  res$state_abundance[is.na(res$state_abundance)] <- 0
  res
}

####################
# Helper: Regulon-level state metrics for one dataset (all state pairs)
####################
compute_regulon_data_for_dataset <- function(
    seurat_obj, state_vec, state_order, abundance_df, auc_mat,
    pairs_for_dataset, ucell_cache, dataset_label
) {
  ## Extract normalised expression
  norm_data <- get_assay_matrix(seurat_obj, "data")
  if (max(norm_data, na.rm = TRUE) == 0) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    norm_data <- get_assay_matrix(seurat_obj, "data")
  }

  ## Common cells with canonical state assignments
  common_cells <- intersect(colnames(norm_data), names(state_vec))
  sv <- as.character(state_vec[common_cells])
  names(sv) <- common_cells
  keep <- sv %in% state_order
  common_cells <- common_cells[keep]
  sv <- sv[keep]
  sample_ids <- setNames(
    as.character(seurat_obj@meta.data[common_cells, "orig.ident"]),
    common_cells
  )
  sample_levels <- sort(unique(sample_ids))

  ## AUC cells
  auc_cells <- intersect(colnames(auc_mat), common_cells)

  ## Deduplicate regulons: prefer non-extended
  regulon_raw <- rownames(auc_mat)
  clean_names <- format_regulon_name(regulon_raw)
  unique_tfs <- unique(clean_names)
  tf_to_reg <- setNames(character(length(unique_tfs)), unique_tfs)
  for (tf in unique_tfs) {
    idx <- which(clean_names == tf)
    non_ext <- idx[!grepl("extended", regulon_raw[idx])]
    tf_to_reg[tf] <- if (length(non_ext) > 0) regulon_raw[non_ext[1]] else regulon_raw[idx[1]]
  }

  ## Build per-sample mean AUCell for every regulon [TFs x samples]
  message("  Computing per-sample mean AUCell for ", length(unique_tfs), " regulons...")
  auc_selected <- as.matrix(auc_mat[tf_to_reg, auc_cells, drop = FALSE])
  rownames(auc_selected) <- names(tf_to_reg)
  regulon_sample_auc <- matrix(NA_real_, nrow = length(unique_tfs), ncol = length(sample_levels),
                               dimnames = list(unique_tfs, sample_levels))
  for (s in sample_levels) {
    s_cells <- auc_cells[sample_ids[auc_cells] == s]
    if (length(s_cells) > 0) {
      regulon_sample_auc[, s] <- rowMeans(auc_selected[, s_cells, drop = FALSE], na.rm = TRUE)
    }
  }

  ## Build per-sample mean gene expression for every TF gene [TFs x samples]
  tfs_in_expr <- intersect(unique_tfs, rownames(norm_data))
  message("  Computing per-sample mean expression for ", length(tfs_in_expr), " TF genes...")
  tf_sample_expr <- matrix(NA_real_, nrow = length(unique_tfs), ncol = length(sample_levels),
                           dimnames = list(unique_tfs, sample_levels))
  if (length(tfs_in_expr) > 0) {
    for (s in sample_levels) {
      s_cells <- common_cells[sample_ids[common_cells] == s]
      if (length(s_cells) > 0) {
        tf_sample_expr[tfs_in_expr, s] <- Matrix::rowMeans(
          norm_data[tfs_in_expr, s_cells, drop = FALSE]
        )
      }
    }
  }

  ## Loop over state pairs
  results_list <- list()
  for (sp in pairs_for_dataset) {
    target_state <- sp$state       # state name in this dataset
    sc_dge_state <- sp$sc_state    # scAtlas state name (for DGE marker lookup)
    label <- sp$label

    message("  State: ", target_state, " (label: ", label, ")")

    ## X-axis: mean AUCell of each regulon in cells of this state
    state_auc_cells <- auc_cells[sv[auc_cells] == target_state]
    if (length(state_auc_cells) < 10) {
      message("    Skipping: only ", length(state_auc_cells), " cells")
      next
    }
    mean_state_auc <- rowMeans(auc_selected[, state_auc_cells, drop = FALSE], na.rm = TRUE)

    ## State abundance per sample
    abund <- abundance_df %>%
      filter(label == target_state) %>%
      select(orig.ident, pct) %>%
      tibble::deframe()

    ## DGE UCell per sample (check cache first)
    if (target_state %in% names(ucell_cache)) {
      dge_ucell <- ucell_cache[[target_state]]
      message("    Using cached DGE UCell")
    } else {
      dge_genes_for_state <- dge_gene_sets[[sc_dge_state]]
      dge_present <- intersect(dge_genes_for_state, rownames(norm_data))
      if (length(dge_present) >= 2) {
        message("    Computing UCell for DGE markers: ", paste(dge_present, collapse = ", "))
        uc <- tryCatch(
          ScoreSignatures_UCell(
            matrix = norm_data[, common_cells, drop = FALSE],
            features = list(DGE = dge_present), name = ""
          ),
          error = function(e) { message("    UCell error: ", e$message); NULL }
        )
        dge_ucell <- if (!is.null(uc)) {
          tapply(uc[, 1], sample_ids[common_cells], mean, na.rm = TRUE)
        } else {
          setNames(rep(NA_real_, length(sample_levels)), sample_levels)
        }
      } else {
        dge_ucell <- setNames(rep(NA_real_, length(sample_levels)), sample_levels)
      }
    }

    ## Compute correlations for each regulon
    common_samp <- intersect(sample_levels, intersect(names(abund), names(dge_ucell)))

    reg_df <- data.frame(
      regulon = unique_tfs,
      mean_state_auc = as.numeric(mean_state_auc[unique_tfs]),
      rho_expr_abundance = NA_real_,
      rho_auc_abundance  = NA_real_,
      rho_expr_dge       = NA_real_,
      rho_auc_dge        = NA_real_,
      dataset     = dataset_label,
      target_state = target_state,
      state_label  = label,
      stringsAsFactors = FALSE
    )

    for (k in seq_along(unique_tfs)) {
      tf <- unique_tfs[k]
      auc_v <- regulon_sample_auc[tf, common_samp]
      reg_df$rho_auc_abundance[k] <- safe_spearman(auc_v, abund[common_samp])
      reg_df$rho_auc_dge[k]       <- safe_spearman(auc_v, dge_ucell[common_samp])
      if (tf %in% tfs_in_expr) {
        expr_v <- tf_sample_expr[tf, common_samp]
        reg_df$rho_expr_abundance[k] <- safe_spearman(expr_v, abund[common_samp])
        reg_df$rho_expr_dge[k]       <- safe_spearman(expr_v, dge_ucell[common_samp])
      }
    }
    results_list[[label]] <- reg_df
  }
  bind_rows(results_list)
}

####################
# Load DGE markers for all overlapping states
####################
message("Loading ranked state markers.")
ranked_markers <- fread("Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv")
dge_gene_sets <- list()

all_sc_states_needed <- unique(c(
  sapply(state_pairs, `[[`, "sc_state"),
  unique(tasks_to_run$target_state)
))
for (ts in all_sc_states_needed) {
  top_genes <- head(ranked_markers[state == ts]$gene, top_n_markers)
  dge_gene_sets[[ts]] <- top_genes
  message("  ", ts, " -> ", paste(top_genes, collapse = ", "))
}

####################
# Process scAtlas
####################
message("\n=== scAtlas ===")
tmdata_sc  <- readRDS("EAC_Ref_epi.rds")
sc_states  <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
sc_abundance <- fread("Metaprogrammes_Results/centred/state_definition/tables/centred_refined_noreg_sample_state_abundance.csv")
sc_auc     <- readRDS(file.path(sc_root, "ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulon_auc.rds"))
sc_auc_mat <- if (inherits(sc_auc, "aucellResults")) getAUC(sc_auc) else as.matrix(sc_auc)

# 1. Heatmap expression
sc_expr <- compute_state_expr(tmdata_sc, sc_states, sc_state_order, "orig.ident")

# 2. Gene-level correlation
sc_results <- list()
for (i in 1:nrow(tasks_to_run)) {
  sc_results[[i]] <- compute_sample_data(
    tmdata_sc, sc_states, sc_state_order, sc_abundance, sc_auc_mat,
    tasks_to_run$gene[i], tasks_to_run$target_state[i], "scAtlas"
  )
}

# 3. Regulon-level correlation
# Build UCell cache from gene-level results
sc_ucell_cache <- list()
for (r in sc_results) {
  ts <- r$target_state[1]
  if (!ts %in% names(sc_ucell_cache)) {
    sc_ucell_cache[[ts]] <- setNames(r$mean_dge_ucell, r$sample)
  }
}

sc_pairs_info <- lapply(state_pairs, function(sp) {
  list(state = sp$sc_state, sc_state = sp$sc_state, label = sp$label)
})
message("Computing regulon-level metrics for scAtlas...")
sc_regulon_data <- compute_regulon_data_for_dataset(
  tmdata_sc, sc_states, sc_state_order, sc_abundance, sc_auc_mat,
  sc_pairs_info, sc_ucell_cache, "scAtlas"
)
message("  scAtlas regulon data: ", nrow(sc_regulon_data), " rows")

rm(tmdata_sc, sc_states, sc_auc, sc_auc_mat)
invisible(gc())

####################
# Process PDO
####################
message("\n=== PDO ===")
tmdata_pdo  <- readRDS(file.path(pdo_root, "PDOs_outs/PDOs_merged.rds"))
pdo_states  <- readRDS(file.path(pdo_root, "PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds"))
pdo_abundance <- fread(file.path(pdo_root, "PDOs_outs/centred_mp_refinement/tables/centred_refined_noreg_sample_state_abundance.csv"))
pdo_auc     <- readRDS(file.path(pdo_root, "PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_regulon_auc.rds"))
pdo_auc_mat <- if (inherits(pdo_auc, "aucellResults")) getAUC(pdo_auc) else as.matrix(pdo_auc)

# 1. Heatmap expression
pdo_expr <- compute_state_expr(tmdata_pdo, pdo_states, pdo_state_order, "orig.ident")

# 2. Gene-level correlation
pdo_results <- list()
for (i in 1:nrow(tasks_to_run)) {
  pdo_results[[i]] <- compute_sample_data(
    tmdata_pdo, pdo_states, pdo_state_order, pdo_abundance, pdo_auc_mat,
    tasks_to_run$gene[i], tasks_to_run$target_state[i], "PDO"
  )
}

# 3. Regulon-level correlation
pdo_ucell_cache <- list()
for (r in pdo_results) {
  ts <- r$target_state[1]
  if (!ts %in% names(pdo_ucell_cache)) {
    pdo_ucell_cache[[ts]] <- setNames(r$mean_dge_ucell, r$sample)
  }
}

pdo_pairs_info <- lapply(state_pairs, function(sp) {
  list(state = sp$pdo_state, sc_state = sp$sc_state, label = sp$label)
})
message("Computing regulon-level metrics for PDO...")
pdo_regulon_data <- compute_regulon_data_for_dataset(
  tmdata_pdo, pdo_states, pdo_state_order, pdo_abundance, pdo_auc_mat,
  pdo_pairs_info, pdo_ucell_cache, "PDO"
)
message("  PDO regulon data: ", nrow(pdo_regulon_data), " rows")

rm(tmdata_pdo, pdo_states, pdo_auc, pdo_auc_mat)
invisible(gc())

####################
# Save data outputs
####################
# Heatmap values
sc_expr  <- ensure_gene_order(sc_expr,  all_genes, sc_state_order)
pdo_expr <- ensure_gene_order(pdo_expr, all_genes, pdo_state_order)
values_df <- data.frame(
  gene = rep(all_genes, times = length(sc_state_order) + length(pdo_state_order)),
  state = c(rep(sc_state_order, each = length(all_genes)),
            rep(pdo_state_order, each = length(all_genes))),
  dataset = c(rep("scAtlas", length(all_genes) * length(sc_state_order)),
              rep("PDO",     length(all_genes) * length(pdo_state_order))),
  expression = c(as.vector(sc_expr), as.vector(pdo_expr)),
  stringsAsFactors = FALSE
)
write.csv(values_df, file.path(out_dir, "Auto_gene_expression_state_heatmap_values.csv"),
          row.names = FALSE)

# Gene correlation values
all_data <- bind_rows(bind_rows(sc_results), bind_rows(pdo_results))
fwrite(all_data, file.path(out_dir, "Auto_expression_correlation_data.csv"))

# Regulon correlation values
regulon_data <- bind_rows(sc_regulon_data, pdo_regulon_data)
fwrite(regulon_data, file.path(out_dir, "Auto_regulon_state_correlation_data.csv"))
message("Saved ", nrow(regulon_data), " regulon-level rows.")

####################
# ========== PLOT 1: HEATMAP ==========
####################
message("\n=== Plotting heatmap ===")
sc_z  <- row_zscore(sc_expr)
pdo_z <- row_zscore(pdo_expr)
gene_group_vec <- rep(names(gene_groups), times = lengths(gene_groups))
names(gene_group_vec) <- all_genes
gene_group_factor <- factor(gene_group_vec, levels = names(gene_groups))

format_expr <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
sc_labels  <- matrix(format_expr(sc_expr),  nrow = nrow(sc_expr),  ncol = ncol(sc_expr))
pdo_labels <- matrix(format_expr(pdo_expr), nrow = nrow(pdo_expr), ncol = ncol(pdo_expr))

all_z <- c(as.vector(sc_z), as.vector(pdo_z))
all_z <- all_z[is.finite(all_z)]
z_limit <- min(max(abs(range(all_z))), 2.5)
col_fun <- colorRamp2(c(-z_limit, 0, z_limit), c("#1D4E89", "#F8F4EC", "#B22222"))

row_ann <- rowAnnotation(Group = gene_group_factor, col = list(Group = group_cols),
                         show_annotation_name = FALSE, simple_anno_size = unit(4, "mm"))
sc_top_ann <- HeatmapAnnotation(
  State = factor(sc_state_order, levels = sc_state_order),
  col = list(State = sc_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(4, "mm"))
pdo_top_ann <- HeatmapAnnotation(
  State = factor(pdo_state_order, levels = pdo_state_order),
  col = list(State = pdo_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(4, "mm"))

ht_sc <- Heatmap(
  sc_z, name = "Row Z-score",
  top_annotation = sc_top_ann, left_annotation = row_ann, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = gene_group_factor,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"), column_names_rot = 45,
  border = TRUE, column_title = "scAtlas (5 states)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sc_labels[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
  },
  heatmap_legend_param = list(title = "Relative\nexpression",
                              title_gp = gpar(fontsize = 10, fontface = "bold"),
                              labels_gp = gpar(fontsize = 9)),
  width = unit(ncol(sc_z) * 1.5, "cm"))

ht_pdo <- Heatmap(
  pdo_z, name = "Row Z-score (PDO)",
  top_annotation = pdo_top_ann, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = gene_group_factor,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"), column_names_rot = 45,
  border = TRUE, column_title = "PDO (6 states)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(pdo_labels[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
  },
  width = unit(ncol(pdo_z) * 1.5, "cm"))

pdf(file.path(out_dir, "Auto_gene_expression_state_heatmap.pdf"),
    width = 20, height = 8, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1,
                                           heights = unit(c(2, 1), c("cm", "null")))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text("Candidate gene expression across epithelial states: scAtlas vs PDO",
          x = unit(0.5, "npc"), y = unit(0.70, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text(paste0("Colour = row Z-scored expression; cell labels = median of ",
                 "sample-level mean log-normalised expression."),
          x = unit(0.5, "npc"), y = unit(0.25, "npc"), gp = gpar(fontsize = 9))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_sc + ht_pdo, newpage = FALSE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     gap = unit(10, "mm"))
popViewport(2)
dev.off()
message("Saved heatmap.")

####################
# ========== PLOT 2: GENE CORRELATION ==========
####################
message("\n=== Plotting gene correlations ===")
make_gene_plot <- function(df, gene_name, target_state, is_negative_control = FALSE) {
  df$dataset <- factor(df$dataset, levels = c("scAtlas", "PDO"))
  dataset_cols <- c("scAtlas" = "#2C7BB6", "PDO" = "#D7191C")
  base_theme <- theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "none")

  make_panel <- function(plot_data, x_var, y_var, x_lab, y_lab, title) {
    pd <- plot_data %>% filter(is.finite(.data[[x_var]]) & is.finite(.data[[y_var]]))
    if (nrow(pd) < 2) {
      dummy_df <- data.frame(dataset = factor(c("scAtlas", "PDO"), levels = c("scAtlas", "PDO")))
      p <- ggplot(dummy_df) +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No data available\n(",
                                ifelse(all(is.na(plot_data[[x_var]])),
                                       paste0(x_lab, " unavailable"),
                                       paste0(y_lab, " unavailable")), ")"),
                 size = 3.5, color = "grey50", fontface = "italic") +
        facet_wrap(~ dataset, nrow = 2, scales = "free") +
        labs(title = title, x = x_lab, y = y_lab) + base_theme +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.border = element_rect(fill = NA, colour = "grey80"))
      return(p)
    }

    cor_df <- pd %>% group_by(dataset) %>% summarise(
      cor_text = {
        valid <- is.finite(.data[[x_var]]) & is.finite(.data[[y_var]])
        if (sum(valid) < 3) "rho = NA" else {
          p_val <- cor.test(.data[[x_var]][valid], .data[[y_var]][valid],
                            method = "spearman", exact = FALSE)$p.value
          paste0("rho = ", sprintf("%.2f", cor(.data[[x_var]][valid],
                 .data[[y_var]][valid], method = "spearman")), ", ",
                 if (p_val < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p_val)))
        }
      }, .groups = "drop")

    label_pos <- pd %>% group_by(dataset) %>% summarise(
      x_pos = min(.data[[x_var]], na.rm = TRUE) +
        0.02 * diff(range(.data[[x_var]], na.rm = TRUE)),
      y_pos = max(.data[[y_var]], na.rm = TRUE) -
        0.02 * diff(range(.data[[y_var]], na.rm = TRUE)),
      .groups = "drop") %>%
      left_join(cor_df, by = "dataset")

    p <- ggplot(pd, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(aes(fill = dataset), shape = 21, size = 2.5, alpha = 0.7,
                 colour = "grey30", stroke = 0.3) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30",
                  linewidth = 0.7, linetype = "dashed", alpha = 0.2) +
      geom_text(data = label_pos, aes(x = x_pos, y = y_pos, label = cor_text),
                hjust = 0, vjust = 1, size = 3, color = "grey20",
                fontface = "italic", inherit.aes = FALSE) +
      scale_fill_manual(values = dataset_cols) +
      facet_wrap(~ dataset, nrow = 2, scales = "free", drop = FALSE) +
      labs(title = title, x = x_lab, y = y_lab) + base_theme

    missing_datasets <- setdiff(c("scAtlas", "PDO"), unique(as.character(pd$dataset)))
    if (length(missing_datasets) > 0) {
      dummy_missing <- data.frame(
        dataset = factor(missing_datasets, levels = c("scAtlas", "PDO")),
        dummy_x = 0, dummy_y = 0)
      p <- p + geom_text(data = dummy_missing,
                         aes(x = dummy_x, y = dummy_y,
                             label = "No data available\n(regulon unavailable)"),
                         inherit.aes = FALSE, size = 3.5, color = "grey50",
                         fontface = "italic")
    }
    p
  }

  p1 <- make_panel(df, "mean_log_expr", "state_abundance",
                   paste0(gene_name, " mean log-norm. expression"),
                   paste0(target_state, "\nstate abundance (%)"),
                   "Gene expr. vs State abundance")
  p2 <- make_panel(df, "mean_regulon_auc", "state_abundance",
                   paste0(gene_name, " regulon AUCell score"),
                   paste0(target_state, "\nstate abundance (%)"),
                   "Regulon AUC vs State abundance")
  p3 <- make_panel(df, "mean_log_expr", "mean_dge_ucell",
                   paste0(gene_name, " mean log-norm. expression"),
                   paste0(target_state, " DGE markers\nUCell enrichment"),
                   "Gene expr. vs DGE marker enrichment")
  p4 <- make_panel(df, "mean_regulon_auc", "mean_dge_ucell",
                   paste0(gene_name, " regulon AUCell score"),
                   paste0(target_state, " DGE markers\nUCell enrichment"),
                   "Regulon AUC vs DGE marker enrichment")

  main_title <- paste0(gene_name, " \u2014 ", target_state, " state correlation (per sample)")
  if (is_negative_control) main_title <- paste0(main_title, " [Negative Control]")

  (p1 | p2 | p3 | p4) + plot_annotation(
    title = main_title,
    subtitle = paste0("DGE markers (top ", top_n_markers, "): ",
                      paste(dge_gene_sets[[target_state]], collapse = ", ")),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35")))
}

for (primary in names(primary_genes)) {
  ts <- primary_genes[[primary]]
  message("  Plotting: ", primary)
  df_prim <- all_data %>% filter(gene == primary, target_state == ts)
  p_prim <- make_gene_plot(df_prim, primary, ts, is_negative_control = FALSE)
  df_neg <- all_data %>% filter(gene == neg_ctrl, target_state == ts)
  p_neg <- make_gene_plot(df_neg, neg_ctrl, ts, is_negative_control = TRUE)

  out_pdf <- file.path(out_dir, paste0("Auto_", primary, "_expression_correlation.pdf"))
  pdf(out_pdf, width = 20, height = 10, useDingbats = FALSE)
  print(p_prim)
  print(p_neg)
  dev.off()
  message("  Saved: ", out_pdf)
}

####################
# ========== PLOT 3: REGULON-LEVEL STATE CORRELATION ==========
####################
message("\n=== Plotting regulon-level state correlations ===")

## Panel definitions (same 4 types as gene-level, same order)
panel_defs <- list(
  list(y_col = "rho_expr_abundance",
       title = "Gene expr. vs State abundance",
       y_lab = expression(rho * "(TF expression, state abundance)")),
  list(y_col = "rho_auc_abundance",
       title = "Regulon AUC vs State abundance",
       y_lab = expression(rho * "(TF regulon AUC, state abundance)")),
  list(y_col = "rho_expr_dge",
       title = "Gene expr. vs DGE marker enrichment",
       y_lab = expression(rho * "(TF expression, DGE enrichment)")),
  list(y_col = "rho_auc_dge",
       title = "Regulon AUC vs DGE marker enrichment",
       y_lab = expression(rho * "(TF regulon AUC, DGE enrichment)"))
)

## Top-10 selection by combined rank of AUCell activity and positive correlation
get_top10 <- function(df, y_col) {
  df_valid <- df %>% filter(is.finite(mean_state_auc), is.finite(.data[[y_col]]))
  if (nrow(df_valid) < 10) return(df_valid$regulon)
  df_valid %>%
    mutate(rank_auc = rank(-mean_state_auc),
           rank_cor = rank(-.data[[y_col]]),
           combined_rank = rank_auc + rank_cor) %>%
    slice_min(combined_rank, n = 10, with_ties = FALSE) %>%
    pull(regulon)
}

## Build one page (one state pair)
make_regulon_page <- function(sc_df, pdo_df, sp) {
  sc_state  <- sp$sc_state
  pdo_state <- sp$pdo_state
  label     <- sp$label
  dge_markers <- dge_gene_sets[[sc_state]]

  base_theme <- theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "none")

  panel_plots <- list()
  for (pd_idx in seq_along(panel_defs)) {
    pdef  <- panel_defs[[pd_idx]]
    y_col <- pdef$y_col

    ## Top 10 for each dataset
    sc_top10  <- get_top10(sc_df, y_col)
    pdo_top10 <- get_top10(pdo_df, y_col)
    shared    <- intersect(sc_top10, pdo_top10)

    ## Prepare combined plot data
    prep <- function(df, top10, ds_label) {
      df %>%
        filter(is.finite(mean_state_auc), is.finite(.data[[y_col]])) %>%
        mutate(
          dataset    = ds_label,
          is_top10   = regulon %in% top10,
          is_shared  = regulon %in% shared,
          label_text = ifelse(regulon %in% top10, regulon, NA_character_),
          point_cat  = case_when(
            regulon %in% shared ~ "Shared top 10",
            regulon %in% top10  ~ "Top 10",
            TRUE                ~ "Other"
          )
        )
    }

    plot_df <- bind_rows(
      prep(sc_df,  sc_top10,  "scAtlas"),
      prep(pdo_df, pdo_top10, "PDO")
    )
    plot_df$dataset   <- factor(plot_df$dataset, levels = c("scAtlas", "PDO"))
    plot_df$point_cat <- factor(plot_df$point_cat,
                                levels = c("Shared top 10", "Top 10", "Other"))

    p <- ggplot(plot_df, aes(x = mean_state_auc, y = .data[[y_col]])) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_point(aes(color = point_cat, size = point_cat), alpha = 0.7) +
      scale_color_manual(
        values = c("Shared top 10" = "#E41A1C", "Top 10" = "#2C7BB6", "Other" = "grey70"),
        drop = FALSE) +
      scale_size_manual(
        values = c("Shared top 10" = 3, "Top 10" = 2.5, "Other" = 1.2),
        drop = FALSE) +
      geom_text_repel(
        data = . %>% filter(!is.na(label_text)),
        aes(label = label_text, color = point_cat),
        size = 2.5, max.overlaps = 20, show.legend = FALSE,
        segment.color = "grey50", segment.size = 0.3, seed = 42) +
      facet_wrap(~ dataset, nrow = 2, scales = "free") +
      labs(title = pdef$title,
           x = "Mean regulon AUCell in state",
           y = pdef$y_lab) +
      base_theme +
      guides(color = "none", size = "none")

    panel_plots[[pd_idx]] <- p
  }

  ## Compose 4-column layout
  combined <- panel_plots[[1]] | panel_plots[[2]] | panel_plots[[3]] | panel_plots[[4]]

  ## Build shared legend manually
  legend_df <- data.frame(
    x = 1:3, y = 1:3,
    cat = factor(c("Shared top 10", "Top 10", "Other"),
                 levels = c("Shared top 10", "Top 10", "Other")))
  p_legend <- ggplot(legend_df, aes(x, y, color = cat, size = cat)) +
    geom_point() +
    scale_color_manual(
      values = c("Shared top 10" = "#E41A1C", "Top 10" = "#2C7BB6", "Other" = "grey70"),
      name = NULL) +
    scale_size_manual(
      values = c("Shared top 10" = 3, "Top 10" = 2.5, "Other" = 1.2),
      name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           size = "none") +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 9))
  shared_legend <- cowplot::get_legend(p_legend)

  combined + plot_annotation(
    title = paste0(label, " \u2014 Regulon activity vs state correlation"),
    subtitle = paste0("scAtlas: ", sc_state, " | PDO: ", pdo_state,
                      "    DGE markers (top ", top_n_markers, "): ",
                      paste(dge_markers, collapse = ", ")),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35")))
}

## Generate 4-page PDF
reg_pdf <- file.path(out_dir, "Auto_regulon_state_correlation.pdf")
pdf(reg_pdf, width = 22, height = 12, useDingbats = FALSE)
for (sp in state_pairs) {
  sc_df  <- regulon_data %>% filter(dataset == "scAtlas", state_label == sp$label)
  pdo_df <- regulon_data %>% filter(dataset == "PDO",     state_label == sp$label)

  if (nrow(sc_df) == 0 && nrow(pdo_df) == 0) {
    message("  Skipping ", sp$label, " (no data)")
    next
  }

  message("  Plotting: ", sp$label)
  p <- tryCatch(make_regulon_page(sc_df, pdo_df, sp), error = function(e) {
    message("    Error: ", e$message)
    NULL
  })
  if (!is.null(p)) print(p)
}
dev.off()
message("Saved: ", reg_pdf)

message("\nDone.")
