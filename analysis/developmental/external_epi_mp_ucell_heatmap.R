####################
# Analysis registry:
#   Status: active
#   Script: analysis/developmental/external_epi_mp_ucell_heatmap.R
#   Methodology: analysis/methodology/developmental/developmental_reference_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# external_epi_mp_ucell_heatmap.R
# Score filtered scATLAS metaprograms (nMP = 19) with UCell in three
# external epithelial datasets:
#   1. Adult oesophagus
#   2. Adult stomach
#   3. Barretts oesophagus
# Outputs mean per-cell-type UCell score heatmaps and score tables.
# Environment: dmtcp
# Enrichment-facing entry point that mirrors the developmental row structure.
####################

library(Seurat)
library(UCell)
library(Matrix)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(clusterProfiler)
library(BiocParallel)
library(gridExtra)
library(SingleCellExperiment)
library(SummarizedExperiment)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# User-tunable parameters
####################
max_cells_per_type <- 5000L
ucell_cores <- 4L
seed_base <- 1234L

if (!is.finite(max_cells_per_type) || max_cells_per_type <= 0) {
  stop("max_cells_per_type must be a finite positive integer.")
}

out_dir <- "Auto_external_epi_mp_ucell"
cache_dir <- file.path(out_dir, "Auto_cache")
summary_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/updates/new_updates/summaries"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Core inputs
####################
mp_path <- "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
cluster_enrich_path <- "cluster_enrich.rds"
developmental_ref_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage"
adult_stomach_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Stomach/data_9_9_annotated_seurat_all_ut.rds"
barretts_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Barretts/alldatahighquality.rds"
oesophagus_meta_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/metadata/EoE_meta.txt"
oesophagus_cell_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE_cell.tsv"
oesophagus_gene_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE_gene.tsv"
oesophagus_mtx_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE.mtx"

####################
# Load and filter MPs using the standard silhouette rule
####################
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

geneNMF.metaprograms <- readRDS(mp_path)
mp_genes <- geneNMF.metaprograms$metaprograms.genes

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp_genes <- mp_genes[!names(mp_genes) %in% bad_mp_names]
}

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_genes)))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_order <- paste0("MP", mp_tree_order)
mp_order <- mp_order[mp_order %in% names(mp_genes)]
mp_genes <- mp_genes[mp_order]
mp_sizes <- vapply(mp_genes, length, integer(1))

# Format column labels with MP# (n=genes) + Description on second line
col_labels <- sapply(mp_order, function(m) {
  desc <- if (m %in% names(mp_descriptions)) mp_descriptions[m] else ""
  if (nchar(desc) > 0) {
    paste0(m, " (n=", mp_sizes[m], ")\n", desc)
  } else {
    paste0(m, " (n=", mp_sizes[m], ")")
  }
})
names(col_labels) <- mp_order

####################
# Reference row orders copied from the developmental marker workflow
####################
oesophagus_labels <- c(
  "Quiescent basal cell",
  "Basal cell (cycling)",
  "Suprabasal",
  "Apical cell"
)
oesophagus_term_map <- setNames(
  paste0(gsub(" ", "_", oesophagus_labels), "_Oesophagus..Adult_Epithelium"),
  oesophagus_labels
)
oesophagus_order <- unname(oesophagus_term_map)

stomach_rename_map <- c(
  "GKN+F" = "GKN+F_PylG_Stomach..Adult_Epithelium",
  "ADH1+GKN1-F" = "ADH1+GKN1-F_PylG_Stomach..Adult_Epithelium",
  "PG/Neck1" = "PG/Neck1_PylG_Stomach..Adult_Epithelium",
  "PG/Neck2" = "PG/Neck2_PylG_Stomach..Adult_Epithelium",
  "NE1" = "NE1_PylG_Stomach..Adult_Epithelium",
  "PC" = "PC_PylG_Stomach..Adult_Epithelium",
  "Chief" = "Chief_FundG_Stomach..Adult_Epithelium",
  "Pr_epi" = "Pr_epi_FundG_Stomach..Adult_Epithelium",
  "NE2" = "NE2_PylG/IntestMeta_Stomach..Adult_Epithelium",
  "Ent" = "Ent_IntestMeta_Stomach..Adult_Epithelium",
  "Gob" = "Gob_IntestMeta_Stomach..Adult_Epithelium"
)
stomach_order <- c(
  "GKN+F_PylG_Stomach..Adult_Epithelium",
  "ADH1+GKN1-F_PylG_Stomach..Adult_Epithelium",
  "PG/Neck1_PylG_Stomach..Adult_Epithelium",
  "PG/Neck2_PylG_Stomach..Adult_Epithelium",
  "NE1_PylG_Stomach..Adult_Epithelium",
  "PC_PylG_Stomach..Adult_Epithelium",
  "Chief_FundG_Stomach..Adult_Epithelium",
  "Pr_epi_FundG_Stomach..Adult_Epithelium",
  "NE2_PylG/IntestMeta_Stomach..Adult_Epithelium",
  "Ent_IntestMeta_Stomach..Adult_Epithelium",
  "Gob_IntestMeta_Stomach..Adult_Epithelium"
)

barrett_term_map <- c(
  "Basal" = "Basal_Normal_Esophagus..Barretts_Oesophagus",
  "Suprabasal" = "Suprabasal_Normal_Esophagus..Barretts_Oesophagus",
  "Suprabasal_Dividing" = "Suprabasal_Dividing_Normal_Esophagus..Barretts_Oesophagus",
  "Intermediate" = "Intermediate_Normal_Esophagus..Barretts_Oesophagus",
  "Superficial" = "Superficial_Normal_Esophagus..Barretts_Oesophagus",
  "Undifferentiated" = "Undifferentiated_Normal_Gastric..Barretts_Oesophagus",
  "Undifferentiated_Dividing" = "Undifferentiated_Dividing_Normal_Gastric..Barretts_Oesophagus",
  "Foveolar_Intermediate" = "Foveolar_Intermediate_Normal_Gastric..Barretts_Oesophagus",
  "Foveolar_differentiated" = "Foveolar_differentiated_Normal_Gastric..Barretts_Oesophagus",
  "Chief" = "Chief_Normal_Gastric..Barretts_Oesophagus",
  "Parietal" = "Parietal_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_GHRL" = "Endocrine_GHRL_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_CHGA" = "Endocrine_CHGA_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_NEUROD1" = "Endocrine_NEUROD1_Normal_Gastric..Barretts_Oesophagus",
  "Columnar_Undifferentiated" = "Columnar_Undifferentiated_Barretts_Esophagus..Barretts_Oesophagus",
  "Columnar_Undifferentiated_Dividing" = "Columnar_Dividing_Barretts_Esophagus..Barretts_Oesophagus",
  "Endocrine_NEUROG3" = "Endocrine_NEUROG3_Barretts_Esophagus..Barretts_Oesophagus",
  "C1" = "Columnar_Intermediate_Barretts_Esophagus..Barretts_Oesophagus",
  "C2" = "Columnar_differentiated_Barretts_Esophagus..Barretts_Oesophagus",
  "Goblet" = "Goblet_Barretts_Esophagus..Barretts_Oesophagus",
  "Duct_Intercalating" = "Duct_Intercalating_Submucosal_Glands..Barretts_Oesophagus",
  "Oncocytes_MUC5B_Low" = "Oncocytes_Submucosal_Glands..Barretts_Oesophagus",
  "Mucous_MUC5B_High" = "Mucous_Submucosal_Glands..Barretts_Oesophagus"
)
barrett_order <- unname(barrett_term_map)

combined_order <- c(oesophagus_order, stomach_order, barrett_order)

####################
# Helper functions
####################
heatmap_cols <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

sample_meta_by_term <- function(meta_df, term_levels, max_cells, seed) {
  meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
  meta_df$term <- factor(meta_df$term, levels = term_levels)
  meta_df <- meta_df[!is.na(meta_df$term), , drop = FALSE]

  set.seed(seed)
  split_meta <- split(meta_df, meta_df$term, drop = TRUE)
  sampled_list <- lapply(split_meta, function(df) {
    if (nrow(df) <= max_cells) {
      return(df)
    }
    df[sample(seq_len(nrow(df)), max_cells, replace = FALSE), , drop = FALSE]
  })

  sampled_meta <- bind_rows(sampled_list)
  sampled_meta$term <- factor(sampled_meta$term, levels = term_levels)
  sampled_meta <- sampled_meta[order(sampled_meta$term, sampled_meta$cell), , drop = FALSE]
  sampled_meta$term <- as.character(sampled_meta$term)
  sampled_meta
}

subset_oesophagus_counts <- function(sampled_meta, cache_tag) {
  cache_path <- file.path(cache_dir, paste0("Auto_oesophagus_subset_", cache_tag, ".rds"))
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  cat("Building cached adult oesophagus epithelial subset:", cache_tag, "\n")

  all_cells <- scan(oesophagus_cell_path, what = character(), quiet = TRUE)
  all_genes <- scan(oesophagus_gene_path, what = character(), quiet = TRUE)
  keep_idx <- match(sampled_meta$cell, all_cells)

  if (anyNA(keep_idx)) {
    missing_cells <- sampled_meta$cell[is.na(keep_idx)]
    stop("Adult_Oesophagus cell barcode mismatch: ", paste(head(missing_cells, 10), collapse = ", "))
  }

  map_df <- data.frame(
    old_idx = keep_idx,
    new_idx = seq_along(keep_idx),
    stringsAsFactors = FALSE
  )
  map_df <- map_df[order(map_df$old_idx), , drop = FALSE]

  map_path <- file.path(cache_dir, paste0("Auto_oesophagus_subset_", cache_tag, "_map.tsv"))
  awk_path <- file.path(cache_dir, paste0("Auto_oesophagus_subset_", cache_tag, "_filter.awk"))
  triplet_path <- file.path(cache_dir, paste0("Auto_oesophagus_subset_", cache_tag, "_triplets.tsv"))
  filtered_mtx_path <- file.path(cache_dir, paste0("Auto_oesophagus_subset_", cache_tag, ".mtx"))

  write.table(
    map_df,
    file = map_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  writeLines(
    c(
      "NR == FNR { map[$1] = $2; next }",
      "FNR <= 2 { next }",
      "($2 in map) { print $1 \" \" map[$2] \" \" $3 }"
    ),
    con = awk_path
  )

  awk_status <- system2(
    "awk",
    args = c("-f", awk_path, map_path, oesophagus_mtx_path),
    stdout = triplet_path
  )
  if (!identical(awk_status, 0L)) {
    stop("Failed to stream-filter the adult oesophagus matrix.")
  }

  nnz_line <- system2("wc", args = c("-l", triplet_path), stdout = TRUE)
  nnz_keep <- as.integer(strsplit(trimws(nnz_line), "\\s+")[[1]][1])
  if (is.na(nnz_keep) || nnz_keep <= 0) {
    stop("No non-zero entries found after adult oesophagus epithelial subsetting.")
  }

  writeLines(
    c(
      "%%MatrixMarket matrix coordinate integer general",
      paste(length(all_genes), nrow(sampled_meta), nnz_keep)
    ),
    con = filtered_mtx_path
  )
  append_ok <- file.append(filtered_mtx_path, triplet_path)
  if (!isTRUE(append_ok)) {
    stop("Failed to append filtered adult oesophagus triplets to Matrix Market header.")
  }

  subset_counts <- readMM(filtered_mtx_path)
  subset_counts <- as(subset_counts, "CsparseMatrix")
  rownames(subset_counts) <- all_genes
  colnames(subset_counts) <- sampled_meta$cell

  cache_obj <- list(
    counts = subset_counts,
    sampled_meta = sampled_meta
  )
  saveRDS(cache_obj, cache_path)

  unlink(c(map_path, awk_path, triplet_path, filtered_mtx_path), force = TRUE)
  cache_obj
}

score_sampled_counts <- function(counts_mat, sampled_meta, full_meta, dataset_source, seed) {
  sampled_meta <- as.data.frame(sampled_meta, stringsAsFactors = FALSE)
  full_meta <- as.data.frame(full_meta, stringsAsFactors = FALSE)

  if (!all(sampled_meta$cell %in% colnames(counts_mat))) {
    stop("Counts matrix is missing sampled cells for ", dataset_source)
  }

  counts_mat <- counts_mat[, sampled_meta$cell, drop = FALSE]
  sampled_meta <- sampled_meta[match(colnames(counts_mat), sampled_meta$cell), , drop = FALSE]
  rownames(sampled_meta) <- sampled_meta$cell

  usable_mp_genes <- lapply(mp_genes, function(genes) {
    intersect(genes, rownames(counts_mat))
  })
  empty_mps <- names(usable_mp_genes)[lengths(usable_mp_genes) == 0]
  if (length(empty_mps) > 0) {
    stop("No overlapping genes for ", dataset_source, " in: ", paste(empty_mps, collapse = ", "))
  }

  set.seed(seed)
  bp_param <- BiocParallel::SerialParam(progressbar = FALSE)
  obj <- CreateSeuratObject(counts = counts_mat, meta.data = sampled_meta)
  obj <- AddModuleScore_UCell(
    obj,
    features = usable_mp_genes,
    slot = "counts",
    BPPARAM = bp_param,
    ncores = 1,
    name = "_UCell"
  )

  ucell_score_cols <- paste0(mp_order, "_UCell")
  if (!all(ucell_score_cols %in% colnames(obj@meta.data))) {
    missing_cols <- ucell_score_cols[!ucell_score_cols %in% colnames(obj@meta.data)]
    stop("UCell score columns missing for ", dataset_source, ": ", paste(missing_cols, collapse = ", "))
  }

  ucell_scores <- obj@meta.data[, ucell_score_cols, drop = FALSE]
  colnames(ucell_scores) <- mp_order
  score_df <- cbind(
    sampled_meta[, c("cell", "term", "dataset_source"), drop = FALSE],
    as.data.frame(ucell_scores, stringsAsFactors = FALSE)
  )

  mean_scores <- score_df %>%
    group_by(term) %>%
    summarise(across(all_of(mp_order), mean), .groups = "drop")

  total_counts <- full_meta %>%
    dplyr::count(term, name = "n_cells_total")
  scored_counts <- sampled_meta %>%
    dplyr::count(term, name = "n_cells_scored")

  summary_df <- full_join(total_counts, scored_counts, by = "term") %>%
    left_join(mean_scores, by = "term") %>%
    mutate(
      dataset_source = dataset_source,
      n_cells_total = ifelse(is.na(n_cells_total), 0L, n_cells_total),
      n_cells_scored = ifelse(is.na(n_cells_scored), 0L, n_cells_scored)
    ) %>%
    select(dataset_source, term, n_cells_total, n_cells_scored, all_of(mp_order))

  list(summary = summary_df, per_cell_scores = score_df)
}

build_custom_annotation_heatmap <- function(cluster_enrich, custom_ref, row_order, cap = 7) {
  ordered_mps <- mp_order[mp_order %in% names(cluster_enrich)]
  terms_use <- as.character(custom_ref$TERM2NAME$term)

  df_list <- lapply(ordered_mps, function(prog) {
    er <- cluster_enrich[[prog]][[custom_ref$ref_name]]
    if (is.null(er)) {
      return(NULL)
    }

    res_tbl <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(res_tbl) || nrow(res_tbl) == 0) {
      return(NULL)
    }

    term_col <- if ("Description" %in% colnames(res_tbl)) {
      res_tbl$Description
    } else {
      res_tbl$ID
    }

    data.frame(
      Program = prog,
      Term = term_col,
      padj = res_tbl$p.adjust,
      Overlap = res_tbl$GeneRatio,
      stringsAsFactors = FALSE
    )
  })

  df <- bind_rows(df_list)
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  final_df <- full_grid %>%
    left_join(df, by = c("Term", "Program")) %>%
    mutate(
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = tidyr::replace_na(Overlap, "")
    )

  mat <- final_df %>%
    select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame()
  rownames(mat) <- mat$Term
  mat <- as.matrix(mat[, ordered_mps, drop = FALSE])

  text_mat <- final_df %>%
    select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame()
  rownames(text_mat) <- text_mat$Term
  text_mat <- as.matrix(text_mat[, ordered_mps, drop = FALSE])

  missing_rows <- setdiff(row_order, rownames(mat))
  if (length(missing_rows) > 0) {
    stop("Missing annotation rows: ", paste(missing_rows, collapse = ", "))
  }

  mat <- matrix(
    as.numeric(mat[row_order, ordered_mps, drop = FALSE]),
    nrow = length(row_order),
    ncol = length(ordered_mps),
    dimnames = list(row_order, ordered_mps)
  )
  text_mat <- text_mat[row_order, ordered_mps, drop = FALSE]

  list(
    mat = mat,
    text_mat = text_mat,
    legend_breaks = c(0, 2, 4, 6),
    legend_labels = c("0", "2", "4", "6+"),
    breaks = seq(0, cap, length.out = length(heatmap_cols) + 1)
  )
}


make_heatmap_obj <- function(score_mat,
                             display_mat,
                             plot_title,
                             gaps_row_idx = NULL,
                             col_fun = NULL,
                             legend_title = "") {
  # If gaps_row is given as an index (e.g. c(4)), we will use row_split
  row_split <- NULL
  if (!is.null(gaps_row_idx) && length(gaps_row_idx) > 0) {
    row_split <- rep(1:(length(gaps_row_idx)+1), times = diff(c(0, gaps_row_idx, nrow(score_mat))))
  }

  Heatmap(score_mat,
          name = legend_title,
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = row_split,
          row_title = NULL,
          column_title = plot_title,
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_rot = 25,
          column_labels = col_labels[colnames(score_mat)],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 9),
          cell_fun = function(j, i, x, y, width, height, fill) {
             if(display_mat[i, j] != "") {
                 grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
             }
          })
}

plot_ucell_heatmap <- function(score_mat, file_stub, plot_title, gaps_row = NULL) {
  ucell_cap <- max(score_mat, na.rm = TRUE)
  col_fun <- colorRamp2(seq(0, ucell_cap, length.out = length(heatmap_cols)), heatmap_cols)
  ucell_display <- ifelse(is.na(score_mat), "", sprintf("%.2f", score_mat))
  
  if (is.null(gaps_row) || length(gaps_row) == 0) {
     row_split <- NULL
  } else {
     row_split <- rep(1:(length(gaps_row)+1), times = diff(c(0, gaps_row, nrow(score_mat))))
  }

  ht <- Heatmap(score_mat,
          name = "Mean UCell",
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = row_split,
          row_title = NULL,
          column_title = plot_title,
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_rot = 25,
          column_labels = col_labels[colnames(score_mat)],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 9),
          cell_fun = function(j, i, x, y, width, height, fill) {
             if(ucell_display[i, j] != "") {
                 grid.text(ucell_display[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
             }
          })

  pdf(file.path(out_dir, paste0(file_stub, ".pdf")), width = 18, height = max(8, nrow(score_mat) * 0.28))
  draw(ht)
  dev.off()
}

plot_annotation_vs_scoring_page <- function(annotation_obj,
                                            scoring_mat,
                                            file_path,
                                            page_title,
                                            gaps_row = NULL) {
  scoring_cap <- max(scoring_mat, na.rm = TRUE)
  col_fun_score <- colorRamp2(seq(0, scoring_cap, length.out = length(heatmap_cols)), heatmap_cols)
  scoring_display <- ifelse(is.na(scoring_mat), "", sprintf("%.2f", scoring_mat))
  
  col_fun_anno <- colorRamp2(seq(0, max(annotation_obj$mat, na.rm=TRUE), length.out=length(heatmap_cols)), heatmap_cols)

  if (is.null(gaps_row) || length(gaps_row) == 0) {
     row_split <- NULL
  } else {
     row_split <- rep(1:(length(gaps_row)+1), times = diff(c(0, gaps_row, nrow(scoring_mat))))
  }

  ht_anno <- Heatmap(annotation_obj$mat,
          name = "Annotation\n-log10(padj)",
          col = col_fun_anno,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = row_split,
          row_title = NULL,
          show_row_names = FALSE,
          column_title = "Enrichment Annotation",
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_rot = 25,
          column_labels = col_labels[colnames(annotation_obj$mat)],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 9),
          cell_fun = function(j, i, x, y, width, height, fill) {
             if(annotation_obj$text_mat[i, j] != "") {
                 grid.text(annotation_obj$text_mat[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
             }
          })

  ht_score <- Heatmap(scoring_mat,
          name = "Reference\nScore",
          col = col_fun_score,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = row_split,
          row_title = NULL,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_max_width = max_text_width(rownames(scoring_mat), gp = gpar(fontsize = 8)),
          column_title = "Reference Scoring (mean UCell)",
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_rot = 25,
          column_labels = col_labels[colnames(scoring_mat)],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 9),
          cell_fun = function(j, i, x, y, width, height, fill) {
             if(scoring_display[i, j] != "") {
                 grid.text(scoring_display[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
             }
          })

  ht_list <- ht_anno + ht_score
  draw(ht_list, column_title = page_title, column_title_gp = gpar(fontsize = 14, fontface = "bold"))
}


####################
# Prepare and score adult oesophagus
####################
cat("Preparing adult oesophagus metadata...\n")
oesophagus_meta <- fread(oesophagus_meta_path)
oesophagus_meta <- oesophagus_meta[NAME != "TYPE"]
oesophagus_meta <- oesophagus_meta %>%
  transmute(
    cell = NAME,
    term = unname(oesophagus_term_map[cell_type_anno]),
    dataset_source = "Adult_Oesophagus"
  ) %>%
  filter(!is.na(term))

oesophagus_sampled <- sample_meta_by_term(
  meta_df = oesophagus_meta,
  term_levels = oesophagus_order,
  max_cells = max_cells_per_type,
  seed = seed_base + 1L
)

oesophagus_cache_tag <- paste0(
  "max", max_cells_per_type,
  "_seed", seed_base + 1L
)
oesophagus_subset <- subset_oesophagus_counts(
  sampled_meta = oesophagus_sampled,
  cache_tag = oesophagus_cache_tag
)

cat("Scoring adult oesophagus MPs with UCell...\n")
oesophagus_scored <- score_sampled_counts(
  counts_mat = oesophagus_subset$counts,
  sampled_meta = oesophagus_subset$sampled_meta,
  full_meta = oesophagus_meta,
  dataset_source = "Adult_Oesophagus",
  seed = seed_base + 11L
)

rm(oesophagus_subset)
gc()

####################
# Prepare and score adult stomach
####################
cat("Preparing adult stomach epithelial subset...\n")
adult_stomach <- readRDS(adult_stomach_path)
stomach_cells <- Cells(adult_stomach)
stomach_meta <- adult_stomach@meta.data
stomach_meta$cell <- stomach_cells
stomach_meta <- stomach_meta %>%
  filter(major_clusters == "epi") %>%
  transmute(
    cell = cell,
    term = unname(stomach_rename_map[as.character(subcluster.v2)]),
    dataset_source = "Adult_Stomach"
  ) %>%
  filter(!is.na(term))

stomach_sampled <- sample_meta_by_term(
  meta_df = stomach_meta,
  term_levels = stomach_order,
  max_cells = max_cells_per_type,
  seed = seed_base + 2L
)

stomach_counts <- GetAssayData(adult_stomach, layer = "counts")
stomach_counts <- stomach_counts[, stomach_sampled$cell, drop = FALSE]

cat("Scoring adult stomach MPs with UCell...\n")
stomach_scored <- score_sampled_counts(
  counts_mat = stomach_counts,
  sampled_meta = stomach_sampled,
  full_meta = stomach_meta,
  dataset_source = "Adult_Stomach",
  seed = seed_base + 12L
)

rm(adult_stomach, stomach_counts)
gc()

####################
# Prepare and score Barretts epithelial subset
####################
cat("Preparing Barretts epithelial subset...\n")
barretts_sce <- readRDS(barretts_path)
barretts_assay_name <- if ("counts" %in% assayNames(barretts_sce)) {
  "counts"
} else {
  assayNames(barretts_sce)[1]
}
barretts_meta <- as.data.frame(colData(barretts_sce), stringsAsFactors = FALSE)
barretts_counts <- assay(barretts_sce, barretts_assay_name)

if (ncol(barretts_counts) != nrow(barretts_meta)) {
  stop("Barretts counts/metadata column mismatch.")
}

barretts_cell_ids <- colnames(barretts_sce)
if (is.null(barretts_cell_ids)) {
  barretts_cell_ids <- colnames(barretts_counts)
}
if (is.null(barretts_cell_ids) || anyNA(barretts_cell_ids) || any(barretts_cell_ids == "") || anyDuplicated(barretts_cell_ids)) {
  barretts_cell_ids <- paste0("Barretts_cell_", seq_len(nrow(barretts_meta)))
}

barretts_meta$cell <- barretts_cell_ids
colnames(barretts_counts) <- barretts_cell_ids
barretts_meta <- barretts_meta %>%
  transmute(
    cell = cell,
    term = unname(barrett_term_map[as.character(cell_type_secondary)]),
    dataset_source = "Barretts"
  ) %>%
  filter(!is.na(term))

barretts_sampled <- sample_meta_by_term(
  meta_df = barretts_meta,
  term_levels = barrett_order,
  max_cells = max_cells_per_type,
  seed = seed_base + 3L
)
barretts_counts <- barretts_counts[, barretts_sampled$cell, drop = FALSE]

cat("Scoring Barretts MPs with UCell...\n")
barretts_scored <- score_sampled_counts(
  counts_mat = barretts_counts,
  sampled_meta = barretts_sampled,
  full_meta = barretts_meta,
  dataset_source = "Barretts",
  seed = seed_base + 13L
)

rm(barretts_sce, barretts_counts)
gc()

####################
# Combine summaries and export
####################
summary_df <- bind_rows(
  oesophagus_scored$summary,
  stomach_scored$summary,
  barretts_scored$summary
)
summary_df$term <- factor(summary_df$term, levels = combined_order)
summary_df <- summary_df %>%
  arrange(term) %>%
  mutate(term = as.character(term))

score_mat <- matrix(
  NA_real_,
  nrow = length(combined_order),
  ncol = length(mp_order),
  dimnames = list(combined_order, mp_order)
)
present_terms <- intersect(summary_df$term, rownames(score_mat))
score_mat[present_terms, mp_order] <- as.matrix(
  summary_df[match(present_terms, summary_df$term), mp_order, drop = FALSE]
)

adult_mat <- score_mat[c(oesophagus_order, stomach_order), mp_order, drop = FALSE]
barretts_mat <- score_mat[barrett_order, mp_order, drop = FALSE]

####################
# Annotation matrices aligned to scoring row order
####################
cluster_enrich <- readRDS(cluster_enrich_path)
adult_ref <- readRDS(file.path(developmental_ref_dir, "enrich_dev_Adult_Epithelium.rds"))
adult_ref$ref_name <- "Adult_Epithelium"
barrett_ref <- readRDS(file.path(developmental_ref_dir, "enrich_dev_Barretts_Oesophagus.rds"))
barrett_ref$ref_name <- "Barretts_Oesophagus"

adult_annotation <- build_custom_annotation_heatmap(
  cluster_enrich = cluster_enrich,
  custom_ref = adult_ref,
  row_order = rownames(adult_mat)
)
barrett_annotation <- build_custom_annotation_heatmap(
  cluster_enrich = cluster_enrich,
  custom_ref = barrett_ref,
  row_order = rownames(barretts_mat)
)

if (!identical(rownames(adult_annotation$mat), rownames(adult_mat))) {
  stop("Adult annotation rows do not match adult scoring rows.")
}
if (!identical(rownames(barrett_annotation$mat), rownames(barretts_mat))) {
  stop("Barretts annotation rows do not match Barretts scoring rows.")
}

write.csv(
  summary_df,
  file = file.path(out_dir, "Auto_external_epi_mp_ucell_summary.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file = file.path(summary_dir, "Auto_external_epi_mp_ucell_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    summary = summary_df,
    combined_matrix = score_mat,
    adult_matrix = adult_mat,
    barretts_matrix = barretts_mat,
    mp_order = mp_order,
    mp_gene_sizes = mp_sizes,
    parameters = list(
      max_cells_per_type = max_cells_per_type,
      ucell_cores = ucell_cores,
      seed_base = seed_base
    )
  ),
  file = file.path(out_dir, "Auto_external_epi_mp_ucell_summary.rds")
)

####################
# Heatmaps
####################
plot_ucell_heatmap(
  score_mat = adult_mat,
  file_stub = "Auto_external_epi_mp_ucell_heatmap_adult_epithelium",
  plot_title = "Adult Epithelium: mean MP UCell score per cell type",
  gaps_row = c(length(oesophagus_order))
)

plot_ucell_heatmap(
  score_mat = barretts_mat,
  file_stub = "Auto_external_epi_mp_ucell_heatmap_barretts_oesophagus",
  plot_title = "Barretts Oesophagus: mean MP UCell score per cell type"
)

plot_ucell_heatmap(
  score_mat = score_mat,
  file_stub = "Auto_external_epi_mp_ucell_heatmap_combined",
  plot_title = "External epithelial datasets: mean MP UCell score per cell type",
  gaps_row = c(
    length(oesophagus_order),
    length(oesophagus_order) + length(stomach_order)
  )
)

####################
# Side-by-side comparison PDF
####################
combined_pdf <- file.path(out_dir, "Auto_external_epi_mp_ucell_annotation_vs_reference.pdf")
pdf(combined_pdf, width = 28, height = 12, onefile = TRUE)
plot_annotation_vs_scoring_page(
  annotation_obj = adult_annotation,
  scoring_mat = adult_mat,
  file_path = combined_pdf,
  page_title = "Adult Epithelium: enrichment annotation vs reference MP scoring",
  gaps_row = c(length(oesophagus_order))
)

plot_annotation_vs_scoring_page(
  annotation_obj = barrett_annotation,
  scoring_mat = barretts_mat,
  file_path = combined_pdf,
  page_title = "Barretts Oesophagus: enrichment annotation vs reference MP scoring"
)
dev.off()

cat("Saved UCell summary tables and heatmaps to:", file.path(getwd(), out_dir), "\n")
cat("Saved lightweight summary CSV to:", file.path(summary_dir, "Auto_external_epi_mp_ucell_summary.csv"), "\n")
cat("Saved annotation-vs-scoring comparison PDF to:", file.path(getwd(), combined_pdf), "\n")
