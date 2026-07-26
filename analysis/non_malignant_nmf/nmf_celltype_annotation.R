####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/nmf_celltype_annotation.R
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
celltype_arg <- args[1]
if (is.na(celltype_arg) || celltype_arg == "") stop("Usage: Rscript nmf_celltype_annotation.R <celltype>")

base_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"

# ---- Cell type → NMF folder & Excel sheet mapping ----
ct_map <- list(
  macrophage  = list(folder = "nmf_macrophage",  sheet = "Macrophages"),
  fibroblast  = list(folder = "nmf_fibroblast",  sheet = "Fibroblasts"),
  endothelial = list(folder = "nmf_endothelial", sheet = "Endothelial"),
  plasma      = list(folder = "nmf_plasma",      sheet = "B_cells"),
  cd4         = list(folder = "nmf_cd4",         sheet = "CD4"),
  cd8         = list(folder = "nmf_cd8",         sheet = "CD8"),
  nk          = list(folder = "nmf_nk",          sheet = NA_character_)
)

ct_info <- ct_map[[celltype_arg]]
if (is.null(ct_info)) stop("Unknown celltype: ", celltype_arg,
                           ". Options: ", paste(names(ct_map), collapse = ", "))

out_dir <- file.path(base_dir, ct_info$folder)
if (!dir.exists(out_dir)) stop("NMF output folder not found: ", out_dir)

# ---- Load MP results ----
mp_path <- file.path(out_dir, "MP_outs_default.rds")
if (!file.exists(mp_path)) stop("MP_outs_default.rds not found at: ", mp_path)

geneNMF.metaprograms <- readRDS(mp_path)
mp_gene_lists  <- geneNMF.metaprograms$metaprograms.genes
mp_assignments <- geneNMF.metaprograms$programs.clusters

cat("Loaded", length(mp_gene_lists), "metaprograms for", celltype_arg, "\n")

# ===========================================================
# 1.  Hallmark reference
# ===========================================================
hallmark_sets     <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

# ===========================================================
# 2.  3CA pan-cancer MP reference (from Immune_NMFs.xlsx)
# ===========================================================
has_3ca <- !is.na(ct_info$sheet)

if (has_3ca) {
  xlsx_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/Immune_NMFs.xlsx"
  mp_df   <- read_excel(xlsx_path, sheet = ct_info$sheet)
  MP_list <- as.list(mp_df)
  MP_list <- lapply(MP_list, function(x) x[!is.na(x)])  # drop NA padding

  mp_term2gene <- data.frame(
    term = rep(names(MP_list), lengths(MP_list)),
    gene = unlist(MP_list),
    row.names = NULL
  )
  mp_term2gene$term <- paste0("3CA_", mp_term2gene$term)
  mp_term2name <- data.frame(
    term = unique(mp_term2gene$term),
    name = unique(mp_term2gene$term)
  )
  cat("Loaded", length(MP_list), "3CA reference MPs from sheet:", ct_info$sheet, "\n")
}

# ===========================================================
# 3.  Enrichment per metaprogram
# ===========================================================
cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {

  genes  <- mp_gene_lists[[mp_name]]
  mp_id  <- as.numeric(gsub("\\D", "", mp_name))
  members <- names(mp_assignments)[mp_assignments == mp_id]

  message(paste0("[", celltype_arg, "] Processing MP: ", mp_name,
                 "  (", length(genes), " genes)"))

  # GO Biological Process
  res_GO <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = "BP", qvalueCutoff = 0.05, readable = TRUE)

  # Hallmark
  res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene,
                    TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)

  results <- list(
    rep_prog = mp_name,
    members  = members,
    genes    = genes,
    GO       = res_GO,
    Hallmark = res_H
  )

  # 3CA pan-cancer MPs (only if sheet exists)
  if (has_3ca) {
    res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene,
                      TERM2NAME = mp_term2name, qvalueCutoff = 0.05)
    results$MPs_3CA <- res_M
  }

  return(results)
})
names(cluster_enrich) <- names(mp_gene_lists)

# Save full enrichment object
saveRDS(cluster_enrich, file.path(out_dir, paste0("cluster_enrich_", celltype_arg, ".rds")))
cat("Enrichment results saved.\n")

# ===========================================================
# 4.  Tree ordering for heatmap columns
# ===========================================================
tree_order       <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order    <- unique(ordered_clusters)

# ===========================================================
# 5.  Heatmap function (adapted from example_anno.R)
# ===========================================================
enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7,
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {

  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)

    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)

    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    if (nrow(r_sig) == 0) return(NULL)

    term <- if ("Description" %in% colnames(r_sig)) r_sig$Description else r_sig$ID
    data.frame(
      Program = prog,
      Term    = term,
      padj    = r_sig$p.adjust,
      Overlap = r_sig$GeneRatio,
      stringsAsFactors = FALSE
    )
  })

  df <- dplyr::bind_rows(df_list)

  if (is.null(df) || nrow(df) == 0) {
    message("No significant results found for: ", element)
    return(invisible(NULL))
  }

  # Select top terms per program
  terms_use <- df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(Program, padj) %>%
    dplyr::group_by(Program) %>%
    dplyr::slice_head(n = top_per_program) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Term) %>%
    dplyr::pull(Term)

  # Cap total terms
  if (length(terms_use) > top_n) {
    terms_use <- df %>%
      dplyr::filter(Term %in% terms_use) %>%
      dplyr::group_by(Term) %>%
      dplyr::summarise(min_p = min(padj), .groups = "drop") %>%
      dplyr::arrange(min_p) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(Term)
  }

  # Column order follows dendrogram tree
  ordered_mps <- paste0("MP", mp_tree_order)
  ordered_mps <- ordered_mps[ordered_mps %in% names(cluster_enrich)]

  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps,
                           stringsAsFactors = FALSE)

  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      score        = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = tidyr::replace_na(Overlap, "")
    )

  # Score matrix
  mat <- final_df %>%
    dplyr::select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>%
    dplyr::select(-Term) %>% as.matrix()

  # Display text matrix
  text_mat <- final_df %>%
    dplyr::select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>%
    dplyr::select(-Term) %>% as.matrix()

  # Align
  mat      <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]

  # Column labels with gene counts
  mp_sizes   <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
  col_labels <- paste0(colnames(mat), "\nn=", mp_sizes)

  # Row ordering by strongest MP assignment
  best_mp   <- colnames(mat)[max.col(mat, ties.method = "first")]
  row_order <- order(match(best_mp, colnames(mat)), -rowSums(mat))
  mat       <- mat[row_order, , drop = FALSE]
  text_mat  <- text_mat[row_order, , drop = FALSE]

  groups   <- colnames(mat)[max.col(mat, ties.method = "first")]
  row_gaps <- which(groups[-length(groups)] != groups[-1])

  # Render
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  pheatmap::pheatmap(mat,
    display_numbers  = text_mat,
    number_color     = "black",
    fontsize_number  = fontsize_row * 1.1,
    labels_col       = col_labels,
    color            = cols,
    breaks           = breaks,
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    gaps_row         = row_gaps,
    border_color     = NA,
    show_colnames    = TRUE,
    angle_col        = 0,
    fontsize_row     = fontsize_row,
    fontsize_col     = fontsize_col,
    main = paste0(celltype_arg, " — ", element, " Enrichment (-log10 padj)")
  )

  return(invisible(mat))
}

# ===========================================================
# 6.  Save heatmaps — 3 PNGs + 1 combined PDF
# ===========================================================
cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

# --- PNG: Hallmark ---
png(file.path(out_dir, paste0("enrich_Hallmark_", celltype_arg, ".png")),
    width = 2000, height = 1750, res = 300)
enrich_heatmap(cluster_enrich, "Hallmark",
               top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()
cat("Saved Hallmark PNG\n")

# --- PNG: GO ---
png(file.path(out_dir, paste0("enrich_GO_", celltype_arg, ".png")),
    width = 2300, height = 2000, res = 300)
enrich_heatmap(cluster_enrich, "GO",
               top_per_program = 6, top_n = 60, cols = cols_palette)
dev.off()
cat("Saved GO PNG\n")

# --- PNG: MPs_3CA (only if sheet available) ---
if (has_3ca) {
  png(file.path(out_dir, paste0("enrich_MPs_3CA_", celltype_arg, ".png")),
      width = 2000, height = 1800, res = 300)
  enrich_heatmap(cluster_enrich, "MPs_3CA",
                 top_per_program = 8, top_n = 80, cols = cols_palette)
  dev.off()
  cat("Saved MPs_3CA PNG\n")
}

# --- Combined PDF ---
pdf(file.path(out_dir, paste0("enrich_combined_", celltype_arg, ".pdf")),
    width = 8, height = 6)
enrich_heatmap(cluster_enrich, "Hallmark",
               top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "GO",
               top_per_program = 6, top_n = 60, cols = cols_palette)
if (has_3ca) {
  enrich_heatmap(cluster_enrich, "MPs_3CA",
                 top_per_program = 8, top_n = 80, cols = cols_palette)
}
dev.off()
cat("Saved combined PDF\n")

cat("\n=== Annotation complete for", celltype_arg, "===\n")

######################################################
# 
# library(grid)
# 
# pdf_out <- file.path(getwd(), "MP_metrics_summary.pdf")
# pdf(pdf_out, width = 8.5, height = 4)  # portrait
# 
# for (ct in names(ct_map)) {
#   folder   <- ct_map[[ct]]$folder
#   rds_path <- file.path(getwd(), folder, "MP_outs_default.rds")
#   
#   grid.newpage()
#   
#   # ---- Load object ----
#   if (!file.exists(rds_path)) {
#     grid.text(
#       paste0("Cell type: ", ct,
#              "\nnumber of samples: NA",
#              "\n\nERROR: MP_outs_default.rds not found:\n", rds_path),
#       x = 0.05, y = 0.97, just = c("left", "top"),
#       gp = gpar(fontsize = 14, fontface = "bold")
#     )
#     next
#   }
#   
#   mp <- tryCatch(readRDS(rds_path), error = function(e) e)
#   if (inherits(mp, "error")) {
#     grid.text(
#       paste0("Cell type: ", ct,
#              "\nnumber of samples: NA",
#              "\n\nERROR reading RDS:\n", conditionMessage(mp)),
#       x = 0.05, y = 0.97, just = c("left", "top"),
#       gp = gpar(fontsize = 14, fontface = "bold")
#     )
#     next
#   }
#   
#   n_samples <- if (!is.null(mp$metaprograms.composition)) ncol(mp$metaprograms.composition) else NA_integer_
#   
#   # ---- Header ----
#   header <- paste0("Cell type: ", ct, "\nnumber of samples: ", n_samples)
#   grid.text(header, x = 0.05, y = 0.97, just = c("left", "top"),
#             gp = gpar(fontsize = 16, fontface = "bold"))
#   
#   metrics <- mp$metaprograms.metrics
#   if (is.null(metrics)) {
#     grid.text("metaprograms.metrics is NULL / not found",
#               x = 0.05, y = 0.90, just = c("left", "top"),
#               gp = gpar(fontsize = 12))
#     next
#   }
#   
#   # ---- Build table with rownames ----
#   metrics_df <- as.data.frame(metrics)
#   metrics_df <- cbind(MP = rownames(metrics_df), metrics_df)
#   rownames(metrics_df) <- NULL
#   
#   # ---- Draw table to fill page ----
#   # ---- Draw table to fill page ----
#   if (requireNamespace("gridExtra", quietly = TRUE)) {
#     
#     tg <- gridExtra::tableGrob(
#       metrics_df,
#       rows = NULL,
#       theme = gridExtra::ttheme_default(
#         base_size = 12,
#         core   = list(fg_params = list(hjust = 0, x = 0.02)),
#         colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))
#       )
#     )
#     
#     # Fill page under header
#     pushViewport(viewport(x = 0.5, y = 0.45, width = 0.95, height = 0.82))
#     grid.draw(tg)
#     popViewport()
#     
#   } else {
#     
#     txt <- paste(capture.output(print(metrics_df)), collapse = "\n")
#     grid.text(txt, x = 0.05, y = 0.88, just = c("left", "top"),
#               gp = gpar(fontsize = 10, fontfamily = "mono"))
#   }
# }
# 
# dev.off()
# cat("Saved PDF:", pdf_out, "\n")
