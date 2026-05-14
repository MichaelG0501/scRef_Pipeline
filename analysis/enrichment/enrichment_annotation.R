####################
# Analysis registry:
#   Status: active
#   Script: analysis/enrichment/enrichment_annotation.R
#   Methodology: analysis/methodology/enrichment/enrichment_external_reference_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Final enrichment annotation script
####################
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

# MP_list <- readRDS("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/MP_list.RDS")
# MP_list <- MP_list$Cancer
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

geneNMF.metaprograms <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp_gene_lists <- geneNMF.metaprograms$metaprograms.genes
mp_assignments <- geneNMF.metaprograms$programs.clusters

# Filter out MPs with silhouette < 0 and NAs
print(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% bad_mp_names]

# Also filter mp_assignments to only keep valid cluster IDs
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_gene_lists)))
mp_assignments <- mp_assignments[mp_assignments %in% valid_cluster_ids & !is.na(mp_assignments)]

individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)

# 2. Map file names to their loaded data
custom_refs <- lapply(custom_files, readRDS)
names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)

# 3. Updated Loop for Metaprogram Enrichment
cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {
  
  genes <- mp_gene_lists[[mp_name]]
  mp_id <- as.numeric(gsub("\\D", "", mp_name)) 
  members <- names(mp_assignments)[mp_assignments == mp_id]
  
  message(paste0("Processing MP: ", mp_name))
  
  # Standard Enrichments
  res_GO <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                     ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
  
  res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene, 
                    TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)
  
  res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene, 
                    TERM2NAME = mp_term2name, qvalueCutoff = 0.05)
  
  # 4. Custom Marker Enrichment Loop
  res_custom_list <- lapply(names(custom_refs), function(ref_name) {
    message(paste0("  -> Running custom enrichment: ", ref_name))
    
    enricher(
      gene          = genes,
      TERM2GENE     = custom_refs[[ref_name]]$TERM2GENE,
      TERM2NAME     = custom_refs[[ref_name]]$TERM2NAME,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    )
  })
  names(res_custom_list) <- names(custom_refs)
  
  base_results <- list(
    rep_prog = mp_name,
    members  = members,
    genes    = genes,
    GO       = res_GO,
    Hallmark = res_H,
    MPs_3CA  = res_M
  )
  
  return(c(base_results, res_custom_list))
})

names(cluster_enrich) <- names(mp_gene_lists)


library(dplyr)
library(tidyr)
library(pheatmap)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
# Remove NAs and keep only valid MPs
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]

enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7, 
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {
  
  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")
  
  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    
    r <- tryCatch(er@result, error = function(e) NULL)
    # Basic structure check
    if (is.null(r) || nrow(r) == 0) return(NULL)
    
    # Filter for significance 
    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    
    data_source <- if(is_custom) r else r_sig
    
    if (nrow(data_source) == 0 && !is_custom) return(NULL)
    
    term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
    data.frame(
      Program = prog, 
      Term = term, 
      padj = data_source$p.adjust, 
      Overlap = data_source$GeneRatio, 
      stringsAsFactors = FALSE
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  if (is.null(df) || nrow(df) == 0) {
    df <- data.frame(
      Program = character(), 
      Term = character(), 
      padj = numeric(), 
      Overlap = character(), 
      stringsAsFactors = FALSE
    )
  }
  
  # 3) Define terms_use (unchanged logic)
  if (is_custom) {

    if (!element %in% names(custom_refs)) {
      message("Custom reference not found for element: ", element)
      return(invisible(NULL))
    }
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    if (nrow(df) == 0) {
      message("No significant results found for: ", element)
      return(invisible(NULL))
    }
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
  }
  
  # 4) Matrix Construction
  ordered_mps <- paste0("MP", mp_tree_order)
  # Keep only MPs that exist in cluster_enrich
  ordered_mps <- ordered_mps[ordered_mps %in% names(cluster_enrich)]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  
  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = if(element %in% c("Hallmark","GO","MPs_3CA") || is_custom) tidyr::replace_na(Overlap, "") else ""
    )
  
  # Create Score Matrix
  mat <- final_df %>%
    dplyr::select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  
  # Create Display Text Matrix
  text_mat <- final_df %>%
    dplyr::select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix() # nolint
  
  # Align columns and rows
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("No matrix content for element: ", element)
    return(invisible(NULL))
  }

  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  mp_sizes <- sapply(colnames(mat), function(x) {
    length(mp_gene_lists[[x]])
  })
  col_labels <- paste0(colnames(mat), "\nn=", mp_sizes)
  
  # 6) Sorting and Gaps logic
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
  
  # 7) Render
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  pheatmap::pheatmap(mat,
                     display_numbers = text_mat,
                     number_color = "black",
                     fontsize_number = fontsize_row * 1.1,
                     labels_col = col_labels,
                     color = cols,
                     breaks = breaks,
                     cluster_rows = cluster_rows_param, 
                     cluster_cols = FALSE,
                     gaps_row = row_gaps,
                     border_color = NA,
                     show_colnames = TRUE,
                     angle_col = 0, # Keeps text horizontal for readability with \n
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     main = paste0(element, " Enrichment (-log10 padj)"))
  
  return(invisible(mat))
}

cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
pdf("New_nMP19_anno.pdf", width = 12, height = 10)
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
cat("Saved combined PDF\n")

# Standard Runs

############

# --- Standard Sets ---

png("enrich_Hallmark.png", width = 3500, height = 1750, res = 300)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_GO.png", width = 4025, height = 2000, res = 300)
enrich_heatmap(cluster_enrich, "GO", top_per_program = 6, top_n = 60, cols = cols_palette)
dev.off()

png("enrich_MPs_3CA.png", width = 3500, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "MPs_3CA", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

# --- Custom Developmental Sets ---

png("enrich_Early_Embryogenesis.png", width = 3850, height = 1500, res = 300)
enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Normal_Development_long.png", width = 5075, height = 3000, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Normal_Development_short.png", width = 5075, height = 3000, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Organogenesis_major.png", width = 3850, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Organogenesis_sub.png", width = 4375, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Adult_Epithelium.png", width = 3938, height = 1900, res = 300)
enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Barretts_Oesophagus.png", width = 3938, height = 1900, res = 300)
enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()
