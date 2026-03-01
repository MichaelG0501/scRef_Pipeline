####################
# Moved from: analysis/enrich_plot.R
# Merged with: 3 functions from analysis/temp_plot_new.R (enrich_heatmap v2, list_to_df, mk_sheet)
# Reorganized as part of analysis/ restructuring
####################
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

save_list <- readRDS("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/enrich_dev.rds")
marker_term2gene <- save_list$TERM2GENE
marker_term2name <- save_list$TERM2NAME

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

mp_gene_lists <- geneNMF.metaprograms$metaprograms.genes
mp_assignments <- geneNMF.metaprograms$programs.clusters

cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {

  genes <- mp_gene_lists[[mp_name]]
  mp_id <- as.numeric(gsub("\\D", "", mp_name)) 
  members <- names(mp_assignments)[mp_assignments == mp_id]
  
  message(paste0("Running enrichment for ", mp_name, " (", length(genes), " genes)"))
  
  res_GO <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  res_H <- enricher(
    gene          = genes,
    TERM2GENE     = hallmark_term2gene,
    TERM2NAME     = hallmark_term2name,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  res_D <- enricher(
    gene          = genes,
    TERM2GENE     = marker_term2gene,
    TERM2NAME     = marker_term2name,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  res_M <- enricher(
    gene          = genes,
    TERM2GENE     = mp_term2gene,
    TERM2NAME     = mp_term2name,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  list(
    rep_prog = mp_name,
    members  = members,
    genes    = genes,
    GO       = res_GO,
    Hallmark = res_H,
    D        = res_D,
    M        = res_M
  )
})

names(cluster_enrich) <- names(mp_gene_lists)



library(patchwork)
library(ggplot2)

# Helper function to generate a plot or a blank placeholder
create_enrich_plot <- function(enrich_obj, title, num_cat = 6) {
  # Check if object exists and has significant results
  if (!is.null(enrich_obj) && nrow(enrich_obj@result[enrich_obj@result$p.adjust < 0.05, ]) > 0) {
    barplot(enrich_obj, showCategory = num_cat) +
      ggtitle(title) +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title  = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11), 
        legend.position = "none",
        plot.title   = element_text(size = 14, face = "bold", hjust = 0.5)
      )
  } else {
    # Return an empty plot object if no results
    ggplot() + theme_void() + ggtitle(paste0(title, "\n(No Sig. Terms)"))
  }
}

aspdf <- FALSE
if (aspdf) {
  pdf("Enrichment_Analysis.pdf", width = 18, height = 10)
}

for (rp in names(cluster_enrich)[c(3,6,7)]) {
  cat("Plotting enrichment for cluster represented by:", rp, "\n")
  enr <- cluster_enrich[[rp]]
  
  p1 <- create_enrich_plot(enr$GO,       "GO Biological Process")
  p2 <- create_enrich_plot(enr$Hallmark, "Hallmark Pathways")
  p3 <- create_enrich_plot(enr$D,        "Developmental Signatures")
  p4 <- create_enrich_plot(enr$M,        "Metaprogram Reference")
  
  # Combine into a layout
  combined <- (p1 + p2 + p3 + p4) + 
    plot_layout(ncol = 2, guides = 'collect') + 
    plot_annotation(
      title = paste("Enrichment Analysis:", rp),
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  print(combined)
}

if (aspdf) {
  dev.off()
  cat("All plots saved to 'Enrichment_Analysis_All_Clusters.pdf'\n")
}


####################### Heatmaps #########################

library(dplyr)
library(tidyr)
library(pheatmap)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)

enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 5,   # NEW: take top 10 from each program
                           top_n = 50,             # cap total terms shown (after union)
                           padj_cut = 0.05,        # used only for scoring mask if you want; NOT required for selection
                           cols = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100),
                           cap = 7, fontsize_row = 7, fontsize_col = 9) {
  
  # 1) Long table: Program, Term, padj
  df <- bind_rows(lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    
    term <- if ("Description" %in% colnames(r)) r$Description else r$ID
    data.frame(Program = prog, Term = term, padj = r$p.adjust, stringsAsFactors = FALSE)
  })) %>%
    filter(!is.na(padj), is.finite(padj), padj > 0)
  
  if (nrow(df) == 0) {
    message("No results for element: ", element)
    return(invisible(NULL))
  }
  
  # 2) NEW term selection:
  #    union of top K terms per program (smallest padj), then trim to top_n by global best padj
  terms_use <- df %>%
    arrange(Program, padj) %>%
    group_by(Program) %>%
    slice_head(n = top_per_program) %>%
    ungroup() %>%
    distinct(Term) %>%
    pull(Term)
  
  # global ranking to optionally cap total displayed terms
  term_rank <- df %>%
    group_by(Term) %>%
    summarise(best_padj = min(padj), .groups = "drop") %>%
    arrange(best_padj)
  
  if (length(terms_use) > top_n) {
    terms_use <- term_rank %>%
      filter(Term %in% terms_use) %>%
      slice_head(n = top_n) %>%
      pull(Term)
  }
  
  # 3) Term x Program matrix of -log10(padj), missing = 0
  mat <- df %>%
    filter(Term %in% terms_use) %>%
    mutate(score = pmin(-log10(padj), cap)) %>%
    group_by(Term, Program) %>%
    summarise(score = max(score), .groups = "drop") %>%
    pivot_wider(names_from = Program, values_from = score, values_fill = 0) %>%
    as.data.frame()
  
  rownames(mat) <- mat$Term
  mat$Term <- NULL
  mat <- as.matrix(mat)
  
  ordered_mps <- paste0("MP", mp_tree_order)[-10]
  mat <- mat[, ordered_mps, drop = FALSE]
  
  # 1) Row Filtering (Subset to desired terms)
  desired <- if(element == "M") unique(mp_term2gene$term) else if(element == "D") unique(marker_term2gene$term) else rownames(mat)
  keep <- intersect(desired, rownames(mat))
  mat <- mat[keep, , drop = FALSE]
  
  # 2) Row Sorting Logic
  if (element == "D") {
    # Sort by Developmental Category (reversed) then significance
    lvls <- c("Adult_Stomach", "Adult_Oesophagus", "Normal_Development", "Organogenesis", "Early_Embryogenesis")
    cats <- factor(sub(".*\\.\\.", "", rownames(mat)), levels = lvls)
    mat  <- mat[order(cats, -rowSums(mat)), , drop = FALSE]
    
    # Define Gaps by Category
    groups <- sub(".*\\.\\.", "", rownames(mat))
    row_gaps <- which(groups[-length(groups)] != groups[-1])
    
  } else if (element == "M") {
    # Cluster all rows (default NMF order) - No Gaps
    row_gaps <- NULL 
    
  } else {
    # GO/Hallmark: Sort by best MP rank in tree
    best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_order <- order(match(best_mp, ordered_mps), rownames(mat))
    mat <- mat[row_order, , drop = FALSE]
    
    # Define Gaps by MP transition
    groups <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_gaps <- which(groups[-length(groups)] != groups[-1])
  }
  
  # 3) Plot
  pheatmap(mat,
           color = cols,
           cluster_rows = (element == "M"), # Only cluster rows for M
           cluster_cols = FALSE,
           gaps_row = row_gaps,
           border_color = NA,
           show_colnames = TRUE,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col,
           main = paste0(element, " enrichment (-log10 padj)"))
  
  invisible(mat)
}

# --------- run all four heatmaps ----------
cols <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
cols <- viridis::magma(100, direction = -1)

enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 10, top_n = 80, cols = cols)
enrich_heatmap(cluster_enrich, "GO",       top_per_program = 10, top_n = 80, cols = cols)
enrich_heatmap(cluster_enrich, "D",        top_per_program = 10, top_n = 80, cols = cols)
enrich_heatmap(cluster_enrich, "M",        top_per_program = 10, top_n = 80, cols = cols)


####################
# Merged from: analysis/temp_plot_new.R
# Functions: enrich_heatmap v2 (most sophisticated), list_to_df, mk_sheet
####################

# enrich_heatmap v2 — sophisticated version with custom developmental RDS support,
# GeneRatio display text, and viridis magma colors
enrich_heatmap_v2 <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7, 
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {
  
  # 1) Detect if this is one of your Custom RDS elements
  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")
  
  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    
    r <- tryCatch(er@result, error = function(e) NULL)
    # Basic structure check
    if (is.null(r) || nrow(r) == 0) return(NULL)
    
    # Filter for significance 
    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    
    # Logic: use significant for colors, but we need structure even if empty
    data_source <- if(is_custom) r else r_sig
    
    # If significant data is missing and it's not custom, skip this program
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
  
  # --- FIX: Ensure df has columns even if it has 0 rows ---
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
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    # For GO/H/M, if df is empty after initialization, we still stop
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
  
  # 4) Matrix Construction (Now safe because df has columns)
  ordered_mps <- paste0("MP", mp_tree_order)[-10]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  
  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      # padj exists now, so pmin won't error
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      # Overlap exists now, so replace_na won't error
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
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  
  # Align columns and rows
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]
  
  mp_sizes <- sapply(colnames(mat), function(x) {
    # If the column name is "MP1", we look up mp_gene_lists[["MP1"]]
    length(mp_gene_lists[[x]])
  })
  col_labels <- paste0(colnames(mat), "\nn=", mp_sizes)
  
  # 6) Sorting and Gaps logic (Keep same as before)
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

# list_to_df — convert named list of gene vectors to padded data.frame
list_to_df <- function(gene_list) {
  max_len <- max(sapply(gene_list, length))
  mat <- sapply(gene_list, function(g)
    c(g, rep(NA, max_len - length(g)))
  )
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# mk_sheet — create a formatted sheet from term2gene data for writexl
mk_sheet <- function(df, sheet_name) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)[, c("term", "gene")]
  sig <- sub("\\.\\..*$", "", df$term)
  genes <- split(df$gene, sig)
  
  # Create matrix and add the cell type names as the first column
  maxn  <- max(lengths(genes))
  mat   <- t(vapply(genes, function(g) c(g, rep(NA, maxn - length(g))), character(maxn)))
  out   <- cbind(rownames(mat), mat)
  
  # Add the Title row at the very top
  out <- rbind(c(sheet_name, rep(NA, ncol(out) - 1)), out)
  
  # Convert to DF and set column names to empty strings to hide "V1, V2..."
  out_df <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out_df) <- rep("", ncol(out_df))
  
  return(out_df)
}
