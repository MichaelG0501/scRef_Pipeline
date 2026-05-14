####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/mp_database_correlation.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/Auto_MP_correlation_v2.R
# Reorganized as part of analysis/ restructuring
####################
####################
# Auto_MP_correlation_v2.R
#
# For each database used in MP annotation (example_anno.R), this script:
#   1. Extracts ALL reference terms that appear in the enrichment results
#      (not just significant ones), matching the same term set used in the
#      overlap enrichment heatmap from example_anno.R
#   2. UCell-scores all epithelial cells using each reference term's gene list
#   3. Computes per-sample Spearman correlation between reference UCell scores
#      and MP UCell scores (Fisher Z-transform → t-test across samples)
#   4. Plots a ComplexHeatmap per database: rows = ref terms, cols = MPs,
#      cells = mean rho with significance stars
#
# Key differences from v1:
#   - GO terms display actual Description names, not GO IDs
#   - ALL terms from each database are included (not just padj < 0.05),
#     so the heatmap layout matches the overlap enrichment plot exactly
#
# Pattern follows compare_pdos_sc.R (lines 314-389).
# Run in dmtcp conda environment.
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ══════════════════════════════════════════════════════════════════════════════
# 0. Load data
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading data...\n")
cluster_enrich <- readRDS("cluster_enrich.rds")
mp_ucell       <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds")       # 75348 cells x 9 MPs
tmdata_all     <- readRDS("EAC_Ref_epi.rds")

# MP gene lists (for reference / naming)
mp_outs <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp_names <- names(mp_outs$metaprograms.genes)  # MP1..MP9

# Filter out MPs with silhouette < 0 (same logic as example_anno.R lines 30-38)
bad_mps <- which(mp_outs$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
cat(sprintf("Filtering out %d MPs with silhouette < 0: %s\n", length(bad_mp_names), paste(bad_mp_names, collapse = ", ")))

# Filter mp_names
mp_names <- mp_names[!mp_names %in% bad_mp_names]
cat(sprintf("Remaining MPs: %s\n", paste(mp_names, collapse = ", ")))

# Filter cluster_enrich (keep only MPs with valid silhouette)
cluster_enrich <- cluster_enrich[names(cluster_enrich) %in% mp_names]
cat(sprintf("Cluster enrich entries: %d\n", length(cluster_enrich)))

# Filter mp_ucell (keep only columns for valid MPs)
mp_ucell <- mp_ucell[, mp_names, drop = FALSE]

# Filter mp_outs$metaprograms.genes
mp_outs$metaprograms.genes <- mp_outs$metaprograms.genes[mp_names]

cat(sprintf("Cells in UCell: %d | Cells in Seurat: %d\n",
            nrow(mp_ucell), ncol(tmdata_all)))

# Align cells: use intersection
common_cells <- intersect(rownames(mp_ucell), colnames(tmdata_all))
cat(sprintf("Common cells: %d\n", length(common_cells)))
mp_ucell   <- mp_ucell[common_cells, ]
tmdata_all <- tmdata_all[, common_cells]

# Load custom reference TERM2NAME for custom databases (needed for term lists)
individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
custom_refs <- lapply(custom_files, readRDS)
names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)
cat(sprintf("Loaded %d custom reference databases\n", length(custom_refs)))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Define databases to process
# ══════════════════════════════════════════════════════════════════════════════
databases <- c("Hallmark", "GO", "MPs_3CA",
               "Adult_Epithelium", "Early_Embryogenesis",
               "Normal_Development_long", "Normal_Development_short",
               "Organogenesis_major", "Organogenesis_sub")

# Which databases are "custom" (developmental) vs standard
custom_dbs <- c("Adult_Epithelium", "Early_Embryogenesis",
                "Normal_Development_long", "Normal_Development_short",
                "Organogenesis_major", "Organogenesis_sub")
standard_dbs <- c("Hallmark", "GO", "MPs_3CA")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Extract ALL terms per database (matching example_anno.R term selection)
#
#    For custom databases: use ALL terms from TERM2NAME
#    For standard databases (GO, Hallmark, MPs_3CA): use the same
#      top_per_program / top_n selection as example_anno.R's enrich_heatmap
#
#    Also builds an ID → Description mapping for GO terms.
# ══════════════════════════════════════════════════════════════════════════════

# Per-database caps mirroring enrich_heatmap settings from example_anno.R
db_caps <- list(
  Hallmark = list(top_per_mp = 8, max_terms = 80),
  GO       = list(top_per_mp = 6, max_terms = 60),
  MPs_3CA  = list(top_per_mp = 8, max_terms = 80)
)

extract_all_terms <- function(cluster_enrich, db_name, is_custom,
                              custom_refs = NULL,
                              top_per_mp = NULL, max_terms = NULL) {
  # -------------------------------------------------------------------
  # Step 1: Build ID → Description mapping (needed for GO)
  # -------------------------------------------------------------------
  id_to_desc <- character(0)
  
  for (mp in names(cluster_enrich)) {
    er <- cluster_enrich[[mp]][[db_name]]
    if (is.null(er)) next
    
    res <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    
    # Build mapping: if Description column exists, map ID -> Description
    if ("Description" %in% colnames(res)) {
      new_map <- setNames(res$Description, res$ID)
      # Only add IDs we haven't seen yet
      new_ids <- setdiff(names(new_map), names(id_to_desc))
      if (length(new_ids) > 0) {
        id_to_desc <- c(id_to_desc, new_map[new_ids])
      }
    }
  }
  
  # -------------------------------------------------------------------
  # Step 2: Determine which terms to use
  # -------------------------------------------------------------------
  if (is_custom) {
    # Custom databases: use ALL terms from TERM2NAME
    terms_use <- as.character(custom_refs[[db_name]]$TERM2NAME$term)
    
  } else {
    # Standard databases: replicate example_anno.R logic
    # Collect all results (significant ones for term selection)
    df_list <- lapply(names(cluster_enrich), function(mp) {
      er <- cluster_enrich[[mp]][[db_name]]
      if (is.null(er)) return(NULL)
      
      res <- tryCatch(er@result, error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) return(NULL)
      
      # Filter significant
      r_sig <- res[which(res$p.adjust < 0.05 & res$p.adjust > 0), ]
      if (nrow(r_sig) == 0) return(NULL)
      
      # Use Description if available, otherwise ID
      term_col <- if ("Description" %in% colnames(r_sig)) r_sig$Description else r_sig$ID
      
      data.frame(
        Program = mp,
        Term = term_col,
        padj = r_sig$p.adjust,
        stringsAsFactors = FALSE
      )
    })
    
    df <- dplyr::bind_rows(df_list)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Same selection logic as example_anno.R enrich_heatmap
    terms_use <- df %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::arrange(Program, padj) %>%
      dplyr::group_by(Program) %>%
      dplyr::slice_head(n = top_per_mp) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Term) %>%
      dplyr::pull(Term)
    
    if (!is.null(max_terms) && length(terms_use) > max_terms) {
      terms_use <- df %>%
        dplyr::filter(Term %in% terms_use) %>%
        dplyr::group_by(Term) %>%
        dplyr::summarise(min_p = min(padj), .groups = "drop") %>%
        dplyr::arrange(min_p) %>%
        dplyr::slice_head(n = max_terms) %>%
        dplyr::pull(Term)
    }
  }
  
  if (length(terms_use) == 0) return(NULL)
  
  # -------------------------------------------------------------------
  # Step 3: For each term, find the gene list from @geneSets
  #         Use Description-based names for the output where applicable
  # -------------------------------------------------------------------
  # We need to map terms_use back to IDs for geneSets lookup
  # For custom dbs: terms_use ARE the IDs (term names used directly)
  # For standard dbs with Description: terms_use are Descriptions, need reverse map
  
  desc_to_id <- setNames(names(id_to_desc), id_to_desc)
  
  gene_lists <- list()
  
  for (mp in names(cluster_enrich)) {
    er <- cluster_enrich[[mp]][[db_name]]
    if (is.null(er)) next
    
    gs <- tryCatch(er@geneSets, error = function(e) NULL)
    if (is.null(gs) || length(gs) == 0) next
    
    for (term in terms_use) {
      if (term %in% names(gene_lists)) next  # already found
      
      # Try direct match first (works for custom dbs and non-GO standard)
      if (term %in% names(gs)) {
        gene_lists[[term]] <- gs[[term]]
        next
      }
      
      # For GO: term is Description, need to look up by ID
      if (term %in% names(desc_to_id)) {
        go_id <- desc_to_id[[term]]
        if (go_id %in% names(gs)) {
          gene_lists[[term]] <- gs[[go_id]]
        }
      }
    }
    
    if (length(gene_lists) == length(terms_use)) break
  }
  
  # Order gene_lists to match terms_use order
  found_terms <- intersect(terms_use, names(gene_lists))
  gene_lists <- gene_lists[found_terms]
  
  cat(sprintf("    Requested %d terms, found gene lists for %d\n",
              length(terms_use), length(gene_lists)))
  
  return(gene_lists)
}

cat("Extracting ALL terms per database (matching example_anno.R layout)...\n")
db_gene_lists <- list()

for (db in databases) {
  is_cust <- db %in% custom_dbs
  caps <- if (is_cust) list(top_per_mp = NULL, max_terms = NULL) else db_caps[[db]]
  
  gl <- extract_all_terms(
    cluster_enrich, db,
    is_custom   = is_cust,
    custom_refs = if (is_cust) custom_refs else NULL,
    top_per_mp  = caps$top_per_mp,
    max_terms   = caps$max_terms
  )
  
  if (!is.null(gl) && length(gl) > 0) {
    db_gene_lists[[db]] <- gl
    cat(sprintf("  %s: %d terms\n", db, length(gl)))
  } else {
    cat(sprintf("  %s: no terms found — skipping\n", db))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. UCell scoring of reference terms
#
#    For each database, score all cells using the reference gene lists.
#    We batch all databases together in a single AddModuleScore_UCell call
#    to avoid repeated passes over the expression matrix.
# ══════════════════════════════════════════════════════════════════════════════
cat("Preparing gene lists for UCell scoring...\n")

# Flatten all gene lists with database prefix to avoid name collisions
all_gene_lists <- list()
term_to_db <- character(0)

for (db in names(db_gene_lists)) {
  for (term in names(db_gene_lists[[db]])) {
    # Create safe name: db__term (double underscore separator)
    safe_name <- paste0(db, "__", make.names(term))
    all_gene_lists[[safe_name]] <- db_gene_lists[[db]][[term]]
    term_to_db[safe_name] <- db
  }
}

cat(sprintf("Total gene sets to score: %d\n", length(all_gene_lists)))

# UCell scoring in batches to manage memory
cat("Running UCell scoring (this may take a while)...\n")

batch_size <- 50
safe_names <- names(all_gene_lists)
n_batches  <- ceiling(length(safe_names) / batch_size)

ref_scores <- data.frame(row.names = common_cells)

for (b in seq_len(n_batches)) {
  idx_start <- (b - 1) * batch_size + 1
  idx_end   <- min(b * batch_size, length(safe_names))
  batch_names <- safe_names[idx_start:idx_end]
  batch_lists <- all_gene_lists[batch_names]

  cat(sprintf("  Batch %d/%d: scoring %d gene sets...\n", b, n_batches, length(batch_lists)))

  tmdata_all <- AddModuleScore_UCell(tmdata_all, features = batch_lists,
                                     ncores = 4, name = "")

  # Extract scores from meta.data
  score_cols <- intersect(batch_names, colnames(tmdata_all@meta.data))
  if (length(score_cols) > 0) {
    ref_scores <- cbind(ref_scores,
                        tmdata_all@meta.data[common_cells, score_cols, drop = FALSE])
    # Clean up meta.data to save memory
    tmdata_all@meta.data <- tmdata_all@meta.data[,
      !(colnames(tmdata_all@meta.data) %in% score_cols), drop = FALSE]
  }
}

cat(sprintf("UCell scoring complete: %d cells x %d reference terms\n",
            nrow(ref_scores), ncol(ref_scores)))

# Save reference UCell scores for reuse
saveRDS(ref_scores, file = "UCell_ref_terms_v2_MP19.rds")
cat("Saved UCell_ref_terms_v2.rds\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Per-sample Spearman cross-correlation
#
#    Following compare_pdos_sc.R (lines 314-389):
#    - mod_mat = MP UCell scores [cells x MPs] → transpose to [MPs x cells]
#    - save    = ref UCell scores [cells x refs] → transpose to [refs x cells]
#    - For each sample: cross-correlate MP scores with ref scores
#    - Fisher Z → t-test → mean_rho + p-values
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing per-sample cross-correlations...\n")

# Transpose to [features x cells] as in compare_pdos_sc.R
mod_mat <- t(as.matrix(mp_ucell))   # [MPs x cells]

# Sample IDs from Seurat metadata
sample_ids <- as.character(tmdata_all$orig.ident[colnames(mod_mat)])
samples    <- unique(sample_ids)
cat(sprintf("Found %d unique samples\n", length(samples)))

# Function: per-sample cross-correlation for one database
compute_cross_cor <- function(ref_score_mat, mod_mat, sample_ids, samples,
                              min_cells = 5) {
  n_ref <- nrow(ref_score_mat)
  n_mp  <- nrow(mod_mat)

  cor_array <- array(NA, dim = c(n_ref, n_mp, length(samples)),
                     dimnames = list(rownames(ref_score_mat),
                                     rownames(mod_mat),
                                     samples))

  for (smp in samples) {
    idx <- which(sample_ids == smp)
    if (length(idx) > min_cells) {
      cor_sample <- cor(
        t(ref_score_mat[, idx, drop = FALSE]),
        t(mod_mat[, idx, drop = FALSE]),
        method = "spearman"
      )
      cor_array[, , smp] <- cor_sample
    }
  }

  # Fisher Z-transform
  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

  # T-test across samples for each (ref, MP) pair
  mean_rho <- matrix(0, n_ref, n_mp,
                     dimnames = list(rownames(ref_score_mat), rownames(mod_mat)))
  p_vals   <- matrix(1, n_ref, n_mp,
                     dimnames = list(rownames(ref_score_mat), rownames(mod_mat)))

  for (i in seq_len(n_ref)) {
    for (j in seq_len(n_mp)) {
      z_scores <- z_array[i, j, ]
      z_scores <- z_scores[!is.na(z_scores)]

      if (length(z_scores) > 1) {
        test_res      <- t.test(z_scores)
        mean_rho[i, j] <- tanh(mean(z_scores))
        p_vals[i, j]   <- test_res$p.value
      }
    }
  }

  list(mean_rho = mean_rho, p_vals = p_vals)
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Run per-database and plot
# ══════════════════════════════════════════════════════════════════════════════
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

all_results <- list()

pdf("Auto_MP_correlation_heatmaps_v2_MP19.pdf", width = 14, height = 10)
for (db in names(db_gene_lists)) {
  cat(sprintf("\n--- Processing database: %s ---\n", db))

  # Get score columns for this database
  db_prefix <- paste0(db, "__")
  db_cols <- grep(paste0("^", db_prefix), colnames(ref_scores), value = TRUE)

  if (length(db_cols) == 0) {
    cat("  No scored terms — skipping\n")
    next
  }

  # Build ref score matrix [refs x cells]
  ref_mat <- t(as.matrix(ref_scores[, db_cols, drop = FALSE]))

  # Clean up row names: remove db prefix, restore original term name
  # The safe_name was db__make.names(term), so strip the prefix
  clean_names <- sub(paste0("^", gsub("([.])", "\\\\\\1", db_prefix)), "",
                     rownames(ref_mat))
  # Undo make.names: replace ".." with " | " and "." with " "
  clean_names <- gsub("\\.\\.", " | ", clean_names)
  clean_names <- gsub("\\.", " ", clean_names)
  rownames(ref_mat) <- clean_names

  cat(sprintf("  %d terms x %d cells\n", nrow(ref_mat), ncol(ref_mat)))

  # Compute cross-correlation
  result <- compute_cross_cor(ref_mat, mod_mat, sample_ids, samples)

  # Adjust p-values across all tests within this database (BH)
  padj <- matrix(p.adjust(as.vector(result$p_vals), method = "BH"),
                 nrow = nrow(result$p_vals),
                 dimnames = dimnames(result$p_vals))

  # Store results for RDS export
  all_results[[db]] <- list(mean_rho = result$mean_rho,
                            p_vals   = result$p_vals,
                            padj     = padj)

  # Dynamic font sizing based on number of rows
  n_rows <- nrow(result$mean_rho)
  row_fontsize <- max(4, min(9, 200 / n_rows))
  cell_fontsize <- max(4, min(8, 150 / n_rows))

  # Plot
  ht <- Heatmap(
    result$mean_rho,
    name = "Mean Rho",
    col  = col_fun,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows    = "ward.D2",
    clustering_method_columns = "ward.D2",
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_title    = "Reference Terms",
    column_title = paste0(db, " — Expression Correlation with MPs"),
    row_names_gp    = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = 10),

    cell_fun = function(j, i, x, y, width, height, fill) {
      p   <- padj[i, j]
      rho <- result$mean_rho[i, j]
      lvl <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      grid.text(sprintf("%.2f\n%s", rho, lvl), x, y,
                gp = gpar(fontsize = cell_fontsize,
                          fontface = "bold"))
    }
  )

  draw(ht, merge_legend = TRUE)
  cat(sprintf("  Plotted %d terms x %d MPs\n", nrow(result$mean_rho), ncol(result$mean_rho)))
}

dev.off()
cat("\nSaved Auto_MP_correlation_heatmaps_v2.pdf\n")

saveRDS(all_results, file = "Auto_MP_correlation_results_v2_MP19.rds")
cat("Saved Auto_MP_correlation_results_v2.rds\n")
cat("Done.\n")
