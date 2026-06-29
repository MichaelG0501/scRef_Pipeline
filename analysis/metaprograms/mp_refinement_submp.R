####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/mp_refinement_submp.R
#   Methodology: analysis/methodology/metaprograms/mp_refinement_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Three-tier MP refinement: keep (sil >= 0.2), remove (sil < 0),
#   split (0 < sil < 0.2).  Sub-splits intermediate MPs via hierarchical
#   clustering on the NMF-program cosine-similarity matrix.  Derives sub-MP
#   gene lists using the exact GeneNMF get_metaprogram_consensus() approach.
#   Produces diagnostic plots, Fisher-Z correlation heatmap, and Jaccard
#   heatmap following analysis/metaprograms/final_mp_correlation.R style.
#
# Inputs:
#   - ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   - ref_outs/geneNMF_outs.rds  (raw multiNMF output)
#   - ref_outs/EAC_Ref_epi.rds   (epithelial Seurat object)
#
# Outputs (ref_outs/Metaprogrammes_Results/mp_refinement/):
#   intermediate/
#     refined_mp_genes.rds          – named list of gene vectors
#     refined_mp_gene_weights.rds   – named list of GeneNMF-style weights
#     refined_mp_assignments.rds    – program -> sub-MP mapping
#     refined_ucell_scores.rds      – UCell scores for refined MPs
#     split_results.rds             – full splitting diagnostic data
#     refined_mp_correlation_matrices.rds
#     refined_mp_jaccard_matrices.rds
#   tables/
#     split_selection_summary.csv
#     refined_mp_gene_sizes.csv
#     refined_mp_correlation_mean_rho.csv
#     refined_mp_jaccard_index.csv
#   figures/
#     mp_splitting_diagnostics.pdf  – silhouette + NMF heatmap per split MP
#     refined_mp_correlation_heatmap.pdf
#     refined_mp_jaccard_heatmap.pdf
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/metaprograms/mp_refinement_submp.R
#
# Conda env: dmtcp (has UCell, ComplexHeatmap, pheatmap, cluster, circlize)
# Resources observed on 7 Jun 2026: ~6 cores, ~30 GB RAM (interactive OK)
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(cluster)
library(RColorBrewer)
library(viridis)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# ============================================================================
# 1. Output directory setup
# ============================================================================

outdir <- "Metaprogrammes_Results/mp_refinement"
for (sub in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}

force_rebuild <- Sys.getenv("SCREF_FORCE_REBUILD", "FALSE") == "TRUE"
replot_only  <- Sys.getenv("SCREF_REPLOT_ONLY",  "FALSE") == "TRUE"
algorithm_version <- "mp_refinement_submp_v4_mean_sil_0.2"

# ============================================================================
# 2. Load data
# ============================================================================

cat("Loading nMP19 metaprograms object...\n")
geneNMF.metaprograms <- readRDS(
  "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)

if (!replot_only) {
  cat("Loading raw NMF programs (geneNMF_outs.rds)...\n")
  geneNMF.programs <- readRDS("geneNMF_outs.rds")
}

# ============================================================================
# 3. GeneNMF helper functions (local re-implementations)
# ============================================================================

normVector <- function(x) {
  s <- sum(x)
  if (s > 0) {
    x <- x / s
  }
  x
}

weightCumul <- function(vector, weight.explained = 0.5) {
  x.sorted <- sort(vector, decreasing = TRUE)
  norm.x   <- normVector(x.sorted)
  cs       <- cumsum(norm.x)
  norm.x[cs < weight.explained]
}

wgtLoad <- function(mat, w) {
  rownorm <- apply(mat, 1, normVector)           # K x G (each column = normalised row)
  spec    <- apply(rownorm, 2, max)               # length G: specificity
  spec.w  <- spec^w
  mat     <- mat * spec.w                         # scale rows by specificity^w
  apply(mat, 2, normVector)                       # normalise each column (pattern)
}

# ============================================================================
# 4. MP triage
# ============================================================================

metrics      <- geneNMF.metaprograms$metaprograms.metrics
sil_scores   <- metrics$silhouette
mp_names_all <- rownames(metrics)

sil_threshold <- 0.2

keep_mps   <- mp_names_all[sil_scores >= sil_threshold]
remove_mps <- mp_names_all[sil_scores < 0]
split_mps  <- mp_names_all[sil_scores > 0 & sil_scores < sil_threshold]

cat("\n=== MP Triage (silhouette threshold:", sil_threshold, ") ===\n")
for (tier in list(list("Keep", keep_mps), list("Remove", remove_mps), list("Split", split_mps))) {
  nms <- tier[[2]]
  cat(sprintf("  %-8s: %s\n", tier[[1]],
              paste(sprintf("%s(%.3f)", nms, sil_scores[match(nms, mp_names_all)]), collapse = ", ")))
}

# ============================================================================
# 5. State-group & colour definitions (from scRef_config.R)
# ============================================================================

state_groups <- list(
  "Classic Proliferative"          = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"                = c("MP13", "MP12"),
  "SMG-like Metaplasia"            = c("MP18", "MP16"),
  "Immune Infiltrating"            = c("MP15")
)

cc_mps <- c("MP1", "MP7", "MP9")

mp_desc_map <- c(
  "MP1"  = "G2M Cell Cycle",        "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",   "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",       "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",      "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)", "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi."
)

group_cols <- c(
  "Classic Proliferative"          = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"                = "#984EA3",
  "SMG-like Metaplasia"            = "#FF7F00",
  "Immune Infiltrating"            = "#377EB8"
)

state_associated_mps   <- unlist(state_groups, use.names = FALSE)
split_mps_for_final    <- intersect(split_mps, state_associated_mps)
keep_mps_for_final     <- intersect(keep_mps, state_associated_mps)

cat("\nState-associated keep :", paste(keep_mps_for_final, collapse = ", "), "\n")
cat("State-associated split:", paste(split_mps_for_final, collapse = ", "), "\n")

# ============================================================================
# 6. Sub-splitting
# ============================================================================

cached_split <- file.path(outdir, "intermediate", "split_results.rds")
split_results <- NULL
split_cache_valid <- FALSE
if (file.exists(cached_split)) {
  split_results_cached <- readRDS(cached_split)
  split_cache_valid <- identical(attr(split_results_cached, "algorithm_version"), algorithm_version)
  if (split_cache_valid) {
    split_results <- split_results_cached
  } else {
    cat("Cached split results use an older algorithm; rebuilding.\n")
  }
  rm(split_results_cached)
}

if (!replot_only && (force_rebuild || !split_cache_valid)) {

  J          <- geneNMF.metaprograms$programs.similarity
  cl_members <- geneNMF.metaprograms$programs.clusters

  split_results <- list()

  for (mp in split_mps) {
    mp_num          <- as.integer(gsub("MP", "", mp))
    member_programs <- names(cl_members)[which(cl_members == mp_num)]
    member_programs <- intersect(member_programs, rownames(J))
    n_progs         <- length(member_programs)

    cat("\n--- Splitting", mp, "(", n_progs, "programs ) ---\n")

    if (n_progs < 4) {
      cat("  Too few programs to split.  Keeping as single MP.\n")
      split_results[[mp]] <- list(skipped = TRUE, reason = "too_few_programs",
                                   member_programs = member_programs)
      next
    }

    J_sub     <- J[member_programs, member_programs]
    Jdist_sub <- as.dist(1 - J_sub)
    tree_sub  <- hclust(Jdist_sub, method = "ward.D2")

    max_k     <- n_progs - 1
    sil_by_k  <- list()
    selected_k <- NA

    for (k in 2:max_k) {
      cl      <- cutree(tree_sub, k = k)
      sil_obj <- silhouette(cl, Jdist_sub)
      current_sil <- mean(sil_obj[, "sil_width"])
      sil_by_k[[as.character(k)]] <- list(
        k        = k,
        mean_sil = current_sil,
        cl       = cl,
        sil_obj  = sil_obj,
        cl_sizes = as.integer(table(cl))
      )
      if (current_sil >= 0.2) {
        selected_k <- k
        break
      }
    }

    mean_sils     <- sapply(sil_by_k, function(x) x$mean_sil)
    
    if (is.na(selected_k)) {
      raw_best_idx  <- which.max(mean_sils)
      selected_k    <- sil_by_k[[raw_best_idx]]$k
    } else {
      raw_best_idx  <- as.character(selected_k)
    }
    
    selected_key  <- as.character(selected_k)
    best_k        <- sil_by_k[[selected_key]]$k
    best_sil      <- sil_by_k[[selected_key]]$mean_sil
    best_cl       <- sil_by_k[[selected_key]]$cl
    raw_best_k    <- best_k

    # Label sub-MPs in dendrogram order (left -> right = a, b, c, ...)
    ordered_progs     <- member_programs[tree_sub$order]
    first_appearance  <- unique(best_cl[ordered_progs])
    letter_map        <- setNames(letters[seq_along(first_appearance)], first_appearance)

    sub_labels        <- paste0(mp, letter_map[as.character(best_cl)])
    names(sub_labels) <- names(best_cl)

    cat("  Raw optimal k =", raw_best_k, "; selected k =", best_k,
        "(mean sil =", round(best_sil, 3), ")\n")
    for (sl in sort(unique(sub_labels))) {
      cat("    ", sl, ":", sum(sub_labels == sl), "programs\n")
    }

    split_results[[mp]] <- list(
      skipped     = FALSE,
      best_k      = best_k,
      best_sil    = best_sil,
      raw_best_k  = raw_best_k,
      raw_best_sil = sil_by_k[[raw_best_idx]]$mean_sil,
      cluster_assignments = best_cl,
      sub_labels  = sub_labels,
      member_programs = member_programs,
      tree        = tree_sub,
      J_sub       = J_sub,
      Jdist_sub   = Jdist_sub,
      sil_by_k    = sil_by_k,
      mean_sils   = mean_sils
    )
  }

  attr(split_results, "algorithm_version") <- algorithm_version
  saveRDS(split_results, cached_split)
  cat("\nSaved split results:", cached_split, "\n")

} else if (split_cache_valid) {
  cat("Loading cached split results...\n")
} else {
  stop("Replot-only mode requested but split cache is missing or outdated. Run without SCREF_REPLOT_ONLY first.")
}

split_selection_summary <- do.call(rbind, lapply(names(split_results), function(mp) {
  res <- split_results[[mp]]
  if (isTRUE(res$skipped)) {
    data.frame(
      original_mp = mp,
      parent_silhouette = unname(metrics[mp, "silhouette"]),
      n_programs = length(res$member_programs),
      raw_optimal_k = NA_integer_,
      selected_k = NA_integer_,
      selected_mean_silhouette = NA_real_,
      reason = res$reason,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      original_mp = mp,
      parent_silhouette = unname(metrics[mp, "silhouette"]),
      n_programs = length(res$member_programs),
      raw_optimal_k = res$raw_best_k,
      selected_k = res$best_k,
      selected_mean_silhouette = res$best_sil,
      reason = "split",
      stringsAsFactors = FALSE
    )
  }
}))
write.csv(split_selection_summary,
          file.path(outdir, "tables", "split_selection_summary.csv"),
          row.names = FALSE)
cat("Saved: tables/split_selection_summary.csv\n")

# ============================================================================
# 7. Sub-MP gene list derivation  (replicates GeneNMF consensus approach)
# ============================================================================

cached_genes <- file.path(outdir, "intermediate", "refined_mp_genes.rds")
cached_gene_weights <- file.path(outdir, "intermediate", "refined_mp_gene_weights.rds")
refined_mp_genes <- NULL
gene_cache_valid <- FALSE
if (file.exists(cached_genes)) {
  refined_mp_genes_cached <- readRDS(cached_genes)
  gene_cache_valid <- identical(attr(refined_mp_genes_cached, "algorithm_version"), algorithm_version)
  if (gene_cache_valid) {
    refined_mp_genes <- refined_mp_genes_cached
  } else {
    cat("Cached refined gene lists use an older algorithm; rebuilding.\n")
  }
  rm(refined_mp_genes_cached)
}

if (!replot_only && (force_rebuild || !gene_cache_valid)) {

  # --- 7a. Weighted loadings --------------------------------------------------
  cat("\nComputing weighted NMF loadings (specificity.weight = 5)...\n")
  nmf.wgt <- lapply(geneNMF.programs, function(model) wgtLoad(model$w, w = 5))

  # --- 7b. Gene x programs table ----------------------------------------------
  cat("Building gene x programs loading table...\n")
  gene.table <- tryCatch({
    processed <- lapply(names(nmf.wgt), function(n) {
      g <- nmf.wgt[[n]]
      colnames(g) <- paste(n, seq_len(ncol(g)), sep = ".")
      g
    })
    Reduce(cbind, processed)
  }, error = function(e) {
    cat("  Standard cbind failed — building with gene-union fallback.\n")
    all_genes <- unique(unlist(lapply(nmf.wgt, rownames)))
    gt <- matrix(0, nrow = length(all_genes), ncol = 0,
                 dimnames = list(all_genes, NULL))
    for (n in names(nmf.wgt)) {
      g <- nmf.wgt[[n]]
      colnames(g) <- paste(n, seq_len(ncol(g)), sep = ".")
      expanded <- matrix(0, nrow = length(all_genes), ncol = ncol(g),
                         dimnames = list(all_genes, colnames(g)))
      expanded[rownames(g), ] <- g
      gt <- cbind(gt, expanded)
    }
    gt
  })
  cat("  Gene table dimensions:", nrow(gene.table), "x", ncol(gene.table), "\n")

  # --- 7c. Per-program gene lists (for confidence filtering) -------------------
  cat("Computing per-program gene lists (weight.explained = 0.9, max = 1000)...\n")
  nmf.genes.single <- {
    result <- lapply(nmf.wgt, function(model) {
      gene.pass <- apply(model, 2, function(x) {
        weightCumul(x, weight.explained = 0.9)
      })
      m <- lapply(gene.pass, function(g) head(g, min(length(g), 1000)))
      isna <- sapply(m, function(x) all(is.na(x)))
      m <- m[!isna]
      names(m) <- seq_len(length(m))
      m
    })
    unlist(result, recursive = FALSE)
  }
  cat("  Total program gene lists:", length(nmf.genes.single), "\n")

  # --- 7d. Consensus gene sets for each sub-MP --------------------------------
  cat("\nDeriving sub-MP gene lists...\n")
  weight_explained <- 0.5
  max_genes        <- 200
  min_confidence   <- 0.5

  submp_genes   <- list()
  submp_weights <- list()

  for (mp in split_mps) {
    res <- split_results[[mp]]
    if (res$skipped) {
      submp_genes[[mp]] <- geneNMF.metaprograms$metaprograms.genes[[mp]]
      cat("  ", mp, ": skipped split, kept original (", length(submp_genes[[mp]]), " genes)\n")
      next
    }

    for (sub_label in sort(unique(res$sub_labels))) {
      progs <- names(res$sub_labels)[res$sub_labels == sub_label]

      # Average loadings per gene (3-SD outlier trimming)
      sub_gt   <- gene.table[, progs, drop = FALSE]
      genes.avg <- apply(as.matrix(sub_gt), 1, function(x) {
        m <- mean(x)
        if (length(x) >= 3) {
          s <- sd(x)
          x.out <- x[x > m - 3 * s & x < m + 3 * s]
        } else {
          x.out <- x
        }
        mean(x.out)
      })
      genes.avg <- sort(genes.avg, decreasing = TRUE)

      # Top genes by cumulative weight
      genes.pass <- weightCumul(genes.avg, weight.explained = weight_explained)

      # Confidence filtering: require gene in > min_confidence fraction of programs
      this <- nmf.genes.single[progs[progs %in% names(nmf.genes.single)]]
      if (length(this) > 0) {
        genes.only <- lapply(this, names)
        genes.sum  <- sort(table(unlist(genes.only)), decreasing = TRUE)
        genes.conf <- genes.sum / length(this)
        genes.conf <- genes.conf[genes.conf > min_confidence]
        genes.pass <- genes.pass[names(genes.pass) %in% names(genes.conf)]
      }

      genes.pass <- head(genes.pass, min(length(genes.pass), max_genes))

      if (length(genes.pass) > 0 && sum(genes.pass) > 0) {
        submp_weights[[sub_label]] <- genes.pass / sum(genes.pass)
        submp_genes[[sub_label]]   <- names(genes.pass)
      } else {
        submp_weights[[sub_label]] <- numeric(0)
        submp_genes[[sub_label]]   <- character(0)
      }
      cat("  ", sub_label, ":", length(submp_genes[[sub_label]]), "genes from",
          length(progs), "programs\n")
    }
  }

  # Combine with kept MP gene lists
  kept_genes       <- geneNMF.metaprograms$metaprograms.genes[keep_mps]
  kept_weights     <- geneNMF.metaprograms$metaprograms.genes.weights[keep_mps]
  refined_mp_genes <- c(kept_genes, submp_genes)
  refined_mp_gene_weights <- c(kept_weights, submp_weights)

  cat("\n=== Refined MP Gene-List Summary ===\n")
  for (nm in names(refined_mp_genes)) {
    cat("  ", nm, ":", length(refined_mp_genes[[nm]]), "genes\n")
  }

  attr(refined_mp_genes, "algorithm_version") <- algorithm_version
  attr(refined_mp_gene_weights, "algorithm_version") <- algorithm_version
  saveRDS(refined_mp_genes, cached_genes)
  saveRDS(refined_mp_gene_weights, cached_gene_weights)
  cat("Saved:", cached_genes, "\n")
  cat("Saved:", cached_gene_weights, "\n")

  # Save program-to-sub-MP mapping
  assign_rows <- list()
  for (mp in split_mps) {
    res <- split_results[[mp]]
    if (res$skipped) {
      assign_rows[[mp]] <- data.frame(
        program = res$member_programs, original_mp = mp, refined_mp = mp,
        stringsAsFactors = FALSE)
    } else {
      assign_rows[[mp]] <- data.frame(
        program = names(res$sub_labels), original_mp = mp,
        refined_mp = as.character(res$sub_labels),
        stringsAsFactors = FALSE)
    }
  }
  assignments <- do.call(rbind, assign_rows)
  rownames(assignments) <- NULL
  saveRDS(assignments, file.path(outdir, "intermediate", "refined_mp_assignments.rds"))
  cat("Saved: intermediate/refined_mp_assignments.rds\n")

} else if (gene_cache_valid) {
  cat("Loading cached gene lists...\n")
} else {
  stop("Replot-only mode requested but refined gene cache is missing or outdated. Run without SCREF_REPLOT_ONLY first.")
}

# Free large objects no longer needed
if (exists("geneNMF.programs"))  rm(geneNMF.programs)
if (exists("gene.table"))        rm(gene.table)
if (exists("nmf.wgt"))           rm(nmf.wgt)
if (exists("nmf.genes.single"))  rm(nmf.genes.single)
invisible(gc())

# ============================================================================
# 8. Build refined MP set for final plots
#    (state-associated only, CC MPs excluded, following final_mp_correlation.R)
# ============================================================================

refined_state_groups <- list()
for (state_name in names(state_groups)) {
  refined_in_state <- character(0)
  for (mp in state_groups[[state_name]]) {
    if (mp %in% keep_mps_for_final) {
      refined_in_state <- c(refined_in_state, mp)
    } else if (mp %in% split_mps_for_final) {
      res <- split_results[[mp]]
      if (res$skipped) {
        refined_in_state <- c(refined_in_state, mp)
      } else {
        refined_in_state <- c(refined_in_state, sort(unique(res$sub_labels)))
      }
    }
    # removed MPs silently excluded
  }
  if (length(refined_in_state) > 0) {
    refined_state_groups[[state_name]] <- refined_in_state
  }
}

refined_mps_ordered <- unlist(refined_state_groups, use.names = FALSE)
cat("\n=== Refined MP Order for Final Plots ===\n")
for (sn in names(refined_state_groups)) {
  cat(" ", sn, ":", paste(refined_state_groups[[sn]], collapse = ", "), "\n")
}

# Subset gene lists to final-plot MPs, drop empty ones
refined_genes_final <- refined_mp_genes[refined_mps_ordered]
zero_gene <- names(refined_genes_final)[sapply(refined_genes_final, length) == 0]
if (length(zero_gene) > 0) {
  cat("Warning: removing MPs with 0 genes:", paste(zero_gene, collapse = ", "), "\n")
  refined_genes_final <- refined_genes_final[sapply(refined_genes_final, length) > 0]
  refined_mps_ordered <- names(refined_genes_final)
  for (sn in names(refined_state_groups))
    refined_state_groups[[sn]] <- setdiff(refined_state_groups[[sn]], zero_gene)
}

# Mapping vectors ----------------------------------------------------------
mp_to_state <- unlist(lapply(names(refined_state_groups), function(gn)
  setNames(rep(gn, length(refined_state_groups[[gn]])), refined_state_groups[[gn]])))

mp_to_parent <- setNames(sub("[a-z]$", "", refined_mps_ordered), refined_mps_ordered)

# Map original MP to its state (for parent-level lookup)
mp_to_state_original <- unlist(lapply(names(state_groups), function(gn)
  setNames(rep(gn, length(state_groups[[gn]])), state_groups[[gn]])))

# Display labels: kept MPs get "MPX Description", sub-MPs get "MPXa Description"
parents_of_ordered <- sub("[a-z]$", "", refined_mps_ordered)
mp_display_labels <- setNames(
  ifelse(parents_of_ordered %in% names(mp_desc_map),
         paste(refined_mps_ordered, mp_desc_map[parents_of_ordered]),
         refined_mps_ordered),
  refined_mps_ordered
)

refined_gene_size_table <- data.frame(
  refined_mp = names(refined_mp_genes),
  parent_mp = sub("[a-z]$", "", names(refined_mp_genes)),
  in_final_plot = names(refined_mp_genes) %in% refined_mps_ordered,
  n_genes = vapply(refined_mp_genes, length, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(refined_gene_size_table,
          file.path(outdir, "tables", "refined_mp_gene_sizes.csv"),
          row.names = FALSE)
cat("Saved: tables/refined_mp_gene_sizes.csv\n")

make_gene_signature <- function(gene_list) {
  paste(vapply(names(gene_list), function(nm) {
    genes <- sort(unique(gene_list[[nm]]))
    paste(nm, length(genes), paste(head(genes, 40), collapse = ","), sep = ":")
  }, character(1)), collapse = "|")
}

# ============================================================================
# 10. UCell scoring
# ============================================================================

cached_ucell <- file.path(outdir, "intermediate", "refined_ucell_scores.rds")
ucell_signature <- make_gene_signature(refined_mp_genes)
refined_ucell <- NULL
ucell_cache_valid <- FALSE
if (file.exists(cached_ucell)) {
  refined_ucell_cached <- readRDS(cached_ucell)
  ucell_cache_valid <- identical(attr(refined_ucell_cached, "gene_signature"), ucell_signature) &&
    all(names(refined_mp_genes) %in% colnames(refined_ucell_cached))
  if (ucell_cache_valid) {
    refined_ucell <- refined_ucell_cached
  } else {
    cat("Cached refined UCell scores are outdated; rebuilding.\n")
  }
  rm(refined_ucell_cached)
}

if (force_rebuild || !ucell_cache_valid) {
  cat("\nLoading epithelial Seurat object...\n")
  tmdata_all <- readRDS("EAC_Ref_epi.rds")

  cat("Computing UCell scores for", length(refined_mp_genes), "refined MPs...\n")
  tmdata_all <- AddModuleScore_UCell(tmdata_all, features = refined_mp_genes,
                                      ncores = 1, name = "")
  refined_ucell <- tmdata_all@meta.data[, names(refined_mp_genes), drop = FALSE]
  attr(refined_ucell, "gene_signature") <- ucell_signature
  saveRDS(refined_ucell, cached_ucell)
  cat("Saved:", cached_ucell, "\n")
} else {
  cat("Loading cached UCell scores...\n")
  refined_ucell <- readRDS(cached_ucell)
  cat("Loading epithelial Seurat object (for sample metadata)...\n")
  tmdata_all <- readRDS("EAC_Ref_epi.rds")
}

# ============================================================================
# 11. Diagnostic plots (multi-page PDF for split MPs)
# ============================================================================

cat("\nGenerating diagnostic plots...\n")

compute_fisher_mean_cor <- function(score_mat, sample_vec, min_cells = 10) {
  score_mat <- as.matrix(score_mat)
  score_mat <- scale(score_mat)
  feature_names <- colnames(score_mat)
  samples <- unique(sample_vec)
  cor_array <- array(NA_real_,
                     dim = c(length(feature_names), length(feature_names), length(samples)),
                     dimnames = list(feature_names, feature_names, samples))

  for (samp in samples) {
    idx <- which(sample_vec == samp)
    if (length(idx) < min_cells) next
    cor_array[, , samp] <- cor(score_mat[idx, , drop = FALSE], method = "spearman")
  }

  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(NA_real_, length(feature_names), length(feature_names),
                     dimnames = list(feature_names, feature_names))
  for (i in seq_along(feature_names)) {
    for (j in seq_along(feature_names)) {
      if (i == j) {
        mean_rho[i, j] <- 1
      } else {
        zs <- z_array[i, j, ]
        zs <- zs[is.finite(zs)]
        mean_rho[i, j] <- if (length(zs) >= 3) tanh(mean(zs)) else NA_real_
      }
    }
  }
  mean_rho
}

original_ucell_path <- "Metaprogrammes_Results/UCell_nMP19_filtered.rds"
original_ucell <- NULL
original_mps <- setdiff(names(geneNMF.metaprograms$metaprograms.genes), remove_mps)
if (file.exists(original_ucell_path)) {
  original_ucell <- readRDS(original_ucell_path)
}
if (is.null(original_ucell) || !all(original_mps %in% colnames(original_ucell))) {
  cat("Computing original MP UCell scores for diagnostic panels...\n")
  original_genes <- geneNMF.metaprograms$metaprograms.genes[original_mps]
  tmdata_all <- AddModuleScore_UCell(tmdata_all, features = original_genes, ncores = 1, name = "")
  original_ucell <- tmdata_all@meta.data[, names(original_genes), drop = FALSE]
}

diagnostic_groups <- c(state_groups, list("Cell Cycle" = cc_mps))
diagnostic_group_cols <- c(group_cols, "Cell Cycle" = "#666666")

pdf(file.path(outdir, "figures", "mp_splitting_diagnostics.pdf"),
    width = 21, height = 8, useDingbats = FALSE)

for (mp in split_mps) {
  res <- split_results[[mp]]
  mp_title <- ifelse(mp %in% names(mp_desc_map), paste(mp, mp_desc_map[mp]), mp)
  if (res$skipped) {
    plot.new()
    text(0.5, 0.5, paste0(mp_title, "\nSkipped: ", res$reason), cex = 2)
    next
  }

  sub_labels_u <- sort(unique(as.character(res$sub_labels)))
  sub_cols <- setNames(colorRampPalette(c("#2c7fb8", "#f03b20"))(length(sub_labels_u)),
                       sub_labels_u)

  par(mfrow = c(1, 3), mar = c(5, 4.5, 4, 2))
  k_vals <- as.integer(names(res$mean_sils))
  sil_vals <- as.numeric(res$mean_sils)
  plot(k_vals, sil_vals,
       type = "b", pch = 19, lwd = 2, col = "#2b6cb0",
       xaxt = "n",
       xlab = "Number of sub-clusters (k)",
       ylab = "Mean silhouette width",
       main = paste0(mp_title, "\nSplit selection"))
  axis(1, at = k_vals)
  abline(h = 0, lty = 2, col = "grey55")
  abline(h = metrics[mp, "silhouette"], lty = 3, col = "#636363", lwd = 1.5)
  abline(v = res$raw_best_k, lty = 2, col = "#e41a1c", lwd = 1.5)
  abline(v = res$best_k, lty = 1, col = "#238b45", lwd = 1.5)
  legend("topleft",
         legend = c(paste0("Raw optimum k=", res$raw_best_k),
                    paste0("Selected k=", res$best_k),
                    paste0("Parent sil=", round(metrics[mp, "silhouette"], 3))),
         col = c("#e41a1c", "#238b45", "#636363"),
         lty = c(2, 1, 3), bty = "n", cex = 0.9)

  sub_sizes <- table(factor(res$sub_labels, levels = sub_labels_u))
  barplot(sub_sizes,
          col = sub_cols[names(sub_sizes)], border = "white",
          main = paste0("Sub-MP program counts\n", sum(sub_sizes), " NMF programs"),
          ylab = "Number of NMF programs", xlab = "Sub-MP", las = 2)

  gene_sizes <- sapply(sub_labels_u, function(lbl) length(refined_mp_genes[[lbl]]))
  barplot(gene_sizes,
          col = sub_cols[names(gene_sizes)], border = "white",
          main = "Sub-MP gene-list sizes",
          ylab = "Number of genes", xlab = "Sub-MP", las = 2)
  par(mfrow = c(1, 1))

  J_sub <- res$J_sub
  ordered_programs <- res$tree$labels[res$tree$order]
  J_ord <- J_sub[ordered_programs, ordered_programs, drop = FALSE]
  program_sub <- factor(res$sub_labels[ordered_programs], levels = sub_labels_u)
  ht_nmf <- Heatmap(J_ord,
    name = "Cosine",
    col = colorRamp2(c(0, 0.5, 1), viridis(3, option = "A", direction = -1)),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = program_sub, column_split = program_sub,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    row_title = NULL, column_title = paste0(mp_title, "\nNMF program similarity"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                                labels_gp = gpar(fontsize = 9)))

  group_name <- names(Filter(function(x) mp %in% x, diagnostic_groups))[1]
  if (length(group_name) == 0 || is.na(group_name)) group_name <- "Unassigned"
  group_members <- if (group_name %in% names(diagnostic_groups)) diagnostic_groups[[group_name]] else mp
  panel_features <- unlist(lapply(group_members, function(parent_mp) {
    if (parent_mp == mp) sub_labels_u else parent_mp
  }), use.names = FALSE)
  panel_features <- panel_features[panel_features %in% c(colnames(refined_ucell), colnames(original_ucell))]

  common_cells <- Reduce(intersect, list(rownames(refined_ucell), rownames(original_ucell), Cells(tmdata_all)))
  panel_scores <- matrix(NA_real_, nrow = length(common_cells), ncol = length(panel_features),
                         dimnames = list(common_cells, panel_features))
  for (feat in panel_features) {
    if (feat %in% colnames(refined_ucell)) {
      panel_scores[, feat] <- refined_ucell[common_cells, feat]
    } else {
      panel_scores[, feat] <- original_ucell[common_cells, feat]
    }
  }
  panel_scores <- panel_scores[, colSums(is.na(panel_scores)) == 0, drop = FALSE]
  panel_samples <- tmdata_all$orig.ident[match(rownames(panel_scores), Cells(tmdata_all))]
  panel_cor <- compute_fisher_mean_cor(panel_scores, panel_samples, min_cells = 10)

  panel_ids <- colnames(panel_cor)
  panel_parents <- sub("[a-z]$", "", panel_ids)
  panel_display <- setNames(
    ifelse(panel_parents %in% names(mp_desc_map),
           paste(panel_ids, mp_desc_map[panel_parents]),
           panel_ids),
    panel_ids
  )
  rownames(panel_cor) <- panel_display[rownames(panel_cor)]
  colnames(panel_cor) <- panel_display[colnames(panel_cor)]
  panel_state <- factor(rep(group_name, length(panel_ids)), levels = names(diagnostic_group_cols))
  panel_parent <- factor(sub("[a-z]$", "", panel_ids), levels = unique(sub("[a-z]$", "", panel_ids)))
  panel_split <- factor(paste(group_name, panel_parent, sep = "::"),
                        levels = unique(paste(group_name, panel_parent, sep = "::")))

  ha_panel_left <- rowAnnotation(
    State = panel_state,
    col = list(State = diagnostic_group_cols),
    show_annotation_name = FALSE, show_legend = FALSE)
  ha_panel_top <- HeatmapAnnotation(
    State = panel_state,
    col = list(State = diagnostic_group_cols),
    show_annotation_name = FALSE, show_legend = FALSE)

  ht_panel <- Heatmap(panel_cor,
    name = "Mean Rho",
    col = colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red")),
    rect_gp = gpar(col = "white", lwd = 0.8),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = panel_split, column_split = panel_split,
    left_annotation = ha_panel_left, top_annotation = ha_panel_top,
    row_title = NULL, column_title = paste0(group_name, "\nSub-MP vs intact group MPs"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_side = "right", column_names_side = "bottom",
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, width, height, fill) {
      rho <- panel_cor[i, j]
      grid.text(ifelse(is.na(rho), "NA", round(rho, 2)), x, y,
                gp = gpar(fontsize = 9))
    },
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                                labels_gp = gpar(fontsize = 9)))

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.52, 0.48), "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_nmf, newpage = FALSE,
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(ht_panel, newpage = FALSE,
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  popViewport()
  popViewport()
}

dev.off()
cat("Saved:", file.path(outdir, "figures", "mp_splitting_diagnostics.pdf"), "\n")

# ============================================================================
# 12. Correlation heatmap  (Fisher-Z per-sample meta-analysis)
# ============================================================================

cat("\nPreparing correlation heatmap...\n")

module_scores <- scale(as.matrix(refined_ucell[, refined_mps_ordered, drop = FALSE]))
mod_mat       <- t(module_scores)           # MPs x Cells

samples_vec <- tmdata_all$orig.ident[match(colnames(mod_mat), Cells(tmdata_all))]
samples_vec[is.na(samples_vec)] <- colnames(mod_mat)[is.na(samples_vec)]
samples <- unique(samples_vec)
mps     <- rownames(mod_mat)
n_mps   <- length(mps)

cat("Computing per-sample Spearman correlations across", length(samples), "samples...\n")

cor_array <- array(NA, dim = c(n_mps, n_mps, length(samples)),
                   dimnames = list(mps, mps, samples))

for (samp in samples) {
  cells_in <- colnames(mod_mat)[samples_vec == samp]
  if (length(cells_in) < 10) next
  cor_array[, , samp] <- cor(t(mod_mat[, cells_in, drop = FALSE]), method = "spearman")
}

z_array  <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
mean_rho <- matrix(0, n_mps, n_mps, dimnames = list(mps, mps))
p_vals   <- matrix(1, n_mps, n_mps, dimnames = list(mps, mps))

for (i in 1:n_mps) for (j in 1:n_mps) {
  if (i == j) { mean_rho[i, j] <- 1; p_vals[i, j] <- 0; next }
  zs <- z_array[i, j, ]; zs <- zs[is.finite(zs)]
  if (length(zs) < 3) { mean_rho[i, j] <- NA; p_vals[i, j] <- NA; next }
  mean_rho[i, j] <- tanh(mean(zs))
  tt <- tryCatch(t.test(zs), error = function(e) NULL)
  p_vals[i, j]   <- if (!is.null(tt)) tt$p.value else NA
}

saveRDS(list(mean_rho = mean_rho, p_values = p_vals, samples = samples),
        file.path(outdir, "intermediate", "refined_mp_correlation_matrices.rds"))
write.csv(mean_rho,
          file.path(outdir, "tables", "refined_mp_correlation_mean_rho.csv"))
cat("Saved: intermediate/refined_mp_correlation_matrices.rds\n")
cat("Saved: tables/refined_mp_correlation_mean_rho.csv\n")

# Rename for display
rownames(mean_rho) <- mp_display_labels[rownames(mean_rho)]
colnames(mean_rho) <- mp_display_labels[colnames(mean_rho)]

# --- Annotation vectors (match final_mp_correlation.R style) ----------------
state_vec   <- factor(mp_to_state[refined_mps_ordered], levels = names(state_groups))
parent_vec  <- factor(mp_to_parent, levels = unique(mp_to_parent))

# Sub-MP colour palette: each parent MP gets a unique hue;
#   sub-MPs within a parent get different luminances of that hue
unique_parents  <- levels(parent_vec)
n_parents       <- length(unique_parents)
parent_hues     <- seq(15, 345, length.out = n_parents + 1)[1:n_parents]

parent_mp_cols <- character(0)
for (i in seq_along(unique_parents)) {
  p       <- unique_parents[i]
  members <- refined_mps_ordered[mp_to_parent == p]
  n_mem   <- length(members)
  if (n_mem == 1) {
    cols <- hcl(parent_hues[i], 70, 60)
  } else {
    luminances <- seq(40, 82, length.out = n_mem)
    cols       <- hcl(parent_hues[i], 70, luminances)
  }
  parent_mp_cols <- c(parent_mp_cols, setNames(cols, members))
}
# Re-key to display labels
names(parent_mp_cols) <- mp_display_labels[names(parent_mp_cols)]

submp_factor <- factor(mp_display_labels[refined_mps_ordered],
                       levels = mp_display_labels[refined_mps_ordered])

# State annotation (colour bar, no legend — legend comes from top)
ha_left <- rowAnnotation(
  State = state_vec,
  col = list(State = group_cols),
  show_annotation_name = FALSE, show_legend = FALSE)

ha_top <- HeatmapAnnotation(
  State    = state_vec,
  `Sub-MP` = submp_factor,
  col = list(State = group_cols, `Sub-MP` = parent_mp_cols),
  show_annotation_name = TRUE, annotation_name_side = "left",
  show_legend = c(TRUE, FALSE))

# Split factor: one level per parent MP, ordered within states
split_factor <- factor(
  paste(mp_to_state[refined_mps_ordered], mp_to_parent, sep = "::"),
  levels = unique(paste(mp_to_state[refined_mps_ordered], mp_to_parent, sep = "::"))
)

# Column titles: state name for first parent per state, empty otherwise
split_to_state <- sub("::.*", "", levels(split_factor))
col_titles <- ifelse(!duplicated(split_to_state), split_to_state, "")

hm_size <- unit(13, "inch")
col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))
cell_fs <- max(6, min(10, round(220 / n_mps)))   # adaptive font size

pdf(file.path(outdir, "figures", "refined_mp_correlation_heatmap.pdf"),
    width = 20, height = 20, useDingbats = FALSE)

ht_cor <- Heatmap(mean_rho,
  name = paste0("Mean Rho\n(", length(samples), " Samples)"),
  col  = col_cor,
  rect_gp = gpar(col = "white", lwd = 0.8),
  cluster_rows = FALSE, cluster_columns = FALSE,
  left_annotation = ha_left, top_annotation = ha_top,
  row_split = split_factor, column_split = split_factor,

  column_title     = col_titles,
  column_title_rot = 20, column_title_side = "top",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  row_title = NULL,

  row_names_side    = "right", column_names_side = "bottom",
  row_names_gp      = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot  = 45,

  width = hm_size, height = hm_size,

  cell_fun = function(j, i, x, y, width, height, fill) {
    p   <- p_vals[i, j]
    rho <- mean_rho[i, j]
    if (is.na(p) || is.na(rho)) {
      grid.text("NA", x, y, gp = gpar(fontsize = cell_fs, col = "grey50"))
    } else {
      stars <- if (p < 0.001) "\n***" else if (p < 0.01) "\n**" else if (p < 0.05) "\n*" else ""
      grid.text(paste0(round(rho, 2), stars), x, y, gp = gpar(fontsize = cell_fs))
    }
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12)))

draw(ht_cor, padding = unit(c(25, 25, 25, 25), "mm"))
dev.off()
cat("Saved:", file.path(outdir, "figures", "refined_mp_correlation_heatmap.pdf"), "\n")

# --- Unsupervised correlation heatmap ---------------------------------------
pdf(file.path(outdir, "figures", "refined_mp_correlation_heatmap_unsupervised.pdf"),
    width = 20, height = 20, useDingbats = FALSE)

ht_cor_unsup <- Heatmap(mean_rho,
  name = paste0("Mean Rho\n(", length(samples), " Samples)"),
  col  = col_cor,
  rect_gp = gpar(col = "white", lwd = 0.8),
  cluster_rows = TRUE, cluster_columns = TRUE,
  left_annotation = ha_left, top_annotation = ha_top,
  
  row_title = NULL, column_title = "Unsupervised clustering",
  column_title_gp  = gpar(fontsize = 16, fontface = "bold"),
  
  row_names_side    = "right", column_names_side = "bottom",
  row_names_gp      = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot  = 45,
  
  width = hm_size, height = hm_size,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    p   <- p_vals[i, j]
    rho <- mean_rho[i, j]
    if (is.na(p) || is.na(rho)) {
      grid.text("NA", x, y, gp = gpar(fontsize = cell_fs, col = "grey50"))
    } else {
      stars <- if (p < 0.001) "\n***" else if (p < 0.01) "\n**" else if (p < 0.05) "\n*" else ""
      grid.text(paste0(round(rho, 2), stars), x, y, gp = gpar(fontsize = cell_fs))
    }
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12)))

draw(ht_cor_unsup, padding = unit(c(25, 25, 25, 25), "mm"))
dev.off()
cat("Saved:", file.path(outdir, "figures", "refined_mp_correlation_heatmap_unsupervised.pdf"), "\n")

# Free Seurat object (no longer needed)
rm(tmdata_all); invisible(gc())

# ============================================================================
# 12. Jaccard heatmap
# ============================================================================

cat("\nComputing Jaccard similarities...\n")

mp_list  <- lapply(refined_genes_final, unique)
mp_jn    <- names(mp_list)
universe <- unique(unlist(mp_list))

jaccard_mat   <- matrix(NA_real_, length(mp_list), length(mp_list),
                        dimnames = list(mp_jn, mp_jn))
overlap_n_mat <- jaccard_mat
pval_mat      <- jaccard_mat

for (i in seq_along(mp_list)) {
  A <- mp_list[[i]]
  for (j in seq_along(mp_list)) {
    B     <- mp_list[[j]]
    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j]   <- if (uni == 0) NA_real_ else inter / uni
    a  <- inter
    b  <- length(setdiff(A, B))
    cc <- length(setdiff(B, A))
    d  <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_
      else fisher.test(matrix(c(a, b, cc, d), nrow = 2),
                       alternative = "greater")$p.value
  }
}

padj_mat <- matrix(p.adjust(as.vector(pval_mat), method = "BH"),
                   nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                   dimnames = dimnames(pval_mat))

stars_mat <- matrix("", nrow(padj_mat), ncol(padj_mat), dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05]  <- "*"
stars_mat[padj_mat < 0.01]  <- "**"
stars_mat[padj_mat < 0.001] <- "***"

display_mat <- matrix(paste0(overlap_n_mat, "\n", stars_mat),
                      nrow = nrow(overlap_n_mat), ncol = ncol(overlap_n_mat),
                      dimnames = dimnames(overlap_n_mat))

saveRDS(list(jaccard = jaccard_mat,
             overlap_n = overlap_n_mat,
             p_values = pval_mat,
             padj = padj_mat),
        file.path(outdir, "intermediate", "refined_mp_jaccard_matrices.rds"))
write.csv(jaccard_mat,
          file.path(outdir, "tables", "refined_mp_jaccard_index.csv"))
cat("Saved: intermediate/refined_mp_jaccard_matrices.rds\n")
cat("Saved: tables/refined_mp_jaccard_index.csv\n")

rownames(jaccard_mat) <- mp_display_labels[rownames(jaccard_mat)]
colnames(jaccard_mat) <- mp_display_labels[colnames(jaccard_mat)]

cell_fs_j <- max(5, min(9, round(180 / n_mps)))

pdf(file.path(outdir, "figures", "refined_mp_jaccard_heatmap.pdf"),
    width = 20, height = 20, useDingbats = FALSE)

ht_jac <- Heatmap(jaccard_mat,
  name = "Jaccard Index",
  col  = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100),
  rect_gp = gpar(col = "grey85", lwd = 0.8),
  cluster_rows = FALSE, cluster_columns = FALSE,
  left_annotation = ha_left, top_annotation = ha_top,
  row_split = split_factor, column_split = split_factor,

  column_title     = col_titles,
  column_title_rot = 20, column_title_side = "top",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  row_title = NULL,

  row_names_side    = "right", column_names_side = "bottom",
  row_names_gp      = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot  = 45,

  width = hm_size, height = hm_size,

  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = cell_fs_j))
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12)))

draw(ht_jac, padding = unit(c(25, 25, 25, 25), "mm"))
dev.off()
cat("Saved:", file.path(outdir, "figures", "refined_mp_jaccard_heatmap.pdf"), "\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== MP Refinement Complete ===\n")
cat("Output dir:", file.path(getwd(), outdir), "\n")
cat("Total refined MPs:", length(refined_mps_ordered), "\n")
cat("  Kept as-is :", sum(refined_mps_ordered %in% keep_mps), "\n")
cat("  From splits:", sum(grepl("[a-z]$", refined_mps_ordered)),
    "sub-MPs from", length(split_mps_for_final), "parents\n")
cat("\nGene list sizes:\n")
for (nm in refined_mps_ordered) {
  cat("  ", sprintf("%-8s", nm), ":", length(refined_genes_final[[nm]]), "genes\n")
}
cat("\nDone.\n")
