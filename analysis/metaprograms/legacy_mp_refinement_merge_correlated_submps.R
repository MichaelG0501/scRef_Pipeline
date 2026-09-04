####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_mp_refinement_merge_correlated_submps.R
#   Methodology: analysis/methodology/metaprograms/mp_refinement_merge_correlated_submps_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Downstream merge layer for analysis/metaprograms/mp_refinement_submp.R.
#   After sub-MP refinement, merge sub-MPs from the same parent MP when each
#   sub-MP has Spearman mean rho > 0.4 with at least 25% of the other sub-MPs
#   from that parent. Re-derive merged MP gene lists from the pooled parent NMF
#   programs using the same GeneNMF-style consensus method as mp_refinement_submp.R.
#
# Inputs:
#   - ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   - ref_outs/geneNMF_outs.rds
#   - ref_outs/EAC_Ref_epi.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/split_results.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_gene_weights.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_ucell_scores.rds
#
# Outputs:
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_gene_weights.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_correlation_matrices.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/merged_refined_mp_merge_decisions.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/merged_refined_mp_gene_sizes.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/merged_refined_mp_correlation_mean_rho.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/merged_refined_MP_genes_summary.xlsx
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_correlation_heatmap_unsupervised_merged.pdf
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_correlation_heatmap_unsupervised_merged_all.pdf
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/merged_refined_mp_enrichment_anno.pdf
#   - updates/new_updates/summaries/mp_refinement_merge_correlated_submps_summary.csv
#
# Cache/replot behavior:
#   SCREF_FORCE_REBUILD=TRUE rebuilds merged gene/UCell caches.
#   SCREF_REPLOT_ONLY=TRUE reuses valid merged caches and regenerates tables/figures.
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/metaprograms/mp_refinement_merge_correlated_submps.R
#
# Conda env: dmtcp
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(openxlsx)


project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
ref_dir <- file.path(project_dir, "ref_outs")
setwd(ref_dir)

outdir <- "Metaprogrammes_Results/mp_refinement"
for (sub in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(project_dir, "updates/new_updates/summaries"),
           recursive = TRUE, showWarnings = FALSE)

force_rebuild <- Sys.getenv("SCREF_FORCE_REBUILD", "FALSE") == "TRUE"
replot_only <- Sys.getenv("SCREF_REPLOT_ONLY", "FALSE") == "TRUE"
cor_threshold <- as.numeric(Sys.getenv("SCREF_SUBMP_MERGE_COR", "0.4"))
fraction_threshold <- as.numeric(Sys.getenv("SCREF_SUBMP_MERGE_FRACTION", "0.25"))
algorithm_version <- paste0(
  "mp_refinement_merge_correlated_submps_v1_cor_",
  cor_threshold, "_frac_", fraction_threshold
)

state_groups <- list(
  "Classic Proliferative"          = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"                = c("MP13", "MP12"),
  "SMG-like Metaplasia"            = c("MP18", "MP16"),
  "Immune Infiltrating"            = c("MP15")
)

mp_desc_map <- c(
  "MP1"  = "G2M Cell Cycle", "MP9" = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition", "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.", "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.", "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair", "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)", "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi."
)

submp_desc_map <- c(
  "MP7+"  = "Fanconi/HR repair progenitor",
  "MP7h"  = "Replication-stress dormant epithelial",
  "MP7j"  = "DNA damage response",
  "MP7r"  = "Stem-like glandular duct progenitor",
  "MP7v"  = "Mucous secretory progenitor",
  "MP10+" = "Metabolic columnar epithelium",
  "MP10e" = "Inflammatory mucous-secretory columnar epithelium",
  "MP8+"  = "Intestinal metaplasia",
  "MP8b"  = "Quiescent glandular-metabolic progenitor",
  "MP8c"  = "NF-κB inflammatory cycling glandular progenitor",
  "MP8e"  = "Cycling intestinal–columnar progenitor",
  "MP12a" = "Enteroendocrine-primed progenitor",
  "MP12b" = "Enteroendocrine differentiation",
  "MP12c" = "Cycling glandular–intestinal progenitor",
  "MP2+"  = "MYC proliferation",
  "MP2v"  = "EMT-V cycling invasive progenitor",
  "MP15a" = "T/NK-like epithelial immune mimicry",
  "MP15b" = "Type I IFN–activated EMT-primed epithelial",
  "MP15c" = "Type II IFN / NF-κB peak inflammatory epithelial",
  "MP5+"  = "Epithelial IFN Resp.",
  "MP16+" = "Secretory Diff. (Gastric)"
)

parent_mp_colors <- c(
  "MP1"  = "#E41A1C", "MP2"  = "#FF7F00", "MP5"  = "#984EA3",
  "MP7"  = "#4DAF4A", "MP8"  = "#377EB8", "MP9"  = "#F781BF",
  "MP10" = "#A65628", "MP12" = "#999999", "MP13" = "#66C2A5",
  "MP14" = "#FC8D62", "MP15" = "#8DA0CB", "MP16" = "#E78AC3",
  "MP17" = "#A6D854", "MP18" = "#FFD92F"
)

group_cols <- c(
  "Classic Proliferative"          = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"                = "#984EA3",
  "SMG-like Metaplasia"            = "#FF7F00",
  "Immune Infiltrating"            = "#377EB8"
)

cat("Loading MP refinement inputs...\n")
geneNMF.metaprograms <- readRDS(
  "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)
split_results <- readRDS(file.path(outdir, "intermediate", "split_results.rds"))
refined_mp_genes <- readRDS(file.path(outdir, "intermediate", "refined_mp_genes.rds"))
refined_mp_gene_weights <- readRDS(file.path(outdir, "intermediate", "refined_mp_gene_weights.rds"))
refined_ucell <- readRDS(file.path(outdir, "intermediate", "refined_ucell_scores.rds"))

metrics <- geneNMF.metaprograms$metaprograms.metrics
sil_scores <- metrics$silhouette
mp_names_all <- rownames(metrics)
sil_threshold <- 0.2

keep_mps <- mp_names_all[sil_scores >= sil_threshold]
remove_mps <- mp_names_all[sil_scores < 0]
split_mps <- mp_names_all[sil_scores > 0 & sil_scores < sil_threshold]
state_associated_mps <- unlist(state_groups, use.names = FALSE)

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

display_label <- function(x) {
  vapply(x, function(xi) {
    if (xi %in% names(submp_desc_map)) return(paste(xi, submp_desc_map[xi]))
    parent <- parent_id(xi)
    if (parent %in% names(mp_desc_map)) return(paste(xi, mp_desc_map[parent]))
    return(xi)
  }, character(1), USE.NAMES = FALSE)
}

make_gene_signature <- function(gene_list) {
  paste(vapply(names(gene_list), function(nm) {
    genes <- sort(unique(gene_list[[nm]]))
    paste(nm, length(genes), paste(head(genes, 40), collapse = ","), sep = ":")
  }, character(1)), collapse = "|")
}

compute_fisher_cor_p <- function(score_mat, sample_vec, min_cells = 10) {
  score_mat <- scale(as.matrix(score_mat))
  feature_names <- colnames(score_mat)
  samples <- unique(sample_vec)
  n_features <- length(feature_names)

  cor_array <- array(NA_real_,
                     dim = c(n_features, n_features, length(samples)),
                     dimnames = list(feature_names, feature_names, samples))

  for (samp in samples) {
    idx <- which(sample_vec == samp)
    if (length(idx) < min_cells) next
    cor_array[, , samp] <- cor(score_mat[idx, , drop = FALSE], method = "spearman")
  }

  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(NA_real_, n_features, n_features,
                     dimnames = list(feature_names, feature_names))
  p_vals <- matrix(NA_real_, n_features, n_features,
                   dimnames = list(feature_names, feature_names))

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

cat("Loading epithelial object for sample labels and merged UCell scoring...\n")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
sample_vec_all <- tmdata_all$orig.ident[match(rownames(refined_ucell), Cells(tmdata_all))]
sample_vec_all[is.na(sample_vec_all)] <- rownames(refined_ucell)[is.na(sample_vec_all)]

cat("Computing full refined sub-MP correlation matrix for merge decisions...\n")
full_refined_cor <- compute_fisher_cor_p(refined_ucell[, names(refined_mp_genes), drop = FALSE],
                                         sample_vec_all, min_cells = 10)
full_mean_rho <- full_refined_cor$mean_rho

merge_decision_rows <- list()
merge_plan <- list()

for (mp in split_mps) {
  res <- split_results[[mp]]
  if (is.null(res)) next

  if (isTRUE(res$skipped)) {
    merge_plan[[mp]] <- list(
      merged_name = NA_character_,
      merge_labels = character(0),
      retained_labels = mp,
      feature_order = mp
    )
    merge_decision_rows[[mp]] <- data.frame(
      original_mp = mp,
      submp = mp,
      n_other_submps = 0L,
      n_cor_gt_threshold = 0L,
      fraction_cor_gt_threshold = NA_real_,
      qualified_for_merge = FALSE,
      merged_feature = mp,
      action = "skipped_split_kept_parent",
      stringsAsFactors = FALSE
    )
    next
  }

  sub_labels <- sort(unique(as.character(res$sub_labels)))
  missing_cor <- setdiff(sub_labels, rownames(full_mean_rho))
  if (length(missing_cor) > 0) {
    stop("Missing refined UCell correlations for: ", paste(missing_cor, collapse = ", "))
  }

  qualified <- setNames(rep(FALSE, length(sub_labels)), sub_labels)
  n_other <- length(sub_labels) - 1L
  n_gt <- setNames(rep(0L, length(sub_labels)), sub_labels)
  frac_gt <- setNames(rep(NA_real_, length(sub_labels)), sub_labels)

  for (submp in sub_labels) {
    other_submps <- setdiff(sub_labels, submp)
    vals <- full_mean_rho[submp, other_submps]
    n_gt[submp] <- sum(is.finite(vals) & vals > cor_threshold)
    frac_gt[submp] <- if (n_other > 0) n_gt[submp] / n_other else NA_real_
    qualified[submp] <- is.finite(frac_gt[submp]) && frac_gt[submp] >= fraction_threshold
  }

  merge_labels <- names(qualified)[qualified]
  if (length(merge_labels) >= 2) {
    merged_name <- paste0(mp, "+")
    retained_labels <- setdiff(sub_labels, merge_labels)
    action <- ifelse(sub_labels %in% merge_labels,
                     "merge_to_parent_plus", "retain_submp")
    merged_feature <- ifelse(sub_labels %in% merge_labels, merged_name, sub_labels)
    feature_order <- c(merged_name, retained_labels)
  } else {
    merged_name <- NA_character_
    merge_labels <- character(0)
    retained_labels <- sub_labels
    action <- rep("retain_submp_no_parent_merge", length(sub_labels))
    merged_feature <- sub_labels
    feature_order <- sub_labels
  }

  merge_plan[[mp]] <- list(
    merged_name = merged_name,
    merge_labels = merge_labels,
    retained_labels = retained_labels,
    feature_order = feature_order
  )

  merge_decision_rows[[mp]] <- data.frame(
    original_mp = mp,
    submp = sub_labels,
    n_other_submps = n_other,
    n_cor_gt_threshold = as.integer(n_gt[sub_labels]),
    fraction_cor_gt_threshold = as.numeric(frac_gt[sub_labels]),
    qualified_for_merge = as.logical(qualified[sub_labels]),
    merged_feature = merged_feature,
    action = action,
    stringsAsFactors = FALSE
  )
}

merge_decisions <- bind_rows(merge_decision_rows)
write.csv(merge_decisions,
          file.path(outdir, "tables", "merged_refined_mp_merge_decisions.csv"),
          row.names = FALSE)
cat("Saved: tables/merged_refined_mp_merge_decisions.csv\n")

normVector <- function(x) {
  s <- sum(x)
  if (s > 0) {
    x <- x / s
  }
  x
}

weightCumul <- function(vector, weight.explained = 0.5) {
  x.sorted <- sort(vector, decreasing = TRUE)
  norm.x <- normVector(x.sorted)
  cs <- cumsum(norm.x)
  norm.x[cs < weight.explained]
}

wgtLoad <- function(mat, w) {
  rownorm <- apply(mat, 1, normVector)
  spec <- apply(rownorm, 2, max)
  spec.w <- spec^w
  mat <- mat * spec.w
  apply(mat, 2, normVector)
}

derive_consensus_signature <- function(progs, gene.table, nmf.genes.single,
                                       weight_explained = 0.5,
                                       max_genes = 200,
                                       min_confidence = 0.5) {
  progs <- progs[progs %in% colnames(gene.table)]
  if (length(progs) == 0) {
    return(list(genes = character(0), weights = numeric(0)))
  }

  sub_gt <- gene.table[, progs, drop = FALSE]
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

  genes.pass <- weightCumul(genes.avg, weight.explained = weight_explained)
  this <- nmf.genes.single[progs[progs %in% names(nmf.genes.single)]]
  if (length(this) > 0) {
    genes.only <- lapply(this, names)
    genes.sum <- sort(table(unlist(genes.only)), decreasing = TRUE)
    genes.conf <- genes.sum / length(this)
    genes.conf <- genes.conf[genes.conf > min_confidence]
    genes.pass <- genes.pass[names(genes.pass) %in% names(genes.conf)]
  }
  genes.pass <- head(genes.pass, min(length(genes.pass), max_genes))

  if (length(genes.pass) > 0 && sum(genes.pass) > 0) {
    list(genes = names(genes.pass), weights = genes.pass / sum(genes.pass))
  } else {
    list(genes = character(0), weights = numeric(0))
  }
}

cached_genes <- file.path(outdir, "intermediate", "merged_refined_mp_genes.rds")
cached_weights <- file.path(outdir, "intermediate", "merged_refined_mp_gene_weights.rds")
cached_assignments <- file.path(outdir, "intermediate", "merged_refined_mp_assignments.rds")

merged_mp_genes <- NULL
merged_mp_gene_weights <- NULL
gene_cache_valid <- FALSE
if (file.exists(cached_genes) && file.exists(cached_weights)) {
  merged_mp_genes_cached <- readRDS(cached_genes)
  gene_cache_valid <- identical(attr(merged_mp_genes_cached, "algorithm_version"),
                                algorithm_version)
  if (gene_cache_valid) {
    merged_mp_genes <- merged_mp_genes_cached
    merged_mp_gene_weights <- readRDS(cached_weights)
  }
  rm(merged_mp_genes_cached)
}

if (!replot_only && (force_rebuild || !gene_cache_valid)) {
  cat("Building merged refined gene lists...\n")
  plus_features <- unlist(lapply(merge_plan, function(x) x$merged_name),
                          use.names = FALSE)
  plus_features <- plus_features[!is.na(plus_features)]

  gene.table <- NULL
  nmf.genes.single <- NULL
  if (length(plus_features) > 0) {
    cat("Loading raw NMF programs and recomputing weighted loadings...\n")
    geneNMF.programs <- readRDS("geneNMF_outs.rds")
    nmf.wgt <- lapply(geneNMF.programs, function(model) wgtLoad(model$w, w = 5))

    gene.table <- tryCatch({
      processed <- lapply(names(nmf.wgt), function(n) {
        g <- nmf.wgt[[n]]
        colnames(g) <- paste(n, seq_len(ncol(g)), sep = ".")
        g
      })
      Reduce(cbind, processed)
    }, error = function(e) {
      cat("  Standard cbind failed; building with gene-union fallback.\n")
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
  }

  merged_mp_genes <- list()
  merged_mp_gene_weights <- list()

  for (mp in keep_mps) {
    if (mp %in% names(refined_mp_genes)) {
      merged_mp_genes[[mp]] <- refined_mp_genes[[mp]]
      merged_mp_gene_weights[[mp]] <- refined_mp_gene_weights[[mp]]
    }
  }

  for (mp in split_mps) {
    plan <- merge_plan[[mp]]
    if (is.null(plan)) next

    for (feature in plan$feature_order) {
      if (grepl("\\+$", feature)) {
        res <- split_results[[mp]]
        progs <- names(res$sub_labels)[as.character(res$sub_labels) %in% plan$merge_labels]
        consensus <- derive_consensus_signature(progs, gene.table, nmf.genes.single)
        merged_mp_genes[[feature]] <- consensus$genes
        merged_mp_gene_weights[[feature]] <- consensus$weights
        cat("  ", feature, ":", length(consensus$genes), "genes from",
            length(progs), "NMF programs\n")
      } else if (feature %in% names(refined_mp_genes)) {
        merged_mp_genes[[feature]] <- refined_mp_genes[[feature]]
        merged_mp_gene_weights[[feature]] <- refined_mp_gene_weights[[feature]]
      }
    }
  }

  attr(merged_mp_genes, "algorithm_version") <- algorithm_version
  attr(merged_mp_gene_weights, "algorithm_version") <- algorithm_version
  saveRDS(merged_mp_genes, cached_genes)
  saveRDS(merged_mp_gene_weights, cached_weights)
  cat("Saved:", cached_genes, "\n")
  cat("Saved:", cached_weights, "\n")

  cl_members <- geneNMF.metaprograms$programs.clusters
  assignment_rows <- list()
  for (mp in setdiff(mp_names_all, remove_mps)) {
    mp_num <- as.integer(gsub("MP", "", mp))
    member_programs <- names(cl_members)[which(cl_members == mp_num)]
    if (mp %in% split_mps && !is.null(split_results[[mp]]) &&
        !isTRUE(split_results[[mp]]$skipped)) {
      res <- split_results[[mp]]
      old_submp <- as.character(res$sub_labels[member_programs])
      plan <- merge_plan[[mp]]
      new_feature <- ifelse(old_submp %in% plan$merge_labels,
                            plan$merged_name, old_submp)
      assignment_rows[[mp]] <- data.frame(
        program = member_programs,
        original_mp = mp,
        refined_submp = old_submp,
        merged_refined_mp = new_feature,
        stringsAsFactors = FALSE
      )
    } else {
      assignment_rows[[mp]] <- data.frame(
        program = member_programs,
        original_mp = mp,
        refined_submp = mp,
        merged_refined_mp = mp,
        stringsAsFactors = FALSE
      )
    }
  }
  merged_assignments <- bind_rows(assignment_rows)
  saveRDS(merged_assignments, cached_assignments)
  cat("Saved:", cached_assignments, "\n")

  if (exists("geneNMF.programs")) rm(geneNMF.programs)
  if (exists("nmf.wgt")) rm(nmf.wgt)
  if (exists("gene.table")) rm(gene.table)
  if (exists("nmf.genes.single")) rm(nmf.genes.single)
  invisible(gc())
} else if (gene_cache_valid) {
  cat("Loading cached merged refined gene lists...\n")
} else {
  stop("Replot-only mode requested but merged gene caches are missing or outdated.")
}

merged_gene_size_table <- data.frame(
  merged_refined_mp = names(merged_mp_genes),
  parent_mp = parent_id(names(merged_mp_genes)),
  is_parent_plus_merge = grepl("\\+$", names(merged_mp_genes)),
  n_genes = vapply(merged_mp_genes, length, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(merged_gene_size_table,
          file.path(outdir, "tables", "merged_refined_mp_gene_sizes.csv"),
          row.names = FALSE)
cat("Saved: tables/merged_refined_mp_gene_sizes.csv\n")

cached_ucell <- file.path(outdir, "intermediate", "merged_refined_ucell_scores.rds")
merged_ucell_signature <- make_gene_signature(merged_mp_genes)
merged_ucell <- NULL
ucell_cache_valid <- FALSE
if (file.exists(cached_ucell)) {
  merged_ucell_cached <- readRDS(cached_ucell)
  ucell_cache_valid <- identical(attr(merged_ucell_cached, "gene_signature"),
                                 merged_ucell_signature) &&
    all(names(merged_mp_genes) %in% colnames(merged_ucell_cached))
  if (ucell_cache_valid) {
    merged_ucell <- merged_ucell_cached
  }
  rm(merged_ucell_cached)
}

if (force_rebuild || !ucell_cache_valid) {
  cat("Building merged refined UCell score matrix...\n")
  new_score_features <- setdiff(names(merged_mp_genes), colnames(refined_ucell))
  if (length(new_score_features) > 0) {
    tmdata_all <- AddModuleScore_UCell(
      tmdata_all,
      features = merged_mp_genes[new_score_features],
      ncores = 1,
      name = ""
    )
    missing_new <- setdiff(new_score_features, colnames(tmdata_all@meta.data))
    if (length(missing_new) > 0) {
      stop("UCell did not return expected merged score columns: ",
           paste(missing_new, collapse = ", "))
    }
  }

  merged_ucell <- matrix(NA_real_, nrow = nrow(refined_ucell),
                         ncol = length(merged_mp_genes),
                         dimnames = list(rownames(refined_ucell),
                                         names(merged_mp_genes)))
  reused_features <- intersect(names(merged_mp_genes), colnames(refined_ucell))
  merged_ucell[, reused_features] <- as.matrix(refined_ucell[, reused_features, drop = FALSE])
  if (length(new_score_features) > 0) {
    merged_ucell[, new_score_features] <- as.matrix(
      tmdata_all@meta.data[rownames(refined_ucell), new_score_features, drop = FALSE]
    )
  }
  merged_ucell <- as.data.frame(merged_ucell, check.names = FALSE)
  attr(merged_ucell, "gene_signature") <- merged_ucell_signature
  saveRDS(merged_ucell, cached_ucell)
  cat("Saved:", cached_ucell, "\n")
} else {
  cat("Loading cached merged refined UCell scores...\n")
}

merged_state_groups <- list()
for (state_name in names(state_groups)) {
  state_features <- character(0)
  for (mp in state_groups[[state_name]]) {
    if (mp %in% keep_mps) {
      state_features <- c(state_features, mp)
    } else if (mp %in% split_mps && !is.null(merge_plan[[mp]])) {
      state_features <- c(state_features, merge_plan[[mp]]$feature_order)
    }
  }
  state_features <- state_features[state_features %in% names(merged_mp_genes)]
  if (length(state_features) > 0) {
    merged_state_groups[[state_name]] <- state_features
  }
}

merged_mps_ordered <- unlist(merged_state_groups, use.names = FALSE)

####################
# Build extended MP groups for MP7 submps and cell cycle MPs
cc_mps <- intersect(c("MP1", "MP9"), names(merged_mp_genes))
mp7_features <- character(0)
if ("MP7" %in% split_mps && !is.null(merge_plan[["MP7"]])) {
  mp7_features <- merge_plan[["MP7"]]$feature_order
} else if ("MP7" %in% names(merged_mp_genes)) {
  mp7_features <- "MP7"
}
mp7_features <- mp7_features[mp7_features %in% names(merged_mp_genes)]
mp7_parent <- mp7_features[!grepl("[a-z]$", mp7_features)]
mp7_submps <- mp7_features[grepl("[a-z]$", mp7_features)]

# Hard-coded user-specified ordering for all MPs (heatmaps + enrichment)
user_mp_order <- c(
  "MP7j", "MP9", "MP1", "MP2+", "MP17", "MP8+", "MP10+", "MP14", "MP5+",
  "MP7r", "MP7v", "MP10e", "MP16+", "MP18",
  "MP8c", "MP15c", "MP12c", "MP2v", "MP8e", "MP12a", "MP13", 
  "MP7+", "MP7h", "MP8b", "MP12b", "MP15a", "MP15b"
)
# Keep only features that actually exist in merged_mp_genes
full_all_mps_ordered <- user_mp_order[user_mp_order %in% names(merged_mp_genes)]
# merged_mps_ordered is the subset used for the state-only heatmap
merged_mps_ordered <- full_all_mps_ordered
####################

zero_gene <- names(merged_mp_genes)[vapply(merged_mp_genes, length, integer(1)) == 0]
if (length(zero_gene) > 0) {
  cat("Warning: removing zero-gene features from final plot: ",
      paste(zero_gene, collapse = ", "), "\n")
  merged_mps_ordered <- setdiff(merged_mps_ordered, zero_gene)
  full_all_mps_ordered <- setdiff(full_all_mps_ordered, zero_gene)
  for (state_name in names(merged_state_groups)) {
    merged_state_groups[[state_name]] <- setdiff(merged_state_groups[[state_name]], zero_gene)
  }
  cc_mps <- setdiff(cc_mps, zero_gene)
  mp7_parent <- setdiff(mp7_parent, zero_gene)
  mp7_submps <- setdiff(mp7_submps, zero_gene)
  mp7_features <- setdiff(mp7_features, zero_gene)
}

cat("Final merged refined MP order:\n")
for (state_name in names(merged_state_groups)) {
  cat("  ", state_name, ": ",
      paste(merged_state_groups[[state_name]], collapse = ", "), "\n", sep = "")
}
cat("  Cell Cycle: ", paste(cc_mps, collapse = ", "), "\n", sep = "")
cat("  DNA Damage Repair (MP7): ", paste(mp7_features, collapse = ", "), "\n", sep = "")

sample_vec_merged <- tmdata_all$orig.ident[match(rownames(merged_ucell), Cells(tmdata_all))]
sample_vec_merged[is.na(sample_vec_merged)] <- rownames(merged_ucell)[is.na(sample_vec_merged)]

####################
# Compute correlation on full set (state MPs + CC + MP7)
cat("Computing merged refined final correlation matrix (full set)...\n")
full_cor <- compute_fisher_cor_p(merged_ucell[, full_all_mps_ordered, drop = FALSE],
                                  sample_vec_merged, min_cells = 10)
full_mean_rho_final <- full_cor$mean_rho
full_p_vals <- full_cor$p_values
samples <- full_cor$samples

# Extract state-only subset for backward-compatible outputs
mean_rho <- full_mean_rho_final[merged_mps_ordered, merged_mps_ordered]
p_vals <- full_p_vals[merged_mps_ordered, merged_mps_ordered]
####################

saveRDS(list(mean_rho = mean_rho, p_values = p_vals, samples = samples,
             merge_decisions = merge_decisions),
        file.path(outdir, "intermediate", "merged_refined_mp_correlation_matrices.rds"))
write.csv(mean_rho,
          file.path(outdir, "tables", "merged_refined_mp_correlation_mean_rho.csv"))
cat("Saved: intermediate/merged_refined_mp_correlation_matrices.rds\n")
cat("Saved: tables/merged_refined_mp_correlation_mean_rho.csv\n")

mp_to_state <- unlist(lapply(names(merged_state_groups), function(state_name) {
  setNames(rep(state_name, length(merged_state_groups[[state_name]])),
           merged_state_groups[[state_name]])
}))

merged_display_labels <- setNames(display_label(merged_mps_ordered), merged_mps_ordered)
display_mean_rho <- mean_rho
rownames(display_mean_rho) <- merged_display_labels[rownames(display_mean_rho)]
colnames(display_mean_rho) <- merged_display_labels[colnames(display_mean_rho)]

parent_vec <- factor(parent_id(merged_mps_ordered),
                     levels = intersect(names(parent_mp_colors), unique(parent_id(merged_mps_ordered))))
parent_anno_cols <- parent_mp_colors[levels(parent_vec)]

ha_left <- rowAnnotation(
  `Parent MP` = parent_vec,
  col = list(`Parent MP` = parent_anno_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ha_top <- HeatmapAnnotation(
  `Parent MP` = parent_vec,
  col = list(`Parent MP` = parent_anno_cols),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  show_legend = TRUE
)

n_mps <- length(merged_mps_ordered)
hm_size <- unit(13, "inch")
col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))
cell_fs <- max(6, min(10, round(220 / n_mps)))

pdf(file.path(outdir, "figures", "refined_mp_correlation_heatmap_unsupervised_merged.pdf"),
    width = 20, height = 20, useDingbats = FALSE)

ht_cor_unsup <- Heatmap(display_mean_rho,
  name = paste0("Mean Rho\n(", length(samples), " Samples)"),
  col = col_cor,
  rect_gp = gpar(col = "white", lwd = 0.8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  left_annotation = ha_left,
  top_annotation = ha_top,
  row_title = NULL,
  column_title = "Unsupervised clustering after correlated sub-MP merge",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot = 45,
  width = hm_size,
  height = hm_size,
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[i, j]
    rho <- mean_rho[i, j]
    if (is.na(p) || is.na(rho)) {
      grid.text("NA", x, y, gp = gpar(fontsize = cell_fs, col = "grey50"))
    } else {
      stars <- if (p < 0.001) "\n***" else if (p < 0.01) "\n**" else if (p < 0.05) "\n*" else ""
      grid.text(paste0(round(rho, 2), stars), x, y, gp = gpar(fontsize = cell_fs))
    }
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12))
)

draw(ht_cor_unsup, padding = unit(c(25, 25, 25, 25), "mm"))
dev.off()
cat("Saved: figures/refined_mp_correlation_heatmap_unsupervised_merged.pdf\n")

####################
# Additional correlation heatmap: All merged MPs + MP7 submps + cell cycle MPs
full_plot_features <- full_all_mps_ordered
if (length(full_plot_features) >= 2) {
  cat("Generating full (All MPs + MP7 + CC) correlation heatmap...\n")
  full_plot_rho <- full_mean_rho_final[full_plot_features, full_plot_features]
  full_plot_pvals <- full_p_vals[full_plot_features, full_plot_features]

  full_plot_display_labels <- setNames(display_label(full_plot_features), full_plot_features)
  full_plot_display_rho <- full_plot_rho
  rownames(full_plot_display_rho) <- full_plot_display_labels[rownames(full_plot_display_rho)]
  colnames(full_plot_display_rho) <- full_plot_display_labels[colnames(full_plot_display_rho)]

  full_parent_vec <- factor(parent_id(full_plot_features),
                            levels = intersect(names(parent_mp_colors), unique(parent_id(full_plot_features))))
  full_parent_anno_cols <- parent_mp_colors[levels(full_parent_vec)]

  ha_left_full_plot <- rowAnnotation(
    `Parent MP` = full_parent_vec,
    col = list(`Parent MP` = full_parent_anno_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  ha_top_full_plot <- HeatmapAnnotation(
    `Parent MP` = full_parent_vec,
    col = list(`Parent MP` = full_parent_anno_cols),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    show_legend = TRUE
  )

  n_full_plot <- length(full_plot_features)
  full_plot_cell_fs <- max(6, min(10, round(220 / n_full_plot)))

  # Name updated based on user request (keep same style, put in single heatmap)
  # But we'll name it something distinct as requested: "same style... but with mp7 related mps and cell cycle mps"
  # Let's save it as refined_mp_correlation_heatmap_unsupervised_merged_all.pdf to not overwrite the original
  # wait, user specifically asked: "I mean the additional .../refined_mp_correlation_heatmap_unsupervised_merged.pdf should be .../refined_mp_correlation_heatmap_unsupervised_merged.pdf but adding MP1 Mp9 and five MP 7 into a single heatmap"
  # I'll save to 'refined_mp_correlation_heatmap_unsupervised_merged_all.pdf'
  pdf(file.path(outdir, "figures", "refined_mp_correlation_heatmap_unsupervised_merged_all.pdf"),
      width = 28, height = 24, useDingbats = FALSE)

  ht_full_plot <- Heatmap(full_plot_display_rho,
    name = paste0("Mean Rho\n(", length(samples), " Samples)"),
    col = col_cor,
    rect_gp = gpar(col = "white", lwd = 0.8),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = ha_left_full_plot,
    top_annotation = ha_top_full_plot,
    row_title = NULL,
    column_title = "Ordered: All Merged sub-MPs + MP7 + Cell Cycle",
    column_title_gp = gpar(fontsize = 18, fontface = "bold"),
    row_names_side = "right",
    column_names_side = "bottom",
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_gp = gpar(fontsize = 16, fontface = "bold"),
    row_names_max_width = unit(20, "cm"),
    column_names_rot = 45,
    width = unit(15, "inch"),
    height = unit(15, "inch"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- full_plot_pvals[i, j]
      rho <- full_plot_rho[i, j]
      if (is.na(p) || is.na(rho)) {
        grid.text("NA", x, y, gp = gpar(fontsize = full_plot_cell_fs, col = "grey50"))
      } else {
        stars <- if (p < 0.001) "\n***" else if (p < 0.01) "\n**" else if (p < 0.05) "\n*" else ""
        grid.text(paste0(round(rho, 2), stars), x, y, gp = gpar(fontsize = full_plot_cell_fs))
      }
    },
    heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                                labels_gp = gpar(fontsize = 12))
  )

  draw(ht_full_plot, padding = unit(c(25, 25, 25, 25), "mm"), heatmap_legend_side = "left", annotation_legend_side = "left")
  dev.off()
  cat("Saved: figures/refined_mp_correlation_heatmap_unsupervised_merged_all.pdf\n")

  png(file.path(outdir, "figures", "refined_mp_correlation_heatmap_unsupervised_merged_all.png"), width = 28, height = 24, units = "in", res = 300)
  draw(ht_full_plot, padding = unit(c(25, 25, 25, 25), "mm"), heatmap_legend_side = "left", annotation_legend_side = "left")
  dev.off()
  cat("Saved: figures/refined_mp_correlation_heatmap_unsupervised_merged_all.png\n")
}
####################

summary_table <- data.frame(
  metric = c(
    "cor_threshold",
    "fraction_threshold",
    "n_parent_plus_features",
    "n_retained_split_submps",
    "n_final_plot_features",
    "n_samples_in_correlation"
  ),
  value = c(
    cor_threshold,
    fraction_threshold,
    sum(grepl("\\+$", names(merged_mp_genes))),
    sum(grepl("[a-z]$", names(merged_mp_genes))),
    length(merged_mps_ordered),
    length(samples)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table,
          file.path(project_dir, "updates/new_updates/summaries",
                    "mp_refinement_merge_correlated_submps_summary.csv"),
          row.names = FALSE)
cat("Saved: updates/new_updates/summaries/mp_refinement_merge_correlated_submps_summary.csv\n")

cat("\n=== Running Enrichment Analysis for Merged/Refined MPs ===\n")
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

####################
# Use full_all_mps_ordered (hard-coded user order) for enrichment
final_col_order <- full_all_mps_ordered
heatmap_gaps_col <- 0
####################

mp_gene_lists <- merged_mp_genes[final_col_order]
####################
# Extend display labels for enrichment to cover all features
merged_display_labels <- setNames(display_label(final_col_order), final_col_order)
####################

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
    genes    = genes,
    GO       = res_GO,
    Hallmark = res_H,
    MPs_3CA  = res_M
  )
  
  return(c(base_results, res_custom_list))
})
names(cluster_enrich) <- names(mp_gene_lists)

enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7, 
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {
  
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
    return(invisible(NULL))
  }
  
  if (is_custom) {
    if (!element %in% names(custom_refs)) return(invisible(NULL))
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    terms_use <- df %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::arrange(Program, padj) %>%
      dplyr::group_by(Program) %>%
      dplyr::slice_head(n = top_per_program) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Term) %>%
      dplyr::pull(Term)
    
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
  
  ordered_mps <- final_col_order
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  
  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = if(element %in% c("Hallmark","GO","MPs_3CA") || is_custom) tidyr::replace_na(Overlap, "") else ""
    )
  
  mat <- final_df %>%
    dplyr::select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  
  text_mat <- final_df %>%
    dplyr::select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()
  
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]

  if (nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))

  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  # Format column labels to be on the same line
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
  
  # Insert clear empty column to separate the groups AFTER row sorting
  if (heatmap_gaps_col > 0 && heatmap_gaps_col < ncol(mat)) {
    mat <- cbind(mat[, 1:heatmap_gaps_col, drop=FALSE], " " = rep(0, nrow(mat)), mat[, (heatmap_gaps_col+1):ncol(mat), drop=FALSE])
    text_mat <- cbind(text_mat[, 1:heatmap_gaps_col, drop=FALSE], " " = rep("", nrow(text_mat)), text_mat[, (heatmap_gaps_col+1):ncol(text_mat), drop=FALSE])
    col_labels <- c(col_labels[1:heatmap_gaps_col], "", col_labels[(heatmap_gaps_col+1):length(col_labels)])
  }
  
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  
  hm_name <- paste0("enrich_", element)
  ht <- ComplexHeatmap::pheatmap(mat,
                     name = hm_name,
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
                     angle_col = "45",
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     main = paste0(element, " Enrichment (-log10 padj)"))
  
  ht@column_names_param$rot <- 35
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 35, 2, 2), "mm"))
                     
  target_mps <- c("MP17", "MP7r", "MP8c", "MP8b")
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
pdf(file.path(outdir, "figures", "merged_refined_mp_enrichment_anno.pdf"), width = 14, height = 12)
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

png(file.path(outdir, "figures", "enrich_merged_Hallmark.png"), width = 3500, height = 2100, res = 300)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_GO.png"), width = 4025, height = 2400, res = 300)
enrich_heatmap(cluster_enrich, "GO", top_per_program = 6, top_n = 60, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_MPs_3CA.png"), width = 3500, height = 2200, res = 300)
enrich_heatmap(cluster_enrich, "MPs_3CA", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Early_Embryogenesis.png"), width = 3850, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Normal_Development_long.png"), width = 5075, height = 3400, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Normal_Development_short.png"), width = 5075, height = 3400, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Organogenesis_major.png"), width = 3850, height = 2200, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Organogenesis_sub.png"), width = 4375, height = 2200, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Adult_Epithelium.png"), width = 3938, height = 2300, res = 300)
enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(outdir, "figures", "enrich_merged_Barretts_Oesophagus.png"), width = 3938, height = 2300, res = 300)
enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

####################
# === External Reference UCell Scoring ===
cat("\n=== UCell Scoring in External References ===\n")
score_external_references_for_heatmap <- function() {
  external_cache <- file.path(outdir, "intermediate", "merged_refined_mp_external_ucell_scores.rds")
  if (file.exists(external_cache) && !force_rebuild) return(readRDS(external_cache))

  clean_counts <- function(mat) {
    rn <- gsub("_", "-", rownames(mat))
    is_dup <- duplicated(rn)
    if (any(is_dup)) {
      mat <- mat[!is_dup, , drop=FALSE]
      rownames(mat) <- rn[!is_dup]
    } else {
      rownames(mat) <- rn
    }
    return(mat)
  }

  summaries <- list()
  
  # 1. Stomach.rds
  stomach_path <- file.path("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds")
  if (file.exists(stomach_path)) {
    cat("Scoring adult stomach annotated Seurat object with external UCell logic.\n")
    stomach_rename <- c(
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
    adult_stomach <- readRDS(stomach_path)
    stomach_meta <- adult_stomach@meta.data
    stomach_meta$cell <- Cells(adult_stomach)
    stomach_meta <- stomach_meta %>%
      dplyr::filter(major_clusters == "epi") %>%
      dplyr::transmute(cell = cell, term = unname(stomach_rename[as.character(subcluster.v2)]), dataset_source = "Adult_Stomach") %>%
      dplyr::filter(!is.na(term))
    
    stomach_sampled <- stomach_meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
    stomach_counts <- GetAssayData(adult_stomach, layer = "counts")[, stomach_sampled$cell, drop=FALSE]
    stomach_counts <- clean_counts(stomach_counts)
    usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(stomach_counts)))
    usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]

    obj <- CreateSeuratObject(counts = stomach_counts, meta.data = as.data.frame(stomach_sampled))
    obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
    score_cols <- paste0(names(usable_mp_genes), "_UCell")
    score_df <- cbind(stomach_sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
    colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
    
    summaries[["Adult_Stomach"]] <- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
    rm(adult_stomach, stomach_counts, obj); invisible(gc())
  }

  # 2. Barretts High Quality
  barretts_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/barretts/alldatahighquality.rds"
  if (file.exists(barretts_path)) {
    cat("Scoring Barretts high-quality SingleCellExperiment with external UCell logic.\n")
    barrett_term_map <- c(
      "Basal" = "Basal_Normal_Esophagus..Barretts_Oesophagus", "Suprabasal" = "Suprabasal_Normal_Esophagus..Barretts_Oesophagus",
      "Suprabasal_Dividing" = "Suprabasal_Dividing_Normal_Esophagus..Barretts_Oesophagus", "Intermediate" = "Intermediate_Normal_Esophagus..Barretts_Oesophagus",
      "Superficial" = "Superficial_Normal_Esophagus..Barretts_Oesophagus", "Undifferentiated" = "Undifferentiated_Normal_Gastric..Barretts_Oesophagus",
      "Undifferentiated_Dividing" = "Undifferentiated_Dividing_Normal_Gastric..Barretts_Oesophagus", "Foveolar_Intermediate" = "Foveolar_Intermediate_Normal_Gastric..Barretts_Oesophagus",
      "Foveolar_differentiated" = "Foveolar_differentiated_Normal_Gastric..Barretts_Oesophagus", "Chief" = "Chief_Normal_Gastric..Barretts_Oesophagus",
      "Parietal" = "Parietal_Normal_Gastric..Barretts_Oesophagus", "Endocrine_GHRL" = "Endocrine_GHRL_Normal_Gastric..Barretts_Oesophagus",
      "Endocrine_CHGA" = "Endocrine_CHGA_Normal_Gastric..Barretts_Oesophagus", "Endocrine_NEUROD1" = "Endocrine_NEUROD1_Normal_Gastric..Barretts_Oesophagus",
      "Columnar_Undifferentiated" = "Columnar_Undifferentiated_Barretts_Esophagus..Barretts_Oesophagus", "Columnar_Undifferentiated_Dividing" = "Columnar_Dividing_Barretts_Esophagus..Barretts_Oesophagus",
      "Endocrine_NEUROG3" = "Endocrine_NEUROG3_Barretts_Esophagus..Barretts_Oesophagus", "C1" = "Columnar_Intermediate_Barretts_Esophagus..Barretts_Oesophagus",
      "C2" = "Columnar_differentiated_Barretts_Esophagus..Barretts_Oesophagus", "Goblet" = "Goblet_Barretts_Esophagus..Barretts_Oesophagus",
      "Duct_Intercalating" = "Duct_Intercalating_Submucosal_Glands..Barretts_Oesophagus", "Oncocytes_MUC5B_Low" = "Oncocytes_Submucosal_Glands..Barretts_Oesophagus",
      "Mucous_MUC5B_High" = "Mucous_Submucosal_Glands..Barretts_Oesophagus"
    )
    barretts_sce <- readRDS(barretts_path)
    if ("SummarizedExperiment" %in% loadedNamespaces()) {
      b_counts <- SummarizedExperiment::assay(barretts_sce, "counts")
      b_counts <- as(b_counts, "CsparseMatrix")
      b_meta <- as.data.frame(SummarizedExperiment::colData(barretts_sce), stringsAsFactors = FALSE)
      b_cells <- colnames(barretts_sce)
      if (is.null(b_cells)) b_cells <- rownames(b_meta)
      if (is.null(b_cells)) b_cells <- paste0("Cell_", seq_len(ncol(b_counts)))
      colnames(b_counts) <- b_cells
      b_meta$cell <- b_cells
      b_meta <- b_meta %>% dplyr::transmute(cell = cell, term = unname(barrett_term_map[as.character(cell_type_secondary)]), dataset_source = "Barretts_HighQuality") %>% dplyr::filter(!is.na(term))
      
      b_sampled <- b_meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
      b_counts <- b_counts[, b_sampled$cell, drop=FALSE]
      b_counts <- clean_counts(b_counts)
      usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(b_counts)))
      usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]

      obj <- CreateSeuratObject(counts = b_counts, meta.data = as.data.frame(b_sampled))
      obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
      score_cols <- paste0(names(usable_mp_genes), "_UCell")
      score_df <- cbind(b_sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
      colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
      
      summaries[["Barretts_HighQuality"]] <- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
    }
    rm(barretts_sce); if(exists("b_counts")) rm(b_counts); if(exists("obj")) rm(obj); invisible(gc())
  }
  
  # 3. Adult Oesophagus
  oesophagus_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/SCP1242"
  oes_mtx <- file.path(oesophagus_dir, "expression", "63f5559ca03dda60911f184e", "EoE-processed.mtx")
  if (file.exists(oes_mtx)) {
    cat("Scoring adult oesophagus MatrixMarket counts with external UCell logic.\n")
    oes_meta_path <- file.path(oesophagus_dir, "metadata", "EoE_meta.txt")
    oes_cell_path <- file.path(oesophagus_dir, "expression", "63f5559ca03dda60911f184e", "EoE_cell_processed.tsv")
    oes_gene_path <- file.path(oesophagus_dir, "expression", "63f5559ca03dda60911f184e", "EoE_gene_processed.tsv")
    if (file.exists(oes_meta_path) && file.exists(oes_cell_path) && file.exists(oes_gene_path)) {
      oesophagus_labels <- c("Quiescent basal cell", "Basal cell (cycling)", "Suprabasal", "Apical cell")
      oesophagus_term_map <- setNames(paste0(gsub(" ", "_", oesophagus_labels), "_Oesophagus..Adult_Epithelium"), oesophagus_labels)
      oes_meta <- data.table::fread(oes_meta_path) %>% dplyr::filter(NAME != "TYPE") %>% dplyr::transmute(cell = as.character(NAME), term = unname(oesophagus_term_map[as.character(cell_type_anno)]), dataset_source = "Adult_Oesophagus") %>% dplyr::filter(!is.na(term))
      oes_sampled <- oes_meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
      
      subset_oesophagus_counts <- function(mtx_path, cell_path, gene_path, sampled_meta) {
        all_cells <- data.table::fread(cell_path, header = FALSE, col.names = "cell")$cell
        all_genes <- data.table::fread(gene_path, header = FALSE, col.names = "gene")$gene
        keep_idx <- match(sampled_meta$cell, all_cells)
        index_map <- data.frame(old_index = keep_idx, new_index = seq_along(keep_idx))
        map_path <- file.path(outdir, "intermediate", "Auto_oesophagus_subset_map.tsv")
        triplet_path <- file.path(outdir, "intermediate", "Auto_oesophagus_subset_triplets.tsv")
        filtered_mtx_path <- file.path(outdir, "intermediate", "Auto_oesophagus_subset.mtx")
        write.table(index_map, map_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        awk_program <- "FNR==NR {map[$1]=$2; next} /^%/ {next} !seen_dims {seen_dims=1; next} NF==3 {if ($2 in map) print $1, map[$2], $3}"
        cmd <- sprintf("awk '%s' %s %s > %s", awk_program, shQuote(map_path), shQuote(mtx_path), shQuote(triplet_path))
        system(cmd)
        n_entries <- as.integer(scan(text = system2("wc", c("-l", triplet_path), stdout = TRUE), what = character(), quiet = TRUE)[1])
        writeLines(c("%%MatrixMarket matrix coordinate integer general", paste(length(all_genes), length(keep_idx), n_entries)), filtered_mtx_path)
        file.append(filtered_mtx_path, triplet_path)
        counts <- Matrix::readMM(filtered_mtx_path)
        counts <- as(counts, "CsparseMatrix")
        rownames(counts) <- all_genes
        colnames(counts) <- sampled_meta$cell
        return(counts)
      }
      
      oes_counts <- subset_oesophagus_counts(oes_mtx, oes_cell_path, oes_gene_path, oes_sampled)
      oes_counts <- clean_counts(oes_counts)
      usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(oes_counts)))
      usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]

      obj <- CreateSeuratObject(counts = oes_counts, meta.data = as.data.frame(oes_sampled))
      obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
      score_cols <- paste0(names(usable_mp_genes), "_UCell")
      score_df <- cbind(oes_sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
      colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
      
      summaries[["Adult_Oesophagus"]] <- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
      rm(oes_counts, obj); invisible(gc())
    }
  }

  # 4. Early Embryogenesis
  early_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/early_embryogenesis/psd.R3.6.em.seurat.ob.rds"
  if (file.exists(early_path)) {
    cat("Scoring Early Embryogenesis...\n")
    early_obj <- readRDS(early_path)
    early_meta <- early_obj@meta.data
    early_meta$cell <- rownames(early_meta)
    early_label <- as.character(early_meta$rename_EML)
    early_label[early_label == "2-4 cell"] <- "Z4cell"
    early_terms <- as.character(custom_refs[["Early_Embryogenesis"]]$TERM2NAME$term)
    early_meta <- early_meta %>%
      dplyr::mutate(term = paste0(gsub(" ", "_", early_label), "..Early_Embryogenesis")) %>%
      dplyr::transmute(cell = cell, term = term, dataset_source = "Early_Embryogenesis_integrated") %>%
      dplyr::filter(term %in% early_terms)
    
    if (nrow(early_meta) > 0) {
      early_sampled <- early_meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
      early_counts <- GetAssayData(early_obj, layer = "counts")[, early_sampled$cell, drop=FALSE]
      early_counts <- clean_counts(early_counts)
      usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(early_counts)))
      usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]
      if(length(usable_mp_genes) > 0) {
        obj <- CreateSeuratObject(counts = early_counts, meta.data = as.data.frame(early_sampled))
        obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
        score_cols <- paste0(names(usable_mp_genes), "_UCell")
        score_df <- cbind(early_sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
        colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
        summaries[["Early_Embryogenesis_integrated"]] <- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
        rm(obj, early_counts)
      }
    }
    rm(early_obj); invisible(gc())
  }
  
  # 5. Normal Development Stomach
  normal_counts_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/Stomach_gene_count.RDS"
  normal_cell_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_cell.RDS"
  normal_gene_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_gene.RDS"
  if (file.exists(normal_counts_path) && file.exists(normal_cell_path) && file.exists(normal_gene_path)) {
    cat("Scoring Normal Development Stomach...\n")
    normal_counts <- readRDS(normal_counts_path)
    normal_counts <- as(normal_counts, "CsparseMatrix")
    normal_cell <- readRDS(normal_cell_path)
    normal_gene <- readRDS(normal_gene_path)
    normal_symbols <- normal_gene$gene_short_name[match(rownames(normal_counts), normal_gene$gene_id)]
    normal_symbols[is.na(normal_symbols) | normal_symbols == ""] <- rownames(normal_counts)[is.na(normal_symbols) | normal_symbols == ""]
    rownames(normal_counts) <- normal_symbols
    
    make_normal_meta <- function(df_cell, count_cells, ref_name, dataset_source) {
      stomach_terms <- grep("_Stomach\\.\\.", as.character(custom_refs[[ref_name]]$TERM2NAME$term), value = TRUE)
      label_to_term <- setNames(stomach_terms, tolower(gsub("[^A-Za-z0-9]+", "", gsub("_", " ", sub("_Stomach\\.\\..*$", "", stomach_terms)))))
      df_cell %>%
        dplyr::filter(sample %in% count_cells, Organ == "Stomach", !is.na(Main_cluster_name)) %>%
        dplyr::transmute(
          cell = as.character(sample),
          term = unname(label_to_term[tolower(gsub("[^A-Za-z0-9]+", "", Main_cluster_name))]),
          dataset_source = dataset_source
        ) %>%
        dplyr::filter(!is.na(term))
    }
    
    normal_long_meta <- make_normal_meta(normal_cell, colnames(normal_counts), "Normal_Development_long", "Normal_Development_Stomach_long")
    normal_short_meta <- make_normal_meta(normal_cell, colnames(normal_counts), "Normal_Development_short", "Normal_Development_Stomach_short")
    
    score_normal <- function(meta, ds) {
      if (nrow(meta) > 0) {
        sampled <- meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
        sub_counts <- normal_counts[, sampled$cell, drop=FALSE]
        sub_counts <- clean_counts(sub_counts)
        usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(sub_counts)))
        usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]
        if(length(usable_mp_genes) > 0) {
          obj <- CreateSeuratObject(counts = sub_counts, meta.data = as.data.frame(sampled))
          obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
          score_cols <- paste0(names(usable_mp_genes), "_UCell")
          score_df <- cbind(sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
          colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
          summaries[[ds]] <<- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
          rm(obj, sub_counts); invisible(gc())
        }
      }
    }
    score_normal(normal_long_meta, "Normal_Development_Stomach_long")
    score_normal(normal_short_meta, "Normal_Development_Stomach_short")
    rm(normal_counts, normal_cell, normal_gene); invisible(gc())
  }
  
  # 6. Organogenesis
  org_meta_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_cell_annotate.txt.gz"
  org_gene_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_gene_annotate.txt.gz"
  org_mtx_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_raw_counts.mtx.gz"
  if (file.exists(org_meta_path) && file.exists(org_gene_path) && file.exists(org_mtx_path)) {
    cat("Scoring Organogenesis...\n")
    org_meta <- data.table::fread(org_meta_path)
    org_genes <- data.table::fread(org_gene_path)
    
    major_terms <- as.character(custom_refs[["Organogenesis_major"]]$TERM2NAME$term)
    sub_terms <- as.character(custom_refs[["Organogenesis_sub"]]$TERM2NAME$term)
    
    org_meta <- org_meta %>%
      dplyr::mutate(
        cell = as.character(cell_id),
        major_term = paste0(gsub(" ", "_", as.character(`developmental system`)), "..Organogenesis_major"),
        sub_base = as.character(final_annotation),
        sub_term = ifelse(
          grepl("epithelium", sub_base, ignore.case = TRUE),
          paste0(gsub(" ", "_", sub_base), "..Organogenesis_sub"),
          paste0(gsub(" ", "_", sub_base), "_Endoderm..Organogenesis_sub")
        )
      )
      
    org_major_meta <- org_meta %>% dplyr::transmute(cell = cell, term = major_term, dataset_source = "Organogenesis_GSE157329_major") %>% dplyr::filter(term %in% major_terms)
    org_sub_meta <- org_meta %>% dplyr::transmute(cell = cell, term = sub_term, dataset_source = "Organogenesis_GSE157329_sub") %>% dplyr::filter(term %in% sub_terms)
    
    if (nrow(org_major_meta) > 0 || nrow(org_sub_meta) > 0) {
      org_counts <- Matrix::readMM(gzfile(org_mtx_path))
      org_counts <- as(org_counts, "CsparseMatrix")
      if (nrow(org_counts) == nrow(org_genes)) {
        rownames(org_counts) <- org_genes$gene_short_name
        colnames(org_counts) <- org_meta$cell
      } else {
        org_counts <- t(org_counts)
        rownames(org_counts) <- org_genes$gene_short_name
        colnames(org_counts) <- org_meta$cell
      }
      
      score_org <- function(meta, ds) {
        if (nrow(meta) > 0) {
          sampled <- meta %>% dplyr::group_by(term) %>% dplyr::slice_sample(n = 5000) %>% dplyr::ungroup()
          sub_counts <- org_counts[, sampled$cell, drop=FALSE]
          sub_counts <- clean_counts(sub_counts)
          usable_mp_genes <- lapply(mp_gene_lists, function(genes) intersect(genes, rownames(sub_counts)))
          usable_mp_genes <- usable_mp_genes[lengths(usable_mp_genes) > 0]
          if(length(usable_mp_genes) > 0) {
            obj <- CreateSeuratObject(counts = sub_counts, meta.data = as.data.frame(sampled))
            obj <- UCell::AddModuleScore_UCell(obj, features = usable_mp_genes, slot = "counts", name = "_UCell", ncores = 1)
            score_cols <- paste0(names(usable_mp_genes), "_UCell")
            score_df <- cbind(sampled[, c("cell", "term", "dataset_source")], obj@meta.data[, score_cols, drop=FALSE])
            colnames(score_df)[match(score_cols, colnames(score_df))] <- names(usable_mp_genes)
            summaries[[ds]] <<- score_df %>% dplyr::group_by(dataset_source, term) %>% dplyr::summarise(n_cells_scored = dplyr::n(), dplyr::across(dplyr::all_of(names(usable_mp_genes)), mean), .groups = "drop")
            rm(obj, sub_counts); invisible(gc())
          }
        }
      }
      score_org(org_major_meta, "Organogenesis_GSE157329_major")
      score_org(org_sub_meta, "Organogenesis_GSE157329_sub")
      rm(org_counts); invisible(gc())
    }
  }

  out <- dplyr::bind_rows(summaries)
  saveRDS(out, external_cache)
  out
}
external_summary <- score_external_references_for_heatmap()

ucell_heatmap <- function(external_summary, element, cols = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100), fontsize_row = 7, fontsize_col = 9) {
  if (!element %in% names(custom_refs)) {
    mat <- matrix(NA, nrow = 1, ncol = length(final_col_order), dimnames = list("No single-cell reference available", final_col_order))
    terms_use <- rownames(mat)
    mp_sizes <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
    col_labels <- paste0(merged_display_labels[colnames(mat)], " (n=", mp_sizes, ")")
    if (heatmap_gaps_col > 0 && heatmap_gaps_col < ncol(mat)) {
      mat <- cbind(mat[, 1:heatmap_gaps_col, drop=FALSE], " " = rep(NA, nrow(mat)), mat[, (heatmap_gaps_col+1):ncol(mat), drop=FALSE])
      col_labels <- c(col_labels[1:heatmap_gaps_col], "", col_labels[(heatmap_gaps_col+1):length(col_labels)])
    }
    text_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
    hm_name <- paste0("ucell_", element)
    ht <- ComplexHeatmap::pheatmap(mat, name = hm_name, display_numbers = text_mat, number_color = "black", fontsize_number = fontsize_row * 1.1, labels_col = col_labels, color = cols, na_col = "#F0F0F0", breaks = seq(0, 1, length.out = length(cols) + 1), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, show_colnames = TRUE, angle_col = "45", fontsize_row = 12, fontsize_col = fontsize_col, main = paste0(element, " Mean UCell in Reference Cells\n(No single-cell data)"))
    ht@column_names_param$rot <- 35
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 35, 2, 2), "mm"))
    
    target_mps <- c("MP17", "MP7r", "MP8c", "MP8b")
    num_slices <- 1
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
  
  terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  df <- if (!is.null(external_summary) && nrow(external_summary) > 0) external_summary %>% dplyr::filter(term %in% terms_use) else data.frame()
  
  mat <- matrix(NA_real_, nrow = length(terms_use), ncol = length(final_col_order), dimnames = list(terms_use, final_col_order))
  
  if (nrow(df) > 0) {
    mat_df <- df %>% dplyr::select(term, dplyr::all_of(intersect(final_col_order, colnames(df)))) %>% as.data.frame()
    rownames(mat_df) <- mat_df$term
    mat_df$term <- NULL
    common_cols <- intersect(final_col_order, colnames(mat_df))
    common_rows <- intersect(terms_use, rownames(mat_df))
    if (length(common_rows) > 0 && length(common_cols) > 0) {
      mat[common_rows, common_cols] <- as.matrix(mat_df[common_rows, common_cols])
    }
    text_mat <- matrix(ifelse(!is.na(mat) & mat > 0, sprintf("%.2f", mat), ""), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    main_title <- paste0(element, " Mean UCell in Reference Cells")
  } else {
    mat[] <- NA
    text_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    main_title <- paste0(element, " Mean UCell in Reference Cells\n(No expression cells available)")
  }
  
  mp_sizes <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
  col_labels <- paste0(merged_display_labels[colnames(mat)], " (n=", mp_sizes, ")")
  
  if (heatmap_gaps_col > 0 && heatmap_gaps_col < ncol(mat)) {
    mat <- cbind(mat[, 1:heatmap_gaps_col, drop=FALSE], " " = rep(NA, nrow(mat)), mat[, (heatmap_gaps_col+1):ncol(mat), drop=FALSE])
    text_mat <- cbind(text_mat[, 1:heatmap_gaps_col, drop=FALSE], " " = rep("", nrow(text_mat)), text_mat[, (heatmap_gaps_col+1):ncol(text_mat), drop=FALSE])
    col_labels <- c(col_labels[1:heatmap_gaps_col], "", col_labels[(heatmap_gaps_col+1):length(col_labels)])
  }
  
  cap <- quantile(mat[mat>0], 0.98, na.rm=TRUE)
  if (is.na(cap) || cap <= 0) cap <- max(mat, na.rm=TRUE)
  if (is.na(cap) || cap <= 0) cap <- 1
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  
  hm_name <- paste0("ucell_", element)
  ht <- ComplexHeatmap::pheatmap(mat,
                     name = hm_name,
                     display_numbers = text_mat,
                     number_color = "black",
                     fontsize_number = fontsize_row * 1.1,
                     labels_col = col_labels,
                     color = cols,
                     na_col = "#F0F0F0",
                     breaks = breaks,
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE,
                     border_color = NA,
                     show_colnames = TRUE,
                     angle_col = "45",
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     main = main_title)
                     
  ht@column_names_param$rot <- 35
  ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 35, 2, 2), "mm"))
  
  target_mps <- c("MP17", "MP7r", "MP8c", "MP8b")
  num_slices <- 1
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

cat("Saving UCell heatmaps to figures/...\n")
pdf(file.path(outdir, "figures", "merged_refined_mp_ucell_reference_anno.pdf"), width = 14, height = 12)
elements_to_plot <- c("Early_Embryogenesis", "Normal_Development_long", "Normal_Development_short", "Organogenesis_major", "Organogenesis_sub", "Adult_Epithelium", "Barretts_Oesophagus")
for (element in elements_to_plot) {
  ucell_heatmap(external_summary, element, cols = cols_palette)
}
dev.off()

for (ref in names(custom_refs)) {
  df <- external_summary %>% dplyr::filter(term %in% as.character(custom_refs[[ref]]$TERM2NAME$term))
  if (nrow(df) > 0) {
    png(file.path(outdir, "figures", paste0("ucell_merged_", ref, ".png")), width = 3938, height = 2300, res = 300)
    ucell_heatmap(external_summary, ref, cols = cols_palette)
    dev.off()
  }
}

####################

cat("\n=== Generating Excel MP Gene Summary ===\n")

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
  "MP12" = "Neuro-responsive Epi."
)

get_desc <- function(mp_name) {
  if (mp_name %in% names(submp_desc_map)) return(submp_desc_map[mp_name])
  if (grepl("[a-z]$", mp_name)) return("")
  parent <- parent_id(mp_name)
  if (parent %in% names(mp_descriptions)) return(mp_descriptions[parent])
  return(paste0(parent, "_unknown"))
}

# full_all_mps_ordered already contains all features including MP7 submps and CC
all_merged_names <- names(merged_mp_genes)

# Use user-specified order for Excel pages too
all_ordered_excel <- full_all_mps_ordered

# Page 1: Split_Separated — split submps first, gap, then merged plus MPs
split_submps_excel <- all_ordered_excel[grepl("[a-z]$", all_ordered_excel)]
merged_plus_mps <- all_ordered_excel[grepl("\\+$", all_ordered_excel)]
page1_order <- c(split_submps_excel, "GAP", merged_plus_mps)

# Page 2: Grouped_By_Parent — group by parent with gaps between, following user order
seen_parents <- character(0)
page2_order <- character(0)
for (feat in all_ordered_excel) {
  p <- parent_id(feat)
  if (!p %in% seen_parents) {
    if (length(seen_parents) > 0) page2_order <- c(page2_order, "GAP")
    seen_parents <- c(seen_parents, p)
  }
  page2_order <- c(page2_order, feat)
}

build_mp_matrix <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)

  get_genes_xl <- function(mp) {
    if (mp == "GAP") return(character(0))
    merged_mp_genes[[mp]]
  }

  max_g <- max(sapply(mp_names_vec, function(x) length(get_genes_xl(x))))
  n_mp <- length(mp_names_vec)
  n_rows <- max_g + 2

  mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
  for (i in seq_along(mp_names_vec)) {
    mp <- mp_names_vec[i]
    if (mp == "GAP") {
      mat[1, i] <- ""
      mat[2, i] <- ""
    } else {
      mat[1, i] <- mp
      mat[2, i] <- get_desc(mp)
      genes <- get_genes_xl(mp)
      if (length(genes) > 0) {
        mat[3:(length(genes)+2), i] <- genes
      }
    }
  }

  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

df_p1 <- build_mp_matrix(page1_order)
df_p2 <- build_mp_matrix(page2_order)

wb <- createWorkbook()
mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
desc_style <- createStyle(fgFill = "#F2F2F2")

add_sheet <- function(wb, sheet_name, df, order_vec) {
  addWorksheet(wb, sheet_name)
  sheet_idx <- length(names(wb))

  if (!is.null(df)) {
    writeData(wb, sheet = sheet_idx, x = df, startCol = 1, startRow = 1, colNames = FALSE)

    for (i in seq_along(order_vec)) {
      if (order_vec[i] == "GAP") {
        setColWidths(wb, sheet_idx, cols = i, widths = 3)
      } else {
        setColWidths(wb, sheet_idx, cols = i, widths = 25)
        addStyle(wb, sheet = sheet_idx, mp_name_style, rows = 1, cols = i, gridExpand = TRUE)
        addStyle(wb, sheet = sheet_idx, desc_style, rows = 2, cols = i, gridExpand = TRUE)
      }
    }
  }
}

add_sheet(wb, "Split_Separated", df_p1, page1_order)
add_sheet(wb, "Grouped_By_Parent", df_p2, page2_order)

xlsx_path <- file.path(outdir, "tables", "merged_refined_MP_genes_summary.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("Saved:", xlsx_path, "\n")
####################

cat("\n=== Correlated Sub-MP Merge Complete ===\n")
cat("Merged parent-plus features:",
    paste(names(merged_mp_genes)[grepl("\\+$", names(merged_mp_genes))],
          collapse = ", "), "\n")
cat("Final state-associated features:", length(merged_mps_ordered), "\n")
cat("Full features (incl. CC + MP7):", length(full_all_mps_ordered), "\n")
cat("Done.\n")
