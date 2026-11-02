#!/usr/bin/env Rscript
####################
# Analysis registry:
# Status: active
# Script: analysis/spatial/xenium_hd_ucell_scores.R
# Description: Score all 17 malignant MPs and the exact saved HD TME signatures
#   on raw Xenium counts using the HD UCell definition and parent restrictions.
# Methodology: analysis/methodology/spatial/xenium_visiumhd_meta_methodology.md
# Inputs: ref_outs/Yates_spatial/source_data/Zenodo_upload/processed/Xenium.h5ad;
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds;
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv;
#   ref_outs/visium_hd_outs/spatial_interactions/intermediate/Auto_visiumhd_tme_mp_gene_sets.rds;
#   ref_outs/xenium_scatlas_states/spatial_neighbourhoods/intermediate/Auto_xenium_metadata.csv.gz.
# Outputs: ref_outs/xenium_scatlas_states/spatial_neighbourhoods/
#   intermediate/Auto_<patient>_ucell_scores.csv.gz;
#   tables/Auto_<patient>_ucell_feature_audit.csv, Auto_signature_coverage.csv;
#   logs/Auto_ucell_session_info.txt.
# Cache/replot: called only for full rebuild; persistent scores support Python
#   contact-only reruns. SCREF_FORCE_REBUILD rebuilds scores; no hidden cache.
# Run: /opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE job_scripts/spatial/spatial_platform_neighbourhoods.sh
# Direct run on compute node: Rscript analysis/spatial/xenium_hd_ucell_scores.R
# Resources/environment: eight CPUs, 128 GB RAM, 24 hours; conda dmtcp.
# Upstream prepare stage writes identity-aligned metadata before this scorer runs.
####################
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(rhdf5)
  library(UCell)
})

####################
wd <- '/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline'
setwd(wd)
out <- file.path(wd, 'ref_outs/xenium_scatlas_states/spatial_neighbourhoods')
h5 <- file.path(wd, 'ref_outs/Yates_spatial/source_data/Zenodo_upload/processed/Xenium.h5ad')
meta <- fread(file.path(out, 'intermediate/Auto_xenium_metadata.csv.gz'))
mp_genes <- readRDS('ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds')
grouping <- fread('ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv')
if (nrow(grouping) != 17L || anyDuplicated(grouping$mp) || anyDuplicated(names(mp_genes)) ||
    any(grouping$mp %in% c("MP2x", "MP11c", "MP18a")) || !setequal(names(mp_genes), grouping$mp)) stop('Current 17-MP contract mismatch')
tme <- readRDS('ref_outs/visium_hd_outs/spatial_interactions/intermediate/Auto_visiumhd_tme_mp_gene_sets.rds')
if (!all(c('gene_sets', 'feature_table') %in% names(tme)) ||
    anyDuplicated(tme$feature_table$feature_id) ||
    !identical(names(tme$gene_sets), tme$feature_table$feature_id)) stop('Malformed HD TME signature object')
sets <- c(setNames(lapply(mp_genes[grouping$mp], function(g) head(unique(g), 100L)), paste0('cancer__', grouping$mp)), tme$gene_sets)
parents <- c(setNames(rep('Carcinoma', 17L), paste0('cancer__', grouping$mp)),
             setNames(as.character(tme$feature_table$annotation_celltype), tme$feature_table$feature_id))
####################
# Align raw matrix axes directly to the exported H5AD identities; row numbers
# alone cannot detect a stale metadata table with a different cell order.
obs_index <- h5readAttributes(h5, 'obs')[['_index']]
var_index <- h5readAttributes(h5, 'var')[['_index']]
cell_ids <- as.character(h5read(h5, paste0('obs/', obs_index)))
genes <- as.character(h5read(h5, paste0('var/', var_index)))
count_attributes <- h5readAttributes(h5, 'layers/counts')
if (!identical(as.character(meta$cell_id), cell_ids) || anyDuplicated(cell_ids) ||
    !identical(as.integer(count_attributes$shape), as.integer(c(nrow(meta), length(genes)))) ||
    count_attributes[['encoding-type']] != 'csr_matrix') stop('Raw count axes disagree with metadata')
if (anyNA(meta$malignant) || !is.logical(meta$malignant) || anyNA(meta$target_type)) stop('Malformed scoring populations')
if (any(lengths(sets) == 0L) || anyNA(parents) || !identical(names(sets), names(parents))) stop('Malformed signature definitions')
####################
if (anyDuplicated(genes)) stop('Duplicate Xenium gene names')
coverage <- data.table(feature_id = names(sets), n_signature = lengths(sets),
                      n_panel_overlap = vapply(sets, function(g) length(intersect(g, genes)), integer(1)))
coverage[, coverage_fraction := n_panel_overlap / n_signature]
coverage[, parent := unname(parents[feature_id])]
fwrite(coverage, file.path(out, 'tables/Auto_signature_coverage.csv'))
usable_sets <- lapply(sets, intersect, y = genes)
if (!identical(meta$row_index, seq_len(nrow(meta)) - 1L)) stop('Metadata is not in H5AD row order')
data_values <- as.numeric(h5read(h5, 'layers/counts/data'))
if (any(!is.finite(data_values)) || any(data_values < 0) || any(data_values != round(data_values))) stop('Xenium counts layer is not raw integer counts')
counts <- sparseMatrix(i = as.integer(h5read(h5, 'layers/counts/indices')) + 1L,
                       p = as.integer(h5read(h5, 'layers/counts/indptr')), x = data_values,
                       dims = as.integer(c(length(genes), nrow(meta))),
                       dimnames = list(genes, as.character(meta$row_index)))
rm(data_values)
validObject(counts)
n_cores <- as.integer(Sys.getenv('SCREF_NCORES', '8'))
if (!is.finite(n_cores) || n_cores < 1L) stop('SCREF_NCORES must be positive')
for (patient in unique(meta$patient)) {
  message(format(Sys.time()), ': UCell ', patient)
  idx <- which(meta$patient == patient)
  sub <- meta[idx]
  score_table <- data.table(row_index = sub$row_index)
  audits <- list()
  for (parent in unique(parents)) {
    feature_ids <- names(parents)[parents == parent]
    local_idx <- which(if (parent == 'Carcinoma') sub$malignant else sub$target_type == parent)
    for (fid in feature_ids) score_table[, (fid) := NA_real_]
    if (length(local_idx) >= 20L) {
      if (any(lengths(usable_sets[feature_ids]) == 0L)) stop('No panel genes for ', patient, ' ', parent)
      scores <- as.matrix(ScoreSignatures_UCell(matrix = counts[, idx[local_idx], drop = FALSE],
                       features = usable_sets[feature_ids], maxRank = 1500L, name = '', ncores = n_cores, force.gc = TRUE))
      if (!all(feature_ids %in% colnames(scores)) || any(!is.finite(scores))) stop('Invalid UCell scores')
      expected_rows <- as.character(sub$row_index[local_idx])
      if (anyDuplicated(rownames(scores)) || !setequal(rownames(scores), expected_rows)) stop('UCell cell identity mismatch')
      scores <- scores[match(expected_rows, rownames(scores)), feature_ids, drop = FALSE]
      for (j in seq_along(feature_ids)) set(score_table, i = local_idx, j = feature_ids[j], value = scores[, j])
      cutoffs <- colMeans(scores) + apply(scores, 2L, sd)
      positives <- colSums(sweep(scores, 2L, cutoffs, '>'))
    } else {
      cutoffs <- setNames(rep(NA_real_, length(feature_ids)), feature_ids)
      positives <- setNames(rep(0L, length(feature_ids)), feature_ids)
    }
    min_positive <- if (parent == 'Carcinoma') 20L else 10L
    a <- coverage[match(feature_ids, feature_id)]
    a[, `:=`(sample = patient, n_parent = length(local_idx), cutoff = unname(cutoffs[feature_ids]),
             n_positive = as.integer(positives[feature_ids]),
             feature_side = fifelse(parent == 'Carcinoma', 'malignant_source', 'tme_target'))]
    a[, eligible := n_parent >= 20L & n_positive >= min_positive]
    a[, skip_reason := fifelse(n_parent < 20L, 'insufficient_parent_observations',
                              fifelse(n_positive < min_positive, 'insufficient_positive_observations', ''))]
    audits[[length(audits) + 1L]] <- a
  }
  fwrite(score_table, file.path(out, 'intermediate', paste0('Auto_', patient, '_ucell_scores.csv.gz')), compress = 'gzip')
  fwrite(rbindlist(audits), file.path(out, 'tables', paste0('Auto_', patient, '_ucell_feature_audit.csv')))
  rm(score_table)
  invisible(gc())
}
writeLines(capture.output(sessionInfo()), file.path(out, 'logs/Auto_ucell_session_info.txt'))
message('Verified raw-count UCell scoring complete')
####################
