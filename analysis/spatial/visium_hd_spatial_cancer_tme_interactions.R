#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: terminal
#   Script: analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R
#   Description: empirical local-neighbourhood validation of scRNA-predicted
#     current malignant MP/state–TME interactions in filtered Visium HD 16 um
#     bins from SUR1231, FFPEA1, and FFPED1 (D1).
#   Methodology:
#     analysis/methodology/spatial/visium_hd_spatial_cancer_tme_interactions_methodology.md
#   Inputs:
#     analysis/shared/scRef_config.R
#     ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#     ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#     ref_outs/visium_hd_outs/malignancy/tables/Auto_<sample>_binned_malignancy.csv.gz
#     ref_outs/visium_hd_outs/state_mapping/tables/Auto_visiumhd_binned_malignant_state_annotations.csv.gz
#     ref_outs/nmf_{fibroblast,endothelial,macrophage,nk,plasma,cd4,cd8}/MP_outs_default.rds
#     ref_outs/non_malignant_mp_correlations/{01,03,05,06}_*/Auto_celltype_correlations_all.csv
#     /rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um/filtered_feature_bc_matrix.h5
#     /rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um/filtered_feature_bc_matrix.h5
#     /rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um/filtered_feature_bc_matrix.h5
#   Outputs under ref_outs/visium_hd_outs/spatial_interactions/:
#     intermediate/Auto_<sample>_spatial_interaction_memberships.rds;
#     intermediate/Auto_visiumhd_spatial_mp_ucell_scores.{rds,csv.gz};
#     intermediate/Auto_visiumhd_tme_mp_gene_sets.rds.
#     tables/Auto_spatial_interactions_{by_sample,by_donor,pooled,scrna_validation}.csv.gz;
#     tables/Auto_spatial_interactions_{recurrence,significant_by_sample,significant_pooled}.csv;
#     tables/Auto_spatial_interaction_{source_feature_audit,target_feature_audit,scrna_candidate_availability,validated_map_summary,parameters,compact_summary}.csv;
#     tables/Auto_paired_section_substitution.csv.gz, Auto_paired_section_effect_difference.csv;
#     tables/Auto_spatial_sample_design.csv, Auto_scatlas_figure_source.csv, Auto_malignant_mp_label_key.csv.
#     figures/Auto_final/{01_scatlas_spatial_evidence,02_supported_tissue_context}.pdf;
#     figures/Auto_supporting/Auto_primary_all_pairs.pdf;
#     figures/Auto_diagnostics/Auto_window_sensitivity.pdf.
#     reports/CURRENT_RESULTS.md;
#     logs/Auto_membership_input_provenance.csv;
#     logs/Auto_visium_hd_spatial_cancer_tme_interactions_{run_summary,session_info}.txt.
#     updates/new_updates/summaries/visium_hd_spatial_cancer_tme_interactions_summary.csv.
#   Additional input: analysis/spatial/spatial_interaction_sample_design.tsv;
#     explicit donor and treatment labels supplied by the user.
#   Cache/replot: SCREF_FORCE_REBUILD=TRUE rebuilds scores and memberships;
#     normal runs reuse matching caches and recompute tests. SCREF_REPLOT_ONLY=TRUE
#     requires current memberships, input fingerprints and per-section tests,
#     and recomputes donor pooling, candidate tests and figures. Flags conflict.
#     SCREF_SPATIAL_PERMUTATIONS defaults to 4999 per tile origin and cannot be lower;
#     SCREF_NCORES controls workers; eight fixed random chunks preserve draws.
#   Run: /opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE job_scripts/spatial/visium_hd_spatial_cancer_tme_interactions.sh
#   Direct run: Rscript analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R
#   Resources/environment: eight CPUs, 160 GB RAM, 24 hours; conda dmtcp.
####################

####################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(scales)
  library(UCell)
  library(digest)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)
source(file.path(wd, "analysis", "shared", "scRef_config.R"))

samples <- c("SUR1231", "FFPEA1", "FFPED1")
sample_display <- c("SUR1231" = "SUR1231 (post)", "FFPEA1" = "A1 (paired post)", "FFPED1" = "D1 (paired pre)")
output_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "spatial_interactions")
output_tiers <- setNames(file.path(output_dir, c("intermediate", "tables", "figures", "logs", "reports")),
                         c("intermediate", "tables", "figures", "logs", "reports"))
for (path in output_tiers) dir.create(path, recursive = TRUE, showWarnings = FALSE)
summary_dir <- file.path(wd, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)


malignancy_dir <- file.path(wd, "ref_outs", "visium_hd_outs", "malignancy", "tables")
state_path <- file.path(
  wd, "ref_outs", "visium_hd_outs", "state_mapping", "tables",
  "Auto_visiumhd_binned_malignant_state_annotations.csv.gz"
)
mp_gene_path <- SCREF_MP_GENES_RDS
mp_grouping_path <- SCREF_MP_GROUPING_CSV

cache_version <- "spatial_ucell_memberships_2026_10_22"
force_rebuild <- tolower(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")) %in% c("true", "t", "1", "yes")
replot_only <- tolower(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")) %in% c("true", "t", "1", "yes")
n_permutations <- as.integer(Sys.getenv("SCREF_SPATIAL_PERMUTATIONS", "4999"))
n_cores <- max(1L, as.integer(Sys.getenv(
  "SCREF_NCORES",
  if (nzchar(Sys.getenv("PBS_JOBID", ""))) Sys.getenv("NCPUS", "8") else "1"
)))
if (!is.finite(n_permutations) || n_permutations < 4999L) stop("SCREF_SPATIAL_PERMUTATIONS must be >= 4999")

neighbour_rings <- c(1L, 2L, 3L)
primary_ring <- 1L
mp_positive_sd_multiplier <- 1
ucell_max_rank <- 1500L
top_genes_per_mp <- 100L
min_source_bins <- 20L
min_whole_target_bins <- 20L
min_parent_bins_for_mp <- 20L
min_mp_positive_bins <- 10L
min_observed_neighbour_slots <- 20L
sample_fdr_cutoff <- 0.05
pooled_fdr_cutoff <- 0.05
recurrence_fdr_cutoff <- 0.10
spatial_block_bins <- 20L
spatial_block_offsets <- rbind(
  c(0L, 0L), c(spatial_block_bins %/% 2L, 0L),
  c(0L, spatial_block_bins %/% 2L),
  c(spatial_block_bins %/% 2L, spatial_block_bins %/% 2L)
)
inference_version <- "donor_aware_origin_robust_support_2026_10_27"
permutation_method <- "max_p_across_shifted_origins"

####################
# Canonical current MP/state metadata
####################
if (!file.exists(mp_grouping_path)) stop("Missing current MP grouping: ", mp_grouping_path)
mp_grouping <- fread(mp_grouping_path)
required_grouping_columns <- c("state", "mp", "description", "plot_label")
if (!all(required_grouping_columns %in% names(mp_grouping))) {
  stop("Current MP grouping lacks: ", paste(setdiff(required_grouping_columns, names(mp_grouping)), collapse = ", "))
}
if (anyDuplicated(mp_grouping$mp)) stop("Current MP grouping contains duplicate MPs")
if (nrow(mp_grouping) != 17L) stop("Current MP grouping must contain exactly 17 MPs; found ", nrow(mp_grouping))

state_crosswalk <- c(
  "Classic proliferation" = "Classic proliferation",
  "Squamous-to-intestinal" = "Squamous-to-intestinal",
  "Basal to intestinal metaplasia" = "Squamous-to-intestinal",
  "Glandular-to-intestinal" = "Glandular-to-intestinal",
  "SMG to intestinal metaplasia" = "Glandular-to-intestinal",
  "Stress-adaptive" = "Stress-adaptive",
  "Stress adaptive" = "Stress-adaptive",
  "Cancer-cell immune mimicry" = "Cancer-cell immune mimicry",
  "Hybrid" = "Hybrid", "Unresolved" = "Unresolved"
)
state_order <- SCREF_PRIMARY_STATE_ORDER
state_colours <- SCREF_STATE_COLOURS[state_order]
cancer_mp_descriptions <- setNames(mp_grouping$description, mp_grouping$mp)
cancer_mp_plot_labels <- setNames(mp_grouping$plot_label, mp_grouping$mp)
cancer_mp_order <- mp_grouping$mp
cancer_excluded_mps <- c("MP2x", "MP11c", "MP18a")
if (any(mp_grouping$mp %in% cancer_excluded_mps)) {
  stop("Excluded MP found in current 17-MP grouping: ", paste(intersect(mp_grouping$mp, cancer_excluded_mps), collapse = ", "))
}
cancer_mp_alias_to_id <- c(
  setNames(names(cancer_mp_descriptions), cancer_mp_descriptions),
  "replication-stress-associated cell cycling" = "MP13+",
  "Stress-reactive columnar epithelium" = "MP6+",
  "Epithelial antiviral interferon response" = "MP11+",
  "Immune-interactive glandular progenitor" = "MP17",
  "Hypoxic inflammatory adaptive plasticity" = "MP12",
  "T/NK-like cancer-cell immune mimicry" = "MP15"
)
####################

celltype_colour_table <- data.frame(
  celltype = c("epithelial", "fibroblast", "endothelial", "macrophage", "mast", "t.cell", "b.cell", "nk.cell", "plasma", "dendritic", "lymph", "erythrocyte", "keratinocyte", "neutrophil", "unresolved", "combined"),
  colour = c("#D73027", "#8C564B", "#1F78B4", "#FF7F00", "#A65628", "#33A02C", "#377EB8", "#984EA3", "#E377C2", "#17BECF", "#6BAED6", "#7F7F7F", "#E6AB02", "#1B9E77", "#BDBDBD", "#555555"),
  stringsAsFactors = FALSE
)
celltype_colours <- setNames(celltype_colour_table$colour, celltype_colour_table$celltype)
celltype_colours <- celltype_colours[names(celltype_colours) != "combined"]
celltype_colours[c("nk", "cd4", "cd8")] <- celltype_colours[c("nk.cell", "t.cell", "t.cell")]

tme_mp_cfg <- data.table(
  compartment = c("fibroblast", "endothelial", "macrophage", "nk", "plasma", "cd4", "cd8"),
  annotation_celltype = c("fibroblast", "endothelial", "macrophage", "nk.cell", "plasma", "t.cell", "t.cell"),
  mp_path = file.path(
    wd, "ref_outs",
    c("nmf_fibroblast", "nmf_endothelial", "nmf_macrophage", "nmf_nk", "nmf_plasma", "nmf_cd4", "nmf_cd8"),
    "MP_outs_default.rds"
  )
)

####################
# scRNA CD4/CD8 predictions are tested against the spatial t-cell compartment
# because the current Visium HD annotation does not resolve those subtypes.
####################
whole_target_spatial_crosswalk <- c(
  "fibroblast" = "fibroblast", "endothelial" = "endothelial",
  "macrophage" = "macrophage", "nk" = "nk.cell", "plasma" = "plasma",
  "cd4" = "t.cell", "cd8" = "t.cell"
)

mode_table <- data.table(
  mode_id = c(
    "01_cancer_mps_vs_tme_mps", "02_cancer_states_vs_tme_mps",
    "03_cancer_mps_vs_whole_celltypes", "04_cancer_states_vs_whole_celltypes"
  ),
  source_kind = c("cancer_mp", "cancer_state", "cancer_mp", "cancer_state"),
  target_kind = c("tme_mp", "tme_mp", "whole_celltype", "whole_celltype"),
  mode_label = c(
    "Malignant MPs vs non-malignant MPs", "Malignant states vs non-malignant MPs",
    "Malignant MPs vs whole cell types", "Malignant states vs whole cell types"
  ),
  scrna_subdir = c(
    "01_cancer_mps_cross_only", "03_cancer_states_cross_only",
    "05_cancer_mps_vs_whole_celltypes", "06_cancer_states_vs_whole_celltypes"
  )
)

run_start <- Sys.time()
run_log <- c(
  paste0("start=", format(run_start, tz = "Europe/London")),
  paste0("cache_version=", cache_version),
  paste0("samples=", paste(samples, collapse = ";")),
  "sample_alias=D1 is canonical manifest sample FFPED1",
  paste0("n_permutations=", n_permutations),
  paste0("n_cores=", n_cores),
  paste0("force_rebuild=", force_rebuild),
  paste0("replot_only=", replot_only)
)

as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

parse_grid_coordinates <- function(barcodes) {
  matched <- regexec("^s_016um_([0-9]+)_([0-9]+)-[0-9]+$", as.character(barcodes))
  pieces <- regmatches(as.character(barcodes), matched)
  valid <- lengths(pieces) == 3L
  row <- col <- rep(NA_integer_, length(barcodes))
  row[valid] <- as.integer(vapply(pieces[valid], `[[`, character(1), 2L))
  col[valid] <- as.integer(vapply(pieces[valid], `[[`, character(1), 3L))
  data.table(grid_row = row, grid_col = col)
}

read_10x_counts <- function(input_dir) {
  path <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(path)) stop("Missing required filtered 10x H5 matrix: ", path)
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) stop("10x H5 lacks the required Gene Expression matrix: ", path)
    counts <- counts[["Gene Expression"]]
  }
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != round(counts@x))) {
    stop("Raw counts must be finite non-negative integers with unique cell/gene identities")
  }
  counts
}

filter_gene_sets_by_silhouette <- function(mp_outs) {
  mp_genes <- mp_outs$metaprograms.genes
  metrics <- mp_outs$metaprograms.metrics
  if (is.null(mp_genes) || is.null(metrics$silhouette) ||
      is.null(names(mp_genes)) || anyDuplicated(names(mp_genes)) ||
      is.null(rownames(metrics)) || anyDuplicated(rownames(metrics)) ||
      !setequal(names(mp_genes), rownames(metrics))) {
    stop("GeneNMF signatures and metric row names must match exactly")
  }
  silhouette <- setNames(metrics$silhouette, rownames(metrics))[names(mp_genes)]
  if (any(!is.finite(silhouette))) stop("Non-finite GeneNMF silhouette")
  mp_genes[silhouette >= 0]
}

load_tme_gene_sets <- function() {
  rows <- list()
  gene_sets <- list()
  for (idx in seq_len(nrow(tme_mp_cfg))) {
    cfg <- tme_mp_cfg[idx]
    if (!file.exists(cfg$mp_path)) stop("Missing GeneNMF object: ", cfg$mp_path)
    message("Loading silhouette-filtered gene sets: ", cfg$compartment)
    mp_outs <- readRDS(cfg$mp_path)
    retained <- filter_gene_sets_by_silhouette(mp_outs)
    rm(mp_outs)
    invisible(gc())
    retained <- lapply(retained, function(genes) head(unique(as.character(genes)), top_genes_per_mp))
    for (mp_name in names(retained)) {
      feature_id <- paste(cfg$compartment, mp_name, sep = "__")
      gene_sets[[feature_id]] <- retained[[mp_name]]
      rows[[length(rows) + 1L]] <- data.table(
        feature_id = feature_id, compartment = cfg$compartment,
        annotation_celltype = cfg$annotation_celltype, mp = mp_name,
        feature_label = paste0(cfg$compartment, " ", mp_name),
        n_signature_genes = length(retained[[mp_name]]), mp_path = cfg$mp_path
      )
    }
  }
  list(gene_sets = gene_sets, feature_table = rbindlist(rows))
}

####################
# The malignant signatures are the exact current 17-MP object, ordered and
# described by the canonical grouping table. Historical coverage failures and
# the explicit MP18a exclusion are rejected if they reappear.
####################
load_cancer_gene_sets <- function() {
  if (!file.exists(mp_gene_path)) stop("Missing current malignant MP genes: ", mp_gene_path)
  mp_genes <- readRDS(mp_gene_path)
  expected <- mp_grouping$mp
  if (!setequal(names(mp_genes), expected)) {
    stop(
      "Current malignant MP gene object/grouping mismatch. Missing: ",
      paste(setdiff(expected, names(mp_genes)), collapse = ", "),
      "; unexpected: ", paste(setdiff(names(mp_genes), expected), collapse = ", ")
    )
  }
  mp_genes <- mp_genes[expected]
  lapply(mp_genes, function(genes) head(unique(as.character(genes)), top_genes_per_mp))
}
####################

####################
# Rank-based MP scoring shared with the single-cell interaction workflow.
####################
score_signature_matrix <- function(counts, gene_sets) {
  usable_sets <- lapply(gene_sets, intersect, y = rownames(counts))
  empty_sets <- names(usable_sets)[lengths(usable_sets) == 0L]
  if (length(empty_sets)) stop("No detected signature genes for: ", paste(empty_sets, collapse = ", "))
  scores <- UCell::ScoreSignatures_UCell(
    matrix = counts,
    features = usable_sets,
    maxRank = ucell_max_rank,
    name = "",
    ncores = n_cores,
    force.gc = TRUE
  )
  scores <- as.matrix(scores[, names(usable_sets), drop = FALSE])
  if (anyDuplicated(rownames(scores)) || !setequal(rownames(scores), colnames(counts)) ||
      any(!is.finite(scores))) stop("Invalid UCell scores or cell identity mismatch")
  scores[colnames(counts), , drop = FALSE]
}

derive_activity_cutoffs <- function(score_matrix) {
  apply(score_matrix, 2L, function(values) {
    spread <- stats::sd(values, na.rm = TRUE)
    if (!is.finite(spread)) stop("Non-finite UCell score spread")
    mean(values, na.rm = TRUE) + mp_positive_sd_multiplier * spread
  })
}
####################

make_sparse_membership <- function(labels, levels) {
  column <- match(as.character(labels), levels)
  keep <- which(!is.na(column))
  Matrix::sparseMatrix(
    i = keep, j = column[keep], x = 1,
    dims = c(length(labels), length(levels)),
    dimnames = list(NULL, levels)
  )
}

make_adjacency <- function(source_meta, target_meta, ring) {
  target_keys <- paste(target_meta$grid_row, target_meta$grid_col, sep = ":")
  if (anyDuplicated(target_keys)) stop("Duplicate target 16 um grid positions")
  offsets <- CJ(dr = seq.int(-ring, ring), dc = seq.int(-ring, ring))
  offsets <- offsets[!(dr == 0L & dc == 0L)]
  edge_i <- edge_j <- vector("list", nrow(offsets))
  for (idx in seq_len(nrow(offsets))) {
    query <- paste(source_meta$grid_row + offsets$dr[idx], source_meta$grid_col + offsets$dc[idx], sep = ":")
    matched <- match(query, target_keys)
    keep <- which(!is.na(matched))
    edge_i[[idx]] <- keep
    edge_j[[idx]] <- matched[keep]
  }
  i <- unlist(edge_i, use.names = FALSE)
  j <- unlist(edge_j, use.names = FALSE)
  Matrix::sparseMatrix(
    i = i, j = j, x = 1,
    dims = c(nrow(source_meta), nrow(target_meta)),
    dimnames = list(source_meta$barcode, target_meta$barcode)
  )
}

safe_bh <- function(p) {
  result <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  result[keep] <- p.adjust(p[keep], method = "BH")
  result
}
####################

####################
# Persistent per-sample feature construction
####################
prepare_sample_membership <- function(sample_name, manifest, mapped_states, tme_reference, cancer_gene_sets) {
  message("Preparing spatial memberships for ", sample_name)
  malignancy_path <- file.path(malignancy_dir, paste0("Auto_", sample_name, "_binned_malignancy.csv.gz"))
  if (!file.exists(malignancy_path)) stop("Missing malignancy table: ", malignancy_path)
  cell_table <- fread(malignancy_path)
  required <- c(
    "barcode", "sample", "Auto_postfilter_celltype", "Auto_postfilter_keep",
    "pxl_col_in_fullres", "pxl_row_in_fullres", "Auto_malignancy"
  )
  missing <- setdiff(required, names(cell_table))
  if (length(missing)) stop(sample_name, " malignancy table lacks: ", paste(missing, collapse = ", "))
  cell_table <- cell_table[sample == sample_name]
  grid <- parse_grid_coordinates(cell_table$barcode)
  cell_table[, `:=`(grid_row = grid$grid_row, grid_col = grid$grid_col)]
  cell_table <- cell_table[as_bool(Auto_postfilter_keep)]
  if (any(!is.finite(cell_table$grid_row)) || any(!is.finite(cell_table$grid_col)) ||
      any(!is.finite(cell_table$pxl_col_in_fullres)) || any(!is.finite(cell_table$pxl_row_in_fullres))) {
    stop("Retained annotation contains invalid spatial coordinates: ", sample_name)
  }
  if (anyDuplicated(cell_table$barcode)) stop("Duplicate barcodes in ", sample_name)

  source_meta <- mapped_states[sample == sample_name]
  if (!nrow(source_meta)) stop("No mapped malignant bins for ", sample_name)
  source_grid <- parse_grid_coordinates(source_meta$barcode)
  source_meta[, `:=`(grid_row = source_grid$grid_row, grid_col = source_grid$grid_col)]
  if (any(!is.finite(source_meta$grid_row)) || any(!is.finite(source_meta$grid_col))) {
    stop("Mapped source contains malformed 16-um coordinates: ", sample_name)
  }
  if (anyDuplicated(source_meta$barcode) || !all(source_meta$barcode %in% cell_table$barcode)) {
    stop("Mapped malignant identities must uniquely match retained annotation")
  }
  source_calls <- cell_table[match(source_meta$barcode, barcode)]
  if (any(source_calls$Auto_postfilter_celltype != "epithelial") ||
      any(!source_calls$Auto_malignancy %in% c("malignant_level_1", "malignant_level_2"))) {
    stop("Mapped source bins must be malignant epithelium")
  }
  if (!nrow(source_meta)) stop("Mapped malignant bins do not overlap retained annotation for ", sample_name)
  setkey(cell_table, barcode)
  source_meta[, `:=`(
    pxl_col_in_fullres = cell_table[barcode, pxl_col_in_fullres],
    pxl_row_in_fullres = cell_table[barcode, pxl_row_in_fullres]
  )]

  raw_states <- as.character(source_meta$Auto_state_B)
  canonical_states <- unname(state_crosswalk[raw_states])
  unknown_states <- unique(raw_states[is.na(canonical_states)])
  if (length(unknown_states)) stop("Unmapped spatial state labels: ", paste(unknown_states, collapse = ", "))
  source_meta[, Auto_state_current := canonical_states]
  state_membership <- make_sparse_membership(canonical_states, state_order)
  rownames(state_membership) <- source_meta$barcode
  state_counts <- Matrix::colSums(state_membership)
  keep_states <- names(state_counts)[state_counts >= min_source_bins]
  state_membership <- state_membership[, keep_states, drop = FALSE]

  target_meta <- cell_table[
    Auto_postfilter_celltype != "epithelial" &
      !is.na(Auto_postfilter_celltype) & Auto_postfilter_celltype != "unresolved"
  ]
  if (length(intersect(source_meta$barcode, target_meta$barcode))) stop("Source and target populations overlap")
  target_meta[, whole_celltype := as.character(Auto_postfilter_celltype)]
  whole_counts <- target_meta[, .N, by = whole_celltype]
  whole_keep <- whole_counts[N >= min_whole_target_bins, whole_celltype]
  whole_membership <- make_sparse_membership(target_meta$whole_celltype, whole_keep)
  whole_features <- data.table(
    target_feature = whole_keep, target_kind = "whole_celltype",
    target_compartment = whole_keep, target_annotation_celltype = whole_keep,
    target_label = whole_keep, permutation_block = "all_non_epithelial",
    n_positive_bins = as.integer(Matrix::colSums(whole_membership))
  )
  colnames(whole_membership) <- paste("whole_celltype", whole_keep, sep = "__")
  whole_features[, target_column := colnames(whole_membership)]

  input_dir <- manifest[sample == sample_name, binned_input][1L]
  if (is.na(input_dir) || !dir.exists(input_dir)) stop("Missing manifest input for ", sample_name)
  counts <- read_10x_counts(input_dir)
  missing_barcodes <- setdiff(c(source_meta$barcode, target_meta$barcode), colnames(counts))
  if (length(missing_barcodes)) stop(sample_name, ": ", length(missing_barcodes), " retained target bins absent from counts")

  ####################
  # UCell malignant-MP activity; state nodes remain assignment based.
  ####################
  cancer_mp_scores <- score_signature_matrix(
    counts[, source_meta$barcode, drop = FALSE], cancer_gene_sets
  )
  cancer_mp_cutoffs <- derive_activity_cutoffs(cancer_mp_scores)
  cancer_mp_positive <- sweep(cancer_mp_scores, 2L, cancer_mp_cutoffs, ">")
  mp_counts <- colSums(cancer_mp_positive)
  keep_cancer_mps <- names(mp_counts)[mp_counts >= min_source_bins]
  cancer_mp_membership <- Matrix::Matrix(
    cancer_mp_positive[, keep_cancer_mps, drop = FALSE], sparse = TRUE
  )
  source_membership <- cbind(cancer_mp_membership, state_membership)
  source_features <- rbindlist(list(
    data.table(
      source_feature = keep_cancer_mps, source_kind = "cancer_mp",
      source_label = unname(cancer_mp_plot_labels[keep_cancer_mps]),
      source_validation_key = keep_cancer_mps,
      source_validation_key_whole = keep_cancer_mps,
      activity_cutoff = as.numeric(cancer_mp_cutoffs[keep_cancer_mps])
    ),
    data.table(
      source_feature = keep_states, source_kind = "cancer_state",
      source_label = keep_states, source_validation_key = keep_states,
      source_validation_key_whole = keep_states, activity_cutoff = NA_real_
    )
  ), fill = TRUE)
  colnames(source_membership) <- paste(source_features$source_kind, source_features$source_feature, sep = "__")
  source_features[, source_column := colnames(source_membership)]
  source_features[, n_positive_bins := as.integer(Matrix::colSums(source_membership))]

  cancer_score_long <- as.data.table(cancer_mp_scores, keep.rownames = "barcode")
  cancer_score_long <- melt(
    cancer_score_long, id.vars = "barcode", variable.name = "feature_id", value.name = "ucell_score"
  )
  cancer_score_long[, `:=`(
    sample = sample_name, feature_side = "malignant_source", compartment = "cancer",
    annotation_celltype = "malignant_epithelial",
    activity_cutoff = as.numeric(cancer_mp_cutoffs[feature_id])
  )]
  cancer_score_long[, positive := ucell_score > activity_cutoff]
  ####################

  score_rows <- list(cancer_score_long)
  score_feature_rows <- list()
  target_mp_membership <- Matrix::Matrix(0, nrow = nrow(target_meta), ncol = 0, sparse = TRUE)
  for (cfg_idx in seq_len(nrow(tme_mp_cfg))) {
    cfg <- tme_mp_cfg[cfg_idx]
    compartment_rows <- tme_reference$feature_table[compartment == cfg$compartment]
    bin_idx <- which(target_meta$whole_celltype == cfg$annotation_celltype)
    if (length(bin_idx) < min_parent_bins_for_mp) {
      score_feature_rows[[length(score_feature_rows) + 1L]] <- compartment_rows[, .(
        target_feature = feature_id, target_kind = "tme_mp", target_compartment = compartment,
        target_annotation_celltype = annotation_celltype, target_label = feature_label,
        permutation_block = annotation_celltype, n_parent_bins = length(bin_idx),
        n_positive_bins = 0L, eligible = FALSE, skip_reason = "insufficient_parent_celltype_bins"
      )]
      next
    }
    barcodes <- target_meta$barcode[bin_idx]
    feature_ids <- compartment_rows$feature_id
    sets <- tme_reference$gene_sets[feature_ids]
    signature_genes <- unique(unlist(sets, use.names = FALSE))
    keep_genes <- intersect(signature_genes, rownames(counts))
    if (!length(keep_genes)) stop("No TME signature genes detected for ", sample_name, " ", cfg$compartment)
    detected_signature_counts <- vapply(sets, function(genes) length(intersect(genes, rownames(counts))), integer(1))
    score_mat <- score_signature_matrix(counts[, barcodes, drop = FALSE], sets)
    activity_cutoffs <- derive_activity_cutoffs(score_mat)
    positive_mat <- Matrix::Matrix(sweep(score_mat, 2L, activity_cutoffs, ">"), sparse = TRUE)
    positive_counts <- Matrix::colSums(positive_mat)
    keep_features <- names(positive_counts)[positive_counts >= min_mp_positive_bins]
    if (length(keep_features)) {
      positive_summary <- summary(positive_mat[, keep_features, drop = FALSE])
      full_block <- Matrix::sparseMatrix(
        i = bin_idx[positive_summary$i],
        j = positive_summary$j,
        x = 1,
        dims = c(nrow(target_meta), length(keep_features)),
        dimnames = list(target_meta$barcode, paste0("tme_mp__", keep_features))
      )
      target_mp_membership <- cbind(target_mp_membership, full_block)
    }
    long_scores <- as.data.table(score_mat, keep.rownames = "barcode")
    long_scores <- melt(long_scores, id.vars = "barcode", variable.name = "feature_id", value.name = "ucell_score")
    long_scores[, `:=`(
      sample = sample_name,
      feature_side = "tme_target",
      compartment = tme_reference$feature_table$compartment[match(feature_id, tme_reference$feature_table$feature_id)],
      annotation_celltype = cfg$annotation_celltype,
      activity_cutoff = as.numeric(activity_cutoffs[feature_id])
    )]
    long_scores[, positive := ucell_score > activity_cutoff]
    score_rows[[length(score_rows) + 1L]] <- long_scores
    score_feature_rows[[length(score_feature_rows) + 1L]] <- compartment_rows[, .(
      target_feature = feature_id, target_kind = "tme_mp", target_compartment = compartment,
      target_annotation_celltype = annotation_celltype, target_label = feature_label,
        permutation_block = annotation_celltype, n_parent_bins = length(bin_idx),
        n_positive_bins = as.integer(positive_counts[feature_id]),
        activity_cutoff = as.numeric(activity_cutoffs[feature_id]),
        n_detected_signature_genes = as.integer(detected_signature_counts[feature_id]),
        eligible = feature_id %in% keep_features,
      skip_reason = ifelse(feature_id %in% keep_features, NA_character_, "insufficient_mp_positive_bins")
    )]
  }
  rm(counts)
  invisible(gc())

  mp_feature_audit <- rbindlist(score_feature_rows, fill = TRUE)
  eligible_mp_features <- mp_feature_audit[eligible %in% TRUE]
  if (ncol(target_mp_membership)) {
    expected_columns <- paste0("tme_mp__", eligible_mp_features$target_feature)
    target_mp_membership <- target_mp_membership[, expected_columns, drop = FALSE]
    eligible_mp_features[, target_column := expected_columns]
  } else {
    eligible_mp_features[, target_column := character()]
  }

  tme_scores <- rbindlist(score_rows, fill = TRUE)
  membership <- list(
    cache_version = cache_version, sample = sample_name,
    source_audit = rbindlist(list(
      data.table(source_kind = "cancer_mp", source_feature = names(mp_counts),
        n_positive_bins = as.integer(mp_counts), eligible = mp_counts >= min_source_bins),
      data.table(source_kind = "cancer_state", source_feature = names(state_counts),
        n_positive_bins = as.integer(state_counts), eligible = state_counts >= min_source_bins)
    ))[, skip_reason := fifelse(eligible, "", "insufficient_positive_bins")],
    whole_target_audit = whole_counts[, .(target_kind = "whole_celltype", target_feature = whole_celltype,
      target_compartment = whole_celltype, target_annotation_celltype = whole_celltype,
      n_positive_bins = N, eligible = N >= min_whole_target_bins,
      skip_reason = fifelse(N >= min_whole_target_bins, "", "insufficient_target_bins"))],
    background = cell_table[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, Auto_postfilter_celltype, Auto_malignancy)],
    source_meta = source_meta[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, Auto_state_B, Auto_state_current)],
    target_meta = target_meta[, .(barcode, pxl_col_in_fullres, pxl_row_in_fullres, grid_row, grid_col, whole_celltype)],
    source_membership = source_membership,
    source_features = source_features,
    target_membership = list(whole_celltype = whole_membership, tme_mp = target_mp_membership),
    target_features = list(whole_celltype = whole_features, tme_mp = eligible_mp_features),
    target_feature_audit = mp_feature_audit
  )
  list(membership = membership, scores = tme_scores)
}

membership_paths <- setNames(
  file.path(output_tiers[["intermediate"]], paste0("Auto_", samples, "_spatial_interaction_memberships.rds")),
  samples
)
mp_score_path <- file.path(output_tiers[["intermediate"]], "Auto_visiumhd_spatial_mp_ucell_scores.rds")
tme_gene_set_path <- file.path(output_tiers[["intermediate"]], "Auto_visiumhd_tme_mp_gene_sets.rds")

if (!file.exists(state_path)) stop("Missing current mapped malignant state/MP table: ", state_path)
manifest <- data.table(
  sample = c("SUR1122", "SUR1231", "FFPEA1", "FFPED1"),
  binned_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"),
  segmented_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs")
)[sample %in% samples]
if (!setequal(manifest$sample, samples)) stop("Manifest lacks one or more requested samples")
mapped_states <- fread(state_path)
if (!setequal(unique(mapped_states$sample), samples)) stop("Mapped state table does not contain exactly the three requested samples")

####################
input_paths <- c(mp_gene_path, mp_grouping_path, state_path, tme_mp_cfg$mp_path,
  file.path(malignancy_dir, paste0("Auto_", samples, "_binned_malignancy.csv.gz")),
  file.path(manifest$binned_input, "filtered_feature_bc_matrix.h5"))
if (any(!file.exists(input_paths))) stop("Missing required input: ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
input_info <- file.info(input_paths)
input_provenance <- data.table(path = input_paths, size_bytes = input_info$size,
  modified_seconds = as.numeric(input_info$mtime))
input_fingerprint <- digest::digest(list(cache_version, input_provenance, ucell_max_rank,
  top_genes_per_mp, mp_positive_sd_multiplier), algo = "sha256")
if (replot_only && force_rebuild) stop("SCREF_REPLOT_ONLY and SCREF_FORCE_REBUILD cannot both be TRUE")
####################
if (replot_only) {
  required_replot <- c(unname(membership_paths), mp_score_path,
                       file.path(output_tiers[["tables"]], "Auto_spatial_interactions_by_sample.csv.gz"),
                       file.path(output_tiers[["tables"]], "Auto_spatial_interactions_pooled.csv.gz"))
  missing <- required_replot[!file.exists(required_replot)]
  if (length(missing)) stop("SCREF_REPLOT_ONLY missing: ", paste(missing, collapse = "; "))
  memberships <- lapply(membership_paths, readRDS)
  if (!all(vapply(memberships, function(x) identical(x$cache_version, cache_version) && identical(x$input_fingerprint, input_fingerprint), logical(1)))) {
    stop("SCREF_REPLOT_ONLY cannot use memberships from an older cache version")
  }
  mp_scores_all <- readRDS(mp_score_path)
  run_log <- c(run_log, "membership_cache=replot_only")
} else {
  need_build <- force_rebuild || !all(file.exists(membership_paths)) || !file.exists(mp_score_path) || !file.exists(tme_gene_set_path)
  if (!need_build) {
    memberships <- lapply(membership_paths, readRDS)
    cache_valid <- all(vapply(memberships, function(x) identical(x$cache_version, cache_version) && identical(x$input_fingerprint, input_fingerprint), logical(1)))
    if (!cache_valid) stop("Stale memberships: explicitly request SCREF_FORCE_REBUILD=TRUE")
  }
  if (need_build) {
    tme_reference <- load_tme_gene_sets()
    cancer_gene_sets <- load_cancer_gene_sets()
    saveRDS(tme_reference, tme_gene_set_path)
    prepared <- lapply(
      samples, prepare_sample_membership,
      manifest = manifest, mapped_states = mapped_states,
      tme_reference = tme_reference, cancer_gene_sets = cancer_gene_sets
    )
    memberships <- lapply(prepared, `[[`, "membership")
    names(memberships) <- samples
    for (sample_name in samples) {
      memberships[[sample_name]]$input_fingerprint <- input_fingerprint
      saveRDS(memberships[[sample_name]], membership_paths[[sample_name]])
    }
    fwrite(input_provenance, file.path(output_tiers[["logs"]], "Auto_membership_input_provenance.csv"))
    mp_scores_all <- rbindlist(lapply(prepared, `[[`, "scores"), fill = TRUE)
    saveRDS(mp_scores_all, mp_score_path)
    fwrite(mp_scores_all, file.path(output_tiers[["intermediate"]], "Auto_visiumhd_spatial_mp_ucell_scores.csv.gz"))
    run_log <- c(run_log, "membership_cache=rebuilt")
  } else {
    names(memberships) <- samples
    mp_scores_all <- readRDS(mp_score_path)
    run_log <- c(run_log, "membership_cache=reused")
  }
}

source_audit <- rbindlist(lapply(memberships, function(x) copy(x$source_audit)[, sample := x$sample]), fill = TRUE)
target_audit <- rbindlist(lapply(memberships, function(x) {
  rbindlist(list(
    copy(x$whole_target_audit),
    copy(x$target_feature_audit)
  ), fill = TRUE)[, sample := x$sample]
}), fill = TRUE)
fwrite(source_audit, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_source_feature_audit.csv"))
fwrite(target_audit, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_target_feature_audit.csv"))
####################

####################
# Canonical plot ordering follows the current grouping table for malignant MPs,
# the shared configuration for states, and GeneNMF object order within each TME
# compartment. These ranks are reused by every interaction comparison.
####################
tme_target_order_key <- unique(target_audit[
  target_kind == "tme_mp", paste(target_compartment, target_feature, sep = "||")
])
whole_target_order <- celltype_colour_table$celltype

source_plot_rank <- function(source_kind, source_feature) {
  ifelse(
    source_kind == "cancer_mp",
    match(source_feature, cancer_mp_order),
    length(cancer_mp_order) + match(source_feature, state_order)
  )
}

target_plot_rank <- function(target_kind, target_compartment, target_feature) {
  ifelse(
    target_kind == "tme_mp",
    match(paste(target_compartment, target_feature, sep = "||"), tme_target_order_key),
    match(target_feature, whole_target_order)
  )
}
####################

####################
# Spatially local, degree-aware neighbourhood permutation tests
####################
make_spatial_permutation_groups <- function(target_meta, target_kind) {
  group_sets <- lapply(seq_len(nrow(spatial_block_offsets)), function(offset_index) {
    row_offset <- spatial_block_offsets[offset_index, 1L]
    col_offset <- spatial_block_offsets[offset_index, 2L]
    tile_row <- floor((target_meta$grid_row + row_offset) / spatial_block_bins)
    tile_col <- floor((target_meta$grid_col + col_offset) / spatial_block_bins)
    tile <- paste(tile_row, tile_col, sep = ":")
    if (target_kind == "tme_mp") tile <- paste(target_meta$whole_celltype, tile, sep = "__")
    split(seq_len(nrow(target_meta)), tile)
  })
  movable_fraction <- vapply(group_sets, function(groups) {
    mean(rep(lengths(groups), lengths(groups)) > 1L)
  }, numeric(1))
  if (any(movable_fraction < 0.5)) {
    stop("Fewer than half of target bins are movable in a local permutation tiling")
  }
  group_sets
}

permutation_summary_chunk <- function(chunk_n, source_weights, target_membership,
                                      permutation_groups, seed, start_iteration) {
  set.seed(seed)
  dims <- c(nrow(source_weights), ncol(target_membership), length(permutation_groups))
  sum_null <- sumsq_null <- array(0, dim = dims)
  greater_equal <- less_equal <- array(0L, dim = dims)
  origin_n <- integer(length(permutation_groups))
  observed <- as.matrix(source_weights %*% target_membership)
  for (iteration in seq_len(chunk_n)) {
    origin <- 1L + ((start_iteration + iteration - 1L) %% length(permutation_groups))
    groups <- permutation_groups[[origin]]
    permuted_rows <- seq_len(nrow(target_membership))
    for (row_idx in groups) {
      if (length(row_idx) > 1L) permuted_rows[row_idx] <- sample(row_idx, length(row_idx), replace = FALSE)
    }
    permuted_counts <- as.matrix(source_weights %*% target_membership[permuted_rows, , drop = FALSE])
    sum_null[, , origin] <- sum_null[, , origin] + permuted_counts
    sumsq_null[, , origin] <- sumsq_null[, , origin] + permuted_counts^2
    greater_equal[, , origin] <- greater_equal[, , origin] + (permuted_counts >= observed)
    less_equal[, , origin] <- less_equal[, , origin] + (permuted_counts <= observed)
    origin_n[origin] <- origin_n[origin] + 1L
  }
  list(n = origin_n, sum = sum_null, sumsq = sumsq_null, ge = greater_equal, le = less_equal)
}

run_permutation_test <- function(membership, ring, target_kind) {
  adjacency <- make_adjacency(membership$source_meta, membership$target_meta, ring)
  source_membership <- membership$source_membership
  target_membership <- membership$target_membership[[target_kind]]
  source_features <- membership$source_features
  target_features <- membership$target_features[[target_kind]]
  if (!ncol(source_membership) || !ncol(target_membership) || !length(adjacency@x)) return(data.table())
  permutation_groups <- make_spatial_permutation_groups(membership$target_meta, target_kind)
  observed <- as.matrix(Matrix::crossprod(source_membership, adjacency %*% target_membership))
  source_slots <- as.numeric(Matrix::crossprod(source_membership, Matrix::rowSums(adjacency)))
  source_counts <- as.numeric(Matrix::colSums(source_membership))
  parent_slots <- sapply(target_features$target_annotation_celltype, function(parent) {
    if (target_kind == "whole_celltype") return(source_slots)
    as.numeric(Matrix::crossprod(source_membership,
      adjacency %*% as.numeric(membership$target_meta$whole_celltype == parent)))
  })
  parent_slots <- matrix(parent_slots, nrow = ncol(source_membership))

  source_weights <- Matrix::crossprod(source_membership, adjacency)
  total_permutations <- n_permutations * length(permutation_groups)
  chunks <- rep(total_permutations %/% 8L, 8L)
  if (total_permutations %% 8L) chunks[seq_len(total_permutations %% 8L)] <- chunks[seq_len(total_permutations %% 8L)] + 1L
  chunk_starts <- c(0L, head(cumsum(chunks), -1L))
  chunks <- chunks[chunks > 0]
  seeds <- 100000L + match(membership$sample, samples) * 10000L + ring * 100L +
    ifelse(target_kind == "whole_celltype", 1L, 2L) + seq_along(chunks)
  chunk_results <- parallel::mclapply(
    seq_along(chunks),
    function(idx) permutation_summary_chunk(
      chunks[idx], source_weights, target_membership,
      permutation_groups, seeds[idx], chunk_starts[idx]
    ),
    mc.cores = min(n_cores, length(chunks)), mc.preschedule = TRUE
  )
  origin_n <- Reduce(`+`, lapply(chunk_results, `[[`, "n"))
  if (any(origin_n != n_permutations)) stop("Incorrect permutation count within a tile origin")
  null_n <- sum(origin_n)
  null_sum <- Reduce(`+`, lapply(chunk_results, `[[`, "sum"))
  null_sumsq <- Reduce(`+`, lapply(chunk_results, `[[`, "sumsq"))
  ge <- Reduce(`+`, lapply(chunk_results, `[[`, "ge"))
  le <- Reduce(`+`, lapply(chunk_results, `[[`, "le"))
  origin_mean <- sweep(null_sum, 3L, origin_n, "/")
  pooled_null_sum <- apply(null_sum, c(1L, 2L), sum)
  null_mean <- pooled_null_sum / null_n
  null_variance <- pmax((apply(null_sumsq, c(1L, 2L), sum) - pooled_null_sum^2 / null_n) / (null_n - 1L), 0)
  null_sd <- sqrt(null_variance)
  z_score <- (observed - null_mean) / null_sd
  z_score[null_sd == 0 & observed == null_mean] <- 0
  z_score[null_sd == 0 & observed > null_mean] <- Inf
  z_score[null_sd == 0 & observed < null_mean] <- -Inf
  p_z_matrix <- matrix(
    pmax(2 * pnorm(-abs(z_score)), .Machine$double.xmin),
    nrow = nrow(observed), ncol = ncol(observed)
  )
  # Each origin defines its own valid within-stratum permutation group. The
  # largest p value requires support under every origin, without treating the
  # dependent tests as independent or mixing incompatible permutation groups.
  origin_p <- array(pmin(1, 2 * pmin(sweep(ge + 1, 3L, origin_n + 1, "/"),
    sweep(le + 1, 3L, origin_n + 1, "/"))), dim = dim(ge))
  origin_signs <- sign(array(observed, dim = dim(origin_mean)) - origin_mean)
  origin_concordant <- apply(origin_signs, c(1L, 2L), function(x) all(x > 0) || all(x < 0))
  p_empirical_matrix <- apply(origin_p, c(1L, 2L), max)
  p_empirical_matrix[!origin_concordant] <- 1

  grid <- CJ(source_index = seq_len(nrow(observed)), target_index = seq_len(ncol(observed)))
  matrix_lookup <- cbind(grid$source_index, grid$target_index)
  grid[, `:=`(
    sample = membership$sample, ring = ring, target_kind = target_kind,
    observed_contacts = as.numeric(observed[matrix_lookup]),
    null_mean_contacts = as.numeric(null_mean[matrix_lookup]),
    null_sd_contacts = as.numeric(null_sd[matrix_lookup]),
    z_score = as.numeric(z_score[matrix_lookup]),
    p_z = as.numeric(p_z_matrix[matrix_lookup]),
    p_empirical = as.numeric(p_empirical_matrix[matrix_lookup]),
    source_neighbour_slots = as.numeric(parent_slots[matrix_lookup]),
    all_tme_neighbour_slots = source_slots[source_index],
    n_source_positive_bins = source_counts[source_index],
    n_target_positive_bins = as.numeric(Matrix::colSums(target_membership))[target_index],
    n_permutations = n_permutations,
    n_permutations_total = null_n,
    permutation_method = permutation_method,
    origin_direction_concordant = as.logical(origin_concordant[matrix_lookup])
  )]
  for (origin in seq_along(origin_n)) {
    grid[, (paste0("null_mean_origin_", origin)) := as.numeric(origin_mean[, , origin][matrix_lookup])]
    grid[, (paste0("tail_ge_origin_", origin)) := as.integer(ge[, , origin][matrix_lookup])]
    grid[, (paste0("tail_le_origin_", origin)) := as.integer(le[, , origin][matrix_lookup])]
    grid[, (paste0("p_origin_", origin)) := as.numeric(origin_p[, , origin][matrix_lookup])]
  }
  grid[, log2_enrichment := log2((observed_contacts + 0.5) / (null_mean_contacts + 0.5))]
  grid <- cbind(
    grid,
    source_features[grid$source_index, .(
      source_feature, source_kind, source_label,
      source_validation_key, source_validation_key_whole
    )],
    target_features[grid$target_index, .(
      target_feature, target_compartment, target_annotation_celltype,
      target_label, permutation_block
    )]
  )
  grid[, mode_id := mode_table$mode_id[match(paste(source_kind, target_kind), paste(mode_table$source_kind, mode_table$target_kind))]]
  grid[, mode_label := mode_table$mode_label[match(mode_id, mode_table$mode_id)]]
  grid[, eligible_pair :=
    n_source_positive_bins >= min_source_bins &
      n_target_positive_bins >= ifelse(target_kind == "whole_celltype", min_whole_target_bins, min_mp_positive_bins) &
      source_neighbour_slots >= min_observed_neighbour_slots]
  grid[eligible_pair == FALSE, `:=`(p_z = NA_real_, p_empirical = NA_real_)]
  grid[, fdr_sample := safe_bh(p_empirical), by = .(sample, mode_id, ring)]
  grid[, direction := fifelse(log2_enrichment > 0, "enriched", fifelse(log2_enrichment < 0, "depleted", "neutral"))]
  grid[]
}

sample_result_path <- file.path(output_tiers[["tables"]], "Auto_spatial_interactions_by_sample.csv.gz")
if (replot_only) {
  sample_results <- fread(sample_result_path)
  if (!all(c("all_tme_neighbour_slots", "permutation_method", "p_origin_4") %in% names(sample_results)) ||
      any(sample_results$permutation_method != permutation_method) ||
      any(sample_results$n_permutations != n_permutations)) stop("Replot requires current per-origin permutation tests")
} else {
  test_rows <- list()
  for (sample_name in samples) {
    membership <- memberships[[sample_name]]
    for (ring in neighbour_rings) {
      for (target_kind in c("whole_celltype", "tme_mp")) {
        message("Permutation test: ", sample_name, " ring ", ring, " ", target_kind)
        test_rows[[length(test_rows) + 1L]] <- run_permutation_test(membership, ring, target_kind)
      }
    }
  }
  sample_results <- rbindlist(test_rows, fill = TRUE)
  setcolorder(sample_results, c(
    "sample", "mode_id", "mode_label", "ring", "source_kind", "source_feature", "source_label",
    "target_kind", "target_compartment", "target_feature", "target_label"
  ))
}
sample_results[is.finite(p_z), p_z := pmax(p_z, .Machine$double.xmin)]
sample_results[, fdr_sample := safe_bh(p_empirical), by = .(sample, mode_id, ring)]
####################

####################
# Cross-sample meta-analysis, recurrence, and scRNA validation
####################
####################
# The biological unit is the donor. A1/D1 are paired pre/post sections.
# Bonferroni merging remains valid without estimating their dependence from one
# pair. Opposite section directions provide no signed donor-level support.
####################
sample_design_path <- file.path(wd, 'analysis/spatial/spatial_interaction_sample_design.tsv')
sample_design <- fread(sample_design_path)[platform == 'visium_hd']
if (!setequal(sample_design$orig.ident, samples) || anyDuplicated(sample_design$orig.ident) ||
    anyNA(sample_design$donor_id)) stop('Invalid Visium HD donor design')
sample_results[, orig.ident := sample]
sample_results[, donor_id := sample_design$donor_id[match(orig.ident, sample_design$orig.ident)]]
sample_results[, treatment_timepoint := sample_design$treatment_timepoint[match(orig.ident, sample_design$orig.ident)]]
fwrite(sample_results, sample_result_path)
fwrite(sample_design, file.path(output_tiers[['tables']], 'Auto_spatial_sample_design.csv'))
interaction_keys <- c('mode_id', 'mode_label', 'ring', 'source_kind', 'source_feature', 'source_label',
  'source_validation_key', 'source_validation_key_whole', 'target_kind', 'target_compartment',
  'target_annotation_celltype', 'target_feature', 'target_label')
if (anyDuplicated(sample_results[, c('sample', interaction_keys), with = FALSE])) {
  stop('Duplicate section-interaction key')
}
aggregate_donors <- function(section_results) {
  section_results[eligible_pair %in% TRUE & is.finite(p_empirical), {
    signs <- sign(log2_enrichment)
    concordant <- all(signs > 0) || all(signs < 0)
    .(n_sections = .N, sections = paste(sample, collapse = ';'),
      donor_log2_enrichment = mean(log2_enrichment),
      donor_p = if (concordant) min(1, .N * min(p_empirical)) else 1,
      section_direction_concordant = concordant,
      observed_contacts = sum(observed_contacts), null_mean_contacts = sum(null_mean_contacts),
      n_permutations = min(n_permutations))
  }, by = c(interaction_keys, 'donor_id')]
}
pool_donors <- function(donor_table) {
  donor_table[, donor_fdr := safe_bh(donor_p), by = .(donor_id, mode_id, ring)]
  result <- donor_table[, {
    donor_z <- sign(donor_log2_enrichment) * qnorm(donor_p / 2, lower.tail = FALSE)
    combined_z <- sum(donor_z) / sqrt(.N)
    .(n_donors_tested = .N, n_sections_tested = sum(n_sections),
      donors_tested = paste(donor_id, collapse = ';'), sections_tested = paste(sections, collapse = ';'),
      pooled_z = combined_z, pooled_p = max(2 * pnorm(-abs(combined_z)), .Machine$double.xmin),
      pooled_log2_enrichment = mean(donor_log2_enrichment),
      contact_weighted_log2_enrichment = log2((sum(observed_contacts) + 0.5) / (sum(null_mean_contacts) + 0.5)),
      n_donors_enriched = sum(donor_log2_enrichment > 0 & section_direction_concordant),
      n_donors_depleted = sum(donor_log2_enrichment < 0 & section_direction_concordant),
      n_donors_fdr_0_10 = sum(donor_fdr < recurrence_fdr_cutoff),
      n_discordant_paired_donors = sum(!section_direction_concordant),
      donor_effects = paste(paste0(donor_id, '=', sprintf('%.4f', donor_log2_enrichment)), collapse = ';'))
  }, by = interaction_keys]
  result[n_donors_tested < 2L, `:=`(pooled_p = NA_real_, pooled_z = NA_real_)]
  result[, pooled_fdr := safe_bh(pooled_p), by = .(mode_id, ring)]
  result[, pooled_direction := fifelse(pooled_log2_enrichment > 0, 'enriched',
    fifelse(pooled_log2_enrichment < 0, 'depleted', 'neutral'))]
  result[, direction_concordant := n_discordant_paired_donors == 0L &
    ((pooled_log2_enrichment > 0 & n_donors_enriched == n_donors_tested) |
     (pooled_log2_enrichment < 0 & n_donors_depleted == n_donors_tested))]
  result
}
donor_results <- aggregate_donors(sample_results)
pooled_results <- pool_donors(donor_results)
fwrite(donor_results, file.path(output_tiers[['tables']], 'Auto_spatial_interactions_by_donor.csv.gz'))

# Each substitution retains SUR1231 (post) and one section of the paired donor.
# These contrasts assess dependence on the paired section, not a treatment effect.
substitution_results <- rbindlist(lapply(c('FFPED1', 'FFPEA1'), function(retained_section) {
  result <- pool_donors(aggregate_donors(sample_results[sample %in% c('SUR1231', retained_section)]))
  result[, paired_section_retained := retained_section]
  result
}))
fwrite(substitution_results, file.path(output_tiers[['tables']], 'Auto_paired_section_substitution.csv.gz'))
paired_effects <- dcast(sample_results[sample %in% c('FFPED1', 'FFPEA1') & eligible_pair %in% TRUE],
  mode_id + ring + source_feature + target_feature ~ sample, value.var = 'log2_enrichment')
paired_effects[, post_minus_pre := FFPEA1 - FFPED1]
fwrite(paired_effects, file.path(output_tiers[['tables']], 'Auto_paired_section_effect_difference.csv'))
recurrence_results <- pooled_results[, .(mode_id, mode_label, ring, source_feature, source_label,
  target_compartment, target_feature, target_label, n_donors_tested, n_sections_tested,
  n_donors_enriched, n_donors_depleted, direction_concordant, pooled_log2_enrichment,
  pooled_p, pooled_fdr, n_donors_fdr_0_10,
  recurrent_same_direction = direction_concordant & n_donors_fdr_0_10 >= 2L,
  pooled_significant = pooled_fdr < pooled_fdr_cutoff & direction_concordant)]
####################

scrna_inputs <- file.path(
  wd, "ref_outs", "non_malignant_mp_correlations", mode_table$scrna_subdir,
  "Auto_celltype_correlations_all.csv"
)
names(scrna_inputs) <- mode_table$mode_id
load_scrna_validation <- function(mode_id) {
  path <- scrna_inputs[[mode_id]]
  if (!file.exists(path)) stop("Missing required scRNA interaction table: ", path)
  input <- fread(path)
  required <- c(
    "compartment1", "compartment2", "mp1_name", "mp1_display", "celltype1_display",
    "mp2_name", "mp2_display", "celltype2_display", "pearson_r",
    "spearman_r", "spearman_p", "spearman_significance", "spearman_sig",
    "shared_sample_n", "edge_id"
  )
  if (!all(required %in% names(input))) {
    stop("scRNA interaction table lacks: ", paste(setdiff(required, names(input)), collapse = ", "))
  }
  input <- input[compartment1 == "cancer" | compartment2 == "cancer"]
  if (!nrow(input)) stop("No cancer-TME rows in required scRNA interaction table: ", path)
  cancer_first <- input$compartment1 == "cancer"
  source_raw <- ifelse(cancer_first, input$mp1_name, input$mp2_name)
  source_kind <- mode_table$source_kind[match(mode_id, mode_table$mode_id)]
  if (source_kind == "cancer_state") {
    canonical <- unname(state_crosswalk[source_raw])
    source_raw[!is.na(canonical)] <- canonical[!is.na(canonical)]
  } else if (mode_id == "03_cancer_mps_vs_whole_celltypes") {
    aliases <- unname(cancer_mp_alias_to_id[source_raw])
    source_raw[!is.na(aliases)] <- aliases[!is.na(aliases)]
  }
  selected_threshold <- if (mode_id == "01_cancer_mps_vs_tme_mps") {
    4
  } else if (mode_id == "02_cancer_states_vs_tme_mps") {
    1.5
  } else {
    -log10(0.05)
  }
  result <- data.table(
    mode_id = mode_id,
    scrna_source_feature = source_raw,
    scrna_source_display = ifelse(cancer_first, input$mp1_display, input$mp2_display),
    scrna_target_compartment = ifelse(cancer_first, input$celltype2_display, input$celltype1_display),
    scrna_target_feature = ifelse(cancer_first, input$mp2_name, input$mp1_name),
    scrna_target_display = ifelse(cancer_first, input$mp2_display, input$mp1_display),
    scrna_spearman_r = input$spearman_r,
    scrna_spearman_p = input$spearman_p,
    scrna_spearman_sig = as_bool(input$spearman_sig),
    scrna_pearson_r = input$pearson_r,
    scrna_selected = as_bool(input$spearman_sig) & is.finite(input$pearson_r) &
      input$spearman_significance >= selected_threshold,
    scrna_selection_threshold = selected_threshold,
    scrna_shared_samples = input$shared_sample_n,
    scrna_edge_id = input$edge_id
  )
  result[, target_compartment := scrna_target_compartment]
  if (mode_id %in% c("03_cancer_mps_vs_whole_celltypes", "04_cancer_states_vs_whole_celltypes")) {
    mapped_target <- unname(whole_target_spatial_crosswalk[result$scrna_target_compartment])
    if (anyNA(mapped_target)) {
      stop(
        "No spatial whole-celltype crosswalk for scRNA compartments: ",
        paste(unique(result$scrna_target_compartment[is.na(mapped_target)]), collapse = ", ")
      )
    }
    result[, target_compartment := mapped_target]
    result[, scrna_target_feature := target_compartment]
  }
  result
}
scrna_validation <- rbindlist(lapply(mode_table$mode_id, load_scrna_validation), fill = TRUE)
if (anyDuplicated(scrna_validation[, .(mode_id, scrna_edge_id)])) stop("Duplicate scATLAS candidate identity")
scrna_validation <- scrna_validation[rep(seq_len(.N), each = length(neighbour_rings))]
scrna_validation[, ring := rep(neighbour_rings, times = .N / length(neighbour_rings))]

validation_base <- copy(pooled_results)
validation_base[, spatial_source_key := ifelse(
  mode_id == "03_cancer_mps_vs_whole_celltypes",
  source_validation_key_whole, source_validation_key
)]
validation_base[, spatial_target_key := ifelse(
  target_kind == "whole_celltype", target_compartment,
  sub("^[^_]+__", "", target_feature)
)]
validation <- merge(
  validation_base,
  scrna_validation,
  by.x = c("mode_id", "ring", "spatial_source_key", "target_compartment", "spatial_target_key"),
  by.y = c("mode_id", "ring", "scrna_source_feature", "target_compartment", "scrna_target_feature"),
  all = TRUE
)
validation[, `:=`(
  scrna_available = is.finite(scrna_spearman_r),
  spatial_available = is.finite(pooled_p) & n_donors_tested >= 2L,
  spatial_scrna_direction_concordant = is.finite(scrna_spearman_r) &
    is.finite(pooled_log2_enrichment) & sign(pooled_log2_enrichment) == sign(scrna_spearman_r),
  same_direction_donor_n = fifelse(
    scrna_spearman_r > 0, n_donors_enriched,
    fifelse(scrna_spearman_r < 0, n_donors_depleted, NA_integer_)
  )
)]
validation[, spatial_unavailable_reason := fifelse(
  spatial_available %in% TRUE, NA_character_,
  fifelse(
    is.finite(pooled_log2_enrichment) & n_donors_tested < 2L,
    "fewer_than_two_eligible_donors", "no_eligible_spatial_test"
  )
)]
####################
validation[scrna_selected %in% TRUE & is.na(ring), ring := primary_ring]
validation[, source_kind := mode_table$source_kind[match(mode_id, mode_table$mode_id)]]
validation[, target_kind := mode_table$target_kind[match(mode_id, mode_table$mode_id)]]
validation[, spatial_test_id := paste(spatial_source_key, target_compartment, spatial_target_key, sep = '||')]
validation[, spatial_validation_fdr := {
  adjusted <- rep(NA_real_, .N)
  selected <- which(scrna_selected %in% TRUE)
  family <- unique(spatial_test_id[selected])
  tested <- selected[is.finite(pooled_p[selected])]
  unique_tested <- tested[!duplicated(spatial_test_id[tested])]
  if (length(unique_tested)) {
    q <- p.adjust(pooled_p[unique_tested], method = 'BH', n = length(family))
    adjusted[tested] <- q[match(spatial_test_id[tested], spatial_test_id[unique_tested])]
  }
  adjusted
}, by = .(mode_id, ring)]
####################
validation[, validation_class := fifelse(
  !scrna_selected %in% TRUE, "Not scRNA-selected",
  fifelse(
    !spatial_available, "Spatially unavailable",
    fifelse(
      spatial_scrna_direction_concordant & same_direction_donor_n >= 2L & n_discordant_paired_donors == 0L & spatial_validation_fdr < pooled_fdr_cutoff,
      "Supported",
      fifelse(spatial_scrna_direction_concordant, "Concordant; not supported", "Direction-discordant")
    )
  )
)]
validation[scrna_selected %in% TRUE & is.na(ring), ring := primary_ring]
validation[, `:=`(
  validation_source_label = fcoalesce(source_label, scrna_source_display, spatial_source_key),
  validation_target_label = fcoalesce(target_label, scrna_target_display, spatial_target_key)
)]

fwrite(pooled_results, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_pooled.csv.gz"))
fwrite(recurrence_results, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_recurrence.csv"))
fwrite(validation, file.path(output_tiers[["tables"]], "Auto_spatial_interactions_scrna_validation.csv.gz"))
fwrite(
  validation[scrna_selected %in% TRUE, .(
    mode_id, ring, spatial_source_key, scrna_target_compartment,
    target_compartment, spatial_target_key,
    scrna_edge_id, scrna_spearman_r, scrna_spearman_p,
    spatial_available, n_donors_tested, donors_tested,
    pooled_log2_enrichment, pooled_p, spatial_validation_fdr,
    spatial_scrna_direction_concordant, same_direction_donor_n,
    spatial_unavailable_reason, validation_class
  )],
  file.path(output_tiers[["tables"]], "Auto_spatial_interaction_scrna_candidate_availability.csv")
)
fwrite(
  sample_results[eligible_pair %in% TRUE & fdr_sample < sample_fdr_cutoff],
  file.path(output_tiers[["tables"]], "Auto_spatial_interactions_significant_by_sample.csv")
)
fwrite(
  pooled_results[pooled_fdr < pooled_fdr_cutoff & direction_concordant %in% TRUE],
  file.path(output_tiers[["tables"]], "Auto_spatial_interactions_significant_pooled.csv")
)
####################

####################
# Export at journal width. Font sizes are specified at the final printed size.
# Every mode keeps its complete comparison on one page.
####################
figure_dirs <- setNames(file.path(output_tiers[['figures']],
  c('Auto_final', 'Auto_supporting', 'Auto_diagnostics')),
  c('final', 'supporting', 'diagnostics'))
for (path in figure_dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
figure_width <- 183 / 25.4
publication_theme <- theme_classic(base_size = 7, base_family = 'sans') +
  theme(axis.text = element_text(size = 7, colour = 'black'),
        axis.title = element_text(size = 7), plot.title = element_text(size = 7, face = 'bold'),
        plot.subtitle = element_text(size = 7), legend.title = element_text(size = 7),
        legend.text = element_text(size = 7), strip.text = element_text(size = 7, face = 'bold'),
        plot.caption = element_text(size = 7, hjust = 0),
        plot.margin = margin(5, 5, 5, 5), legend.position = 'bottom')
state_short <- setNames(c('Classic proliferation', 'Squamous-to-intestinal',
  'Glandular-to-intestinal', 'Stress-adaptive', 'Immune mimicry'), state_order)
source_short <- function(x) {
  result <- as.character(x)
  matched <- result %in% names(state_short)
  result[matched] <- unname(state_short[result[matched]])
  result
}
target_short <- function(x) gsub('__', ' ', x, fixed = TRUE)
empty_plot <- function(title, detail) {
  ggplot() + annotate('text', x = 0, y = 0, label = detail, size = 7 / ggplot2::.pt) +
    theme_void(base_size = 7) + labs(title = title) + publication_theme
}
colour_effect <- function(limit = 1) {
  scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B',
    midpoint = 0, limits = c(-limit, limit), oob = scales::squish,
    na.value = '#D9D9D9', name = 'log2 enrichment',
    breaks = c(-limit, 0, limit), labels = c(paste0('<= ', -limit), '0', paste0('>= ', limit)))
}
mode_titles <- setNames(c('Malignant MPs / TME programmes', 'Malignant states / TME programmes',
  'Malignant MPs / TME cell types', 'Malignant states / TME cell types'), mode_table$mode_id)

# Each candidate retains its original scATLAS target identity. CD4/CD8 whole-type
# candidates explicitly share the spatial T-cell test and never imply subtype validation.
primary_candidates <- copy(validation[ring == primary_ring & scrna_selected %in% TRUE])
primary_candidates[, candidate_label := paste(source_short(spatial_source_key),
  ifelse(target_kind == 'tme_mp', paste(scrna_target_compartment, spatial_target_key),
    ifelse(scrna_target_compartment %in% c('cd4', 'cd8'),
      paste0(toupper(scrna_target_compartment), ' / spatial T cells'), scrna_target_compartment)), sep = ' : ')]
primary_candidates[, source_order := source_plot_rank(source_kind, spatial_source_key)]
setorder(primary_candidates, mode_id, source_order, scrna_target_compartment, spatial_target_key)
primary_candidates[, candidate_number := seq_len(.N), by = mode_id]
primary_candidates[, evidence_id := paste(mode_id, candidate_number, sep = '__')]
fwrite(primary_candidates, file.path(output_tiers[['tables']], 'Auto_scatlas_figure_source.csv'))

pdf(file.path(figure_dirs[['final']], '01_scatlas_spatial_evidence.pdf'),
    width = figure_width, height = 190 / 25.4, onefile = TRUE, useDingbats = FALSE)
for (mode in mode_table$mode_id) {
  page <- copy(primary_candidates[mode_id == mode])
  if (!nrow(page)) {
    print(empty_plot(mode_titles[[mode]], 'No prespecified scATLAS candidates'))
    next
  }
  page[, row_label := factor(evidence_id, levels = rev(evidence_id))]
  row_labels <- setNames(page$candidate_label, page$evidence_id)
  section_rows <- rbindlist(lapply(seq_len(nrow(page)), function(i) {
    candidate <- page[i]
    found <- sample_results[mode_id == mode & ring == primary_ring &
      source_feature == candidate$spatial_source_key & target_compartment == candidate$target_compartment &
      if (candidate$target_kind == 'tme_mp') sub('^[^_]+__', '', target_feature) == candidate$spatial_target_key
      else target_feature == candidate$spatial_target_key]
    found[, evidence_id := candidate$evidence_id]
    found
  }), fill = TRUE)
  section_rows[, row_label := factor(evidence_id, levels = levels(page$row_label))]
  section_rows[, section := factor(sample, levels = c('SUR1231', 'FFPED1', 'FFPEA1'))]
  page[, evidence := factor('Spatial log2 enrichment',
    levels = c('scATLAS Spearman rho', 'Spatial log2 enrichment'))]
  section_rows[, evidence := factor('Spatial log2 enrichment', levels = levels(page$evidence))]
  scrna_points <- copy(page)[, evidence := factor('scATLAS Spearman rho', levels = levels(page$evidence))]
  scrna_range <- data.table(evidence = factor('scATLAS Spearman rho', levels = levels(page$evidence)),
    value = c(-1, 1), row_label = page$row_label[1L])
  combined <- ggplot(page, aes(pooled_log2_enrichment, row_label)) +
    geom_vline(xintercept = 0, colour = 'grey75', linewidth = 0.3) +
    geom_blank(data = scrna_range, aes(x = value)) +
    geom_point(data = scrna_points, aes(x = scrna_spearman_r), size = 1.6, colour = '#444444') +
    geom_point(data = section_rows[eligible_pair %in% TRUE],
      aes(x = log2_enrichment, colour = section, shape = section),
      position = position_dodge(width = 0.55), size = 1.7, stroke = 0.5) +
    geom_point(data = page[spatial_available %in% TRUE], shape = 23,
      aes(fill = validation_class == 'Supported'), size = 2, colour = 'black', stroke = 0.5) +
    geom_text(data = page[!spatial_available %in% TRUE], x = 0, label = 'x', size = 7 / ggplot2::.pt) +
    scale_colour_manual(values = c(SUR1231 = '#0072B2', FFPED1 = '#D55E00', FFPEA1 = '#777777'),
      labels = sample_display[c('SUR1231', 'FFPED1', 'FFPEA1')], drop = FALSE, name = NULL) +
    scale_shape_manual(values = c(SUR1231 = 16, FFPED1 = 17, FFPEA1 = 1),
      labels = sample_display[c('SUR1231', 'FFPED1', 'FFPEA1')], drop = FALSE, name = NULL) +
    scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'white'),
      labels = c('FALSE' = 'Donor mean', 'TRUE' = 'scATLAS-supported'), name = NULL) +
    scale_y_discrete(labels = row_labels, limits = levels(page$row_label), drop = FALSE) +
    facet_grid(~evidence, scales = 'free_x') +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    labs(x = NULL, y = NULL, title = mode_titles[[mode]],
      subtitle = 'All sections contribute; A1/D1 share one donor. Section effects have no estimated confidence intervals.',
      caption = 'Black diamond: candidate-family FDR < 0.05 and both donors match scATLAS direction.\nx: fewer than two eligible donors. Whole T-cell evidence cannot distinguish CD4 from CD8.') +
    publication_theme + theme(legend.box = 'vertical')
  print(combined)
}
dev.off()

# Complete effect matrices retain non-significant tests. Grey means unavailable;
# a black ring marks FDR support with agreement across independent donors.
pdf(file.path(figure_dirs[['supporting']], 'Auto_primary_all_pairs.pdf'),
    width = figure_width, height = 170 / 25.4, onefile = TRUE, useDingbats = FALSE)
for (mode in mode_table$mode_id) {
  page <- copy(pooled_results[mode_id == mode & ring == primary_ring])
  source_levels <- if (mode_table[mode_id == mode, source_kind] == 'cancer_mp') cancer_mp_order else state_order
  target_levels <- unique(page$target_feature[order(target_plot_rank(page$target_kind, page$target_compartment, page$target_feature))])
  complete <- CJ(source_feature = source_levels, target_feature = target_levels, sorted = FALSE)
  page <- merge(complete, page, by = c('source_feature', 'target_feature'), all.x = TRUE)
  page[, source_axis := factor(source_feature, levels = rev(source_levels))]
  page[, target_axis := factor(target_feature, levels = target_levels)]
  page[n_donors_tested < 2 | is.na(n_donors_tested), pooled_log2_enrichment := NA_real_]
  panel <- ggplot(page, aes(target_axis, source_axis)) +
    geom_tile(aes(fill = pooled_log2_enrichment), colour = 'white', linewidth = 0.15) +
    geom_point(data = page[pooled_fdr < 0.05 & direction_concordant %in% TRUE],
      shape = 1, size = 1.3, stroke = 0.4) + colour_effect(1) +
    scale_x_discrete(labels = target_short, drop = FALSE) +
    scale_y_discrete(labels = source_short, drop = FALSE) +
    labs(title = mode_titles[[mode]], subtitle = 'Primary independent-donor mean; ring: FDR < 0.05 and matching directions',
      caption = 'Grey: fewer than two eligible donors. Colour is clipped at +/-1; exact effects remain in source tables.',
      x = NULL, y = NULL) + publication_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(panel)
}
dev.off()

# Maps are tissue context, not evidence of signalling. Use the 16-um lattice for
# physical axes and a scale bar; do not invent a high-resolution contact image.
make_spatial_plot <- function(candidate) {
  points <- rbindlist(lapply(c('SUR1231', 'FFPED1', 'FFPEA1'), function(sample_name) {
    membership <- memberships[[sample_name]]
    source_column <- paste(candidate$source_kind, candidate$source_feature, sep = '__')
    target_column <- paste(candidate$target_kind, candidate$target_feature, sep = '__')
    if (!source_column %in% colnames(membership$source_membership) ||
        !target_column %in% colnames(membership$target_membership[[candidate$target_kind]])) {
      stop('Supported map lacks required membership: ', sample_name)
    }
    background <- copy(membership$background)[, .(barcode, grid_row, grid_col)]
    background[, population := 'Other retained bins']
    source_ids <- membership$source_meta$barcode[as.numeric(membership$source_membership[, source_column]) > 0]
    target_ids <- membership$target_meta$barcode[as.numeric(membership$target_membership[[candidate$target_kind]][, target_column]) > 0]
    background[barcode %in% target_ids, population := 'TME target']
    background[barcode %in% source_ids, population := 'Malignant source']
    background[, `:=`(x_um = (grid_col - min(grid_col)) * 16,
      y_um = (grid_row - min(grid_row)) * 16,
      section = factor(sample_name, levels = c('SUR1231', 'FFPED1', 'FFPEA1')))]
    background
  }))
  points[, population := factor(population, levels = c('Other retained bins', 'TME target', 'Malignant source'))]
  setorder(points, population)
  bar_y <- max(points$y_um) + 160
  bars <- unique(points[, .(section)])[, `:=`(x_um = 0, y_um = bar_y)]
  ggplot(points, aes(x_um, y_um)) +
    geom_point(aes(colour = population), size = 0.23, shape = 15) +
    geom_segment(data = bars, aes(xend = 500, yend = y_um), linewidth = 0.6) +
    geom_text(data = bars, aes(x = 250, y = y_um + 240), label = '500 um', size = 7 / ggplot2::.pt) +
    scale_colour_manual(values = c('Other retained bins' = '#E5E5E5', 'TME target' = '#0072B2',
      'Malignant source' = '#D55E00'), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    facet_wrap(~section, nrow = 1, labeller = as_labeller(sample_display)) +
    scale_y_reverse() + coord_fixed() +
    labs(title = paste(source_short(candidate$source_feature), target_short(candidate$target_feature), sep = ' : '),
      subtitle = paste0('scATLAS rho = ', sprintf('%.2f', candidate$scrna_spearman_r),
        '; candidate-family FDR = ', format(candidate$spatial_validation_fdr, digits = 2)),
      caption = 'A1/D1 share one donor. Identical physical scales across panels. Proximity does not establish signalling.',
      x = NULL, y = NULL) + publication_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
}

validated_map_rows <- unique(validation[ring == primary_ring & validation_class == 'Supported'],
  by = c('mode_id', 'source_feature', 'target_feature'))
fwrite(validated_map_rows, file.path(output_tiers[['tables']], 'Auto_spatial_interaction_validated_map_summary.csv'))
pdf(file.path(figure_dirs[['final']], '02_supported_tissue_context.pdf'),
    width = figure_width, height = 100 / 25.4, onefile = TRUE, useDingbats = FALSE)
if (!nrow(validated_map_rows)) {
  print(empty_plot('scATLAS-supported tissue context', 'No candidate meets the independent-donor support rule'))
} else {
  for (i in seq_len(nrow(validated_map_rows))) {
    candidate <- validated_map_rows[i]
    print(make_spatial_plot(candidate))
  }
}
dev.off()

# Ring effects and pre/post substitution are sensitivity evidence, not additional
# biological replicates. Full statistics remain independently reproducible.
sensitivity <- merge(pooled_results[ring == primary_ring,
  .(mode_id, source_feature, target_feature, primary_effect = pooled_log2_enrichment)],
  pooled_results[ring != primary_ring], by = c('mode_id', 'source_feature', 'target_feature'))
sensitivity_plot <- ggplot(sensitivity, aes(primary_effect, pooled_log2_enrichment, colour = factor(ring))) +
  geom_hline(yintercept = 0, colour = 'grey80', linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = 'grey80', linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.5) +
  facet_wrap(~factor(mode_id, levels = mode_table$mode_id, labels = unname(mode_titles)), ncol = 2) +
  scale_colour_manual(values = c('2' = '#E69F00', '3' = '#0072B2'), labels = c('32 um', '48 um'), name = NULL) +
  labs(x = 'Primary 16-um window log2 enrichment', y = 'Broader window log2 enrichment') + publication_theme

ggsave(file.path(figure_dirs[['diagnostics']], 'Auto_window_sensitivity.pdf'), sensitivity_plot,
  width = figure_width, height = 135 / 25.4, device = 'pdf', useDingbats = FALSE)
fwrite(mp_grouping, file.path(output_tiers[['tables']], 'Auto_malignant_mp_label_key.csv'))
writeLines(c('# Current Visium HD spatial evidence', '',
  'All three sections contribute to two donors. A1 and D1 are a pre/post pair combined before pooling.', '',
  '- [scATLAS candidate evidence](../figures/Auto_final/01_scatlas_spatial_evidence.pdf): every selected candidate; individual sections and primary donor mean.',
  '- [Supported tissue context](../figures/Auto_final/02_supported_tissue_context.pdf): supported pairs only; physical scale bars.',
  '- [All-pair effects](../figures/Auto_supporting/Auto_primary_all_pairs.pdf): complete primary effects and explicit unavailable cells.',
  '- [Window sensitivity](../figures/Auto_diagnostics/Auto_window_sensitivity.pdf): cumulative 32/48-um windows against the primary 16-um window.',
  '- [Candidate source data](../tables/Auto_scatlas_figure_source.csv) and [MP descriptions](../tables/Auto_malignant_mp_label_key.csv).', '',
  'PDFs use a 183-mm canvas and 7-pt text. Final figures require the accompanying statistical limitations.',
  'Loose historical PDFs outside the three purpose folders are superseded and excluded from the current figure set.'),
  file.path(output_tiers[['reports']], 'CURRENT_RESULTS.md'))
####################

parameters <- data.table(
  parameter = c(
    "permutation_method", "n_permutations_total", "source_slot_rule", "inference_version", "input_fingerprint", "cache_version", "samples", "sample_alias", "bin_size_um", "primary_ring",
    "sensitivity_rings", "n_permutations", "mp_scoring", "mp_positive_rule", "ucell_max_rank",
    "top_genes_per_mp", "min_source_bins", "min_whole_target_bins",
    "min_parent_bins_for_mp", "min_mp_positive_bins", "min_observed_neighbour_slots",
    "spatial_permutation_block_bins", "spatial_permutation_block_um", "spatial_permutation_origins",
    "sample_fdr_cutoff", "pooled_fdr_cutoff", "sample_p_method", "pooled_method",
    "candidate_direction", "malignant_mp_plot_order", "state_plot_order",
    "whole_celltype_null", "tme_mp_null", "state_crosswalk"
  ),
  value = c(
    permutation_method, as.character(n_permutations * nrow(spatial_block_offsets)),
    "target_parent_for_tme_mp;all_tme_for_whole_celltype", inference_version, input_fingerprint, cache_version, paste(samples, collapse = ";"), "D1=FFPED1", "16", as.character(primary_ring),
    paste(neighbour_rings, collapse = ";"), as.character(n_permutations), "UCell on raw counts",
    paste0("sample/parent mean + ", mp_positive_sd_multiplier, " SD"), as.character(ucell_max_rank),
    as.character(top_genes_per_mp), as.character(min_source_bins), as.character(min_whole_target_bins),
    as.character(min_parent_bins_for_mp), as.character(min_mp_positive_bins), as.character(min_observed_neighbour_slots),
    as.character(spatial_block_bins), as.character(16L * spatial_block_bins),
    paste(apply(spatial_block_offsets, 1L, paste, collapse = ":"), collapse = ";"),
    as.character(sample_fdr_cutoff), as.character(pooled_fdr_cutoff), "two-sided empirical permutation p",
    "Bonferroni minimum within concordant donor; equal-weight signed Stouffer across donors",
    "scRNA Spearman rho", paste(cancer_mp_order, collapse = ";"), paste(state_order, collapse = ";"),
    "permute whole-celltype labels within shifted local tiles and sample",
    "permute joint MP-positive matrix within parent celltype, shifted local tile, and sample",
    paste(paste(names(state_crosswalk), state_crosswalk, sep = "->"), collapse = ";")
  )
)
fwrite(parameters, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_parameters.csv"))

compact_summary <- rbindlist(lapply(seq_len(nrow(mode_table)), function(idx) {
  mode_row <- mode_table[idx]
  sample_primary <- sample_results[mode_id == mode_row$mode_id & ring == primary_ring & eligible_pair %in% TRUE]
  pooled_primary <- pooled_results[mode_id == mode_row$mode_id & ring == primary_ring]
  data.table(
    mode_id = mode_row$mode_id, mode_label = mode_row$mode_label,
    n_sample_pair_tests = nrow(sample_primary),
    n_sample_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff, na.rm = TRUE),
    n_sample_enriched_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff & sample_primary$log2_enrichment > 0, na.rm = TRUE),
    n_sample_depleted_fdr_hits = sum(sample_primary$fdr_sample < sample_fdr_cutoff & sample_primary$log2_enrichment < 0, na.rm = TRUE),
    n_pooled_pair_tests = sum(is.finite(pooled_primary$pooled_p)),
    n_pooled_fdr_hits = sum(pooled_primary$pooled_fdr < pooled_fdr_cutoff, na.rm = TRUE),
    n_pooled_fdr_direction_concordant = sum(
      pooled_primary$pooled_fdr < pooled_fdr_cutoff & pooled_primary$direction_concordant,
      na.rm = TRUE
    ),
    n_recurrent_same_direction = sum(pooled_primary$direction_concordant & pooled_primary$n_donors_fdr_0_10 >= 2L, na.rm = TRUE),
    n_scatlas_supported_candidates = validation[
      mode_id == mode_row$mode_id & ring == primary_ring & validation_class == "Supported", .N
    ]
  )
}))
fwrite(compact_summary, file.path(summary_dir, "visium_hd_spatial_cancer_tme_interactions_summary.csv"))
fwrite(compact_summary, file.path(output_tiers[["tables"]], "Auto_spatial_interaction_compact_summary.csv"))

####################
# Keep the current index interpretable without relying on session audit notes.
result_lines <- c(
  '', '## Results and interpretation', '',
  '| Comparison | Eligible two-donor tests | Donor-concordant pooled FDR < 0.05 | scATLAS-supported candidates |',
  '|---|---:|---:|---:|',
  vapply(seq_len(nrow(compact_summary)), function(i) {
    row <- compact_summary[i]
    sprintf('| %s | %d | %d | %d |', row$mode_label, row$n_pooled_pair_tests,
      row$n_pooled_fdr_direction_concordant, row$n_scatlas_supported_candidates)
  }, character(1)), '',
  sprintf('The primary test counts adjacency edges within one 16-um grid step, including diagonal neighbours. Each of four shifted 320-um tile origins has %s permutations. The section test takes the largest origin-specific two-sided p value and requires matching nonzero effect directions across origins. TME-programme permutations also condition on the parent cell type. These nulls preserve local composition, but not arbitrary fine-scale spatial autocorrelation.', format(n_permutations, big.mark = ',')), '',
  'Both A1 and D1 contribute to their shared donor. Concordant section evidence is combined by a Bonferroni minimum-p rule; discordant paired directions contribute zero signed evidence. Equal-donor signed Stouffer pooling then combines the two independent donors. Plotted donor means are descriptive effects without estimated confidence intervals.', '',
  'Candidate support uses a separate BH family of unique scATLAS-selected spatial hypotheses and requires matching directions in both donors. Whole T-cell evidence cannot resolve CD4 from CD8. The source table retains unavailable candidates and the upstream correlation direction.', '',
  'The two-donor design limits generalisation. The single pre/post pair cannot establish a treatment effect. Bin proximity does not establish membrane contact, ligand-receptor signalling or causation.', '',
  '- [Primary numerical summary](../tables/Auto_spatial_interaction_compact_summary.csv)',
  '- [All section tests](../tables/Auto_spatial_interactions_by_sample.csv.gz)',
  '- [Donor evidence](../tables/Auto_spatial_interactions_by_donor.csv.gz)',
  '- [Paired-section substitution sensitivity](../tables/Auto_paired_section_substitution.csv.gz)'
)
cat(paste(result_lines, collapse = '\n'), '\n',
    file = file.path(output_tiers[['reports']], 'CURRENT_RESULTS.md'), append = TRUE)
####################

run_log <- c(
  run_log,
  paste0("n_sample_results=", nrow(sample_results)),
  paste0("n_pooled_results=", nrow(pooled_results)),
  paste0("n_sample_primary_fdr_hits=", sum(sample_results$ring == primary_ring & sample_results$fdr_sample < sample_fdr_cutoff, na.rm = TRUE)),
  paste0("n_pooled_primary_fdr_hits=", sum(pooled_results$ring == primary_ring & pooled_results$pooled_fdr < pooled_fdr_cutoff, na.rm = TRUE)),
  paste0("end=", format(Sys.time(), tz = "Europe/London")),
  "status=complete"
)
writeLines(run_log, file.path(output_tiers[["logs"]], "Auto_visium_hd_spatial_cancer_tme_interactions_run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_tiers[["logs"]], "Auto_visium_hd_spatial_cancer_tme_interactions_session_info.txt"))
####################
