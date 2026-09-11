####################
# Legacy interactive Visium HD annotation audit
#
# Status: legacy interactive audit
# Script: analysis/spatial/legacy_visiumhd/legacy_visium_hd_annotation_full_audit.R
# Description: Run and inspect one Visium HD sample and representation through
#   every production annotation stage without user-defined helper functions.
# Methodology:
#   analysis/methodology/spatial/visium_hd_final_annotation_methodology.md
# Inputs:
#   analysis/spatial/visium_hd_samples.tsv
#   ref_outs/EAC_Ref_merged.rds when rebuilding RCTD
#   ref_outs/visium_hd_outs/rctd/tables/Auto_<sample>_binned_rctd_annotations.csv.gz
# Outputs:
#   Objects and plots remain in the interactive R session by default.
#   Set WRITE_AUDIT_OUTPUTS <- TRUE to write under
#   ref_outs/visium_hd_outs/interactive_audit/.
# Run:
#   Open in RStudio using the dmtcp environment and execute one numbered block
#   at a time. RCTD reconstruction and whole-transcriptome DGE are heavy.
#
# This file deliberately contains no user-defined functions. The repeated code
# makes every transformation, threshold, and decision visible for auditing.
####################

####################
# 0. Packages and working directory
####################
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(leidenbase)
  library(Matrix)
  library(patchwork)
  library(Seurat)
  library(spacexr)
})

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)

####################
# 1. Select one sample and expose every adjustable parameter
####################

# Required selection.
SELECTED_SAMPLE <- getOption("visium.audit.sample", "SUR1231")
SELECTED_METHOD <- getOption(
  "visium.audit.method",
  "segmented"
) # "segmented" or "binned"

# Execution controls.
RUN_RCTD_FROM_SCRATCH <- FALSE
RUN_WHOLE_TRANSCRIPTOME_DGE <- FALSE
OPEN_INTERACTIVE_TABLES <- getOption(
  "visium.audit.open_tables",
  TRUE
)
WRITE_AUDIT_OUTPUTS <- getOption(
  "visium.audit.write_outputs",
  FALSE
)

# RCTD parameters. These matter only for SELECTED_METHOD == "binned".
RCTD_REFERENCE_PATH <- file.path(WD, "ref_outs", "EAC_Ref_merged.rds")
RCTD_MIN_UMIS <- 100L
RCTD_MAX_CORES <- 7L
RCTD_EPITHELIAL_REFERENCE_CAP <- 6000L
RCTD_OTHER_REFERENCE_CAP <- 2000L
RCTD_REFERENCE_SEED <- 666L
RCTD_DOUBLET_MODE <- "doublet"

# Shared observation QC.
ANNOTATION_MIN_COUNTS <- 100L
ANNOTATION_MAX_MT_PERCENT <- 15

# Normalisation and annotation graph.
NORMALISATION_SCALE_FACTOR <- 10000
N_HIGHLY_VARIABLE_GENES <- 3000L
ANNOTATION_N_PCS <- 50L
ANNOTATION_NEIGHBOURS <- 15L
ANNOTATION_NEIGHBOUR_METRIC <- "cosine"
ANNOTATION_LEIDEN_RESOLUTION <- 6
ANNOTATION_LEIDEN_PARTITION <- "RBConfigurationVertexPartition"
ANNOTATION_LEIDEN_ITERATIONS <- 2L
ANNOTATION_RANDOM_SEED <- 0L

# Marker evidence and local refinement.
MARKER_MIN_CLUSTER_CELLS <- 20L
MARKER_MIN_DETECTION_INCREASE <- 0.02
MARKER_REQUIRED_LONG_PANEL <- 2L
MARKER_REQUIRED_SHORT_PANEL <- 1L
SHORT_PANEL_MAX_GENES <- 2L
AMBIGUOUS_REFINEMENT_RESOLUTION <- 1

# Per-observation coexpression filtering.
COEXPRESSION_MIN_DETECTED_GENES <- 2L

# Clear legacy-style UMAP used only for display, never for annotation calls.
UMAP_DISPLAY_MIN_COUNTS <- 200L
UMAP_MIN_CELLS_PER_GENE <- 10L
UMAP_N_HIGHLY_VARIABLE_GENES <- 3000L
UMAP_N_PCS <- 50L
UMAP_NEIGHBOURS <- 15L
UMAP_NEIGHBOUR_METRIC <- "cosine"
UMAP_SCALE_MAX <- 10
UMAP_MIN_DIST <- 0.5
UMAP_SPREAD <- 1
UMAP_RANDOM_SEED <- 0L

# Score-only alternative audit. These settings do not alter the production
# annotation above; they let the alternative be assessed before replacement.
SCORE_ONLY_RESOLUTION_GRID <- c(0.5, 1, 2, 3, 4, 6)
SCORE_ONLY_SELECTED_RESOLUTION <- getOption(
  "visium.audit.score_resolution",
  6
)
SCORE_ONLY_MIN_CLUSTER_CELLS <- 5L
SCORE_ONLY_MARKER_MIN_POSITIVE_CELLS <- 2L
SCORE_ONLY_PROTECTED_MIN_STANDARDIZED_SCORE <- 0.25
SCORE_ONLY_REQUIRED_LONG_PANEL <- 2L
SCORE_ONLY_REQUIRED_SHORT_PANEL <- 1L

# Plot controls.
AUDIT_CELL_TYPE <- "t.cell"
SPATIAL_POINT_SIZE <- 0.35
UMAP_POINT_SIZE <- 0.45
PLOT_ALPHA <- 0.8
LEGEND_POINT_SIZE <- 4

markers <- list(
  erythrocyte = c("HBA1", "HBA2", "HBB"),
  keratinocyte = c("FLG", "IVL"),
  lymph = c("CCL21"),
  neutrophil = c("CTSG", "ELANE", "MPO", "AZU1"),
  endothelial = c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5"),
  epithelial = c("KRT7", "MUC1", "KRT19", "EPCAM"),
  fibroblast = c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"),
  b.cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  plasma = c("MZB1", "JCHAIN", "DERL3"),
  dendritic = c("CLEC10A", "CCR7", "CD86"),
  macrophage = c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"),
  mast = c("MS4A2", "CPA3", "TPSB2", "TPSAB1"),
  nk.cell = c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"),
  t.cell = c("CD3E", "CD3D", "CD2", "CD3G")
)

protected_types <- setdiff(names(markers), c("epithelial", "fibroblast"))
related_groups <- list(
  c("macrophage", "mast", "dendritic"),
  c("t.cell", "nk.cell", "dendritic"),
  c("b.cell", "plasma")
)

celltype_colours <- c(
  epithelial = "#D73027",
  fibroblast = "#8C564B",
  endothelial = "#1F78B4",
  macrophage = "#FF7F00",
  mast = "#A65628",
  t.cell = "#33A02C",
  b.cell = "#377EB8",
  nk.cell = "#984EA3",
  plasma = "#E377C2",
  dendritic = "#17BECF",
  lymph = "#6BAED6",
  erythrocyte = "#7F7F7F",
  keratinocyte = "#E6AB02",
  neutrophil = "#1B9E77",
  unresolved = "#BDBDBD"
)

parameter_audit <- data.frame(
  parameter = c(
    "SELECTED_SAMPLE", "SELECTED_METHOD", "RCTD_MIN_UMIS",
    "ANNOTATION_MIN_COUNTS", "ANNOTATION_MAX_MT_PERCENT",
    "NORMALISATION_SCALE_FACTOR", "N_HIGHLY_VARIABLE_GENES",
    "ANNOTATION_N_PCS", "ANNOTATION_NEIGHBOURS",
    "ANNOTATION_NEIGHBOUR_METRIC", "ANNOTATION_LEIDEN_RESOLUTION",
    "ANNOTATION_LEIDEN_PARTITION", "ANNOTATION_LEIDEN_ITERATIONS",
    "MARKER_MIN_CLUSTER_CELLS", "MARKER_MIN_DETECTION_INCREASE",
    "MARKER_REQUIRED_LONG_PANEL", "MARKER_REQUIRED_SHORT_PANEL",
    "AMBIGUOUS_REFINEMENT_RESOLUTION",
    "COEXPRESSION_MIN_DETECTED_GENES", "UMAP_DISPLAY_MIN_COUNTS",
    "UMAP_MIN_CELLS_PER_GENE", "UMAP_N_HIGHLY_VARIABLE_GENES",
    "UMAP_N_PCS", "UMAP_NEIGHBOURS", "UMAP_NEIGHBOUR_METRIC",
    "UMAP_SCALE_MAX", "UMAP_MIN_DIST", "UMAP_SPREAD",
    "SCORE_ONLY_RESOLUTION_GRID", "SCORE_ONLY_SELECTED_RESOLUTION",
    "SCORE_ONLY_MIN_CLUSTER_CELLS",
    "SCORE_ONLY_MARKER_MIN_POSITIVE_CELLS",
    "SCORE_ONLY_PROTECTED_MIN_STANDARDIZED_SCORE",
    "SCORE_ONLY_REQUIRED_LONG_PANEL",
    "SCORE_ONLY_REQUIRED_SHORT_PANEL"
  ),
  value = as.character(c(
    SELECTED_SAMPLE, SELECTED_METHOD, RCTD_MIN_UMIS,
    ANNOTATION_MIN_COUNTS, ANNOTATION_MAX_MT_PERCENT,
    NORMALISATION_SCALE_FACTOR, N_HIGHLY_VARIABLE_GENES,
    ANNOTATION_N_PCS, ANNOTATION_NEIGHBOURS,
    ANNOTATION_NEIGHBOUR_METRIC, ANNOTATION_LEIDEN_RESOLUTION,
    ANNOTATION_LEIDEN_PARTITION, ANNOTATION_LEIDEN_ITERATIONS,
    MARKER_MIN_CLUSTER_CELLS, MARKER_MIN_DETECTION_INCREASE,
    MARKER_REQUIRED_LONG_PANEL, MARKER_REQUIRED_SHORT_PANEL,
    AMBIGUOUS_REFINEMENT_RESOLUTION,
    COEXPRESSION_MIN_DETECTED_GENES, UMAP_DISPLAY_MIN_COUNTS,
    UMAP_MIN_CELLS_PER_GENE, UMAP_N_HIGHLY_VARIABLE_GENES,
    UMAP_N_PCS, UMAP_NEIGHBOURS, UMAP_NEIGHBOUR_METRIC,
    UMAP_SCALE_MAX, UMAP_MIN_DIST, UMAP_SPREAD,
    paste(SCORE_ONLY_RESOLUTION_GRID, collapse = ";"),
    SCORE_ONLY_SELECTED_RESOLUTION, SCORE_ONLY_MIN_CLUSTER_CELLS,
    SCORE_ONLY_MARKER_MIN_POSITIVE_CELLS,
    SCORE_ONLY_PROTECTED_MIN_STANDARDIZED_SCORE,
    SCORE_ONLY_REQUIRED_LONG_PANEL, SCORE_ONLY_REQUIRED_SHORT_PANEL
  )),
  stringsAsFactors = FALSE
)
print(parameter_audit, row.names = FALSE)

####################
# 2. Resolve the selected Space Ranger input
####################
manifest <- data.frame(
  sample = c("SUR1122", "SUR1231", "FFPEA1", "FFPED1"),
  binned_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"),
  segmented_input = c("/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs", "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs"),
  stringsAsFactors = FALSE
)

if (!SELECTED_SAMPLE %in% manifest$sample) {
  stop("SELECTED_SAMPLE is absent from hardcoded manifest")
}
if (!SELECTED_METHOD %in% c("binned", "segmented")) {
  stop("SELECTED_METHOD must be 'binned' or 'segmented'")
}

selected_manifest <- manifest[manifest$sample == SELECTED_SAMPLE, , drop = FALSE]
input_dir <- if (SELECTED_METHOD == "binned") {
  selected_manifest$binned_input[[1]]
} else {
  selected_manifest$segmented_input[[1]]
}

counts_h5 <- if (SELECTED_METHOD == "binned") {
  file.path(input_dir, "filtered_feature_bc_matrix.h5")
} else {
  file.path(input_dir, "filtered_feature_cell_matrix.h5")
}
if (!file.exists(counts_h5)) stop("Missing expression matrix: ", counts_h5)

counts <- Read10X_h5(counts_h5)
if (is.list(counts)) {
  counts <- if ("Gene Expression" %in% names(counts)) {
    counts[["Gene Expression"]]
  } else {
    counts[[1]]
  }
}
counts <- Matrix(counts, sparse = TRUE)
cat("Raw matrix:", nrow(counts), "genes x", ncol(counts), "observations\n")

####################
# 3. Read spatial coordinates for inspection
####################
if (SELECTED_METHOD == "binned") {
  positions_path <- file.path(input_dir, "spatial", "tissue_positions.parquet")
  if (!file.exists(positions_path)) stop("Missing positions: ", positions_path)
  spatial_coordinates <- as.data.frame(read_parquet(positions_path))
  barcode_column <- intersect(c("barcode", "Barcode"), colnames(spatial_coordinates))
  x_column <- intersect(
    c("pxl_col_in_fullres", "array_col"),
    colnames(spatial_coordinates)
  )
  y_column <- intersect(
    c("pxl_row_in_fullres", "array_row"),
    colnames(spatial_coordinates)
  )
  if (length(barcode_column) != 1L || length(x_column) == 0L ||
      length(y_column) == 0L) {
    stop("Unsupported binned positions columns")
  }
  spatial_coordinates <- data.frame(
    barcode = as.character(spatial_coordinates[[barcode_column[[1]]]]),
    pxl_col_in_fullres = as.numeric(spatial_coordinates[[x_column[[1]]]]),
    pxl_row_in_fullres = as.numeric(spatial_coordinates[[y_column[[1]]]]),
    stringsAsFactors = FALSE
  )
} else {
  # The dmtcp R environment does not contain sf. Coordinates are read from the
  # production annotation table, where they were calculated directly from the
  # Space Ranger cell-segmentation GeoJSON centroids. Expression and every
  # annotation decision below are recomputed independently from raw counts.
  production_annotation_path <- file.path(
    WD, "ref_outs", "visium_hd_outs", "tables",
    paste0("Auto_", SELECTED_SAMPLE, "_segmented_cell_annotations.csv.gz")
  )
  if (!file.exists(production_annotation_path)) {
    stop("Segmented centroid cache is missing: ", production_annotation_path)
  }
  production_coordinates <- fread(
    production_annotation_path,
    select = c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres"),
    data.table = FALSE
  )
  spatial_coordinates <- production_coordinates
  rm(production_coordinates)
}

rownames(spatial_coordinates) <- spatial_coordinates$barcode

####################
# 4. Optional binned RCTD reconstruction, or load the finalized RCTD calls
####################
rctd_result <- NULL
rctd_object <- NULL

if (SELECTED_METHOD == "binned" && RUN_RCTD_FROM_SCRATCH) {
  set.seed(RCTD_REFERENCE_SEED)
  reference_object <- readRDS(RCTD_REFERENCE_PATH)
  reference_label_column <- intersect(
    c("celltype_update", "celltype_manual", "celltype_group"),
    colnames(reference_object@meta.data)
  )
  if (length(reference_label_column) == 0L) {
    stop("No recognized cell-type column in RCTD reference")
  }
  reference_label_column <- reference_label_column[[1]]
  reference_labels <- as.character(
    reference_object@meta.data[[reference_label_column]]
  )
  names(reference_labels) <- rownames(reference_object@meta.data)
  reference_labels[reference_labels %in% c("t.cell", "nk.cell")] <- "t_nk.cell"
  valid_reference <- !is.na(reference_labels) &
    nzchar(reference_labels) &
    !reference_labels %in% c("unresolved", "unresolved_inconsistent")
  reference_labels <- reference_labels[valid_reference]

  sampled_reference_cells <- character()
  for (reference_type in sort(unique(reference_labels))) {
    candidate_cells <- names(reference_labels)[reference_labels == reference_type]
    reference_cap <- if (reference_type == "epithelial") {
      RCTD_EPITHELIAL_REFERENCE_CAP
    } else {
      RCTD_OTHER_REFERENCE_CAP
    }
    sampled_reference_cells <- c(
      sampled_reference_cells,
      sample(candidate_cells, min(length(candidate_cells), reference_cap))
    )
  }

  reference_counts <- LayerData(
    reference_object,
    assay = "RNA",
    layer = "counts"
  )
  reference_counts <- round(
    reference_counts[, sampled_reference_cells, drop = FALSE]
  )
  reference_cell_types <- factor(reference_labels[sampled_reference_cells])
  names(reference_cell_types) <- sampled_reference_cells
  rctd_reference <- Reference(
    reference_counts[, names(reference_cell_types), drop = FALSE],
    reference_cell_types,
    Matrix::colSums(reference_counts)
  )

  binned_umis <- Matrix::colSums(counts)
  rctd_keep <- binned_umis >= RCTD_MIN_UMIS
  rctd_counts <- counts[, rctd_keep, drop = FALSE]
  binned_umis <- binned_umis[rctd_keep]
  shared_barcodes <- intersect(colnames(rctd_counts), spatial_coordinates$barcode)
  rctd_counts <- rctd_counts[, shared_barcodes, drop = FALSE]
  binned_umis <- binned_umis[shared_barcodes]
  rctd_coordinates <- data.frame(
    x = spatial_coordinates[shared_barcodes, "pxl_col_in_fullres"],
    y = spatial_coordinates[shared_barcodes, "pxl_row_in_fullres"],
    row.names = shared_barcodes
  )

  rctd_puck <- SpatialRNA(rctd_coordinates, rctd_counts, binned_umis)
  rctd_object <- create.RCTD(
    rctd_puck,
    rctd_reference,
    max_cores = RCTD_MAX_CORES,
    test_mode = FALSE,
    UMI_min = RCTD_MIN_UMIS
  )
  rctd_object <- run.RCTD(
    rctd_object,
    doublet_mode = RCTD_DOUBLET_MODE
  )
  rctd_result <- as.data.frame(rctd_object@results$results_df)
  rctd_result$barcode <- rownames(rctd_result)
}

if (SELECTED_METHOD == "binned" && !RUN_RCTD_FROM_SCRATCH) {
  rctd_table_path <- file.path(
    WD, "ref_outs", "visium_hd_outs", "rctd", "tables",
    paste0("Auto_", SELECTED_SAMPLE, "_binned_rctd_annotations.csv.gz")
  )
  if (!file.exists(rctd_table_path)) {
    stop("Missing finalized RCTD table: ", rctd_table_path)
  }
  rctd_result <- fread(rctd_table_path, data.table = FALSE)
}

if (SELECTED_METHOD == "binned") {
  print(sort(table(rctd_result$spot_class, useNA = "ifany"), decreasing = TRUE))
  rctd_singlets <- as.character(rctd_result$spot_class) == "singlet"
  annotation_barcodes <- intersect(
    as.character(rctd_result$barcode[rctd_singlets]),
    colnames(counts)
  )
} else {
  annotation_barcodes <- colnames(counts)
}

####################
# 5. Construct the selected observation object and apply shared QC
####################
audit_object <- CreateSeuratObject(
  counts = counts[, annotation_barcodes, drop = FALSE],
  project = paste(SELECTED_SAMPLE, SELECTED_METHOD, sep = "_"),
  min.cells = 0,
  min.features = 0
)
audit_object[["percent.mt"]] <- PercentageFeatureSet(
  audit_object,
  pattern = "^MT-"
)

qc_before <- audit_object@meta.data
qc_before$barcode <- rownames(qc_before)

if (SELECTED_METHOD == "segmented") {
  qc_keep <- audit_object$nCount_RNA >= ANNOTATION_MIN_COUNTS &
    audit_object$percent.mt <= ANNOTATION_MAX_MT_PERCENT
  audit_object <- subset(audit_object, cells = colnames(audit_object)[qc_keep])
}

if (ncol(audit_object) < MARKER_MIN_CLUSTER_CELLS) {
  stop("Too few observations remain for cluster annotation")
}

raw_counts <- LayerData(audit_object, assay = "RNA", layer = "counts")
cat("Annotation universe:", ncol(audit_object), "observations\n")
print(summary(audit_object$nCount_RNA))
print(summary(audit_object$percent.mt))

####################
# 6. Log1p CP10K normalisation and raw normalized marker scores
####################
audit_object <- NormalizeData(
  audit_object,
  normalization.method = "LogNormalize",
  scale.factor = NORMALISATION_SCALE_FACTOR,
  verbose = FALSE
)
normalised_expression <- LayerData(
  audit_object,
  assay = "RNA",
  layer = "data"
)

marker_scores <- matrix(
  0,
  nrow = ncol(audit_object),
  ncol = length(markers),
  dimnames = list(colnames(audit_object), names(markers))
)

available_markers <- list()
for (cell_type in names(markers)) {
  available_genes <- intersect(markers[[cell_type]], rownames(audit_object))
  available_markers[[cell_type]] <- available_genes
  if (length(available_genes) > 0L) {
    marker_scores[, cell_type] <- Matrix::colMeans(
      normalised_expression[available_genes, , drop = FALSE]
    )
  }
  audit_object[[paste0(cell_type, "_score")]] <- marker_scores[, cell_type]
}

marker_score_summary <- as.data.frame(marker_scores) |>
  data.table::as.data.table(keep.rownames = "barcode")

if (OPEN_INTERACTIVE_TABLES) View(marker_score_summary)

####################
# 7. Annotation graph: HVGs, centred but unscaled PCA, cosine neighbours,
#    and deliberate resolution-6 overclustering
####################
set.seed(ANNOTATION_RANDOM_SEED)
audit_object <- FindVariableFeatures(
  audit_object,
  selection.method = "vst",
  nfeatures = min(N_HIGHLY_VARIABLE_GENES, nrow(audit_object)),
  verbose = FALSE
)
annotation_hvgs <- VariableFeatures(audit_object)

# do.scale = FALSE matches the production Scanpy PCA, which mean-centres but
# does not unit-variance scale the annotation HVGs.
audit_object <- ScaleData(
  audit_object,
  features = annotation_hvgs,
  do.center = TRUE,
  do.scale = FALSE,
  verbose = FALSE
)
annotation_npcs_used <- min(
  ANNOTATION_N_PCS,
  length(annotation_hvgs) - 1L,
  ncol(audit_object) - 1L
)
audit_object <- RunPCA(
  audit_object,
  features = annotation_hvgs,
  npcs = annotation_npcs_used,
  reduction.name = "pca_annotation",
  seed.use = ANNOTATION_RANDOM_SEED,
  verbose = FALSE
)
audit_object <- FindNeighbors(
  audit_object,
  reduction = "pca_annotation",
  dims = seq_len(min(40L, annotation_npcs_used)),
  k.param = min(ANNOTATION_NEIGHBOURS, ncol(audit_object) - 1L),
  annoy.metric = ANNOTATION_NEIGHBOUR_METRIC,
  graph.name = c("annotation_nn", "annotation_snn"),
  verbose = FALSE
)
annotation_igraph <- graph_from_adjacency_matrix(
  audit_object[["annotation_snn"]],
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
annotation_partition <- leiden_find_partition(
  igraph = annotation_igraph,
  partition_type = ANNOTATION_LEIDEN_PARTITION,
  edge_weights = E(annotation_igraph)$weight,
  seed = ANNOTATION_RANDOM_SEED + 1L,
  resolution_parameter = ANNOTATION_LEIDEN_RESOLUTION,
  num_iter = ANNOTATION_LEIDEN_ITERATIONS,
  verbose = FALSE
)
annotation_membership <- annotation_partition$membership - 1L
names(annotation_membership) <- V(annotation_igraph)$name
audit_object$Auto_manual_cluster <- as.character(
  annotation_membership[colnames(audit_object)]
)
audit_object$Auto_manual_cluster_parent <- audit_object$Auto_manual_cluster
audit_object$Auto_manual_cluster_refined <- FALSE

cluster_sizes_initial <- sort(
  table(audit_object$Auto_manual_cluster),
  decreasing = TRUE
)
print(cluster_sizes_initial)

####################
# 8. First-pass marker-by-marker cluster-versus-rest evidence
####################
all_marker_genes <- sort(unique(unlist(available_markers, use.names = FALSE)))
marker_counts <- raw_counts[all_marker_genes, , drop = FALSE]
observation_totals <- pmax(Matrix::colSums(raw_counts), 1)
marker_cp10k <- t(t(marker_counts) / observation_totals) *
  NORMALISATION_SCALE_FACTOR
marker_detected <- marker_counts > 0
cluster_vector <- as.character(audit_object$Auto_manual_cluster)

gene_evidence_rows <- list()
gene_evidence_index <- 0L

for (cluster_id in sort(unique(cluster_vector))) {
  in_cluster <- cluster_vector == cluster_id
  cluster_n <- sum(in_cluster)
  rest_n <- sum(!in_cluster)
  cluster_detected <- Matrix::rowSums(marker_detected[, in_cluster, drop = FALSE])
  rest_detected <- Matrix::rowSums(marker_detected[, !in_cluster, drop = FALSE])
  mean_cluster <- Matrix::rowMeans(marker_cp10k[, in_cluster, drop = FALSE])
  mean_rest <- Matrix::rowMeans(marker_cp10k[, !in_cluster, drop = FALSE])
  pct_cluster <- cluster_detected / cluster_n
  pct_rest <- rest_detected / rest_n
  marker_p_values <- phyper(
    cluster_detected - 1,
    Matrix::rowSums(marker_detected),
    ncol(marker_detected) - Matrix::rowSums(marker_detected),
    cluster_n,
    lower.tail = FALSE
  )
  marker_adjusted <- p.adjust(marker_p_values, method = "BH")
  names(marker_adjusted) <- all_marker_genes

  for (cell_type in names(available_markers)) {
    for (gene in available_markers[[cell_type]]) {
      gene_evidence_index <- gene_evidence_index + 1L
      gene_evidence_rows[[gene_evidence_index]] <- data.frame(
        Auto_manual_cluster = cluster_id,
        cell_type = cell_type,
        gene = gene,
        cluster_n = cluster_n,
        mean_cp10k_cluster = as.numeric(mean_cluster[[gene]]),
        mean_cp10k_rest = as.numeric(mean_rest[[gene]]),
        log2fc = log2(
          (as.numeric(mean_cluster[[gene]]) + 1) /
            (as.numeric(mean_rest[[gene]]) + 1)
        ),
        pct_cluster = as.numeric(pct_cluster[[gene]]),
        pct_rest = as.numeric(pct_rest[[gene]]),
        pct_delta = as.numeric(pct_cluster[[gene]] - pct_rest[[gene]]),
        p_value = as.numeric(marker_p_values[[gene]]),
        p_adjusted = as.numeric(marker_adjusted[[gene]]),
        supported = cluster_n >= MARKER_MIN_CLUSTER_CELLS &&
          as.numeric(pct_cluster[[gene]] - pct_rest[[gene]]) >=
            MARKER_MIN_DETECTION_INCREASE &&
          as.numeric(mean_cluster[[gene]]) > as.numeric(mean_rest[[gene]]),
        stringsAsFactors = FALSE
      )
    }
  }
}
gene_evidence_initial <- rbindlist(gene_evidence_rows, fill = TRUE)

cluster_score_table <- as.data.frame(marker_scores)
cluster_score_table$Auto_manual_cluster <- cluster_vector
cluster_score_means <- aggregate(
  cluster_score_table[, names(markers), drop = FALSE],
  by = list(Auto_manual_cluster = cluster_score_table$Auto_manual_cluster),
  FUN = mean
)
rownames(cluster_score_means) <- cluster_score_means$Auto_manual_cluster

type_evidence_rows <- list()
type_evidence_index <- 0L
for (cluster_id in sort(unique(gene_evidence_initial$Auto_manual_cluster))) {
  for (cell_type in names(markers)) {
    evidence_subset <- gene_evidence_initial[
      gene_evidence_initial$Auto_manual_cluster == cluster_id &
        gene_evidence_initial$cell_type == cell_type,
    ]
    if (nrow(evidence_subset) == 0L) next
    supported_subset <- evidence_subset[evidence_subset$supported, ]
    required_markers <- if (nrow(evidence_subset) <= SHORT_PANEL_MAX_GENES) {
      MARKER_REQUIRED_SHORT_PANEL
    } else {
      MARKER_REQUIRED_LONG_PANEL
    }
    type_evidence_index <- type_evidence_index + 1L
    type_evidence_rows[[type_evidence_index]] <- data.frame(
      Auto_manual_cluster = cluster_id,
      cell_type = cell_type,
      cluster_n = evidence_subset$cluster_n[[1]],
      n_markers_available = nrow(evidence_subset),
      n_markers_supported = nrow(supported_subset),
      required_markers = required_markers,
      marker_fraction_supported = nrow(supported_subset) / nrow(evidence_subset),
      marker_median_log2fc = if (nrow(supported_subset) > 0L) {
        median(supported_subset$log2fc)
      } else {
        NA_real_
      },
      marker_max_padjusted = if (nrow(supported_subset) > 0L) {
        max(supported_subset$p_adjusted)
      } else {
        NA_real_
      },
      supported_markers = paste(supported_subset$gene, collapse = ";"),
      passes_marker_evidence = nrow(supported_subset) >= required_markers,
      cluster_normalised_marker_score = cluster_score_means[
        cluster_id, cell_type
      ],
      stringsAsFactors = FALSE
    )
  }
}
type_evidence_initial <- rbindlist(type_evidence_rows, fill = TRUE)

if (OPEN_INTERACTIVE_TABLES) {
  View(gene_evidence_initial)
  View(type_evidence_initial)
}

####################
# 9. Identify clusters supporting multiple protected lineages
####################
protected_passing_initial <- type_evidence_initial[
  passes_marker_evidence & cell_type %in% protected_types
]
protected_counts_initial <- protected_passing_initial[
  , .(n_supported_protected_types = uniqueN(cell_type)),
  by = Auto_manual_cluster
]
ambiguous_clusters <- protected_counts_initial[
  n_supported_protected_types > 1,
  Auto_manual_cluster
]
print(protected_passing_initial[
  Auto_manual_cluster %in% ambiguous_clusters
])

####################
# 10. One low-resolution local split for ambiguous protected-lineage clusters
####################
for (parent_cluster in ambiguous_clusters) {
  parent_cells <- colnames(audit_object)[
    audit_object$Auto_manual_cluster == parent_cluster
  ]
  if (length(parent_cells) < 2L * MARKER_MIN_CLUSTER_CELLS) next

  refinement_object <- subset(audit_object, cells = parent_cells)
  refinement_object <- FindVariableFeatures(
    refinement_object,
    selection.method = "vst",
    nfeatures = min(N_HIGHLY_VARIABLE_GENES, nrow(refinement_object)),
    verbose = FALSE
  )
  refinement_hvgs <- VariableFeatures(refinement_object)
  refinement_object <- ScaleData(
    refinement_object,
    features = refinement_hvgs,
    do.center = TRUE,
    do.scale = FALSE,
    verbose = FALSE
  )
  refinement_npcs <- min(
    ANNOTATION_N_PCS,
    length(refinement_hvgs) - 1L,
    ncol(refinement_object) - 1L
  )
  if (refinement_npcs < 2L) next
  refinement_object <- RunPCA(
    refinement_object,
    features = refinement_hvgs,
    npcs = refinement_npcs,
    reduction.name = "pca_refinement",
    seed.use = ANNOTATION_RANDOM_SEED,
    verbose = FALSE
  )
  refinement_object <- FindNeighbors(
    refinement_object,
    reduction = "pca_refinement",
    dims = seq_len(min(40L, refinement_npcs)),
    k.param = min(ANNOTATION_NEIGHBOURS, ncol(refinement_object) - 1L),
    annoy.metric = ANNOTATION_NEIGHBOUR_METRIC,
    graph.name = c("refinement_nn", "refinement_snn"),
    verbose = FALSE
  )
  refinement_igraph <- graph_from_adjacency_matrix(
    refinement_object[["refinement_snn"]],
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  refinement_partition <- leiden_find_partition(
    igraph = refinement_igraph,
    partition_type = ANNOTATION_LEIDEN_PARTITION,
    edge_weights = E(refinement_igraph)$weight,
    seed = ANNOTATION_RANDOM_SEED + 1L,
    resolution_parameter = AMBIGUOUS_REFINEMENT_RESOLUTION,
    num_iter = ANNOTATION_LEIDEN_ITERATIONS,
    verbose = FALSE
  )
  refinement_membership <- refinement_partition$membership - 1L
  names(refinement_membership) <- V(refinement_igraph)$name
  refined_labels <- paste0(
    parent_cluster, ".",
    as.character(refinement_membership[colnames(refinement_object)])
  )
  names(refined_labels) <- colnames(refinement_object)
  audit_object$Auto_manual_cluster[parent_cells] <- refined_labels[parent_cells]
  audit_object$Auto_manual_cluster_refined[parent_cells] <- TRUE
}

cluster_vector <- as.character(audit_object$Auto_manual_cluster)
print(sort(table(cluster_vector), decreasing = TRUE))

####################
# 11. Recalculate marker evidence after local refinement
####################
gene_evidence_rows <- list()
gene_evidence_index <- 0L

for (cluster_id in sort(unique(cluster_vector))) {
  in_cluster <- cluster_vector == cluster_id
  cluster_n <- sum(in_cluster)
  rest_n <- sum(!in_cluster)
  cluster_detected <- Matrix::rowSums(marker_detected[, in_cluster, drop = FALSE])
  rest_detected <- Matrix::rowSums(marker_detected[, !in_cluster, drop = FALSE])
  mean_cluster <- Matrix::rowMeans(marker_cp10k[, in_cluster, drop = FALSE])
  mean_rest <- Matrix::rowMeans(marker_cp10k[, !in_cluster, drop = FALSE])
  pct_cluster <- cluster_detected / cluster_n
  pct_rest <- rest_detected / rest_n
  marker_p_values <- phyper(
    cluster_detected - 1,
    Matrix::rowSums(marker_detected),
    ncol(marker_detected) - Matrix::rowSums(marker_detected),
    cluster_n,
    lower.tail = FALSE
  )
  marker_adjusted <- p.adjust(marker_p_values, method = "BH")
  names(marker_adjusted) <- all_marker_genes

  for (cell_type in names(available_markers)) {
    for (gene in available_markers[[cell_type]]) {
      gene_evidence_index <- gene_evidence_index + 1L
      gene_evidence_rows[[gene_evidence_index]] <- data.frame(
        Auto_manual_cluster = cluster_id,
        cell_type = cell_type,
        gene = gene,
        cluster_n = cluster_n,
        mean_cp10k_cluster = as.numeric(mean_cluster[[gene]]),
        mean_cp10k_rest = as.numeric(mean_rest[[gene]]),
        log2fc = log2(
          (as.numeric(mean_cluster[[gene]]) + 1) /
            (as.numeric(mean_rest[[gene]]) + 1)
        ),
        pct_cluster = as.numeric(pct_cluster[[gene]]),
        pct_rest = as.numeric(pct_rest[[gene]]),
        pct_delta = as.numeric(pct_cluster[[gene]] - pct_rest[[gene]]),
        p_value = as.numeric(marker_p_values[[gene]]),
        p_adjusted = as.numeric(marker_adjusted[[gene]]),
        supported = cluster_n >= MARKER_MIN_CLUSTER_CELLS &&
          as.numeric(pct_cluster[[gene]] - pct_rest[[gene]]) >=
            MARKER_MIN_DETECTION_INCREASE &&
          as.numeric(mean_cluster[[gene]]) > as.numeric(mean_rest[[gene]]),
        stringsAsFactors = FALSE
      )
    }
  }
}
gene_evidence <- rbindlist(gene_evidence_rows, fill = TRUE)

cluster_score_table$Auto_manual_cluster <- cluster_vector
cluster_score_means <- aggregate(
  cluster_score_table[, names(markers), drop = FALSE],
  by = list(Auto_manual_cluster = cluster_score_table$Auto_manual_cluster),
  FUN = mean
)
rownames(cluster_score_means) <- cluster_score_means$Auto_manual_cluster

type_evidence_rows <- list()
type_evidence_index <- 0L
for (cluster_id in sort(unique(gene_evidence$Auto_manual_cluster))) {
  for (cell_type in names(markers)) {
    current_cell_type <- cell_type
    evidence_subset <- gene_evidence[
      Auto_manual_cluster == cluster_id & cell_type == current_cell_type
    ]
    if (nrow(evidence_subset) == 0L) next
    supported_subset <- evidence_subset[supported == TRUE]
    required_markers <- if (nrow(evidence_subset) <= SHORT_PANEL_MAX_GENES) {
      MARKER_REQUIRED_SHORT_PANEL
    } else {
      MARKER_REQUIRED_LONG_PANEL
    }
    type_evidence_index <- type_evidence_index + 1L
    type_evidence_rows[[type_evidence_index]] <- data.frame(
      Auto_manual_cluster = cluster_id,
      cell_type = cell_type,
      cluster_n = evidence_subset$cluster_n[[1]],
      n_markers_available = nrow(evidence_subset),
      n_markers_supported = nrow(supported_subset),
      required_markers = required_markers,
      marker_fraction_supported = nrow(supported_subset) / nrow(evidence_subset),
      marker_median_log2fc = if (nrow(supported_subset) > 0L) {
        median(supported_subset$log2fc)
      } else {
        NA_real_
      },
      marker_max_padjusted = if (nrow(supported_subset) > 0L) {
        max(supported_subset$p_adjusted)
      } else {
        NA_real_
      },
      supported_markers = paste(supported_subset$gene, collapse = ";"),
      passes_marker_evidence = nrow(supported_subset) >= required_markers,
      cluster_normalised_marker_score = cluster_score_means[
        cluster_id, cell_type
      ],
      stringsAsFactors = FALSE
    )
  }
}
type_evidence <- rbindlist(type_evidence_rows, fill = TRUE)

####################
# 12. Apply the hierarchical cell-type assignment explicitly
####################
cluster_assignment_rows <- list()
cluster_assignment_index <- 0L

for (cluster_id in sort(unique(type_evidence$Auto_manual_cluster))) {
  cluster_evidence <- type_evidence[Auto_manual_cluster == cluster_id]
  passing <- cluster_evidence[passes_marker_evidence == TRUE]
  protected_passing <- passing[cell_type %in% protected_types]

  if (nrow(protected_passing) > 0L) {
    candidates <- protected_passing
    evidence_stage <- "non_epithelial_non_fibroblast_dge"
  } else if (nrow(passing[cell_type == "fibroblast"]) > 0L) {
    candidates <- passing[cell_type == "fibroblast"]
    evidence_stage <- "fibroblast_dge"
  } else {
    candidates <- cluster_evidence[
      cell_type %in% c("epithelial", "fibroblast")
    ]
    evidence_stage <- "epithelial_fibroblast_normalised_score"
  }

  if (grepl("_dge$", evidence_stage)) {
    candidate_order <- order(
      -candidates$n_markers_supported,
      -candidates$marker_fraction_supported,
      -candidates$marker_median_log2fc,
      -candidates$cluster_normalised_marker_score,
      candidates$cell_type,
      na.last = TRUE
    )
  } else {
    candidate_order <- order(
      -candidates$cluster_normalised_marker_score,
      candidates$cell_type
    )
  }
  selected_evidence <- candidates[candidate_order[[1]]]

  cluster_assignment_index <- cluster_assignment_index + 1L
  cluster_assignment_rows[[cluster_assignment_index]] <- data.frame(
    Auto_manual_cluster = cluster_id,
    Auto_annotation_celltype = selected_evidence$cell_type,
    Auto_annotation_evidence_stage = evidence_stage,
    Auto_annotation_supported_markers = selected_evidence$supported_markers,
    Auto_annotation_n_markers_supported =
      selected_evidence$n_markers_supported,
    Auto_annotation_required_markers = selected_evidence$required_markers,
    Auto_annotation_marker_median_log2fc =
      selected_evidence$marker_median_log2fc,
    Auto_annotation_marker_max_padjusted =
      selected_evidence$marker_max_padjusted,
    Auto_annotation_cluster_normalised_score =
      selected_evidence$cluster_normalised_marker_score,
    stringsAsFactors = FALSE
  )
}
cluster_assignments <- rbindlist(cluster_assignment_rows, fill = TRUE)

assignment_lookup <- setNames(
  cluster_assignments$Auto_annotation_celltype,
  cluster_assignments$Auto_manual_cluster
)
celltype_assignment_values <- assignment_lookup[
  audit_object$Auto_manual_cluster
]
names(celltype_assignment_values) <- colnames(audit_object)
audit_object$Auto_annotation_celltype <- celltype_assignment_values
audit_object$Auto_annotation_celltype_pre_filter <-
  audit_object$Auto_annotation_celltype

for (assignment_column in setdiff(
  colnames(cluster_assignments),
  c("Auto_manual_cluster", "Auto_annotation_celltype")
)) {
  assignment_values <- setNames(
    cluster_assignments[[assignment_column]],
    cluster_assignments$Auto_manual_cluster
  )
  cell_assignment_values <- assignment_values[
    audit_object$Auto_manual_cluster
  ]
  names(cell_assignment_values) <- colnames(audit_object)
  audit_object[[assignment_column]] <- cell_assignment_values
}

print(sort(table(audit_object$Auto_annotation_celltype), decreasing = TRUE))

####################
# 13. Per-observation incompatible coexpression filtering
####################
active_panel_matrix <- matrix(
  FALSE,
  nrow = ncol(audit_object),
  ncol = length(markers),
  dimnames = list(colnames(audit_object), names(markers))
)
for (cell_type in names(available_markers)) {
  genes <- available_markers[[cell_type]]
  if (length(genes) >= COEXPRESSION_MIN_DETECTED_GENES) {
    active_panel_matrix[, cell_type] <- Matrix::colSums(
      raw_counts[genes, , drop = FALSE] > 0
    ) >= COEXPRESSION_MIN_DETECTED_GENES
  }
}

active_marker_types <- character(ncol(audit_object))
coexpression_status <- rep("compatible", ncol(audit_object))

for (cell_index in seq_len(ncol(audit_object))) {
  active_types <- colnames(active_panel_matrix)[active_panel_matrix[cell_index, ]]
  active_marker_types[[cell_index]] <- paste(active_types, collapse = "|")
  incompatible <- FALSE

  if (length(active_types) >= 2L) {
    for (first_index in seq_len(length(active_types) - 1L)) {
      for (second_index in seq.int(first_index + 1L, length(active_types))) {
        first_type <- active_types[[first_index]]
        second_type <- active_types[[second_index]]
        pair_allowed <- first_type %in% c("epithelial", "fibroblast") ||
          second_type %in% c("epithelial", "fibroblast")
        for (related_group in related_groups) {
          if (all(c(first_type, second_type) %in% related_group)) {
            pair_allowed <- TRUE
          }
        }
        if (!pair_allowed) incompatible <- TRUE
      }
    }
  }
  if (incompatible) coexpression_status[[cell_index]] <-
    "incompatible_coexpression"
}

names(active_marker_types) <- colnames(audit_object)
names(coexpression_status) <- colnames(audit_object)
audit_object$Auto_annotation_active_marker_types <- active_marker_types
audit_object$Auto_annotation_coexpression <- coexpression_status
audit_object$Auto_annotation_pass_doublet_filter <-
  coexpression_status == "compatible"

incompatible_cells <- colnames(audit_object)[
  coexpression_status == "incompatible_coexpression"
]
if (length(incompatible_cells) > 0L) {
  replacement_labels <- active_marker_types[incompatible_cells]
  replacement_labels[!nzchar(replacement_labels)] <- "unresolved"
  audit_object$Auto_annotation_celltype[incompatible_cells] <-
    replacement_labels
}
audit_object$Auto_annotation_keep_epithelial <-
  audit_object$Auto_annotation_pass_doublet_filter &
  audit_object$Auto_annotation_celltype == "epithelial"

cat("Coexpression filter:\n")
print(table(audit_object$Auto_annotation_coexpression))
cat("Final exact labels among passing observations:\n")
print(sort(
  table(audit_object$Auto_annotation_celltype[
    audit_object$Auto_annotation_pass_doublet_filter
  ]),
  decreasing = TRUE
))

####################
# 14. Optional whole-transcriptome cluster DGE
####################
whole_transcriptome_dge <- NULL
if (RUN_WHOLE_TRANSCRIPTOME_DGE) {
  Idents(audit_object) <- audit_object$Auto_manual_cluster
  rankable_clusters <- names(table(Idents(audit_object)))[
    table(Idents(audit_object)) >= 2L
  ]
  rankable_cells <- colnames(audit_object)[
    as.character(Idents(audit_object)) %in% rankable_clusters
  ]
  rankable_object <- subset(audit_object, cells = rankable_cells)
  Idents(rankable_object) <- rankable_object$Auto_manual_cluster
  whole_transcriptome_dge <- FindAllMarkers(
    rankable_object,
    only.pos = TRUE,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1,
    verbose = FALSE
  )
  whole_transcriptome_dge <- as.data.table(whole_transcriptome_dge)
  setorder(whole_transcriptome_dge, cluster, -avg_log2FC, p_val_adj)
  whole_transcriptome_dge <- whole_transcriptome_dge[
    , head(.SD, 100L),
    by = cluster
  ]
  if (OPEN_INTERACTIVE_TABLES) View(whole_transcriptome_dge)
}

####################
# 15. Clear legacy-style display UMAP, independent of annotation clustering
####################
display_cells <- colnames(audit_object)[
  audit_object$nCount_RNA >= UMAP_DISPLAY_MIN_COUNTS &
    audit_object$percent.mt <= ANNOTATION_MAX_MT_PERCENT
]
if (length(display_cells) < MARKER_MIN_CLUSTER_CELLS) {
  warning("Too few cells pass display QC; using all annotated observations")
  display_cells <- colnames(audit_object)
}

display_object <- subset(audit_object, cells = display_cells)
display_gene_detection <- Matrix::rowSums(
  LayerData(display_object, assay = "RNA", layer = "counts") > 0
)
display_genes <- names(display_gene_detection)[
  display_gene_detection >= UMAP_MIN_CELLS_PER_GENE
]
display_object <- subset(display_object, features = display_genes)
display_object <- FindVariableFeatures(
  display_object,
  selection.method = "vst",
  nfeatures = min(UMAP_N_HIGHLY_VARIABLE_GENES, nrow(display_object)),
  verbose = FALSE
)
display_hvgs <- VariableFeatures(display_object)
display_object <- ScaleData(
  display_object,
  features = display_hvgs,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  scale.max = UMAP_SCALE_MAX,
  verbose = FALSE
)
display_npcs_used <- min(
  UMAP_N_PCS,
  length(display_hvgs) - 1L,
  ncol(display_object) - 1L
)
display_object <- RunPCA(
  display_object,
  features = display_hvgs,
  npcs = display_npcs_used,
  reduction.name = "pca_display",
  seed.use = UMAP_RANDOM_SEED,
  verbose = FALSE
)
display_object <- FindNeighbors(
  display_object,
  reduction = "pca_display",
  dims = seq_len(min(40L, display_npcs_used)),
  k.param = min(UMAP_NEIGHBOURS, ncol(display_object) - 1L),
  annoy.metric = UMAP_NEIGHBOUR_METRIC,
  graph.name = c("display_nn", "display_snn"),
  verbose = FALSE
)
display_object <- RunUMAP(
  display_object,
  reduction = "pca_display",
  dims = seq_len(min(40L, display_npcs_used)),
  reduction.name = "umap_display",
  n.neighbors = min(UMAP_NEIGHBOURS, ncol(display_object) - 1L),
  min.dist = UMAP_MIN_DIST,
  spread = UMAP_SPREAD,
  metric = UMAP_NEIGHBOUR_METRIC,
  seed.use = UMAP_RANDOM_SEED,
  verbose = FALSE
)

display_embedding <- Embeddings(display_object, "umap_display")
display_object$Auto_annotation_celltype <-
  audit_object$Auto_annotation_celltype[colnames(display_object)]
display_object$Auto_annotation_pass_doublet_filter <-
  audit_object$Auto_annotation_pass_doublet_filter[colnames(display_object)]

####################
# 16. Spatial, UMAP, marker-score, and marker-evidence audit plots
####################
annotation_table <- audit_object@meta.data
annotation_table$barcode <- rownames(annotation_table)
annotation_table <- merge(
  annotation_table,
  spatial_coordinates,
  by = "barcode",
  all.x = TRUE,
  sort = FALSE
)

plot_labels <- sort(unique(annotation_table$Auto_annotation_celltype))
missing_colours <- setdiff(plot_labels, names(celltype_colours))
if (length(missing_colours) > 0L) {
  celltype_colours <- c(
    celltype_colours,
    setNames(rep("#636363", length(missing_colours)), missing_colours)
  )
}

spatial_plot <- ggplot(
  annotation_table[
    annotation_table$Auto_annotation_pass_doublet_filter, ,
    drop = FALSE
  ],
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    colour = Auto_annotation_celltype
  )
) +
  geom_point(size = SPATIAL_POINT_SIZE, alpha = PLOT_ALPHA) +
  scale_colour_manual(values = celltype_colours, drop = FALSE) +
  scale_y_reverse() +
  coord_equal() +
  labs(
    title = paste(SELECTED_SAMPLE, SELECTED_METHOD, "spatial"),
    colour = "Cell type"
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.55, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = LEGEND_POINT_SIZE)))

umap_plot_data <- data.frame(
  barcode = rownames(display_embedding),
  UMAP_1 = display_embedding[, 1],
  UMAP_2 = display_embedding[, 2],
  Auto_annotation_celltype =
    audit_object$Auto_annotation_celltype[rownames(display_embedding)],
  Auto_annotation_pass_doublet_filter =
    audit_object$Auto_annotation_pass_doublet_filter[rownames(display_embedding)],
  stringsAsFactors = FALSE
)

umap_plot <- ggplot(
  umap_plot_data[
    umap_plot_data$Auto_annotation_pass_doublet_filter, ,
    drop = FALSE
  ],
  aes(x = UMAP_1, y = UMAP_2, colour = Auto_annotation_celltype)
) +
  geom_point(size = UMAP_POINT_SIZE, alpha = PLOT_ALPHA) +
  scale_colour_manual(values = celltype_colours, drop = FALSE) +
  coord_equal() +
  labs(
    title = paste(SELECTED_SAMPLE, "legacy-style display UMAP"),
    colour = "Cell type"
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.55, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = LEGEND_POINT_SIZE)))

selected_marker_score <- marker_scores[, AUDIT_CELL_TYPE]
marker_distribution_data <- data.frame(
  barcode = rownames(marker_scores),
  marker_score = selected_marker_score,
  cluster = audit_object$Auto_manual_cluster[rownames(marker_scores)],
  final_cell_type =
    audit_object$Auto_annotation_celltype[rownames(marker_scores)],
  stringsAsFactors = FALSE
)

marker_density_plot <- ggplot(
  marker_distribution_data,
  aes(x = marker_score, colour = final_cell_type)
) +
  geom_density(linewidth = 0.8, adjust = 0.8) +
  scale_colour_manual(values = celltype_colours, drop = FALSE) +
  labs(
    title = paste(AUDIT_CELL_TYPE, "normalized marker-score distribution"),
    x = "Mean log1p CP10K marker score",
    y = "Density",
    colour = "Final cell type"
  ) +
  theme_classic(base_size = 13)

selected_marker_evidence <- gene_evidence[
  cell_type == AUDIT_CELL_TYPE
]
marker_evidence_plot <- ggplot(
  selected_marker_evidence,
  aes(
    x = reorder(
      paste(Auto_manual_cluster, gene, sep = " | "),
      pct_delta
    ),
    y = pct_delta,
    fill = supported
  )
) +
  geom_col() +
  geom_hline(
    yintercept = MARKER_MIN_DETECTION_INCREASE,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#238B45", `FALSE` = "#BDBDBD")) +
  labs(
    title = paste(AUDIT_CELL_TYPE, "cluster marker evidence"),
    x = "Cluster | marker",
    y = "Detection increase over rest",
    fill = "Supported"
  ) +
  theme_classic(base_size = 11)

print(spatial_plot)
print(umap_plot)
print(marker_density_plot)
print(marker_evidence_plot)

####################
# 17. Score-only alternative: resolution, score, support, and call diagnostics
#
# Candidate identity is the highest cluster-average within-sample standardized
# marker-panel score. Marker detection is reported as an independent quality
# flag, not used to manufacture a higher relative score. No DGE or local
# refinement is used in this branch.
####################
display_snn_igraph <- graph_from_adjacency_matrix(
  display_object[["display_snn"]],
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

resolution_memberships <- data.frame(
  barcode = colnames(display_object),
  stringsAsFactors = FALSE
)
resolution_cluster_sizes <- list()
resolution_cluster_size_index <- 0L

for (candidate_resolution in SCORE_ONLY_RESOLUTION_GRID) {
  candidate_partition <- leiden_find_partition(
    igraph = display_snn_igraph,
    partition_type = ANNOTATION_LEIDEN_PARTITION,
    edge_weights = E(display_snn_igraph)$weight,
    seed = ANNOTATION_RANDOM_SEED + 1L,
    resolution_parameter = candidate_resolution,
    num_iter = ANNOTATION_LEIDEN_ITERATIONS,
    verbose = FALSE
  )
  candidate_membership <- candidate_partition$membership - 1L
  names(candidate_membership) <- V(display_snn_igraph)$name
  resolution_column <- paste0("resolution_", candidate_resolution)
  resolution_memberships[[resolution_column]] <- as.character(
    candidate_membership[resolution_memberships$barcode]
  )
  candidate_sizes <- as.data.frame(
    table(resolution_memberships[[resolution_column]]),
    stringsAsFactors = FALSE
  )
  colnames(candidate_sizes) <- c("cluster", "cluster_n")
  candidate_sizes$resolution <- candidate_resolution
  resolution_cluster_size_index <- resolution_cluster_size_index + 1L
  resolution_cluster_sizes[[resolution_cluster_size_index]] <- candidate_sizes
}
resolution_cluster_sizes <- rbindlist(
  resolution_cluster_sizes,
  fill = TRUE
)

resolution_plot_data <- merge(
  data.frame(
    barcode = rownames(display_embedding),
    UMAP_1 = display_embedding[, 1],
    UMAP_2 = display_embedding[, 2],
    stringsAsFactors = FALSE
  ),
  resolution_memberships,
  by = "barcode",
  all.x = TRUE,
  sort = FALSE
)
resolution_plot_data <- melt(
  as.data.table(resolution_plot_data),
  id.vars = c("barcode", "UMAP_1", "UMAP_2"),
  variable.name = "resolution",
  value.name = "cluster"
)
resolution_plot_data$resolution <- sub(
  "^resolution_",
  "Resolution ",
  resolution_plot_data$resolution
)

resolution_umap_plot <- ggplot(
  resolution_plot_data,
  aes(x = UMAP_1, y = UMAP_2, colour = cluster)
) +
  geom_point(size = 0.18, alpha = 0.75) +
  facet_wrap(~resolution, ncol = 3) +
  coord_equal() +
  guides(colour = "none") +
  labs(title = "Clear UMAP: cluster resolution comparison") +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14)
  )

resolution_size_plot <- ggplot(
  resolution_cluster_sizes,
  aes(
    x = factor(resolution, levels = SCORE_ONLY_RESOLUTION_GRID),
    y = cluster_n
  )
) +
  geom_boxplot(outlier.size = 1, width = 0.7) +
  geom_hline(
    yintercept = SCORE_ONLY_MIN_CLUSTER_CELLS,
    linetype = "dashed",
    colour = "#D73027",
    linewidth = 0.7
  ) +
  scale_y_log10() +
  labs(
    title = "Cluster-size distribution by resolution",
    x = "Leiden resolution",
    y = "Cells per cluster (log10)"
  ) +
  theme_classic(base_size = 13)

selected_resolution_column <- paste0(
  "resolution_",
  SCORE_ONLY_SELECTED_RESOLUTION
)
if (!selected_resolution_column %in% colnames(resolution_memberships)) {
  stop("SCORE_ONLY_SELECTED_RESOLUTION must occur in SCORE_ONLY_RESOLUTION_GRID")
}
score_only_cluster <- resolution_memberships[[selected_resolution_column]]
names(score_only_cluster) <- resolution_memberships$barcode

display_marker_scores <- marker_scores[
  resolution_memberships$barcode,
  ,
  drop = FALSE
]
display_score_means <- colMeans(display_marker_scores)
display_score_sds <- apply(display_marker_scores, 2, sd)
display_score_sds[!is.finite(display_score_sds) | display_score_sds == 0] <- 1
display_standardized_scores <- sweep(
  display_marker_scores,
  2,
  display_score_means,
  "-"
)
display_standardized_scores <- sweep(
  display_standardized_scores,
  2,
  display_score_sds,
  "/"
)

score_only_cell_table <- data.frame(
  barcode = rownames(display_marker_scores),
  score_only_cluster = score_only_cluster[rownames(display_marker_scores)],
  display_marker_scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
score_only_z_cell_table <- data.frame(
  barcode = rownames(display_standardized_scores),
  score_only_cluster =
    score_only_cluster[rownames(display_standardized_scores)],
  display_standardized_scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

score_only_cluster_raw <- aggregate(
  score_only_cell_table[, names(markers), drop = FALSE],
  by = list(score_only_cluster = score_only_cell_table$score_only_cluster),
  FUN = mean
)
score_only_cluster_z <- aggregate(
  score_only_z_cell_table[, names(markers), drop = FALSE],
  by = list(score_only_cluster = score_only_z_cell_table$score_only_cluster),
  FUN = mean
)
rownames(score_only_cluster_raw) <- score_only_cluster_raw$score_only_cluster
rownames(score_only_cluster_z) <- score_only_cluster_z$score_only_cluster

score_only_raw_long <- melt(
  as.data.table(score_only_cluster_raw),
  id.vars = "score_only_cluster",
  variable.name = "cell_type",
  value.name = "mean_log1p_cp10k_score"
)
score_only_z_long <- melt(
  as.data.table(score_only_cluster_z),
  id.vars = "score_only_cluster",
  variable.name = "cell_type",
  value.name = "mean_standardized_score"
)

raw_score_heatmap <- ggplot(
  score_only_raw_long,
  aes(
    x = score_only_cluster,
    y = cell_type,
    fill = mean_log1p_cp10k_score
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#2166AC") +
  labs(
    title = "Cluster-average raw normalized marker scores",
    x = "Cluster",
    y = NULL,
    fill = "Mean log1p\nCP10K"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

z_score_limit <- max(
  abs(score_only_z_long$mean_standardized_score),
  na.rm = TRUE
)
standardized_score_heatmap <- ggplot(
  score_only_z_long,
  aes(
    x = score_only_cluster,
    y = cell_type,
    fill = mean_standardized_score
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-z_score_limit, z_score_limit)
  ) +
  labs(
    title = "Cluster-average within-sample standardized marker scores",
    x = "Cluster",
    y = NULL,
    fill = "Mean z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

display_raw_counts <- raw_counts[
  ,
  resolution_memberships$barcode,
  drop = FALSE
]
display_normalised_expression <- normalised_expression[
  ,
  resolution_memberships$barcode,
  drop = FALSE
]

score_only_marker_rows <- list()
score_only_marker_index <- 0L
for (cluster_id in sort(unique(score_only_cluster))) {
  cluster_cells <- names(score_only_cluster)[score_only_cluster == cluster_id]
  cluster_n <- length(cluster_cells)
  for (cell_type in names(available_markers)) {
    for (gene in available_markers[[cell_type]]) {
      score_only_marker_index <- score_only_marker_index + 1L
      score_only_marker_rows[[score_only_marker_index]] <- data.frame(
        score_only_cluster = cluster_id,
        cluster_n = cluster_n,
        cell_type = cell_type,
        gene = gene,
        positive_cells = sum(
          display_raw_counts[gene, cluster_cells, drop = TRUE] > 0
        ),
        detection_fraction = mean(
          display_raw_counts[gene, cluster_cells, drop = TRUE] > 0
        ),
        mean_log1p_cp10k = mean(
          display_normalised_expression[gene, cluster_cells, drop = TRUE]
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
score_only_marker_audit <- rbindlist(score_only_marker_rows, fill = TRUE)
score_only_marker_audit[
  ,
  marker_detected_for_audit :=
    positive_cells >= SCORE_ONLY_MARKER_MIN_POSITIVE_CELLS
]

score_only_panel_audit <- score_only_marker_audit[
  ,
  .(
    cluster_n = first(cluster_n),
    n_markers_available = .N,
    n_markers_detected = sum(marker_detected_for_audit),
    detected_markers = paste(
      gene[marker_detected_for_audit],
      collapse = ";"
    ),
    minimum_positive_cells = as.integer(
      if (any(marker_detected_for_audit)) {
        min(positive_cells[marker_detected_for_audit])
      } else {
        0L
      }
    ),
    mean_gene_detection_fraction = mean(detection_fraction),
    mean_gene_log1p_cp10k = mean(mean_log1p_cp10k)
  ),
  by = .(score_only_cluster, cell_type)
]
score_only_panel_audit[
  ,
  required_markers := ifelse(
    n_markers_available <= SHORT_PANEL_MAX_GENES,
    SCORE_ONLY_REQUIRED_SHORT_PANEL,
    SCORE_ONLY_REQUIRED_LONG_PANEL
  )
]
score_only_panel_audit[
  ,
  passes_absolute_marker_audit :=
    cluster_n >= SCORE_ONLY_MIN_CLUSTER_CELLS &
    n_markers_detected >= required_markers
]
score_only_panel_audit <- merge(
  score_only_panel_audit,
  score_only_z_long,
  by = c("score_only_cluster", "cell_type"),
  all.x = TRUE
)
score_only_panel_audit <- merge(
  score_only_panel_audit,
  score_only_raw_long,
  by = c("score_only_cluster", "cell_type"),
  all.x = TRUE
)

score_only_call_rows <- list()
score_only_call_index <- 0L
for (cluster_id in sort(unique(score_only_panel_audit$score_only_cluster))) {
  cluster_panels <- score_only_panel_audit[
    score_only_cluster == cluster_id
  ]
  setorder(
    cluster_panels,
    -mean_standardized_score,
    -mean_log1p_cp10k_score,
    cell_type
  )
  top_panel <- cluster_panels[1]
  second_panel <- cluster_panels[2]
  supported_panels <- cluster_panels[
    passes_absolute_marker_audit == TRUE
  ]
  supported_protected <- supported_panels[
    cell_type %in% protected_types &
      mean_standardized_score >=
        SCORE_ONLY_PROTECTED_MIN_STANDARDIZED_SCORE
  ]
  if (nrow(supported_protected) > 0L) {
    setorder(
      supported_protected,
      -mean_standardized_score,
      -mean_log1p_cp10k_score,
      cell_type
    )
    supported_top_panel <- supported_protected[1]
    supported_assignment_stage <- "protected_marker_score"
  } else {
    supported_epithelial_fibroblast <- supported_panels[
      cell_type %in% c("epithelial", "fibroblast")
    ]
    if (nrow(supported_epithelial_fibroblast) == 0L) {
      supported_epithelial_fibroblast <- cluster_panels[
        cell_type %in% c("epithelial", "fibroblast") &
          mean_log1p_cp10k_score > 0
      ]
      supported_assignment_stage <-
        "epithelial_fibroblast_positive_raw_fallback"
    } else {
      supported_assignment_stage <-
        "epithelial_fibroblast_marker_score"
    }
    setorder(
      supported_epithelial_fibroblast,
      -mean_log1p_cp10k_score,
      -mean_standardized_score,
      cell_type
    )
    supported_top_panel <- supported_epithelial_fibroblast[1]
  }
  if (nrow(supported_top_panel) > 0L) {
    supported_score_only_celltype <- supported_top_panel$cell_type
    supported_top_standardized_score <-
      supported_top_panel$mean_standardized_score
    supported_top_raw_score <-
      supported_top_panel$mean_log1p_cp10k_score
    supported_detected_markers <- supported_top_panel$detected_markers
  } else {
    supported_score_only_celltype <- "unresolved"
    supported_assignment_stage <- "unresolved_no_replicated_marker_support"
    supported_top_standardized_score <- NA_real_
    supported_top_raw_score <- NA_real_
    supported_detected_markers <- ""
  }
  score_only_call_index <- score_only_call_index + 1L
  score_only_call_rows[[score_only_call_index]] <- data.frame(
    score_only_cluster = cluster_id,
    cluster_n = top_panel$cluster_n,
    score_only_celltype = top_panel$cell_type,
    top_standardized_score = top_panel$mean_standardized_score,
    second_celltype = second_panel$cell_type,
    second_standardized_score = second_panel$mean_standardized_score,
    top_second_gap =
      top_panel$mean_standardized_score -
      second_panel$mean_standardized_score,
    top_raw_score = top_panel$mean_log1p_cp10k_score,
    n_markers_detected = top_panel$n_markers_detected,
    required_markers = top_panel$required_markers,
    detected_markers = top_panel$detected_markers,
    passes_absolute_marker_audit =
      top_panel$passes_absolute_marker_audit,
    supported_score_only_celltype = supported_score_only_celltype,
    supported_top_standardized_score =
      supported_top_standardized_score,
    supported_top_raw_score = supported_top_raw_score,
    supported_detected_markers = supported_detected_markers,
    supported_assignment_stage = supported_assignment_stage,
    stringsAsFactors = FALSE
  )
}
score_only_cluster_calls <- rbindlist(score_only_call_rows, fill = TRUE)

score_only_cluster_order <- score_only_cluster_calls[
  order(
    supported_score_only_celltype,
    -cluster_n,
    score_only_cluster
  )
]
score_only_cluster_order[
  ,
  cluster_display_label := paste0(
    "C", score_only_cluster,
    " | ", supported_score_only_celltype,
    " | n=", cluster_n
  )
]
score_only_display_lookup <- setNames(
  score_only_cluster_order$cluster_display_label,
  score_only_cluster_order$score_only_cluster
)
score_only_display_levels <-
  score_only_cluster_order$cluster_display_label

score_only_raw_plot_data <- copy(score_only_raw_long)
score_only_raw_plot_data[
  ,
  cluster_display_label :=
    score_only_display_lookup[score_only_cluster]
]
score_only_raw_plot_data[
  ,
  cluster_display_label := factor(
    cluster_display_label,
    levels = score_only_display_levels
  )
]
raw_score_heatmap <- ggplot(
  score_only_raw_plot_data,
  aes(
    x = cluster_display_label,
    y = cell_type,
    fill = mean_log1p_cp10k_score
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#2166AC") +
  labs(
    title = "Cluster-average raw normalized marker scores",
    x = "Cluster | final assignment | cell count",
    y = NULL,
    fill = "Mean log1p\nCP10K"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

score_only_z_plot_data <- copy(score_only_z_long)
score_only_z_plot_data[
  ,
  cluster_display_label :=
    score_only_display_lookup[score_only_cluster]
]
score_only_z_plot_data[
  ,
  cluster_display_label := factor(
    cluster_display_label,
    levels = score_only_display_levels
  )
]
standardized_score_heatmap <- ggplot(
  score_only_z_plot_data,
  aes(
    x = cluster_display_label,
    y = cell_type,
    fill = mean_standardized_score
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-z_score_limit, z_score_limit)
  ) +
  labs(
    title = "Cluster-average within-sample standardized marker scores",
    x = "Cluster | final assignment | cell count",
    y = NULL,
    fill = "Mean z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

score_only_call_lookup <- setNames(
  score_only_cluster_calls$score_only_celltype,
  score_only_cluster_calls$score_only_cluster
)
score_only_gap_lookup <- setNames(
  score_only_cluster_calls$top_second_gap,
  score_only_cluster_calls$score_only_cluster
)
score_only_support_lookup <- setNames(
  score_only_cluster_calls$supported_score_only_celltype != "unresolved",
  score_only_cluster_calls$score_only_cluster
)
supported_score_only_call_lookup <- setNames(
  score_only_cluster_calls$supported_score_only_celltype,
  score_only_cluster_calls$score_only_cluster
)

score_only_plot_data <- data.frame(
  barcode = rownames(display_embedding),
  UMAP_1 = display_embedding[, 1],
  UMAP_2 = display_embedding[, 2],
  score_only_cluster = score_only_cluster[rownames(display_embedding)],
  score_only_celltype = score_only_call_lookup[
    score_only_cluster[rownames(display_embedding)]
  ],
  supported_score_only_celltype = supported_score_only_call_lookup[
    score_only_cluster[rownames(display_embedding)]
  ],
  top_second_gap = score_only_gap_lookup[
    score_only_cluster[rownames(display_embedding)]
  ],
  passes_absolute_marker_audit = score_only_support_lookup[
    score_only_cluster[rownames(display_embedding)]
  ],
  stringsAsFactors = FALSE
)
score_only_cluster_centroids <- aggregate(
  score_only_plot_data[, c("UMAP_1", "UMAP_2")],
  by = list(score_only_cluster = score_only_plot_data$score_only_cluster),
  FUN = median
)

score_only_labels <- sort(unique(score_only_plot_data$score_only_celltype))
score_only_labels <- sort(unique(c(
  score_only_labels,
  score_only_plot_data$supported_score_only_celltype
)))
score_only_missing_colours <- setdiff(
  score_only_labels,
  names(celltype_colours)
)
if (length(score_only_missing_colours) > 0L) {
  celltype_colours <- c(
    celltype_colours,
    setNames(
      rep("#636363", length(score_only_missing_colours)),
      score_only_missing_colours
    )
  )
}

score_only_umap_plot <- ggplot(
  score_only_plot_data,
  aes(x = UMAP_1, y = UMAP_2, colour = score_only_celltype)
) +
  geom_point(size = UMAP_POINT_SIZE, alpha = PLOT_ALPHA) +
  geom_label(
    data = score_only_cluster_centroids,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = score_only_cluster
    ),
    inherit.aes = FALSE,
    size = 2.5,
    label.size = 0,
    label.padding = unit(0.08, "lines"),
    fill = scales::alpha("white", 0.75)
  ) +
  scale_colour_manual(values = celltype_colours, drop = FALSE) +
  coord_equal() +
  labs(
    title = paste0(
      "Score-only candidate calls | resolution ",
      SCORE_ONLY_SELECTED_RESOLUTION
    ),
    colour = "Top score"
  ) +
  theme_void(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = LEGEND_POINT_SIZE)))

supported_score_only_umap_plot <- ggplot(
  score_only_plot_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    colour = supported_score_only_celltype
  )
) +
  geom_point(size = UMAP_POINT_SIZE, alpha = PLOT_ALPHA) +
  geom_label(
    data = score_only_cluster_centroids,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = score_only_cluster
    ),
    inherit.aes = FALSE,
    size = 2.5,
    label.size = 0,
    label.padding = unit(0.08, "lines"),
    fill = scales::alpha("white", 0.75)
  ) +
  scale_colour_manual(values = celltype_colours, drop = FALSE) +
  coord_equal() +
  labs(
    title = "Supported score-only calls; unsupported clusters unresolved",
    colour = "Supported top score"
  ) +
  theme_void(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = LEGEND_POINT_SIZE)))

score_gap_umap_plot <- ggplot(
  score_only_plot_data,
  aes(x = UMAP_1, y = UMAP_2, colour = top_second_gap)
) +
  geom_point(size = UMAP_POINT_SIZE, alpha = PLOT_ALPHA) +
  scale_colour_gradient(low = "#F7FBFF", high = "#08306B") +
  coord_equal() +
  labs(
    title = "Top-versus-second standardized score gap",
    colour = "Score gap"
  ) +
  theme_void(base_size = 13)

absolute_support_umap_plot <- ggplot(
  score_only_plot_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    colour = passes_absolute_marker_audit
  )
) +
  geom_point(size = UMAP_POINT_SIZE, alpha = PLOT_ALPHA) +
  scale_colour_manual(
    values = c(`TRUE` = "#238B45", `FALSE` = "#D73027"),
    labels = c(`TRUE` = "Supported", `FALSE` = "Weak")
  ) +
  coord_equal() +
  labs(
    title = "Absolute marker-support audit of top score",
    colour = NULL
  ) +
  theme_void(base_size = 13)

panel_support_plot <- ggplot(
  score_only_panel_audit,
  aes(
    x = score_only_cluster,
    y = cell_type,
    size = n_markers_detected,
    colour = mean_log1p_cp10k_score
  )
) +
  geom_point() +
  scale_colour_gradient(low = "#F7FBFF", high = "#08519C") +
  scale_size_continuous(range = c(0.5, 6)) +
  labs(
    title = "Panel support: raw score and number of detected markers",
    x = "Cluster",
    y = NULL,
    size = "Markers\n>= detection floor",
    colour = "Mean log1p\nCP10K score"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

selected_type_marker_audit <- score_only_marker_audit[
  cell_type == AUDIT_CELL_TYPE
]
selected_type_marker_dotplot <- ggplot(
  selected_type_marker_audit,
  aes(
    x = score_only_cluster,
    y = gene,
    size = detection_fraction,
    colour = mean_log1p_cp10k
  )
) +
  geom_point() +
  scale_size_continuous(range = c(0.3, 7)) +
  scale_colour_gradient(low = "#F7FBFF", high = "#CB181D") +
  labs(
    title = paste(AUDIT_CELL_TYPE, "individual-marker audit"),
    x = "Cluster",
    y = NULL,
    size = "Detection\nfraction",
    colour = "Mean log1p\nCP10K"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

selected_type_feature_rows <- list()
selected_type_feature_index <- 0L
for (gene in available_markers[[AUDIT_CELL_TYPE]]) {
  selected_type_feature_index <- selected_type_feature_index + 1L
  selected_type_feature_rows[[selected_type_feature_index]] <- data.frame(
    barcode = resolution_memberships$barcode,
    UMAP_1 = display_embedding[resolution_memberships$barcode, 1],
    UMAP_2 = display_embedding[resolution_memberships$barcode, 2],
    gene = gene,
    expression = as.numeric(
      display_normalised_expression[
        gene,
        resolution_memberships$barcode,
        drop = TRUE
      ]
    ),
    stringsAsFactors = FALSE
  )
}
selected_type_feature_data <- rbindlist(
  selected_type_feature_rows,
  fill = TRUE
)
selected_type_feature_umap <- ggplot(
  selected_type_feature_data,
  aes(x = UMAP_1, y = UMAP_2, colour = expression)
) +
  geom_point(size = 0.25, alpha = 0.8) +
  facet_wrap(~gene, ncol = 2) +
  scale_colour_gradient(low = "#F0F0F0", high = "#CB181D") +
  coord_equal() +
  labs(
    title = paste(AUDIT_CELL_TYPE, "individual marker expression"),
    colour = "log1p\nCP10K"
  ) +
  theme_void(base_size = 12)

selected_type_spatial_data <- annotation_table[
  match(rownames(marker_scores), annotation_table$barcode),
  c("barcode", "pxl_col_in_fullres", "pxl_row_in_fullres"),
  drop = FALSE
]
selected_type_spatial_data$marker_score <- marker_scores[
  selected_type_spatial_data$barcode,
  AUDIT_CELL_TYPE
]
selected_type_spatial_score_plot <- ggplot(
  selected_type_spatial_data,
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    colour = marker_score
  )
) +
  geom_point(size = SPATIAL_POINT_SIZE, alpha = PLOT_ALPHA) +
  scale_colour_gradient(low = "#F0F0F0", high = "#CB181D") +
  scale_y_reverse() +
  coord_equal() +
  labs(
    title = paste(AUDIT_CELL_TYPE, "raw normalized spatial score"),
    colour = "Mean log1p\nCP10K"
  ) +
  theme_void(base_size = 13)

selected_type_cluster_distribution <- data.frame(
  barcode = resolution_memberships$barcode,
  score_only_cluster = score_only_cluster[resolution_memberships$barcode],
  marker_score = marker_scores[
    resolution_memberships$barcode,
    AUDIT_CELL_TYPE
  ],
  score_only_celltype = score_only_call_lookup[
    score_only_cluster[resolution_memberships$barcode]
  ],
  stringsAsFactors = FALSE
)
selected_type_cluster_boxplot <- ggplot(
  selected_type_cluster_distribution,
  aes(
    x = reorder(score_only_cluster, marker_score, FUN = median),
    y = marker_score,
    fill = score_only_celltype
  )
) +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.35) +
  scale_fill_manual(values = celltype_colours, drop = FALSE) +
  labs(
    title = paste(AUDIT_CELL_TYPE, "raw score distribution within clusters"),
    x = "Cluster ordered by median score",
    y = "Mean log1p CP10K marker score",
    fill = "Top score call"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

current_display_labels <- audit_object$Auto_annotation_celltype[
  resolution_memberships$barcode
]
score_only_transition_table <- as.data.table(data.frame(
  current_annotation = current_display_labels,
  score_only_annotation = score_only_call_lookup[
    score_only_cluster[resolution_memberships$barcode]
  ],
  stringsAsFactors = FALSE
))[
  ,
  .N,
  by = .(current_annotation, score_only_annotation)
][order(-N)]

supported_score_only_transition_table <- as.data.table(data.frame(
  current_annotation = current_display_labels,
  supported_score_only_annotation = supported_score_only_call_lookup[
    score_only_cluster[resolution_memberships$barcode]
  ],
  stringsAsFactors = FALSE
))[
  ,
  .N,
  by = .(current_annotation, supported_score_only_annotation)
][order(-N)]

print(resolution_umap_plot)
print(resolution_size_plot)
print(raw_score_heatmap)
print(standardized_score_heatmap)
print(score_only_umap_plot)
print(supported_score_only_umap_plot)
print(score_gap_umap_plot)
print(absolute_support_umap_plot)
print(panel_support_plot)
print(selected_type_marker_dotplot)
print(selected_type_feature_umap)
print(selected_type_spatial_score_plot)
print(selected_type_cluster_boxplot)
print(score_only_cluster_calls)
print(score_only_transition_table)
print(supported_score_only_transition_table)

####################
# 18. Optional audit-only output; production outputs are never overwritten
####################
if (WRITE_AUDIT_OUTPUTS) {
  audit_output_dir <- file.path(
    WD, "ref_outs", "visium_hd_outs", "interactive_audit",
    paste(SELECTED_SAMPLE, SELECTED_METHOD, sep = "_")
  )
  dir.create(audit_output_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    parameter_audit,
    file.path(audit_output_dir, "annotation_parameters.csv")
  )
  fwrite(
    annotation_table,
    file.path(audit_output_dir, "cell_annotations.csv.gz"),
    compress = "gzip"
  )
  fwrite(
    gene_evidence,
    file.path(audit_output_dir, "cluster_marker_evidence.csv")
  )
  fwrite(
    type_evidence,
    file.path(audit_output_dir, "cluster_celltype_evidence.csv")
  )
  fwrite(
    cluster_assignments,
    file.path(audit_output_dir, "cluster_assignments.csv")
  )
  fwrite(
    resolution_cluster_sizes,
    file.path(audit_output_dir, "score_only_resolution_cluster_sizes.csv")
  )
  fwrite(
    score_only_marker_audit,
    file.path(audit_output_dir, "score_only_marker_audit.csv")
  )
  fwrite(
    score_only_panel_audit,
    file.path(audit_output_dir, "score_only_panel_audit.csv")
  )
  fwrite(
    score_only_cluster_calls,
    file.path(audit_output_dir, "score_only_cluster_calls.csv")
  )
  fwrite(
    score_only_transition_table,
    file.path(audit_output_dir, "current_vs_score_only_transitions.csv")
  )
  fwrite(
    supported_score_only_transition_table,
    file.path(
      audit_output_dir,
      "current_vs_supported_score_only_transitions.csv"
    )
  )
  if (!is.null(whole_transcriptome_dge)) {
    fwrite(
      whole_transcriptome_dge,
      file.path(audit_output_dir, "cluster_top_dge.csv.gz"),
      compress = "gzip"
    )
  }
  ggsave(
    file.path(audit_output_dir, "spatial_annotation.pdf"),
    spatial_plot,
    width = 11,
    height = 9,
    useDingbats = FALSE
  )
  ggsave(
    file.path(audit_output_dir, "umap_annotation.pdf"),
    umap_plot,
    width = 11,
    height = 9,
    useDingbats = FALSE
  )
  ggsave(
    file.path(audit_output_dir, "marker_score_and_evidence.pdf"),
    marker_density_plot / marker_evidence_plot,
    width = 14,
    height = 14,
    useDingbats = FALSE
  )
  ggsave(
    file.path(audit_output_dir, "score_only_resolution_diagnostics.pdf"),
    resolution_umap_plot / resolution_size_plot,
    width = 16,
    height = 14,
    useDingbats = FALSE
  )
  ggsave(
    file.path(audit_output_dir, "score_only_score_heatmaps.pdf"),
    raw_score_heatmap / standardized_score_heatmap,
    width = max(18, 0.48 * nrow(score_only_cluster_calls)),
    height = 14,
    useDingbats = FALSE
  )
  ggsave(
    file.path(audit_output_dir, "score_only_call_diagnostics.pdf"),
    score_only_umap_plot + supported_score_only_umap_plot +
      score_gap_umap_plot + absolute_support_umap_plot +
      panel_support_plot,
    width = 18,
    height = 20,
    useDingbats = FALSE
  )
  ggsave(
    file.path(
      audit_output_dir,
      paste0("score_only_", AUDIT_CELL_TYPE, "_diagnostics.pdf")
    ),
    selected_type_feature_umap /
      (selected_type_marker_dotplot | selected_type_spatial_score_plot) /
      selected_type_cluster_boxplot,
    width = 18,
    height = 20,
    useDingbats = FALSE
  )
}

cat(
  "\nAudit complete. Inspect these main objects:\n",
  "  audit_object\n",
  "  parameter_audit\n",
  "  marker_scores\n",
  "  gene_evidence_initial / type_evidence_initial\n",
  "  ambiguous_clusters\n",
  "  gene_evidence / type_evidence\n",
  "  cluster_assignments\n",
  "  spatial_plot / umap_plot\n",
  "  marker_density_plot / marker_evidence_plot\n",
  "  resolution_umap_plot / resolution_size_plot\n",
  "  raw_score_heatmap / standardized_score_heatmap\n",
  "  score_only_cluster_calls / score_only_transition_table\n",
  "  supported_score_only_transition_table\n",
  "  score_only_umap_plot / supported_score_only_umap_plot\n",
  "  score_gap_umap_plot\n",
  "  absolute_support_umap_plot / panel_support_plot\n",
  "  selected_type_feature_umap / selected_type_marker_dotplot\n",
  "  selected_type_spatial_score_plot / selected_type_cluster_boxplot\n",
  sep = ""
)
####################
