####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/final_mp_scenic.R
#   Methodology: analysis/methodology/cell_states/final_mp_scenic_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# final_mp_scenic.R
#
# Final-MP-focused SCENIC workflow for OAC epithelial cells.
# Uses the curated final MP panel:
#   - the final 17 centred refined MPs, including three cell-cycle-associated MPs
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#
# Output:
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_selected_cells.csv
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_selected_cells.pdf
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_assignment_summary.csv
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_gene_sets.rds
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_gene_membership.csv
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulon_auc.rds
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_rss.rds
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulon_heatmap.pdf
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_network.pdf
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_network_edges.csv
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulon_targets.csv
#   updates/new_updates/summaries/Auto_final_mp_scenic_summary.csv
#
# Usage:
#   Rscript analysis/cell_states/final_mp_scenic.R
#   Rscript analysis/cell_states/final_mp_scenic.R prepare_only=true
#   Rscript analysis/cell_states/final_mp_scenic.R db_dir=/path/to/cistarget
####################

####################
# Additional live outputs:
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulons_by_mp.xlsx
#   ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulons_by_state.xlsx
# Each workbook contains an all-regulon overview and ranked per-MP/per-state
# sheets with SCENIC RSS, mean AUCell activity, and regulon-target summaries.
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(data.table)
library(scales)
library(igraph)
library(ggraph)
library(tidygraph)
library(grid)
library(doParallel)
library(foreach)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(x[1])) {
    return(y)
  }
  x[1]
}

parse_args <- function(args) {
  if (length(args) == 0) {
    return(list())
  }
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    out[[key]] <- value
  }
  out
}

to_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "1", "yes", "y")
}



format_regulon_name <- function(x) {
  x <- gsub("_extended$", "", x)
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  x <- gsub(" \\([0-9]+ genes\\)$", "", x)
  x
}

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(
      across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

get_assay_matrix <- function(seurat_obj, slot_name = c("counts", "data")) {
  slot_name <- match.arg(slot_name)

  mat <- tryCatch(
    GetAssayData(seurat_obj, assay = "RNA", slot = slot_name),
    error = function(e) NULL
  )
  if (!is.null(mat)) return(mat)

  mat <- tryCatch(
    LayerData(seurat_obj, assay = "RNA", layer = slot_name),
    error = function(e) NULL
  )
  if (!is.null(mat)) return(mat)

  mat <- tryCatch(
    LayerData(seurat_obj[["RNA"]], layer = slot_name),
    error = function(e) NULL
  )
  if (!is.null(mat)) return(mat)

  assay_obj <- seurat_obj@assays$RNA
  mat <- tryCatch(slot(assay_obj, slot_name), error = function(e) NULL)
  if (!is.null(mat)) return(mat)

  stop("Unable to retrieve RNA ", slot_name, " matrix from Seurat object.")
}

extract_regulon_targets <- function(x) {
  if (requireNamespace("GSEABase", quietly = TRUE) && methods::is(x, "GeneSet")) {
    return(unique(GSEABase::geneIds(x)))
  }
  if (is.character(x)) {
    return(unique(x))
  }
  if (is.list(x) && !is.null(x$gene)) {
    return(unique(as.character(x$gene)))
  }
  if (!is.null(names(x))) {
    return(unique(names(x)))
  }
  unique(as.character(x))
}

detect_db_files <- function(db_dir) {
  if (!dir.exists(db_dir)) {
    stop(
      "SCENIC database directory not found: ", db_dir, "\n",
      "Set SCENIC_DB_DIR or pass db_dir=/path/to/cistarget."
    )
  }

  db_files <- list.files(db_dir, pattern = "\\.feather$", full.names = FALSE)
  db_files <- db_files[grepl("hg38|refseq-r80|hgnc", db_files, ignore.case = TRUE)]
  preferred_db_files <- db_files[grepl("mc9nr|refseq-r80", db_files, ignore.case = TRUE)]
  if (length(preferred_db_files) > 0) {
    db_files <- preferred_db_files
  }
  if (length(db_files) == 0) {
    stop(
      "No human cisTarget feather databases found in ", db_dir, ".\n",
      "Expected files such as hg38/refseq-r80 500bp and 10kb rankings."
    )
  }

  picked <- unique(c(
    db_files[grepl("500bp", db_files, ignore.case = TRUE)][1],
    db_files[grepl("10kb", db_files, ignore.case = TRUE)][1],
    db_files
  ))
  picked <- picked[!is.na(picked)]
  picked
}

patch_scenic_annotation_lookup <- function() {
  scenic_ns <- asNamespace("SCENIC")
  original_fun <- get("getDbAnnotations", envir = scenic_ns)
  patched_fun <- original_fun

  body(patched_fun) <- quote({
    dbAnnotFiles <- scenicOptions@settings$db_annotFiles
    if (!is.null(dbAnnotFiles)) {
      motifAnnotations <- NULL
      for (annotPath in dbAnnotFiles) {
        motifAnnot <- data.table::fread(annotPath)
        motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
        colnames(motifAnnot)[1] <- "motif"
        levels(motifAnnot$annotationSource) <- c(
          levels(motifAnnot$annotationSource),
          c(
            "directAnnotation",
            "inferredBy_Orthology",
            "inferredBy_MotifSimilarity",
            "inferredBy_MotifSimilarity_n_Orthology"
          )
        )
        motifAnnotations <- rbind(motifAnnotations, motifAnnot)
      }
    } else {
      if (is.na(getDatasetInfo(scenicOptions, "org"))) {
        stop("Please provide an organism (scenicOptions@inputDatasetInfo$org).")
      }
      org <- getDatasetInfo(scenicOptions, "org")
      if (is.na(org)) {
        stop("Please provide an organism (scenicOptions@inputDatasetInfo$org).")
      }
      if (!org %in% c("hgnc", "mgi", "dmel")) {
        stop("Organism not recognized (scenicOptions@inputDatasetInfo$org).")
      }

      if (org == "hgnc") motifAnnotName <- "motifAnnotations_hgnc"
      if (org == "mgi") motifAnnotName <- "motifAnnotations_mgi"
      if (org == "dmel") motifAnnotName <- "motifAnnotations_dmel"

      if (!is.null(scenicOptions@settings$db_mcVersion)) {
        if (scenicOptions@settings$db_mcVersion == "v8") {
          motifAnnotName <- paste0(motifAnnotName, "_v8")
        }
      }

      annot_env <- new.env(parent = baseenv())
      data(list = motifAnnotName, package = "RcisTarget", envir = annot_env, verbose = FALSE)
      if (!exists(motifAnnotName, envir = annot_env, inherits = FALSE)) {
        v9_name <- paste0(motifAnnotName, "_v9")
        data(list = v9_name, package = "RcisTarget", envir = annot_env, verbose = FALSE)
        if (exists(v9_name, envir = annot_env, inherits = FALSE)) {
          assign(motifAnnotName, get(v9_name, envir = annot_env), envir = annot_env)
        }
      }
      motifAnnotations <- get(motifAnnotName, envir = annot_env, inherits = FALSE)
    }
    return(motifAnnotations)
  })

  unlockBinding("getDbAnnotations", scenic_ns)
  assign("getDbAnnotations", patched_fun, envir = scenic_ns)
  lockBinding("getDbAnnotations", scenic_ns)
  invisible(TRUE)
}

patch_scenic_gene_filtering <- function() {
  scenic_ns <- asNamespace("SCENIC")
  original_fun <- get("geneFiltering", envir = scenic_ns)
  patched_fun <- original_fun

  body(patched_fun) <- quote({
    outFile_genesKept <- NULL
    dbFilePath <- NULL

    if (class(scenicOptions) == "ScenicOptions") {
      dbFilePath <- getDatabases(scenicOptions)[[1]]
      outFile_genesKept <- getIntName(scenicOptions, "genesKept")
    } else {
      dbFilePath <- scenicOptions[["dbFilePath"]]
      outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
    }

    if (is.null(dbFilePath)) stop("dbFilePath")
    if (is.data.frame(exprMat)) {
      supportedClasses <- paste(
        gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")),
        collapse = ", "
      )
      supportedClasses <- gsub("-method", "", supportedClasses)
      stop(
        "'exprMat' should be one of the following classes: ",
        supportedClasses,
        "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)"
      )
    }
    if (any(table(rownames(exprMat)) > 1)) {
      stop("The rownames (gene id/name) in the expression matrix should be unique.")
    }

    if (inherits(exprMat, "Matrix") || inherits(exprMat, "sparseMatrix")) {
      nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = TRUE)
      nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = TRUE)
    } else {
      nCountsPerGene <- rowSums(exprMat, na.rm = TRUE)
      nCellsPerGene <- rowSums(exprMat > 0, na.rm = TRUE)
    }

    message("Maximum value in the expression matrix: ", max(exprMat, na.rm = TRUE))
    message(
      "Ratio of detected vs non-detected: ",
      signif(sum(exprMat > 0, na.rm = TRUE) / sum(exprMat == 0, na.rm = TRUE), 2)
    )
    message("Number of counts (in the dataset units) per gene:")
    print(summary(nCountsPerGene))
    message("Number of cells in which each gene is detected:")
    print(summary(nCellsPerGene))
    message("\nNumber of genes left after applying the following filters (sequential):")

    genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
    message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", minCountsPerGene)
    nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
    genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
    message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", minSamples, " cells")

    library(RcisTarget)
    motifRankings <- importRankings(dbFilePath)
    genesInDatabase <- colnames(getRanking(motifRankings))
    genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
    message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
    genesKept <- genesLeft_minCells_inDatabases

    if (!is.null(outFile_genesKept)) {
      saveRDS(genesKept, file = outFile_genesKept)
      if (getSettings(scenicOptions, "verbose")) {
        message("Gene list saved in ", outFile_genesKept)
      }
    }

    return(genesKept)
  })

  unlockBinding("geneFiltering", scenic_ns)
  assign("geneFiltering", patched_fun, envir = scenic_ns)
  lockBinding("geneFiltering", scenic_ns)
  invisible(TRUE)
}

scenic_gene_filtering_sparse <- function(
    exprMat,
    scenicOptions,
    minCountsPerGene,
    minSamples
) {
  outFile_genesKept <- NULL
  dbFilePath <- NULL

  if (class(scenicOptions) == "ScenicOptions") {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  } else {
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }

  if (is.null(dbFilePath)) stop("dbFilePath")
  if (is.data.frame(exprMat)) {
    supportedClasses <- paste(
      gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")),
      collapse = ", "
    )
    supportedClasses <- gsub("-method", "", supportedClasses)
    stop(
      "'exprMat' should be one of the following classes: ",
      supportedClasses,
      "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)"
    )
  }
  if (any(table(rownames(exprMat)) > 1)) {
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  }

  if (inherits(exprMat, "Matrix") || inherits(exprMat, "sparseMatrix")) {
    nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = TRUE)
  } else {
    nCountsPerGene <- rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- rowSums(exprMat > 0, na.rm = TRUE)
  }

  message("Maximum value in the expression matrix: ", max(exprMat, na.rm = TRUE))
  message(
    "Ratio of detected vs non-detected: ",
    signif(sum(exprMat > 0, na.rm = TRUE) / sum(exprMat == 0, na.rm = TRUE), 2)
  )
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  message("\nNumber of genes left after applying the following filters (sequential):")

  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", minCountsPerGene)
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", minSamples, " cells")

  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath)
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  genesKept <- genesLeft_minCells_inDatabases

  if (!is.null(outFile_genesKept)) {
    saveRDS(genesKept, file = outFile_genesKept)
    if (getSettings(scenicOptions, "verbose")) {
      message("Gene list saved in ", outFile_genesKept)
    }
  }

  genesKept
}



arg_list <- parse_args(commandArgs(trailingOnly = TRUE))
prepare_only <- to_flag(arg_list[["prepare_only"]], default = FALSE)
n_cores <- as.integer(arg_list[["n_cores"]] %||% "8")
cells_per_mp <- as.integer(arg_list[["cells_per_mp"]] %||% "350")
min_cells_per_mp <- as.integer(arg_list[["min_cells_per_mp"]] %||% "80")
min_best_z <- as.numeric(arg_list[["min_best_z"]] %||% "0.6")
min_gap_z <- as.numeric(arg_list[["min_gap_z"]] %||% "0.15")
top_genes_per_set <- as.integer(arg_list[["top_genes_per_set"]] %||% "100")
top_regulons_per_mp <- as.integer(arg_list[["top_regulons_per_mp"]] %||% "8")
db_dir <- arg_list[["db_dir"]] %||% Sys.getenv("SCENIC_DB_DIR", unset = "")
if (!nzchar(db_dir)) {
  db_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr"
}
db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)

out_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/final_mp_scenic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ephemeral_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/final_mp_scenic"
dir.create(ephemeral_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(ephemeral_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Final MP definitions
####################
sc_mps <- c(
  "MP1", "MP5", "MP13+",
  "MP2+",
  "MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+",
  "MP8+", "MP8b", "MP16", "MP18b", "MP17",
  "MP12",
  "MP15"
)

sc_mp_descriptions <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)

mp_group_map <- c(
  "MP1" = "Cell cycle",
  "MP5" = "Cell cycle",
  "MP13+" = "Cell cycle",
  "MP2+" = "Classic proliferation",
  "MP14" = "Squamous-to-intestinal",
  "MP3+" = "Squamous-to-intestinal",
  "MP6+" = "Squamous-to-intestinal",
  "MP11+" = "Squamous-to-intestinal",
  "MP9+" = "Squamous-to-intestinal",
  "MP10+" = "Squamous-to-intestinal",
  "MP8+" = "Glandular-to-intestinal",
  "MP8b" = "Glandular-to-intestinal",
  "MP16" = "Glandular-to-intestinal",
  "MP18b" = "Glandular-to-intestinal",
  "MP17" = "Glandular-to-intestinal",
  "MP12" = "Stress-adaptive",
  "MP15" = "Cancer-cell immune mimicry"
)

group_cols <- c(
  "Cell cycle" = "#6B7280",
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

state_level_order <- c(
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)

state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

####################
# Load inputs
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds")
mp_gene_sets <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")

final_mp_ids <- sc_mps
display_label_map <- setNames(paste(sc_mps, sc_mp_descriptions[sc_mps]), sc_mps)

####################
# Score assembly and cell selection
####################
selection_cache <- file.path(cache_dir, "mp_cell_selection_final17.rds")

if (file.exists(selection_cache)) {
  message("Loading cached MP cell selection.")
  cached_sel <- readRDS(selection_cache)
  common_cells <- cached_sel$common_cells
  score_mat <- cached_sel$score_mat
  z_mat <- cached_sel$z_mat
  assignment_df <- cached_sel$assignment_df
  eligible_df <- cached_sel$eligible_df
  selected_df <- cached_sel$selected_df
} else {
  common_cells <- Reduce(
    intersect,
    list(
      Cells(tmdata_all),
      names(final_states),
      rownames(ucell_scores)
    )
  )
  if (length(common_cells) == 0) {
    stop("No overlapping cells across required inputs.")
  }

  tmdata_all <- tmdata_all[, common_cells]
  final_states <- final_states[common_cells]
  ucell_scores <- ucell_scores[common_cells, , drop = FALSE]

  missing_mps <- setdiff(final_mp_ids, colnames(ucell_scores))
  if (length(missing_mps) > 0) {
    stop(
      "Selected final MPs are missing from merged_refined_ucell_scores.rds: ",
      paste(missing_mps, collapse = ", ")
    )
  }

  score_mat <- as.matrix(ucell_scores[, final_mp_ids, drop = FALSE])

  sample_var <- tmdata_all$orig.ident
  study_var <- tmdata_all$study
  names(sample_var) <- Cells(tmdata_all)
  names(study_var) <- Cells(tmdata_all)
  z_mat <- z_normalise(score_mat, sample_var, study_var)

  best_idx <- max.col(z_mat, ties.method = "first")
  best_id <- colnames(z_mat)[best_idx]
  best_z <- z_mat[cbind(seq_len(nrow(z_mat)), best_idx)]
  sorted_scores <- t(apply(z_mat, 1, sort, decreasing = TRUE))
  gap_z <- sorted_scores[, 1] - sorted_scores[, 2]

  assignment_df <- data.frame(
    cell = rownames(z_mat),
    final_mp_id = best_id,
    final_mp_label = display_label_map[best_id],
    mp_group = mp_group_map[best_id],
    best_z = as.numeric(best_z),
    gap_z = as.numeric(gap_z),
    orig.ident = as.character(tmdata_all$orig.ident[rownames(z_mat)]),
    study = as.character(tmdata_all$study[rownames(z_mat)]),
    final_state = as.character(final_states[rownames(z_mat)]),
    stringsAsFactors = FALSE
  )

  eligible_df <- assignment_df %>%
    filter(best_z >= min_best_z, gap_z >= min_gap_z)

  eligible_counts <- eligible_df %>%
    count(final_mp_id, final_mp_label, mp_group, sort = TRUE, name = "n_eligible")

  eligible_mps <- eligible_counts %>%
    filter(n_eligible >= min_cells_per_mp) %>%
    pull(final_mp_id)

  if (length(eligible_mps) < 3) {
    stop(
      "Fewer than 3 MPs passed the selection filters. ",
      "Try lower thresholds or run with prepare_only=true to inspect the summary."
    )
  }

  selected_df <- eligible_df %>%
    filter(final_mp_id %in% eligible_mps) %>%
    group_by(final_mp_id) %>%
    arrange(desc(best_z), desc(gap_z), .by_group = TRUE) %>%
    slice_head(n = cells_per_mp) %>%
    ungroup()
  
  message("Caching MP cell selection.")
  saveRDS(list(
    common_cells = common_cells,
    score_mat = score_mat,
    z_mat = z_mat,
    assignment_df = assignment_df,
    eligible_df = eligible_df,
    selected_df = selected_df
  ), selection_cache)
}

selected_df$final_mp_label <- factor(
  selected_df$final_mp_label,
  levels = display_label_map[final_mp_ids[final_mp_ids %in% selected_df$final_mp_id]]
)
selected_df$mp_group <- factor(selected_df$mp_group, levels = names(group_cols))

selected_counts <- selected_df %>%
  count(final_mp_id, final_mp_label, mp_group, sort = FALSE, name = "n_selected")

write.csv(
  selected_df,
  file.path(out_dir, "Auto_final_mp_scenic_selected_cells.csv"),
  row.names = FALSE
)

assignment_summary <- assignment_df %>%
  group_by(final_mp_id, final_mp_label, mp_group) %>%
  summarise(
    n_cells = n(),
    n_above_threshold = sum(best_z >= min_best_z & gap_z >= min_gap_z),
    mean_best_z = mean(best_z, na.rm = TRUE),
    mean_gap_z = mean(gap_z, na.rm = TRUE),
    n_samples = n_distinct(orig.ident),
    n_studies = n_distinct(study),
    .groups = "drop"
  ) %>%
  arrange(match(final_mp_id, final_mp_ids))

write.csv(
  assignment_summary,
  file.path(out_dir, "Auto_final_mp_scenic_assignment_summary.csv"),
  row.names = FALSE
)

p_selected <- ggplot(selected_counts, aes(x = final_mp_label, y = n_selected, fill = mp_group)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(
    title = "Final MP-selected cells for SCENIC",
    subtitle = paste0(
      "best_z >= ", min_best_z,
      ", gap_z >= ", min_gap_z,
      ", up to ", cells_per_mp, " cells per MP"
    ),
    x = NULL,
    y = "Selected cells",
    fill = "MP group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "Auto_final_mp_scenic_selected_cells.pdf"),
  p_selected,
  width = 12,
  height = 8
)

####################
# Final MP gene sets
####################
sc_mp_gene_sets <- mp_gene_sets[sc_mps]
sc_mp_gene_sets <- lapply(sc_mp_gene_sets, function(x) unique(x[!is.na(x) & nzchar(x)]))
sc_mp_gene_sets <- lapply(sc_mp_gene_sets, function(x) head(x, top_genes_per_set))

final_gene_sets <- sc_mp_gene_sets
names(final_gene_sets) <- display_label_map[names(final_gene_sets)]

saveRDS(
  final_gene_sets,
  file.path(out_dir, "Auto_final_mp_scenic_gene_sets.rds")
)

gene_membership_df <- bind_rows(lapply(names(final_gene_sets), function(mp_label) {
  genes <- final_gene_sets[[mp_label]]
  data.frame(
    final_mp_label = mp_label,
    gene = genes,
    rank = seq_along(genes),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  gene_membership_df,
  file.path(out_dir, "Auto_final_mp_scenic_gene_membership.csv"),
  row.names = FALSE
)

if (prepare_only) {
  summary_dir <- file.path(
    "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
    "updates", "new_updates", "summaries"
  )
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- data.frame(
    mode = "prepare_only",
    n_common_cells = length(common_cells),
    n_selected_cells = nrow(selected_df),
    n_selected_mps = length(unique(selected_df$final_mp_id)),
    db_dir = db_dir,
    stringsAsFactors = FALSE
  )
  write.csv(
    summary_df,
    file.path(summary_dir, "Auto_final_mp_scenic_summary.csv"),
    row.names = FALSE
  )
  message("Prepare-only mode complete. Selection and gene-set outputs saved in ", out_dir)
  quit(save = "no")
}

####################
# SCENIC dependency checks
####################
required_pkgs <- c("SCENIC", "AUCell", "RcisTarget", "GENIE3", "doRNG", "doMC")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required SCENIC packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Run this script with prepare_only=true to validate selection without SCENIC."
  )
}

suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(GENIE3)
})

patch_scenic_annotation_lookup()
patch_scenic_gene_filtering()
db_files <- detect_db_files(db_dir)

####################
# Expression matrix for SCENIC
####################
counts_mat <- get_assay_matrix(tmdata_all, "counts")
counts_mat <- counts_mat[, selected_df$cell, drop = FALSE]
if (!inherits(counts_mat, "dgCMatrix")) {
  counts_mat <- as(counts_mat, "dgCMatrix")
}

expr_genes <- unique(unlist(final_gene_sets, use.names = FALSE))
expr_genes <- intersect(expr_genes, rownames(counts_mat))

old_wd <- getwd()
setwd(ephemeral_dir)
on.exit(setwd(old_wd), add = TRUE)

scenicOptions <- initializeScenic(
  org = "hgnc",
  dbDir = db_dir,
  dbs = db_files,
  datasetTitle = "Auto_final_mp_scenic",
  nCores = n_cores
)

min_counts_per_gene <- max(3 * 0.01 * ncol(counts_mat), 20)
min_samples <- max(0.01 * ncol(counts_mat), 20)
genes_kept_path <- file.path("int", "1.1_genesKept.Rds")
####################
# Reuse saved SCENIC intermediates when available so reruns restart from the
# last missing stage instead of repeating the full network build.
####################
if (file.exists(genes_kept_path)) {
  message("Reusing existing genesKept intermediate: ", genes_kept_path)
  genes_kept <- readRDS(genes_kept_path)
} else {
  genes_kept <- scenic_gene_filtering_sparse(
    counts_mat,
    scenicOptions = scenicOptions,
    minCountsPerGene = min_counts_per_gene,
    minSamples = min_samples
  )
}
expr_mat_filtered <- counts_mat[genes_kept, , drop = FALSE]

db_tfs <- tryCatch(getDbTfs(scenicOptions), error = function(e) character(0))
focus_genes <- unique(c(
  intersect(expr_genes, rownames(expr_mat_filtered)),
  intersect(db_tfs, rownames(expr_mat_filtered))
))

if (length(focus_genes) >= 500) {
  expr_mat_use <- expr_mat_filtered[focus_genes, , drop = FALSE]
} else {
  expr_mat_use <- expr_mat_filtered
}
if (!is.matrix(expr_mat_use)) {
  expr_mat_use <- as.matrix(expr_mat_use)
}

corr_mat_path <- file.path("int", "1.2_corrMat.Rds")
if (file.exists(corr_mat_path)) {
  message("Reusing existing correlation matrix: ", corr_mat_path)
} else {
  runCorrelation(expr_mat_use, scenicOptions)
}

genie3_linklist_path <- file.path("int", "1.4_GENIE3_linkList.Rds")
if (file.exists(genie3_linklist_path)) {
  message("Reusing existing GENIE3 network: ", genie3_linklist_path)
} else {
  runGenie3(expr_mat_use, scenicOptions)
}

tf_modules_path <- file.path("int", "1.6_tfModules_asDF.Rds")
if (file.exists(tf_modules_path)) {
  message("Reusing existing TF modules: ", tf_modules_path)
} else {
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
}

regulons_path <- file.path("int", "2.6_regulons_asGeneSet.Rds")
regulons_mat_path <- file.path("int", "2.6_regulons_asIncidMat.Rds")
if (file.exists(regulons_path) && file.exists(regulons_mat_path)) {
  message("Reusing existing regulons: ", regulons_path)
} else {
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
}

regulon_auc_path <- file.path("int", "3.4_regulonAUC.Rds")
if (file.exists(regulon_auc_path)) {
  message("Reusing existing regulon AUC: ", regulon_auc_path)
} else {
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat = counts_mat)
}

regulon_auc <- loadInt(scenicOptions, "aucell_regulonAUC")
regulons <- loadInt(scenicOptions, "regulons")
auc_mat <- getAUC(regulon_auc)

cell_label_map <- setNames(as.character(selected_df$final_mp_label), selected_df$cell)
cell_label_map <- cell_label_map[colnames(auc_mat)]

mean_auc_mat <- sapply(split(names(cell_label_map), cell_label_map), function(cells) {
  rowMeans(auc_mat[, cells, drop = FALSE], na.rm = TRUE)
})
mean_auc_mat <- as.matrix(mean_auc_mat)

rss_mat <- tryCatch(
  calcRSS(AUC = auc_mat, cellAnnotation = cell_label_map),
  error = function(e) NULL
)
if (is.null(rss_mat)) {
  rss_mat <- mean_auc_mat
}
rss_mat <- as.matrix(rss_mat)

setwd(out_dir)

saveRDS(regulon_auc, "Auto_final_mp_scenic_regulon_auc.rds")
saveRDS(rss_mat, "Auto_final_mp_scenic_rss.rds")

####################
# Regulon heatmap
####################
mp_label_order <- levels(selected_df$final_mp_label)
mean_auc_mat <- mean_auc_mat[, mp_label_order, drop = FALSE]
rss_mat <- rss_mat[, mp_label_order, drop = FALSE]

top_regulons <- unique(unlist(lapply(mp_label_order, function(mp_label) {
  vals <- sort(rss_mat[, mp_label], decreasing = TRUE)
  names(vals)[seq_len(min(top_regulons_per_mp, length(vals)))]
})))
top_regulons <- top_regulons[!is.na(top_regulons)]

plot_rss_mat <- rss_mat[top_regulons, mp_label_order, drop = FALSE]
plot_rss_scaled <- t(scale(t(plot_rss_mat)))
plot_rss_scaled[!is.finite(plot_rss_scaled)] <- 0

column_groups <- mp_group_map[final_mp_ids]
column_groups <- setNames(column_groups, display_label_map[final_mp_ids])
column_groups <- factor(column_groups[mp_label_order], levels = names(group_cols))

ha_cols <- HeatmapAnnotation(
  MP_group = column_groups,
  col = list(MP_group = group_cols),
  show_annotation_name = FALSE
)

rownames(plot_rss_scaled) <- format_regulon_name(rownames(plot_rss_scaled))
rss_col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf("Auto_final_mp_scenic_regulon_heatmap.pdf", width = 16, height = 10, useDingbats = FALSE)
draw(
  Heatmap(
    plot_rss_scaled,
    name = "Scaled RSS",
    col = rss_col_fun,
    top_annotation = ha_cols,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    heatmap_legend_param = list(title = "Scaled RSS")
  ),
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grid.text(
  "SCENIC regulon specificity across final MPs",
  x = unit(4, "mm"),
  y = unit(1, "npc") - unit(4, "mm"),
  just = c("left", "top"),
  gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

####################
# MP clustering heatmap (Top 100 regulons)
####################
# Select top 100 regulons globally based on max RSS across MPs
global_rss_top100 <- names(sort(apply(rss_mat[, mp_label_order, drop = FALSE], 1, max), decreasing = TRUE))
global_rss_top100 <- global_rss_top100[1:max(100, length(global_rss_top100))]
global_rss_top100 <- global_rss_top100[!is.na(global_rss_top100)]

plot_rss_clust_mat <- rss_mat[global_rss_top100, mp_label_order, drop = FALSE]
plot_rss_clust_scaled <- t(scale(t(plot_rss_clust_mat)))
plot_rss_clust_scaled[!is.finite(plot_rss_clust_scaled)] <- 0
rownames(plot_rss_clust_scaled) <- format_regulon_name(rownames(plot_rss_clust_scaled))

pdf("Auto_final_mp_scenic_mp_clustering_heatmap.pdf", width = 16, height = 12, useDingbats = FALSE)
draw(
  Heatmap(
    plot_rss_clust_scaled,
    name = "Scaled RSS",
    col = rss_col_fun,
    top_annotation = ha_cols,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    heatmap_legend_param = list(title = "Scaled RSS")
  ),
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grid.text(
  "MP clustering by top 100 SCENIC regulon specificity",
  x = unit(4, "mm"),
  y = unit(1, "npc") - unit(4, "mm"),
  just = c("left", "top"),
  gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

####################
# Regulon network
####################
edge_df <- bind_rows(lapply(mp_label_order, function(mp_label) {
  vals <- sort(rss_mat[, mp_label], decreasing = TRUE)
  keep <- names(vals)[seq_len(min(top_regulons_per_mp, length(vals)))]
  data.frame(
    regulon = keep,
    regulon_label = format_regulon_name(keep),
    mp_label = mp_label,
    weight = as.numeric(vals[keep]),
    stringsAsFactors = FALSE
  )
})) %>%
  distinct(regulon_label, mp_label, .keep_all = TRUE) %>%
  filter(is.finite(weight), weight > 0)

write.csv(edge_df, "Auto_final_mp_scenic_network_edges.csv", row.names = FALSE)

node_df <- data.frame(
  name = unique(c(edge_df$regulon_label, edge_df$mp_label)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    node_type = ifelse(name %in% mp_label_order, "MP", "Regulon"),
    mp_group = ifelse(node_type == "MP", as.character(column_groups[name]), "Regulon")
  )

network_graph <- tbl_graph(
  nodes = node_df,
  edges = edge_df %>% transmute(from = regulon_label, to = mp_label, weight = weight),
  directed = FALSE
)

network_fill <- c(group_cols, Regulon = "grey35")

pdf("Auto_final_mp_scenic_network.pdf", width = 16, height = 10, useDingbats = FALSE)
print(
  ggraph(network_graph, layout = "stress") +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "grey70") +
    scale_edge_width(range = c(0.4, 2.2)) +
    scale_edge_alpha(range = c(0.3, 0.9)) +
    geom_node_point(
      aes(fill = mp_group, shape = node_type),
      size = 5,
      colour = "black",
      stroke = 0.3
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3
    ) +
    scale_shape_manual(values = c(MP = 21, Regulon = 22)) +
    scale_fill_manual(values = network_fill, drop = FALSE) +
    theme_void(base_size = 12) +
    labs(
      title = "SCENIC final-MP regulatory network"
    ) +
    guides(edge_width = "none", edge_alpha = "none", shape = "none", fill = "none")
)
dev.off()

####################
# Top-100 mean-AUC MP network (shared regulons across MPs)
####################
n_top_global <- 100
min_per_category <- 5

# Guarantee minimum 5 per MP, fill remaining from global top-100
mp_per_cat_regs <- lapply(mp_label_order, function(mp_label) {
  vals <- sort(mean_auc_mat[, mp_label], decreasing = TRUE)
  names(vals)[seq_len(min(min_per_category, length(vals)))]
})
guaranteed_regs <- unique(unlist(mp_per_cat_regs))
guaranteed_regs <- guaranteed_regs[!is.na(guaranteed_regs)]

global_top <- names(sort(apply(mean_auc_mat[, mp_label_order, drop = FALSE], 1, max), decreasing = TRUE))
remaining <- setdiff(global_top, guaranteed_regs)
n_fill <- max(0, n_top_global - length(guaranteed_regs))
auc_top_regulons <- unique(c(guaranteed_regs, head(remaining, n_fill)))

# Build edge list: connect each selected regulon to every MP where it ranks
# in that MP's per-regulon top list (i.e. above the median of top regulons)
auc_edge_df <- bind_rows(lapply(mp_label_order, function(mp_label) {
  vals <- mean_auc_mat[auc_top_regulons, mp_label]
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(NULL)
  threshold <- quantile(vals, 0.5, na.rm = TRUE)
  keep <- names(vals)[vals >= threshold]
  if (length(keep) == 0) return(NULL)
  data.frame(
    regulon_label = format_regulon_name(keep),
    mp_label = mp_label,
    weight = as.numeric(vals[keep]),
    stringsAsFactors = FALSE
  )
})) %>% distinct(regulon_label, mp_label, .keep_all = TRUE)

auc_node_df <- data.frame(
  name = unique(c(auc_edge_df$regulon_label, auc_edge_df$mp_label)),
  stringsAsFactors = FALSE
) %>% mutate(
  node_type = ifelse(name %in% mp_label_order, "MP", "Regulon"),
  mp_group = ifelse(node_type == "MP", as.character(column_groups[name]), "Regulon")
)

# Flag regulons connected to multiple MPs
reg_degree <- auc_edge_df %>% count(regulon_label, name = "n_mps")
auc_node_df <- auc_node_df %>%
  left_join(reg_degree, by = c("name" = "regulon_label")) %>%
  mutate(is_shared = !is.na(n_mps) & n_mps > 1)

auc_network_graph <- tbl_graph(
  nodes = auc_node_df,
  edges = auc_edge_df %>% transmute(from = regulon_label, to = mp_label, weight = weight),
  directed = FALSE
)

auc_network_fill <- c(group_cols, Regulon = "grey35")
pdf("Auto_final_mp_scenic_network_top100auc.pdf", width = 20, height = 14, useDingbats = FALSE)
print(
  ggraph(auc_network_graph, layout = "stress") +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "grey70") +
    scale_edge_width(range = c(0.3, 2.0)) +
    scale_edge_alpha(range = c(0.2, 0.8)) +
    geom_node_point(
      aes(fill = mp_group, shape = node_type, size = ifelse(node_type == "MP", 7, ifelse(is_shared, 5, 3.5))),
      colour = "black", stroke = 0.3
    ) +
    scale_size_identity() +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 30) +
    scale_shape_manual(values = c(MP = 21, Regulon = 22)) +
    scale_fill_manual(values = auc_network_fill, drop = FALSE) +
    theme_void(base_size = 12) +
    labs(title = "Top regulons by mean AUC across MPs") +
    guides(edge_width = "none", edge_alpha = "none", shape = "none", fill = "none")
)
dev.off()

####################
# Regulon target summary
####################
mp_gene_sets_by_id <- sc_mp_gene_sets
regulon_target_df <- bind_rows(lapply(seq_len(nrow(edge_df)), function(i) {
  reg_name <- edge_df$regulon[i]
  reg_targets <- extract_regulon_targets(regulons[[reg_name]])
  mp_label <- edge_df$mp_label[i]
  mp_id <- names(display_label_map)[match(mp_label, display_label_map)]
  overlap_n <- sum(reg_targets %in% mp_gene_sets_by_id[[mp_id]])
  data.frame(
    mp_label = mp_label,
    regulon = reg_name,
    regulon_label = edge_df$regulon_label[i],
    rss_weight = edge_df$weight[i],
    n_targets = length(reg_targets),
    overlap_with_mp_genes = overlap_n,
    targets_preview = paste(head(reg_targets, 30), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  regulon_target_df,
  "Auto_final_mp_scenic_regulon_targets.csv",
  row.names = FALSE
)

####################
# Excel summary of scATLAS SCENIC regulons by final MP
####################
regulon_target_df <- regulon_target_df %>%
  rowwise() %>%
  mutate(
    .targets = list(extract_regulon_targets(regulons[[format_regulon_name(regulon)]])),
    n_targets = length(.targets[[1]]),
    overlap_with_mp_genes = sum(.targets[[1]] %in% mp_gene_sets_by_id[[names(display_label_map)[match(mp_label, display_label_map)]]]),
    targets_preview = paste(head(.targets[[1]], 30), collapse = "; ")
  ) %>%
  ungroup() %>%
  select(-.targets)
write.csv(regulon_target_df, "Auto_final_mp_scenic_regulon_targets.csv", row.names = FALSE)

write_scenic_regulon_workbook <- function(
    specificity_mat,
    activity_mat,
    label_order,
    output_file,
    analysis_level,
    mp_gene_sets_by_label = NULL
) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("The openxlsx package is required to write SCENIC regulon Excel summaries.")
  }

  if (identical(analysis_level, "final MP")) {
    label_order <- intersect(display_label_map[sc_mps], label_order)
  }

  specificity_mat <- as.matrix(specificity_mat[, label_order, drop = FALSE])
  activity_mat <- as.matrix(activity_mat[, label_order, drop = FALSE])
  regulon_ids <- intersect(rownames(specificity_mat), rownames(activity_mat))
  if (length(regulon_ids) == 0) {
    stop("No overlapping regulons are available for the ", analysis_level, " Excel summary.")
  }
  specificity_mat <- specificity_mat[regulon_ids, , drop = FALSE]
  activity_mat <- activity_mat[regulon_ids, , drop = FALSE]

  get_regulon_targets_by_id <- function(regulon_id) {
    regulon_key <- format_regulon_name(regulon_id)
    extract_regulon_targets(regulons[[regulon_key]])
  }
  target_df <- bind_rows(lapply(regulon_ids, function(regulon_id) {
    targets <- get_regulon_targets_by_id(regulon_id)
    data.frame(
      regulon_id = regulon_id,
      n_targets = length(targets),
      targets_preview = paste(head(targets, 100), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))

  overview_df <- target_df %>% select(regulon_id, n_targets, targets_preview)
  for (label in label_order) {
    overview_df[[label]] <- specificity_mat[regulon_ids, label]
  }

  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#1F4E78",
    fontColour = "#FFFFFF",
    halign = "center",
    valign = "center",
    wrapText = TRUE
  )
  title_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fontSize = 14,
    fontColour = "#1F1F1F"
  )
  numeric_style <- openxlsx::createStyle(numFmt = "0.0000")
  text_style <- openxlsx::createStyle(valign = "top", wrapText = TRUE)

  overview_sheet <- "All regulons"
  openxlsx::addWorksheet(wb, overview_sheet)
  openxlsx::writeData(
    wb,
    overview_sheet,
    paste0("scATLAS SCENIC regulons by ", analysis_level),
    startRow = 1,
    startCol = 1
  )
  openxlsx::addStyle(wb, overview_sheet, title_style, rows = 1, cols = 1)
  openxlsx::writeData(wb, overview_sheet, overview_df, startRow = 3, headerStyle = header_style)
  openxlsx::addStyle(
    wb, overview_sheet, numeric_style,
    rows = 4:(nrow(overview_df) + 3),
    cols = 4:ncol(overview_df), gridExpand = TRUE, stack = TRUE
  )
  openxlsx::addStyle(
    wb, overview_sheet, text_style,
    rows = 4:(nrow(overview_df) + 3), cols = 1:3, gridExpand = TRUE, stack = TRUE
  )
  openxlsx::setColWidths(wb, overview_sheet, cols = 1, widths = 28)
  openxlsx::setColWidths(wb, overview_sheet, cols = 2, widths = 12)
  openxlsx::setColWidths(wb, overview_sheet, cols = 3, widths = 60)
  openxlsx::setColWidths(wb, overview_sheet, cols = 4:ncol(overview_df), widths = 14)
  openxlsx::freezePane(wb, overview_sheet, firstActiveRow = 4, firstActiveCol = 3)

  for (cols in list(4:ncol(overview_df))) {
    values <- as.numeric(as.matrix(overview_df[, cols, drop = FALSE]))
    values <- values[is.finite(values)]
    if (length(values) > 0) {
      limits <- as.numeric(stats::quantile(values, c(0.05, 0.5, 0.95), na.rm = TRUE))
      if (length(unique(limits)) == 3) {
        openxlsx::conditionalFormatting(
          wb, overview_sheet, cols = cols, rows = 4:(nrow(overview_df) + 3),
          type = "colourScale", style = c("#FFFFFF", "#F4B183", "#C00000"), rule = limits
        )
      }
    }
  }

  for (index in seq_along(label_order)) {
    label <- label_order[index]
    sheet_name <- substr(paste0(sprintf("%02d", index), " ", label), 1, 31)
    sheet_name <- gsub("/", "-", sheet_name, fixed = TRUE)
    sheet_df <- data.frame(
      regulon_id = regulon_ids,
      rss = specificity_mat[regulon_ids, label],
      rss_rank = rank(-specificity_mat[regulon_ids, label], ties.method = "min", na.last = "keep"),
      mean_auc = activity_mat[regulon_ids, label],
      n_targets = target_df$n_targets[match(regulon_ids, target_df$regulon_id)],
      targets_preview = target_df$targets_preview[match(regulon_ids, target_df$regulon_id)],
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(rss), desc(mean_auc), regulon_id)

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, paste0(label, " — ", analysis_level, " regulons"), startRow = 1)
    openxlsx::addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
    openxlsx::writeData(wb, sheet_name, sheet_df, startRow = 3, headerStyle = header_style)
    openxlsx::addStyle(
      wb, sheet_name, numeric_style, rows = 4:(nrow(sheet_df) + 3),
      cols = c(2, 4), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      wb, sheet_name, text_style, rows = 4:(nrow(sheet_df) + 3),
      cols = c(1, 5:6), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::conditionalFormatting(
      wb, sheet_name, cols = 2, rows = 4:(nrow(sheet_df) + 3),
      type = "colourScale", style = c("#FFFFFF", "#F4B183", "#C00000")
    )
    openxlsx::conditionalFormatting(
      wb, sheet_name, cols = 4, rows = 4:(nrow(sheet_df) + 3),
      type = "colourScale", style = c("#FFFFFF", "#9DC3E6", "#2F5597")
    )
    openxlsx::setColWidths(wb, sheet_name, cols = 1, widths = 28)
    openxlsx::setColWidths(wb, sheet_name, cols = 2:5, widths = 13)
    openxlsx::setColWidths(wb, sheet_name, cols = 6, widths = 60)
    openxlsx::freezePane(wb, sheet_name, firstActiveRow = 4, firstActiveCol = 3)
  }

  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  message("Saved ", analysis_level, " SCENIC regulon Excel summary: ", output_file)
}

mp_gene_sets_by_label <- setNames(sc_mp_gene_sets, display_label_map[names(sc_mp_gene_sets)])
write_scenic_regulon_workbook(
  specificity_mat = rss_mat,
  activity_mat = mean_auc_mat,
  label_order = intersect(display_label_map[sc_mps], mp_label_order),
  output_file = "Auto_final_mp_scenic_regulons_by_mp.xlsx",
  analysis_level = "final MP",
  mp_gene_sets_by_label = mp_gene_sets_by_label
)
####################

####################
# State-level regulon summaries
####################
state_df <- selected_df %>%
  filter(final_state %in% state_level_order) %>%
  mutate(final_state = factor(final_state, levels = state_level_order))

if (nrow(state_df) > 0 && dplyr::n_distinct(state_df$final_state) >= 2) {
  state_cells <- as.character(state_df$cell)
  state_label_map <- setNames(as.character(state_df$final_state), state_cells)
  state_auc_mat <- auc_mat[, state_cells, drop = FALSE]
  state_label_map <- state_label_map[colnames(state_auc_mat)]
  state_label_order <- levels(droplevels(state_df$final_state))
  state_label_order <- state_label_order[state_label_order %in% unique(as.character(state_df$final_state))]

  state_mean_auc_mat <- sapply(split(names(state_label_map), state_label_map), function(cells) {
    rowMeans(state_auc_mat[, cells, drop = FALSE], na.rm = TRUE)
  })
  state_mean_auc_mat <- as.matrix(state_mean_auc_mat)
  state_mean_auc_mat <- state_mean_auc_mat[, state_label_order, drop = FALSE]

  state_rss_mat <- tryCatch(
    calcRSS(AUC = state_auc_mat, cellAnnotation = state_label_map),
    error = function(e) NULL
  )
  if (is.null(state_rss_mat)) {
    state_rss_mat <- state_mean_auc_mat
  }
  state_rss_mat <- as.matrix(state_rss_mat)
  state_rss_mat <- state_rss_mat[, state_label_order, drop = FALSE]

  saveRDS(state_rss_mat, "Auto_final_mp_scenic_state_rss.rds")

  state_top_regulons <- unique(unlist(lapply(state_label_order, function(state_label) {
    vals <- sort(state_rss_mat[, state_label], decreasing = TRUE)
    names(vals)[seq_len(min(top_regulons_per_mp, length(vals)))]
  })))
  state_top_regulons <- state_top_regulons[!is.na(state_top_regulons)]

  state_plot_rss_mat <- state_rss_mat[state_top_regulons, state_label_order, drop = FALSE]
  state_plot_rss_scaled <- t(scale(t(state_plot_rss_mat)))
  state_plot_rss_scaled[!is.finite(state_plot_rss_scaled)] <- 0
  rownames(state_plot_rss_scaled) <- format_regulon_name(rownames(state_plot_rss_scaled))

  state_annotation <- HeatmapAnnotation(
    State = factor(state_label_order, levels = state_level_order),
    col = list(State = state_cols),
    show_annotation_name = FALSE
  )

  pdf(
    "Auto_final_mp_scenic_state_regulon_heatmap.pdf",
    width = 13,
    height = 10,
    useDingbats = FALSE
  )
  draw(
    Heatmap(
      state_plot_rss_scaled,
      name = "Scaled RSS",
      col = rss_col_fun,
      top_annotation = state_annotation,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled RSS")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(
    "SCENIC regulon specificity across final states",
    x = unit(4, "mm"),
    y = unit(1, "npc") - unit(4, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()

  state_edge_df <- bind_rows(lapply(state_label_order, function(state_label) {
    vals <- sort(state_rss_mat[, state_label], decreasing = TRUE)
    keep <- names(vals)[seq_len(min(top_regulons_per_mp, length(vals)))]
    data.frame(
      regulon = keep,
      regulon_label = format_regulon_name(keep),
      state_label = state_label,
      weight = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) %>%
    distinct(regulon_label, state_label, .keep_all = TRUE) %>%
    filter(is.finite(weight), weight > 0)

  write.csv(
    state_edge_df,
    "Auto_final_mp_scenic_state_network_edges.csv",
    row.names = FALSE
  )

  ####################
  # Excel summary of scATLAS SCENIC regulons by final state
  ####################
  write_scenic_regulon_workbook(
    specificity_mat = state_rss_mat,
    activity_mat = state_mean_auc_mat,
    label_order = state_label_order,
    output_file = "Auto_final_mp_scenic_regulons_by_state.xlsx",
    analysis_level = "final state"
  )
  ####################

  state_node_df <- data.frame(
    name = unique(c(state_edge_df$regulon_label, state_edge_df$state_label)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      node_type = ifelse(name %in% state_label_order, "State", "Regulon"),
      state_group = ifelse(node_type == "State", name, "Regulon")
    )

  state_network_graph <- tbl_graph(
    nodes = state_node_df,
    edges = state_edge_df %>% transmute(from = regulon_label, to = state_label, weight = weight),
    directed = FALSE
  )

  state_network_fill <- c(state_cols, Regulon = "grey35")

  pdf("Auto_final_mp_scenic_state_network.pdf", width = 15, height = 10, useDingbats = FALSE)
  print(
    ggraph(state_network_graph, layout = "stress") +
      geom_edge_link(aes(width = weight, alpha = weight), colour = "grey70") +
      scale_edge_width(range = c(0.4, 2.2)) +
      scale_edge_alpha(range = c(0.3, 0.9)) +
      geom_node_point(
        aes(fill = state_group, shape = node_type),
        size = 5,
        colour = "black",
        stroke = 0.3
      ) +
      geom_node_text(
        aes(label = name),
        repel = TRUE,
        size = 3
      ) +
      scale_shape_manual(values = c(State = 21, Regulon = 22)) +
      scale_fill_manual(values = state_network_fill, drop = FALSE) +
      theme_void(base_size = 12) +
      labs(
        title = "SCENIC final-state regulatory network"
      ) +
      guides(edge_width = "none", edge_alpha = "none", shape = "none", fill = "none")
  )
  dev.off()

  ####################
  # Top-100 mean-AUC state network (shared regulons across states)
  ####################
  state_per_cat_regs <- lapply(state_label_order, function(st) {
    vals <- sort(state_mean_auc_mat[, st], decreasing = TRUE)
    names(vals)[seq_len(min(min_per_category, length(vals)))]
  })
  state_guaranteed_regs <- unique(unlist(state_per_cat_regs))
  state_guaranteed_regs <- state_guaranteed_regs[!is.na(state_guaranteed_regs)]

  state_global_top <- names(sort(apply(state_mean_auc_mat[, state_label_order, drop = FALSE], 1, max), decreasing = TRUE))
  state_remaining <- setdiff(state_global_top, state_guaranteed_regs)
  state_n_fill <- max(0, n_top_global - length(state_guaranteed_regs))
  state_auc_top_regulons <- unique(c(state_guaranteed_regs, head(state_remaining, state_n_fill)))

  state_auc_edge_df <- bind_rows(lapply(state_label_order, function(st) {
    vals <- state_mean_auc_mat[state_auc_top_regulons, st]
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) == 0) return(NULL)
    threshold <- quantile(vals, 0.5, na.rm = TRUE)
    keep <- names(vals)[vals >= threshold]
    if (length(keep) == 0) return(NULL)
    data.frame(
      regulon_label = format_regulon_name(keep),
      state_label = st,
      weight = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) %>% distinct(regulon_label, state_label, .keep_all = TRUE)

  state_auc_node_df <- data.frame(
    name = unique(c(state_auc_edge_df$regulon_label, state_auc_edge_df$state_label)),
    stringsAsFactors = FALSE
  ) %>% mutate(
    node_type = ifelse(name %in% state_label_order, "State", "Regulon"),
    state_group = ifelse(node_type == "State", name, "Regulon")
  )

  state_reg_degree <- state_auc_edge_df %>% count(regulon_label, name = "n_states")
  state_auc_node_df <- state_auc_node_df %>%
    left_join(state_reg_degree, by = c("name" = "regulon_label")) %>%
    mutate(is_shared = !is.na(n_states) & n_states > 1)

  state_auc_network_graph <- tbl_graph(
    nodes = state_auc_node_df,
    edges = state_auc_edge_df %>% transmute(from = regulon_label, to = state_label, weight = weight),
    directed = FALSE
  )

  state_auc_network_fill <- c(state_cols, Regulon = "grey35")
  pdf("Auto_final_mp_scenic_state_network_top100auc.pdf", width = 18, height = 12, useDingbats = FALSE)
  print(
    ggraph(state_auc_network_graph, layout = "stress") +
      geom_edge_link(aes(width = weight, alpha = weight), colour = "grey70") +
      scale_edge_width(range = c(0.3, 2.0)) +
      scale_edge_alpha(range = c(0.2, 0.8)) +
      geom_node_point(
        aes(fill = state_group, shape = node_type, size = ifelse(node_type == "State", 7, ifelse(is_shared, 5, 3.5))),
        colour = "black", stroke = 0.3
      ) +
      scale_size_identity() +
      geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 30) +
      scale_shape_manual(values = c(State = 21, Regulon = 22)) +
      scale_fill_manual(values = state_auc_network_fill, drop = FALSE) +
      theme_void(base_size = 12) +
      labs(title = "Top regulons by mean AUC across states") +
      guides(edge_width = "none", edge_alpha = "none", shape = "none", fill = "none")
  )
  dev.off()
}

####################
# Publish all terminal SCENIC objects and figures from the ephemeral working
# directory into live storage. The SCENIC int/ cache remains ephemeral only.
####################
publish_files <- list.files(
  ephemeral_dir,
  pattern = "^Auto_final_mp_scenic_.*\\.(rds|RDS|csv|pdf|png|xlsx)$",
  full.names = TRUE
)
if (length(publish_files) > 0) {
  copied <- file.copy(publish_files, file.path(out_dir, basename(publish_files)), overwrite = TRUE)
  if (!all(copied)) stop("Failed to publish one or more final SCENIC outputs to live storage")
}
####################

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_df <- data.frame(
  mode = "scenic",
  n_common_cells = length(common_cells),
  n_selected_cells = nrow(selected_df),
  n_selected_mps = length(unique(selected_df$final_mp_id)),
  n_selected_states = dplyr::n_distinct(selected_df$final_state[selected_df$final_state %in% state_level_order]),
  n_regulons = nrow(auc_mat),
  db_files = paste(db_files, collapse = " | "),
  stringsAsFactors = FALSE
)
write.csv(
  summary_df,
  file.path(summary_dir, "Auto_final_mp_scenic_summary.csv"),
  row.names = FALSE
)

message("Saved final MP SCENIC outputs in ", out_dir, "; ephemeral computation cache: ", ephemeral_dir)
