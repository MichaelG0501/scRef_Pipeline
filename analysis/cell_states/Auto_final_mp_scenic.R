####################
# Auto_final_mp_scenic.R
#
# Final-MP-focused SCENIC workflow for OAC epithelial cells.
# Uses the curated final MP panel:
#   - retained scATLAS MPs, including the three CC-associated MPs
#   - the three retained pan-cancer 3CA MPs from unresolved relabeling
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_final_states.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/UCell_3CA_MPs.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
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
#   Rscript analysis/cell_states/Auto_final_mp_scenic.R
#   Rscript analysis/cell_states/Auto_final_mp_scenic.R prepare_only=true
#   Rscript analysis/cell_states/Auto_final_mp_scenic.R db_dir=/path/to/cistarget
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

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

format_3ca_label <- function(x) {
  x <- clean_3ca_name(x)
  x <- sub("^3CA_mp_([0-9]+) ", "3CA MP\\1 ", x)
  x
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

load_retained_3ca <- function(coverage_path, state_path, ucell_3ca) {
  if (file.exists(coverage_path)) {
    coverage_df <- read.csv(coverage_path, stringsAsFactors = FALSE)
    retained <- coverage_df %>%
      filter(
        n_samples >= 50,
        n_studies >= 6,
        pct_cells >= 1
      ) %>%
      arrange(desc(n_cells)) %>%
      pull(mp_label)
    retained <- retained[retained %in% colnames(ucell_3ca)]
    if (length(retained) > 0) {
      return(unique(retained))
    }
  }

  if (!file.exists(state_path)) {
    stop("Could not find 3CA coverage CSV or fallback state file.")
  }

  state_B <- readRDS(state_path)
  common_cells <- intersect(names(state_B), rownames(ucell_3ca))
  state_B <- state_B[common_cells]
  unresolved_cells <- names(state_B)[state_B == "Unresolved"]
  if (length(unresolved_cells) == 0) {
    stop("Fallback unresolved-cell 3CA retention failed: no unresolved cells available.")
  }

  cc_fixed <- c(
    "X3CA_mp_1.Cell.Cycle...G2.M",
    "X3CA_mp_2.Cell.Cycle...G1.S",
    "X3CA_mp_3.Cell.Cylce.HMG.rich",
    "X3CA_mp_4.Chromatin",
    "X3CA_mp_5.Cell.cycle.single.nucleus"
  )
  unresolved_3ca <- ucell_3ca[unresolved_cells, setdiff(colnames(ucell_3ca), cc_fixed), drop = FALSE]
  top_3ca_mp <- colnames(unresolved_3ca)[max.col(unresolved_3ca, ties.method = "first")]
  mp_counts <- sort(table(top_3ca_mp), decreasing = TRUE)
  names(mp_counts)[seq_len(min(3, length(mp_counts)))]
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
  db_dir <- file.path(getwd(), "final_mp_scenic", "cistarget_databases")
}
db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)

out_dir <- "final_mp_scenic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Final MP definitions
####################
sc_mps <- c(
  "MP1", "MP7", "MP9",
  "MP2",
  "MP17", "MP14", "MP5", "MP10", "MP8",
  "MP13", "MP12",
  "MP18", "MP16",
  "MP15"
)

sc_mp_descriptions <- c(
  "MP1" = "G2M Cell Cycle",
  "MP7" = "DNA Damage Repair",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epithelium",
  "MP5" = "Epithelial IFN Response",
  "MP10" = "Columnar Differentiation",
  "MP8" = "Intestinal Differentiation",
  "MP13" = "Hypoxic Inflammatory Epithelium",
  "MP12" = "Neuro-responsive Epithelium",
  "MP18" = "Secretory Differentiation (Intestinal)",
  "MP16" = "Secretory Differentiation (Gastric)",
  "MP15" = "Immune Infiltration"
)

mp_group_map <- c(
  "MP1" = "Cell Cycle",
  "MP7" = "Cell Cycle",
  "MP9" = "Cell Cycle",
  "MP2" = "Classic Proliferative",
  "MP17" = "Basal to Intestinal Metaplasia",
  "MP14" = "Basal to Intestinal Metaplasia",
  "MP5" = "Basal to Intestinal Metaplasia",
  "MP10" = "Basal to Intestinal Metaplasia",
  "MP8" = "Basal to Intestinal Metaplasia",
  "MP13" = "Stress-adaptive",
  "MP12" = "Stress-adaptive",
  "MP18" = "SMG-like Metaplasia",
  "MP16" = "SMG-like Metaplasia",
  "MP15" = "Immune Infiltrating"
)

group_cols <- c(
  "Cell Cycle" = "#D4AF37",
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "Pan-cancer 3CA" = "#6A3D9A"
)

state_level_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

####################
# Load inputs
####################
tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states <- readRDS("Auto_final_states.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

three_ca_df <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
three_ca_list <- as.list(three_ca_df)
three_ca_list <- lapply(three_ca_list, function(x) unique(x[!is.na(x) & nzchar(x)]))
names(three_ca_list) <- make.names(sub("^MP", "3CA_mp_", names(three_ca_list)))

coverage_path <- "task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv"
retained_3ca <- load_retained_3ca(
  coverage_path = coverage_path,
  state_path = "Auto_topmp_v2_noreg_states_B.rds",
  ucell_3ca = ucell_3ca
)
retained_3ca <- retained_3ca[retained_3ca %in% names(three_ca_list)]
if (length(retained_3ca) == 0) {
  stop("No retained 3CA MPs could be resolved from coverage results.")
}

final_mp_ids <- c(sc_mps, retained_3ca)
display_label_map <- c(
  setNames(paste(sc_mps, sc_mp_descriptions[sc_mps]), sc_mps),
  setNames(vapply(retained_3ca, format_3ca_label, character(1)), retained_3ca)
)
mp_group_map <- c(
  mp_group_map,
  setNames(rep("Pan-cancer 3CA", length(retained_3ca)), retained_3ca)
)

####################
# Score assembly and cell selection
####################
selection_cache <- file.path(cache_dir, "mp_cell_selection.rds")

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
      rownames(ucell_scores),
      rownames(ucell_3ca)
    )
  )
  if (length(common_cells) == 0) {
    stop("No overlapping cells across required inputs.")
  }

  tmdata_all <- tmdata_all[, common_cells]
  final_states <- final_states[common_cells]
  ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
  ucell_3ca <- ucell_3ca[common_cells, retained_3ca, drop = FALSE]

  missing_sc <- setdiff(sc_mps, colnames(ucell_scores))
  if (length(missing_sc) > 0) {
    stop(
      "Selected final scATLAS MPs are missing from UCell_nMP19_filtered.rds: ",
      paste(missing_sc, collapse = ", ")
    )
  }

  score_mat <- cbind(
    as.matrix(ucell_scores[, sc_mps, drop = FALSE]),
    as.matrix(ucell_3ca[, retained_3ca, drop = FALSE])
  )
  score_mat <- score_mat[, final_mp_ids, drop = FALSE]

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
sc_mp_gene_sets <- geneNMF.metaprograms$metaprograms.genes[sc_mps]
sc_mp_gene_sets <- lapply(sc_mp_gene_sets, function(x) unique(x[!is.na(x) & nzchar(x)]))
sc_mp_gene_sets <- lapply(sc_mp_gene_sets, function(x) head(x, top_genes_per_set))

three_ca_gene_sets <- three_ca_list[retained_3ca]
three_ca_gene_sets <- lapply(three_ca_gene_sets, function(x) head(x, top_genes_per_set))

final_gene_sets <- c(sc_mp_gene_sets, three_ca_gene_sets)
final_gene_sets <- final_gene_sets[final_mp_ids]
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
    "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
    "updates", "new_updates", "summaries"
  )
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- data.frame(
    mode = "prepare_only",
    n_common_cells = length(common_cells),
    n_selected_cells = nrow(selected_df),
    n_selected_mps = length(unique(selected_df$final_mp_id)),
    retained_3ca = paste(format_3ca_label(retained_3ca), collapse = " | "),
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
setwd(out_dir)
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
mp_gene_sets_by_id <- c(sc_mp_gene_sets, three_ca_gene_sets)
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

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_df <- data.frame(
  mode = "scenic",
  n_common_cells = length(common_cells),
  n_selected_cells = nrow(selected_df),
  n_selected_mps = length(unique(selected_df$final_mp_id)),
  n_selected_states = dplyr::n_distinct(selected_df$final_state[selected_df$final_state %in% state_level_order]),
  retained_3ca = paste(format_3ca_label(retained_3ca), collapse = " | "),
  n_regulons = nrow(auc_mat),
  db_files = paste(db_files, collapse = " | "),
  stringsAsFactors = FALSE
)
write.csv(
  summary_df,
  file.path(summary_dir, "Auto_final_mp_scenic_summary.csv"),
  row.names = FALSE
)

message("Saved final MP SCENIC outputs in ", file.path(getwd()))
