####################
# Auto_mp_cross_celltype_correlations.R
# Build a cross-celltype MP co-occurrence network from the full
# EAC_Ref_merged_strict.rds atlas, annotate positive edges with
# ligand-receptor support, and export reproducible summary tables.
#
# Inputs:
#   ref_outs/EAC_Ref_merged_strict.rds
#   ref_outs/meta_full_epi.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/nmf_*/MP_outs_default.rds
#   ref_outs/nmf_*/UCell_default.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx
#
# Main outputs:
#   ref_outs/non_malignant_mp_correlations/
#   ref_outs/non_malignant_mp_correlations/cache/
#   updates/new_updates/summaries/Auto_mp_cross_celltype_correlations_summary.csv
#
# Notes:
#   - Uses the complete atlas for sample/cell matching and LR expression ranking.
#   - Cancer cells are the malignant epithelial subset defined in meta_full_epi.rds.
#   - LR filtering keeps only "literature supported" and "putative" pairs.
#   - Large intermediate objects are cached as RDS files and reused when present.
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(Matrix)
library(scales)
library(tibble)
library(purrr)
library(grid)
library(openxlsx)

####################
# Paths, parameters, and display metadata
####################

resolve_project_dir <- function() {
  candidate_dirs <- c(
    "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
  )
  keep_dirs <- candidate_dirs[file.exists(candidate_dirs)]
  if (length(keep_dirs) == 0) {
    stop("Could not resolve project directory")
  }
  keep_dirs[1]
}

project_dir <- resolve_project_dir()
ref_dir <- file.path(project_dir, "ref_outs")
setwd(ref_dir)

out_dir <- "non_malignant_mp_correlations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(project_dir, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
ucell_cutoff <- if (length(args) >= 1 && nzchar(args[1])) as.numeric(args[1]) else 0.25

if (!is.finite(ucell_cutoff) || ucell_cutoff <= 0 || ucell_cutoff >= 1) {
  stop("UCell cutoff must be a numeric value between 0 and 1")
}

cache_version <- "2026-04-17_v5"
force_rebuild <- toupper(Sys.getenv("AUTO_MPXCELL_FORCE_REBUILD", "FALSE")) %in% c("TRUE", "1", "YES")
pair_evidence_allowed <- c("literature supported", "putative")
min_positive_samples <- 5
min_shared_samples_per_study <- 10
min_pair_samples <- 3
positive_sig_cutoff <- 4
negative_sig_cutoff <- -log10(0.05)
top_ranked_genes <- 4000
cutoff_grid <- sort(unique(c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, ucell_cutoff)))
annotated_positive_edges_n <- 18
state_marker_top_n <- 100
final_cancer_state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)
final_cancer_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

celltype_display_order <- c("cancer", "fibroblast", "endothelial", "cd8", "cd4", "macrophage", "nk", "plasma")
focal_page_order <- celltype_display_order[-length(celltype_display_order)]

publication_colour_lookup <- c(
  "cancer" = "#984EA3",
  "fibroblast" = "#FF7F00",
  "endothelial" = "#4DAF4A",
  "cd8" = "#56B4E9",
  "cd4" = "#00CED1",
  "macrophage" = "#A65628",
  "nk" = "#FFD700",
  "plasma" = "#999999"
)

cancer_state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)
cancer_cc_mps <- c("MP1", "MP7", "MP9")

celltype_cfg <- tibble(
  compartment = c("cancer", "fibroblast", "endothelial", "cd8", "cd4", "macrophage", "nk", "plasma"),
  display = c("cancer", "fibroblast", "endothelial", "cd8", "cd4", "macrophage", "nk", "plasma"),
  atlas_celltype = c("epithelial", "fibroblast", "endothelial", "t.cell", "t.cell", "macrophage", "nk.cell", "plasma"),
  subtype = c(NA, NA, NA, "cd8", "cd4", NA, NA, NA),
  mp_path = c(
    file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds"),
    file.path("nmf_fibroblast", "MP_outs_default.rds"),
    file.path("nmf_endothelial", "MP_outs_default.rds"),
    file.path("nmf_cd8", "MP_outs_default.rds"),
    file.path("nmf_cd4", "MP_outs_default.rds"),
    file.path("nmf_macrophage", "MP_outs_default.rds"),
    file.path("nmf_nk", "MP_outs_default.rds"),
    file.path("nmf_plasma", "MP_outs_default.rds")
  ),
  score_path = c(
    file.path("Metaprogrammes_Results", "UCell_nMP19_filtered.rds"),
    file.path("nmf_fibroblast", "UCell_default.rds"),
    file.path("nmf_endothelial", "UCell_default.rds"),
    file.path("nmf_cd8", "UCell_default.rds"),
    file.path("nmf_cd4", "UCell_default.rds"),
    file.path("nmf_macrophage", "UCell_default.rds"),
    file.path("nmf_nk", "UCell_default.rds"),
    file.path("nmf_plasma", "UCell_default.rds")
  ),
  plot_colour = unname(publication_colour_lookup[c("cancer", "fibroblast", "endothelial", "cd8", "cd4", "macrophage", "nk", "plasma")])
)

cancer_mp_descriptions <- c(
  "MP1" = "G2M Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP5" = "Epithelial IFN Resp.",
  "MP7" = "DNA Damage Repair",
  "MP8" = "Intestinal Diff.",
  "MP9" = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)"
)

pair_order_lookup <- setNames(seq_along(celltype_display_order), celltype_display_order)
plot_colour_lookup <- publication_colour_lookup

####################
# Core helpers
####################

derive_study <- function(sample_ids) {
  sample_ids <- as.character(sample_ids)
  vapply(
    strsplit(sample_ids, "_", fixed = TRUE),
    function(bits) {
      bits <- bits[bits != ""]
      if (length(bits) >= 2) {
        paste(bits[1:2], collapse = "_")
      } else if (length(bits) == 1) {
        bits[1]
      } else {
        NA_character_
      }
    },
    character(1)
  )
}

safe_cor_test <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || stats::var(x) == 0 || stats::var(y) == 0) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }

  test_res <- tryCatch(
    {
      if (identical(method, "spearman")) {
        stats::cor.test(x, y, method = method, exact = FALSE)
      } else {
        stats::cor.test(x, y, method = method)
      }
    },
    error = function(e) NULL
  )

  if (is.null(test_res)) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }

  list(
    estimate = unname(test_res$estimate),
    p.value = test_res$p.value
  )
}

sort_mp_names <- function(mp_names) {
  mp_names[order(as.numeric(gsub("\\D", "", mp_names)))]
}

order_cancer_mps <- function(mp_outs, keep_mps) {
  keep_mps <- unique(keep_mps)
  state_ordered_mps <- unlist(cancer_state_groups, use.names = FALSE)

  orig_order <- NULL
  if (!is.null(mp_outs$programs.tree$order) && !is.null(mp_outs$programs.clusters)) {
    orig_tree_order <- mp_outs$programs.tree$order
    orig_clusters <- mp_outs$programs.clusters[orig_tree_order]
    orig_order <- paste0("MP", unique(orig_clusters))
  }
  if (is.null(orig_order) || length(orig_order) == 0) {
    orig_order <- sort_mp_names(keep_mps)
  }

  orig_order <- orig_order[orig_order %in% keep_mps]
  unique(c(
    orig_order[orig_order %in% cancer_cc_mps],
    state_ordered_mps[state_ordered_mps %in% keep_mps],
    orig_order
  ))
}

focal_direction_binary <- function(focal_is_driver, match_mode) {
  ifelse(
    (focal_is_driver & match_mode == "driver_top_ligand_target_mp_receptor") |
      (!focal_is_driver & match_mode == "driver_top_receptor_target_mp_ligand"),
    1L,
    0L
  )
}

focal_direction_label <- function(direction_binary) {
  ifelse(direction_binary == 1L, "ligand", "receptor")
}

get_focal_pairs <- function(focal, available_pairs) {
  focal_idx <- match(focal, celltype_display_order)
  partner_levels <- celltype_display_order[seq.int(focal_idx + 1, length(celltype_display_order))]
  pair_levels <- paste(focal, partner_levels, sep = " vs ")
  pair_levels[pair_levels %in% available_pairs]
}

get_all_focal_pairs <- function(focal, available_pairs) {
  pair_levels <- vapply(
    setdiff(celltype_display_order, focal),
    function(partner) {
      bits <- c(focal, partner)
      bits <- bits[order(pair_order_lookup[bits])]
      paste(bits, collapse = " vs ")
    },
    character(1)
  )
  unique(pair_levels[pair_levels %in% available_pairs])
}

####################
# Generalized focal-pair ordering used by both the bubble plot and the
# LR workbook exports. When include_within is TRUE, the focal cell type
# is shown against itself first, followed by all cross-celltype pairs in
# the fixed display order.
####################
get_expected_focal_pairs <- function(focal, include_within = FALSE, available_pairs = NULL) {
  pair_levels <- character(0)
  if (isTRUE(include_within)) {
    pair_levels <- c(pair_levels, paste(focal, focal, sep = " vs "))
  }
  pair_levels <- c(
    pair_levels,
    vapply(
      setdiff(celltype_display_order, focal),
      function(partner) {
        bits <- c(focal, partner)
        bits <- bits[order(pair_order_lookup[bits])]
        paste(bits, collapse = " vs ")
      },
      character(1)
    )
  )
  pair_levels <- unique(pair_levels)
  if (!is.null(available_pairs)) {
    pair_levels <- pair_levels[pair_levels %in% available_pairs]
  }
  pair_levels
}

get_partner_celltype <- function(pair_bits, focal) {
  if (length(pair_bits) == 2 && identical(pair_bits[1], pair_bits[2])) {
    focal
  } else {
    setdiff(pair_bits, focal)
  }
}

make_sample_meta <- function(meta_df) {
  meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
  meta_df <- tibble::rownames_to_column(meta_df, "cell")

  study_vec <- if ("study" %in% colnames(meta_df)) as.character(meta_df$study) else rep(NA_character_, nrow(meta_df))
  study_vec[is.na(study_vec) | study_vec == ""] <- derive_study(meta_df$orig.ident[is.na(study_vec) | study_vec == ""])

  tibble(
    cell = meta_df$cell,
    sample = as.character(meta_df$orig.ident),
    study = study_vec
  )
}

filter_gene_sets_by_silhouette <- function(mp_outs) {
  mp_genes <- mp_outs$metaprograms.genes
  sil <- mp_outs$metaprograms.metrics$silhouette
  sil_names <- paste0("MP", seq_along(sil))
  names(sil) <- sil_names
  keep_names <- sil_names[!is.na(sil) & sil >= 0]
  mp_genes[intersect(names(mp_genes), keep_names)]
}

get_assay_data_safely <- function(seurat_obj, assay = "RNA") {
  tryCatch(
    GetAssayData(seurat_obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(seurat_obj, assay = assay, slot = "data")
  )
}

format_cancer_mp_label <- function(mp_name) {
  if (mp_name %in% names(cancer_mp_descriptions)) {
    unname(cancer_mp_descriptions[mp_name])
  } else {
    mp_name
  }
}

format_mp_display <- function(compartment, mp_name) {
  if (identical(compartment, "cancer")) {
    format_cancer_mp_label(mp_name)
  } else {
    mp_name
  }
}

format_node_label <- function(compartment, display_name, mp_name) {
  paste(format_mp_display(compartment, mp_name), display_name)
}

make_node_id <- function(compartment, feature_name) {
  paste0(
    compartment,
    "__",
    if (identical(compartment, "cancer") && feature_name %in% final_cancer_state_order) {
      make.names(feature_name)
    } else {
      feature_name
    }
  )
}

calc_adjusted_scores <- function(score_df, sample_meta, mp_names, cutoff) {
  common_cells <- intersect(rownames(score_df), sample_meta$cell)
  if (length(common_cells) == 0) {
    stop("No overlapping cells between score matrix and metadata")
  }

  sample_meta <- sample_meta[match(common_cells, sample_meta$cell), , drop = FALSE]
  score_mat <- as.matrix(score_df[common_cells, mp_names, drop = FALSE])
  positive_mat <- score_mat > cutoff

  sample_factor <- factor(sample_meta$sample, levels = unique(sample_meta$sample))
  sample_counts <- as.numeric(table(sample_factor))
  names(sample_counts) <- levels(sample_factor)

  adj_sum <- rowsum(positive_mat * 1, group = sample_factor, reorder = FALSE)
  adj_pct <- sweep(adj_sum, 1, sample_counts[rownames(adj_sum)], "/") * 100
  adj_pct <- as.data.frame(adj_pct, stringsAsFactors = FALSE)
  adj_pct$sample <- rownames(adj_pct)
  adj_pct$study <- derive_study(adj_pct$sample)
  adj_pct$n_cells <- sample_counts[adj_pct$sample]

  coverage_counts <- colSums(as.matrix(adj_pct[, mp_names, drop = FALSE]) > 0, na.rm = TRUE)

  list(
    adjusted_scores = adj_pct,
    coverage_counts = coverage_counts,
    sample_counts = sample_counts
  )
}

####################
# Assignment-based adjusted scores for finalized cancer-state labels.
# These cells are already discretely labeled, so positivity is defined
# as state membership rather than a UCell threshold.
####################
calc_adjusted_scores_from_positive_mat <- function(positive_mat, sample_meta, mp_names) {
  common_cells <- intersect(rownames(positive_mat), sample_meta$cell)
  if (length(common_cells) == 0) {
    stop("No overlapping cells between positive matrix and metadata")
  }

  sample_meta <- sample_meta[match(common_cells, sample_meta$cell), , drop = FALSE]
  positive_mat <- as.matrix(positive_mat[common_cells, mp_names, drop = FALSE]) > 0

  sample_factor <- factor(sample_meta$sample, levels = unique(sample_meta$sample))
  sample_counts <- as.numeric(table(sample_factor))
  names(sample_counts) <- levels(sample_factor)

  adj_sum <- rowsum(positive_mat * 1, group = sample_factor, reorder = FALSE)
  adj_pct <- sweep(adj_sum, 1, sample_counts[rownames(adj_sum)], "/") * 100
  adj_pct <- as.data.frame(adj_pct, stringsAsFactors = FALSE)
  adj_pct$sample <- rownames(adj_pct)
  adj_pct$study <- derive_study(adj_pct$sample)
  adj_pct$n_cells <- sample_counts[adj_pct$sample]

  coverage_counts <- colSums(as.matrix(adj_pct[, mp_names, drop = FALSE]) > 0, na.rm = TRUE)

  list(
    adjusted_scores = adj_pct,
    coverage_counts = coverage_counts,
    sample_counts = sample_counts
  )
}

calc_cutoff_sensitivity <- function(score_df, sample_meta, mp_names, cutoffs, compartment, display_name) {
  common_cells <- intersect(rownames(score_df), sample_meta$cell)
  if (length(common_cells) == 0) {
    return(tibble())
  }

  sample_meta <- sample_meta[match(common_cells, sample_meta$cell), , drop = FALSE]
  score_mat <- as.matrix(score_df[common_cells, mp_names, drop = FALSE])
  sample_factor <- factor(sample_meta$sample, levels = unique(sample_meta$sample))
  denominator_samples <- length(levels(sample_factor))

  bind_rows(lapply(cutoffs, function(cutoff) {
    positive_mat <- score_mat > cutoff
    sample_positive <- rowsum(positive_mat * 1, group = sample_factor, reorder = FALSE)
    coverage_n <- colSums(sample_positive > 0, na.rm = TRUE)
    positive_cell_fraction <- colMeans(positive_mat, na.rm = TRUE) * 100

    tibble(
      compartment = compartment,
      celltype_display = display_name,
      cutoff = cutoff,
      mp_name = mp_names,
      mp_display = vapply(mp_names, function(mp) format_mp_display(compartment, mp), character(1)),
      mp_index = as.numeric(gsub("\\D", "", mp_names)),
      node_label = vapply(mp_names, function(mp) format_node_label(compartment, display_name, mp), character(1)),
      denominator_samples = denominator_samples,
      sample_coverage_n = as.numeric(coverage_n[mp_names]),
      sample_coverage_pct = 100 * as.numeric(coverage_n[mp_names]) / denominator_samples,
      positive_cell_fraction_pct = as.numeric(positive_cell_fraction[mp_names]),
      coverage_pass = as.numeric(coverage_n[mp_names]) > min_positive_samples
    )
  }))
}

coerce_support_present <- function(x) {
  if (is.logical(x)) {
    !is.na(x) & x
  } else {
    x_chr <- trimws(as.character(x))
    !is.na(x_chr) & x_chr != "" & toupper(x_chr) != "FALSE"
  }
}

clean_lr_columns <- function(lr_df) {
  lr_df <- as.data.frame(lr_df, stringsAsFactors = FALSE)
  colnames(lr_df) <- gsub("\\.+", "_", colnames(lr_df))
  lower_names <- tolower(colnames(lr_df))

  pair_idx <- match("pair_name", lower_names)
  ligand_idx <- match("ligand_approvedsymbol", lower_names)
  receptor_idx <- match("receptor_approvedsymbol", lower_names)
  source_idx <- match("pair_source", lower_names)
  evidence_idx <- match("pair_evidence", lower_names)

  if (is.na(ligand_idx)) ligand_idx <- grep("^ligand", lower_names)[1]
  if (is.na(receptor_idx)) receptor_idx <- grep("^receptor", lower_names)[1]

  if (is.na(ligand_idx) || is.na(receptor_idx)) {
    char_cols <- which(vapply(lr_df, function(x) is.character(x) || is.factor(x), logical(1)))
    if (length(char_cols) < 2) {
      stop("Could not identify ligand and receptor columns in LR table")
    }
    ligand_idx <- char_cols[1]
    receptor_idx <- char_cols[2]
  }

  support_defs <- c(
    "DLRP" = "dlrp",
    "HPMR" = "hpmr",
    "IUPHAR" = "iuphar",
    "HPRD" = "hprd",
    "STRING.binding" = "string_binding",
    "STRING.experiment" = "string_experiment",
    "PMID.Manual" = "pmid_manual"
  )

  support_present <- lapply(names(support_defs), function(nm) {
    idx <- match(support_defs[[nm]], lower_names)
    if (is.na(idx)) {
      rep(FALSE, nrow(lr_df))
    } else {
      coerce_support_present(lr_df[[idx]])
    }
  })
  names(support_present) <- names(support_defs)
  support_mat <- as.data.frame(support_present, stringsAsFactors = FALSE)

  support_sources <- vapply(seq_len(nrow(support_mat)), function(i) {
    active <- names(support_mat)[which(as.logical(unlist(support_mat[i, , drop = TRUE])))]
    if (length(active) == 0) "" else paste(active, collapse = "; ")
  }, character(1))

  out <- tibble(
    pair_name = if (!is.na(pair_idx)) as.character(lr_df[[pair_idx]]) else NA_character_,
    ligand = toupper(trimws(gsub(" ", "", as.character(lr_df[[ligand_idx]])))),
    receptor = toupper(trimws(gsub(" ", "", as.character(lr_df[[receptor_idx]])))),
    pair_source = if (!is.na(source_idx)) as.character(lr_df[[source_idx]]) else NA_character_,
    pair_evidence = if (!is.na(evidence_idx)) as.character(lr_df[[evidence_idx]]) else NA_character_,
    pair_support_n = rowSums(as.data.frame(lapply(support_mat, as.integer), stringsAsFactors = FALSE)),
    pair_support_sources = support_sources
  )

  out %>%
    filter(!is.na(ligand), !is.na(receptor), ligand != "", receptor != "") %>%
    mutate(
      pair_source = ifelse(is.na(pair_source) | pair_source == "", "unknown", pair_source),
      pair_evidence = ifelse(is.na(pair_evidence) | pair_evidence == "", "unknown", pair_evidence)
    ) %>%
    distinct()
}

load_ramilowski_reference <- function(project_dir) {
  env_path <- Sys.getenv("RAMILOWSKI_LR_PATH", unset = "")
  candidate_paths <- c(
    env_path,
    "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx",
    file.path(project_dir, "ref_outs", "references", "ncomms8866-s3.xlsx"),
    file.path(project_dir, "ref_outs", "references", "Ramilowski_2015_Supplementary_Data_2.xlsx"),
    file.path(project_dir, "analysis", "references", "ncomms8866-s3.xlsx"),
    file.path(project_dir, "analysis", "references", "Ramilowski_2015_Supplementary_Data_2.xlsx")
  )
  candidate_paths <- candidate_paths[nzchar(candidate_paths)]

  for (candidate in candidate_paths) {
    if (!file.exists(candidate)) next

    ext <- tolower(tools::file_ext(candidate))
    lr_df <- NULL
    sheet_used <- NA_character_

    lr_df <- tryCatch(
      {
        if (ext %in% c("xlsx", "xls")) {
          sheets <- readxl::excel_sheets(candidate)
          sheet_used <- if ("All.Pairs" %in% sheets) "All.Pairs" else sheets[1]
          readxl::read_excel(candidate, sheet = sheet_used)
        } else if (ext == "csv") {
          read.csv(candidate, stringsAsFactors = FALSE, check.names = FALSE)
        } else if (ext %in% c("tsv", "txt")) {
          read.delim(candidate, stringsAsFactors = FALSE, check.names = FALSE)
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )

    if (!is.null(lr_df)) {
      cleaned_all <- tryCatch(clean_lr_columns(lr_df), error = function(e) NULL)
      if (!is.null(cleaned_all) && nrow(cleaned_all) > 0) {
        cleaned_keep <- cleaned_all %>%
          filter(pair_evidence %in% pair_evidence_allowed)

        return(list(
          source = candidate,
          sheet = sheet_used,
          pairs = cleaned_keep,
          n_pairs_all = nrow(cleaned_all),
          n_pairs_retained = nrow(cleaned_keep),
          n_pairs_excluded = nrow(cleaned_all) - nrow(cleaned_keep)
        ))
      }
    }
  }

  NULL
}

####################
# State-based cancer definition support. The finalized cancer-state
# labels come from Auto_final_states.rds, while LR target gene sets are
# derived from the ranked recurrent state-marker table.
####################
load_cancer_state_reference <- function() {
  state_path <- "Auto_final_states.rds"
  marker_candidates <- c(
    file.path("Auto_six_state_markers", "Auto_six_state_markers_ranked.csv"),
    file.path("Auto_six_state_markers", "Auto_six_state_markers_final.csv")
  )
  marker_path <- marker_candidates[file.exists(marker_candidates)][1]

  if (!file.exists(state_path)) {
    stop("Missing finalized cancer states: ", state_path)
  }
  if (is.na(marker_path) || !nzchar(marker_path)) {
    stop("Missing ranked cancer-state marker table under ref_outs/Auto_six_state_markers/")
  }

  state_labels <- readRDS(state_path)
  state_label_names <- names(state_labels)
  state_labels <- as.character(unname(state_labels))
  names(state_labels) <- state_label_names
  state_labels <- state_labels[!is.na(state_labels) & nzchar(state_labels)]

  marker_df <- read.csv(marker_path, stringsAsFactors = FALSE)
  if (!all(c("state", "gene") %in% colnames(marker_df))) {
    stop("Cancer-state marker table does not contain required columns: state, gene")
  }

  state_keep <- final_cancer_state_order[final_cancer_state_order %in% unique(state_labels)]
  marker_df <- marker_df %>%
    filter(state %in% state_keep)

  if ("best_state_match" %in% colnames(marker_df)) {
    marker_df <- marker_df %>%
      filter(is.na(best_state_match) | best_state_match)
  }

  marker_df <- marker_df %>%
    mutate(
      ranking_score_sort = if ("ranking_score" %in% colnames(marker_df)) -ranking_score else NA_real_,
      reproducibility_score_sort = if ("reproducibility_score" %in% colnames(marker_df)) -reproducibility_score else NA_real_,
      specificity_gap_sort = if ("specificity_gap" %in% colnames(marker_df)) -specificity_gap else NA_real_,
      median_log2FC_hit_sort = if ("median_log2FC_hit" %in% colnames(marker_df)) -median_log2FC_hit else NA_real_
    ) %>%
    arrange(
      state,
      ranking_score_sort,
      reproducibility_score_sort,
      specificity_gap_sort,
      median_log2FC_hit_sort,
      gene
    )

  state_gene_sets <- marker_df %>%
    group_by(state) %>%
    slice_head(n = state_marker_top_n) %>%
    summarise(gene_set = list(unique(toupper(gene))), .groups = "drop")

  gene_sets <- setNames(state_gene_sets$gene_set, state_gene_sets$state)
  gene_sets <- gene_sets[state_keep]

  list(
    state_labels = state_labels,
    state_order = state_keep,
    state_cols = final_cancer_state_cols[state_keep],
    marker_path = marker_path,
    gene_sets = gene_sets
  )
}

extract_top_genes <- function(expr_mat, cells_use, top_n = 4000) {
  cells_use <- intersect(cells_use, colnames(expr_mat))
  if (length(cells_use) == 0) {
    return(character(0))
  }

  gene_means <- Matrix::rowMeans(expr_mat[, cells_use, drop = FALSE])
  gene_means <- sort(gene_means, decreasing = TRUE)
  toupper(names(gene_means)[seq_len(min(top_n, length(gene_means)))])
}

summarise_edge_lr_labels <- function(lr_pairs_df) {
  if (nrow(lr_pairs_df) == 0) {
    return(tibble())
  }

  lr_pairs_df %>%
    arrange(desc(pair_evidence == "literature supported"), desc(pair_support_n), desc(spearman_significance)) %>%
    mutate(pair_label = ifelse(
      !is.na(pair_name) & pair_name != "",
      pair_name,
      paste(ligand, receptor, sep = " -> ")
    )) %>%
    group_by(edge_id) %>%
    summarise(
      n_lr_pairs = n_distinct(pair_label),
      n_literature_supported = sum(pair_evidence == "literature supported", na.rm = TRUE),
      n_putative = sum(pair_evidence == "putative", na.rm = TRUE),
      max_pair_support_n = max(pair_support_n, na.rm = TRUE),
      top_lr_label = paste(head(unique(pair_label), 3), collapse = "; "),
      .groups = "drop"
    )
}

make_group_layout <- function(node_df) {
  node_df <- node_df %>%
    mutate(celltype_display = factor(celltype_display, levels = celltype_display_order)) %>%
    arrange(celltype_display, mp_plot_order, node_label)

  group_list <- split(node_df, node_df$celltype_display, drop = TRUE)
  group_names <- names(group_list)
  group_angles <- seq(0, 2 * pi, length.out = length(group_names) + 1)[-length(group_names) - 1]
  group_radius <- 6
  local_radius <- 1.5

  bind_rows(lapply(seq_along(group_list), function(i) {
    grp_df <- group_list[[i]]
    grp_n <- nrow(grp_df)
    if (grp_n == 1) {
      local_angles <- 0
    } else {
      local_angles <- seq(-pi / 5, pi / 5, length.out = grp_n)
    }

    tibble(
      node_id = grp_df$node_id,
      x = group_radius * cos(group_angles[i]) + local_radius * cos(group_angles[i] + local_angles),
      y = group_radius * sin(group_angles[i]) + local_radius * sin(group_angles[i] + local_angles)
    )
  }))
}

make_empty_pair_plot <- function(pair_name) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No eligible data", size = 4.5, color = "grey45") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_size = 10) +
    labs(title = pair_name)
}

make_pair_bubble_plot <- function(plot_df, pair_name, x_nodes, y_nodes) {
  if (nrow(plot_df) == 0) {
    return(make_empty_pair_plot(pair_name))
  }

  plot_df <- plot_df %>%
    filter(
      is.finite(spearman_significance),
      is.finite(pearson_r),
      !is.na(spearman_sig)
    )

  if (nrow(plot_df) == 0) {
    return(make_empty_pair_plot(pair_name))
  }

  plot_df <- plot_df %>%
    mutate(
      node1_label = factor(as.character(node1_label), levels = x_nodes),
      node2_label = factor(as.character(node2_label), levels = y_nodes)
    )

  ggplot(plot_df, aes(x = node1_label, y = node2_label)) +
    geom_point(
      aes(size = spearman_significance, fill = pearson_r, alpha = spearman_sig),
      shape = 21,
      color = "black",
      stroke = 0.2
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey90") +
    scale_size_continuous(range = c(0.7, 7.5), limits = c(0, NA)) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.20)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = pair_name,
      x = NULL,
      y = NULL,
      fill = "Pearson r",
      size = "-log10\nSpearman p",
      alpha = "Spearman\np < 0.05"
    )
}

compute_layout <- function(edge_df, node_df) {
  if (nrow(edge_df) == 0 || !requireNamespace("igraph", quietly = TRUE)) {
    return(make_group_layout(node_df))
  }

  graph_df <- edge_df %>%
    transmute(
      from = node1_id,
      to = node2_id,
      weight = pmax(spearman_significance, 0.1)
    )

  graph_obj <- igraph::graph_from_data_frame(
    d = graph_df,
    directed = FALSE,
    vertices = node_df$node_id
  )

  set.seed(123)
  lay <- igraph::layout_with_fr(
    graph_obj,
    weights = igraph::E(graph_obj)$weight,
    niter = 5000
  )

  tibble(
    node_id = igraph::V(graph_obj)$name,
    x = lay[, 1],
    y = lay[, 2]
  )
}

plot_network <- function(edge_df, node_df, title_text, subtitle_text, edge_low, edge_high, annotation_df = NULL) {
  if (nrow(edge_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No edges passed the threshold", size = 6) +
        theme_void() +
        labs(title = title_text, subtitle = subtitle_text)
    )
  }

  node_keep <- unique(c(edge_df$node1_id, edge_df$node2_id))
  node_plot_df <- node_df %>%
    filter(node_id %in% node_keep) %>%
    distinct(node_id, .keep_all = TRUE)

  layout_df <- compute_layout(edge_df, node_plot_df) %>%
    mutate(
      x = x * 1.35,
      y = y * 1.20
    )

  edge_plot_df <- edge_df %>%
    left_join(layout_df, by = c("node1_id" = "node_id")) %>%
    rename(x = x, y = y) %>%
    left_join(layout_df, by = c("node2_id" = "node_id"), suffix = c("", "_end")) %>%
    rename(xend = x_end, yend = y_end)

  node_plot_df <- node_plot_df %>%
    left_join(layout_df, by = "node_id")

  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = node_plot_df,
      aes(x = x, y = y, label = node_label),
      size = 3.1,
      box.padding = 0.35,
      point.padding = 0.20,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 123,
      segment.alpha = 0.5
    )
  } else {
    geom_text(
      data = node_plot_df,
      aes(x = x, y = y, label = node_label),
      size = 3.0,
      nudge_y = 0.18
    )
  }

  base_plot <- ggplot() +
    geom_segment(
      data = edge_plot_df,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        color = spearman_significance,
        linewidth = spearman_significance
      ),
      alpha = 0.80,
      lineend = "round"
    ) +
    geom_point(
      data = node_plot_df,
      aes(x = x, y = y, size = coverage_pct, fill = celltype_display),
      shape = 21,
      color = "white",
      stroke = 0.7
    ) +
    label_layer +
    scale_fill_manual(values = plot_colour_lookup) +
    scale_color_gradient(low = edge_low, high = edge_high) +
    scale_size(range = c(4, 11)) +
    scale_linewidth(range = c(0.4, 2.4)) +
    guides(
      fill = guide_legend(override.aes = list(size = 8, shape = 21)),
      size = guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "white", stroke = 0.7)),
      color = guide_colorbar(barheight = unit(50, "pt")),
      linewidth = guide_legend()
    ) +
    coord_equal(clip = "off") +
    theme_void(base_size = 15) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 13, color = "grey25"),
      plot.margin = margin(20, 120, 20, 120)
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      fill = "Node colour",
      size = "Positive sample\ncoverage (%)",
      color = "Interaction\nsignificance\n(-log10 p)",
      linewidth = "Interaction\nsignificance\n(-log10 p)"
    )

  if (is.null(annotation_df) || nrow(annotation_df) == 0) {
    return(base_plot)
  }

  annotation_plot_df <- edge_plot_df %>%
    inner_join(annotation_df, by = "edge_id") %>%
    mutate(
      label_x = (x + xend) / 2,
      label_y = (y + yend) / 2
    )

  base_plot +
    geom_label(
      data = annotation_plot_df,
      aes(x = label_x, y = label_y, label = edge_note),
      size = 2.3,
      alpha = 0.92,
      linewidth = 0.15,
      fill = alpha("white", 0.92),
      label.padding = unit(0.09, "lines")
    )
}

####################
# Plot subtitle helper so state-based cancer modes explicitly state that
# cancer-state positivity comes from finalized labels rather than UCell.
####################
describe_mode_positivity <- function(mode_cfg) {
  if (identical(mode_cfg$cancer_definition, "state")) {
    paste0(
      "Cancer states use finalized cell-state labels; all non-cancer compartments use UCell > ",
      ucell_cutoff
    )
  } else {
    paste0("All compartments use UCell > ", ucell_cutoff)
  }
}

####################
# Compact all-node interaction map. Significant positive interactions
# are drawn once in the lower triangle; dot size gives the percentage of
# eligible samples where both nodes are present and fill gives Spearman rho.
####################
prepare_interaction_dotmap_data <- function(edge_df, node_df) {
  if (nrow(edge_df) == 0) {
    return(tibble())
  }
  if (!"n_lr_pairs" %in% colnames(edge_df)) {
    edge_df$n_lr_pairs <- NA_real_
  }
  if (!"top_lr_label" %in% colnames(edge_df)) {
    edge_df$top_lr_label <- NA_character_
  }
  if (!"co_positive_sample_pct" %in% colnames(edge_df)) {
    edge_df$co_positive_sample_pct <- NA_real_
  }
  if ("edge_label.x" %in% colnames(edge_df)) {
    edge_df$edge_label <- edge_df$edge_label.x
  }
  if ("celltype_pair.x" %in% colnames(edge_df)) {
    edge_df$celltype_pair <- edge_df$celltype_pair.x
  }
  if ("pair_scope.x" %in% colnames(edge_df)) {
    edge_df$pair_scope <- edge_df$pair_scope.x
  }

  node_order_df <- node_df %>%
    arrange(celltype_order, mp_plot_order, node_label) %>%
    mutate(node_order = row_number()) %>%
    select(node_id, node_label, celltype_display, node_order)

  edge_df %>%
    left_join(node_order_df, by = c("node1_id" = "node_id")) %>%
    rename(
      node1_order = node_order,
      node1_plot_label = node_label,
      node1_celltype_plot = celltype_display
    ) %>%
    left_join(node_order_df, by = c("node2_id" = "node_id"), suffix = c("", "_node2")) %>%
    rename(
      node2_order = node_order,
      node2_plot_label = node_label,
      node2_celltype_plot = celltype_display
    ) %>%
    filter(!is.na(node1_order), !is.na(node2_order)) %>%
    mutate(
      x_label = ifelse(node1_order <= node2_order, node1_plot_label, node2_plot_label),
      y_label = ifelse(node1_order <= node2_order, node2_plot_label, node1_plot_label),
      x_order = pmin(node1_order, node2_order),
      y_order = pmax(node1_order, node2_order),
      lr_supported = !is.na(n_lr_pairs) & n_lr_pairs > 0,
      co_positive_sample_pct = ifelse(
        is.na(co_positive_sample_pct),
        NA_real_,
        pmax(0, pmin(100, co_positive_sample_pct))
      )
    ) %>%
    select(
      edge_id,
      edge_label,
      celltype_pair,
      pair_scope,
      compartment1,
      compartment2,
      node1_label,
      node2_label,
      mp1_name,
      mp2_name,
      shared_sample_n,
      co_positive_sample_n,
      co_positive_sample_pct,
      pearson_r,
      spearman_r,
      spearman_p,
      spearman_significance,
      n_lr_pairs,
      top_lr_label,
      x_label,
      y_label,
      x_order,
      y_order,
      lr_supported
    ) %>%
    mutate(
      edge_label = as.character(edge_label),
      celltype_pair = as.character(celltype_pair),
      pair_scope = as.character(pair_scope),
      node1_label = as.character(node1_label),
      node2_label = as.character(node2_label),
      mp1_name = as.character(mp1_name),
      mp2_name = as.character(mp2_name),
      top_lr_label = as.character(top_lr_label)
    )
}

plot_interaction_dotmap <- function(edge_df, node_df, title_text, subtitle_text) {
  ordered_nodes <- node_df %>%
    arrange(celltype_order, mp_plot_order, node_label) %>%
    mutate(node_order = row_number())
  node_levels <- ordered_nodes$node_label
  separator_df <- ordered_nodes %>%
    count(celltype_display, celltype_order, name = "n") %>%
    arrange(celltype_order) %>%
    mutate(separator = cumsum(n) + 0.5) %>%
    filter(separator < max(separator))

  plot_df <- prepare_interaction_dotmap_data(edge_df, node_df)
  if (nrow(plot_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No significant positive interactions", size = 5) +
        theme_void(base_size = 12) +
        labs(title = title_text, subtitle = subtitle_text)
    )
  }

  plot_df <- plot_df %>%
    mutate(
      x_label = factor(x_label, levels = node_levels),
      y_label = factor(y_label, levels = node_levels)
    )

  ggplot(plot_df, aes(x = x_label, y = y_label)) +
    geom_vline(data = separator_df, aes(xintercept = separator), inherit.aes = FALSE, color = "grey88", linewidth = 0.35) +
    geom_hline(data = separator_df, aes(yintercept = separator), inherit.aes = FALSE, color = "grey88", linewidth = 0.35) +
    geom_point(
      aes(size = co_positive_sample_pct, fill = spearman_r),
      shape = 21,
      color = "grey45",
      stroke = 0.25,
      alpha = 0.92
    ) +
    geom_point(
      data = plot_df %>% filter(lr_supported),
      aes(size = co_positive_sample_pct, fill = spearman_r),
      shape = 21,
      color = "black",
      stroke = 0.95,
      alpha = 0.98
    ) +
    geom_point(
      data = plot_df %>% filter(lr_supported),
      aes(size = co_positive_sample_pct),
      shape = 4,
      color = "black",
      stroke = 0.55,
      alpha = 0.85
    ) +
    scale_x_discrete(drop = FALSE, position = "top") +
    scale_y_discrete(drop = FALSE, limits = rev(node_levels)) +
    scale_fill_gradient(
      low = "#F7F7F7",
      high = "#B2182B",
      limits = c(0, if(nrow(plot_df) > 0) max(plot_df$spearman_r, na.rm = TRUE) else 1),
      oob = squish,
      name = "Spearman\nrho"
    ) +
    scale_size_continuous(
      range = c(1.2, 8.5),
      limits = c(0, 100),
      breaks = c(10, 25, 50, 75),
      name = "Co-positive\nsamples (%)"
    ) +
    guides(
      fill = guide_colorbar(barheight = unit(60, "pt")),
      size = guide_legend(override.aes = list(shape = 21, fill = "grey70", color = "grey35"))
    ) +
    coord_fixed(clip = "off") +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 0, vjust = 0.5, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 15, color = "grey25"),
      plot.margin = margin(18, 28, 18, 18)
    ) +
    labs(
      title = title_text
    )
}

prepare_focal_interaction_dotmap_data <- function(edge_df, node_df, focal, include_within = FALSE) {
  if (nrow(edge_df) == 0) {
    return(tibble())
  }
  if (!"n_lr_pairs" %in% colnames(edge_df)) {
    edge_df$n_lr_pairs <- NA_real_
  }
  if (!"top_lr_label" %in% colnames(edge_df)) {
    edge_df$top_lr_label <- NA_character_
  }
  if (!"co_positive_sample_pct" %in% colnames(edge_df)) {
    edge_df$co_positive_sample_pct <- NA_real_
  }
  if ("edge_label.x" %in% colnames(edge_df)) {
    edge_df$edge_label <- edge_df$edge_label.x
  }
  if ("celltype_pair.x" %in% colnames(edge_df)) {
    edge_df$celltype_pair <- edge_df$celltype_pair.x
  }
  if ("pair_scope.x" %in% colnames(edge_df)) {
    edge_df$pair_scope <- edge_df$pair_scope.x
  }

  node_lookup <- node_df %>%
    arrange(celltype_order, mp_plot_order, node_label) %>%
    mutate(node_order = row_number()) %>%
    select(node_id, node_label, celltype_display, celltype_order, mp_plot_order, node_order)

  edge_df %>%
    left_join(node_lookup, by = c("node1_id" = "node_id")) %>%
    rename(
      node1_plot_label = node_label,
      node1_celltype_plot = celltype_display,
      node1_celltype_order = celltype_order,
      node1_mp_plot_order = mp_plot_order,
      node1_order = node_order
    ) %>%
    left_join(node_lookup, by = c("node2_id" = "node_id"), suffix = c("", "_node2")) %>%
    rename(
      node2_plot_label = node_label,
      node2_celltype_plot = celltype_display,
      node2_celltype_order = celltype_order,
      node2_mp_plot_order = mp_plot_order,
      node2_order = node_order
    ) %>%
    filter(node1_celltype_plot == focal | node2_celltype_plot == focal) %>%
    filter(include_within | node1_celltype_plot != node2_celltype_plot) %>%
    mutate(
      focal_label = ifelse(node1_celltype_plot == focal, node1_plot_label, node2_plot_label),
      partner_label = ifelse(node1_celltype_plot == focal, node2_plot_label, node1_plot_label),
      partner_celltype = ifelse(node1_celltype_plot == focal, node2_celltype_plot, node1_celltype_plot),
      partner_celltype_order = ifelse(node1_celltype_plot == focal, node2_celltype_order, node1_celltype_order),
      partner_mp_plot_order = ifelse(node1_celltype_plot == focal, node2_mp_plot_order, node1_mp_plot_order),
      lr_supported = !is.na(n_lr_pairs) & n_lr_pairs > 0,
      co_positive_sample_pct = ifelse(
        is.na(co_positive_sample_pct),
        NA_real_,
        pmax(0, pmin(100, co_positive_sample_pct))
      )
    ) %>%
    mutate(focal_celltype = focal) %>%
    select(
      focal_celltype,
      partner_celltype,
      edge_id,
      edge_label,
      celltype_pair,
      pair_scope,
      focal_label,
      partner_label,
      node1_label,
      node2_label,
      mp1_name,
      mp2_name,
      shared_sample_n,
      co_positive_sample_n,
      co_positive_sample_pct,
      pearson_r,
      spearman_r,
      spearman_p,
      spearman_significance,
      n_lr_pairs,
      top_lr_label,
      lr_supported,
      partner_celltype_order,
      partner_mp_plot_order
    ) %>%
    mutate(
      focal_celltype = as.character(focal_celltype),
      partner_celltype = as.character(partner_celltype),
      edge_label = as.character(edge_label),
      celltype_pair = as.character(celltype_pair),
      pair_scope = as.character(pair_scope),
      focal_label = as.character(focal_label),
      partner_label = as.character(partner_label),
      node1_label = as.character(node1_label),
      node2_label = as.character(node2_label),
      mp1_name = as.character(mp1_name),
      mp2_name = as.character(mp2_name),
      top_lr_label = as.character(top_lr_label)
    )
}

make_focal_interaction_grid <- function(node_df, focal, include_within = FALSE) {
  focal_nodes <- node_df %>%
    filter(celltype_display == focal) %>%
    arrange(mp_plot_order, node_label) %>%
    transmute(focal_label = node_label, focal_order = row_number())

  partner_nodes <- node_df %>%
    filter(include_within | celltype_display != focal) %>%
    arrange(celltype_order, mp_plot_order, node_label) %>%
    transmute(
      partner_label = node_label,
      partner_celltype = celltype_display,
      partner_celltype_order = celltype_order,
      partner_mp_plot_order = mp_plot_order,
      partner_order = row_number()
    )
  partner_celltype_levels <- celltype_display_order[celltype_display_order %in% unique(partner_nodes$partner_celltype)]

  expand.grid(
    focal_label = focal_nodes$focal_label,
    partner_label = partner_nodes$partner_label,
    stringsAsFactors = FALSE
  ) %>%
    left_join(focal_nodes, by = "focal_label") %>%
    left_join(partner_nodes, by = "partner_label") %>%
    mutate(
      focal_label = factor(focal_label, levels = focal_nodes$focal_label),
      partner_label = factor(partner_label, levels = partner_nodes$partner_label),
      partner_celltype = factor(partner_celltype, levels = partner_celltype_levels)
    )
}

plot_focal_interaction_dotmap <- function(edge_df, node_df, focal, include_within = FALSE) {
  focal_grid <- make_focal_interaction_grid(node_df, focal, include_within = include_within)
  focal_levels <- levels(focal_grid$focal_label)
  partner_levels <- levels(focal_grid$partner_label)
  partner_celltype_levels <- levels(focal_grid$partner_celltype)

  plot_df <- prepare_focal_interaction_dotmap_data(edge_df, node_df, focal, include_within = include_within) %>%
    mutate(
      focal_label = factor(focal_label, levels = focal_levels),
      partner_label = factor(partner_label, levels = partner_levels),
      partner_celltype = factor(partner_celltype, levels = partner_celltype_levels)
    )

  base_plot <- ggplot() +
    geom_tile(
      data = focal_grid,
      aes(x = partner_label, y = focal_label),
      fill = "#FBFBF8",
      color = "#ECE8DE",
      linewidth = 0.18
    ) +
    facet_grid(. ~ partner_celltype, scales = "free_x", space = "free_x", drop = FALSE) +
    scale_y_discrete(limits = rev(focal_levels), drop = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 9.5, color = "grey15"),
      axis.text.y = element_text(size = 10.5, color = "grey10"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      strip.background = element_rect(fill = alpha(plot_colour_lookup[focal], 0.12), color = NA),
      strip.text.x = element_text(face = "bold", size = 11, color = "grey10"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 10)),
      plot.caption = element_blank(),
      plot.margin = margin(16, 26, 42, 18)
    ) +
    labs(
      title = paste0(focal, " interaction map")
    )

  if (nrow(plot_df) == 0) {
    return(
      base_plot +
        annotate("text", x = 1, y = 1, label = "No significant positive interactions", size = 4, color = "grey45")
    )
  }

  dot_plot <- base_plot +
    geom_point(
      data = plot_df,
      aes(x = partner_label, y = focal_label, size = co_positive_sample_pct, fill = spearman_r),
      shape = 21,
      color = "grey45",
      stroke = 0.28,
      alpha = 0.93
    ) +
    scale_fill_gradient(
      low = "#F8F4ED",
      high = "#B7212E",
      limits = c(0, if(nrow(plot_df) > 0) max(plot_df$spearman_r, na.rm = TRUE) else 1),
      oob = squish,
      name = "Spearman\nrho"
    ) +
    scale_size_continuous(
      range = c(1.3, 8.6),
      limits = c(0, 100),
      breaks = c(10, 25, 50, 75),
      name = "Co-positive\nsamples (%)"
    ) +
    guides(
      fill = guide_colorbar(barheight = unit(55, "pt"), order = 1),
      size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey70", color = "grey35"))
    )

  lr_plot_df <- plot_df %>% filter(lr_supported)
  if (nrow(lr_plot_df) == 0) {
    return(dot_plot)
  }

  dot_plot +
    geom_point(
      data = lr_plot_df,
      aes(x = partner_label, y = focal_label, size = co_positive_sample_pct, fill = spearman_r),
      shape = 21,
      color = "black",
      stroke = 0.95,
      alpha = 0.98
    ) +
    geom_point(
      data = lr_plot_df,
      aes(x = partner_label, y = focal_label, shape = "Supported"),
      size = 2.15,
      color = "black",
      stroke = 0.55
    ) +
    scale_shape_manual(name = "LR support", values = c("Supported" = 4), drop = FALSE) +
    guides(shape = guide_legend(order = 3, override.aes = list(size = 3.2, color = "black")))
}

write_focal_interaction_dotmap <- function(edge_df, node_df, out_path, include_within = FALSE) {
  pdf(out_path, width = 22, height = 10.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  for (focal in celltype_display_order) {
    print(plot_focal_interaction_dotmap(edge_df, node_df, focal, include_within = include_within))
  }
}

add_co_positive_support <- function(edge_df, compartment_data) {
  if (nrow(edge_df) == 0) {
    return(edge_df)
  }

  co_n <- numeric(nrow(edge_df))
  co_pct <- numeric(nrow(edge_df))
  co_n[] <- NA_real_
  co_pct[] <- NA_real_

  for (edge_idx in seq_len(nrow(edge_df))) {
    edge_row <- edge_df[edge_idx, , drop = FALSE]
    comp1 <- compartment_data[[edge_row$compartment1]]
    comp2 <- compartment_data[[edge_row$compartment2]]
    if (is.null(comp1) || is.null(comp2)) next

    eligible_studies <- strsplit(edge_row$eligible_studies, ";", fixed = TRUE)[[1]]
    eligible_studies <- eligible_studies[nzchar(eligible_studies)]
    if (length(eligible_studies) == 0) next

    samples1 <- comp1$adjusted_scores$sample[comp1$adjusted_scores$study %in% eligible_studies]
    samples2 <- comp2$adjusted_scores$sample[comp2$adjusted_scores$study %in% eligible_studies]
    eligible_samples <- intersect(samples1, samples2)
    if (length(eligible_samples) == 0) next

    score1 <- comp1$adjusted_scores[match(eligible_samples, comp1$adjusted_scores$sample), , drop = FALSE]
    score2 <- comp2$adjusted_scores[match(eligible_samples, comp2$adjusted_scores$sample), , drop = FALSE]
    keep <- !is.na(score1$sample) & !is.na(score2$sample)
    score1 <- score1[keep, , drop = FALSE]
    score2 <- score2[keep, , drop = FALSE]
    if (nrow(score1) == 0) next

    co_n[edge_idx] <- sum(score1[[edge_row$mp1_name]] > 0 & score2[[edge_row$mp2_name]] > 0, na.rm = TRUE)
    co_pct[edge_idx] <- 100 * co_n[edge_idx] / nrow(score1)
  }

  edge_df$co_positive_sample_n <- co_n
  edge_df$co_positive_sample_pct <- co_pct
  edge_df
}

load_or_build_cache <- function(cache_filename, build_fun, step_name, cache_dir_use = cache_dir) {
  dir.create(cache_dir_use, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(cache_dir_use, cache_filename)

  if (!force_rebuild && file.exists(cache_path)) {
    cache_obj <- readRDS(cache_path)
    if (is.list(cache_obj) && identical(cache_obj$cache_version, cache_version)) {
      message("Loading cached ", step_name, " from ", cache_path)
      return(cache_obj$data)
    }
    message("Ignoring stale cache for ", step_name, ": ", cache_path)
  }

  message("Computing ", step_name)
  result <- build_fun()
  saveRDS(
    list(
      cache_version = cache_version,
      created_at = as.character(Sys.time()),
      data = result
    ),
    cache_path
  )
  result
}

prepare_adjusted_score_export <- function(compartment_name, adjusted_scores_df, mp_names) {
  export_df <- adjusted_scores_df %>%
    select(sample, study, n_cells, all_of(mp_names))

  display_map <- setNames(
    vapply(mp_names, function(mp) format_mp_display(compartment_name, mp), character(1)),
    mp_names
  )

  colnames(export_df)[match(mp_names, colnames(export_df))] <- unname(display_map[mp_names])
  export_df
}

make_safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  x <- trimws(x, which = "both")
  substr(x, 1, 31)
}

excel_styles <- list(
  title = createStyle(fontSize = 15, textDecoration = "bold", fontColour = "#1F2937"),
  subtitle = createStyle(fontSize = 10, fontColour = "#4B5563"),
  header = createStyle(
    textDecoration = "bold",
    fgFill = "#1F4E78",
    fontColour = "#FFFFFF",
    halign = "center",
    border = "Bottom",
    borderColour = "#1F4E78"
  ),
  literature = createStyle(fgFill = "#D9EAD3"),
  putative = createStyle(fgFill = "#FFF2CC"),
  group = createStyle(fgFill = "#EAF2F8"),
  body = createStyle(valign = "top"),
  number = createStyle(numFmt = "0.000"),
  integer = createStyle(numFmt = "0"),
  wrap = createStyle(wrapText = TRUE)
)

write_styled_sheet <- function(wb, sheet_name, data_df, title_text, subtitle_text = NULL, tab_colour = NULL, manual_widths = NULL) {
  addWorksheet(wb, sheetName = sheet_name, tabColour = tab_colour, gridLines = FALSE)
  writeData(wb, sheet = sheet_name, x = title_text, startRow = 1, startCol = 1)
  addStyle(wb, sheet = sheet_name, style = excel_styles$title, rows = 1, cols = 1, stack = TRUE)
  table_start_row <- 3
  if (!is.null(subtitle_text) && nzchar(subtitle_text)) {
    writeData(wb, sheet = sheet_name, x = subtitle_text, startRow = 2, startCol = 1)
    addStyle(wb, sheet = sheet_name, style = excel_styles$subtitle, rows = 2, cols = 1, stack = TRUE)
    table_start_row <- 4
  }

  if (nrow(data_df) == 0) {
    writeData(wb, sheet = sheet_name, x = "No rows available for this view.", startRow = table_start_row, startCol = 1)
    return(invisible(NULL))
  }

  writeDataTable(
    wb,
    sheet = sheet_name,
    x = data_df,
    startRow = table_start_row,
    startCol = 1,
    tableStyle = "TableStyleMedium2",
    withFilter = TRUE
  )

  data_rows <- (table_start_row + 1):(nrow(data_df) + table_start_row)

  addStyle(
    wb,
    sheet = sheet_name,
    style = excel_styles$header,
    rows = table_start_row,
    cols = seq_len(ncol(data_df)),
    gridExpand = TRUE,
    stack = TRUE
  )
  addStyle(
    wb,
    sheet = sheet_name,
    style = excel_styles$body,
    rows = data_rows,
    cols = seq_len(ncol(data_df)),
    gridExpand = TRUE,
    stack = TRUE
  )

  if ("Pair evidence" %in% colnames(data_df)) {
    evidence_col <- match("Pair evidence", colnames(data_df))
    evidence_rows <- data_rows
    lit_rows <- evidence_rows[data_df[["Pair evidence"]] == "literature supported"]
    put_rows <- evidence_rows[data_df[["Pair evidence"]] == "putative"]
    if (length(lit_rows) > 0) {
      addStyle(wb, sheet = sheet_name, style = excel_styles$literature, rows = lit_rows, cols = evidence_col, gridExpand = TRUE, stack = TRUE)
    }
    if (length(put_rows) > 0) {
      addStyle(wb, sheet = sheet_name, style = excel_styles$putative, rows = put_rows, cols = evidence_col, gridExpand = TRUE, stack = TRUE)
    }
  }

  numeric_cols <- intersect(c("Pair support n", "Pearson r", "Spearman r", "-log10 Spearman p", "Driver positive cells", "Target positive cells", "Focal positive cells", "Partner positive cells"), colnames(data_df))
  for (nm in numeric_cols) {
    style_use <- if (nm %in% c("Pair support n", "Driver positive cells", "Target positive cells", "Focal positive cells", "Partner positive cells")) excel_styles$integer else excel_styles$number
    addStyle(
      wb,
      sheet = sheet_name,
      style = style_use,
      rows = data_rows,
      cols = match(nm, colnames(data_df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  wrap_cols <- intersect(c("Edge label", "LR pair"), colnames(data_df))
  for (nm in wrap_cols) {
    addStyle(
      wb,
      sheet = sheet_name,
      style = excel_styles$wrap,
      rows = data_rows,
      cols = match(nm, colnames(data_df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  if ("Partner cell type" %in% colnames(data_df)) {
    addStyle(
      wb,
      sheet = sheet_name,
      style = excel_styles$group,
      rows = data_rows,
      cols = match("Partner cell type", colnames(data_df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  if ("-log10 Spearman p" %in% colnames(data_df) && nrow(data_df) >= 2) {
    significance_col <- match("-log10 Spearman p", colnames(data_df))
    conditionalFormatting(
      wb,
      sheet = sheet_name,
      cols = significance_col,
      rows = data_rows,
      style = c("#FEE0D2", "#FC9272", "#DE2D26"),
      type = "colourScale"
    )
  }

  freezePane(wb, sheet = sheet_name, firstActiveRow = table_start_row + 1, firstActiveCol = 2)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data_df)), widths = "auto")
  if (!is.null(manual_widths)) {
    for (nm in names(manual_widths)) {
      if (nm %in% colnames(data_df)) {
        setColWidths(wb, sheet = sheet_name, cols = match(nm, colnames(data_df)), widths = manual_widths[[nm]])
      }
    }
  }
}

prepare_lr_export_table <- function(lr_pairs_df) {
  if (nrow(lr_pairs_df) == 0) {
    return(tibble())
  }

  lr_pairs_df %>%
    mutate(
      pair_label = ifelse(!is.na(pair_name) & pair_name != "", pair_name, paste(ligand, receptor, sep = " -> ")),
      evidence_rank = ifelse(pair_evidence == "literature supported", 1, 2),
      celltype1_order = pair_order_lookup[celltype1_display],
      celltype2_order = pair_order_lookup[celltype2_display]
    ) %>%
    arrange(celltype1_order, celltype2_order, desc(spearman_significance), desc(pair_support_n), evidence_rank, pair_label)
}

write_lr_workbooks <- function(lr_pairs_df, lr_edge_summary, out_dir, include_within = FALSE) {
  if (nrow(lr_pairs_df) == 0) {
    return(character(0))
  }

  lr_export <- prepare_lr_export_table(lr_pairs_df)

  edge_overview <- lr_edge_summary %>%
    mutate(
      celltype1_display = sub(" vs .*", "", celltype_pair),
      celltype2_display = sub(".* vs ", "", celltype_pair),
      celltype1_order = pair_order_lookup[celltype1_display],
      celltype2_order = pair_order_lookup[celltype2_display]
    ) %>%
    arrange(celltype1_order, celltype2_order, desc(n_lr_pairs), desc(edge_label))

  focal_workbook_path <- file.path(out_dir, "Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx")

  wb_focal <- createWorkbook(creator = "Codex")
  write_styled_sheet(
    wb = wb_focal,
    sheet_name = "Overview",
    data_df = edge_overview %>%
      select(edge_label, celltype_pair, n_lr_pairs, top_lr_label) %>%
      rename(
        "Edge label" = edge_label,
        "Celltype pair" = celltype_pair,
        "LR pairs" = n_lr_pairs,
        "Top example pairs" = top_lr_label
      ),
    title_text = "LR overview"
  )

  for (focal in celltype_display_order) {
    focal_pairs <- get_expected_focal_pairs(focal, include_within = include_within)
    focal_df <- lr_export %>%
      filter((celltype1_display == focal | celltype2_display == focal), celltype_pair %in% focal_pairs) %>%
      mutate(
        partner_celltype = as.character(ifelse(celltype1_display == focal, celltype2_display, celltype1_display)),
        partner_order = pair_order_lookup[partner_celltype],
        focal_node = as.character(ifelse(celltype1_display == focal, node1_label, node2_label)),
        partner_node = as.character(ifelse(celltype1_display == focal, node2_label, node1_label)),
        focal_is_driver = driver_node == focal_node,
        direction_binary = focal_direction_binary(focal_is_driver, match_mode),
        direction_label = as.character(focal_direction_label(direction_binary)),
        focal_positive_cells = ifelse(focal_is_driver, driver_positive_cells, target_positive_cells),
        partner_positive_cells = ifelse(focal_is_driver, target_positive_cells, driver_positive_cells)
      ) %>%
      arrange(partner_order, desc(spearman_significance), desc(pair_support_n), evidence_rank, pair_label) %>%
      transmute(
        "Partner cell type" = partner_celltype,
        "Edge label" = edge_label,
        "Focal node" = focal_node,
        "Partner node" = partner_node,
        "Direction" = direction_label,
        "Ligand" = ligand,
        "Receptor" = receptor,
        "LR pair" = pair_label,
        "Pair evidence" = pair_evidence,
        "Pearson r" = pearson_r,
        "Spearman r" = spearman_r,
        "-log10 Spearman p" = spearman_significance,
        "Focal positive cells" = focal_positive_cells,
        "Partner positive cells" = partner_positive_cells
      ) %>%
      distinct()

    present_pairs <- unique(lr_export$celltype_pair[lr_export$celltype1_display == focal | lr_export$celltype2_display == focal])
    missing_pairs <- setdiff(focal_pairs, present_pairs)
    if (length(missing_pairs) > 0) {
      placeholder_df <- bind_rows(lapply(missing_pairs, function(pair_name) {
        partner_celltype <- sub(".* vs ", "", pair_name)
        tibble(
          "Partner cell type" = partner_celltype,
          "Edge label" = pair_name,
          "Focal node" = "",
          "Partner node" = "",
          "Direction" = "",
          "Ligand" = "",
          "Receptor" = "",
          "LR pair" = "No retained LR pair",
          "Pair evidence" = "",
          "Pearson r" = NA_real_,
          "Spearman r" = NA_real_,
          "-log10 Spearman p" = NA_real_,
          "Focal positive cells" = NA_integer_,
          "Partner positive cells" = NA_integer_
        )
      }))

      focal_df <- bind_rows(focal_df, placeholder_df) %>%
        mutate(partner_order = pair_order_lookup[`Partner cell type`]) %>%
        arrange(partner_order, desc(`-log10 Spearman p`), desc(`Pearson r`)) %>%
        select(-partner_order)
    }

    write_styled_sheet(
      wb = wb_focal,
      sheet_name = make_safe_sheet_name(focal),
      data_df = focal_df,
      title_text = focal,
      tab_colour = plot_colour_lookup[[focal]],
      manual_widths = c("Partner cell type" = 20)
    )
  }

  saveWorkbook(wb_focal, focal_workbook_path, overwrite = TRUE)

  focal_workbook_path
}

####################
# Multi-mode analysis runner
####################

####################
# Load the atlas, finalized cancer states, and LR reference once, then
# run four separate analysis modes under distinct output subfolders:
# 1) cancer MPs, cross-celltype only
# 2) cancer MPs, cross + within-celltype
# 3) cancer final states, cross-celltype only
# 4) cancer final states, cross + within-celltype
####################
message("Loading complete atlas, epithelial malignancy metadata, finalized cancer states, and LR reference")

merged_obj_mode <- readRDS("EAC_Ref_merged_strict.rds")
merged_meta_mode <- merged_obj_mode@meta.data
merged_cells_mode <- rownames(merged_meta_mode)
merged_expr_mode_all <- get_assay_data_safely(merged_obj_mode, assay = "RNA")
rm(merged_obj_mode)
invisible(gc())
epi_meta_mode <- readRDS("meta_full_epi.rds")
cancer_state_ref_mode <- load_cancer_state_reference()
lr_ref_mode <- load_ramilowski_reference(project_dir)

malignant_cells_mode <- rownames(epi_meta_mode)[epi_meta_mode$malignancy %in% c("malignant_level_1", "malignant_level_2")]
full_sample_lookup_mode <- make_sample_meta(merged_meta_mode) %>%
  distinct(sample, study)

analysis_modes <- tibble(
  analysis_id = c(
    "cancer_mps_cross_only",
    "cancer_mps_cross_and_within",
    "cancer_states_cross_only",
    "cancer_states_cross_and_within"
  ),
  analysis_label = c(
    "Cancer MPs with cross-celltype correlations",
    "Cancer MPs with cross- and within-celltype correlations",
    "Cancer final states with cross-celltype correlations",
    "Cancer final states with cross- and within-celltype correlations"
  ),
  cancer_definition = c("mp", "mp", "state", "state"),
  include_within = c(FALSE, TRUE, FALSE, TRUE),
  out_subdir = c(
    "01_cancer_mps_cross_only",
    "02_cancer_mps_cross_and_within",
    "03_cancer_states_cross_only",
    "04_cancer_states_cross_and_within"
  )
)

run_analysis_mode <- function(mode_cfg) {
  mode_cfg <- as.list(mode_cfg)
  analysis_out_dir <- file.path(out_dir, mode_cfg$out_subdir)
  analysis_cache_dir <- file.path(analysis_out_dir, "cache")
  dir.create(analysis_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(analysis_cache_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running analysis mode: ", mode_cfg$analysis_id)

  ####################
  # Mode-specific compartment builder
  ####################
  build_compartment_data_mode <- function(cfg_row) {
    message("Preparing compartment: ", cfg_row$compartment, " [", mode_cfg$analysis_id, "]")

    if (cfg_row$compartment == "cancer" && identical(mode_cfg$cancer_definition, "state")) {
      state_labels <- cancer_state_ref_mode$state_labels
      cells_use <- intersect(names(state_labels), merged_cells_mode)
      cells_use <- intersect(cells_use, malignant_cells_mode)
      state_labels <- state_labels[cells_use]
      state_keep <- cancer_state_ref_mode$state_order[cancer_state_ref_mode$state_order %in% unique(state_labels)]
      cells_use <- names(state_labels)[state_labels %in% state_keep]
      state_labels <- state_labels[cells_use]
      if (length(state_keep) == 0 || length(cells_use) == 0) {
        stop("No finalized cancer-state cells available for state-based mode")
      }

      score_mat <- vapply(
        state_keep,
        function(state_name) as.numeric(state_labels == state_name),
        numeric(length(state_labels))
      )
      rownames(score_mat) <- cells_use
      colnames(score_mat) <- state_keep
      score_df <- as.data.frame(score_mat, stringsAsFactors = FALSE)
      sample_meta <- make_sample_meta(merged_meta_mode[cells_use, , drop = FALSE])
      ####################
      # Finalized cancer states are categorical labels, so adjusted scores
      # are defined directly from state membership without a UCell cutoff.
      ####################
      adj_bits <- calc_adjusted_scores_from_positive_mat(score_df, sample_meta, state_keep)
      cutoff_tbl <- tibble()

      node_stats <- tibble(
        compartment = cfg_row$compartment,
        celltype_display = cfg_row$display,
        celltype_order = unname(pair_order_lookup[cfg_row$display]),
        node_id = vapply(state_keep, function(x) make_node_id(cfg_row$compartment, x), character(1)),
        node_label = vapply(state_keep, function(x) format_node_label(cfg_row$compartment, cfg_row$display, x), character(1)),
        mp_name = state_keep,
        mp_display = state_keep,
        mp_index = seq_along(state_keep),
        mp_plot_order = seq_along(state_keep),
        sample_coverage_n = as.numeric(adj_bits$coverage_counts[state_keep]),
        denominator_samples = nrow(adj_bits$adjusted_scores),
        coverage_pct = 100 * as.numeric(adj_bits$coverage_counts[state_keep]) / nrow(adj_bits$adjusted_scores),
        coverage_pass = as.numeric(adj_bits$coverage_counts[state_keep]) > min_positive_samples,
        n_cells = nrow(score_df),
        plot_colour = cfg_row$plot_colour
      )

      return(list(
        compartment = cfg_row$compartment,
        display = cfg_row$display,
        atlas_celltype = cfg_row$atlas_celltype,
        subtype = cfg_row$subtype,
        score_path = NA_character_,
        mp_path = "Auto_final_states.rds",
        score_matrix = score_df,
        mp_names = state_keep,
        mp_display = setNames(state_keep, state_keep),
        gene_sets = cancer_state_ref_mode$gene_sets[state_keep],
        positive_rule = "assigned_state",
        cells_use = cells_use,
        sample_meta = sample_meta,
        adjusted_scores = adj_bits$adjusted_scores,
        node_stats = node_stats,
        cutoff_sensitivity = cutoff_tbl
      ))
    }

    mp_outs <- readRDS(cfg_row$mp_path)
    mp_genes <- filter_gene_sets_by_silhouette(mp_outs)
    score_df <- readRDS(cfg_row$score_path)
    keep_mps <- intersect(names(mp_genes), colnames(score_df))
    keep_mps <- if (cfg_row$compartment == "cancer") order_cancer_mps(mp_outs, keep_mps) else sort_mp_names(keep_mps)
    if (length(keep_mps) == 0) {
      stop("No silhouette-retained MPs found for ", cfg_row$compartment)
    }

    cells_use <- intersect(rownames(score_df), merged_cells_mode)
    if (!is.na(cfg_row$atlas_celltype) && "celltype_update" %in% colnames(merged_meta_mode)) {
      atlas_cells <- rownames(merged_meta_mode)[merged_meta_mode$celltype_update == cfg_row$atlas_celltype]
      cells_use <- intersect(cells_use, atlas_cells)
    }
    if (cfg_row$compartment == "cancer") {
      cells_use <- intersect(cells_use, malignant_cells_mode)
    }
    if (length(cells_use) == 0) {
      stop("No cells available for compartment ", cfg_row$compartment)
    }

    score_df <- score_df[cells_use, keep_mps, drop = FALSE]
    sample_meta <- make_sample_meta(merged_meta_mode[cells_use, , drop = FALSE])
    adj_bits <- calc_adjusted_scores(score_df, sample_meta, keep_mps, cutoff = ucell_cutoff)
    cutoff_tbl <- calc_cutoff_sensitivity(
      score_df = score_df,
      sample_meta = sample_meta,
      mp_names = keep_mps,
      cutoffs = cutoff_grid,
      compartment = cfg_row$compartment,
      display_name = cfg_row$display
    ) %>%
      mutate(
        celltype_order = unname(pair_order_lookup[cfg_row$display]),
        mp_plot_order = match(mp_name, keep_mps)
      )

    mp_display <- vapply(keep_mps, function(mp) format_mp_display(cfg_row$compartment, mp), character(1))
    node_stats <- tibble(
      compartment = cfg_row$compartment,
      celltype_display = cfg_row$display,
      celltype_order = unname(pair_order_lookup[cfg_row$display]),
      node_id = vapply(keep_mps, function(x) make_node_id(cfg_row$compartment, x), character(1)),
      node_label = vapply(keep_mps, function(mp) format_node_label(cfg_row$compartment, cfg_row$display, mp), character(1)),
      mp_name = keep_mps,
      mp_display = mp_display,
      mp_index = as.numeric(gsub("\\D", "", keep_mps)),
      mp_plot_order = seq_along(keep_mps),
      sample_coverage_n = as.numeric(adj_bits$coverage_counts[keep_mps]),
      denominator_samples = nrow(adj_bits$adjusted_scores),
      coverage_pct = 100 * as.numeric(adj_bits$coverage_counts[keep_mps]) / nrow(adj_bits$adjusted_scores),
      coverage_pass = as.numeric(adj_bits$coverage_counts[keep_mps]) > min_positive_samples,
      n_cells = nrow(score_df),
      plot_colour = cfg_row$plot_colour
    )

    list(
      compartment = cfg_row$compartment,
      display = cfg_row$display,
      atlas_celltype = cfg_row$atlas_celltype,
      subtype = cfg_row$subtype,
      score_path = cfg_row$score_path,
      mp_path = cfg_row$mp_path,
      score_matrix = NULL,
      mp_names = keep_mps,
      mp_display = setNames(mp_display, keep_mps),
      gene_sets = mp_genes[keep_mps],
      positive_rule = "ucell_cutoff",
      cells_use = cells_use,
      sample_meta = sample_meta,
      adjusted_scores = adj_bits$adjusted_scores,
      node_stats = node_stats,
      cutoff_sensitivity = cutoff_tbl
    )
  }

  ####################
  # Step 1 cache: compartment-level adjusted scores
  ####################
  step1_cache_mode <- load_or_build_cache(
    cache_filename = "Auto_celltype_step1_compartment_cache.rds",
    step_name = paste0("compartment adjusted scores [", mode_cfg$analysis_id, "]"),
    cache_dir_use = analysis_cache_dir,
    build_fun = function() {
      compartment_data <- setNames(
        lapply(seq_len(nrow(celltype_cfg)), function(i) build_compartment_data_mode(celltype_cfg[i, , drop = FALSE])),
        celltype_cfg$compartment
      )
      node_summary <- bind_rows(lapply(compartment_data, `[[`, "node_stats")) %>%
        arrange(celltype_order, mp_plot_order)
      cutoff_summary <- bind_rows(lapply(compartment_data, `[[`, "cutoff_sensitivity")) %>%
        arrange(celltype_order, mp_plot_order, cutoff)
      list(
        compartment_data = compartment_data,
        node_summary = node_summary,
        cutoff_summary = cutoff_summary
      )
    }
  )

  compartment_data <- step1_cache_mode$compartment_data
  node_summary <- step1_cache_mode$node_summary
  cutoff_summary <- step1_cache_mode$cutoff_summary

  # Mode-specific significance threshold
  mode_positive_sig_cutoff <- if (identical(mode_cfg$cancer_definition, "state")) 1.5 else positive_sig_cutoff

  ####################
  # Save adjusted-score tables and cutoff diagnostics
  ####################
  for (compartment_name in names(compartment_data)) {
    write.csv(
      prepare_adjusted_score_export(compartment_name, compartment_data[[compartment_name]]$adjusted_scores, compartment_data[[compartment_name]]$mp_names),
      file.path(analysis_out_dir, paste0("Auto_adjusted_scores_", compartment_name, ".csv")),
      row.names = FALSE
    )
  }
  write.csv(node_summary, file.path(analysis_out_dir, "Auto_celltype_node_summary.csv"), row.names = FALSE)
  write.csv(cutoff_summary, file.path(analysis_out_dir, "Auto_celltype_cutoff_sensitivity.csv"), row.names = FALSE)

  if (nrow(cutoff_summary) > 0) {
    cutoff_subtitle_text <- if (identical(mode_cfg$cancer_definition, "state")) {
      paste0(
        "Dashed line = active UCell cutoff (",
        ucell_cutoff,
        "); cancer-state labels are assignment-based and excluded"
      )
    } else {
      paste0("Dashed line = active cutoff (", ucell_cutoff, ")")
    }

    cutoff_count_plot <- cutoff_summary %>%
      group_by(celltype_display, cutoff) %>%
      summarise(n_nodes_passing = sum(coverage_pass), .groups = "drop") %>%
      ggplot(aes(x = cutoff, y = n_nodes_passing, color = celltype_display)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      geom_vline(xintercept = ucell_cutoff, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = plot_colour_lookup) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank()) +
      labs(
        title = paste0("Positive-node coverage across cutoff choices: ", mode_cfg$analysis_label),
        subtitle = cutoff_subtitle_text,
        x = "UCell positivity cutoff",
        y = paste0("Number of nodes with > ", min_positive_samples, " positive samples"),
        color = "Cell type"
      )

    cutoff_box_plot <- cutoff_summary %>%
      ggplot(aes(x = factor(cutoff), y = sample_coverage_n, fill = celltype_display)) +
      geom_boxplot(outlier.size = 0.7) +
      geom_hline(yintercept = min_positive_samples, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = plot_colour_lookup) +
      facet_wrap(~celltype_display, scales = "free_y", ncol = 4) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(), legend.position = "none") +
      labs(
        title = "Per-node sample coverage across cutoff choices",
        x = "UCell positivity cutoff",
        y = "Positive samples per node"
      )

    ggsave(
      file.path(analysis_out_dir, "Auto_celltype_cutoff_sensitivity.pdf"),
      cutoff_count_plot / cutoff_box_plot + plot_layout(heights = c(1, 1)),
      width = 14,
      height = 11
    )
  }

  ####################
  # Step 2 cache: pairwise correlations
  ####################
  step2_cache_mode <- load_or_build_cache(
    cache_filename = "Auto_celltype_step2_correlation_cache.rds",
    step_name = paste0("pairwise celltype correlations [", mode_cfg$analysis_id, "]"),
    cache_dir_use = analysis_cache_dir,
    build_fun = function() {
      pair_results <- list()
      pair_study_rows <- list()
      pair_grid <- combn(names(compartment_data), 2, simplify = FALSE)
      if (isTRUE(mode_cfg$include_within)) {
        pair_grid <- c(lapply(names(compartment_data), function(comp_name) c(comp_name, comp_name)), pair_grid)
      }

      for (pair_idx in seq_along(pair_grid)) {
        comp_a_name <- pair_grid[[pair_idx]][1]
        comp_b_name <- pair_grid[[pair_idx]][2]
        same_compartment <- identical(comp_a_name, comp_b_name)
        comp_a <- compartment_data[[comp_a_name]]
        comp_b <- compartment_data[[comp_b_name]]
        nodes_a <- comp_a$node_stats %>% filter(coverage_pass)
        nodes_b <- comp_b$node_stats %>% filter(coverage_pass)
        if (same_compartment && nrow(nodes_a) < 2) next
        if (!same_compartment && (nrow(nodes_a) == 0 || nrow(nodes_b) == 0)) next

        shared_samples <- if (same_compartment) {
          comp_a$adjusted_scores$sample
        } else {
          intersect(comp_a$adjusted_scores$sample, comp_b$adjusted_scores$sample)
        }
        if (length(shared_samples) == 0) next

        shared_info <- tibble(sample = shared_samples) %>%
          left_join(full_sample_lookup_mode, by = "sample") %>%
          mutate(study = ifelse(is.na(study) | study == "", derive_study(sample), study))

        study_df <- shared_info %>%
          count(study, name = "shared_samples") %>%
          mutate(
            compartment1 = comp_a_name,
            compartment2 = comp_b_name,
            celltype1_display = comp_a$display,
            celltype2_display = comp_b$display,
            pair_scope = ifelse(same_compartment, "within", "cross"),
            eligible = shared_samples >= min_shared_samples_per_study
          )

        eligible_studies <- study_df$study[study_df$eligible]
        eligible_samples <- shared_info$sample[shared_info$study %in% eligible_studies]
        pair_study_rows[[length(pair_study_rows) + 1]] <- study_df
        if (length(eligible_samples) < min_pair_samples) next

        score_a <- comp_a$adjusted_scores[match(eligible_samples, comp_a$adjusted_scores$sample), , drop = FALSE]
        score_b <- if (same_compartment) score_a else comp_b$adjusted_scores[match(eligible_samples, comp_b$adjusted_scores$sample), , drop = FALSE]
        keep_pair <- !is.na(score_a$sample) & !is.na(score_b$sample)
        score_a <- score_a[keep_pair, , drop = FALSE]
        score_b <- score_b[keep_pair, , drop = FALSE]
        if (nrow(score_a) < min_pair_samples || nrow(score_b) < min_pair_samples) next

        if (same_compartment) {
          feature_pairs <- combn(nodes_a$mp_name, 2, simplify = FALSE)
        } else {
          feature_grid <- expand.grid(mp_a = nodes_a$mp_name, mp_b = nodes_b$mp_name, stringsAsFactors = FALSE)
          feature_pairs <- split(feature_grid, seq_len(nrow(feature_grid)))
        }
        if (length(feature_pairs) == 0) next

        for (feature_idx in seq_along(feature_pairs)) {
          if (same_compartment) {
            mp_a <- feature_pairs[[feature_idx]][1]
            mp_b <- feature_pairs[[feature_idx]][2]
          } else {
            mp_a <- feature_pairs[[feature_idx]]$mp_a
            mp_b <- feature_pairs[[feature_idx]]$mp_b
          }

          pear_bits <- safe_cor_test(score_a[[mp_a]], score_b[[mp_b]], method = "pearson")
          spear_bits <- safe_cor_test(score_a[[mp_a]], score_b[[mp_b]], method = "spearman")
          co_positive_sample_n <- sum(score_a[[mp_a]] > 0 & score_b[[mp_b]] > 0, na.rm = TRUE)
          co_positive_sample_pct <- 100 * co_positive_sample_n / nrow(score_a)

          pair_results[[length(pair_results) + 1]] <- tibble(
            edge_id = paste(comp_a_name, mp_a, comp_b_name, mp_b, sep = "__"),
            compartment1 = comp_a_name,
            compartment2 = comp_b_name,
            celltype1_display = comp_a$display,
            celltype2_display = comp_b$display,
            pair_scope = ifelse(same_compartment, "within", "cross"),
            node1_id = make_node_id(comp_a_name, mp_a),
            node2_id = make_node_id(comp_b_name, mp_b),
            node1_label = comp_a$node_stats$node_label[match(mp_a, comp_a$node_stats$mp_name)],
            node2_label = comp_b$node_stats$node_label[match(mp_b, comp_b$node_stats$mp_name)],
            mp1_name = mp_a,
            mp2_name = mp_b,
            mp1_display = comp_a$node_stats$mp_display[match(mp_a, comp_a$node_stats$mp_name)],
            mp2_display = comp_b$node_stats$mp_display[match(mp_b, comp_b$node_stats$mp_name)],
            shared_sample_n = nrow(score_a),
            co_positive_sample_n = co_positive_sample_n,
            co_positive_sample_pct = co_positive_sample_pct,
            eligible_studies = paste(sort(eligible_studies), collapse = ";"),
            pearson_r = pear_bits$estimate,
            pearson_p = pear_bits$p.value,
            spearman_r = spear_bits$estimate,
            spearman_p = spear_bits$p.value
          )
        }
      }

      result_df <- bind_rows(pair_results)
      study_summary <- bind_rows(pair_study_rows)
      if (nrow(result_df) == 0) {
        stop("No pairwise correlations were generated for mode ", mode_cfg$analysis_id)
      }

      result_df <- result_df %>%
        mutate(
          spearman_significance = -log10(pmax(spearman_p, .Machine$double.xmin)),
          pearson_direction = case_when(
            is.na(pearson_r) ~ NA_character_,
            pearson_r > 0 ~ "positive",
            pearson_r < 0 ~ "negative",
            TRUE ~ "zero"
          ),
          spearman_sig = !is.na(spearman_p) & spearman_p < 0.05,
          celltype_pair = paste(celltype1_display, "vs", celltype2_display),
          edge_label = paste(node1_label, node2_label, sep = " <-> ")
        )

      list(
        result_df = result_df,
        study_summary = study_summary,
        positive_edges = result_df %>%
          filter(spearman_sig, !is.na(pearson_r), pearson_r > 0, spearman_significance >= mode_positive_sig_cutoff) %>%
          arrange(desc(spearman_significance), desc(pearson_r)),
        negative_edges = result_df %>%
          filter(spearman_sig, !is.na(pearson_r), pearson_r < 0, spearman_significance >= negative_sig_cutoff) %>%
          arrange(desc(spearman_significance), pearson_r)
      )
    }
  )

  result_df <- step2_cache_mode$result_df
  study_summary <- step2_cache_mode$study_summary
  
  # Re-filter edges from result_df in case thresholds changed since cache creation
  positive_edges <- result_df %>%
    filter(spearman_sig, !is.na(pearson_r), pearson_r > 0, spearman_significance >= mode_positive_sig_cutoff) %>%
    arrange(desc(spearman_significance), desc(pearson_r))
  negative_edges <- result_df %>%
    filter(spearman_sig, !is.na(pearson_r), pearson_r < 0, spearman_significance >= negative_sig_cutoff) %>%
    arrange(desc(spearman_significance), pearson_r)

  if (!"edge_label" %in% colnames(result_df)) {
    result_df <- result_df %>% mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
  }
  if (!"edge_label" %in% colnames(positive_edges) && nrow(positive_edges) > 0) {
    positive_edges <- positive_edges %>% mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
  }
  if (!"edge_label" %in% colnames(negative_edges) && nrow(negative_edges) > 0) {
    negative_edges <- negative_edges %>% mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
  }
  result_df <- add_co_positive_support(result_df, compartment_data)
  positive_edges <- add_co_positive_support(positive_edges, compartment_data)
  negative_edges <- add_co_positive_support(negative_edges, compartment_data)

  write.csv(study_summary, file.path(analysis_out_dir, "Auto_celltype_shared_sample_summary.csv"), row.names = FALSE)
  write.csv(result_df, file.path(analysis_out_dir, "Auto_celltype_correlations_all.csv"), row.names = FALSE)
  write.csv(positive_edges, file.path(analysis_out_dir, "Auto_celltype_correlations_positive.csv"), row.names = FALSE)
  write.csv(negative_edges, file.path(analysis_out_dir, "Auto_celltype_correlations_negative.csv"), row.names = FALSE)

  ####################
  # Correlation overview plots
  ####################
  pdf(file.path(analysis_out_dir, "Auto_celltype_correlation_bubble.pdf"), width = 16, height = 12, onefile = TRUE)
  for (focal in celltype_display_order) {
    focal_pairs <- get_expected_focal_pairs(focal, include_within = mode_cfg$include_within)
    pair_plots <- lapply(focal_pairs, function(pair_name) {
      pair_bits <- strsplit(pair_name, " vs ", fixed = TRUE)[[1]]
      partner <- get_partner_celltype(pair_bits, focal)
      pair_df <- result_df %>% filter(celltype_pair == pair_name)
      if (nrow(pair_df) > 0 && pair_bits[1] != focal) {
        pair_df <- pair_df %>%
          transmute(
            edge_id = edge_id,
            node1_label_plot = node2_label,
            node2_label_plot = node1_label,
            spearman_significance = spearman_significance,
            pearson_r = pearson_r,
            spearman_sig = spearman_sig
          ) %>%
          rename(node1_label = node1_label_plot, node2_label = node2_label_plot)
      }
      x_nodes <- node_summary %>% filter(celltype_display == focal) %>% arrange(celltype_order, mp_plot_order) %>% pull(node_label)
      y_nodes <- node_summary %>% filter(celltype_display == partner) %>% arrange(celltype_order, mp_plot_order) %>% pull(node_label)
      make_pair_bubble_plot(pair_df, paste(focal, partner, sep = " vs "), x_nodes, y_nodes)
    })

    ncol_page <- 3
    pad_n <- (ncol_page - (length(pair_plots) %% ncol_page)) %% ncol_page
    if (pad_n > 0) {
      pair_plots <- c(pair_plots, rep(list(patchwork::plot_spacer()), pad_n))
    }
    print(
      patchwork::wrap_plots(pair_plots, ncol = ncol_page, guides = "collect") +
        patchwork::plot_annotation(title = paste("Bubble plot:", focal))
    )
  }
  dev.off()

  ggsave(
    file.path(analysis_out_dir, "Auto_celltype_positive_network.pdf"),
    plot_network(
      edge_df = positive_edges,
      node_df = node_summary,
      title_text = paste("Positive correlations:", mode_cfg$analysis_label),
      subtitle_text = paste0(
        "Edges shown when Pearson > 0, Spearman p < 0.05, and -log10(p) >= ",
        mode_positive_sig_cutoff,
        "; ",
        describe_mode_positivity(mode_cfg)
      ),
      edge_low = "#F9D7CF",
      edge_high = "#8F0000"
    ),
    width = 15,
    height = 11
  )
  ggsave(
    file.path(analysis_out_dir, "Auto_celltype_negative_network.pdf"),
    plot_network(
      edge_df = negative_edges,
      node_df = node_summary,
      title_text = paste("Negative correlations:", mode_cfg$analysis_label),
      subtitle_text = "Edges shown when Pearson < 0 and Spearman p < 0.05",
      edge_low = "#DCEAF7",
      edge_high = "#08519C"
    ),
    width = 15,
    height = 11
  )

  ####################
  # Ligand-receptor annotation
  ####################
  lr_status <- tibble(
    analysis_id = mode_cfg$analysis_id,
    status = if (is.null(lr_ref_mode)) "missing" else "loaded",
    source = if (is.null(lr_ref_mode)) NA_character_ else lr_ref_mode$source,
    sheet = if (is.null(lr_ref_mode)) NA_character_ else lr_ref_mode$sheet,
    cancer_definition = mode_cfg$cancer_definition,
    include_within = mode_cfg$include_within,
    retained_evidence = paste(pair_evidence_allowed, collapse = "; "),
    n_pairs_all = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_all,
    n_pairs_retained = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_retained,
    n_pairs_excluded = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_excluded,
    cancer_state_marker_path = if (identical(mode_cfg$cancer_definition, "state")) cancer_state_ref_mode$marker_path else NA_character_,
    cancer_state_marker_top_n = if (identical(mode_cfg$cancer_definition, "state")) state_marker_top_n else NA_integer_
  )
  write.csv(lr_status, file.path(analysis_out_dir, "Auto_celltype_ligand_receptor_status.csv"), row.names = FALSE)

  score_cache_mode <- list()
  get_compartment_scores_mode <- function(compartment_name) {
    if (!is.null(score_cache_mode[[compartment_name]])) {
      return(score_cache_mode[[compartment_name]])
    }
    comp <- compartment_data[[compartment_name]]
    score_df <- if (!is.null(comp$score_matrix)) comp$score_matrix else readRDS(comp$score_path)
    keep_cells <- intersect(rownames(score_df), comp$cells_use)
    score_df <- score_df[keep_cells, comp$mp_names, drop = FALSE]
    score_cache_mode[[compartment_name]] <<- score_df
    score_df
  }
  get_positive_cells_mode <- function(compartment_name, mp_name, shared_samples) {
    comp <- compartment_data[[compartment_name]]
    score_df <- get_compartment_scores_mode(compartment_name)
    candidate_cells <- comp$sample_meta$cell[comp$sample_meta$sample %in% shared_samples]
    candidate_cells <- intersect(candidate_cells, rownames(score_df))
    if (length(candidate_cells) == 0) return(character(0))
    if (identical(comp$positive_rule, "assigned_state")) {
      candidate_cells[score_df[candidate_cells, mp_name] > 0]
    } else {
      candidate_cells[score_df[candidate_cells, mp_name] > ucell_cutoff]
    }
  }
  match_driver_target_lr_mode <- function(driver_top_genes, target_mp_genes, edge_row, driver_label, target_label) {
    if (length(driver_top_genes) == 0 || length(target_mp_genes) == 0 || is.null(lr_ref_mode)) {
      return(tibble())
    }
    driver_top_genes <- unique(toupper(driver_top_genes))
    target_mp_genes <- unique(toupper(target_mp_genes))
    bind_rows(
      lr_ref_mode$pairs %>%
        filter(ligand %in% driver_top_genes, receptor %in% target_mp_genes) %>%
        mutate(match_mode = "driver_top_ligand_target_mp_receptor"),
      lr_ref_mode$pairs %>%
        filter(receptor %in% driver_top_genes, ligand %in% target_mp_genes) %>%
        mutate(match_mode = "driver_top_receptor_target_mp_ligand")
    ) %>%
      mutate(
        edge_id = edge_row$edge_id,
        edge_label = edge_row$edge_label,
        direction = paste(driver_label, "->", target_label),
        driver_node = driver_label,
        target_node = target_label,
        compartment1 = edge_row$compartment1,
        compartment2 = edge_row$compartment2,
        celltype1_display = edge_row$celltype1_display,
        celltype2_display = edge_row$celltype2_display,
        node1_label = edge_row$node1_label,
        node2_label = edge_row$node2_label,
        mp1_name = edge_row$mp1_name,
        mp2_name = edge_row$mp2_name,
        mp1_display = edge_row$mp1_display,
        mp2_display = edge_row$mp2_display,
        celltype_pair = edge_row$celltype_pair,
        pair_scope = edge_row$pair_scope,
        pearson_r = edge_row$pearson_r,
        spearman_r = edge_row$spearman_r,
        spearman_significance = edge_row$spearman_significance
      )
  }

  step3_cache_mode <- load_or_build_cache(
    cache_filename = "Auto_celltype_step3_lr_cache.rds",
    step_name = paste0("ligand-receptor annotation [", mode_cfg$analysis_id, "]"),
    cache_dir_use = analysis_cache_dir,
    build_fun = function() {
      lr_pairs_df <- tibble()
      lr_edge_summary <- tibble()
      positive_edges_lr <- positive_edges
      if (is.null(lr_ref_mode) || nrow(positive_edges) == 0) {
        return(list(lr_pairs_df = lr_pairs_df, lr_edge_summary = lr_edge_summary, positive_edges_lr = positive_edges_lr))
      }

      lr_rows <- list()
      for (edge_idx in seq_len(nrow(positive_edges))) {
        edge_row <- positive_edges[edge_idx, , drop = FALSE]
        shared_studies <- strsplit(edge_row$eligible_studies, ";", fixed = TRUE)[[1]]
        shared_studies <- shared_studies[nzchar(shared_studies)]
        comp1_scores <- compartment_data[[edge_row$compartment1]]$adjusted_scores
        comp2_scores <- compartment_data[[edge_row$compartment2]]$adjusted_scores
        eligible_samples <- intersect(
          comp1_scores$sample[comp1_scores$study %in% shared_studies],
          comp2_scores$sample[comp2_scores$study %in% shared_studies]
        )
        if (length(eligible_samples) == 0) next

        node1_cells <- get_positive_cells_mode(edge_row$compartment1, edge_row$mp1_name, eligible_samples)
        node2_cells <- get_positive_cells_mode(edge_row$compartment2, edge_row$mp2_name, eligible_samples)
        if (length(node1_cells) == 0 || length(node2_cells) == 0) next

        node1_top <- extract_top_genes(merged_expr_mode_all, node1_cells, top_n = top_ranked_genes)
        node2_top <- extract_top_genes(merged_expr_mode_all, node2_cells, top_n = top_ranked_genes)
        node1_genes <- toupper(unique(compartment_data[[edge_row$compartment1]]$gene_sets[[edge_row$mp1_name]]))
        node2_genes <- toupper(unique(compartment_data[[edge_row$compartment2]]$gene_sets[[edge_row$mp2_name]]))

        lr_rows[[length(lr_rows) + 1]] <- bind_rows(
          match_driver_target_lr_mode(node1_top, node2_genes, edge_row, edge_row$node1_label, edge_row$node2_label) %>%
            mutate(driver_positive_cells = length(node1_cells), target_positive_cells = length(node2_cells)),
          match_driver_target_lr_mode(node2_top, node1_genes, edge_row, edge_row$node2_label, edge_row$node1_label) %>%
            mutate(driver_positive_cells = length(node2_cells), target_positive_cells = length(node1_cells))
        )
      }

      lr_pairs_df <- bind_rows(lr_rows)
      if (nrow(lr_pairs_df) > 0) {
        lr_pairs_df <- lr_pairs_df %>%
          distinct(edge_id, direction, ligand, receptor, match_mode, .keep_all = TRUE)
        lr_edge_summary <- lr_pairs_df %>%
          group_by(edge_id, edge_label, celltype_pair, pair_scope) %>%
          summarise(n_lr_rows = n(), n_lr_directions = n_distinct(direction), .groups = "drop") %>%
          left_join(summarise_edge_lr_labels(lr_pairs_df), by = "edge_id") %>%
          arrange(desc(n_literature_supported), desc(n_lr_pairs), desc(max_pair_support_n), desc(n_lr_rows))
        positive_edges_lr <- positive_edges %>%
          left_join(lr_edge_summary, by = "edge_id") %>%
          arrange(desc(spearman_significance), desc(n_lr_pairs), desc(n_literature_supported))
      }
      list(lr_pairs_df = lr_pairs_df, lr_edge_summary = lr_edge_summary, positive_edges_lr = positive_edges_lr)
    }
  )

  lr_pairs_df <- step3_cache_mode$lr_pairs_df
  lr_edge_summary <- step3_cache_mode$lr_edge_summary
  
  # Re-sync positive_edges_lr with potentially updated positive_edges (from threshold changes)
  # This ensures new edges (sig 1.5-4.0) appear as dots even if Step 3 cache only had sig > 4.0.
  positive_edges_lr <- positive_edges %>%
    left_join(lr_edge_summary, by = "edge_id") %>%
    arrange(desc(spearman_significance), desc(n_lr_pairs), desc(n_literature_supported))

  positive_edges_lr <- add_co_positive_support(positive_edges_lr, compartment_data)

  if (nrow(lr_pairs_df) > 0) {
    write.csv(lr_pairs_df, file.path(analysis_out_dir, "Auto_celltype_ligand_receptor_pairs.csv"), row.names = FALSE)
    write.csv(lr_edge_summary, file.path(analysis_out_dir, "Auto_celltype_ligand_receptor_edge_summary.csv"), row.names = FALSE)
    write.csv(positive_edges_lr, file.path(analysis_out_dir, "Auto_celltype_positive_edges_lr_annotated.csv"), row.names = FALSE)
  }

  ####################
  # Focal-celltype interaction dotmap. Each PDF page focuses on one
  # cell type, with its MPs/states as rows and partner MPs/states as
  # column groups.
  ####################
  interaction_dotmap_df <- bind_rows(lapply(
    celltype_display_order,
    function(focal) prepare_focal_interaction_dotmap_data(
      positive_edges_lr,
      node_summary,
      focal = focal,
      include_within = mode_cfg$include_within
    )
  ))
  write.csv(
    interaction_dotmap_df,
    file.path(analysis_out_dir, "Auto_celltype_interaction_dotmap_data.csv"),
    row.names = FALSE
  )
  write_focal_interaction_dotmap(
    edge_df = positive_edges_lr,
    node_df = node_summary,
    file.path(analysis_out_dir, "Auto_celltype_interaction_dotmap.pdf"),
    include_within = mode_cfg$include_within
  )

  excel_output_paths <- character(0)
  if (nrow(lr_edge_summary) > 0) {
    edge_annotations <- positive_edges_lr %>%
      filter(!is.na(n_lr_pairs), n_lr_pairs > 0) %>%
      slice_head(n = annotated_positive_edges_n) %>%
      transmute(edge_id, edge_note = top_lr_label)

    ggsave(
      file.path(analysis_out_dir, "Auto_celltype_positive_network_lr_annotated.pdf"),
      plot_network(
        edge_df = positive_edges_lr,
        node_df = node_summary,
        title_text = paste("Positive correlations with LR support:", mode_cfg$analysis_label),
        subtitle_text = paste0(
          "Only literature-supported and putative LR pairs retained; ",
          describe_mode_positivity(mode_cfg)
        ),
        edge_low = "#F9D7CF",
        edge_high = "#8F0000",
        annotation_df = edge_annotations
      ),
      width = 16,
      height = 12
    )

    lr_edge_barplot <- lr_edge_summary %>%
      slice_head(n = 25) %>%
      mutate(edge_label = factor(edge_label, levels = rev(edge_label))) %>%
      ggplot(aes(x = n_lr_pairs, y = edge_label, fill = n_literature_supported)) +
      geom_col() +
      scale_fill_gradient(low = "#FEE0D2", high = "#A50F15") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16)
      ) +
      labs(
        title = "Top positive edges by LR support",
        subtitle = "Counts shown after removing all excluded LR pairs",
        x = "Number of candidate LR pairs",
        y = NULL,
        fill = "Literature\nsupported"
      )

    lr_pair_barplot <- lr_pairs_df %>%
      mutate(pair_label = ifelse(!is.na(pair_name) & pair_name != "", pair_name, paste(ligand, receptor, sep = " -> "))) %>%
      count(pair_label, pair_evidence, sort = TRUE, name = "n_edges_supported") %>%
      slice_head(n = 25) %>%
      mutate(pair_label = factor(pair_label, levels = rev(unique(pair_label)))) %>%
      ggplot(aes(x = n_edges_supported, y = pair_label, fill = pair_evidence)) +
      geom_col(position = "stack") +
      scale_fill_manual(values = c("literature supported" = "#66A61E", "putative" = "#E6AB02")) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16)
      ) +
      labs(
        title = "Most recurrent retained LR pairs",
        x = "Positive edges supported",
        y = NULL,
        fill = "Pair evidence"
      )

    ggsave(
      file.path(analysis_out_dir, "Auto_celltype_ligand_receptor_summary.pdf"),
      lr_edge_barplot / lr_pair_barplot + plot_layout(heights = c(1, 1)),
      width = 12,
      height = 12
    )

    excel_output_paths <- write_lr_workbooks(
      lr_pairs_df = lr_pairs_df,
      lr_edge_summary = lr_edge_summary,
      out_dir = analysis_out_dir,
      include_within = mode_cfg$include_within
    )
  }

  cancer_edge_df <- result_df %>%
    filter(compartment1 == "cancer" | compartment2 == "cancer")
  cancer_positive_edge_df <- cancer_edge_df %>%
    filter(spearman_sig, !is.na(pearson_r), pearson_r > 0)

  ####################
  # Cancer-TME interaction summary Excel (TME-centric with Gene Lists)
  ####################
  if (nrow(cancer_positive_edge_df) > 0) {
    message("Generating TME-centric Cancer-TME interaction Excel")
    
    # 1. Prepare ordering information
    celltype_preferred_order <- c("fibroblast", "endothelial", "cd8", "cd4", "macrophage", "nk", "plasma")
    cancer_mp_order <- compartment_data$cancer$node_stats$mp_name
    
    # 2. Add co-positive support if not already present
    if (!"co_positive_sample_pct" %in% colnames(cancer_positive_edge_df)) {
        cancer_positive_edge_df <- add_co_positive_support(cancer_positive_edge_df, compartment_data)
    }

    # 3. Process interactions and group by TME MP
    tme_interaction_groups <- cancer_positive_edge_df %>%
      mutate(
        cancer_mp = ifelse(compartment1 == "cancer", mp1_name, mp2_name),
        cancer_display = ifelse(compartment1 == "cancer", mp1_display, mp2_display),
        tme_compartment = ifelse(compartment1 == "cancer", compartment2, compartment1),
        tme_mp = ifelse(compartment1 == "cancer", mp2_name, mp1_name),
        tme_display = ifelse(compartment1 == "cancer", mp2_display, mp1_display),
        log10p = spearman_significance,
        support_pct = co_positive_sample_pct,
        cancer_order = match(cancer_mp, cancer_mp_order)
      ) %>%
      group_by(tme_compartment, tme_mp, tme_display) %>%
      nest() %>%
      mutate(
        ct_order = match(tme_compartment, celltype_preferred_order),
        mp_num = as.numeric(gsub("\\D", "", tme_mp))
      ) %>%
      arrange(ct_order, mp_num)

    # 4. Create Excel Workbook
    interaction_wb <- createWorkbook()
    addWorksheet(interaction_wb, "TME-Cancer Interactions")

    # Styles
    tmeHeaderStyle <- createStyle(fontSize = 12, fontColour = "#FFFFFF", fgFill = "#2F5597", halign = "center", textDecoration = "bold")
    cancerNameStyle <- createStyle(fontSize = 10, fgFill = "#E9E9E9", halign = "left", textDecoration = "bold")
    cancerStatStyle <- createStyle(fontSize = 9, halign = "left", textDecoration = "italic")
    geneStyle <- createStyle(fontSize = 10, halign = "left")

    # 5. Fill columns (One per TME MP)
    curr_col <- 1
    for (i in seq_len(nrow(tme_interaction_groups))) {
      group <- tme_interaction_groups[i, ]
      # Order partners by Cancer MP tree order
      partners <- group$data[[1]] %>% arrange(cancer_order)
      
      # Row 1: TME MP Name
      writeData(interaction_wb, sheet = 1, x = paste0(group$tme_compartment, " ", group$tme_display), startRow = 1, startCol = curr_col)
      addStyle(interaction_wb, sheet = 1, style = tmeHeaderStyle, rows = 1, cols = curr_col)
      
      # Rows 2+: Cancer MPs
      curr_row <- 2
      for (j in seq_len(nrow(partners))) {
        p <- partners[j, ]
        # Cancer MP Name
        writeData(interaction_wb, sheet = 1, x = p$cancer_display, startRow = curr_row, startCol = curr_col)
        addStyle(interaction_wb, sheet = 1, style = cancerNameStyle, rows = curr_row, cols = curr_col)
        curr_row <- curr_row + 1
        
        # Stats & Support
        stat_text <- sprintf("-log10(p): %.2f, Support: %.1f%%", p$log10p, p$support_pct)
        writeData(interaction_wb, sheet = 1, x = stat_text, startRow = curr_row, startCol = curr_col)
        addStyle(interaction_wb, sheet = 1, style = cancerStatStyle, rows = curr_row, cols = curr_col)
        curr_row <- curr_row + 1
      }
      
      # TME-MP Gene List
      genes <- compartment_data[[group$tme_compartment]]$gene_sets[[group$tme_mp]]
      if (is.null(genes) || length(genes) == 0) {
          # Fallback to original RDS if cache is incomplete
          mp_path_orig <- file.path("ref_outs", compartment_data[[group$tme_compartment]]$mp_path)
          if (file.exists(mp_path_orig)) {
              mp_outs_orig <- readRDS(mp_path_orig)
              if ("metaprograms.genes" %in% names(mp_outs_orig)) {
                  genes <- mp_outs_orig$metaprograms.genes[[group$tme_mp]]
              } else if ("programs.genes" %in% names(mp_outs_orig)) {
                  genes <- mp_outs_orig$programs.genes[[group$tme_mp]]
              }
          }
      }
      if (is.null(genes)) genes <- "No genes found"
      
      # Write genes (one per row)
      writeData(interaction_wb, sheet = 1, x = genes, startRow = curr_row, startCol = curr_col)
      addStyle(interaction_wb, sheet = 1, style = geneStyle, rows = curr_row:(curr_row + length(genes) - 1), cols = curr_col)
      
      setColWidths(interaction_wb, sheet = 1, cols = curr_col, widths = 30)
      curr_col <- curr_col + 1
    }

    interaction_excel_path <- file.path(analysis_out_dir, "Auto_cancer_tme_interactions_TME_centric.xlsx")
    saveWorkbook(interaction_wb, interaction_excel_path, overwrite = TRUE)
    excel_output_paths <- c(excel_output_paths, interaction_excel_path)
  }

  mode_summary <- tibble(
    analysis_id = mode_cfg$analysis_id,
    analysis_label = mode_cfg$analysis_label,
    output_dir = analysis_out_dir,
    cancer_definition = mode_cfg$cancer_definition,
    include_within = mode_cfg$include_within,
    cache_version = cache_version,
    ucell_cutoff = ucell_cutoff,
    n_compartments = length(compartment_data),
    n_nodes_total = nrow(node_summary),
    n_nodes_passing_coverage = sum(node_summary$coverage_pass),
    n_pairwise_tests = nrow(result_df),
    n_cancer_pairwise_tests = nrow(cancer_edge_df),
    n_positive_edges = nrow(positive_edges),
    n_negative_edges = nrow(negative_edges),
    n_cancer_positive_edges_p05 = nrow(cancer_positive_edge_df),
    n_cancer_positive_edges_sig = sum(cancer_positive_edge_df$spearman_significance >= mode_positive_sig_cutoff, na.rm = TRUE),
    max_cancer_positive_spearman_significance = if (nrow(cancer_positive_edge_df) == 0) NA_real_ else max(cancer_positive_edge_df$spearman_significance, na.rm = TRUE),
    n_lr_rows = nrow(lr_pairs_df),
    n_lr_supported_edges = if (nrow(lr_edge_summary) == 0) 0 else sum(lr_edge_summary$n_lr_pairs > 0),
    n_lr_pairs_all_reference = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_all,
    n_lr_pairs_retained_reference = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_retained,
    n_lr_pairs_excluded_reference = if (is.null(lr_ref_mode)) NA_integer_ else lr_ref_mode$n_pairs_excluded,
    cancer_positive_rule = if (identical(mode_cfg$cancer_definition, "state")) "assigned_state_label" else paste0("ucell_gt_", ucell_cutoff),
    excel_outputs = if (length(excel_output_paths) == 0) "" else paste(excel_output_paths, collapse = "; ")
  )
  write.csv(mode_summary, file.path(analysis_out_dir, "Auto_celltype_mode_summary.csv"), row.names = FALSE)
  invisible(gc())
  message("Saved outputs to: ", file.path(getwd(), analysis_out_dir))
  mode_summary
}

mode_summaries <- bind_rows(lapply(seq_len(nrow(analysis_modes)), function(i) run_analysis_mode(analysis_modes[i, , drop = FALSE])))
write.csv(mode_summaries, file.path(summary_dir, "Auto_mp_cross_celltype_correlations_summary.csv"), row.names = FALSE)
message("Saved outputs for all modes under: ", file.path(getwd(), out_dir))

if (FALSE) {

####################
# Load the atlas and malignancy metadata
####################

message("Loading complete atlas and epithelial malignancy metadata")

merged_obj <- readRDS("EAC_Ref_merged_strict.rds")
merged_meta <- merged_obj@meta.data
merged_cells <- rownames(merged_meta)
epi_meta <- readRDS("meta_full_epi.rds")

malignant_cells <- rownames(epi_meta)[epi_meta$malignancy %in% c("malignant_level_1", "malignant_level_2")]
full_sample_lookup <- make_sample_meta(merged_meta) %>%
  distinct(sample, study)

####################
# Step 1 cache: compartment-level adjusted scores
####################

build_compartment_data <- function(cfg_row) {
  message("Preparing compartment: ", cfg_row$compartment)

  mp_outs <- readRDS(cfg_row$mp_path)
  mp_genes <- filter_gene_sets_by_silhouette(mp_outs)
  score_df <- readRDS(cfg_row$score_path)

  keep_mps <- intersect(names(mp_genes), colnames(score_df))
  keep_mps <- if (cfg_row$compartment == "cancer") {
    order_cancer_mps(mp_outs, keep_mps)
  } else {
    sort_mp_names(keep_mps)
  }
  if (length(keep_mps) == 0) {
    stop("No silhouette-retained MPs found for ", cfg_row$compartment)
  }

  cells_use <- intersect(rownames(score_df), merged_cells)

  if (!is.na(cfg_row$atlas_celltype) && "celltype_update" %in% colnames(merged_meta)) {
    atlas_cells <- rownames(merged_meta)[merged_meta$celltype_update == cfg_row$atlas_celltype]
    cells_use <- intersect(cells_use, atlas_cells)
  }

  if (cfg_row$compartment == "cancer") {
    cells_use <- intersect(cells_use, malignant_cells)
  }

  if (length(cells_use) == 0) {
    stop("No cells available for compartment ", cfg_row$compartment)
  }

  score_df <- score_df[cells_use, keep_mps, drop = FALSE]
  sample_meta <- make_sample_meta(merged_meta[cells_use, , drop = FALSE])
  adj_bits <- calc_adjusted_scores(score_df, sample_meta, keep_mps, cutoff = ucell_cutoff)
  cutoff_tbl <- calc_cutoff_sensitivity(
    score_df = score_df,
    sample_meta = sample_meta,
    mp_names = keep_mps,
    cutoffs = cutoff_grid,
    compartment = cfg_row$compartment,
    display_name = cfg_row$display
  ) %>%
    mutate(
      celltype_order = unname(pair_order_lookup[cfg_row$display]),
      mp_plot_order = match(mp_name, keep_mps)
    )

  mp_display <- vapply(keep_mps, function(mp) format_mp_display(cfg_row$compartment, mp), character(1))
  node_label <- vapply(keep_mps, function(mp) format_node_label(cfg_row$compartment, cfg_row$display, mp), character(1))

  node_stats <- tibble(
    compartment = cfg_row$compartment,
    celltype_display = cfg_row$display,
    celltype_order = unname(pair_order_lookup[cfg_row$display]),
    node_id = paste0(cfg_row$compartment, "__", keep_mps),
    node_label = node_label,
    mp_name = keep_mps,
    mp_display = mp_display,
    mp_index = as.numeric(gsub("\\D", "", keep_mps)),
    mp_plot_order = seq_along(keep_mps),
    sample_coverage_n = as.numeric(adj_bits$coverage_counts[keep_mps]),
    denominator_samples = nrow(adj_bits$adjusted_scores),
    coverage_pct = 100 * as.numeric(adj_bits$coverage_counts[keep_mps]) / nrow(adj_bits$adjusted_scores),
    coverage_pass = as.numeric(adj_bits$coverage_counts[keep_mps]) > min_positive_samples,
    n_cells = nrow(score_df),
    plot_colour = cfg_row$plot_colour
  )

  list(
    compartment = cfg_row$compartment,
    display = cfg_row$display,
    atlas_celltype = cfg_row$atlas_celltype,
    subtype = cfg_row$subtype,
    score_path = cfg_row$score_path,
    mp_path = cfg_row$mp_path,
    mp_names = keep_mps,
    mp_display = setNames(mp_display, keep_mps),
    gene_sets = mp_genes[keep_mps],
    cells_use = cells_use,
    sample_meta = sample_meta,
    adjusted_scores = adj_bits$adjusted_scores,
    node_stats = node_stats,
    cutoff_sensitivity = cutoff_tbl
  )
}

step1_cache <- load_or_build_cache(
  cache_filename = "Auto_cross_celltype_step1_compartment_cache.rds",
  step_name = "compartment adjusted scores",
  build_fun = function() {
    compartment_data <- setNames(
      lapply(seq_len(nrow(celltype_cfg)), function(i) build_compartment_data(celltype_cfg[i, , drop = FALSE])),
      celltype_cfg$compartment
    )

    node_summary <- bind_rows(lapply(compartment_data, `[[`, "node_stats")) %>%
      arrange(celltype_order, mp_plot_order)
    cutoff_summary <- bind_rows(lapply(compartment_data, `[[`, "cutoff_sensitivity")) %>%
      arrange(celltype_order, mp_plot_order, cutoff)

    list(
      compartment_data = compartment_data,
      node_summary = node_summary,
      cutoff_summary = cutoff_summary
    )
  }
)

compartment_data <- step1_cache$compartment_data
node_summary <- step1_cache$node_summary
cutoff_summary <- step1_cache$cutoff_summary

####################
# Save adjusted-score tables and cutoff diagnostics
####################

message("Saving adjusted-score and threshold diagnostics")

for (compartment_name in names(compartment_data)) {
  adj_scores <- compartment_data[[compartment_name]]$adjusted_scores
  mp_names <- compartment_data[[compartment_name]]$mp_names
  write.csv(
    prepare_adjusted_score_export(compartment_name, adj_scores, mp_names),
    file.path(out_dir, paste0("Auto_adjusted_scores_", compartment_name, ".csv")),
    row.names = FALSE
  )
}

write.csv(
  node_summary,
  file.path(out_dir, "Auto_cross_celltype_node_summary.csv"),
  row.names = FALSE
)
write.csv(
  cutoff_summary,
  file.path(out_dir, "Auto_cross_celltype_cutoff_sensitivity.csv"),
  row.names = FALSE
)

cutoff_count_plot <- cutoff_summary %>%
  group_by(celltype_display, cutoff) %>%
  summarise(n_mps_passing = sum(coverage_pass), .groups = "drop") %>%
  ggplot(aes(x = cutoff, y = n_mps_passing, color = celltype_display)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = ucell_cutoff, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = plot_colour_lookup) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 17)
  ) +
  labs(
    title = "MP retention is highly sensitive to the UCell positivity cutoff",
    subtitle = paste0("Dashed line = active cutoff (", ucell_cutoff, "); 0.5 is sparse for these UCell score ranges"),
    x = "UCell positivity cutoff",
    y = paste0("Number of MPs with > ", min_positive_samples, " positive samples"),
    color = "Cell type"
  )

cutoff_box_plot <- cutoff_summary %>%
  ggplot(aes(x = factor(cutoff), y = sample_coverage_n, fill = celltype_display)) +
  geom_boxplot(outlier.size = 0.7) +
  geom_hline(yintercept = min_positive_samples, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = plot_colour_lookup) +
  facet_wrap(~celltype_display, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  ) +
  labs(
    title = "Per-MP sample coverage across cutoff choices",
    x = "UCell positivity cutoff",
    y = "Positive samples per MP"
  )

ggsave(
  file.path(out_dir, "Auto_cross_celltype_cutoff_sensitivity.pdf"),
  cutoff_count_plot / cutoff_box_plot + plot_layout(heights = c(1, 1)),
  width = 14,
  height = 11
)

####################
# Step 2 cache: pairwise cross-celltype correlations
####################

step2_cache <- load_or_build_cache(
  cache_filename = "Auto_cross_celltype_step2_correlation_cache.rds",
  step_name = "pairwise cross-celltype correlations",
  build_fun = function() {
    pair_results <- list()
    pair_study_rows <- list()
    pair_grid <- combn(names(compartment_data), 2, simplify = FALSE)

    for (pair_idx in seq_along(pair_grid)) {
      comp_a_name <- pair_grid[[pair_idx]][1]
      comp_b_name <- pair_grid[[pair_idx]][2]
      comp_a <- compartment_data[[comp_a_name]]
      comp_b <- compartment_data[[comp_b_name]]

      nodes_a <- comp_a$node_stats %>% filter(coverage_pass)
      nodes_b <- comp_b$node_stats %>% filter(coverage_pass)
      if (nrow(nodes_a) == 0 || nrow(nodes_b) == 0) next

      shared_samples <- intersect(comp_a$adjusted_scores$sample, comp_b$adjusted_scores$sample)
      if (length(shared_samples) == 0) next

      shared_info <- tibble(sample = shared_samples) %>%
        left_join(full_sample_lookup, by = "sample") %>%
        mutate(study = ifelse(is.na(study) | study == "", derive_study(sample), study))

      study_df <- shared_info %>%
        count(study, name = "shared_samples") %>%
        mutate(
          compartment1 = comp_a_name,
          compartment2 = comp_b_name,
          celltype1_display = comp_a$display,
          celltype2_display = comp_b$display,
          eligible = shared_samples >= min_shared_samples_per_study
        )

      eligible_studies <- study_df$study[study_df$eligible]
      eligible_samples <- shared_info$sample[shared_info$study %in% eligible_studies]
      pair_study_rows[[length(pair_study_rows) + 1]] <- study_df

      if (length(eligible_samples) < min_pair_samples) next

      score_a <- comp_a$adjusted_scores[match(eligible_samples, comp_a$adjusted_scores$sample), , drop = FALSE]
      score_b <- comp_b$adjusted_scores[match(eligible_samples, comp_b$adjusted_scores$sample), , drop = FALSE]
      keep_pair <- !is.na(score_a$sample) & !is.na(score_b$sample)
      score_a <- score_a[keep_pair, , drop = FALSE]
      score_b <- score_b[keep_pair, , drop = FALSE]

      if (nrow(score_a) < min_pair_samples || nrow(score_b) < min_pair_samples) next

      for (mp_a in nodes_a$mp_name) {
        x <- score_a[[mp_a]]
        for (mp_b in nodes_b$mp_name) {
          y <- score_b[[mp_b]]
          pear_bits <- safe_cor_test(x, y, method = "pearson")
          spear_bits <- safe_cor_test(x, y, method = "spearman")

          pair_results[[length(pair_results) + 1]] <- tibble(
            edge_id = paste(comp_a_name, mp_a, comp_b_name, mp_b, sep = "__"),
            compartment1 = comp_a_name,
            compartment2 = comp_b_name,
            celltype1_display = comp_a$display,
            celltype2_display = comp_b$display,
            node1_id = paste0(comp_a_name, "__", mp_a),
            node2_id = paste0(comp_b_name, "__", mp_b),
            node1_label = format_node_label(comp_a_name, comp_a$display, mp_a),
            node2_label = format_node_label(comp_b_name, comp_b$display, mp_b),
            mp1_name = mp_a,
            mp2_name = mp_b,
            mp1_display = format_mp_display(comp_a_name, mp_a),
            mp2_display = format_mp_display(comp_b_name, mp_b),
            shared_sample_n = nrow(score_a),
            eligible_studies = paste(sort(eligible_studies), collapse = ";"),
            pearson_r = pear_bits$estimate,
            pearson_p = pear_bits$p.value,
            spearman_r = spear_bits$estimate,
            spearman_p = spear_bits$p.value
          )
        }
      }
    }

    result_df <- bind_rows(pair_results)
    study_summary <- bind_rows(pair_study_rows)

    if (nrow(result_df) == 0) {
      stop("No pairwise cross-celltype correlations were generated")
    }

    result_df <- result_df %>%
      mutate(
        spearman_significance = -log10(pmax(spearman_p, .Machine$double.xmin)),
        pearson_direction = case_when(
          is.na(pearson_r) ~ NA_character_,
          pearson_r > 0 ~ "positive",
          pearson_r < 0 ~ "negative",
          TRUE ~ "zero"
        ),
        spearman_sig = !is.na(spearman_p) & spearman_p < 0.05,
        celltype_pair = paste(celltype1_display, "vs", celltype2_display)
      )

    ####################
    # Add a stable edge label once at the correlation-table stage so every
    # downstream LR/export step can rely on the same field.
    ####################
    result_df <- result_df %>%
      mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))

    # Use the threshold defined in the mode scope
    positive_edges <- result_df %>%
      filter(
        spearman_sig,
        !is.na(pearson_r),
        pearson_r > 0,
        spearman_significance >= mode_positive_sig_cutoff
      ) %>%
      arrange(desc(spearman_significance), desc(pearson_r))

    negative_edges <- result_df %>%
      filter(
        spearman_sig,
        !is.na(pearson_r),
        pearson_r < 0,
        spearman_significance >= negative_sig_cutoff
      ) %>%
      arrange(desc(spearman_significance), pearson_r)

    list(
      result_df = result_df,
      study_summary = study_summary,
      positive_edges = positive_edges,
      negative_edges = negative_edges
    )
  }
)

result_df <- step2_cache$result_df
study_summary <- step2_cache$study_summary
positive_edges <- step2_cache$positive_edges
negative_edges <- step2_cache$negative_edges

####################
# Backfill edge_label for any cached step-2 objects that were created before
# the label column was added, so LR annotation remains cache-safe.
####################
if (!"edge_label" %in% colnames(result_df)) {
  result_df <- result_df %>%
    mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
}
if (!"edge_label" %in% colnames(positive_edges) && nrow(positive_edges) > 0) {
  positive_edges <- positive_edges %>%
    mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
}
if (!"edge_label" %in% colnames(negative_edges) && nrow(negative_edges) > 0) {
  negative_edges <- negative_edges %>%
    mutate(edge_label = paste(node1_label, node2_label, sep = " <-> "))
}

write.csv(
  study_summary,
  file.path(out_dir, "Auto_cross_celltype_shared_sample_summary.csv"),
  row.names = FALSE
)
write.csv(
  result_df,
  file.path(out_dir, "Auto_cross_celltype_correlations_all.csv"),
  row.names = FALSE
)
write.csv(
  positive_edges,
  file.path(out_dir, "Auto_cross_celltype_correlations_positive.csv"),
  row.names = FALSE
)
write.csv(
  negative_edges,
  file.path(out_dir, "Auto_cross_celltype_correlations_negative.csv"),
  row.names = FALSE
)

####################
# Correlation overview plots
####################

message("Saving cross-celltype correlation overview plots")

####################
# Split the bubble overview into focal-celltype pages and render each
# comparison as a fixed-size panel so partially filled pages do not stretch.
####################
pdf(file.path(out_dir, "Auto_cross_celltype_correlation_bubble.pdf"), width = 16, height = 12, onefile = TRUE)
for (focal in celltype_display_order) {
  focal_pairs <- get_all_focal_pairs(focal, apply(combn(celltype_display_order, 2), 2, paste, collapse = " vs "))

  pair_plots <- lapply(focal_pairs, function(pair_name) {
    pair_bits <- strsplit(pair_name, " vs ", fixed = TRUE)[[1]]
    pair_df <- result_df %>%
      filter(celltype_pair == pair_name)
    partner <- setdiff(pair_bits, focal)
    display_pair_name <- paste(focal, partner, sep = " vs ")

    if (nrow(pair_df) > 0 && pair_bits[1] != focal) {
      pair_df <- pair_df %>%
        transmute(
          edge_id = edge_id,
          node1_label_plot = node2_label,
          node2_label_plot = node1_label,
          spearman_significance = spearman_significance,
          pearson_r = pearson_r,
          spearman_sig = spearman_sig
        ) %>%
        rename(
          node1_label = node1_label_plot,
          node2_label = node2_label_plot
        )
    }

    x_nodes <- node_summary %>%
      filter(celltype_display == focal) %>%
      arrange(celltype_order, mp_plot_order) %>%
      pull(node_label)
    y_nodes <- node_summary %>%
      filter(celltype_display == partner) %>%
      arrange(celltype_order, mp_plot_order) %>%
      pull(node_label)

    make_pair_bubble_plot(pair_df, display_pair_name, x_nodes, y_nodes)
  })

  ncol_page <- 3
  pad_n <- (ncol_page - (length(pair_plots) %% ncol_page)) %% ncol_page
  if (pad_n > 0) {
    pair_plots <- c(pair_plots, rep(list(patchwork::plot_spacer()), pad_n))
  }

  page_plot <- patchwork::wrap_plots(pair_plots, ncol = ncol_page, guides = "collect") +
    patchwork::plot_annotation(title = paste("Bubble plot:", focal))

  print(page_plot)
}
dev.off()

positive_network_plot <- plot_network(
  edge_df = positive_edges,
  node_df = node_summary,
  title_text = "Positive cross-celltype MP correlations",
  subtitle_text = paste0(
    "Edges shown when Pearson > 0, Spearman p < 0.05, and -log10(p) >= ",
    positive_sig_cutoff
  ),
  edge_low = "#F9D7CF",
  edge_high = "#8F0000"
)

negative_network_plot <- plot_network(
  edge_df = negative_edges,
  node_df = node_summary,
  title_text = "Negative cross-celltype MP correlations",
  subtitle_text = "Edges shown when Pearson < 0 and Spearman p < 0.05",
  edge_low = "#DCEAF7",
  edge_high = "#08519C"
)

ggsave(
  file.path(out_dir, "Auto_cross_celltype_positive_network.pdf"),
  positive_network_plot,
  width = 15,
  height = 11
)
ggsave(
  file.path(out_dir, "Auto_cross_celltype_negative_network.pdf"),
  negative_network_plot,
  width = 15,
  height = 11
)

####################
# Ligand-receptor support
####################

message("Loading ligand-receptor reference and annotating positive edges")

lr_ref <- load_ramilowski_reference(project_dir)
lr_status <- tibble(
  status = if (is.null(lr_ref)) "missing" else "loaded",
  source = if (is.null(lr_ref)) NA_character_ else lr_ref$source,
  sheet = if (is.null(lr_ref)) NA_character_ else lr_ref$sheet,
  retained_evidence = paste(pair_evidence_allowed, collapse = "; "),
  n_pairs_all = if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_all,
  n_pairs_retained = if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_retained,
  n_pairs_excluded = if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_excluded
)

score_cache <- list()

get_compartment_scores <- function(compartment_name) {
  if (!is.null(score_cache[[compartment_name]])) {
    return(score_cache[[compartment_name]])
  }

  comp <- compartment_data[[compartment_name]]
  score_df <- readRDS(comp$score_path)
  keep_cells <- intersect(rownames(score_df), comp$cells_use)
  score_df <- score_df[keep_cells, comp$mp_names, drop = FALSE]
  score_cache[[compartment_name]] <<- score_df
  score_df
}

get_positive_cells <- function(compartment_name, mp_name, shared_samples) {
  comp <- compartment_data[[compartment_name]]
  score_df <- get_compartment_scores(compartment_name)
  sample_meta <- comp$sample_meta
  candidate_cells <- sample_meta$cell[sample_meta$sample %in% shared_samples]
  candidate_cells <- intersect(candidate_cells, rownames(score_df))
  if (length(candidate_cells) == 0) {
    return(character(0))
  }
  candidate_cells[score_df[candidate_cells, mp_name] > ucell_cutoff]
}

match_driver_target_lr <- function(driver_top_genes, target_mp_genes, edge_row, driver_label, target_label) {
  if (length(driver_top_genes) == 0 || length(target_mp_genes) == 0 || is.null(lr_ref)) {
    return(tibble())
  }

  driver_top_genes <- unique(toupper(driver_top_genes))
  target_mp_genes <- unique(toupper(target_mp_genes))
  lr_pairs <- lr_ref$pairs

  driver_as_ligand <- lr_pairs %>%
    filter(ligand %in% driver_top_genes, receptor %in% target_mp_genes) %>%
    mutate(match_mode = "driver_top_ligand_target_mp_receptor")

  driver_as_receptor <- lr_pairs %>%
    filter(receptor %in% driver_top_genes, ligand %in% target_mp_genes) %>%
    mutate(match_mode = "driver_top_receptor_target_mp_ligand")

  bind_rows(driver_as_ligand, driver_as_receptor) %>%
    mutate(
      edge_id = edge_row$edge_id,
      edge_label = edge_row$edge_label,
      direction = paste(driver_label, "->", target_label),
      driver_node = driver_label,
      target_node = target_label,
      compartment1 = edge_row$compartment1,
      compartment2 = edge_row$compartment2,
      celltype1_display = edge_row$celltype1_display,
      celltype2_display = edge_row$celltype2_display,
      node1_label = edge_row$node1_label,
      node2_label = edge_row$node2_label,
      mp1_name = edge_row$mp1_name,
      mp2_name = edge_row$mp2_name,
      mp1_display = edge_row$mp1_display,
      mp2_display = edge_row$mp2_display,
      celltype_pair = edge_row$celltype_pair,
      pearson_r = edge_row$pearson_r,
      spearman_r = edge_row$spearman_r,
      spearman_significance = edge_row$spearman_significance
    )
}

step3_cache <- load_or_build_cache(
  cache_filename = "Auto_cross_celltype_step3_lr_cache.rds",
  step_name = "ligand-receptor annotation",
  build_fun = function() {
    lr_pairs_df <- tibble()
    lr_edge_summary <- tibble()
    positive_edges_lr <- positive_edges

    if (is.null(lr_ref) || nrow(positive_edges) == 0) {
      return(list(
        lr_pairs_df = lr_pairs_df,
        lr_edge_summary = lr_edge_summary,
        positive_edges_lr = positive_edges_lr
      ))
    }

    message("Extracting top genes from positive cells for LR support")
    merged_expr <- get_assay_data_safely(merged_obj, assay = "RNA")
    lr_rows <- list()

    for (edge_idx in seq_len(nrow(positive_edges))) {
      edge_row <- positive_edges[edge_idx, , drop = FALSE]
      shared_studies <- strsplit(edge_row$eligible_studies, ";", fixed = TRUE)[[1]]
      shared_studies <- shared_studies[nzchar(shared_studies)]

      comp1_scores <- compartment_data[[edge_row$compartment1]]$adjusted_scores
      comp2_scores <- compartment_data[[edge_row$compartment2]]$adjusted_scores
      eligible_samples <- intersect(
        comp1_scores$sample[comp1_scores$study %in% shared_studies],
        comp2_scores$sample[comp2_scores$study %in% shared_studies]
      )
      if (length(eligible_samples) == 0) next

      node1_cells <- get_positive_cells(edge_row$compartment1, edge_row$mp1_name, eligible_samples)
      node2_cells <- get_positive_cells(edge_row$compartment2, edge_row$mp2_name, eligible_samples)
      if (length(node1_cells) == 0 || length(node2_cells) == 0) next

      node1_top <- extract_top_genes(merged_expr, node1_cells, top_n = top_ranked_genes)
      node2_top <- extract_top_genes(merged_expr, node2_cells, top_n = top_ranked_genes)

      node1_mp_genes <- toupper(unique(compartment_data[[edge_row$compartment1]]$gene_sets[[edge_row$mp1_name]]))
      node2_mp_genes <- toupper(unique(compartment_data[[edge_row$compartment2]]$gene_sets[[edge_row$mp2_name]]))

      node1_to_node2 <- match_driver_target_lr(
        driver_top_genes = node1_top,
        target_mp_genes = node2_mp_genes,
        edge_row = edge_row,
        driver_label = edge_row$node1_label,
        target_label = edge_row$node2_label
      ) %>%
        mutate(
          driver_positive_cells = length(node1_cells),
          target_positive_cells = length(node2_cells)
        )

      node2_to_node1 <- match_driver_target_lr(
        driver_top_genes = node2_top,
        target_mp_genes = node1_mp_genes,
        edge_row = edge_row,
        driver_label = edge_row$node2_label,
        target_label = edge_row$node1_label
      ) %>%
        mutate(
          driver_positive_cells = length(node2_cells),
          target_positive_cells = length(node1_cells)
        )

      lr_rows[[length(lr_rows) + 1]] <- bind_rows(node1_to_node2, node2_to_node1)
    }

    lr_pairs_df <- bind_rows(lr_rows)

    if (nrow(lr_pairs_df) > 0) {
      lr_pairs_df <- lr_pairs_df %>%
        distinct(
          edge_id,
          direction,
          ligand,
          receptor,
          match_mode,
          .keep_all = TRUE
        )

      lr_edge_summary <- lr_pairs_df %>%
        group_by(edge_id, edge_label, celltype_pair) %>%
        summarise(
          n_lr_rows = n(),
          n_lr_directions = n_distinct(direction),
          .groups = "drop"
        ) %>%
        left_join(summarise_edge_lr_labels(lr_pairs_df), by = "edge_id") %>%
        arrange(desc(n_literature_supported), desc(n_lr_pairs), desc(max_pair_support_n), desc(n_lr_rows))

      positive_edges_lr <- positive_edges %>%
        left_join(lr_edge_summary, by = "edge_id") %>%
        arrange(desc(spearman_significance), desc(n_lr_pairs), desc(n_literature_supported))
    }

    list(
      lr_pairs_df = lr_pairs_df,
      lr_edge_summary = lr_edge_summary,
      positive_edges_lr = positive_edges_lr
    )
  }
)

lr_pairs_df <- step3_cache$lr_pairs_df
lr_edge_summary <- step3_cache$lr_edge_summary
positive_edges_lr <- step3_cache$positive_edges_lr

write.csv(
  lr_status,
  file.path(out_dir, "Auto_cross_celltype_ligand_receptor_status.csv"),
  row.names = FALSE
)

if (nrow(lr_pairs_df) > 0) {
  write.csv(
    lr_pairs_df,
    file.path(out_dir, "Auto_cross_celltype_ligand_receptor_pairs.csv"),
    row.names = FALSE
  )
  write.csv(
    lr_edge_summary,
    file.path(out_dir, "Auto_cross_celltype_ligand_receptor_edge_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    positive_edges_lr,
    file.path(out_dir, "Auto_cross_celltype_positive_edges_lr_annotated.csv"),
    row.names = FALSE
  )
}

####################
# LR-annotated positive network and Excel exports
####################

excel_output_paths <- character(0)

if (nrow(lr_edge_summary) > 0) {
  edge_annotations <- positive_edges_lr %>%
    filter(!is.na(n_lr_pairs), n_lr_pairs > 0) %>%
    slice_head(n = annotated_positive_edges_n) %>%
    transmute(edge_id, edge_note = top_lr_label)

  positive_lr_network_plot <- plot_network(
    edge_df = positive_edges_lr,
    node_df = node_summary,
    title_text = "Positive cross-celltype MP correlations with LR support",
    subtitle_text = paste0(
      "Only literature-supported and putative LR pairs retained; positivity = UCell > ",
      ucell_cutoff
    ),
    edge_low = "#F9D7CF",
    edge_high = "#8F0000",
    annotation_df = edge_annotations
  )

  lr_edge_barplot <- lr_edge_summary %>%
    slice_head(n = 25) %>%
    mutate(edge_label = factor(edge_label, levels = rev(edge_label))) %>%
    ggplot(aes(x = n_lr_pairs, y = edge_label, fill = n_literature_supported)) +
    geom_col() +
    scale_fill_gradient(low = "#FEE0D2", high = "#A50F15") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16)
    ) +
    labs(
      title = "Top positive edges by LR support",
      subtitle = "Counts shown after removing all excluded LR pairs",
      x = "Number of candidate LR pairs",
      y = NULL,
      fill = "Literature\nsupported"
    )

  lr_pair_barplot <- lr_pairs_df %>%
    mutate(pair_label = ifelse(!is.na(pair_name) & pair_name != "", pair_name, paste(ligand, receptor, sep = " -> "))) %>%
    count(pair_label, pair_evidence, sort = TRUE, name = "n_edges_supported") %>%
    slice_head(n = 25) %>%
    mutate(pair_label = factor(pair_label, levels = rev(unique(pair_label)))) %>%
    ggplot(aes(x = n_edges_supported, y = pair_label, fill = pair_evidence)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("literature supported" = "#66A61E", "putative" = "#E6AB02")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16)
    ) +
    labs(
      title = "Most recurrent retained LR pairs",
      x = "Positive edges supported",
      y = NULL,
      fill = "Pair evidence"
    )

  ggsave(
    file.path(out_dir, "Auto_cross_celltype_positive_network_lr_annotated.pdf"),
    positive_lr_network_plot,
    width = 16,
    height = 12
  )
  ggsave(
    file.path(out_dir, "Auto_cross_celltype_ligand_receptor_summary.pdf"),
    lr_edge_barplot / lr_pair_barplot + plot_layout(heights = c(1, 1)),
    width = 12,
    height = 12
  )

  excel_output_paths <- write_lr_workbooks(
    lr_pairs_df = lr_pairs_df,
    lr_edge_summary = lr_edge_summary,
    out_dir = out_dir
  )
}

####################
# Final summary
####################

####################
# Keep the summary table strictly character-typed so mixed numeric/path rows
# can be combined without vctrs type errors at the end of the run.
####################
summary_row <- function(metric, value) {
  tibble(metric = metric, value = as.character(value))
}

summary_df <- bind_rows(
  summary_row("cache_version", cache_version),
  summary_row("ucell_cutoff", ucell_cutoff),
  summary_row("n_compartments", length(compartment_data)),
  summary_row("n_nodes_total", nrow(node_summary)),
  summary_row("n_nodes_passing_coverage", sum(node_summary$coverage_pass)),
  summary_row("n_pairwise_tests", nrow(result_df)),
  summary_row("n_positive_edges", nrow(positive_edges)),
  summary_row("n_negative_edges", nrow(negative_edges)),
  summary_row("n_lr_rows", nrow(lr_pairs_df)),
  summary_row("n_lr_supported_edges", if (nrow(lr_edge_summary) == 0) 0 else sum(lr_edge_summary$n_lr_pairs > 0)),
  summary_row("n_lr_pairs_all_reference", if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_all),
  summary_row("n_lr_pairs_retained_reference", if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_retained),
  summary_row("n_lr_pairs_excluded_reference", if (is.null(lr_ref)) NA_integer_ else lr_ref$n_pairs_excluded),
  summary_row("excel_outputs", if (length(excel_output_paths) == 0) "" else paste(excel_output_paths, collapse = "; "))
)

write.csv(
  summary_df,
  file.path(summary_dir, "Auto_mp_cross_celltype_correlations_summary.csv"),
  row.names = FALSE
)

message("Saved outputs to: ", file.path(getwd(), out_dir))
}
