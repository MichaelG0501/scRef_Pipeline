####################
# Analysis registry:
#   Status: terminal
#   Script: analysis/cell_states/final_mp_scenic_parse_overlap.R
#   Methodology: analysis/methodology/cell_states/scenic_parse_overlap_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Compare within-dataset, regulon-enrichment signatures of
#     scRef MPs/final states with balanced Parse treatment timepoints.
#   Inputs:
#     ref_outs/final_mp_scenic/Auto_final_mp_scenic_rss.rds
#     ref_outs/final_mp_scenic/Auto_final_mp_scenic_state_rss.rds
#     /rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/
#       parse_outs/cell_states/timepoint_scenic/intermediate/
#       combined_timepoint_balanced2600_rss.rds (preferred) or
#       combined_timepoint_rss.rds (fallback)
#   Outputs:
#     intermediate/: canonical RSS and pairwise similarity matrices
#     tables/: MP/state × timepoint evidence and leading-regulon tables
#     figures/: separate primary-similarity and top-regulon-overlap heatmaps,
#       plus per-regulon evidence-profile PDFs for each MP/state's closest
#       timepoint
#     logs/: run metadata
#     updates/new_updates/summaries/: compact best-match summary
#   Cache/replot: no costly computation; all tables/figures are regenerated
#     directly from the saved RSS matrices.
#   Run: Rscript analysis/cell_states/final_mp_scenic_parse_overlap.R
#   Env: dmtcp
####################

####################
# Cross-dataset regulon-enrichment comparison
#
# RSS values are only comparable *within* the SCENIC run from which they were
# calculated. Consequently, this script does not correlate raw RSS values or
# min-max scale a joint scRef/Parse matrix. For every canonical TF regulon it
# instead z-scores RSS across labels within each dataset, then compares the
# resulting enrichment signatures. A positive cosine therefore means that the
# same TF regulons are selectively enriched in an MP/state and a timepoint
# relative to their own respective comparison labels.
####################

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

####################
# Configuration and utility functions
####################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(x[1])) return(y)
  x[1]
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

canonical_tf <- function(x) {
  x <- trimws(x)
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  x <- gsub(" \\([0-9]+ genes\\)$", "", x)
  gsub("_extended$", "", x)
}

canonicalise_rss <- function(rss_mat) {
  raw_names <- rownames(rss_mat)
  tf_names <- canonical_tf(raw_names)
  is_extended <- grepl("_extended", raw_names, fixed = TRUE)
  selected_index <- unlist(lapply(split(seq_along(raw_names), tf_names), function(indices) {
    indices[order(is_extended[indices], raw_names[indices])][1]
  }), use.names = FALSE)
  selected_index <- selected_index[order(tf_names[selected_index])]
  out <- rss_mat[selected_index, , drop = FALSE]
  rownames(out) <- tf_names[selected_index]
  list(
    rss = out,
    selected_raw_regulon = setNames(raw_names[selected_index], rownames(out)),
    n_collapsed = length(raw_names) - nrow(out)
  )
}

row_zscore <- function(mat) {
  out <- t(apply(mat, 1, function(values) {
    value_sd <- stats::sd(values, na.rm = TRUE)
    if (!is.finite(value_sd) || value_sd == 0) return(rep(0, length(values)))
    as.numeric((values - mean(values, na.rm = TRUE)) / value_sd)
  }))
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

top_enriched_regulons <- function(z_mat, label, n) {
  values <- z_mat[, label]
  values <- values[is.finite(values) & values > 0]
  values <- sort(values, decreasing = TRUE)
  names(values)[seq_len(min(n, length(values)))]
}

positive_cosine <- function(x, y) {
  x <- pmax(x, 0)
  y <- pmax(y, 0)
  denominator <- sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2))
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  sum(x * y) / denominator
}

weighted_jaccard <- function(first, second) {
  first_weights <- setNames(rev(seq_along(first)) / length(first), first)
  second_weights <- setNames(rev(seq_along(second)) / length(second), second)
  features <- union(names(first_weights), names(second_weights))
  if (length(features) == 0) return(NA_real_)
  first_values <- first_weights[features]
  second_values <- second_weights[features]
  first_values[is.na(first_values)] <- 0
  second_values[is.na(second_values)] <- 0
  sum(pmin(first_values, second_values)) / sum(pmax(first_values, second_values))
}

format_leading_regulons <- function(first_z, second_z, max_n = 8) {
  concordant <- intersect(names(first_z)[first_z > 0], names(second_z)[second_z > 0])
  if (length(concordant) == 0) return("")
  contribution <- first_z[concordant] * second_z[concordant]
  paste(names(sort(contribution, decreasing = TRUE))[seq_len(min(max_n, length(contribution)))], collapse = "; ")
}

make_similarity_table <- function(sc_z, parse_z, entities, timepoints, top_n, entity_name) {
  sc_top <- lapply(entities, function(entity) top_enriched_regulons(sc_z, entity, top_n))
  names(sc_top) <- entities
  parse_top <- lapply(timepoints, function(timepoint) top_enriched_regulons(parse_z, timepoint, top_n))
  names(parse_top) <- timepoints

  detail <- bind_rows(lapply(entities, function(entity) {
    bind_rows(lapply(timepoints, function(timepoint) {
      first_z <- sc_z[, entity]
      second_z <- parse_z[, timepoint]
      shared_top <- intersect(sc_top[[entity]], parse_top[[timepoint]])
      data.frame(
        entity = entity,
        timepoint = timepoint,
        n_shared_regulons = length(intersect(names(first_z), names(second_z))),
        enrichment_cosine = positive_cosine(first_z, second_z),
        signed_spearman = suppressWarnings(stats::cor(first_z, second_z, method = "spearman")),
        top_regulon_weighted_jaccard = weighted_jaccard(sc_top[[entity]], parse_top[[timepoint]]),
        n_shared_top_regulons = length(shared_top),
        shared_top_regulons = paste(shared_top, collapse = "; "),
        leading_concordant_regulons = format_leading_regulons(first_z, second_z),
        stringsAsFactors = FALSE
      )
    }))
  }))
  names(detail)[names(detail) == "entity"] <- entity_name
  list(detail = detail, sc_top = sc_top, parse_top = parse_top)
}

table_to_matrix <- function(detail, row_name, value, row_order, column_order) {
  out <- matrix(NA_real_, nrow = length(row_order), ncol = length(column_order),
                dimnames = list(row_order, column_order))
  for (i in seq_len(nrow(detail))) {
    out[detail[[row_name]][i], detail$timepoint[i]] <- detail[[value]][i]
  }
  out
}

plot_similarity_heatmap <- function(detail, row_name, row_order, timepoints, output_file,
                                    title, row_annotation = NULL, row_split = NULL,
                                    row_fontsize = 9) {
  cosine_mat <- table_to_matrix(detail, row_name, "enrichment_cosine", row_order, timepoints)
  best_columns <- apply(cosine_mat, 1, function(values) which.max(values))
  cosine_colours <- colorRamp2(c(0, 0.25, 0.5, 0.75),
                                c("#FFFFFF", "#DCEAF7", "#77ADD7", "#1F5A94"))
  timepoint_annotation <- HeatmapAnnotation(
    Timepoint = factor(timepoints, levels = timepoints),
    col = list(Timepoint = c(
      "T0" = "#0072B2", "T1" = "#E69F00", "T2" = "#009E73",
      "T4" = "#D55E00", "R4" = "#CC79A7", "eR4" = "#56B4E9"
    )[timepoints]),
    show_annotation_name = FALSE
  )
  cell_labels <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cosine_mat[i, j]), x, y,
              gp = gpar(fontsize = if (length(row_order) > 10) 8 else 10,
                        col = ifelse(cosine_mat[i, j] >= 0.52, "white", "black")))
    if (j == best_columns[i]) {
      grid.rect(x, y, width * 0.91, height * 0.84,
                gp = gpar(fill = NA, col = "#111111", lwd = 1.6))
    }
  }
  pdf(output_file, width = 14, height = if (length(row_order) > 10) 14 else 8,
      useDingbats = FALSE)
  draw(
    Heatmap(
      cosine_mat,
      name = "Enrichment\ncosine",
      col = cosine_colours,
      top_annotation = timepoint_annotation,
      left_annotation = row_annotation,
      row_split = row_split,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = row_fontsize),
      column_names_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = cell_labels,
      heatmap_legend_param = list(title = "Regulon\nenrichment\nsimilarity")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(title, x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text("Darker = more similar selectively enriched TF-regulons. Black outline = closest timepoint for that row.",
            x = unit(4, "mm"), y = unit(1, "npc") - unit(10, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 8))
  dev.off()
  invisible(list(cosine = cosine_mat, best_columns = best_columns))
}

plot_top_regulon_overlap_heatmap <- function(detail, row_name, row_order, timepoints, output_file,
                                              title, row_annotation = NULL, row_split = NULL,
                                              row_fontsize = 9) {
  overlap_mat <- table_to_matrix(detail, row_name, "top_regulon_weighted_jaccard", row_order, timepoints)
  best_columns <- apply(overlap_mat, 1, function(values) which.max(values))
  overlap_colours <- colorRamp2(c(0, 0.05, 0.12, 0.22), c("#FFFFFF", "#FEE8C8", "#FDAE6B", "#E6550D"))
  timepoint_annotation <- HeatmapAnnotation(
    Timepoint = factor(timepoints, levels = timepoints),
    col = list(Timepoint = c(
      "T0" = "#0072B2", "T1" = "#E69F00", "T2" = "#009E73",
      "T4" = "#D55E00", "R4" = "#CC79A7", "eR4" = "#56B4E9"
    )[timepoints]),
    show_annotation_name = FALSE
  )
  cell_labels <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", overlap_mat[i, j]), x, y,
              gp = gpar(fontsize = if (length(row_order) > 10) 8 else 10,
                        col = ifelse(overlap_mat[i, j] >= 0.14, "white", "black")))
    if (j == best_columns[i]) {
      grid.rect(x, y, width * 0.91, height * 0.84,
                gp = gpar(fill = NA, col = "#111111", lwd = 1.6))
    }
  }
  pdf(output_file, width = 14, height = if (length(row_order) > 10) 14 else 8,
      useDingbats = FALSE)
  draw(
    Heatmap(
      overlap_mat,
      name = "Weighted\nJaccard",
      col = overlap_colours,
      top_annotation = timepoint_annotation,
      left_annotation = row_annotation,
      row_split = row_split,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = row_fontsize),
      column_names_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_rot = 0,
      cell_fun = cell_labels,
      heatmap_legend_param = list(title = "Top-regulon\nweighted Jaccard")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(title, x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text("Supporting evidence only: overlap of the top enriched TF regulons. It is not used to select the primary closest timepoint. Black outline = greatest overlap in that row.",
            x = unit(4, "mm"), y = unit(1, "npc") - unit(10, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 8))
  dev.off()
  invisible(list(jaccard = overlap_mat, best_columns = best_columns))
}

plot_regulon_evidence_profiles <- function(sc_z, parse_z, best_table, entity_name,
                                           timepoints, output_file, profile_n) {
  z_colours <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#FFFFFF", "#B2182B"))
  pdf(output_file, width = 12, height = 8, onefile = TRUE, useDingbats = FALSE)
  for (entity in best_table[[entity_name]]) {
    pair <- best_table[best_table[[entity_name]] == entity, , drop = FALSE]
    best_timepoint <- pair$best_timepoint[1]
    first_z <- sc_z[, entity]
    second_z <- parse_z[, best_timepoint]
    concordant <- intersect(names(first_z)[first_z > 0], names(second_z)[second_z > 0])
    contribution <- first_z[concordant] * second_z[concordant]
    selected_tfs <- names(sort(contribution, decreasing = TRUE))[seq_len(min(profile_n, length(contribution)))]
    if (length(selected_tfs) == 0) {
      plot.new()
      title(main = paste(entity, "has no positively concordant canonical TF regulons."))
      next
    }
    profile_mat <- cbind(`scRef entity` = first_z[selected_tfs], parse_z[selected_tfs, timepoints, drop = FALSE])
    best_column <- match(best_timepoint, colnames(profile_mat))
    column_annotation <- HeatmapAnnotation(
      Dataset = factor(c("scRef", rep("Parse", length(timepoints))), levels = c("scRef", "Parse")),
      col = list(Dataset = c("scRef" = "#5B5B5B", "Parse" = "#4C78A8")),
      show_annotation_name = FALSE
    )
    evidence_cell_fun <- function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", profile_mat[i, j]), x, y,
                gp = gpar(fontsize = 9, col = ifelse(abs(profile_mat[i, j]) > 1.35, "white", "black")))
      if (j == best_column) {
        grid.rect(x, y, width * 0.92, height * 0.88,
                  gp = gpar(fill = NA, col = "#111111", lwd = 1.7))
      }
    }
    draw(
      Heatmap(
        profile_mat,
        name = "Within-run\nRSS z-score",
        col = z_colours,
        top_annotation = column_annotation,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 11, fontface = "bold"),
        column_names_gp = gpar(fontsize = 11, fontface = "bold"),
        cell_fun = evidence_cell_fun,
        heatmap_legend_param = list(title = "Within-run RSS z-score")
      ),
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
    grid.text(paste0(entity_name, ": ", entity, " | closest Parse timepoint: ", best_timepoint,
                     " (enrichment cosine = ", sprintf("%.2f", pair$enrichment_cosine[1]), ")"),
              x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
              just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.text("Rows are the leading TF regulons enriched in both the scRef entity and its closest timepoint. Values are RSS z-scores calculated separately within each SCENIC run; red = selectively enriched. The black outline marks the matched timepoint.",
              x = unit(4, "mm"), y = unit(1, "npc") - unit(10, "mm"),
              just = c("left", "top"), gp = gpar(fontsize = 8))
  }
  dev.off()
}

####################
# Load and canonicalise within-run RSS matrices
####################
arg_list <- parse_args(commandArgs(trailingOnly = TRUE))
top_n <- as.integer(arg_list[["top_n"]] %||% "20")
profile_n <- as.integer(arg_list[["profile_n"]] %||% "10")
if (!is.finite(top_n) || top_n < 2) stop("top_n must be an integer of at least 2.")
if (!is.finite(profile_n) || profile_n < 3) stop("profile_n must be an integer of at least 3.")

out_base <- "final_mp_scenic/parse_overlap"
tables_dir <- file.path(out_base, "tables")
figures_dir <- file.path(out_base, "figures")
intermediate_dir <- file.path(out_base, "intermediate")
logs_dir <- file.path(out_base, "logs")
summary_dir <- "../updates/new_updates/summaries"
for (directory in c(tables_dir, figures_dir, intermediate_dir, logs_dir, summary_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

sc_mp_rss_path <- "final_mp_scenic/Auto_final_mp_scenic_rss.rds"
sc_state_rss_path <- "final_mp_scenic/Auto_final_mp_scenic_state_rss.rds"
parse_base <- "/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/parse_outs/cell_states/timepoint_scenic/intermediate"
parse_balanced_rss_path <- file.path(parse_base, "combined_timepoint_balanced2600_rss.rds")
parse_full_rss_path <- file.path(parse_base, "combined_timepoint_rss.rds")
parse_rss_path <- if (file.exists(parse_balanced_rss_path)) parse_balanced_rss_path else parse_full_rss_path

if (!all(file.exists(c(sc_mp_rss_path, sc_state_rss_path, parse_rss_path)))) {
  stop("Required RSS input is missing. Expected: ", paste(c(sc_mp_rss_path, sc_state_rss_path, parse_rss_path), collapse = "; "))
}

message("Loading scRef MP, state, and Parse RSS matrices...")
sc_mp_raw <- as.matrix(readRDS(sc_mp_rss_path))
sc_state_raw <- as.matrix(readRDS(sc_state_rss_path))
parse_raw <- as.matrix(readRDS(parse_rss_path))

sc_mp_canonical <- canonicalise_rss(sc_mp_raw)
sc_state_canonical <- canonicalise_rss(sc_state_raw)
parse_canonical <- canonicalise_rss(parse_raw)

parse_timepoints <- c("T0", "T1", "T2", "T4", "R4", "eR4")
parse_timepoints <- intersect(parse_timepoints, colnames(parse_canonical$rss))
if (length(parse_timepoints) < 2) stop("Fewer than two expected Parse timepoints were found in the RSS matrix.")

sc_mps <- c("MP1", "MP5", "MP13+", "MP2+", "MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+",
            "MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP12", "MP15")
sc_mp_descriptions <- c(
  "MP1" = "G2/M cell cycle", "MP5" = "G1/S cell cycle", "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis", "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium", "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response", "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia", "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia", "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation", "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)
mp_group_map <- c(
  "MP1" = "Cell cycle", "MP5" = "Cell cycle", "MP13+" = "Cell cycle", "MP2+" = "Classic proliferation",
  "MP14" = "Basal to intestinal metaplasia", "MP3+" = "Basal to intestinal metaplasia",
  "MP6+" = "Basal to intestinal metaplasia", "MP11+" = "Basal to intestinal metaplasia",
  "MP9+" = "Basal to intestinal metaplasia", "MP10+" = "Basal to intestinal metaplasia",
  "MP8+" = "SMG to intestinal metaplasia", "MP8b" = "SMG to intestinal metaplasia",
  "MP16" = "SMG to intestinal metaplasia", "MP18b" = "SMG to intestinal metaplasia",
  "MP17" = "SMG to intestinal metaplasia",
  "MP12" = "Stress adaptive", "MP15" = "Cancer-cell immune mimicry"
)
mp_group_cols <- c(
  "Cell cycle" = "#6B7280", "Classic proliferation" = "#E41A1C",
  "Basal to intestinal metaplasia" = "#4DAF4A", "SMG to intestinal metaplasia" = "#FF7F00",
  "Stress adaptive" = "#984EA3", "Cancer-cell immune mimicry" = "#377EB8"
)
state_order <- c("Classic proliferation", "Basal to intestinal metaplasia", "SMG to intestinal metaplasia",
                 "Stress adaptive", "Cancer-cell immune mimicry")

mp_labels <- setNames(paste(sc_mps, sc_mp_descriptions[sc_mps]), sc_mps)
mp_columns <- unname(mp_labels[mp_labels %in% colnames(sc_mp_canonical$rss)])
state_columns <- intersect(state_order, colnames(sc_state_canonical$rss))
if (length(mp_columns) < 2 || length(state_columns) < 2) {
  stop("Could not recover the expected final-MP/state RSS columns.")
}

shared_mp_tfs <- intersect(rownames(sc_mp_canonical$rss), rownames(parse_canonical$rss))
shared_state_tfs <- intersect(rownames(sc_state_canonical$rss), rownames(parse_canonical$rss))
if (length(shared_mp_tfs) < 10 || length(shared_state_tfs) < 10) {
  stop("Too few canonical TF regulons are shared between scRef and Parse RSS matrices.")
}

message("Comparing ", length(shared_mp_tfs), " shared canonical TF regulons for MPs and ",
        length(shared_state_tfs), " for states using ", basename(parse_rss_path), ".")
mp_sc_z <- row_zscore(sc_mp_canonical$rss[shared_mp_tfs, mp_columns, drop = FALSE])
state_sc_z <- row_zscore(sc_state_canonical$rss[shared_state_tfs, state_columns, drop = FALSE])
mp_parse_z <- row_zscore(parse_canonical$rss[shared_mp_tfs, parse_timepoints, drop = FALSE])
state_parse_z <- row_zscore(parse_canonical$rss[shared_state_tfs, parse_timepoints, drop = FALSE])

####################
# Pairwise evidence: enrichment concordance and leading TF regulons
####################
mp_results <- make_similarity_table(mp_sc_z, mp_parse_z, mp_columns, parse_timepoints, top_n, "mp")
state_results <- make_similarity_table(state_sc_z, state_parse_z, state_columns, parse_timepoints, top_n, "state")
mp_detail <- mp_results$detail
state_detail <- state_results$detail

write.csv(mp_detail, file.path(tables_dir, "mp_timepoint_regulon_enrichment_similarity.csv"), row.names = FALSE)
write.csv(state_detail, file.path(tables_dir, "state_timepoint_regulon_enrichment_similarity.csv"), row.names = FALSE)

mp_best <- mp_detail %>%
  group_by(mp) %>%
  slice_max(enrichment_cosine, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(mp, best_timepoint = timepoint, enrichment_cosine, signed_spearman,
            top_regulon_weighted_jaccard, n_shared_top_regulons, leading_concordant_regulons)
state_best <- state_detail %>%
  group_by(state) %>%
  slice_max(enrichment_cosine, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(state, best_timepoint = timepoint, enrichment_cosine, signed_spearman,
            top_regulon_weighted_jaccard, n_shared_top_regulons, leading_concordant_regulons)
write.csv(mp_best, file.path(tables_dir, "mp_best_matching_timepoint_by_regulon_enrichment.csv"), row.names = FALSE)
write.csv(state_best, file.path(tables_dir, "state_best_matching_timepoint_by_regulon_enrichment.csv"), row.names = FALSE)

mp_leading <- mp_detail %>%
  select(mp, timepoint, enrichment_cosine, signed_spearman, top_regulon_weighted_jaccard,
         n_shared_top_regulons, shared_top_regulons, leading_concordant_regulons) %>%
  arrange(mp, desc(enrichment_cosine))
state_leading <- state_detail %>%
  select(state, timepoint, enrichment_cosine, signed_spearman, top_regulon_weighted_jaccard,
         n_shared_top_regulons, shared_top_regulons, leading_concordant_regulons) %>%
  arrange(state, desc(enrichment_cosine))
write.csv(mp_leading, file.path(tables_dir, "mp_timepoint_leading_concordant_regulons.csv"), row.names = FALSE)
write.csv(state_leading, file.path(tables_dir, "state_timepoint_leading_concordant_regulons.csv"), row.names = FALSE)

####################
# Presentation heatmaps with timepoint order fixed to the treatment course
####################
mp_ids <- names(mp_labels)[match(mp_columns, mp_labels)]
mp_groups <- unname(mp_group_map[mp_ids])
names(mp_groups) <- mp_columns
mp_row_annotation <- rowAnnotation(
  `MP group` = factor(mp_groups, levels = names(mp_group_cols)),
  col = list(`MP group` = mp_group_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 9)
)

mp_plot_matrices <- plot_similarity_heatmap(
  detail = mp_detail,
  row_name = "mp",
  row_order = mp_columns,
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "mp_timepoint_regulon_enrichment_similarity_heatmap.pdf"),
  title = "Which Parse treatment timepoints share each scRef MP's selectively enriched regulons?",
  row_annotation = mp_row_annotation,
  row_split = factor(mp_groups, levels = names(mp_group_cols)),
  row_fontsize = 8
)
state_plot_matrices <- plot_similarity_heatmap(
  detail = state_detail,
  row_name = "state",
  row_order = state_columns,
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "state_timepoint_regulon_enrichment_similarity_heatmap.pdf"),
  title = "Which Parse treatment timepoints share each scRef state's selectively enriched regulons?",
  row_fontsize = 10
)

mp_jaccard_matrices <- plot_top_regulon_overlap_heatmap(
  detail = mp_detail,
  row_name = "mp",
  row_order = mp_columns,
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "mp_timepoint_top_regulon_jaccard_heatmap.pdf"),
  title = "Supporting evidence: top-TF-regulon overlap between scRef MPs and Parse timepoints",
  row_annotation = mp_row_annotation,
  row_split = factor(mp_groups, levels = names(mp_group_cols)),
  row_fontsize = 8
)
state_jaccard_matrices <- plot_top_regulon_overlap_heatmap(
  detail = state_detail,
  row_name = "state",
  row_order = state_columns,
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "state_timepoint_top_regulon_jaccard_heatmap.pdf"),
  title = "Supporting evidence: top-TF-regulon overlap between scRef states and Parse timepoints",
  row_fontsize = 10
)

plot_regulon_evidence_profiles(
  sc_z = mp_sc_z,
  parse_z = mp_parse_z,
  best_table = mp_best,
  entity_name = "mp",
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "mp_timepoint_regulon_match_evidence_profiles.pdf"),
  profile_n = profile_n
)
plot_regulon_evidence_profiles(
  sc_z = state_sc_z,
  parse_z = state_parse_z,
  best_table = state_best,
  entity_name = "state",
  timepoints = parse_timepoints,
  output_file = file.path(figures_dir, "state_timepoint_regulon_match_evidence_profiles.pdf"),
  profile_n = profile_n
)

####################
# Replot object, concise update summary, and provenance log
####################
saveRDS(
  list(
    method = "within-run RSS row z-score; positive cosine and signed Spearman across shared canonical direct TF regulons",
    parse_rss_input = parse_rss_path,
    top_n = top_n,
    profile_n = profile_n,
    mp_shared_canonical_tfs = shared_mp_tfs,
    state_shared_canonical_tfs = shared_state_tfs,
    mp_sc_z = mp_sc_z,
    state_sc_z = state_sc_z,
    mp_parse_z = mp_parse_z,
    state_parse_z = state_parse_z,
    mp_detail = mp_detail,
    state_detail = state_detail,
    mp_heatmap_matrices = mp_plot_matrices,
    state_heatmap_matrices = state_plot_matrices,
    mp_jaccard_matrices = mp_jaccard_matrices,
    state_jaccard_matrices = state_jaccard_matrices
  ),
  file.path(intermediate_dir, "regulon_enrichment_similarity_intermediate.rds")
)

summary_table <- bind_rows(
  mp_best %>% transmute(level = "MP", entity = mp, best_timepoint, enrichment_cosine,
                         signed_spearman, top_regulon_weighted_jaccard, n_shared_top_regulons,
                         leading_concordant_regulons),
  state_best %>% transmute(level = "State", entity = state, best_timepoint, enrichment_cosine,
                            signed_spearman, top_regulon_weighted_jaccard, n_shared_top_regulons,
                            leading_concordant_regulons)
)
write.csv(summary_table, file.path(summary_dir, "final_mp_scenic_parse_overlap_summary.csv"), row.names = FALSE)

run_log <- c(
  paste("timestamp", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), sep = ": "),
  paste("scRef MP RSS", normalizePath(sc_mp_rss_path), sep = ": "),
  paste("scRef state RSS", normalizePath(sc_state_rss_path), sep = ": "),
  paste("Parse RSS", parse_rss_path, sep = ": "),
  paste("Parse balance", if (identical(parse_rss_path, parse_balanced_rss_path)) "balanced 2600 cells/timepoint" else "unbalanced fallback", sep = ": "),
  paste("canonical TFs shared for MPs", length(shared_mp_tfs), sep = ": "),
  paste("canonical TFs shared for states", length(shared_state_tfs), sep = ": "),
  paste("top_n", top_n, sep = ": "),
  paste("profile_n", profile_n, sep = ": "),
  "primary interpretation: high positive cosine means matching within-dataset TF-regulon enrichment, not equality of raw RSS or target-gene networks."
)
writeLines(run_log, file.path(logs_dir, "final_mp_scenic_parse_overlap_run.txt"))

message("Done. Primary figures: ", normalizePath(figures_dir))
message("The old raw-RSS/min-max combined score is intentionally not regenerated.")
