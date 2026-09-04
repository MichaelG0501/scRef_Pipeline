####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_refined_mp_split_correlation_ordered_heatmap.R
#   Methodology: analysis/methodology/metaprograms/refined_mp_split_correlation_ordered_heatmap_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Program-resolution refined MP correlation heatmap aligned to
#   refined_mp_nmf_ordered_heatmap.R. Finalized merged MP blocks keep the same
#   heights/widths as in the ordered NMF heatmap, while each block is internally
#   split by its pre-merge refined sub-MPs.
#
# Inputs:
#   - ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_ucell_scores.rds
#   - ref_outs/EAC_Ref_epi.rds
#
# Outputs:
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.pdf
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.png
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_split_correlation_ordered_final_blocks.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_split_correlation_ordered_sub_blocks.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_split_correlation_ordered_matrix.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_split_display_correlation_matrices.rds
#   - updates/new_updates/summaries/refined_mp_split_correlation_ordered_heatmap_summary.csv
#
# Cache/replot behavior:
#   SCREF_FORCE_REBUILD=TRUE recomputes split-feature Fisher-Z correlations
#   from refined UCell scores. Otherwise valid cached correlations are reused.
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R
#
# Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
ref_dir <- file.path(project_dir, "ref_outs")
setwd(ref_dir)

outdir <- file.path("Metaprogrammes_Results", "mp_refinement")
for (subdir in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
}
summary_dir <- file.path(project_dir, "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

force_rebuild <- Sys.getenv("SCREF_FORCE_REBUILD", "FALSE") == "TRUE"

nmf_file <- file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")
assignment_file <- file.path(outdir, "intermediate", "merged_refined_mp_assignments.rds")
gene_file <- file.path(outdir, "intermediate", "merged_refined_mp_genes.rds")
refined_ucell_file <- file.path(outdir, "intermediate", "refined_ucell_scores.rds")
epi_file <- "EAC_Ref_epi.rds"

required_inputs <- c(nmf_file, assignment_file, gene_file, refined_ucell_file, epi_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required input(s). Run mp_refinement_submp.R and ",
    "mp_refinement_merge_correlated_submps.R first: ",
    paste(missing_inputs, collapse = ", ")
  )
}

cat("Loading finalized refined MP assignments and NMF layout...\n")
geneNMF.metaprograms <- readRDS(nmf_file)
merged_assignments <- readRDS(assignment_file)
merged_mp_genes <- readRDS(gene_file)

required_assignment_cols <- c("program", "original_mp", "refined_submp", "merged_refined_mp")
missing_cols <- setdiff(required_assignment_cols, colnames(merged_assignments))
if (length(missing_cols) > 0) {
  stop("Assignment table is missing column(s): ", paste(missing_cols, collapse = ", "))
}

program_similarity <- geneNMF.metaprograms$programs.similarity
original_tree_order <- rownames(program_similarity)[geneNMF.metaprograms$programs.tree$order]

merged_assignments <- merged_assignments |>
  mutate(
    program = as.character(program),
    original_mp = as.character(original_mp),
    refined_submp = as.character(refined_submp),
    merged_refined_mp = as.character(merged_refined_mp)
  ) |>
  filter(program %in% rownames(program_similarity))

if (anyDuplicated(merged_assignments$program) > 0) {
  duplicated_programs <- unique(merged_assignments$program[duplicated(merged_assignments$program)])
  stop("Programs have duplicate finalized refined MP assignments: ",
       paste(head(duplicated_programs, 20), collapse = ", "))
}

####################
# Strict finalized MP order shared with refined_mp_nmf_ordered_heatmap.R.
####################
strict_refined_mp_order <- c(
  "MP7j", "MP9", "MP1", "MP2+", "MP17", "MP8+", "MP10+", "MP14", "MP5+",
  "MP7r", "MP7v", "MP10e", "MP16+", "MP18",
  "MP8c", "MP15c", "MP12c", "MP2v", "MP8e", "MP12a", "MP13",
  "MP7+", "MP7h", "MP8b", "MP12b", "MP15a", "MP15b"
)
####################

finalized_mps <- unique(merged_assignments$merged_refined_mp)
missing_from_strict_order <- setdiff(finalized_mps, strict_refined_mp_order)
if (length(missing_from_strict_order) > 0) {
  stop(
    "Finalized refined MP(s) are absent from strict_refined_mp_order: ",
    paste(missing_from_strict_order, collapse = ", "),
    ". Update the order explicitly before plotting."
  )
}

ordered_features <- strict_refined_mp_order[strict_refined_mp_order %in% finalized_mps]
program_rank <- setNames(seq_along(original_tree_order), original_tree_order)

ordered_programs <- unlist(lapply(ordered_features, function(feature) {
  feature_rows <- merged_assignments[merged_assignments$merged_refined_mp == feature, , drop = FALSE]
  split_order <- unique(feature_rows$refined_submp[order(program_rank[feature_rows$program])])
  unlist(lapply(split_order, function(split_feature) {
    progs <- feature_rows$program[feature_rows$refined_submp == split_feature]
    progs[order(program_rank[progs])]
  }), use.names = FALSE)
}), use.names = FALSE)

missing_programs <- setdiff(merged_assignments$program, ordered_programs)
if (length(missing_programs) > 0) {
  stop("Assigned program(s) were not placed in the ordered heatmap: ",
       paste(head(missing_programs, 20), collapse = ", "))
}

ordered_assignments <- merged_assignments[match(ordered_programs, merged_assignments$program), , drop = FALSE]
final_feature_by_program <- setNames(ordered_assignments$merged_refined_mp, ordered_programs)
split_feature_by_program <- setNames(ordered_assignments$refined_submp, ordered_programs)

final_runs <- rle(unname(final_feature_by_program))
final_end <- cumsum(final_runs$lengths)
final_start <- final_end - final_runs$lengths + 1L

split_runs <- rle(paste(unname(final_feature_by_program), unname(split_feature_by_program), sep = "||"))
split_end <- cumsum(split_runs$lengths)
split_start <- split_end - split_runs$lengths + 1L
split_parts <- strsplit(split_runs$values, "\\|\\|")

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

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

display_label <- function(x) {
  vapply(x, function(xi) {
    if (xi %in% names(submp_desc_map)) {
      return(paste0(xi, ": ", submp_desc_map[[xi]]))
    }
    parent <- parent_id(xi)
    if (parent %in% names(mp_desc_map)) {
      return(paste0(xi, ": ", mp_desc_map[[parent]]))
    }
    xi
  }, character(1), USE.NAMES = FALSE)
}

final_block_table <- data.frame(
  refined_mp = final_runs$values,
  parent_mp = parent_id(final_runs$values),
  strict_order_rank = match(final_runs$values, strict_refined_mp_order),
  start = final_start,
  end = final_end,
  n_programs = final_runs$lengths,
  n_split_blocks = vapply(final_runs$values, function(x) {
    length(unique(ordered_assignments$refined_submp[ordered_assignments$merged_refined_mp == x]))
  }, integer(1)),
  n_genes = vapply(final_runs$values, function(x) length(merged_mp_genes[[x]]), integer(1)),
  label = display_label(final_runs$values),
  stringsAsFactors = FALSE
)

split_block_table <- data.frame(
  finalized_mp = vapply(split_parts, `[`, character(1), 1),
  split_mp = vapply(split_parts, `[`, character(1), 2),
  parent_mp = parent_id(vapply(split_parts, `[`, character(1), 2)),
  start = split_start,
  end = split_end,
  n_programs = split_runs$lengths,
  stringsAsFactors = FALSE
)

write.csv(final_block_table,
          file.path(outdir, "tables", "refined_mp_split_correlation_ordered_final_blocks.csv"),
          row.names = FALSE)
write.csv(split_block_table,
          file.path(outdir, "tables", "refined_mp_split_correlation_ordered_sub_blocks.csv"),
          row.names = FALSE)

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

display_split_features <- unique(unname(split_feature_by_program))
cached_cor <- file.path(outdir, "intermediate", "refined_mp_split_display_correlation_matrices.rds")
split_cor <- NULL
cor_cache_valid <- FALSE
if (file.exists(cached_cor) && !force_rebuild) {
  split_cor_cached <- readRDS(cached_cor)
  cor_cache_valid <- all(display_split_features %in% rownames(split_cor_cached$mean_rho)) &&
    all(display_split_features %in% colnames(split_cor_cached$mean_rho))
  if (cor_cache_valid) {
    split_cor <- split_cor_cached
  }
  rm(split_cor_cached)
}

if (is.null(split_cor)) {
  cat("Computing split/full feature correlations from refined UCell scores...\n")
  refined_ucell <- readRDS(refined_ucell_file)
  missing_features <- setdiff(display_split_features, colnames(refined_ucell))
  if (length(missing_features) > 0) {
    stop(
      "refined_ucell_scores.rds is missing split/full feature(s) needed for plotting: ",
      paste(missing_features, collapse = ", ")
    )
  }
  tmdata_all <- readRDS(epi_file)
  sample_vec <- tmdata_all@meta.data$orig.ident[match(rownames(refined_ucell), rownames(tmdata_all@meta.data))]
  sample_vec[is.na(sample_vec)] <- rownames(refined_ucell)[is.na(sample_vec)]
  split_cor <- compute_fisher_cor_p(refined_ucell[, display_split_features, drop = FALSE],
                                    sample_vec, min_cells = 10)
  split_cor$display_split_features <- display_split_features
  saveRDS(split_cor, cached_cor)
  cat("Saved:", cached_cor, "\n")
} else {
  cat("Using cached split/full feature correlations:", cached_cor, "\n")
}

mean_rho_features <- split_cor$mean_rho[display_split_features, display_split_features, drop = FALSE]
split_idx <- match(unname(split_feature_by_program), display_split_features)
rho_ordered <- mean_rho_features[split_idx, split_idx, drop = FALSE]
rownames(rho_ordered) <- ordered_programs
colnames(rho_ordered) <- ordered_programs

saveRDS(
  list(
    rho_program_resolution = rho_ordered,
    ordered_programs = ordered_programs,
    final_feature_by_program = final_feature_by_program,
    split_feature_by_program = split_feature_by_program,
    final_blocks = final_block_table,
    split_blocks = split_block_table,
    feature_correlation = split_cor
  ),
  file.path(outdir, "intermediate", "refined_mp_split_correlation_ordered_matrix.rds")
)

cat("Generating ordered split-MP correlation heatmap...\n")
n_prog <- nrow(rho_ordered)
col_fun <- colorRamp2(c(-0.4, 0, 0.4), c("#2166AC", "#FFFFFF", "#B2182B"))

ht <- Heatmap(
  rho_ordered,
  name = "Mean rho",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = paste0("Split-MP correlation aligned to finalized refined MP blocks (n = ", n_prog, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 4,
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    direction = "horizontal",
    legend_width = unit(8, "cm")
  ),
  width = unit(205, "mm"),
  height = unit(205, "mm")
)

####################
# Match the ordered-NMF label separation and final-block callout style.
####################
label_break_before <- c("MP17", "MP7r", "MP8c", "MP8b", "MP12b", "MP15a", "MP15b")

adjust_label_positions <- function(raw_y, labels, base_gap = 0.020,
                                   break_gap = 0.046, low = 0.018, high = 0.982) {
  n <- length(raw_y)
  label_y <- pmax(low, pmin(high, raw_y))
  gap_before <- rep(base_gap, n)
  gap_before[labels %in% label_break_before] <- break_gap
  gap_before[1] <- 0

  for (i in seq(2L, n)) {
    label_y[i] <- min(label_y[i], label_y[i - 1L] - gap_before[i])
  }
  if (label_y[n] < low) {
    label_y <- label_y + (low - label_y[n])
  }
  for (i in seq(n - 1L, 1L)) {
    label_y[i] <- max(label_y[i], label_y[i + 1L] + gap_before[i + 1L])
  }
  if (label_y[1] > high) {
    label_y <- label_y - (label_y[1] - high)
  }
  if (label_y[n] < low) {
    scaled_gaps <- gap_before[-1L]
    scaled_gaps <- scaled_gaps * ((high - low) / sum(scaled_gaps))
    label_y <- c(high, high - cumsum(scaled_gaps))
  }

  pmax(low, pmin(high, label_y))
}
####################

draw_heatmap <- function() {
  draw(ht, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 165), "mm"))
  decorate_heatmap_body("Mean rho", {
    for (rr in seq_len(nrow(split_block_table))) {
      s <- split_block_table$start[rr]
      e <- split_block_table$end[rr]
      centre <- (s + e - 1) / (2 * n_prog)
      extent <- (e - s + 1) / n_prog
      grid.rect(
        x = unit(centre, "npc"),
        y = unit(1 - centre, "npc"),
        width = unit(extent, "npc"),
        height = unit(extent, "npc"),
        gp = gpar(fill = NA, col = "#444444", lwd = 1.0, lty = "solid")
      )
    }

    raw_y <- 1 - ((final_block_table$start + final_block_table$end - 1) / (2 * n_prog))
    label_y <- adjust_label_positions(raw_y, final_block_table$refined_mp)

    for (rr in seq_len(nrow(final_block_table))) {
      s <- final_block_table$start[rr]
      e <- final_block_table$end[rr]
      centre <- (s + e - 1) / (2 * n_prog)
      extent <- (e - s + 1) / n_prog
      box_y <- 1 - centre
      ly <- label_y[rr]

      grid.rect(
        x = unit(centre, "npc"),
        y = unit(box_y, "npc"),
        width = unit(extent, "npc"),
        height = unit(extent, "npc"),
        gp = gpar(fill = NA, col = "#111111", lwd = 2.2, lty = "dotted")
      )
      grid.lines(
        x = unit(c(centre + extent / 2, 1.008), "npc"),
        y = unit(c(box_y, ly), "npc"),
        gp = gpar(col = "#111111", lwd = 0.55, lty = "dotted")
      )
      grid.text(
        final_block_table$label[rr],
        x = unit(1.012, "npc"),
        y = unit(ly, "npc"),
        just = "left",
        gp = gpar(fontsize = 8.8, col = "#111111", lineheight = 0.86)
      )
    }
  })
}

pdf_path <- file.path(outdir, "figures", "refined_mp_split_correlation_ordered_heatmap.pdf")
png_path <- file.path(outdir, "figures", "refined_mp_split_correlation_ordered_heatmap.png")

cairo_pdf(pdf_path, width = 15.0, height = 9.8)
draw_heatmap()
dev.off()

png(png_path, width = 15.0, height = 9.8, units = "in", res = 450, bg = "white", type = "cairo")
draw_heatmap()
dev.off()

summary_table <- data.frame(
  metric = c(
    "n_ordered_programs",
    "n_finalized_refined_mps",
    "n_split_display_blocks",
    "n_split_display_features",
    "min_programs_per_finalized_mp",
    "max_programs_per_finalized_mp",
    "pdf_output",
    "png_output"
  ),
  value = c(
    n_prog,
    nrow(final_block_table),
    nrow(split_block_table),
    length(display_split_features),
    min(final_block_table$n_programs),
    max(final_block_table$n_programs),
    pdf_path,
    png_path
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table,
          file.path(summary_dir, "refined_mp_split_correlation_ordered_heatmap_summary.csv"),
          row.names = FALSE)

cat("Saved:", pdf_path, "\n")
cat("Saved:", png_path, "\n")
cat("Saved: tables/refined_mp_split_correlation_ordered_final_blocks.csv\n")
cat("Saved: tables/refined_mp_split_correlation_ordered_sub_blocks.csv\n")
cat("Saved: intermediate/refined_mp_split_correlation_ordered_matrix.rds\n")
