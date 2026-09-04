####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_refined_mp_nmf_ordered_heatmap.R
#   Methodology: analysis/methodology/metaprograms/refined_mp_nmf_ordered_heatmap_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Replots the GeneNMF program similarity heatmap after MP refinement, ordering
#   NMF programs by the finalized merged refined MP order from
#   mp_refinement_merge_correlated_submps.R. Dotted borders mark the finalized
#   refined MP blocks, not the original nMP19 clusters.
#
# Inputs:
#   - ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_genes.rds
#
# Outputs:
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.pdf
#   - ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.png
#   - ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_nmf_ordered_blocks.csv
#   - ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_nmf_ordered_similarity.rds
#   - updates/new_updates/summaries/refined_mp_nmf_ordered_heatmap_summary.csv
#
# Cache/replot behavior:
#   This is a plot-only script. It always rebuilds the ordered similarity matrix
#   from the finalized merged refinement assignment cache.
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R
#
# Conda env: dmtcp
####################

suppressPackageStartupMessages({
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

nmf_file <- file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")
assignment_file <- file.path(outdir, "intermediate", "merged_refined_mp_assignments.rds")
gene_file <- file.path(outdir, "intermediate", "merged_refined_mp_genes.rds")

required_inputs <- c(nmf_file, assignment_file, gene_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required input(s). Run mp_refinement_submp.R and ",
    "mp_refinement_merge_correlated_submps.R first: ",
    paste(missing_inputs, collapse = ", ")
  )
}

cat("Loading finalized refined MP inputs...\n")
geneNMF.metaprograms <- readRDS(nmf_file)
merged_assignments <- readRDS(assignment_file)
merged_mp_genes <- readRDS(gene_file)

required_assignment_cols <- c("program", "original_mp", "merged_refined_mp")
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
    merged_refined_mp = as.character(merged_refined_mp)
  ) |>
  filter(program %in% rownames(program_similarity))

if (anyDuplicated(merged_assignments$program) > 0) {
  duplicated_programs <- unique(merged_assignments$program[duplicated(merged_assignments$program)])
  stop("Programs have duplicate finalized refined MP assignments: ",
       paste(head(duplicated_programs, 20), collapse = ", "))
}

####################
# Strict finalized MP order from mp_refinement_merge_correlated_submps.R.
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
if (!identical(sort(ordered_features), sort(finalized_mps))) {
  stop("The strict order does not resolve to exactly the finalized refined MP set.")
}

program_rank <- setNames(seq_along(original_tree_order), original_tree_order)
ordered_programs <- unlist(lapply(ordered_features, function(feature) {
  progs <- merged_assignments$program[merged_assignments$merged_refined_mp == feature]
  progs[order(program_rank[progs])]
}), use.names = FALSE)

missing_programs <- setdiff(merged_assignments$program, ordered_programs)
if (length(missing_programs) > 0) {
  stop("Assigned program(s) were not placed in the ordered heatmap: ",
       paste(head(missing_programs, 20), collapse = ", "))
}

sim_ordered <- program_similarity[ordered_programs, ordered_programs, drop = FALSE]
feature_by_program <- setNames(
  merged_assignments$merged_refined_mp[match(ordered_programs, merged_assignments$program)],
  ordered_programs
)

runs <- rle(unname(feature_by_program))
run_end <- cumsum(runs$lengths)
run_start <- run_end - runs$lengths + 1L

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

block_table <- data.frame(
  refined_mp = runs$values,
  parent_mp = parent_id(runs$values),
  strict_order_rank = match(runs$values, strict_refined_mp_order),
  start = run_start,
  end = run_end,
  n_programs = runs$lengths,
  n_genes = vapply(runs$values, function(x) length(merged_mp_genes[[x]]), integer(1)),
  label = display_label(runs$values),
  stringsAsFactors = FALSE
)

write.csv(block_table,
          file.path(outdir, "tables", "refined_mp_nmf_ordered_blocks.csv"),
          row.names = FALSE)
saveRDS(
  list(
    similarity = sim_ordered,
    ordered_programs = ordered_programs,
    feature_by_program = feature_by_program,
    blocks = block_table,
    strict_refined_mp_order = strict_refined_mp_order
  ),
  file.path(outdir, "intermediate", "refined_mp_nmf_ordered_similarity.rds")
)

cat("Generating ordered refined MP NMF heatmap...\n")
n_prog <- nrow(sim_ordered)
col_fun <- colorRamp2(c(0, 0.16, 0.42), c("#FFF7F3", "#FB6A4A", "#67000D"))

ht <- Heatmap(
  sim_ordered,
  name = "Jaccard",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = paste0("NMF Programmes ordered by finalized refined MPs (n = ", n_prog, ")"),
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
# Add deliberate label breathing room at biological/order breakpoints.
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

draw_heatmap <- function() {
  draw(ht, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 165), "mm"))
  decorate_heatmap_body("Jaccard", {
    raw_y <- 1 - ((block_table$start + block_table$end - 1) / (2 * n_prog))
    label_y <- adjust_label_positions(raw_y, block_table$refined_mp)

    for (rr in seq_len(nrow(block_table))) {
      s <- block_table$start[rr]
      e <- block_table$end[rr]
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
        block_table$label[rr],
        x = unit(1.012, "npc"),
        y = unit(ly, "npc"),
        just = "left",
        gp = gpar(fontsize = 8.8, col = "#111111", lineheight = 0.86)
      )
    }
  })
}

pdf_path <- file.path(outdir, "figures", "refined_mp_nmf_ordered_heatmap.pdf")
png_path <- file.path(outdir, "figures", "refined_mp_nmf_ordered_heatmap.png")

cairo_pdf(pdf_path, width = 15.0, height = 9.8)
draw_heatmap()
dev.off()

png(png_path, width = 15.0, height = 9.8, units = "in", res = 450, bg = "white")
draw_heatmap()
dev.off()

summary_table <- data.frame(
  metric = c(
    "n_ordered_programs",
    "n_finalized_refined_mps",
    "strict_order_complete",
    "min_programs_per_refined_mp",
    "max_programs_per_refined_mp",
    "pdf_output",
    "png_output"
  ),
  value = c(
    n_prog,
    nrow(block_table),
    identical(sort(ordered_features), sort(finalized_mps)),
    min(block_table$n_programs),
    max(block_table$n_programs),
    pdf_path,
    png_path
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table,
          file.path(summary_dir, "refined_mp_nmf_ordered_heatmap_summary.csv"),
          row.names = FALSE)

cat("Saved:", pdf_path, "\n")
cat("Saved:", png_path, "\n")
cat("Saved: tables/refined_mp_nmf_ordered_blocks.csv\n")
cat("Saved: intermediate/refined_mp_nmf_ordered_similarity.rds\n")
