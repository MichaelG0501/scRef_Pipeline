####################
# Retained for provenance. Current entry points are
# analysis/clinical/Auto_bayesprism_bulk_deconvolution.R and
# analysis/clinical/Auto_malignant_bulk_survival.R.
# This historical workflow is not part of the current execution order.
# Its reference misclassified state-unassigned malignant cells as normal.
####################
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_Auto_bayesprism_malignant_state_deconvolution.R
#   Description: Runs BayesPrism on QC-pass OCCAMS raw counts using the full
#     scATLAS reference, models malignant cells as one cell type with the five
#     current malignant states plus Unresolved and Hybrid, and extracts tumour
#     fraction, malignant-state fractions, and tumour-only MP activity.
#   Methodology: analysis/methodology/clinical/bayesprism_malignant_state_deconvolution_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Inputs:
#   - ref_outs/EAC_Ref_merged.rds
#   - ref_outs/meta_full_epi.rds
#   - ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   - ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz
#   - ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   - ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz
#   - BayesPrism R package installed in dmtcp
#
# Outputs:
#   - intermediate/Auto_OCCAMS_bayesprism_result.rds
#   - intermediate/Auto_OCCAMS_bayesprism_tumour_expression.rds
#   - tables/Auto_OCCAMS_bayesprism_state_fractions.csv
#   - tables/Auto_OCCAMS_bayesprism_mp_scores.csv
#   - tables/Auto_OCCAMS_bayesprism_gene_coverage.csv
#   - tables/Auto_OCCAMS_bayesprism_reference_labels.csv
#   - reports/Auto_OCCAMS_bayesprism_run_report.txt
#   - updates/new_updates/summaries/Auto_OCCAMS_bayesprism_summary.csv
#
# Cache/replot: SCREF_FORCE_REBUILD=TRUE rebuilds the BayesPrism run and all
#   extraction outputs. Plot-only reruns are not applicable because this script
#   creates no figures. SCREF_RESUME_EXTRACTION=TRUE reuses the fitted RDS
#   and reference-label audit to regenerate scores, fractions, and reports.
# Run: /opt/pbs/bin/qsub analysis/clinical/legacy_Auto_bayesprism_malignant_state_deconvolution.sh
# Conda env: dmtcp
####################

library(BayesPrism)
library(data.table)
library(Matrix)
library(Seurat)
library(GSVA)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

set.seed(123)
N_CORES <- 8L

out_root <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "bayesprism")
out_dirs <- setNames(
  file.path(out_root, c("intermediate", "tables", "logs", "reports")),
  c("intermediate", "tables", "logs", "reports")
)
invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(SCREF_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(
  SCREF_REF_OUTS_DIR, "OCCAMS", "counts", "OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz"
)
qc_path <- file.path(
  SCREF_REF_OUTS_DIR, "OCCAMS", "tables", "OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv"
)
gtf_path <- file.path(
  SCREF_REF_OUTS_DIR, "OCCAMS", "source_data", "GENCODE_v19_GRCh37",
  "gencode.v19.annotation.gtf.gz"
)
merged_path <- file.path(SCREF_REF_OUTS_DIR, "EAC_Ref_merged.rds")
meta_epi_path <- SCREF_META_FULL_EPI_RDS
state_path <- SCREF_STATE_NOREG_RDS
mp_path <- SCREF_MP_GENES_RDS
grouping_path <- SCREF_MP_GROUPING_CSV

result_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_bayesprism_result.rds")
tumour_expression_path <- file.path(
  out_dirs[["intermediate"]], "Auto_OCCAMS_bayesprism_tumour_expression.rds"
)
fraction_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_bayesprism_state_fractions.csv")
mp_score_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_bayesprism_mp_scores.csv")
coverage_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_bayesprism_gene_coverage.csv")
label_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_bayesprism_reference_labels.csv")
report_path <- file.path(out_dirs[["reports"]], "Auto_OCCAMS_bayesprism_run_report.txt")
summary_path <- file.path(
  SCREF_SUMMARY_DIR, "Auto_OCCAMS_bayesprism_summary.csv"
)

force_rebuild <- identical(toupper(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")), "TRUE")
####################
resume_extraction <- identical(toupper(Sys.getenv("SCREF_RESUME_EXTRACTION", "FALSE")), "TRUE")
if (resume_extraction && force_rebuild) stop("Resume and force rebuild are mutually exclusive")
if (resume_extraction && !all(file.exists(c(result_path, label_path)))) {
  stop("Resume requires the fitted BayesPrism RDS and reference-label audit")
}
####################
output_paths <- c(
  result_path, tumour_expression_path, fraction_path, mp_score_path,
  coverage_path, label_path, report_path, summary_path
)
if (!force_rebuild && !resume_extraction && all(file.exists(output_paths))) {
  message("Complete BayesPrism outputs exist; set SCREF_FORCE_REBUILD=TRUE to rebuild.")
  quit(save = "no", status = 0)
}
if (!force_rebuild && !resume_extraction && any(file.exists(output_paths))) {
  stop("Incomplete BayesPrism cache. Set SCREF_FORCE_REBUILD=TRUE to rebuild.")
}

required_inputs <- c(
  counts_path, qc_path, gtf_path, merged_path, meta_epi_path, state_path,
  mp_path, grouping_path
)
if (any(!file.exists(required_inputs))) {
  stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

####################
# Load and validate canonical current inputs.
####################
qc <- fread(qc_path)
if (!all(c("Subject_ID", "Annotation_Compatible_Below_50_Flag") %in% colnames(qc))) {
  stop("OCCAMS QC table lacks required subject and annotation fields")
}
keep_subjects <- qc[Annotation_Compatible_Below_50_Flag == "no", Subject_ID]
if (length(keep_subjects) != 281L || anyDuplicated(keep_subjects)) {
  stop("Expected exactly 281 unique OCCAMS QC-pass subjects")
}

counts_dt <- fread(counts_path)
if (ncol(counts_dt) < 2L || colnames(counts_dt)[1L] != "Geneid") {
  stop("Unexpected OCCAMS count-matrix header")
}
gene_ids <- counts_dt[[1L]]
if (anyDuplicated(gene_ids)) stop("OCCAMS gene IDs are duplicated")
subject_ids <- colnames(counts_dt)[-1L]
missing_subjects <- setdiff(keep_subjects, subject_ids)
extra_subjects <- setdiff(subject_ids, keep_subjects)
if (length(missing_subjects) > 0L || length(extra_subjects) != 1L) {
  stop("OCCAMS count matrix must contain all 281 QC-pass subjects and exactly one excluded subject")
}
bulk_counts <- as.matrix(counts_dt[, -1L, with = FALSE])
rownames(bulk_counts) <- gene_ids
bulk_counts <- bulk_counts[, keep_subjects, drop = FALSE]
bulk_counts <- t(bulk_counts)
storage.mode(bulk_counts) <- "numeric"
if (any(!is.finite(bulk_counts)) || any(bulk_counts < 0)) {
  stop("OCCAMS counts contain invalid values")
}

gtf <- fread(gtf_path, sep = "\t", header = FALSE, comment.char = "#", quote = "")
if (ncol(gtf) < 9L) stop("GENCODE GTF has fewer than nine columns")
gene_rows <- gtf[V3 == "gene"]
gene_map <- setNames(
  sub('.*gene_name "([^"]+)".*', "\\1", gene_rows$V9),
  sub('.*gene_id "([^"]+)".*', "\\1", gene_rows$V9)
)
gene_map <- gene_map[nzchar(names(gene_map)) & nzchar(gene_map)]
if (length(gene_map) < 10000L) stop("GENCODE gene map is unexpectedly small")
gene_map_unversioned <- setNames(unname(gene_map), sub("\\..*", "", names(gene_map)))
gene_map <- c(gene_map, gene_map_unversioned[!names(gene_map_unversioned) %in% names(gene_map)])
gene_symbols_known <- unique(unname(gene_map))
canonical_symbols <- function(ids) {
  mapped <- unname(gene_map[ids])
  symbol_ids <- is.na(mapped) & ids %in% gene_symbols_known
  mapped[symbol_ids] <- ids[symbol_ids]
  mapped
}
collapse_gene_columns <- function(matrix, gene_symbols) {
  keep <- !is.na(gene_symbols) & nzchar(gene_symbols)
  if (!any(keep)) stop("No GENCODE-mapped genes remain")
  t(rowsum(t(matrix[, keep, drop = FALSE]), group = gene_symbols[keep], reorder = FALSE))
}

####################
if (!resume_extraction) {
####################
merged <- readRDS(merged_path)
meta_full_epi <- readRDS(meta_epi_path)
state_vector <- readRDS(state_path)
mp_genes <- readRDS(mp_path)
mp_grouping <- fread(grouping_path)
required_grouping_columns <- c("state", "mp", "description", "plot_label")
if (!all(required_grouping_columns %in% names(mp_grouping))) {
  stop("Current MP grouping lacks required columns")
}
if (nrow(mp_grouping) != 17L || anyDuplicated(mp_grouping$mp)) {
  stop("Current MP grouping must contain exactly 17 unique MPs")
}

if (!"celltype_update" %in% colnames(merged[[]])) {
  stop("EAC_Ref_merged.rds lacks celltype_update")
}
if (!"orig.ident" %in% colnames(merged[[]])) {
  stop("EAC_Ref_merged.rds lacks orig.ident")
}
if (!"RNA" %in% names(merged@assays)) {
  stop("EAC_Ref_merged.rds lacks RNA assay")
}
reference_counts <- GetAssayData(merged, assay = "RNA", layer = "counts")
if (!inherits(reference_counts, "dgCMatrix")) {
  stop("Reference RNA counts are not a dgCMatrix")
}

####################
# Assign one tumour cell type with explicit malignant states.
####################
primary_states <- SCREF_PRIMARY_STATE_ORDER
technical_states <- SCREF_TECHNICAL_STATE_ORDER
malignant_states <- c(primary_states, technical_states)

cell_ids <- colnames(reference_counts)
cell_type_raw <- as.character(merged$celltype_update)
names(cell_type_raw) <- cell_ids
orig_ident <- as.character(merged$orig.ident)
names(orig_ident) <- cell_ids

if (anyNA(cell_type_raw) || any(!nzchar(cell_type_raw))) {
  stop("Reference contains missing celltype_update labels")
}
exclude_cells <- cell_type_raw == "unresolved_inconsistent"
if (any(exclude_cells)) {
  reference_counts <- reference_counts[, !exclude_cells, drop = FALSE]
  cell_type_raw <- cell_type_raw[!exclude_cells]
  orig_ident <- orig_ident[!exclude_cells]
  cell_ids <- colnames(reference_counts)
}

epi_cells <- cell_ids[cell_type_raw[cell_ids] == "epithelial"]
if (is.null(names(state_vector))) stop("Current state vector has no cell names")
state_vector <- setNames(as.character(state_vector), names(state_vector))
meta_malignant_cells <- intersect(epi_cells, rownames(meta_full_epi))
malignant_cells <- intersect(meta_malignant_cells, names(state_vector))
if (length(malignant_cells) == 0L) stop("No malignant cells matched the current state vector")
excluded_without_state <- setdiff(meta_malignant_cells, malignant_cells)
if (length(excluded_without_state) > 0L) {
  message("Excluding epithelial reference cells without current state labels: ", length(excluded_without_state))
}
malignant_state_values <- state_vector[malignant_cells]
if (anyNA(malignant_state_values) || any(!malignant_state_values %in% malignant_states)) {
  stop("Current malignant state labels contain missing or unexpected values")
}

cell_type_labels <- gsub("\\.", "_", cell_type_raw)
nonmalignant_epithelial <- setdiff(epi_cells, malignant_cells)
cell_type_labels[match(nonmalignant_epithelial, cell_ids)] <- "Non_malignant_epithelial"
cell_type_labels[match(malignant_cells, cell_ids)] <- "Malignant"
cell_state_labels <- cell_type_labels
cell_state_labels[match(malignant_cells, cell_ids)] <- malignant_state_values
names(cell_type_labels) <- cell_ids
names(cell_state_labels) <- cell_ids
names(orig_ident) <- cell_ids

state_counts <- table(cell_state_labels)
if (any(state_counts < 20L)) {
  stop("Every BayesPrism cell state must contain at least 20 cells")
}
if (!all(malignant_states %in% cell_state_labels[malignant_cells])) {
  stop("Not all five malignant states plus technical states are represented")
}

reference_counts <- reference_counts[, cell_ids, drop = FALSE]
reference_by_state <- do.call(
  rbind,
  lapply(unique(cell_state_labels), function(state_name) {
    Matrix::rowSums(reference_counts[, cell_state_labels == state_name, drop = FALSE])
  })
)
rownames(reference_by_state) <- unique(cell_state_labels)
colnames(reference_by_state) <- rownames(reference_counts)
reference_by_state <- as.matrix(reference_by_state)
storage.mode(reference_by_state) <- "numeric"

if (any(!is.finite(reference_by_state)) || any(reference_by_state < 0)) {
  stop("Reference pseudobulk counts contain invalid values")
}

####################
# Align raw bulk and reference genes without a marker-only restriction.
####################
reference_by_state <- collapse_gene_columns(
  reference_by_state, canonical_symbols(colnames(reference_by_state))
)
bulk_counts <- collapse_gene_columns(
  bulk_counts, canonical_symbols(colnames(bulk_counts))
)
shared_genes <- intersect(colnames(reference_by_state), colnames(bulk_counts))
if (length(shared_genes) < 1000L) {
  stop("Fewer than 1,000 GENCODE-symbol genes shared by bulk and reference")
}
reference_by_state <- reference_by_state[, shared_genes, drop = FALSE]
bulk_counts <- bulk_counts[, shared_genes, drop = FALSE]

reference_state_labels <- rownames(reference_by_state)
reference_type_labels <- ifelse(
  reference_state_labels %in% malignant_states,
  "Malignant",
  reference_state_labels
)
if (anyDuplicated(reference_state_labels)) {
  stop("BayesPrism GEP state labels are not unique")
}

label_table <- data.frame(
  cell_type = reference_type_labels,
  cell_state = reference_state_labels,
  n_cells = as.integer(state_counts[reference_state_labels]),
  stringsAsFactors = FALSE
)
write.csv(label_table, label_path, row.names = FALSE)

####################
# Run BayesPrism with malignant-cell update enabled.
####################
prism <- new.prism(
  reference = reference_by_state,
  input.type = "GEP",
  cell.type.labels = reference_type_labels,
  cell.state.labels = reference_state_labels,
  key = "Malignant",
  mixture = bulk_counts
)
bp_result <- run.prism(prism = prism, n.cores = N_CORES)
saveRDS(bp_result, result_path)
####################
} else {
  message("Resuming extraction from fitted BayesPrism result")
  bp_result <- readRDS(result_path)
  label_table <- read.csv(label_path, check.names = FALSE)
  mp_genes <- readRDS(mp_path)
  mp_grouping <- fread(grouping_path)
  primary_states <- SCREF_PRIMARY_STATE_ORDER
  technical_states <- SCREF_TECHNICAL_STATE_ORDER
  malignant_states <- c(primary_states, technical_states)
}
if (nrow(mp_grouping) != 17L || anyDuplicated(mp_grouping$mp)) {
  stop("Expected exactly 17 unique current MPs")
}
####################

####################
# Extract fractions and tumour-only expression.
####################
state_fraction <- get.fraction(
  bp = bp_result, which.theta = "first", state.or.type = "state"
)
type_fraction_initial <- get.fraction(
  bp = bp_result, which.theta = "first", state.or.type = "type"
)
type_fraction_final <- get.fraction(
  bp = bp_result, which.theta = "final", state.or.type = "type"
)

if (!all(malignant_states %in% colnames(state_fraction))) {
  stop("BayesPrism state fractions lack current malignant states")
}
if (!"Malignant" %in% colnames(type_fraction_initial)) {
  stop("BayesPrism initial type fractions lack Malignant")
}
if (!"Malignant" %in% colnames(type_fraction_final)) {
  stop("BayesPrism final type fractions lack Malignant")
}

malignant_fraction_initial <- rowSums(state_fraction[, malignant_states, drop = FALSE])
malignant_fraction_final <- type_fraction_final[, "Malignant"]
####################
if (!setequal(rownames(state_fraction), keep_subjects) || nrow(state_fraction) != 281L ||
    anyDuplicated(rownames(state_fraction))) stop("Posterior subjects do not match QC-pass cohort")
if (!identical(rownames(state_fraction), rownames(type_fraction_final))) stop("Posterior sample order mismatch")
if (any(!is.finite(state_fraction)) || any(state_fraction < 0) ||
    max(abs(rowSums(state_fraction) - 1)) > 1e-6) stop("Invalid state posterior fractions")
if (any(!is.finite(type_fraction_final)) || any(type_fraction_final < 0) ||
    max(abs(rowSums(type_fraction_final) - 1)) > 1e-6) stop("Invalid final type posterior fractions")
####################
if (any(!is.finite(malignant_fraction_initial)) || any(malignant_fraction_initial <= 0)) {
  stop("BayesPrism returned non-positive malignant fractions")
}
fraction_df <- data.frame(
  sample_id = rownames(state_fraction),
  malignant_fraction_initial = unname(malignant_fraction_initial),
  malignant_fraction_final = unname(malignant_fraction_final),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (state_name in malignant_states) {
  fraction_df[[paste0("fraction_bulk__", state_name)]] <- state_fraction[, state_name]
  fraction_df[[paste0("fraction_tumour__", state_name)]] <-
    state_fraction[, state_name] / malignant_fraction_initial
}
write.csv(fraction_df, fraction_path, row.names = FALSE)

tumour_expression <- get.exp(
  bp = bp_result, state.or.type = "type", cell.name = "Malignant"
)
saveRDS(tumour_expression, tumour_expression_path)

####################
# Score current MPs in tumour-only expression.
####################
required_mps <- mp_grouping$mp
if (!all(required_mps %in% names(mp_genes))) {
  stop("Current MP gene object lacks one or more grouped MPs")
}
mp_genes <- mp_genes[required_mps]
mp_genes <- lapply(mp_genes, unique)

gene_symbols <- canonical_symbols(colnames(tumour_expression))
if (anyNA(gene_symbols) || any(!nzchar(gene_symbols))) {
  stop("GENCODE gene map lacks a symbol for tumour expression genes")
}
keep_symbol <- rep(TRUE, length(gene_symbols))
tumour_expression_symbol <- t(rowsum(
  t(tumour_expression[, keep_symbol, drop = FALSE]),
  group = gene_symbols[keep_symbol],
  reorder = FALSE
))
####################
if (!identical(rownames(tumour_expression_symbol), rownames(state_fraction)) ||
    any(!is.finite(tumour_expression_symbol)) || any(tumour_expression_symbol < 0)) {
  stop("Invalid tumour expression or sample order")
}
tumour_library_sizes <- rowSums(tumour_expression_symbol)
if (any(tumour_library_sizes <= 0)) stop("Tumour expression has empty samples")
tumour_expression_log <- log2(sweep(tumour_expression_symbol, 1, tumour_library_sizes, "/") * 1e6 + 1)
####################

gene_coverage <- data.frame(
  mp = names(mp_genes),
  genes_defined = vapply(mp_genes, length, integer(1)),
  genes_present = vapply(mp_genes, function(genes) sum(genes %in% colnames(tumour_expression_log)), integer(1)),
  stringsAsFactors = FALSE
)
if (any(gene_coverage$genes_present < 5L)) {
  stop("At least one current MP has fewer than five tumour-only genes")
}
write.csv(gene_coverage, coverage_path, row.names = FALSE)

mp_scores <- gsva(
  t(tumour_expression_log),
  mp_genes,
  method = "gsva",
  kcdf = "Gaussian",
  verbose = FALSE
)
mp_score_df <- data.frame(
  sample_id = colnames(mp_scores),
  as.data.frame(t(mp_scores), check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)
####################
if (!identical(colnames(mp_scores), rownames(state_fraction)) ||
    !identical(rownames(mp_scores), required_mps) || any(!is.finite(mp_scores))) {
  stop("GSVA output must contain 17 finite MP scores for each retained subject")
}
####################
write.csv(mp_score_df, mp_score_path, row.names = FALSE)

run_summary <- data.frame(
  cohort = "OCCAMS",
  subjects = nrow(state_fraction),
  reference_states = nrow(label_table),
  malignant_cells = sum(label_table$n_cells[label_table$cell_type == "Malignant"]),
  shared_genes = ncol(tumour_expression),
  current_mps = nrow(mp_scores),
  primary_malignant_states = length(primary_states),
  technical_malignant_states = length(technical_states),
  bayesprism_cores = N_CORES,
  stringsAsFactors = FALSE
)
write.csv(run_summary, summary_path, row.names = FALSE)

report_lines <- c(
  "BayesPrism malignant-state deconvolution",
  "result=PASS",
  paste0("subjects=", nrow(state_fraction)),
  paste0("reference_states=", nrow(label_table)),
  paste0("malignant_cells=", run_summary$malignant_cells),
  paste0("inference_genes=", ncol(tumour_expression)),
  paste0("current_mps=", nrow(mp_scores)),
  paste0("primary_malignant_states=", paste(primary_states, collapse = " | ")),
  paste0("technical_malignant_states=", paste(technical_states, collapse = " | ")),
  "state_fraction_source=initial BayesPrism state posterior",
  "malignant_fraction_source=initial and final BayesPrism type posterior",
  "mp_activity_source=initial BayesPrism malignant-type expression",
  "mp_activity_transform=log2(CPM + 1), Gaussian GSVA with genes in rows",
  paste0("malignant_fraction_final_min_median_max=", paste(signif(c(min(malignant_fraction_final), median(malignant_fraction_final), max(malignant_fraction_final)), 5), collapse = " | ")),
  paste0("min_mp_gene_coverage=", min(gene_coverage$genes_present / gene_coverage$genes_defined)),
  "validation=281 unique QC-pass subjects; finite 17-MP scores; fractions sum to one",
  "limitation=BayesPrism fractions describe RNA/read contribution, not validated histological cell counts",
  "limitation=state-less epithelial cells were retained as Non_malignant_epithelial in the fitted reference; tumour estimates require biological validation",
  paste0("outputs=", out_root)
)
writeLines(report_lines, report_path)
message("BayesPrism malignant-state deconvolution complete: ", out_root)
