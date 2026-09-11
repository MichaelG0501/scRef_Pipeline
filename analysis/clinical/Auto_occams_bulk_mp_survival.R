####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/Auto_occams_bulk_mp_survival.R
#   Description: Builds QC-pass OCCAMS gene-symbol logCPM, scores the current
#     17 centred refined MPs and five non-cell-cycle state unions with GSVA,
#     resolves subject-level metadata, and fits overall-survival Cox models.
#   Methodology: analysis/methodology/clinical/Auto_occams_bulk_clinical_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Inputs:
#   - ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz
#   - ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   - ref_outs/OCCAMS/source_data/OCCAMS_metadata.csv
#   - ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
# Outputs:
#   intermediate/:
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_mp_scores.rds
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_state_scores.rds
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_dge_scores.rds
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_survival_model_data.rds
#   tables/:
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_metadata_field_audit.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_metadata_conflicts.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_gene_set_coverage.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_survival_cox_results.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_survival_optimal_cut_results.csv
#   figures/:
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_survival_volcano.pdf
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_survival_km_selected.pdf
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_survival_optimal_cut_volcano_km.pdf
#     - ref_outs/OCCAMS/clinical/figures/Auto_OCCAMS_survival_optimal_cut_km_selected.pdf
#   logs/:
#     - ref_outs/OCCAMS/clinical/logs/Auto_occams_bulk_mp_survival.stdout.log
#     - ref_outs/OCCAMS/clinical/logs/Auto_occams_bulk_mp_survival.stderr.log
#   reports/ and summaries/:
#     - ref_outs/OCCAMS/clinical/reports/Auto_OCCAMS_survival_run_summary.txt
#     - updates/new_updates/summaries/Auto_OCCAMS_survival_summary.csv
# Cache/replot behavior:
#   SCREF_FORCE_REBUILD=TRUE rebuilds normalized expression and GSVA caches.
#   Existing valid caches are reused otherwise. SCREF_REPLOT_ONLY=TRUE requires
#   persistent model data/results and redraws figures without counts or GSVA.
# Run command:
#   /opt/pbs/bin/qsub analysis/clinical/Auto_occams_bulk_mp_survival.sh
# Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(GSVA)
library(survival)
library(survminer)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

project_dir <- SCREF_PROJECT_DIR
out_root <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "clinical")
out_dirs <- setNames(file.path(out_root, c("intermediate", "tables", "figures", "logs", "reports")),
                     c("intermediate", "tables", "figures", "logs", "reports"))
invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(SCREF_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "counts", "OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz")
qc_path <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "tables", "OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv")
metadata_path <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "source_data", "OCCAMS_metadata.csv")
gtf_path <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "source_data", "GENCODE_v19_GRCh37", "gencode.v19.annotation.gtf.gz")
mp_path <- SCREF_MP_GENES_RDS
grouping_path <- SCREF_MP_GROUPING_CSV

logcpm_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds")
mp_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_mp_scores.rds")
state_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_state_scores.rds")
dge_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_dge_scores.rds")
model_data_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_survival_model_data.rds")
subject_metadata_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_qc_pass_subject_metadata.csv")
field_audit_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_metadata_field_audit.csv")
conflict_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_metadata_conflicts.csv")
coverage_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_gene_set_coverage.csv")
cox_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_survival_cox_results.csv")
optimal_cut_results_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_survival_optimal_cut_results.csv")
volcano_path <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_survival_volcano.pdf")
km_path <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_survival_km_selected.pdf")
optimal_cut_path <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_survival_optimal_cut_volcano_km.pdf")
optimal_cut_km_path <- file.path(out_dirs[["figures"]], "Auto_OCCAMS_survival_optimal_cut_km_selected.pdf")
report_path <- file.path(out_dirs[["reports"]], "Auto_OCCAMS_survival_run_summary.txt")
summary_path <- file.path(SCREF_SUMMARY_DIR, "Auto_OCCAMS_survival_summary.csv")

force_rebuild <- identical(toupper(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")), "TRUE")
replot_only <- identical(toupper(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")), "TRUE")

required_inputs <- c(counts_path, qc_path, metadata_path, gtf_path, mp_path, grouping_path)
if (!replot_only && any(!file.exists(required_inputs))) {
  stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

clean_value <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | tolower(x) %in% c("", "na", "n/a", "unknown", "not_recorded", "not recorded")
  x[bad] <- NA_character_
  x
}

consensus_value <- function(x) {
  values <- unique(stats::na.omit(clean_value(x)))
  if (length(values) == 1) values[[1]] else NA_character_
}

normalise_t_stage <- function(x) {
  x <- tolower(clean_value(x))
  out <- rep(NA_character_, length(x))
  out[grepl("^t0$|no_evidence", x)] <- "T0"
  out[grepl("^t1|lamina|mucosae|submucosa", x)] <- "T1"
  out[grepl("^t2$|muscularis_propria", x)] <- "T2"
  out[grepl("^t3$|adventitia", x)] <- "T3"
  out[grepl("^t4|adjacent_structures|pleura|pericardium|diaphragm", x)] <- "T4"
  out
}

normalise_site <- function(x) {
  x <- tolower(clean_value(x))
  out <- rep(NA_character_, length(x))
  out[grepl("goj|siewert", x)] <- "GOJ"
  out[grepl("oesoph", x)] <- "Oesophagus"
  out[grepl("lower third|middle third|upper third", x)] <- "Oesophagus"
  out[grepl("gastric|stomach|cardia|body", x)] <- "Stomach"
  out
}

normalise_m_stage <- function(x) {
  x <- toupper(clean_value(x))
  ifelse(x %in% c("M0", "M1"), x, NA_character_)
}

normalise_n_stage <- function(x) {
  x <- toupper(clean_value(x))
  ifelse(x %in% c("N0", "N1", "N2", "N3"), x, NA_character_)
}

collapse_occams_metadata <- function(metadata, keep_ids) {
  metadata <- metadata[metadata[["SHA IDs"]] %in% keep_ids, , drop = FALSE]
  split_rows <- split(metadata, metadata[["SHA IDs"]])
  raw_cols <- setdiff(colnames(metadata), "SHA IDs")

  audit <- bind_rows(lapply(raw_cols, function(column_name) {
    subject_values <- lapply(split_rows, function(d) unique(stats::na.omit(clean_value(d[[column_name]]))))
    data.frame(
      source_column = column_name,
      n_subjects_available = sum(lengths(subject_values) > 0),
      n_subjects_discordant = sum(lengths(subject_values) > 1),
      n_unique_nonmissing_values = length(unique(unlist(subject_values, use.names = FALSE))),
      stringsAsFactors = FALSE
    )
  }))

  conflicts <- bind_rows(lapply(names(split_rows), function(subject_id) {
    d <- split_rows[[subject_id]]
    bind_rows(lapply(raw_cols, function(column_name) {
      values <- unique(stats::na.omit(clean_value(d[[column_name]])))
      if (length(values) <= 1) return(NULL)
      data.frame(subject_id = subject_id, source_column = column_name,
                 values = paste(sort(values), collapse = " | "), stringsAsFactors = FALSE)
    }))
  }))

  collapsed <- bind_rows(lapply(names(split_rows), function(subject_id) {
    d <- split_rows[[subject_id]]
    cv <- function(column_name) consensus_value(d[[column_name]])
    cvn <- function(column_name, normaliser) consensus_value(normaliser(d[[column_name]]))
    deceased_days <- suppressWarnings(as.numeric(cv("Deceased Survival Days")))
    last_known_days <- suppressWarnings(as.numeric(cv("Last Known Survival Days")))
    os_event <- if (!is.na(deceased_days)) 1L else if (!is.na(last_known_days)) 0L else NA_integer_
    os_days <- if (!is.na(deceased_days)) deceased_days else last_known_days

    recurrence_raw <- cv("Has Original Disease Reoccurred")
    recurrence_tokens <- if (is.na(recurrence_raw)) {
      character(0)
    } else {
      tolower(trimws(unlist(strsplit(recurrence_raw, ",", fixed = TRUE))))
    }
    recurrence <- if (any(recurrence_tokens == "yes")) "Yes" else if (any(recurrence_tokens == "no")) "No" else NA_character_
    mandard <- cv("RP Mandard Score For Response")
    mandard_binary <- if (!is.na(mandard) && mandard %in% c("TRG1", "TRG2")) "TRG1-2" else if (!is.na(mandard) && mandard %in% c("TRG3", "TRG4", "TRG5")) "TRG3-5" else NA_character_
    treatment_response <- cv("TR Response")
    response_binary <- if (!is.na(treatment_response) && treatment_response %in% c("CR", "PR")) "CR/PR" else if (!is.na(treatment_response) && treatment_response %in% c("SD", "PD")) "SD/PD" else NA_character_
    grade <- cv("RP Tumour Grading Differentiation Status")
    grade_group <- if (!is.na(grade) && grade %in% c("well", "moderate_to_well", "moderate")) "Well/moderate" else if (!is.na(grade) && grade %in% c("moderate_to_poor", "poor")) "Moderate/poor" else NA_character_
    performance <- cv("TP Performance Status")
    performance_group <- if (!is.na(performance) && performance == "0") "0" else if (!is.na(performance) && performance %in% c("1", "2")) "1-2" else NA_character_
    comorbidities <- cv("TP Comorbidities")
    comorbidity_group <- if (is.na(comorbidities)) NA_character_ else if (tolower(comorbidities) == "no comorbidities") "None recorded" else "Any recorded"
    barrett_raw <- cv("RP Baretts Adjacent To Tumour Microscopic Dysplasia")
    barrett_dysplasia <- if (is.na(barrett_raw)) NA_character_ else if (grepl("low grade|high grade", tolower(barrett_raw))) "Dysplasia present" else if (tolower(barrett_raw) == "no dysplasia") "No dysplasia" else NA_character_
    positive_nodes <- suppressWarnings(as.numeric(cv("RP NStage Number Of Positive Nodes")))
    positive_node_group <- if (is.na(positive_nodes)) NA_character_ else if (positive_nodes == 0) "0" else if (positive_nodes <= 3) "1-3" else ">=4"
    procedure <- cv("ST Procedure")
    procedure_group <- if (is.na(procedure)) NA_character_ else if (grepl("Ivor-Lewis", procedure, fixed = TRUE)) "Ivor-Lewis oesophagectomy" else if (grepl("Oesophagectomy|MIO", procedure)) "Other oesophagectomy" else if (grepl("Gastrectomy", procedure)) "Gastrectomy" else "Other operation"

    data.frame(
      subject_id = subject_id,
      metadata_row_count = nrow(d),
      gender = cv("DI Patient Gender"),
      age_at_diagnosis = suppressWarnings(as.numeric(cv("DI Age At Diagnosis"))),
      age_group = ifelse(suppressWarnings(as.numeric(cv("DI Age At Diagnosis"))) > 60, ">60", "<=60"),
      diagnostic_histology = cv("RD Histology"),
      os_days = os_days,
      os_event = os_event,
      deceased_survival_days = deceased_days,
      last_known_survival_days = last_known_days,
      fe_end_point = cv("FE End Point"),
      ep_end_point = cv("EP End Point"),
      diagnostic_combined_tumour_site = cvn("RD Combined Tumour Site", normalise_site),
      diagnostic_oesophagus_site = cvn("RD Oesophagus Site", normalise_site),
      pathology_location = cvn("RP Location", normalise_site),
      pretreatment_combined_tumour_site = cvn("PS Pre Treatment Combined Tumour Site", normalise_site),
      pretreatment_oesophagus_site = cvn("PS Pre Treatment Oesophagus Site", normalise_site),
      diagnostic_siewert_classification = cv("RD Siewert Classification"),
      pretreatment_siewert_classification = cv("PS Siewert Classification"),
      pathology_siewert_classification = cv("RP Siewert Classification"),
      pretreatment_t_stage = cvn("PS TStage Primary Tumour Final Pretreatment Staging", normalise_t_stage),
      pretreatment_n_stage_tnm7 = cvn("PS NStage Primary Tumour Final Pretreatment Staging TNM7", normalise_n_stage),
      pretreatment_n_stage_tnm6 = cvn("PS NStage Primary Tumour Final Pretreatment Staging TNM6", normalise_n_stage),
      pretreatment_m_stage = cvn("PS MStage Primary Tumour Final Pretreatment Staging", normalise_m_stage),
      treatment_plan_performance_status = performance,
      treatment_plan_performance_status_group = performance_group,
      treatment_plan_intent = cv("TP Treatment Intent"),
      treatment_plan_curative_modality = cv("TP Curative Treatment Modality"),
      surgery_main = cv("ST Main Surgery"),
      surgery_procedure_group = procedure_group,
      pathology_t_stage = cvn("RP TStage Primary Tumour", normalise_t_stage),
      pathology_n_stage_tnm7 = cvn("RP Nstage RP TNM7", normalise_n_stage),
      pathology_n_stage_tnm6 = cvn("RP Nstage RP TNM6", normalise_n_stage),
      pathology_m_stage = cvn("RP MStage Distant Metastasis", normalise_m_stage),
      pathology_positive_node_count = positive_nodes,
      pathology_positive_node_group = positive_node_group,
      pathology_grade = grade,
      pathology_grade_group = grade_group,
      pathology_signet_ring_cells = cv("RP Signet Ring Cells Present"),
      pathology_barrett_adjacent_dysplasia = barrett_dysplasia,
      pathology_neoadjuvant_therapy = cv("RP History Of Neo Adjuvant Therapy"),
      pathology_mandard_score = mandard,
      pathology_mandard_response_group = mandard_binary,
      treatment_response = treatment_response,
      treatment_response_group = response_binary,
      recurrence_ever_recorded = recurrence,
      treatment_plan_comorbidities = comorbidities,
      treatment_plan_comorbidity_group = comorbidity_group,
      stringsAsFactors = FALSE
    )
  }))

  list(metadata = collapsed, audit = audit, conflicts = conflicts)
}

extract_gene_map <- function(path) {
  command <- paste("zgrep -P", shQuote("\\tgene\\t"), shQuote(path))
  gtf <- data.table::fread(cmd = command, sep = "\t", header = FALSE, select = c(3, 9),
                           showProgress = FALSE)
  data.table::setnames(gtf, c("feature", "attributes"))
  gtf <- gtf[feature == "gene"]
  gtf[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
  gtf[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attributes)]
  unique(gtf[, .(gene_id, gene_name)])
}

run_gsva <- function(expression, gene_sets) {
  filtered <- lapply(gene_sets, function(x) intersect(unique(as.character(x)), rownames(expression)))
  filtered <- filtered[lengths(filtered) >= 5]
  if (length(filtered) == 0) stop("No gene set retained at >=5 genes")
  scores <- GSVA::gsva(expression, filtered, method = "gsva", kcdf = "Gaussian",
                       min.sz = 5, max.sz = Inf, verbose = FALSE)
  list(scores = scores, filtered_sets = filtered)
}

build_caches <- function() {
  qc <- read.csv(qc_path, stringsAsFactors = FALSE, check.names = FALSE)
  keep_ids <- qc$Subject_ID[tolower(qc$Annotation_Compatible_Below_50_Flag) == "no"]
  if (length(keep_ids) != 281 || anyDuplicated(keep_ids)) {
    stop("Expected exactly 281 unique QC-pass subjects; found ", length(unique(keep_ids)))
  }

  metadata_raw <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE,
                           fileEncoding = "latin1")
  collapsed <- collapse_occams_metadata(metadata_raw, keep_ids)
  if (!setequal(collapsed$metadata$subject_id, keep_ids)) stop("QC-pass metadata coverage is incomplete")
  write.csv(collapsed$metadata, subject_metadata_path, row.names = FALSE, na = "")
  write.csv(collapsed$audit, field_audit_path, row.names = FALSE, na = "")
  write.csv(collapsed$conflicts, conflict_path, row.names = FALSE, na = "")

  counts_dt <- data.table::fread(counts_path, check.names = FALSE, showProgress = TRUE)
  gene_ids <- counts_dt[[1]]
  counts <- as.matrix(counts_dt[, -1, with = FALSE])
  rownames(counts) <- gene_ids
  storage.mode(counts) <- "numeric"
  if (!setequal(colnames(counts), qc$Subject_ID) || ncol(counts) != 282) {
    stop("Count matrix does not match the validated 282-subject QC table")
  }
  counts <- counts[, keep_ids, drop = FALSE]

  gene_map <- extract_gene_map(gtf_path)
  symbol <- gene_map$gene_name[match(rownames(counts), gene_map$gene_id)]
  keep_gene <- !is.na(symbol) & nzchar(symbol)
  counts_by_symbol <- rowsum(counts[keep_gene, , drop = FALSE], symbol[keep_gene], reorder = FALSE)
  min_samples <- ceiling(0.10 * ncol(counts_by_symbol))
  expressed <- rowSums(edgeR::cpm(counts_by_symbol) >= 1) >= min_samples
  counts_by_symbol <- counts_by_symbol[expressed, , drop = FALSE]
  dge <- edgeR::DGEList(counts = counts_by_symbol)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  saveRDS(logcpm, logcpm_path, compress = "xz")

  mp_genes <- readRDS(mp_path)
  grouping <- read.csv(grouping_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_genes <- mp_genes[grouping$mp]
  if (length(mp_genes) != 17 || any(vapply(mp_genes, is.null, logical(1)))) {
    stop("Current centred refined MP input did not resolve to the expected 17 MPs")
  }
  state_grouping <- grouping[grouping$state != "Cell cycle", , drop = FALSE]
  state_genes <- lapply(split(state_grouping$mp, state_grouping$state), function(mps) {
    unique(unlist(mp_genes[mps], use.names = FALSE))
  })
  
  ranked_markers_path <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv")
  if (!file.exists(ranked_markers_path)) stop("Missing ranked markers for DGE")
  ranked_markers <- read.csv(ranked_markers_path, stringsAsFactors = FALSE)
  state_dges <- lapply(split(ranked_markers, ranked_markers$state), function(df) {
    head(df$gene, 20)
  })
  names(state_dges) <- paste0(names(state_dges), " (DGE)")

  mp_result <- run_gsva(logcpm, mp_genes)
  state_result <- run_gsva(logcpm, state_genes)
  dge_result <- run_gsva(logcpm, state_dges)
  
  mp_scores <- t(mp_result$scores)
  state_scores <- t(state_result$scores)
  dge_scores <- t(dge_result$scores)
  
  saveRDS(mp_scores, mp_score_path, compress = "xz")
  saveRDS(state_scores, state_score_path, compress = "xz")
  saveRDS(dge_scores, dge_score_path, compress = "xz")

  mp_matched <- lengths(mp_result$filtered_sets)[names(mp_genes)]
  state_matched <- lengths(state_result$filtered_sets)[names(state_genes)]
  dge_matched <- lengths(dge_result$filtered_sets)[names(state_dges)]
  mp_matched[is.na(mp_matched)] <- 0L
  state_matched[is.na(state_matched)] <- 0L
  dge_matched[is.na(dge_matched)] <- 0L
  
  coverage <- bind_rows(
    data.frame(feature_type = "MP", feature = names(mp_genes),
               source_genes = lengths(mp_genes), matched_genes = unname(mp_matched), stringsAsFactors = FALSE),
    data.frame(feature_type = "State", feature = names(state_genes),
               source_genes = lengths(state_genes), matched_genes = unname(state_matched), stringsAsFactors = FALSE),
    data.frame(feature_type = "State DGE", feature = names(state_dges),
               source_genes = lengths(state_dges), matched_genes = unname(dge_matched), stringsAsFactors = FALSE)
  ) %>% mutate(coverage_percent = 100 * matched_genes / source_genes)
  write.csv(coverage, coverage_path, row.names = FALSE)
  list(metadata = collapsed$metadata, mp_scores = mp_scores, state_scores = state_scores, dge_scores = dge_scores,
       n_expression_genes = nrow(logcpm))
}

make_model_data <- function(metadata, mp_scores, state_scores, dge_scores) {
  score_df <- cbind(as.data.frame(mp_scores, check.names = FALSE),
                    as.data.frame(state_scores, check.names = FALSE),
                    as.data.frame(dge_scores, check.names = FALSE))
  score_df$subject_id <- rownames(score_df)
  metadata %>% inner_join(score_df, by = "subject_id")
}

fit_one_cox <- function(data, feature, feature_type, split_method) {
  d <- data[!is.na(data$os_days) & !is.na(data$os_event) & !is.na(data[[feature]]) & data$os_days > 0, , drop = FALSE]
  if (nrow(d) < 20 || sum(d$os_event) < 10 || stats::sd(d[[feature]]) == 0) return(NULL)
  score_sd <- stats::sd(d[[feature]])
  threshold_low <- NA_real_
  threshold_high <- NA_real_
  if (split_method == "continuous_per_sd") {
    d$model_value <- as.numeric(scale(d[[feature]]))
  } else if (split_method == "median") {
    threshold_low <- stats::median(d[[feature]])
    threshold_high <- threshold_low
    d$model_value <- factor(ifelse(d[[feature]] > threshold_high, "High", "Low"), levels = c("Low", "High"))
  } else {
    quantiles <- stats::quantile(d[[feature]], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    threshold_low <- quantiles[[1]]
    threshold_high <- quantiles[[2]]
    d <- d[d[[feature]] <= threshold_low | d[[feature]] >= threshold_high, , drop = FALSE]
    d$model_value <- factor(ifelse(d[[feature]] >= threshold_high, "High", "Low"), levels = c("Low", "High"))
  }
  fit <- try(survival::coxph(survival::Surv(os_days, os_event) ~ model_value, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  ss <- summary(fit)
  ci <- ss$conf.int[1, ]
  data.frame(
    cohort = "OCCAMS QC-pass", endpoint = "Overall survival", time_field = "deceased_survival_days else last_known_survival_days",
    event_field = "deceased_survival_days present", feature_type = feature_type, feature = feature,
    split_method = split_method, score_scaling = ifelse(split_method == "continuous_per_sd", "one cohort SD", "High versus Low"),
    covariates = "none", n = fit$n, events = fit$nevent, score_sd = score_sd,
    threshold_low = threshold_low, threshold_high = threshold_high,
    hazard_ratio = unname(ci[["exp(coef)"]]), ci_lower_95 = unname(ci[["lower .95"]]),
    ci_upper_95 = unname(ci[["upper .95"]]), p_value = ss$coefficients[1, "Pr(>|z|)"], stringsAsFactors = FALSE
  )
}

fit_all_cox <- function(model_data, mp_features, state_features, dge_features) {
  features <- c(mp_features, state_features, dge_features)
  types <- c(rep("MP", length(mp_features)), rep("State", length(state_features)), rep("State DGE", length(dge_features)))
  splits <- c("continuous_per_sd", "median", "q1_vs_q4")
  results <- bind_rows(lapply(seq_along(features), function(i) {
    bind_rows(lapply(splits, function(split_method) fit_one_cox(model_data, features[[i]], types[[i]], split_method)))
  }))
  results %>% group_by(feature_type, split_method) %>% mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
}

plot_volcanoes <- function(cox_results) {
  plot_data <- cox_results %>% filter(!is.na(feature_type)) %>%
    mutate(log2_hr = log2(hazard_ratio), neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
           significant_fdr = p_adj_bh < 0.05,
           label = ifelse(feature_type == "MP", paste(feature, SCREF_MP_DESCRIPTIONS[feature]), feature))
           
  make_volcano_panel <- function(d, split_name) {
    if (nrow(d) == 0) return(NULL)
    ggplot(d, aes(log2_hr, neg_log10_p)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
      geom_point(aes(colour = significant_fdr), size = 2.8, alpha = 0.9) +
      ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 100,
                               box.padding = 0.5, point.padding = 0.2, fontface = "bold") +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick3"), guide = "none") +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Split: ", split_name), x = "log2(HR)", y = "-log10(p)")
  }
  
  grDevices::cairo_pdf(volcano_path, width = 18, height = 8, onefile = TRUE)
  
  for (ft in unique(plot_data$feature_type)) {
    panels <- list()
    for (split_name in c("continuous_per_sd", "median", "q1_vs_q4")) {
      d <- plot_data %>% filter(feature_type == ft, split_method == split_name)
      panels[[split_name]] <- make_volcano_panel(d, split_name)
    }
    
    for (sm in names(panels)) {
      if (is.null(panels[[sm]])) {
        panels[[sm]] <- ggplot() + theme_void() + annotate("text", x=0, y=0, label="No model available") + labs(title = sm)
      }
    }
    
    page <- gridExtra::arrangeGrob(
      panels[["continuous_per_sd"]], panels[["median"]], panels[["q1_vs_q4"]],
      ncol = 3,
      top = grid::textGrob(paste("OCCAMS QC-pass overall survival:", ft), gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  dev.off()
}

plot_selected_km <- function(model_data, cox_results) {
  best_p <- cox_results %>%
    filter(!is.na(feature_type), split_method %in% c("median", "q1_vs_q4")) %>%
    group_by(feature_type, feature) %>%
    slice_min(p_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(p_value)
  
  bad_features <- best_p %>% filter(hazard_ratio > 1) %>% head(5) %>% pull(feature)
  good_features <- best_p %>% filter(hazard_ratio < 1) %>% head(5) %>% pull(feature)
  
  selected <- unique(c(bad_features, good_features))
  selection_rule <- "Top 5 bad and top 5 good survival features by absolute minimum p-value across categorical splits (median, q1_vs_q4)"
  
  grDevices::cairo_pdf(km_path, width = 10, height = 8, onefile = TRUE)
  for (feature in selected) {
    d <- model_data[!is.na(model_data$os_days) & !is.na(model_data$os_event) & !is.na(model_data[[feature]]) & model_data$os_days > 0, , drop = FALSE]
    
    # Use the split method that yielded the best categorical p-value
    best_split <- best_p %>% filter(feature == !!feature) %>% pull(split_method) %>% .[1]
    
    if (best_split == "q1_vs_q4") {
      quantiles <- stats::quantile(d[[feature]], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      d <- d[d[[feature]] <= quantiles[1] | d[[feature]] >= quantiles[2], , drop = FALSE]
      d$score_group <- factor(ifelse(d[[feature]] >= quantiles[2], "High", "Low"), levels = c("Low", "High"))
      plot_title_suffix <- "(Q1 vs Q4)"
    } else {
      threshold <- median(d[[feature]])
      d$score_group <- factor(ifelse(d[[feature]] > threshold, "High", "Low"), levels = c("Low", "High"))
      plot_title_suffix <- "(median split)"
    }
    
    fit <- survival::survfit(survival::Surv(os_days, os_event) ~ score_group, data = d)
    
    ftype <- best_p %>% filter(feature == !!feature) %>% pull(feature_type) %>% .[1]
    title_label <- feature
    if (ftype == "MP" && feature %in% names(SCREF_MP_DESCRIPTIONS)) {
      title_label <- paste(feature, SCREF_MP_DESCRIPTIONS[feature])
    }
    
    p <- survminer::ggsurvplot(
      fit, data = d, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
      palette = c("#377EB8", "#E41A1C"), 
      xlab = "Time (Days)", ylab = "Overall Survival Probability",
      title = paste(title_label, "\n", plot_title_suffix),
      legend.labs = c("Low", "High"),
      ggtheme = theme_minimal()
    )
    print(p)
  }
  dev.off()
  data.frame(feature = selected, selection_rule = selection_rule, stringsAsFactors = FALSE)
}

plot_optimal_cut_volcano_km <- function(model_data, feature_types, output_path, km_output_path) {
  d_all <- model_data
  specs <- data.frame(feature = unlist(feature_types, use.names = FALSE),
                      feature_type = rep(names(feature_types), lengths(feature_types)),
                      stringsAsFactors = FALSE)
  results <- dplyr::bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    feature <- specs$feature[i]
    d <- d_all[!is.na(d_all$os_days) & !is.na(d_all$os_event) &
                 !is.na(d_all[[feature]]) & d_all$os_days > 0, , drop = FALSE]
    if (nrow(d) < 50) return(NULL)
    cuts <- unique(stats::quantile(d[[feature]], seq(0.20, 0.80, by = 0.05),
                                   na.rm = TRUE, names = FALSE))
    candidates <- dplyr::bind_rows(lapply(cuts, function(cutpoint) {
      grp <- ifelse(d[[feature]] >= cutpoint, 1, 0)
      if (sum(grp) == 0 || sum(grp) == nrow(d)) return(NULL)
      fit <- try(survival::coxph(survival::Surv(os_days, os_event) ~ grp, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      ss <- summary(fit)
      data.frame(cutpoint = cutpoint, hazard_ratio = ss$conf.int[1, "exp(coef)"],
                 p_value = ss$coefficients[1, "Pr(>|z|)"], n = fit$n,
                 events = fit$nevent, stringsAsFactors = FALSE)
    }))
    if (nrow(candidates) == 0) return(NULL)
    best <- candidates[order(candidates$p_value, candidates$cutpoint), ][1, ]
    cbind(specs[i, , drop = FALSE], best)
  }))
  if (nrow(results) == 0) stop("No valid optimal-cut survival models available")
  results <- results %>% mutate(
    log2_hr = log2(hazard_ratio),
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    nominal_p_lt_0_05 = p_value < 0.05,
    label = ifelse(feature_type == "MP" & feature %in% names(SCREF_MP_DESCRIPTIONS),
                   paste(feature, SCREF_MP_DESCRIPTIONS[feature]), feature)
  )
  make_panel <- function(d, title) {
    ggplot(d, aes(log2_hr, neg_log10_p)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
      geom_point(aes(colour = p_value < 0.05), size = 2.8, alpha = 0.9) +
      ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 100) +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick3"),
                           guide = "none") +
      theme_minimal(base_size = 12) +
      labs(title = title, x = "log2(HR)", y = "-log10(p-value)")
  }
  panels <- lapply(names(feature_types), function(ft) {
    make_panel(results[results$feature_type == ft, , drop = FALSE],
               paste0(ft, " (Optimal Cut)"))
  })
  grDevices::cairo_pdf(output_path, width = 18, height = 8, onefile = TRUE)
  grid::grid.newpage()
  grid::grid.draw(gridExtra::arrangeGrob(grobs = panels, ncol = 3,
    top = grid::textGrob("OCCAMS Survival (Optimal Cutpoint; univariable Cox model)",
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))))

  dev.off()

  grDevices::cairo_pdf(km_output_path, width = 10, height = 8, onefile = TRUE)
  first_km_page <- TRUE
  selected <- results %>% arrange(p_value) %>%
    group_by(hazard_ratio > 1) %>% slice_head(n = 5) %>% ungroup() %>%
    pull(feature) %>% unique()
  for (feature in selected) {
    d <- d_all[!is.na(d_all$os_days) & !is.na(d_all$os_event) &
                 !is.na(d_all[[feature]]) & d_all$os_days > 0, , drop = FALSE]
    result_index <- match(feature, results[["feature"]])
    cutpoint <- results[["cutpoint"]][result_index]
    selected_p <- results[["p_value"]][result_index]
    d$grp <- ifelse(d[[feature]] >= cutpoint, 1, 0)
    fit <- survival::survfit(survival::Surv(os_days, os_event) ~ grp, data = d)
    title_label <- if (feature %in% names(SCREF_MP_DESCRIPTIONS))
      paste(feature, SCREF_MP_DESCRIPTIONS[feature]) else feature
    km_plot <- survminer::ggsurvplot(
      fit, data = d, risk.table = TRUE, pval = paste0("Optimal-cut Cox Wald p = ", format.pval(selected_p, digits = 3)), conf.int = FALSE,
      palette = c("#377EB8", "#E41A1C"), xlab = "Time (Days)",
      ylab = "Overall Survival Probability",
      title = paste0(title_label, "\nOptimal cut = ", signif(cutpoint, 4)),
      legend.labs = c("Low", "High"), ggtheme = theme_minimal()
    )
    print(km_plot, newpage = !first_km_page)
    first_km_page <- FALSE
  }
  dev.off()
  results
}

if (replot_only) {
  required_cache <- c(model_data_path, cox_path, mp_score_path, state_score_path, dge_score_path, logcpm_path)
  if (any(!file.exists(required_cache))) stop("SCREF_REPLOT_ONLY requires: ", paste(required_cache, collapse = ", "))
  model_data <- readRDS(model_data_path)
  cox_results <- read.csv(cox_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_scores <- readRDS(mp_score_path)
  state_scores <- readRDS(state_score_path)
  dge_scores <- readRDS(dge_score_path)
  n_expression_genes <- nrow(readRDS(logcpm_path))
} else {
  cache_paths <- c(logcpm_path, mp_score_path, state_score_path, subject_metadata_path)
  if (force_rebuild || any(!file.exists(cache_paths))) {
    built <- build_caches()
    metadata <- built$metadata
    mp_scores <- built$mp_scores
    state_scores <- built$state_scores
    dge_scores <- built$dge_scores
    n_expression_genes <- built$n_expression_genes
  } else {
    metadata <- read.csv(subject_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
    mp_scores <- readRDS(mp_score_path)
    state_scores <- readRDS(state_score_path)
    dge_scores <- readRDS(dge_score_path)
    n_expression_genes <- nrow(readRDS(logcpm_path))
  }
  if (nrow(metadata) != 281 || nrow(mp_scores) != 281 || nrow(state_scores) != 281 || nrow(dge_scores) != 281) {
    stop("All downstream objects must contain exactly 281 QC-pass subjects")
  }
  model_data <- make_model_data(metadata, mp_scores, state_scores, dge_scores)
  if (nrow(model_data) != 281) stop("Model data join did not retain exactly 281 subjects")
  saveRDS(model_data, model_data_path, compress = "xz")
  cox_results <- fit_all_cox(model_data, colnames(mp_scores), colnames(state_scores), colnames(dge_scores))
  write.csv(cox_results, cox_path, row.names = FALSE, na = "")
}

plot_volcanoes(cox_results)
km_selection <- plot_selected_km(model_data, cox_results)
optimal_cut_results <- plot_optimal_cut_volcano_km(model_data, list(MP = colnames(mp_scores), State = colnames(state_scores), `State DGE` = colnames(dge_scores)), optimal_cut_path, optimal_cut_km_path)
write.csv(optimal_cut_results, optimal_cut_results_path, row.names = FALSE, na = "")

summary_df <- cox_results %>%
  group_by(feature_type, split_method) %>%
  summarise(n_features = n_distinct(feature), n_models = n(), n_subjects_max = max(n),
            n_events_max = max(events), n_nominal_p_lt_0_05 = sum(p_value < 0.05),
            n_bh_fdr_lt_0_05 = sum(p_adj_bh < 0.05), .groups = "drop")
write.csv(summary_df, summary_path, row.names = FALSE)

report_lines <- c(
  "OCCAMS current-centred bulk survival run",
  paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "qc_rule=Annotation_Compatible_Below_50_Flag == no",
  paste0("qc_pass_subjects=", nrow(model_data)),
  paste0("survival_evaluable_subjects=", sum(!is.na(model_data$os_days) & !is.na(model_data$os_event) & model_data$os_days > 0)),
  paste0("survival_events=", sum(model_data$os_event == 1, na.rm = TRUE)),
  paste0("survival_censored=", sum(model_data$os_event == 0, na.rm = TRUE)),
  paste0("expression_genes=", n_expression_genes),
  paste0("mp_features=", sum(cox_results$feature_type == "MP" & cox_results$split_method == "continuous_per_sd")),
  paste0("state_features=", sum(cox_results$feature_type == "State" & cox_results$split_method == "continuous_per_sd")),
  paste0("km_selected=", paste(km_selection$feature, collapse = ",")),
  paste0("result=PASS")
)
writeLines(report_lines, report_path)
message("Completed OCCAMS QC-pass bulk MP/state survival workflow.")
