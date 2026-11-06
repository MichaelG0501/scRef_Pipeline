####################
# Analysis registry
# Status: active
# Script: analysis/clinical/Auto_malignant_bulk_survival.R
# Description: Matched bulk and malignant-expression survival associations in
#   OCCAMS and TCGA EAC, with within-malignant state RNA-fraction associations.
# Methodology: analysis/methodology/clinical/malignant_bulk_survival_methodology.md
# Exact inputs relative to the live project root:
#   analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx
#     Patient_Master_388: authoritative patient, treatment, specimen, AJCC stage,
#     age, sex, OS, Recurrence_Ever and Estimated_Recurrence_Time_Days_From_Diagnosis_Surrogate fields. Alive/Not recorded dead is treated as censored;
#     this source label does not establish independently confirmed vital status.
#   analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx
#     EAC_RNA88_Final: authoritative 88 EAC subjects, exact primary RNA sample,
#     AJCC8 stage and treatment-pathway/surgery classifications.
#   ref_outs/TCGA/surgery_clinical_associations/reference/clinical_cbioportal.tsv
#     Patient ID, Diagnosis Age, Sex, Overall Survival (Months), Overall Survival
#     Status supply fields absent from the September TCGA workbook. OS months
#     are multiplied by 30.4375; only explicit status codes define events.
#   ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#     Exact sample/patient mapping and historical OS concordance audit only.
#   ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   ref_outs/clinical_selected_survival/source_data/TCGA_BCR_endpoint_candidates.csv
#     Explicit historical new tumour event/death endpoint, including new primary; diagnosis origin.
#   ref_outs/clinical_selected_survival/tables/cohort_membership.csv
#     Validation only: cohort definitions come from current workbooks above.
#   ref_outs/bulk_malignant_deconvolution/intermediate/{OCCAMS,TCGA}_counts.rds
#   ref_outs/bulk_malignant_deconvolution/intermediate/{OCCAMS,TCGA}_extracted.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv
#   analysis/shared/scRef_config.R
# Outputs under ref_outs/clinical_malignant_survival/:
#   README.md: pointer to the current result index and methodology.
#   intermediate/four_method_cache.rds; tables/four_method_{cox_results,cutpoint_grid,concordance}.csv.
#   Final PDFs: one-page four-method overviews; supporting 07/08 have four forest pages each.
#   Supporting 05_continuous_os_adjustment.pdf and 06_continuous_rfs_adjustment.pdf:
#     two pages each (MPs then states), unadjusted versus stage adjusted.
#     Supporting 02/05/06 use shared significance symbols and panel n/event counts.
#   intermediate/occams_context_cache.rds: five OCCAMS selections, four OS methods.
#   tables/occams_context_{cox_results,cutpoint_grid,model_data,membership,figure_index}.csv.
#   Supporting 02: 12 pages; specimen, pathway and combined comparisons, one method/page; OS/RFS panels.
#   intermediate/rfs_analysis_cache.rds: matched exploratory RFS analyses.
#   tables/rfs_{input_manifest,patient_metadata,model_data,cox_results,concordance_results}.csv;
#     tables/endpoint_availability.csv. RFS sources are explicit, with no endpoint substitution.
#   intermediate/analysis_cache.rds: complete reproducible scores/model/plot data.
#   tables/: input_manifest, patient_metadata, cohort_membership, cohort_summary,
#     metadata_source_audit, gene_sets, gene_coverage, scoring_gene_universe, feature_definitions,
#     scores_and_fractions, model_data, cox_results, concordance_results,
#     bulk_malignant_comparison, km_source_data, km_risk_tables, stage_km_source,
#     stage_km_risk_tables, stage_cox_results, stage_model_summary and figure_index CSV files.
#   figures/final/01_malignant_mp_survival.pdf (17 MPs per page; four methods)
#   figures/final/02_malignant_state_survival.pdf (15 state features per page; four methods)
#   figures/supporting/{01_bulk_vs_malignant,02_neoadjuvant_context,
#     03_prespecified_state_km,04_pathological_stage_km}.pdf
#   figures/diagnostics/01_model_diagnostics.pdf
#   reports/{current_results.md,validation.csv}; logs/package_versions.csv
#   updates/new_updates/summaries/malignant_bulk_survival.csv (<100 KB).
# Cache: hashes, statistical function bodies, settings and package versions must
#   match. Changed inputs or model code stop;
#   SCREF_FORCE_REBUILD=TRUE explicitly rebuilds. SCREF_REPLOT_ONLY=TRUE requires
#   the matching cache. Existing current output filenames are updated in place.
#   SCREF_SURVIVAL_THREADS controls GSVA workers (default 4). The exploratory
#   optimal-cut comparison searches 13 fixed thresholds with multiplicity correction.
# Direct allocated-node run: Rscript analysis/clinical/Auto_malignant_bulk_survival.R
# PBS: /opt/pbs/bin/qsub job_scripts/clinical/Auto_malignant_bulk_survival.sh
# Environment dmtcp; 4 CPUs, 24 GB RAM, 02:00:00. PBS logs in project root.
# Plot-only PBS run: add -v SCREF_REPLOT_ONLY=TRUE and resource overrides
#   -l select=1:ncpus=1:mem=8gb -l walltime=00:30:00 to the qsub command.
# Surgery groups are exact existing 68/72 subject lists, not matched individuals.
# Treatment histories and endpoint time origins are not asserted equivalent.
####################

library(readxl)
library(GSVA)
library(survival)
library(splines)
library(ggplot2)
library(patchwork)

project_dir <- '/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline'
setwd(project_dir)
source('analysis/shared/scRef_config.R')
if (!nzchar(Sys.getenv('PBS_JOBID'))) stop('Run within a PBS allocation')
force <- identical(toupper(Sys.getenv('SCREF_FORCE_REBUILD', 'FALSE')), 'TRUE')
replot <- identical(toupper(Sys.getenv('SCREF_REPLOT_ONLY', 'FALSE')), 'TRUE')
if (force && replot) stop('Force rebuild and replot are mutually exclusive')
threads <- as.integer(Sys.getenv('SCREF_SURVIVAL_THREADS', '4'))
stopifnot(is.finite(threads), threads >= 1L)
out <- file.path(project_dir, 'ref_outs/clinical_malignant_survival')
for (d in c('intermediate', 'tables', 'reports', 'logs', 'figures/final',
            'figures/supporting', 'figures/diagnostics'))
  dir.create(file.path(out, d), recursive=TRUE, showWarnings=FALSE)
inputs <- c(
  occams_workbook='analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx',
  tcga_workbook='analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx',
  tcga_clinical='ref_outs/TCGA/surgery_clinical_associations/reference/clinical_cbioportal.tsv',
  tcga_mapping='ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds',
  occams_qc='ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv',
  prior_membership='ref_outs/clinical_selected_survival/tables/cohort_membership.csv',
  OCCAMS_counts='ref_outs/bulk_malignant_deconvolution/intermediate/OCCAMS_counts.rds',
  TCGA_counts='ref_outs/bulk_malignant_deconvolution/intermediate/TCGA_counts.rds',
  OCCAMS_extracted='ref_outs/bulk_malignant_deconvolution/intermediate/OCCAMS_extracted.rds',
  TCGA_extracted='ref_outs/bulk_malignant_deconvolution/intermediate/TCGA_extracted.rds',
  mp_genes=SCREF_MP_GENES_RDS, mp_grouping=SCREF_MP_GROUPING_CSV,
  state_markers='ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv',
  config='analysis/shared/scRef_config.R')
if (any(!file.exists(inputs))) stop('Missing required inputs: ', paste(inputs[!file.exists(inputs)], collapse=', '))
manifest <- data.frame(input=names(inputs), path=normalizePath(inputs), md5=unname(tools::md5sum(inputs)))
####################
rfs_inputs <- c(tcga_events='ref_outs/clinical_selected_survival/source_data/TCGA_BCR_endpoint_candidates.csv')
if(!all(file.exists(rfs_inputs))) stop('Missing required RFS source')
####################
packages <- c('GSVA', 'survival', 'splines', 'readxl', 'ggplot2', 'patchwork')
versions <- data.frame(package=c('R', packages), version=c(as.character(getRversion()),
  vapply(packages, function(p) as.character(packageVersion(p)), character(1))))
settings <- list(version=3L, gsva='Gaussian log2(CPM+1), same retained gene universe',
  parent_n=c(OCCAMS=281L, TCGA=88L), min_n=20L, min_events=10L,
  min_events_per_coefficient=5L, min_factor_deaths=2L,
  nonlinear_min_events=30L, nonlinear_events_per_coefficient=10L,
####################
  nonlinear_max_leverage=.5,
####################
  main_family_size=32L, logratio_family_size=5L, state_marker_n=20L,
  score_scale='one SD within endpoint-eligible cohort, fixed across adjustments',
  fraction_scale='10 percentage points within all eight malignant states',
  logratio_scale='doubling state/geometric mean of other seven states')
cache_path <- file.path(out, 'intermediate/analysis_cache.rds')
primary_cohorts <- c('occams_surgery', 'tcga_yes')
cohort_labels <- c(occams_surgery='OCCAMS surgery pathway', tcga_yes='TCGA surgery Yes',
####################
  tcga_all='TCGA surgery pathway: all patients',
####################
  occams_neo='OCCAMS neoadjuvant: all specimens',
  occams_neo_biopsy='OCCAMS neoadjuvant: biopsy', occams_neo_surgery='OCCAMS neoadjuvant: surgery')
state_names <- SCREF_PRIMARY_STATE_ORDER
all_states <- c(state_names, SCREF_TECHNICAL_STATE_ORDER, 'Malignant state unassigned')
table_out <- function(x, name) write.csv(x, file.path(out, 'tables', paste0(name, '.csv')), row.names=FALSE, na='')
positive <- function(x) is.finite(x) & x > 0
stage_value <- function(x) {
  z <- trimws(as.character(x)); z[!z %in% as.character(1:4)] <- NA_character_
  z
}
sex_value <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[!z %in% c('male', 'female')] <- NA_character_
  z
}

read_metadata <- function() {
  o <- as.data.frame(read_excel(inputs[['occams_workbook']], sheet='Patient_Master_388'))
  t <- as.data.frame(read_excel(inputs[['tcga_workbook']], sheet='EAC_RNA88_Final'))
  qc <- read.csv(inputs[['occams_qc']], stringsAsFactors=FALSE)
  cb <- read.delim(inputs[['tcga_clinical']], check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))
####################
  cb <- cb[cb[['Study ID']]=='esca_tcga_gdc' & cb[['Sample Type']]=='Primary Tumor',,drop=FALSE]
####################
  mapping <- readRDS(inputs[['tcga_mapping']])
  stopifnot(nrow(o)==388L, nrow(t)==88L, !anyDuplicated(o$SHA_ID), !anyDuplicated(t$Patient_ID),
    !anyDuplicated(cb[['Patient ID']]), !anyDuplicated(mapping$sample_barcode))
  qc_ids <- qc$Subject_ID[qc$Annotation_Compatible_Below_50_Flag=='no']
  stopifnot(length(qc_ids)==281L, all(qc_ids %in% o$SHA_ID))
  o <- o[match(qc_ids, o$SHA_ID), , drop=FALSE]
  excluded <- c('26d51bb9b7c9974b5e753fc3520ffcfa2abdbd132509a887f4a5bc5630aa05d9',
    '72089d94e1a29cc0b7df4fb8f9694343cda2c6e71ae1083fc266b5b36b62323d',
    '9ea99e16eaaf486784cb0cd7e1c846a4dcc7bbaecfccc9a04d7853093f2cf775')
  stopifnot(all(excluded %in% o$SHA_ID), all(o$Vital_Status %in% c('Dead', 'Alive/Not recorded dead')))
  surgery <- o$Treatment_Intent_Final=='curative' & o$Curative_Treatment_Modality_Final=='surgery' &
    o$RNA_Metadata_Matched=='Yes' & !o$SHA_ID %in% excluded
  neo <- o$Treatment_Intent_Final=='curative' & o$Curative_Treatment_Modality_Final=='neoadjuvant (chemo + CRT)' &
    o$RNA_Metadata_Matched=='Yes' & o$RNA_specimen_source %in% c('biopsy', 'surgery') & !o$SHA_ID %in% excluded
  surgery[is.na(surgery)] <- FALSE; neo[is.na(neo)] <- FALSE
  om <- data.frame(subject_id=o$SHA_ID, sample_id=o$SHA_ID, dataset='OCCAMS',
    OS_time=as.numeric(o$Survival_Days), OS_event=as.integer(o$Vital_Status=='Dead'),
    age=as.numeric(o$Age_At_Diagnosis), sex=sex_value(o$Gender),
    pstage=stage_value(o$Pathological_Stage_AJCC), cstage=stage_value(o$Clinical_Stage_AJCC),
    specimen=as.character(o$RNA_specimen_source), treatment=as.character(o$Curative_Treatment_Modality_Final),
    treatment_context=as.character(o$RNA_treatment_status), surgery_selected=surgery, neo_selected=neo,
    prior_exclusion=o$SHA_ID %in% excluded, vital_status_source=o$Vital_Status, stringsAsFactors=FALSE)
  j <- match(t$RNA_primary_sample_ID, mapping$sample_barcode)
  k <- match(t$Patient_ID, cb[['Patient ID']])
  stopifnot(!anyNA(j), !anyNA(k), all(mapping$case_barcode[j]==t$Patient_ID),
    all(mapping$HistologyGroup[j]=='EAC'), all(mapping$sample_type_code[j]=='01'))
  status <- cb[['Overall Survival Status']][k]
  stopifnot(all(is.na(status) | status %in% c('0:LIVING', '1:DECEASED')))
  event <- as.integer(status=='1:DECEASED')
  time <- as.numeric(cb[['Overall Survival (Months)']][k])*30.4375
  tm <- data.frame(subject_id=t$Patient_ID, sample_id=t$RNA_primary_sample_ID, dataset='TCGA',
    OS_time=time, OS_event=event, age=as.numeric(cb[['Diagnosis Age']][k]), sex=sex_value(cb[['Sex']][k]),
    pstage=stage_value(t$Pathological_Stage_AJCC8), cstage=stage_value(t$Clinical_Stage_AJCC8),
    specimen='primary tumour', treatment=as.character(t$Flowchart_Pathway_Final),
    treatment_context=as.character(t$History_neoadjuvant_for_TCGA_specimen_raw),
    surgery_selected=t$Flowchart_Pathway_Final=='Surgery group' & t$Flowchart_Surgery_Status_Final=='Yes',
    neo_selected=FALSE, prior_exclusion=FALSE, vital_status_source=status, stringsAsFactors=FALSE)
  tm$surgery_selected[is.na(tm$surgery_selected)] <- FALSE
  meta <- rbind(om, tm)
  stopifnot(!anyDuplicated(meta$subject_id), !anyDuplicated(meta$sample_id))
  meta$OS_eligible <- positive(meta$OS_time) & meta$OS_event %in% c(0L, 1L)
  meta$OS_exclusion <- ifelse(!positive(meta$OS_time), 'missing or nonpositive OS interval',
    ifelse(!meta$OS_event %in% c(0L, 1L), 'missing explicit vital status', ''))
  cohorts <- list(occams_surgery=om$subject_id[om$surgery_selected], tcga_yes=tm$subject_id[tm$surgery_selected],
    occams_neo=om$subject_id[om$neo_selected],
    occams_neo_biopsy=om$subject_id[om$neo_selected & om$specimen=='biopsy'],
    occams_neo_surgery=om$subject_id[om$neo_selected & om$specimen=='surgery'])
  stopifnot(identical(as.integer(lengths(cohorts)), c(68L, 72L, 157L, 77L, 80L)))
####################
  cohorts$tcga_all <- tm$subject_id[tm$treatment=='Surgery group']
  stopifnot(length(cohorts$tcga_all)==86L)
####################
  previous <- read.csv(inputs[['prior_membership']], stringsAsFactors=FALSE)
  for (key in names(cohorts)) stopifnot(setequal(cohorts[[key]], previous$subject_id[previous$cohort==key]))
  audit <- data.frame(subject_id=tm$subject_id, sample_id=tm$sample_id, workbook_stage=tm$pstage,
    cbio_OS_days=time, mapped_OS_days=as.numeric(mapping$OS_time[j]),
    OS_days_difference=time-as.numeric(mapping$OS_time[j]), cbio_OS_event=event,
    mapped_OS_event=as.numeric(mapping$OS_event[j]), stringsAsFactors=FALSE)
  list(meta=meta, cohorts=cohorts, audit=audit)
}

####################
read_rfs_metadata <- function(metadata) {
  o <- as.data.frame(read_excel(inputs[['occams_workbook']],sheet='Patient_Master_388'))
  e <- read.csv(rfs_inputs[['tcga_events']],stringsAsFactors=FALSE,na.strings='')
  stopifnot(!anyDuplicated(o$SHA_ID),!anyDuplicated(e$Patient_ID))
  meta <- metadata
  meta$RFS_time <- meta$RFS_event <- NA_real_
  meta$RFS_exclusion <- meta$RFS_definition <- ''
  blank <- function(x) is.na(x) | trimws(as.character(x))==''
  for (i in seq_len(nrow(meta))) {
    if(meta$dataset[i]=='OCCAMS') {
      j <- match(meta$subject_id[i],o$SHA_ID)
      stopifnot(!is.na(j))
      os_time <- as.numeric(o$Survival_Days[j])
      os_event <- as.integer(o$Vital_Status[j]=='Dead')
      recurrence <- as.numeric(o$Estimated_Recurrence_Time_Days_From_Diagnosis_Surrogate[j])
      rec_yes <- !is.na(o$Recurrence_Ever[j]) && o$Recurrence_Ever[j]=='Yes'
      rec_no <- !is.na(o$Recurrence_Ever[j]) && o$Recurrence_Ever[j]=='No'
      meta$RFS_definition[i] <- 'Exploratory recurrence/death; assumed common first-therapy origin'
      if(!positive(os_time)) meta$RFS_exclusion[i] <- 'missing/nonpositive survival interval'
      else if(rec_yes && !positive(recurrence)) meta$RFS_exclusion[i] <- 'recorded recurrence without positive timed interval'
      else if(positive(recurrence) && recurrence>os_time) meta$RFS_exclusion[i] <- 'recurrence interval exceeds survival/follow-up interval'
      else if(!rec_yes && !rec_no && os_event==0L) meta$RFS_exclusion[i] <- 'living without recurrence ascertainment'
      else if(rec_yes) {
        meta$RFS_time[i] <- if(os_event==1L) min(recurrence,os_time) else recurrence
        meta$RFS_event[i] <- 1L
      } else {
        meta$RFS_time[i] <- os_time
        meta$RFS_event[i] <- os_event
      }
    } else {
      j <- match(meta$subject_id[i],e$Patient_ID)
      stopifnot(!is.na(j))
      meta$RFS_definition[i] <- 'Original new tumour event or death; includes new primary; diagnosis origin'
      if(!blank(e$exclusion_reason[j])) meta$RFS_exclusion[i] <- e$exclusion_reason[j]
      else if(blank(e$tumour_event_indicator[j]) && blank(e$tumour_event_days[j]) && blank(e$death_days[j]))
        meta$RFS_exclusion[i] <- 'living without tumour-event ascertainment'
      if(meta$RFS_exclusion[i]=='') {
        meta$RFS_time[i] <- as.numeric(e$candidate_rfs_days[j])
        meta$RFS_event[i] <- as.numeric(e$candidate_rfs_event[j])
      }
    }
  }
  meta$RFS_eligible <- positive(meta$RFS_time) & meta$RFS_event %in% c(0L,1L)
  invalid <- !meta$RFS_eligible & meta$RFS_exclusion==''
  meta$RFS_exclusion[invalid] <- 'missing/nonpositive RFS interval or event status'
  meta
}
####################

make_features <- function() {
  grouping <- read.csv(inputs[['mp_grouping']], stringsAsFactors=FALSE)
  genes <- readRDS(inputs[['mp_genes']])[grouping$mp]
  stopifnot(length(genes)==17L, !any(vapply(genes, is.null, logical(1))),
    !any(c('MP2x', 'MP11c', 'MP18a') %in% names(genes)), identical(names(SCREF_STATE_GROUPS), state_names))
  unions <- lapply(SCREF_STATE_GROUPS, function(mps) unique(unlist(genes[mps], use.names=FALSE)))
  names(unions) <- paste0(state_names, ' (MP union)')
  markers <- read.csv(inputs[['state_markers']], stringsAsFactors=FALSE)
  dge <- lapply(state_names, function(s) head(markers$gene[markers$state==s], settings$state_marker_n))
  names(dge) <- paste0(state_names, ' (DGE)')
  stopifnot(all(lengths(dge)==20L), all(vapply(dge, function(x) !anyDuplicated(x), logical(1))))
  sets <- c(genes, unions, dge)
  definitions <- data.frame(feature=c(names(sets), paste0(state_names, ' (fraction)'), paste0(state_names, ' (logratio)')),
    family=c(rep('MP', 17L), rep('State MP union', 5L), rep('State DGE', 5L), rep('State fraction', 5L), rep('State logratio', 5L)),
    state=c(grouping$state, rep(state_names, 4L)),
    label=c(paste(grouping$mp, grouping$description, sep=' | '), rep(state_names, 4L)),
    stringsAsFactors=FALSE)
  definitions$expected_direction <- ifelse(definitions$state %in% c('Cell cycle', 'Classic proliferation', 'Stress-adaptive'),
    'higher hazard hypothesised', ifelse(definitions$state %in% c('Squamous-to-intestinal', 'Glandular-to-intestinal'),
      'lower hazard hypothesised', 'no prespecified direction'))
  list(sets=sets, definitions=definitions)
}

score_expression <- function(metadata, features) {
  counts <- setNames(lapply(c('OCCAMS', 'TCGA'), function(d) readRDS(inputs[[paste0(d, '_counts')]])), c('OCCAMS', 'TCGA'))
  extracted <- setNames(lapply(c('OCCAMS', 'TCGA'), function(d) readRDS(inputs[[paste0(d, '_extracted')]])), c('OCCAMS', 'TCGA'))
  common <- Reduce(intersect, c(lapply(counts, rownames), lapply(extracted, function(x) colnames(x$tumour_expression))))
  stopifnot(length(common)>1000L, !anyDuplicated(common))
  # Remove a gene only if constant in any complete parent compartment. The same
  # universe is then used by all four GSVA calls, independently of outcomes.
  expression <- list()
  for (dataset in names(counts)) {
    ids <- metadata$sample_id[metadata$dataset==dataset]
    stopifnot(length(ids)==settings$parent_n[[dataset]], setequal(colnames(counts[[dataset]]), ids),
      setequal(rownames(extracted[[dataset]]$tumour_expression), ids),
      setequal(rownames(extracted[[dataset]]$within), ids))
    for (source in c('bulk', 'malignant')) {
      x <- if (source=='bulk') counts[[dataset]][common, ids, drop=FALSE] else
        t(extracted[[dataset]]$tumour_expression[ids, common, drop=FALSE])
      stopifnot(all(is.finite(x)), all(x>=0), all(colSums(x)>0))
      expression[[paste(dataset, source, sep='__')]] <- log2(sweep(x, 2L, colSums(x), '/')*1e6+1)
    }
  }
  variable <- Reduce(`&`, lapply(expression, function(x) apply(x, 1L, function(y) diff(range(y))>0)))
  common <- common[variable]
  sets <- lapply(features$sets, function(g) intersect(unique(g), common))
  stopifnot(length(common)>1000L, all(lengths(sets)>=5L))
  scores <- list()
  for (key in names(expression)) {
    message('Gaussian GSVA: ', key, ', n=', ncol(expression[[key]]))
    x <- GSVA::gsva(expression[[key]][common, , drop=FALSE], sets, method='gsva', kcdf='Gaussian',
      min.sz=5L, max.sz=Inf, mx.diff=TRUE, abs.ranking=FALSE, verbose=FALSE, parallel.sz=threads)
    stopifnot(identical(rownames(x), names(sets)), all(is.finite(x)))
    scores[[key]] <- t(x)
  }
  fractions <- list()
  for (dataset in names(extracted)) {
    e <- extracted[[dataset]]; ids <- metadata$sample_id[metadata$dataset==dataset]
    w <- e$within[ids, all_states, drop=FALSE]
    stopifnot(all(is.finite(w)), all(w>0), max(abs(rowSums(w)-1))<1e-6)
    lr <- sapply(state_names, function(s) log2(w[, s])-rowMeans(log2(w[, setdiff(all_states, s), drop=FALSE])))
    colnames(lr) <- paste0(state_names, ' (logratio)')
    p <- w[, state_names, drop=FALSE]; colnames(p) <- paste0(state_names, ' (fraction)')
    fractions[[dataset]] <- cbind(p, lr)
  }
  coverage <- data.frame(feature=names(sets), source_genes=lengths(features$sets), retained_genes=lengths(sets),
    retained_fraction=lengths(sets)/lengths(features$sets), universe_genes=length(common))
  list(scores=scores, fractions=fractions, coverage=coverage, common_genes=common,
    raw_fractions=lapply(extracted, `[[`, 'fractions'))
}

build_model_data <- function(metadata, cohorts, features, scored) {
  pieces <- list(); i <- 0L
  for (cohort in names(cohorts)) {
    meta <- metadata[match(cohorts[[cohort]], metadata$subject_id), , drop=FALSE]
    dataset <- unique(meta$dataset); stopifnot(length(dataset)==1L)
####################
    f <- scored$raw_fractions[[dataset]]
    meta$malignant_fraction <- f$malignant_final[match(meta$sample_id,f$sample_id)]
    stopifnot(all(is.finite(meta$malignant_fraction)),all(meta$malignant_fraction>=0 & meta$malignant_fraction<=1))
####################
    for (source in c('bulk', 'malignant')) {
      x <- scored$scores[[paste(dataset, source, sep='__')]][meta$sample_id, , drop=FALSE]
      if (source=='malignant') x <- cbind(x, scored$fractions[[dataset]][meta$sample_id, , drop=FALSE])
      for (feature in colnames(x)) {
        value <- x[, feature]
        family <- features$definitions$family[match(feature, features$definitions$feature)]
        scale <- if (family=='State fraction') .1 else if (family=='State logratio') 1 else sd(value[meta$OS_eligible])
        if (!is.finite(scale) || scale<=0) stop('Invalid prespecified predictor scale: ', cohort, '/', source, '/', feature)
        i <- i+1L
        pieces[[i]] <- data.frame(cohort=cohort, source=source, feature=feature, family=family,
          meta, value=value, scale=scale, predictor=value/scale, stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind, pieces)
}

fit_one <- function(d, adjustment) {
  result <- data.frame(cohort=d$cohort[1], dataset=d$dataset[1], source=d$source[1],
    feature=d$feature[1], family=d$family[1], adjustment=adjustment,
    n_selected=nrow(d), n_endpoint=sum(d$OS_eligible), n=NA_integer_, events=NA_integer_,
    n_excluded_endpoint=sum(!d$OS_eligible), n_excluded_covariates=NA_integer_, parameters=NA_integer_,
    scale=d$scale[1], status='estimated', warning='', logHR=NA_real_, SE=NA_real_, HR=NA_real_,
    lower_CI=NA_real_, upper_CI=NA_real_, P_value=NA_real_, PH_feature_P=NA_real_, PH_global_P=NA_real_,
    linearity_P=NA_real_, linearity_max_leverage=NA_real_, linearity_status='not eligible', stringsAsFactors=FALSE)
  covariates <- switch(adjustment, unadjusted=character(0), stage_complete_unadjusted=character(0),
    stage='pstage', stage_purity=c('pstage','malignant_fraction'), age_sex_stage=c('age', 'sex', 'pstage'))
  d <- d[d$OS_eligible, , drop=FALSE]
####################
  if (adjustment=='stage_complete_unadjusted') d <- d[!is.na(d$pstage),,drop=FALSE]
####################
  if (length(covariates)) d <- d[complete.cases(d[, covariates, drop=FALSE]), , drop=FALSE]
  result$n <- nrow(d); result$events <- sum(d$OS_event)
  result$n_excluded_covariates <- result$n_endpoint-result$n
  fail <- function(reason) { result$status <- reason; result }
  if (nrow(d)<settings$min_n) return(fail('not estimated: fewer than 20 complete subjects'))
  if (result$events<settings$min_events) return(fail('not estimated: fewer than 10 deaths'))
  if (sd(d$predictor)==0) return(fail('not estimated: constant predictor'))
  d$pstage <- factor(d$pstage, levels=sort(unique(d$pstage[!is.na(d$pstage)])))
  d$sex <- factor(d$sex, levels=sort(unique(d$sex[!is.na(d$sex)])))
  for (variable in intersect(covariates, c('sex', 'pstage'))) {
    if (nlevels(d[[variable]])<2L) return(fail(paste('not estimated: single', variable, 'category')))
    if (any(tapply(d$OS_event, d[[variable]], sum)<settings$min_factor_deaths))
      return(fail(paste('not estimated: fewer than 2 deaths in a', variable, 'category')))
  }
  formula <- reformulate(c('predictor', covariates), response='Surv(OS_time, OS_event)')
  design <- model.matrix(formula, d)[, -1L, drop=FALSE]
  result$parameters <- ncol(design)
  if (qr(design)$rank<ncol(design)) return(fail('not estimated: rank-deficient design'))
  if (result$events<settings$min_events_per_coefficient*ncol(design))
    return(fail('not estimated: fewer than 5 deaths per coefficient'))
  warnings <- character(0)
  fit <- withCallingHandlers(coxph(formula, data=d, ties='efron', x=TRUE,
    control=coxph.control(iter.max=50L)), warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') })
  if (length(warnings)) {
    result$warning <- paste(unique(warnings), collapse='; ')
    stop('Cox warning for ',result$cohort,' / ',result$source,' / ',result$feature,' / ',adjustment,': ',result$warning)
  }
  s <- summary(fit)
  stopifnot('predictor' %in% rownames(s$coefficients))
  values <- c(coef(fit)['predictor'], s$coefficients['predictor', 'se(coef)'], s$coefficients['predictor', 'Pr(>|z|)'])
  if (any(!is.finite(values))) stop('Nonfinite Cox coefficient, SE or P: ',result$feature)
  result$logHR <- unname(coef(fit)['predictor']); result$SE <- s$coefficients['predictor', 'se(coef)']
  result$HR <- s$conf.int['predictor', 'exp(coef)']; result$lower_CI <- s$conf.int['predictor', 'lower .95']
  result$upper_CI <- s$conf.int['predictor', 'upper .95']; result$P_value <- s$coefficients['predictor', 'Pr(>|z|)']
  ph <- cox.zph(fit, transform='km', terms=TRUE)$table
  result$PH_feature_P <- ph['predictor', 'p']; result$PH_global_P <- ph['GLOBAL', 'p']
  spline_parameters <- ncol(design)+2L
  if (result$events>=settings$nonlinear_min_events &&
      result$events>=settings$nonlinear_events_per_coefficient*spline_parameters && length(unique(d$predictor))>=10L) {
    spline_formula <- reformulate(c('splines::ns(predictor, df=3)', covariates), response='Surv(OS_time, OS_event)')
####################
    # A spline must have support across patients, rather than a basis direction
    # determined mostly by one observation. This check uses the design only.
    spline_design <- model.matrix(spline_formula,d)
    spline_qr <- qr(spline_design)
    if (spline_qr$rank<ncol(spline_design)) {
      result$linearity_status <- 'not eligible: rank-deficient spline design'
      return(result)
    }
    result$linearity_max_leverage <- max(rowSums(qr.Q(spline_qr)^2))
    if (result$linearity_max_leverage>=settings$nonlinear_max_leverage) {
      result$linearity_status <- 'not eligible: spline design leverage at least 0.5'
      return(result)
    }
####################
    spline_warnings <- character(0)
    flexible <- withCallingHandlers(coxph(spline_formula, data=d, ties='efron', control=coxph.control(iter.max=50L)),
      warning=function(w) { spline_warnings <<- c(spline_warnings, conditionMessage(w)); invokeRestart('muffleWarning') })
    if (length(spline_warnings)) stop('Spline model warning for ',result$cohort,' / ',result$source,' / ',
      result$feature,' / ',adjustment,': ',paste(spline_warnings,collapse='; ')) else {
      result$linearity_P <- pchisq(max(0, 2*(as.numeric(logLik(flexible))-as.numeric(logLik(fit)))), df=2L, lower.tail=FALSE)
      result$linearity_status <- 'estimated: natural cubic spline df3 versus linear'
    }
  }
  result
}

run_models <- function(model_data) {
  keys <- interaction(model_data$cohort, model_data$source, model_data$feature, drop=TRUE, lex.order=TRUE)
  results <- do.call(rbind, lapply(split(model_data, keys), function(d)
    do.call(rbind, lapply(c('unadjusted', 'stage_complete_unadjusted', 'stage', 'stage_purity', 'age_sex_stage'), function(a) fit_one(d, a)))))
  results$P_BH <- NA_real_; results$PH_feature_BH <- NA_real_
  family_key <- interaction(results$cohort, results$source, results$adjustment,
    ifelse(results$family=='State logratio', 'logratio', 'main'), drop=TRUE)
  for (idx in split(seq_len(nrow(results)), family_key)) {
    n <- if (all(results$family[idx]=='State logratio')) 5L else 32L
    results$P_BH[idx] <- p.adjust(results$P_value[idx], method='BH', n=n)
    results$PH_feature_BH[idx] <- p.adjust(results$PH_feature_P[idx], method='BH', n=n)
  }
  results
}

cross_cohort <- function(results, definitions) {
  a <- results[results$cohort=='occams_surgery', ]; b <- results[results$cohort=='tcga_yes', ]
  columns <- c('source', 'feature', 'family', 'adjustment')
  x <- merge(a[, c(columns, 'n', 'events', 'status', 'logHR', 'SE', 'HR', 'P_value', 'P_BH')],
    b[, c(columns, 'n', 'events', 'status', 'logHR', 'SE', 'HR', 'P_value', 'P_BH')], by=columns, suffixes=c('_OCCAMS', '_TCGA'), all=TRUE)
  x$concordant <- with(x, is.finite(logHR_OCCAMS) & is.finite(logHR_TCGA) & logHR_OCCAMS*logHR_TCGA>0)
  available <- with(x, is.finite(P_value_OCCAMS) & is.finite(P_value_TCGA))
  x$conjunction_P <- NA_real_
  x$conjunction_P[available] <- ifelse(x$concordant[available],
    pmax(x$P_value_OCCAMS[available], x$P_value_TCGA[available]), 1)
  x$conjunction_BH <- NA_real_
  groups <- interaction(x$source, x$adjustment, ifelse(x$family=='State logratio', 'logratio', 'main'), drop=TRUE)
  for (idx in split(seq_len(nrow(x)), groups))
    x$conjunction_BH[idx] <- p.adjust(x$conjunction_P[idx], 'BH', n=if (all(x$family[idx]=='State logratio')) 5L else 32L)
  x$nominal_replication <- x$concordant & !is.na(x$conjunction_P) & x$conjunction_P<.05
  x$FDR_replication <- x$concordant & !is.na(x$conjunction_BH) & x$conjunction_BH<.05
  x$direction <- ifelse(!available, 'not estimable', ifelse(!x$concordant, 'discordant',
    ifelse(x$logHR_OCCAMS>0, 'higher hazard in both', 'lower hazard in both')))
####################
  x$expected_direction <- definitions$expected_direction[match(x$feature,definitions$feature)]
  stopifnot(!anyNA(x$expected_direction))
  x$matches_hypothesis <- ifelse(x$expected_direction=='no prespecified direction',NA,
    x$concordant & ((x$expected_direction=='higher hazard hypothesised' & x$logHR_OCCAMS>0) |
      (x$expected_direction=='lower hazard hypothesised' & x$logHR_OCCAMS<0)))
####################
  x
}

make_km <- function(model_data) {
  curves <- risks <- list(); i <- 0L
  keep <- model_data$cohort %in% primary_cohorts & model_data$source=='malignant' &
    model_data$family %in% c('State MP union', 'State fraction') & model_data$OS_eligible
  d <- model_data[keep, ]
  for (z in split(d, interaction(d$cohort, d$feature, drop=TRUE))) {
    cut <- median(z$value); z$group <- factor(ifelse(z$value>cut, 'Above median', 'At or below median'),
      levels=c('At or below median', 'Above median'))
    if (any(table(z$group)<5L)) stop('Median KM group below five subjects')
    fit <- survfit(Surv(OS_time/365.25, OS_event) ~ group, data=z, conf.type='log-log')
    strata <- rep(names(fit$strata), fit$strata)
    i <- i+1L
    curves[[i]] <- data.frame(cohort=z$cohort[1], feature=z$feature[1], family=z$family[1], cutpoint=cut,
      time=fit$time, survival=fit$surv, lower=fit$lower, upper=fit$upper, censored=fit$n.censor,
      group=sub('group=', '', strata), stringsAsFactors=FALSE)
    times <- seq(0, floor(max(z$OS_time)/365.25), by=2)
    s <- summary(fit, times=times, extend=TRUE)
    risks[[i]] <- data.frame(cohort=z$cohort[1], feature=z$feature[1], family=z$family[1], cutpoint=cut,
      time=s$time, n_risk=s$n.risk, group=sub('group=', '', s$strata), stringsAsFactors=FALSE)
  }
  list(curves=do.call(rbind, curves), risks=do.call(rbind, risks))
}

make_stage_km <- function(metadata, cohorts) {
  curves <- risks <- list()
  for (cohort in primary_cohorts) {
    d <- metadata[metadata$subject_id %in% cohorts[[cohort]] & metadata$OS_eligible & !is.na(metadata$pstage), ]
    d$stage <- factor(d$pstage, levels=as.character(1:4))
    fit <- survfit(Surv(OS_time/365.25, OS_event) ~ stage, data=d, conf.type='log-log')
    curves[[cohort]] <- data.frame(cohort=cohort, time=fit$time, survival=fit$surv,
      lower=fit$lower, upper=fit$upper, censored=fit$n.censor, stage=sub('stage=', '', rep(names(fit$strata), fit$strata)))
    s <- summary(fit, times=seq(0, floor(max(d$OS_time)/365.25), by=2), extend=TRUE)
    risks[[cohort]] <- data.frame(cohort=cohort, time=s$time, n_risk=s$n.risk, stage=sub('stage=', '', s$strata))
  }
  list(curves=do.call(rbind, curves), risks=do.call(rbind, risks))
}

####################
fit_stage_models <- function(metadata,cohorts) {
  coefficients <- summaries <- list()
  for (cohort in primary_cohorts) {
    d <- metadata[metadata$subject_id %in% cohorts[[cohort]] & metadata$OS_eligible & !is.na(metadata$pstage), ]
    d$pstage <- factor(d$pstage,levels=as.character(1:4))
    stopifnot(nrow(d)>=settings$min_n,sum(d$OS_event)>=3L*settings$min_events_per_coefficient,
      all(tapply(d$OS_event,d$pstage,sum)>=settings$min_factor_deaths))
    fit <- withCallingHandlers(coxph(Surv(OS_time,OS_event)~pstage,data=d,ties='efron',x=TRUE),
      warning=function(w) stop('Stage-only Cox model: ',conditionMessage(w)))
    s <- summary(fit)
    stopifnot(all(is.finite(s$coefficients)),all(is.finite(s$conf.int)))
    coefficients[[cohort]] <- data.frame(cohort=cohort,stage=as.character(2:4),reference_stage='1',
      n=nrow(d),events=sum(d$OS_event),HR=s$conf.int[,'exp(coef)'],
      lower_CI=s$conf.int[,'lower .95'],upper_CI=s$conf.int[,'upper .95'],
      P_value=s$coefficients[,'Pr(>|z|)'],row.names=NULL)
    summaries[[cohort]] <- data.frame(cohort=cohort,n=nrow(d),events=sum(d$OS_event),
      global_LR_P=unname(s$logtest['pvalue']),
      global_PH_P=cox.zph(fit,transform='km')$table['GLOBAL','p'])
  }
  list(coefficients=do.call(rbind,coefficients),summary=do.call(rbind,summaries))
}
####################
####################
analysis_code <- lapply(list(positive,stage_value,sex_value,read_metadata,make_features,score_expression,build_model_data,
  fit_one,run_models,cross_cohort,make_km,make_stage_km,fit_stage_models),function(f) paste(deparse(body(f)),collapse='\n'))
####################
if (file.exists(cache_path) && !force) {
  cache <- readRDS(cache_path)
  if (!identical(cache$manifest, manifest) || !identical(cache$settings, settings) || !identical(cache$versions, versions) ||
      !identical(cache$analysis_code,analysis_code))
    stop('Cache provenance changed; explicitly rebuild with SCREF_FORCE_REBUILD=TRUE')
  message('Using verified current analysis cache')
} else {
  if (replot) stop('Replot requires the current matching analysis cache')
  clinical <- read_metadata(); features <- make_features()
  scored <- score_expression(clinical$meta, features)
  model_data <- build_model_data(clinical$meta, clinical$cohorts, features, scored)
  results <- run_models(model_data)
  cache <- list(manifest=manifest, settings=settings, versions=versions, analysis_code=analysis_code, clinical=clinical,
    features=features, scored=scored, model_data=model_data, results=results,
    concordance=cross_cohort(results,features$definitions), km=make_km(model_data), stage_km=make_stage_km(clinical$meta, clinical$cohorts),
    stage_models=fit_stage_models(clinical$meta,clinical$cohorts))
  saveRDS(cache, cache_path)
}

table_out(cache$manifest, 'input_manifest'); table_out(cache$clinical$meta, 'patient_metadata')
table_out(cache$clinical$audit, 'metadata_source_audit'); table_out(cache$features$definitions, 'feature_definitions')
table_out(cache$scored$coverage, 'gene_coverage')
####################
table_out(data.frame(gene=cache$scored$common_genes),'scoring_gene_universe')
####################
gene_sets <- do.call(rbind, lapply(names(cache$features$sets), function(f) data.frame(feature=f,
  gene=cache$features$sets[[f]], retained=cache$features$sets[[f]] %in% cache$scored$common_genes)))
table_out(gene_sets, 'gene_sets')
membership <- do.call(rbind, lapply(names(cache$clinical$cohorts), function(c)
  data.frame(cohort=c, subject_id=cache$clinical$cohorts[[c]])))
table_out(membership, 'cohort_membership')
cohort_summary <- do.call(rbind, lapply(names(cache$clinical$cohorts), function(c) {
  d <- cache$clinical$meta[match(cache$clinical$cohorts[[c]], cache$clinical$meta$subject_id), ]
  data.frame(cohort=c, n_selected=nrow(d), n_OS=sum(d$OS_eligible), OS_events=sum(d$OS_event[d$OS_eligible]),
    n_pathological_stage=sum(d$OS_eligible & !is.na(d$pstage)),
    n_age_sex_stage=sum(d$OS_eligible & complete.cases(d[, c('age', 'sex', 'pstage')])))}))
table_out(cohort_summary, 'cohort_summary')
all_scores <- do.call(rbind, lapply(names(cache$scored$scores), function(key) {
  x <- cache$scored$scores[[key]]
  do.call(rbind, lapply(colnames(x), function(f) data.frame(parent=key, sample_id=rownames(x), feature=f, value=x[, f])))}))
all_fractions <- do.call(rbind, lapply(names(cache$scored$fractions), function(key) {
  x <- cache$scored$fractions[[key]]
  do.call(rbind, lapply(colnames(x), function(f) data.frame(parent=paste0(key, '__malignant'), sample_id=rownames(x), feature=f, value=x[, f])))}))
table_out(rbind(all_scores, all_fractions), 'scores_and_fractions')
table_out(cache$model_data, 'model_data'); table_out(cache$results, 'cox_results')
table_out(cache$concordance, 'concordance_results'); table_out(cache$km$curves, 'km_source_data')
table_out(cache$km$risks, 'km_risk_tables'); table_out(cache$stage_km$curves, 'stage_km_source')
table_out(cache$stage_km$risks, 'stage_km_risk_tables')
####################
table_out(cache$stage_models$coefficients,'stage_cox_results')
table_out(cache$stage_models$summary,'stage_model_summary')
####################
bulk <- cache$results[cache$results$source=='bulk', ]; malignant <- cache$results[cache$results$source=='malignant', ]
comparison <- merge(bulk, malignant, by=c('cohort', 'dataset', 'feature', 'family', 'adjustment'), suffixes=c('_bulk', '_malignant'))
comparison$direction_changed <- with(comparison, is.finite(logHR_bulk) & is.finite(logHR_malignant) & logHR_bulk*logHR_malignant<0)
comparison$delta_logHR <- comparison$logHR_malignant-comparison$logHR_bulk
table_out(comparison, 'bulk_malignant_comparison')
write.csv(versions, file.path(out, 'logs/package_versions.csv'), row.names=FALSE)

theme_set(theme_classic(base_size=12, base_family='sans')+
  theme(plot.title=element_text(face='bold', size=14), plot.subtitle=element_text(size=10.5),
    strip.background=element_blank(), strip.text=element_text(face='bold', size=11),
    axis.text=element_text(colour='black', size=10), legend.position='bottom',
    legend.title=element_blank(), plot.title.position='plot',plot.caption.position='plot',
    plot.margin=margin(9, 12, 9, 9)))
cohort_colours <- c(OCCAMS='#0072B2', TCGA='#D55E00')
family_labels <- c('MP'='Metaprogrammes', 'State MP union'='MP-union\nstate scores',
  'State DGE'='DGE\nstate scores', 'State fraction'='Within-malignant\nstate fractions')
format_feature <- function(feature) {
  i <- match(feature, cache$features$definitions$feature)
  label <- cache$features$definitions$label[i]
  vapply(label, function(x) paste(strwrap(x, width=47), collapse='\n'), character(1))
}

####################
interval_coordinates <- function(d) {
  # A common log2 axis makes cohort comparisons legible. Arrowheads identify
  # intervals extending past the display; complete estimates remain in CSVs.
  axis_limit <- 3
  d$effect <- pmax(-axis_limit,pmin(axis_limit,d$logHR/log(2)))
  d$left <- pmax(-axis_limit,(d$logHR-qnorm(.975)*d$SE)/log(2))
  d$right <- pmin(axis_limit,(d$logHR+qnorm(.975)*d$SE)/log(2))
  d$left <- pmin(d$left,axis_limit)
  d$right <- pmax(d$right,-axis_limit)
  d$left_overflow <- (d$logHR-qnorm(.975)*d$SE)/log(2) < -axis_limit
  d$right_overflow <- (d$logHR+qnorm(.975)*d$SE)/log(2)>axis_limit
  d
}
forest_axis <- function() scale_x_continuous(limits=c(-3.15,3.15),breaks=-3:3,
  labels=c('0.125','0.25','0.5','1','2','4','8'))
interval_overflow <- function(d) list(
  geom_text(data=d,aes(x=-3,label=ifelse(left_overflow,'\u2039','')),
    position=position_dodge(.6),size=5,show.legend=FALSE),
  geom_text(data=d,aes(x=3,label=ifelse(right_overflow,'\u203a','')),
    position=position_dodge(.6),size=5,show.legend=FALSE))
####################

####################
# Counts are taken from the displayed models. Ranges reflect feature-specific
# quartile sample sizes, rather than implying a single common denominator.
count_range <- function(x) {
  stopifnot(length(x)>0L,all(is.finite(x)))
  limits <- range(x)
  if(limits[1]==limits[2]) as.character(limits[1]) else paste(limits,collapse='–')
}
panel_count_label <- function(d,variable,levels,labels) {
  paste(vapply(seq_along(levels),function(i) {
    z <- d[as.character(d[[variable]])==levels[i],]
    paste0(labels[i],': n=',count_range(z$n),', events=',count_range(z$events))
  },character(1)),collapse='\n')
}
####################
forest_plot <- function(d, feature_order, title, subtitle, facet='adjustment') {
  d$feature_plot <- factor(d$feature, levels=rev(feature_order))
####################
  d$family <- factor(d$family,levels=names(family_labels))
####################
  d$adjustment <- factor(d$adjustment, levels=c('unadjusted', 'stage', 'age_sex_stage'),
    labels=c('Unadjusted', 'Pathological stage adjusted', 'Age, sex and stage adjusted'))
####################
  panel_names <- c('Unadjusted'='Unadjusted','Pathological stage adjusted'='Stage adjusted',
    'Age, sex and stage adjusted'='Age, sex and stage adjusted')
  strip_labels <- setNames(vapply(levels(droplevels(d$adjustment)),function(a)
    paste(panel_names[[a]],panel_count_label(d[d$adjustment==a,],'dataset',
      c('OCCAMS','TCGA'),c('OCCAMS','TCGA')),sep='\n'),character(1)),levels(droplevels(d$adjustment)))
  d$significance <- factor(ifelse(d$P_BH<.05,'FDR',ifelse(d$P_value<.05,'Nominal','NS')),
    levels=c('NS','Nominal','FDR'))
####################
  estimated <- d[d$status=='estimated', ]
  if (!nrow(estimated)) stop('No estimable models for requested forest panel')
####################
  estimated <- interval_coordinates(estimated)
####################
  dodge <- position_dodge(width=.6)
  p <- ggplot(d, aes(y=feature_plot, colour=dataset, group=dataset))+
####################
    geom_blank(aes(x=0))+
####################
    geom_vline(xintercept=0, linewidth=.4, colour='grey55', linetype=2)+
    geom_errorbar(data=estimated, aes(x=effect,xmin=left,xmax=right), orientation='y', width=.16,
      linewidth=.5, position=dodge)+
    geom_point(data=estimated,aes(x=effect,shape=significance),size=2.3,stroke=.6,position=dodge,show.legend=TRUE)+
    forest_axis()+
    scale_y_discrete(labels=format_feature)+
    scale_colour_manual(values=cohort_colours,name=NULL)+
    scale_shape_manual(values=c(NS=1,Nominal=16,FDR=18),drop=FALSE,name='Evidence')+
    labs(x='Hazard ratio (95% CI)', y=NULL, title=title)+
    theme(panel.grid.major.x=element_line(colour='grey92', linewidth=.25), panel.spacing.x=grid::unit(1, 'lines'),
      strip.text.x=element_text(size=10.5,face='plain',lineheight=1.1))
####################
  p <- p+geom_text(data=estimated,aes(x=-3,label=ifelse(left_overflow,'\u2039','')),
    position=dodge,size=5,show.legend=FALSE)+
    geom_text(data=estimated,aes(x=3,label=ifelse(right_overflow,'\u203a','')),
      position=dodge,size=5,show.legend=FALSE)
####################
  if (facet=='adjustment') p <- p+facet_grid(.~adjustment,labeller=labeller(adjustment=strip_labels))
  if (facet=='family_adjustment') p <- p+facet_grid(family~adjustment, scales='free_y', space='free_y',
    labeller=labeller(family=family_labels,adjustment=strip_labels))+theme(strip.text.y=element_text(angle=0,size=10))
  p
}

main <- cache$results[cache$results$cohort %in% primary_cohorts & cache$results$source=='malignant' &
  cache$results$adjustment %in% c('unadjusted', 'stage') & cache$results$family!='State logratio', ]
mp_order <- cache$features$definitions$feature[cache$features$definitions$family=='MP']
state_order <- cache$features$definitions$feature[cache$features$definitions$family %in% c('State MP union', 'State DGE', 'State fraction')]

paired <- cache$results[cache$results$cohort %in% primary_cohorts & cache$results$adjustment=='unadjusted' &
  cache$results$family %in% c('MP', 'State MP union', 'State DGE'), ]
paired$feature_plot <- factor(paired$feature, levels=rev(c(mp_order, state_order[!grepl('fraction', state_order)])))
paired$source <- factor(paired$source, levels=c('bulk', 'malignant'), labels=c('Whole bulk', 'Malignant expression'))
####################
paired$family <- factor(paired$family,levels=names(family_labels))
####################
paired_estimated <- interval_coordinates(paired[paired$status=='estimated', ])
p <- ggplot(paired, aes(y=feature_plot, colour=source, shape=source))+
####################
  geom_blank(aes(x=0))+
####################
  geom_vline(xintercept=0, linetype=2, colour='grey60', linewidth=.4)+
  geom_errorbar(data=paired_estimated,aes(x=effect,xmin=left,xmax=right),orientation='y',width=.12,position=position_dodge(.6),linewidth=.4)+
  geom_point(data=paired_estimated,aes(x=effect),position=position_dodge(.6),size=1.9)+interval_overflow(paired_estimated)+
  facet_grid(family~dataset,scales='free_y',space='free_y',labeller=labeller(family=family_labels))+forest_axis()+
  scale_y_discrete(labels=format_feature)+scale_colour_manual(values=c('Whole bulk'='#999999', 'Malignant expression'='#0072B2'))+
  labs(x='Unadjusted hazard ratio per score SD (95% CI)', y=NULL, title='Paired bulk and malignant-expression associations',
    subtitle='Identical patients, retained genes, log2(CPM+1) transformation and Gaussian GSVA in both compartments.')+
  labs(caption='Boundary points and arrowheads indicate estimates or CIs beyond the axis. Complete estimates are tabulated.')+
  theme(strip.text.y=element_text(angle=0,size=10),plot.caption=element_text(size=9,hjust=0))
ggsave(file.path(out, 'figures/supporting/01_bulk_vs_malignant.pdf'), p, width=12, height=11, device=cairo_pdf)

km_colours <- c('At or below median'='#0072B2', 'Above median'='#D55E00')
####################
km_step_bands <- function(d,group_column) {
  do.call(rbind,lapply(split(d,d[[group_column]]),function(z) {
    z <- z[order(z$time),,drop=FALSE]
    previous <- z[-nrow(z),,drop=FALSE]
    previous$time <- z$time[-1L]
    expanded <- rbind(previous,z)
    expanded <- expanded[order(expanded$time),,drop=FALSE]
    expanded[is.finite(expanded$lower)&is.finite(expanded$upper),,drop=FALSE]
  }))
}
####################
km_panel <- function(feature, cohort) {
  d <- cache$km$curves[cache$km$curves$feature==feature & cache$km$curves$cohort==cohort, ]
  r <- cache$km$risks[cache$km$risks$feature==feature & cache$km$risks$cohort==cohort, ]
  start <- d[!duplicated(d$group), ]; start$time <- 0; start$survival <- start$lower <- start$upper <- 1; start$censored <- 0
  d <- rbind(start, d)
  end <- max(d$time)
  time_breaks <- sort(unique(r$time))
  p <- ggplot(d, aes(time, survival, colour=group, group=group))+
####################
    geom_ribbon(data=km_step_bands(d,'group'),aes(ymin=lower,ymax=upper,fill=group),
      alpha=.13,colour=NA,show.legend=FALSE)+scale_fill_manual(values=km_colours)+
####################
    geom_step(linewidth=.65)+geom_point(data=d[d$censored>0, ], shape=3, size=1.3, stroke=.35)+
    scale_colour_manual(values=km_colours)+scale_y_continuous(limits=c(0, 1), labels=scales::label_percent(accuracy=1))+
    scale_x_continuous(limits=c(0, end),breaks=time_breaks)+labs(x=NULL, y='Overall survival', title=cohort_labels[[cohort]],
      subtitle=if (grepl('fraction', feature, fixed=TRUE)) 'Within-malignant RNA fraction' else 'Malignant MP-union score')+
    theme(legend.position='none', plot.title=element_text(size=12))
  risk <- ggplot(r, aes(time, factor(group, levels=c('Above median', 'At or below median')), label=n_risk, colour=group))+
    geom_text(size=3.2, show.legend=FALSE)+scale_colour_manual(values=km_colours)+
    scale_x_continuous(limits=c(0, end),breaks=time_breaks)+labs(x='Years from recorded OS origin', y=NULL)+
    theme_void(base_size=10)+theme(axis.text.y=element_text(size=9), axis.text.x=element_text(), axis.title.x=element_text(),
      plot.margin=margin(0, 12, 5, 9))
  p/risk+plot_layout(heights=c(4, 1))
}
cairo_pdf(file.path(out, 'figures/supporting/03_prespecified_state_km.pdf'), width=12, height=9, onefile=TRUE)
for (state in state_names) {
  panels <- lapply(c(' (MP union)', ' (fraction)'), function(suffix)
    lapply(primary_cohorts, function(c) km_panel(paste0(state, suffix), c)))
  p <- wrap_plots(unlist(panels, recursive=FALSE), ncol=2)+
    plot_annotation(title=state, subtitle='Fixed cohort median split. Blue: at/below median. Orange: above median. Shading: 95% confidence intervals.\nNumbers below each panel are at risk; primary inference uses continuous Cox models.')
  print(p)
}
dev.off()

stage_colours <- c('1'='#0072B2', '2'='#009E73', '3'='#E69F00', '4'='#D55E00')
stage_panels <- lapply(primary_cohorts, function(cohort) {
  d <- cache$stage_km$curves[cache$stage_km$curves$cohort==cohort, ]; r <- cache$stage_km$risks[cache$stage_km$risks$cohort==cohort, ]
####################
  stage_summary <- cache$stage_models$summary[cache$stage_models$summary$cohort==cohort, ]
####################
  start <- d[!duplicated(d$stage), ]; start$time <- 0; start$survival <- start$lower <- start$upper <- 1; start$censored <- 0
  d <- rbind(start, d); end <- max(d$time)
  time_breaks <- sort(unique(r$time))
  p <- ggplot(d, aes(time, survival, colour=stage, group=stage))+
####################
    geom_ribbon(data=km_step_bands(d,'stage'),aes(ymin=lower,ymax=upper,fill=stage),
      alpha=.09,colour=NA,show.legend=FALSE)+scale_fill_manual(values=stage_colours)+
####################
    geom_step(linewidth=.65)+
    geom_point(data=d[d$censored>0, ], shape=3, size=1.5, stroke=.35)+scale_colour_manual(values=stage_colours, name='Pathological stage')+
    scale_x_continuous(limits=c(0, end),breaks=time_breaks)+scale_y_continuous(limits=c(0, 1), labels=scales::label_percent())+
    labs(x=NULL, y='Overall survival', title=cohort_labels[[cohort]],
      subtitle=sprintf('n = %d; deaths = %d; overall stage LR P = %.3g',
        stage_summary$n,stage_summary$events,stage_summary$global_LR_P))
  risk <- ggplot(r, aes(time, factor(stage, levels=as.character(4:1)), label=n_risk, colour=stage))+
    geom_text(size=3.5, show.legend=FALSE)+scale_colour_manual(values=stage_colours)+scale_x_continuous(limits=c(0, end),breaks=time_breaks)+
    labs(x='Years from recorded OS origin', y='Stage')+theme_void()+
    theme(axis.text.y=element_text(), axis.text.x=element_text(), axis.title.x=element_text(), axis.title.y=element_text())
  p/risk+plot_layout(heights=c(4, 1.2))
})
p <- wrap_plots(stage_panels, ncol=2)+plot_annotation(title='Pathological stage in the primary surgery cohorts',
  subtitle='Current workbook AJCC categories; unknown stage excluded. Curves are descriptive; no matching of individual patients.')
ggsave(file.path(out, 'figures/supporting/04_pathological_stage_km.pdf'), p, width=12, height=6.5, device=cairo_pdf)

diagnostics <- main; diagnostics$feature_plot <- factor(diagnostics$feature, levels=rev(c(mp_order, state_order)))
####################
diagnostics$family <- factor(diagnostics$family,levels=names(family_labels))
####################
p <- ggplot(diagnostics, aes(-log10(PH_feature_P), feature_plot, colour=dataset, shape=dataset))+
  geom_blank(aes(x=0))+
  geom_vline(xintercept=-log10(.05), linetype=2, colour='grey55')+
  geom_point(data=diagnostics[diagnostics$status=='estimated', ],position=position_dodge(.5), size=2)+
  facet_grid(family~adjustment,scales='free_y',space='free_y',labeller=labeller(family=family_labels))+
  scale_y_discrete(labels=format_feature)+scale_colour_manual(values=cohort_colours)+
  labs(x='−log10(feature proportional-hazards test P)', y=NULL, title='Proportional-hazards diagnostics',
    subtitle='Dashed line: nominal P = 0.05. Complete PH and eligible spline tests are in cox_results.csv; no automatic model replacement.')+
  theme(strip.text.y=element_text(angle=0,size=10))
ggsave(file.path(out, 'figures/diagnostics/01_model_diagnostics.pdf'), p, width=12, height=11, device=cairo_pdf)

####################
# Four complementary displays requested by the user. Binary analyses are
# exploratory; the stage-adjusted continuous model remains primary.
binary_model <- function(d,method,prob=NA_real_) {
  d <- d[d$OS_eligible, ]
  r <- data.frame(cohort=d$cohort[1],dataset=d$dataset[1],source='malignant',feature=d$feature[1],family=d$family[1],
    adjustment=method,n=nrow(d),events=sum(d$OS_event),low_n=NA_integer_,high_n=NA_integer_,cutpoint=NA_real_,
    q1=NA_real_,q3=NA_real_,probability=prob,status='estimated',logHR=NA_real_,SE=NA_real_,HR=NA_real_,
    lower_CI=NA_real_,upper_CI=NA_real_,P_raw=NA_real_,P_value=NA_real_)
  if(method=='median') {r$cutpoint<-median(d$value);high<-d$value>r$cutpoint}
  if(method=='q1q4') {
    q<-quantile(d$value,c(.25,.75),names=FALSE);r$q1<-q[1];r$q3<-q[2]
    if(q[1]>=q[2]) {r$status<-'not estimated: tied quartiles';return(r)}
    d<-d[d$value<=q[1]|d$value>=q[2], ];high<-d$value>=q[2]
  }
  if(method=='optimal_cut') {r$cutpoint<-unname(quantile(d$value,prob));high<-d$value>=r$cutpoint}
  d$high<-as.integer(high);r$n<-nrow(d);r$events<-sum(d$OS_event);r$low_n<-sum(!high);r$high_n<-sum(high)
  if(r$n<20L || r$events<5L || min(r$low_n,r$high_n)<5L || any(tapply(d$OS_event,d$high,sum)==0)) {
    r$status<-'not estimated: insufficient group/event support';return(r)
  }
  fit<-withCallingHandlers(coxph(Surv(OS_time,OS_event)~high,data=d,ties='efron',control=coxph.control(iter.max=50)),
    warning=function(w) stop(r$cohort,' / ',r$feature,' / ',method,': ',conditionMessage(w)))
  s<-summary(fit);r$logHR<-unname(coef(fit));r$SE<-s$coefficients[1,'se(coef)']
  r$HR<-s$conf.int[1,'exp(coef)'];r$lower_CI<-s$conf.int[1,'lower .95'];r$upper_CI<-s$conf.int[1,'upper .95']
  r$P_raw<-r$P_value<-s$coefficients[1,'Pr(>|z|)']
  stopifnot(all(is.finite(c(r$logHR,r$SE,r$P_value))))
  r
}
four_method_path <- file.path(out,'intermediate/four_method_cache.rds')
four_contract <- list(parent=unname(tools::md5sum(cache_path)),grid=seq(.20,.80,.05),version=1L,
  binary_code=paste(deparse(body(binary_model)),collapse='\n'))
if(file.exists(four_method_path) && !force) {
  four<-readRDS(four_method_path)
  if(!identical(four$contract,four_contract))stop('Four-method cache changed; explicit rebuild required')
} else {
  if(replot)stop('First four-method run requires analysis mode')
  selected<-cache$model_data[cache$model_data$cohort %in% primary_cohorts & cache$model_data$source=='malignant' & cache$model_data$family!='State logratio', ]
  binary<-grids<-list()
  for(d in split(selected,interaction(selected$cohort,selected$feature,drop=TRUE))) {
    binary[[length(binary)+1L]]<-binary_model(d,'median')
    binary[[length(binary)+1L]]<-binary_model(d,'q1q4')
    grid<-do.call(rbind,lapply(four_contract$grid,function(p)binary_model(d,'optimal_cut',p)))
    grids[[length(grids)+1L]]<-grid
    available<-which(grid$status=='estimated')
    if(!length(available))stop('No estimable optimal-cut candidate: ',d$feature[1])
    best<-grid[available[which.min(grid$P_raw[available])], ]
    best$P_value<-min(1,length(four_contract$grid)*best$P_raw)
    binary[[length(binary)+1L]]<-best
  }
  binary<-do.call(rbind,binary)
  continuous<-cache$results[cache$results$cohort %in% primary_cohorts & cache$results$source=='malignant' &
    cache$results$adjustment=='unadjusted' & cache$results$family!='State logratio', ]
  continuous$adjustment<-'continuous';continuous$P_raw<-continuous$P_value
  for(k in setdiff(names(binary),names(continuous)))continuous[[k]]<-NA_real_
  results<-rbind(continuous[,names(binary)],binary)
  results$P_BH<-NA_real_
  for(i in split(seq_len(nrow(results)),interaction(results$cohort,results$adjustment)))results$P_BH[i]<-p.adjust(results$P_value[i],'BH',n=32L)
  four<-list(contract=four_contract,results=results,grid=do.call(rbind,grids),concordance=cross_cohort(results,cache$features$definitions))
  saveRDS(four,four_method_path)
}
table_out(four$results,'four_method_cox_results');table_out(four$grid,'four_method_cutpoint_grid')
table_out(four$concordance,'four_method_concordance')
state_palette<-c('Cell cycle'='#666666',SCREF_STATE_COLOURS)
method_titles<-c(continuous='Continuous Cox model',median='Median split',q1q4='Q1 versus Q4',optimal_cut='Optimal cutpoint (exploratory)')
method_details<-c(continuous='All eligible patients. Expression: HR per 1 SD; fractions: HR per 10 percentage points.',
  median='All eligible patients. HR compares above-median with at-or-below-median values.',
  q1q4='HR compares the upper quartile with the lower quartile. Middle 50% excluded.',
  optimal_cut='HR compares high versus low at the minimum-P cutpoint; 13 candidates from the 20th to 80th percentile.')
format_p<-function(x)ifelse(is.na(x),'NA',ifelse(x<.001,'<0.001',sprintf('%.3f',x)))
four_method_figure<-function(features,path,title) {
  defs<-cache$features$definitions[match(features,cache$features$definitions$feature), ]
  defs$y<-rev(seq_len(nrow(defs)));defs$colour<-unname(state_palette[defs$state])
  defs$display<-ifelse(defs$family=='MP',defs$label,paste0(defs$label,'\n',defs$family))
  defs$display<-vapply(defs$display,function(x)paste(strwrap(x,40),collapse='\n'),character(1))
  n<-nrow(defs)
  base<-theme_void(base_size=11)+theme(plot.title=element_text(size=12,face='bold',hjust=.5),plot.margin=margin(4,3,4,3))
  yscale<-function()scale_y_continuous(limits=c(.4,n+.6),expand=c(0,0))
  cairo_pdf(path,width=12,height=11,onefile=TRUE)
  for(method in names(method_titles)) {
    labels<-ggplot(defs,aes(1,y,label=display,colour=state))+geom_text(hjust=1,size=3.35,lineheight=.95)+
      scale_colour_manual(values=state_palette,name=NULL)+guides(colour=guide_legend(nrow=2,override.aes=list(label='■',size=4)))+
      scale_x_continuous(limits=c(0,1.02))+yscale()+labs(title='MP / state')+base
    panels<-list(labels)
    for(cohort in primary_cohorts) {
      d<-four$results[four$results$cohort==cohort & four$results$adjustment==method & four$results$feature %in% features, ]
      d$y<-defs$y[match(d$feature,defs$feature)];d$colour<-defs$colour[match(d$feature,defs$feature)]
      d$significance<-ifelse(d$status!='estimated','Not estimable',ifelse(d$P_BH<.05,'FDR significant',ifelse(d$P_value<.05,'Nominal only','Not significant')))
      d$statistics<-ifelse(d$status!='estimated','Not estimable',paste0(sprintf('HR %.3g [%.3g, %.3g]',d$HR,d$lower_CI,d$upper_CI),
        '\nP',ifelse(method=='optimal_cut','*',''),'=',format_p(d$P_value),'  q=',format_p(d$P_BH),'\n',d$significance))
      d$point_shape<-ifelse(d$P_BH<.05,18,ifelse(d$P_value<.05,16,1))
      z<-interval_coordinates(d[d$status=='estimated', ])
      forest<-ggplot(d,aes(y=y,colour=colour))+geom_vline(xintercept=0,linetype=2,colour='grey50')+
        geom_errorbar(data=z,aes(x=effect,xmin=left,xmax=right),orientation='y',width=.12,linewidth=.55)+
        geom_point(data=z,aes(x=effect,shape=point_shape),size=2.4)+scale_shape_identity()+scale_colour_identity()+
        geom_text(data=z,aes(x=-3,label=ifelse(left_overflow,'<','')),size=3.5)+
        geom_text(data=z,aes(x=3,label=ifelse(right_overflow,'>','')),size=3.5)+
        scale_x_continuous(limits=c(-3.15,3.15),breaks=c(-2,0,2),labels=c('0.25','1','4'))+yscale()+labs(title=paste0(ifelse(cohort=='occams_surgery','OCCAMS','TCGA'),'\nn=',paste(unique(range(d$n)),collapse='–')),x='Hazard ratio')+base+
        theme(axis.text.x=element_text(size=9),axis.title.x=element_text(size=10),axis.line.x=element_line())
      stats<-ggplot(d,aes(0,y,label=statistics))+geom_text(hjust=0,size=3.05,lineheight=.95)+
        scale_x_continuous(limits=c(0,1))+yscale()+labs(title='HR [95% CI]; significance')+base+theme(plot.title=element_text(size=10))
      panels<-c(panels,list(forest,stats))
    }
    p<-wrap_plots(panels,nrow=1,widths=c(3.5,2,2.3,2,2.3),guides='collect')+plot_annotation(
      title=paste(title,'|',method_titles[[method]]),subtitle=paste0(method_details[[method]],'\nUnadjusted comparisons; stage-adjusted continuous inference is retained separately.'),
      caption=paste0('HR < 1: lower hazard / better survival. HR > 1: higher hazard / poorer survival. Bars show 95% CI; boundary symbols indicate values beyond the axis.\n',
        'Open circle: not significant (P >= 0.05). Filled circle: nominal only (P < 0.05, q >= 0.05). Diamond: FDR significant (q < 0.05 across 32 features per method).\n',
        ifelse(method=='optimal_cut','P* corrects 13 searches; q is BH FDR. Selected HR/CI are descriptive; CIs do not account for cutpoint selection.','P is two-sided Wald P; q is within-cohort BH FDR. Binary methods are exploratory; comparing four methods does not establish replication.')),
      theme=theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=11),plot.caption=element_text(size=9,hjust=0),plot.margin=margin(9,9,9,9)))
    print(p)
  }
  dev.off()
}
four_method_figure(mp_order,file.path(out,'figures/supporting/07_four_method_mp_forests.pdf'),'Malignant MP survival')
four_method_figure(state_order,file.path(out,'figures/supporting/08_four_method_state_forests.pdf'),'Malignant state survival')
####################

####################
four_method_overview <- function(features,path,title) {
  defs<-cache$features$definitions[match(features,cache$features$definitions$feature), ]
  defs$y<-rev(seq_len(nrow(defs)))
  defs$display<-ifelse(defs$family=='MP',defs$label,paste(defs$label,sub('State ','',defs$family),sep=' | '))
  defs$display<-vapply(defs$display,function(x)paste(strwrap(x,38),collapse='\n'),character(1))
  n<-nrow(defs)
  d<-four$results[four$results$feature %in% features, ]
  d$y<-defs$y[match(d$feature,defs$feature)]
  d$dataset<-factor(d$dataset,levels=c('OCCAMS','TCGA'))
  d$method<-factor(d$adjustment,levels=names(method_titles),labels=c('Continuous','Median split','Q1 vs Q4','Optimal cut*'))
  d$significance<-ifelse(d$status!='estimated','Not estimated',ifelse(d$P_BH<.05,'FDR',ifelse(d$P_value<.05,'Nominal','NS')))
  d$display<-ifelse(d$status!='estimated','Not\nestimated',paste0(sprintf('%.3g',d$HR),' (',d$significance,')\nP',ifelse(d$adjustment=='optimal_cut','*',''),' ',format_p(d$P_value),
    '\nq ',format_p(d$P_BH)))
  d$trend<-ifelse(d$status!='estimated','Not estimated',ifelse(d$logHR>0,'Higher hazard / poorer survival','Lower hazard / better survival'))
  labels<-ggplot(defs,aes(1,y,label=display,colour=state))+geom_text(hjust=1,size=3.25,lineheight=.95)+
    scale_colour_manual(values=state_palette,name=NULL)+guides(colour=guide_legend(nrow=2,override.aes=list(label='■',size=4)))+
    scale_x_continuous(limits=c(0,1.02))+scale_y_continuous(limits=c(.4,n+.6),expand=c(0,0))+
    theme_void()+theme(legend.position='bottom',legend.text=element_text(size=10),plot.margin=margin(27,5,5,5))
  values<-ggplot(d,aes(dataset,y,fill=trend))+geom_tile(width=.96,height=.96,colour='white')+
    geom_text(aes(label=display),size=3.05,lineheight=.98)+
    facet_grid(.~method)+scale_y_continuous(limits=c(.4,n+.6),expand=c(0,0))+
    scale_x_discrete(position='top',expand=expansion(add=.5))+
    scale_fill_manual(values=c('Higher hazard / poorer survival'='#FBE2DD','Lower hazard / better survival'='#DFEDF7','Not estimated'='#EEEEEE'),name=NULL)+
    theme_minimal(base_size=11)+theme(panel.grid=element_blank(),axis.title=element_blank(),axis.text.y=element_blank(),
      axis.text.x=element_text(size=10,colour='black'),strip.text=element_text(size=11,face='bold'),legend.position='bottom',
      legend.text=element_text(size=10),plot.margin=margin(5,5,5,5),panel.spacing.x=grid::unit(.45,'lines'))
  p<-wrap_plots(labels,values,nrow=1,widths=c(3.25,8.75),guides='collect')+plot_annotation(
    title=title,subtitle='Unadjusted surgery-cohort comparisons: OCCAMS and TCGA. HR > 1 = poorer survival; HR < 1 = better survival.\nContinuous: per score SD / 10 percentage points of fraction. Binary methods: high versus low; Q1/Q4 excludes the middle half.',
    caption='Each cell: HR (significance), P, q. NS: P >= 0.05. Nominal: P < 0.05 but q >= 0.05. FDR: q < 0.05 across 32 features.\n*Optimal cut: P* corrects for 13 searched thresholds; q is calculated after that correction. Raw minimum P is tabulated separately.\nFull 95% CIs, group sizes and cutpoints are in supporting forests 07/08 and source tables. Stage-adjusted continuous models remain primary.',
    theme=theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10.5),plot.caption=element_text(size=9,hjust=0)))
  ggsave(path,p,width=12,height=11,device=cairo_pdf)
}
four_method_overview(mp_order,file.path(out,'figures/final/01_malignant_mp_survival.pdf'),'Malignant MPs: four survival analyses')
four_method_overview(state_order,file.path(out,'figures/final/02_malignant_state_survival.pdf'),'Malignant states: four survival analyses')
####################

####################
# Reuse the verified parent scores; the combined clinical selection is an
# additional subset, not a new GSVA reference population.
context_path <- file.path(out,'intermediate/occams_context_cache.rds')
context_cohorts <- c('occams_combined','occams_neo','occams_neo_biopsy','occams_neo_surgery','occams_surgery')
context_contract <- list(parent=unname(tools::md5sum(cache_path)),four=four_contract,
  cohorts=context_cohorts,version=1L)
if (file.exists(context_path) && !force) {
  context <- readRDS(context_path)
  if (!identical(context$contract,context_contract)) stop('OCCAMS context cache changed; explicit rebuild required')
} else {
  if (replot) stop('First OCCAMS context run requires analysis mode')
  combined_ids <- unique(c(cache$clinical$cohorts$occams_neo,cache$clinical$cohorts$occams_surgery))
  previous <- read.csv(inputs[['prior_membership']],stringsAsFactors=FALSE)
  stopifnot(length(combined_ids)==225L,
    setequal(combined_ids,previous$subject_id[previous$cohort=='occams_combined']))
  combined <- build_model_data(cache$clinical$meta,list(occams_combined=combined_ids),cache$features,cache$scored)
  selected <- rbind(cache$model_data[cache$model_data$cohort %in% context_cohorts,],combined)
  selected <- selected[selected$source=='malignant' & selected$family!='State logratio',]
  continuous <- cache$results[cache$results$cohort %in% context_cohorts & cache$results$source=='malignant' &
    cache$results$adjustment=='unadjusted' & cache$results$family!='State logratio',]
  additional <- selected[selected$cohort=='occams_combined',]
  additional_results <- do.call(rbind,lapply(split(additional,additional$feature),function(d)fit_one(d,'unadjusted')))
  for (column in setdiff(names(continuous),names(additional_results))) additional_results[[column]] <- NA_real_
  continuous <- rbind(continuous,additional_results[,names(continuous)])
  binary <- grids <- list()
  for (d in split(selected,interaction(selected$cohort,selected$feature,drop=TRUE))) {
    binary[[length(binary)+1L]] <- binary_model(d,'median')
    binary[[length(binary)+1L]] <- binary_model(d,'q1q4')
    grid <- do.call(rbind,lapply(four_contract$grid,function(p)binary_model(d,'optimal_cut',p)))
    grids[[length(grids)+1L]] <- grid
    available <- which(grid$status=='estimated')
    if (!length(available)) stop('No estimable optimal-cut candidate: ',d$cohort[1],' / ',d$feature[1])
    best <- grid[available[which.min(grid$P_raw[available])],]
    best$P_value <- min(1,length(four_contract$grid)*best$P_raw)
    binary[[length(binary)+1L]] <- best
  }
  binary <- do.call(rbind,binary)
  continuous$adjustment <- 'continuous'; continuous$P_raw <- continuous$P_value
  for (column in setdiff(names(binary),names(continuous))) continuous[[column]] <- NA_real_
  results <- rbind(continuous[,names(binary)],binary)
  results$P_BH <- NA_real_
  for (i in split(seq_len(nrow(results)),interaction(results$cohort,results$adjustment)))
    results$P_BH[i] <- p.adjust(results$P_value[i],'BH',n=32L)
  stopifnot(nrow(results)==640L)
  context <- list(contract=context_contract,results=results,grid=do.call(rbind,grids),model_data=selected,
    membership=unique(selected[,c('cohort','subject_id','sample_id','OS_eligible')]))
  saveRDS(context,context_path)
}
####################
# RFS is the exploratory endpoint used by the earlier selected-cohort workflow.
# Reuse the fixed malignant scores; endpoint eligibility determines score SDs.
rfs_path <- file.path(out,'intermediate/rfs_analysis_cache.rds')
rfs_manifest <- data.frame(input=names(rfs_inputs),path=normalizePath(rfs_inputs),
  md5=unname(tools::md5sum(rfs_inputs)))
rfs_cohorts <- cache$clinical$cohorts
rfs_cohorts$occams_combined <- unique(c(rfs_cohorts$occams_neo,rfs_cohorts$occams_surgery))
####################
# Use a fixed tighter convergence tolerance for every binary RFS fit.
# Near-null coefficients can otherwise trigger the Cox infinite-coefficient heuristic.
rfs_binary_model <- function(d,method,prob=NA_real_) {
  d <- d[d$RFS_eligible, ]
  r <- data.frame(cohort=d$cohort[1],dataset=d$dataset[1],source='malignant',feature=d$feature[1],family=d$family[1],
    adjustment=method,n=nrow(d),events=sum(d$RFS_event),low_n=NA_integer_,high_n=NA_integer_,cutpoint=NA_real_,
    q1=NA_real_,q3=NA_real_,probability=prob,status='estimated',logHR=NA_real_,SE=NA_real_,HR=NA_real_,
    lower_CI=NA_real_,upper_CI=NA_real_,P_raw=NA_real_,P_value=NA_real_)
  if(method=='median') {r$cutpoint<-median(d$value);high<-d$value>r$cutpoint}
  if(method=='q1q4') {
    q<-quantile(d$value,c(.25,.75),names=FALSE);r$q1<-q[1];r$q3<-q[2]
    if(q[1]>=q[2]) {r$status<-'not estimated: tied quartiles';return(r)}
    d<-d[d$value<=q[1]|d$value>=q[2], ];high<-d$value>=q[2]
  }
  if(method=='optimal_cut') {r$cutpoint<-unname(quantile(d$value,prob));high<-d$value>=r$cutpoint}
  d$high<-as.integer(high);r$n<-nrow(d);r$events<-sum(d$RFS_event);r$low_n<-sum(!high);r$high_n<-sum(high)
  if(r$n<20L || r$events<5L || min(r$low_n,r$high_n)<5L || any(tapply(d$RFS_event,d$high,sum)==0)) {
    r$status<-'not estimated: insufficient group/event support';return(r)
  }
  fit<-withCallingHandlers(coxph(Surv(RFS_time,RFS_event)~high,data=d,ties='efron',control=coxph.control(iter.max=50,eps=1e-12,toler.chol=1e-14)),
    warning=function(w) stop(r$cohort,' / ',r$feature,' / ',method,': ',conditionMessage(w)))
  s<-summary(fit);r$logHR<-unname(coef(fit));r$SE<-s$coefficients[1,'se(coef)']
  r$HR<-s$conf.int[1,'exp(coef)'];r$lower_CI<-s$conf.int[1,'lower .95'];r$upper_CI<-s$conf.int[1,'upper .95']
  r$P_raw<-r$P_value<-s$coefficients[1,'Pr(>|z|)']
  stopifnot(all(is.finite(c(r$logHR,r$SE,r$P_value))))
  r
}
####################
build_rfs_analysis <- function() {
  metadata <- read_rfs_metadata(cache$clinical$meta[cache$clinical$meta$subject_id %in%
    unique(unlist(rfs_cohorts)),])
  endpoint_metadata <- metadata
  endpoint_metadata$OS_time <- metadata$RFS_time
  endpoint_metadata$OS_event <- metadata$RFS_event
  endpoint_metadata$OS_eligible <- metadata$RFS_eligible
  endpoint_metadata$OS_exclusion <- metadata$RFS_exclusion
  model_cohorts <- rfs_cohorts[unique(c(context_cohorts,primary_cohorts))]
  selected <- build_model_data(endpoint_metadata,model_cohorts,cache$features,cache$scored)
  selected <- selected[selected$source=='malignant' & selected$family!='State logratio',]
  continuous <- do.call(rbind,lapply(split(selected,interaction(selected$cohort,selected$feature,drop=TRUE)),function(d) {
    adjustments <- if(d$cohort[1] %in% primary_cohorts) c('unadjusted','stage') else 'unadjusted'
    do.call(rbind,lapply(adjustments,function(a)fit_one(d,a)))
  }))
  continuous$status <- gsub('deaths','events',continuous$status,fixed=TRUE)
  continuous$P_BH <- NA_real_
  for(i in split(seq_len(nrow(continuous)),interaction(continuous$cohort,continuous$adjustment,drop=TRUE)))
    continuous$P_BH[i] <- p.adjust(continuous$P_value[i],'BH',n=32L)
  binary <- grids <- list()
  for(d in split(selected[selected$cohort %in% context_cohorts,],
    interaction(selected$cohort[selected$cohort %in% context_cohorts],
      selected$feature[selected$cohort %in% context_cohorts],drop=TRUE))) {
    binary[[length(binary)+1L]] <- rfs_binary_model(d,'median')
    binary[[length(binary)+1L]] <- rfs_binary_model(d,'q1q4')
    grid <- do.call(rbind,lapply(four_contract$grid,function(p)rfs_binary_model(d,'optimal_cut',p)))
    grids[[length(grids)+1L]] <- grid
    available <- which(grid$status=='estimated')
    if(!length(available)) stop('No estimable RFS optimal-cut candidate: ',d$cohort[1],' / ',d$feature[1])
    best <- grid[available[which.min(grid$P_raw[available])],]
    best$P_value <- min(1,length(four_contract$grid)*best$P_raw)
    binary[[length(binary)+1L]] <- best
  }
  binary <- do.call(rbind,binary)
  context_continuous <- continuous[continuous$cohort %in% context_cohorts & continuous$adjustment=='unadjusted',]
  context_continuous$adjustment <- 'continuous'
  context_continuous$P_raw <- context_continuous$P_value
  for(column in setdiff(names(binary),names(context_continuous))) context_continuous[[column]] <- NA_real_
  context_results <- rbind(context_continuous[,names(binary)],binary)
  context_results$P_BH <- NA_real_
  for(i in split(seq_len(nrow(context_results)),interaction(context_results$cohort,context_results$adjustment,drop=TRUE)))
    context_results$P_BH[i] <- p.adjust(context_results$P_value[i],'BH',n=32L)
  stopifnot(nrow(context_results)==640L,nrow(continuous)==256L)
  list(metadata=metadata,model_data=selected,continuous=continuous,
    context_results=context_results,grid=do.call(rbind,grids),
    concordance=cross_cohort(continuous[continuous$cohort %in% primary_cohorts,],cache$features$definitions))
}
rfs_contract <- list(parent=unname(tools::md5sum(cache_path)),inputs=rfs_manifest,
  four=four_contract,cohorts=rfs_cohorts,version=1L,
  endpoint_code=paste(deparse(body(read_rfs_metadata)),collapse='\n'),
  analysis_code=paste(deparse(body(build_rfs_analysis)),collapse='\n'),
  binary_code=paste(deparse(body(rfs_binary_model)),collapse='\n'))
if(file.exists(rfs_path) && !force) {
  rfs <- readRDS(rfs_path)
  if(!identical(rfs$contract,rfs_contract)) stop('RFS cache changed; explicit rebuild required')
} else {
  if(replot) stop('First RFS run requires analysis mode')
  rfs <- build_rfs_analysis()
  rfs$contract <- rfs_contract
  saveRDS(rfs,rfs_path)
}
endpoint_source <- function(d,endpoint) {
  names(d)[match(c('OS_time','OS_event','OS_eligible','OS_exclusion'),names(d))] <-
    c('endpoint_time','endpoint_event','endpoint_eligible','endpoint_exclusion')
  d$endpoint <- endpoint
  d
}
table_out(rfs_manifest,'rfs_input_manifest')
table_out(rfs$metadata,'rfs_patient_metadata')
table_out(endpoint_source(rfs$model_data,'RFS'),'rfs_model_data')
table_out(transform(rfs$continuous,endpoint='RFS'),'rfs_cox_results')
table_out(transform(rfs$concordance,endpoint='RFS'),'rfs_concordance_results')
endpoint_summary <- do.call(rbind,lapply(names(rfs_cohorts),function(key) {
  d <- rfs$metadata[match(rfs_cohorts[[key]],rfs$metadata$subject_id),]
  do.call(rbind,lapply(c('OS','RFS'),function(endpoint) {
    eligible <- d[[paste0(endpoint,'_eligible')]]
    events <- d[[paste0(endpoint,'_event')]]
    stage_complete <- eligible & !is.na(d$pstage)
    data.frame(cohort=key,endpoint=endpoint,n_selected=nrow(d),n=sum(eligible),
      events=sum(events[eligible]),n_stage=sum(stage_complete),events_stage=sum(events[stage_complete]))
  }))
}))
table_out(endpoint_summary,'endpoint_availability')
context_results <- rbind(transform(context$results,endpoint='OS'),transform(rfs$context_results,endpoint='RFS'))
context_grid <- rbind(transform(context$grid,endpoint='OS'),transform(rfs$grid,endpoint='RFS'))
os_context_source <- endpoint_source(context$model_data,'OS')
rfs_context_source <- endpoint_source(rfs$model_data[rfs$model_data$cohort %in% context_cohorts,],'RFS')
# Export a common patient-level schema; RFS-specific provenance is in rfs_patient_metadata.
context_source <- rbind(os_context_source,rfs_context_source[,names(os_context_source)])
table_out(context_results,'occams_context_cox_results')
table_out(context_grid,'occams_context_cutpoint_grid')
table_out(context_source,'occams_context_model_data')
table_out(unique(context_source[,c('cohort','subject_id','sample_id','endpoint','endpoint_eligible')]),
  'occams_context_membership')

for(endpoint in c('OS','RFS')) {
  d <- if(endpoint=='OS') main else rfs$continuous[rfs$continuous$cohort %in% primary_cohorts,]
  filename <- if(endpoint=='OS') '05_continuous_os_adjustment.pdf' else '06_continuous_rfs_adjustment.pdf'
  cairo_pdf(file.path(out,'figures/supporting',filename),width=12,height=11,onefile=TRUE)
  print(forest_plot(d[d$family=='MP',],mp_order,paste(endpoint,'| MPs'),NULL))
  print(forest_plot(d[d$family!='MP',],state_order,paste(endpoint,'| States'),NULL,'family_adjustment'))
  dev.off()
}
####################

context_panels <- list(
  specimens=list(title='Neoadjuvant specimens',cohorts=c('occams_neo_biopsy','occams_neo_surgery'),
    labels=c('Biopsy','Surgery'),colours=c('#009E73','#CC79A7')),
  pathways=list(title='OCCAMS treatment pathways',cohorts=c('occams_neo','occams_surgery'),
    labels=c('Neoadjuvant','Surgery pathway'),colours=c('#0072B2','#D55E00')),
  combined=list(title='OCCAMS combined pathways',cohorts='occams_combined',
    labels='Combined',colours='#333333'))
context_method_labels <- c(continuous='Continuous',median='Median',q1q4='Q1 vs Q4',optimal_cut='Optimal cut')
context_page_index <- list()
cairo_pdf(file.path(out,'figures/supporting/02_neoadjuvant_context.pdf'),width=12,height=11,onefile=TRUE)
for (comparison in names(context_panels)) {
  panel <- context_panels[[comparison]]
  for (method in names(context_method_labels)) {
    d <- context_results[context_results$cohort %in% panel$cohorts & context_results$adjustment==method,]
    d$endpoint <- factor(d$endpoint,levels=c('OS','RFS'))
    d$feature_plot <- factor(d$feature,levels=rev(c(mp_order,state_order)))
    d$family <- factor(d$family,levels=names(family_labels))
    d$cohort_plot <- factor(d$cohort,levels=panel$cohorts,labels=panel$labels)
    d$significance <- factor(ifelse(d$P_BH<.05,'FDR',ifelse(d$P_value<.05,'Nominal','NS')),
      levels=c('NS','Nominal','FDR'))
    estimated <- interval_coordinates(d[d$status=='estimated',])
    endpoint_labels <- setNames(vapply(c('OS','RFS'),function(endpoint)
      paste(endpoint,panel_count_label(d[d$endpoint==endpoint,],'cohort',panel$cohorts,panel$labels),
        sep='\n'),character(1)),c('OS','RFS'))
    p <- ggplot(d,aes(y=feature_plot,colour=cohort_plot,group=cohort_plot))+
      geom_blank(aes(x=0))+
      geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.4)+
      geom_errorbar(data=estimated,aes(x=effect,xmin=left,xmax=right),orientation='y',
        width=.12,position=position_dodge(.6),linewidth=.5)+
      geom_point(data=estimated,aes(x=effect,shape=significance),position=position_dodge(.6),size=2.3,show.legend=TRUE)+
      interval_overflow(estimated)+
      facet_grid(family~endpoint,scales='free_y',space='free_y',
        labeller=labeller(family=family_labels,endpoint=endpoint_labels))+
      forest_axis()+scale_y_discrete(labels=format_feature)+
      scale_colour_manual(values=setNames(panel$colours,panel$labels),name=NULL)+
      scale_shape_manual(values=c(NS=1,Nominal=16,FDR=18),drop=FALSE,name='Evidence')+
      labs(x=if(method=='continuous')'Hazard ratio per SD / 10 percentage points' else 'Hazard ratio: high vs low',
        y=NULL,title=paste(panel$title,'|',context_method_labels[[method]]))+
      theme(strip.text.y=element_text(angle=0,size=10),plot.title=element_text(size=16,face='bold'),
        legend.text=element_text(size=11),axis.text.y=element_text(size=10),
        strip.text.x=element_text(size=10.5,face='plain',lineheight=1.1))
    print(p)
    context_page_index[[length(context_page_index)+1L]] <- data.frame(page=length(context_page_index)+1L,
      comparison=comparison,method=method,cohorts=paste(panel$cohorts,collapse=' | '),endpoints='OS | RFS',features=32L)
  }
}
dev.off()
context_page_index <- do.call(rbind,context_page_index)
table_out(context_page_index,'occams_context_figure_index')
####################

figure_index <- data.frame(role=c(rep('final', 2), rep('supporting', 4), 'diagnostics'),
  path=c('figures/final/01_malignant_mp_survival.pdf', 'figures/final/02_malignant_state_survival.pdf',
    'figures/supporting/01_bulk_vs_malignant.pdf', 'figures/supporting/02_neoadjuvant_context.pdf',
    'figures/supporting/03_prespecified_state_km.pdf', 'figures/supporting/04_pathological_stage_km.pdf',
    'figures/diagnostics/01_model_diagnostics.pdf'),
  pages=c(1L, 1L, 1L, nrow(context_page_index), 5L, 1L, 1L),
  purpose=c('All 17 malignant MP estimates in both surgery cohorts', 'All 15 state score/fraction estimates in both surgery cohorts',
    'Paired bulk-versus-malignant estimates on identical patients', 'OCCAMS specimens, treatment pathways and combined selection; OS/RFS panels and four methods on separate pages',
    'Prespecified median state curves with risk counts', 'Current pathological stage curves and risk counts', 'PH diagnostic overview'))
figure_index <- rbind(figure_index,data.frame(role='supporting',path=c('figures/supporting/05_continuous_os_adjustment.pdf','figures/supporting/06_continuous_rfs_adjustment.pdf'),pages=2L,purpose=c('OS: MPs then states, unadjusted and stage adjusted','Exploratory RFS: MPs then states, unadjusted and stage adjusted')))
figure_index <- rbind(figure_index,data.frame(role='supporting',path=c('figures/supporting/07_four_method_mp_forests.pdf','figures/supporting/08_four_method_state_forests.pdf'),pages=4L,purpose=c('Four-method MP estimates and confidence intervals','Four-method state estimates and confidence intervals')))
table_out(figure_index, 'figure_index')
stopifnot(all(file.exists(file.path(out, figure_index$path))), all(file.info(file.path(out, figure_index$path))$size>1000))
primary_replication <- cache$concordance[cache$concordance$source=='malignant' & cache$concordance$family!='State logratio', ]
summary <- do.call(rbind, lapply(c('unadjusted', 'stage_complete_unadjusted', 'stage', 'stage_purity', 'age_sex_stage'), function(a) {
  x <- primary_replication[primary_replication$adjustment==a, ]
  data.frame(adjustment=a, features=32L, jointly_estimable=sum(is.finite(x$conjunction_P)),
    concordant=sum(x$concordant), nominal_replication=sum(x$nominal_replication), FDR_replication=sum(x$FDR_replication),
####################
    nominal_expected_direction=sum(x$nominal_replication & x$matches_hypothesis,na.rm=TRUE),
    FDR_expected_direction=sum(x$FDR_replication & x$matches_hypothesis,na.rm=TRUE),
####################
    FDR_higher_hazard=sum(x$FDR_replication & x$direction=='higher hazard in both'),
    FDR_lower_hazard=sum(x$FDR_replication & x$direction=='lower hazard in both'))}))
validation <- data.frame(check=c('OCCAMS parent', 'TCGA parent', 'OCCAMS surgery membership', 'TCGA surgery membership',
  'main features per malignant cohort/model', 'minimum retained signature genes', 'final figure pages'),
  observed=c(281L, 88L, 68L, 72L, 32L, min(cache$scored$coverage$retained_genes), 2L),
  result='pass')
stopifnot(all(table(main$cohort, main$adjustment)==32L))
write.csv(validation, file.path(out, 'reports/validation.csv'), row.names=FALSE)
dir.create(SCREF_SUMMARY_DIR, recursive=TRUE, showWarnings=FALSE)
write.csv(summary, file.path(SCREF_SUMMARY_DIR, 'malignant_bulk_survival.csv'), row.names=FALSE)
####################
stage_support <- primary_replication[primary_replication$adjustment=='stage' & primary_replication$nominal_replication, ]
stage_support <- stage_support[order(stage_support$conjunction_BH,stage_support$feature), ]
support_text <- if (!nrow(stage_support)) {
  'No feature has a concordant nominal P < 0.05 in both cohorts in the primary stage-adjusted analysis.'
} else {
  c('Features with concordant nominal P < 0.05 in both cohorts in the primary stage-adjusted analysis:', '',
    vapply(seq_len(nrow(stage_support)),function(i) {
      z <- stage_support[i, ]
      sprintf('- %s: %s; OCCAMS HR %.3g (P %.3g), TCGA HR %.3g (P %.3g); cross-cohort BH q %.3g. %s',
        z$feature,z$direction,z$HR_OCCAMS,z$P_value_OCCAMS,z$HR_TCGA,z$P_value_TCGA,z$conjunction_BH,
        if (z$FDR_replication) 'Meets the prespecified cross-cohort FDR threshold.' else 'Does not meet the cross-cohort FDR threshold.')
    },character(1)))
}
####################
####################
summary_labels <- c(unadjusted='Unadjusted',stage_complete_unadjusted='Unadjusted, stage-complete cases',
  stage='Stage adjusted (primary)',stage_purity='Stage and inferred purity adjusted',age_sex_stage='Age, sex and stage adjusted')
summary_markdown <- c('| Model | Same direction / 32 | Nominal P < 0.05 in both | Cross-cohort BH q < 0.05 |',
  '|---|---:|---:|---:|',vapply(seq_len(nrow(summary)),function(i) {
    z <- summary[i, ]
    sprintf('| %s | %d | %d | %d |',summary_labels[[z$adjustment]],z$concordant,z$nominal_replication,z$FDR_replication)
  },character(1)))
unadjusted_support <- primary_replication[primary_replication$adjustment=='unadjusted' & primary_replication$nominal_replication, ]
unadjusted_text <- if (!nrow(unadjusted_support)) 'No unadjusted feature has concordant nominal P < 0.05 in both cohorts.' else
  c('Unadjusted findings are secondary and do not replace the primary stage-adjusted result:', '',
    vapply(seq_len(nrow(unadjusted_support)),function(i) {
      z <- unadjusted_support[i, ]
      sprintf('- %s: %s; OCCAMS HR %.3g (P %.3g), TCGA HR %.3g (P %.3g); cross-cohort BH q %.3g.',
        z$feature,z$direction,z$HR_OCCAMS,z$P_value_OCCAMS,z$HR_TCGA,z$P_value_TCGA,z$conjunction_BH)
    },character(1)), '',
    'The stress-adaptive MP-union signature is exactly MP12. Those two labels are the same test, not independent supporting findings. The fixed 32-feature correction retains these correlated and duplicate definitions conservatively.')
####################
####################
paired_expression_text <- c('Paired comparison restricted to the same 27 expression features:', '',
  '| Model | Expression source | Same direction / 27 | Nominal P < 0.05 in both | Cross-cohort BH q < 0.05 |',
  '|---|---|---:|---:|---:|')
for (adjustment in c('unadjusted','stage')) for (source in c('bulk','malignant')) {
  z <- cache$concordance[cache$concordance$adjustment==adjustment & cache$concordance$source==source &
    cache$concordance$family %in% c('MP','State MP union','State DGE'), ]
  stopifnot(nrow(z)==27L)
  paired_expression_text <- c(paired_expression_text,sprintf('| %s | %s | %d | %d | %d |',
    summary_labels[[adjustment]],ifelse(source=='bulk','Whole bulk','Malignant expression'),
    sum(z$concordant),sum(z$nominal_replication),sum(z$FDR_replication)))
}
####################
####################
four_method_text <- c('Exploratory unadjusted OS comparisons across the four displayed methods:', '',
  '| Method | Same direction / 32 | Nominal P < 0.05 in both | Cross-cohort BH q < 0.05 |',
  '|---|---:|---:|---:|')
for (method in names(method_titles)) {
  z <- four$concordance[four$concordance$adjustment==method, ]
  four_method_text <- c(four_method_text, sprintf('| %s | %d | %d | %d |',
    method_titles[[method]], sum(z$concordant), sum(z$nominal_replication), sum(z$FDR_replication)))
}
four_method_text <- c(four_method_text, '',
  'These method-specific q values do not correct for choosing among four methods. No method gives FDR-supported cross-cohort replication. Cohort-specific significance therefore must not be interpreted as replicated evidence.', '',
  '[Four-method estimates](../tables/four_method_cox_results.csv) | [Four-method concordance](../tables/four_method_concordance.csv) | [Cutpoint search grid](../tables/four_method_cutpoint_grid.csv)')
####################
####################
rfs_text <- c('Exploratory RFS in the primary surgery cohorts:', '',
  '| Model | Same direction / 32 | Nominal P < 0.05 in both | Cross-cohort BH q < 0.05 |',
  '|---|---:|---:|---:|')
for(adjustment in c('unadjusted','stage')) {
  z <- rfs$concordance[rfs$concordance$adjustment==adjustment,]
  rfs_text <- c(rfs_text,sprintf('| %s | %d | %d | %d |',adjustment,
    sum(z$concordant),sum(z$nominal_replication),sum(z$FDR_replication)))
}
rfs_support <- rfs$concordance[rfs$concordance$nominal_replication,]
rfs_text <- c(rfs_text,'',
  if(!nrow(rfs_support)) 'No RFS feature has nominal significance in both cohorts.' else
    vapply(seq_len(nrow(rfs_support)),function(i) {
      z <- rfs_support[i,]
      sprintf('- %s (%s): %s; OCCAMS P %.3g, TCGA P %.3g; cross-cohort BH q %.3g.',
        z$feature,z$adjustment,z$direction,z$P_value_OCCAMS,z$P_value_TCGA,z$conjunction_BH)
    },character(1)), '',
  'Nominal associations do not establish FDR-supported replication. MP12 and the stress-adaptive MP union are identical tests, not independent findings. RFS remains exploratory because the historical endpoint definitions and origins differ across cohorts.',
  '', '[RFS cross-cohort results](../tables/rfs_concordance_results.csv)')
####################
report <- c('# Current malignant survival analysis', '',
  'Primary comparison: the exact OCCAMS curative-intent surgery-pathway 68 and TCGA surgery-Yes 72 subject lists. Endpoint and covariate exclusions are explicit in the tables. These are comparable clinical subsets, not individually matched patients or proven equivalent treatment histories.', '',
  'OCCAMS metadata come from the 25 September 2026 workbook. TCGA stage, surgery and sample selection come from the 18 September 2026 workbook; OS, age and sex come from the local cBioPortal clinical table because the workbook does not contain them. Exact GDC sample mapping is independently checked. Current patient lists match the previous selected-cohort analysis.', '',
  'Both expression compartments use the same genes, log2(CPM+1) transformation, Gaussian GSVA and complete parent populations (281 OCCAMS, 88 TCGA). The 32-feature primary family comprises 17 MPs, five MP-union state scores, five DGE state scores, and five within-malignant state RNA fractions. Fractions retain technical and unassigned states in their denominator. State logratios are a separate five-feature sensitivity family.', '',
  'The primary model is adjusted for categorical pathological stage. Continuous Cox models use one score SD or 10 percentage points of state fraction. Unadjusted, stage-complete unadjusted, stage/purity and age/sex/stage models are prespecified sensitivity analyses. Feasibility rules can leave an adjusted estimate unreported; the script never replaces it with a simpler model. Unexpected fitting warnings stop execution. Complete model status, exclusions, PH tests and eligible spline tests are retained.', '',
  'Multiplicity: BH correction across all 32 primary features within each cohort, compartment and adjustment. Concordant cross-cohort conjunction P is the larger two-sided cohort P; discordant estimates receive P=1. BH correction across 32 conjunction tests defines FDR replication. Expected biological directions are hypotheses, never filters or selection rules.', '',
  summary_markdown, '',
####################
  support_text, '',unadjusted_text, '',paired_expression_text, '',four_method_text, '',rfs_text, '',
####################
  'Each final PDF has one page comparing four explicitly labelled methods: continuous, median, Q1/Q4 and optimal cutpoint. State colours, numerical P/q values and explicit significance labels accompany every feature. Optimal-cutpoint P values are corrected for all 13 searched thresholds before BH correction. Supporting 05 contains OS and supporting 06 exploratory RFS, each with MPs and states on separate pages and unadjusted/stage-adjusted panels. The OS stage-adjusted models remain primary. The final overviews show all features; supporting forests 07/08 provide the corresponding 95% confidence intervals. No signature is selected because its direction matches an expectation. Descriptive median Kaplan–Meier curves do not provide the primary significance test. Bulk-versus-malignant differences are paired descriptive comparisons; their effect difference is not tested as if independent. Both cohorts were examined previously: concordant findings provide retrospective cross-cohort support, not untouched external validation.', '',
  'Limitations: limited deaths constrain adjustment and spline diagnostics. Broad pathological stages can retain residual confounding and differ after neoadjuvant treatment. OCCAMS uses its recorded survival interval; TCGA uses the cBioPortal diagnosis-based OS interval, so identical time origins are not asserted. OCCAMS Alive/Not recorded dead is censored as encoded by its source. BayesPrism fractions are RNA contributions, not validated cell proportions, and reference uncertainty is not propagated into Cox standard errors. Many inferred fractions are near zero, with high relative conditional posterior uncertainty; numerical model estimates do not establish accurately measured state abundance. Normal epithelial contamination and similar malignant-state programs remain limitations. Survival associations are prognostic observations, not causal validation.', '',
  '[Final MP figure](../figures/final/01_malignant_mp_survival.pdf) | [Final state figure](../figures/final/02_malignant_state_survival.pdf)',
  '[Figure index](../tables/figure_index.csv) | [Cox results](../tables/cox_results.csv) | [Cross-cohort concordance](../tables/concordance_results.csv)',
  'Exploratory RFS reproduces the earlier selected-cohort endpoint: OCCAMS recurrence/death under the previously approved common-first-therapy assumption; TCGA new tumour event/death from diagnosis, including new primary. These are not harmonized PFS endpoints. Missing or inconsistent recurrence records are excluded explicitly. Supporting 02 places OS and RFS side by side on every method/comparison page. Supporting 02/05/06 use open circles for NS, solid circles for nominal P<0.05 and diamonds for BH q<0.05. Panel headers report cohort n and events; ranges denote feature-dependent counts.', '',
  '[Endpoint availability](../tables/endpoint_availability.csv) | [RFS results](../tables/rfs_cox_results.csv) | [RFS definitions and exclusions](../tables/rfs_patient_metadata.csv)',
  '[OS adjustment](../figures/supporting/05_continuous_os_adjustment.pdf) | [RFS adjustment](../figures/supporting/06_continuous_rfs_adjustment.pdf)',
  '[OCCAMS subgroup survival](../figures/supporting/02_neoadjuvant_context.pdf) | [Subgroup page index](../tables/occams_context_figure_index.csv)',
  '[Cohort counts](../tables/cohort_summary.csv) | [Metadata source audit](../tables/metadata_source_audit.csv) | [Validation](validation.csv)',
  '[State-fraction resolution](../../bulk_malignant_deconvolution/tables/state_fraction_resolution.csv)',
  '[Methodology](../../../analysis/methodology/clinical/malignant_bulk_survival_methodology.md)')
writeLines(report, file.path(out, 'reports/current_results.md'))
####################
writeLines(c('# Current malignant survival analysis','',
  'Start with the [current results and figure index](reports/current_results.md).',
  'The selected final figures are in `figures/final/`; complete model and source tables are in `tables/`.',
  'Methods: [malignant-compartment survival](../../analysis/methodology/clinical/malignant_bulk_survival_methodology.md).'),
  file.path(out,'README.md'))
####################
message('Complete: ', out)
