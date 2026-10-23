####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/survival_selected_cohorts.R
#   Description: Selected-cohort OS/exploratory RFS stage KM panels, final-17-MP
#     and dual-state-signature OS volcanoes, and neoadjuvant specimen boxplots,
#     with stable current outputs, archived revisions, and reproducible caches.
#   Methodology: analysis/methodology/clinical/survival_selected_cohorts_methodology.md
#   Map: analysis/ANALYSIS_MAP.md.
# Exact inputs, relative to the live project root:
#   analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx
#   analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx
#   ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds
#   ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#   ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds
#   ref_outs/clinical_selected_survival/source_data/TCGA_BCR_endpoint_candidates.csv
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
# Outputs under ref_outs/clinical_selected_survival/:
#   intermediate/<signature>_<run_id>.rds: reproducible patient/score/model cache.
#   intermediate/cache_pointer.txt
#   tables/: patient data, cohort/endpoint eligibility, GSVA coverage (all 27 signatures),
#     MP scores/descriptions, Cox estimates, optimal-cut grid, KM source/risk tables,
#     panel inventory, source hashes, and package versions.
#   figures/selected_cohort_stage_survival.pdf (7 pages)
#   figures/selected_cohort_mp_os_volcanoes.pdf (6 pages)
#   reports/run_summary.txt; PBS stdout/stderr directly in project root.
#   updates/new_updates/summaries/selected_cohort_survival_<run_id>.csv (<100 KB).
# Cache/replot: SCREF_FORCE_REBUILD=TRUE builds a new immutable cache;
#   SCREF_REPLOT_ONLY=TRUE requires a matching cache and performs plotting only.
#   SCREF_SELECTED_SURVIVAL_RUN_ID labels the revision; R archives prior output tiers.
#   SCREF_SELECTED_SURVIVAL_THREADS controls GSVA workers (default 4).
#   OCCAMS RFS uses the user's explicitly approved common-first-therapy assumption.
# Run on an allocated compute node: Rscript analysis/clinical/survival_selected_cohorts.R
# Optional PBS launcher: /opt/pbs/bin/qsub job_scripts/clinical/survival_selected_cohorts.sh
# Resources: 4 CPUs, 24 GB RAM, 02:00:00; scheduler stdout/stderr in project root.
# Additional input: ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv
# Additional outputs: figures/occams_neoadjuvant_biopsy_vs_surgery_boxplots.pdf
#   (2 pages: MPs, then both state signature families), and
#   figures/selected_cohort_mp_and_state_os_volcanoes.pdf (6 pages).
#   tables/state_signature_genes.csv, all_feature_gsva_scores.csv,
#   all_feature_os_cox_results.csv, all_feature_os_optimal_cut_grid.csv,
#   all_feature_volcano_label_layout.csv, specimen_boxplot_source_data.csv,
#   specimen_boxplot_statistics.csv. State scores are bulk GSVA, not proportions.
# R owns output archiving, caching, analysis and plotting; no Python/audit dependency.
# Conda environment: dmtcp; installed GSVA 1.46.0 API, no alternate analysis inputs.
####################

library(readxl)
library(GSVA)
library(survival)
library(ggplot2)
library(gridExtra)
library(grid)

project_dir <- '/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline'
source(file.path(project_dir, 'analysis/shared/scRef_config.R'))
setwd(project_dir)
out_root <- file.path(project_dir, 'ref_outs/clinical_selected_survival')
run_id <- Sys.getenv('SCREF_SELECTED_SURVIVAL_RUN_ID', paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid()))
if (!grepl('^[A-Za-z0-9_.-]+$', run_id)) stop('Invalid run ID')
run_dir <- out_root
dirs <- setNames(file.path(run_dir, c('intermediate','tables','figures','reports','logs')), c('intermediate','tables','figures','reports','logs'))
invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))
cache_dir <- file.path(out_root, 'intermediate')
dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
force_rebuild <- identical(toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE')), 'TRUE')
replot_only <- identical(toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE')), 'TRUE')
if (force_rebuild && replot_only) stop('FORCE_REBUILD and REPLOT_ONLY cannot both be TRUE')
threads <- as.integer(Sys.getenv('SCREF_SELECTED_SURVIVAL_THREADS', '4'))
if (!is.finite(threads) || threads < 1) stop('Invalid thread count')
input_paths <- c(
  occams_workbook='analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx',
  tcga_workbook='analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx',
  qc='ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv',
  occams_expr='ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds',
  tcga_meta='ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds',
  tcga_expr='ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds',
  tcga_events='ref_outs/clinical_selected_survival/source_data/TCGA_BCR_endpoint_candidates.csv',
  mp_genes=SCREF_MP_GENES_RDS, mp_grouping=SCREF_MP_GROUPING_CSV,
  state_markers='ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv'
)
if (any(!file.exists(input_paths))) stop('Missing required input: ', paste(input_paths[!file.exists(input_paths)], collapse=', '))
input_hashes <- data.frame(input=names(input_paths), path=normalizePath(input_paths), md5=unname(tools::md5sum(input_paths)), stringsAsFactors=FALSE)
package_versions <- data.frame(package=c('R','GSVA','survival','ggplot2','readxl','gridExtra'),
  version=c(as.character(getRversion()), vapply(c('GSVA','survival','ggplot2','readxl','gridExtra'), function(p) as.character(packageVersion(p)), character(1))))
settings <- list(version='selected_survival_v4_two_state_signatures', occams_neo_typed=157L, occams_surgery=68L, tcga_all=86L, tcga_yes=72L,
  occams_rfs='exploratory common first-recorded-therapy origin explicitly approved by user',
  tcga_rfs='original new tumour event or death; includes new primary to reproduce earlier cohort',
  tcga_expression='raw TPM, Gaussian GSVA, matching survival_integrated_analysis.R',
  gsva_parent_cohorts=c('OCCAMS typed neoadjuvant 157','OCCAMS matched surgery pathway 68','TCGA surgery pathway 86'),
  cut_probabilities=seq(.20,.80,.05), min_model_n=20L, min_events=5L, min_binary_group=5L)
####################
# R owns revision preservation so direct Rscript runs need no shell/Python setup.
archive_dir <- file.path(out_root,'archive',paste0('before_',run_id))
if (dir.exists(archive_dir)) stop('Revision archive exists; choose a new run ID: ',run_id)
archive_output <- function(path,destination) {
  if (file.exists(destination)) stop('Archive destination exists: ',destination)
  dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
  if (!file.rename(path,destination)) stop('Cannot archive output: ',path)
}
for (tier in c('figures','tables','reports','validation','intermediate')) {
  path <- file.path(out_root,tier)
  if (dir.exists(path) && length(list.files(path,all.files=TRUE,no..=TRUE)))
    archive_output(path,file.path(archive_dir,tier))
  if (tier!='validation') dir.create(path,recursive=TRUE,showWarnings=FALSE)
}
job_id <- sub('\\..*$','',Sys.getenv('PBS_JOBID',''))
old_logs <- list.files(dirs[['logs']],full.names=TRUE)
if (nzchar(job_id)) old_logs <- old_logs[!grepl(job_id,basename(old_logs),fixed=TRUE)]
for (path in old_logs) archive_output(path,file.path(archive_dir,'logs',basename(path)))
old_summaries <- list.files(file.path(project_dir,'updates/new_updates/summaries'),
  pattern='^selected_cohort_survival_.*\\.csv$',full.names=TRUE)
for (path in old_summaries) archive_output(path,file.path(archive_dir,'compact_summaries',basename(path)))
####################
signature_file <- file.path(dirs[['intermediate']], 'signature_input.rds')
saveRDS(list(settings=settings, inputs=input_hashes, versions=package_versions), signature_file)
signature <- unname(tools::md5sum(signature_file))
cache_matches <- sort(list.files(out_root, pattern=paste0('^',signature,'_.*\\.rds$'), full.names=TRUE, recursive=TRUE))
write_table <- function(x, name) {
  path <- file.path(dirs[['tables']], name)
  if (file.exists(path)) stop('Refusing to overwrite: ', path)
  write.csv(x, path, row.names=FALSE, na='')
}
stage_value <- function(x) {
  z <- trimws(tolower(as.character(x)))
  z[is.na(z) | !z %in% as.character(1:4)] <- 'unknown'
  z
}
blank <- function(x) is.na(x) | trimws(as.character(x)) == ''
positive <- function(x) is.finite(x) & x > 0
excluded_occams <- c(
  '26d51bb9b7c9974b5e753fc3520ffcfa2abdbd132509a887f4a5bc5630aa05d9',
  '72089d94e1a29cc0b7df4fb8f9694343cda2c6e71ae1083fc266b5b36b62323d',
  '9ea99e16eaaf486784cb0cd7e1c846a4dcc7bbaecfccc9a04d7853093f2cf775')

build_patient_data <- function() {
  o <- as.data.frame(read_excel(input_paths[['occams_workbook']], sheet='Patient_Master_388'))
  t <- as.data.frame(read_excel(input_paths[['tcga_workbook']], sheet='EAC_RNA88_Final'))
  qc <- read.csv(input_paths[['qc']], stringsAsFactors=FALSE)
  stopifnot(nrow(o)==388L, !anyDuplicated(o$SHA_ID), nrow(t)==88L, !anyDuplicated(t$Patient_ID))
  stopifnot(all(excluded_occams %in% o$SHA_ID), !any(blank(o$Clinical_Stage_AJCC)), !any(blank(o$Pathological_Stage_AJCC)),
            !any(blank(t$Clinical_Stage_AJCC8)), !any(blank(t$Pathological_Stage_AJCC8)))
  qc_ids <- qc$Subject_ID[qc$Annotation_Compatible_Below_50_Flag == 'no']
  stopifnot(length(qc_ids)==281L, all(excluded_occams %in% qc_ids))
  keep <- o$SHA_ID %in% qc_ids & !o$SHA_ID %in% excluded_occams
  o <- o[keep, , drop=FALSE]
  stopifnot(nrow(o)==278L, all(o$Vital_Status %in% c('Dead','Alive/Not recorded dead')))
  raw_neo <- !is.na(o$Curative_Treatment_Modality_Final) & o$Treatment_Intent_Final == 'curative' & o$Curative_Treatment_Modality_Final == 'neoadjuvant (chemo + CRT)'
  raw_surgery <- !is.na(o$Curative_Treatment_Modality_Final) & o$Treatment_Intent_Final == 'curative' & o$Curative_Treatment_Modality_Final == 'surgery'
  typed_neo <- raw_neo & o$RNA_specimen_source %in% c('biopsy','surgery') & o$RNA_Metadata_Matched == 'Yes'
  matched_surgery <- raw_surgery & o$RNA_Metadata_Matched == 'Yes'
  stopifnot(sum(matched_surgery)==68L)
  stopifnot(sum(raw_neo)==165L, sum(raw_surgery)==72L, sum(typed_neo)==157L,
    sum(typed_neo & o$RNA_specimen_source == 'biopsy', na.rm=TRUE)==77L,
    sum(typed_neo & o$RNA_specimen_source == 'surgery', na.rm=TRUE)==80L)
  os_time <- as.numeric(o$Survival_Days)
  os_event <- as.integer(o$Vital_Status == 'Dead')
  rec_time <- as.numeric(o$Estimated_Recurrence_Time_Days_From_Diagnosis_Surrogate)
  rfs_time <- rfs_event <- rep(NA_real_, nrow(o))
  rfs_reason <- rep('', nrow(o))
  for (i in seq_len(nrow(o))) {
    rec_yes <- !is.na(o$Recurrence_Ever[i]) && o$Recurrence_Ever[i] == 'Yes'
    rec_no <- !is.na(o$Recurrence_Ever[i]) && o$Recurrence_Ever[i] == 'No'
    if (!positive(os_time[i])) rfs_reason[i] <- 'missing/nonpositive survival interval'
    else if (rec_yes && !positive(rec_time[i])) rfs_reason[i] <- 'recorded recurrence without positive timed interval'
    else if (positive(rec_time[i]) && rec_time[i] > os_time[i]) rfs_reason[i] <- 'recurrence interval exceeds survival/follow-up interval'
    else if (!rec_yes && !rec_no && os_event[i] == 0L) rfs_reason[i] <- 'living without recurrence ascertainment'
    else if (rec_yes) {
      rfs_time[i] <- if (os_event[i]==1L) min(rec_time[i],os_time[i]) else rec_time[i]
      rfs_event[i] <- 1L
    } else {
      rfs_time[i] <- os_time[i]
      rfs_event[i] <- os_event[i]
    }
  }
  om <- data.frame(subject_id=o$SHA_ID, sample_id=o$SHA_ID, dataset='OCCAMS',
    cstage=stage_value(o$Clinical_Stage_AJCC), pstage=stage_value(o$Pathological_Stage_AJCC),
    source=as.character(o$RNA_specimen_source), pathway=as.character(o$Curative_Treatment_Modality_Final),
    surgery_status='', neo_selected=typed_neo, surgery_selected=matched_surgery,
    OS_time=os_time, OS_event=os_event, RFS_time=rfs_time, RFS_event=rfs_event,
    OS_exclusion=ifelse(positive(os_time),'','missing/nonpositive survival interval'), RFS_exclusion=rfs_reason,
    RFS_definition='Exploratory recurrence/death; assumed common first-therapy origin', stringsAsFactors=FALSE)
  t <- t[t$Flowchart_Pathway_Final == 'Surgery group', , drop=FALSE]
  stopifnot(nrow(t)==86L, sum(t$Flowchart_Surgery_Status_Final == 'Yes')==72L,
    sum(t$Flowchart_Surgery_Status_Final == 'No')==12L, sum(t$Flowchart_Surgery_Status_Final == 'Unresolved')==2L)
  tm <- readRDS(input_paths[['tcga_meta']])
  tm <- tm[tm$sample_type_code == '01' & tm$HistologyGroup == 'EAC', , drop=FALSE]
  j <- match(t$RNA_primary_sample_ID, tm$sample_barcode)
  stopifnot(!anyNA(j), all(tm$case_barcode[j] == t$Patient_ID))
  events <- read.csv(input_paths[['tcga_events']], stringsAsFactors=FALSE, na.strings='')
  e <- events[match(t$Patient_ID, events$Patient_ID), , drop=FALSE]
  stopifnot(!anyNA(e$Patient_ID))
  has_indicator <- !blank(e$tumour_event_indicator)
  has_event <- !blank(e$tumour_event_days)
  has_death <- !blank(e$death_days)
  rfs_reason <- ifelse(!blank(e$exclusion_reason), e$exclusion_reason,
    ifelse(!has_indicator & !has_event & !has_death, 'living without tumour-event ascertainment', ''))
  rfs_time <- as.numeric(e$candidate_rfs_days)
  rfs_event <- as.numeric(e$candidate_rfs_event)
  rfs_time[rfs_reason != ''] <- NA_real_
  rfs_event[rfs_reason != ''] <- NA_real_
  tmeta <- data.frame(subject_id=t$Patient_ID, sample_id=t$RNA_primary_sample_ID, dataset='TCGA',
    cstage=stage_value(t$Clinical_Stage_AJCC8), pstage=stage_value(t$Pathological_Stage_AJCC8),
    source='primary tumour', pathway='Surgery group', surgery_status=t$Flowchart_Surgery_Status_Final,
    neo_selected=FALSE, surgery_selected=FALSE,
    OS_time=as.numeric(tm$OS_time[j]), OS_event=as.numeric(tm$OS_event[j]), RFS_time=rfs_time, RFS_event=rfs_event,
    OS_exclusion=ifelse(positive(tm$OS_time[j]),'','missing/nonpositive survival interval'), RFS_exclusion=rfs_reason,
    RFS_definition='Original new tumour event or death; includes one new primary; diagnosis origin', stringsAsFactors=FALSE)
  all <- rbind(om,tmeta)
  stopifnot(!anyDuplicated(all$subject_id), !any(all$subject_id %in% excluded_occams))
  all
}

cohort_indices <- function(meta) {
  is_o <- meta$dataset == 'OCCAMS'
  is_t <- meta$dataset == 'TCGA'
  list(occams_combined=which(is_o & (meta$neo_selected | meta$surgery_selected)),
    occams_neo=which(is_o & meta$neo_selected),
    occams_neo_biopsy=which(is_o & meta$neo_selected & meta$source=='biopsy'),
    occams_neo_surgery=which(is_o & meta$neo_selected & meta$source=='surgery'),
    occams_surgery=which(is_o & meta$surgery_selected),
    tcga_all=which(is_t), tcga_yes=which(is_t & meta$surgery_status=='Yes'),
    tcga_no=which(is_t & meta$surgery_status=='No'), tcga_unresolved=which(is_t & meta$surgery_status=='Unresolved'))
}
cohort_labels <- c(occams_combined='OCCAMS neoadjuvant + surgery pathways', occams_neo='OCCAMS neoadjuvant: biopsy + surgery',
  occams_neo_biopsy='OCCAMS neoadjuvant: biopsy specimens', occams_neo_surgery='OCCAMS neoadjuvant: surgery specimens',
  occams_surgery='OCCAMS surgery pathway', tcga_all='TCGA surgery pathway: all patients', tcga_yes='TCGA surgery pathway: surgery Yes',
  tcga_no='TCGA surgery pathway: surgery No', tcga_unresolved='TCGA surgery pathway: surgery unresolved')
mp_cohorts <- c('occams_neo','occams_neo_biopsy','occams_neo_surgery','occams_surgery','tcga_all','tcga_yes')
split_methods <- c('continuous','median','q1q4','optimal_cut')

score_parent <- function(expr, sample_ids, mp_genes, dataset) {
  stopifnot(is.matrix(expr), !anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)), all(sample_ids %in% colnames(expr)))
  expr <- expr[, sample_ids, drop=FALSE]
  if (any(!is.finite(expr))) stop('Nonfinite expression in ',dataset)
  sets <- lapply(mp_genes,function(g) intersect(unique(g),rownames(expr)))
  if (any(lengths(sets)<5L)) stop('Insufficient signature gene coverage in ',dataset)
  message('Computing selected-parent Gaussian GSVA: ',dataset,' n=',ncol(expr))
  score <- GSVA::gsva(expr,sets,method='gsva',kcdf='Gaussian',verbose=FALSE,parallel.sz=threads)
  if (!identical(rownames(score),names(mp_genes)) || !setequal(colnames(score),sample_ids) || any(!is.finite(score))) stop('GSVA output validation failed')
  coverage <- data.frame(dataset=dataset,mp=names(mp_genes),source_genes=lengths(mp_genes),matched_genes=lengths(sets),stringsAsFactors=FALSE)
  list(scores=t(score[,sample_ids,drop=FALSE]),coverage=coverage)
}

cox_one <- function(d, feature, method, prob=NA_real_) {
  if (!feature %in% names(d)) stop('Missing exact signature score column: ',feature)
  d <- d[positive(d$OS_time) & d$OS_event %in% c(0,1) & is.finite(d[[feature]]), , drop=FALSE]
  result <- data.frame(feature=feature,method=method,status='eligible',n=nrow(d),events=sum(d$OS_event),
    low_n=NA_integer_,high_n=NA_integer_,cutpoint=NA_real_,q1=NA_real_,q3=NA_real_,quantile_prob=prob,
    logHR=NA_real_,SE=NA_real_,HR=NA_real_,lower_CI=NA_real_,upper_CI=NA_real_,P_value=NA_real_,stringsAsFactors=FALSE)
  fail <- function(reason) {result$status <- reason;result}
  if (nrow(d)<settings$min_model_n) return(fail('ineligible_n_below_20'))
  if (sum(d$OS_event)<settings$min_events) return(fail('ineligible_fewer_than_5_deaths'))
  if (!is.finite(sd(d[[feature]])) || sd(d[[feature]])==0) return(fail('ineligible_zero_variance'))
  if (method=='continuous') d$split_val <- as.numeric(scale(d[[feature]]))
  else {
    if (method=='median') {
      cut <- median(d[[feature]])
      group <- ifelse(d[[feature]]>cut,'High','Low')
      result$cutpoint <- cut
    } else if (method=='q1q4') {
      q <- quantile(d[[feature]],c(.25,.75),names=FALSE)
      result$q1 <- q[1];result$q3 <- q[2]
      if (q[1]>=q[2]) return(fail('ineligible_tied_quartile_thresholds'))
      d <- d[d[[feature]]<=q[1] | d[[feature]]>=q[2], , drop=FALSE]
      group <- ifelse(d[[feature]]>=q[2],'High','Low')
    } else if (method=='optimal_cut') {
      cut <- quantile(d[[feature]],prob,names=FALSE)
      result$cutpoint <- cut
      group <- ifelse(d[[feature]]>=cut,'High','Low')
    } else stop('Unknown method')
    d$split_val <- factor(group,levels=c('Low','High'))
    result$n <- nrow(d);result$events <- sum(d$OS_event)
    result$low_n <- sum(group=='Low');result$high_n <- sum(group=='High')
    if (nrow(d)<settings$min_model_n) return(fail('ineligible_n_below_20'))
    if (sum(d$OS_event)<settings$min_events) return(fail('ineligible_fewer_than_5_deaths'))
    if (min(table(d$split_val))<settings$min_binary_group) return(fail('ineligible_group_below_5'))
    if (any(tapply(d$OS_event,d$split_val,sum)==0)) return(fail('ineligible_zero_deaths_in_group'))
  }
  fit <- withCallingHandlers(coxph(Surv(OS_time,OS_event) ~ split_val,data=d,control=coxph.control(iter.max=50)),
    warning=function(w) stop('Cox fit warning for ',feature,' / ',method,': ',conditionMessage(w)))
  ss <- summary(fit)
  values <- c(ss$coefficients[1,'coef'],ss$coefficients[1,'se(coef)'],ss$coefficients[1,'exp(coef)'],ss$conf.int[1,'lower .95'],ss$conf.int[1,'upper .95'],ss$coefficients[1,'Pr(>|z|)'])
  if (any(!is.finite(values))) stop('Nonfinite fitted Cox estimate: ',feature,' / ',method)
  result[1,c('logHR','SE','HR','lower_CI','upper_CI','P_value')] <- values
  result$status <- 'fitted'
  result
}

if (length(cache_matches)>0L && !force_rebuild) {
  cache_path <- tail(cache_matches,1)
  cache <- readRDS(cache_path)
  stopifnot(identical(cache$signature,signature), identical(cache$settings,settings))
  canonical_cache <- file.path(cache_dir,basename(cache_path))
  if (normalizePath(cache_path) != normalizePath(canonical_cache,mustWork=FALSE)) {
    stopifnot(file.copy(cache_path,canonical_cache,overwrite=FALSE))
    cache_path <- canonical_cache
  }
  message('Using validated analysis cache: ',cache_path)
} else {
  if (replot_only) stop('REPLOT_ONLY requested but no matching complete cache exists')
  meta <- build_patient_data()
  cohorts <- cohort_indices(meta)
  stopifnot(identical(as.integer(lengths(cohorts)),c(225L,157L,77L,80L,68L,86L,72L,12L,2L)))
  grouping <- read.csv(input_paths[['mp_grouping']],stringsAsFactors=FALSE)
  mp_genes <- readRDS(input_paths[['mp_genes']])[grouping$mp]
  stopifnot(length(mp_genes)==17L, !any(vapply(mp_genes,is.null,logical(1))), !any(c('MP2x','MP11c','MP18a') %in% names(mp_genes)))
####################
  # Match the integrated-survival definitions exactly; DGE rows retain source order.
  state_groups <- list('Classic proliferation'=c('MP2+'),
    'Squamous-to-intestinal'=c('MP14','MP3+','MP6+','MP11+','MP9+','MP10+'),
    'Glandular-to-intestinal'=c('MP18b','MP16','MP17','MP8b','MP8+'),
    'Stress-adaptive'='MP12', 'Cancer-cell immune mimicry'='MP15')
  stopifnot(all(unlist(state_groups) %in% names(mp_genes)))
  state_genes <- lapply(state_groups,function(mps) unique(unlist(mp_genes[mps],use.names=FALSE)))
  names(state_genes) <- paste0(names(state_groups),' (MP union)')
  markers <- read.csv(input_paths[['state_markers']],stringsAsFactors=FALSE)
  stopifnot(all(names(state_groups) %in% markers$state))
  dge_genes <- lapply(names(state_groups),function(s) head(markers$gene[markers$state==s],20L))
  names(dge_genes) <- paste0(names(state_groups),' (DGE)')
  stopifnot(all(lengths(dge_genes)==20L),all(lengths(state_genes)>=5L))
  feature_sets <- c(mp_genes,state_genes,dge_genes)
  feature_family <- setNames(c(rep('MP',17),rep('State MP union',5),rep('State DGE',5)),names(feature_sets))
  # Separate calls mirror the integrated script's three signature-family GSVA calls.
  score_all <- function(expr,ids,dataset) {
    parts <- lapply(list(mp_genes,state_genes,dge_genes),function(sets) score_parent(expr,ids,sets,dataset))
    list(scores=do.call(cbind,lapply(parts,`[[`,'scores')),
      coverage=do.call(rbind,lapply(parts,`[[`,'coverage')))
  }
####################
  oe <- readRDS(input_paths[['occams_expr']])
  te <- readRDS(input_paths[['tcga_expr']])
  so <- score_all(oe,meta$sample_id[cohorts$occams_neo],'OCCAMS neoadjuvant typed 157')
  st <- score_all(te,meta$sample_id[cohorts$tcga_all],'TCGA surgery pathway 86')
####################
  # Score the newly requested surgery pathway on its own selected parent cohort.
  # Neoadjuvant and TCGA scoring populations remain unchanged.
  so_surgery <- score_all(oe,meta$sample_id[cohorts$occams_surgery],'OCCAMS matched surgery pathway 68')
  scores <- rbind(so$scores,so_surgery$scores,st$scores)
####################
  results <- grids <- list()
  for (cohort in mp_cohorts) {
    d <- meta[cohorts[[cohort]], , drop=FALSE]
    j <- match(d$sample_id,rownames(scores));stopifnot(!anyNA(j))
    d <- cbind(d,as.data.frame(scores[j,,drop=FALSE],optional=TRUE))
    for (feature in names(feature_sets)) {
      for (method in split_methods) {
        if (method=='optimal_cut') {
          candidates <- do.call(rbind,lapply(settings$cut_probabilities,function(p) cox_one(d,feature,method,p)))
          candidates$cohort <- cohort
          grids[[paste(cohort,feature)]] <- candidates
          ok <- candidates$status=='fitted'
          chosen <- if (any(ok)) candidates[which(ok)[which.min(candidates$P_value[ok])],,drop=FALSE] else candidates[1,,drop=FALSE]
          if (!any(ok)) chosen$status <- 'ineligible_no_estimable_optimal_cut'
          r <- chosen
        } else r <- cox_one(d,feature,method)
        r$cohort <- cohort
        results[[paste(cohort,feature,method)]] <- r
      }
    }
  }
  cox_results <- do.call(rbind,results);rownames(cox_results)<-NULL
  cox_results$feature_family <- unname(feature_family[cox_results$feature])
  cox_results$FDR_BH_all_features <- NA_real_
  cox_results$FDR_BH <- NA_real_
  for (key in unique(paste(cox_results$cohort,cox_results$method))) {
    ix <- paste(cox_results$cohort,cox_results$method)==key & cox_results$status=='fitted'
    cox_results$FDR_BH_all_features[ix] <- p.adjust(cox_results$P_value[ix],method='BH')
    for (family in unique(feature_family)) {
      j <- ix & cox_results$feature_family==family
      cox_results$FDR_BH[j] <- p.adjust(cox_results$P_value[j],method='BH')
    }
  }
  cox_results$endpoint <- 'OS'
  cox_results$covariates <- 'none'
  cox_results$optimal_cut_inference <- ifelse(cox_results$method=='optimal_cut','exploratory minimum raw Wald p; not cut-search adjusted','prespecified split')
  cache <- list(signature=signature,settings=settings,meta=meta,cohorts=cohorts,mp_genes=mp_genes,scores=scores,
    feature_sets=feature_sets,feature_family=feature_family,
    coverage=rbind(so$coverage,so_surgery$coverage,st$coverage),cox_results=cox_results,optimal_grid=do.call(rbind,grids))
  cache_path <- file.path(cache_dir,paste0(signature,'_',run_id,'.rds'))
  if (file.exists(cache_path)) stop('Refusing to overwrite a cache')
  saveRDS(cache,cache_path,compress='xz')
}
meta <- cache$meta;cohorts <- cache$cohorts;cox_results <- cache$cox_results
writeLines(cache_path,file.path(dirs[['intermediate']],'cache_pointer.txt'))
write_table(input_hashes,'input_manifest.csv');write_table(package_versions,'package_versions.csv')
write_table(meta,'patient_endpoint_data.csv')
write_table(data.frame(sample_id=rownames(cache$scores),cache$scores[,names(cache$mp_genes)],check.names=FALSE),'mp_gsva_scores.csv')
write_table(data.frame(mp=names(cache$mp_genes),description=unname(SCREF_MP_DESCRIPTIONS[names(cache$mp_genes)])),'mp_descriptions.csv')
write_table(cache$coverage,'gsva_gene_coverage.csv');write_table(cox_results[cox_results$feature_family=='MP',],'mp_os_cox_results.csv')
write_table(cache$optimal_grid[cache$optimal_grid$feature %in% names(cache$mp_genes),],'mp_os_optimal_cut_grid.csv')
####################
write_table(data.frame(sample_id=rownames(cache$scores),cache$scores,check.names=FALSE),'all_feature_gsva_scores.csv')
write_table(cox_results,'all_feature_os_cox_results.csv')
write_table(cache$optimal_grid,'all_feature_os_optimal_cut_grid.csv')
write_table(do.call(rbind,lapply(names(cache$feature_sets)[cache$feature_family!='MP'],function(f)
  data.frame(feature=f,family=cache$feature_family[[f]],gene=cache$feature_sets[[f]],rank=seq_along(cache$feature_sets[[f]])))),
  'state_signature_genes.csv')
####################
members <- do.call(rbind,lapply(names(cohorts),function(key) data.frame(cohort=key,subject_id=meta$subject_id[cohorts[[key]]])))
write_table(members,'cohort_membership.csv')
cohort_summary <- do.call(rbind,lapply(names(cohorts),function(key) {
  d<-meta[cohorts[[key]],,drop=FALSE]
  data.frame(cohort=key,n=nrow(d),OS_n=sum(positive(d$OS_time)),OS_events=sum(d$OS_event[positive(d$OS_time)]),
    RFS_n=sum(positive(d$RFS_time)),RFS_events=sum(d$RFS_event[positive(d$RFS_time)]),
    clinical_known=sum(d$cstage %in% as.character(1:4)),pathological_known=sum(d$pstage %in% as.character(1:4)))
}))
write_table(cohort_summary,'cohort_endpoint_summary.csv')

stage_colours <- c('1'='#0072B2','2'='#009E73','3'='#D55E00','4'='#8E44AD')
km_sources <- km_risks <- km_panels <- list()
panel_number <- 0L
theme_readable <- theme_bw(base_size=19) + theme(plot.title=element_text(size=23,face='bold'),
  plot.subtitle=element_text(size=17,lineheight=1.05),axis.title=element_text(size=20),axis.text=element_text(size=17),
  legend.title=element_blank(),legend.text=element_text(size=17),legend.position='top',
  legend.key.width=unit(.7,'cm'),legend.key.spacing.x=unit(.5,'cm'),
  plot.margin=margin(10,18,8,12),panel.grid.minor=element_blank())

km_panel <- function(d,endpoint,stage_column,cohort_key,page) {
  panel_number <<- panel_number+1L
  time <- d[[paste0(endpoint,'_time')]];event <- d[[paste0(endpoint,'_event')]]
  known <- d[[stage_column]] %in% as.character(1:4)
  endpoint_ok <- positive(time) & event %in% c(0,1)
  keep <- known & endpoint_ok
  z <- data.frame(subject_id=d$subject_id[keep],time=time[keep],event=event[keep],stage=factor(d[[stage_column]][keep],levels=1:4))
  z$stage <- droplevels(z$stage)
  stage_label <- if (stage_column=='cstage') 'Clinical stage' else 'Pathological stage'
  endpoint_label <- endpoint
  title <- paste(endpoint_label,stage_label,sep=' | ')
  p_value <- NA_real_
  if (nrow(z)>=20L && sum(z$event)>=5L && nlevels(z$stage)>=2L) {
    test <- survdiff(Surv(time,event) ~ stage,data=z)
    p_value <- pchisq(test$chisq,length(test$n)-1L,lower.tail=FALSE)
  }
  km_panels[[as.character(panel_number)]] <<- data.frame(pdf_page=page,panel=panel_number,cohort=cohort_key,endpoint=endpoint,stage=stage_column,
    cohort_n=nrow(d),evaluable_n=nrow(z),events=sum(z$event),unknown_stage=sum(!known),endpoint_missing=sum(!endpoint_ok),
    logrank_p=p_value,low_coverage=nrow(z)<20L)
  if (nrow(z)==0L) return(ggplotGrob(ggplot()+theme_void()+annotate('text',x=0,y=0,label=paste(title,'No evaluable patients',sep='\n'),size=9)))
  fit <- survfit(Surv(time/365.25,event) ~ stage,data=z)
  ss <- summary(fit,censored=TRUE)
  curve_stage <- if (nlevels(z$stage)==1L) rep(levels(z$stage),length(ss$time)) else sub('^stage=','',as.character(ss$strata))
  curve <- data.frame(years=ss$time,survival=ss$surv,n_event=ss$n.event,n_censor=ss$n.censor,stage=curve_stage)
  start <- data.frame(years=0,survival=1,n_event=0,n_censor=0,stage=levels(z$stage))
  curve <- rbind(start,curve);curve <- curve[order(curve$stage,curve$years),]
  curve$cohort <- cohort_key;curve$endpoint <- endpoint;curve$stage_column <- stage_column;curve$pdf_page <- page
  km_sources[[as.character(panel_number)]] <<- curve
  x_max <- max(10,ceiling(max(c(d$OS_time,d$RFS_time),na.rm=TRUE)/365.25/2)*2)
  breaks <- seq(0,x_max,2)
  rs <- summary(fit,times=breaks,extend=TRUE)
  risk_stage <- if (nlevels(z$stage)==1L) rep(levels(z$stage),length(rs$time)) else sub('^stage=','',as.character(rs$strata))
  risk <- data.frame(years=rs$time,n_risk=rs$n.risk,stage=risk_stage,cohort=cohort_key,endpoint=endpoint,stage_column=stage_column,pdf_page=page)
  km_risks[[as.character(panel_number)]] <<- risk
  n_stage <- table(z$stage)
  legend_labels <- setNames(paste0('Stage ',names(n_stage),' (N=',as.integer(n_stage),')'),names(n_stage))
  if (is.finite(p_value)) {
    p_text <- format.pval(p_value,digits=3,eps=1e-16)
    p_label <- paste0('Global log-rank p ',if (startsWith(p_text,'<')) '' else '= ',p_text)
  } else p_label <- 'Log-rank not reported: sparse data'
  subtitle <- paste0('n=',nrow(z),' | Events = ',sum(z$event),' | unknown removed = ',sum(!known))
  x_label <- 'Time (years)'
  y_label <- if (endpoint=='OS') 'Overall survival probability' else paste(endpoint_label,'probability')
  p <- ggplot(curve,aes(years,survival,colour=stage,group=stage))+geom_step(linewidth=1.25)+
    geom_point(data=curve[curve$n_censor>0,,drop=FALSE],shape=3,size=2.5,stroke=.8)+
    scale_colour_manual(values=stage_colours,breaks=names(n_stage),labels=legend_labels)+
    scale_x_continuous(breaks=breaks,limits=c(0,x_max),expand=expansion(mult=c(.015,.025)))+
    scale_y_continuous(limits=c(0,1.02),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.015)))+
    labs(title=title,subtitle=subtitle,x=x_label,y=y_label)+theme_readable+
    guides(colour=guide_legend(nrow=1,byrow=TRUE))+
    annotate('text',x=x_max*.98,y=.96,label=p_label,hjust=1,vjust=1,size=5.3,colour='black')
  risk$stage <- factor(risk$stage,levels=rev(names(n_stage)))
  pr <- ggplot(risk,aes(years,stage,label=n_risk,colour=as.character(stage)))+
    geom_text(size=6.2,show.legend=FALSE)+scale_colour_manual(values=stage_colours)+
    scale_x_continuous(breaks=breaks,limits=c(0,x_max),expand=expansion(mult=c(.015,.025)))+
    scale_y_discrete(labels=function(x) paste('Stage',x))+labs(title='Number at risk',x=NULL,y=NULL)+
    theme_classic(base_size=17)+theme(plot.title=element_text(size=18,face='bold'),axis.text=element_text(size=17),
      axis.text.y=element_text(margin=margin(r=18)),
      axis.line=element_blank(),axis.ticks=element_blank(),plot.margin=margin(0,18,5,12))
  g1 <- ggplotGrob(p);g2 <- ggplotGrob(pr)
  widths <- unit.pmax(g1$widths,g2$widths);g1$widths<-widths;g2$widths<-widths
  arrangeGrob(g1,g2,ncol=1,heights=c(.77,.23))
}

page_grob <- function(grobs,title,ncol,nrow) {
  body <- arrangeGrob(grobs=grobs,ncol=ncol,nrow=nrow)
  arrangeGrob(textGrob(title,gp=gpar(fontsize=29,fontface='bold')),body,
    ncol=1,heights=c(.065,.935))
}

stage_plan <- list(list(cohort='occams_combined',stages='cstage'),list(cohort='occams_neo',stages=c('cstage','pstage')),
  list(cohort='occams_neo_biopsy',stages=c('cstage','pstage')),list(cohort='occams_neo_surgery',stages=c('cstage','pstage')),
  list(cohort='occams_surgery',stages=c('cstage','pstage')),list(cohort='tcga_all',stages=c('cstage','pstage')),
  list(cohort='tcga_yes',stages=c('cstage','pstage')))
stage_pages <- list()
for (page in seq_along(stage_plan)) {
  key <- stage_plan[[page]]$cohort;d<-meta[cohorts[[key]],,drop=FALSE]
  grobs <- list()
  for (endpoint in c('OS','RFS')) for (stage in stage_plan[[page]]$stages) grobs[[paste(endpoint,stage)]] <- km_panel(d,endpoint,stage,key,page)
  stage_pages[[page]] <- page_grob(grobs,paste0(cohort_labels[[key]],' | N=',nrow(d)),
    ncol=2,nrow=if(page==1L) 1 else 2)
}
km_panel_table <- do.call(rbind,km_panels)
km_risk_table <- do.call(rbind,km_risks)
write_table(do.call(rbind,km_sources),'km_curve_source_data.csv');write_table(km_risk_table,'km_number_at_risk.csv')
write_table(km_panel_table,'km_panel_inventory.csv')

####################
# Full MP descriptions occupy fixed left/right label columns. This avoids label
# collisions without dropping MPs or moving plotted HR/p-value coordinates.
method_labels <- c(continuous='Continuous',median='Median split',q1q4='Q1/Q4 split',optimal_cut='Optimal cut')
volcano_collections <- layout_collections <- list()
for (display in c('mp','all')) {
volcano_pages <- list()
label_layouts <- list()
for (page in seq_along(mp_cohorts)) {
  key <- mp_cohorts[page]
  all_fitted <- cox_results[cox_results$cohort==key & cox_results$status=='fitted' &
    (display=='all' | cox_results$feature_family=='MP'),,drop=FALSE]
  x_limit <- max(1,max(abs(log2(all_fitted$HR)))*1.22)
  y_limit <- max(2,max(-log10(pmax(all_fitted$P_value,.Machine$double.xmin)))*1.15)
  plots <- list()
  for (method in split_methods) {
    d <- all_fitted[all_fitted$method==method,,drop=FALSE]
    if (nrow(d)!=if (display=='mp') 17L else 27L) stop('Missing fitted feature models: ',key,' / ',method)
    d$significant <- d$P_value<.05
    d$x <- log2(d$HR); d$y <- -log10(pmax(d$P_value,.Machine$double.xmin))
    d$full_name <- paste(d$feature,unname(SCREF_MP_DESCRIPTIONS[d$feature]))
    is_state <- d$feature_family!='MP'
    d$full_name[is_state] <- d$feature[is_state]
    if (anyNA(d$full_name)) stop('Missing MP description')
    d$label <- vapply(d$full_name,function(x) paste(strwrap(x,width=26),collapse='\n'),character(1))
    order_x <- order(d$x,d$y)
    d$side <- 'right';d$side[order_x[seq_len(floor(nrow(d)/2))]] <- 'left'
    d$label_y <- NA_real_
    for (side in c('left','right')) {
      ix <- which(d$side==side); ix <- ix[order(d$y[ix],d$x[ix])]
      # Allocate nonoverlapping multiline slots instead of stochastic label repulsion.
      line_units <- lengths(strsplit(d$label[ix],'\n',fixed=TRUE)) + 1.2
      centres <- (cumsum(line_units)-line_units/2)/sum(line_units)
      d$label_y[ix] <- (.035+.93*centres)*y_limit
    }
    d$label_x <- ifelse(d$side=='left',-1.4,1.4)*x_limit
    d$label_hjust <- ifelse(d$side=='left',1,0)
    d$cohort <- key
    label_layouts[[paste(key,method)]] <- d[,c('cohort','method','feature','full_name','label','side','label_y','label_x')]
    count_range <- function(x) if (length(unique(x))==1L) as.character(x[1]) else paste(range(x),collapse='-')
    subtitle <- paste0('Overall Survival | N=',count_range(d$n),' | Events = ',count_range(d$events))
    plots[[method]] <- ggplot(d,aes(x,y))+
      annotate('segment',x=-x_limit,xend=x_limit,y=-log10(.05),yend=-log10(.05),linetype='dashed',colour='grey50',linewidth=.6)+
      annotate('segment',x=0,xend=0,y=0,yend=y_limit,linetype='dashed',colour='grey50',linewidth=.6)+
      geom_segment(aes(xend=label_x,yend=label_y),colour='grey70',linewidth=.3)+
      geom_point(aes(colour=significant,shape=feature_family),size=3.7)+
      scale_shape_manual(values=c('MP'=16,'State MP union'=17,'State DGE'=15),name='Signature family')+
      geom_text(aes(x=label_x,y=label_y,label=label,hjust=label_hjust),size=if(display=='all') 4.8 else 5.8,lineheight=.95,colour='black')+
      scale_colour_manual(values=c('FALSE'='#777777','TRUE'='#B22222'),guide='none')+
      scale_x_continuous(breaks=pretty(c(-x_limit,x_limit),n=3))+
      coord_cartesian(xlim=c(-4*x_limit,4*x_limit),ylim=c(0,y_limit),clip='on')+
      labs(title=method_labels[[method]],subtitle=subtitle,x='log2(hazard ratio)',y='-log10(raw Wald p)')+theme_readable+
      theme(panel.grid=element_blank(),plot.margin=margin(12,26,12,44))
  }
####################
  page_title <- if (key=='occams_surgery') 'OCCAMS surgery pathway | RNA matched 68/72' else
    paste0(cohort_labels[[key]],' | N=',length(cohorts[[key]]))
  volcano_pages[[page]] <- page_grob(lapply(plots,ggplotGrob),page_title,2,2)
####################
}
volcano_collections[[display]] <- volcano_pages
layout_collections[[display]] <- do.call(rbind,label_layouts)
}
volcano_pages <- volcano_collections$mp
write_table(layout_collections$mp,'volcano_label_layout.csv')
write_table(layout_collections$all,'all_feature_volcano_label_layout.csv')
# Quartile membership varies by MP; exact per-MP deaths remain auditable.
quartile_members <- list()
for (key in mp_cohorts) {
  d <- meta[cohorts[[key]],,drop=FALSE]
  d <- d[positive(d$OS_time) & d$OS_event %in% c(0,1),,drop=FALSE]
  for (feature in names(cache$mp_genes)) {
    scores <- cache$scores[match(d$sample_id,rownames(cache$scores)),feature]
    q <- quantile(scores,c(.25,.75),names=FALSE)
    keep <- scores<=q[1] | scores>=q[2]
    quartile_members[[paste(key,feature)]] <- data.frame(cohort=key,mp=feature,subject_id=d$subject_id[keep],
      score=scores[keep],group=ifelse(scores[keep]>=q[2],'High','Low'),OS_event=d$OS_event[keep])
  }
}
write_table(do.call(rbind,quartile_members),'q1q4_patient_membership.csv')
####################

####################
# Independent-patient specimen comparison; both groups use the same 157-parent GSVA.
neo <- meta[cohorts$occams_neo,,drop=FALSE]
box_source <- do.call(rbind,lapply(names(cache$feature_sets),function(f)
  data.frame(subject_id=neo$subject_id,group=neo$source,feature=f,
    family=cache$feature_family[[f]],value=cache$scores[match(neo$sample_id,rownames(cache$scores)),f])))
box_stats <- do.call(rbind,lapply(split(box_source,box_source$feature),function(d) {
  data.frame(feature=d$feature[1],family=d$family[1],biopsy_n=sum(d$group=='biopsy'),
    surgery_n=sum(d$group=='surgery'),biopsy_median=median(d$value[d$group=='biopsy']),
    surgery_median=median(d$value[d$group=='surgery']),
    p_value=wilcox.test(value~group,data=d,paired=FALSE,exact=FALSE)$p.value)
}))
box_stats$FDR_BH <- NA_real_
for (family in unique(box_stats$family)) {
  ix <- box_stats$family==family
  box_stats$FDR_BH[ix] <- p.adjust(box_stats$p_value[ix],method='BH')
}
write_table(box_source,'specimen_boxplot_source_data.csv')
write_table(box_stats,'specimen_boxplot_statistics.csv')
boxplots <- lapply(c('MP','States'),function(kind) {
  features <- names(cache$feature_sets)[if(kind=='MP') cache$feature_family=='MP' else cache$feature_family!='MP']
  d <- box_source[box_source$feature %in% features,]
  # Shared five-state x axis on the state page; rows distinguish signature definitions.
  d$axis_label <- if(kind=='MP') d$feature else sub(' \\((MP union|DGE)\\)$','',d$feature)
  axis_order <- unique(d$axis_label)
  d$axis_label <- factor(d$axis_label,levels=axis_order)
  d$group <- factor(d$group,levels=c('biopsy','surgery'))
  d$family <- factor(d$family,levels=c('MP','State MP union','State DGE'))
  ann <- box_stats[match(features,box_stats$feature),]
  ann$axis_label <- factor(if(kind=='MP') ann$feature else sub(' \\((MP union|DGE)\\)$','',ann$feature),levels=axis_order)
  ann$y_pos <- vapply(features,function(f) max(d$value[d$feature==f]),numeric(1))+.08*diff(range(d$value))
  ann$label <- ifelse(ann$FDR_BH<.001,'***',ifelse(ann$FDR_BH<.01,'**',ifelse(ann$FDR_BH<.05,'*','')))
  labels <- if(kind=='MP') setNames(paste(features,SCREF_MP_DESCRIPTIONS[features],sep=' | '),features) else setNames(axis_order,axis_order)
  p <- ggplot(d,aes(axis_label,value,fill=group,colour=group))+
    geom_boxplot(position=position_dodge(width=.75),width=.6,outlier.shape=NA,alpha=.8,linewidth=.4,colour='black')+
    geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75,seed=1701),alpha=.7,size=1,stroke=0,show.legend=FALSE)+
    geom_text(data=ann[ann$label!='',],aes(x=axis_label,y=y_pos,label=label),inherit.aes=FALSE,size=3,fontface='bold')+
    scale_fill_manual(values=c(biopsy='#0072B2',surgery='#D55E00'),labels=c('Biopsy (n=77)','Surgery (n=80)'))+
    scale_colour_manual(values=c(biopsy='#0072B2',surgery='#D55E00'),guide='none')+
    scale_x_discrete(labels=labels)+scale_y_continuous(expand=expansion(mult=c(.02,.18)))+
    labs(title=paste0('OCCAMS neoadjuvant | Biopsy versus surgery specimens | ',kind),
      subtitle='Independent patients; two-sided Wilcoxon tests. Stars: BH FDR <0.05, <0.01, <0.001 within signature family.',
      x=NULL,y='GSVA enrichment score',fill='RNA specimen')+
    coord_cartesian(clip='off')+theme_classic(base_size=12)+
    theme(plot.title=element_text(face='bold',size=15),plot.subtitle=element_text(size=10,colour='grey35'),
      axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.line.x=element_blank(),
      legend.position='top',legend.title=element_text(face='bold'),plot.margin=margin(12,16,12,12))
  if(kind=='States') p <- p+facet_wrap(~family,ncol=1,scales='free_y')
  ggplotGrob(p)
})
####################
write_pdf <- function(pages,path,width=24,height=17) {
  if (file.exists(path)) stop('Refusing to overwrite PDF: ',path)
  cairo_pdf(path,width=width,height=height,family='sans',onefile=TRUE)
  for (g in pages) {grid.newpage();grid.draw(g)}
  dev.off()
  if (!file.exists(path) || file.info(path)$size<1000) stop('PDF not produced: ',path)
}
write_pdf(stage_pages,file.path(dirs[['figures']],'selected_cohort_stage_survival.pdf'))
write_pdf(volcano_pages,file.path(dirs[['figures']],'selected_cohort_mp_os_volcanoes.pdf'))
write_pdf(volcano_collections$all,file.path(dirs[['figures']],'selected_cohort_mp_and_state_os_volcanoes.pdf'),width=32,height=26)
write_pdf(boxplots,file.path(dirs[['figures']],'occams_neoadjuvant_biopsy_vs_surgery_boxplots.pdf'),width=18,height=12)
stopifnot(length(stage_pages)==7L,length(volcano_pages)==6L,nrow(cox_results)==648L,length(boxplots)==2L)
summary_dir <- file.path(project_dir,'updates/new_updates/summaries')
dir.create(summary_dir,recursive=TRUE,showWarnings=FALSE)
summary_path <- file.path(summary_dir,paste0('selected_cohort_survival_',run_id,'.csv'))
if (file.exists(summary_path)) stop('Refusing to overwrite compact summary')
write.csv(cohort_summary,summary_path,row.names=FALSE)
stopifnot(file.info(summary_path)$size<100000)
report <- c('STATUS: PASS',paste('Run ID:',run_id),paste('Analysis cache:',cache_path),
  'PDF 1: 7 stage-survival pages; PDF 2: 6 OS MP-volcano pages.',
  'PDF 3: 6 OS volcano pages with all 17 MPs and both five-state signature types; all 27 features labelled per panel.',
  'PDF 4: 2 OCCAMS neoadjuvant specimen boxplot pages: MPs, then MP-union and DGE states in separate rows.',
  'Specimen tests: unpaired two-sided Wilcoxon, BH within each signature family; not paired treatment effects.',
  'OCCAMS: typed neoadjuvant 157 (biopsy 77; surgery specimens 80), matched surgery pathway 68, combined 225.',
  'TCGA: surgery pathway 86; surgery Yes 72, No 12, unresolved 2. Stage/RFS known/evaluable Yes subset 64.',
  'OCCAMS RFS: exploratory; user approved assuming a common first-recorded-therapy origin for recurrence and survival intervals.',
  'TCGA RFS: preserves earlier cohort/definition including one New Primary Tumor event; time origin is initial pathologic diagnosis.',
  'Three exact OCCAMS exclusions applied. Unknown stage excluded per panel; low-coverage KM panels are explicitly descriptive.',
  'MP OS: Gaussian GSVA as in original integrated script; parent-cohort scores reused in subgroups; continuous, median, Q1/Q4, within-cohort optimal cut.',
  'Cox models are univariable; optimal-cut raw p-values are exploratory, not search-adjusted.',
  paste('Fitted Cox rows:',sum(cox_results$status=='fitted'),'of',nrow(cox_results)),
  paste('Figures:',dirs[['figures']]),paste('Tables:',dirs[['tables']]),paste('Compact summary:',summary_path))
writeLines(report,file.path(dirs[['reports']],'run_summary.txt'))
message(paste(report,collapse='\n'))
message('SELECTED_COHORT_SURVIVAL_COMPLETE')
