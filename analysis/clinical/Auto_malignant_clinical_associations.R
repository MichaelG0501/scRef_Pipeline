####################
# Analysis registry
# Status: active
# Script: analysis/clinical/Auto_malignant_clinical_associations.R
# Description: Deconvoluted MP/state clinical associations in five OCCAMS and
#   two TCGA clinical selections. One patient is the statistical unit.
#   Age is <=50 versus >50. Clinical categories require at least three patients,
#   as requested; all retained categories enter the rank tests and FDR annotation.
#   Final site compares Oesophagus + GOJ I versus GOJ II/III.
#   Mandard is shown as TRG1-5 and TRG1-3 versus TRG4-5.
# Methodology: analysis/methodology/clinical/malignant_clinical_associations_methodology.md
# Metadata authority (multiple sources are required; no missing-value substitution):
#   analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx
#     Patient_Master_388: all covered OCCAMS clinical and cohort fields.
#   ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv
#     Comorbidity only, because the current workbook lacks this field.
#   analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx
#     EAC_RNA88_Final: TCGA staging, selection, sample IDs and treatment fields.
#   ref_outs/TCGA/surgery_clinical_associations/reference/clinical_cbioportal.tsv
#     TCGA demographic/exposure fields explicitly listed in tcga_spec below.
# Other exact live inputs:
#   ref_outs/clinical_malignant_survival/intermediate/analysis_cache.rds
#     Fixed 281/88-parent malignant GSVA scores and within-malignant RNA fractions.
#   ref_outs/clinical_selected_survival/tables/cohort_membership.csv (identity audit only)
#   ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   analysis/shared/scRef_config.R
# Outputs under ref_outs/clinical_malignant_associations/:
#   intermediate/clinical_associations.rds (full replot data and source contract)
#   tables/{input_manifest,subject_metadata,cohort_membership,variable_inventory,
#     association_statistics,plot_source_data,selected_scores,feature_definitions,
#     metadata_source_audit,figure_index}.csv
#   figures/Auto_final/{01_occams_surgery,02_tcga_yes}_{MP,State}.pdf
#     All eligible variables together for each surgery cohort; one variable per page.
#   figures/Auto_supporting/{occams_combined,occams_neo,occams_neo_biopsy,
#     occams_neo_surgery,tcga_all}_{MP,State}.pdf
#     All eligible variables in each remaining cohort; no duplicated final pages.
#   reports/{current_results.md,validation.csv}; logs/package_versions.csv
#   updates/new_updates/summaries/malignant_clinical_associations.csv (<100 KB)
# Additional output: figures/Auto_supporting/occams_neoadjuvant_biopsy_vs_surgery_boxplots.pdf;
#   tables/specimen_boxplot_{source_data,statistics}.csv. Two pages: MPs and state GSVA.
#   SCREF_SPECIMEN_ONLY=TRUE with SCREF_REPLOT_ONLY=TRUE exports this PDF only.
# Cache: changed analytical code/input/package contracts stop. Explicit
#   SCREF_FORCE_REBUILD=TRUE rebuilds; SCREF_REPLOT_ONLY=TRUE requires a matching cache.
#   Plot changes alone do not recompute scores or statistical tests. No downloads.
# Run: /opt/pbs/bin/qsub job_scripts/clinical/Auto_malignant_clinical_associations.sh
# Direct allocated-node run: Rscript analysis/clinical/Auto_malignant_clinical_associations.R
# Environment: dmtcp; 1 CPU, 8 GB, 01:00:00; PBS logs in the project root.
####################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

source('analysis/shared/scRef_config.R')
setwd(SCREF_PROJECT_DIR)
if (!nzchar(Sys.getenv('PBS_JOBID'))) stop('Run within a PBS allocation')
out <- file.path(SCREF_REF_OUTS_DIR, 'clinical_malignant_associations')
for (tier in c('intermediate','tables','reports','logs','figures/Auto_final','figures/Auto_supporting'))
  dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
inputs <- c(
  occams_workbook='analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx',
  occams_comorbidity='ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv',
  tcga_workbook='analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx',
  tcga_cbio='ref_outs/TCGA/surgery_clinical_associations/reference/clinical_cbioportal.tsv',
  score_cache='ref_outs/clinical_malignant_survival/intermediate/analysis_cache.rds',
  membership='ref_outs/clinical_selected_survival/tables/cohort_membership.csv',
  qc='ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv',
  grouping=SCREF_MP_GROUPING_CSV,config='analysis/shared/scRef_config.R')
if (!all(file.exists(inputs))) stop('Missing required clinical-association inputs')
manifest <- data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
force <- identical(toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE')),'TRUE')
replot <- identical(toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE')),'TRUE')
if (force && replot) stop('Force and replot modes are mutually exclusive')
minimum_n <- 3L
workbook_path <- inputs[['occams_workbook']]
metadata_path <- inputs[['occams_comorbidity']]
expected_cohort_n <- c(225L,157L,77L,80L,68L)
cohort_labels <- c(occams_combined='OCCAMS combined',
  occams_neo='OCCAMS neoadjuvant',
  occams_neo_biopsy='OCCAMS neoadjuvant biopsy',
  occams_neo_surgery='OCCAMS neoadjuvant surgery',
  occams_surgery='OCCAMS surgery pathway',
  tcga_all='TCGA surgery pathway',tcga_yes='TCGA surgery Yes')

variable_config <- data.frame(
  variable = c(
    "gender", "age_group", "diagnostic_combined_tumour_site", "pretreatment_combined_tumour_site",
    "pathology_location", "diagnostic_histology", "diagnostic_siewert_classification",
    "pretreatment_siewert_classification", "pathology_siewert_classification",
    "pretreatment_t_stage", "pretreatment_n_stage_tnm7", "pretreatment_m_stage",
    "treatment_plan_performance_status_group", "treatment_plan_intent",
    "treatment_plan_curative_modality", "surgery_main", "treatment_plan_comorbidity_group",
    "pathology_t_stage", "pathology_n_stage_tnm7", "pathology_m_stage",
    "pathology_positive_node_group", "pathology_grade_group", "pathology_signet_ring_cells",
    "pathology_barrett_adjacent_dysplasia", "pathology_neoadjuvant_therapy",
    "surgery_procedure_group", "pathology_mandard_score", "pathology_mandard_response_group",
    "treatment_response_group", "recurrence_ever_recorded"
  ),
  label = c(
    "Gender", "Age at diagnosis (>50)", "Diagnostic combined tumour site",
    "Final pretreatment combined tumour site", "Resection pathology location", "Diagnostic histology",
    "Diagnostic Siewert classification", "Final pretreatment Siewert classification",
    "Resection pathology Siewert classification", "Final pretreatment T stage",
    "Final pretreatment N stage (TNM7)", "Final pretreatment M stage", "Treatment-plan performance status",
    "Treatment-plan intent", "Treatment-plan curative modality", "Main surgery",
    "Treatment-plan comorbidity burden", "Resection pathology T stage",
    "Resection pathology N stage (TNM7)", "Resection pathology M stage",
    "Resection positive-node burden", "Resection pathology grade", "Resection signet-ring cells",
    "Resection adjacent Barrett-associated dysplasia", "Recorded neoadjuvant therapy at resection",
    "Operation type", "Mandard grade", "Mandard TRG1-3 vs TRG4-5",
    "Recorded treatment response grouping", "Disease recurrence ever recorded"
  ),
  category = c(
    rep("Demographic", 2), "Diagnostic", "Pretreatment", "Resection pathology",
    "Diagnostic", "Diagnostic", "Pretreatment", "Resection pathology",
    rep("Pretreatment", 3), rep("Treatment plan", 3), "Surgery", "Treatment plan",
    rep("Resection pathology", 8), "Surgery", rep("Resection pathology", 2),
    "Treatment response", "Longitudinal follow-up"
  ),
  rationale = c(
    "Demographic variable", "Demographic/prognostic variable", "RD anatomical site only",
    "PS anatomical site only", "RP anatomical location only", "RD disease class only",
    "RD GOJ anatomical class only", "PS GOJ anatomical class only", "RP GOJ anatomical class only",
    "PS pretreatment local extent", "PS pretreatment nodal extent using TNM7 only",
    "PS pretreatment metastatic status", "TP recorded fitness", "TP recorded treatment intent",
    "TP recorded curative modality", "ST recorded surgery", "TP recorded comorbidity context",
    "RP resection local extent", "RP resection nodal extent using TNM7 only",
    "RP resection metastatic status", "RP positive-node count grouping", "RP differentiation",
    "RP histological phenotype", "RP dysplasia field context", "RP neoadjuvant exposure record",
    "ST surgical approach", "RP pathological response", "RP grouped pathological response",
    "TR recorded response", "Any Yes recorded in the single recurrence source column"
  ),
  source_column = c(
    "DI Patient Gender", "DI Age At Diagnosis", "RD Combined Tumour Site",
    "PS Pre Treatment Combined Tumour Site", "RP Location", "RD Histology",
    "RD Siewert Classification", "PS Siewert Classification", "RP Siewert Classification",
    "PS TStage Primary Tumour Final Pretreatment Staging",
    "PS NStage Primary Tumour Final Pretreatment Staging TNM7",
    "PS MStage Primary Tumour Final Pretreatment Staging", "TP Performance Status",
    "TP Treatment Intent", "TP Curative Treatment Modality", "ST Main Surgery", "TP Comorbidities",
    "RP TStage Primary Tumour", "RP Nstage RP TNM7", "RP MStage Distant Metastasis",
    "RP NStage Number Of Positive Nodes", "RP Tumour Grading Differentiation Status",
    "RP Signet Ring Cells Present", "RP Baretts Adjacent To Tumour Microscopic Dysplasia",
    "RP History Of Neo Adjuvant Therapy", "ST Procedure", "RP Mandard Score For Response",
    "RP Mandard Score For Response", "TR Response", "Has Original Disease Reoccurred"
  ),
  timepoint = c(
    "DI", "DI", "RD", "PS", "RP", "RD", "RD", "PS", "RP", "PS", "PS", "PS",
    "TP", "TP", "TP", "ST", "TP", "RP", "RP", "RP", "RP", "RP", "RP", "RP",
    "RP", "ST", "RP", "RP", "TR", "longitudinal source field"
  ),
  stringsAsFactors = FALSE
)

####################
# A source is chosen explicitly, never selected by file modification time or by
# filling missing curated values from the superseded metadata.
modern_columns <- c(gender="Gender", age_group="Age_At_Diagnosis",
  diagnostic_combined_tumour_site="RD_Combined_Tumour_Site", pretreatment_combined_tumour_site="PS_Combined_Tumour_Site",
  pathology_location="RP_Location", diagnostic_histology="RD_Histology", diagnostic_siewert_classification="RD_Siewert",
  pretreatment_siewert_classification="PS_Siewert", pathology_siewert_classification="RP_Siewert",
  pretreatment_t_stage="cT", pretreatment_n_stage_tnm7="cN", pretreatment_m_stage="cM",
  treatment_plan_performance_status_group="TP_Performance_Status_All", treatment_plan_intent="Treatment_Intent_Final",
  treatment_plan_curative_modality="Curative_Treatment_Modality_Final", surgery_main="Surgery_Ever",
  pathology_t_stage="pT", pathology_n_stage_tnm7="pN", pathology_m_stage="pM",
  pathology_positive_node_group="pN_Positive_Nodes", pathology_grade_group="RP_Grade",
  pathology_signet_ring_cells="RP_Signet_Ring_Cells", pathology_barrett_adjacent_dysplasia="RP_Barretts_Adjacent_Dysplasia",
  pathology_neoadjuvant_therapy="RP_Neoadjuvant_History", surgery_procedure_group="Surgery_Procedures_All",
  pathology_mandard_score="RP_Mandard_Score", pathology_mandard_response_group="RP_Mandard_Score",
  treatment_response_group="Curative_Neoadjuvant_Clinical_Response", recurrence_ever_recorded="Recurrence_Ever")
ix <- match(names(modern_columns), variable_config$variable)
variable_config$source_column[ix] <- unname(modern_columns)
variable_config$label[variable_config$variable=="pretreatment_n_stage_tnm7"] <- "Harmonized pretreatment N stage"
variable_config$label[variable_config$variable=="pathology_n_stage_tnm7"] <- "Harmonized pathology N stage"
variable_config$label[variable_config$variable=="treatment_plan_intent"] <- "Final initial treatment intent"
variable_config$label[variable_config$variable=="treatment_plan_curative_modality"] <- "Final initial curative treatment pathway"
variable_config$label[variable_config$variable=="surgery_main"] <- "Surgery ever recorded"
variable_config$label[variable_config$variable=="treatment_response_group"] <- "Curative neoadjuvant clinical response"
variable_config$rationale[ix] <- paste0("Corrected patient master: ", unname(modern_columns))
variable_config$timepoint[variable_config$variable=="treatment_response_group"] <- "curative neoadjuvant"
variable_config <- bind_rows(variable_config, data.frame(
  variable=c("clinical_stage","pathological_stage","rna_specimen_source","rna_treatment_status","final_site","final_siewert","neoadjuvant_type"),
  label=c("Clinical AJCC stage","Pathological AJCC stage","RNA specimen source","RNA treatment status","Curated final anatomical site","Curated final Siewert classification","Curative neoadjuvant treatment type"),
  category=c("Pretreatment","Resection pathology","RNA specimen","RNA specimen","Curated diagnosis","Curated diagnosis","Treatment plan"),
  rationale="Explicit corrected patient-master field; no older-value substitution",
  source_column=c("Clinical_Stage_AJCC","Pathological_Stage_AJCC","RNA_specimen_source","RNA_treatment_status","Site_Final","Siewert_Final","Curative_Neoadjuvant_Type"),
  timepoint=c("clinical staging","pathological staging","RNA collection","RNA collection","curated final","curated final","curative neoadjuvant")))
####################
variable_config <- bind_rows(variable_config,data.frame(
  variable="final_site_siewert_group",label="Oesophagus + GOJ I vs GOJ II/III",
  category="Curated diagnosis",rationale="Requested combined final anatomical site and Siewert grouping",
  source_column="Site_Final + Siewert_Final",timepoint="curated final"))
####################
variable_config$source_file <- ifelse(variable_config$variable=="treatment_plan_comorbidity_group",metadata_path,workbook_path)
variable_config$source_sheet <- ifelse(variable_config$variable=="treatment_plan_comorbidity_group","not applicable","Patient_Master_388")

clean_curated <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c("", "na", "n/a", "unknown", "not_recorded", "not recorded", "tx", "nx", "mx")] <- NA_character_
  x
}
consensus_tokens <- function(x) {
  vapply(x,function(v) {
    tokens <- unique(na.omit(clean_curated(strsplit(ifelse(is.na(v),"",v),"[,;|]")[[1]])))
    if(length(tokens)==1L) tokens else NA_character_
  },character(1))
}
resolve_metadata <- function(older, qc_ids) {
  w <- as.data.frame(read_excel(workbook_path,sheet="Patient_Master_388"))
  if(nrow(w)!=388L || anyDuplicated(w$SHA_ID) || !all(qc_ids %in% w$SHA_ID)) stop("Invalid corrected patient master")
  w <- w[match(qc_ids,w$SHA_ID),,drop=FALSE]
  older <- older[match(qc_ids,older$subject_id),,drop=FALSE]
  m <- data.frame(subject_id=qc_ids,stringsAsFactors=FALSE)
  for(i in seq_len(nrow(variable_config))) {
    variable <- variable_config$variable[i]; column <- variable_config$source_column[i]
    if(variable=="final_site_siewert_group") next
    if(variable=="treatment_plan_comorbidity_group") m[[variable]] <- clean_curated(older[[variable]])
    else {
      if(!column %in% names(w)) stop("Missing curated column: ",column)
      m[[variable]] <- consensus_tokens(clean_curated(w[[column]]))
    }
  }
  ####################
  known_siewert <- m$final_siewert %in% c("1","2","3")
  if(any(known_siewert & (is.na(m$final_site) | m$final_site!="GOJ")))
    stop("Recorded final Siewert class outside the current GOJ metadata contract")
  m$final_site_siewert_group <- ifelse(
    m$final_site %in% c("Oesophagus lower","Oesophagus middle") |
      (m$final_site %in% "GOJ" & m$final_siewert %in% "1"),"Oesophagus + GOJ I",
    ifelse(m$final_site %in% "GOJ" & m$final_siewert %in% c("2","3"),"GOJ II/III",NA_character_))
  ####################
  age <- suppressWarnings(as.numeric(m$age_group))
  m$age_group <- ifelse(is.finite(age),ifelse(age>50,">50","<=50"),NA_character_)
  m$treatment_plan_performance_status_group <- consensus_tokens(m$treatment_plan_performance_status_group)
  m$treatment_plan_performance_status_group <- ifelse(m$treatment_plan_performance_status_group=="0","0",
    ifelse(m$treatment_plan_performance_status_group %in% c("1","2"),"1-2",NA_character_))
  n <- suppressWarnings(as.numeric(m$pathology_positive_node_group))
  m$pathology_positive_node_group <- ifelse(is.finite(n) & n>=0,ifelse(n==0,"0",ifelse(n<=3,"1-3",">=4")),NA_character_)
  grade <- tolower(consensus_tokens(m$pathology_grade_group))
  m$pathology_grade_group <- ifelse(grade %in% c("well","moderate_to_well","moderate"),"Well/moderate",
    ifelse(grade %in% c("moderate_to_poor","poor"),"Moderate/poor",NA_character_))
  barrett <- tolower(consensus_tokens(m$pathology_barrett_adjacent_dysplasia))
  m$pathology_barrett_adjacent_dysplasia <- ifelse(grepl("low grade|high grade",barrett),"Dysplasia present",
    ifelse(barrett=="no dysplasia","No dysplasia",NA_character_))
  trg <- toupper(consensus_tokens(m$pathology_mandard_score))
  trg[!trg %in% paste0("TRG",1:5)] <- NA_character_
  m$pathology_mandard_score <- trg
  m$pathology_mandard_response_group <- ifelse(trg %in% c("TRG1","TRG2","TRG3"),"TRG1-3",
    ifelse(trg %in% c("TRG4","TRG5"),"TRG4-5",NA_character_))
  response <- toupper(consensus_tokens(m$treatment_response_group))
  m$treatment_response_group <- ifelse(response %in% c("CR","PR"),"CR/PR",ifelse(response %in% c("SD","PD"),"SD/PD",NA_character_))
  procedure <- consensus_tokens(m$surgery_procedure_group)
  m$surgery_procedure_group <- ifelse(is.na(procedure),NA_character_,ifelse(grepl("Ivor.Lewis",procedure),"Ivor-Lewis oesophagectomy",
    ifelse(grepl("oesophagectomy|MIO",procedure,ignore.case=TRUE),"Other oesophagectomy",
      ifelse(grepl("gastrectomy",procedure,ignore.case=TRUE),"Gastrectomy","Other operation"))))
  for(variable in c("clinical_stage","pathological_stage")) m[[variable]][!m[[variable]] %in% as.character(1:4)] <- NA_character_
  for(variable in c("pretreatment_t_stage","pathology_t_stage")) {
    value <- toupper(consensus_tokens(m[[variable]]))
    m[[variable]] <- ifelse(grepl("^T[0-4]",value),substr(value,1,2),NA_character_)
  }
  for(variable in c("pretreatment_n_stage_tnm7","pathology_n_stage_tnm7")) {
    value <- toupper(consensus_tokens(m[[variable]]))
    m[[variable]] <- ifelse(value %in% paste0("N",0:3),value,NA_character_)
  }
  excluded <- c("26d51bb9b7c9974b5e753fc3520ffcfa2abdbd132509a887f4a5bc5630aa05d9",
    "72089d94e1a29cc0b7df4fb8f9694343cda2c6e71ae1083fc266b5b36b62323d",
    "9ea99e16eaaf486784cb0cd7e1c846a4dcc7bbaecfccc9a04d7853093f2cf775")
  eligible <- !m$subject_id %in% excluded
  matched <- !is.na(w$RNA_Metadata_Matched) & w$RNA_Metadata_Matched=="Yes"
  curative <- !is.na(w$Treatment_Intent_Final) & w$Treatment_Intent_Final=="curative"
  neo <- eligible & matched & curative & !is.na(w$Curative_Treatment_Modality_Final) &
    w$Curative_Treatment_Modality_Final=="neoadjuvant (chemo + CRT)" & w$RNA_specimen_source %in% c("biopsy","surgery")
  surgery <- eligible & matched & curative & !is.na(w$Curative_Treatment_Modality_Final) & w$Curative_Treatment_Modality_Final=="surgery"
  cohorts <- list(occams_combined=which(neo|surgery),occams_neo=which(neo),
    occams_neo_biopsy=which(neo & w$RNA_specimen_source=="biopsy"),
    occams_neo_surgery=which(neo & w$RNA_specimen_source=="surgery"),occams_surgery=which(surgery))
  if(!identical(as.integer(lengths(cohorts)),expected_cohort_n)) stop("Cohorts differ from selected-survival contract")
  audit <- variable_config[,c('variable','label','source_column','source_file','source_sheet')]
  list(metadata=m,cohorts=cohorts,audit=audit)
}
####################

clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c('','na','n/a','unknown','not reported','not available','not applicable','[not available]','gx','nx','mx','tx')] <- NA_character_
  x
}
yes_no <- function(x) { z<-tolower(clean(x)); ifelse(z %in% c('true','yes'),'Yes',ifelse(z %in% c('false','no'),'No',NA_character_)) }
stage <- function(x) {z<-clean(x); ifelse(z %in% as.character(1:4),z,NA_character_)}
tn <- function(x,letter) {z<-toupper(clean(x)); ifelse(grepl(paste0('^',letter,'[0-4]'),z),substr(z,1,2),NA_character_)}
tcga_spec <- data.frame(variable=c('sex','age_group','grade','pathological_stage','clinical_stage','pT','pN','cT','cN','pM','positive_nodes','location','alcohol','smoking','prior_malignancy','prior_treatment','race','ethnicity','preoperative_treatment'),
  label=c('Sex','Age at diagnosis (>50 years)','Histological grade','Pathological AJCC8 stage','Clinical AJCC8 stage','Pathological T stage','Pathological nodal involvement','Clinical T stage','Clinical N stage','Pathological M classification','Positive lymph-node burden','Tumour location','Recorded alcohol history','Smoking exposure (pack-years)','Prior malignancy','Recorded prior treatment','Recorded race','Recorded ethnicity','Recorded preoperative treatment'),
  source=c('cbio','cbio',rep('workbook',10),rep('cbio',6),'workbook'),
  field=c('Sex','Diagnosis Age','Grade_AJCC8_used','Pathological_Stage_AJCC8','Clinical_Stage_AJCC8','pT_AJCC8_used','pN_AJCC8_used','cT_AJCC8_used','cN_AJCC8_used','M_AJCC8_used_for_path','Positive_nodes_HE','Tumour_location','Alcohol History Documented','Person Cigarette Smoking History Pack Year Value','Prior Malignancy','Prior Treatment','Race Category','Ethnicity Category','Treatment_prior_to_surgery_raw'),
  role=c(rep('Auto_final',4),rep('Auto_supporting',15)),stringsAsFactors=FALSE)

read_tcga_metadata <- function(parent_metadata) {
  w <- as.data.frame(read_excel(inputs[['tcga_workbook']],sheet='EAC_RNA88_Final'))
  cb <- read.delim(inputs[['tcga_cbio']],check.names=FALSE,stringsAsFactors=FALSE)
  cb <- cb[cb[['Study ID']]=='esca_tcga_gdc' & cb[['Sample Type']]=='Primary Tumor',,drop=FALSE]
  stopifnot(nrow(w)==88L,!anyDuplicated(w$Patient_ID),!anyDuplicated(cb[['Patient ID']]))
  j <- match(w$Patient_ID,cb[['Patient ID']])
  k <- match(w$Patient_ID,parent_metadata$subject_id)
  stopifnot(!anyNA(j),!anyNA(k),identical(w$RNA_primary_sample_ID,parent_metadata$sample_id[k]))
  cb <- cb[j,,drop=FALSE]
  m <- data.frame(subject_id=w$Patient_ID,sample_id=w$RNA_primary_sample_ID)
  for (i in seq_len(nrow(tcga_spec))) {
    column <- tcga_spec$field[i]
    source <- if (tcga_spec$source[i]=='workbook') w else cb
    if (!column %in% names(source)) stop('Missing TCGA source column: ',column)
    m[[tcga_spec$variable[i]]] <- clean(source[[column]])
  }
  age <- suppressWarnings(as.numeric(m$age_group)); m$age_group <- ifelse(is.finite(age),ifelse(age>50,'>50','<=50'),NA_character_)
  m$grade <- ifelse(m$grade %in% c('G1','G2'),'G1/G2',ifelse(m$grade %in% c('G3','G4'),'G3/G4',NA_character_))
  for(v in c('clinical_stage','pathological_stage')) m[[v]] <- stage(m[[v]])
  for(v in c('pT','cT')) m[[v]] <- tn(m[[v]],'T')
  m$cN <- tn(m$cN,'N'); n <- tn(m$pN,'N'); m$pN <- ifelse(n=='N0','N0',ifelse(n %in% c('N1','N2','N3'),'N+',NA_character_))
  m$pM[!m$pM %in% c('M0','M1')] <- NA_character_
  n <- suppressWarnings(as.numeric(m$positive_nodes)); m$positive_nodes <- ifelse(is.finite(n)&n>=0,ifelse(n==0,'0',ifelse(n<=3,'1-3','>=4')),NA_character_)
  n <- suppressWarnings(as.numeric(m$smoking)); m$smoking <- ifelse(is.finite(n)&n>=0,ifelse(n==0,'0',ifelse(n<=20,'>0-20','>20')),NA_character_)
  for(v in c('alcohol','prior_malignancy','prior_treatment')) m[[v]] <- yes_no(m[[v]])

  surgery <- !is.na(w$Flowchart_Pathway_Final) & w$Flowchart_Pathway_Final=='Surgery group'
  yes <- surgery & !is.na(w$Flowchart_Surgery_Status_Final) & w$Flowchart_Surgery_Status_Final=='Yes'
  cohorts <- list(tcga_all=which(surgery),tcga_yes=which(yes))
  stopifnot(identical(as.integer(lengths(cohorts)),c(86L,72L)))
  audit <- data.frame(variable=tcga_spec$variable,label=tcga_spec$label,source_column=tcga_spec$field,
    source_file=ifelse(tcga_spec$source=='workbook',inputs[['tcga_workbook']],inputs[['tcga_cbio']]),
    source_sheet=ifelse(tcga_spec$source=='workbook','EAC_RNA88_Final','not applicable'))
  list(metadata=m,cohorts=cohorts,audit=audit)
}

build_associations <- function() {
  parent <- readRDS(inputs[['score_cache']])
  stopifnot(all(file.exists(parent$manifest$path)),
    identical(unname(tools::md5sum(parent$manifest$path)),parent$manifest$md5))
  definitions <- parent$features$definitions
  definitions <- definitions[definitions$family!='State logratio',,drop=FALSE]
  stopifnot(nrow(definitions)==32L,sum(definitions$family=='MP')==17L)
  qc <- read.csv(inputs[['qc']],stringsAsFactors=FALSE)
  qc_ids <- qc$Subject_ID[qc$Annotation_Compatible_Below_50_Flag=='no']
  stopifnot(length(qc_ids)==281L,!anyDuplicated(qc_ids))
  old <- read.csv(inputs[['occams_comorbidity']],check.names=FALSE)
  stopifnot(all(qc_ids %in% old$subject_id),!anyDuplicated(old$subject_id))
  occams <- resolve_metadata(old,qc_ids)
  occams$metadata$sample_id <- occams$metadata$subject_id
  tcga <- read_tcga_metadata(parent$clinical$meta)
  resolved <- list(OCCAMS=occams,TCGA=tcga)
  old_membership <- read.csv(inputs[['membership']])
  score_parts <- metadata_parts <- membership_parts <- source_parts <- list()
  for (dataset in names(resolved)) {
    obj <- resolved[[dataset]]
    meta <- obj$metadata
    expression <- parent$scored$scores[[paste0(dataset,'__malignant')]]
    fractions <- parent$scored$fractions[[dataset]]
    x <- cbind(expression,fractions)[meta$sample_id,definitions$feature,drop=FALSE]
    stopifnot(nrow(x)==nrow(meta),all(is.finite(x)),!anyDuplicated(meta$subject_id))
    score_parts[[dataset]] <- data.frame(subject_id=meta$subject_id,sample_id=meta$sample_id,dataset=dataset,
      x,check.names=FALSE) %>% pivot_longer(all_of(definitions$feature),names_to='feature',values_to='score')
    metadata_parts[[dataset]] <- mutate(meta,dataset=dataset)
    source_parts[[dataset]] <- mutate(obj$audit,dataset=dataset)
    for (cohort in names(obj$cohorts)) {
      ids <- meta$subject_id[obj$cohorts[[cohort]]]
      stopifnot(setequal(ids,old_membership$subject_id[old_membership$cohort==cohort]))
      membership_parts[[cohort]] <- data.frame(cohort=cohort,subject_id=ids,dataset=dataset)
    }
  }
  scores <- bind_rows(score_parts) %>% left_join(definitions,by='feature')
  membership <- bind_rows(membership_parts)
  inventory <- plot_parts <- list()
  for (dataset in names(resolved)) {
    obj <- resolved[[dataset]]
    spec <- if (dataset=='OCCAMS') data.frame(variable=variable_config$variable,label=variable_config$label) else tcga_spec[,c('variable','label')]
    for (cohort in names(obj$cohorts)) {
      m <- obj$metadata[obj$cohorts[[cohort]],,drop=FALSE]
      for (i in seq_len(nrow(spec))) {
        variable <- spec$variable[i]
        values <- m[[variable]]
        counts <- table(values)
        ####################
        # Apply the user-requested minimum of three patients per clinical category.
        retained <- names(counts)[counts>=minimum_n]
        eligible <- length(retained)>=2L
        test_eligible <- eligible
        ####################
        priority <- cohort %in% c('occams_surgery','tcga_yes')
        key <- paste(cohort,variable,sep='__')
        inventory[[key]] <- data.frame(cohort=cohort,dataset=dataset,variable=variable,label=spec$label[i],
          cohort_n=nrow(m),n_available=sum(!is.na(values)),n_missing=sum(is.na(values)),
          levels_and_counts=paste(names(counts),as.integer(counts),sep='=',collapse=' | '),
          retained_levels=paste(retained,collapse=' | '),n_retained=sum(counts[retained]),
          excluded_small_levels=paste(names(counts)[counts<minimum_n],collapse=' | '),
          n_excluded_small=sum(counts[counts<minimum_n]),
          eligible=eligible,test_eligible=test_eligible,
          reason=if(eligible) '' else 'Fewer than two categories with at least 3 patients',
          test_reason=if(!eligible) 'Fewer than two categories with at least 3 patients' else '',
          role=if(priority) 'Auto_final' else 'Auto_supporting')
        if (eligible) {
          groups <- data.frame(subject_id=m$subject_id,group=values) %>% filter(group %in% retained)
          plot_parts[[key]] <- inner_join(groups,scores[scores$dataset==dataset,],by='subject_id') %>%
            mutate(cohort=cohort,variable=variable)
        }
      }
    }
  }
  plot_data <- bind_rows(plot_parts)
  statistics <- plot_data %>% group_by(cohort,variable,feature,family) %>% group_modify(~{
    d <- .x
    ####################
    counts <- table(d$group)
    group_summary <- paste(vapply(split(d,d$group),function(x)paste0(x$group[1],': n=',nrow(x),
      ', median=',signif(median(x$score),5)),character(1)),collapse=' | ')
    ####################
    if (diff(range(d$score))==0) return(data.frame(n=nrow(d),groups=n_distinct(d$group),
      test='not estimated',status='constant score',p=NA_real_,group_summary=group_summary))
    k <- n_distinct(d$group)
    test <- withCallingHandlers(if (k==2L) wilcox.test(score~group,data=d,exact=FALSE) else
      kruskal.test(score~group,data=d),warning=function(w)stop(conditionMessage(w)))
    stopifnot(is.finite(test$p.value))
    data.frame(n=nrow(d),groups=k,test=if(k==2L)'Wilcoxon rank-sum' else 'Kruskal-Wallis',
      status='estimated',p=test$p.value,
      group_summary=group_summary)
  }) %>% ungroup() %>% group_by(cohort,variable) %>% mutate(q=p.adjust(p,'BH',n=32L)) %>%
    group_by(cohort,variable,family) %>% mutate(q_family=p.adjust(p,'BH',n=if(first(family)=='MP')17L else 5L)) %>% ungroup()
  statistics$significance <- ifelse(statistics$status!='estimated','Not estimated',
    ifelse(statistics$q<.05,'FDR',ifelse(statistics$p<.05,'Nominal','NS')))
  stopifnot(nrow(statistics)==sum(bind_rows(inventory)$eligible)*32L)
  list(metadata=bind_rows(metadata_parts),membership=membership,scores=scores,definitions=definitions,
    audit=bind_rows(source_parts),inventory=bind_rows(inventory),plot_data=plot_data,statistics=statistics)
}

packages <- c('readxl','dplyr','tidyr','ggplot2')
versions <- data.frame(package=c('R',packages),version=c(as.character(getRversion()),
  vapply(packages,function(p)as.character(packageVersion(p)),character(1))))
contract <- list(version=5L,inputs=manifest,versions=versions,minimum_n=minimum_n,
  occams_config=variable_config,tcga_config=tcga_spec,
  code=lapply(list(clean_curated,consensus_tokens,resolve_metadata,clean,yes_no,stage,tn,
    read_tcga_metadata,build_associations),function(f)paste(deparse(body(f)),collapse='\n')))
cache_path <- file.path(out,'intermediate/clinical_associations.rds')
if (file.exists(cache_path) && !force) {
  objects <- readRDS(cache_path)
  if (!identical(objects$contract,contract)) stop('Clinical-association cache changed; explicit rebuild required')
} else {
  if (replot) stop('Replot requires the matching clinical-association cache')
  objects <- build_associations()
  objects$contract <- contract
  saveRDS(objects,cache_path)
}
write_table <- function(x,name) write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
write_table(manifest,'input_manifest');write_table(objects$metadata,'subject_metadata')
write_table(objects$membership,'cohort_membership');write_table(objects$inventory,'variable_inventory')
write_table(objects$statistics,'association_statistics');write_table(objects$plot_data,'plot_source_data')
write_table(objects$scores,'selected_scores');write_table(objects$definitions,'feature_definitions')
write_table(objects$audit,'metadata_source_audit')
write.csv(versions,file.path(out,'logs/package_versions.csv'),row.names=FALSE)

####################
specimen_only <- identical(toupper(Sys.getenv('SCREF_SPECIMEN_ONLY','FALSE')),'TRUE')
if (specimen_only && !replot) stop('SCREF_SPECIMEN_ONLY requires SCREF_REPLOT_ONLY=TRUE')
figure_index <- if (specimen_only) list(read.csv(file.path(out,'tables/figure_index.csv'))) else list()
####################
for (cohort in if(specimen_only) character() else names(cohort_labels)) {
  role <- if(cohort %in% c('occams_surgery','tcga_yes')) 'Auto_final' else 'Auto_supporting'
  variables <- objects$inventory$variable[objects$inventory$cohort==cohort & objects$inventory$eligible]
  for (type in c('MP','State')) {
    filename <- paste0(if(role=='Auto_final') if(cohort=='occams_surgery') '01_' else '02_' else '',
      cohort,'_',type,'.pdf')
    relative_path <- file.path('figures',role,filename)
    cairo_pdf(file.path(out,relative_path),width=if(type=='MP')18 else 16,
      height=if(type=='MP')10 else 12,onefile=TRUE)
    for (page in seq_along(variables)) {
      variable <- variables[page]
      features <- objects$definitions$feature[if(type=='MP')objects$definitions$family=='MP' else objects$definitions$family!='MP']
      d <- objects$plot_data %>% filter(.data$cohort==.env$cohort,.data$variable==.env$variable,feature %in% features)
      d$display_score <- ifelse(d$family=='State fraction',100*d$score,d$score)
      labels <- if(type=='MP') setNames(vapply(features,function(f)
        paste(c(f,strwrap(SCREF_MP_DESCRIPTIONS[[f]],width=18)),collapse='\n'),character(1)),features) else
        setNames(objects$definitions$state[match(features,objects$definitions$feature)],features)
      d$x <- factor(unname(labels[d$feature]),levels=unique(unname(labels[features])))
      d$family <- factor(d$family,levels=c('MP','State MP union','State DGE','State fraction'))
      ####################
      groups <- sort(unique(d$group))
      d$group <- factor(d$group,levels=groups)
      counts <- d %>% distinct(subject_id,group) %>% count(group,.drop=FALSE)
      ####################
      palette <- setNames(scales::hue_pal()(length(groups)),groups)
      legend_labels <- setNames(paste0(counts$group,' (n=',counts$n,')'),as.character(counts$group))
      ann <- d %>% group_by(feature,x,family) %>% summarise(
        y=max(display_score)+.10*diff(range(display_score))+.03,.groups='drop') %>%
        left_join(objects$statistics %>% filter(.data$cohort==.env$cohort,.data$variable==.env$variable),
          by=c('feature','family'))
      ann$star <- ifelse(is.na(ann$q_family),'',ifelse(ann$q_family<.001,'***',
        ifelse(ann$q_family<.01,'**',ifelse(ann$q_family<.05,'*',''))))
      ann$family <- factor(ann$family,levels=levels(d$family))
      info <- objects$inventory[objects$inventory$cohort==cohort & objects$inventory$variable==variable,]
      p <- ggplot(d,aes(x,display_score,fill=group,colour=group))+
        geom_boxplot(position=position_dodge(.75),width=.6,outlier.shape=NA,alpha=.8,
          linewidth=.4,colour='black')+
        geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75,seed=1701),
          size=1,alpha=.65,stroke=0,show.legend=FALSE)+
        geom_text(data=ann[ann$star!='',],aes(x=x,y=y,label=star),inherit.aes=FALSE,size=4.5,fontface='bold')+
        scale_fill_manual(values=palette,labels=legend_labels,drop=FALSE)+
        scale_colour_manual(values=palette,guide='none',drop=FALSE)+
        scale_y_continuous(expand=expansion(mult=c(.04,.16)))+
        labs(title=paste(cohort_labels[[cohort]],'|',info$label),x=NULL,
          y=if(type=='MP')'GSVA enrichment score' else 'Score / fraction (%)',fill=NULL)+
        theme_classic(base_size=14)+theme(plot.title=element_text(size=18,face='bold'),
          axis.text=element_text(size=12.5,colour='black'),legend.position='top',
          legend.text=element_text(size=14),strip.text=element_text(size=14),
          plot.margin=margin(12,16,12,12))
      if(type=='MP') p <- p+theme(axis.text.x=element_text(size=12.5,angle=45,hjust=1,vjust=1))
      if(type=='State') p <- p+facet_wrap(~family,ncol=1,scales='free_y',labeller=labeller(family=c(
        'State MP union'='MP-union score','State DGE'='DGE score','State fraction'='Malignant RNA (%)')))+
        scale_x_discrete(labels=function(x)vapply(x,function(s)paste(strwrap(s,width=18),collapse='\n'),character(1)))
      print(p)
      figure_index[[length(figure_index)+1L]] <- data.frame(role=role,path=relative_path,page=page,
        cohort=cohort,variable=variable,feature_type=type,features=length(features),n=info$n_retained)
    }
    dev.off()
  }
}
####################
# The requested specimen comparison reuses the malignant parent scores and tests.
# Patients are independent; the two specimen groups are not paired measurements.
specimen_path <- 'figures/Auto_supporting/occams_neoadjuvant_biopsy_vs_surgery_boxplots.pdf'
specimen_source <- objects$plot_data %>% filter(cohort=='occams_neo',
  variable=='rna_specimen_source',family!='State fraction')
specimen_stats <- objects$statistics %>% filter(cohort=='occams_neo',
  variable=='rna_specimen_source',family!='State fraction')
counts <- specimen_source %>% distinct(subject_id,group) %>% count(group)
stopifnot(nrow(counts)==2L,counts$n[counts$group=='biopsy']==77L,
  counts$n[counts$group=='surgery']==80L,nrow(specimen_stats)==27L)
write_table(specimen_source,'specimen_boxplot_source_data')
write_table(specimen_stats,'specimen_boxplot_statistics')
cairo_pdf(file.path(out,specimen_path),width=18,height=12,onefile=TRUE)
for (kind in c('MP','States')) {
  features <- objects$definitions$feature[if(kind=='MP') objects$definitions$family=='MP' else
    objects$definitions$family %in% c('State MP union','State DGE')]
  d <- specimen_source[specimen_source$feature %in% features,]
  d$axis_label <- if(kind=='MP') d$feature else d$state
  axis_order <- if(kind=='MP') features else unique(objects$definitions$state[match(features,objects$definitions$feature)])
  d$axis_label <- factor(d$axis_label,levels=axis_order)
  d$group <- factor(d$group,levels=c('biopsy','surgery'))
  d$family <- factor(d$family,levels=c('MP','State MP union','State DGE'))
  ann <- d %>% group_by(feature,axis_label,family) %>% summarise(
    y_pos=max(score)+.08*diff(range(d$score)),.groups='drop') %>%
    left_join(specimen_stats,by=c('feature','family'))
  ann$family <- factor(ann$family,levels=levels(d$family))
  ann$star <- ifelse(is.na(ann$q_family),'',ifelse(ann$q_family<.001,'***',
    ifelse(ann$q_family<.01,'**',ifelse(ann$q_family<.05,'*',''))))
  labels <- if(kind=='MP') setNames(paste(features,SCREF_MP_DESCRIPTIONS[features],sep=' | '),features) else
    setNames(axis_order,axis_order)
  p <- ggplot(d,aes(axis_label,score,fill=group,colour=group))+
    geom_boxplot(position=position_dodge(.75),width=.6,outlier.shape=NA,alpha=.8,linewidth=.4,colour='black')+
    geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75,seed=1701),
      alpha=.7,size=1,stroke=0,show.legend=FALSE)+
    geom_text(data=ann[ann$star!='',],aes(x=axis_label,y=y_pos,label=star),inherit.aes=FALSE,size=4,fontface='bold')+
    scale_fill_manual(values=c(biopsy='#0072B2',surgery='#D55E00'),labels=c('Biopsy (n=77)','Surgery (n=80)'))+
    scale_colour_manual(values=c(biopsy='#0072B2',surgery='#D55E00'),guide='none')+
    scale_x_discrete(labels=labels)+scale_y_continuous(expand=expansion(mult=c(.02,.18)))+
    labs(title=paste0('OCCAMS neoadjuvant | Biopsy versus surgery specimens | ',kind),
      x=NULL,y='Deconvoluted malignant GSVA score',fill='RNA specimen')+
    theme_classic(base_size=16)+theme(plot.title=element_text(face='bold',size=18),
      axis.text.x=element_text(angle=45,hjust=1,vjust=1,colour='black'),
      axis.line.x=element_blank(),legend.position='top',legend.title=element_text(face='bold'),
      plot.margin=margin(12,16,12,12))
  if(kind=='States') p <- p+facet_wrap(~family,ncol=1,scales='free_y')
  print(p)
}
dev.off()
figure_index <- lapply(figure_index,function(x)x[x$path!=specimen_path,,drop=FALSE])
figure_index[[length(figure_index)+1L]] <- data.frame(role='Auto_supporting',path=specimen_path,
  page=1:2,cohort='occams_neo',variable='rna_specimen_source',feature_type=c('MP','State'),features=c(17L,10L),n=157L)
####################
figure_index <- bind_rows(figure_index)
write_table(figure_index,'figure_index')
stopifnot(all(file.info(file.path(out,unique(figure_index$path)))$size>1000))
summary <- objects$statistics %>% group_by(cohort,family) %>% summarise(tests=n(),
  estimated=sum(status=='estimated'),not_estimated=sum(status!='estimated'),
  nominal=sum(p<.05,na.rm=TRUE),FDR_family=sum(q_family<.05,na.rm=TRUE),FDR_32=sum(q<.05,na.rm=TRUE),.groups='drop')
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE)
summary_path <- file.path(SCREF_SUMMARY_DIR,'malignant_clinical_associations.csv')
write.csv(summary,summary_path,row.names=FALSE)
stopifnot(file.info(summary_path)$size<100000)
validation <- data.frame(check=c('Exact historical cohort membership','Parent score coverage',
  'All 32 features per eligible cohort/variable','All indexed figures exist'),passed=TRUE,
  detail=c(paste(names(cohort_labels),vapply(names(cohort_labels),function(key)
    sum(objects$membership$cohort==key),integer(1)),sep='=',collapse='; '),
    '281 OCCAMS; 88 TCGA; current malignant scores and fractions',
    paste(nrow(objects$statistics),'feature comparisons'),paste(length(unique(figure_index$path)),'PDFs;',nrow(figure_index),'pages')))
write.csv(validation,file.path(out,'reports/validation.csv'),row.names=FALSE)
files <- unique(figure_index$path)
writeLines(c('# Malignant-compartment clinical associations','',
  'Five historical OCCAMS selections and the TCGA full surgery pathway/surgery-Yes selections use exact historical patient membership. Clinical analyses do not exclude patients for missing survival endpoints.',
  'The September 25 OCCAMS and September 18 TCGA workbooks supply current curated fields. OCCAMS comorbidity and the explicitly listed TCGA cBioPortal fields require separate sources; the source audit lists each field. Missing curated fields are never filled from older values.',
  'Scores are reused from the verified malignant-survival cache: Gaussian GSVA on common-gene log2(CPM+1), with fixed parent populations of 281 and 88. No subgroup rescaling or gene-set rescoring is performed. Five state fractions are within-malignant RNA contributions, not measured cell proportions.',
  'Age is split at <=50 versus >50 in both datasets. OCCAMS final Siewert is shown as separate I/II/III classes. The combined anatomical comparison is Oesophagus + GOJ I versus GOJ II/III, using Site_Final and Siewert_Final; GOJ cases without a class remain missing. Final Siewert classes are recorded only for GOJ cases in the current workbook; no Siewert class is inferred for oesophageal location records.',
  'The 32 features comprise 17 MPs, five state MP unions, five DGE signatures and five fractions. Two groups use an unpaired Wilcoxon test; more groups use Kruskal-Wallis. Clinical categories with fewer than three patients are excluded from both plots and tests, as requested. At least two retained categories are required. All retained categories enter the rank test. Constant scores remain unestimated. No test silently drops a displayed category. Unknown/conflicting categories remain missing. Mandard is shown both as TRG1-5 and TRG1-3 versus TRG4-5; only grades meeting the three-patient threshold are displayed; original counts and exclusions remain in the inventory.',
  'Figure stars use q_family, the original BH correction within cohort, clinical variable and signature family (17 MPs or five state measurements). One, two and three stars indicate q_family below 0.05, 0.01 and 0.001. No star means either the threshold is not met or the comparison is not estimated; the statistics table distinguishes these cases. Exact P, q_family and the broader 32-feature q remain in tables. Neither correction controls searching across all variables or overlapping cohorts; these are exploratory associations.',
  'Each cohort has one MP PDF and one State PDF containing every eligible clinical variable, with one variable per page. The two surgery cohorts are in Auto_final and the five remaining cohorts are in Auto_supporting; demographic and other clinical variables are not separated. All 17 MPs retain full descriptions on the shared x axis. State pages contain MP-union, DGE and RNA-fraction rows. Colours identify clinical groups, points are small, and only FDR stars are annotated. There are no explanatory subtitles, footnotes or exact P/q labels.',
  'The dedicated neoadjuvant biopsy-versus-surgery PDF contains two pages: 17 malignant MPs, then MP-union and DGE state scores. It compares 77 biopsy with 80 surgery patients using the same deconvoluted scores and family-adjusted tests as the RNA specimen source comparison; no subtitles are shown.',
  'Treatment response, recurrence and postoperative variables describe their recorded timepoint; these associations are not baseline prediction tests or causal treatment effects. Deconvolution uncertainty is not propagated into rank-test P values. Sparse inferred states remain uncertain.',
  '',paste0('- [',basename(files),'](../',files,')'),'',
  '[Figure/page index](../tables/figure_index.csv) | [Clinical eligibility](../tables/variable_inventory.csv)',
  '[Full statistics](../tables/association_statistics.csv) | [Metadata sources](../tables/metadata_source_audit.csv)',
  '[Validation](validation.csv) | [Methodology](../../../analysis/methodology/clinical/malignant_clinical_associations_methodology.md)'),
  file.path(out,'reports/current_results.md'))
writeLines(c('# Malignant clinical associations','',
  'Start with [current results and figures](reports/current_results.md).'),file.path(out,'README.md'))
message('Complete: ',out,'; ',nrow(figure_index),' pages')
