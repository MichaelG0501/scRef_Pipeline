####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/Auto_occams_bulk_clinical_associations.R
#   Description: Tests current centred refined OCCAMS bulk MP/state GSVA scores
#     across corrected clinical strata in five selected OCCAMS cohorts.
#   Methodology: analysis/methodology/clinical/Auto_occams_bulk_clinical_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
# Inputs:
#   - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_mp_scores.rds
#   - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_state_scores.rds
#   - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_dge_scores.rds
#   - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv
#     (only comorbidity remains authoritative here; all covered fields use workbook)
#   - analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx
#     sheet Patient_Master_388 (explicit latest corrected source, not mtime inference)
#   - ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
# Outputs:
#   intermediate/:
#     - ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_clinical_association_plot_data.rds
#   tables/:
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_variable_inventory.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_stats.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_subject_metadata.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_metadata_audit.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_cohort_membership.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_pdf_inventory.csv
#     - ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_clinical_association_input_manifest.csv
#   figures/:
#     - ref_outs/OCCAMS/clinical/figures/selected_cohorts/Auto_OCCAMS_<cohort>_clinical_scores_<MP|State>.pdf
#       cohorts: occams_combined, occams_neo, occams_neo_biopsy,
#       occams_neo_surgery, occams_surgery (five PDFs per type; one variable/page).
#       State pages contain both MP-union and top-20-DGE signatures in separate rows.
#   Scheduler logs: <project_root>/<job_id>.OU and .ER.
#   Previous owned outputs: ref_outs/OCCAMS/clinical/archive/clinical_associations/before_<run_id>/
#   reports/ and summaries/:
#     - ref_outs/OCCAMS/clinical/reports/Auto_OCCAMS_clinical_association_run_summary.txt
#     - updates/new_updates/summaries/Auto_OCCAMS_clinical_association_summary.csv
# Cache/replot behavior:
#   SCREF_FORCE_REBUILD=TRUE rebuilds persistent plot/statistics data.
#   SCREF_REPLOT_ONLY=TRUE redraws ten PDFs from a matching persistent plot-data RDS.
#   Caches validate hashes of all inputs/script, package versions, and settings;
#   an input change automatically invalidates the cache. Both flags together fail.
# Run command:
#   /opt/pbs/bin/qsub job_scripts/clinical/Auto_occams_bulk_clinical_associations.sh
#   Allocated compute node: Rscript analysis/clinical/Auto_occams_bulk_clinical_associations.R
# Conda env: dmtcp; PBS resources: 2 CPUs, 12 GB, 01:00:00.
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(readxl)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

out_root <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "clinical")
out_dirs <- setNames(file.path(out_root, c("intermediate", "tables", "figures", "logs", "reports")),
                     c("intermediate", "tables", "figures", "logs", "reports"))
invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(SCREF_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

mp_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_mp_scores.rds")
state_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_state_scores.rds")
metadata_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_qc_pass_subject_metadata.csv")
grouping_path <- SCREF_MP_GROUPING_CSV
plot_data_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_clinical_association_plot_data.rds")
inventory_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_variable_inventory.csv")
stats_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_stats.csv")
####################
figure_dir <- file.path(out_dirs[["figures"]], "selected_cohorts")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
workbook_path <- file.path(SCREF_PROJECT_DIR, "analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx")
qc_path <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv")
dge_score_path <- file.path(out_dirs[["intermediate"]], "Auto_OCCAMS_qc_pass_dge_scores.rds")
resolved_metadata_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_subject_metadata.csv")
metadata_audit_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_metadata_audit.csv")
membership_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_cohort_membership.csv")
pdf_inventory_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_pdf_inventory.csv")
manifest_path <- file.path(out_dirs[["tables"]], "Auto_OCCAMS_clinical_association_input_manifest.csv")
cohort_labels <- c(occams_combined = "Neoadjuvant + surgery pathways", occams_neo = "Neoadjuvant: all typed specimens",
  occams_neo_biopsy = "Neoadjuvant: biopsy specimens", occams_neo_surgery = "Neoadjuvant: surgery specimens",
  occams_surgery = "Surgery pathway (RNA matched 68/72)")
expected_cohort_n <- c(225L,157L,77L,80L,68L)
run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
####################
report_path <- file.path(out_dirs[["reports"]], "Auto_OCCAMS_clinical_association_run_summary.txt")
summary_path <- file.path(SCREF_SUMMARY_DIR, "Auto_OCCAMS_clinical_association_summary.csv")

force_rebuild <- identical(toupper(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")), "TRUE")
replot_only <- identical(toupper(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")), "TRUE")
minimum_level_n <- 10L
if (force_rebuild && replot_only) stop("FORCE_REBUILD and REPLOT_ONLY cannot both be TRUE")

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
    "Gender", "Age at diagnosis (>60)", "Diagnostic combined tumour site",
    "Final pretreatment combined tumour site", "Resection pathology location", "Diagnostic histology",
    "Diagnostic Siewert classification", "Final pretreatment Siewert classification",
    "Resection pathology Siewert classification", "Final pretreatment T stage",
    "Final pretreatment N stage (TNM7)", "Final pretreatment M stage", "Treatment-plan performance status",
    "Treatment-plan intent", "Treatment-plan curative modality", "Main surgery",
    "Treatment-plan comorbidity burden", "Resection pathology T stage",
    "Resection pathology N stage (TNM7)", "Resection pathology M stage",
    "Resection positive-node burden", "Resection pathology grade", "Resection signet-ring cells",
    "Resection adjacent Barrett-associated dysplasia", "Recorded neoadjuvant therapy at resection",
    "Operation type", "Mandard tumour regression grade", "Mandard response grouping",
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
    if(variable=="treatment_plan_comorbidity_group") m[[variable]] <- clean_curated(older[[variable]])
    else {
      if(!column %in% names(w)) stop("Missing curated column: ",column)
      m[[variable]] <- consensus_tokens(clean_curated(w[[column]]))
    }
  }
  age <- suppressWarnings(as.numeric(m$age_group))
  m$age_group <- ifelse(is.finite(age),ifelse(age>60,">60","<=60"),NA_character_)
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
  m$pathology_mandard_response_group <- ifelse(trg %in% c("TRG1","TRG2"),"TRG1-2",
    ifelse(trg %in% c("TRG3","TRG4","TRG5"),"TRG3-5",NA_character_))
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
  audit <- bind_rows(lapply(variable_config$variable,function(v) {
    old <- if(v %in% names(older)) clean_curated(older[[v]]) else rep(NA_character_,nrow(m))
    value <- m[[v]]
    data.frame(variable=v,n_current_available=sum(!is.na(value)),n_current_missing=sum(is.na(value)),
      n_older_available=sum(!is.na(old)),n_changed=sum(ifelse(is.na(old),"<missing>",old)!=ifelse(is.na(value),"<missing>",value)))
  })) %>% left_join(variable_config,by="variable")
  list(metadata=m,cohorts=cohorts,audit=audit)
}
####################

feature_long <- function(score_matrix, feature_type) {
  x <- as.data.frame(score_matrix, check.names = FALSE)
  x$subject_id <- rownames(x)
  x %>% pivot_longer(-subject_id, names_to = "feature", values_to = "score") %>%
    mutate(feature_type = feature_type)
}

prepare_variable_data <- function(metadata, scores_long, variable_name) {
  x <- metadata %>%
    select(subject_id, group = all_of(variable_name)) %>%
    mutate(group = trimws(as.character(group)), group = na_if(group, "")) %>%
    filter(!is.na(group))
  level_counts <- x %>% count(group, name = "level_n")
  retained <- level_counts %>% filter(level_n >= minimum_level_n) %>% pull(group)
  x <- x %>% filter(group %in% retained)
  if (n_distinct(x$group) < 2) return(NULL)
  x %>% inner_join(scores_long, by = "subject_id") %>%
    mutate(clinical_variable = variable_name)
}

compute_stats <- function(plot_data) {
  plot_data %>%
    group_by(cohort, clinical_variable, feature_type, feature) %>%
    group_modify(~ {
      d <- .x %>% filter(!is.na(group), !is.na(score))
      groups <- sort(unique(d$group))
      group_summary <- d %>% group_by(group) %>%
        summarise(n = n_distinct(subject_id), median = median(score), mean = mean(score), .groups = "drop")
      test_name <- if (length(groups) == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"
      test <- if (length(groups) == 2) wilcox.test(score ~ group, data = d, exact = FALSE) else kruskal.test(score ~ group, data = d)
      if(!is.finite(test$p.value)) stop("Nonfinite clinical association test")
      data.frame(
        test = test_name, n_groups = length(groups), n_subjects = n_distinct(d$subject_id),
        p_value = test$p.value,
        group_summary = paste0(group_summary$group, " (n=", group_summary$n,
                               ", median=", sprintf("%.3f", group_summary$median), ")", collapse = " | "),
        stringsAsFactors = FALSE
      )
    }) %>% ungroup() %>%
    group_by(cohort, clinical_variable, feature_type) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>% ungroup()
}

build_plot_objects <- function() {
  required <- c(mp_score_path, state_score_path, dge_score_path, metadata_path, grouping_path, workbook_path, qc_path)
  if (any(!file.exists(required))) stop("Missing upstream OCCAMS survival output(s): ", paste(required[!file.exists(required)], collapse = ", "))
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_scores <- readRDS(mp_score_path)
  state_scores <- readRDS(state_score_path)
  dge_scores <- readRDS(dge_score_path)
  if (nrow(metadata) != 281 || nrow(mp_scores) != 281 || nrow(state_scores) != 281) {
    stop("Expected exactly 281 QC-pass subjects in metadata, MP scores, and state scores")
  }
  if (!setequal(metadata$subject_id, rownames(mp_scores)) || !setequal(metadata$subject_id, rownames(state_scores))) {
    stop("Subject identities differ between metadata and score matrices")
  }
  qc <- read.csv(qc_path,stringsAsFactors=FALSE)
  qc_ids <- qc$Subject_ID[qc$Annotation_Compatible_Below_50_Flag=="no"]
  if(length(qc_ids)!=281L || anyDuplicated(qc_ids) || !setequal(qc_ids,metadata$subject_id) ||
     !setequal(qc_ids,rownames(dge_scores)) || anyDuplicated(metadata$subject_id)) stop("QC-pass identity mismatch")
  resolved <- resolve_metadata(metadata,qc_ids)
  metadata <- resolved$metadata

  grouping <- read.csv(grouping_path, stringsAsFactors = FALSE, check.names = FALSE)
  mp_labels <- setNames(grouping$plot_label, grouping$mp)
  if(ncol(mp_scores)!=17L || !setequal(colnames(mp_scores),grouping$mp) || ncol(state_scores)!=5L || ncol(dge_scores)!=5L ||
     !setequal(colnames(state_scores),SCREF_PRIMARY_STATE_ORDER) || !setequal(colnames(dge_scores),paste0(SCREF_PRIMARY_STATE_ORDER," (DGE)")) ||
     any(!is.finite(mp_scores)) || any(!is.finite(state_scores)) || any(!is.finite(dge_scores))) stop("Invalid 17-MP/five-state score contract")
  score_long <- bind_rows(feature_long(mp_scores, "MP"), feature_long(state_scores, "State MP union"),feature_long(dge_scores,"State DGE"))
  plot_data <- bind_rows(lapply(names(resolved$cohorts),function(key) {
    m <- metadata[resolved$cohorts[[key]],,drop=FALSE]
    bind_rows(lapply(variable_config$variable, function(variable_name) prepare_variable_data(m, score_long, variable_name))) %>% mutate(cohort=key)
  })) %>%
    left_join(variable_config, by = c("clinical_variable" = "variable")) %>%
    mutate(feature_label = ifelse(feature_type == "MP", unname(mp_labels[feature]), feature))

  inventory <- bind_rows(lapply(names(resolved$cohorts),function(key) {
    cohort_metadata <- metadata[resolved$cohorts[[key]],,drop=FALSE]
    tested_variables <- unique(plot_data$clinical_variable[plot_data$cohort==key])
    bind_rows(lapply(variable_config$variable, function(variable_name) {
    values <- cohort_metadata[[variable_name]]
    values <- trimws(as.character(values))
    values[values == ""] <- NA_character_
    counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
    data.frame(
      cohort=key,cohort_n=nrow(cohort_metadata),variable = variable_name,
      n_available = sum(!is.na(values)), n_missing = sum(is.na(values)),
      n_levels_observed = length(counts), n_levels_at_least_10 = sum(counts >= minimum_level_n),
      levels_and_counts = paste0(names(counts), "=", as.integer(counts), collapse = " | "),
      included_in_tests = variable_name %in% tested_variables,
      exclusion_rule = ifelse(variable_name %in% tested_variables, "", "fewer than two levels with n>=10"),
      stringsAsFactors = FALSE
    )
    }))
  })) %>% left_join(variable_config, by = "variable")

  stats <- compute_stats(plot_data)

  list(plot_data = plot_data, stats = stats,
       inventory = inventory, mp_order = grouping$mp,metadata=metadata,metadata_audit=resolved$audit,
       membership=bind_rows(lapply(names(resolved$cohorts),function(key) data.frame(cohort=key,subject_id=metadata$subject_id[resolved$cohorts[[key]]]))),
       state_order = c(intersect(SCREF_PRIMARY_STATE_ORDER, colnames(state_scores)),paste0(SCREF_PRIMARY_STATE_ORDER," (DGE)")))
}

plot_boxplots <- function(objects) {
  pdf_inventory <- list()
  for(cohort_key in names(cohort_labels)) {
  for (type_name in c("MP", "State")) {
    out_path <- file.path(figure_dir,paste0("Auto_OCCAMS_",cohort_key,"_clinical_scores_",type_name,".pdf"))
    grDevices::cairo_pdf(out_path, width = 18, height = if(type_name=="MP") 10 else 12, onefile = TRUE)
    page_count <- 0L
    for (variable_name in unique(objects$plot_data$clinical_variable[objects$plot_data$cohort==cohort_key])) {
      variable_data <- objects$plot_data %>% filter(cohort==cohort_key,clinical_variable == variable_name)
      d <- variable_data %>% filter(if(type_name=="MP") feature_type=="MP" else feature_type!="MP")
      if (nrow(d) == 0) next
      variable_label <- unique(variable_data$label)
      
      feature_order <- if (type_name == "MP") objects$mp_order else objects$state_order
      feature_order <- feature_order[feature_order %in% d$feature]
      d$feature_label <- factor(d$feature_label,levels=unique(d$feature_label[match(feature_order,d$feature)]))
      if(type_name=="State") d$feature_label <- factor(sub(" \\(DGE\\)$","",as.character(d$feature_label)),levels=SCREF_PRIMARY_STATE_ORDER)
      d$feature_type <- factor(d$feature_type,levels=c("MP","State MP union","State DGE"))
      d$group <- as.factor(d$group)
      
      legend_counts <- d %>% distinct(subject_id, group) %>% count(group, name = "n_samples")
      legend_labels <- setNames(
        paste0(legend_counts$group, " (n=", legend_counts$n_samples, ")"),
        as.character(legend_counts$group)
      )
      
      group_levels <- levels(d$group)
      palette <- scales::hue_pal()(length(group_levels))
      names(palette) <- group_levels
      
      stats_subset <- objects$stats %>% filter(cohort==cohort_key,clinical_variable == variable_name)
      annot_df <- d %>%
        group_by(feature, feature_label,feature_type) %>%
        summarise(y_pos = max(score, na.rm = TRUE) + max(0.08 * diff(range(score, na.rm=TRUE)), 0.02), .groups = "drop") %>%
        left_join(stats_subset, by = c("feature", "feature_type")) %>%
        mutate(sig_label = case_when(
          is.na(p_adj_bh) ~ "",
          p_adj_bh < 0.001 ~ "***",
          p_adj_bh < 0.01 ~ "**",
          p_adj_bh < 0.05 ~ "*",
          TRUE ~ "ns"
        )) %>%
        mutate(sig_label = ifelse(p_adj_bh < 0.05, sig_label, ""))
      annot_df$feature_type <- factor(annot_df$feature_type,levels=levels(d$feature_type))
      
      p <- ggplot(d, aes(x = feature_label, y = score, fill = group, color = group)) +
        geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, linewidth = 0.4, color = "black") +
        geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,seed=1701), alpha = 0.7, size = 1.0, stroke = 0, show.legend = FALSE) +
        geom_text(data = annot_df %>% filter(sig_label != ""), aes(x = feature_label, y = y_pos, label = sig_label), inherit.aes = FALSE, size = 4, fontface = "bold") +
        scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
        scale_color_manual(values = palette, guide = "none", drop = FALSE) +
        labs(title = paste0("OCCAMS | ",cohort_labels[[cohort_key]]," | ",variable_label, " - ", type_name, " scores"),
             subtitle = "Subject-level GSVA; levels require n>=10. Stars: BH FDR <0.05, <0.01, <0.001 within cohort/variable/signature family.",
             x = NULL, y = "GSVA score", fill = "Clinical group") +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 10, colour = "grey35"),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.line.x = element_blank(),
              legend.position = "top",
              legend.title = element_text(face = "bold"),
              plot.margin = margin(12, 16, 12, 12))
      if(type_name=="State") p <- p+facet_wrap(~feature_type,ncol=1,scales="free_y")
      print(p)
      page_count <- page_count+1L
    }
    dev.off()
    if(page_count<1L || file.info(out_path)$size<1000) stop("Missing/empty clinical-score PDF: ",out_path)
    pdf_inventory[[paste(cohort_key,type_name)]] <- data.frame(cohort=cohort_key,feature_type=type_name,
      cohort_n=expected_cohort_n[match(cohort_key,names(cohort_labels))],pages=page_count,path=out_path)
  }
  }
  bind_rows(pdf_inventory)
}



####################
required_inputs <- c(mp_score_path,state_score_path,dge_score_path,metadata_path,grouping_path,workbook_path,qc_path,
  file.path(SCREF_PROJECT_DIR,"analysis/clinical/Auto_occams_bulk_clinical_associations.R"))
if(any(!file.exists(required_inputs))) stop("Missing required input: ",paste(required_inputs[!file.exists(required_inputs)],collapse=", "))
manifest <- data.frame(path=required_inputs,md5=unname(tools::md5sum(required_inputs)))
cache_contract <- list(version="corrected_master_five_cohorts_v1",inputs=manifest,variables=variable_config,
  minimum_level_n=minimum_level_n,cohort_n=expected_cohort_n,
  packages=vapply(c("dplyr","tidyr","ggplot2","scales","readxl"),function(p) as.character(packageVersion(p)),character(1)))
objects <- if(file.exists(plot_data_path) && !force_rebuild) readRDS(plot_data_path) else NULL
cache_valid <- !is.null(objects) && identical(objects$cache_contract,cache_contract)
if(replot_only && !cache_valid) stop("REPLOT_ONLY requires a matching corrected-metadata five-cohort cache")
if(!cache_valid) objects <- build_plot_objects()
# Archive only outputs owned by this script; upstream scores/metadata and other
# OCCAMS clinical figures are never moved or overwritten.
archive_dir <- file.path(out_root,"archive/clinical_associations",paste0("before_",run_id))
owned_paths <- c(plot_data_path,inventory_path,stats_path,report_path,summary_path,
  resolved_metadata_path,metadata_audit_path,membership_path,pdf_inventory_path,manifest_path,
  file.path(out_dirs[["figures"]],c("Auto_OCCAMS_clinical_score_boxplots_MP.pdf","Auto_OCCAMS_clinical_score_boxplots_State.pdf")),figure_dir)
for(path in owned_paths) {
  if(file.exists(path) && (!dir.exists(path) || length(list.files(path)))) {
    parent <- if(startsWith(path,SCREF_SUMMARY_DIR)) "compact_summaries" else basename(dirname(path))
    destination <- file.path(archive_dir,parent,basename(path))
    dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
    if(file.exists(destination) || !file.rename(path,destination)) stop("Could not preserve previous output: ",path)
  }
}
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
objects$cache_contract <- cache_contract
saveRDS(objects,plot_data_path,compress="xz")
write.csv(objects$inventory,inventory_path,row.names=FALSE,na="")
write.csv(objects$stats,stats_path,row.names=FALSE,na="")
write.csv(objects$metadata,resolved_metadata_path,row.names=FALSE,na="")
write.csv(objects$metadata_audit,metadata_audit_path,row.names=FALSE,na="")
write.csv(objects$membership,membership_path,row.names=FALSE,na="")
write.csv(manifest,manifest_path,row.names=FALSE)
pdf_inventory <- plot_boxplots(objects)
if(nrow(pdf_inventory)!=10L || any(table(pdf_inventory$feature_type)!=5L)) stop("Expected five MP and five State PDFs")
write.csv(pdf_inventory,pdf_inventory_path,row.names=FALSE)
####################

summary_df <- objects$stats %>% group_by(cohort, feature_type) %>%
  summarise(qc_pass_subjects = 281L, n_clinical_variables = n_distinct(clinical_variable),
            n_tests = n(), n_nominal_p_lt_0_05 = sum(p_value < 0.05, na.rm = TRUE),
            n_bh_fdr_lt_0_05 = sum(p_adj_bh < 0.05, na.rm = TRUE), .groups = "drop")
write.csv(summary_df, summary_path, row.names = FALSE)
if(file.info(summary_path)$size>=100000) stop("Compact summary exceeds 100 KB")
writeLines(c(
  "OCCAMS current-centred bulk clinical association run",
  paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "qc_pass_subjects=281",
  paste0("metadata_workbook=",workbook_path,"; sheet=Patient_Master_388"),
  "cohort_sizes=combined 225; neoadjuvant 157; neoadjuvant biopsy 77; neoadjuvant surgery specimens 80; surgery pathway 68",
  "scores=shared 281-subject Gaussian GSVA; no subgroup-specific rescoring; MP-union and DGE state signatures included",
  "legacy_metadata_use=comorbidity only; absent from corrected workbook; no missing corrected values filled from older metadata",
  paste0("pdf_count=",nrow(pdf_inventory),"; output_folder=",figure_dir),
  paste0("previous_outputs_preserved=",archive_dir),
  paste0("minimum_subjects_per_displayed_level=", minimum_level_n),
  paste0("clinical_variables_tested=", n_distinct(objects$plot_data$clinical_variable)),
  paste0("association_tests=", nrow(objects$stats)),
  "result=PASS"
), report_path)
message("Completed OCCAMS QC-pass clinical association workflow.")
