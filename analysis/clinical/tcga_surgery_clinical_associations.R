####################
# Analysis registry
# Status: active
# Script: analysis/clinical/tcga_surgery_clinical_associations.R
# Description: Selected TCGA EAC surgery-Yes clinical associations for 17 MPs and
#   five MP-union / five DGE state signatures. Patient is the statistical unit.
# Methodology: analysis/methodology/clinical/tcga_surgery_clinical_associations_methodology.md
# Inputs relative to project root:
#   analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx (EAC_RNA88_Final)
#   ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#   ref_outs/clinical_selected_survival/tables/{all_feature_gsva_scores,cohort_membership,input_manifest}.csv
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv
#   analysis/shared/scRef_config.R
# External import: /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/clinical_cbioportal.tsv
#   (copied into the canonical reference/ tier; contributes fields absent from the corrected workbook).
# Stable outputs: ref_outs/TCGA/surgery_clinical_associations/
#   reference/clinical_cbioportal.tsv
#   intermediate/clinical_associations.rds
#   tables/{subject_metadata,variable_inventory,association_statistics,plot_source_data,
#     selected_scores,metadata_source_audit,input_manifest,figure_index}.csv
#   figures/Auto_final/{01_MP_priority_clinical,02_State_priority_clinical}.pdf
#   figures/Auto_supporting/{01_MP_other_clinical,02_State_other_clinical}.pdf
#   reports/{current_results.md,run_summary.txt}
#   updates/new_updates/summaries/tcga_surgery_clinical_associations.csv (<100 KB)
# Cache: input/script/package hashes invalidate changed inputs; FORCE_REBUILD rebuilds
#   metadata/tests only. SCREF_REPLOT_ONLY=TRUE requires a matching complete cache.
# Scores reuse the selected-survival 86-patient Gaussian raw-TPM GSVA parent; no rescoring.
# Run: /opt/pbs/bin/qsub job_scripts/clinical/tcga_surgery_clinical_associations.sh
# Direct allocated-node command: Rscript analysis/clinical/tcga_surgery_clinical_associations.R
# Environment: dmtcp; PBS 2 CPUs, 8 GB, 01:00:00, -koed; scheduler logs in project root.
####################
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
source('/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R')
if(!nzchar(Sys.getenv('PBS_JOBID'))) stop('Run on a PBS compute node')
setwd(SCREF_PROJECT_DIR)
out <- file.path(SCREF_REF_OUTS_DIR,'TCGA/surgery_clinical_associations')
for(d in c('reference','intermediate','tables','figures/Auto_final','figures/Auto_supporting','reports'))
  dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
external <- '/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/clinical_cbioportal.tsv'
snapshot <- file.path(out,'reference/clinical_cbioportal.tsv')
if(!file.exists(external)) stop('Required cBioPortal source is absent')
if(!file.exists(snapshot) || tools::md5sum(snapshot)!=tools::md5sum(external)) {
  if(file.exists(snapshot)) {
    archive <- file.path(SCREF_PROJECT_DIR,'Auto_session_checks/tcga_surgery_clinical/Auto_archive',
      paste0('clinical_cbioportal_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.tsv'))
    dir.create(dirname(archive),recursive=TRUE,showWarnings=FALSE)
    if(file.exists(archive) || !file.rename(snapshot,archive)) stop('Cannot preserve old clinical snapshot')
  }
  stopifnot(file.copy(external,snapshot,overwrite=FALSE))
}
inputs <- c(workbook='analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx',
  reconstruction='ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds',
  scores='ref_outs/clinical_selected_survival/tables/all_feature_gsva_scores.csv',
  membership='ref_outs/clinical_selected_survival/tables/cohort_membership.csv',
  score_manifest='ref_outs/clinical_selected_survival/tables/input_manifest.csv',
  mp_genes=SCREF_MP_GENES_RDS,grouping=SCREF_MP_GROUPING_CSV,
  markers='ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv',
  cbio=snapshot,config='analysis/shared/scRef_config.R',script='analysis/clinical/tcga_surgery_clinical_associations.R')
stopifnot(all(file.exists(inputs)))
manifest <- data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
force <- toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))=='TRUE'
replot <- toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))=='TRUE'
if(force && replot) stop('FORCE_REBUILD and REPLOT_ONLY are mutually exclusive')
minimum_n <- 10L
contract <- list(version=1L,inputs=manifest,minimum_n=minimum_n,
  packages=vapply(c('readxl','dplyr','tidyr','ggplot2'),function(p) as.character(packageVersion(p)),character(1)))
clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c('','na','n/a','unknown','not reported','not available','not applicable','[not available]','gx','nx','mx','tx')] <- NA_character_
  x
}
yes_no <- function(x) { z<-tolower(clean(x)); ifelse(z %in% c('true','yes'),'Yes',ifelse(z %in% c('false','no'),'No',NA_character_)) }
stage <- function(x) {z<-clean(x); ifelse(z %in% as.character(1:4),z,NA_character_)}
tn <- function(x,letter) {z<-toupper(clean(x)); ifelse(grepl(paste0('^',letter,'[0-4]'),z),substr(z,1,2),NA_character_)}
spec <- data.frame(variable=c('sex','age_group','grade','pathological_stage','clinical_stage','pT','pN','cT','cN','pM','positive_nodes','location','alcohol','smoking','prior_malignancy','prior_treatment','race','ethnicity','preoperative_treatment'),
  label=c('Sex','Age at diagnosis (>60 years)','Histological grade','Pathological AJCC8 stage','Clinical AJCC8 stage','Pathological T stage','Pathological nodal involvement','Clinical T stage','Clinical N stage','Pathological M classification','Positive lymph-node burden','Tumour location','Recorded alcohol history','Smoking exposure (pack-years)','Prior malignancy','Recorded prior treatment','Recorded race','Recorded ethnicity','Recorded preoperative treatment'),
  source=c('cbio','cbio',rep('workbook',10),rep('cbio',6),'workbook'),
  field=c('Sex','Diagnosis Age','Grade_AJCC8_used','Pathological_Stage_AJCC8','Clinical_Stage_AJCC8','pT_AJCC8_used','pN_AJCC8_used','cT_AJCC8_used','cN_AJCC8_used','M_AJCC8_used_for_path','Positive_nodes_HE','Tumour_location','Alcohol History Documented','Person Cigarette Smoking History Pack Year Value','Prior Malignancy','Prior Treatment','Race Category','Ethnicity Category','Treatment_prior_to_surgery_raw'),
  role=c(rep('Auto_final',4),rep('Auto_supporting',15)),stringsAsFactors=FALSE)
build <- function() {
  w <- as.data.frame(read_excel(inputs[['workbook']],sheet='EAC_RNA88_Final'))
  stopifnot(nrow(w)==88L,!anyDuplicated(w$Patient_ID))
  w <- w[w$Flowchart_Pathway_Final=='Surgery group' & w$Flowchart_Surgery_Status_Final=='Yes',,drop=FALSE]
  stopifnot(nrow(w)==72L,!anyNA(w$Patient_ID),!anyDuplicated(w$RNA_primary_sample_ID))
  members <- read.csv(inputs[['membership']])
  stopifnot(setequal(w$Patient_ID,members$subject_id[members$cohort=='tcga_yes']))
  recon <- readRDS(inputs[['reconstruction']])
  stopifnot(!anyDuplicated(recon$sample_barcode))
  r <- recon[match(w$RNA_primary_sample_ID,recon$sample_barcode),,drop=FALSE]
  stopifnot(!anyNA(r$case_barcode),all(r$case_barcode==w$Patient_ID),all(r$HistologyGroup=='EAC'),all(r$sample_type_code=='01'))
  cb <- read.delim(snapshot,check.names=FALSE,stringsAsFactors=FALSE)
  # Clinical records identify patients explicitly. Require one primary-tumour
  # record per selected patient; preserve each source's distinct sample ID.
  cb <- cb[cb[['Study ID']]=='esca_tcga_gdc' & cb[['Sample Type']]=='Primary Tumor' & cb[['Patient ID']] %in% w$Patient_ID,,drop=FALSE]
  stopifnot(nrow(cb)==72L,!anyDuplicated(cb[['Patient ID']]),all(w$Patient_ID %in% cb[['Patient ID']]))
  cb <- cb[match(w$Patient_ID,cb[['Patient ID']]),,drop=FALSE]
  old_manifest <- read.csv(inputs[['score_manifest']])
  for(pair in list(c('tcga_workbook','workbook'),c('mp_genes','mp_genes'),c('mp_grouping','grouping'),c('state_markers','markers'))) {
    recorded <- old_manifest$md5[old_manifest$input==pair[1]]
    if(length(recorded)!=1L || recorded!=unname(tools::md5sum(inputs[[pair[2]]]))) stop('Stale upstream score provenance: ',pair[1])
  }
  grouping <- read.csv(inputs[['grouping']])
  features <- c(grouping$mp,paste0(SCREF_PRIMARY_STATE_ORDER,' (MP union)'),paste0(SCREF_PRIMARY_STATE_ORDER,' (DGE)'))
  family <- setNames(c(rep('MP',17),rep('State MP union',5),rep('State DGE',5)),features)
  scores <- read.csv(inputs[['scores']],check.names=FALSE)
  stopifnot(!anyDuplicated(scores$sample_id),all(features %in% names(scores)),length(grouping$mp)==17L)
  scores <- scores[match(w$RNA_primary_sample_ID,scores$sample_id),c('sample_id',features),drop=FALSE]
  stopifnot(identical(scores$sample_id,w$RNA_primary_sample_ID),all(is.finite(as.matrix(scores[,features]))))
  m <- data.frame(subject_id=w$Patient_ID,sample_id=w$RNA_primary_sample_ID,cbio_sample_id=cb[['Sample ID']])
  for(i in seq_len(nrow(spec))) m[[spec$variable[i]]] <- clean(if(spec$source[i]=='workbook') w[[spec$field[i]]] else cb[[spec$field[i]]])
  age <- suppressWarnings(as.numeric(m$age_group)); m$age_group <- ifelse(is.finite(age),ifelse(age>60,'>60','<=60'),NA_character_)
  m$grade <- ifelse(m$grade %in% c('G1','G2'),'G1/G2',ifelse(m$grade %in% c('G3','G4'),'G3/G4',NA_character_))
  for(v in c('clinical_stage','pathological_stage')) m[[v]] <- stage(m[[v]])
  for(v in c('pT','cT')) m[[v]] <- tn(m[[v]],'T')
  m$cN <- tn(m$cN,'N'); n <- tn(m$pN,'N'); m$pN <- ifelse(n=='N0','N0',ifelse(n %in% c('N1','N2','N3'),'N+',NA_character_))
  m$pM[!m$pM %in% c('M0','M1')] <- NA_character_
  n <- suppressWarnings(as.numeric(m$positive_nodes)); m$positive_nodes <- ifelse(is.finite(n)&n>=0,ifelse(n==0,'0',ifelse(n<=3,'1-3','>=4')),NA_character_)
  n <- suppressWarnings(as.numeric(m$smoking)); m$smoking <- ifelse(is.finite(n)&n>=0,ifelse(n==0,'0',ifelse(n<=20,'>0-20','>20')),NA_character_)
  for(v in c('alcohol','prior_malignancy','prior_treatment')) m[[v]] <- yes_no(m[[v]])
  # Independent source comparisons audit disagreements; curated staging always wins.
  audit <- bind_rows(lapply(list(c('sex','Gender'),c('age_group','Age_at_diagnosis'),c('pathological_stage','Stage_Simple'),c('grade','Grade')),function(pair) {
    data.frame(subject_id=m$subject_id,variable=pair[1],selected_value=as.character(m[[pair[1]]]),
      reconstructed_value=as.character(r[[pair[2]]]),note='Source comparison; recoding/edition differences are not silently coalesced')
  }))
  long <- scores %>% mutate(subject_id=m$subject_id) %>% pivot_longer(all_of(features),names_to='feature',values_to='score') %>%
    mutate(family=unname(family[feature]))
  inventory <- plots <- list()
  for(i in seq_len(nrow(spec))) {
    v <- spec$variable[i]; tab <- table(m[[v]]); retained <- names(tab)[tab>=minimum_n]
    eligible <- length(retained)>=2L
    inventory[[v]] <- cbind(spec[i,],data.frame(n_available=sum(!is.na(m[[v]])),n_missing=sum(is.na(m[[v]])),
      levels_and_counts=paste(names(tab),as.integer(tab),sep='=',collapse=' | '),retained_levels=paste(retained,collapse=' | '),
      n_retained=sum(tab[retained]),eligible=eligible,reason=if(eligible) '' else 'Fewer than two levels with at least 10 subjects'))
    if(eligible) plots[[v]] <- long %>% inner_join(data.frame(subject_id=m$subject_id,group=m[[v]]),by='subject_id') %>%
      filter(group %in% retained) %>% mutate(variable=v)
  }
  plot_data <- bind_rows(plots)
  statistics <- plot_data %>% group_by(variable,family,feature) %>% group_modify(~{
    d <- .x; k <- n_distinct(d$group)
    test <- if(k==2L) wilcox.test(score~group,data=d,exact=FALSE) else kruskal.test(score~group,data=d)
    if(!is.finite(test$p.value)) stop('Nonfinite association test')
    data.frame(n=n_distinct(d$subject_id),groups=k,test=if(k==2L) 'Wilcoxon rank-sum' else 'Kruskal-Wallis',p=test$p.value,
      group_summary=paste(vapply(split(d,d$group),function(x) paste0(x$group[1],': n=',nrow(x),', median=',signif(median(x$score),4)),character(1)),collapse=' | '))
  }) %>% ungroup() %>% group_by(variable,family) %>% mutate(FDR=p.adjust(p,'BH')) %>% ungroup()
  list(metadata=m,scores=scores,inventory=bind_rows(inventory),statistics=statistics,plot_data=plot_data,audit=audit,
    features=features,family=family,contract=contract)
}
cache_path <- file.path(out,'intermediate/clinical_associations.rds')
objects <- if(file.exists(cache_path)&&!force) readRDS(cache_path) else NULL
valid <- !is.null(objects)&&identical(objects$contract,contract)
if(replot&&!valid) stop('No matching complete replot cache')
if(!valid) {objects<-build();saveRDS(objects,cache_path,compress='xz')}
write_table <- function(x,name) write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
write_table(manifest,'input_manifest');write_table(objects$metadata,'subject_metadata');write_table(objects$scores,'selected_scores')
write_table(objects$inventory,'variable_inventory');write_table(objects$statistics,'association_statistics')
write_table(objects$plot_data,'plot_source_data');write_table(objects$audit,'metadata_source_audit')
index <- list()
for(role in c('Auto_final','Auto_supporting')) for(type in c('MP','State')) {
  filename <- paste0(if(type=='MP') '01_MP_' else '02_State_',if(role=='Auto_final') 'priority_clinical.pdf' else 'other_clinical.pdf')
  variables <- objects$inventory$variable[objects$inventory$eligible & objects$inventory$role==role]
  if(!length(variables)) stop('No eligible variables for ',role)
  path <- file.path(out,'figures',role,filename)
  cairo_pdf(path,width=if(type=='MP') 18 else 16,height=10,onefile=TRUE)
  page <- 0L
  for(v in variables) {
    groups <- if(type=='MP') list(objects$features[objects$family=='MP']) else list(objects$features[objects$family!='MP'])
    for(features in groups) {
      d <- objects$plot_data %>% filter(variable==v,feature %in% features)
      labels <- if(type=='MP') setNames(vapply(features,function(f) paste(c(f,strwrap(SCREF_MP_DESCRIPTIONS[[f]],width=18)),collapse='\n'),character(1)),features) else
        setNames(sub(' \\((MP union|DGE)\\)$','',features),features)
      d$x <- unname(labels[d$feature]); ordering <- unique(unname(labels[features]))
      d$x <- factor(d$x,levels=ordering);d$group <- factor(d$group);d$family <- factor(d$family,levels=c('MP','State MP union','State DGE'))
      ann <- d %>% group_by(feature,x,family) %>% summarise(y=max(score)+.10*diff(range(score))+.03,.groups='drop') %>%
        left_join(objects$statistics %>% filter(variable==v),by=c('feature','family')) %>%
        mutate(star=ifelse(FDR<.001,'***',ifelse(FDR<.01,'**',ifelse(FDR<.05,'*',''))))
      ann$family <- factor(ann$family,levels=levels(d$family))
      counts <- d %>% distinct(subject_id,group) %>% count(group)
      palette <- setNames(scales::hue_pal()(nrow(counts)),as.character(counts$group))
      legend <- setNames(paste0(counts$group,' (n=',counts$n,')'),as.character(counts$group))
      title <- spec$label[match(v,spec$variable)]
      p <- ggplot(d,aes(x,score,fill=group,colour=group))+
        geom_boxplot(position=position_dodge(.75),width=.6,outlier.shape=NA,alpha=.8,linewidth=.5,colour='black')+
        geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75,seed=1701),size=1.3,alpha=.65,stroke=0,show.legend=FALSE)+
        geom_text(data=ann[ann$star!='',],aes(x=x,y=y,label=star),inherit.aes=FALSE,size=6,fontface='bold')+
        scale_fill_manual(values=palette,labels=legend)+scale_colour_manual(values=palette,guide='none')+
        scale_y_continuous(expand=expansion(mult=c(.04,.16)))+
        labs(title=paste('TCGA EAC surgery Yes |',title),subtitle=paste0('Cohort N=72; displayed n=',n_distinct(d$subject_id),'. Stars: BH FDR <0.05 / <0.01 / <0.001 within signature family.'),
          x=NULL,y='GSVA enrichment score',fill=NULL)+theme_classic(base_size=19)+
        theme(plot.title=element_text(size=23,face='bold'),plot.subtitle=element_text(size=14),axis.text.x=element_text(size=14),
          legend.position='top',legend.text=element_text(size=17),strip.text=element_text(size=17,face='bold'),plot.margin=margin(14,16,14,16))
      if(type=='MP') p <- p+theme(axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=1))
      if(type=='State') p <- p+facet_wrap(~family,ncol=1,scales='free_y')+
        scale_x_discrete(labels=function(x) vapply(x,function(z) paste(strwrap(z,width=18),collapse='\n'),character(1)))
      print(p);page<-page+1L
      index[[paste(role,type,page)]] <- data.frame(role=role,file=paste0('figures/',role,'/',filename),page=page,variable=v,
        feature_type=type,features=paste(features,collapse=' | '),displayed_n=n_distinct(d$subject_id))
    }
  }
  dev.off();stopifnot(file.info(path)$size>1000)
}
index <- bind_rows(index);write_table(index,'figure_index')
summary <- objects$statistics %>% group_by(variable,family) %>% summarise(tests=n(),FDR_lt_05=sum(FDR<.05),.groups='drop')
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE)
summary_path <- file.path(SCREF_SUMMARY_DIR,'tcga_surgery_clinical_associations.csv')
write.csv(summary,summary_path,row.names=FALSE);stopifnot(file.info(summary_path)$size<100000)
files <- unique(index$file)
writeLines(c('# TCGA surgery-Yes clinical associations','',
  'The cohort contains 72 primary EAC patients with surgery status Yes in the corrected September 18 workbook.',
  'Priority figures cover sex, age, grade and pathological stage. Supporting figures cover other eligible clinical factors.',
  'Each clinical variable occupies one MP page containing all 17 MPs and one state page containing both state definitions.',
  '',paste0('- [',basename(files),'](../',files,')'),'',
  '[Clinical variable eligibility and source fields](../tables/variable_inventory.csv)',
  '[Complete tests](../tables/association_statistics.csv)',
  '[Figure/page index](../tables/figure_index.csv)','',
  'Scores retain the selected-survival 86-patient parent GSVA definition. Associations are unadjusted and exploratory.',
  'Unknown values are excluded. Levels require at least 10 patients. Bulk signals include nonmalignant cells.',
  'Corrected staging is authoritative; sparse clinical stage and site categories are not imputed.'),file.path(out,'reports/current_results.md'))
report <- c('STATUS: PASS','cohort=TCGA EAC surgery Yes; n=72',paste0('eligible_variables=',sum(objects$inventory$eligible)),
  paste0('tests=',nrow(objects$statistics)),paste0('pdfs=',length(files),'; pages=',nrow(index)),
  'metadata=corrected AJCC8 workbook plus cBioPortal fields absent from workbook; exact reconstruction sample mapping',
  'scores=17 MPs + 5 MP-union states + 5 top-20-DGE states; selected-survival parent 86; Gaussian raw-TPM GSVA',
  'statistics=unpaired Wilcoxon/Kruskal-Wallis; BH within clinical variable and signature family; minimum level n=10')
writeLines(report,file.path(out,'reports/run_summary.txt'));message(paste(report,collapse='\n'))
