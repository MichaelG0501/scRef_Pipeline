####################
# Status: active
# Script: analysis/GSE_Krish/gse_malignant_clinical_associations.R
# Description: GSE207526 malignant MP/state associations with supplied sex, age <50 versus >=50, and broad baseline stage;
#   all, surgery-only and neoadjuvant selections. Mandard analyses only within neoadjuvant.
# Methodology: analysis/methodology/GSE_Krish/gse_krish_methodology.md.
# Inputs: ref_outs/GSE_Krish/intermediate/survival_cache.rds; tables/clinical_metadata_clean.csv;
#   analysis/shared/scRef_config.R. Scores are reused exactly; no score recalculation.
# Outputs under ref_outs/GSE_Krish/: intermediate/clinical_cache.rds;
#   tables/{clinical_variable_inventory,clinical_statistics,clinical_plot_source,mandard_ordinal_statistics,
#   clinical_figure_index,clinical_subject_metadata}.csv; figures/Auto_final/04_all_MP_clinical.pdf through 09_neo_State_clinical.pdf;
#   reports/clinical_results.md; updates/new_updates/summaries/GSE_Krish_clinical.csv (<100 KB).
# Cache: metadata/score hashes and statistical code must match. SCREF_FORCE_REBUILD=TRUE rebuilds;
#   SCREF_REPLOT_ONLY=TRUE renders matching cache. Categories require >=3 patients; missing not imputed.
# Run: Rscript analysis/GSE_Krish/gse_malignant_clinical_associations.R
# PBS: job_scripts/GSE_Krish/gse_malignant_clinical_associations.sh; dmtcp; 1 CPU, 8 GB, 1 hour.
####################
library(dplyr)
library(ggplot2)
source('analysis/shared/scRef_config.R')
if(!nzchar(Sys.getenv('PBS_JOBID')))stop('Run within PBS')
out<-file.path(SCREF_REF_OUTS_DIR,'GSE_Krish')
force<-toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))=='TRUE'
replot<-toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))=='TRUE';stopifnot(!(force&&replot))
parent_path<-file.path(out,'intermediate/survival_cache.rds');metadata_path<-file.path(out,'tables/clinical_metadata_clean.csv')
stopifnot(file.exists(parent_path),file.exists(metadata_path))
parent<-readRDS(parent_path);defs<-parent$scored$definitions
cohort_labels<-c(all='All cohort',surgery='Surgery only',neo='Neoadjuvant')
variables<-c(sex='Sex',age_group='Age at diagnosis (<50 vs ≥50 years)',baseline_stage='Baseline AJCC stage',
 treatment_group='Mandard-defined treatment group',
 mandard_score='Mandard grade',mandard_response='Mandard TRG1–3 versus TRG4–5')
table_out<-function(x,name)write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')

build_clinical<-function(){
 meta<-read.csv(metadata_path,stringsAsFactors=FALSE,na.strings=c('','NA'))
 stopifnot(identical(meta$patient_id,parent$scored$meta$patient_id))
 ####################
 # The supplied five-year intervals classify <50 versus >=50 without inferring exact ages.
 bounds<-strsplit(meta$age_band,'-',fixed=TRUE)
 stopifnot(all(lengths(bounds)==2L))
 lower_age<-as.numeric(vapply(bounds,`[`,character(1),1L))
 upper_age<-as.numeric(vapply(bounds,`[`,character(1),2L))
 stopifnot(all(is.finite(lower_age)),all(is.finite(upper_age)),all(lower_age<=upper_age),
  all(upper_age<50 | lower_age>=50))
 meta$age_group<-ifelse(upper_age<50,'<50','≥50')
 ####################
 ####################
 # Preserve the source text in the clinical CSV, but exclude unstageable records from stage tests.
 meta$baseline_stage_detailed[!grepl('^stage (I|II|III|IV)[ABC]?$',meta$baseline_stage_detailed)]<-NA_character_
 ####################
 meta$baseline_stage<-sub('^stage (IV|III|II|I).*','\\1',meta$baseline_stage_detailed)
 meta$baseline_stage[!meta$baseline_stage%in%c('I','II','III','IV')]<-NA_character_
 meta$mandard_score<-as.character(meta$mandard_score)
 inv<-plots<-stats<-list();i<-0L
 for(cohort in names(cohort_labels)){
  m<-meta[meta$patient_id%in%parent$models$cohorts[[cohort]],]
  for(variable in names(variables)){
   i<-i+1L;counts<-table(m[[variable]],useNA='no');retained<-names(counts)[counts>=3L]
   applicable<-!(variable%in%c('mandard_score','mandard_response')&&cohort!='neo')
   applicable<-applicable&&!(variable=='treatment_group'&&cohort!='all')
   eligible<-applicable&&length(retained)>=2L
   inv[[i]]<-data.frame(cohort=cohort,variable=variable,label=variables[[variable]],selected=nrow(m),
    missing=sum(is.na(m[[variable]])),category_counts=paste(names(counts),counts,sep=':',collapse=' | '),
    retained_categories=paste(retained,collapse=' | '),n_retained=if(applicable)sum(counts[retained]) else 0L,
    eligible=eligible,reason=if(!applicable)'not applicable by prespecified cohort definition' else if(!eligible)'fewer than two categories with >=3 patients' else 'eligible')
   if(!eligible)next
   m<-m[order(m$patient_id),,drop=FALSE]
   selected<-m[m[[variable]]%in%retained,,drop=FALSE]
   scores<-parent$scored$long[parent$scored$long$source=='malignant'&parent$scored$long$patient_id%in%selected$patient_id,]
   scores$group<-selected[[variable]][match(scores$patient_id,selected$patient_id)]
   scores$cohort<-cohort;scores$variable<-variable;plots[[i]]<-scores
   stats[[i]]<-do.call(rbind,lapply(split(scores,scores$feature),function(z){
    p<-NA_real_;test<-if(length(unique(z$group))==2L)'Wilcoxon' else 'Kruskal-Wallis'
    status<-'estimated'
    if(diff(range(z$value))==0)status<-'not estimated: constant score' else {
     p<-if(test=='Wilcoxon')wilcox.test(value~group,data=z,exact=FALSE)$p.value else kruskal.test(value~group,data=z)$p.value
    }
    data.frame(cohort=cohort,variable=variable,feature=z$feature[1],family=z$family[1],n=nrow(z),test=test,status=status,p=p)
   }))
  }
 }
 stats<-do.call(rbind,stats);stats$q_family<-stats$q_32<-NA_real_
 for(j in split(seq_len(nrow(stats)),interaction(stats$cohort,stats$variable,stats$family,drop=TRUE)))
  stats$q_family[j]<-p.adjust(stats$p[j],'BH',n=if(stats$family[j[1]]=='MP')17L else 5L)
 for(j in split(seq_len(nrow(stats)),interaction(stats$cohort,stats$variable,drop=TRUE)))stats$q_32[j]<-p.adjust(stats$p[j],'BH',n=32L)
 neo<-meta[meta$treatment_group=='Neoadjuvant',]
 ordinal<-do.call(rbind,lapply(split(parent$scored$long[parent$scored$long$source=='malignant'&parent$scored$long$patient_id%in%neo$patient_id,],
  parent$scored$long$feature[parent$scored$long$source=='malignant'&parent$scored$long$patient_id%in%neo$patient_id]),function(z){
   grade<-as.numeric(neo$mandard_score[match(z$patient_id,neo$patient_id)])
   stopifnot(nrow(z)==68L,all(is.finite(grade)))
   t<-cor.test(z$value,grade,method='spearman',exact=FALSE)
   data.frame(cohort='neo',feature=z$feature[1],family=z$family[1],n=nrow(z),rho=unname(t$estimate),p=t$p.value)
  }))
 ordinal$q_family<-NA_real_
 for(j in split(seq_len(nrow(ordinal)),ordinal$family))ordinal$q_family[j]<-p.adjust(ordinal$p[j],'BH',n=if(ordinal$family[j[1]]=='MP')17L else 5L)
 list(metadata=meta,inventory=do.call(rbind,inv),plot_data=do.call(rbind,plots),statistics=stats,ordinal=ordinal)
}
contract<-list(parent=unname(tools::md5sum(parent_path)),metadata=unname(tools::md5sum(metadata_path)),
 variables=variables,min_category_n=3L,code=paste(deparse(body(build_clinical)),collapse='\n'),
 versions=vapply(c('R'=NA_character_,'dplyr','ggplot2'),function(p)if(is.na(p))as.character(getRversion()) else as.character(packageVersion(p)),character(1)))
path<-file.path(out,'intermediate/clinical_cache.rds')
if(file.exists(path)&&!force){saved<-readRDS(path);if(!identical(saved$contract,contract))stop('Stale clinical cache: explicitly rebuild');objects<-saved$objects
}else{if(replot)stop('Replot requires matching cache');objects<-build_clinical();saveRDS(list(contract=contract,objects=objects),path)}
table_out(objects$inventory,'clinical_variable_inventory');table_out(objects$statistics,'clinical_statistics')
table_out(objects$plot_data,'clinical_plot_source');table_out(objects$ordinal,'mandard_ordinal_statistics')
####################
table_out(objects$metadata,'clinical_subject_metadata')
####################
figure_index<-list()
for(cohort in names(cohort_labels))for(type in c('MP','State')){
 eligible<-objects$inventory$variable[objects$inventory$cohort==cohort&objects$inventory$eligible]
 stopifnot(length(eligible)>0L)
 ####################
 number<-3L+(match(cohort,names(cohort_labels))-1L)*2L+match(type,c('MP','State'))
 file<-sprintf('figures/Auto_final/%02d_%s_%s_clinical.pdf',number,cohort,type)
 ####################
 cairo_pdf(file.path(out,file),width=if(type=='MP')14 else 13,height=if(type=='MP')10 else 12,onefile=TRUE)
 for(page in seq_along(eligible)){
  variable<-eligible[page]
  features<-defs$feature[if(type=='MP')defs$family=='MP' else defs$family!='MP']
  d<-objects$plot_data[objects$plot_data$cohort==cohort&objects$plot_data$variable==variable&objects$plot_data$feature%in%features,]
  d$display_score<-ifelse(d$family=='State fraction',100*d$value,d$value)
  labels<-if(type=='MP')setNames(vapply(features,function(f)paste(c(f,strwrap(SCREF_MP_DESCRIPTIONS[[f]],width=24)),collapse='\n'),character(1)),features) else
   setNames(defs$state[match(features,defs$feature)],features)
  d$x<-factor(unname(labels[d$feature]),levels=unique(unname(labels[features])))
  d$family<-factor(d$family,levels=c('MP','State MP union','State DGE','State fraction'))
  groups<-sort(unique(d$group));d$group<-factor(d$group,levels=groups)
  counts<-d%>%distinct(patient_id,group)%>%count(group,.drop=FALSE)
  palette<-setNames(scales::hue_pal()(length(groups)),groups)
  legend_labels<-setNames(paste0(counts$group,' (n=',counts$n,')'),as.character(counts$group))
  ann<-d%>%group_by(feature,x,family)%>%summarise(y=max(display_score)+.10*diff(range(display_score))+.03,.groups='drop')%>%
   left_join(objects$statistics[objects$statistics$cohort==cohort&objects$statistics$variable==variable,],by=c('feature','family'))
  ann$star<-ifelse(is.na(ann$q_family),'',ifelse(ann$q_family<.001,'***',ifelse(ann$q_family<.01,'**',ifelse(ann$q_family<.05,'*',''))))
  ####################
  ann$family<-factor(ann$family,levels=levels(d$family))
  ####################
  p<-ggplot(d,aes(x,display_score,fill=group,colour=group))+
   geom_boxplot(position=position_dodge(.75),width=.6,outlier.shape=NA,alpha=.8,linewidth=.4,colour='black')+
   geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75,seed=1701),size=1,alpha=.65,stroke=0,show.legend=FALSE)+
   geom_text(data=ann[ann$star!='',],aes(x=x,y=y,label=star),inherit.aes=FALSE,size=4.5,fontface='bold')+
   scale_fill_manual(values=palette,labels=legend_labels,drop=FALSE)+scale_colour_manual(values=palette,guide='none',drop=FALSE)+
   scale_y_continuous(expand=expansion(mult=c(.04,.16)))+
   labs(title=paste(cohort_labels[[cohort]],'|',variables[[variable]]),x=NULL,y=if(type=='MP')'GSVA enrichment score' else 'Score / fraction (%)',fill=NULL)+
   theme_classic(base_size=14)+theme(plot.title=element_text(size=18,face='bold'),axis.text=element_text(size=14.5,colour='black'),
    legend.position='top',legend.text=element_text(size=14),strip.text=element_text(size=14),plot.margin=margin(12,16,12,12))
  if(type=='MP')p<-p+theme(axis.text.x=element_text(size=14.5,angle=90,hjust=1,vjust=.5))
  if(type=='State')p<-p+facet_wrap(~family,ncol=1,scales='free_y',labeller=labeller(family=c(
   'State MP union'='MP-union score','State DGE'='DGE score','State fraction'='Malignant RNA (%)')))+
   scale_x_discrete(labels=function(x)vapply(x,function(s)paste(strwrap(gsub('-to-','-to- ',s,fixed=TRUE),width=18),collapse='\n'),character(1)))
  print(p);figure_index[[length(figure_index)+1L]]<-data.frame(path=file,page=page,cohort=cohort,variable=variable,type=type,features=length(features))
 }
 dev.off()
}
table_out(do.call(rbind,figure_index),'clinical_figure_index')
summary<-objects$statistics%>%group_by(cohort,variable,family)%>%summarise(tests=n(),nominal=sum(p<.05,na.rm=TRUE),FDR=sum(q_family<.05,na.rm=TRUE),.groups='drop')
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE);write.csv(summary,file.path(SCREF_SUMMARY_DIR,'GSE_Krish_clinical.csv'),row.names=FALSE)
writeLines(c('# GSE207526 malignant-expression clinical associations',
 'All 110 patients are retained independently of OS eligibility. The Mandard-defined strata contain 68 neoadjuvant and 42 surgery-only patients.',
 'Sex, five-year age bands and detailed baseline AJCC stage are supplied in the original metadata. Broad stage collapses only the documented Roman-numeral stage; no AJCC edition is inferred.',
 'Age comparisons use <50 versus ≥50, derived exactly from the supplied five-year intervals: all cohort 9/101; neoadjuvant 8/60; surgery-only 1/41. Exact ages are not inferred. The surgery-only age comparison is ineligible under the minimum-three-per-category rule.',
 'Detailed AJCC substages and five-year age-band comparisons are excluded from tests and figures. Broad stage II/III/IV remains.',
 'Mandard grade and TRG1–3 versus TRG4–5 are analysed only within neoadjuvant. They are not repeated in the combined or surgery-only groups.',
 'Each categorical comparison retains categories with at least three patients, and requires at least two such categories. Excluded category counts and missingness are tabulated.',
 'One TRG5 patient is excluded from the five-grade boxplot/test by this threshold, but remains in the binary TRG comparison and the 68-patient ordinal Spearman analysis.',
 'Two-group tests use two-sided unpaired Wilcoxon with normal approximation and tie handling. Multi-group tests use Kruskal-Wallis.',
 'Figure stars use BH within cohort/variable/signature family (17 MPs or five state features); broader 32-feature q values are also tabulated.',
 'Ordinal Mandard Spearman tests include all 68 graded patients; BH is calculated within each signature family. These complement the categorical exploratory comparisons.',
 'One clinical variable per page; all 17 MPs remain together. State pages contain MP-union, DGE and RNA-fraction rows. No page splits are added for font size.',
 'Age is not converted to exact years or an arbitrary midpoint. Missing clinical values are never imputed. Treatment group comparisons are descriptive user-defined proxy comparisons.',
 '[Clinical statistics](../tables/clinical_statistics.csv) | [Variable inventory](../tables/clinical_variable_inventory.csv)',
 '[Ordinal Mandard statistics](../tables/mandard_ordinal_statistics.csv) | [Figure index](../tables/clinical_figure_index.csv)'),file.path(out,'reports/clinical_results.md'))
message('CLINICAL_VALIDATED: ',nrow(objects$statistics),' categorical and ',nrow(objects$ordinal),' ordinal tests')
