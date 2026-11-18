####################
# Status: active
# Script: analysis/GEO/geo_microarray_quality.R
# Description: GSE19417 malignant-fraction, posterior-uncertainty and eight-state composition QC;
#   same three-panel design as the shared bulk-malignant deconvolution quality PDF.
# Methodology: analysis/methodology/GEO/geo_microarray_methodology.md.
# Inputs: ref_outs/GEO/intermediate/GSE_extracted.rds;
#   ref_outs/GEO/tables/clinical_metadata.csv; analysis/shared/scRef_config.R.
# Outputs under ref_outs/GEO/: figures/Auto_diagnostics/01_deconvolution_quality.pdf;
#   tables/{state_composition,composition_malignant_fractions,state_posterior_cv_summary,
#   state_fraction_resolution,quality_figure_index,quality_input_manifest}.csv;
#   reports/state_composition_results.md; logs/quality_packages.csv;
#   updates/new_updates/summaries/GEO_composition.csv (<100 KB).
# Cache/replot: deterministic rendering from the persistent extracted posterior; never refits BayesPrism.
#   SCREF_REPLOT_ONLY=TRUE and SCREF_FORCE_REBUILD=TRUE both refresh these derived outputs.
# Run: Rscript analysis/GEO/geo_microarray_quality.R
# PBS: job_scripts/GEO/geo_microarray_quality.sh; dmtcp; 1 CPU, 8 GB, 30 minutes.
####################
library(ggplot2)
library(patchwork)
source('analysis/shared/scRef_config.R')
if(!nzchar(Sys.getenv('PBS_JOBID')))stop('Run within PBS')
out<-file.path(SCREF_REF_OUTS_DIR,'GEO')
for(tier in c('figures/Auto_diagnostics','tables','reports','logs'))dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
inputs<-c(extracted=file.path(out,'intermediate/GSE_extracted.rds'),clinical=file.path(out,'tables/clinical_metadata.csv'),
 config='analysis/shared/scRef_config.R')
stopifnot(all(file.exists(inputs)))
x<-readRDS(inputs[['extracted']]);meta<-read.csv(inputs[['clinical']],stringsAsFactors=FALSE)
states<-c(SCREF_PRIMARY_STATE_ORDER,SCREF_TECHNICAL_STATE_ORDER,'Malignant state unassigned')
stopifnot(nrow(x$within)==52L,all(states%in%colnames(x$within)),
 setequal(rownames(x$within),meta$patient_id),all(is.finite(x$within)),all(x$within>=0),
 max(abs(rowSums(x$within)-1))<1e-6)
fractions<-x$fractions[match(rownames(x$within),x$fractions$sample_id),]
stopifnot(identical(fractions$sample_id,rownames(x$within)),all(is.finite(x$state_cv)),all(x$state_cv>=0),
 all(fractions$malignant_initial>=0&fractions$malignant_initial<=1),
 all(fractions$malignant_final>=0&fractions$malignant_final<=1))
fractions$treatment_group<-'All EAC specimens'
composition<-do.call(rbind,lapply(states,function(s)data.frame(sample_id=rownames(x$within),
 sample_rank=rank(fractions$malignant_final,ties.method='first'),state=s,fraction=x$within[,s],
 treatment_group=fractions$treatment_group)))
composition$state<-factor(composition$state,levels=states)
cv_summary<-do.call(rbind,lapply(states,function(s)data.frame(state=s,median_cv=median(x$state_cv[,s]),
 percentile95_cv=unname(quantile(x$state_cv[,s],.95)),maximum_cv=max(x$state_cv[,s]))))
resolution<-do.call(rbind,lapply(states,function(s)data.frame(state=s,n=nrow(x$within),
 median_within_fraction=median(x$within[,s]),samples_below_point_one_percent=sum(x$within[,s]<.001),
 median_initial_state_cv=median(x$state_cv[,s]),samples_cv_above_point_five=sum(x$state_cv[,s]>.5))))
table_out<-function(d,name)write.csv(d,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
table_out(composition,'state_composition');table_out(fractions,'composition_malignant_fractions')
table_out(cv_summary,'state_posterior_cv_summary');table_out(resolution,'state_fraction_resolution')
table_out(data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs))),'quality_input_manifest')
colours<-c(SCREF_STATE_COLOURS,'Malignant state unassigned'='#626262')
state_label<-function(v)vapply(v,function(s)paste(strwrap(gsub('-to-','-to- ',s,fixed=TRUE),width=22),collapse='\n'),character(1))
theme_set(theme_classic(base_size=14)+theme(legend.position='bottom',strip.background=element_blank(),
 plot.title=element_text(face='bold',size=16),plot.title.position='plot',axis.text=element_text(colour='black',size=12),
 legend.text=element_text(size=13),plot.margin=margin(10,12,10,10)))
p1<-ggplot(fractions,aes(malignant_initial,malignant_final))+
 geom_abline(slope=1,intercept=0,colour='grey65',linetype=2)+geom_point(size=1.8,alpha=.65,colour='#0072B2')+
 coord_fixed(xlim=c(0,1),ylim=c(0,1))+scale_x_continuous(breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+
 labs(x='Initial malignant RNA fraction',y='Updated malignant RNA fraction',title='a  Malignant compartment',subtitle='GSE19417; n=52')
cv_summary$state<-factor(cv_summary$state,levels=rev(states))
p3<-ggplot(cv_summary,aes(median_cv,state))+
 geom_errorbar(aes(xmin=median_cv,xmax=percentile95_cv),orientation='y',width=0,linewidth=.6,colour='#0072B2')+
 geom_point(size=2.3,colour='#0072B2')+scale_y_discrete(labels=state_label)+
 labs(x='Initial state-fraction posterior CV',y=NULL,title='b  Conditional state uncertainty',
 subtitle='Conditional on artificial intensity scale')+
 theme(plot.subtitle=element_text(size=11),axis.text.y=element_text(size=12))
p2<-ggplot(composition,aes(sample_rank,fraction,fill=state))+geom_col(width=1)+
 scale_fill_manual(values=colours,drop=FALSE)+scale_y_continuous(labels=scales::label_percent(),expand=c(0,0))+
 scale_x_continuous(expand=c(0,0))+guides(fill=guide_legend(ncol=2,byrow=TRUE))+
 labs(x='Samples ordered by updated malignant fraction',y='Within-malignant RNA fraction',fill=NULL,title='c  State composition')
p<-(p1|p3)/p2+plot_layout(heights=c(1,1.2))
file<-'figures/Auto_diagnostics/01_deconvolution_quality.pdf'
ggsave(file.path(out,file),p,width=14,height=12,device=cairo_pdf)
table_out(data.frame(path=file,page=1L,purpose='Malignant compartment, conditional uncertainty and 52-sample eight-state composition'),
 'quality_figure_index')
write.csv(data.frame(package=c('R','ggplot2','patchwork'),version=c(as.character(getRversion()),
 as.character(packageVersion('ggplot2')),as.character(packageVersion('patchwork')))),file.path(out,'logs/quality_packages.csv'),row.names=FALSE)
summary<-data.frame(samples=52L,genes=ncol(x$tumour_expression),median_updated_malignant_fraction=median(fractions$malignant_final),
 median_unassigned_within_malignant=median(x$within[,'Malignant state unassigned']),max_closure_error=max(abs(rowSums(x$within)-1)))
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE);write.csv(summary,file.path(SCREF_SUMMARY_DIR,'GEO_composition.csv'),row.names=FALSE)
writeLines(c('# GSE19417 state composition',
 'The single-page QC figure follows the shared deconvolution quality layout: initial versus updated malignant RNA fractions, conditional posterior uncertainty, and stacked state composition.',
 'All 52 EAC resection specimens are ordered by updated malignant fraction; tie order follows the extracted sample order. No new fitting or state assignment is performed.',
 'Composition uses the initial posterior, normalized across all eight malignant labels. The five biological states, Unresolved, Hybrid and Malignant state unassigned remain in the denominator.',
 'Fractions describe inferred RNA contributions, not validated cell counts. Posterior CV measures conditional sampling uncertainty rather than reference uncertainty.',
 'The resolution table identifies fractions below 0.1% and posterior CV above 0.5 descriptively; these thresholds do not filter any patient or state.',
 '[State-composition quality figure](../figures/Auto_diagnostics/01_deconvolution_quality.pdf)',
 '[Composition source data](../tables/state_composition.csv) | [Posterior CV](../tables/state_posterior_cv_summary.csv) | [Resolution](../tables/state_fraction_resolution.csv)'),
 file.path(out,'reports/state_composition_results.md'))
message('STATE_COMPOSITION_VALIDATED: 52 patients, eight malignant labels')
