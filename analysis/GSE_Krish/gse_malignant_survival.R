####################
# Status: active
# Script: analysis/GSE_Krish/gse_malignant_survival.R
# Description: OS associations for current 17 MPs and five-state signatures/fractions in GSE207526;
#   whole cohort and user-defined surgery-only/neoadjuvant strata; paired bulk/malignant comparison.
# Methodology: analysis/methodology/GSE_Krish/gse_krish_methodology.md.
# Inputs: ref_outs/GSE_Krish/intermediate/{GSE_counts,GSE_extracted}.rds;
#   tables/clinical_metadata_clean.csv; canonical centred MP genes/grouping and ranked state markers;
#   analysis/shared/scRef_config.R.
# Outputs under ref_outs/GSE_Krish/: intermediate/survival_cache.rds;
#   tables/{survival_input_manifest,feature_definitions,scoring_gene_coverage,scoring_gene_universe,
#   gene_sets,scored_features,model_data,cox_results,cutpoint_grid,bulk_malignant_comparison,
#   survival_figure_index}.csv; logs/survival_packages.csv;
#   figures/Auto_final/{01_all_cohort_survival,02_treatment_group_survival,03_bulk_vs_malignant}.pdf;
#   reports/survival_results.md; updates/new_updates/summaries/GSE_Krish_survival.csv (<100 KB).
# Cache: hashes/settings/statistical function bodies must match. SCREF_FORCE_REBUILD=TRUE rebuilds;
#   SCREF_REPLOT_ONLY=TRUE redraws from matching persistent cache. Gaussian GSVA uses all 110
#   biopsies in each compartment once; survival missingness does not change scoring cohorts.
# Run: Rscript analysis/GSE_Krish/gse_malignant_survival.R
# PBS: job_scripts/GSE_Krish/gse_malignant_survival.sh; dmtcp; 4 CPUs, 24 GB, 2 hours.
####################
library(GSVA)
library(survival)
library(ggplot2)
source('analysis/shared/scRef_config.R')
if(!nzchar(Sys.getenv('PBS_JOBID')))stop('Run within PBS')
force<-toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))=='TRUE'
replot<-toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))=='TRUE'
stopifnot(!(force && replot))
threads<-as.integer(Sys.getenv('SCREF_SURVIVAL_THREADS','4'))
out<-file.path(SCREF_REF_OUTS_DIR,'GSE_Krish')
for(tier in c('intermediate','tables','logs','reports','figures/Auto_final'))dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
inputs<-c(counts=file.path(out,'intermediate/GSE_counts.rds'),extracted=file.path(out,'intermediate/GSE_extracted.rds'),
 clinical=file.path(out,'tables/clinical_metadata_clean.csv'),genes=SCREF_MP_GENES_RDS,grouping=SCREF_MP_GROUPING_CSV,
 markers='ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv',config='analysis/shared/scRef_config.R')
stopifnot(all(file.exists(inputs)))
manifest<-data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
settings<-list(version=1L,min_n=20L,min_events=10L,binary_min_events=5L,min_group=5L,
 optimal_grid=seq(.2,.8,.05),multiplicity=32L,gsva='Gaussian log2(CPM+1); common nonconstant genes',seed=22062026L)
packages<-c('GSVA','survival','ggplot2')
versions<-data.frame(package=c('R',packages),version=c(as.character(getRversion()),vapply(packages,function(p)as.character(packageVersion(p)),character(1))))
table_out<-function(x,name)write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
cohort_labels<-c(all='All cohort',surgery='Surgery only',neo='Neoadjuvant')
family_labels<-c(MP='Metaprogrammes','State MP union'='MP-union\nstate scores','State DGE'='DGE\nstate scores',
 'State fraction'='Within-malignant\nRNA fractions')
method_labels<-c(continuous='Continuous',median='Median',q1q4='Q1 vs Q4',optimal_cut='Optimal cut (exploratory)')

make_scores<-function(){
 meta<-read.csv(inputs[['clinical']],stringsAsFactors=FALSE,na.strings=c('','NA'))
 stopifnot(nrow(meta)==110L,!anyDuplicated(meta$patient_id))
 meta$OS_time<-meta$os_days;meta$OS_event<-meta$os_event
 meta$OS_eligible<-is.finite(meta$OS_time)&meta$OS_time>0&meta$OS_event %in% c(0,1)
 grouping<-read.csv(inputs[['grouping']],stringsAsFactors=FALSE)
 genes<-readRDS(inputs[['genes']])[grouping$mp]
 stopifnot(length(genes)==17L,!any(c('MP2x','MP11c','MP18a')%in%names(genes)))
 states<-SCREF_PRIMARY_STATE_ORDER
 unions<-lapply(SCREF_STATE_GROUPS,function(mps)unique(unlist(genes[mps],use.names=FALSE)))
 names(unions)<-paste0(states,' (MP union)')
 markers<-read.csv(inputs[['markers']],stringsAsFactors=FALSE)
 dge<-lapply(states,function(s)head(markers$gene[markers$state==s],20L));names(dge)<-paste0(states,' (DGE)')
 stopifnot(all(lengths(dge)==20L),all(vapply(dge,function(x)!anyDuplicated(x),logical(1))))
 sets<-c(genes,unions,dge)
 defs<-data.frame(feature=c(names(sets),paste0(states,' (fraction)')),
  family=c(rep('MP',17),rep('State MP union',5),rep('State DGE',5),rep('State fraction',5)),
  state=c(grouping$state,rep(states,3)),label=c(paste(grouping$mp,grouping$description,sep=' | '),rep(states,3)))
 counts<-readRDS(inputs[['counts']]);e<-readRDS(inputs[['extracted']])
 stopifnot(setequal(colnames(counts),meta$patient_id),setequal(rownames(e$tumour_expression),meta$patient_id))
 common<-intersect(rownames(counts),colnames(e$tumour_expression));stopifnot(length(common)>1000L)
 matrices<-list(bulk=counts[common,meta$patient_id,drop=FALSE],malignant=t(e$tumour_expression[meta$patient_id,common,drop=FALSE]))
 matrices<-lapply(matrices,function(x){stopifnot(all(is.finite(x)),all(x>=0),all(colSums(x)>0));log2(sweep(x,2,colSums(x),'/')*1e6+1)})
 variable<-Reduce(`&`,lapply(matrices,function(x)apply(x,1,function(y)diff(range(y))>0)))
 common<-common[variable];used<-lapply(sets,function(g)intersect(unique(g),common))
 stopifnot(length(common)>1000L,all(lengths(used)>=5L))
 scores<-lapply(matrices,function(x){
  a<-GSVA::gsva(x[common,,drop=FALSE],used,method='gsva',kcdf='Gaussian',min.sz=5L,
   max.sz=Inf,mx.diff=TRUE,abs.ranking=FALSE,verbose=FALSE,parallel.sz=threads)
  stopifnot(identical(rownames(a),names(used)),all(is.finite(a)));t(a)})
 fractions<-e$within[meta$patient_id,states,drop=FALSE]
 stopifnot(all(is.finite(fractions)),all(fractions>=0),all(fractions<=1))
 colnames(fractions)<-paste0(states,' (fraction)');scores$malignant<-cbind(scores$malignant,fractions)
 long<-do.call(rbind,lapply(names(scores),function(src)do.call(rbind,lapply(colnames(scores[[src]]),function(f)
  data.frame(patient_id=rownames(scores[[src]]),source=src,feature=f,family=defs$family[match(f,defs$feature)],value=scores[[src]][,f])))))
 list(meta=meta,definitions=defs,sets=sets,scores=scores,long=long,common=common,
  coverage=data.frame(feature=names(used),defined=lengths(sets),retained=lengths(used)))
}

fit_model<-function(d,method,prob=NA_real_){
 d<-d[d$OS_eligible,,drop=FALSE]
 r<-data.frame(cohort=d$cohort[1],source=d$source[1],feature=d$feature[1],family=d$family[1],method=method,
  n=nrow(d),events=sum(d$OS_event),low_n=NA_integer_,high_n=NA_integer_,cutpoint=NA_real_,
  q1=NA_real_,q3=NA_real_,probability=prob,scale=NA_real_,status='estimated',logHR=NA_real_,SE=NA_real_,
  HR=NA_real_,lower_CI=NA_real_,upper_CI=NA_real_,P_raw=NA_real_,P_value=NA_real_,PH_P=NA_real_)
 fail<-function(reason){r$status<-paste('not estimated:',reason);r}
 if(method=='continuous'){
  r$scale<-if(d$family[1]=='State fraction').1 else sd(d$value)
  if(!is.finite(r$scale)||r$scale<=0)return(fail('constant predictor'))
  d$predictor<-d$value/r$scale
 }else{
  if(method=='median'){r$cutpoint<-median(d$value);high<-d$value>r$cutpoint}
  if(method=='q1q4'){
   q<-quantile(d$value,c(.25,.75),names=FALSE);r$q1<-q[1];r$q3<-q[2]
   if(q[1]>=q[2])return(fail('tied quartiles'))
   d<-d[d$value<=q[1]|d$value>=q[2],,drop=FALSE];high<-d$value>=q[2]
  }
  if(method=='optimal_cut'){r$cutpoint<-unname(quantile(d$value,prob));high<-d$value>=r$cutpoint}
  d$predictor<-as.integer(high);r$n<-nrow(d);r$events<-sum(d$OS_event);r$low_n<-sum(!high);r$high_n<-sum(high)
  if(min(r$low_n,r$high_n)<settings$min_group||length(unique(high))<2L||any(tapply(d$OS_event,high,sum)==0))
   return(fail('insufficient group/event support'))
 }
 if(r$n<settings$min_n||r$events<if(method=='continuous')settings$min_events else settings$binary_min_events)
  return(fail('insufficient patients/events'))
 fit<-withCallingHandlers(coxph(Surv(OS_time,OS_event)~predictor,data=d,ties='efron',x=TRUE,control=coxph.control(iter.max=50L)),
  warning=function(w)stop(r$cohort,'/',r$source,'/',r$feature,'/',method,': ',conditionMessage(w)))
 s<-summary(fit);r$logHR<-unname(coef(fit));r$SE<-s$coefficients[1,'se(coef)']
 r$HR<-s$conf.int[1,'exp(coef)'];r$lower_CI<-s$conf.int[1,'lower .95'];r$upper_CI<-s$conf.int[1,'upper .95']
 r$P_raw<-r$P_value<-s$coefficients[1,'Pr(>|z|)']
 stopifnot(all(is.finite(c(r$logHR,r$SE,r$P_value))))
 if(method=='continuous')r$PH_P<-cox.zph(fit,transform='km')$table['predictor','p']
 r
}

run_models<-function(scored){
 cohorts<-list(all=scored$meta$patient_id,surgery=scored$meta$patient_id[scored$meta$treatment_group=='Surgery only'],
  neo=scored$meta$patient_id[scored$meta$treatment_group=='Neoadjuvant'])
 stopifnot(identical(lengths(cohorts),c(all=110L,surgery=42L,neo=68L)))
 data<-do.call(rbind,lapply(names(cohorts),function(k){
  x<-scored$long[scored$long$patient_id%in%cohorts[[k]],];j<-match(x$patient_id,scored$meta$patient_id)
  data.frame(cohort=k,x,OS_time=scored$meta$OS_time[j],OS_event=scored$meta$OS_event[j],OS_eligible=scored$meta$OS_eligible[j])}))
 results<-grids<-list();i<-0L
 for(d in split(data,interaction(data$cohort,data$source,data$feature,drop=TRUE))){
  i<-i+1L
  grid<-do.call(rbind,lapply(settings$optimal_grid,function(p)fit_model(d,'optimal_cut',p)))
  eligible<-which(grid$status=='estimated')
  if(length(eligible)){
   best<-grid[eligible[which.min(grid$P_raw[eligible])],,drop=FALSE]
   best$P_value<-min(1,best$P_raw*nrow(grid))
  }else{best<-grid[1,,drop=FALSE];best$status<-'not estimated: no supported cutpoint'}
  results[[i]]<-rbind(fit_model(d,'continuous'),fit_model(d,'median'),fit_model(d,'q1q4'),best);grids[[i]]<-grid
 }
 results<-do.call(rbind,results);results$P_BH<-NA_real_
 for(j in split(seq_len(nrow(results)),interaction(results$cohort,results$source,results$method,drop=TRUE)))
  results$P_BH[j]<-p.adjust(results$P_value[j],'BH',n=settings$multiplicity)
 list(data=data,results=results,grid=do.call(rbind,grids),cohorts=cohorts)
}

code<-lapply(list(make_scores,fit_model,run_models),function(f)paste(deparse(body(f)),collapse='\n'))
contract<-list(manifest=manifest,settings=settings,versions=versions,code=code)
cache_path<-file.path(out,'intermediate/survival_cache.rds')
if(file.exists(cache_path)&&!force){cache<-readRDS(cache_path);if(!identical(cache$contract,contract))stop('Stale survival cache: explicitly rebuild')
}else{if(replot)stop('Replot requires matching cache');scored<-make_scores();models<-run_models(scored);cache<-list(contract=contract,scored=scored,models=models);saveRDS(cache,cache_path)}
scored<-cache$scored;results<-cache$models$results
stopifnot(nrow(results)==3L*(32L+27L)*4L)
table_out(manifest,'survival_input_manifest');table_out(scored$definitions,'feature_definitions')
table_out(scored$coverage,'scoring_gene_coverage');table_out(data.frame(gene=scored$common),'scoring_gene_universe')
table_out(do.call(rbind,lapply(names(scored$sets),function(f)data.frame(feature=f,gene=scored$sets[[f]],retained=scored$sets[[f]]%in%scored$common))),'gene_sets')
table_out(scored$long,'scored_features');table_out(cache$models$data,'model_data');table_out(results,'cox_results');table_out(cache$models$grid,'cutpoint_grid')
a<-results[results$source=='bulk',];b<-results[results$source=='malignant'&results$family!='State fraction',]
paired<-merge(a,b,by=c('cohort','feature','family','method'),suffixes=c('_bulk','_malignant'))
stopifnot(nrow(paired)==3L*27L*4L,all(paired$n_bulk==paired$n_malignant|paired$method!='continuous'))
paired$delta_logHR<-paired$logHR_malignant-paired$logHR_bulk;paired$direction_changed<-paired$logHR_bulk*paired$logHR_malignant<0
# Binary cutpoints are compartment-specific, exactly as in the comparator workflow; continuous comparisons use identical patients.
table_out(paired,'bulk_malignant_comparison');write.csv(versions,file.path(out,'logs/survival_packages.csv'),row.names=FALSE)

format_feature<-function(f){labels<-scored$definitions$label[match(f,scored$definitions$feature)];vapply(labels,function(z)paste(strwrap(z,width=55),collapse='\n'),character(1))}
prepare_plot<-function(d){
 d$feature_plot<-factor(d$feature,levels=rev(scored$definitions$feature));d$family<-factor(d$family,levels=names(family_labels))
 d$significance<-factor(ifelse(is.na(d$P_BH),'NS',ifelse(d$P_BH<.05,'FDR',ifelse(d$P_value<.05,'Nominal','NS'))),levels=c('NS','Nominal','FDR'))
 d$effect<-pmax(-3,pmin(3,d$logHR/log(2)));d$left<-pmax(-3,pmin(3,(d$logHR-qnorm(.975)*d$SE)/log(2)))
 d$right<-pmax(-3,pmin(3,(d$logHR+qnorm(.975)*d$SE)/log(2)))
 d$overflow_left<-(d$logHR-qnorm(.975)*d$SE)/log(2)< -3;d$overflow_right<-(d$logHR+qnorm(.975)*d$SE)/log(2)>3
 d
}
forest_base<-function(d,paired=FALSE){
 e<-d[d$status=='estimated',]
 ####################
 shape_scale<-if(paired)scale_shape_manual(values=c('Whole bulk'=16,'Malignant expression'=17),drop=FALSE,name=NULL) else
  scale_shape_manual(values=c(NS=1,Nominal=16,FDR=18),drop=FALSE,name='Evidence')
 ####################
 ggplot(d,aes(y=feature_plot,colour=display_group,group=display_group))+
  geom_blank(aes(x=0))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.4)+
  geom_errorbar(data=e,aes(x=effect,xmin=left,xmax=right),orientation='y',width=.12,position=position_dodge(.6),linewidth=.5)+
  geom_point(data=e,aes(x=effect,shape=if(paired)display_group else significance),position=position_dodge(.6),size=2.5,show.legend=TRUE)+
  geom_text(data=e,aes(x=-3,label=ifelse(overflow_left,'\u2039','')),position=position_dodge(.6),size=5,show.legend=FALSE)+
  geom_text(data=e,aes(x=3,label=ifelse(overflow_right,'\u203a','')),position=position_dodge(.6),size=5,show.legend=FALSE)+
  scale_x_continuous(limits=c(-3.15,3.15),breaks=if(paired)c(-2,0,2) else -3:3,
   labels=if(paired)c('0.25','1','4') else c('.125','.25','.5','1','2','4','8'))+
  scale_y_discrete(labels=format_feature)+shape_scale+
  theme_classic(base_size=14)+theme(legend.position='bottom',legend.title=element_blank(),axis.text.y=element_text(size=13.5,colour='black'),
   axis.text.x=element_text(size=13,colour='black'),strip.background=element_blank(),strip.text.y=element_text(angle=0,size=13),
   strip.text.x=element_text(size=13),plot.title=element_text(size=18,face='bold'),plot.title.position='plot',
   plot.margin=margin(12,14,12,12))
}
figure_index<-list()
for(kind in c('all','groups')){
 file<-paste0('figures/Auto_final/',if(kind=='all')'01_all_cohort_survival.pdf' else '02_treatment_group_survival.pdf')
 cairo_pdf(file.path(out,file),width=13,height=14,onefile=TRUE)
 for(method in names(method_labels)){
  keys<-if(kind=='all')'all' else c('surgery','neo')
  d<-prepare_plot(results[results$source=='malignant'&results$cohort%in%keys&results$method==method,])
  counts<-vapply(keys,function(k){z<-d[d$cohort==k,];paste0(cohort_labels[[k]],': n=',paste(unique(range(z$n)),collapse='–'),', events=',paste(unique(range(z$events)),collapse='–'))},character(1))
  d$display_group<-factor(d$cohort,levels=keys,labels=unname(cohort_labels[keys]))
  p<-forest_base(d)+facet_grid(family~.,scales='free_y',space='free_y',labeller=labeller(family=family_labels))+
   scale_colour_manual(values=if(kind=='all')c('All cohort'='#333333') else c('Surgery only'='#D55E00',Neoadjuvant='#0072B2'))+
   labs(x=if(method=='continuous')'Hazard ratio per SD / 10 percentage points' else 'Hazard ratio: high vs low',y=NULL,
    title=paste(if(kind=='all')'GSE207526 all cohort' else 'GSE207526 treatment groups','|',method_labels[[method]]),
    subtitle=paste(counts,collapse='; '))
  print(p);figure_index[[length(figure_index)+1L]]<-data.frame(path=file,page=match(method,names(method_labels)),method=method,cohorts=paste(keys,collapse=' | '),features=32L)
 }
 dev.off()
}
file<-'figures/Auto_final/03_bulk_vs_malignant.pdf'
cairo_pdf(file.path(out,file),width=14,height=13,onefile=TRUE)
d<-prepare_plot(results[results$method=='continuous'&results$family!='State fraction',])
d$cohort_plot<-factor(d$cohort,levels=names(cohort_labels),labels=unname(cohort_labels))
d$display_group<-factor(d$source,levels=c('bulk','malignant'),labels=c('Whole bulk','Malignant expression'))
p<-forest_base(d,paired=TRUE)+facet_grid(family~cohort_plot,scales='free_y',space='free_y',labeller=labeller(family=family_labels))+
 scale_colour_manual(values=c('Whole bulk'='#999999','Malignant expression'='#0072B2'),name=NULL)+
 labs(x='Unadjusted hazard ratio per score SD (95% CI)',y=NULL,title='Paired bulk and malignant-expression associations')
print(p);dev.off()
figure_index[[length(figure_index)+1L]]<-data.frame(path=file,page=1L,method='continuous',cohorts='all | surgery | neo',features=27L)
table_out(do.call(rbind,figure_index),'survival_figure_index')
summary<-aggregate(cbind(nominal=as.integer(results$P_value<.05),FDR=as.integer(results$P_BH<.05)),
 list(cohort=results$cohort,source=results$source,method=results$method),sum,na.rm=TRUE)
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE);write.csv(summary,file.path(SCREF_SUMMARY_DIR,'GSE_Krish_survival.csv'),row.names=FALSE)
writeLines(c('# GSE207526 malignant-expression survival',
 'All 110 tumour biopsies were scored; 107 have usable OS status: 68 neoadjuvant (22 deaths) and 39 surgery-only (24 deaths).',
 'Treatment groups are user-defined from Mandard availability, not confirmed treatment records. Biopsies were collected before therapy.',
 'All comparisons are unadjusted exploratory associations. Age is only available in bands; baseline AJCC edition and OS origin are not specified.',
 'Continuous Cox HRs are per subgroup-specific score SD; state fractions use 10 percentage points. Efron ties. PH tests are tabulated.',
 'Median, Q1/Q4 and optimal-cut models follow the comparator support thresholds. Unsupported tests are explicitly not estimated.',
 'Optimal cutpoints search 13 fixed quantiles; P=min(1,13*minimum Wald P), followed by BH across the fixed 32-feature family within cohort/source/method.',
 'Selected cutpoint HR/CI are exploratory and do not account for selection. Methods and overlapping subgroups are not independent validations.',
 'Whole-bulk and malignant scores use identical retained genes, log2(CPM+1), and Gaussian GSVA on all 110 specimens. No rescoring by treatment group.',
 'Bulk/malignant continuous models use identical patients; differences are paired descriptive comparisons, not independent tests.',
 'Deconvolution uncertainty is not propagated. State fractions describe RNA contributions, not measured cell fractions.',
 '[All cohort](../figures/Auto_final/01_all_cohort_survival.pdf)',
 '[Treatment groups](../figures/Auto_final/02_treatment_group_survival.pdf)',
 '[Bulk versus malignant](../figures/Auto_final/03_bulk_vs_malignant.pdf)',
 '[Cox results](../tables/cox_results.csv) | [Cutpoint candidates](../tables/cutpoint_grid.csv)',
 '[Model source data](../tables/model_data.csv) | [Gene coverage](../tables/scoring_gene_coverage.csv)'),file.path(out,'reports/survival_results.md'))
message('SURVIVAL_VALIDATED: ',nrow(results),' model rows; 9 PDF pages')
