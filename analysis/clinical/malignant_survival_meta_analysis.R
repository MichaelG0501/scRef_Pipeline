####################
# Status: active
# Script: analysis/clinical/malignant_survival_meta_analysis.R
# Description: All-patient malignant OS associations and configured-dataset meta-analysis;
#   no surgery, treatment or specimen subgroup restriction. One patient per dataset.
# Methodology: analysis/methodology/clinical/malignant_survival_meta_analysis_methodology.md
# Exact live inputs:
#   ref_outs/clinical_malignant_survival/tables/{patient_metadata,scores_and_fractions,
#     feature_definitions,input_manifest}.csv
#   ref_outs/GSE_Krish/tables/{model_data,feature_definitions,survival_input_manifest}.csv
#   With SCREF_INCLUDE_GEO=TRUE: ref_outs/GEO/tables/{model_data,scoring_input_manifest}.csv.
# Metadata authority: OCCAMS Patient_Master_388 in the 2026-09-25 corrected workbook;
#   TCGA EAC_RNA88_Final in the 2026-09-18 workbook plus explicit cBioPortal OS fields;
#   GSE207526_GSE112926_meta.txt, curated by analysis/GSE_Krish/geo_oac_bulk_metadata_prep.py.
#   Exact source paths and hashes are validated through the two producer manifests.
#   Missing endpoint fields are never substituted. All recorded OS origins are preserved.
# Outputs: ref_outs/clinical_malignant_meta/intermediate/meta_analysis.rds;
#   tables/{input_manifest,dataset_summary,model_data,cohort_estimates,cutpoint_grid,
#     meta_estimates,leave_one_dataset_out,feature_definitions,figure_index}.csv;
#   figures/Auto_final/{01_mp_meta_analysis,02_state_meta_analysis}.pdf (three pages each; 20 x 12 inches);
#   figures/Auto_supporting/01_optimal_cut_exploratory.pdf (two pages);
#   reports/{current_results.md,validation.csv}; logs/package_versions.csv;
#   updates/new_updates/summaries/malignant_survival_meta.csv (<100 KB).
# Cache: input hashes, statistical function bodies and package versions must match.
#   SCREF_INCLUDE_GEO=TRUE explicitly adds GSE19417; all its inputs must exist.
#   Default TRUE includes GEO; FALSE explicitly selects RNA-seq-only analysis.
#   SCREF_FORCE_REBUILD=TRUE rebuilds; SCREF_REPLOT_ONLY=TRUE requires matching cache.
# Run: Rscript analysis/clinical/malignant_survival_meta_analysis.R
# PBS: job_scripts/clinical/malignant_survival_meta_analysis.sh; dmtcp; 1 CPU/8 GB/1 hour.
####################
library(survival)
library(ggplot2)
source('analysis/shared/scRef_config.R')
stopifnot(nzchar(Sys.getenv('PBS_JOBID')))
force <- toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))=='TRUE'
replot <- toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))=='TRUE'
####################
include_geo <- toupper(Sys.getenv('SCREF_INCLUDE_GEO','TRUE'))=='TRUE'
####################
stopifnot(!(force && replot))
out <- 'ref_outs/clinical_malignant_meta'
for(tier in c('intermediate','tables','reports','logs','figures/Auto_final','figures/Auto_supporting'))
  dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
inputs <- c(metadata='ref_outs/clinical_malignant_survival/tables/patient_metadata.csv',
 scores='ref_outs/clinical_malignant_survival/tables/scores_and_fractions.csv',
 definitions='ref_outs/clinical_malignant_survival/tables/feature_definitions.csv',
 manifest='ref_outs/clinical_malignant_survival/tables/input_manifest.csv',
 gse_data='ref_outs/GSE_Krish/tables/model_data.csv',
 gse_definitions='ref_outs/GSE_Krish/tables/feature_definitions.csv',
 gse_manifest='ref_outs/GSE_Krish/tables/survival_input_manifest.csv')
####################
if(include_geo)inputs <- c(inputs,geo_data='ref_outs/GEO/tables/model_data.csv',
 geo_manifest='ref_outs/GEO/tables/scoring_input_manifest.csv')
####################
stopifnot(all(file.exists(inputs)))
read_table <- function(p) read.csv(p,stringsAsFactors=FALSE,na.strings=c('','NA'))
manifest <- data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
for(key in c('manifest','gse_manifest',if(include_geo)'geo_manifest')) {
 m <- read_table(inputs[[key]])
 stopifnot(all(file.exists(m$path)),identical(unname(tools::md5sum(m$path)),m$md5))
}
versions <- data.frame(package=c('R','survival','ggplot2'),version=c(as.character(getRversion()),
 as.character(packageVersion('survival')),as.character(packageVersion('ggplot2'))))
settings <- list(parent_n=c(OCCAMS=281L,TCGA=88L,GSE207526=110L),
 include_geo=include_geo,min_n=20L,min_events=10L,min_group=5L,optimal_grid=seq(.2,.8,.05),multiplicity=32L,
 meta='REML; modified Hartung-Knapp, t with k-1 df; all configured datasets required',version=1L)
####################
if(include_geo)settings$parent_n <- c(settings$parent_n,GSE19417=52L)
datasets <- names(settings$parent_n)
####################
methods <- c(continuous='Continuous',median='Median',q1q4='Q1 vs Q4')
write_table <- function(x,name) write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')

prepare_data <- function() {
 meta <- read_table(inputs[['metadata']]);scores <- read_table(inputs[['scores']])
 scores <- scores[grepl('__malignant$',scores$parent),]
 scores$dataset <- sub('__malignant$','',scores$parent)
 stopifnot(!anyDuplicated(meta[c('dataset','sample_id')]))
 j <- match(paste(scores$dataset,scores$sample_id),paste(meta$dataset,meta$sample_id));stopifnot(!anyNA(j))
 d <- data.frame(dataset=scores$dataset,patient_id=meta$subject_id[j],feature=scores$feature,
 value=scores$value,OS_time=meta$OS_time[j],OS_event=meta$OS_event[j])
 g <- read_table(inputs[['gse_data']]);g <- g[g$cohort=='all' & g$source=='malignant',]
 d <- rbind(d,data.frame(dataset='GSE207526',g[c('patient_id','feature','value','OS_time','OS_event')]))
 ####################
 if(include_geo) {
  geo <- read_table(inputs[['geo_data']]);stopifnot(all(geo$source=='malignant'))
  d <- rbind(d,data.frame(dataset='GSE19417',geo[c('patient_id','feature','value','OS_time','OS_event')]))
 }
 ####################
 defs <- read_table(inputs[['definitions']]);gd <- read_table(inputs[['gse_definitions']])
 defs <- defs[defs$feature %in% gd$feature,];d <- d[d$feature %in% defs$feature,]
 stopifnot(nrow(defs)==32L,setequal(defs$feature,gd$feature),
   identical(defs$family,gd$family[match(defs$feature,gd$feature)]),all(is.finite(d$value)),
   !anyDuplicated(d[c('dataset','patient_id','feature')]))
 d$family <- defs$family[match(d$feature,defs$feature)]
 stopifnot(!anyNA(d$family))
 d$eligible <- is.finite(d$OS_time)&d$OS_time>0&d$OS_event %in% c(0,1)
 summary <- do.call(rbind,lapply(datasets,function(k) {
  z <- d[d$dataset==k & d$feature==defs$feature[1],]
  stopifnot(nrow(z)==settings$parent_n[[k]])
  data.frame(dataset=k,n_selected=nrow(z),n_OS=sum(z$eligible),events=sum(z$OS_event[z$eligible]))
 }))
 list(data=d,definitions=defs,summary=summary)
}
fit_cox <- function(d,method,prob=NA_real_) {
 d <- d[d$eligible,]
 r <- data.frame(dataset=d$dataset[1],feature=d$feature[1],family=d$family[1],method=method,
 n=nrow(d),events=sum(d$OS_event),scale=NA_real_,low_n=NA_integer_,high_n=NA_integer_,
 cutpoint=NA_real_,q1=NA_real_,q3=NA_real_,probability=prob,status='estimated',
 logHR=NA_real_,SE=NA_real_,HR=NA_real_,lower_CI=NA_real_,upper_CI=NA_real_,
 P_raw=NA_real_,P_value=NA_real_,PH_P=NA_real_)
 fail <- function(reason) {r$status <- paste('not estimated:',reason);r}
 if(method=='continuous') {
  r$scale <- if(d$family[1]=='State fraction') .1 else sd(d$value)
  if(!is.finite(r$scale)||r$scale<=0)return(fail('constant predictor'))
  d$predictor <- d$value/r$scale
 } else {
  if(method=='median') {r$cutpoint <- median(d$value);high <- d$value>r$cutpoint}
  if(method=='q1q4') {
   q <- quantile(d$value,c(.25,.75),names=FALSE);r$q1 <- q[1];r$q3 <- q[2]
   if(q[1]>=q[2])return(fail('tied quartiles'))
   d <- d[d$value<=q[1]|d$value>=q[2],];high <- d$value>=q[2]
  }
  if(method=='optimal_cut') {r$cutpoint <- unname(quantile(d$value,prob));high <- d$value>=r$cutpoint}
  d$predictor <- as.integer(high);r$n <- nrow(d);r$events <- sum(d$OS_event)
  r$low_n <- sum(!high);r$high_n <- sum(high)
  if(min(r$low_n,r$high_n)<settings$min_group||length(unique(high))<2L||any(tapply(d$OS_event,high,sum)==0))
   return(fail('insufficient group/event support'))
 }
 if(r$n<settings$min_n||r$events<if(method=='continuous')settings$min_events else 5L)
  return(fail('insufficient patients/events'))
 f <- withCallingHandlers(coxph(Surv(OS_time,OS_event)~predictor,data=d,ties='efron',x=TRUE,
   control=coxph.control(eps=1e-12,toler.chol=1e-14,iter.max=50L)),warning=function(w)stop(conditionMessage(w)))
 s <- summary(f);r$logHR <- unname(coef(f));r$SE <- s$coefficients[1,'se(coef)']
 r$HR <- exp(r$logHR);r$lower_CI <- s$conf.int[1,'lower .95'];r$upper_CI <- s$conf.int[1,'upper .95']
 r$P_raw <- r$P_value <- s$coefficients[1,'Pr(>|z|)']
 stopifnot(all(is.finite(c(r$logHR,r$SE,r$P_value))),r$SE>0)
 if(method=='continuous')r$PH_P <- cox.zph(f)$table['predictor','p']
 r
}
fit_meta <- function(d,required=datasets) {
 stopifnot(nrow(d)==length(required),setequal(d$dataset,required))
 r <- data.frame(feature=d$feature[1],family=d$family[1],method=d$method[1],
 k=length(required),n=sum(d$n),events=sum(d$events),status='estimated',logHR=NA_real_,SE=NA_real_,
 HR=NA_real_,lower_CI=NA_real_,upper_CI=NA_real_,P_value=NA_real_,tau2=NA_real_,
 Q=NA_real_,Q_P=NA_real_,I2=NA_real_,all_same_direction=FALSE,all_nominal=FALSE)
 if(any(d$status!='estimated')) {r$status <- 'not estimated: required dataset unsupported';return(r)}
 y <- d$logHR;v <- d$SE^2;k <- length(y)
 objective <- function(tau) {w <- 1/(v+tau);mu <- sum(w*y)/sum(w)
  .5*(sum(log(v+tau))+log(sum(w))+sum(w*(y-mu)^2))}
 bound <- max(1,100*var(y)+max(v))
 opt <- optimize(objective,c(0,bound),tol=1e-12)
 stopifnot(opt$minimum<bound*.99)
 tau <- if(objective(0)<=opt$objective)0 else opt$minimum
 w <- 1/(v+tau);mu <- sum(w*y)/sum(w)
 se <- sqrt(max(1,sum(w*(y-mu)^2)/(k-1))/sum(w))
 critical <- qt(.975,k-1);wf <- 1/v;Q <- sum(wf*(y-sum(wf*y)/sum(wf))^2)
 r$logHR <- mu;r$SE <- se;r$HR <- exp(mu);r$lower_CI <- exp(mu-critical*se);r$upper_CI <- exp(mu+critical*se)
 r$P_value <- 2*pt(-abs(mu/se),k-1);r$tau2 <- tau;r$Q <- Q;r$Q_P <- pchisq(Q,k-1,lower.tail=FALSE)
 typical_v <- (k-1)/(sum(wf)-sum(wf^2)/sum(wf));r$I2 <- 100*tau/(tau+typical_v)
 r$all_same_direction <- all(y>0)||all(y<0);r$all_nominal <- all(d$P_value<.05)
 r
}
run_analysis <- function() {
 a <- prepare_data();fits <- grids <- list();i <- 0L
 for(d in split(a$data,interaction(a$data$dataset,a$data$feature,drop=TRUE))) {
  i <- i+1L;grid <- do.call(rbind,lapply(settings$optimal_grid,function(p)fit_cox(d,'optimal_cut',p)))
  ok <- which(grid$status=='estimated')
  if(length(ok)) {best <- grid[ok[which.min(grid$P_raw[ok])],];best$P_value <- min(1,13*best$P_raw)}
  else {best <- grid[1,];best$status <- 'not estimated: no supported cutpoint'}
  fits[[i]] <- rbind(do.call(rbind,lapply(names(methods),function(m)fit_cox(d,m))),best);grids[[i]] <- grid
 }
 fits <- do.call(rbind,fits);fits$P_BH <- NA_real_
 for(j in split(seq_len(nrow(fits)),interaction(fits$dataset,fits$method,drop=TRUE)))
  fits$P_BH[j] <- p.adjust(fits$P_value[j],'BH',n=32L)
 primary <- fits[fits$method %in% names(methods),]
 groups <- split(primary,interaction(primary$feature,primary$method,drop=TRUE))
 pooled <- do.call(rbind,lapply(groups,fit_meta));pooled$P_BH <- NA_real_
 for(j in split(seq_len(nrow(pooled)),pooled$method)) pooled$P_BH[j] <- p.adjust(pooled$P_value[j],'BH',n=32L)
 loo <- do.call(rbind,lapply(datasets,function(k)do.call(rbind,lapply(groups,function(d) {
  r <- fit_meta(d[d$dataset!=k,],setdiff(datasets,k));r$excluded_dataset <- k;r
 }))))
 list(prepared=a,fits=fits,grid=do.call(rbind,grids),pooled=pooled,loo=loo)
}
contract <- list(manifest=manifest,settings=settings,versions=versions,
 code=lapply(list(prepare_data,fit_cox,fit_meta,run_analysis),function(f)paste(deparse(body(f)),collapse='\n')))
cache_path <- file.path(out,'intermediate/meta_analysis.rds')
if(file.exists(cache_path)&&!force) {cache <- readRDS(cache_path);stopifnot(identical(cache$contract,contract))} else {if(replot)stop('Matching cache required');cache <- c(list(contract=contract),run_analysis());saveRDS(cache,cache_path)}
write_table(manifest,'input_manifest');write_table(cache$prepared$summary,'dataset_summary')
write_table(cache$prepared$data,'model_data');write_table(cache$prepared$definitions,'feature_definitions')
write_table(cache$fits,'cohort_estimates');write_table(cache$grid,'cutpoint_grid')
write_table(cache$pooled,'meta_estimates');write_table(cache$loo,'leave_one_dataset_out')
write.csv(versions,file.path(out,'logs/package_versions.csv'),row.names=FALSE)

colours <- c(OCCAMS='#0072B2',TCGA='#D55E00',GSE207526='#009E73',GSE19417='#CC79A7',Pooled='#222222')
families <- c(MP='Metaprogrammes','State MP union'='MP-union\nstate scores','State DGE'='DGE\nstate scores',
 'State fraction'='Within-malignant\nRNA fractions')
defs <- cache$prepared$definitions
format_feature <- function(x) vapply(defs$label[match(x,defs$feature)],function(z)paste(strwrap(z,47),collapse='\n'),character(1))
count_label <- function(x) paste(unique(range(x)),collapse='–')
plot_forest <- function(features,method,pooled=TRUE) {
 d <- cache$fits[cache$fits$feature %in% features & cache$fits$method==method,]
 if(pooled) {
  z <- cache$pooled[cache$pooled$feature %in% features & cache$pooled$method==method,]
  z$dataset <- 'Pooled';cols <- c('dataset','feature','family','method','n','events','status','logHR','SE','HR','lower_CI','upper_CI','P_value','P_BH')
  d <- rbind(d[cols],z[cols])
 }
 labels <- setNames(vapply(unique(d$dataset),function(k) {z <- d[d$dataset==k,]
  paste0(k,'\nn=',count_label(z$n),', events=',count_label(z$events))},character(1)),unique(d$dataset))
 d$panel <- factor(d$dataset,levels=names(colours),labels=labels[names(colours)])
 d$feature_plot <- factor(d$feature,levels=rev(features));d$family <- factor(d$family,levels=names(families))
 d$evidence <- factor(ifelse(d$P_BH<.05,'FDR',ifelse(d$P_value<.05,'Nominal','NS')),levels=c('NS','Nominal','FDR'))
 ####################
 # One common zoom is used for every dataset and pooled column on each page.
 # Expression/fraction continuous models share HR 0.5-2; binary contrasts share
 # HR 0.25-4. Extreme intervals retain arrows and complete source-table values.
 d$axis_span <- if(method=='continuous')1 else 2
 e <- d[d$status=='estimated',]
 e$effect <- pmax(-e$axis_span,pmin(e$axis_span,log2(e$HR)))
 e$left <- pmax(-e$axis_span,pmin(e$axis_span,log2(e$lower_CI)))
 e$right <- pmax(-e$axis_span,pmin(e$axis_span,log2(e$upper_CI)))
 row_labels <- families
 if(method=='continuous')row_labels <- c(MP='Metaprogrammes\nper 1 SD',
  'State MP union'='MP-union scores\nper 1 SD','State DGE'='DGE scores\nper 1 SD',
  'State fraction'='RNA fractions\nper 10 percentage points')
 ####################
 ggplot(d,aes(y=feature_plot,colour=dataset)) +
  geom_blank(aes(x=-axis_span))+geom_blank(aes(x=axis_span))+
  geom_vline(xintercept=0,linetype=2,colour='grey45',linewidth=.5)+
  geom_errorbar(data=e,aes(xmin=left,xmax=right,x=effect),orientation='y',width=.12,linewidth=.5)+
  geom_point(data=e,aes(x=effect,shape=evidence),size=2.4,show.legend=TRUE)+
  geom_text(data=e[log2(e$lower_CI) < -e$axis_span,],aes(x=-axis_span,label='\u2039'),size=5,show.legend=FALSE)+
  geom_text(data=e[log2(e$upper_CI)>e$axis_span,],aes(x=axis_span,label='\u203a'),size=5,show.legend=FALSE)+
  geom_text(data=d[d$status!='estimated',],aes(x=0,label='NE'),colour='grey50',size=3,show.legend=FALSE)+
  facet_grid(family~panel,scales='free_y',space='free_y',drop=TRUE,labeller=labeller(family=row_labels))+
  scale_colour_manual(values=colours,guide='none')+scale_shape_manual(values=c(NS=1,Nominal=16,FDR=18),drop=FALSE,name=NULL)+
  scale_y_discrete(labels=format_feature)+
  scale_x_continuous(expand=expansion(mult=.02),
   breaks=function(limits) {span<-max(abs(limits))/1.04;c(-span,0,span)},
   labels=function(x)format(signif(2^x,3),trim=TRUE))+
  labs(x=if(method=='continuous') {if(all(features %in% defs$feature[defs$family=='MP']))'Hazard ratio per 1 SD (95% CI)' else 'Hazard ratio (95% CI)'} else 'Hazard ratio: high vs low (95% CI)',
   y=NULL,title=paste(if(all(features %in% defs$feature[defs$family=='MP']))'MPs' else 'States','|',
    if(method=='optimal_cut')'Optimal cut' else methods[[method]]))+
  theme_classic(base_size=14)+theme(legend.position='bottom',strip.background=element_blank(),
   panel.spacing.x=grid::unit(36,'pt'),
   strip.text.x=element_text(size=12),strip.text.y=element_text(size=12,angle=0),axis.text.y=element_text(size=12,colour='black'),
   axis.text.x=element_text(size=12,colour='black'),plot.title=element_text(size=18,face='bold'),plot.margin=margin(10,12,10,10))
}
index <- list()
for(kind in c('MP','State')) {
 features <- defs$feature[if(kind=='MP')defs$family=='MP' else defs$family!='MP']
 file <- paste0('figures/Auto_final/',if(kind=='MP')'01_mp' else '02_state','_meta_analysis.pdf')
 cairo_pdf(file.path(out,file),width=20,height=12,onefile=TRUE)
 for(m in names(methods))print(plot_forest(features,m))
 dev.off();index[[length(index)+1L]] <- data.frame(path=file,pages=3L,purpose=paste(kind,'all-patient',length(datasets),'dataset OS meta-analysis'))
}
file <- 'figures/Auto_supporting/01_optimal_cut_exploratory.pdf'
cairo_pdf(file.path(out,file),width=20,height=12,onefile=TRUE)
for(kind in c('MP','State'))print(plot_forest(defs$feature[if(kind=='MP')defs$family=='MP' else defs$family!='MP'],'optimal_cut',FALSE))
dev.off();index[[3]] <- data.frame(path=file,pages=2L,purpose='Dataset-specific searched cutpoints; no pooled estimate')
write_table(do.call(rbind,index),'figure_index')
checks <- data.frame(check=c('all parent patients represented','cohort fits','96 meta rows','leave-one-out rows'),
 passed=c(nrow(cache$prepared$data)==sum(settings$parent_n)*32L,nrow(cache$fits)==length(datasets)*128L,
 nrow(cache$pooled)==96L,nrow(cache$loo)==length(datasets)*96L))
stopifnot(all(checks$passed));write.csv(checks,file.path(out,'reports/validation.csv'),row.names=FALSE)
summary <- cache$pooled[c('feature','family','method','status','HR','P_value','P_BH','I2','all_same_direction','all_nominal')]
dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE)
write.csv(summary,file.path(SCREF_SUMMARY_DIR,'malignant_survival_meta.csv'),row.names=FALSE)
writeLines(c('# All-patient malignant OS meta-analysis',
 paste('Parent samples:',paste(names(settings$parent_n),settings$parent_n,collapse='; '),'. No treatment, surgery or specimen subgroup restrictions are applied.'),
 paste(apply(cache$prepared$summary,1,function(x)paste(x,collapse=' | ')),collapse='\n'),
 'Cohort-specific Cox estimates use Efron ties; continuous expression HRs are per dataset SD, fractions per 10 percentage points. Scores are reused without rescoring.',
 paste('REML random-effects pooling uses modified Hartung-Knapp standard errors bounded below by the ordinary random-effects SE; t tests use',length(datasets)-1L,'degrees of freedom. All configured datasets are required for each pooled estimate.'),
 'Expression HRs are per one dataset-specific score SD; fraction HRs are per 10 percentage points. These units are not equivalent. All datasets and pooled estimates share the same axis on each page: HR 0.5-2 for continuous models, HR 0.25-4 for binary models. Boundary arrows mark intervals extending outside these display limits.',
 'Three pages per final PDF: continuous (primary), median and Q1/Q4. Open circles: NS; solid circles: nominal P<0.05; diamonds: BH q<0.05 across 32 features per method. NE means an explicitly unsupported model. Boundary arrows identify intervals extending beyond the display.',
 'Optimal-cut results are cohort-specific and exploratory, with a 13-threshold P correction before BH. Selected HRs/SEs remain biased by selection and are not meta-analysed.',
 'A pooled association is not independent validation. Direction agreement and all-cohort nominal significance are tabulated separately. MP12 and the stress MP union can duplicate the same signature.',
 'Few datasets are available; heterogeneity is imprecise. OS origins, treatment exposure and specimen timing differ. Mixing these patients follows the requested all-patient estimand; unadjusted pooled associations do not control clinical confounding.',
 if(include_geo)'GSE19417 is included as exploratory microarray evidence from tumour-channel intensities; it is not RNA-seq. Leave-one-out results excluding GSE19417 give the RNA-seq-only sensitivity comparison. Microarray posterior CV is conditional on an artificial mixture scale.' else 'Current figures use only the three RNA-seq datasets. GEO inclusion requires an explicit SCREF_INCLUDE_GEO=TRUE run with all required prepared inputs; no automatic cohort substitution is permitted.',
 '[MPs](../figures/Auto_final/01_mp_meta_analysis.pdf) | [States](../figures/Auto_final/02_state_meta_analysis.pdf)',
 '[Pooled statistics](../tables/meta_estimates.csv) | [Dataset estimates](../tables/cohort_estimates.csv) | [Sensitivity](../tables/leave_one_dataset_out.csv)'),file.path(out,'reports/current_results.md'))
message('META_COMPLETE: ',nrow(cache$pooled),' pooled specifications; ',sum(settings$parent_n),' parent patients; 8 PDF pages')
