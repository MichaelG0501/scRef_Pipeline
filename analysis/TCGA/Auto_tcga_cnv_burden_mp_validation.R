####################
# Analysis registry:
# Status: active
# Script: analysis/TCGA/Auto_tcga_cnv_burden_mp_validation.R
# Description: within-tumour single-cell InferCNA burden/MP contrasts and between-tumour TCGA validation.
# Methodology: analysis/methodology/TCGA/tcga_cnv_burden_mp_validation_methodology.md
# Inputs: ref_outs/Auto_malignant_subclone_mp/threshold_1/Auto_malignant_subclone_cells.csv;
#   ref_outs/meta_full_epi.rds;
#   ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_sample_state_eligibility.csv;
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds;
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds;
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv;
#   ref_outs/by_samples/<first audited sample>/<sample>_outs.rds (InferCNA score reproduction);
#   ref_outs/TCGA/cna_gene_event_validation/intermediate/Auto_targeted_cna_results.rds;
#   ref_outs/TCGA/cna_gene_event_validation/reports/Auto_manual_infercna_provenance.csv;
#   ref_outs/TCGA/cna_recurrent_event_validation/intermediate/source_data/esca_tcga_gdc_segments.seg.
# Outputs: ref_outs/TCGA/cnv_burden_mp_validation/
#   intermediate/Auto_burden_sources.rds; Auto_burden_results.rds; Auto_burden_figures.rds;
#   tables/Auto_cell_burden_mp.csv; Auto_sample_burden_groups.csv; Auto_sample_mp_contrasts.csv;
#     Auto_tcga_burden_mp.csv; Auto_burden_effect_tests.csv; Auto_burden_continuous_tests.csv;
#     Auto_burden_sensitivity.csv; Auto_burden_summary.csv; Auto_burden_boxplot_tests.csv;
#     Auto_burden_tertile_groups.csv; Auto_burden_tertile_sample_means.csv; Auto_burden_tertile_tests.csv;
#   figures/Auto_burden_mp_cross_cohort_effects.pdf;
#     Auto_burden_mp_continuous_rectangles.pdf (three pages: trends, pooled scATLAS boxes, TCGA boxes);
#   reports/Auto_burden_validation.txt; logs/Auto_burden_run_summary.rds/.txt;
#     logs/Auto_burden_replot_run_summary.rds/.txt;
#   updates/new_updates/summaries/Auto_burden_summary.csv.
# Cache: SCREF_FORCE_REBUILD=TRUE rebuilds fingerprinted source cache.
#   SCREF_REPLOT_ONLY=TRUE renders saved grobs without loading cells or fitting tests.
# Run: /opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_cnv_burden_mp_validation.sh; conda: dmtcp.
####################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)
source("analysis/shared/scRef_config.R")
source("analysis/shared/scRef_helpers.R")
source("analysis/TCGA/Auto_tcga_workflow_helpers.R")
tcga_execution_guard()
set.seed(42)
script_path <- file.path(SCREF_ANALYSIS_DIR,"TCGA","Auto_tcga_cnv_burden_mp_validation.R")
out_dir <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cnv_burden_mp_validation")
tiers <- ensure_output_tiers(out_dir)
figure_path <- file.path(tiers[["intermediate"]],"Auto_burden_figures.rds")
summary_path <- file.path(tiers[["tables"]],"Auto_burden_summary.csv")
if (tcga_flag("SCREF_REPLOT_ONLY")) {
  tcga_replot(figure_path,summary_path,file.path(tiers[["logs"]],"Auto_burden_replot_run_summary.rds"))
  quit(save="no")
}
cell_path <- file.path(SCREF_REF_OUTS_DIR,"Auto_malignant_subclone_mp","threshold_1","Auto_malignant_subclone_cells.csv")
parent_path <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cna_gene_event_validation","intermediate","Auto_targeted_cna_results.rds")
provenance_path <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cna_gene_event_validation","reports","Auto_manual_infercna_provenance.csv")
segment_path <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cna_recurrent_event_validation","intermediate","source_data","esca_tcga_gdc_segments.seg")
mp_sets <- tcga_validate_gene_sets(readRDS(SCREF_MP_GENES_RDS))
mp_order <- names(mp_sets)
mp_labels <- setNames(label_mps(mp_order),mp_order)
study_path <- file.path(SCREF_REF_OUTS_DIR,"Metaprogrammes_Results","centred","state_markers","Auto_five_state_sample_state_eligibility.csv")
inputs <- c(cell_path,study_path,parent_path,provenance_path,segment_path,SCREF_META_FULL_EPI_RDS,SCREF_MP_UCELL_RDS,SCREF_MP_GENES_RDS)
tcga_require_files(inputs)
cache <- tcga_cache(file.path(tiers[["intermediate"]],"Auto_burden_sources.rds"),inputs,function() {
  manual <- read.csv(cell_path,check.names=FALSE)
  provenance <- read.csv(provenance_path)
  stopifnot(all(provenance$source=="manual InferCNA; not Numbat"),all(provenance$abs_difference<1e-10),
            !anyDuplicated(manual$cell),all(is.finite(manual$cna_signal)),all(manual$cna_signal>=0))
  meta <- readRDS(SCREF_META_FULL_EPI_RDS)
  scores <- readRDS(SCREF_MP_UCELL_RDS)
  stopifnot("orig.ident" %in% names(meta),all(manual$cell %in% rownames(meta)),
            all(manual$cell %in% rownames(scores)),all(mp_order %in% colnames(scores)))
  meta <- meta[manual$cell,,drop=FALSE]
  study_map <- read.csv(study_path) |> select(sample,study) |> distinct()
  stopifnot(!anyDuplicated(study_map$sample),all(manual$sample %in% study_map$sample))
  study <- as.character(study_map$study[match(manual$sample,study_map$sample)])
  stopifnot(identical(as.character(meta$orig.ident),as.character(manual$sample)),
            !anyNA(study),all(nzchar(study)))
  score_mat <- as.matrix(scores[manual$cell,mp_order,drop=FALSE])
  stopifnot(all(is.finite(score_mat)))
  cells <- data.frame(cell=manual$cell,sample=manual$sample,study=study,burden=manual$cna_signal,
                      score_mat,check.names=FALSE)
  first_sample <- sort(unique(cells$sample))[1]
  outs_path <- file.path(SCREF_REF_OUTS_DIR,"by_samples",first_sample,paste0(first_sample,"_outs.rds"))
  tcga_require_files(outs_path)
  outs <- readRDS(outs_path)
  ids <- head(cells$cell[cells$sample==first_sample],20)
  observed <- colMeans(outs[,ids,drop=FALSE]^2)
  max_error <- max(abs(observed-cells$burden[match(ids,cells$cell)]))
  stopifnot(max_error<1e-10)
  rm(outs,meta,scores,score_mat)
  parent <- readRDS(parent_path)
  bulk <- parent$tcga_arm |> select(sample_barcode,sample_key,seg_sample_id,all_of(mp_order)) |> distinct()
  stopifnot(!anyDuplicated(bulk$sample_key),!anyDuplicated(bulk$seg_sample_id))
  seg <- fread(segment_path,data.table=FALSE)
  names(seg) <- make.names(names(seg))
  stopifnot(all(c("ID","chrom","loc.start","loc.end","seg.mean") %in% names(seg)))
  seg <- seg |> transmute(seg_sample_id=as.character(ID),chr=gsub("^chr","",as.character(chrom)),
                         start=as.numeric(loc.start),end=as.numeric(loc.end),signal=as.numeric(seg.mean)) |>
    filter(seg_sample_id %in% bulk$seg_sample_id,chr %in% as.character(1:22)) |>
    arrange(seg_sample_id,chr,start)
  stopifnot(all(is.finite(seg$signal)),all(seg$end>=seg$start),all(seg$start>=0))
  overlap <- seg |> group_by(seg_sample_id,chr) |> summarise(bad=any(start<lag(end),na.rm=TRUE),.groups="drop")
  if (any(overlap$bad)) stop("Overlapping autosomal TCGA segments; do not double-count burden.")
  # The two audited shared endpoints belong once to the earlier segment.
  seg <- seg |> group_by(seg_sample_id,chr) |>
    mutate(shared_endpoint=!is.na(lag(end))&start==lag(end),start=start+as.integer(shared_endpoint)) |> ungroup()
  seg$bp <- seg$end-seg$start+1
  stopifnot(all(seg$bp>0))
  burden <- seg |> group_by(seg_sample_id) |> summarise(covered_bp=sum(bp),
    fga=sum(bp[abs(signal)>=0.20])/sum(bp),fga_010=sum(bp[abs(signal)>=0.10])/sum(bp),
    fga_030=sum(bp[abs(signal)>=0.30])/sum(bp),rms=sqrt(weighted.mean(signal^2,bp)),.groups="drop")
  bulk <- bulk |> left_join(burden,by="seg_sample_id")
  stopifnot(!anyNA(bulk$fga),all(bulk$fga>=0&bulk$fga<=1))
  list(cells=cells,bulk=bulk,score_audit=data.frame(sample=first_sample,n_cells=length(ids),max_error=max_error,input=outs_path))
})
cells <- cache$value$cells
bulk <- cache$value$bulk
sample_groups <- cells |> group_by(sample,study) |> summarise(n_cells=n(),cutoff=median(burden),
  mean_burden=mean(burden),low_n=sum(burden<=median(burden)),high_n=sum(burden>median(burden)),.groups="drop") |>
  mutate(eligible=low_n>=20L&high_n>=20L,reason=ifelse(eligible,"included","fewer than 20 cells in a median-defined group"))
stopifnot(!anyDuplicated(sample_groups$sample))
cells <- cells |> left_join(sample_groups |> select(sample,cutoff,eligible),by="sample") |>
  mutate(group=factor(ifelse(burden>cutoff,"High","Low"),levels=c("Low","High")))
long <- cells |> filter(eligible) |> pivot_longer(all_of(mp_order),names_to="mp",values_to="score") |>
  group_by(sample,mp) |> mutate(z=(score-mean(score))/sd(score)) |> ungroup()
means <- long |> group_by(sample,study,mp,group) |> summarise(raw=mean(score),z=mean(z),n_cells=n(),.groups="drop")
contrasts <- means |> select(sample,study,mp,group,z,raw) |>
  pivot_wider(names_from=group,values_from=c(z,raw)) |> mutate(effect=z_High-z_Low,raw_difference=raw_High-raw_Low)
bulk$group <- factor(ifelse(bulk$fga>median(bulk$fga),"High","Low"),levels=c("Low","High"))
bulk_long <- bulk |> pivot_longer(all_of(mp_order),names_to="mp",values_to="score") |>
  group_by(mp) |> mutate(z=(score-mean(score))/sd(score)) |> ungroup()
boot_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)<5L) return(c(low=NA_real_,high=NA_real_))
  quantile(replicate(1000,mean(sample(x,length(x),replace=TRUE))),c(.025,.975),names=FALSE)
}
safe_signed <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)<5L || !any(x!=0)) return(NA_real_)
  suppressWarnings(wilcox.test(x,mu=0,exact=FALSE)$p.value)
}
spearman <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if (sum(ok)<5L || length(unique(x[ok]))<2L || length(unique(y[ok]))<2L) return(c(rho=NA_real_,p=NA_real_))
  test <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  c(rho=unname(test$estimate),p=test$p.value)
}
studies <- sort(unique(sample_groups$study))
sc_effects <- bind_rows(lapply(c("All studies",studies),function(scope) {
  df <- if (scope=="All studies") contrasts else contrasts[contrasts$study==scope,]
  bind_rows(lapply(mp_order,function(mp) {
    x <- df$effect[df$mp==mp];x <- x[is.finite(x)];ci <- boot_mean(x)
    data.frame(dataset="scATLAS within sample",scope=scope,mp=mp,n=length(x),
      effect=if(length(x))mean(x) else NA_real_,low=ci[1],high=ci[2],p_value=safe_signed(x))
  }))
}))
tcga_effects <- bind_rows(lapply(mp_order,function(mp) {
  df <- bulk_long[bulk_long$mp==mp,]
  lo <- df$z[df$group=="Low"];hi <- df$z[df$group=="High"]
  stopifnot(length(lo)>=5L,length(hi)>=5L,all(is.finite(c(lo,hi))))
  ci <- quantile(replicate(1000,mean(sample(hi,length(hi),replace=TRUE))-mean(sample(lo,length(lo),replace=TRUE))),c(.025,.975),names=FALSE)
  data.frame(dataset="TCGA between samples",scope="TCGA EAC",mp=mp,n=nrow(df),effect=mean(hi)-mean(lo),
    low=ci[1],high=ci[2],p_value=suppressWarnings(wilcox.test(hi,lo,exact=FALSE)$p.value))
}))
effects <- bind_rows(sc_effects,tcga_effects) |> group_by(dataset,scope) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
cell_rhos <- long |> group_by(sample,study,mp) |> summarise(rho=spearman(burden,score)[["rho"]],.groups="drop")
continuous <- bind_rows(lapply(c("All studies",studies),function(scope) {
  df <- if (scope=="All studies") cell_rhos else cell_rhos[cell_rhos$study==scope,]
  bind_rows(lapply(mp_order,function(mp) {
    x <- df$rho[df$mp==mp];x <- x[is.finite(x)]
    data.frame(dataset="scATLAS within sample",scope=scope,mp=mp,n=length(x),
               rho=if(length(x))mean(x) else NA_real_,p_value=safe_signed(x))
  }))
}),lapply(mp_order,function(mp) {
  stat <- spearman(bulk$fga,bulk[[mp]])
  data.frame(dataset="TCGA between samples",scope="TCGA EAC",mp=mp,n=nrow(bulk),rho=stat[["rho"]],p_value=stat[["p"]])
})) |> group_by(dataset,scope) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
sensitivity <- bind_rows(lapply(c("fga_010","fga","fga_030","rms"),function(metric) {
  bind_rows(lapply(mp_order,function(mp) {
    stat <- spearman(bulk[[metric]],bulk[[mp]])
    data.frame(dataset="TCGA between samples",metric=metric,mp=mp,rho=stat[["rho"]],p_value=stat[["p"]])
  }))
})) |> group_by(metric) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
####################
# Boxplot inference uses paired raw sample-group means, independently of standardized forest tests.
boxplot_tests <- bind_rows(lapply(mp_order,function(mp) {
  df <- contrasts |> filter(.data$mp==.env$mp)
  data.frame(dataset="scATLAS",mp=mp,n=nrow(df),mean_difference=mean(df$raw_difference),
    test="paired Wilcoxon signed-rank on raw sample means",p_value=safe_signed(df$raw_difference))
}),lapply(mp_order,function(mp) {
  df <- bulk_long |> filter(.data$mp==.env$mp)
  hi <- df$score[df$group=="High"];lo <- df$score[df$group=="Low"]
  data.frame(dataset="TCGA EAC",mp=mp,n=nrow(df),mean_difference=mean(hi)-mean(lo),
    test="unpaired Wilcoxon rank-sum",p_value=suppressWarnings(wilcox.test(hi,lo,exact=FALSE)$p.value))
})) |> group_by(dataset) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
# Outer-tertile sensitivity: predeclared quantiles, ties kept together, middle excluded.
tertile_groups <- cells |> group_by(sample,study) |> summarise(q1=quantile(burden,1/3,type=7),q2=quantile(burden,2/3,type=7),
  low_n=sum(burden<=q1),middle_n=sum(burden>q1&burden<=q2),high_n=sum(burden>q2),.groups="drop") |>
  mutate(eligible=low_n>=20L&high_n>=20L)
tertile_means <- long |> left_join(tertile_groups |> select(sample,q1,q2,tertile_eligible=eligible),by="sample") |>
  filter(tertile_eligible) |> mutate(tertile_group=ifelse(burden<=q1,"Low",ifelse(burden>q2,"High","Middle"))) |>
  filter(tertile_group!="Middle") |> group_by(sample,study,mp,tertile_group) |>
  summarise(raw=mean(score),n_cells=n(),.groups="drop")
tertile_contrasts <- tertile_means |> select(sample,study,mp,tertile_group,raw) |>
  pivot_wider(names_from=tertile_group,values_from=raw) |> mutate(difference=High-Low)
tcga_tertile_cuts <- quantile(bulk$fga,c(1/3,2/3),type=7,names=FALSE)
bulk_tertile <- bulk_long |> mutate(tertile_group=ifelse(fga<=tcga_tertile_cuts[1],"Low",ifelse(fga>tcga_tertile_cuts[2],"High","Middle")))
tertile_tests <- bind_rows(lapply(c("All studies",studies),function(scope) {
  df <- if(scope=="All studies") tertile_contrasts else tertile_contrasts[tertile_contrasts$study==scope,]
  bind_rows(lapply(mp_order,function(mp) {
    x <- df$difference[df$mp==mp]
    data.frame(dataset="scATLAS",scope=scope,mp=mp,n=length(x),mean_difference=if(length(x))mean(x) else NA_real_,
      p_value=safe_signed(x))
  }))
}),lapply(mp_order,function(mp) {
  df <- bulk_tertile |> filter(.data$mp==.env$mp)
  hi <- df$score[df$tertile_group=="High"];lo <- df$score[df$tertile_group=="Low"]
  stopifnot(length(hi)>=5L,length(lo)>=5L)
  data.frame(dataset="TCGA EAC",scope="TCGA EAC",mp=mp,n=length(hi)+length(lo),mean_difference=mean(hi)-mean(lo),
    p_value=suppressWarnings(wilcox.test(hi,lo,exact=FALSE)$p.value))
})) |> group_by(dataset,scope) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
tertile_groups_export <- bind_rows(mutate(tertile_groups,dataset="scATLAS"),
 data.frame(sample="TCGA EAC cohort",study="TCGA EAC",q1=tcga_tertile_cuts[1],q2=tcga_tertile_cuts[2],
   low_n=sum(bulk$fga<=tcga_tertile_cuts[1]),middle_n=sum(bulk$fga>tcga_tertile_cuts[1]&bulk$fga<=tcga_tertile_cuts[2]),
   high_n=sum(bulk$fga>tcga_tertile_cuts[2]),eligible=TRUE,dataset="TCGA EAC"))
write.csv(boxplot_tests,file.path(tiers[["tables"]],"Auto_burden_boxplot_tests.csv"),row.names=FALSE)
write.csv(tertile_groups_export,file.path(tiers[["tables"]],"Auto_burden_tertile_groups.csv"),row.names=FALSE)
write.csv(tertile_means,file.path(tiers[["tables"]],"Auto_burden_tertile_sample_means.csv"),row.names=FALSE)
write.csv(tertile_tests,file.path(tiers[["tables"]],"Auto_burden_tertile_tests.csv"),row.names=FALSE)
####################
for (name in c("cells","sample_groups","contrasts","bulk","effects","continuous","sensitivity")) {
  filename <- c(cells="Auto_cell_burden_mp.csv",sample_groups="Auto_sample_burden_groups.csv",contrasts="Auto_sample_mp_contrasts.csv",bulk="Auto_tcga_burden_mp.csv",effects="Auto_burden_effect_tests.csv",continuous="Auto_burden_continuous_tests.csv",sensitivity="Auto_burden_sensitivity.csv")[[name]]
  write.csv(get(name),file.path(tiers[["tables"]],filename),row.names=FALSE)
}
results <- list(sample_groups=sample_groups,means=means,contrasts=contrasts,bulk_long=bulk_long,effects=effects,
                cell_rhos=cell_rhos,continuous=continuous,sensitivity=sensitivity,score_audit=cache$value$score_audit,
                boxplot_tests=boxplot_tests,tertile_groups=tertile_groups_export,tertile_means=tertile_means,tertile_tests=tertile_tests)
saveRDS(results,file.path(tiers[["intermediate"]],"Auto_burden_results.rds"))
forest <- function(df,title) {
  df <- df |> mutate(mp_label=factor(mp_labels[mp],levels=rev(mp_labels[mp_order])),star=ifelse(is.finite(fdr)&fdr<.05,"*",""))
  ggplot(df,aes(effect,mp_label,colour=dataset)) + geom_vline(xintercept=0,colour="grey55",linetype=2) +
    geom_segment(aes(x=low,xend=high,yend=mp_label),linewidth=1,na.rm=TRUE) + geom_point(size=3,na.rm=TRUE) +
    geom_text(aes(label=star),nudge_y=.23,size=6,show.legend=FALSE) + facet_wrap(~scope,nrow=1) +
    scale_colour_manual(values=c("scATLAS within sample"="#2166AC","TCGA between samples"="#B2182B")) +
    labs(title=title,x="Standardized MP difference: high minus low burden (95% bootstrap CI)",y=NULL) +
    theme_classic(base_size=19)+theme(axis.text.y=element_text(size=14,colour="black"),strip.text=element_text(face="bold",size=19),
      legend.position="bottom",legend.title=element_blank(),plot.title=element_text(face="bold",size=23),plot.margin=margin(18,24,18,18))
}
tcga_save_figure("Auto_burden_mp_cross_cohort_effects.pdf",list(forest(effects |> filter(scope %in% c("All studies","TCGA EAC")),
 "CNV burden versus MPs: within-tumour scATLAS and between-tumour TCGA")),width=24,height=12)
####################
# Separate pooled scATLAS, individual studies and TCGA with explicit facet gutters and borders.
heat <- continuous |> mutate(mp_label=factor(mp_labels[mp],levels=rev(mp_labels[mp_order])),
  scope=factor(scope,levels=c("All studies",studies,"TCGA EAC")),
  section=factor(ifelse(scope=="All studies","Pooled scATLAS",ifelse(scope=="TCGA EAC","TCGA EAC","Individual studies")),
    levels=c("Pooled scATLAS","Individual studies","TCGA EAC")),
  star=ifelse(is.finite(fdr)&fdr<.05,"*",""))
p <- ggplot(heat,aes(scope,mp_label,fill=rho))+geom_tile(colour="grey85",linewidth=.4)+geom_text(aes(label=star),size=7)+
  facet_grid(.~section,scales="free_x",space="free_x",
    labeller=labeller(section=c("Pooled scATLAS"="Pooled\nscATLAS","Individual studies"="Individual studies","TCGA EAC"="TCGA EAC")))+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-.5,.5),oob=scales::squish,na.value="white",name="Spearman rho")+
  labs(title="Continuous burden–MP trends: equal-sample mean rho in scATLAS; tumour rho in TCGA",x=NULL,y=NULL)+
  theme_classic(base_size=20)+theme(axis.text.x=element_text(angle=45,hjust=1,size=15),axis.text.y=element_text(size=14),
    strip.text=element_text(size=20,face="bold"),strip.background=element_rect(fill="grey95",colour="black"),
    panel.spacing.x=grid::unit(1.3,"lines"),panel.border=element_rect(fill=NA,colour="black",linewidth=1.1),axis.line=element_blank(),
    plot.title=element_text(face="bold",size=22),plot.margin=margin(18,30,18,20))
# Clinical-style layout: all 17 MPs along x, dodged low/high boxes and sample-level dots.
clinical_boxplot <- function(df,stats,title,subtitle,y_label) {
  df <- df |> mutate(mp=factor(mp,levels=mp_order),group=factor(group,levels=c("Low","High")))
  counts <- df |> distinct(sample,group) |> count(group,name="n")
  legend_labels <- setNames(paste0(counts$group," (n=",counts$n,")"),counts$group)
  y_span <- diff(range(df$value));y_offset <- max(.01,.08*y_span)
  annotations <- df |> group_by(mp) |> summarise(y=max(value)+y_offset,.groups="drop") |>
    left_join(stats |> select(mp,fdr),by="mp") |>
    mutate(label=case_when(fdr<.001~"***",fdr<.01~"**",fdr<.05~"*",TRUE~""))
  axis_labels <- setNames(paste0(mp_order,"\n",unname(SCREF_MP_DESCRIPTIONS[mp_order])),mp_order)
  ggplot(df,aes(mp,value,fill=group,colour=group))+
    geom_boxplot(position=position_dodge(width=.75),width=.6,outlier.shape=NA,alpha=.8,linewidth=.6,colour="black")+
    geom_point(position=position_jitterdodge(jitter.width=.15,dodge.width=.75),alpha=.7,size=1.5,stroke=0,show.legend=FALSE)+
    geom_text(data=annotations |> filter(label!=""),aes(x=mp,y=y,label=label),inherit.aes=FALSE,size=7,fontface="bold")+
    scale_fill_manual(values=c(Low="#2166AC",High="#B2182B"),labels=legend_labels,drop=FALSE)+
    scale_colour_manual(values=c(Low="#2166AC",High="#B2182B"),guide="none",drop=FALSE)+
    scale_x_discrete(labels=axis_labels)+scale_y_continuous(expand=expansion(mult=c(.02,.18)))+
    labs(title=title,subtitle=subtitle,x=NULL,y=y_label,fill="CNV burden")+
    theme_classic(base_size=20)+theme(plot.title=element_text(face="bold",size=24),
      plot.subtitle=element_text(size=16,colour="grey35"),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=14,colour="black"),
      axis.line.x=element_blank(),legend.position="top",legend.title=element_text(face="bold"),plot.margin=margin(18,28,18,18))
}
sc_boxes <- means |> transmute(sample,mp,group,value=raw)
tcga_boxes <- bulk_long |> transmute(sample=sample_barcode,mp,group,value=score)
p_sc <- clinical_boxplot(sc_boxes,boxplot_tests |> filter(dataset=="scATLAS"),
 "All studies: scATLAS MP activity by within-sample CNV burden",
 "Median split within each sample; dots = sample/group means; paired Wilcoxon signed-rank; stars = BH FDR across 17 MPs.",
 "Mean sample-group UCell score")
p_tcga <- clinical_boxplot(tcga_boxes,boxplot_tests |> filter(dataset=="TCGA EAC"),
 "TCGA EAC: MP activity by tumour CNV burden",
 "Cohort median FGA split; dots = individual tumours; Wilcoxon rank-sum; stars = BH FDR across 17 MPs.","Bulk GSVA MP score")
results$display_sources <- list(heat=heat,sc_boxes=sc_boxes,tcga_boxes=tcga_boxes)
saveRDS(results,file.path(tiers[["intermediate"]],"Auto_burden_results.rds"))
tcga_save_figure("Auto_burden_mp_continuous_rectangles.pdf",list(p,p_sc,p_tcga),width=25,height=12)
####################
saveRDS(list(figures=tcga_figure_bundle),figure_path)
summary <- data.frame(metric=c("manual_cells","manual_samples","included_samples","excluded_samples","studies","tcga_eac_samples", "tcga_fga_cutoff","tcga_low_n","tcga_high_n","infercna_score_max_error","sc_pooled_fdr_mps","tcga_fdr_mps"),
 value=c(nrow(cells),nrow(sample_groups),sum(sample_groups$eligible),sum(!sample_groups$eligible),length(studies),nrow(bulk),median(bulk$fga),sum(bulk$group=="Low"),sum(bulk$group=="High"),cache$value$score_audit$max_error,
 sum(effects$fdr[effects$scope=="All studies"]<.05,na.rm=TRUE),sum(effects$fdr[effects$scope=="TCGA EAC"]<.05,na.rm=TRUE)))
write.csv(summary,summary_path,row.names=FALSE)
summary_copy <- tcga_live_path(file.path(SCREF_SUMMARY_DIR,"Auto_burden_summary.csv"))
write.csv(summary,summary_copy,row.names=FALSE)
stopifnot(file.info(summary_copy)$size<100000,nrow(contrasts)==17L*sum(sample_groups$eligible),
 all(cells$group[cells$burden>cells$cutoff]=="High"),all(cells$group[cells$burden<=cells$cutoff]=="Low"),
 nrow(effects)==17L*(length(studies)+2L),!anyDuplicated(effects[c("dataset","scope","mp")]))
writeLines(c("PASS: manual InferCNA score independently reproduced from raw matrix; orig.ident exact; current 17 MPs.",
 "PASS: median split within sample; at least 20 cells in each included group; ties kept together.",
 "PASS: TCGA autosomal segment FGA at absolute log2 ratio >=0.20; no overlapping segments; median cohort split.",
 "PASS: sample-level inference; 1000 bootstrap resamples of samples; full 17-MP BH families.",
 "Caution: within-tumour scATLAS and between-tumour TCGA are distinct contrasts and distinct burden/MP assays.",
 "Caution: relative high/low are median-defined, not validated biological CNA thresholds; single-cell RNA signal can reflect technical variation.",
 "Caution: sample units use orig.ident; repeated samples from a patient, if present, are not automatically collapsed."),
 file.path(tiers[["reports"]],"Auto_burden_validation.txt"))
run <- start_run_summary(script_path,inputs,c(figure_path,summary_path),list(sc_burden="InferCNA mean squared gene signal",boxplot_test="paired raw sample means in scATLAS; unpaired TCGA",tertile_sensitivity="type-7 1/3 and 2/3; middle excluded",tcga_burden="autosomal segment FGA abs log2>=0.20",minimum_cells_per_group=20,bootstrap_replicates=1000,seed=42))
run <- add_cache_status(run,cache$label,cache$path,cache$reused)
write_run_summary(run,file.path(tiers[["logs"]],"Auto_burden_run_summary.rds"))
message("CNV burden analysis complete: ",out_dir)
####################
