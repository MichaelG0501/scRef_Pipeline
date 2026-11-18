####################
# Status: active
# Script: analysis/GEO/geo_microarray_scores.R
# Description: Fixed-parent malignant GSVA and state RNA fractions for all 52 GSE19417 EAC specimens.
# Methodology: analysis/methodology/GEO/geo_microarray_methodology.md
# Inputs: ref_outs/GEO/intermediate/{GSE_counts,GSE_extracted}.rds; tables/clinical_metadata.csv;
#   canonical current MP genes/grouping, ranked five-state markers; analysis/shared/scRef_config.R.
# Outputs: ref_outs/GEO/intermediate/scoring_cache.rds; tables/{scoring_input_manifest,
#   scoring_gene_coverage,scored_features,model_data,feature_definitions,gene_sets}.csv;
#   reports/scoring_validation.txt. No plots.
# Cache: input/settings/package/code contract; SCREF_FORCE_REBUILD=TRUE rebuilds; no replot mode.
# Run/PBS: Rscript analysis/GEO/geo_microarray_scores.R; job_scripts/GEO/geo_microarray_scores.sh;
#   dmtcp; 4 CPUs, 24 GB, 1 hour.
####################
library(GSVA)
source('analysis/shared/scRef_config.R')
stopifnot(nzchar(Sys.getenv('PBS_JOBID')),toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))!='TRUE')
out <- 'ref_outs/GEO'
inputs <- c(mixture=file.path(out,'intermediate/GSE_counts.rds'),extracted=file.path(out,'intermediate/GSE_extracted.rds'),
 clinical=file.path(out,'tables/clinical_metadata.csv'),genes=SCREF_MP_GENES_RDS,grouping=SCREF_MP_GROUPING_CSV,
 markers='ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv',config='analysis/shared/scRef_config.R')
stopifnot(all(file.exists(inputs)))
manifest <- data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
make_scores <- function() {
 meta <- read.csv(inputs[['clinical']],stringsAsFactors=FALSE)
 grouping <- read.csv(inputs[['grouping']],stringsAsFactors=FALSE)
 genes <- readRDS(inputs[['genes']])[grouping$mp];states <- SCREF_PRIMARY_STATE_ORDER
 stopifnot(nrow(meta)==52L,!anyDuplicated(meta$patient_id),length(genes)==17L)
 unions <- lapply(SCREF_STATE_GROUPS,function(mps)unique(unlist(genes[mps],use.names=FALSE)))
 names(unions) <- paste0(states,' (MP union)')
 markers <- read.csv(inputs[['markers']],stringsAsFactors=FALSE)
 dge <- lapply(states,function(s)head(markers$gene[markers$state==s],20L));names(dge) <- paste0(states,' (DGE)')
 stopifnot(all(lengths(dge)==20L))
 sets <- c(genes,unions,dge)
 defs <- data.frame(feature=c(names(sets),paste0(states,' (fraction)')),
  family=c(rep('MP',17),rep('State MP union',5),rep('State DGE',5),rep('State fraction',5)),
  state=c(grouping$state,rep(states,3)),label=c(paste(grouping$mp,grouping$description,sep=' | '),rep(states,3)))
 mixture <- readRDS(inputs[['mixture']]);e <- readRDS(inputs[['extracted']])
 stopifnot(setequal(meta$patient_id,rownames(e$tumour_expression)))
 x <- t(e$tumour_expression[meta$patient_id,,drop=FALSE])
 stopifnot(all(is.finite(x)),all(x>=0),all(colSums(x)>0))
 x <- log2(sweep(x,2,colSums(x),'/')*1e6+1)
 common <- rownames(x)[apply(x,1,function(z)diff(range(z))>0)]
 used <- lapply(sets,function(g)intersect(unique(g),common));stopifnot(length(common)>1000L,all(lengths(used)>=5L))
 scores <- t(GSVA::gsva(x[common,,drop=FALSE],used,method='gsva',kcdf='Gaussian',min.sz=5L,
  max.sz=Inf,mx.diff=TRUE,abs.ranking=FALSE,verbose=FALSE,parallel.sz=4L))
 f <- e$within[meta$patient_id,states,drop=FALSE];colnames(f) <- paste0(states,' (fraction)')
 scores <- cbind(scores,f);stopifnot(ncol(scores)==32L,all(is.finite(scores)))
 data <- do.call(rbind,lapply(colnames(scores),function(feature)data.frame(patient_id=meta$patient_id,
  source='malignant',feature=feature,value=scores[,feature],OS_time=meta$OS_time,OS_event=meta$OS_event)))
 list(data=data,definitions=defs,coverage=data.frame(feature=names(used),defined=lengths(sets),retained=lengths(used)),
  gene_sets=do.call(rbind,lapply(names(sets),function(f)data.frame(feature=f,gene=sets[[f]],retained=sets[[f]]%in%common))))
}
contract <- list(manifest=manifest,version=1L,GSVA=as.character(packageVersion('GSVA')),code=deparse(body(make_scores)))
p <- file.path(out,'intermediate/scoring_cache.rds')
if(file.exists(p)&&toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))!='TRUE') {
 cache <- readRDS(p);stopifnot(identical(cache$contract,contract))
} else {cache <- list(contract=contract,result=make_scores());saveRDS(cache,p)}
write.csv(manifest,file.path(out,'tables/scoring_input_manifest.csv'),row.names=FALSE)
for(name in c('data','definitions','coverage','gene_sets')) {
 file <- c(data='model_data',definitions='feature_definitions',coverage='scoring_gene_coverage',gene_sets='gene_sets')[[name]]
 write.csv(cache$result[[name]],file.path(out,'tables',paste0(file,'.csv')),row.names=FALSE,na='')
}
write.csv(cache$result$data,file.path(out,'tables/scored_features.csv'),row.names=FALSE,na='')
writeLines(c('52 EAC specimens; fixed-parent scoring; no treatment subgroup restriction.',
 '27 expression features plus five within-malignant RNA fractions.',
 'Microarray deconvolution and cross-platform survival synthesis remain exploratory.'),file.path(out,'reports/scoring_validation.txt'))
message('GEO_SCORES_COMPLETE: 52 patients, 32 features')
