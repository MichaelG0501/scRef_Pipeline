####################
# Status: active
# Script: analysis/GSE_Krish/gse_bayesprism_deconvolution.R
# Description: BayesPrism malignant expression and eight-state RNA fractions for 110 GSE207526 tumour biopsies.
# Methodology: analysis/methodology/GSE_Krish/gse_krish_methodology.md.
# Inputs: ref_outs/GSE_Krish/reference/GSE207526/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt.gz;
#   ref_outs/GSE_Krish/tables/{clinical_metadata_clean,geo_sample_mapping}.csv;
#   ref_outs/bulk_malignant_deconvolution/reference/reference_bundle.rds (corrected scATLAS reference);
#   analysis/shared/scRef_config.R; current centred MP genes/grouping specified there.
# Outputs under ref_outs/GSE_Krish/: intermediate/{GSE_counts,GSE_fit,GSE_extracted}.rds;
#   tables/{deconvolution_input_manifest,inference_gene_filter,gene_coverage,GSE_fractions,
#   GSE_state_posterior_cv,excluded_expression_rows}.csv; reports/deconvolution_validation.txt;
#   logs/deconvolution_packages.csv. All persisted matrices remain in live.
# Cache: input/settings fingerprint must match; SCREF_FORCE_REBUILD=TRUE rebuilds.
#   No plotting; SCREF_REPLOT_ONLY=TRUE rejected. Seed 22062026, SCREF_DECONV_CORES default 16.
# Run: Rscript analysis/GSE_Krish/gse_bayesprism_deconvolution.R
# PBS: job_scripts/GSE_Krish/gse_bayesprism_deconvolution.sh; dmtcp; 16 CPUs, 96 GB, 6 hours.
####################
library(BayesPrism)
library(data.table)
source('analysis/shared/scRef_config.R')
if (!nzchar(Sys.getenv('PBS_JOBID'))) stop('Run within a PBS allocation')
if (toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE'))=='TRUE') stop('No plotting in this step')
force <- toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'))=='TRUE'
cores <- as.integer(Sys.getenv('SCREF_DECONV_CORES','16'))
stopifnot(is.finite(cores),cores>=1L)
out <- file.path(SCREF_REF_OUTS_DIR,'GSE_Krish')
for (tier in c('intermediate','tables','reports','logs')) dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
inputs <- c(matrix=file.path(out,'reference/GSE207526/GSE207526_110.EAC.and.10.Normal.for.GSEA.txt.gz'),
  clinical=file.path(out,'tables/clinical_metadata_clean.csv'),mapping=file.path(out,'tables/geo_sample_mapping.csv'),
  reference='ref_outs/bulk_malignant_deconvolution/reference/reference_bundle.rds',
  genes=SCREF_MP_GENES_RDS,grouping=SCREF_MP_GROUPING_CSV,config='analysis/shared/scRef_config.R')
stopifnot(all(file.exists(inputs)))
manifest <- data.frame(input=names(inputs),path=normalizePath(inputs),md5=unname(tools::md5sum(inputs)))
settings <- list(version=1L,seed=22062026L,outlier_cut=.01,outlier_fraction=.1,
  reference='existing corrected scATLAS pseudobulk; no new reference classification',
  BayesPrism=as.character(packageVersion('BayesPrism')))
contract <- list(manifest=manifest,settings=settings)
fit_path <- file.path(out,'intermediate/GSE_fit.rds')
meta <- read.csv(inputs[['clinical']],stringsAsFactors=FALSE)
stopifnot(nrow(meta)==110L,!anyDuplicated(meta$patient_id))
bundle <- readRDS(inputs[['reference']])
stopifnot(identical(bundle$settings$normal_epithelial_rule,
 'outside canonical malignant set; cna_non_malignant AND non_malignant_clus'),
 bundle$settings$version==2L)
# Validate the biological reference sources inside the compute allocation.
stopifnot(all(file.exists(bundle$manifest$path)),
 identical(unname(tools::md5sum(bundle$manifest$path)),bundle$manifest$md5))
stopifnot(identical(file.exists(bundle$epithelial_sources$path),bundle$epithelial_sources$present))
present <- bundle$epithelial_sources$present
stopifnot(identical(unname(tools::md5sum(bundle$epithelial_sources$path[present])),bundle$epithelial_sources$md5[present]))
x <- as.data.frame(fread(cmd=paste('gzip -dc',shQuote(inputs[['matrix']])),header=FALSE,skip=1,quote='"'))
header <- scan(gzfile(inputs[['matrix']]),what=character(),sep='\t',nlines=1,quiet=TRUE,quote='"')
stopifnot(ncol(x)==length(header)+1L,length(header)==120L,all(meta$patient_id %in% header))
names(x) <- c('gene',header)
keep <- !is.na(x$gene) & nzchar(x$gene)
write.csv(data.frame(source_row=which(!keep),gene=x$gene[!keep],reason='missing gene symbol'),
 file.path(out,'tables/excluded_expression_rows.csv'),row.names=FALSE)
stopifnot(!anyDuplicated(x$gene[keep]))
counts <- as.matrix(x[keep,meta$patient_id,drop=FALSE]);storage.mode(counts)<-'double';rownames(counts)<-x$gene[keep]
rm(x);gc()
stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==round(counts)),all(colSums(counts)>0))
common <- intersect(colnames(bundle$reference),rownames(counts))
stopifnot(length(common)>1000L)
counts <- counts[common,,drop=FALSE]
positive <- rowSums(counts)>0
outlier <- rowMeans(sweep(counts,2,colSums(counts),'/')>settings$outlier_cut)>settings$outlier_fraction
filter <- data.frame(gene=common,positive=positive,outlier=outlier,retained=positive & !outlier)
genes <- common[filter$retained];counts<-counts[genes,,drop=FALSE]
stopifnot(length(genes)>1000L)
write.csv(filter,file.path(out,'tables/inference_gene_filter.csv'),row.names=FALSE)
saveRDS(counts,file.path(out,'intermediate/GSE_counts.rds'))
if(file.exists(fit_path) && !force) {
 saved<-readRDS(fit_path)
 if(!identical(saved$contract,contract))stop('Stale fit; explicit SCREF_FORCE_REBUILD=TRUE required')
 fit<-saved$fit
} else {
 set.seed(settings$seed)
 prism<-new.prism(reference=bundle$reference[,genes,drop=FALSE],mixture=t(counts),input.type='GEP',
  cell.type.labels=bundle$labels$cell_type,cell.state.labels=bundle$labels$cell_state,
  key='Malignant',outlier.cut=1,outlier.fraction=1)
 stopifnot(identical(colnames(prism@mixture),genes))
 fit<-run.prism(prism=prism,n.cores=cores,gibbs.control=list(seed=settings$seed))
 saveRDS(list(contract=contract,fit=fit),fit_path)
}
first<-get.fraction(fit,which.theta='first',state.or.type='state')
initial<-get.fraction(fit,which.theta='first',state.or.type='type')
final<-get.fraction(fit,which.theta='final',state.or.type='type')
expression<-get.exp(fit,state.or.type='type',cell.name='Malignant')
states<-c(SCREF_PRIMARY_STATE_ORDER,SCREF_TECHNICAL_STATE_ORDER,'Malignant state unassigned')
stopifnot(all(states %in% colnames(first)),identical(rownames(first),rownames(expression)),
 setequal(rownames(expression),meta$patient_id),identical(colnames(expression),genes),
 all(is.finite(expression)),all(expression>=0),all(rowSums(expression)>0),
 max(abs(rowSums(first)-1))<1e-6,max(abs(rowSums(final)-1))<1e-6)
malignant<-rowSums(first[,states,drop=FALSE]);stopifnot(all(malignant>0))
within<-first[,states,drop=FALSE]/malignant
stopifnot(max(abs(malignant-initial[,'Malignant']))<1e-6,max(abs(rowSums(within)-1))<1e-6)
cv<-fit@posterior.initial.cellState@theta.cv
stopifnot(identical(dimnames(cv),dimnames(first)),all(is.finite(cv)),all(cv>=0))
fraction<-data.frame(sample_id=rownames(first),malignant_initial=malignant,malignant_final=final[,'Malignant'])
for(s in states){fraction[[paste0('bulk__',s)]]<-first[,s];fraction[[paste0('within__',s)]]<-within[,s]}
extracted<-list(cohort='GSE207526',fractions=fraction,within=within,bulk=first,state_cv=cv,
 type_initial=initial,type_final=final,tumour_expression=expression,contract=contract,reference_labels=bundle$labels)
saveRDS(extracted,file.path(out,'intermediate/GSE_extracted.rds'))
write.csv(fraction,file.path(out,'tables/GSE_fractions.csv'),row.names=FALSE)
write.csv(data.frame(sample_id=rownames(cv),cv,check.names=FALSE),file.path(out,'tables/GSE_state_posterior_cv.csv'),row.names=FALSE)
grouping<-read.csv(inputs[['grouping']]);sets<-readRDS(inputs[['genes']])[grouping$mp]
coverage<-data.frame(mp=names(sets),defined=lengths(sets),inferred=vapply(sets,function(g)sum(unique(g)%in%genes),integer(1)))
stopifnot(nrow(coverage)==17L,all(coverage$inferred>=5L))
write.csv(coverage,file.path(out,'tables/gene_coverage.csv'),row.names=FALSE)
write.csv(manifest,file.path(out,'tables/deconvolution_input_manifest.csv'),row.names=FALSE)
write.csv(data.frame(package=c('R','BayesPrism','data.table'),version=c(as.character(getRversion()),
 as.character(packageVersion('BayesPrism')),as.character(packageVersion('data.table')))),file.path(out,'logs/deconvolution_packages.csv'),row.names=FALSE)
writeLines(c(paste('Validated subjects:',nrow(expression)),paste('Inference genes:',length(genes)),
 'Original biopsies were collected before treatment. Treatment strata are user-defined from Mandard availability.',
 'Fractions estimate RNA contributions, not directly measured cell proportions.'),file.path(out,'reports/deconvolution_validation.txt'))
message('DECONVOLUTION_VALIDATED: ',nrow(expression),' patients, ',length(genes),' genes')
