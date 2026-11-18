####################
# Analysis registry
# Status: active
# Script: analysis/developmental/prepare_developmental_nmf_references.R
# Description: Prepare annotated, non-integrated Seurat RNA objects for developmental NMF/reference work.
# Methodology: analysis/methodology/developmental/developmental_nmf_reference_preparation.md
# Inputs: SCREF_DEVELOPMENTAL_SOURCE (default /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental),
#   original RDS/H5AD/MTX and source annotations mapped by source_dirs below.
# Outputs: <bundle>/<study>.rds for all seven references; cell/gene metadata under Original_data/<study>/.
#   Documentation and reports are written only to project ref_outs/developmental_reference/reports/.
#   ref_outs/developmental_reference/reports/nmf_reference_preparation_summary.csv
# Cache: existing validated preparation-version objects reused unless SCREF_FORCE_REBUILD=TRUE.
#   SCREF_REPLOT_ONLY: not applicable; no plots or NMF are generated.
# Run: /opt/pbs/bin/qsub job_scripts/developmental/prepare_developmental_nmf_references.sh
# Environment: dmtcp; PBS 4 CPUs, 180 GB, 06:00:00. No batch correction or gene filtering.
####################
library(Seurat)
library(Matrix)
library(data.table)
library(hdf5r)
library(SingleCellExperiment)
library(SummarizedExperiment)

project <- '/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline'
bundle <- Sys.getenv('SCREF_DEVELOPMENTAL_BUNDLE','/rds/general/user/el1224/home/OAC_data/Developmental')
####################
source_root <- Sys.getenv('SCREF_DEVELOPMENTAL_SOURCE','/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental')
source_dirs <- c(Yu_2021_developing_endoderm='raw_data/yu_2021',Zhao_early_embryogenesis='raw_data/early_embryogenesis',Xu_organogenesis_GSE157329='raw_data/organogenesis',Cao_fetal_stomach='raw_data/normal_development',Ding_adult_oesophagus='raw_data/adult_oesophagus/SCP1242/expression/63f53992d91a88956d36dc4f',Tsubosaka_adult_stomach='.',Nowicki_Osuch_Zhuang_Barretts_cardia_glands='raw_data/barretts')
pending_studies <- character() # Seven final references; limitations embedded in each object.
read_h5_metadata <- function(group) {
 keys<-group$ls(recursive=FALSE)$name
 index<-as.character(group[['_index']][])
 values<-lapply(setdiff(keys,'_index'),function(key) {
  node<-group[[key]]
  if(inherits(node,'H5Group')) {
   stopifnot(all(c('codes','categories') %in% node$ls(recursive=FALSE)$name))
   codes<-as.integer(node[['codes']][]);cats<-node[['categories']][]
   out<-cats[pmax(codes+1L,1L)];out[codes<0L]<-NA;out
  } else node[]
 })
 names(values)<-setdiff(keys,'_index')
 frame<-as.data.frame(values,stringsAsFactors=FALSE,check.names=FALSE)
 rownames(frame)<-index
 list(frame=frame,index=index)
}
####################
force <- identical(Sys.getenv('SCREF_FORCE_REBUILD','FALSE'),'TRUE')
version <- 'annotated_rna_reference_1'
folders <- c('Yu_2021_developing_endoderm','Zhao_early_embryogenesis','Xu_organogenesis_GSE157329','Cao_fetal_stomach','Ding_adult_oesophagus','Tsubosaka_adult_stomach','Nowicki_Osuch_Zhuang_Barretts_cardia_glands')
write_gz <- function(x,path) { con<-gzfile(path,'wt');on.exit(close(con));write.csv(x,con,row.names=FALSE) }
read_table <- function(path) as.data.frame(fread(path,data.table=FALSE,check.names=FALSE))
clean_missing <- function(x) { x<-as.character(x);x[is.na(x)|x==''|x %in% c('NA','nan','None')]<-NA_character_;x }
normalize_genes <- function(m,ids,symbols) {
 if(length(ids)!=nrow(m)||length(symbols)!=nrow(m)) stop('Gene mapping dimension mismatch')
 symbols<-clean_missing(symbols);unmapped<-is.na(symbols);symbols[unmapped]<-ids[unmapped]
 symbols<-gsub('_','-',symbols,fixed=TRUE)
 if(anyNA(symbols)||any(symbols=='')) stop('Missing gene identifiers')
 unique_symbols<-unique(symbols)
 if(anyDuplicated(symbols)) {
  if(current_processed) stop('Duplicated genes in processed log expression; summation is not justified')
  collapse<-sparseMatrix(i=match(symbols,unique_symbols),j=seq_along(symbols),x=1,dims=c(length(unique_symbols),length(symbols)))
  m<-collapse %*% m
 }
 rownames(m)<-unique_symbols
 list(matrix=as(m,'dgCMatrix'),mapping=data.frame(source_gene_id=ids,source_symbol=symbols,output_gene=symbols,unmapped_symbol=unmapped))
}
validate <- function(obj) {
 stopifnot(inherits(obj,'Seurat'),identical(obj@misc$reference_preparation$version,version))
 md<-obj@meta.data;a<-obj@assays[['RNA']]
 stopifnot(nrow(md)==ncol(a@data),identical(rownames(md),colnames(a@data)),!anyDuplicated(rownames(md)),nrow(a@data)>0)
 stopifnot(all(is.finite(a@data@x)),all(a@data@x>=0))
 stopifnot(all(c('source_dataset','sample_id','donor_id','orig.ident','cell_type_source','age_label','age_weeks_source','age_days_source','sample_id_status') %in% colnames(md)))
 invisible(TRUE)
}
summaries<-list()
for(folder in folders) {
 cat('\nPreparing',folder,format(Sys.time()),'\n');flush.console()
 ####################
 d<-file.path(source_root,source_dirs[[folder]])
 artifact_dir<-file.path(bundle,'Original_data',folder)
 dir.create(artifact_dir,recursive=TRUE,showWarnings=FALSE)
 outfile<-if(folder %in% pending_studies)file.path(artifact_dir,'nmf_reference.rds')else file.path(bundle,paste0(folder,'.rds'))
 ####################
 if(file.exists(outfile)&&!force) {
  obj<-readRDS(outfile);validate(obj)
 } else {
  current_processed<-FALSE;donor<-sample<-age<-weeks<-days<-tissue<-annotation<-NULL
  sample_status<-'explicit source sample identifier';donor_status<-'explicit source donor identifier'
  if(folder=='Yu_2021_developing_endoderm') {
   h<-H5File$new(file.path(d,'developing_human_atlas.h5ad'),'r');x<-h[['X']]
   obs<-read_h5_metadata(h[['obs']]);md<-obs$frame;md$cell_id<-obs$index
   var<-read_h5_metadata(h[['var']]);gm<-var$frame;gm$gene_id<-var$index
   shape<-as.integer(x$attr_open('shape')$read())
   stopifnot(identical(as.character(x$attr_open('encoding-type')$read()),'csr_matrix'))
   m<-sparseMatrix(i=as.integer(x[['indices']][])+1L,p=as.integer(x[['indptr']][]),x=as.numeric(x[['data']][]),dims=rev(shape))
   h$close_all();colnames(m)<-md$cell_id;ids<-gm$gene_id;symbols<-gm$feature_name
   sample<-md$Sample;donor<-md$donor_id;age<-md$Age_week;weeks<-suppressWarnings(as.numeric(age));tissue<-md$annotated_tissue;annotation<-md$author_cell_type
   input_status<-'Integer counts in H5AD X; LogNormalize 10000 applied'
  } else if(folder=='Zhao_early_embryogenesis') {
   source<-readRDS(file.path(d,'psd.R3.6.em.seurat.ob.rds'));md<-source@meta.data
   md$cell_id<-rownames(md);m<-source@assays[['RNA']]@counts;ids<-symbols<-rownames(m)
   sample<-rep('pooled_atlas',nrow(md));donor<-rep(NA_character_,nrow(md));age<-md$devTime;tissue<-md$EML;annotation<-md$rename_EML
   sample_status<-'pooled atlas analysis unit; not an independent donor; source sample identities unresolved'
   donor_status<-'unresolved';input_status<-'Source integer RNA counts; LogNormalize 10000 applied';rm(source)
  } else if(folder=='Xu_organogenesis_GSE157329') {
   md<-read_table(file.path(d,'GSE157329_cell_annotate.txt.gz'));gm<-read_table(file.path(d,'GSE157329_gene_annotate.txt.gz'))
   m<-as(readMM(gzfile(file.path(d,'GSE157329_raw_counts.mtx.gz'))),'dgCMatrix')
   stopifnot(nrow(m)==nrow(gm),ncol(m)==nrow(md));colnames(m)<-md$cell_id;ids<-gm$gene_id;symbols<-gm$gene_short_name
   sample<-md$sample;donor<-md$embryo;age<-md$stage;tissue<-md$dissection_part;annotation<-md$final_annotation
   input_status<-'Source integer count matrix; LogNormalize 10000 applied'
  } else if(folder=='Cao_fetal_stomach') {
   md<-readRDS(file.path(d,'df_cell.RDS'));md$cell_id<-as.character(md$sample);gm<-readRDS(file.path(d,'df_gene.RDS'))
   m<-readRDS(file.path(d,'Stomach_gene_count.RDS'));md<-md[match(colnames(m),md$cell_id),,drop=FALSE]
   ids<-rownames(m);symbols<-gm$gene_short_name[match(ids,gm$gene_id)]
   donor<-md$Fetus_id;sample<-md$Fetus_id;sample_status<-'donor-level stomach unit: Fetus_id; source sample field is a cell ID'
   age<-md$Development_day;days<-suppressWarnings(as.numeric(age));tissue<-md$Organ;annotation<-md$Main_cluster_name
   input_status<-'Source integer count matrix; LogNormalize 10000 applied'
  } else if(folder=='Ding_adult_oesophagus') {
   md<-read_table(file.path(source_root,'raw_data/adult_oesophagus/SCP1242/metadata/EoE_meta.txt'));md<-md[md$NAME!='TYPE',,drop=FALSE];md$cell_id<-md$NAME
   genes<-fread(file.path(d,'EoE_gene.tsv'),header=FALSE,data.table=FALSE)[[1]]
   cells<-fread(file.path(d,'EoE_cell.tsv'),header=FALSE,data.table=FALSE)[[1]]
   m<-as(readMM(file.path(d,'EoE.mtx')),'dgCMatrix');stopifnot(nrow(m)==length(genes),ncol(m)==length(cells))
   colnames(m)<-cells;md<-md[match(cells,md$cell_id),,drop=FALSE];ids<-symbols<-genes
   sample<-md$biosample_id;donor<-md$donor_id;age<-rep('Adult (study designation)',nrow(md));tissue<-md$organ;annotation<-md$cell_type_anno
   input_status<-'Source integer count matrix; LogNormalize 10000 applied'
  } else if(folder=='Tsubosaka_adult_stomach') {
   source<-readRDS(file.path(d,'Stomach.rds'));md<-source@meta.data;md$cell_id<-rownames(md)
   m<-source@assays[['RNA']]@counts;ids<-symbols<-rownames(m)
   sample<-md$orig.ident;donor<-md$patient;age<-rep('Adult (study designation)',nrow(md));tissue<-rep('Stomach; exact region requires source curation',nrow(md));annotation<-md$subcluster.v2
   input_status<-'Fractional processed/decontaminated count-like RNA counts; LogNormalize 10000 applied; original provenance retained'
   rm(source)
  } else if(folder=='Nowicki_Osuch_Zhuang_Barretts_cardia_glands') {
   source<-readRDS(file.path(d,'alldatahighquality.rds'));md<-as.data.frame(colData(source));md$cell_id<-rownames(md)
   md$cell_id_basis<-'Source colData row identifiers; original barcode names not supplied in this object'
   m<-as(assay(source,'logcounts'),'dgCMatrix');colnames(m)<-md$cell_id;ids<-rownames(source);symbols<-as.character(rowData(source)$Symbol)
   sample<-md$Sample;donor<-md$Patient;age<-rep('Adult (study designation)',nrow(md));tissue<-md$Tissue;annotation<-md$cell_type_secondary
   current_processed<-TRUE;input_status<-'Source logcounts retained unchanged; exact original transform unverified; no raw counts supplied'
   rm(source)
  }
  stopifnot(!anyNA(md$cell_id),!anyDuplicated(md$cell_id),ncol(m)==nrow(md),identical(as.character(md$cell_id),colnames(m)))
  stopifnot(all(is.finite(m@x)),all(m@x>=0))
  mapped<-normalize_genes(m,ids,symbols);m<-mapped$matrix
  write_gz(mapped$mapping,file.path(artifact_dir,'nmf_reference_gene_mapping.csv.gz'));rm(mapped)
  rownames(md)<-md$cell_id
  if('orig.ident' %in% names(md))md$source_orig_ident<-md$orig.ident
  md$source_dataset<-folder;md$sample_id<-clean_missing(sample);md$donor_id<-clean_missing(donor)
  md$sample_id_status<-sample_status;md$donor_id_status<-donor_status
  md$orig.ident<-ifelse(is.na(md$sample_id),NA_character_,paste(folder,md$sample_id,sep='::'))
  md$cell_type_source<-clean_missing(annotation);md$tissue_source<-clean_missing(tissue);md$age_label<-clean_missing(age)
  md$age_weeks_source<-if(is.null(weeks))rep(NA_real_,nrow(md))else weeks
  md$age_days_source<-if(is.null(days))rep(NA_real_,nrow(md))else days
  md$age_basis<-'Source age label retained; no cross-study conversion or imputation'
  options(Seurat.object.assay.version='v3')
  if(current_processed) {
   a<-CreateAssayObject(data=m);obj<-CreateSeuratObject(counts=a,meta.data=md,project=folder)
  } else {
   obj<-CreateSeuratObject(counts=m,meta.data=md,project=folder,min.cells=0,min.features=0)
   obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
  }
  obj@misc$reference_preparation<-list(version=version,source_folder=d,input_status=input_status,assay='RNA',nmf_layer='data',sample_id_status=sample_status,donor_id_status=donor_status,
   readiness=if(current_processed)'prepared using supplied nonnegative logcounts; exact source transform undocumented; no renormalization'else if(folder=='Zhao_early_embryogenesis')'prepared as one pooled atlas unit; no donor-level replication; curate lineage before NMF'else 'prepared; select lineage and eligible samples before NMF',
   integrated=FALSE,scaled=FALSE,all_source_cells_retained=TRUE)
  validate(obj)
  saveRDS(obj,outfile,compress=FALSE)
  write_gz(obj@meta.data,file.path(artifact_dir,'nmf_reference_metadata.csv.gz'))
  rm(list=intersect(c('m','md','gm','ids','symbols'),ls()));invisible(gc())
 }
 meta<-obj@meta.data;prep<-obj@misc$reference_preparation
 summaries[[folder]]<-data.frame(dataset=folder,rds=outfile,cells=nrow(meta),genes=nrow(obj@assays[['RNA']]@data),sample_units=length(unique(na.omit(meta$orig.ident))),donors=length(unique(na.omit(meta$donor_id))),sample_id_status=prep$sample_id_status,input_status=prep$input_status,readiness=prep$readiness,bytes=file.info(outfile)$size)
 cat('SAVED',outfile,'cells',nrow(meta),'\n');flush.console();rm(obj);invisible(gc())
}
summary<-do.call(rbind,summaries)
write.csv(summary,file.path(project,'ref_outs/developmental_reference/reports/nmf_reference_preparation_summary.csv'),row.names=FALSE)
cat('Seven annotated reference RDS objects prepared successfully\n')
