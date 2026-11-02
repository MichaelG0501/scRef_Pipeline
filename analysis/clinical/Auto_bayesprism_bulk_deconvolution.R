####################
# Analysis registry
# Status: active
# Script: analysis/clinical/Auto_bayesprism_bulk_deconvolution.R
# Description: One corrected scATLAS reference for OCCAMS and TCGA EAC BayesPrism.
#   Malignancy comes from upstream calls, independently of current state coverage.
# Methodology: analysis/methodology/clinical/malignant_bulk_survival_methodology.md
# Inputs relative to the live project root:
#   ref_outs/{EAC_Ref_merged,meta_full_epi}.rds
#   ref_outs/by_samples/<orig.ident>/<orig.ident>_epi.rds: existing InferCNA
#     cell and cluster classifications for epithelial cells outside the canonical
#     malignant metadata. Exact paths/presence/hashes are recorded in
#     tables/epithelial_source_manifest.csv. Missing calls are excluded, never
#     interpreted as normal and never replaced with another file or classifier.
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz
#   ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv
#   ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz
#   analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx
#     sheet EAC_RNA88_Final: authoritative EAC patient and primary-sample selection.
#   ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#     contributes explicit GDC file/sample/patient IDs and pinned file checksums only.
#   ref_outs/TCGA/esca_gdc_reconstruction/raw/gdc_files/<file_id>/<file_name>
#     canonical live copies of GDC GENCODE-v36 STAR unstranded counts; no downloads.
#   analysis/shared/scRef_config.R
# Clinical survival metadata are read separately by Auto_malignant_bulk_survival.R.
# Outputs: ref_outs/bulk_malignant_deconvolution/
#   README.md: pointer to the current deconvolution result index.
#   reference/reference_bundle.rds; tables/{input_manifest,reference_labels,
#     reference_cell_audit,reference_gene_filter,tcga_count_manifest,
#     epithelial_source_manifest,reference_donor_coverage}.csv
#   tables/{malignant_reference_similarity,state_posterior_cv_summary,state_fraction_resolution}.csv
#   intermediate/{OCCAMS,TCGA}_{counts,fit,extracted}.rds
#   tables/{OCCAMS,TCGA}_{fractions,state_posterior_cv,gene_coverage,inference_gene_filter}.csv
#   figures/Auto_diagnostics/01_deconvolution_quality.pdf
#   reports/{current_results.md,validation.csv}; logs/package_versions.csv
#   updates/new_updates/summaries/bulk_malignant_deconvolution.csv (<100 KB)
# Cache: preparation and fits require matching hashes/settings. Changed inputs
#   stop; SCREF_FORCE_REBUILD=TRUE explicitly rebuilds the selected task.
#   SCREF_REPLOT_ONLY=TRUE is supported only for the qc task using extracted RDS.
# Tasks, in order: prepare; OCCAMS and TCGA (independent); qc.
# Direct allocated-node run: Rscript analysis/clinical/Auto_bayesprism_bulk_deconvolution.R <task>
# PBS: /opt/pbs/bin/qsub -v SCREF_DECONV_TASK=<task> job_scripts/clinical/Auto_bayesprism_bulk_deconvolution.sh
# Environment dmtcp; fit 16 CPUs/96 GB/6 h; prepare 4 CPUs/96 GB/1 h; qc 1 CPU/8 GB/30 min.
# All PBS stdout/stderr stay in the project root. No temporary tool is a dependency.
####################
library(BayesPrism)
library(Matrix)
library(Seurat)
library(data.table)
library(readxl)
library(ggplot2)
library(patchwork)

project_dir <- '/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline'
setwd(project_dir)
source('analysis/shared/scRef_config.R')
if (!nzchar(Sys.getenv('PBS_JOBID'))) stop('Run within a PBS allocation')
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L || !args[1] %in% c('prepare','OCCAMS','TCGA','qc')) stop('Supply prepare, OCCAMS, TCGA, or qc')
task <- args[1]
force <- identical(toupper(Sys.getenv('SCREF_FORCE_REBUILD','FALSE')),'TRUE')
replot <- identical(toupper(Sys.getenv('SCREF_REPLOT_ONLY','FALSE')),'TRUE')
if (force && replot) stop('Force rebuild and replot are mutually exclusive')
if (replot && task!='qc') stop('Replot applies only to qc')
cores <- as.integer(Sys.getenv('SCREF_DECONV_CORES','16'))
if (!is.finite(cores) || cores<1L) stop('Invalid core count')
out <- file.path(project_dir,'ref_outs/bulk_malignant_deconvolution')
for (d in c('reference','intermediate','tables','reports','logs','figures/Auto_diagnostics'))
  dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
inputs <- c(merged='ref_outs/EAC_Ref_merged.rds',epithelial='ref_outs/meta_full_epi.rds',
  states=SCREF_STATE_NOREG_RDS,mp_genes=SCREF_MP_GENES_RDS,grouping=SCREF_MP_GROUPING_CSV,
  occams_counts='ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz',
  occams_qc='ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv',
  gtf='ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz',
  tcga_workbook='analysis/clinical/TCGA_ESCA_EAC_AJCC8_MASTER_corrected_2026-09-18.xlsx',
  tcga_mapping='ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds',
  config='analysis/shared/scRef_config.R')
stopifnot(all(file.exists(inputs)))
settings <- list(version=2L,malignant_calls=c('malignant_level_1','malignant_level_2'),
  normal_epithelial_rule='outside canonical malignant set; cna_non_malignant AND non_malignant_clus',
  min_cells_gene=20L,min_cells_state=20L,outlier_cut=.01,outlier_fraction=.1,
  excluded_gene_classes='mitochondrial; RPL/RPS ribosomal; GENCODE-v19 X/Y',seed=22062026L,
  bayesprism=as.character(packageVersion('BayesPrism')))
states <- c(SCREF_PRIMARY_STATE_ORDER,SCREF_TECHNICAL_STATE_ORDER,'Malignant state unassigned')
bundle_path <- file.path(out,'reference/reference_bundle.rds')
table_out <- function(x,name) write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
fingerprint <- function(paths) data.frame(input=names(paths),path=normalizePath(paths),md5=unname(tools::md5sum(paths)))
####################
epithelial_manifest <- function(paths) {
  present <- file.exists(paths)
  checksum <- rep(NA_character_,length(paths))
  checksum[present] <- unname(tools::md5sum(paths[present]))
  data.frame(path=paths,present=present,md5=checksum,stringsAsFactors=FALSE)
}
####################
count_path <- function(cohort) file.path(out,'intermediate',paste0(cohort,'_counts.rds'))
extract_path <- function(cohort) file.path(out,'intermediate',paste0(cohort,'_extracted.rds'))
collapse_symbols <- function(x,symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  if (!any(keep)) stop('No mapped symbols')
  rowsum(x[keep,,drop=FALSE],group=symbols[keep],reorder=FALSE)
}
check_counts <- function(x) {
  stopifnot(is.matrix(x),!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),
    all(is.finite(x)),all(x>=0),all(x==round(x)),all(colSums(x)>0))
}

if (task=='prepare') {
  manifest <- fingerprint(inputs)
  if (file.exists(bundle_path) && !force) {
    old <- readRDS(bundle_path)
    if (!identical(old$manifest,manifest) || !identical(old$settings,settings)) stop('Stale reference; force rebuild required')
####################
    stopifnot(all(file.exists(old$counts_manifest$path)),
      identical(unname(tools::md5sum(old$counts_manifest$path)),old$counts_manifest$md5))
    stopifnot(identical(old$epithelial_sources,epithelial_manifest(old$epithelial_sources$path)))
####################
    stopifnot(all(file.exists(vapply(c('OCCAMS','TCGA'),count_path,character(1)))))
    message('Verified existing reference and count inputs');quit(save='no')
  }
  message('Preparing canonical raw-count matrices')
  qc <- fread(inputs[['occams_qc']])
  occams_ids <- qc$Subject_ID[qc$Annotation_Compatible_Below_50_Flag=='no']
  stopifnot(length(occams_ids)==281L,!anyDuplicated(occams_ids))
  gtf <- fread(inputs[['gtf']],sep='\t',header=FALSE,skip='chr',quote='',select=c(1,3,9))
  gtf <- gtf[gtf[[2]]=='gene',]
  annotation <- data.frame(ensembl=sub('.*gene_id "([^"]+)".*','\\1',gtf[[3]]),
    symbol=sub('.*gene_name "([^"]+)".*','\\1',gtf[[3]]),chromosome=gtf[[1]])
  annotation$ensembl <- sub('\\..*$','',annotation$ensembl)
  stopifnot(!anyDuplicated(annotation$ensembl),nrow(annotation)>10000L)
  symbol_map <- setNames(annotation$symbol,annotation$ensembl)
  map_symbols <- function(ids) {
    known <- ids %in% annotation$symbol
    result <- unname(symbol_map[sub('\\..*$','',ids)])
    result[known] <- ids[known]
    result
  }
  raw <- fread(inputs[['occams_counts']])
  stopifnot(names(raw)[1]=='Geneid',all(occams_ids %in% names(raw)))
  occams <- collapse_symbols(as.matrix(raw[,..occams_ids]),map_symbols(raw[[1]]))
  rm(raw);gc()
  workbook <- as.data.frame(read_excel(inputs[['tcga_workbook']],sheet='EAC_RNA88_Final'))
  metadata <- readRDS(inputs[['tcga_mapping']])
  stopifnot(nrow(workbook)==88L,!anyDuplicated(workbook$Patient_ID),!anyDuplicated(metadata$sample_barcode))
  metadata <- metadata[match(workbook$RNA_primary_sample_ID,metadata$sample_barcode),,drop=FALSE]
  stopifnot(!anyNA(metadata$file_id),identical(metadata$case_barcode,workbook$Patient_ID),
    all(metadata$HistologyGroup=='EAC'),all(metadata$sample_type_code=='01'))
  paths <- file.path(project_dir,'ref_outs/TCGA/esca_gdc_reconstruction/raw/gdc_files',metadata$file_id,metadata$file_name)
  stopifnot(all(file.exists(paths)),identical(unname(tools::md5sum(paths)),metadata$md5sum))
  count_manifest <- data.frame(subject_id=metadata$case_barcode,sample_id=metadata$sample_barcode,
    file_id=metadata$file_id,path=paths,md5=metadata$md5sum)
  tcga_parts <- lapply(paths,function(path) {
    x <- fread(path,skip=1L)
    stopifnot(all(c('gene_id','gene_name','unstranded') %in% names(x)))
    x <- x[grepl('^ENSG',gene_id),]
    list(genes=x$gene_id,symbols=x$gene_name,counts=x$unstranded)
  })
  stopifnot(all(vapply(tcga_parts,function(x) identical(x$genes,tcga_parts[[1]]$genes) &&
    identical(x$symbols,tcga_parts[[1]]$symbols),logical(1))))
  tcga <- do.call(cbind,lapply(tcga_parts,`[[`,'counts'))
  colnames(tcga) <- metadata$sample_barcode
  tcga <- collapse_symbols(tcga,tcga_parts[[1]]$symbols)
  rm(tcga_parts);gc()
  check_counts(occams);check_counts(tcga)

  message('Aggregating scATLAS reference with independent malignancy calls')
  atlas <- readRDS(inputs[['merged']])
  epithelial <- readRDS(inputs[['epithelial']])
  state_vector <- readRDS(inputs[['states']])
  meta <- atlas[[]]
  counts <- GetAssayData(atlas,assay='RNA',layer='counts')
  stopifnot(inherits(counts,'dgCMatrix'),identical(colnames(counts),rownames(meta)),
    all(c('orig.ident','celltype_update') %in% names(meta)),
    all(c('orig.ident','malignancy') %in% names(epithelial)),!is.null(names(state_vector)),
    !anyDuplicated(names(state_vector)),all(state_vector %in% states))
  types <- gsub('\\.','_',as.character(meta$celltype_update))
  stopifnot(!anyNA(types),all(nzchar(types)),!anyNA(meta$orig.ident))
  malignant_ids <- rownames(epithelial)[epithelial$malignancy %in% settings$malignant_calls]
  stopifnot(all(names(state_vector) %in% malignant_ids),all(malignant_ids %in% rownames(meta)))
  is_malignant <- rownames(meta) %in% malignant_ids
  stopifnot(all(meta$celltype_update[is_malignant]=='epithelial'))
####################
  # Normal epithelial reference cells require both independent CNA evidence
  # and a nonmalignant cluster call. Ambiguous or unavailable calls are ineligible.
  extra_ids <- rownames(meta)[meta$celltype_update=='epithelial' & !is_malignant]
  extra_samples <- unique(as.character(meta$orig.ident[match(extra_ids,rownames(meta))]))
  source_paths <- file.path(project_dir,'ref_outs/by_samples',extra_samples,paste0(extra_samples,'_epi.rds'))
  normal_sources <- epithelial_manifest(source_paths)
  calls <- data.frame(cell_id=extra_ids,classification=NA_character_,malignant_clus=NA_character_,
    source_status='source unavailable',stringsAsFactors=FALSE)
  for (k in seq_along(extra_samples)) {
    ids <- extra_ids[meta$orig.ident[match(extra_ids,rownames(meta))]==extra_samples[k]]
    j <- match(ids,calls$cell_id)
    if (!normal_sources$present[k]) next
    source_meta <- readRDS(source_paths[k])[[]]
    stopifnot(all(c('orig.ident','classification','malignant_clus') %in% names(source_meta)),
      all(source_meta$orig.ident==extra_samples[k]))
    m <- match(ids,rownames(source_meta))
    calls$classification[j] <- as.character(source_meta$classification[m])
    calls$malignant_clus[j] <- as.character(source_meta$malignant_clus[m])
    calls$source_status[j] <- ifelse(is.na(m),'cell excluded by upstream QC','classified')
    rm(source_meta);gc()
  }
  normal_ids <- calls$cell_id[calls$classification %in% 'cna_non_malignant' &
    calls$malignant_clus %in% 'non_malignant_clus']
  if (length(normal_ids)<settings$min_cells_state) stop('Insufficient independently supported normal epithelial reference cells')
  uncertain_ids <- setdiff(extra_ids,normal_ids)
  types[rownames(meta) %in% normal_ids] <- 'Non_malignant_epithelial'
  types[rownames(meta) %in% uncertain_ids] <- 'Epithelial malignancy unresolved'
####################
  types[is_malignant] <- 'Malignant'
  labels <- types
  labels[is_malignant] <- 'Malignant state unassigned'
  idx <- match(names(state_vector),rownames(meta))
  labels[idx] <- as.character(state_vector)
  keep_cells <- meta$celltype_update!='unresolved_inconsistent' & !rownames(meta) %in% uncertain_ids
  stopifnot(all(is_malignant<=keep_cells))
  cell_audit <- data.frame(cell_id=rownames(meta),orig.ident=meta$orig.ident,
    original_cell_type=meta$celltype_update,malignancy=epithelial$malignancy[match(rownames(meta),rownames(epithelial))],
    reference_type=types,reference_state=labels,included=keep_cells)
####################
  cell_audit$cna_classification <- calls$classification[match(cell_audit$cell_id,calls$cell_id)]
  cell_audit$cna_cluster <- calls$malignant_clus[match(cell_audit$cell_id,calls$cell_id)]
  cell_audit$additional_epithelial_source <- calls$source_status[match(cell_audit$cell_id,calls$cell_id)]
  cell_audit$exclusion_reason <- ifelse(keep_cells,'',ifelse(meta$celltype_update=='unresolved_inconsistent',
    'unresolved cell type','epithelial malignancy not independently resolved'))
  donor_coverage <- aggregate(list(n_cells=rep(1L,sum(keep_cells))),
    list(orig.ident=meta$orig.ident[keep_cells],cell_state=labels[keep_cells]),sum)
####################
  counts <- counts[,keep_cells,drop=FALSE]
  labels <- labels[keep_cells];types <- types[keep_cells]
  levels <- unique(labels)
  design <- sparseMatrix(i=seq_along(labels),j=match(labels,levels),x=1,
    dims=c(length(labels),length(levels)),dimnames=list(colnames(counts),levels))
  expressed_cells <- Matrix::rowSums(counts>0)
  reference <- as.matrix(counts %*% design)
  symbols <- map_symbols(rownames(reference))
  gene_audit <- data.frame(gene=rownames(reference),symbol=symbols,expressed_cells=expressed_cells)
  valid <- !is.na(symbols) & expressed_cells>=settings$min_cells_gene
  reference <- collapse_symbols(reference[valid,,drop=FALSE],symbols[valid])
  reference <- t(reference)
  reference_types <- vapply(levels,function(s) {
    x <- unique(types[labels==s]);stopifnot(length(x)==1L);x
  },character(1))
  label_table <- data.frame(cell_state=levels,cell_type=unname(reference_types),n_cells=as.integer(table(factor(labels,levels=levels))))
  stopifnot(all(label_table$n_cells>=settings$min_cells_state),all(states %in% levels))
  common <- Reduce(intersect,list(colnames(reference),rownames(occams),rownames(tcga)))
  exclude <- grepl('^MT-|^RPL[0-9]|^RPS[0-9]',common) |
    common %in% annotation$symbol[annotation$chromosome %in% c('chrX','chrY','X','Y')]
  common <- common[!exclude]
  common <- common[colSums(reference[,common,drop=FALSE])>0 &
    rowSums(occams[common,,drop=FALSE])>0 & rowSums(tcga[common,,drop=FALSE])>0]
  stopifnot(length(common)>1000L)
  gene_audit$included <- gene_audit$symbol %in% common & valid
  reference <- reference[,common,drop=FALSE]
  occams <- occams[common,,drop=FALSE];tcga <- tcga[common,,drop=FALSE]
  counts_manifest <- fingerprint(setNames(paths,metadata$sample_barcode))
  bundle <- list(reference=reference,labels=label_table,genes=common,settings=settings,
    manifest=manifest,counts_manifest=counts_manifest,tcga_mapping=count_manifest,epithelial_sources=normal_sources)
  saveRDS(occams,count_path('OCCAMS'));saveRDS(tcga,count_path('TCGA'))
  saveRDS(bundle,bundle_path)
  table_out(manifest,'input_manifest');table_out(label_table,'reference_labels')
  table_out(cell_audit,'reference_cell_audit');table_out(gene_audit,'reference_gene_filter')
  table_out(count_manifest,'tcga_count_manifest')
####################
  table_out(normal_sources,'epithelial_source_manifest');table_out(donor_coverage,'reference_donor_coverage')
####################
  message('Reference ready: ',length(common),' genes; ',sum(is_malignant),' malignant cells')
}

if (task %in% c('OCCAMS','TCGA')) {
  stopifnot(file.exists(bundle_path),file.exists(count_path(task)))
  bundle <- readRDS(bundle_path)
  stopifnot(identical(bundle$settings,settings))
####################
  if (!identical(bundle$manifest,fingerprint(inputs))) stop('Upstream reference inputs changed; rerun prepare')
  stopifnot(identical(bundle$epithelial_sources,epithelial_manifest(bundle$epithelial_sources$path)))
  stopifnot(all(file.exists(bundle$counts_manifest$path)),
    identical(unname(tools::md5sum(bundle$counts_manifest$path)),bundle$counts_manifest$md5))
  # Apply BayesPrism's default bulk-outlier rule to both cohorts and remove the
  # union before either fit. This keeps the inferred gene universe identical.
  all_counts <- setNames(lapply(c('OCCAMS','TCGA'),function(c) readRDS(count_path(c))),c('OCCAMS','TCGA'))
  stopifnot(identical(rownames(all_counts$OCCAMS),rownames(all_counts$TCGA)))
  outliers <- lapply(all_counts,function(x) {
    check_counts(x)
    rowMeans(sweep(x,2,colSums(x),'/')>settings$outlier_cut)>settings$outlier_fraction
  })
  gene_filter <- data.frame(gene=rownames(all_counts$OCCAMS),
    outlier_occams=outliers$OCCAMS,outlier_tcga=outliers$TCGA,
    retained=!(outliers$OCCAMS|outliers$TCGA))
  genes <- gene_filter$gene[gene_filter$retained]
  stopifnot(length(genes)>1000L)
  counts <- all_counts[[task]][genes,,drop=FALSE]
  rm(all_counts);gc()
  table_out(gene_filter,paste0(task,'_inference_gene_filter'))
  contract <- list(settings=settings,inference_filter='union of cohort-specific default outlier flags',
    inputs=fingerprint(c(reference=bundle_path,OCCAMS_counts=count_path('OCCAMS'),TCGA_counts=count_path('TCGA'))))
####################
  fit_path <- file.path(out,'intermediate',paste0(task,'_fit.rds'))
  if (file.exists(fit_path) && !force) {
    saved <- readRDS(fit_path)
    if (!identical(saved$contract,contract)) stop('Stale fit; force rebuild required')
    fit <- saved$fit
  } else {
    set.seed(settings$seed)
    prism <- new.prism(reference=bundle$reference[,genes,drop=FALSE],mixture=t(counts),input.type='GEP',
      cell.type.labels=bundle$labels$cell_type,cell.state.labels=bundle$labels$cell_state,
      key='Malignant',outlier.cut=1,outlier.fraction=1)
####################
    stopifnot(identical(colnames(prism@mixture),genes))
####################
    fit <- run.prism(prism=prism,n.cores=cores,gibbs.control=list(seed=settings$seed))
    saveRDS(list(contract=contract,fit=fit),fit_path)
  }
  first <- get.fraction(fit,which.theta='first',state.or.type='state')
  initial <- get.fraction(fit,which.theta='first',state.or.type='type')
  final <- get.fraction(fit,which.theta='final',state.or.type='type')
  expression <- get.exp(fit,state.or.type='type',cell.name='Malignant')
####################
  stopifnot(identical(colnames(expression),genes))
####################
  stopifnot(all(states %in% colnames(first)),identical(rownames(first),rownames(final)),
    identical(rownames(first),rownames(expression)),all(is.finite(expression)),all(expression>=0),
    all(rowSums(expression)>0),all(is.finite(first)),all(first>=0),
    max(abs(rowSums(first)-1))<1e-6,max(abs(rowSums(final)-1))<1e-6)
  malignant <- rowSums(first[,states,drop=FALSE])
  stopifnot(all(malignant>0),max(abs(malignant-initial[,'Malignant']))<1e-6)
  within <- first[,states,drop=FALSE]/malignant
  fraction <- data.frame(sample_id=rownames(first),malignant_initial=malignant,
    malignant_final=final[,'Malignant'],check.names=FALSE)
  for (s in states) {
    fraction[[paste0('bulk__',s)]] <- first[,s]
    fraction[[paste0('within__',s)]] <- within[,s]
  }
  cv <- fit@posterior.initial.cellState@theta.cv
  stopifnot(identical(dimnames(cv),dimnames(first)),all(is.finite(cv)),all(cv>=0))
  extracted <- list(cohort=task,fractions=fraction,within=within,bulk=first,
    state_cv=cv,type_initial=initial,type_final=final,tumour_expression=expression,
    contract=contract,reference_labels=bundle$labels)
  saveRDS(extracted,extract_path(task))
  table_out(fraction,paste0(task,'_fractions'))
  table_out(data.frame(sample_id=rownames(cv),cv,check.names=FALSE),paste0(task,'_state_posterior_cv'))
  grouping <- read.csv(inputs[['grouping']]);sets <- readRDS(inputs[['mp_genes']])[grouping$mp]
  coverage <- data.frame(mp=names(sets),defined=lengths(sets),
    inferred=vapply(sets,function(g) sum(unique(g)%in%colnames(expression)),integer(1)))
  stopifnot(nrow(coverage)==17L,all(coverage$inferred>=5L))
  table_out(coverage,paste0(task,'_gene_coverage'))
  message(task,' fitted and extracted: ',nrow(expression),' subjects; ',ncol(expression),' genes')
}

if (task=='qc') {
  objects <- setNames(lapply(c('OCCAMS','TCGA'),function(c) readRDS(extract_path(c))),c('OCCAMS','TCGA'))
####################
  stopifnot(identical(colnames(objects$OCCAMS$tumour_expression),colnames(objects$TCGA$tumour_expression)))
  reference <- readRDS(bundle_path)$reference[states,colnames(objects$OCCAMS$tumour_expression),drop=FALSE]
  similarity <- cor(t(reference),method='spearman')
  similarity_table <- as.data.frame(as.table(similarity),stringsAsFactors=FALSE)
  names(similarity_table) <- c('state_one','state_two','spearman_rho')
  table_out(similarity_table,'malignant_reference_similarity')
  cv_summary <- rbindlist(lapply(objects,function(x) rbindlist(lapply(states,function(s)
    data.frame(cohort=x$cohort,state=s,median_cv=median(x$state_cv[,s]),
      percentile95_cv=unname(quantile(x$state_cv[,s],.95)),maximum_cv=max(x$state_cv[,s]))))))
  table_out(cv_summary,'state_posterior_cv_summary')
####################
  resolution <- rbindlist(lapply(objects,function(x) rbindlist(lapply(states,function(s)
    data.frame(cohort=x$cohort,state=s,n=nrow(x$within),
      median_within_fraction=median(x$within[,s]),
      samples_below_point_one_percent=sum(x$within[,s]<.001),
      median_initial_state_cv=median(x$state_cv[,s]),
      samples_cv_above_point_five=sum(x$state_cv[,s]>.5))))))
  table_out(resolution,'state_fraction_resolution')
####################
####################
  fractions <- rbindlist(lapply(objects,function(x) transform(x$fractions,cohort=x$cohort)))
  composition <- rbindlist(lapply(objects,function(x) {
    d <- data.table(sample_id=rownames(x$within),x$within,check.names=FALSE)
    d[,sample_rank:=rank(x$fractions$malignant_final,ties.method='first')]
    z <- melt(d,id.vars=c('sample_id','sample_rank'),variable.name='state',value.name='fraction')
    z[,cohort:=x$cohort];z
  }))
  composition$state <- factor(composition$state,levels=states)
  colours <- c(SCREF_STATE_COLOURS,'Malignant state unassigned'='#626262')
  theme_set(theme_classic(base_size=11)+theme(legend.position='bottom',strip.background=element_blank(),
    strip.text=element_text(face='bold'),plot.title=element_text(face='bold'),plot.margin=margin(8,8,8,8)))
  p1 <- ggplot(fractions,aes(malignant_initial,malignant_final))+
    geom_abline(slope=1,intercept=0,colour='grey65',linetype=2)+geom_point(size=1.2,alpha=.6,colour='#0072B2')+
    facet_wrap(~cohort)+coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    scale_x_continuous(breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+
    theme(panel.spacing.x=grid::unit(1.2,'lines'))+
    labs(x='Initial malignant RNA fraction',y='Updated malignant RNA fraction',title='a  Malignant compartment')
  p2 <- ggplot(composition,aes(sample_rank,fraction,fill=state))+geom_col(width=1)+
    facet_wrap(~cohort,ncol=1,scales='free_x',strip.position='right')+scale_fill_manual(values=colours,drop=FALSE)+
    scale_y_continuous(labels=scales::label_percent(),expand=c(0,0))+
    guides(fill=guide_legend(ncol=2,byrow=TRUE))+labs(x='Samples ordered by updated malignant fraction',
      y='Within-malignant RNA fraction',fill=NULL,title='c  State composition')
####################
  cv_summary$state <- factor(cv_summary$state,levels=rev(states))
  p3 <- ggplot(cv_summary,aes(median_cv,state,colour=cohort,shape=cohort))+
    geom_errorbar(aes(xmin=median_cv,xmax=percentile95_cv),orientation='y',width=0,linewidth=.55,position=position_dodge(.5))+
    geom_point(size=2,position=position_dodge(.5))+
    scale_colour_manual(values=c(OCCAMS='#0072B2',TCGA='#D55E00'))+
    scale_y_discrete(labels=function(x) vapply(x,function(s) paste(strwrap(s,22),collapse='\n'),character(1)))+
    labs(x='Initial state-fraction posterior CV',y=NULL,colour=NULL,shape=NULL,
      title='b  Conditional state uncertainty',subtitle='Point: median; line endpoint: 95th percentile across samples')+
    theme(axis.text.y=element_text(size=10),plot.subtitle=element_text(size=9))
  quality_plot <- (p1|p3)/p2+plot_layout(heights=c(1,1.15))
  ggsave(file.path(out,'figures/Auto_diagnostics/01_deconvolution_quality.pdf'),quality_plot,
    width=12,height=11,device=cairo_pdf)
####################
  validation <- rbindlist(lapply(objects,function(x) data.frame(cohort=x$cohort,n=nrow(x$within),
    genes=ncol(x$tumour_expression),malignant_median=median(x$fractions$malignant_final),
    state_unassigned_median=median(x$within[,'Malignant state unassigned']),
    max_closure_error=max(abs(rowSums(x$within)-1)),max_posterior_cv=max(x$state_cv))))
  stopifnot(identical(validation$n,c(281L,88L)),all(validation$max_closure_error<1e-6))
  write.csv(validation,file.path(out,'reports/validation.csv'),row.names=FALSE)
  dir.create(SCREF_SUMMARY_DIR,recursive=TRUE,showWarnings=FALSE)
  write.csv(validation,file.path(SCREF_SUMMARY_DIR,'bulk_malignant_deconvolution.csv'),row.names=FALSE)
  writeLines(c('# Current malignant deconvolution',
    'Both cohorts use the same corrected scATLAS reference and preprocessing. Malignancy is independent of state availability.',
    'All upstream malignant cells remain malignant. Cells without a current state are explicitly state-unassigned.',
    'Fractions represent RNA contributions over retained genes, not validated cell counts. State fractions use the initial posterior; total malignant fractions also use the updated posterior.',
    'Normal epithelial cells require both cna_non_malignant and non_malignant_clus calls, outside the canonical malignant metadata. Cells without concordant calls are excluded. These computational labels are not independent histological validation.',
    'Related state profiles, donor imbalance and reference uncertainty remain biological limitations; the similarity table quantifies reference-profile overlap.',
    'Posterior CV describes conditional sampling uncertainty, not reference uncertainty or biological validation.',
####################
    'Several initial state fractions have high relative posterior uncertainty and many within-malignant fractions are below 0.1%. The resolution table reports these observations without treating near-zero inferred values as accurately measured biological variation. State-fraction survival associations require this qualification; numerical posterior closure is not biological validation.',
    '[State-fraction resolution](../tables/state_fraction_resolution.csv) | [Posterior CV](../tables/state_posterior_cv_summary.csv)',
####################
    '[Quality figure](../figures/Auto_diagnostics/01_deconvolution_quality.pdf) | [Validation](validation.csv) | [Reference labels](../tables/reference_labels.csv)'),
    file.path(out,'reports/current_results.md'))
####################
  writeLines(c('# Current OCCAMS and TCGA deconvolution','',
    'Start with the [current results and quality checks](reports/current_results.md).',
    'Corrected reference, fitted posteriors and extracted quantities remain in `reference/` and `intermediate/`.',
    'Methods: [malignant-compartment workflow](../../analysis/methodology/clinical/malignant_bulk_survival_methodology.md).'),
    file.path(out,'README.md'))
####################
}
write.csv(data.frame(package=c('R','BayesPrism','Seurat','Matrix'),
  version=c(as.character(getRversion()),vapply(c('BayesPrism','Seurat','Matrix'),function(p) as.character(packageVersion(p)),character(1)))),
  file.path(out,'logs/package_versions.csv'),row.names=FALSE)
message('Task complete: ',task)
