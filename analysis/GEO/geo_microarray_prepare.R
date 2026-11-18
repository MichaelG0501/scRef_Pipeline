####################
# Status: active
# Script: analysis/GEO/geo_microarray_prepare.R
# Description: Recover tumour-channel intensities for all 52 histologically confirmed EACs in GSE19417.
# Methodology: analysis/methodology/GEO/geo_microarray_methodology.md
# Inputs: ref_outs/GEO/tables/download_manifest.csv; downloaded GPL4372 annotation and all sample files;
#   ref_outs/geo_survival/Auto_GSE19417_meta.csv (public clinical endpoint/histology mapping).
# Outputs: ref_outs/GEO/intermediate/GEO_linear_expression.rds;
#   tables/{clinical_metadata,probe_mapping,array_qc,preparation_manifest}.csv; reports/preparation.txt.
# Cache: deterministic preparation; explicit current outputs updated. No fitting or plot-only mode.
# Run/PBS: Rscript analysis/GEO/geo_microarray_prepare.R; job_scripts/GEO/geo_microarray_prepare.sh;
#   dmtcp, 1 CPU, 8 GB, 30 minutes.
####################
library(limma)
library(data.table)
stopifnot(nzchar(Sys.getenv('PBS_JOBID')))
out <- 'ref_outs/GEO'
for(tier in c('intermediate','tables','reports'))dir.create(file.path(out,tier),recursive=TRUE,showWarnings=FALSE)
m <- read.csv('ref_outs/geo_survival/Auto_GSE19417_meta.csv',stringsAsFactors=FALSE)
manifest <- read.csv(file.path(out,'tables/download_manifest.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(m)==76L,nrow(manifest)==78L,!anyDuplicated(m$sample_geo_accession))
meta <- m[m$HistologyGroup=='EAC',];stopifnot(nrow(meta)==52L,!anyDuplicated(meta$patient_id))
stopifnot(all(meta$source_name_ch1=='Stratagene Universal Human Reference RNA'),all(meta$label_ch2=='Cy5'))
paths <- manifest$path[match(meta$sample_geo_accession,manifest$sample)]
stopifnot(!anyNA(paths),all(file.exists(paths)))
tables <- lapply(paths,function(p)as.data.frame(fread(cmd=paste('gzip -dc',shQuote(p)))))
probe_ids <- tables[[1]]$ID_REF
stopifnot(!anyDuplicated(probe_ids),all(vapply(tables,function(t)identical(t$ID_REF,probe_ids),logical(1))))
# INTENSITY2 is the recorded tumour channel, not the exponentiated ratio.
x <- do.call(cbind,lapply(tables,function(t)t$INTENSITY2));storage.mode(x)<-'double'
complete_probe <- apply(x,1,function(z)all(is.finite(z)&z>=0))
stopifnot(sum(complete_probe)>1000L)
x <- x[complete_probe,,drop=FALSE]
stopifnot(all(colSums(x)>0))
colnames(x) <- meta$sample_geo_accession;rownames(x) <- as.character(probe_ids[complete_probe])
####################
ratio_diagnostics <- t(vapply(tables,function(t) {
  ok <- is.finite(t$INTENSITY1)&is.finite(t$INTENSITY2)&is.finite(t$RATIO)&t$INTENSITY1>0&t$INTENSITY2>0
  raw_ratio <- log10(t$INTENSITY2[ok]/t$INTENSITY1[ok])
  error <- abs(raw_ratio-t$RATIO[ok])
  c(median_difference=median(error),p99_difference=unname(quantile(error,.99)),
    maximum_difference=max(error),channel_ratio_correlation=cor(raw_ratio,t$RATIO[ok]))
},numeric(4)))
qc <- data.frame(excluded_incomplete_probes=sum(!complete_probe),sample_id=colnames(x),
 zero_fraction=colMeans(x==0),total_intensity=colSums(x),ratio_diagnostics)
# GEO channel provenance defines tumour INTENSITY2. Modelled RATIO need not equal
# the simple ratio of reported intensities; retain those differences as diagnostics.
stopifnot(all(is.finite(qc$channel_ratio_correlation)),all(qc$channel_ratio_correlation>0))
####################
# Quantile normalisation corrects sample intensity distributions, without logging ratios or inventing counts.
x <- normalizeBetweenArrays(x,method='quantile')
p <- manifest$path[manifest$sample=='platform'];lines <- readLines(gzfile(p),warn=FALSE)
a <- match('!platform_table_begin',lines);b <- match('!platform_table_end',lines)
stopifnot(is.finite(a),is.finite(b),b>a)
annotation <- read.delim(text=paste(lines[(a+1L):(b-1L)],collapse='\n'),check.names=FALSE,quote='"',comment.char='')
symbol <- trimws(annotation[['Gene symbol']][match(rownames(x),as.character(annotation$ID))])
keep <- !is.na(symbol)&nzchar(symbol)&!grepl('///|;|,',symbol)
map <- data.frame(probe_id=rownames(x),gene_symbol=symbol,retained=keep)
groups <- split(which(keep),symbol[keep]);linear <- t(vapply(groups,function(j)apply(x[j,,drop=FALSE],2,median),numeric(ncol(x))))
colnames(linear) <- meta$sample_geo_accession
stopifnot(nrow(linear)>1000L,all(is.finite(linear)),all(linear>=0))
clinical <- data.frame(patient_id=meta$sample_geo_accession,original_patient_id=meta$patient_id,
 OS_time=meta$OS_time,OS_event=meta$OS_event,histology=meta$HistologyGroup)
saveRDS(linear,file.path(out,'intermediate/GEO_linear_expression.rds'))
write.csv(clinical,file.path(out,'tables/clinical_metadata.csv'),row.names=FALSE,na='')
write.csv(map,file.path(out,'tables/probe_mapping.csv'),row.names=FALSE,na='')
write.csv(qc,file.path(out,'tables/array_qc.csv'),row.names=FALSE)
inputs <- c(paths,p,'ref_outs/geo_survival/Auto_GSE19417_meta.csv',file.path(out,'tables/download_manifest.csv'))
write.csv(data.frame(path=normalizePath(inputs),md5=unname(tools::md5sum(inputs))),file.path(out,'tables/preparation_manifest.csv'),row.names=FALSE)
writeLines(c('52 EAC specimens; gastric and squamous cancers excluded by explicit histology.',
 'Tumour INTENSITY2 identified by GEO channel provenance; positive concordance with recorded RATIO verified, full discrepancies tabulated.',
 'Linear quantile normalisation across all 52 EAC samples; median collapse of unambiguous probe symbols.',
 'Microarray intensities are not sequencing counts. BayesPrism use remains exploratory.'),file.path(out,'reports/preparation.txt'))
