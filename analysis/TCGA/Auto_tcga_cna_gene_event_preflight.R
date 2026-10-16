####################
# Analysis registry:
#   Status: active
#   Script: analysis/TCGA/Auto_tcga_cna_gene_event_preflight.R
#   Description: verify manual InferCNA subclone lineage and catalogue TCGA gene CNA profiles.
#   Methodology: not required (deterministic provenance and API catalogue checks; no association fitting).
#   Inputs: InferCNA.R; analysis/cnv/cnv_malignant_subclone_mp_heatmap.R;
#     analysis/cnv/cna_subclone_expression_correlation.R;
#     ref_outs/Auto_malignant_subclone_mp/threshold_1/Auto_malignant_subclone_cells.csv;
#     ref_outs/Auto_cna_subclone_expression/threshold_1/rds/Auto_cna_subclone_expression_results.rds;
#     ref_outs/by_samples/<sample>/<sample>_outs.rds and <sample>_epi.rds;
#     ref_outs/meta_full_epi.rds;
#     /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt;
#     https://www.cbioportal.org/api/studies/esca_tcga_gdc/molecular-profiles
#   Outputs: ref_outs/TCGA/cna_gene_event_validation/reports/Auto_manual_infercna_provenance.csv;
#     reports/Auto_input_inventory.txt; intermediate/Auto_cbio_molecular_profiles.rds;
#     tables/Auto_cbio_molecular_profiles.csv; logs/Auto_preflight_run_summary.rds/.txt.
#   Cache/replot: network catalogue is queried afresh; existing files are retained.
#   Run: /opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_cna_gene_event_preflight.sh; conda: dmtcp.
####################
library(data.table)
library(dplyr)
library(httr2)
library(Seurat)
source("analysis/shared/scRef_config.R")
source("analysis/shared/scRef_helpers.R")
source("analysis/TCGA/Auto_tcga_workflow_helpers.R")
tcga_execution_guard()
out_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "cna_gene_event_validation")
tiers <- ensure_output_tiers(out_dir)
cell_path <- file.path(SCREF_REF_OUTS_DIR, "Auto_malignant_subclone_mp", "threshold_1", "Auto_malignant_subclone_cells.csv")
result_path <- file.path(SCREF_REF_OUTS_DIR, "Auto_cna_subclone_expression", "threshold_1", "rds", "Auto_cna_subclone_expression_results.rds")
cells <- fread(cell_path, data.table = FALSE)
results <- readRDS(result_path)
source_paths <- c("InferCNA.R", "analysis/cnv/cnv_malignant_subclone_mp_heatmap.R", "analysis/cnv/cna_subclone_expression_correlation.R")
source_text <- lapply(source_paths, readLines)
stopifnot(any(grepl("outs <- infercna\\(", source_text[[1]])),
          any(grepl("subclone <- infer_cna_subclones", source_text[[2]], fixed = TRUE)),
          !any(grepl("numbat", unlist(source_text), ignore.case = TRUE)))
sample_id <- sort(unique(cells$sample))[1]
outs_path <- file.path(SCREF_REF_OUTS_DIR, "by_samples", sample_id, paste0(sample_id, "_outs.rds"))
outs <- readRDS(outs_path)
# Reproduce the cached arm means from the exact raw InferCNA matrix for one complete sample.
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
go <- read.table(gene_order_path, header = FALSE, col.names = c("gene", "chr", "start", "end"), stringsAsFactors = FALSE)
go <- go[!duplicated(go$gene), ]
centromere <- c(chr1=121700000,chr2=91800000,chr3=87900000,chr4=50600000,chr5=48400000,chr6=61000000,chr7=59900000,chr8=45600000,chr9=49000000,chr10=40200000,chr11=53400000,chr12=35500000,chr13=17700000,chr14=17200000,chr15=19000000,chr16=36800000,chr17=25100000,chr18=18500000,chr19=26200000,chr20=28100000,chr21=12000000,chr22=15000000,chrX=61000000)
go <- go[go$gene %in% rownames(outs) & go$chr %in% names(centromere), ]
go$arm <- paste0(go$chr, ifelse(go$start <= centromere[go$chr], "p", "q"))
sample_cells <- cells[cells$sample == sample_id, ]
check_rows <- list()
for (clone in unique(sample_cells$subclone)) {
  ids <- sample_cells$cell[sample_cells$subclone == clone]
  stopifnot(all(ids %in% colnames(outs)))
  for (arm in unique(go$arm)) {
    genes <- go$gene[go$arm == arm]
    value <- mean(outs[genes, ids, drop = FALSE], na.rm = TRUE)
    cached <- results$arm_long
    cached <- cached[cached$sample == sample_id & cached$subclone == clone & cached$arm == arm, ]
    stopifnot(nrow(cached) == 1L)
    check_rows[[length(check_rows) + 1L]] <- data.frame(sample=sample_id, subclone=clone, arm=arm,
      source="manual InferCNA; not Numbat", observed=value, cached=cached$arm_mean,
      abs_difference=abs(value-cached$arm_mean), input=outs_path)
  }
}
checks <- bind_rows(check_rows)
stopifnot(all(checks$abs_difference < 1e-10))
write.csv(checks, file.path(tiers[["reports"]], "Auto_manual_infercna_provenance.csv"), row.names = FALSE)
# Audit every assigned cell against canonical orig.ident without barcode parsing.
meta <- readRDS(SCREF_META_FULL_EPI_RDS)
stopifnot(all(cells$cell %in% rownames(meta)),
          identical(as.character(meta[cells$cell, "orig.ident"]), as.character(cells$sample)))
paths <- file.path(SCREF_REF_OUTS_DIR,"by_samples",unique(cells$sample),paste0(unique(cells$sample),"_epi.rds"))
tcga_require_files(paths)
epi <- readRDS(paths[1])
writeLines(c(paste("manual cells:", nrow(cells)),paste("samples:",length(unique(cells$sample))),
             paste("cached subclones:",nrow(results$features)),paste("gene source:",gene_order_path),
             paste("expression files GB:",round(sum(file.info(paths)$size)/1024^3,3)),
             paste("largest expression RDS GB:",round(max(file.info(paths)$size)/1024^3,3)),
             capture.output(print(Assays(epi))),capture.output(print(Layers(epi[["RNA"]]))),
             capture.output(str(results,max.level=1))),file.path(tiers[["reports"]],"Auto_input_inventory.txt"))
profiles_url <- "https://www.cbioportal.org/api/studies/esca_tcga_gdc/molecular-profiles"
response <- request(profiles_url) |> req_timeout(120) |> req_retry(max_tries=3) |> req_perform()
profiles <- resp_body_json(response, simplifyVector=TRUE)
saveRDS(profiles,file.path(tiers[["intermediate"]],"Auto_cbio_molecular_profiles.rds"))
write.csv(profiles,file.path(tiers[["tables"]],"Auto_cbio_molecular_profiles.csv"),row.names=FALSE)
run <- start_run_summary("analysis/TCGA/Auto_tcga_cna_gene_event_preflight.R", c(source_paths,cell_path,result_path,outs_path,profiles_url),
                        file.path(tiers[["reports"]],"Auto_manual_infercna_provenance.csv"), list(provenance_pass=TRUE,max_arm_error=max(checks$abs_difference)))
write_run_summary(run,file.path(tiers[["logs"]],"Auto_preflight_run_summary.rds"))
message("Preflight PASS: manual InferCNA arm means reproduce exactly; orig.ident verified; molecular profile catalogue saved.")
