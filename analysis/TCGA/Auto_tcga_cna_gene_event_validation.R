####################
# Analysis registry:
#   Status: active
#   Script: analysis/TCGA/Auto_tcga_cna_gene_event_validation.R
#   Description: manual-InferCNA arm/MP and curated gene-CNA/MP comparisons with TCGA.
#   Methodology: analysis/methodology/TCGA/tcga_cna_gene_event_validation_methodology.md
#   Inputs (all persistent paths inside the live project):
#     analysis/TCGA/OAC_CNV.xlsx; analysis/TCGA/41588_2018_331_MOESM3_ESM.xlsx;
#     ref_outs/TCGA/cna_recurrent_event_validation/tables/Auto_oac_cnv_curated_reparsed.csv;
#     ref_outs/TCGA/cna_recurrent_event_validation/intermediate/Auto_tcga_cna_plot_source_data.rds;
#     ref_outs/Auto_cna_subclone_expression/threshold_1/rds/Auto_cna_subclone_expression_results.rds;
#     ref_outs/Auto_malignant_subclone_mp/threshold_1/Auto_malignant_subclone_cells.csv;
#     ref_outs/by_samples/<orig.ident>/<orig.ident>_epi.rds (RNA counts);
#     ref_outs/meta_full_epi.rds;
#     ref_outs/TCGA/cna_recurrent_event_validation/tables/Auto_tcga_arm_cna_long.csv;
#     ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds;
#     ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv;
#     ref_outs/TCGA/cna_gene_event_validation/reports/Auto_manual_infercna_provenance.csv;
#     ref_outs/TCGA/cna_gene_event_validation/intermediate/Auto_cbio_molecular_profiles.rds;
#     API: https://www.cbioportal.org/api (esca_tcga_gdc_cna, gene-level GISTIC).
#   Outputs: ref_outs/TCGA/cna_gene_event_validation/
#     intermediate/Auto_tcga_gene_cna.rds; Auto_sc_manual_subclone_gene_expression.rds;
#       Auto_targeted_cna_results.rds; Auto_targeted_cna_figures.rds;
#       source_data/Auto_cbio_genes.rds; Auto_cbio_cna_sample_ids.rds; Auto_cbio_gene_cna_response.rds;
#     tables/Auto_selected_oac_genes.csv; Auto_tcga_gene_cna_long.csv;
#       Auto_gene_measurement_coverage.csv; Auto_arm_mp_tests.csv; Auto_gene_mp_tests.csv;
#       Auto_arm_pair_ranking.csv; Auto_gene_pair_ranking.csv; Auto_top_arm_pairs.csv;
#       Auto_top_gene_pairs.csv; Auto_sc_sample_level_sensitivity.csv; Auto_targeted_cna_summary.csv;
#       Auto_gene_rectangle_annotations.csv; Auto_gene_pair_boxplot_source.csv;
#       Auto_gene_expression_tertile_thresholds.csv;
#     figures/Auto_oac_gene_mp_spearman_rectangles.pdf;
#       Auto_top_arm_mp_paired_boxplots.pdf; Auto_top_arm_mp_paired_scatterplots.pdf;
#       Auto_top_oac_gene_mp_paired_plots.pdf;
#     reports/Auto_analysis_validation.txt; logs/Auto_targeted_cna_run_summary.rds/.txt;
#       logs/Auto_targeted_cna_replot_run_summary.rds/.txt;
#     updates/new_updates/summaries/Auto_targeted_cna_summary.csv.
#   Cache/replot: SCREF_FORCE_REBUILD=TRUE rebuilds derived caches; API snapshots remain pinned.
#     SCREF_REPLOT_ONLY=TRUE renders the saved figures without loading expression or recomputing tests.
#     SCREF_REFRESH_TCGA_CNA_API=TRUE explicitly refreshes the source snapshots.
#   Run: /opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_cna_gene_event_validation.sh; conda: dmtcp.
####################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(Matrix)
library(httr2)
library(jsonlite)
library(ggpubr)
library(gridExtra)
library(grid)
source("analysis/shared/scRef_config.R")
source("analysis/shared/scRef_helpers.R")
source("analysis/TCGA/Auto_tcga_workflow_helpers.R")
tcga_execution_guard()
set.seed(42)
script_path <- file.path(SCREF_ANALYSIS_DIR,"TCGA","Auto_tcga_cna_gene_event_validation.R")
out_dir <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cna_gene_event_validation")
tiers <- ensure_output_tiers(out_dir)
source_dir <- file.path(tiers[["intermediate"]],"source_data")
dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)
figure_path <- file.path(tiers[["intermediate"]],"Auto_targeted_cna_figures.rds")
summary_path <- file.path(tiers[["tables"]],"Auto_targeted_cna_summary.csv")
if (tcga_flag("SCREF_REPLOT_ONLY")) {
  tcga_replot(figure_path,summary_path,file.path(tiers[["logs"]],"Auto_targeted_cna_replot_run_summary.rds"))
  quit(save="no")
}
parent_dir <- file.path(SCREF_REF_OUTS_DIR,"TCGA","cna_recurrent_event_validation")
sc_dir <- file.path(SCREF_REF_OUTS_DIR,"Auto_cna_subclone_expression","threshold_1")
cell_path <- file.path(SCREF_REF_OUTS_DIR,"Auto_malignant_subclone_mp","threshold_1","Auto_malignant_subclone_cells.csv")
sc_path <- file.path(sc_dir,"rds","Auto_cna_subclone_expression_results.rds")
parent_source_path <- file.path(parent_dir,"intermediate","Auto_tcga_cna_plot_source_data.rds")
provenance_path <- file.path(tiers[["reports"]],"Auto_manual_infercna_provenance.csv")
profile_path <- file.path(tiers[["intermediate"]],"Auto_cbio_molecular_profiles.rds")
oac_path <- file.path(parent_dir,"tables","Auto_oac_cnv_curated_reparsed.csv")
tcga_require_files(c(cell_path,sc_path,parent_source_path,provenance_path,profile_path,oac_path))
provenance <- read.csv(provenance_path)
if (!nrow(provenance) || any(provenance$abs_difference >= 1e-10) || any(provenance$source != "manual InferCNA; not Numbat")) stop("Manual InferCNA provenance gate failed.")
mp_sets <- tcga_validate_gene_sets(readRDS(SCREF_MP_GENES_RDS))
mp_order <- names(mp_sets)
mp_labels <- setNames(label_mps(mp_order),mp_order)
raw_oac <- read.csv(oac_path,check.names=FALSE)
aliases <- c(HER2="ERBB2",P16="CDKN2A",P14ARF="CDKN2A",PARK2="PRKN",RASSF1A="RASSF1")
selected_genes <- raw_oac |>
  mutate(gene=ifelse(gene %in% names(aliases),unname(aliases[gene]),gene)) |>
  arrange(direction,rank,excel_row) |> distinct(direction,gene,.keep_all=TRUE) |>
  group_by(direction) |> slice_head(n=10) |> ungroup() |>
  mutate(direction=factor(direction,levels=c("gain","loss"))) |> arrange(direction,rank,excel_row) |>
  mutate(event_label=paste(ifelse(direction=="gain","Gain","Loss"),gene),event_order=row_number())
stopifnot(nrow(selected_genes)==20L,all(table(selected_genes$direction)==10L))
write.csv(selected_genes,file.path(tiers[["tables"]],"Auto_selected_oac_genes.csv"),row.names=FALSE)
profiles <- readRDS(profile_path)
profile_id <- "esca_tcga_gdc_cna"
profile <- profiles[profiles$molecularProfileId==profile_id,]
if (nrow(profile)!=1L || profile$datatype!="DISCRETE" || profile$molecularAlterationType!="COPY_NUMBER_ALTERATION") stop("Required gene-level GISTIC profile is unavailable; no silent method substitution.")
api_json <- function(path,body=NULL) {
  req <- request(paste0("https://www.cbioportal.org/api",path)) |> req_timeout(120) |> req_retry(max_tries=3)
  if (!is.null(body)) req <- req_body_json(req,body,auto_unbox=TRUE)
  resp_body_json(req_perform(req),simplifyVector=TRUE)
}
source_snapshot <- function(filename,fun) {
  path <- file.path(source_dir,filename)
  if (file.exists(path) && !tcga_flag("SCREF_REFRESH_TCGA_CNA_API")) return(readRDS(path))
  value <- fun()
  saveRDS(value,path)
  value
}
gene_map <- source_snapshot("Auto_cbio_genes.rds",function() {
  bind_rows(lapply(unique(selected_genes$gene),function(gene) api_json(paste0("/genes/",gene))))
})
if (!all(unique(selected_genes$gene) %in% gene_map$hugoGeneSymbol)) stop("Unmapped selected HGNC genes.")
cna_sample_ids <- source_snapshot("Auto_cbio_cna_sample_ids.rds",function() api_json("/sample-lists/esca_tcga_gdc_cna/sample-ids"))
gene_response <- source_snapshot("Auto_cbio_gene_cna_response.rds",function() {
  api_json(paste0("/molecular-profiles/",profile_id,"/molecular-data/fetch?projection=DETAILED"),
           list(entrezGeneIds=as.integer(gene_map$entrezGeneId),sampleIds=as.character(cna_sample_ids)))
})
if (!all(c("sampleId","entrezGeneId","value") %in% names(gene_response))) stop("Unexpected molecular-data response schema.")
if (!any(gene_response$value==0)) stop("No explicit neutral calls returned; missing calls must not be imputed as zero.")
parent <- readRDS(parent_source_path)
sc_results <- readRDS(sc_path)
features <- as.data.frame(sc_results$features)
manual_cells <- fread(cell_path,data.table=FALSE)
counts_by_clone <- manual_cells |> count(sample,subclone,name="manual_n")
features <- features |> left_join(counts_by_clone,by=c("sample","subclone"))
if (anyNA(features$manual_n) || any(features$n_cells!=features$manual_n)) stop("Cached subclone membership does not match manual InferCNA cell table.")
meta_epi <- readRDS(SCREF_META_FULL_EPI_RDS)
if (!all(manual_cells$cell %in% rownames(meta_epi)) ||
    !identical(as.character(meta_epi[manual_cells$cell,"orig.ident"]),as.character(manual_cells$sample))) stop("orig.ident provenance mismatch.")
rm(meta_epi)
features <- features |> filter(subclone!="Unresolved")
tcga_meta <- parent$tcga_scores |> select(sample_barcode,sample_key,case_barcode,all_of(mp_order))
if (anyDuplicated(tcga_meta$sample_key)) stop("Ambiguous RNA-to-CNA TCGA sample keys.")
tcga_gene <- gene_response |> as_tibble() |>
  transmute(cna_sample_id=sampleId,entrezGeneId=as.integer(entrezGeneId),cna=as.numeric(value),sample_key=substr(sampleId,1,15)) |>
  left_join(gene_map |> select(entrezGeneId,gene=hugoGeneSymbol),by="entrezGeneId") |>
  inner_join(tcga_meta,by="sample_key")
if (any(!tcga_gene$cna %in% -2:2) || anyDuplicated(tcga_gene[c("sample_key","gene")])) stop("Invalid or duplicated gene-level CNA values.")
write.csv(tcga_gene,file.path(tiers[["tables"]],"Auto_tcga_gene_cna_long.csv"),row.names=FALSE)
saveRDS(list(data=tcga_gene,profile=profile,genes=gene_map),file.path(tiers[["intermediate"]],"Auto_tcga_gene_cna.rds"))
# Store source counts-derived activity and exact orig.ident/subclone membership for all replotting.
expression_paths <- file.path(SCREF_REF_OUTS_DIR,"by_samples",unique(manual_cells$sample),paste0(unique(manual_cells$sample),"_epi.rds"))
tcga_require_files(expression_paths)
expression_cache <- tcga_cache(file.path(tiers[["intermediate"]],"Auto_sc_manual_subclone_gene_expression.rds"),
  inputs=c(cell_path,expression_paths,oac_path,SCREF_MP_GENES_RDS), build_fun=function() {
    rows <- list()
    for (sample_id in unique(manual_cells$sample)) {
      message("RNA expression: ",sample_id)
      epi <- readRDS(file.path(SCREF_REF_OUTS_DIR,"by_samples",sample_id,paste0(sample_id,"_epi.rds")))
      sample_cells <- manual_cells[manual_cells$sample==sample_id & manual_cells$subclone!="Unresolved",]
      counts <- GetAssayData(epi,assay="RNA",layer="counts")
      if (!all(sample_cells$cell %in% colnames(counts))) stop("Manual cells missing from RNA counts: ",sample_id)
      if (!identical(as.character(epi@meta.data[sample_cells$cell,"orig.ident"]),as.character(sample_cells$sample))) stop("RNA orig.ident mismatch.")
      counts <- counts[,sample_cells$cell,drop=FALSE]
      depth <- Matrix::colSums(counts)
      if (any(depth<=0) || any(!is.finite(depth))) stop("Invalid RNA library sizes.")
      present <- intersect(unique(selected_genes$gene),rownames(counts))
      expr <- log1p(sweep(as.matrix(counts[present,,drop=FALSE]),2,depth,"/")*10000)
      for (clone in unique(sample_cells$subclone)) {
        ids <- sample_cells$cell[sample_cells$subclone==clone]
        means <- rowMeans(expr[,ids,drop=FALSE])
        rows[[length(rows)+1L]] <- data.frame(sample=sample_id,subclone=clone,subclone_id=paste(sample_id,clone,sep="::"),
          gene=unique(selected_genes$gene),expression=unname(means[unique(selected_genes$gene)]),n_cells=length(ids))
      }
      rm(epi,counts,expr)
      gc(verbose=FALSE)
    }
    bind_rows(rows)
  })
sc_gene <- expression_cache$value |> inner_join(features |> select(subclone_id,all_of(paste0("mp__",mp_order))),by="subclone_id")
coverage <- bind_rows(sc_gene |> group_by(gene) |> summarise(dataset="scATLAS RNA",n_total=n(),n_measured=sum(is.finite(expression)),.groups="drop"),
                     tcga_gene |> group_by(gene) |> summarise(dataset="TCGA gene GISTIC",n_total=n(),n_measured=sum(is.finite(cna)),.groups="drop"))
write.csv(coverage,file.path(tiers[["tables"]],"Auto_gene_measurement_coverage.csv"),row.names=FALSE)
spearman <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if (sum(ok)<5L || length(unique(x[ok]))<2L || length(unique(y[ok]))<2L) return(data.frame(n=sum(ok),rho=NA_real_,p_value=NA_real_))
  test <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  data.frame(n=sum(ok),rho=unname(test$estimate),p_value=test$p.value)
}
arm_ids <- intersect(rownames(sc_results$arm_matrix),features$subclone_id)
sc_arm <- as.data.frame(sc_results$arm_matrix[arm_ids,,drop=FALSE])
sc_arm$subclone_id <- rownames(sc_arm)
sc_arm <- sc_arm |> inner_join(features |> select(subclone_id,sample,n_cells,all_of(paste0("mp__",mp_order))),by="subclone_id")
tcga_arm <- read.csv(file.path(parent_dir,"tables","Auto_tcga_arm_cna_long.csv"),check.names=FALSE)
if (anyDuplicated(tcga_arm[c("sample_barcode","arm")])) stop("Duplicated TCGA sample-arm rows.")
tcga_arm <- tcga_arm |> inner_join(tcga_meta,by=c("sample_barcode","sample_key","case_barcode"))
arm_order <- colnames(sc_results$arm_matrix)
arm_tests <- bind_rows(lapply(arm_order,function(arm) {
  bind_rows(lapply(mp_order,function(mp) {
    sc <- spearman(sc_arm[[arm]],sc_arm[[paste0("mp__",mp)]])
    tc <- tcga_arm[tcga_arm$arm==arm,]
    bulk <- spearman(tc$arm_mean,tc[[mp]])
    bind_rows(cbind(data.frame(dataset="scATLAS manual InferCNA",arm=arm,mp=mp),sc),
              cbind(data.frame(dataset="TCGA segment CNA",arm=arm,mp=mp),bulk))
  }))
})) |> group_by(dataset) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
gene_tests <- bind_rows(lapply(unique(selected_genes$gene),function(gene) {
  sc <- sc_gene[sc_gene$gene==gene,]
  tc <- tcga_gene[tcga_gene$gene==gene,]
  bind_rows(lapply(mp_order,function(mp) {
    bind_rows(cbind(data.frame(dataset="scATLAS gene RNA expression",gene=gene,mp=mp),spearman(sc$expression,sc[[paste0("mp__",mp)]])),
              cbind(data.frame(dataset="TCGA gene CNA (GISTIC)",gene=gene,mp=mp),spearman(tc$cna,tc[[mp]])))
  }))
})) |> group_by(dataset) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup() |>
  mutate(gene_in_mp=mapply(function(gene,mp) gene %in% mp_sets[[mp]],gene,mp))
####################
# Require concordance; equally weight ranks of weaker-dataset significance and effect strength.
rank_pairs <- function(tests,key,sc_label,tc_label) {
  sc <- tests |> filter(dataset==sc_label) |> select(all_of(key),mp,sc_rho=rho,sc_p=p_value,sc_fdr=fdr,sc_n=n)
  tc <- tests |> filter(dataset==tc_label) |> select(all_of(key),mp,tcga_rho=rho,tcga_p=p_value,tcga_fdr=fdr,tcga_n=n)
  inner_join(sc,tc,by=c(key,"mp")) |>
    mutate(same_direction=is.finite(sc_rho)&is.finite(tcga_rho)&sign(sc_rho)==sign(tcga_rho)&sc_rho!=0,
      testable=is.finite(sc_rho)&is.finite(tcga_rho)&is.finite(sc_p)&is.finite(tcga_p),
      both_fdr=is.finite(sc_fdr)&is.finite(tcga_fdr)&sc_fdr<0.05&tcga_fdr<0.05,
      either_fdr=(is.finite(sc_fdr)&sc_fdr<0.05)|(is.finite(tcga_fdr)&tcga_fdr<0.05),
      significance_tier=ifelse(both_fdr,2L,ifelse(either_fdr,1L,0L)),
      max_nominal_p=pmax(sc_p,tcga_p),min_nominal_p=pmin(sc_p,tcga_p),
      min_abs_rho=pmin(abs(sc_rho),abs(tcga_rho)),eligible=testable&same_direction,
      significance_rank=rank(ifelse(eligible,max_nominal_p,NA_real_),ties.method="min",na.last="keep"),
      strength_rank=rank(ifelse(eligible,-min_abs_rho,NA_real_),ties.method="min",na.last="keep"),
      combined_rank_score=significance_rank+strength_rank) |>
    arrange(desc(eligible),desc(significance_tier),combined_rank_score,max_nominal_p,desc(min_abs_rho),min_nominal_p,.data[[key]],match(mp,mp_order)) |>
    mutate(joint_rank=row_number())
}
arm_ranking <- rank_pairs(arm_tests,"arm","scATLAS manual InferCNA","TCGA segment CNA")
gene_ranking <- rank_pairs(gene_tests,"gene","scATLAS gene RNA expression","TCGA gene CNA (GISTIC)")
top_arm <- arm_ranking |> filter(eligible) |> slice_head(n=10)
top_gene <- gene_ranking |> filter(eligible) |> slice_head(n=10)
####################
# Sample-aggregated sensitivity keeps subclones from one tumour from being counted as independent samples.
sensitivity <- bind_rows(lapply(unique(selected_genes$gene),function(gene) {
  d <- sc_gene |> filter(.data$gene==.env$gene) |> group_by(sample) |>
    summarise(expression=weighted.mean(expression,n_cells,na.rm=TRUE),
      across(all_of(paste0("mp__",mp_order)),~weighted.mean(.x,n_cells)),.groups="drop")
  bind_rows(lapply(mp_order,function(mp)cbind(data.frame(analysis="gene RNA",feature=gene,mp=mp),spearman(d$expression,d[[paste0("mp__",mp)]]))))
}), lapply(arm_order,function(arm) {
  d <- sc_arm |> group_by(sample) |> summarise(value=weighted.mean(.data[[arm]],n_cells),
    across(all_of(paste0("mp__",mp_order)),~weighted.mean(.x,n_cells)),.groups="drop")
  bind_rows(lapply(mp_order,function(mp)cbind(data.frame(analysis="arm CNA",feature=arm,mp=mp),spearman(d$value,d[[paste0("mp__",mp)]]))))
})) |> group_by(analysis) |> mutate(fdr=p.adjust(p_value,"BH")) |> ungroup()
for (name in c("arm_tests","gene_tests","arm_ranking","gene_ranking","top_arm","top_gene","sensitivity")) {
  filename <- c(arm_tests="Auto_arm_mp_tests.csv",gene_tests="Auto_gene_mp_tests.csv",arm_ranking="Auto_arm_pair_ranking.csv",gene_ranking="Auto_gene_pair_ranking.csv",top_arm="Auto_top_arm_pairs.csv",top_gene="Auto_top_gene_pairs.csv",sensitivity="Auto_sc_sample_level_sensitivity.csv")[[name]]
  write.csv(get(name),file.path(tiers[["tables"]],filename),row.names=FALSE)
}
results <- list(selected_genes=selected_genes,gene_map=gene_map,profile=profile,sc_gene=sc_gene,tcga_gene=tcga_gene,
  sc_arm=sc_arm,tcga_arm=tcga_arm,arm_tests=arm_tests,gene_tests=gene_tests,arm_ranking=arm_ranking,gene_ranking=gene_ranking,
  top_arm=top_arm,top_gene=top_gene,sensitivity=sensitivity,mp_sets=mp_sets)
saveRDS(results,file.path(tiers[["intermediate"]],"Auto_targeted_cna_results.rds"))
####################
# Match the reference all-arm rectangles: identical tile geometry, fonts and saturated rho scale.
build_targeted_tile_plot <- function(plot_df, title_text, fill_title,
                                    x_text_size = 15, y_text_size = 14,
                                    base_size = 20) {
  max_abs <- max(c(0,abs(plot_df$effect_for_fill)), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  if (!"star_label" %in% colnames(plot_df)) {
      plot_df <- plot_df |>
        mutate(
          star_label = ifelse(!is.na(.data$p_adj) & .data$p_adj < 0.05, "\u2605", "")
        )
  }

  ggplot(plot_df, aes(x = .data$arm, y = .data$mp_label)) +
    geom_tile(aes(fill = .data$effect_for_fill), colour = "grey85", linewidth = 0.4, width = 1, height = 1) +
    geom_text(
      data = plot_df |> filter(.data$star_label != ""),
      aes(label = .data$star_label), colour = "black", size = 8, vjust = 0.4, hjust = 0.5
    ) +
    facet_grid(. ~ dataset) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-0.25, 0.25), oob = scales::squish, na.value = "white", 
      name = fill_title
    ) +
    labs(title = title_text, x = "Recurrent OAC gene events", y = NULL,
         caption = NULL) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = x_text_size, colour = "black"),
      axis.text.y = element_text(size = y_text_size, colour = "black"),
      axis.title.x = element_text(size = base_size + 1, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = base_size + 1, face = "bold", margin = margin(r = 8)),
      strip.text = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      panel.spacing.x = grid::unit(1.5, "lines"),
      plot.margin = margin(16, 32, 16, 20),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank()
    )
}

gene_plot <- gene_tests |> left_join(selected_genes |> select(gene,direction,event_label,event_order),by="gene") |>
  left_join(gene_ranking |> select(gene,mp,same_direction,tcga_fdr),by=c("gene","mp")) |>
  mutate(mp_label=factor(mp_labels[mp],levels=mp_labels[mp_order]),
         dataset=factor(dataset,levels=c("scATLAS gene RNA expression","TCGA gene CNA (GISTIC)")),
         tcga_significant=is.finite(tcga_fdr)&tcga_fdr<0.05,
         star_label=ifelse(tcga_significant&same_direction,"\u2605",""))
rectangle_sources <- list();rectangle_pages <- list()
for (direction_label in c("gain","loss")) {
  for (view in c("all trends","TCGA significant only")) {
    df <- gene_plot |> filter(as.character(direction)==direction_label) |>
      mutate(arm=factor(event_label,levels=selected_genes$event_label[as.character(selected_genes$direction)==direction_label]),
        view=.env$view,effect_for_fill=ifelse(.env$view=="all trends"|tcga_significant,rho,NA_real_),
        plot_caption=paste0("Stars: TCGA BH FDR < 0.05 and matching nonzero rho direction in both datasets. ",
          ifelse(.env$view=="all trends","All correlations shown.","Both panels show only TCGA BH FDR < 0.05 pairs; blank pages remain plotted.")))
    rectangle_sources[[length(rectangle_sources)+1L]] <- df
    rectangle_pages[[length(rectangle_pages)+1L]] <- build_targeted_tile_plot(df,
      paste("Top 10 recurrent OAC",direction_label,"genes versus MP activity:",view),"Spearman rho")
  }
}
rectangle_annotations <- bind_rows(rectangle_sources)
write.csv(rectangle_annotations,file.path(tiers[["tables"]],"Auto_gene_rectangle_annotations.csv"),row.names=FALSE)
tcga_save_figure("Auto_oac_gene_mp_spearman_rectangles.pdf",rectangle_pages,width=25,height=10)
####################
paired_page <- function(left,right,title) {
  force(left);force(right);force(title)
  function() gridExtra::arrangeGrob(left,right,ncol=2,
    top=grid::textGrob(str_wrap(title,width=90),gp=grid::gpar(fontsize=18,fontface="bold")))
}
stat_label <- function(row,dataset) {
  prefix <- if(dataset=="sc") "sc_" else "tcga_"
  paste0("rho = ",sprintf("%.2f",row[[paste0(prefix,"rho")]]),
    "; P = ",format.pval(row[[paste0(prefix,"p")]],digits=2),
    "; BH FDR = ",format.pval(row[[paste0(prefix,"fdr")]],digits=2),
    "; n = ",row[[paste0(prefix,"n")]])
}
empty_page <- function(label) {
  ggplot() + annotate("text",x=0,y=0,label=str_wrap(label,width=80),size=6) +
    xlim(-1,1)+ylim(-1,1)+theme_void()
}
  ####################
  box <- function(df,title,statistic) {
    counts <- table(df$group)
    present <- names(counts)[counts>0]
    comparisons <- if (length(present)>=2L) combn(present,2,simplify=FALSE) else list()
    plot <- ggplot(df,aes(group,y,fill=group)) + geom_boxplot(outlier.shape=NA,width=0.62,linewidth=0.8,alpha=0.88)
    if ("n_cells" %in% names(df)) {
      plot <- plot + geom_point(aes(size=n_cells),position=position_jitter(width=0.14,height=0),alpha=0.42,colour="black") +
        scale_size_continuous(range=c(1.8,6))
    } else {
      plot <- plot + geom_point(position=position_jitter(width=0.14,height=0),alpha=0.42,size=3,colour="black")
    }
    plot + ggpubr::stat_compare_means(comparisons=comparisons,method="wilcox.test",label="p.signif",
      symnum.args=list(cutpoints=c(0,0.0001,0.001,0.01,0.05,Inf),symbols=c("****","***","**","*","NS")),
      size=7,vjust=-0.2,step.increase=0.12,tip.length=0.02) +
      scale_fill_manual(values=c(Loss="#2166AC",Neutral="grey72",Gain="#B2182B",Low="#2166AC",Medium="grey72",High="#B2182B",
        "Deep loss"="#053061",Amplified="#67001F"),drop=FALSE) +
      scale_x_discrete(drop=FALSE,labels=setNames(paste0(names(counts),"\n(n=",counts,")"),names(counts))) +
      labs(title=title,subtitle=statistic,x=NULL,y="Mean MP activity") +
      theme_classic(base_size=22)+theme(legend.position="none",plot.margin=margin(14,24,14,24),
        plot.title=element_text(face="bold",size=26),plot.subtitle=element_text(size=16),
        axis.text=element_text(size=22,colour="black"),axis.title=element_text(size=24,face="bold",colour="black"))
  }
  ####################
arm_box_pages <- list();arm_scatter_pages <- list()
for (i in seq_len(nrow(top_arm))) {
  row <- top_arm[i,];arm <- row$arm;mp <- row$mp
  sc <- data.frame(x=sc_arm[[arm]],y=sc_arm[[paste0("mp__",mp)]],n_cells=sc_arm$n_cells)
  tc <- tcga_arm |> filter(.data$arm==.env$arm) |> transmute(x=arm_mean,y=.data[[mp]])
  sc$group <- factor(ifelse(sc$x>=0.05,"Gain",ifelse(sc$x<=-0.05,"Loss","Neutral")),levels=c("Loss","Neutral","Gain"))
  tc$group <- factor(ifelse(tc$x>=parent$selected_threshold,"Gain",ifelse(tc$x<=-parent$selected_threshold,"Loss","Neutral")),levels=c("Loss","Neutral","Gain"))
  title <- paste0(i,". ",arm," versus ",mp_labels[mp])
  arm_box_pages[[i]] <- paired_page(box(sc,"scATLAS manual InferCNA",stat_label(row,"sc")),
    box(tc,"TCGA EAC segment CNA",stat_label(row,"tcga")),title)
  scatter <- function(df,title,statistic,x_label) {
    ggplot(df,aes(x,y)) +geom_point(size=2.5,alpha=0.65,colour="#285E79")+
      labs(title=title,subtitle=statistic,x=x_label,y="Mean MP activity")+
      theme_classic(base_size=16)+theme(plot.margin=margin(14,18,14,18),plot.title=element_text(face="bold",size=19))
  }
  arm_scatter_pages[[i]] <- paired_page(scatter(sc,"scATLAS manual InferCNA",stat_label(row,"sc"),"Mean InferCNA arm signal"),
    scatter(tc,"TCGA EAC segment CNA",stat_label(row,"tcga"),"Mean segment log2 ratio"),title)
}
if (!length(arm_box_pages)) {
  arm_box_pages <- list(empty_page("No arm/MP pair has finite correlation tests in both datasets."))
  arm_scatter_pages <- arm_box_pages
}
tcga_save_figure("Auto_top_arm_mp_paired_boxplots.pdf",arm_box_pages,width=22,height=10)
tcga_save_figure("Auto_top_arm_mp_paired_scatterplots.pdf",arm_scatter_pages,width=22,height=10)
####################
gene_pair_pages <- list();gene_box_sources <- list();gene_bin_thresholds <- list()
for (i in seq_len(nrow(top_gene))) {
  row <- top_gene[i,];gene <- row$gene;mp <- row$mp
  sc <- sc_gene |> filter(.data$gene==.env$gene,is.finite(expression),is.finite(.data[[paste0("mp__",mp)]])) |>
    transmute(sample,subclone_id,x=expression,y=.data[[paste0("mp__",mp)]],n_cells)
  tc <- tcga_gene |> filter(.data$gene==.env$gene,is.finite(cna),is.finite(.data[[mp]])) |>
    transmute(sample=sample_barcode,x=cna,y=.data[[mp]])
  cuts <- quantile(sc$x,c(1/3,2/3),type=7,names=FALSE)
  # Equal-valued expression observations always remain in the same bin, including tied quantiles.
  sc$group <- factor(ifelse(sc$x<=cuts[1],"Low",ifelse(sc$x<=cuts[2],"Medium","High")),levels=c("Low","Medium","High"))
  tc$group <- factor(tc$x,levels=-2:2,labels=c("Deep loss","Loss","Neutral","Gain","Amplified"))
  gene_bin_thresholds[[i]] <- data.frame(gene=gene,mp=mp,lower_tertile=cuts[1],upper_tertile=cuts[2],
    low_n=sum(sc$group=="Low"),medium_n=sum(sc$group=="Medium"),high_n=sum(sc$group=="High"))
  gene_box_sources[[i]] <- bind_rows(mutate(sc,dataset="scATLAS gene RNA expression",gene=.env$gene,mp=.env$mp),
                                    mutate(tc,dataset="TCGA gene CNA (GISTIC)",gene=.env$gene,mp=.env$mp))
  left <- box(sc,"scATLAS: gene RNA",stat_label(row,"sc")) + labs(x=paste(gene,"RNA expression tertiles"),y="Subclone mean MP activity")
  right <- box(tc,"TCGA: gene CNA",stat_label(row,"tcga")) + labs(x=paste(gene,"GISTIC call"),y="Bulk GSVA MP activity")
  overlap <- if (gene %in% mp_sets[[mp]]) " [gene is part of this MP signature]" else ""
  gene_pair_pages[[i]] <- paired_page(left,right,paste0(i,". ",gene," versus ",mp_labels[mp],overlap))
}
if (!length(gene_pair_pages)) gene_pair_pages <- list(empty_page("No curated gene/MP pair has finite same-direction correlations in both datasets."))
gene_box_source <- bind_rows(gene_box_sources)
gene_bin_thresholds <- bind_rows(gene_bin_thresholds)
write.csv(gene_box_source,file.path(tiers[["tables"]],"Auto_gene_pair_boxplot_source.csv"),row.names=FALSE)
write.csv(gene_bin_thresholds,file.path(tiers[["tables"]],"Auto_gene_expression_tertile_thresholds.csv"),row.names=FALSE)
results$rectangle_annotations <- rectangle_annotations
results$gene_box_source <- gene_box_source
results$gene_bin_thresholds <- gene_bin_thresholds
saveRDS(results,file.path(tiers[["intermediate"]],"Auto_targeted_cna_results.rds"))
tcga_save_figure("Auto_top_oac_gene_mp_paired_plots.pdf",gene_pair_pages,width=22,height=10)
####################
saveRDS(list(figures=tcga_figure_bundle),figure_path)
summary <- data.frame(metric=c("manual_cells","manual_samples","manual_subclones","tcga_arm_samples","tcga_gene_samples", "selected_gain_genes","selected_loss_genes","testable_arm_pairs","testable_gene_pairs","same_direction_arm_pairs","same_direction_gene_pairs","plotted_arm_pairs","plotted_gene_pairs","gene_cna_profile"),
  value=c(nrow(manual_cells),length(unique(manual_cells$sample)),nrow(features),length(unique(tcga_arm$sample_barcode)),length(unique(tcga_gene$sample_barcode)),10,10,
    sum(arm_ranking$testable),sum(gene_ranking$testable),sum(arm_ranking$same_direction),sum(gene_ranking$same_direction),nrow(top_arm),nrow(top_gene),profile_id))
write.csv(summary,summary_path,row.names=FALSE)
summary_copy <- tcga_live_path(file.path(SCREF_SUMMARY_DIR,"Auto_targeted_cna_summary.csv"))
write.csv(summary,summary_copy,row.names=FALSE)
stopifnot(file.info(summary_copy)$size<100000,nrow(arm_tests)==2L*length(arm_order)*17L,
          nrow(gene_tests)==2L*20L*17L,!anyDuplicated(arm_tests[c("dataset","arm","mp")]),
          !anyDuplicated(gene_tests[c("dataset","gene","mp")]),all(top_arm$eligible),all(top_gene$eligible),
          nrow(top_arm)==min(10L,sum(arm_ranking$eligible)),nrow(top_gene)==min(10L,sum(gene_ranking$eligible)))
writeLines(c("PASS: manual InferCNA provenance gate, orig.ident and subclone cell counts verified.",
  "PASS: exact current 17 MPs; 10 gain and 10 loss genes; explicit observed gene-level GISTIC calls.",
  "PASS: no missing call imputed as neutral; absent RNA genes retained as unmeasured.",
  "PASS: full-family BH correction, unique pair rows and user-approved selection rule.",
  paste("Testable/plotted arm pairs:",sum(arm_ranking$testable),nrow(top_arm)),
  paste("Testable/plotted gene pairs:",sum(gene_ranking$testable),nrow(top_gene)),
  "scATLAS gene predictor is RNA expression; TCGA gene predictor is DNA-based copy number.",
  "Subclones are nested within patients; sample-level sensitivity is provided.",
  "Gene-in-MP membership is flagged because it can create score-expression coupling."),
  file.path(tiers[["reports"]],"Auto_analysis_validation.txt"))
run <- start_run_summary(script_path,c(cell_path,sc_path,parent_source_path,oac_path,provenance_path,profile_path,expression_paths),
  c(figure_path,summary_path,file.path(tiers[["intermediate"]],"Auto_targeted_cna_results.rds")),
  list(cna_profile=profile_id,top_genes_per_direction=10,pair_ranking="same direction required; both/either FDR support first, then sum of ranks of max nominal P and negative min absolute rho",
       star_rule="TCGA BH FDR <0.05 and same nonzero direction",gene_expression_bins="type-7 tertiles; ties kept together",top_pairs=10,
       sample_unit="scATLAS manual InferCNA subclone; TCGA EAC primary tumour",gene_sc_predictor="mean log1p CP10K RNA"))
run <- add_cache_status(run,expression_cache$label,expression_cache$path,expression_cache$reused)
write_run_summary(run,file.path(tiers[["logs"]],"Auto_targeted_cna_run_summary.rds"))
message("Targeted gene/arm CNA comparison complete: ",out_dir)
