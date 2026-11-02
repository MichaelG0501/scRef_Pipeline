####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/export_scatlas_tme_signatures.R
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Export silhouette-filtered ranked genes and summary for the
#     non-malignant GeneNMF MP panels (fibroblast, endothelial, macrophage, nk,
#     plasma, cd4, cd8) for spatial TME MP scoring.
#   Inputs:
#     - ref_outs/nmf_{fibroblast,endothelial,macrophage,nk,plasma,cd4,cd8}/MP_outs_default.rds
#   Outputs: <output_dir>/Auto_scATLAS_tme_mp_gene_ranked.csv,
#     <output_dir>/Auto_scATLAS_tme_mp_summary.csv
#   Cache/replot: deterministic table export.
#   Run: called by map_scatlas_states_xenium.sh or Rscript
#   Conda env: dmtcp
####################

args <- commandArgs(trailingOnly = TRUE)

wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
output_dir <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(wd, "ref_outs/xenium_scatlas_states")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tme_mp_cfg <- data.frame(
  compartment = c("fibroblast", "endothelial", "macrophage", "nk", "plasma", "cd4", "cd8"),
  annotation_celltype = c("fibroblast", "endothelial", "macrophage", "nk.cell", "plasma", "t.cell", "t.cell"),
  mp_path = file.path(
    wd, "ref_outs",
    c("nmf_fibroblast", "nmf_endothelial", "nmf_macrophage", "nmf_nk", "nmf_plasma", "nmf_cd4", "nmf_cd8"),
    "MP_outs_default.rds"
  ),
  stringsAsFactors = FALSE
)

filter_gene_sets_by_silhouette <- function(mp_outs) {
  mp_genes <- mp_outs$metaprograms.genes
  metrics <- mp_outs$metaprograms.metrics
  if (is.null(mp_genes) || is.null(metrics$silhouette) ||
      is.null(names(mp_genes)) || anyDuplicated(names(mp_genes)) ||
      is.null(rownames(metrics)) || anyDuplicated(rownames(metrics)) ||
      !setequal(names(mp_genes), rownames(metrics))) {
    stop("GeneNMF signatures and metric row names must match exactly")
  }
  sil <- setNames(metrics$silhouette, rownames(metrics))[names(mp_genes)]
  if (any(!is.finite(sil))) stop("Non-finite GeneNMF silhouette")
  keep <- sil >= 0
  list(genes = mp_genes[keep], sil = sil[keep])
}

ranked_rows <- list()
summary_rows <- list()

for (i in seq_len(nrow(tme_mp_cfg))) {
  comp <- tme_mp_cfg$compartment[i]
  celltype <- tme_mp_cfg$annotation_celltype[i]
  path <- tme_mp_cfg$mp_path[i]
  
  if (!file.exists(path)) {
    stop("Missing GeneNMF object: ", path)
  }
  
  mp_outs <- readRDS(path)
  filtered <- filter_gene_sets_by_silhouette(mp_outs)
  retained_genes <- filtered$genes
  retained_sil <- filtered$sil
  
  for (mp_name in names(retained_genes)) {
    genes <- unique(as.character(retained_genes[[mp_name]]))
    feature_id <- paste(comp, mp_name, sep = "__")
    sil_val <- as.numeric(retained_sil[mp_name])
    
    ranked_rows[[length(ranked_rows) + 1L]] <- data.frame(
      compartment = comp,
      annotation_celltype = celltype,
      mp = mp_name,
      feature_id = feature_id,
      gene = genes,
      rank = seq_along(genes),
      silhouette = sil_val,
      stringsAsFactors = FALSE
    )
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = comp,
      annotation_celltype = celltype,
      mp = mp_name,
      feature_id = feature_id,
      n_genes = length(genes),
      silhouette = sil_val,
      stringsAsFactors = FALSE
    )
  }
}

ranked_df <- do.call(rbind, ranked_rows)
summary_df <- do.call(rbind, summary_rows)

write.csv(
  ranked_df,
  file.path(output_dir, "Auto_scATLAS_tme_mp_gene_ranked.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  summary_df,
  file.path(output_dir, "Auto_scATLAS_tme_mp_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Exported TME MP signatures to: ", output_dir)
message("Total TME MPs: ", nrow(summary_df), " across ", length(unique(summary_df$compartment)), " compartments")
