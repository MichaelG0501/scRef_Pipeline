####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/export_scatlas_visium_signatures.R
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1) {
  args[1]
} else {
  "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_visium_initial"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

geneNMF_path <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"

geneNMF.metaprograms <- readRDS(geneNMF_path)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

mp_descriptions <- c(
  "MP1" = "G2M Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP5" = "Epithelial IFN Resp.",
  "MP7" = "DNA Damage Repair",
  "MP8" = "Intestinal Diff.",
  "MP9" = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)"
)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15")
)

cc_mps <- c("MP1", "MP7", "MP9")

ranked_rows <- do.call(
  rbind,
  lapply(names(mp.genes), function(mp_name) {
    data.frame(
      mp = mp_name,
      gene = mp.genes[[mp_name]],
      rank = seq_along(mp.genes[[mp_name]]),
      description = unname(mp_descriptions[mp_name]),
      is_cc = mp_name %in% cc_mps,
      stringsAsFactors = FALSE
    )
  })
)

state_rows <- do.call(
  rbind,
  lapply(names(state_groups), function(state_name) {
    data.frame(
      state = state_name,
      mp = state_groups[[state_name]],
      mp_description = unname(mp_descriptions[state_groups[[state_name]]]),
      stringsAsFactors = FALSE
    )
  })
)

signature_summary <- aggregate(
  gene ~ mp + description + is_cc,
  data = ranked_rows,
  FUN = length
)
colnames(signature_summary)[colnames(signature_summary) == "gene"] <- "n_genes"

orig_tree_order <- geneNMF.metaprograms$programs.tree$order
orig_clusters <- geneNMF.metaprograms$programs.clusters[orig_tree_order]
mp_tree_order <- paste0("MP", unique(orig_clusters))
mp_tree_order <- mp_tree_order[mp_tree_order %in% names(mp.genes)]
state_ordered_mps <- unlist(state_groups, use.names = FALSE)
mp_order <- unique(c(
  mp_tree_order[mp_tree_order %in% cc_mps],
  state_ordered_mps[state_ordered_mps %in% names(mp.genes)],
  mp_tree_order
))

mp_order_df <- data.frame(
  mp = mp_order,
  description = unname(mp_descriptions[mp_order]),
  is_cc = mp_order %in% cc_mps,
  plot_order = seq_along(mp_order),
  stringsAsFactors = FALSE
)

write.csv(
  ranked_rows,
  file.path(output_dir, "Auto_scATLAS_mp_gene_ranked.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  signature_summary,
  file.path(output_dir, "Auto_scATLAS_mp_signature_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  state_rows,
  file.path(output_dir, "Auto_scATLAS_state_groups.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  mp_order_df,
  file.path(output_dir, "Auto_scATLAS_mp_order.csv"),
  row.names = FALSE,
  quote = TRUE
)

writeLines(
  c(
    paste0("geneNMF_path=", geneNMF_path),
    paste0("negative_silhouette_filtered=", paste(paste0("MP", bad_mps), collapse = ",")),
    paste0("cc_mps=", paste(cc_mps, collapse = ",")),
    paste0("mp_order=", paste(mp_order, collapse = ",")),
    "default_top_n_for_spatial=100"
  ),
  con = file.path(output_dir, "Auto_scATLAS_signature_metadata.txt")
)
