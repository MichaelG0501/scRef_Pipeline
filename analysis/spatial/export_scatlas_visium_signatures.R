####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/export_scatlas_visium_signatures.R
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Export ranked genes, full descriptions, state groups, and plot
#     order for the final centred 17-MP panel used by Visium/Xenium mapping.
#   Inputs:
#     - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#     - ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   Outputs: <output_dir>/Auto_scATLAS_{mp_gene_ranked,mp_signature_summary,state_groups,mp_order}.csv
#     and <output_dir>/Auto_scATLAS_signature_metadata.txt.
#   Cache/replot: deterministic table export; overwrites the selected live output directory.
#   Run: called by map_scatlas_states_visium.sh or map_scatlas_states_xenium.sh.
#   Conda env: dmtcp
####################

args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1) {
  args[1]
} else {
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/spatial_signatures/visium_initial"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

geneNMF_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds"
grouping_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv"
mp.genes <- readRDS(geneNMF_path)
grouping <- read.csv(grouping_path, check.names = FALSE)

mp_descriptions <- setNames(grouping$description, grouping$mp)

state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

cc_mps <- c("MP1", "MP5", "MP13+")

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

mp_tree_order <- grouping$mp[grouping$mp %in% names(mp.genes)]
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
    "excluded_mps=MP2x,MP11c,MP18a",
    paste0("cc_mps=", paste(cc_mps, collapse = ",")),
    paste0("mp_order=", paste(mp_order, collapse = ",")),
    "default_top_n_for_spatial=100"
  ),
  con = file.path(output_dir, "Auto_scATLAS_signature_metadata.txt")
)
