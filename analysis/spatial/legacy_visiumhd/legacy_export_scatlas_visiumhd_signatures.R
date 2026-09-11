#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_export_scatlas_visiumhd_signatures.R
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1) {
  args[1]
} else {
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/spatial/visium_hd_outs"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

refined_genes_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds"

merged_mp_genes <- readRDS(refined_genes_path)

cc_mps <- c("MP1", "MP5", "MP13+")
state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP2x"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)
excluded_mps <- c("MP11c", "MP18a")

mp_desc_map <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP2x" = "Wnt-active glandular stem/progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry",
  "MP11c" = "Excluded",
  "MP18a" = "Excluded"
)

# Filter for only MPs that have defined genes
valid_mps <- names(merged_mp_genes)[vapply(merged_mp_genes, length, integer(1)) > 0]

ranked_rows <- do.call(
  rbind,
  lapply(valid_mps, function(mp_name) {
    desc <- unname(mp_desc_map[mp_name])
    if (is.na(desc)) desc <- mp_name
    data.frame(
      mp = mp_name,
      gene = merged_mp_genes[[mp_name]],
      rank = seq_along(merged_mp_genes[[mp_name]]),
      description = desc,
      is_cc = mp_name %in% cc_mps,
      stringsAsFactors = FALSE
    )
  })
)

state_rows <- do.call(
  rbind,
  lapply(names(state_groups), function(state_name) {
    mps <- state_groups[[state_name]]
    desc <- unname(mp_desc_map[mps])
    desc[is.na(desc)] <- mps[is.na(desc)]
    data.frame(
      state = state_name,
      mp = mps,
      mp_description = desc,
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

# generate logical plot order for MPs
all_available_mps <- valid_mps
mp_numbers <- as.integer(sub("^MP", "", sub("[a-z\\+]*$", "", all_available_mps)))
mp_order <- all_available_mps[order(mp_numbers, all_available_mps)]

mp_order_df <- data.frame(
  mp = mp_order,
  description = unname(mp_desc_map[mp_order]),
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
    paste0("refined_genes_path=", refined_genes_path),
    paste0("cc_mps=", paste(cc_mps, collapse = ",")),
    paste0("mp_order=", paste(mp_order, collapse = ",")),
    "default_top_n_for_spatial=100"
  ),
  con = file.path(output_dir, "Auto_scATLAS_signature_metadata.txt")
)
