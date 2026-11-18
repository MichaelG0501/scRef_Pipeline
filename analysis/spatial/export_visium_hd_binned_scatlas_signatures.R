#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/export_visium_hd_binned_scatlas_signatures.R
#   Description: Export current centred-refined scATLAS MP signatures and
#     Approach-B state groups for malignant Visium HD mapping in both
#     representations. Current labels come from analysis/shared/scRef_config.R.
#     First argument selects the exact live output directory; the segmented
#     launcher uses state_mapping_segmented/intermediate/signatures/.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_state_mapping_methodology.md
#   Inputs:
#     ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/
#       merged_refined_mp_genes.rds
#   Outputs:
#     intermediate/signatures/: ranked MP genes, state groups, MP order,
#       signature summary, and provenance metadata under
#       ref_outs/visium_hd_outs/state_mapping/
#   Cache/replot: inexpensive; always rebuilt from the live MP gene object.
#   Run: Rscript analysis/spatial/export_visium_hd_binned_scatlas_signatures.R
#   Environment: dmtcp
####################

####################
args <- commandArgs(trailingOnly = TRUE)
wd <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(wd)

output_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(wd, "ref_outs", "visium_hd_outs", "state_mapping", "intermediate", "signatures")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genes_path <- file.path(
  wd, "ref_outs", "Metaprogrammes_Results", "centred", "mp_refinement",
  "intermediate", "merged_refined_mp_genes.rds"
)
state_definition_path <- file.path(
  wd, "analysis", "metaprograms", "centred",
  "06_centred_refined_state_definition_noreg.R"
)
if (!file.exists(genes_path)) stop("Missing centred-refined MP genes: ", genes_path)
if (!file.exists(state_definition_path)) stop("Missing centred-refined state definition: ", state_definition_path)
mp_genes <- readRDS(genes_path)

####################
# Use only the current central state/description contract for both spatial
# representations. No legacy state labels or descriptions are exported.
####################
source(file.path(wd, "analysis/shared/scRef_config.R"))
cc_mps <- SCREF_CC_MPS
excluded_mps <- c("MP2x", "MP11c", "MP18a")
state_groups <- SCREF_STATE_GROUPS
mp_descriptions <- SCREF_MP_DESCRIPTIONS
if (any(names(mp_genes) %in% excluded_mps)) stop("Excluded MP in current refined genes")
excluded_mps <- character(0)
####################

required_mps <- unique(c(cc_mps, unlist(state_groups, use.names = FALSE), excluded_mps))
missing_mps <- setdiff(required_mps, names(mp_genes))
if (length(missing_mps)) stop("Centred-refined MP object lacks: ", paste(missing_mps, collapse = ", "))
mp_genes <- mp_genes[required_mps]
if (any(vapply(mp_genes, length, integer(1)) == 0L)) stop("One or more MP signatures are empty")

ranked_rows <- do.call(rbind, lapply(names(mp_genes), function(mp_name) {
  data.frame(
    mp = mp_name,
    gene = as.character(mp_genes[[mp_name]]),
    rank = seq_along(mp_genes[[mp_name]]),
    description = unname(mp_descriptions[[mp_name]]),
    is_cc = mp_name %in% cc_mps,
    is_state_defining = mp_name %in% unlist(state_groups, use.names = FALSE),
    is_excluded = mp_name %in% excluded_mps,
    stringsAsFactors = FALSE
  )
}))

state_rows <- do.call(rbind, lapply(names(state_groups), function(state_name) {
  mps <- state_groups[[state_name]]
  data.frame(
    state = state_name,
    mp = mps,
    mp_description = unname(mp_descriptions[mps]),
    stringsAsFactors = FALSE
  )
}))

signature_summary <- aggregate(
  gene ~ mp + description + is_cc + is_state_defining + is_excluded,
  data = ranked_rows,
  FUN = length
)
names(signature_summary)[names(signature_summary) == "gene"] <- "n_genes"
mp_order <- c(cc_mps, unlist(state_groups, use.names = FALSE), excluded_mps)
mp_order_table <- data.frame(
  mp = mp_order,
  description = unname(mp_descriptions[mp_order]),
  is_cc = mp_order %in% cc_mps,
  is_state_defining = mp_order %in% unlist(state_groups, use.names = FALSE),
  is_excluded = mp_order %in% excluded_mps,
  plot_order = seq_along(mp_order),
  stringsAsFactors = FALSE
)

write.csv(ranked_rows, file.path(output_dir, "Auto_scATLAS_centred_refined_mp_gene_ranked.csv"), row.names = FALSE)
write.csv(signature_summary, file.path(output_dir, "Auto_scATLAS_centred_refined_mp_signature_summary.csv"), row.names = FALSE)
write.csv(state_rows, file.path(output_dir, "Auto_scATLAS_centred_refined_state_groups.csv"), row.names = FALSE)
write.csv(mp_order_table, file.path(output_dir, "Auto_scATLAS_centred_refined_mp_order.csv"), row.names = FALSE)
writeLines(
  c(
    paste0("source=", genes_path),
    paste0("source_md5=", unname(tools::md5sum(genes_path))),
    paste0("state_definition_source=", state_definition_path),
    paste0("state_definition_source_md5=", unname(tools::md5sum(state_definition_path))),
    "negative_silhouette_filter=completed upstream before centred MP refinement",
    paste0("cc_mps=", paste(cc_mps, collapse = ";")),
    paste0("excluded_mps=", paste(excluded_mps, collapse = ";")),
    paste0("state_groups=", paste(names(state_groups), collapse = ";")),
    "spatial_top_n=100",
    "state_threshold=0.5",
    "hybrid_gap=0.3"
  ),
  file.path(output_dir, "Auto_scATLAS_centred_refined_signature_metadata.txt")
)
####################
