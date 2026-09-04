####################
# Analysis registry:
#   Status: active
#   Script: analysis/shared/scRef_config.R
#   Methodology: not required (configuration only)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# scRef_config.R
#
# Central configuration for downstream scRef analysis scripts.
#
# Purpose:
#   Keep canonical paths, state definitions, plotting defaults, and output-tier
#   conventions in one place. New scripts may source this file directly, or copy
#   the relevant constants only when a fully self-contained HPC job is required.
#
# Output:
#   None. This file defines constants only.
####################

SCREF_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
SCREF_REF_OUTS_DIR <- file.path(SCREF_PROJECT_DIR, "ref_outs")
SCREF_ANALYSIS_DIR <- file.path(SCREF_PROJECT_DIR, "analysis")
SCREF_UPDATES_DIR <- file.path(SCREF_PROJECT_DIR, "updates", "new_updates")
SCREF_SUMMARY_DIR <- file.path(SCREF_UPDATES_DIR, "summaries")
SCREF_EPHEMERAL_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
SCREF_EPHEMERAL_REF_OUTS_DIR <- file.path(SCREF_EPHEMERAL_PROJECT_DIR, "ref_outs")

SCREF_PREFERRED_STATE_METHOD <- "centred refined noreg, Approach B"
SCREF_CENTRED_MP_DIR <- file.path(
  SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "centred", "mp_refinement"
)
SCREF_CENTRED_STATE_DIR <- file.path(
  SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "centred", "state_definition"
)
SCREF_STATE_NOREG_RDS <- file.path(
  SCREF_CENTRED_STATE_DIR, "intermediate", "centred_refined_noreg_states.rds"
)
SCREF_STATE_NOREG_MP_ADJ_RDS <- file.path(
  SCREF_CENTRED_STATE_DIR, "intermediate", "centred_refined_noreg_mp_adj.rds"
)
SCREF_STATE_NOREG_GROUP_MAX_RDS <- file.path(
  SCREF_CENTRED_STATE_DIR, "intermediate", "centred_refined_noreg_group_max.rds"
)
# Compatibility alias for scripts that historically distinguished final from noreg states.
SCREF_FINAL_STATE_RDS <- SCREF_STATE_NOREG_RDS

SCREF_EPI_RDS <- file.path(SCREF_REF_OUTS_DIR, "EAC_Ref_epi.rds")
SCREF_META_FULL_EPI_RDS <- file.path(SCREF_REF_OUTS_DIR, "meta_full_epi.rds")
SCREF_LEGACY_MP_OBJECT_RDS <- file.path(
  SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "centred", "geneNMF_metaprograms_nMP_19.rds"
)
SCREF_MP_GENES_RDS <- file.path(
  SCREF_CENTRED_MP_DIR, "intermediate", "merged_refined_mp_genes.rds"
)
SCREF_MP_UCELL_RDS <- file.path(
  SCREF_CENTRED_MP_DIR, "intermediate", "merged_refined_ucell_scores.rds"
)
SCREF_MP_GROUPING_CSV <- file.path(
  SCREF_CENTRED_MP_DIR, "tables", "centred_refined_mp_state_grouping.csv"
)
SCREF_3CA_UCELL_RDS <- file.path(SCREF_REF_OUTS_DIR, "UCell_3CA_MPs.rds")
SCREF_CELL_CYCLE_GENE_CSV <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
SCREF_CLINICAL_XLSX <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx"

SCREF_OUTPUT_TIERS <- c("intermediate", "tables", "figures", "logs", "reports")

SCREF_DEFAULT_PLOT <- list(
  dpi = 300,
  pdf_width = 12,
  pdf_height = 8,
  slide_width = 13.333,
  slide_height = 7.5,
  heatmap_width = 16,
  heatmap_height = 10,
  base_size = 14,
  axis_text_size = 12,
  legend_text_size = 12,
  legend_title_size = 13,
  point_size = 2.8,
  label_size = 4
)

SCREF_MP_DESCRIPTIONS <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "Single-nucleus cell cycle",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Inflammatory-reactive columnar epithelium",
  "MP11+" = "Epithelial type I interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP18b" = "Mucous-secretory differentiation",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP17" = "MHC-II glandular progenitor",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP12" = "Hypoxic-inflammatory adaptive plasticity",
  "MP15" = "Cancer-cell immune mimicry"
)

SCREF_CC_MPS <- c("MP1", "MP5", "MP13+")
SCREF_STATE_GROUPS <- list(
  "Classic proliferation" = c("MP2+"),
  "Squamous-to-intestinal" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "Glandular-to-intestinal" = c("MP18b", "MP16", "MP17", "MP8b", "MP8+"),
  "Stress-adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

SCREF_PRIMARY_STATE_ORDER <- c(
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)

SCREF_FINAL_EXTRA_STATE_ORDER <- character(0)
SCREF_TECHNICAL_STATE_ORDER <- c("Unresolved", "Hybrid")

SCREF_STATE_COLOURS <- c(
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

SCREF_STATE_THRESHOLDS <- list(
  approach_b_unresolved_min_z = 0.5,
  hybrid_gap_threshold = 0.3,
  min_cells_per_state_for_pseudotime = 20,
  min_total_cells_for_pseudotime = 80,
  max_cells_per_state_heatmap = 500,
  mp_silhouette_min = 0
)

SCREF_METADATA_COLUMNS <- list(
  sample = "orig.ident",
  study = "study",
  malignancy = "malignancy",
  celltype = "celltype_update",
  treatment = "Treatment",
  response = "Clinical response"
)
