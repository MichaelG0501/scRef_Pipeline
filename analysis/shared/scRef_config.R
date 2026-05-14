####################
# Analysis registry:
#   Status: active
#   Script: analysis/shared/scRef_config.R
#   Methodology: analysis/methodology/README.md
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
# Methodology:
#   analysis/methodology/README.md
#
# Output:
#   None. This file defines constants only.
####################

SCREF_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
SCREF_REF_OUTS_DIR <- file.path(SCREF_PROJECT_DIR, "ref_outs")
SCREF_ANALYSIS_DIR <- file.path(SCREF_PROJECT_DIR, "analysis")
SCREF_UPDATES_DIR <- file.path(SCREF_PROJECT_DIR, "updates", "new_updates")
SCREF_SUMMARY_DIR <- file.path(SCREF_UPDATES_DIR, "summaries")

SCREF_PREFERRED_STATE_METHOD <- "Approach B, noreg"
SCREF_STATE_NOREG_RDS <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_states_B.rds")
SCREF_STATE_NOREG_MP_ADJ_RDS <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_mp_adj.rds")
SCREF_STATE_NOREG_GROUP_MAX_RDS <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_group_max.rds")
SCREF_FINAL_STATE_RDS <- file.path(SCREF_REF_OUTS_DIR, "Auto_final_states.rds")

SCREF_EPI_RDS <- file.path(SCREF_REF_OUTS_DIR, "EAC_Ref_epi.rds")
SCREF_META_FULL_EPI_RDS <- file.path(SCREF_REF_OUTS_DIR, "meta_full_epi.rds")
SCREF_MP_OBJECT_RDS <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")
SCREF_MP_UCELL_RDS <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "UCell_nMP19_filtered.rds")
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
  "MP1" = "G2M Cell Cycle",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8" = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7" = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi."
)

SCREF_CC_MPS <- c("MP1", "MP7", "MP9")
SCREF_STATE_GROUPS <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15")
)

SCREF_PRIMARY_STATE_ORDER <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating"
)

SCREF_FINAL_EXTRA_STATE_ORDER <- c("3CA_EMT_and_Protein_maturation")
SCREF_TECHNICAL_STATE_ORDER <- c("Unresolved", "Hybrid")

SCREF_STATE_COLOURS <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
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
