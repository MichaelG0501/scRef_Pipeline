# Cell-State Workflow Methodology

This document describes shared methodology for current centred-state downstream scripts under `analysis/cell_states/`. State construction itself is documented in `analysis/methodology/metaprograms/centred_refined_state_definition_noreg_methodology.md`.

## Current State Decision

The selected state definition is the centred refined, noreg-only Approach B implementation. The canonical inputs are:

- `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds`
- `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds`
- `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_group_max.rds`

The five primary states are `Classic proliferation`, `Basal to intestinal metaplasia`, `SMG to intestinal metaplasia`, `Stress adaptive`, and `Cancer-cell immune mimicry`. Technical labels are `Unresolved` and `Hybrid`. No current downstream workflow should fall back to `Auto_final_states.rds`, `Auto_topmp_v2_*`, or legacy nMP=19 UCell objects.

## Terminal State Figures

`state_mp_sample_abundance.R`, `top_diverse_sample_state_umap.R`, `basal_smg_marker_mp_dissection_heatmap.R`, and `centred_refined_state_marker_discovery.R` read the current centred state vector and current 17-MP grouping. Abundance denominators are all aligned epithelial cells for the relevant sample; MP labels use the descriptions in `centred_refined_mp_state_grouping.csv`.

## Pseudotime And Hybrid Analyses

Pseudotime scripts use Monocle3 on sample-level subsets. Unless a script states a different biological question, the basal state is the root and only the five primary states enter trajectory fitting. Per-sample pseudotime vectors, fitted objects, projection tables, and distance matrices are persistent live intermediates so plotting changes do not require refitting Monocle. Sample eligibility thresholds are script-specific and must be stated in the header or a dedicated methodology file.

Hybrid pair labels are reconstructed from the two highest current group scores among cells labelled `Hybrid`. The node layout compares distance matrices produced by `pseudotime_state_distance_matrix.R`; MDS layout selection and fit diagnostics are descriptive and do not imply lineage direction.

## Legacy Scripts

The `legacy_*` scripts document the uncentred nMP=19, reg/noreg sensitivity, 3CA unresolved-relabel, manual, clustered, or old top-MP alternatives. They are retained for provenance only and must not produce inputs consumed by current scripts.
