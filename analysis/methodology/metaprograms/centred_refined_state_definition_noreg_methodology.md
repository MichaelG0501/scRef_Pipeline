# Centred Refined State Definition Noreg Methodology

This script defines per-cell states from the centred refined metaprogramme scores. It is a noreg-only analogue of `analysis/cell_states/state_definition_approach_b_reg_noreg.R`: centred refined MP UCell scores are sample-centred and study-scaled directly, with no cell-cycle regression branch and no cell-cycle MPs used for state definition.

The input score matrix is `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`. The cell-cycle MPs `MP1`, `MP5`, and `MP13+` are retained as plotted heatmap rows only. The state-defining groups are:

- Classic proliferation: `MP2+`
- Basal to intestinal metaplasia: `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, `MP10+`
- SMG to intestinal metaplasia: `MP8+`, `MP8b`, `MP16`, `MP18b`, `MP17`, `MP2x`
- Stress adaptive: `MP12`
- Cancer-cell immune mimicry: `MP15`

`MP11c` and `MP18a` are also retained as plotted heatmap rows in an `Excluded` group but are not used to assign states. For each cell, group scores are the maximum adjusted score across MPs in each non-cell-cycle state group. Cells with maximum group score below `0.5` are labelled `Unresolved`; non-unresolved cells with a top-minus-second group score gap below `0.3` are labelled `Hybrid`.

Outputs are written under `ref_outs/Metaprogrammes_Results/centred/state_definition/` using `intermediate/`, `tables/`, `figures/`, and `logs/` tiers. The script also writes a compact summary to `updates/new_updates/summaries/centred_refined_noreg_state_definition_summary.csv`.

Figures mirror the existing Approach B and sample-abundance workflows: a sampled per-cell adjusted MP heatmap with state/CNA/cell-cycle/study annotations, study and overall state proportions with pie chart, cell-cycle score boxplot by state, and per-sample state abundance pages sorted by state diversity and by study.
