# Centred Refined State Definition Noreg Methodology

This document describes `analysis/metaprograms/centred/06_centred_refined_state_definition_noreg.R`, the canonical state-definition workflow. Centred refined MP UCell scores are sample-centred and study-scaled directly, with no cell-cycle regression branch and no cell-cycle MPs used for state definition.

The input score matrix is `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`. The cell-cycle MPs `MP1`, `MP5`, and `MP13+` are retained as plotted heatmap rows only. The state-defining groups are:

- Classic proliferation: `MP2+`
- Basal to intestinal metaplasia: `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, `MP10+`
- SMG to intestinal metaplasia: `MP8+`, `MP8b`, `MP16`, `MP18b`, `MP17`
- Stress adaptive: `MP12`
- Cancer-cell immune mimicry: `MP15`

The final panel contains 17 MPs: three cell-cycle rows and 14 state-defining rows. Step 04 removes `MP2x` and `MP11c` because each is represented in fewer than three samples, and removes `MP18a` as the documented explicit quality/coverage exclusion. They are therefore absent from this script rather than plotted as an excluded group.

For each cell, every MP is centred within `orig.ident`, then scaled by study so study-specific dispersion does not dominate. A state-group score is the maximum adjusted score across MPs in that non-cell-cycle group. Cells with maximum group score below `0.5` are labelled `Unresolved`; among the remainder, cells whose top-minus-second group score gap is below `0.3` are labelled `Hybrid`. These are the fixed Approach B thresholds used by the current workflow.

Outputs are written under `ref_outs/Metaprogrammes_Results/centred/state_definition/` using `intermediate/`, `tables/`, `figures/`, and `logs/` tiers. The script also writes a compact summary to `updates/new_updates/summaries/centred_refined_noreg_state_definition_summary.csv`.

Figures mirror the existing Approach B and sample-abundance workflows: a sampled per-cell adjusted MP heatmap with state/CNA/cell-cycle/study annotations, study and overall state proportions with pie chart, cell-cycle score boxplot by state, and per-sample state abundance pages sorted by state diversity and by study.

The script is a lightweight deterministic rebuild from persistent live inputs. It does not rerun NMF, refinement, or UCell scoring. `meta_full_epi.rds` is optional and adds CNA annotations only; its absence does not change state assignment.
