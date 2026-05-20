# Auto CNA Subclone Expression Visuals V2 Methodology

## Purpose

This v2 analysis is a presentation-focused revision of the CNA subclone versus expression/state workflow. It keeps the original v1 discovery outputs intact and reads `ref_outs/Auto_cna_subclone_expression/rds/Auto_cna_subclone_expression_results.rds`. The older `analysis/cnv/cna_subclone_expression_correlation.R` script is therefore not a delete-candidate: it remains the upstream builder for the v1 result RDS that this v2 script consumes.

## Changes From V1

- The consensus CNA heatmap now hides subclone row names and suppresses row-split cluster-name labels on the left. Cluster identity is shown only as a left annotation colour bar.
- Recurrent CNA event statistics are recomputed directly from the event chromosome arm and direction, avoiding the v1 overcounting issue in event-positive/event-negative subclone counts.
- Recurrent event visualisations now use the recurrent events plus the next two ranked CNA events by sample/subclone recurrence. The current run plots six events.
- Recurrent event visualisations now cover:
  - all filtered metaprogrammes,
  - the six real state categories: Classic Proliferative, Basal to Intestinal Metaplasia, SMG-like Metaplasia, Stress-adaptive, Immune Infiltrating, and 3CA EMT/Protein maturation,
  - QC/CNA metrics including `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `cc_score`, `cs_score`, CNA signal/correlation, CNA burden, and gained/lost arm counts.
- Hybrid and Unresolved are excluded from state association plots as requested.
- Significance is encoded with both point size (`-log10(FDR)`) and stars (`*`, `**`, `***`) using BH FDR within each feature group.
- Recurrent-event significance now uses the cohort-level Wilcoxon comparison of event-positive subclones against event-negative subclones. This is appropriate as a screening association for MP UCell means and state proportions because they are computed with the same scoring/state definitions across all subclones; it is not a covariate-adjusted mixed model and can still reflect study/sample composition.
- Recurrent-event distributions are shown as event-positive versus event-negative boxplots using standardized feature values, so MPs, states, and QC/CNA metrics can be displayed without mixing incompatible native axes.
- The chr8q/MYC plot now compares three CNA groups: 8q loss, no 8q CNA, and 8q gain.
- Largest-subclone plots use boxplots of standardized deltas, not a single median point or native-unit deltas on a shared axis, so RNA counts, states, MPs, and QC metrics are not visually forced onto an inappropriate scale.
- Pairwise CNA distance plots now include all MPs, the six states, and QC/CNA metrics.

## Outputs

All outputs are written under `ref_outs/Auto_cna_subclone_expression_v2/`.

- `figures/Auto_v2_cna_consensus_heatmap_no_row_labels.pdf`
- `figures/Auto_v2_recurrent_cna_event_associations_all_features.pdf`
- `figures/Auto_v2_recurrent_cna_event_boxplots_all_features.pdf`
- `figures/Auto_v2_recurrent_cna_event_per_sample_deltas.pdf`
- `figures/Auto_v2_gain_chr8q_myc_mp_per_sample.pdf`
- `figures/Auto_v2_largest_subclone_effects_all_features.pdf`
- `figures/Auto_v2_pairwise_cna_distance_all_features.pdf`
- `tables/Auto_v2_recurrent_cna_event_subclone_presence.csv`
- `tables/Auto_v2_recurrent_cna_event_feature_values.csv`
- `tables/Auto_v2_recurrent_cna_event_feature_tests.csv`
- `tables/Auto_v2_recurrent_cna_event_per_sample_feature_deltas.csv`
- `tables/Auto_v2_largest_subclone_feature_tests.csv`
- `tables/Auto_v2_pairwise_cna_distance_all_feature_tests.csv`
- `tables/Auto_v2_run_summary.csv`
- `rds/Auto_v2_visualisation_results.rds`
- `updates/new_updates/summaries/Auto_cna_subclone_expression_v2_summary.csv`

## Run Command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cnv/Auto_cna_subclone_expression_visuals_v2.R
```

## Validation

The v2 event-presence table is generated from one row per subclone and event after filtering to the event arm. The May 18 rerun analysed 68 samples and 135 subclones, plotted six CNA events, and gave `gain_chr8q` as 22 event-positive subclones and 113 event-negative subclones, with 12 same-sample paired comparisons retained only as a secondary sensitivity table.
