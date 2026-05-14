# Summary Methodology

`cross_sample_summary.R` summarizes per-sample pipeline outputs from `ref_outs/by_samples/`.

Expected inputs:

- per-sample post-QC RDS files
- annotation outputs where available
- epithelial/malignancy outputs where available

Outputs are terminal summary tables and figures. They should not be used as state-definition inputs.

`legacy_summary_qc_plots.R` is an older summary plotting script retained for reference only.
