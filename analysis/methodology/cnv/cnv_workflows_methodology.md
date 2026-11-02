# CNV Workflow Methodology

This document covers CNV and malignant subclone scripts under `analysis/cnv/`.

## Core CNV Scripts

`cnv_profiling.R`, `legacy_cnv_subsetting.R`, and `cnv_plotting.R` form the older core CNV workflow; the reference-spiking subsetting script is legacy and is not an active dependency:

1. infer CNA profiles for epithelial cells by sample
2. subset cells by CNA signal
3. create publication-style CNA heatmaps

These scripts share the package pattern `data.table`, `dplyr`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `Seurat`, and `infercna`.

## Malignant Subclone MP Heatmap

`cnv_malignant_subclone_mp_heatmap.R` is the active cohort-level CNA subclone versus MP/state heterogeneity workflow.

Implemented subclone logic:

1. Load per-sample InferCNA outputs and epithelial RDS objects.
2. Select malignant epithelial cells.
3. Compute CNA signal and top CNA-signal genes.
4. Run Louvain kNN clustering with `k = 15`.
5. Convert gene-level CNA to chromosome-arm summaries.
6. Call arm-level gains/losses at `+/-0.10`.
7. Drop provisional clusters below 20 cells or 5 percent sample fraction.
8. Merge clusters with same-shape strength-scaled profiles.
9. Retain robust subclones when arm-level CNA patterns are distinct.
10. Plot CNA heatmaps, MP/state distributions, and subclone statistics.

No silhouette-based subclone selection is used.

## Untracked CNV Scripts

`cna_subclone_expression_correlation.R`, `cna_subclone_expression_correlation.sh`, and `mp_chromosomal_mapping.R` are currently untracked and should remain unstaged unless the user explicitly requests adding them. Their methodology lives in `analysis/methodology/cnv/cna_subclone_expression_correlation_methodology.md`, also currently unstaged.

####################
## Manual InferCNA subclone MP effect estimability

`cnv_malignant_subclone_mp_heatmap.R` uses the default R-scaled within-sample
MAD of retained-cell MP scores. A nonfinite or nonpositive MAD makes the
standardized median difference unestimable (`NA`) and the significant flag
`FALSE`; epsilon replacement is forbidden. Export `sample_mad` and
`effect_estimable`. Retain raw means, Wilcoxon P values, within-sample BH
correction and the requested absolute-effect cutoff. No gene-detection filter
or arbitrary positive-MAD floor is added to these MP-score comparisons.
`SKIP_CNA_LOOP=TRUE` rejects caches lacking the estimability fields or using
significance flags inconsistent with this policy and the requested threshold.
This source policy applies to future recomputation; previously saved
threshold results must not be assumed to have been recomputed by a source edit.
####################
