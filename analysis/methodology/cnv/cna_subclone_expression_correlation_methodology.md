# CNA Subclone Expression Correlation Methodology

## Purpose

This analysis asks whether CNA subclone structure explains malignant epithelial expression states in the scATLAS. It starts from the completed malignant CNA subclone calls in `ref_outs/Auto_malignant_subclone_mp/` and adds arm-level CNA profiles, clone dominance testing, pairwise CNA-expression distances, consensus CNA pattern discovery, OAC-specific external CNA annotation, and presentation-quality visualisations with corrected event-counting statistics.

This is a single self-contained script that replaces the previous two-script workflow (`cna_subclone_expression_correlation.R` + `Auto_cna_subclone_expression_visuals_v2.R`).

## Inputs

- `Auto_malignant_subclone_cells.csv`: malignant cells with sample, subclone, top MP, final state, CNA signal, RNA QC, `cc_score`, and `cs_score`.
- `Auto_malignant_subclone_summary.csv`: analysed/skipped sample status and subclone count.
- `Auto_malignant_subclone_mp_subclone_tests.csv`: mean MP z-score per sample/subclone/MP.
- `by_samples/<sample>/<sample>_outs.rds`: inferCNA gene-by-cell CNA matrix.
- `OAC_CNV.xlsx`: manually curated OAC CNA gain/loss genes, frequencies, pathways, and clinical relevance.
- `41588_2018_331_MOESM3_ESM.xlsx`: OCCAMS supplementary GISTIC peaks and CNA driver sheets from Frankell et al., Nature Genetics 2019.
- Optional live annotation from cBioPortal TCGA PanCancer Atlas EAC GISTIC calls.

## Arm-Level CNA Profiles

For every analysed sample, the script reloads the inferCNA matrix and restricts it to malignant cells present in the completed subclone table. Genes are mapped to hg38 chromosome arms using `hg38_gencode_v27.txt` and fixed centromere coordinates.

For each sample/subclone/arm:

- `arm_mean` is the mean inferCNA value over all genes on that arm and all cells in the subclone.
- `cna_call = 1` for `arm_mean >= 0.10`, `-1` for `arm_mean <= -0.10`, otherwise `0`.
- CNA burden is summarised as mean absolute arm signal, arm-signal SD, number/fraction of gained arms, and number/fraction of lost arms.

Missing arm values are treated as neutral (`0`) for downstream distances and consensus clustering.

## Dominant Clone Analysis

Within each sample, the largest subclone is compared against the second largest using an exact one-sided binomial test on cell counts. A sample is marked `significant_dominant` only when largest subclone fraction >= 0.50, largest-minus-second gap >= 0.15, and BH-adjusted binomial P < 0.05.

For multi-subclone samples, each expression/state/QC feature is tested as largest clone minus the cell-count-weighted mean of the remaining subclones. Largest-subclone effects use boxplots of standardized deltas so RNA counts, states, MPs, and QC metrics are not visually forced onto an incompatible scale.

## Pairwise Subclone Distance

For every pair of subclones within a sample, the script computes mean absolute arm-CNA distance, Euclidean arm-CNA distance, fraction of arms with discordant CNA calls, MP Euclidean and mean absolute distance, state-proportion L1 and Jensen-Shannon distance, and absolute deltas for key features. Associations are tested with Spearman correlations, sample-centered Spearman, and sample fixed-effect linear models.

## Cross-Sample Consensus CNA Patterns

Subclones are clustered using Ward hierarchical clustering on continuous arm-level CNA profiles. The number of clusters is selected from k=2 to k=6 by average silhouette width. The consensus heatmap hides subclone row names and suppresses row-split cluster-name labels; cluster identity is shown only as a left annotation colour bar.

## Recurrent Event Analysis

Each arm-level gain/loss is summarised by fraction of subclones and samples with the event. An event is marked recurrent at >= 15% of samples and >= 3 samples. Event presence is computed directly from arm and direction, avoiding the v1 overcounting issue. Visualisations cover the recurrent events plus the next ranked CNA events by sample/subclone recurrence.

Recurrent event visualisations cover:
- all filtered metaprogrammes
- the six real state categories (excluding Hybrid and Unresolved)
- QC/CNA metrics including nCount_RNA, nFeature_RNA, percent.mt, cc_score, cs_score, CNA signal/correlation, CNA burden, and gained/lost arm counts

Significance is encoded with both point size (-log10 FDR) and stars using BH FDR within each feature group. The cohort-level Wilcoxon comparison of event-positive versus event-negative subclones is the primary test. The chr8q/MYC plot compares three CNA groups: 8q loss, no 8q CNA, and 8q gain.

## External CNA Annotation

- `OAC_CNV.xlsx` contributes curated OAC genes, broad frequencies, pathways, and therapeutic relevance.
- OCCAMS GISTIC peak sheets are transposed into peak-by-row tables and mapped to chromosome arms.
- OCCAMS high-confidence and candidate driver sheets contribute gene-level CNA driver labels.
- TCGA cBioPortal is queried for gene-level GISTIC frequencies. If unavailable, the script continues with local annotations.

## Outputs

All outputs are written under `ref_outs/Auto_cna_subclone_expression/`.

### Tables
- `Auto_subclone_arm_cna_long.csv`: long arm-level CNA means and calls.
- `Auto_subclone_feature_summary_with_clusters.csv`: per-subclone features.
- `Auto_dominant_clone_feature_tests.csv`: largest-clone versus rest feature tests.
- `Auto_pairwise_cna_expression_distance_tests.csv`: CNA distance vs expression divergence.
- `Auto_recurrent_cna_events_annotated.csv`: recurrent arm events with annotations.
- `Auto_recurrent_cna_event_feature_tests.csv`: recurrent event feature associations.
- `Auto_recurrent_cna_event_subclone_presence.csv`: per-subclone event presence.
- `Auto_recurrent_cna_event_feature_values.csv`: subclone feature values by event.
- `Auto_recurrent_cna_event_per_sample_feature_deltas.csv`: per-sample paired deltas.
- `Auto_largest_subclone_feature_tests.csv`: standardized largest-subclone tests.
- `Auto_pairwise_cna_distance_all_feature_tests.csv`: pairwise CNA-distance tests.
- `Auto_run_summary.csv`

### Figures
- `Auto_cna_consensus_heatmap.pdf`: arm-level CNA consensus heatmap.
- `Auto_recurrent_cna_event_associations_all_features.pdf`: event association dot plots.
- `Auto_recurrent_cna_event_boxplots_all_features.pdf`: event boxplots.
- `Auto_recurrent_cna_event_per_sample_deltas.pdf`: per-sample event deltas.
- `Auto_gain_chr8q_myc_mp_per_sample.pdf`: chr8q CNA vs MYC MP.
- `Auto_largest_subclone_effects_all_features.pdf`: largest-subclone boxplots.
- `Auto_pairwise_cna_distance_all_features.pdf`: pairwise CNA-distance dot plots.
- `Auto_per_sample_heatmap_recurrent_events.pdf`: per-sample cell-level CNA heatmaps.

### RDS
- `Auto_cna_subclone_expression_results.rds`: full intermediate cache.

## Cache and Replot

The script saves a full intermediate RDS after computation. Set `SCREF_REPLOT_ONLY=TRUE` to skip computation and regenerate plots from the cached RDS.

## Run Command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cnv/cna_subclone_expression_correlation.R
```

## References

- cBioPortal REST API: https://docs.cbioportal.org/web-api-and-clients/
- cBioPortal GISTIC values: https://docs.cbioportal.org/file-formats/
- Frankell et al. 2019 OCCAMS EAC study: https://www.nature.com/articles/s41588-018-0331-5
