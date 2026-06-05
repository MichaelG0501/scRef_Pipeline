####################
# TCGA CNA Recurrent Event Validation Methodology
####################

## Purpose

`analysis/TCGA/tcga_cna_recurrent_event_validation.R` validates whether chromosome-arm recurrent CNA events found in the scRef malignant subclone workflow show matching MP/state expression associations in TCGA EAC bulk RNA-seq. It also discovers recurrent arm events directly in TCGA segment data and checks whether significant TCGA event-expression associations correspond to scRef recurrent events.

## Inputs

- TCGA GDC segment file: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg`
- TCGA reconstructed RNA-seq metadata: `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds`
- TCGA reconstructed TPM matrix or existing GSVA cache:
  - `ref_outs/TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores.rds`
  - fallback: `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds`
- scRef recurrent event summary: `ref_outs/Auto_cna_subclone_expression/tables/Auto_v2_recomputed_recurrent_cna_event_summary.csv`
- scRef recurrent event feature tests: `ref_outs/Auto_cna_subclone_expression/tables/Auto_v2_recurrent_cna_event_feature_tests.csv`
- CNA annotation workbooks:
  - `ref_outs/OAC_CNV.xlsx`
  - `ref_outs/41588_2018_331_MOESM3_ESM.xlsx`

## Workbook Interpretation

`OAC_CNV.xlsx` is treated as a curated OAC recurrent CNA table. It has separate gain and loss sections with ranked cytobands, genes, approximate frequencies, pathways, and clinical relevance. These rows form the primary external reference set for TCGA arm-call threshold selection.

`41588_2018_331_MOESM3_ESM.xlsx` is parsed by sheet role:

- ST1/ST2 are OCCAMS GISTIC peak sheets. They are transposed peak tables with cytoband, q values, wide peak boundaries, and genes in the peak.
- ST3/ST4 are deletion driver gene sheets.
- ST5/ST6 are amplification driver gene sheets.

The driver sheets are parsed by locating the real `Gene`/`hgnc_symbol` header row rather than assuming a fixed skip count. Gene symbols are normalized for matching while preserving raw symbols in output. This keeps ST5 row 25, `ENSG00000136997 / MYC*`, as a high-confidence `MYC` amplification driver on `chr8q`.

## TCGA Arm-Level CNA Calls

Segment means are converted to arm-level values by length-weighted overlap with hg38 chromosome-arm intervals. Segments crossing a centromere contribute separately to p and q arms. TCGA segment IDs such as `TCGA-2H-A9GF-01` are matched to RNA-seq sample barcodes such as `TCGA-2H-A9GF-01A` using the first 15 TCGA barcode characters.

The main cohort is TCGA EAC primary tumors with both segment CNA and RNA-seq MP/state scores.

## Threshold Selection

The script scans absolute arm mean thresholds from 0.05 to 0.30. For each threshold it calls arm gains/losses, marks recurrent events at at least 15% of TCGA EAC primary samples and at least 3 samples, then compares the recurrent event set with the curated OAC events from `OAC_CNV.xlsx`.

The selected threshold is the threshold with the highest F1 score against the curated OAC event set, then highest Jaccard index, then closest distance to 0.10 as a tie-breaker. Full threshold sensitivity is written so the choice remains auditable.

## Expression Association Tests

TCGA MP and state scores are reused from the gender-validation GSVA cache when present. If the cache is absent, the script reconstructs the same silhouette-filtered MP gene sets and final state gene sets, then runs GSVA on `log2(TPM + 1)`.

For each event-feature pair:

- patients are split into event-present and event-absent groups;
- feature values are TCGA bulk GSVA scores;
- Wilcoxon rank-sum tests compare event-present versus event-absent samples;
- Cliff's delta and median/mean deltas provide effect size;
- BH FDR is computed within feature type and globally.

Two event universes are tested:

1. The top 8 scRef recurrent events from `Auto_v2_recomputed_recurrent_cna_event_summary.csv`, ranked as in `cna_subclone_expression_correlation.R` by recurrence status, sample fraction, and subclone fraction. The dotplot shows these 8 events; the boxplot uses the top 4, matching the scRef recurrent-event plotting convention.
2. TCGA-discovered recurrent events at the selected TCGA threshold.

For TCGA-discovered events, the output includes whether the same event was recurrent in scRef and the scRef sample recurrence fraction.

The validation figures follow the presentation styling of `analysis/cnv/cna_subclone_expression_correlation.R`: MP and state pages are separated, MP/state order is fixed to the scRef plotting order, dotplots use median standardized event delta and `-log10(FDR)`, and boxplots show standardized TCGA GSVA scores.

The rectangular MP/CNA validation layer adds two direct scRef-vs-TCGA comparisons. For all chromosome arms, the script uses the continuous Spearman approach from `analysis/cnv/Auto_mp18_cna_investigation.R`: each arm-level CNA mean is correlated with each MP score, separately in scRef subclones and TCGA bulk samples, with BH FDR correction across the full MP-by-arm family in each dataset. For recurrent events, the script keeps the existing event-positive versus event-negative Wilcoxon framework from `analysis/cnv/cna_subclone_expression_correlation.R` and plots the standardized event delta for every recurrent event/MP pair. Rectangular heatmaps use blue-red direction coloring, asterisks for dataset-level FDR significance, gold outlines for same-trend scRef/TCGA concordance, and magenta plus-marked outlines for same-trend pairs significant in both datasets. The all-arm layer writes both a significant-only plot and an all-trends plot.

## Outputs

Outputs are written under `ref_outs/TCGA/cna_recurrent_event_validation/`.

### Tables

- `Auto_oac_cnv_curated_reparsed.csv`
- `Auto_occams_gistic_peaks_reparsed.csv`
- `Auto_occams_cna_driver_genes_reparsed.csv`
- `Auto_cnv_annotation_long_reparsed.csv`
- `Auto_cnv_event_annotation_summary_reparsed.csv`
- `Auto_tcga_arm_cna_long.csv`
- `Auto_tcga_arm_cna_calls_selected_threshold.csv`
- `Auto_tcga_cna_event_threshold_sensitivity.csv`
- `Auto_tcga_cna_threshold_optimization.csv`
- `Auto_tcga_recurrent_cna_events.csv`
- `Auto_scRef_recurrent_events_tcga_feature_values.csv`
- `Auto_scRef_recurrent_events_tcga_feature_tests.csv`
- `Auto_tcga_discovered_recurrent_event_feature_values.csv`
- `Auto_tcga_discovered_recurrent_event_feature_tests.csv`
- `Auto_tcga_discovered_significant_event_overlap_scRef.csv`
- `Auto_scRef_tcga_event_feature_concordance.csv`
- `Auto_scRef_all_arm_mp_spearman_tests.csv`
- `Auto_tcga_all_arm_mp_spearman_tests.csv`
- `Auto_scRef_tcga_all_arm_mp_spearman_concordance.csv`
- `Auto_scRef_tcga_recurrent_event_mp_concordance_rectangles.csv`
- `Auto_tcga_cna_rectangle_heatmap_summary.csv`
- `Auto_tcga_cna_validation_conclusion.csv`
- `Auto_tcga_cna_recurrent_event_validation_summary.csv`

### Figures

- `Auto_tcga_cna_threshold_optimization.pdf`
- `Auto_tcga_cna_event_association_dotplots.pdf`
- `Auto_tcga_cna_event_boxplots.pdf`
- `Auto_scRef_tcga_all_arm_mp_spearman_rectangles.pdf`
- `Auto_scRef_tcga_all_arm_mp_spearman_rectangles_all_trends.pdf`
- `Auto_scRef_tcga_recurrent_event_mp_rectangles.pdf`

### Logs and Summaries

- `logs/Auto_tcga_cna_recurrent_event_validation_run_summary.rds`
- `logs/Auto_tcga_cna_recurrent_event_validation_run_summary.txt`
- `updates/new_updates/summaries/Auto_tcga_cna_recurrent_event_validation_summary.csv`

## Run Command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/TCGA/tcga_cna_recurrent_event_validation.R
```
