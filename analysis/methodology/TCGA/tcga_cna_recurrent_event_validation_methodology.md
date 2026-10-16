####################
# TCGA CNA Recurrent Event Validation Methodology
####################

## Purpose

`analysis/TCGA/tcga_cna_recurrent_event_validation.R` validates whether chromosome-arm recurrent CNA events found in the scRef malignant subclone workflow show matching MP/state expression associations in TCGA EAC bulk RNA-seq. It also discovers recurrent arm events directly in TCGA segment data and checks whether significant TCGA event-expression associations correspond to scRef recurrent events.

## Inputs

- TCGA GDC segment file: `ref_outs/TCGA/cna_recurrent_event_validation/intermediate/source_data/esca_tcga_gdc_segments.seg` (initial import: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg`)
- TCGA reconstructed RNA-seq metadata: `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds`
- Required TCGA expression provenance and shared GSVA cache:
  - `ref_outs/TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores_centred17.rds`
  - `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds`
- scRef recurrent event summary: `ref_outs/Auto_cna_subclone_expression/threshold_1/tables/Auto_v2_recomputed_recurrent_cna_event_summary.csv`
- scRef recurrent event feature tests: `ref_outs/Auto_cna_subclone_expression/threshold_1/tables/Auto_v2_recurrent_cna_event_feature_tests.csv`
- CNA annotation workbooks:
  - `analysis/TCGA/OAC_CNV.xlsx`
  - `analysis/TCGA/41588_2018_331_MOESM3_ESM.xlsx`

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

TCGA MP and state scores are consumed only from the versioned final-17 gender-validation GSVA cache. The canonical TPM matrix, MP genes, grouping CSV and shared configuration are fingerprinted and checked against its producer contract. Missing, stale or unversioned scores stop; the gender workflow must be rebuilt first. This script does not calculate a substitute GSVA cache or read a mixture TXT as an alternative matrix. MP names, current gene-set identities, score dimensions and sample coverage must match. Legacy nMP19, excluded MPs and appended 3CA sets are not accepted.

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
/opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_cna_recurrent_event_validation.sh
```


## Persistent source data and plot-only reruns

The canonical root is `ref_outs/TCGA/cna_recurrent_event_validation/` inside the live project. All persistent source data, intermediate objects, tables, figures and logs remain live. The segment input is imported once into `intermediate/source_data/`; the original external file is retained. Later runs use the project snapshot. `SCREF_FORCE_REBUILD=TRUE` rebuilds this script's arm-call cache, not the independently produced gender GSVA cache or input snapshot. The arm cache fingerprints the segment file and CNA script. Missing required workbooks or scRef result tables stop immediately; comparative panels are not silently skipped.

Additional required input is `ref_outs/Auto_cna_subclone_expression/threshold_1/rds/Auto_cna_subclone_expression_results.rds`. It supplies arm matrices, existing scRef cluster assignments and MP features. The cache's `sample` identifiers must originate from upstream `orig.ident`; the TCGA-only 15-character matching rule is never applied to scRef barcodes. Historical state-column spellings in the scRef table are mapped explicitly to current `Squamous-to-intestinal`, `Glandular-to-intestinal` and `Stress-adaptive` labels for concordance joins; no old state-vector object is loaded.

The final source bundle is `intermediate/Auto_tcga_cna_plot_source_data.rds`. It retains TCGA sample scores, TCGA/scRef arm matrices, cluster assignments, chr8q/MYC plot values, feature metadata, event test/value tables and the selected threshold. TCGA clusters are also exported as `tables/Auto_tcga_cna_cluster_assignments.csv`. The complete figure bundle is `intermediate/Auto_tcga_cna_figure_bundle.rds`; it contains rendered ggplot grobs and the native paired ComplexHeatmap object for all successfully produced PDFs, including the combined heatmap and chr8q/MYC comparison. This preserves plotting inputs that previously existed only in memory.

`SCREF_REPLOT_ONLY=TRUE` renders that bundle without reading raw segments, loading scRef analytical caches, recomputing GSVA/tests or reclustering. It also refreshes the compact summary copy and writes `logs/Auto_tcga_cna_replot_run_summary.rds/.txt`. The first replot requires a successful updated full run. Plot-style changes require a full run to regenerate the grobs. Force rebuild and replot-only are mutually exclusive. Both plotting workflows enforce a compact summary smaller than 100 KB. Missing replot artifacts stop immediately.

The wrapper requests 2 CPUs, 8 GB and 4 hours in `dmtcp`, with explicit name and `#PBS -koed`. The R entry point requires a PBS allocation. Run reconstruction, then gender validation, then CNA validation; verify each producer's expected outputs before starting its consumer.

## Clustering, filtering and limitations

The paired heatmap uses the existing scRef CNA cluster labels and Ward.D2 row order. TCGA clustering uses Euclidean arm-mean distances with Ward.D2 linkage. Candidate cluster counts are 2 through min(6, n-1); solutions containing singleton clusters are excluded. The highest mean silhouette wins. Fewer than four samples, a silhouette error or no finite admissible solution stops; there is no arbitrary replacement cluster count. The plot is a hierarchical-cluster comparison, despite the historical `consensus` filename; it does not perform resampled consensus clustering.

The chr8q comparison preserves the existing scRef illustrative call threshold of absolute arm mean 0.05 and uses the selected TCGA threshold for bulk calls. Available group pairs use Wilcoxon tests and nominal significance symbols. Missing categories do not generate fabricated group comparisons. The subplot's scRef unit is a subclone and the TCGA unit is a sample, so scRef tests do not account for correlation among subclones from the same sample. The main recurrent-event tests require at least three event-positive and three event-negative observations. Continuous Spearman tests require at least five finite pairs and variation in both quantities. Undefined tests remain missing; unexpected test failures stop.

Threshold optimization uses the EAC primary cohort with segments before restricting event tests to RNA-score overlap. Calibration against curated OAC recurrence is descriptive and not an independently calibrated error-rate estimate. The code retains historical neutral-zero values for uncovered arm intervals and non-finite cached arm values; coverage columns must therefore be inspected before interpreting arm-level results. The hg38 coordinate assumption must be verified against segment provenance before adopting a new segment snapshot. The discovery lookup reports overlap with the selected top-eight scRef event universe, not every possible scRef recurrent event. Association or absence of statistical significance does not establish a causal relation or biological independence.

## Validation and adoption

R and shell syntax, all required `dmtcp` package imports, the Gaussian GSVA producer interface, persistent-path enforcement, cache-staleness detection, 17-MP/group validation and pure figure replot are checked using small PBS fixtures under `ref_outs/TCGA/contract_validation/<job_id>/`. Full biological reruns require restored canonical annotation workbooks, valid upstream scores, readable scRef result caches and verified registered outputs. A zero PBS exit status alone is insufficient. Existing historical outputs and scripts are retained; updating this implementation does not make old results newly validated.


The October 2026 contract validation completed as PBS jobs `4373148.pbs-7` and `4373174.pbs-7`, both state `F`, `Exit_status = 0`. The final run also exercised the actual STAR parser and TPM reconstruction on synthetic inputs, including rejection of alternate TPM schemas and negative TPM values. Expected fixture PDF, expression cache, replot logs and result report were verified under `ref_outs/TCGA/contract_validation/4373174.pbs-7/`. The installed GSVA version was 1.46.0. The test harness is retained as `temp/Auto_tcga_contract_validation.sh`; it does not execute production analyses or replace production outputs.


The restored annotation workbooks are read directly from `analysis/TCGA/`, as supplied by the user. The production rerun reuses the existing canonical reconstructed TPM/metadata and rebuilds gender GSVA before CNA validation. The production validation harness is `temp/Auto_tcga_production_validation.sh`; it verifies registered output files, sample/feature coverage and saved-grob replot, then writes `reports/Auto_output_validation.csv`, `reports/Auto_figure_page_counts.csv`, PDF metadata and PNG previews under each workflow root.


Rendered ggplot grobs are measured on the final Cairo PDF device at the intended page dimensions before caching. This avoids font-metric differences from an implicit default device. Composite chr8q plots are constructed on that device as well. The paired heatmap stores the native ComplexHeatmap object, retaining its fixed row orders and disabled clustering, and renders directly on a 20-by-10-inch Cairo page with 12 mm padding. This avoids off-screen capture drift. Dataset-specific scRef and TCGA cluster legends retain their independent palettes. Plot-only reruns perform layout rendering only; they do not recluster. Rectangular heatmaps omit the redundant y-axis title, split the Spearman legend title over two lines and retain explicit margins. Long event-association and event-boxplot titles wrap at 70 characters while preserving their font size. Production PDF validation checks every text span against page bounds in addition to file/page counts; numerical analyses are unchanged by these layout fixes.


Production regeneration on 15 October 2026 completed with gender job `4373277.pbs-7`, CNA job `4373346.pbs-7` and validation job `4373347.pbs-7`, all state `F` and exit status zero. The existing reconstructed TPM and metadata were reused. Validation verified 53 registered nonempty files, nine PDFs across 15 pages, text within page bounds and standalone cached replot. Separate CSV checks verified unique observational units and 87 TCGA samples per CNA event-feature comparison. Exact outputs and run provenance are recorded in `ref_outs/TCGA/cna_recurrent_event_validation/reports/Auto_production_run_summary.txt`; observation checks are in `reports/Auto_observational_unit_validation.csv` under the same CNA root. The canonical CNA arm threshold is 0.12 for this run.
