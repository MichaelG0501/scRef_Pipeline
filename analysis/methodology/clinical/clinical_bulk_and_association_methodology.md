# Clinical, Bulk, And Association Methodology

This document covers clinical association plotting and TCGA/GEO bulk survival workflows under `analysis/clinical/`.

## Clinical Association Plotting

Final clinical association plots should use:

- `ref_outs/meta_full_epi.rds`
- `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds`
- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`
- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx`, sheet 3

`clinical_association_final_boxplots.R` renders sample-level MP and state boxplots and writes Wilcoxon/Kruskal-Wallis summaries with BH adjustment within clinical variable and feature type. `clinical_association_final_stacked.R` renders current-state composition by clinical group. `clinical_association_mp_ucell_plots.R` summarizes sample-mean current UCell scores before group-level plotting; no per-cell test is used as a substitute for sample replication.

Clinical variables are normalized before plotting:

- `orig.ident` is reconstructed from `Author`, `Year`, and `Sample Name`.
- `Age_Group` is `>60` versus `<=60`.
- `Gender` is converted to `Female`/`Male`.
- `Treatment` is ordered as `Tx-naive`/`Post` while preserving the repository's existing spelling in source data.
- `Clinical response` is converted from `R` to `Responder` and all other observed coded responses to `Nonresponder`.

## TCGA Survival

`analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R` is the current TCGA-ESCA survival workflow. It consumes the reconstructed gene-symbol TPM matrix and reconstructed clinical metadata, requires exact sample-barcode agreement, and restricts analysis to unique primary-tumour (`sample_type_code == 01`) EAC samples. Overall-survival time is the reconstructed cBioPortal OS time in days and the event is death. Rows require positive time and an event code of 0 or 1. The September 2026 audit identified 88 unique evaluable cases and 45 deaths, with no duplicate cases or samples and no discordance between vital status and OS event.

Expression is transformed as `log2(TPM + 1)` before Gaussian-kernel GSVA. The MP input must be the exact current 17 centred refined MPs. State-union gene sets come from `SCREF_STATE_GROUPS`; state-marker sets use the first 20 genes per state from the current ranked marker table. Every requested gene set must retain at least five genes or the workflow stops. No alternate input, histology inference, or legacy MP fallback is permitted.

The prespecified survival models are all univariable Cox proportional-hazards models with no clinical, purity, node, stage, age, or treatment covariates. Continuous scores are standardized within the EAC cohort and reported per one-SD increase. Categorical models use a cohort median split or compare the upper and lower quartiles while excluding the middle half. Results record the model formula, `covariates = none`, scaling, hazard ratio, 95% confidence interval, raw p-value, sample count, event count, and BH adjustment within feature type and split method.

The exploratory optimal-cut analysis tests score quantiles from 20% to 80% in 5% increments and selects the smallest Cox Wald p-value for each feature. Its displayed and saved p-values are unadjusted minima across searched cuts and are not confirmatory p-values. The KM pages display that exact selected Cox Wald p-value rather than recomputing a differently labelled test. KM features are the five lowest-p protective and five lowest-p adverse associations. The optimal-cut volcano is 18 by 8 inches with base font 12, label size 2.8, and title size 14; the separate KM file is 10 by 8 inches.

Persistent outputs include the model-ready data, standard Cox table, optimal-cut table, volcanoes, KM curves, compact summary, and run report. Bulk GSVA can reflect both tumour and non-tumour expression and must not be interpreted as tumour-cell-specific activity. Optimal-cut findings require validation in an independent cohort.

## GEO And Cross-Platform Survival

`geo_survival_data_prep.R` downloads and prepares GEO cohorts. Downloads are cached under `ref_outs/geo_survival/raw/`; platform probes are collapsed to gene symbols using a highest-variance probe-per-gene rule.

`bulk_tcga_geo_qc.R` harmonizes TCGA RNA-seq and GEO microarray expression:

1. Intersect shared genes.
2. Transform TCGA as `log2(TPM + 1)`.
3. Keep supplied GEO processed log-scale matrix.
4. Standardize each dataset by gene-level z-scores.
5. Use PCA, expression-strength metrics, and histology checks to mark samples as retained or removed.

`bulk_tcga_geo_integrated_survival.R` recomputes the current centred MP/state scores on the harmonized expression matrix and runs dataset-aware Cox models, including direction and interaction summaries. `bulk_tcga_geo_meta_survival.R` fits dataset-specific Cox models and combines log hazard ratios by random-effects meta-analysis. `bulk_tcga_geo_feature_presence.R` is a coverage/QC visualization, not an inferential survival result. The GEO-only script is retained for cohort-specific checking but uses the same current signatures.

## Output And Replot Policy

All inputs, model-ready score matrices, Cox/meta-analysis result tables, figures, and compact summaries are written to live. Raw GEO downloads are also retained in live because they are required to reproduce the prepared expression matrices. Volcano/KM style plots should be reproducible from persistent tables without recomputing GSVA or gene-set scores. Network downloads and GSVA/Cox work must run through PBS.

<!-- #################### -->
## OCCAMS Bulk Cohort

The current OCCAMS workflow is implemented by `Auto_occams_bulk_mp_survival.R` and `Auto_occams_bulk_clinical_associations.R`. It uses exactly the 281 subjects passing the documented GRCh37 annotation-compatibility QC flag, the current 17 centred refined MPs, and five non-cell-cycle state unions. Detailed endpoint construction, repeated-row handling, clinical-variable selection, minimum group sizes, tests, cache semantics, and limitations are in `Auto_occams_bulk_clinical_methodology.md`.
<!-- #################### -->
