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

`analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R` is the current TCGA-ESCA MP survival workflow. Historical reg/noreg and QC-filtered nMP=19 scripts carry `legacy_` filenames. Current survival workflows use the final 17 merged refined gene lists and five centred state groups. Models may compare:

- continuous score
- median split
- q1 versus q4

All model tables must record cohort definition, event/time fields, score scaling, covariates, split rule, hazard ratio, confidence interval, raw p-value, and multiplicity adjustment. State gene sets are unions of the current refined MP genes in each state; cell-cycle MPs do not define states.

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
