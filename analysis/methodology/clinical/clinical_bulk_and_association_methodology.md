# Clinical, Bulk, And Association Methodology

This document covers clinical association plotting and TCGA/GEO bulk survival workflows under `analysis/clinical/`.

## Clinical Association Plotting

Final clinical association plots should use:

- `ref_outs/meta_full_epi.rds`
- `ref_outs/Auto_final_states.rds`
- `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- `ref_outs/UCell_3CA_MPs.rds` when retained 3CA states are shown
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx`, sheet 3

`clinical_association_final_boxplots.R` renders sample-level MP and state boxplots and writes statistical summaries. `clinical_association_final_stacked.R` renders stacked state composition plots and is currently untracked. Future cleanup should merge both plotting modes into one tracked script with separate output sections for boxplots and stacked bars.

Clinical variables are normalized before plotting:

- `orig.ident` is reconstructed from `Author`, `Year`, and `Sample Name`.
- `Age_Group` is `>60` versus `<=60`.
- `Gender` is converted to `Female`/`Male`.
- `Treatment` is ordered as `Tx-naive`/`Post` while preserving the repository's existing spelling in source data.
- `Clinical response` is converted from `R` to `Responder` and all other observed coded responses to `Nonresponder`.

## TCGA Survival

`tcga_mp_state_survival_reg_noreg.R` is the reference TCGA script. It compares MP/state gene-set scoring strategies across split methods:

- continuous score
- median split
- q1 versus q4

It uses MP gene sets from the nMP19 object, DGE-derived state gene sets where available, and state labels from the current final or noreg state vectors. New TCGA scripts should default to noreg/final labels and avoid `state_temp.rds`.

`tcga_mp_state_survival_qc_filtered.R` repeats this workflow on TCGA samples retained by cross-platform QC.

## GEO And Cross-Platform Survival

`geo_survival_data_prep.R` downloads and prepares GEO cohorts. Downloads are cached under `ref_outs/geo_survival/raw/`; platform probes are collapsed to gene symbols using a highest-variance probe-per-gene rule.

`bulk_tcga_geo_qc.R` harmonizes TCGA RNA-seq and GEO microarray expression:

1. Intersect shared genes.
2. Transform TCGA as `log2(TPM + 1)`.
3. Keep supplied GEO processed log-scale matrix.
4. Standardize each dataset by gene-level z-scores.
5. Use PCA, expression-strength metrics, and histology checks to mark samples as retained or removed.

`bulk_tcga_geo_integrated_survival.R` recomputes MP/state scores on the harmonized expression matrix and runs dataset-aware Cox models, including direction and interaction summaries.

## Output And Replot Policy

Survival scripts should save model-ready score matrices and Cox result tables before plotting. Volcano/KM style plots should be reproducible from those tables without recomputing GSVA or gene-set scores.
