# Metaprogram Methodology

This document covers metaprogram scoring and comparison scripts under `analysis/metaprograms/`.

## nMP19 UCell Scoring

`mp_ucell_scoring.R` scores the nMP19 geneNMF metaprograms on `EAC_Ref_epi.rds`.

Implemented rules:

1. Load `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`.
2. Remove MPs with silhouette below zero.
3. Score retained MP gene lists with UCell.
4. Save filtered scores to `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`.
5. Produce self-similarity and heatmap diagnostics.

This filtered UCell matrix is required by the current state definition, final state relabeling, MP correlation plots, clinical association scripts, and SCENIC workflow.

## 3CA MP UCell Scoring

`mp_3ca_ucell_scoring.R` regenerates 3CA pan-cancer MP scores.

Inputs:

- `ref_outs/EAC_Ref_epi.rds`
- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`

Outputs:

- `ref_outs/UCell_3CA_MPs.rds`
- `ref_outs/Auto_3ca_mp_ucell_summary.csv`

The old `delete_mp_external_scoring_superseded.R` script should not be used because it contains a partially commented older workflow and is superseded by the current scorer.

## MP Correlation And Cross-Dataset Comparisons

`final_mp_correlation.R` is the preferred final MP correlation/Jaccard figure script. It uses the current state-associated MP order and excludes cell-cycle MPs unless the figure explicitly needs them.

`mp_pdo_sc_crossdataset_correlation.R` compares PDO and scRef metaprograms. It is terminal and should not write state labels for downstream use.

`legacy_mp_correlation_sc_kmeans_state_temp.R` and `legacy_mp_correlation_pdo_kmeans_state_temp.R` are retained only to document earlier k-means state-like labels. New scripts must not depend on `state_temp.rds`.
