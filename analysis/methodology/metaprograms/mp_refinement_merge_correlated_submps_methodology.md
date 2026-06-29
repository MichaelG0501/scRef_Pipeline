# MP Refinement Correlated Sub-MP Merge Methodology

## Script

`analysis/metaprograms/mp_refinement_merge_correlated_submps.R`

## Purpose

This script is a downstream layer for `analysis/metaprograms/mp_refinement_submp.R`. It keeps the original sub-MP refinement outputs intact, then merges sub-MPs that appear redundant by within-parent UCell correlation.

## Inputs

- `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `ref_outs/geneNMF_outs.rds`
- `ref_outs/EAC_Ref_epi.rds`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/split_results.rds`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_genes.rds`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_gene_weights.rds`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_ucell_scores.rds`

## Merge Rule

The script computes per-sample Spearman correlations between refined sub-MP UCell scores and combines them by Fisher-Z averaging, matching `mp_refinement_submp.R`.

For each split parent MP, each sub-MP is assessed against all other sub-MPs from the same parent. A sub-MP qualifies for parent-level merge when:

- its correlation is `> 0.4` with at least 25% of the other sub-MPs from the same parent.

If at least two sub-MPs from a parent qualify, they are merged into one parent-plus feature named like `MP2+`. Non-qualifying sub-MPs are retained with their original sub-MP names and original gene lists.

## Merged Gene-List Derivation

For each parent-plus feature, all NMF programs assigned to the qualified sub-MPs are pooled. The gene list and weights are then recalculated using the same GeneNMF-style consensus implementation in `mp_refinement_submp.R`:

1. Weight NMF loadings with `specificity.weight = 5`.
2. Normalize vectors by vector sum using the local `normVector()` implementation.
3. Average weighted loadings across pooled NMF programs with 3-SD outlier trimming.
4. Select genes with `weight.explained = 0.5`.
5. Keep genes present in `> 0.5` of the pooled per-program gene lists.
6. Cap the merged signature at 200 genes.

## Outputs

All main outputs are written under `ref_outs/Metaprogrammes_Results/mp_refinement/`.

- `intermediate/merged_refined_mp_genes.rds`
- `intermediate/merged_refined_mp_gene_weights.rds`
- `intermediate/merged_refined_mp_assignments.rds`
- `intermediate/merged_refined_ucell_scores.rds`
- `intermediate/merged_refined_mp_correlation_matrices.rds`
- `tables/merged_refined_mp_merge_decisions.csv`
- `tables/merged_refined_mp_gene_sizes.csv`
- `tables/merged_refined_mp_correlation_mean_rho.csv`
- `figures/refined_mp_correlation_heatmap_unsupervised_merged.pdf`
- `updates/new_updates/summaries/mp_refinement_merge_correlated_submps_summary.csv`

## Cache Behavior

`SCREF_FORCE_REBUILD=TRUE` rebuilds merged gene and UCell caches. `SCREF_REPLOT_ONLY=TRUE` reuses valid merged caches and regenerates tables and the final heatmap.

## Run Command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/metaprograms/mp_refinement_merge_correlated_submps.R
```
