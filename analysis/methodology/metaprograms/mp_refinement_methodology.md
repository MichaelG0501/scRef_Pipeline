# MP Refinement And Sub-MP Splitting Methodology

## Script

`analysis/metaprograms/mp_refinement_submp.R`

## Purpose

Refine the selected nMP19 scATLAS metaprogrammes by keeping high-silhouette MPs intact, removing negative-silhouette MPs, and splitting low-positive-silhouette MPs into sub-MPs for exploratory interpretation and terminal heatmaps.

This workflow is optional and does not replace the canonical `Metaprogrammes_Results/UCell_nMP19_filtered.rds` scoring unless explicitly adopted downstream.

## Inputs

- `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `ref_outs/geneNMF_outs.rds`
- `ref_outs/EAC_Ref_epi.rds`
- Optional cache: `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`

## Triage Rules

- Keep MPs intact when `metaprograms.metrics$silhouette >= 0.2`.
- Remove MPs when silhouette is `< 0`, matching `mp_ucell_scoring.R`.
- Split MPs when silhouette is `> 0` and `< 0.2`.
- Final state-associated plots exclude cell-cycle MPs and follow the same state grouping/order as `final_mp_correlation.R`.

## Split Selection

For each split-eligible MP, the script subsets `programs.similarity` to the NMF programs assigned to that parent MP, converts cosine similarity to distance (`1 - similarity`), and performs Ward hierarchical clustering.

The script incrementally evaluates split for `k = 2..10`. It stops early when the mean silhouette width exceeds `0.1` and uses that `k`. If no `k` exceeds `0.1`, it defaults to the `k` with the highest mean silhouette width.

Sub-MP labels are assigned in dendrogram order as `MP#a`, `MP#b`, etc.

## Gene-List Derivation

Sub-MP gene lists reproduce the GeneNMF consensus approach used by `getMetaPrograms()`:

1. Weight NMF loadings with `specificity.weight = 5`.
2. Use GeneNMF-style `normVector()`, which normalizes by vector sum, not Euclidean/L2 norm.
3. Build a gene-by-program loading table with program names matching the original GeneNMF naming convention.
4. For each sub-MP, average weighted loadings across assigned programs with the same 3-SD outlier trimming used by GeneNMF.
5. Select genes by `weight.explained = 0.5`, filter by per-program presence confidence `> 0.5`, and cap at `max.genes = 200`.

The script saves both gene vectors and normalized GeneNMF-style weights.

## Outputs

All outputs are under `ref_outs/Metaprogrammes_Results/mp_refinement/`.

- `intermediate/split_results.rds`
- `intermediate/refined_mp_genes.rds`
- `intermediate/refined_mp_gene_weights.rds`
- `intermediate/refined_mp_assignments.rds`
- `intermediate/refined_ucell_scores.rds`
- `intermediate/refined_mp_correlation_matrices.rds`
- `intermediate/refined_mp_jaccard_matrices.rds`
- `tables/split_selection_summary.csv`
- `tables/refined_mp_gene_sizes.csv`
- `tables/refined_mp_correlation_mean_rho.csv`
- `tables/refined_mp_jaccard_index.csv`
- `figures/mp_splitting_diagnostics.pdf`
- `figures/refined_mp_correlation_heatmap.pdf`
- `figures/refined_mp_correlation_heatmap_unsupervised.pdf`
- `figures/refined_mp_jaccard_heatmap.pdf`

## Cache Behavior

`SCREF_FORCE_REBUILD=TRUE` rebuilds split, gene-list, and UCell caches. `SCREF_REPLOT_ONLY=TRUE` reuses valid caches and regenerates tables/figures only.

Caches include an internal algorithm signature so older intermediates are rebuilt instead of silently reused after method changes.

## Run Command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/metaprograms/mp_refinement_submp.R
```
