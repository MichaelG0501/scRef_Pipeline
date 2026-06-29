####################
# Methodology: refined MP ordered NMF heatmap
####################

`analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R` is a plot-only
terminal script for visualising the original GeneNMF programme similarity
matrix after finalized MP refinement.

The script reads the nMP19 GeneNMF metaprogram object and the finalized merged
refinement assignment table from
`ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds`.
That assignment table maps every original NMF programme to its finalized merged
refined MP after the two-step refinement workflow:

1. `analysis/metaprograms/mp_refinement_submp.R` splits intermediate-silhouette
   parent MPs into sub-MPs using hierarchical clustering of the programme
   similarity matrix.
2. `analysis/metaprograms/mp_refinement_merge_correlated_submps.R` merges
   correlated sub-MPs from the same parent and writes the finalized
   `merged_refined_mp` assignments.

Programs are ordered by the finalized refined MP order used in
`mp_refinement_merge_correlated_submps.R`:
`MP7j, MP9, MP1, MP2+, MP17, MP8+, MP10+, MP14, MP5+, MP7r, MP7v, MP10e,
MP16+, MP18, MP8c, MP15c, MP12c, MP2v, MP8e, MP12a, MP13, MP7+, MP7h, MP8b,
MP12b, MP15a, MP15b`.

Within each finalized refined MP block, programmes retain their original
GeneNMF dendrogram order. The heatmap therefore preserves local programme
structure while enforcing the finalized biological MP order globally. Dotted
diagonal borders and callout labels are drawn from contiguous runs of the
finalized `merged_refined_mp` labels, so the borders represent the final refined
MPs rather than original nMP19 clusters.

The script fails if a finalized refined MP is absent from the strict order. This
prevents future refinement changes from being silently appended in an arbitrary
position.

Outputs:

- `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.pdf`
- `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.png`
- `ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_nmf_ordered_blocks.csv`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_nmf_ordered_similarity.rds`
- `updates/new_updates/summaries/refined_mp_nmf_ordered_heatmap_summary.csv`

Run command:

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R
```
