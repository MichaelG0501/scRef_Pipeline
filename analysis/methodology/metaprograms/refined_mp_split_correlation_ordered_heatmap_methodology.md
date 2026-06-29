####################
# Methodology: split MP correlation ordered heatmap
####################

`analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R` creates a
program-resolution correlation heatmap aligned to
`analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R`.

The plot uses the same finalized refined MP order as the ordered NMF heatmap:
`MP7j, MP9, MP1, MP2+, MP17, MP8+, MP10+, MP14, MP5+, MP7r, MP7v, MP10e,
MP16+, MP18, MP8c, MP15c, MP12c, MP2v, MP8e, MP12a, MP13, MP7+, MP7h, MP8b,
MP12b, MP15a, MP15b`.

Rows and columns are expanded to the original NMF-program resolution. This keeps
the width and height of every finalized MP block identical to the ordered NMF
heatmap. Within each finalized MP block, programs are grouped by their pre-merge
`refined_submp` label from
`ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds`.
For example, single/full MPs such as `MP7j` and `MP9` occupy one internal split
block, while a merged MP such as `MP2+` is internally subdivided into all of its
constituent pre-merge `MP2*` sub-MPs.

The heatmap colour for each program-resolution pixel is the Fisher-Z averaged
per-sample Spearman correlation between the corresponding split/full MP UCell
features. Correlations are computed from
`ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_ucell_scores.rds`
using sample labels from `ref_outs/EAC_Ref_epi.rds`, then cached as
`refined_mp_split_display_correlation_matrices.rds`.

Finalized MP blocks are drawn with thick dotted diagonal borders matching the
ordered NMF heatmap. Internal pre-merge split-MP blocks are drawn with thin grey
solid borders. The external labels and deliberate label gaps follow the ordered
NMF heatmap style, including extra separation before `MP17`, `MP7r`, `MP8c`,
`MP8b`, `MP12b`, `MP15a`, and `MP15b`.

Outputs:

- `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.pdf`
- `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.png`
- `ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_split_correlation_ordered_final_blocks.csv`
- `ref_outs/Metaprogrammes_Results/mp_refinement/tables/refined_mp_split_correlation_ordered_sub_blocks.csv`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_split_correlation_ordered_matrix.rds`
- `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/refined_mp_split_display_correlation_matrices.rds`
- `updates/new_updates/summaries/refined_mp_split_correlation_ordered_heatmap_summary.csv`

Run command:

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R
```
