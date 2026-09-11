# SCENIC Parse regulon-enrichment comparison methodology

## Purpose

`analysis/cell_states/final_mp_scenic_parse_overlap.R` answers a directional biological question: for each scRef final MP or final state, which Parse treatment timepoint has the most similar *selectively enriched TF-regulon signature*? It does not require the MP/state and timepoint labels to represent the same biology. For example, a Stress adaptive state can be closest to an acute-response timepoint such as T1 when their selectively enriched regulons agree.

## Inputs

- `ref_outs/final_mp_scenic/Auto_final_mp_scenic_rss.rds`: scRef regulon RSS by final MP.
- `ref_outs/final_mp_scenic/Auto_final_mp_scenic_state_rss.rds`: scRef regulon RSS by final state.
- Parse balanced timepoint RSS, preferred: `/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/parse_outs/cell_states/timepoint_scenic/intermediate/combined_timepoint_balanced2600_rss.rds`.
- Parse full timepoint RSS fallback: `/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/parse_outs/cell_states/timepoint_scenic/intermediate/combined_timepoint_rss.rds`.

The Parse balanced RSS is preferred because it recomputes timepoint RSS after equal downsampling to 2,600 cells per timepoint. The Parse workflow has one combined SCENIC network for all six timepoints, so its within-Parse comparison is internally consistent.

## Method

1. Regulan names are reduced to a canonical TF name by removing the gene-count suffix and `_extended`. If both direct and extended regulons exist for a TF, the direct regulon is retained. This avoids double weighting one TF and avoids the old script's ambiguous name-to-row mapping.
2. For each canonical TF separately, RSS is z-scored across the labels within its own dataset:
   - across scRef MPs for the MP analysis;
   - across scRef final states for the state analysis;
   - across Parse timepoints for the Parse comparison.
3. Every MP/state × timepoint pair is then assessed over the canonical TFs shared by the two runs. The primary quantity is the cosine similarity of the positive within-run z-score vectors. High values mean the same TF regulons are selectively enriched in both labels relative to the other labels in their respective datasets.
4. Two supporting quantities are displayed alongside the primary score:
   - signed Spearman correlation across all z-scored TF profiles, which captures agreement including relative depletion;
   - weighted Jaccard overlap of the top `top_n` positively enriched TF regulons (default 20), which gives greater weight to the highest-ranked shared TFs.
5. The evidence tables list both the literal shared top TFs and the leading concordant TFs, ranked by the product of the two positive z-scores. These are the regulators to inspect when judging a proposed temporal association.

## Why raw RSS correlation and target overlap are not used as the primary similarity

RSS is a label-specific score whose scale and background distribution depend on the cells, labels, and independently inferred SCENIC network. A raw RSS correlation between scRef and Parse therefore mixes biological similarity with run-specific scale and label-composition effects. Likewise, Parse target sets come from one combined network and do not change from one timepoint to another; target-gene overlap can be useful TF-level context but is not evidence that one timepoint is closer than another. The old raw-RSS/min-max combined score and its target-overlap contribution are intentionally retired.

## Outputs

All main outputs are beneath `ref_outs/final_mp_scenic/parse_overlap/`:

- `tables/mp_timepoint_regulon_enrichment_similarity.csv` and `tables/state_timepoint_regulon_enrichment_similarity.csv`: all pairwise metrics plus interpretable TF drivers.
- `tables/*_best_matching_timepoint_by_regulon_enrichment.csv`: primary-cosine best timepoint for every MP/state.
- `tables/*_leading_concordant_regulons.csv`: pairs ranked within every MP/state with leading TFs.
- `figures/mp_timepoint_regulon_enrichment_similarity_heatmap.pdf` and `figures/state_timepoint_regulon_enrichment_similarity_heatmap.pdf`: intentionally simple fixed T0 → eR4 treatment-course heatmaps. Fill and the sole cell number are the positive enrichment cosine; a black outline marks the closest timepoint in every row. Supporting metrics remain in the CSV rather than crowding the figure.
- `figures/mp_timepoint_top_regulon_jaccard_heatmap.pdf` and `figures/state_timepoint_top_regulon_jaccard_heatmap.pdf`: separate supporting-evidence figures. They show only weighted Jaccard overlap of the top enriched TF regulons. This metric is not combined with, and does not choose, the primary enrichment-cosine match.
- `figures/mp_timepoint_regulon_match_evidence_profiles.pdf` and `figures/state_timepoint_regulon_match_evidence_profiles.pdf`: one page per MP/state demonstrating the match calculation. Each page shows its leading concordantly enriched TF-regulons across the scRef entity and every Parse timepoint. Values are within-run RSS z-scores and the black outline marks the selected closest timepoint.
- `intermediate/regulon_enrichment_similarity_intermediate.rds`: z-score matrices and all pairwise results for replotting.
- `logs/final_mp_scenic_parse_overlap_run.txt`: input and parameter provenance.
- `updates/new_updates/summaries/final_mp_scenic_parse_overlap_summary.csv`: compact, login-node-readable best-match summary.

## Run command

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cell_states/final_mp_scenic_parse_overlap.R top_n=20
```
