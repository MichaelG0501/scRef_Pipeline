# Progress Update — 11 May 2026

**Project:** scRef_Pipeline / Parse_Pipeline / PDOs_Pipeline

---

## Context
This period focused on Parse T2T4-high pattern filtering, assessing Parse RNA velocity, and systematically comparing different MP activity scoring methods (full list vs cumulative 70% vs weighted rank) to see how they impact state definitions and compositions.

## Summary of Work
* Filtered Parse T2T4-high patterns and enriched selected MPs with 3CA and Hallmarks.
* Analyzed Parse RNA velocity across samples (grid, stream, direction summaries, and sample networks).
* Investigated PDOs pairwise pattern filtering using high-resolution mean/median pair trends.
* Evaluated different approaches for scoring MP activities (Full gene list vs cumulative 70% vs weighted rank).
* Compared final state compositions and pairwise concordance heatmaps based on the different scoring methodologies.

## Key Outputs
| File/Output | Description |
|---|---|
| `Auto_parse_highres_T2T4_mean_median_trends_nMP117_selected.pdf` | Parse T2T4-high pattern filtering trends. |
| `Auto_velocity_grid_by_sample.png` | Parse RNA velocity grid plots. |
| `Auto_pdo_flot_highres_mean_median_pair_trends_nMP156_selected.pdf` | PDOs pairwise pattern filtering trends. |
| `Auto_activity_pairwise_scatter.pdf` | Comparison of MP scoring state definition methods. |
| `Auto_state_pairwise_concordance_heatmaps.pdf` | Final states pairwise concordance heatmaps for PDOs. |

## Decisions Made
1. Incorporated a weighted rank / cumulative threshold approach for MP scoring to check for state stability and composition changes.
2. Relied on pairwise concordance to determine the robustness of state definitions across different scoring metrics.

## Next Steps
- [ ] scATLAS identify subclones CNA patterns and states/MP correlation
- [ ] PDOs RNA velocity analysis
- [ ] PDOs FLOT matched pairs geneNMF analysis + trend filtering
- [ ] PDOs Numbat subclone analysis for validation
- [ ] scATLAS chemical inhibitors screening
