# Progress Update — 26 March 2026

**Project:** scRef Pipeline — OAC Single-Cell Reference Atlas

---

## Context

This cycle focused on consolidating pseudotime trajectory analysis, finalising MP/state ordering consistency in correlation and Jaccard heatmaps, updating survival volcano plots to per-MP resolution, adding statistical annotations to EAC vs ESCC comparisons, and extending PDO cross-validation with clinical abundance plots.

## Summary of Work

- **Pseudotime trajectories** (Part A) — unified Monocle3 pipeline producing per-sample trajectory reports.
- **Pairwise state distance** — expanded with a summary of **five methodologies** (Directed Pseudotime Mean/Median, Principal Graph Geodesic Centroid/Medoid, and UMAP Euclidean Centroid).
- **PDO sample abundance** — finalised stacked bar abundance plots for both MPs and states with clinical top annotation bars.
- **Comparative Similarity heatmaps** — presented MP Spearman Correlation and Jaccard similarity **side-by-side** for direct comparison.
- **Final outputs update** — per-cell heatmap and cell cycle boxplots updated to the latest `task4_unresolved_relabel` versions.
- **Survival volcano (per-MP)** — updated to show individual MPs, continuous Cox volcano across 4 panels.
- **EAC vs ESCC delta** — updated comparison showing Wilcoxon rank-sum significance.
- **PDO vs scRef expression** — presented cross-dataset per-MP expression correlations **side-by-side** (PDO-based vs. scATLAS-based).

## Key Outputs

| File | Description |
|---|---|
| `ref_outs/Auto_pseudotime_trajectory_summary_reports_partA.pdf` | Monocle3 pseudotime trajectory reports (12 samples) |
| `ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_distance_heatmap_all_methods.pdf` | Pairwise state distance heatmap (5 methods) |
| `ref_outs/task6_hybrid_pairwise_distance/Auto_task6_hybrid_pairwise_nodeplot_all_methods.pdf` | Network nodeplot of state distances |
| `PDOs_outs/sample_abundance/Auto_sample_abundance_pdo.pdf` | PDO abundance (MPs + states + clinical) |
| `ref_outs/Metaprogrammes_Results/Auto_final_mp_correlation_heatmap.pdf` | State-ordered MP Spearman correlation |
| `ref_outs/Metaprogrammes_Results/Auto_final_mp_jaccard_heatmap.pdf` | State-ordered MP Jaccard similarity |
| `ref_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_heatmap.pdf` | Final per-cell heatmap (post-relabel update) |
| `ref_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_cc_boxplot.pdf` | Cell cycle boxplot (6 states, post-relabel update) |
| `ref_outs/task2_survival/Auto_task2_survival_volcano_continuous.pdf` | Per-MP continuous Cox volcano plots |
| `ref_outs/task8_tcga_eac_escc_compare/Auto_task8_tcga_eac_escc_compare_plots.pdf` | EAC vs ESCC delta with Wilcoxon statistics |
| `PDOs_outs/Auto_PDO_mp_crossref_correlation_meta.pdf` | Per-MP expression correlation (in PDOs) |
| `PDOs_outs/Auto_PDO_mp_crossref_in_sc_correlation_meta.pdf` | Per-MP expression correlation (in scATLAS) |

## Decisions Made

1. Pseudotime trajectories consolidated into a single unified script producing multi-panel PDF reports per sample
2. MP ordering in correlation and Jaccard heatmaps now matches state grouping for visual consistency
3. Survival volcano plots now show individual MPs rather than merged state-level groupings
4. EAC vs ESCC delta plot now includes per-sample Wilcoxon rank-sum test statistics

## Next Steps

- [ ] SCENIC regulon analysis on final MPs (pending doMC fix)
- [ ] Finalise PPTX conversion and upload after content confirmation
- [ ] Extend pseudotime analysis to Part B subsets if needed
- [ ] PDO matched pair temporal analysis refinement
