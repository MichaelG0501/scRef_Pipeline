# Centred Refined MP Ordered Heatmaps Methodology

This script is a terminal plotting layer for the centred GeneNMF refinement workflow. It does not rerun NMF, UCell scoring, sub-MP splitting, or correlated sub-MP merging. It reads the finalized centred refinement caches and produces three ordered heatmaps.

The NMF clustering heatmap uses `geneNMF_metaprograms_nMP_<optimal>.rds`, where `<optimal>` is read from `ref_outs/Metaprogrammes_Results/centred/optimal_nMP.rds`. The script extracts the GeneNMF program-program Jaccard similarity matrix and orders programs first by the requested finalized centred MP order, then by the original GeneNMF hierarchical tree order within each MP. Dotted diagonal boxes mark each finalized centred refined MP block. The label shown next to each block is the requested `MP + description` string.

The MP correlation heatmap uses `merged_refined_ucell_scores.rds` from the centred refinement output directory. UCell scores are z-scaled across cells, then Spearman correlations are computed independently within each `orig.ident` sample with at least 10 cells. Sample-level correlations are Fisher-Z transformed, averaged, and back-transformed to mean rho values. A one-sample t-test over finite sample-level Fisher-Z values provides the displayed significance stars. The annotation bars use the requested state grouping: classic proliferation red, basal-to-intestinal green, SMG-to-intestinal orange, stress-adaptive purple, immune mimicry blue, and cell-cycle grey.

The program-resolution refined MP correlation heatmap uses `refined_ucell_scores.rds` before correlated sub-MP merging. It computes Fisher-Z averaged per-sample Spearman correlations across all 90 pre-merge/full refined features, then expands the feature-level correlation matrix back to NMF-program resolution using `merged_refined_mp_assignments.rds`. Final centred refined MP blocks follow the requested MP order, while each final block is internally split by its `refined_submp` program assignments. Solid internal boxes mark pre-merge split features and dotted diagonal boxes mark final centred refined MP blocks, matching the style of `analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R`.

The MP-level correlation heatmap places column labels at a 30-degree bottom angle, corresponding to a 60-degree clockwise rotation from the previous vertical label position.

The plotting order is exactly:

- Cell cycle: `MP1`, `MP5`, `MP13+`
- Classic proliferation: `MP2+`
- Basal to intestinal metaplasia: `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, `MP10+`
- SMG to intestinal metaplasia: `MP8+`, `MP8b`, `MP16`, `MP18b`, `MP17`
- Stress adaptive: `MP12`
- Cancer-cell immune mimicry: `MP15`
The plotted panel therefore has 17 MPs. `MP2x` and `MP11c` are removed upstream for coverage below three samples, while `MP18a` is the documented explicit exclusion; none is carried into these terminal heatmaps.

Outputs are written under `ref_outs/Metaprogrammes_Results/centred/mp_refinement/` using the standard `figures/`, `tables/`, and `intermediate/` tiers. A lightweight run summary is written to `updates/new_updates/summaries/centred_refined_mp_ordered_heatmaps_summary.csv`.
