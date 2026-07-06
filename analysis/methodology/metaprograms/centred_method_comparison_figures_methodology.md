####################
# Centred Method Comparison Figures Methodology
####################

`analysis/metaprograms/centred_method_comparison_figures.R` compares the historical uncentred GeneNMF workflow with the centred GeneNMF workflow. It is a terminal plotting workflow and does not modify the upstream GeneNMF, rank-selection, or refinement scripts.

The uncentred initial nMP19 object and filtered UCell matrix are read from the historical read-only cache under `/rds/general/project/tumourheterogeneity1/live/temp_save/`. For uncentred rank-selection diagnostics, the script caches the historical nMP19 tree/similarity object under `ref_outs/Metaprogrammes_Results/centred_comparison/intermediate/` and cuts that tree at k = 8-30, matching the diagnostic calculation that uses `programs.similarity` and `programs.tree`.

Initial MPs are filtered by silhouette score, retaining MPs with silhouette >= 0. Programmes with missing metaprogram assignment are excluded from the unsupervised programme-similarity heatmaps, matching the previous MPNA-removal plotting logic. Rank-selection diagnostics recompute average silhouette width and within-cluster sum of squares for nMP 8-30 and mark the kneedle-style inflection/elbow.

3CA enrichment uses the 3CA MP gene sets from `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv` and `clusterProfiler::enricher()`. Heatmap values are capped `-log10(adjusted p)` scores.

Expression-score correlations are Spearman correlations of UCell scores. For initial centred MPs, UCell scores are computed sample-by-sample from `ref_outs/by_samples/<sample>/<sample>_epi_f.rds` for the same cells present in the historical uncentred UCell matrix and cached in the comparison intermediate directory. Within-method correlation heatmaps use per-sample Fisher-Z averaging. Cross-method comparison heatmaps use all matched cells.

Label transfer for centred MPs is based on gene-list Jaccard similarity to uncentred MPs. When the best Jaccard is at least 0.25, the centred MP keeps its own MP number but receives the matched uncentred description. Below that threshold, the centred MP label is left unchanged. The transfer tables are written under `tables/`.
