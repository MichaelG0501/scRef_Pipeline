# Centred GeneNMF and Refined MP Methodology

This document covers the canonical centred metaprogram workflow in `analysis/metaprograms/centred/01_centred_geneNMF.R`, `02_nmf_rank_selection_diagnostics.R`, `03_mp_refinement_submp.R`, and `04_mp_refinement_merge_correlated_submps.R`.

## Centred GeneNMF programs

Step 01 loads each per-sample epithelial object ending in `_epi_f.rds`, retains `malignant_level_1` and `malignant_level_2` cells, and excludes samples with fewer than 10 malignant cells. `GeneNMF::multiNMF()` is run with `assay="RNA"`, ranks `k=4:9`, `min.exp=0.05`, and `center=TRUE`. GeneNMF subtracts each gene mean and truncates negative centred values to zero before non-negative factorisation. Candidate metaprogram solutions are extracted for nMP 8 through 30 with cosine similarity, specificity weight 5, cumulative weight explained 0.5, and minimum confidence 0.5. The rank diagnostic currently selects nMP=19 from the silhouette-curve knee.

## Rank selection diagnostic

Step 02 evaluates every available nMP solution from 8 through 30 using the same program-program cosine similarity matrix and hierarchical tree stored by GeneNMF. It converts similarity to distance as `1 - cosine similarity` and cuts the tree at the candidate nMP. Average silhouette width measures how much closer each NMF program is to its assigned MP than to neighbouring MPs. As an independent compactness diagnostic, within-cluster sum of squares is the sum of squared pairwise cosine distances divided by twice the number of programs in each cluster, summed across clusters.

For each curve, nMP and the diagnostic value are independently scaled to [0,1]. The diagnostic knee is the point with the greatest perpendicular distance from the straight line joining the first and last curve points. The silhouette knee is the prespecified selector because silhouette jointly reflects within-MP cohesion and separation from alternative MPs; the WSS knee is reported as a sensitivity diagnostic and does not override it. All candidate nMP objects must be present, because a missing candidate would interrupt the continuous curve and make the geometric knee incomparable. The current complete candidate set selects nMP=19. Negative-silhouette MPs are excluded from the initial enrichment plots because they are less similar to their assigned cluster than to an alternative cluster; the later refinement steps apply the more detailed retain/split/remove rules below.

## Parent-MP triage and splitting

Step 03 converts `sampleCoverage` to the number of contributing samples and applies three simultaneous quality criteria:

- retain without splitting when silhouette is at least 0.2, coverage is at least 3 samples, and the consensus has more than 5 genes;
- split when silhouette is strictly between 0 and 0.2 and the same coverage/gene criteria pass;
- remove when silhouette is below 0, coverage is below 3 samples, or the consensus has at most 5 genes.

For each split candidate with at least four NMF programs, the program cosine-similarity submatrix is converted to distance `1 - similarity` and clustered with Ward.D2. Candidate cuts range from k=2 to one fewer than the number of programs. The first k whose mean silhouette reaches 0.2 is selected; if no cut reaches 0.2, the cut with maximum mean silhouette is used. Candidates with fewer than four programs remain unsplit. Sub-MP letters follow left-to-right dendrogram order.

Consensus signatures reproduce the GeneNMF weighting logic: program loadings are specificity-weighted with exponent 5, genes within 50% cumulative average weight are retained, genes must occur in more than half the contributing programs, and at most 200 genes are returned. UCell scores are calculated on `EAC_Ref_epi.rds`. Correlations are computed per `orig.ident` sample with at least 10 cells, Fisher-Z averaged across samples, and back-transformed. Jaccard overlap supplies an orthogonal gene-list diagnostic.

`split_results.rds` is a required input to step 04, so it is saved persistently in the live output tree and duplicated to ephemeral only as a rebuild cache. `SCREF_FORCE_REBUILD=TRUE` invalidates caches; `SCREF_REPLOT_ONLY=TRUE` requires a compatible cache.

## Correlated sub-MP merging and final filter

Step 04 considers sub-MPs only within their original parent MP. A sub-MP qualifies for merging when its Fisher-Z mean Spearman rho exceeds 0.4 with at least 25% of the other sub-MPs in that parent. A parent-plus feature is created only when at least two sub-MPs qualify. Non-qualifying sub-MPs remain separate. The merged consensus is rebuilt from the pooled contributing NMF programs rather than taking a union of previously truncated gene lists.

After merging, final refined MPs must be represented by programs from at least 3 samples and contain at least 5 genes. `MP18a` is additionally excluded as the recorded low-coverage/quality decision. In the current result, `MP2x` and `MP11c` fail the sample-coverage criterion and `MP18a` is explicitly removed, producing a 17-MP panel.

Enrichment uses BH-adjusted q-value 0.05 for Hallmark, GO Biological Process, 3CA MPs, and the persistent developmental/adult epithelial reference collections. The complete enrichment objects are saved to live as `intermediate/cluster_enrich_centred.rds`; heatmaps and summary tables are derived from this object.

## Persistent outputs and reruns

All gene lists, weights, assignments, UCell matrices, correlations, enrichment results, tables, and final figures are written under `ref_outs/Metaprogrammes_Results/centred/mp_refinement/`. These files are required for downstream state definition or replotting and must remain in live storage. Only duplicate, fully rebuildable caches may be ephemeral. These steps are HPC workflows and must be run through PBS, not on a login node.
