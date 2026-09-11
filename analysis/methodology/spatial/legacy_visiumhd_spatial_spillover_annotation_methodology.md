# Visium HD spatial spillover-corrected annotation

## Purpose

This workflow adds a fourth annotation method, `spatial`, for Space Ranger
segmented cells. It uses the same marker panels as the manual segmented and
custom-bin annotations, but corrects marker counts for spatially predicted
external RNA before assigning cell types. It does not use the hierarchical
non-epithelial/fibroblast/epithelial decision order. Cell identities are
assigned globally at the expression-cluster level.

## Inputs and filtering

For each sample, the script reads `filtered_feature_cell_matrix.h5` and the
centroids of `cell_segmentations.geojson`. Cells must have at least 200 total
UMIs and at most 15% mitochondrial UMIs. The full transcriptome is CP10K/log1p
normalised, depth and mitochondrial fraction are regressed from highly
variable genes, and a single Leiden clustering is calculated for the whole
sample. This clustering is independent of the spillover parameters and is
used to pool transcriptionally similar cells, including spatially isolated
members of a lineage.

## Spatial kernel

Cell centroids form an irregular point cloud, so the method does not rasterise
the cells. The median positive nearest-neighbour centroid distance, `d_NN`, is
the sample-specific length unit. Candidate Gaussian scales are
`sigma = {0.75, 1.25, 2.0} * d_NN`; neighbours beyond `3 * sigma` are omitted.
The diagonal is zero, so a cell never contributes to its own expected
background.

For donor cell `j` and recipient cell `i`, the unnormalised weight is

`w_ji = exp(-distance(i,j)^2 / (2*sigma^2))`.

Each donor row is normalised to sum to one over its observed recipients. This
donor normalisation gives the contamination fraction a direct interpretation:
if `rho = 0.10`, up to 10% of the donor marker molecules are spatially
redistributed among nearby recipients according to the Gaussian weights.
Boundary cells do not lose mass merely because part of a raster kernel would
fall outside the tissue.

## Global-expression calibration

Calibration is performed on raw UMI counts for every marker gene. For marker
gene `g`, the expected external count in recipient `i` is

`B_ig = rho * sum_j W_ji * X_jg`,

and the corrected count is

`X_corrected_ig = max(0, X_ig - B_ig)`.

No additional multiplication by the tissue-wide mean expression is applied.
The donor counts `X_jg` already calibrate the background to global abundance:
an abundant transcript such as `COL1A1` contributes more background molecules
than a rare transcript because more observed donor UMIs enter the convolution.
Multiplying by a global mean again would double-count abundance and
over-correct fibroblast/epithelial markers.

Candidate contamination fractions are 0.05, 0.10, 0.20, and 0.30. Every
`sigma/rho` pair is evaluated without using existing annotation labels:

1. **Residual neighbour dependence.** For each cell type, positive low-marker recipient
   cells are the lower 75% of the positive raw marker-score distribution. The absolute
   Spearman correlation between their corrected score and predicted external
   background is measured. Residual spillover remains when this is high.
2. **Source preservation.** For each cell type, the upper 10% of raw
   marker-scoring cells are source anchors. Their corrected/raw mean score
   ratio must remain high. A penalty begins below 0.75.
3. **Zero inflation.** A penalty begins if correction increases the median
   zero-score fraction by more than 0.20 across marker panels.
4. **Cluster separation.** Corrected cell-type scores are standardised within
   the sample. The median cluster-level gap between the highest and second
   highest cell-type scores rewards parameters that improve distinguishable
   marker profiles.

The objective minimised is:

`residual dependence + 2*max(0, 0.75-retention) + 1.5*max(0, zero increase-0.20) - 0.1*cluster margin`.

Thus, a dense field with raw fibroblast score 5 can receive an expected
background of 1 and retain a corrected score of 4, while an isolated score-1
cell with no neighbouring donor signal has background zero and remains 1.
The source-preservation penalty prevents choosing parameters that erase a
uniform genuine compartment.

## Cluster annotation and ambiguity

Corrected marker counts are divided by each cell's original total UMI count,
converted to log1p CP10K, and averaged within each marker panel. Per-cell panel
scores are standardised across the sample so panels with different baseline
expression are comparable. Leiden-cluster means are then calculated for all
cell types simultaneously.

Cluster-local marker support is reported as a high-confidence evidence tier: at least two genes (or all genes for panels shorter than two) must
have at least 5% detection, at least two percentage points greater detection
than the rest of the tissue, and higher mean corrected expression than the
rest of the tissue. This prevents a single abundant spillover transcript from being mistaken for high-confidence evidence. Global cluster annotation can still use a corrected panel profile without this tier when its standardised score passes the calibrated minimum; this fallback is required for spatially isolated or tissue-wide compartments that cannot show cluster-versus-rest enrichment.

Minimum standardised score values `{0, 0.25, 0.5}` and top-to-competitor gaps
`{0.15, 0.30, 0.50}` are assessed by 50 within-cluster bootstrap resamples.
The selected pair maximises label-set stability plus a coverage term weighted at 0.75
and penalty for unnecessary multi-label calls. Eligible types within the
selected gap are retained as `typeA|typeB`, but competitors beyond the top
corrected profile must pass the cluster-local marker-support tier. Clusters
with no eligible corrected type are `unresolved`. This is a global rule, not
a hierarchy.

## Downstream use

Only exact `epithelial` calls enter InferCNA. InferCNA uses the two most
abundant qualifying reference types among endothelial, macrophage, and
fibroblast, with at least 20 cells per type. Malignant levels 1 and 2 then enter
the centred scATLAS state mapping. The workflow exports per-sample annotation
and malignancy spatial maps, UMAPs, CNA scatter plots, and a four-method CNA
comparison.

## Audit outputs

The complete calibration grid, selected parameters, threshold-bootstrap grid,
cluster/type marker evidence, final annotations, sparse selected kernel, and a
diagnostic PDF are saved under `ref_outs/visium_hd_outs/`. These files permit
review of FFPED1 `t.cell`, `nk.cell`, and `macrophage` recovery without relying
on the final plot alone.

Pages 1 and 2 of each spillover diagnostic PDF are controlled comparisons.
The left panel applies the final corrected annotation pipeline to uncorrected
marker counts; the right applies it to corrected marker counts. QC cells,
whole-transcriptome Leiden clusters and UMAP, marker panels, minimum z-score,
ambiguity gap, and marker-support rules are identical. Therefore, label
changes on these pages are attributable only to marker-count subtraction.
Per-sample transition tables and changed-cell percentages are also exported.

`SCREF_RUN_MODE=downstream_only` reuses completed spatial InferCNA tables to
recover state mapping, spatial/UMAP diagnostics, and the four-method comparison
after a downstream-only failure.

## Final calibration (22 July 2026)

All three samples selected the most local tested kernel (`sigma = 0.75` times
the median nearest-neighbour distance) and `rho = 0.30`. High-source signal
retention was 0.981 for SUR1231, 0.994 for FFPEA1, and 0.967 for FFPED1;
residual background dependence was 0.216, 0.240, and 0.333, respectively.
The selected `rho` is the prespecified upper guard of the search range, so it
should be interpreted as the strongest permitted conservative correction, not
as proof that the physical contamination fraction is exactly 30% or that a
more aggressive subtraction would be justified. Extending beyond this guard
was deliberately avoided because it would extrapolate beyond the calibrated
source-preservation regime.

The bootstrap selected minimum standardised score zero in every sample. Zero
is an interpretable minimum: the winning cluster profile must be above that
marker panel's tissue-wide mean after correction. Selected ambiguity gaps were
0.50 for SUR1231/FFPEA1 and 0.15 for FFPED1, with non-top competitors still
required to pass marker support.

## Raw-versus-corrected annotation audit (23 July 2026)

Holding QC cells, whole-transcriptome clusters, UMAP coordinates, marker
panels, minimum z-score, ambiguity gap, and support rules fixed, correction
changed 0/23,368 SUR1231 labels, 1,258/75,668 FFPEA1 labels (1.66%), and
0/22,791 FFPED1 labels. The FFPEA1 changes comprise one complete cluster from
epithelial to lymph. Median positive-score reductions for fibroblast were
4.0%, 5.0%, and 8.5% in SUR1231, FFPEA1, and FFPED1; the corresponding upper
decile reductions were 17.5%, 23.4%, and 34.0%.

This shows that the selected correction changes marker magnitude but has
little effect on cluster-level identities. It does not demonstrate that the
visible halo is removed. The selected contamination fraction also lies at the
upper search boundary. Before increasing subtraction in production, evaluate
an extended strength/scale sensitivity grid using shared colour limits,
source-compartment preservation, halo-versus-distance reduction, and label
transitions. Raising `rho` alone may leave standardized cluster rankings
unchanged, because clustering is calculated from the uncorrected whole
transcriptome and annotation pools corrected scores over each fixed cluster.

## References

- Ni et al. SpotClean, *Nature Communications* (2022):
  https://www.nature.com/articles/s41467-022-30587-y
- DeSpotX spatial contamination model:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC13192984/
