####################
# Visium HD Binned Post-annotation Filtering and Malignancy

## Scope

This workflow operates only on the final 16 um binned annotations for `SUR1231`,
`FFPEA1`, and `FFPED1`. `SUR1122` and every segmented-cell result are excluded.
It does not alter the upstream annotation tables.

Run order:

1. `analysis/spatial/visium_hd_binned_post_annotation_filter.R`
2. `analysis/spatial/visium_hd_binned_malignancy.R`

The PBS wrapper `analysis/spatial/run_visium_hd_binned_malignancy.sh` runs both
steps in this order.

## Post-annotation Filter

The input marker scores are the fixed-panel scores already saved by the final
annotation pipeline. Each score is the mean `log1p(CP10K)` expression across the
marker genes for that cell type. No marker list, normalization, cluster, or
annotation is recomputed.

For a combined annotation such as `endothelial|fibroblast`, the first label is
the principal annotation and is used as the home lineage. The complete original
label remains in `Auto_annotation_celltype`.

The filter has two checks:

1. **Home-marker expression.** The home score must be strictly greater than
   zero. This is the least stringent meaningful expression gate: one or more
   genes from the assigned marker panel must have a nonzero count. It is more
   appropriate than a fixed score of 0.25 at this stage because 16 um bins have
   variable and sometimes low UMI depth.
2. **Unrelated-lineage coexpression.** A bin fails only when the strongest
   unrelated nonstructural competitor has a score of at least 1.0 and is at
   least as high as the home score. Epithelial and fibroblast scores are never
   considered off-lineage competitors because their abundant transcripts are
   the principal expected spillover sources. Closely related pairs are also
   ignored: T/NK, B/plasma, macrophage/dendritic, and
   keratinocyte/epithelial.

This is deliberately permissive. It removes bins lacking direct support for
their assigned type and clear unrelated coexpression, but does not recluster,
relabel doublets, or create an unresolved rescue step. Every decision, score,
competitor, threshold, and reason is written to the filtered per-bin table.

The cell-type count report shows one page per sample. Each page contains
separate before- and after-filtering vertical bar panels with actual bin counts;
it does not use stacked bars. Pages are landscape, the y axis is the observation
count, and one fixed y-axis scale is used across the left and right panels.

Two additional count reports compare cell segmentation on the left with 16 um
bins on the right. One PDF uses pre-filter annotations and the other applies
the exact same permissive expression/coexpression filter to both methods.
These comparison-only segmented filter calls do not enter malignancy analysis,
which remains strictly binned. Within each sample page, both methods use
vertical bars, share the same count y-axis, and show the actual count above
every bar.
Cell types with zero observations in both methods on a sample page are omitted;
a zero in only one method remains visible for direct comparison.

Before-filter comparison counts preserve the exact finalized annotation label,
including pipe-delimited combined labels such as `erythrocyte|fibroblast`. This
is the same category definition used by the pre-filter annotation-diagnostic
legend. After-filter counts use `Auto_postfilter_celltype`, because filtering
validates the principal first label as the home lineage; these counts therefore
match the post-filter diagnostic legend. No combined-label observation is added
to a single-lineage bar in a pre-filter figure.

All Visium HD cell-type bars, spatial maps, UMAPs, and interactive audit plots
read `analysis/shared/visium_hd_celltype_colours.tsv`. Single-lineage colours
therefore exactly match the production annotation diagnostics. Pipe-delimited
combined labels use the shared neutral grey `#555555`.

The canonical annotation diagnostics in
`ref_outs/visium_hd_outs/figures/annotation_diagnostics/` remain the unfiltered
annotation-stage outputs. Matched post-filter spatial and assignment-UMAP
diagnostics are written separately under its `after_filtering/` subdirectory.
`visium_hd_postfilter_diagnostics.py` calls the same production plotting
function used for the pre-filter annotation diagnostics. It therefore preserves
the exact `19 x 9` layout, spatial/UMAP aspect handling, colour map, alpha,
right-side count legend, rasterization, typography, and point-size formula.
For direct comparison, point size is calculated from the original pre-filter
observation count for that sample. Only retained bins are drawn.

## InferCNA Inputs and References

Only post-filter epithelial bins are malignancy targets. Candidate normal
references start from post-filter endothelial, macrophage, and fibroblast bins,
but a bin is eligible only when all of the following independent evidence
agrees: the finalized manual annotation is a single exact label, the RCTD call
is a singlet, and the RCTD first type equals that manual label. The annotation
itself is not changed by this reference-only gate.

Exactly the two most abundant concordant candidate types with at least 15 bins
are selected. The minimum of 15 retains the independently concordant SUR1231
macrophage group (16 bins), comparable to the 21 macrophages used by the legacy
RCTD binned run. If two types are not available, the sample is recorded with
counts and the workflow stops after attempting the remaining samples.

Every qualifying bin from each selected group is retained. The installed
InferCNA `refCorrect()` implementation calculates a separate mean profile for
each reference group and corrects against the range between group means;
therefore an abundant fibroblast group does not numerically weight the
endothelial/macrophage group. Retaining all qualifying bins estimates each
group mean more precisely than the previous random cap. The selected barcodes
are saved in live storage.

Only epithelial targets and selected reference bins enter InferCNA. Counts are
converted to CPM, restricted to genes recognized by the hg38 gene-order input,
and passed to `infercna(..., isLog = FALSE)`. The expensive full InferCNA matrix
is a rebuildable cache in ephemeral storage. All per-bin CNA metrics,
thresholds, calls, selected reference barcodes, signature genes, tables, and
figures needed for audit and replotting are retained in live storage.

## Malignancy Calling

`cnaScatterPlot(..., refCells = ref_barcodes)` supplies CNA signal and CNA
correlation while excluding the normal references from the average CNA profile
used for correlation. The legacy `excludeFromAvg` argument was invalid for the
installed InferCNA version and was silently forwarded as a plotting parameter;
it is not retained. Both thresholds are
estimated from the pooled selected reference bins as reference mean plus one
standard deviation. This retains the legacy two-metric evidence rule while
avoiding the previously over-stringent two-standard-deviation threshold.

- **CNA malignant:** both CNA signal and CNA correlation exceed their thresholds.
- **CNA unresolved:** exactly one metric exceeds its threshold.
- **CNA non-malignant:** neither metric exceeds its threshold.

CNA-malignant epithelial bins are `malignant_level_1`.

The secondary cancer signature follows the scATLAS `Malignancy.R` procedure.
The union signature is read from `ref_outs/cancer_signatures.txt`; within each
sample, the 50 genes with the highest mean `log1p(CP10K)` expression among
post-filter epithelial targets are retained. A bin score is the mean across
those genes, and a score of at least 1 is positive. Only CNA-unresolved bins are
rescued as `malignant_level_2`. A signature-positive bin with both CNA metrics
below threshold is not called malignant, because the signature alone is not
specific enough to override normal CNA evidence.

Keratinocytes are not malignancy targets. The legacy segmented-to-binned
keratinocyte projection, nearest-bin rescue, cross-method reference sharing,
and later keratinocyte profile overrides are intentionally omitted: those were
corrective experiments for older RCTD-only epithelial labels and are not
appropriate for the current direct binned annotation.

## Parameters and Reproducibility

The final parameter tables are written beside the result tables. Set
`SCREF_FORCE_REBUILD=TRUE` to ignore compatible InferCNA caches. Set
`SCREF_REPLOT_ONLY=TRUE` to regenerate diagnostic figures from the live per-bin
malignancy tables without reading count matrices or ephemeral caches.

Final outputs are under:

- `ref_outs/visium_hd_outs/post_annotation_filter/`
- `ref_outs/visium_hd_outs/malignancy/`

Each malignancy result has an InferCNA scatter plot and a spatial malignancy
map in the multi-page diagnostic PDF. The spatial map is additionally saved as
a standalone presentation-readable PDF and PNG for each sample. Spatial points
use size `0.45` and alpha `0.72`, with background/normal points drawn before
malignant levels, to reduce occlusion in dense 16 um regions.


####################
The combined malignancy count figure uses every post-filter epithelial target
in each of the three samples as its denominator. It shows one separate
horizontal stacked panel per sample and prints the absolute bin count for every
segment. The four final states are non-malignant, unresolved, malignant level
1, and malignant level 2. Unresolved bins remain explicit because one CNA
metric exceeds its reference threshold; they are not silently counted as
non-malignant. Colours exactly match the InferCNA and spatial malignancy plots.
The underlying count and percentage table is saved beside the per-bin tables.
####################

Compact summaries are also written to `updates/new_updates/summaries/`.

## Reference-selection Audit

`visium_hd_binned_reference_selection_audit.R` compares the previous random
post-filter reference, the final manual/RCTD-concordant reference, and the
legacy RCTD reference. The final reference thresholds are close to the legacy
binned thresholds and substantially tighter than the previous current run:

| Sample | Reference definition | Signal threshold | Correlation threshold |
|---|---:|---:|---:|
| SUR1231 | Previous random-balanced | 0.002565 | 0.2463 |
| SUR1231 | Final concordant | 0.000119 | 0.0172 |
| SUR1231 | Legacy RCTD | 0.000229 | -0.0047 |
| FFPEA1 | Previous random-balanced | 0.000653 | 0.0615 |
| FFPEA1 | Final concordant | 0.000255 | 0.0246 |
| FFPEA1 | Legacy RCTD | 0.000292 | 0.0312 |
| FFPED1 | Previous random-balanced | 0.003604 | 0.1335 |
| FFPED1 | Final concordant | 0.000764 | -0.0374 |
| FFPED1 | Legacy RCTD | 0.000880 | -0.0257 |

This shows that the legacy binned advantage was principally caused by a purer
normal reference definition rather than by the epithelial annotation itself.
####################
