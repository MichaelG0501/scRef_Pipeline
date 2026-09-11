# Final Visium HD cell-type annotation

####################
## Scope

The final workflow contains two representations and one shared annotation
algorithm:

1. `binned`: Space Ranger 16 um bins. Every bin retained by the 100-UMI RCTD
   run enters score-based annotation; RCTD singlet/doublet calls are preserved
   as metadata and are not an exclusion gate.
2. `segmented`: Space Ranger cell-segmented expression. Segmented cells pass
   the same score-based annotation after count and mitochondrial QC.

The sample manifest is `analysis/spatial/visium_hd_samples.tsv`. It contains
`SUR1122`, `SUR1231`, `FFPEA1`, and `FFPED1`.

## RCTD metadata

`visium_hd_rctd_doublet_detection.R` uses raw binned UMI counts, a 100-UMI
minimum, doublet mode, and `ref_outs/EAC_Ref_merged.rds`. The reference removes
unresolved cells, combines T and NK cells only for RCTD fitting, and caps
reference sampling at 2,000 cells per non-epithelial type and 6,000 epithelial
cells. RCTD classification remains strictly `singlet`, `doublet_certain`,
`doublet_uncertain`, or `reject`, but all classes enter the custom expression
graph and annotation. The RCTD class and first/second types remain separate
metadata and never determine the custom cell-type label.

## Shared preprocessing and graph

Segmented cells require at least 100 UMIs and less than 15% mitochondrial counts.
Binned observations have already passed the RCTD 100-UMI input gate. Counts are
normalized to 10,000 per observation and transformed by `log1p`.

One expression graph is built per sample and representation:

- genes detected in at least 10 observations;
- 3,000 Seurat-style highly variable genes;
- regression of total counts and mitochondrial percentage;
- scaling capped at 10;
- up to 50 principal components, with up to 40 used;
- 15 cosine neighbours;
- Leiden resolution 10 when at least 1,000 observations are analyzed, otherwise
  resolution 6;
- UMAP `min_dist=0.5`, `spread=1`, random seed 0. The interactive R audit
  embeds the exact Leiden SNN with deterministic random initialization because
  spectral initialization is undefined for disconnected tissue graphs in the
  installed UMAP/scikit-learn versions.

Leiden clustering and the displayed UMAP use the same neighbour graph. There
is no separate presentation-only embedding. This guarantees that the UMAP
used to inspect assignments represents the exact graph on which clusters were
defined. The size-aware high resolution avoids diluting small marker-positive
populations in the large datasets without reducing the median SUR1122 cluster
size below the already sparse value of five. Clusters with the same final
identity are visually merged by their common annotation colour.

The resolution rule was selected from an FFPEA1 segmented audit at resolutions
6, 8, 10, and 12. At resolution 6, old C92 contained 3,311 observations and
only 298 (9.0%) exceeded T-cell score 1; old C78 contained 1,352 observations
and only 128 (9.5%) exceeded keratinocyte score 1. Thus the reported 24.8% for
C92 was the proportion of all high-scoring observations contributed by C92,
not the within-cluster percentage. Lowering the raw threshold would have
incorrectly relabelled both complete mixed clusters.

Resolution 10 was the lowest tested value that separated both populations. It
produced a 365-observation T-cell cluster with standardized score 1.410 and
raw score 0.267, and a 361-observation keratinocyte cluster with standardized
score 3.586 and raw score 0.323. It generated no clusters below 20 observations
in FFPEA1. Across the other datasets with at least 1,000 observations,
resolution 10 generated no clusters below five. In SUR1122, however, it
reduced median cluster size from five to three, so SUR1122 retains resolution
6. Audit results are stored under
`ref_outs/visium_hd_outs/resolution_audit/`.

## Marker scores

Marker panels are unchanged from the prior finalized scripts. Each
observation receives a raw normalized score equal to the mean log1p-CP10K
expression of the panel genes. Each panel is then standardized across all
analyzed observations in that sample and representation. Both raw and
standardized observation scores are averaged within every Leiden cluster.
The panel mean and standard deviation are fitted exactly once, before any
targeted refinement. Child clusters reuse those original per-observation
standardized values; neither the number nor composition of child clusters can
change the scale.

Marker-detection support is exported as a diagnostic only. It does not gate
assignment or targeted refinement. DGE, hierarchical lineage protection, and
per-cell coexpression filtering are not used.

## Two-step assignment

Assignment is strictly score based:

1. Rank every cell type within a cluster by decreasing cluster-average
   standardized marker score.
2. Starting at rank 1, test the corresponding cluster-average raw normalized
   score against `> 0.25`.
3. Select the first ranked type that passes. Continue to rank 2, rank 3, and
   so on whenever a higher standardized candidate fails its raw threshold.
4. Assign `unresolved` only if every panel has raw normalized score `<= 0.25`.

The raw threshold prevents a relative standardized winner with negligible
absolute transcript support from determining identity. Unlike the previous
lineage-specific z-score floor, it applies identically to every cell type and
therefore does not force a fixed fraction of a dominant lineage to fail.
`0.25` is the current explicit operational threshold requested for auditing;
it is exposed in both the production parameter table and the interactive
script rather than being hidden in a fallback. It is a cluster average of
per-observation mean log1p-CP10K scores, not a direct count of 0.25 marker
reads per cell. This distinction is why partition refinement, rather than a
lower threshold, is used for a marker-positive minority inside a mixed
cluster.

## Targeted expression-graph refinement

For each initial cluster, the ranked assignment traversal records only the
cell types rejected before the first raw-score-passing type. For example, if
T cell and plasma rank first and second but both fail raw `> 0.25`, followed by
a passing fibroblast score, T cell and plasma are the two recorded rejected
candidates. Types ranked below the first passing assignment are not refinement
candidates. If no type passes, every ranked type is recorded.

For each recorded candidate, the pipeline counts observations in the parent
cluster with individual raw normalized marker score `> 1.0`. A parent is
eligible when any recorded candidate has at least 20 such observations. The
fixed count of 20 is retained as a conservative guard against reclustering for
a few isolated high-score observations; the former 8% prevalence, replicated
marker-support, structural-lineage, child-size, and capture-fraction gates are
removed.

Each eligible parent receives one locally recomputed whole-transcriptome graph
using the same normalization, HVG selection, regression, PCA, and neighbour
settings as the global graph. Marker scores trigger refinement but are not
included in this graph. Leiden is run once at resolution `1.0`. Every resulting
local community is retained as a separate provenance-labelled child, such as
`C120_a`, `C120_b`, and `C120_c`; communities are not merged according to marker
scores.

After all eligible parents have been split, cluster-average raw and
standardized marker scores are calculated for every child. Only these cluster
means are new. The per-observation standardized scores are not recalculated:
they remain on the scale fitted from the original full sample before any split,
so children are evaluated relative to the same marker distributions as the
major clusters.
The normal two-step assignment is then repeated: rank by standardized score
and select the first type with raw score `> 0.25`. This is the terminal pass;
children are never recursively reclustered. A split may therefore change a
cell-type call, retain the parent's call, or remain unresolved.

This procedure has one intentional limitation: 20 high-scoring observations
can still represent spatial spillover rather than a real population. The
whole-transcriptome local graph provides the independent partitioning evidence,
while the terminal raw-score validation prevents a child from being assigned
solely because it triggered refinement. The refinement audit records rejected
candidates, per-candidate high-score counts, eligible candidates, child IDs,
and terminal status.

## Structural doublet labels

Epithelial and fibroblast are treated as structural types only for reporting
an ambiguous cluster, not for changing score rank. If epithelial or fibroblast
is the selected primary type and the highest-ranked non-structural type also
passes the raw threshold, their standardized-score difference is calculated.
When the difference is at most `0.25`, the output is
`nonstructural|structural`, for example `macrophage|fibroblast`. If macrophage
itself is the highest valid standardized candidate, the call remains
`macrophage`. No normalization is altered to force a preferred label.

Pipe labels are retained in the output for manual review. They are not removed
by a cell-level coexpression filter. RCTD class does not exclude binned
observations.

## Shared Cell-type Colours

Production annotation diagnostics and downstream Visium HD cell-type plots read
`analysis/shared/visium_hd_celltype_colours.tsv`. This prevents independent R
and Python palettes from assigning different colours to the same lineage.
Pipe-delimited combined labels use the shared `combined` colour (`#555555`).

## Diagnostics and interactive audit

For every sample and representation, production writes:

- `tables/Auto_<sample>_<method>_cell_annotations.csv.gz`;
- cluster score calls with full rank order, initial winner, selected rank, raw
  and standardized scores, doublet partner, and score gap;
- a targeted-refinement audit containing rejected candidates, their individual
  high-score counts, parent eligibility, fixed local resolution, and child IDs;
- cluster/cell-type score evidence and marker detection diagnostics;
- the exact assignment UMAP coordinates;
- spatial/UMAP annotation figures with large points and legend keys;
- a per-cell-type cluster-assignment UMAP PDF;
- paginated raw and standardized cluster score heatmaps labelled by cluster,
  final assignment, and observation count.

`visium_hd_annotation.R` is a short, function-free, sequential audit for one
selected sample and method. It exposes global Leiden resolution, raw score,
structural-doublet gap, the individual score-1 trigger, the minimum of 20 high
observations, and local resolution 1.0. It builds one internally consistent
graph for both assignment and display, prints parent/child refinement and
cluster calls, and produces a cluster-numbered UMAP, spatial map, and score
heatmaps. The former comprehensive
audit is retained as
`analysis/spatial/legacy_visiumhd/legacy_visium_hd_annotation_full_audit.R`.
Native `leidenbase` clusters the Seurat SNN. Seurat's graph UMAP is connected
through `reticulate` to the existing jupyter Python environment. It uses
deterministic random initialization to handle disconnected SNN components.
Because the editable audit uses Seurat/leidenbase whereas production uses
Scanpy, numeric cluster IDs and exact partitions can differ; final pipeline
labels and IDs must be read from the canonical Python outputs.

## Execution

Submit `analysis/spatial/run_visium_hd_annotation.sh` through PBS. The wrapper
uses live logging, reuses completed RCTD results where possible, and rebuilds
both annotation representations. `SCREF_REUSE_COMPLETE=TRUE` reuses complete
annotation tables; combine it with comma-separated `SCREF_FORCE_SAMPLES` to
rebuild selected samples only.
####################
