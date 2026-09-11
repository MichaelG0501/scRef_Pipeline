####################
# Visium HD malignant-state/MP–TME spatial interaction methodology

## Scope

This terminal workflow tests whether malignant centred-refined states or MPs
occur next to particular non-malignant cell types or non-malignant GeneNMF MPs
more or less often than expected from their within-section abundance and
spatial opportunity. It uses only filtered 16 µm bins from `SUR1231`, `FFPEA1`,
and `FFPED1`. `FFPED1` is the canonical manifest identifier for the requested
sample D1, whose Space Ranger input folder is named `D1`. The three samples are
never joined into a common spatial graph.

The analysis script is
`analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R`; its PBS wrapper
is `analysis/spatial/visium_hd_spatial_cancer_tme_interactions.sh`.

## Four analysis modes

The workflow mirrors the biologically relevant cancer–TME branches of the
single-cell interaction workflow:

1. `01_cancer_mps_vs_tme_mps`: malignant centred-refined MPs versus
   silhouette-retained non-malignant GeneNMF MPs.
2. `02_cancer_states_vs_tme_mps`: malignant centred-refined states versus
   silhouette-retained non-malignant GeneNMF MPs.
3. `03_cancer_mps_vs_whole_celltypes`: malignant centred-refined MPs versus
   whole final annotated non-epithelial cell types.
4. `04_cancer_states_vs_whole_celltypes`: malignant centred-refined states
   versus whole final annotated non-epithelial cell types.

All four modes are cancer-to-TME tests. Within-cancer and within-TME edges are
not introduced because they are outside the requested validation question.

## Inputs and gates

### Malignant source bins

The source table is
`ref_outs/visium_hd_outs/state_mapping/tables/Auto_visiumhd_binned_malignant_state_annotations.csv.gz`.
That upstream table already restricts bins to final epithelial targets with
`Auto_malignancy` equal to `malignant_level_1` or `malignant_level_2`. It
contains the centred-refined state assignment and all adjusted malignant MP
scores.

The spatial mapper retains historical display labels, so this workflow applies
an explicit label-only crosswalk to the current centred-refined names:

- `Classic Proliferative` -> `Classic proliferation`
- `Basal to Intestinal Metaplasia` -> `Basal to intestinal metaplasia`
- `SMG-like Metaplasia` -> `SMG to intestinal metaplasia`
- `Stress-adaptive` -> `Stress adaptive`
- `Immune Infiltrating` -> `Cancer-cell immune mimicry`

MP memberships and state compositions are unchanged. `Hybrid` and `Unresolved`
are not tested as biological state nodes. Malignant MP positivity is adjusted
MP score greater than `0.5`. `MP11c` and `MP18a` are excluded, matching the
active single-cell cancer-MP interaction workflow; cell-cycle MPs remain
auditable malignant MP nodes.

### TME target bins

Final filtered cell-type calls and coordinates come from each live binned
malignancy table under
`ref_outs/visium_hd_outs/malignancy/tables/Auto_<sample>_binned_malignancy.csv.gz`.
Only `Auto_postfilter_keep == TRUE` bins are used. Whole-cell-type targets are
exact `Auto_postfilter_celltype` calls excluding `epithelial` and `unresolved`.
This permits additional spatially observed cell types beyond the seven
GeneNMF compartments when they pass the minimum-bin requirement.

### Non-malignant MP definitions and scoring

The GeneNMF definitions come from:

- `nmf_fibroblast/MP_outs_default.rds`
- `nmf_endothelial/MP_outs_default.rds`
- `nmf_macrophage/MP_outs_default.rds`
- `nmf_nk/MP_outs_default.rds`
- `nmf_plasma/MP_outs_default.rds`
- `nmf_cd4/MP_outs_default.rds`
- `nmf_cd8/MP_outs_default.rds`

Before scoring, every MP with silhouette below zero is removed. This is the
mandatory repository-wide MP silhouette filter. Up to the first 100 ranked
genes are retained per MP.

Each MP is scored only in bins carrying its parent annotation. CD4 and CD8 MPs
are both scored within final `t.cell` bins, providing subtype/program-level
resolution without inventing a whole-bin CD4/CD8 annotation. Per sample and
parent cell type, counts are normalized to CP10K and transformed with `log1p`.
Every available signature gene is standardized across those bins, the MP score
is the unweighted mean of its gene z-scores, and each MP score is standardized
again across the same sample/parent-bin population. A bin is MP-positive when
the adjusted score is greater than `0.5`.

The resulting score table, filtered gene sets, and membership matrices are
saved persistently in live storage because they are required to reproduce the
plots without rescoring Space Ranger matrices.

## Spatial neighbourhood definition

The 16 µm grid row and column are parsed from the Space Ranger barcode, for
example `s_016um_00107_00066-1`. The primary graph uses queen adjacency:
all eight possible bins at Chebyshev distance one from each malignant source
bin. This is a true immediate-bin contact test and cannot bridge gaps as a
fixed-k nearest-neighbour graph can.

Sensitivity analyses repeat the full test at cumulative Chebyshev rings two
and three, corresponding to neighbourhood extents of approximately 32 and
48 µm. These are sensitivity modes; ring one is the primary result.

## Sample-specific permutation inference

For a source feature and target feature, the observed statistic is the number
of directed source-to-target adjacency edges. The malignant source labels and
the tissue graph remain fixed.

The target null differs by target family:

- Whole-cell-type mode: final cell-type labels are permuted jointly across all
  retained non-epithelial target bins in that sample. This preserves every
  cell-type count and the exact source degrees/tissue mask.
- TME-MP mode: the joint MP-positive membership matrix is permuted only among
  bins of the same annotated parent cell type in that sample. This preserves
  the parent cell type's spatial distribution, every MP-positive count, and
  within-compartment MP overlap/correlation while testing MP-specific spatial
  placement.

The default is 499 permutations and can be changed with
`SCREF_SPATIAL_PERMUTATIONS`. Permutation sums and sums of squares give a null
mean, null standard deviation, and z-score. Two-sided normal-tail p-values are
used for BH correction because their resolution is not capped at
`1 / (permutations + 1)` across the large interaction family. Two-sided
empirical permutation p-values are also exported as calibration diagnostics.
BH correction is applied independently within sample, biological mode, and
neighbour ring.

Pairs are skipped rather than treated as negative evidence when they fail any
minimum:

- 20 source-positive malignant bins;
- 20 whole-cell-type target bins;
- 20 parent cell-type bins before MP scoring;
- 10 MP-positive target bins;
- 20 observed source-to-TME neighbour slots.

Thus rare populations in one sample do not prevent the other samples from
being tested.

## Pooled three-sample evidence

The pooled analysis never creates cross-sample neighbours. For each pair tested
in at least two samples, the independently derived sample z-scores are combined
with weighted Stouffer meta-analysis, using the square root of the number of
source-positive bins as the weight. Pooled log2 enrichment is calculated from
summed observed contacts relative to summed sample-null means. BH correction is
then applied within biological mode and ring.

Recurrence is reported separately. A recurrent same-direction interaction has
the same enrichment direction in all tested samples and sample FDR below 0.10
in at least two samples. This prevents a strong pooled result from being
misdescribed as recurrent when driven by only one section.

## Validation against the scRNA interaction workflow

Spatial pairs are matched to the corresponding live scRNA modes under
`ref_outs/non_malignant_mp_correlations/`:

- mode 01 -> `01_cancer_mps_cross_only`
- mode 02 -> `03_cancer_states_cross_only`
- mode 03 -> `05_cancer_mps_vs_whole_celltypes`
- mode 04 -> `06_cancer_states_vs_whole_celltypes`

The exported validation table records availability, scRNA Spearman direction
and p-value, spatial pooled direction/FDR, and whether directions agree. A
spatially enriched edge validates recurrent proximity, not molecular signaling
direction. Ligand-receptor causality requires the separate expression/LR
evidence already produced by the scRNA workflow.

## Outputs

All critical outputs are persistent under
`ref_outs/visium_hd_outs/spatial_interactions/`:

- `intermediate/`: TME MP scores/gene sets and per-sample replot-ready spatial
  membership RDS files.
- `tables/`: all sample tests, pooled tests, significant subsets, recurrence,
  radius sensitivity, source/target eligibility, map selection, parameters,
  and scRNA validation.
- `figures/Auto_spatial_interaction_dotmap_with_neg.pdf`: interaction dot-map
  pages with the three independently tested samples side by side.
- `figures/Auto_spatial_interaction_dotmap_with_neg_pooled.pdf`: pooled
  three-sample dot maps.
- `figures/Auto_0*_spatial_maps.pdf`: tissue maps for significant interactions,
  or explicitly labelled exploratory best-ranked pairs when a mode/sample has
  no FDR hit.
- `figures/Auto_spatial_interaction_recurrence.pdf`: recurrent/pooled summary.
- `figures/Auto_spatial_interaction_neighbourhood_sensitivity.pdf`: 16 versus
  32/48 µm robustness.
- `figures/Auto_spatial_interaction_scrna_validation.pdf`: spatial-versus-scRNA
  direction concordance.
- `logs/`: run parameters, cache status, timestamps, and session information.

A compact result summary is also written to
`updates/new_updates/summaries/visium_hd_spatial_cancer_tme_interactions_summary.csv`.

Set `SCREF_FORCE_REBUILD=TRUE` to ignore persistent memberships/tests. Set
`SCREF_REPLOT_ONLY=TRUE` to rebuild figures from the live memberships and
tables without rereading the large Space Ranger matrices.
####################
