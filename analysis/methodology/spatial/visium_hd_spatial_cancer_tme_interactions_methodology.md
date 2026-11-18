# Visium HD cancer–TME spatial association analysis

## Current 16-um-bin production contract

The user selected the previous verified 16 µm bin analysis as the current HD
representation. Its saved section tests, pooled tables, figures and reports
are restored from `spatial_interactions/reference/binned/` by copy. No previous
P values, null draws or significance thresholds are replaced. The restoration
checks the section-test SHA-256. Current HD/Xenium/meta evidence is regenerated
from these original section tests and the unchanged Xenium tests.

The maintained R producer again reads binned malignancy under `malignancy/`,
binned states under `state_mapping/`, and the exact `square_016um` H5 matrices.
It parses native 16 µm grid barcodes and links offsets up to one, two or three
grid steps, excluding the same position. It does not use segmented expression
or centroid geometry. All remaining scoring, eligibility and permutation
settings below describe the current binned estimand. The five-page comparison
retains the user-requested equal specimen weighting and explicit independence
assumption; original donor-aware tables remain separately labelled.

`SCREF_FORCE_REBUILD=TRUE,SCREF_MEMBERSHIPS_ONLY=TRUE` reconstructs binned UCell
and membership caches without repeating spatial permutations. It verifies the
original input fingerprint and writes
`reports/Auto_binned_cache_restoration_validation.csv`. The normal rebuild and
replot entries remain available through the matching existing PBS wrapper.

Segmented HD results, caches and frozen producer source are preserved by copy
under `spatial_interactions/reference/segmented/`. Its upstream segmented
malignancy/state outputs remain live. They do not supply current HD inference.
Cross-platform and directional segmented result references are retained under
their respective `reference/segmented_hd/` tiers. No files are deleted or moved.
Restoration provenance is recorded in `logs/Auto_binned_restoration.json`.

The four ranked interaction workbooks are maintained by
`spatial_platform_neighbourhoods.py` under `spatial_cross_platform/tables/`;
see the shared methodology for dataset-specific significance and sheet layout.

Bins can contain multiple cells or mixed RNA. The existing permissive
postfilter does not establish physical singlets; adjacency is bin-level
neighbourhood enrichment, not validated cell contact or molecular signalling.

## Purpose and scope

This terminal workflow asks whether cancer–TME interactions identified by
`analysis/non_malignant_nmf/mp_cross_celltype_correlations.R` recur as local
spatial enrichment in Visium HD. The executable is
`analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R`; the PBS wrapper
is `job_scripts/spatial/visium_hd_spatial_cancer_tme_interactions.sh`.

The analysis uses filtered 16 µm bins from all three supplied sections. Their
relationship and treatment status were confirmed by the dataset owner:

| Section | Donor | Treatment status |
|---|---|---|
| FFPED1 | Shared FFPE donor | Pretreatment |
| FFPEA1 | Shared FFPE donor | Post-treatment |
| SUR1231 | Separate SUR1231 donor | Post-treatment |

Sections retain separate graphs and permutation tests. FFPED1 and FFPEA1 are
repeated observations of one donor, not independent biological replicates.
All three sections contribute to the primary pool through two independent
donors. Both donors are independent of the six Xenium donors used downstream.
One paired donor cannot support a population-level treatment-effect analysis.

Four modes are retained:

1. malignant MPs versus non-malignant MPs;
2. malignant states versus non-malignant MPs;
3. malignant MPs versus whole non-epithelial cell types;
4. malignant states versus whole non-epithelial cell types.

The scRNA-selected edges form a separate candidate family. All other spatial
pairs are exploratory. Correlation between scRNA-seq samples and proximity
within tissue are different estimands. Agreement provides spatial support for
a candidate, not reproduction of its scRNA correlation coefficient or proof
of cell–cell communication.

## Current malignant MP and state contract

The script sources `analysis/shared/scRef_config.R` and reads:

- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`;
- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`;
- `ref_outs/visium_hd_outs/state_mapping/tables/Auto_visiumhd_binned_malignant_state_annotations.csv.gz`.

The grouping table must contain exactly the current 17 MPs and must match the
MP gene object exactly. `MP2x`, `MP11c`, and `MP18a` are explicitly rejected if
they reappear. MP names, descriptions, exact plotting order, state groups, and
Approach-B labels come from the current shared configuration/grouping table;
they are not redefined in this workflow.

The spatial state table was produced with the same five MP group memberships
but retains three earlier display labels. The following label-only crosswalk is
applied before analysis:

- `Basal to intestinal metaplasia` → `Squamous-to-intestinal`;
- `SMG to intestinal metaplasia` → `Glandular-to-intestinal`;
- `Stress adaptive` → `Stress-adaptive`.

`Classic proliferation` and `Cancer-cell immune mimicry` are unchanged.
`Hybrid` and `Unresolved` remain technical assignments and are not biological
source nodes. An unknown state label stops the run.

Malignant source bins must occur in the mapped state table and the final
post-filter annotation table. The mapped table is already restricted to
epithelial bins with `Auto_malignancy` equal to `malignant_level_1` or
`malignant_level_2`.

## Non-malignant targets and MP filtering

Whole-cell-type targets are final `Auto_postfilter_celltype` calls among bins
with `Auto_postfilter_keep == TRUE`, excluding `epithelial` and `unresolved`.

Non-malignant MP definitions come from the current `MP_outs_default.rds`
objects for fibroblast, endothelial, macrophage, NK, plasma, CD4, and CD8.
Only MPs with GeneNMF silhouette at least zero are retained, and the first 100
ranked genes form each spatial signature. The script stops for a missing or
malformed MP object; it does not substitute another object.

Each MP is scored only in bins assigned to its parent cell type. CD4 and CD8
programs are both scored in final `t.cell` bins because the binned annotation
does not provide a validated CD4/CD8 split. These are therefore T-cell program
labels, not claims that a bin is a purified CD4 or CD8 cell.

## MP scoring and positivity

Both malignant and non-malignant MPs are rescored from raw Space Ranger counts
with UCell (`maxRank = 1500`), using all measured genes as the ranking background.
Malignant signatures also use the first 100 unique ranked genes. Signature
genes absent from the matrix are intersected out explicitly. Every scored
signature must retain at least one gene. Scores must be finite and align
exactly with bin identities. Rank normalization does not eliminate sparsity,
panel-coverage or annotation effects.

This makes the score definition compatible with the upstream single-cell
interaction workflow, which also uses UCell. The numerical UCell scale is
platform-specific, so a scRNA-derived numeric cutoff is not transferred to
Visium HD. Within each spatial sample and eligible parent population, an MP is
positive when its UCell score exceeds:

`mean UCell + 1 × SD(UCell)`

For malignant MPs the population is mapped malignant bins; for a TME MP it is
the relevant parent-cell-type bins. States remain assignment based. Raw UCell
scores, fitted cutoffs, positivity, and feature eligibility are persisted for
audit and replotting.

Minimum eligibility is 20 positive malignant source bins, 20 whole-cell-type
target bins, 20 parent bins before TME-MP scoring, 10 positive TME-MP bins, and
20 source neighbourhood slots. For TME-MP tests, slots count neighbours of
the relevant parent; for whole-type tests, they count all TME neighbours.
No lower-threshold fallback is used.

## Spatial neighbourhood

The primary graph links a malignant source to eligible TME targets at Chebyshev
grid distance at most one, excluding the same grid position. Ring one includes
up to eight neighbouring 16 µm bins, with centre offsets up to ±16 µm per axis.
Rings two and three are cumulative windows at ±32 and ±48 µm. Maximum Euclidean
centre separations are approximately 22.63, 45.25 and 67.88 µm.

The graph does not force a neighbour count or interpolate missing bins. It is
a proximity graph: diagonal or broader-window edges do not establish a
continuous tissue path or membrane contact. Ring one is primary; broader rings
are separate sensitivity analyses, not additional opportunities to choose a
significant radius.

## Local permutation inference

The observed statistic is the directed count of source-positive to
target-positive adjacency edges. Source labels, the tissue mask, coordinates,
and source degrees remain fixed.

Targets are permuted within 20 × 20 grid-bin tiles
(320 × 320 µm) inside each sample. Four tile origins—`(0,0)`, `(10,0)`,
`(0,10)`, and `(10,10)` bins—define four separate conditional randomisation
tests. Each origin retains its own null counts and empirical p-value.

- Whole-cell-type labels are shuffled jointly within each local tile. This
  preserves tile-specific cell-type counts.
- The complete TME-MP positivity matrix is shuffled jointly within parent cell
  type and local tile. This preserves parent-cell-type location, each MP count,
  and within-bin MP overlap while testing program-specific placement.

At least half of target bins must be movable under every shifted tiling or the
run stops. The default and minimum are 4,999 permutations per origin, or 19,996
draws in total for each section, target kind and ring. Eight fixed seeded
chunks make the draws independent of the worker count. Each null assumes
exchangeability of target memberships within its strata. It preserves broad
local composition, but not arbitrary finer spatial autocorrelation.

With observed edge count `O`, permutation counts `O_b` and `B` permutations,
the two-sided empirical p-value is
`min(1, 2 * min((1 + count(O_b >= O))/(B + 1), (1 + count(O_b <= O))/(B + 1)))`
within each origin. The section p-value is the maximum of the four origin
p-values. Directions must also be nonzero and agree relative to every
origin-specific null mean; otherwise the section p-value is one. This
conservative conjunction requires evidence across all four boundary choices.
It does not assume independence among origins. The tables retain every
origin's null mean, upper/lower tail counts and p-value.

Permutations from different origin-specific groups are not pooled to compute a
test p-value: such a union need not itself be a permutation group. This follows
the group-structure requirement for random permutation tests
([Hemerik and Goeman, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6405018/)).
Exchangeability within the chosen spatial strata remains a scientific assumption.

The descriptive effect is `log2((O + 0.5)/(mean_expected + 0.5))`, where
`mean_expected` averages the four origin-specific expected counts. Zero denotes the local null
expectation; the 0.5 continuity correction is fixed. BH correction is applied
within section, mode and ring.
Normal-tail p-values derived from permutation z-scores are retained only as
diagnostics and never drive FDR calls.

## Donor-level evidence

Within each interaction and ring, the descriptive donor effect is the arithmetic
mean of eligible section log2 enrichments. One eligible section contributes
its empirical p-value and effect direction. For `m` eligible sections of one
donor, concordant nonzero directions give
`p_donor = min(1, m * min(p_section))`. Discordant or zero section directions
give `p_donor = 1` and signed evidence zero. The donor remains represented;
its descriptive effects are retained. Ineligible sections remain explicitly
unavailable and do not receive substitute evidence.

The Bonferroni minimum-p rule is conservative under arbitrary within-donor
dependence when the joint section null holds. It tests whether a donor contains
section-level evidence, not whether every section is significant or whether
the mean effect differs from zero. Direction consistency is an additional
conservative requirement. This retains both paired sections without counting
them as independent donors.

For each donor, `z_d = sign_d * normal_upper_quantile(p_donor/2)`. Equal-donor
Stouffer pooling uses `Z = sum(z_d)/sqrt(n_donors)` and a two-sided normal-tail
p-value. At least two eligible donors are required. BH correction is within
mode and ring. The mean donor effect and the separately reported
contact-weighted effect are descriptive; neither is a variance-weighted effect
estimate or an effect confidence interval.

An exploratory pooled call requires FDR below 0.05 and concordant nonzero
direction across eligible donors. Donor recurrence remains descriptive and
cannot turn repeated sections into independent replication. With two donors,
results provide limited evidence of reproducibility, not population prevalence.

## Comparison with scATLAS candidates

The exact `Auto_celltype_correlations_all.csv` tables from the four matching
upstream modes are required. Missing inputs stop the run. Earlier description
aliases are converted to current MP IDs before joining, preventing a display
label change from breaking an MP match.

The current Visium HD annotation resolves `t.cell`, not CD4 and CD8 separately.
For whole-cell-type validation only, upstream CD4 and CD8 candidates are tested
against the same spatial `t.cell` compartment; upstream NK is tested against
spatial `nk.cell`. The candidate-availability table retains both the original
scRNA compartment and mapped spatial compartment. This crosswalk is not used
for TME-MP targets, where CD4 and CD8 programs remain distinct within their
shared spatial t-cell parent compartment.

An upstream edge is prespecified as follows, matching the upstream selection:

- malignant MP–TME MP: Spearman p below 0.05, finite Pearson effect, and
  `-log10(Spearman p) >= 4`;
- malignant state–TME MP: the same rule with threshold `1.5`;
- whole-cell-type modes: Spearman p below 0.05 and a finite Pearson effect.

The finite Pearson value is retained only as an upstream availability check.
Candidate direction is defined consistently by the upstream Spearman
correlation, which is also the effect plotted on the validation x-axis.

Within mode and ring, candidate-family BH uses all unique selected spatial
hypotheses. CD4/CD8 whole-type rows that map to the same source–`t.cell` test
contribute one hypothesis. Original candidate rows remain available for
interpretation, but neither establishes subtype-specific validation.
Unavailable selected hypotheses remain in the multiplicity denominator and
retain missing test results in the outputs.

Spatial support requires candidate-family FDR below 0.05, at least two eligible
donors, consistent donor directions and agreement with the upstream Spearman
direction. Other results distinguish concordant evidence without support,
discordant evidence, spatially unavailable hypotheses and non-selected pairs.
Raw upstream selection p-values do not establish that the upstream study
controlled FDR. Candidate-family correction concerns the spatial tests.

## Figure interpretation

Figures distinguish effect magnitude, availability and significance. Blank or
explicitly unavailable entries are not zero effects. Current MP/state labels
and a common zero-centred enrichment scale make comparisons interpretable.
Individual section evidence remains visible alongside donor pooling.

Whole-type and within-parent MP effects answer different questions. Whole-
fibroblast depletion means fewer fibroblast neighbours than expected from local
non-epithelial composition. A positive fibroblast-MP effect means enrichment
of that programme relative to locally available fibroblast bins. Their signs
can differ without contradiction.

Main figures use editable vector text and are assessed at their intended
publication size. [Nature's artwork guide](https://www.nature.com/documents/nature-final-artwork.pdf)
specifies 183 mm for two-column figures and 5–7 pt ordinary text. The layout
aims for 7 pt labels at that size without adding pages merely to enlarge text.
Spatial maps retain a physical scale bar. Source tables remain live.

## Outputs and cache behavior

Persistent outputs are under
`ref_outs/visium_hd_outs/spatial_interactions/`:

- `intermediate/`: current-version per-sample memberships, filtered gene sets,
  and malignant/TME UCell scores with cutoffs;
- `tables/`: section, donor, pooled, recurrence, eligibility, candidate,
  map-summary, parameter, figure-source, paired-section sensitivity and compact-summary tables;
- `figures/Auto_final/01_scatlas_spatial_evidence.pdf`: four pages, one comparison
  per page, showing every selected scATLAS candidate and individual section effects;
- `figures/Auto_final/02_supported_tissue_context.pdf`: tissue maps for supported
  candidates, with all three sections on each page and a shared physical scale;
- `figures/Auto_supporting/Auto_primary_all_pairs.pdf`: complete primary effects
  and donor-concordant FDR calls, one page per comparison;
- `figures/Auto_diagnostics/Auto_window_sensitivity.pdf`: broader-window effects
  against the primary window, with all four comparisons on one page;
- `reports/CURRENT_RESULTS.md`: the current results, limitations and figure index;
- `logs/`: run summary and R session information.

The current index identifies the selected figure set. Production reports and
downstream code use these current paths.

The plotting workflow also writes
`updates/new_updates/summaries/visium_hd_spatial_cancer_tme_interactions_summary.csv`.

`SCREF_FORCE_REBUILD=TRUE` rebuilds scores, memberships, tests, tables, and
figures. `SCREF_REPLOT_ONLY=TRUE` requires current-version persistent
memberships, scores and section result tables. It regenerates donor inference,
candidate summaries and figures without repeating UCell or permutations.
Missing or incompatible inputs stop immediately.

## Limitations and validation checks

Visium HD bins can contain mixed RNA and do not prove direct ligand–receptor
communication. Parent-cell-type errors can propagate into TME-MP scores.
Local tiling controls broad regional composition but cannot guarantee complete
removal of fine-scale spatial autocorrelation. Ring sensitivity, empirical
section tests, equal donor weighting, direction recurrence and the explicit
candidate family limit overinterpretation. One paired donor does not establish
a treatment effect.

Recent benchmarking emphasizes that edge enrichment and the fraction of source
cells involved describe different tissue properties. The current statistic is
edge enrichment, not a contact probability or directed signalling score
([Schiller et al., 2026](https://doi.org/10.1038/s41467-026-71699-z)).

Every production run must verify:

1. exact 17-MP agreement between grouping and gene objects;
2. successful UCell scoring and non-empty signatures;
3. eligibility and skip-reason tables, especially T-cell/NK coverage;
4. empirical rather than normal-tail FDR columns;
5. paired-donor identity, candidate availability and donor-level direction concordance;
6. all three ring outputs and current cache version;
7. canonical malignant MP order in every comparison figure;
8. correct candidate-family multiplicity, figure-source data and physical map scale bars;
9. compact summary, run log, session information, and readable PDFs.

## Run

Submit from the live project root:

    /opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE job_scripts/spatial/visium_hd_spatial_cancer_tme_interactions.sh

The wrapper uses eight CPUs, 160 GB RAM, a 24-hour walltime, `#PBS -koed`, and
the required `dmtcp` environment. PBS stdout/stderr go directly to the project
root. The job must be monitored to terminal PBS
status and the expected persistent outputs must be checked after completion.

## Shared comparison report

Pages 2–5 now call the original `plot_focal_interaction_dotmap()` in the
maintained scATLAS R producer, preserving its exact palette, physical size
range, fonts, margins, grid and facet layout. Spatial effect and FDR replace
the original correlation and co-positive measurements; no statistical values
are relabelled as different quantities. See the shared methodology for the
live plot spec, source tables and current PDF assembly components.

The requested `figures/Auto_spatial_interaction_four_mode_significant_comparison.pdf` is regenerated from specimen-pooled comparison tables, treating all three sections separately with equal weights under the explicitly requested specimen-independence assumption. The primary donor-aware tables produced earlier in this R workflow remain separately available. It has a common four-panel overview followed by four tiled detail pages, matching the scATLAS detail-page dimensions and shared cell-type order. The HD producer calls the maintained Python renderer after scientific exports. Selection is primary-window specimen-pooled BH FDR <0.05 with at least two eligible specimens. This does not establish biological independence of A1/D1; their shared patient identity is retained. Rendering methods, other platform criteria and exact execution are documented in [the shared figure methodology](xenium_visiumhd_meta_methodology.md#shared-five-page-comparison-reports).
