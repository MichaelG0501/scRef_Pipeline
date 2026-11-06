# Xenium and Visium HD contact recomputation and meta-analysis

## Scope and execution

This workflow measures local cancer–TME spatial enrichment in six Xenium donors
and combines it with Visium HD evidence from three sections of two donors.
All 17 current malignant MPs and the five biological states are included.
Ordinary Visium does not contribute interaction statistics.

The paper's 16 Xenium tissue samples are represented within six released patient
coordinate fields: P4_A/B, P8, P10_A/B/C, P11_A/B/C, P12_A/B/C/D and P13_A/B/C.
The authors' visualization notebook extracts these regions with coordinate
windows. All 16 regions contain cells in the local release, and all 558,318
released cell identities match the current pipeline metadata. Thus six denotes
independent donors, not a selection of six of the 16 tissue samples. Current
scoring and contact inference operate within patient fields; the release has no
separate tissue-region column. The visualization rectangles can overlap and
must not be treated as a validated disjoint tissue segmentation. Coverage and
provenance are recorded in the Xenium output's
`reports/Auto_paper_sample_reconciliation.md` and
`reports/Auto_paper_xenium_region_coverage.csv`.

`analysis/spatial/Auto_xenium_neighbourhoods.sh` is the essential compute-node
orchestrator. It prepares metadata with `spatial_platform_neighbourhoods.py`,
scores raw counts with `xenium_hd_ucell_scores.R`, then runs contact inference,
meta-analysis and figures in Python. Submit from the project root with:

```bash
/opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE job_scripts/spatial/spatial_platform_neighbourhoods.sh
```

The launcher requests eight CPUs, 128 GB RAM and 24 hours. It uses `#PBS -koed`
and writes scheduler stdout/stderr directly in the project root. R UCell runs
in `dmtcp`; Python runs in the existing `bidcell_temp` environment. Scientific
computation runs only on allocated compute nodes. No script depends on the
reconstructable local launcher or on a temporary audit script.

## Exact inputs and identities

- `ref_outs/Yates_spatial/source_data/Zenodo_upload/processed/Xenium.h5ad`:
  raw CSR `layers/counts`, 4,999 panel genes, cell identities, patient,
  original cell types and SOPA image-pixel spatial coordinates. The authors'
  SOPA aggregation exports polygon centroids in image coordinates, without a
  pixel-to-micrometre conversion.
- `ref_outs/Yates_spatial/source_data/Zenodo_upload/raw/Xenium/<patient>/cell_boundaries.parquet`:
  raw XOA boundary vertices used to check the declared pixel calibration.
- `ref_outs/xenium_scatlas_states/Auto_xenium_cell_annotations.csv.gz`:
  existing malignant Approach B state assignments. The H5AD observation index
  is matched exactly; patient and cell type must also agree. The published
  `patient` field supplies explicit sample identity and is represented as
  `orig.ident`; identity is never inferred from a trimmed barcode.
- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
  and `ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`:
  exact current 17-MP genes, labels and order. MP1, MP5 and MP13+ are included as
  activity sources, matching HD. MP2x, MP11c and MP18a remain excluded.
- `ref_outs/visium_hd_outs/spatial_interactions/intermediate/Auto_visiumhd_tme_mp_gene_sets.rds`:
  exact saved, silhouette-filtered, top-100 HD TME signatures and parent labels.
- `ref_outs/visium_hd_outs/spatial_interactions/tables/Auto_spatial_interactions_by_sample.csv.gz`
  and `Auto_spatial_interaction_parameters.csv`: existing HD per-section
  results, all four modes and three rings; HD pooled results are not inputs.

The count layer must contain finite non-negative integer values. Its CSR
arrays are represented as a genes-by-cells CSC matrix in R. The Matrix
constructor sorts the source CSR's unsorted gene indices without changing
counts or cell order. Raw count dimensions and metadata order must agree.

## Approved target matching

| Xenium labels | Shared parent |
|---|---|
| Fibroblast, Inflammatory CAF, Adipose CAF | fibroblast |
| Endothelial | endothelial |
| TAM1, TAM2 | macrophage |
| T, Treg | t.cell |
| B | b.cell |
| DC | dendritic |
| Mast | mast |

Adipocyte, Muscle and Nerve/adrenal retain explicit Xenium-only labels. They
remain nuisance target categories in local whole-type permutations, but cannot
enter cross-platform tests without an HD counterpart. Carcinoma is the source
population. Epithelial and Epithelial (Pre-cancerous) are excluded from TME
targets, matching HD's epithelial exclusion. Unrecognised labels stop the run.
HD NK, plasma and other unmatched classes remain HD-only. The approved B-cell
mapping does not create a separate Xenium plasma compartment.

## UCell activity and state sources

Use R UCell in `dmtcp`, with the same installed version used by the HD workflow,
`maxRank=1500`, raw counts, first 100 ranked genes per malignant MP, and exact
saved HD TME signatures. Scoring ranks all 4,999 genes in the count matrix;
only signature genes missing from the panel are removed, as in HD's signature
intersection. Coverage is audited for every signature. A tested parent with a
signature having no panel genes stops the run; there is no alternative scorer.

All 17 malignant MPs are scored only within Carcinoma cells. Non-malignant MPs
are scored only within their matched parent population. CD4 and CD8 MPs are
both scored in the combined T + Treg population. These are separate programme
features, not mutually exclusive CD4/CD8 cell identities. Co-positive cells
are allowed. No NK or plasma labels are inferred from programme scores.

An MP is positive when its score is strictly greater than the patient/parent
mean plus one sample standard deviation. Minimums match HD: 20 source-positive
observations; 20 whole-type target observations; 20 parent observations for TME
scoring; 10 TME-MP-positive observations; and 20 source neighbourhood slots per
tested pair. For TME programmes, these slots must contact the programme's
parent population; whole-cell-type tests use all non-epithelial TME slots. Audits retain unavailable and ineligible features.

The five current malignant state assignments are retained from the existing
mapper, as HD also uses its upstream state mapping. UCell activity is not used
to redefine states. Hybrid and Unresolved are excluded as biological sources.

## Spatial graph and null model

Convert the authors' SOPA centroid coordinates from image pixels to micrometres
using 0.2125 micrometres per pixel. The authors' preprocessing is documented in
[their SOPA notebook](https://github.com/vanallenlab/EAC-multiome/blob/main/code/python/notebooks/spatial-transcriptomics/4.%20SOPA-xenium-cellpose.ipynb);
[SOPA exports polygon centroids](https://github.com/gustaveroussy/sopa/blob/master/sopa/aggregation/table.py).
The physical pixel scale and XOA boundary vertex units are documented in the
[10x Xenium format specification](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html).
For each supplied sample, verify 4,096 raw XOA boundary coordinate values:
at least 95% must align within 0.05 pixels of the 0.2125-micrometre pixel grid.
Failure stops preparation. Persist raw pixel centroids, converted coordinates,
scale, units and a per-sample calibration audit. The check corroborates the
pixel-grid scale. The centroid coordinate convention is inferred from the
authors' image-based Cellpose/SOPA workflow and SOPA's direct polygon-centroid
export; boundary-grid alignment alone does not establish the H5AD units. Ranks and UCell scores do not depend on spatial units.

Use the exact physical equivalent of the HD queen-window condition:
`max(abs(source_x-target_x), abs(source_y-target_y)) <= 16*ring` micrometres.
Construct sparse malignant-to-TME adjacency directly with a Chebyshev KD-tree
query on calibrated cell centroids; do not bin Xenium cells. Ring one is the
primary ±16-micrometre square window. Rings two and three are cumulative ±32
and ±48-micrometre square windows, matching the HD grid ring footprints. Their
maximum Euclidean separations are 22.63, 45.25 and 67.88 micrometres. Duplicate
coordinates stop the run. Sources and targets are disjoint, so no self-edges
occur. The adjacency is binary, with no repeated counting of the same pair.

The centre-distance condition and axis-aligned physical footprints match HD;
cell and bin densities and biological measurement units still differ. A shared window measures local proximity, not membrane
contact or demonstrated signalling.

Sources, degrees, coordinates and target positions remain fixed. Whole-type
membership rows are jointly shuffled within each patient and shifted
320-micrometre coordinate blocks. TME-MP membership rows are jointly shuffled
within parent type and those same blocks, preserving within-cell MP overlap.
Four block origins (0,0), (160,0), (0,160), (160,160) define separate permutation
tests. Each null assumes exchangeability of target memberships inside its strata.
It preserves broad local composition and within-cell MP overlap, but not
arbitrary fine-scale spatial autocorrelation. At least 50% of target
observations must be movable in every origin, or the run stops. The default
and minimum are 4,999 permutations per origin, matching the HD run: 19,996
draws in total for each patient, target kind and ring.
Eight fixed, independently seeded chunks are scheduled over the requested
workers, so changing the worker count does not change permutation draws.

Observed contacts are source-membership transpose times adjacency times
target-membership. Null means and SDs, tail counts and two-sided add-one
empirical p-values are calculated from permutations. With observed count `O`,
permuted counts `O_b` and `B` permutations, the p-value is
`min(1, 2 * min((1 + count(O_b >= O))/(B + 1), (1 + count(O_b <= O))/(B + 1)))`
within each origin. The section p-value is the maximum across origins, with
matching nonzero effect directions relative to all four origin-specific null
means required; otherwise it is one. This conservative conjunction requires
support across boundary choices without assuming origin independence. Every
origin's expected count, upper/lower tail counts and p-value is persisted.

The origin-specific shuffles each form a permutation group. Mixing their draws
into one p-value would not generally preserve that group structure
([Hemerik and Goeman, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6405018/)).
The effect uses `log2((observed+0.5)/(null_mean+0.5))`, with `null_mean` averaging
the four origin-specific expected counts. It is descriptive, rather than a
test based on a mixture of null draws. Ineligible pairs receive missing p-values.
BH adjustment is within patient, mode and ring, exactly as in HD. Four modes
are malignant MPs/states against TME MPs/whole cell types.

## Biological unit and donor-level meta-analysis

The dataset owner confirmed that FFPED1 is pretreatment and FFPEA1 is
post-treatment from the same donor. SUR1231 is post-treatment from another
donor. All three HD sections contribute to the primary pool, but the paired
sections do not count as independent replicates. The six Xenium patients,
P4, P8, P10, P11, P12 and P13, are distinct donors and independent of both HD
donors. Their treatment statuses are not assumed from the HD metadata.
The primary pool therefore has at most eight independent donors.

Use eligible empirical p-values from Xenium and the HD per-section table.
Platform-qualified section and donor IDs prevent collisions. Duplicate
interaction/section keys or unexpected identities stop the run. Missing
compartments remain unavailable, with missing effects and test results.
They are never assigned a zero effect or p=1 merely to fill a cross-platform
comparison.

Within each mode, ring and source–target pair, first collapse eligible sections
of the same donor. The descriptive effect is their mean log2 enrichment.
One eligible section contributes its empirical p-value and effect direction.
For `m` eligible sections, concordant nonzero effect directions give
`p_donor = min(1, m * min(p_section))`. Discordant or zero section directions
give `p_donor = 1` and signed evidence zero. The donor and its descriptive
section effects remain in the analysis.

This Bonferroni minimum-p rule is conservative under arbitrary within-donor
dependence for the joint section null. It tests whether the donor contains
section-level evidence; it is neither a test of the mean effect nor a claim
that every section is independently significant. The direction filter is an
additional conservative restriction. Both paired sections affect that rule
and the descriptive donor mean.

For independent donors, compute
`z_d = sign_d * normal_upper_quantile(p_donor/2)` and
`Z = sum(z_d)/sqrt(n_donors)`. Convert `Z` to a two-sided normal-tail p-value
and apply BH within mode and ring. Donors receive equal weight, regardless
of cell/bin count or number of sections. The empirical permutation resolution
is retained. A pooled test requires at least two eligible independent donors.
Ring one is primary; rings two and three are separate sensitivity families.

All eligible shared pairs are exploratory hypotheses. This pool is separate
from the HD scATLAS-candidate family. Mean donor log2 enrichment is descriptive. Because the platforms measure
cells and bins with different gene coverage, neither is a common physical
contact effect estimate or an effect confidence interval.

## Cross-platform replication and sensitivity

A pooled FDR below 0.05 establishes combined evidence, not replication in both
technologies. A single platform can drive pooled significance even when the
other has weak or opposing evidence. Therefore the workflow reports separate
pooled and replicated result sets.

For a replication test, both platforms must contain at least two eligible
independent donors. Compute each platform's donor-level signed-Stouffer p-value.
Require matching nonzero platform Z directions and agreement between each
platform Z and its descriptive mean-effect direction. When these agree, set
`p_replication = max(p_HD, p_Xenium)`; discordant directions receive
`p_replication = 1`. The maximum-p conjunction tests the null that at least
one platform lacks evidence. It is BH-adjusted within mode and ring. A
replicated call requires replication FDR below 0.05 and direction agreement.
If either platform lacks enough donors, replication is unavailable, not failed.

The workflow also reports within-platform BH results, equal-platform-weight
Stouffer sensitivity, and leave-one-donor-out sensitivity. Equal-platform
weighting combines the two normalized platform Z statistics with equal weights;
it does not change the primary equal-donor analysis. Leave-one-donor-out drops
all sections of the omitted donor, never only one repeated section. A result
with fewer than two donors remaining in a platform cannot establish replication
in that sensitivity analysis.

Two paired-section substitution analyses retain SUR1231 and every Xenium donor
while representing the FFPE donor by FFPED1 alone or FFPEA1 alone. A separate
paired difference table reports post-minus-pre section effect differences.
These are descriptive sensitivity checks from one paired donor. They do not
estimate a population treatment effect, and pooling all sections does not
remove possible treatment or cohort differences.

## Persistent outputs, cache and plotting

Xenium output root: `ref_outs/xenium_scatlas_states/spatial_neighbourhoods/`.
`intermediate/` contains identity-aligned metadata, per-patient UCell scores,
source/target membership matrices and spatial graph/adjacency arrays.
`tables/` contains crosswalk counts, signature coverage, eligibility audits,
per-patient contacts, combined contacts and graph diagnostics. `logs/` contains
R session information and exact contact parameters. All are canonical live
outputs needed for audit, replot or contact-only reproduction.

Meta output root: `ref_outs/spatial_cross_platform/`. Current tables are:

- `Auto_sample_donor_design.csv` and `Auto_target_label_crosswalk.csv`;
- `Auto_harmonised_interactions_by_sample.csv.gz` and
  `Auto_harmonised_interactions_by_donor.csv.gz`;
- `Auto_cross_platform_availability.csv`, `Auto_within_platform_meta.csv`
  and `Auto_cross_platform_meta.csv`;
- `Auto_leave_one_donor_out.csv.gz`, `Auto_paired_hd_section_differences.csv`
  and `Auto_paired_section_substitution.csv.gz`;
- `Auto_cross_platform_pooled_primary.csv` and
  `Auto_cross_platform_replicated_primary.csv`.

`figures/Auto_final/01_cross_platform_contact_enrichment.pdf` shows the four
source–target modes, one intended comparison per page. Each malignant MP label
includes its full current description from the centred refinement grouping
table. TME programmes retain their compartment-specific identifiers; the
malignant MP descriptions must not be assigned to identically numbered TME MPs.
Rows without any pooled estimate are omitted, including unmatched NK and
plasma targets. This display filter does not change the testing families or
the complete availability tables. Partially unavailable rows remain, with
missing cells shown in grey. Colour represents descriptive mean donor log2
enrichment, clipped at ±2. A side legend distinguishes filled circles
(replication FDR < 0.05) from open circles (pooled FDR < 0.05 without replication
FDR < 0.05). The figure contains comparison titles and no explanatory subtitles
or footnotes. HD scATLAS-candidate support belongs to a separate testing family
and is displayed in the HD final evidence figure, linked from the current
results index. Additional
section/donor evidence is in `figures/Auto_supporting/`; coverage and eligibility
are in `figures/Auto_diagnostics/`. Individual-effect plots do not fabricate
confidence intervals. Text remains editable and layouts are checked at 183 mm
width, following [Nature's artwork guidance](https://www.nature.com/documents/nature-final-artwork.pdf).

The current index is `reports/Auto_current_results.md`. The current scientific
report and compact result table are `reports/Auto_cross_platform_run_report.txt`
and `reports/Auto_cross_platform_result_summary.csv`. `logs/` records input
provenance. The compact update is
`updates/new_updates/summaries/Auto_xenium_visiumhd_meta_summary.csv`.
Xenium-only donor pooling is saved in
`ref_outs/xenium_scatlas_states/spatial_neighbourhoods/tables/Auto_xenium_interactions_pooled.csv`;
its selected figure is
`figures/Auto_final/01_xenium_contact_enrichment.pdf` relative to that root.

The orchestrator requires exactly one of `SCREF_FORCE_REBUILD`,
`SCREF_CONTACT_ONLY`, `SCREF_META_ONLY` or `SCREF_REPLOT_ONLY` to be `TRUE`. Full PBS runs
recompute metadata, UCell scores, contacts and pooling. `SCREF_REPLOT_ONLY=TRUE`
regenerates figures from persistent meta/sample tables without rescoring.
For a scoped update of the four-page cross-platform final PDF and its current
index, set `SCREF_FIGURE_SCOPE=cross-platform-final` together with
`SCREF_REPLOT_ONLY=TRUE`. The default figure scope is `all`. The corresponding
Python option is `--stage replot --figure-scope cross-platform-final`; a
restricted scope is rejected for analytical stages. A scoped replot needs
one CPU, 8 GB RAM and 15 minutes:

```bash
/opt/pbs/bin/qsub -l select=1:ncpus=1:mem=8gb -l walltime=00:15:00 -v SCREF_REPLOT_ONLY=TRUE,SCREF_FIGURE_SCOPE=cross-platform-final job_scripts/spatial/spatial_platform_neighbourhoods.sh
```

The Python `--stage contacts` command is an explicit contact-only rerun using
persisted identity-aligned scores; it is not an automatic stale-cache fallback.
`SCREF_CONTACT_ONLY=TRUE` in the PBS wrapper repeats validation, metadata
preparation/calibration and contacts using persisted audited UCell scores. It
stops after contact inference so the two platforms can finish independently.
`SCREF_META_ONLY=TRUE` pools completed HD/Xenium section tests and exports figures.
Contact-only reruns are used for spatial corrections that do not change scoring inputs.
Score reuse stops if counts, annotations, MP definitions, grouping or the R
scorer are newer than any persisted score or feature-audit table. Replotting
checks the inference version and recorded input sizes/modification times;
editing plot code alone does not invalidate the numerical inference contract.
Current meta-analysis rejects an uncalibrated or older-version contact table.

For a permutation-only method update, reuse the audited UCell scores with
`SCREF_CONTACT_ONLY=TRUE`. After both platform section-test jobs have completed,
submit the same launcher with `SCREF_META_ONLY=TRUE`. Use PBS `afterok`
dependencies on both producers so pooling cannot read partially updated inputs.
Full rebuilds still execute all stages in order within one job.

## Validation and limits

The PBS run checks BH adjustment, donor collapse, signed-Stouffer cancellation and amplification,
permutation invariance, merged-label contact conservation, identity joins,
count integrity, 17-MP completeness, coverage, positivity audit consistency,
null mobility, donor/section identities, replication conjunction and duplicate meta keys. Completion requires PBS exit 0 and all
expected tables, reports, figures and persistent score/membership files.

UCell values remain affected by targeted-panel gene selection and per-cell
sparsity. Panel overlap is not equivalent to detectable expression. Thresholds
are fitted independently by patient/parent; rare or absent populations remain
unavailable. Cell-centroid proximity and HD bin adjacency are different spatial
units despite their matched physical windows. Edge enrichment is not the
fraction of source cells that contact a target, nor evidence of signalling.
Recent benchmarking demonstrates that these spatial summaries can differ
([Schiller et al., 2026](https://doi.org/10.1038/s41467-026-71699-z)). Local
shuffles do not remove arbitrary finer spatial autocorrelation. Gene coverage,
cell annotation, treatment, tissue preparation and cohort differences remain
potential sources of heterogeneity. Only two HD donors are available.
Use the conjunction-based replication field and donor sensitivity tables for
claims of cross-platform reproducibility; pooled significance alone is
insufficient.
