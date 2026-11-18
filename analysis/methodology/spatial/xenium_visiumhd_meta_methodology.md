# Xenium and Visium HD contact recomputation and meta-analysis

## Current HD 16-um-bin inputs

The user restored the previous verified 16 µm binned HD analysis. Current
meta-analysis requires `bin_size_um=16` in the HD parameters and uses the
original saved section tests. Xenium retains its calibrated segmented-cell
centroids and unchanged section tests. The shared 16/32/48 µm physical windows
and conditional permutation settings are retained, but the measurement units
remain bins versus cells. Segmented HD results are preserved only as reference
outputs and do not supply current inference. The exact restoration/cache
contract is in the [HD methodology](visium_hd_spatial_cancer_tme_interactions_methodology.md#current-16-um-bin-production-contract).

## Scope and execution

### Plain Excel exports

Four workbooks under `ref_outs/spatial_cross_platform/tables/` contain one
worksheet per mode:

- `Auto_01_scATLAS_TME_interactions.xlsx`
- `Auto_02_VisiumHD_TME_interactions.xlsx`
- `Auto_03_Xenium_TME_interactions.xlsx`
- `Auto_04_pooled_TME_interactions.xlsx`

Sheets are plain header/data tables: no titles, explanatory rows, colours,
custom cell styles, ranks, filters or duplicate display flags. Basic column
widths keep the full malignant MP descriptions readable. Rows retain the
bubble-plot cell-type/programme order within the significant block and the
remaining block, separated by three blank rows. TME programme IDs are included
inside the TME cell-type value, so each interaction remains identifiable
without a separate programme column.

The single `Significant (bubble rule)` column matches the corresponding bubble
selection exactly. For scATLAS this is the original exported positive/negative
selection, using mode-specific nominal-P/display criteria. BH FDR remains a
numeric diagnostic; it does not redefine this plot flag. The original plot
criteria require finite nonzero Pearson direction and Spearman P <0.05, with
additional -log10(P) >=4 for MP/TME-MP modes and >=1.5 for state/TME-MP modes;
whole-cell-type modes have no additional display cutoff. The exported display
source is authoritative, including preserved original selections. This fixes
the earlier mismatch between a BH-based workbook flag and nominal-P bubbles,
without altering the plots or their statistics.

scATLAS exports only the flag, malignant feature, TME target, direction,
Spearman rho/P, BH FDR and N samples. Spatial exports retain mean log2
enrichment, pooled P/FDR, N samples and each individual specimen's log2
enrichment. Per-specimen P/FDR, patient counts, co-positive percentages and
study-name lists are omitted as requested. Pooled sheets also contain separate
HD and Xenium mean/P/FDR/sample counts and combined mean/P/FDR. Combined
statistics do not change the approved pooled flag, which still requires
native FDR <0.05 in both spatial platforms and matching nonzero means.

Every meta run updates the workbooks. `--stage workbooks` and the existing
launcher's `SCREF_WORKBOOKS_ONLY=TRUE` export only matching current tables;
input-provenance mismatches stop the export. Installed openpyxl is the sole
writer. Missing values remain blank. Complete source inputs and exact rules
are recorded in `logs/Auto_interaction_workbook_provenance.json`; sheet counts
are in `reports/Auto_interaction_workbook_summary.csv`, with a compact update
under `updates/new_updates/summaries/`. The external workbook index holds
method details, keeping the Excel files data-only. Independent checks compare
every exported value and significance flag with the original source tables
and figure selections, verify ordering and separators, and reject any styled
cells or extra content.

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

## Shared five-page comparison reports

The HD, Xenium and cross-platform comparison PDFs now use one shared renderer
in `spatial_platform_neighbourhoods.py`. Page 1 is a less compact 610 × 430 mm
four-panel overview, with a shared full MP/state row axis and wider plotting panels to prevent crowded target labels. Pages 2–5 are one mode each, in numeric mode order, at
exactly 22 × 10.5 inches (558.8 × 266.7 mm), matching page 1 of the scATLAS
`Auto_celltype_interaction_dotmap_with_neg.pdf`. Spatial detail pages call
`plot_focal_interaction_dotmap()` in the original maintained
`mp_cross_celltype_correlations.R` producer, through its
`SCREF_DOTMAP_RENDER_SPEC` plot-only entry. They preserve its exact colours,
22 × 10.5-inch canvas, fonts, tile/grid stroke, 48-degree labels, facet spacing,
strip fill, margins, legend layout, point alpha/stroke and 1.3–8.6-mm point-size
range. Only the measured quantities, scale labels and comparison title change. Spatial dot size
represents FDR evidence, rather than borrowing the scATLAS co-positive-sample
percentage; spatial fill represents mean log2 enrichment, not Spearman rho.
No subtitles or LR-support marks are invented for spatial data.

All plots use the upstream scATLAS parent order: fibroblast, endothelial, CD8,
CD4, macrophage, NK, plasma. A whole T-cell spatial target occupies the T-cell
position without implying a resolved CD4/CD8 subtype. Other spatial targets
follow those shared parents. MP numbers are naturally ordered within parent.
Detail pages show all source MPs/states and the available target universe,
including blank tiles for interactions failing the figure selection.

### Specimen-based comparison inference (user-requested update)

Each HD section is treated as a separate statistical specimen: FFPED1, FFPEA1
and SUR1231 receive equal weights, matching the per-specimen treatment of the
six Xenium samples. Each retains its independently calculated local-permutation
P value and descriptive log2 enrichment. The three HD section graphs and null
tests are not recomputed. At least two eligible specimens are required for a
pooled test. Convert each section P to signed normal evidence and combine
`Z = sum(sign(effect) * normal_upper_quantile(P/2)) / sqrt(n_specimens)`.
Apply a two-sided normal tail and BH over every mode/window's finite tests.
Mean effects are arithmetic means of eligible specimen effects.

This calculation explicitly **assumes statistical independence of specimen
tests** for Stouffer pooling. The user requested this analysis after being told
that A1/D1 share a patient and that different treatment conditions do not prove
independence. Therefore it is specimen-based, conditional spatial evidence;
it does not turn three specimens into three independent biological patients.
Actual patient IDs are preserved. No estimated paired correction or mixed model
is applied in this comparison calculation. The original donor-aware scientific
tables remain available for sensitivity and do not supply the displayed means
or within-HD FDR in these reports.

Current comparison inputs/results are
`tables/Auto_specimen_within_platform_meta.csv` and
`tables/Auto_specimen_cross_platform_meta.csv` under spatial_cross_platform,
plus `tables/Auto_specimen_pooled_interactions.csv` under the HD root.
`logs/Auto_specimen_pooling_provenance.json` records the raw section inputs,
unit and independence assumption. Source tables retain `n_statistical_units`,
`n_donors_tested`, specimen IDs and biological patient IDs as distinct fields.

### Selection rules

- **HD:** primary-window specimen-pooled BH FDR < 0.05 and at least two eligible
  specimens. All three sections receive equal weights when eligible.
- **Xenium:** primary-window specimen-pooled BH FDR < 0.05 and at least two
  specimens. Because there is one field per independent patient, this matches
  the existing equal-donor Xenium calculation.
- **Cross-platform:** BH FDR < 0.05 separately in the three-specimen HD pool
  and the six-specimen Xenium pool, with matching nonzero platform mean-effect
  directions. Each FDR is adjusted over its complete mode/window family. Dot
  size uses the larger platform FDR. Displayed combined effects give every
  eligible specimen equal weight. This is an intersection of separately
  analysed datasets, not a new joint replicability-FDR guarantee.

Overview fill is blue–white–red around zero, clipped symmetrically at the larger of
0.25 or the 95th percentile absolute primary-window significant effect across both native spatial analyses. The same limit is shared by all three spatial comparison reports. Size uses capped
−log10 FDR (maximum 12); the cross-platform legend names the maximum platform
FDR explicitly. Means remain descriptive and have no estimated confidence
intervals. Displaying only significant hypotheses does not remove cohort,
coverage or annotation uncertainty.

Detail fill uses the original scATLAS palette (`#08519C`, `#F8F4ED`, `#B7212E`)
and the same shared spatial effect limits. Size retains the original physical
range while replacing co-positive percentage with capped −log10 FDR, as
confirmed by the user. No spatial effect or FDR value is described as a
Spearman coefficient or co-positive percentage.

Output locations:

- HD: `ref_outs/visium_hd_outs/spatial_interactions/figures/Auto_spatial_interaction_four_mode_significant_comparison.pdf`.
- Xenium: existing `figures/Auto_final/01_xenium_contact_enrichment.pdf`.
- Cross-platform: existing `figures/Auto_final/01_cross_platform_contact_enrichment.pdf`.
- Each root retains `tables/Auto_significant_comparison_figure_source.csv` and
  `reports/Auto_significant_comparison_summary.csv` plus its compact update.
- Spatial detail replot inputs are `tables/Auto_comparison_detail_{nodes,edges}.csv`
  and `logs/Auto_detail_plot_spec.json`. Current assembly components are
  `intermediate/Auto_comparison_{overview,details}.pdf`. They remain internal
  reproducibility inputs rather than additional selected figures.

The formerly historical HD PDF at the requested path is regenerated from the
current specimen-pooled comparison tables under the explicitly documented
independence assumption. Original section tests are reused without alteration. The original R HD producer calls the shared renderer after its
scientific tables have been written, using the fixed bidcell_temp Python binary.

### scATLAS summary

`ref_outs/non_malignant_mp_correlations/Auto_four_mode_interaction_summary.pdf`
is a five-page report directly under the requested parent folder. Its overview
uses the union of original positive cancer-focal selections from each mode's
`Auto_celltype_interaction_dotmap_data.csv` and cancer-focal edges from
`Auto_celltype_correlations_negative.csv`. Consequently it preserves the
upstream nominal Spearman/display rules; it does not silently substitute the
separate FDR-based scATLAS/HD directional-screen rule. Colour is Spearman rho,
size is co-positive sample percentage. Pages 2–5 are the original cancer-focal
first pages extracted without alteration from modes 01, 03, 05 and 06
`Auto_celltype_interaction_dotmap_with_neg.pdf`. Their colour scales, point
sizes, fonts, legends, LR annotations and every drawing are preserved exactly;
they are not reconstructed by the spatial renderer. Poppler `pdfseparate` and
`pdfunite` in the existing bidcell_temp environment extract and assemble these
pages. Live `intermediate/` keeps the current assembly components. Its exact source and counts are under
the parent `tables/` and `reports/` tiers. The maintained scATLAS R producer
calls the same renderer after exporting all mode tables. The canonical alias
reader in `scatlas_hd_xenium_directional.py` is a maintained dependency.

### Execution and validation

The existing matching launcher supports `SCREF_COMPARISON_FIGURES_ONLY=TRUE`
for the shared four-report redraw. Submit with one CPU, 8 GB and 30 minutes:

```bash
/opt/pbs/bin/qsub -l select=1:ncpus=1:mem=8gb -l walltime=00:30:00 -v SCREF_COMPARISON_FIGURES_ONLY=TRUE job_scripts/spatial/spatial_platform_neighbourhoods.sh
```

Other analytical/replot modes retain their previous launcher behavior.
`--stage hd-figures` and `--stage atlas-figures` are the R producers' explicit
plot-only callbacks. `--stage comparison-figures` regenerates all four reports.
No scores or spatial permutations are recalculated by these commands. A
separate audit reconstructs the cross-platform FDR intersection, checks native
HD/Xenium significance, verifies five pages in each report, checks exact detail
page dimensions, renders every page and checks text bounds. It also compares
each appended scATLAS page against its original PDF page to verify identical
rendering, text and drawing geometry. Figure reports
and indexes identify the current selection rule; older full numerical results
remain reproducible independently of these display changes.

The exact non-malignant spatial signatures are the `gene_sets` member of
`ref_outs/visium_hd_outs/spatial_interactions/intermediate/Auto_visiumhd_tme_mp_gene_sets.rds`.
They retain GeneNMF silhouette >= 0 and the first 100 unique ranked genes per
programme. Xenium intersects these same lists with its panel. Original complete
lists are `$metaprograms.genes` in `ref_outs/nmf_<compartment>/MP_outs_default.rds`,
for fibroblast, endothelial, macrophage, nk, plasma, cd4 and cd8. Reading a
single spatial signature in R is reproducible directly from the canonical RDS:

```r
tme <- readRDS("ref_outs/visium_hd_outs/spatial_interactions/intermediate/Auto_visiumhd_tme_mp_gene_sets.rds")
tme$gene_sets[["fibroblast__MP9"]]
```

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
Set `SCREF_FIGURE_SCOPE=final-comparisons` with `SCREF_REPLOT_ONLY=TRUE` to
update both final bubble PDFs, their source/summary tables and indexes without
regenerating supporting figures. For a scoped update of the one-page cross-platform final PDF and its current
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
