# scATLAS significance with concordant Visium HD and Xenium trends

## Purpose and execution

The maintained producer is `analysis/spatial/scatlas_hd_xenium_directional.py`.
It updates `ref_outs/spatial_scatlas_hd_xenium_directional/` in place. The final
hit rule, requested on 16 November 2026, requires scATLAS significance and mean
effect direction agreement in both spatial platforms. Neither spatial assay
must be statistically significant. The precise interpretation is
**scATLAS-significant, spatially concordant interactions**.

The user retained the previously agreed scATLAS FDR <0.05 definition. This
secondary exploratory rule was revised after examining earlier results; the
final list has no formal spatial replication FDR guarantee. The four comparison
modes are separate families, not an experiment-wide FDR family. Sign agreement
can occur by chance or with effects arbitrarily close to zero.

The matching ignored launcher `job_scripts/spatial/scatlas_hd_xenium_directional.sh`
uses one CPU, 8 GB and 30 minutes, with mandatory dmtcp initialization followed
by the existing bidcell_temp Python environment. This table-level analysis
reuses persisted UCell scores and permutation inference.

## Inputs and hypotheses

Four full correlation tables are read from
`ref_outs/non_malignant_mp_correlations/<mode>/Auto_celltype_correlations_all.csv`:
`01_cancer_mps_cross_only`, `03_cancer_states_cross_only`,
`05_cancer_mps_vs_whole_celltypes`, and `06_cancer_states_vs_whole_celltypes`.
Every edge with exactly one malignant endpoint enters its respective mode.
Direction comes from Spearman rho; BH is calculated over all finite Spearman
P values within that complete mode before any spatial filtering. No old network
display threshold or Pearson gate is used. The denominator is retained per row.

Current malignant MP IDs, order and full labels come from the centred refined
17-MP grouping table. Explicit historical description aliases map to current
IDs; unknown labels stop execution. TME programme IDs are compartment-specific;
malignant descriptions must not be assigned to identically numbered TME MPs.
Whole CD4/CD8 atlas hypotheses share one spatial T-cell test. Both hypotheses
remain in the scATLAS family; summary counts also report unique spatial tests.
The spatial evidence cannot establish whole CD4/CD8 subtype specificity.

Visium HD effects come from
`ref_outs/visium_hd_outs/spatial_interactions/tables/Auto_spatial_interactions_pooled.csv.gz`.
Only the primary ring-one ±16 µm axis-distance window is used. The three sections
retain two independent donor weights: FFPEA1 and FFPED1 are paired observations
of one donor; SUR1231 is the other. The descriptive pooled effect is the mean
of eligible donor effects. Donor effects average the eligible sections within
each donor. No favourable radius or individual section is selected.

Xenium means come from `ref_outs/spatial_cross_platform/tables/Auto_within_platform_meta.csv`.
The donor table `Auto_harmonised_interactions_by_donor.csv.gz` independently
reconstructs the Xenium means. Input freshness is checked against
`ref_outs/spatial_cross_platform/logs/Auto_input_provenance.json`. Changed inputs
or an old inference version stop the run. The original HD and Xenium
methodologies define UCell positivity, eligibility and local permutation nulls.
This producer does not repeat those analyses.

## Exact final-hit rule

1. scATLAS within-mode Spearman BH FDR is strictly below 0.05.
2. HD has a finite mean log2 enrichment from at least two eligible donors and
   its nonzero sign matches the scATLAS Spearman sign.
3. Xenium has a finite mean log2 enrichment from at least two eligible donors
   and its nonzero sign matches that same scATLAS sign.

These three conditions define `passes_directional_rule`. There is no spatial
P-value, FDR, Z-sign, minimum effect-size, unanimous-donor or donor-majority gate.
A descriptive mean can differ in sign from signed-Stouffer evidence; the table
retains this diagnostic without using it to override the requested mean trend.
Missing or zero effects fail the trend requirement. Missing parent labels and
fewer than two eligible donors are explicitly unavailable.

The table retains HD and Xenium P/FDR values and per-donor heterogeneity.
`hd_confirmed` remains a stricter diagnostic requiring HD FDR <0.05 and matching
HD Z direction. Its export `hd_confirmed_pairs.csv` is supplementary; it does
not define the final hit list. All final hits are in the existing
`xenium_directionally_consistent_pairs.csv` and `figure_source.csv` paths.

## Interpretation and coverage

scATLAS sample-level co-occurrence and local spatial enrichment are different
measurements. Agreement supports a biological hypothesis but does not reproduce
a correlation coefficient or demonstrate physical contact or signalling. HD
bins and Xenium cells differ in measurement scale. Treatment and cohort
heterogeneity remain. Only two HD donors are available.

The Xenium matrix header supplies the actual gene count: 4,999 in this release.
Signature coverage comes from the current `Auto_signature_coverage.csv`, with
separate malignant source and TME target coverage fractions. State assignments
have no artificial single-signature coverage value. The six Xenium donor fields
contain the paper's 16 tissue regions. NK/plasma lack matched Xenium parent
annotations. Low coverage may affect sensitivity but cannot be assumed to
explain a nonsignificant result.

## Stable outputs and replotting

- `tables/all_pairs.csv`: complete hypotheses, effects, P/FDR, availability,
  direction checks, final decisions, subtype caveats and signature coverage.
- `tables/xenium_directionally_consistent_pairs.csv`: final-hit list.
- `tables/hd_confirmed_pairs.csv`: stricter HD-significance diagnostic subset.
- `tables/donor_evidence.csv`: spatial donor data for all scATLAS-significant pairs.
- `tables/selection_summary.csv`: attrition and unique spatial-test counts.
- `tables/signature_coverage.csv` and `tables/figure_source.csv`: replot sources.
- `figures/Auto_final/01_scatlas_hd_xenium_evidence.pdf`: two landscape pages per
  populated mode, at 300 mm width with 12 pt ordinary text (7.32 pt if reduced
  to 183 mm). Enrichment hits precede depletion hits, with canonical MP/target
  order within each direction. Red denotes enrichment and blue depletion.
  Page one separates scATLAS rho from spatial donor-mean log2 enrichment.
  Labels are black, with one line per interaction. Fib./Mac./Endo. abbreviate
  target parents, epi. abbreviates epithelium, and IFN-I abbreviates type I
  interferon. Full names remain alongside `pair_axis_label` in figure_source.csv.
  Page two repeats the scATLAS Spearman-rho panel alongside separate HD and
  Xenium sample panels, all with the same row order. Small grey
  points show eligible individual sample effects on each pair's row. Larger
  red/blue diamonds show the exact platform means used on page one. HD means
  retain two donor weights, averaging FFPEA1/FFPED1 within the paired donor;
  they are not arithmetic means of all three plotted section points. No
  uncertainty intervals are estimated. Missing sample effects are omitted;
  their availability and values remain in `tables/sample_figure_source.csv`.
- `reports/CURRENT_RESULTS.md`, `reports/validation.csv` and
  `logs/input_manifest.json`: current interpretation, checks and provenance.
- `updates/new_updates/summaries/scatlas_hd_xenium_directional.csv`: compact summary.

Run from the project root:

```bash
/opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE job_scripts/spatial/scatlas_hd_xenium_directional.sh
```

`SCREF_FORCE_REBUILD=TRUE` updates this screen from persistent inputs.
`SCREF_REPLOT_ONLY=TRUE` requires matching input sizes, modification times and
method version. Exactly one flag must be TRUE; there is no automatic fallback.
No new analysis script or output root is required for this strategy update.

Validation covers upstream provenance, donor identities, unique hypotheses,
HD BH values as a diagnostic, independently reconstructed Xenium means, source
freshness during execution, independent R BH/decision checks and final PDF
content/bounds at export width and its 183 mm equivalent. PBS exit zero and verified outputs are required.
