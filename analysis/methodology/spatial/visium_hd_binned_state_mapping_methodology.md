# Visium HD Binned Malignant scATLAS State Mapping

## Scope

This workflow maps the current centred-refined scATLAS metaprograms and
Approach-B states onto 16 um Visium HD bins for `SUR1231`, `FFPEA1`, and
`FFPED1`. `SUR1122` and all segmented-cell outputs are excluded.

State mapping is restricted before scoring to bins satisfying both:

1. `is_epithelial_target == TRUE` in the finalized binned malignancy table;
2. `Auto_malignancy` is `malignant_level_1` or `malignant_level_2`.

Thus normal epithelial, CNA-non-malignant, unresolved epithelial, reference,
non-epithelial, and post-filter-rejected bins cannot affect expression
standardization, MP scores, state thresholds, maps, or colocalisation.

## Run Order

1. `export_visium_hd_binned_scatlas_signatures.R`
2. `visium_hd_binned_state_mapping.py`
3. `visium_hd_binned_state_colocalisation.R`

The PBS wrapper `run_visium_hd_binned_states.sh` runs these in order. Set
`SCREF_REPLOT_ONLY=TRUE` to reuse the live mapped per-bin table and rebuild maps,
abundance summaries, and colocalisation without reading Space Ranger matrices.

## Signature Export

The exporter reads the persistent centred-refined gene lists from:

`ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`

Negative-silhouette parent MPs were removed upstream during centred MP triage;
the refined/merged object is the current downstream gene-set object. The export
retains the top 100 ranked genes per MP for spatial scoring. Signature CSVs and
metadata are saved in live storage under
`ref_outs/visium_hd_outs/state_mapping/intermediate/signatures/`.

The state-defining groups exactly match the centred-refined noreg state script:

- **Classic proliferation:** `MP2+`
- **Basal to intestinal metaplasia:** `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, `MP10+`
- **SMG to intestinal metaplasia:** `MP8+`, `MP8b`, `MP16`, `MP18b`, `MP17`, `MP2x`
- **Stress adaptive:** `MP12`
- **Cancer-cell immune mimicry:** `MP15`

These labels are retained verbatim from
`analysis/metaprograms/centred/06_centred_refined_state_definition_noreg.R`;
legacy display aliases such as `SMG-like Metaplasia` and `Immune Infiltrating`
are not applied. The signature metadata records MD5 hashes for both that state
definition script and the merged refined MP gene object, making the exact
inputs used by each spatial export auditable.

Cell-cycle MPs `MP1`, `MP5`, and `MP13+` are scored for audit but never define
states. `MP11c` and `MP18a` remain excluded from state definition.

## MP Scoring

Only malignant level-1/level-2 epithelial count rows are loaded. Genes common
to all three Space Ranger matrices are retained so every sample uses the same
signature features.

For each bin, counts are normalized to CP10K and transformed with `log1p`. For
each sample separately, every gene is standardized across its malignant bins:

`z_gene = (log1p_CP10K - sample_gene_mean) / sample_gene_sd`

An MP raw score is the unweighted mean of its available top-100 gene z-scores.
This is the sparse-matrix implementation of the legacy Visium HD equation. MP
scores are then sample-centered and divided by the pooled standard deviation
across malignant bins from the three spatial samples. Signature coverage is
exported for audit.

The important change from the legacy run is the scale population: the legacy
script scored all epithelial observations and applied malignancy afterward.
The finalized workflow fits every spatial expression and MP scale using only
the explicitly requested malignant level-1/level-2 epithelial bins.

## State Assignment

For each state, its score is the maximum adjusted score among its defining MPs.
The highest group is selected, then the canonical centred-refined Approach-B
rules are applied:

- `Unresolved` if the maximum state-group score is below `0.5`;
- `Hybrid` if the top-minus-second state-group gap is below `0.3`;
- otherwise assign the highest state group.

The top-MP diagnostic is restricted to state-defining non-cell-cycle MPs and
uses the same `0.5` minimum. All raw MP scores, adjusted MP scores, group scores,
state calls, gaps, malignancy evidence, and spatial coordinates are persisted
in the live per-bin mapping table.

## Colocalisation

The primary raw colocalisation metric reproduces the legacy calculation.
Within each sample, `Hybrid` and `Unresolved` bins are first excluded. For each
remaining biological-state bin, `FNN::get.knn()` identifies its six nearest
biological-state malignant bins in full-resolution pixel coordinates. The raw
score is:

`same_neighbor_score = number of six neighbours with the focal state / 6`

There is deliberately no distance cutoff in this reproduction. Consequently,
the six-neighbour metric can bridge a spatial gap when malignant bins are
sparse, and common states have a higher chance baseline. To keep those
limitations auditable without changing the reproduced result, the workflow
also saves:

- mean and maximum neighbour distance per bin;
- the focal state's within-sample abundance;
- `same_neighbor_excess`, defined as raw score minus that abundance;
- `same_neighbor_ratio`, defined as raw score divided by that abundance.

The raw same-neighbour score remains the primary colocalisation figure. The
abundance-adjusted excess answers a complementary question: whether local
same-state frequency is greater than the frequency expected from that state's
sample-wide abundance. It is included as a separate report page and standalone
figure, while the unadjusted six-neighbour fraction remains the primary
colocalisation result. Both standalone plots pool per-bin values across samples
and use identical boxplot/jitter styling and dimensions; the abundance
subtraction itself is still calculated separately within each sample.

Spatial state and top-MP maps draw all other post-filter retained bins first in
light grey to show tissue context. Those background bins do not enter MP
scoring, state assignment, abundance, or colocalisation. Malignant state points
are drawn on top with a minimum scatter size of `1.8` rather than `1.2`.

## Outputs

All outputs are persistent under `ref_outs/visium_hd_outs/state_mapping/`:

- `intermediate/signatures/`: exported MP/state definitions;
- `tables/`: per-bin scores/states, signature coverage, parameters, abundance,
  per-bin colocalisation, and sample/pooled summaries;
- `figures/`: per-sample state/top-MP maps, abundance plots, colocalisation
  plots, and a multi-page spatial summary;
- `logs/`: mapping and colocalisation run summaries.

Compact abundance and colocalisation summaries are also written to
`updates/new_updates/summaries/`.
