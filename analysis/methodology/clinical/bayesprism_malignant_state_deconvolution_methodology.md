# Historical OCCAMS BayesPrism malignant-state deconvolution

This document describes the retained OCCAMS-only trial, not the current workflow. Its state-coverage rule incorrectly assigned 21,171 malignant cells without current state labels to the nonmalignant epithelial reference. Use the [corrected common-reference workflow](malignant_bulk_survival_methodology.md). The historical producers now have `legacy_` filenames and their results are excluded from current survival inference.

## Scope

This is the first BayesPrism trial and uses OCCAMS only. It estimates:

- total malignant fraction of each bulk sample;
- abundance of the five current malignant states in bulk and within the malignant compartment;
- malignant-compartment expression inferred by BayesPrism;
- current metaprogram activity within that inferred malignant compartment.

The five primary states are `Classic proliferation`, `Squamous-to-intestinal`,
`Glandular-to-intestinal`, `Stress-adaptive`, and `Cancer-cell immune mimicry`.
`Unresolved` and `Hybrid` remain explicit technical malignant states. They are
reported separately and included in total malignant fraction.

## Inputs and cohort restriction

The script uses the canonical GRCh37 OCCAMS matrix
`ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz` and retains
exactly the 281 subjects with
`Annotation_Compatible_Below_50_Flag == no` in
`ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv`.
The count matrix contains the 281 retained subjects plus one documented QC-failing subject; the script drops that extra subject.

The single-cell reference is `ref_outs/EAC_Ref_merged.rds`. Malignant cells are
the intersection of `celltype_update == epithelial`, cell IDs in
`ref_outs/meta_full_epi.rds`, and names of the current state vector. Epithelial
reference cells without current state labels were retained in the fitted trial
as `Non_malignant_epithelial`. This label is not independent evidence of
nonmalignancy and can bias malignant fractions downward. The fitted result is
therefore a computational trial requiring biological validation. Cells labelled
`unresolved_inconsistent` are excluded.
`orig.ident` is required as the reference sample identity field.
Current state labels come only from the centred refined no-reg state vector.

## Reference construction

Raw RNA counts are aggregated into one pseudobulk gene-expression profile per
current cell state. All nonmalignant cell types remain separate reference
profiles. All malignant cells share BayesPrism cell type `Malignant`; their
five primary, `Unresolved`, and `Hybrid` labels are BayesPrism cell states.
Each state must contain at least 20 cells.

The reference and OCCAMS matrix are mapped to GENCODE v19 symbols, recognizing
versioned and unversioned Ensembl IDs and existing symbols. Duplicate symbols
are summed. No marker-only restriction is applied; at least 1,000 shared symbols
are required. BayesPrism then applies its default mixture-outlier filtering.

BayesPrism receives raw-count GEP profiles with `key = "Malignant"`.
It estimates initial state fractions, initial type fractions, final type
fractions, and malignant-type expression. BayesPrism does not provide a final
state posterior in this GEP configuration; state abundances therefore use its
initial state posterior. Total malignant abundance is reported from both
initial state sums and final malignant-type posterior.

## Malignant MP activity

The inferred malignant-type expression is mapped through the pinned GENCODE
v19 GRCh37 gene annotation. Duplicate gene symbols are summed. Expression is
normalized within each sample to CPM and transformed as `log2(CPM + 1)`, then current refined MP gene sets are scored
with Gaussian-kernel GSVA. The script uses the current 17 MPs and stops if any
MP has fewer than five genes present.

GSVA receives genes in rows and subjects in columns. Its output has MPs in
rows and subjects in columns and is transposed for the sample-by-MP CSV.
Validation checks exact subject identity/order, 17 finite MP scores, nonnegative
expression, positive tumour library sizes, and posterior fractions summing to
one within `1e-6`. These checks establish computational validity, not biological
accuracy. Fractions represent RNA/read contributions, not direct cell counts.

## Outputs and validation

Outputs are written under `ref_outs/OCCAMS/bayesprism/`:

- `intermediate/Auto_OCCAMS_bayesprism_result.rds`: fitted BayesPrism object;
- `intermediate/Auto_OCCAMS_bayesprism_tumour_expression.rds`: inferred malignant expression;
- `tables/Auto_OCCAMS_bayesprism_state_fractions.csv`: bulk and tumour-only state fractions;
- `tables/Auto_OCCAMS_bayesprism_mp_scores.csv`: tumour-only MP scores;
- `tables/Auto_OCCAMS_bayesprism_gene_coverage.csv`: MP gene coverage;
- `tables/Auto_OCCAMS_bayesprism_reference_labels.csv`: reference label audit;
- `reports/Auto_OCCAMS_bayesprism_run_report.txt`: compact run audit;
- `updates/new_updates/summaries/Auto_OCCAMS_bayesprism_summary.csv`: compact summary.

The script stops on missing inputs, wrong subject count, subject mismatch,
invalid counts, missing state labels, unexpected states, underrepresented
states, insufficient gene overlap, duplicate BayesPrism labels, or inadequate
MP coverage. `SCREF_FORCE_REBUILD=TRUE` is required to rebuild an incomplete or
complete cache. No fallback deconvolution method is used.

`SCREF_RESUME_EXTRACTION=TRUE` requires the fitted result RDS and label audit,
skips reference loading and sampling, and regenerates extraction outputs and
the report. It is mutually exclusive with `SCREF_FORCE_REBUILD=TRUE`. The PBS
wrapper uses `set -euo pipefail` so an R failure produces a failed PBS exit.

## Visual QC

`analysis/clinical/legacy_Auto_occams_bayesprism_qc_plots.R` reads the persistent fraction,
MP score, coverage and QC CSVs and renders a four-panel PDF/PNG dashboard under
`figures/Auto_OCCAMS_bayesprism_qc_overview.{pdf,png}`. Panels show initial versus
updated malignant RNA fractions, MP gene coverage, within-tumour state
composition, and tumour-only MP activity. Composition and heatmap share sample
ordering by updated malignant fraction. MP scores are standardized within each
MP across subjects; colour saturation is capped at +/-2.5 SD, while saved CSVs
retain the full values. No clustering, significance test, or accuracy threshold
is applied. The dashed scatter line is identity, not a target imposed on the
reference update. Initial and updated estimates are not independent validation.

The plot script verifies 281 unique QC-pass subjects, the exact 17 MPs, finite
scores, nonnegative fractions bounded by one, and tumour-state closure and
initial bulk-state closure within `1e-6`. It writes sample order, full heatmap
values and state means/medians to `tables/Auto_OCCAMS_bayesprism_qc_*.csv`, a
report to `reports/Auto_OCCAMS_bayesprism_visual_qc.txt`, and a compact summary
to `updates/new_updates/summaries/Auto_OCCAMS_bayesprism_visual_qc.csv`.
`SCREF_FORCE_REBUILD=TRUE` or `SCREF_REPLOT_ONLY=TRUE` regenerates these outputs
from persistent CSVs. The dashboard explicitly records the reference-label
limitation and distinguishes numerical integrity from biological accuracy.
