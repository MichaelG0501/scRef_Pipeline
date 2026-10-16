# Selected-cohort survival analysis

## Scope and source contract

`analysis/clinical/survival_selected_cohorts.R` creates the two PDFs specified on 15 October 2026. It preserves the original integrated-survival script and all existing work. Patient, not cell, is the statistical unit. Identifiers are full OCCAMS `SHA_ID` values and TCGA workbook `Patient_ID`/`RNA_primary_sample_ID` values. There is no barcode-trimming sample inference or alternative-input fallback.

Exact inputs are listed in the executable registry. The September 25 OCCAMS patient-master workbook and September 18 TCGA RNA88 workbook are read directly for current treatment grouping, specimen type, and harmonized clinical/pathological stage. OCCAMS RNA eligibility uses `Annotation_Compatible_Below_50_Flag == no` in the canonical GRCh37 QC table. TCGA OS uses the reconstructed primary-EAC metadata; its RFS source values come from the prior checksum-verified original GDC BCR XML audit, not the sparsely populated cBioPortal DFS fields.

Every analytical input is hashed into a run manifest; source identities and cohort counts are checked before analysis. The earlier 64-patient model-data and risk-table CSVs are historical verification records, not rerun dependencies. The standalone R script does not call audit scripts, render PDF previews, install software, or require Python/PyMuPDF. The processed TCGA tumour-event CSV remains essential clinical endpoint source data.

## Explicit OCCAMS exclusions and selected cohorts

The three user-specified exclusions are always applied by exact subject identifier:

- `26d51bb9b7c9974b5e753fc3520ffcfa2abdbd132509a887f4a5bc5630aa05d9`: organ-unknown endoscopy.
- `72089d94e1a29cc0b7df4fb8f9694343cda2c6e71ae1083fc266b5b36b62323d`: stomach palliative.
- `9ea99e16eaaf486784cb0cd7e1c846a4dcc7bbaecfccc9a04d7853093f2cf775`: stomach surgery.

This leaves 278 of the original 281 QC-pass RNA subjects in the audit data. Initial-curative neoadjuvant and surgery-pathway RNA cohorts contain 165 and 72. The user explicitly selected only the 157 neoadjuvant subjects with typed, matched RNA specimen metadata: 77 biopsy and 80 surgery specimens. The eight untyped neoadjuvant RNA subjects remain in the patient audit, but enter none of the requested plots or MP scores. The surgery pathway applies the same RNA_Metadata_Matched == Yes rule, retaining 68 of 72 RNA subjects. All 157 matched neoadjuvant records have known specimen type. Of the 68 matched surgery-pathway records, 47 are surgery specimens, 18 biopsy, and three explicitly unknown; the three unknown-source matched records remain eligible because this pathway has no requested specimen-specific split. The combined OCCAMS stage cohort therefore has 225 subjects (157 + 68). The previous mixed-selection count of 229 was inconsistent and is superseded. Specimen source denotes the RNA specimen; it is not substituted with eventual operation status or treatment timing.

TCGA uses exactly the 86-patient workbook `Surgery group`, comprising surgery `Yes` 72, `No` 12, and `Unresolved` 2. The two explicit neoadjuvant-plus-surgery cases are outside this group. The 86-person pathway must not be relabelled as 86 confirmed untreated-before-surgery cases: 83 preoperative-treatment fields are unrecorded and three explicitly report no treatment. The cohort itself is not reduced to a stage-specific or endpoint-specific complete-case set. Those exclusions occur separately within each analysis panel.

## OS and RFS construction

OCCAMS OS time is `Survival_Days`; the event is 1 for `Dead`, 0 for the source label `Alive/Not recorded dead`, matching the original script's convention. An absent/nonpositive time is ineligible. The audit preserves the source's imperfect censoring label and does not claim confirmed living status for all non-death records. Its original time origin remains undocumented, and remains documented here; plotted axes use Time (years).

The user explicitly authorized an exploratory OCCAMS RFS assumption: recurrence and death/follow-up intervals are treated as sharing a first-recorded-therapy origin. Existing `Estimated_Recurrence_Time_Days_From_Diagnosis_Surrogate` values were verified in the source audit to equal recurrence date minus first recorded therapy date. Despite the column name, these are not measured diagnosis-origin recurrence times. Exploratory RFS is calculated as follows:

1. Missing/nonpositive survival time is excluded.
2. Explicit recurrence `Yes` without a positive timed recurrence is excluded, including deceased patients, because the earlier recurrence time is unknown.
3. A recurrence interval beyond the reported survival/follow-up interval is excluded as inconsistent rather than truncated or repaired.
4. A living/non-death record with unknown recurrence status is excluded. An observed death can contribute the composite death event when recurrence status is unrecorded; an earlier undocumented recurrence remains a limitation.
5. Explicit timed recurrence contributes the recurrence event; for a deceased subject with both times, the earlier recurrence/death interval is used. Explicit `No` or an unrecorded recurrence with observed death uses the reported survival interval and death indicator.

At the user’s request, panels are labelled **RFS** and report actual evaluable sample/event counts without explanatory footers. The common-origin assumption remains explicit in this methodology, settings, patient audit, and run report. This is a recurrence/death association endpoint under an approved assumption, not a verified native clinical PFS or proof that every neoadjuvant patient achieved a disease-free baseline. Diagnosis/therapy offsets, incomplete progression ascertainment, and specimen timing prevent a stronger interpretation.

TCGA RFS preserves the exact earlier source definition and patients: earliest recorded new tumour event or death, with latest recorded positive follow-up for censoring. A living patient without recorded tumour-event ascertainment is excluded. Original source-audit invalid-time exclusions remain applied. One event is `New Primary Tumor` (TCGA-L5-A43E), retained to reproduce the previously supplied cohort/definition rather than silently changing its endpoint. The run report and patient audit state this explicitly. Time zero is initial pathologic diagnosis, as specified by the official ESCA follow-up form, not automatically surgery or first treatment. Sources: https://gdc.cancer.gov/system/files/public/file/Esophageal%20Followup%20Form.pdf .

The known-pathological-stage surgery-yes RFS analysis must reproduce exactly 64 patients, 42 events, stage counts 10/14/28/12, log-rank p = 0.001081798772, and the supplied numbers at risk at years 0, 2, 4, 6, 8, and 10. The five unknown-stage and three living-without-ascertainment exclusions are preserved. Other TCGA panels may have different evaluable N while retaining parent cohort sizes 86 and 72.

## Stage KM panels and page order

Clinical and pathological stage remain separate. Only explicit harmonized stages 1–4 enter a stage comparison; `unknown` remains nonempty in the patient audit and is excluded per panel, never imputed or counted as a biological stage. Kaplan–Meier curves use `survival::survfit`; censored observations are marked, and stage-specific starting sample counts, global evaluable N, events, unknown-stage counts, and risk tables are shown. Log-rank p-values use `survival::survdiff` and are reported only with at least 20 evaluable patients, five events, and two stage groups. Sparse KM panels still show their observed curves and risk counts, with a descriptive-only statement instead of an unstable p-value. This particularly affects TCGA surgery-yes clinical stage, known for only two of 72 subjects.

PDF 1 contains exactly seven 24-by-17-inch pages:

1. Combined OCCAMS typed-neoadjuvant plus surgery pathways (225): clinical-stage OS and exploratory RFS side by side.
2. OCCAMS typed neoadjuvant (157): clinical/pathological OS across the top, clinical/pathological exploratory RFS below.
3. OCCAMS neoadjuvant biopsy specimens (77): same four-panel arrangement.
4. OCCAMS neoadjuvant surgery specimens (80): same four-panel arrangement.
5. OCCAMS matched surgery pathway (68): same four-panel arrangement.
6. TCGA surgery pathway (86): clinical/pathological OS above, clinical/pathological RFS below.
7. TCGA surgery-yes subgroup (72): same four-panel arrangement. Surgery-no and unresolved patients are audited separately but are not added to this requested surgery-yes page.

All follow-up is shown. Time axes extend to at least ten years, with two-year ticks, and expand for longer observed follow-up. Curve source tables and every number-at-risk value are persistent live outputs.

## MP expression/scoring and OS Cox models

Only the current refined 17-MP collection is used, ordered by the canonical grouping CSV. MP2x, MP11c, and MP18a are prohibited. Cell-cycle MPs MP1, MP5, and MP13+ remain testable because the request covers all 17 MPs. Each gene set requires at least five matched genes, and every gene set must yield finite scores; missing genes/objects do not invoke a substitute source.

OCCAMS uses the existing QC-pass TMM-normalized gene-symbol logCPM matrix, whose expression filtering/normalization was defined on the upstream 281-subject cohort. GSVA is recomputed using only the selected 157 neoadjuvant subjects. TCGA uses the persistent gene-symbol TPM matrix restricted to the 86 selected primary-tumour cases. Gaussian-kernel GSVA on raw TPM deliberately matches `survival_integrated_analysis.R`; it is not silently changed to the log2(TPM+1) transformation used by another TCGA workflow. The installed GSVA 1.46.0 API is used directly. The surgery-pathway 68 are independently scored as their own Gaussian-GSVA parent population using the same normalized OCCAMS logCPM input. The existing 157-person neoadjuvant and 86-person TCGA scoring populations remain unchanged. Parent-cohort scores are reused for the 77/80 OCCAMS specimen subgroups and TCGA surgery-yes subgroup, avoiding a change in the score definition when drawing a subgroup plot.

PDF 2 contains exactly six pages: OCCAMS typed neoadjuvant combined (157), its biopsy subgroup (77), its surgery-specimen subgroup (80), OCCAMS matched surgery pathway (68 of 72 QC-pass pathway subjects), TCGA surgery pathway (86), and TCGA surgery-yes (72). The added surgery-pathway page has 67 OS-evaluable subjects and 43 deaths; one matched-RNA subject lacks a positive survival interval. It is distinct from neoadjuvant surgery specimens. Each page contains four **OS-only** univariable Cox volcanoes in a 2-by-2 layout:

- Continuous GSVA score standardized to one SD within the endpoint-evaluable displayed cohort.
- Median split, High strictly greater than the median versus Low at or below it.
- Q1/Q4, retaining scores at/below the 25th percentile or at/above the 75th percentile, comparing High with Low.
- Within-cohort optimal cut: score quantiles 20%–80% in 5% steps, High at/above the cut versus Low below it. Each MP independently selects the smallest fitted raw Wald p-value; ties use the first quantile. Discovery and projection across datasets are not performed because the user selected within-cohort cuts.

Models require at least 20 evaluable patients, at least five deaths, nonzero score variance, and, for binary models, at least five patients and at least one death per group. A non-estimable candidate is recorded with its explicit eligibility reason; no different model or endpoint is substituted. Unexpected fitting errors, convergence warnings, and nonfinite coefficients stop the run. All 17 MP/method rows are retained in the output, including eligibility status. The default model is unadjusted `Surv(OS_time, OS_event) ~ split_val`, without stage, purity, age, or treatment covariates.

Outputs include hazard ratios, 95% CIs, log-HRs, SEs, raw Wald p-values, BH FDR within cohort/method, sample/event/group counts, cut values, selected quantiles, and all 13 optimal-cut candidates per feature. Volcanoes show log2(HR) against minus-log10(raw Wald p), with the 0.05 line and raw-p significance colouring. The minimum optimal-cut p-value is exploratory and is not corrected for searching thresholds; BH across MPs does not correct that search. Methods, model tables, and run reports state this limitation; the user requested no explanatory figure footers. No significance result is labelled independently prognostic or causal.

## Storage, caching, replot, and validation

All persistent files remain under the live project. The current output root is `ref_outs/clinical_selected_survival/`, with `figures/`, `tables/`, `intermediate/`, `reports/`, `logs/`, and `source_data/`. Historical `validation/` previews are archived on rerun; new previews are not generated. Only the two current PDFs occupy `figures/`. Original source-audit CSVs are copied into `source_data/`; historical workbooks-audit/XML-audit trees and superseded runs are preserved in `archive/prior_structure_20261015_145606/`. No files are deleted. A source-data JSON manifest records moves and copies. `validation_dependencies` is a compatibility symlink to `source_data/software`, preserving unrelated TCGA validation jobs’ existing dependency path.

Before subsequent runs, the R script itself archives nonempty current generated output tiers to `archive/before_<run_id>/`, never overwriting an earlier revision. Current job-specific logs stay in `logs/`; previous logs and compact summaries are archived with the previous revision. Analysis caches are keyed by input hashes, package versions, and settings; a matching cache can be copied from an archived revision. `SCREF_FORCE_REBUILD=TRUE` recomputes scores/models. `SCREF_REPLOT_ONLY=TRUE` requires a complete matching cache. Both flags cannot be true. A compact summary is written under `updates/new_updates/summaries/`.

Plot titles omit PDF/page numbers; panel subtitles contain only the requested endpoint/count fields. Stage panels show `n=... | Events = ... | unknown removed = ...`; log-rank results appear within the curve area. Every stage time axis uses years with two-year intervals, exceeding the user’s minimum six-month tick spacing. Full MP names occupy fixed left/right label columns, with connector lines to unchanged HR/p coordinates. All 17 names appear in every MP panel. Q1/Q4 patient membership varies with each MP’s score ranking, so death counts legitimately vary across MPs. The displayed event range is the minimum–maximum of the exact per-MP model counts; `q1q4_patient_membership.csv` makes each count auditable.

The original one-time MP16/MP18b direction and PDF-layout checks are complete. They are not part of reruns. The retired validator and original complex wrapper remain as `analysis/clinical/legacy_validate_selected_survival_pdfs.py` and `analysis/clinical/legacy_survival_selected_cohorts.sh` for provenance; the old audit R script was already absent when this simplification began. No current executable calls them. Historical audit results remain in live archives.

Run the R file directly on an allocated dmtcp compute node: `Rscript analysis/clinical/survival_selected_cohorts.R`. The optional 17-line `analysis/clinical/survival_selected_cohorts.sh` only requests four CPUs, 24 GB and two hours, activates dmtcp, and runs that command. It contains no workflow logic and is not an R dependency. PBS remains necessary for running this analysis from a login node. There is no commit operation in this workflow; the R script can be staged on its own.

Routine runs retain normal input/cohort, model-estimability, expected output-count and file-existence guards inside R. They produce seven stage pages, six volcano pages, 26 KM panel records and 408 Cox result rows. They do not repeat historical reference-cohort comparisons, independent model refits, PDF text/typography checks, or preview rendering. Statistical definitions, cohort membership, model estimates and figure styling are unchanged by this operational simplification.
