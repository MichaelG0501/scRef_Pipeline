# OCCAMS Bulk Survival and Clinical Association Methodology

## Scope and inputs

### Corrected-metadata clinical-association revision (17 October 2026)

The clinical-association script now reads `analysis/clinical/OCCAMS_clinical_updated_RNA_source_status_FINAL_2026-09-25.xlsx`, sheet `Patient_Master_388`, directly. This explicit corrected patient master supersedes the older resolved CSV for every covered clinical field. File modification times are not used to select authority. Future curated releases require an explicit path/column contract update rather than silently choosing a different file. The older `Auto_OCCAMS_qc_pass_subject_metadata.csv` remains an input only for treatment-plan comorbidity, which is absent from the workbook, and for a change audit. Missing or conflicting corrected values are not filled with older values. The upstream survival producer and its older metadata/endpoint logic are unchanged by this request; historical statements below about that producer do not describe the revised association metadata authority.

The association audit retains all 281 QC-pass subjects, identified by exact SHA IDs, while plotting precisely the five OCCAMS selections used in `selected_cohort_stage_survival.pdf`: combined matched neoadjuvant plus matched surgery pathways (225), typed matched neoadjuvant (157), its biopsy specimens (77), its surgery specimens (80), and matched surgery pathway (68). The same three exact organ/site exclusions documented in `survival_selected_cohorts_methodology.md` are applied to these selections. Neoadjuvant selection requires initial curative intent, final modality `neoadjuvant (chemo + CRT)`, RNA metadata matched `Yes`, and biopsy/surgery specimen source. The surgery pathway requires initial curative intent, final modality `surgery`, and matched RNA metadata, without dropping its three unknown-source specimens. No OS/RFS availability filter is applied to clinical-score plots. These are overlapping cohort views, not independent validation datasets.

Source mappings are stored in the inventory and metadata audit. Demographic and timepoint-specific RD/PS/RP fields use their corresponding corrected workbook columns; no cross-timepoint anatomical substitution is made. Treatment intent and modality use the corrected initial-pathway `Treatment_Intent_Final` and `Curative_Treatment_Modality_Final`, not the all-record treatment history. Surgery means `Surgery_Ever`; response uses `Curative_Neoadjuvant_Clinical_Response`. Pretreatment and pathological nodal variables use curated `cN` and `pN`. Their historical internal names ending `_tnm7` are retained for compatibility, but titles explicitly describe harmonized N stage and no longer claim these are untouched TNM7 source fields. The workbook's `Clinical_Stage_AJCC` and `Pathological_Stage_AJCC` are additional separate comparisons; only explicit 1–4 are retained. Corrected final site/Siewert, RNA specimen type/treatment status, and neoadjuvant type are also inventoried and tested when eligible. Source-specific missing site fields remain missing even when a final curated site is available.

Blank, unknown, not-recorded, and Tx/Nx/Mx entries are missing. Comma/semicolon/pipe-separated repeated entries require one unique nonmissing token; discordant histories are not reduced to their first value or treated as new biological categories. Age is split at >60 years. Performance status uses 0 versus 1–2; grade uses well/moderate versus moderate/poor; positive nodes use 0, 1–3, >=4; Mandard uses TRG1–2 versus TRG3–5; clinical response uses CR/PR versus SD/PD. Surgery-procedure groups use the single concordant corrected procedure string. T-stage substages are reduced to T0–T4. Missing values and sparse categories remain visible in cohort-specific inventory counts and exclusion reasons.

The association script reuses the upstream 281-subject Gaussian-GSVA MP, MP-union state, and top-20-DGE state matrices without re-estimating scores separately in each selection. All subsets therefore share one scoring scale. This differs deliberately from selected-survival's parent-cohort-specific scoring and does not claim numerical equivalence to those scores. Expected dimensions are 281 by 17, 281 by five, and 281 by five, with finite values and exact QC identities. State signatures must match the current five state names. Cell-cycle MPs remain on the MP page but not in the upstream state unions.

Within each cohort and clinical variable, a level requires at least 10 unique subjects and at least two levels must remain. Eligibility is recalculated independently per cohort; constant pathway/specimen variables are not assigned artificial contrasts. Each signature is tested by two-sided unpaired Wilcoxon rank-sum (`exact=FALSE`) for two groups or Kruskal–Wallis for more groups. BH adjustment is within cohort, clinical variable, and signature family (MP, State MP union, State DGE). Stars represent the adjusted values themselves: <0.05, <0.01, <0.001. Test errors or nonfinite results stop the run without another analysis method. Associations remain exploratory and unadjusted for clinical confounding; bulk scores mix malignant and nonmalignant expression. Pathology/response/recurrence are post-baseline variables, not validated treatment effects.

Figures are now in `ref_outs/OCCAMS/clinical/figures/selected_cohorts/`: five cohort PDFs for MPs and five for states, named `Auto_OCCAMS_<cohort>_clinical_scores_<MP|State>.pdf`. Each PDF has one page per eligible clinical variable. MP pages display all 17 MPs; State pages contain separate rows for the five MP-union and five DGE signatures. Box/dodge/jitter/theme conventions are retained, with fixed jitter seed, sample-count legends, cohort titles, and GSVA score rather than percentage axes. No empty PDF is accepted. `Auto_OCCAMS_clinical_association_pdf_inventory.csv` records exact page counts and cohort sizes.

The script writes its own corrected 281-subject association metadata, per-field change/source audit, cohort membership, input hash manifest, variable inventory, test results, and complete plot/cache RDS under the existing live `tables/` and `intermediate/` tiers. It does not overwrite upstream metadata or score caches. All previous outputs owned by this association script, including the old two general-cohort PDFs, are moved recoverably into `archive/clinical_associations/before_<run_id>/`. Unrelated OCCAMS clinical outputs are untouched. A compact cohort/signature-family summary remains under `updates/new_updates/summaries/` and is checked below 100 KB.

Cache contracts include all input/script hashes, package versions, source mappings, cohort counts, and the minimum level size. Changes automatically invalidate stale caches. `SCREF_FORCE_REBUILD=TRUE` rebuilds association data/statistics; `SCREF_REPLOT_ONLY=TRUE` requires a matching complete cache and redraws the ten PDFs without recalculating tests. Both flags together fail. PBS launcher: `job_scripts/clinical/Auto_occams_bulk_clinical_associations.sh`, requesting two CPUs, 12 GB and one hour in dmtcp. Submit from the project root; scheduler stdout/stderr go directly to that root with `-koed`. Validation requires exit status 0, PASS report, five PDFs per type, nonempty files, PDF page counts matching inventory, exact cohort membership/counts, and inspection of representative MP and both-state pages.

The completed revision used analysis job `4390787.pbs-7` and independent one-time PDF/membership check `4390791.pbs-7`, both with exit status 0. All ten PDFs passed; 240 pages were checked for all 17 MP labels or both state-family rows. Cohort membership matched selected-survival exactly. Eligible variables/pages per PDF were 30, 28, 22, 21, and 19 for combined, neoadjuvant, biopsy, neoadjuvant surgery specimens, and surgery pathway respectively. There were 3,240 association tests. Representative combined-cohort MP and State pages were visually inspected. One-time check results and four previews are in `figures/selected_cohorts/validation/`; they are not executable dependencies of the association script. Previous association outputs are preserved in `archive/clinical_associations/before_20261017_142302_1934856/`.

This workflow analyses the corrected GRCh37 OCCAMS gene-count matrix with the current centred refined 17-MP panel. The statistical unit is one OCCAMS subject, identified by the unmodified SHA ID. It never reads BAMs and it never uses the historical GRCh38-annotated count file.

The inclusion rule is exactly `Annotation_Compatible_Below_50_Flag == no` in `OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv`. This retains 281 of 282 matrix subjects. The excluded subject remains in the canonical source matrix; filtering occurs only in these downstream scripts.

## Metadata resolution and endpoints

The source metadata contain repeated rows for some SHA IDs. Resolution is strictly within one named source column: documented encodings are normalized within each row, and a subject value is retained only when all non-missing normalized values from that same column agree. Discordant repeated-row values become missing and are written to `Auto_OCCAMS_metadata_conflicts.csv`. First-row selection is forbidden.

No clinical column is filled from another column. The prefixes are retained as distinct contexts: `DI` demographic, `RD` diagnostic, `PS` final pretreatment staging, `TP` treatment-plan, `ST` surgery, `TR` treatment response, and `RP` resection pathology. Diagnostic, pretreatment, and resection tumour site and Siewert variables are therefore separate. TNM6 and TNM7 N-stage fields are also separate; the association workflow tests the TNM7 fields only because TNM6 coverage is too sparse, while retaining TNM6 values in the subject metadata for provenance.

Overall survival is the one explicit endpoint construction requiring two source columns:

- event = 1 and time = `Deceased Survival Days` when a concordant death time is present;
- event = 0 and time = `Last Known Survival Days` otherwise;
- a subject missing both times is not included in a survival model.

`FE End Point` and `EP End Point` are retained for audit but do not override an explicit death-time field. In the 281-subject QC-pass cohort, this yields 280 evaluable subjects, 189 deaths, and 91 censored observations. Survival time is in days.

The recurrence source column is explicitly longitudinal within a single field. `recurrence_ever_recorded` is Yes when that source field contains any `Yes`, No when it contains at least one `No` and no `Yes`, and missing otherwise. It is labelled exploratory and is never combined with recurrence date or site columns.

## Expression normalization and scoring

GENCODE v19 gene IDs are mapped to gene symbols using the same pinned GRCh37 GTF used by featureCounts. Counts mapping to the same symbol are summed. A gene is retained at CPM >= 1 in at least 10% of the 281 subjects. Library sizes are TMM-normalized with edgeR and transformed to log2 CPM with prior count 1.

GSVA with Gaussian kernel scores all 17 current centred refined MPs. State signatures are unions of constituent current MP genes in the five biological state groups. The three cell-cycle MPs (MP1, MP5, and MP13+) remain individually testable but are excluded from state unions. MP2x, MP11c, and MP18a are absent from the current input and are not reintroduced. Gene-set coverage is reported.

## Survival analysis

For every MP, state-union, and state-marker score, univariable Cox proportional-hazards models are fitted using:

1. a continuous score standardized to one cohort standard deviation;
2. high versus low at the cohort median;
3. upper versus lower quartile, excluding the middle two quartiles.

The exploratory optimal-cut display searches score quantiles from 20% through 80% in 5% increments and remains univariable. It does not adjust for pathology N stage, pretreatment N stage, or CIBERSORTx malignant fraction. Its minimum raw p-value is descriptive and is not interpreted as a cutpoint-search-adjusted confirmatory p-value.

Each model result records cohort, endpoint, time/event construction, scaling, split thresholds, sample/event counts, hazard ratio, 95% confidence interval, raw p-value, and BH FDR. BH adjustment is performed within feature type and split method for the prespecified models. These are association models, not causal or independently prognostic models.

## Clinical variables and association tests

Every tested variable is mapped to one exact source column in `Auto_OCCAMS_clinical_variable_inventory.csv`, together with its source prefix/timepoint. RD, PS, TP, ST, TR, and RP variables are displayed separately and are never used as substitutes for one another. Main staging comparisons use PS TNM7 pretreatment N stage and RP TNM7 pathology N stage as separate analyses.

Derived groupings are limited to transparent recoding within one source column: performance status 0 versus 1-2; grade well/moderate versus moderate/poor; positive nodes 0, 1-3, or >=4; Mandard TRG1-2 versus TRG3-5; and treatment response CR/PR versus SD/PD. `Mx`, `Nx`, `Tx`, unknown, and not-recorded values are missing rather than biological levels. Adjacent Barrett-associated dysplasia is present only for explicit low- or high-grade entries and absent only for explicit no dysplasia. Operation groups are derived solely from `ST Procedure`.

Only levels with at least 10 subjects are displayed or tested, and a variable must retain at least two such levels. MP and state scores are tested at subject level with Wilcoxon rank-sum for two retained groups and Kruskal-Wallis for more than two. BH adjustment is within clinical variable and feature type. No cell is treated as an independent replicate. Treatment response, resection pathology, Mandard response, and recurrence are post-baseline exploratory associations and are labelled by their context.

## Cache semantics, validation, and limitations

All normalized matrices, scores, source-column-specific metadata, model data, statistical tables, plot data, figures, inventories, and summaries are persistent in live storage. `SCREF_FORCE_REBUILD=TRUE` recreates computational caches. `SCREF_REPLOT_ONLY=TRUE` redraws figures from persistent model/plot data without rerunning normalization or GSVA.

Both scripts hard-stop unless metadata and score objects contain exactly 281 QC-pass subjects with matching identities. Validation requires successful PBS exit, PASS run reports, non-empty PDFs/tables, 17 MP and five state scores, 280 survival-evaluable subjects with 189 events, and manual inspection of representative figure pages. Bulk GSVA reflects mixed tumour and microenvironment expression and must not be interpreted as tumour-cell-specific activity without orthogonal validation. Sparse or conflicting clinical fields are reported rather than imputed.
