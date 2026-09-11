# OCCAMS Bulk Survival and Clinical Association Methodology

## Scope and inputs

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
