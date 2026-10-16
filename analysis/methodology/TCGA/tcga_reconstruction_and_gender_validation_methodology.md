####################
# TCGA ESCA reconstruction and scRef gender-validation methodology
####################

## Scope

This methodology covers the active scripts under `analysis/TCGA/`:

- `tcga_esca_reconstruct_data.R`
- `tcga_gender_state_mp_validation.R`

The goal is to rebuild the TCGA-ESCA bulk RNA-seq expression matrix and clinical metadata without relying on the deleted historical TCGA directory, then test whether scRef sex-associated metaprogram and state signals show concordant directions in TCGA EAC bulk RNA-seq.

## TCGA reconstruction

`tcga_esca_reconstruct_data.R` queries the public GDC API for open-access `TCGA-ESCA` RNA-seq files with:

- data category: `Transcriptome Profiling`
- data type: `Gene Expression Quantification`
- workflow type: `STAR - Counts`
- experimental strategy: `RNA-Seq`

The script writes a fresh GDC manifest and file metadata table under `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/source_data/`, downloads each STAR-count TSV into `intermediate/source_data/gdc_files/<file_id>/`, and verifies size and MD5 checksum before using the file. Existing verified files are reused. Set `SCREF_TCGA_SKIP_DOWNLOAD=TRUE` to process only existing downloads. Set `SCREF_TCGA_OVERWRITE_BAD=TRUE` only when a failed partial download should be replaced.

Clinical metadata come only from the supplied cBioPortal TSV for study `esca_tcga_gdc`. The initial import source is `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/clinical_cbioportal.tsv`; the canonical project snapshot is `ref_outs/TCGA/esca_gdc_reconstruction/intermediate/source_data/clinical_cbioportal.tsv`. The workflow does not call the cBioPortal API or generate patient/sample long/wide API exports. GDC diagnosis-level clinical T/N/M, stage, grade and location are queried separately and cached as `intermediate/source_data/Auto_gdc_diagnoses_query.rds`, with a readable `tables/Auto_gdc_diagnoses.csv` export. The processed metadata table joins GDC file/sample metadata with cBioPortal clinical fields and derives:

- `sample_barcode`, `case_barcode`, and `sample_type_code`
- broad GDC histology field `type`
- `HistologyGroup` as `EAC`, `ESCC`, or `Other`
- `Gender`, `OS_time`, `OS_event`, stage, grade, age, mutation burden, FGA, TMB, smoking, and selected sample annotations where available

TPM values are read from the STAR-count `tpm_unstranded` column. Rows with GDC synthetic count categories (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`) are removed. Gene IDs are version-stripped, gene symbols are taken from the STAR-count `gene_name` column, and duplicated symbols are summed to produce a gene-symbol-by-sample TPM matrix.

Primary outputs are written under `ref_outs/TCGA/esca_gdc_reconstruction/`. Compatibility copies are also written for older downstream scripts:

- `ref_outs/tcga_esca_meta.rds`
- `ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`

## scRef to TCGA gender validation

`tcga_gender_state_mp_validation.R` compares sample-level scRef OAC associations with TCGA EAC bulk associations.

For scRef:

1. Load `meta_full_epi.rds`, the centred-refined noreg state vector, final-17 centred UCell scores, final-17 MP genes, and clinical metadata from `Concise_Summary_EAC_Ref.xlsx`.
2. Recode clinical `Gender` to `Female`/`Male`.
3. Compute sample-level MP activity as the mean UCell score per `orig.ident`.
4. Compute sample-level state abundance as the percentage of epithelial cells in each final state, excluding `Unresolved` and `Hybrid`.

For TCGA:

1. Load reconstructed TPM and metadata.
2. Restrict the main validation cohort to primary-tumour EAC samples (`sample_type_code == "01"`, `HistologyGroup == "EAC"`).
3. Transform expression as `log2(TPM + 1)`.
4. Score scRef MPs and final-state gene sets using GSVA with Gaussian kernel.
5. Cache GSVA scores in `ref_outs/TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores_centred17.rds`.

MP gene sets are read directly from the final centred-refined 17-MP object. State gene sets are unions of the five current noreg Approach B groups in `analysis/shared/scRef_config.R`. No legacy nMP19 or appended 3CA gene sets enter the active validation.

## Statistics and concordance

Within each platform and feature, Female-vs-Male association is summarized using:

- female and male sample counts
- median feature value in each group
- median difference: Female minus Male
- Wilcoxon rank-sum p value
- BH-adjusted p value within feature type
- Cliff's delta as a direction-preserving effect size

Concordance is called for each feature by comparing the sign of scRef Cliff's delta with the sign of TCGA EAC Cliff's delta. The terminal figure contains:

1. a scRef-vs-TCGA effect-size scatter plot with concordant and discordant directions,
2. TCGA EAC MP boxplots for all 17 current MPs,
3. TCGA EAC state boxplots,
4. a direction bar plot showing TCGA effect sizes with scRef effect-size overlay.

The compact machine-readable summary is written to:

`updates/new_updates/summaries/Auto_tcga_gender_scRef_concordance_summary.csv`


## Storage, execution and rebuild contract

All persistent products are canonical inside `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/`. Both workflows create `intermediate/`, `tables/`, `figures/`, `logs/` and `reports/` under their registered `ref_outs/TCGA/` roots. Reconstruction stores original STAR-count downloads, the GDC manifest, query responses, download-verification status and the clinical TSV snapshot under `intermediate/source_data/`. Metadata CSV and the gene-key CSV belong in `tables/`; their RDS counterparts belong in `intermediate/`. Existing historical `raw/` and intermediate CSV files are retained, but the updated producer does not read or refresh them. The two documented compatibility copies remain live because existing downstream scripts consume them.

The reconstruction stores `intermediate/Auto_tcga_esca_expression_source.rds`, containing the gene key and gene-symbol TPM matrix. It reuses the pinned GDC file query and diagnosis query unless `SCREF_FORCE_REBUILD=TRUE`; verified downloaded files remain reusable even during a force rebuild. The expression cache fingerprints the manifest metadata, clinical snapshot and every contributing STAR-count file. Unversioned or stale computational caches stop with an explicit rebuild instruction. No input-file, schema or analysis-method substitution is permitted. STAR input requires `gene_id`, `gene_name`, `gene_type` and `tpm_unstranded`; malformed, negative or non-finite TPM values and duplicate gene IDs stop the run. A missing canonical clinical snapshot is imported from its documented external source exactly once. To change the snapshot later, explicitly replace the canonical input; force rebuild does not silently change source snapshots.

`SCREF_TCGA_SKIP_DOWNLOAD=TRUE` forbids downloading missing or invalid STAR-count files. `SCREF_TCGA_OVERWRITE_BAD=TRUE` explicitly permits replacing a failed-checksum download. `SCREF_TCGA_DOWNLOAD_TIMEOUT` sets the download timeout in seconds (default 1800). `SCREF_TCGA_MAX_FILES` is a positive debug limit; results produced with that limit are incomplete and must not be interpreted as full-cohort reconstruction. Reconstruction has no figures and rejects `SCREF_REPLOT_ONLY=TRUE`.

Gender validation consumes only canonical reconstruction RDS files. It imports the clinical workbook into `gender_validation/intermediate/source_data/Concise_Summary_EAC_Ref.xlsx` once and uses that persistent snapshot. It validates the exact 17 MP names and the five-state grouping CSV against shared configuration. MP2x, MP11c and MP18a are forbidden. Every MP/state set must retain at least five measured genes. Expression must remain finite after `log2(TPM + 1)`. GSVA-cache fingerprints cover the TPM matrix, MP genes, grouping CSV and shared configuration; gene-set identities, score dimensions and sample coverage are also validated. `SCREF_FORCE_REBUILD=TRUE` rebuilds scores. It does not change the five-state definition or the 0.5 unresolved / 0.3 hybrid thresholds in the upstream state producer.

The complete gender plot source is `intermediate/Auto_tcga_gender_plot_source_data.rds`, containing per-sample TCGA scores, long-form TCGA values, scRef sample summaries, statistics, concordance and gene sets. `intermediate/Auto_tcga_gender_figure_bundle.rds` stores rendered figure grobs and the PNG plot. `SCREF_REPLOT_ONLY=TRUE` loads only that bundle and the compact result summary; it never loads cell-level objects, recomputes tests or runs GSVA. A successful updated full run is required before the first replot. Changing plotting style requires a full run to rebuild the bundle. Replot regenerates the recorded PDF and PNG and writes a separate `logs/Auto_tcga_gender_replot_run_summary.rds/.txt`. Force rebuild and replot-only cannot both be enabled. The four-page PDF uses 20-by-11-inch pages; the scatter PNG remains 13-by-8.5 inches at 300 dpi. Summary CSV size is checked to remain below 100 KB.

Run through PBS, in this order:

```bash
/opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_esca_reconstruct_data.sh
# Wait for completion and verify metadata, TPM, gene key, staging workbook and run summary.
/opt/pbs/bin/qsub analysis/TCGA/Auto_tcga_gender_state_mp_validation.sh
```

Use `/opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE <wrapper>` for an explicit rebuild, or `-v SCREF_REPLOT_ONLY=TRUE` for gender figure regeneration. Both entry points require a PBS allocation, use `dmtcp`, initialize conda with the mandated hook and write live PBS logs with `#PBS -koed`. Reconstruction requests 2 CPUs, 24 GB and 6 hours; gender validation requests 4 CPUs, 12 GB and 4 hours. These are single jobs; they do not launch per-sample jobs.

## Statistical unit, filtering and limitations

The scRef unit is `orig.ident`, never a cell or a barcode-derived sample. Each sample contributes one mean UCell value per MP. State percentages use only cells assigned to one of the five biological states as the denominator; technical states are excluded. TCGA uses reconstructed RNA sample barcodes, with primary EAC as the main cohort and primary ESCC/combined ESCA as secondary comparisons. Gender tests require at least two samples of each sex; otherwise the test and effect estimate are unavailable. BH adjustment is separate for MP and state features within each cohort. Clinical `orig.ident` must be unique before joining.

The comparison is observational and unadjusted. scRef state proportions and TCGA state-union GSVA scores measure different quantities; concordance concerns effect direction, not equality of magnitudes. The clinical source can cover fewer subjects than the full ESCA query. The AJCC8 workbook is an annotation template restricted to supplied cBioPortal sample IDs, not a calculated AJCC8 stage assignment. GDC diagnosis selection retains the first diagnosis matching esophageal/gastroesophageal origin; absent fields remain missing. Metadata coalescing is the documented field-level source hierarchy, not an alternative analytical workflow.

## Validation

Validation covers R/shell syntax, loading all required packages in `dmtcp`, the installed Gaussian GSVA interface, storage-boundary rejection, cache build/reuse/staleness, exact 17-MP/group enforcement and figure regeneration from saved grobs. Validation fixtures run through PBS and persist separately under `ref_outs/TCGA/contract_validation/<job_id>/`. Syntax and fixture checks do not establish biological results. A production rerun must additionally verify sample alignment, finite matrices, expected feature counts, all registered exports and a readable multipage PDF before adopting new outputs.


The October 2026 contract validation completed as PBS jobs `4373148.pbs-7` and `4373174.pbs-7`, both state `F`, `Exit_status = 0`. The final run also exercised the actual STAR parser and TPM reconstruction on synthetic inputs, including rejection of alternate TPM schemas and negative TPM values. Expected fixture PDF, expression cache, replot logs and result report were verified under `ref_outs/TCGA/contract_validation/4373174.pbs-7/`. The installed GSVA version was 1.46.0. The test harness is retained as `temp/Auto_tcga_contract_validation.sh`; it does not execute production analyses or replace production outputs.


Production figure exports measure ggplot grobs on the final Cairo PDF device before caching. The production validation harness `temp/Auto_tcga_production_validation.sh` checks the four-page PDF, required outputs, exact MP columns, saved-grob replot and text bounds using the pinned project-local PyMuPDF 1.26.7 package. Validation inventories, page metadata and previews are retained under `gender_validation/reports/`.


Production regeneration on 15 October 2026 completed with gender job `4373277.pbs-7`, CNA job `4373346.pbs-7` and validation job `4373347.pbs-7`, all state `F` and exit status zero. The existing reconstructed TPM and metadata were reused. Validation verified 53 registered nonempty files, nine PDFs across 15 pages, text within page bounds and standalone cached replot. Separate CSV checks verified unique observational units and 87 TCGA samples per CNA event-feature comparison. Exact outputs and run provenance are recorded in `ref_outs/TCGA/cna_recurrent_event_validation/reports/Auto_production_run_summary.txt`; observation checks are in `reports/Auto_observational_unit_validation.csv` under the same CNA root. The canonical CNA arm threshold is 0.12 for this run.
