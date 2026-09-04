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

The script writes a fresh GDC manifest and file metadata table under `ref_outs/TCGA/esca_gdc_reconstruction/raw/`, downloads each STAR-count TSV into `raw/gdc_files/<file_id>/`, and verifies size and MD5 checksum before using the file. Existing verified files are reused. Set `SCREF_TCGA_SKIP_DOWNLOAD=TRUE` to process only existing downloads. Set `SCREF_TCGA_OVERWRITE_BAD=TRUE` only when a failed partial download should be replaced.

Clinical metadata are pulled from the public cBioPortal API study `esca_tcga_gdc`. Patient and sample clinical attributes are saved in both long and wide forms. The processed metadata table joins GDC file/sample metadata with cBioPortal clinical fields and derives:

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
2. TCGA EAC MP boxplots for the top scRef sex-associated MPs,
3. TCGA EAC state boxplots,
4. a direction bar plot showing TCGA effect sizes with scRef effect-size overlay.

The compact machine-readable summary is written to:

`updates/new_updates/summaries/Auto_tcga_gender_scRef_concordance_summary.csv`
