# Progress Update — 30 March 2026

**Project:** OAC Single-Cell Reference Pipeline (scRef_Pipeline)

---

## Context

This cycle focused on three main areas: (1) external validation of scATLAS metaprograms in adult stomach and Barrett's oesophagus epithelial references, (2) cross-platform bulk survival analysis integrating TCGA and GEO GSE19417 data with QC-based sample filtering, and (3) PDO differential expression analysis comparing treated vs untreated organoids across cell states.

## Summary of Work

### scATLAS External Reference Validation
- Scored scATLAS nMP=19 metaprograms on external epithelial datasets (adult stomach, Barrett's oesophagus)
- Generated annotation vs reference heatmaps showing MP activity patterns in non-malignant developmental contexts
- Validated that MPs map meaningfully to known epithelial differentiation trajectories

### Bulk Survival Analysis
- Developed cross-platform QC workflow harmonising TCGA RNA-seq and GEO GSE19417 microarray
- Applied standardisation (log2(TPM+1) for TCGA, z-score normalisation per dataset)
- QC filtering removed low-expression outliers and histology-inconsistent samples
- Ran survival volcanoes on:
  - TCGA unfiltered vs QC-filtered (comparison)
  - GEO GSE19417 EAC-only cohort
  - Pooled TCGA + GEO with dataset covariate

### Clinical Associations
- Generated sample-level state composition (stacked barplots sorted by diversity)
- Clinical association boxplots for states and MPs across treatment/histology variables

### PDO DGE Analysis
- Completed Type1 and Type2 differential expression analyses on matched PDO pairs
- Type1: DEGs per state across all PDOs
- Type2: DEGs by state comparing treated vs untreated within matched pairs
- Generated Venn diagrams, UpSet plots, and combined DEG heatmaps

## Key Outputs

| Output | Path | Description |
|--------|------|-------------|
| External MP validation (stomach) | `ref_outs/Auto_external_epi_mp_ucell/...annotation_vs_reference.pdf` (p2) | MP activity in adult stomach epithelium |
| External MP validation (Barrett's) | `ref_outs/Auto_external_epi_mp_ucell/...annotation_vs_reference.pdf` (p4) | MP activity in Barrett's oesophagus |
| Sample abundance | `ref_outs/task3_sample_abundance/Auto_task3_sample_abundance.pdf` (p7) | State proportions per sample |
| Clinical boxplots (states) | `ref_outs/Auto_clinical_assoc_state_boxplots_final.pdf` (p3) | State associations with clinical variables |
| Clinical boxplots (MPs) | `ref_outs/Auto_clinical_assoc_mp_boxplots_final.pdf` (p3) | MP associations with clinical variables |
| Bulk QC review | `ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_review.pdf` (p1) | TCGA + GEO harmonisation QC |
| TCGA survival (unfiltered) | `ref_outs/task2_survival/...volcano...pdf` (p1) | Pre-QC TCGA survival volcano |
| TCGA survival (filtered) | `ref_outs/task2_filtered_survival/...volcano...pdf` (p1) | Post-QC TCGA survival volcano |
| GEO survival | `ref_outs/geo_survival/geo_task2_survival/...volcano...pdf` (p1) | GSE19417 EAC-only survival |
| Cross-platform survival | `ref_outs/bulk_crossplatform/survival/...volcano...pdf` (p5) | Pooled TCGA + GEO survival |
| PDO hybrid states | `PDOs_Pipeline/PDOs_outs/Auto_PDO_hybrid_proportion_noreg.pdf` | PDO state proportions |
| PDO hybrid network | `PDOs_Pipeline/PDOs_outs/hybrid_pairwise/Auto_PDO_hybrid_pairwise_nodeplot_noreg.pdf` | Pairwise hybrid transitions |
| PDO Type2 volcano | `PDOs_Pipeline/PDOs_outs/DGE_matched_analysis/Type2_volcano_plots_by_state.pdf` (p2) | Treated vs untreated DEGs by state |
| PDO DEG overlap | `PDOs_Pipeline/PDOs_outs/DGE_matched_analysis/Type1_Venn_up.png`, `Type1_UpSet_up.pdf` | State DEG overlap |
| PDO DEG heatmap | `PDOs_Pipeline/PDOs_outs/DGE_matched_analysis/Combined_DEG_heatmaps.pdf` (p3) | Combined DEG expression |

## Decisions Made

1. **Cross-platform harmonisation**: Used shared genes only, with dataset-specific transforms (log2 TPM+1 for TCGA, keep log-scale for GEO), followed by row-wise z-score standardisation.

2. **QC strategy**: Applied binary Keep/Remove decisions based on expression-strength metrics and histology consistency in PCA space.

3. **Pooled survival model**: Included dataset as covariate in Cox regression to account for platform batch effects while testing MP/state hazard ratios.

4. **PDO DGE analysis types**:
   - Type1: Per-state DEGs aggregated across all PDOs
   - Type2: Within-state DEGs comparing treated vs untreated matched pairs

## Next Steps

- [ ] Validate cross-platform survival findings in additional external cohorts if available
- [ ] Investigate direction heterogeneity between TCGA and GEO for key MPs
- [ ] Complete SCENIC regulon analysis on final states
- [ ] Generate final publication figures with consistent styling
