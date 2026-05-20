# Progress Update — 07 May 2026

**Project:** scRef_Pipeline / Parse_Pipeline / PDOs_Pipeline

---

## Context
This update centered heavily on the Parse pipeline, specifically inferring cell trajectories, pseudotime, and CNV comparisons between T0 and PDO. We also filtered geneNMF metaprograms based on their post-treatment trends and validated scATLAS states spatially using the Visium dataset.

## Summary of Work
* Computed cell trajectory and pseudotime using Monocle3, inferring sample distance from T0.
* Filtered geneNMF programmes (nMP117) by expected trend after treatment (mean/median trends).
* Annotated selected geneNMF programmes using 3CA pan-cancer and Hallmarks enrichment.
* Conducted Parse samples CNV inference and compared profiles (PDO vs T0 and all 9 samples).
* Extended CNV subclone analysis to the PDOs pipeline, analyzing cohorts and QC metrics.
* Spatially validated scATLAS malignant cell states using Visium tumor bins and distribution summaries.

## Key Outputs
| File/Output | Description |
|---|---|
| `Auto_parse_sample_pseudotime_combined.pdf` | Parse Cell Trajectory / Pseudotime from Monocle3. |
| `Auto_parse_highres_mean_median_trends_nMP117_selected.pdf` | Metaprogram trends after treatment. |
| `Auto_parse_cnv_heatmap_PDO_vs_T0_Carroll_2023_v3.pdf` | Parse CNV Profile comparing T0 vs PDO. |
| `Auto_PDO_cnv_subclone_mp_cohort_summary.pdf` | PDOs CNV subclonal analysis summary. |
| `Auto_visium_spatial_states.pdf` | scATLAS states spatial validation (Visium) final assignment. |

## Decisions Made
1. Utilized post-treatment mean/median trends to select and filter high-resolution metaprograms.
2. Assessed Visium spatial state scoring to map scATLAS states accurately.

## Next Steps
- [ ] RNA velocity parse (implementing from loom files)
- [ ] States concordance analysis for different MP gene list selection
- [ ] scATLAS identify subclones CNA patterns and states/MP correlation
- [ ] PDOs RNA velocity analysis
- [ ] PDOs FLOT matched pairs geneNMF analysis + trend filtering
