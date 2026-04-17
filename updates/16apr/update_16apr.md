# Progress Update — 16 April 2026

**Projects:** `snSeq_Pipeline` and `scRef_Pipeline`

---

## Context

This cycle focused on two main areas: (1) rebuilding the snRNA-seq workflow from raw counts through expression filtering, inferCNA, malignancy labelling, and scATLAS state mapping, and (2) finalising cross-celltype metaprogram interaction analysis with ligand-receptor annotation in the scATLAS reference atlas.

## Summary of Work

### snSeq_Pipeline
- Rebuilt stages 1 to 6 directly from raw 10x inputs across 77 samples and 3 technologies.
- Step 1 QC uses fixed thresholds after CPM plus log normalisation: mitochondrial RNA `< 15%`, detected genes `>= 200`, and mean housekeeping expression `>= 0.5`.
- Step 2 performs per-sample PCA, UMAP, and graph clustering, then assigns initial cluster labels using top-ranked markers from the 3CA-derived marker workbook.
- Step 3 derives robust OAC-specific annotation markers and refines cluster labels using weighted marker scores plus canonical-marker rescue.
- Step 4 removes weak-lineage and mixed-lineage cells using matched-lineage marker support and loose co-expression filtering, leaving 167,982 final singlets from 428,892 post-QC cells; 99,144 retained cells are epithelial.
- inferCNA uses batch-matched non-epithelial references (`cynthia_sn`, `multiomes`, `gemx`) with no global fallback. Across 41 completed samples, 99,144 epithelial cells were evaluated and 55,664 malignant epithelial cells were retained.
- snRNA-seq state mapping uses the 14 scRef-retained `nMP = 19` metaprograms, the `noreg` workflow only, and Approach B top-MP group assignment on 55,664 malignant epithelial cells.
- Before unresolved relabeling, the dominant state labels were Hybrid (27.9%), Unresolved (21.7%), Basal to Intestinal Metaplasia (16.0%), Stress-adaptive (12.9%), and Classic Proliferative (8.1%).
- Unresolved relabeling using retained 3CA pan-cancer metaprograms reassigned 1,168 cells.

### scRef_Pipeline
- Finalised cross-celltype metaprogram co-occurrence analysis on the full atlas using sample-level adjusted scores: the percentage of MP-positive cells per sample within each compartment.
- Cancer MPs were scored only in malignant epithelial cells; non-malignant MPs were scored only within their matched TME compartments.
- Applied the standard silhouette filter, required non-trivial sample coverage, and tested only celltype pairs with adequate matched-sample overlap.
- Built positive and negative cross-celltype networks from Pearson and Spearman correlation results.
- Added ligand-receptor support only to retained positive edges using `Ligand_Receptor_Pairs.xlsx`, top 4,000 ranked genes per node, and evidence restricted to `literature supported` or `putative`.

## Key Outputs

| Output | Path | Description |
|--------|------|-------------|
| Annotation review PDF | `snSeq_Pipeline/sn_outs/snseq_annotation_featureplots_single_F_post_N2_biopsy.pdf` | Well-resolved annotation example |
| Annotation review PDF | `snSeq_Pipeline/sn_outs/snseq_annotation_featureplots_single_I_post_N1_biopsy.pdf` | Weak / messy annotation example |
| Pre/post QC heatmaps | `snSeq_Pipeline/sn_outs/Auto_QC_snSeq_prefilter.png`, `snSeq_Pipeline/sn_outs/Auto_QC_snSeq_final.png` | Cohort heatmaps before and after expression filtering |
| snRNA-seq summary plots | `snSeq_Pipeline/analysis/summary/*.png` | Retention, composition, and UMAP overview panels |
| Malignancy overview | `snSeq_Pipeline/sn_outs/Auto_malignancy_scatter_comparison.png` | inferCNA scatter comparison across batches |
| Malignancy validation | `snSeq_Pipeline/sn_outs/analysis/infercna/validation_cancer_signatures_vs_cna.pdf` | Cancer signature support for CNA-based labelling |
| State mapping outputs | `snSeq_Pipeline/sn_outs/Auto_topmp_v2_noreg_*.pdf` | Heatmap, proportion, and CC-score plots |
| Bubble overview | `ref_outs/non_malignant_mp_correlations/Auto_cross_celltype_correlation_bubble.pdf` | Cross-celltype MP co-occurrence by focal compartment |
| Correlation networks | `ref_outs/non_malignant_mp_correlations/Auto_cross_celltype_positive_network.pdf`, `...negative_network.pdf` | Positive and negative cross-celltype MP networks |
| LR-annotated network | `ref_outs/non_malignant_mp_correlations/Auto_cross_celltype_positive_network_lr_annotated.pdf` | Positive network annotated with candidate LR pairs |
| LR summary panel | `ref_outs/non_malignant_mp_correlations/LR_Summary.png` | Summary of candidate ligand-receptor interactions |

## Other Ongoing Tasks

- Finalise unresolved snRNA-seq cell-state relabeling with retained 3CA pan-cancer metaprograms, then inspect the two samples with matched PDO partners.
- Compare the alternative presentation frameworks and keep the clearest version.
- Explore membrane and surfaceome databases for downstream receptor-focused filtering.
- Download and analyse the OSCC PDO dataset.
- Perform CNA subclone analysis.
