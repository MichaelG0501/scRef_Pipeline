# Progress Update — 09 April 2026

**Project:** OAC Single-Cell Reference Pipeline (scRef_Pipeline)

---

## Context

This cycle focused on three areas: (1) cross-platform meta-analysis survival validation integrating TCGA and GEO EAC cohorts using random-effects models, (2) SCENIC regulon analysis for both scATLAS and PDO datasets to identify transcription factor networks governing metaprograms and cell states, and (3) derivation of robust state-specific marker genes using a per-sample recurrence-based workflow.

## Summary of Work

### Cross-Platform Meta-Analysis Survival
- Completed random-effects meta-analysis (DerSimonian-Laird) across TCGA RNA-seq and GEO GSE19417 microarray for all MPs and cell states
- Tested continuous, median-split, and quartile-split Cox models for each feature
- Generated EAC-only volcano plots comparing reference-based vs DGE-derived gene set scoring approaches
- Computed heterogeneity statistics (I², τ², Cochran's Q) per feature

### SCENIC Regulon Analysis
- Ran pySCENIC on scATLAS epithelial cells assigned to final MPs and cell states
- Ran pySCENIC on PDO cells assigned to PDO-specific MPs and cell states
- Generated regulon specificity scores (RSS), top-100 AUC heatmaps, and network visualizations for both datasets
- Implemented intermediate file caching to avoid redundant SCENIC recomputation
- Removed legends and subtitles from network plots for cleaner presentation

### State Marker Gene Derivation
- Developed per-sample recurrence-based marker gene workflow for the 6 finalized epithelial states
- Three-step pipeline: pooled candidate screening → per-sample Wilcoxon DGE → cross-sample/study recurrence filtering
- Generated z-scored heatmap of top 15 markers per state
- Classic Proliferative yielded 959 markers; Basal to Intestinal Metaplasia 142; other states fewer due to stricter specificity requirements

## Key Outputs

| Output | Path | Description |
|--------|------|-------------|
| Meta-analysis volcano | `ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_volcano_eac_only.pdf` | Random-effects meta-analysis survival volcano (MP + State, 3 split methods) |
| Meta-analysis results | `ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_random_effects_results.csv` | Full random-effects output (HR, CI, I², τ²) |
| Direction summary | `ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_direction_summary.csv` | TCGA vs GEO directional comparison |
| scATLAS SCENIC network (MP) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_network.pdf` | TF regulon network by metaprogram |
| scATLAS SCENIC network (State) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_state_network.pdf` | TF regulon network by cell state |
| scATLAS regulon heatmap (MP) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_regulon_heatmap.pdf` | RSS heatmap by metaprogram |
| scATLAS regulon heatmap (State) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_state_regulon_heatmap.pdf` | RSS heatmap by cell state |
| scATLAS top-100 AUC (MP) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_network_top100auc.pdf` | Top regulons AUC heatmap by MP |
| scATLAS top-100 AUC (State) | `ref_outs/final_mp_scenic/Auto_final_mp_scenic_state_network_top100auc.pdf` | Top regulons AUC heatmap by state |
| PDO SCENIC network (MP) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_network.pdf` | PDO TF network by MP |
| PDO SCENIC network (State) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_state_network.pdf` | PDO TF network by state |
| PDO regulon heatmap (MP) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_regulon_heatmap.pdf` | PDO RSS heatmap by MP |
| PDO regulon heatmap (State) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_state_regulon_heatmap.pdf` | PDO RSS heatmap by state |
| PDO top-100 AUC (MP) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_network_top100auc.pdf` | PDO top regulons AUC by MP |
| PDO top-100 AUC (State) | `PDOs_Pipeline/PDOs_outs/final_mp_scenic/Auto_PDO_final_mp_scenic_state_network_top100auc.pdf` | PDO top regulons AUC by state |
| Marker methodology | `analysis/methodology/Auto_six_state_marker_methodology.md` | Per-sample recurrence workflow documentation |
| Marker heatmap | `ref_outs/Auto_six_state_markers/Auto_six_state_marker_heatmap.pdf` | Top 15 markers per state (z-scored) |
| Marker summary | `ref_outs/Auto_six_state_markers/Auto_six_state_marker_summary.csv` | All aggregated marker statistics |

## Decisions Made

1. **Meta-analysis approach**: Used DerSimonian-Laird random-effects model to account for expected heterogeneity between RNA-seq (TCGA) and microarray (GEO) platforms. Reported I² and Cochran's Q alongside meta HR.

2. **SCENIC intermediate caching**: Implemented file-based caching for adjacency matrices, AUC matrices, and regulon objects to enable iterative visualization refinements without rerunning the full ~8h SCENIC pipeline.

3. **Marker gene strategy**: Chose per-sample recurrence over pooled DGE to avoid large-sample dominance. Genes must be significantly upregulated in ≥20% of eligible samples and ≥35% of studies to qualify as markers.

4. **Marker specificity rule**: In addition to recurrence, each marker must have the highest median expression in its target state compared to all other states (specificity gap > 0).

## Next Steps

- [ ] TME ligand-receptor interaction analysis: database selection, filtering, and visualization
- [ ] snRNA-seq malignancy classification and scATLAS state mapping
- [ ] Present snRNA-seq results systematically once all components are complete
- [ ] Finalize publication figures with consistent styling across all analyses
