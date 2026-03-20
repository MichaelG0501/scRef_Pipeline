# Progress Update — 19 March 2026

**Project:** scRef_Pipeline and PDOs_Pipeline

---

## Context
This cycle focused on refining state assignments, performing downstream survival and comparative analyses (e.g., CIBERSORTx vs TCGA bulk), and exploring specific biological relationships like Barrett's vs Intestinal Metaplasia. Pan-cancer MP correlations and PDO state concordance were also analyzed to validate findings across datasets.

## Summary of Work
- Updated topMP unresolved state assignments with per-cell heatmaps.
- Computed self MP correlations including the new unresolved MPs, visualizing Jaccard similarity and per-sample correlation.
- Compared CIBERSORTx deconvolution activities to raw whole TCGA expressions.
- Conducted four-in-one GSVA survival analyses for Metaprograms and Cell States, comparing malignant subsets, whole cohorts, and DGE methods.
- Analyzed the Spearman correlation of proportions between Barrett's Oesophagus and Intestinal Metaplasia.
- Compared MP expression between OAC and OSCC cohorts from TCGA.
- Performed detailed pseudotime distribution mappings per sample for both states and MPs (specifically Basal to Intestinal Metaplasia).
- Visualized Pan-cancer MP correlations (crossdata bar and scatter plots) and state concordance (barplot and heatmap) for PDOs.
- Generated MP and State abundance summary profiles per sample.
- Formatted hybrid cell mappings into a per-cell heatmap.

## Key Outputs
| File | Description |
|---|---|
| Auto_task4_unresolved_relabel_heatmap.pdf | Unresolved state relabeling heatmap. |
| Auto_task5_nMP19_jaccard_with_pancancer_heatmap.pdf | Jaccard self similarity including unresolved MPs. |
| Auto_task5_nMP19_correlation_heatmap_persample.pdf | Per-sample MP correlation heatmap. |
| mp_activity_comparison_heatmap_ordered.pdf | CIBERSORTx deconvolution vs whole TCGA comparison. |
| Auto_survival_tcga_mp_state_volcano_methods_reg_noreg.pdf | GSVA survival analysis for Metaprograms (four-in-one). |
| Auto_task2_survival_mp_state_volcano_methods_noreg.pdf | GSVA survival analysis for Cell States (four-in-one). |
| Auto_task7_basal_smg_relationship_plots.pdf | Barretts vs Intestinal metaplasia proportion spearman. |
| Auto_task8_tcga_eac_escc_compare_plots.pdf | OAC vs OSCC MP expression comparison. |
| Auto_task1_partA_top12_pseudotime_states.pdf | Pseudotime state-colored plots per sample. |
| Auto_task1_partB_Basal_to_Intestinal_Metaplasia...pdf | Pseudotime MP-colored plots for Basal metaplasia per sample. |
| Auto_PDO_mp_correlation_crossdata_bar.pdf | Pan-cancer MPs correlation barplot (PDOs). |
| Auto_PDO_mp_correlation_crossdata_scatter.pdf | Pan-cancer MPs correlation scatterplot (PDOs). |
| Auto_PDO_concordance_barplot.pdf | PDO state concordance barplot. |
| Auto_PDO_concordance_heatmap.pdf | PDO state concordance heatmap. |
| Auto_task3_sample_abundance.pdf | Sample abundance for MPs and states. |
| Auto_task6_hybrid_heatmap.pdf | Hybrid per cell heatmap. |
| update_19mar.tex | Weekly LaTeX presentation source. |

## Decisions Made
1. Selected key comparisons highlighting methods consistency (GSVA vs DGE and malignant vs whole).
2. Grouped PDO pan-cancer correlations and state concordance onto shared side-by-side slides to conserve presentation space while delivering comprehensive context.
3. Focused explicitly on 'Basal to Intestinal Metaplasia' pseudotime states for targeted narrative flow.

## Next Steps
- [ ] Review slides and confirm content for PPTX conversion.
- [ ] Incorporate pipeline tweaks or updated color palettes if requested.
