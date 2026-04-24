# Progress Update — 23 April 2026

**Project:** scRef_Pipeline / snSeq_Pipeline / PDOs_Pipeline

---

## Context
This cycle focused on consolidating a cross-project results deck around TME interaction evidence, cancer-linked shared gene relationships, unresolved pan-cancer assignment in snRNA-seq, matched snRNA-seq/PDO comparisons, treated-vs-untreated PDO response summaries, and marker-check simulation outputs for state identification.

## Summary of Work
- Structured the deck in a fixed seven-part backbone so results flow from TME interaction evidence to cross-platform PDO/snRNA-seq interpretation and simulation-based marker validation.
- Selected the four TME interaction dot-map outputs (one per mode) and restricted to page 1 to keep the comparison focused on significant positive interaction context.
- Added the updated shared TME-MP cancer relationship screenshot and paired it with the unresolved pan-cancer snRNA-seq barplot output.
- Included the updated matched snRNA-seq/PDO comparison pages, the full six-page treated-vs-untreated PDO summary PDF, and the marker-check simulation sequence (time-series → scatter log2FC → cross-correlation → marker detection dotplot).
- Prepared the final one-page ongoing-task slide using the same concise sub-bullet style as prior weekly updates.

## Key Outputs
| Output File | Description |
|---|---|
| `ref_outs/non_malignant_mp_correlations/01_cancer_mps_cross_only/Auto_celltype_interaction_dotmap.pdf` | TME interaction map (mode 1), page 1 used for presentation. |
| `ref_outs/non_malignant_mp_correlations/02_cancer_mps_cross_and_within/Auto_celltype_interaction_dotmap.pdf` | TME interaction map (mode 2), page 1 used for presentation. |
| `ref_outs/non_malignant_mp_correlations/03_cancer_states_cross_only/Auto_celltype_interaction_dotmap.pdf` | TME interaction map (mode 3), page 1 used for presentation. |
| `ref_outs/non_malignant_mp_correlations/04_cancer_states_cross_and_within/Auto_celltype_interaction_dotmap.pdf` | TME interaction map (mode 4), page 1 used for presentation. |
| `ref_outs/non_malignant_mp_correlations/tme_mps_gene_list.png` | Significant relationship with cancer screenshot figure. |
| `sn_outs/Auto_topmp_v2_noreg_unresolved_pan_cancer_barplot.pdf` | snRNA-seq unresolved pan-cancer assignment bar plot (bar plot only included). |
| `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.pdf` | snRNA-seq/PDO matched pair comparison (selected pages). |
| `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_presentation_final.pdf` | Six-page treated-vs-untreated PDO summary used as six separate slides. |
| `PDOs_outs/Auto_marker_selection_simulation/Auto_time_series_simulation_lineplots.pdf` | Time-series marker simulation line plots (page 1 and page 2). |
| `PDOs_outs/Auto_marker_selection_simulation/Auto_qpcr_scatter_log2fc.pdf` | Marker simulation scatter log2 fold-change summary. |
| `PDOs_outs/Auto_marker_selection_simulation/Auto_qpcr_cross_correlation_heatmap.pdf` | Marker simulation cross-correlation heatmap summary. |
| `Auto_marker_detection_by_state_dotplot.pdf` | Marker detection by state dotplot used as final simulation panel. |

## Decisions Made
1. **Section ordering:** Prioritized biological narrative continuity (TME evidence first, then cancer-linked screenshot, then pan-cancer barplot, then cross-platform matching and treatment response, then simulation diagnostics).
2. **Slide density control:** Used one figure per slide for all heavy visual outputs to preserve legibility and keep sentence style consistent with prior updates.
3. **Marker simulation ordering:** Used time-series first (both pages), then scatter log2FC, then cross-correlation, then marker detection dotplot to progress from dynamics to relationship-level and detection-level validation.

## Next Steps
- [ ] Confirm final slide text/titles and whether any panel requires page substitution before PPTX conversion.
- [ ] Convert the approved PDF slide deck to PPTX and upload the final file to Google Drive.
- [ ] Reuse the new reusable backbone-prompt template for the next update cycle and only swap asset paths/page rules.
