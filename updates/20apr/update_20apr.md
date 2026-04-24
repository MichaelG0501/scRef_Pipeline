# Progress Update — 20 April 2026

**Project:** snSeq_Pipeline / PDOs_Pipeline

---

## Context
This update focuses on finalizing the resolution of unresolved cell states through pan-cancer subclassifications, evaluating state sample abundance, assessing matched pair temporal cell state distributions between snSeq and PDOs, characterizing state-specific surface markers in PDOs, and validating cross-dataset pseudobulk correlation with external 3CA references.

## Summary of Work
- Characterized the unresolved cell population by establishing 8 new state labels using a 3CA pan-cancer top-MP relabeling approach.
- Plotted sample abundance distributions categorized by top MPs and study origin, detailing dataset diversity.
- Executed matched-pair analysis between PDO and snSeq modalities to visualize cell state transitions.
- Generated expression dotplots mapping the top five prioritized surface markers for the five primary states in PDO datasets.
- Performed cross-dataset 3CA pseudobulk correlation to confirm conservation of state-specific metaprograms in PDO models.

## Key Outputs
| Output File | Description |
|---|---|
| `Auto_task4_unresolved_relabel_heatmap.pdf` | Heatmap of the newly subclassified unresolved states using pan-cancer MPs. |
| `Auto_task4_unresolved_relabel_proportion.pdf` | Relative proportion distributions of the newly assigned unresolved sub-states. |
| `Auto_task4_unresolved_relabel_cc_boxplot.pdf` | Boxplots comparing cell cycle scores across the subclassified unresolved states. |
| `sn_sample_abundance.pdf` | Breakdown of sample abundance and state representation sorted by top MPs and study. |
| `pdo_sn_matched.pdf` | Matched pair temporal comparison plots illustrating state distribution dynamics in PDOs vs snSeq. |
| `Auto_five_state_surface_marker_dotplot.pdf` | Dotplot visualizing the expression intensity and frequency of prioritized surface markers across five PDO states. |
| `Auto_3CA_pseudobulk_correlation_crossdata.pdf` | Correlation matrix comparing PDO pseudobulk profiles with 3CA external references. |

## Decisions Made
1. **Unresolved Subclassification Strategy:** Used 3CA pan-cancer MP reference matching to relabel unresolved cells without forcing them into strict scRef state definitions, allowing identification of unique but coherent signatures.
2. **Temporal Context:** Positioned the matched-pair analysis to highlight individual sample variations directly rather than relying on aggregated condition metrics.
3. **Marker Visualization Setup:** Leveraged expression dotplots to convey both the prevalence and expression level of candidate surface markers concurrently for PDO sorting applications.

## Next Steps
- [ ] Finalize and upload the slides update to the lab shared directory/Google Drive.
- [ ] Review state plasticity indicated by the matched-pair temporal dynamics in the PDO pipeline.
- [ ] Prepare functional characterization experiments using the prioritized five-state surface markers.
