# Progress Update — 23 March 2026

**Project:** scRef_Pipeline / OAC Reference Analysis

---

## Context
This cycle focused on finalizing cell state definitions through pan-cancer (3CA) relabeling of unresolved epithelial cells and performing robust survival and pseudotime analysis across these finalized states. We implemented batch-corrected pseudotime trajectories to ensure biological state transitions are consistent across studies.

## Summary of Work
- **Finalized States**: Updated `states_unresolved_relabel.R` to relabel "Unresolved" cells by assigning them to their top 3CA metaprograms with sufficient coverage (>50 samples, >6 studies). 
- **Expanded States**: Integrated 3CA metaprograms (e.g., Epithelial-1, Epithelial-5, EMT-and-Protein-Maturation) into the final state classification, resulting in a refined 8-state model.
- **TCGA Survival**: Implemented multi-method survival analysis comparing Continuous Cox, Median Split, and Quantile Split (Q1 vs Q4) for both metaprograms and cell states.
- **Batch-Corrected Pseudotime**: Developed `Auto_pseudotime_batch_correction.R` using Harmony and scVI to verify state trajectories are not driven by batch effects.
- **Cross-Cohort Comparison**: Analyzed survival and state proportions across EAC and ESCC cohorts, identifying distinct metaprogram and state enrichment.
- **DGE Analysis**: Performed DGE overlap analysis between defined states and Hallmark/Reactome pathways to validate state biology.

## Key Outputs
| File | Description |
| :--- | :--- |
| `states_unresolved_relabel.R` | Relabeling unresolved cells using 3CA MPs and finalizing state nomenclature. |
| `Auto_final_percell_heatmap.R` | High-resolution heatmap of finalized states across 8,000 cells. |
| `Auto_pseudotime_batch_correction.R` | Pseudotime analysis with Harmony/scVI integration to handle batch effects. |
| `Auto_topmp_v2_noreg_states_B.rds` | Updated state assignments with 3CA relabelled subclasses. |
| `Auto_task2_survival_*.pdf` | Paired continuous, median, and quantile survival comparison plots. |
| `Auto_task2_dge_overlap_analysis.pdf` | Side-by-side hallmark/pathway enrichment for state DEGs. |
| `Auto_task8_tcga_eac_escc_compare_plots.pdf` | Side-by-side comparison of state activity in EAC vs ESCC cohorts. |

## Decisions Made
1. **Unresolved Relabelling**: Chose a 3CA-based approach for unresolved cells to provide higher resolution without introducing purely study-specific clusters.
2. **Batch Correction in Pseudotime**: Opted for both Harmony and scVI to confirm that trajectory findings (root=Basal/Classic Prolif) are biologically robust.
3. **State Merges**: Merged `3CA_mp_30 Respiration 1` into 'Classic Proliferative' and combined `3CA_mp_12/17` into '3CA_EMT_and_Protein_maturation' to maintain parsimony and focus on clinical relevance.

## Next Steps
- [ ] Incorporate clinical variable associations (Treatment response, Stage) for the finalized 8-state model.
- [ ] Validate specific hybrid cell configurations in matched PDO pairs.
- [ ] Prepare finalized figures for publication-ready UMAP layouts.
- [ ] Finalize the mastersheet of all sample metadata across cohorts.
