# Progress Update — 2 March 2026

**Project:** OAC scATLAS — Single-cell Reference Pipeline  
**Cycle:** Since commit `734e773` ("Added weekly updates flow")

---

## Context

This cycle focused on two areas: (1) quality-checking MP15 (T/NK infiltration program) for potential doublet/contamination artefacts in the epithelial compartment, and (2) implementing and comparing two systematic approaches for assigning cell states from metaprogram scores — Louvain clustering on MP UCell scores vs top-MP residual assignment. Both approaches were run on 75,348 OAC epithelial cells across all studies.

## Summary of Work

- **MP15 contamination investigation**: Examined whether epithelial cells with high MP15 (T/NK infiltration) scores co-express non-epithelial markers. Generated per-sample heatmaps for the top 3 MP15-contributing samples (Alcindor SRR27335939, Alcindor SRR27335940, Carroll PreTx) and a pan-cohort comparison of MP15-high vs MP15-low cells using z-normalised UCell scores.
- **Single-sample QC co-expression heatmap**: Generated for Alcindor_2025_SRR27335939 to visualise epithelial/T/NK marker co-expression patterns directly.
- **Cluster-based state assignment**: Louvain clustering (res 0.5/0.8/1.0) on PCA of cell-cycle-regressed, Z-normalised MP UCell scores. Identified 4 states at res 0.5: Classical/Metabolic (C1), Ciliated/Stressed (C2), Cycling (C3), Basal/EMT (C4). Generated mean heatmap, single-cell heatmap, silhouette, UMAP, and composition visualisations.
- **Top-MP state assignment**: Each cell assigned to its dominant non-cell-cycle MP by highest Z-score (threshold 0.5 → Unresolved). Identified 5 named states + Unresolved. Generated mean heatmap, single-cell heatmap, violin, and proportion plots.
- **Systematic comparison**: Confusion matrix, ARI (0.22), NMI (0.34), bootstrap stability (100× 80% subsample), Cramér's V study-bias, DEG analysis + fgsea Hallmark enrichment for biological coherence.

## Key Outputs

### New Scripts

| File | Purpose |
|------|---------|
| `analysis/Auto_MP15_investigation.R` | MP15 contamination check: z-normalised UCell scores, per-sample and pan-cohort heatmaps |
| `analysis/cell_states/Auto_states_cluster.R` (+`.sh`) | Louvain clustering on CC-regressed MP scores (PCA → FindNeighbors → FindClusters), 5 visualisations |
| `analysis/cell_states/Auto_states_topmp.R` (+`.sh`) | Top-MP residual state assignment (max Z-score, threshold 0.5), 4 visualisations |
| `analysis/cell_states/Auto_states_comparison.R` (+`.sh`) | Quantitative comparison: confusion matrix, ARI/NMI, bootstrap, Cramér's V, DEG coherence |

### New Data Outputs (in `ref_outs/`)

| File | Description |
|------|-------------|
| `Auto_cluster_states.rds` | Named vector: barcode → cluster state |
| `Auto_cluster_umap_embeddings.rds` | UMAP coordinates from cluster analysis |
| `Auto_cluster_mp_adj.rds` | CC-regressed MP score matrix (cluster) |
| `Auto_topmp_states.rds` | Named vector: barcode → top-MP state |
| `Auto_topmp_mp_adj.rds` | CC-regressed MP score matrix (top-MP) |
| `Auto_states_comparison_summary.csv` | Comparison metrics table |

### New Figures (in `ref_outs/`)

| Figure | Description |
|--------|-------------|
| `Auto_MP15_top3_*.png` (×3) | Per-sample MP15 gene expression in top 3 contributing samples |
| `Auto_MP15_allSamples_high_vs_low.png` | Pan-cohort MP15-high vs MP15-low comparison |
| `Auto_QC_Alcindor_2025_SRR27335939.png` | Single-sample co-expression QC heatmap |
| `Auto_cluster_states_mean_heatmap.pdf` | Mean Z-score MP heatmap per cluster |
| `Auto_cluster_states_singlecell_heatmap.pdf` | Single-cell level heatmap |
| `Auto_cluster_states_silhouette.pdf` | Silhouette scores per cluster |
| `Auto_cluster_states_umap.pdf` | UMAP coloured by cluster states |
| `Auto_cluster_states_composition.pdf` | Per-sample state composition |
| `Auto_topmp_states_mean_heatmap.pdf` | Mean Z-score MP heatmap per top-MP state |
| `Auto_topmp_states_singlecell_heatmap.pdf` | Single-cell level heatmap |
| `Auto_topmp_states_violin.pdf` | Violin plots of MP scores per state |
| `Auto_topmp_states_proportion.pdf` | Per-sample state proportions |
| `Auto_states_comparison_umap.pdf` | Side-by-side UMAP comparison |
| `Auto_states_confusion.pdf` | Confusion matrix between approaches |
| `Auto_states_studybias.pdf` | Cramér's V study confounding analysis |
| `Auto_states_coherence_cluster.pdf` | Hallmark pathway enrichment of cluster DEGs |

### Modified Files

| File | Change |
|------|--------|
| `AGENTS.md` | Updated with documentation for all new analysis scripts |

## Decisions Made

1. **Cell-cycle regression before state assignment**: Both approaches regress out G1/S and G2/M MP scores before computing states to avoid cycling cells dominating the analysis.
2. **Z-normalisation of UCell scores**: Per-MP Z-normalisation ensures equal weighting across MPs with different score distributions.
3. **Top-MP threshold of 0.5**: Cells with max Z-score < 0.5 labelled "Unresolved" to avoid forcing weak assignments.
4. **Resolution 0.5 for Louvain**: Selected as the primary clustering resolution (0.8 and 1.0 also computed) — yields 4 biologically interpretable states.
5. **Bootstrap at 80% subsampling**: 100 iterations of 80% random subsampling to assess partition stability.

## Comparison Results Summary

| Metric | Cluster | Top-MP |
|--------|---------|--------|
| States | 11 | 12 |
| ARI | 0.221 | 0.221 |
| NMI | 0.345 | 0.345 |
| Bootstrap ARI | 0.412 | **1.000** |
| Cramér's V | 0.065 | **0.055** |
| DEGs (top 50/state) | 150 | 50 |
| Hallmark pathways (NES>1.5) | 13 | — |

Top-MP approach is perfectly stable and slightly less study-confounded. Cluster approach finds more DEGs and enriched pathways but is less reproducible under subsampling.

## Next Steps

- [ ] Decide on final state assignment approach based on comparison metrics and biological interpretability
- [ ] Address MP15: determine whether MP15-high cells should be filtered, flagged, or retained with caveats
- [ ] Run downstream survival and clinical analysis with the selected cell states
- [ ] Extend non-malignant TME NMF analysis across remaining cell types

---

*Generated: 2 March 2026*
