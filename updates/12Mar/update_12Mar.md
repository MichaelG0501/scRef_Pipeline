# Progress Update — 12 March 2026

**Project:** scRef Pipeline — OAC single-cell RNA-seq public reference dataset analysis

---

## Context

This cycle focused on finalising Approach B cell state assignment with two modes (CC-regressed and non-regressed), running survival and clinical association analyses on TCGA-ESCA, two sensitivity analyses for survival parameters, updating the enrichment annotation pipeline with a Barretts reference, and completing an initial GeneNMF characterisation of the PDO dataset.

---

## Summary of Work

- **Approach B state assignment (reg vs noreg):** Unified script comparing CC-regressed and non-regressed MP Z-score state calls. States are highly consistent across modes; CC regression does not materially change dominant state identities. Classic\_Prolif cells carry the highest CC scores as expected.
- **Pairwise hybrid node plot:** Pairwise network of real states connected by hybrid-cell edges. Barretts–EMT is the dominant hybrid edge (~4.1%). Triplet hybrids excluded from edge rendering.
- **Unresolved cell pan-cancer subclassification:** Unresolved cells (Z < 0.5 threshold) scored against 3CA pan-cancer MPs and subclassified into 7 subclasses in both reg and noreg modes.
- **TCGA survival analysis:** Cox PH volcano plots for EAC (and ESCC) across both reg/noreg modes. One EAC state reaches nominal significance; effect directions consistent across modes.
- **Clinical associations:** Full clinical variable association analysis (treatment, stage, sample type, etc.) across both reg/noreg modes; combined PDF (10 pages).
- **Sensitivity — gene-weight threshold:** Three MP gene-selection strategies tested (min-length top-N, weight ≥ min 10th-weight, all genes). Cox volcano results are stable across all three strategies.
- **Sensitivity — TCGA sample type filter:** Compared all sample_type_code vs primary tumour only (code 01). Effect sizes and significance preserved; primary-tumour restriction confirmed appropriate.
- **Enrichment annotation — Barretts reference added:** `Barretts_Oesophagus` gene set integrated into `enrichment_annotation.R`. Adult Epithelium terms reordered by MP. Pending full rerun with updated `enrich_dev.rds`.
- **PDO GeneNMF:** GeneNMF run on PDO lines (nMP range 4–35). Optimal nMP selected by silhouette + WSS. MPs annotated against 10 databases. Jaccard self-similarity and per-sample Spearman correlation computed against scRef MPs.

---

## Key Outputs

| File | Location | Description |
|:-----|:---------|:------------|
| `Auto_states_topmpB_reg_noreg.R` | `analysis/cell_states/` | Unified Approach B state assignment, reg + noreg modes |
| `Auto_topmp_v2_reg_states_B.rds` | `ref_outs/` | Per-cell state calls — reg mode |
| `Auto_topmp_v2_noreg_states_B.rds` | `ref_outs/` | Per-cell state calls — noreg mode |
| `Auto_topmp_v2_reg_noreg_heatmap_B_cconly.pdf` | `ref_outs/` | Per-cell MP Z-score heatmap (reg p1, noreg p2) |
| `Auto_topmp_v2_reg_noreg_proportion_B_withpie.pdf` | `ref_outs/` | State proportions stacked bar + pie |
| `Auto_topmp_v2_reg_noreg_ccscore_boxplot_B.pdf` | `ref_outs/` | CC score by state boxplot |
| `Auto_states_hybrid_pairwise_nodeplot.R` | `analysis/cell_states/` | Pairwise hybrid network script |
| `Auto_topmp_v2_hybrid_pairwise_nodeplot_reg_noreg.pdf` | `ref_outs/` | Pairwise hybrid node plot |
| `Auto_states_unresolved_pan_cancer_reg_noreg.R` | `analysis/cell_states/` | Unresolved cell pan-cancer subclassification |
| `Auto_topmp_v2_reg_noreg_unresolved_heatmap.pdf` | `ref_outs/` | Unresolved subclass heatmap (reg + noreg) |
| `Auto_survival_clinical_mps_v2_reg_noreg.R` | `analysis/clinical/` | TCGA survival + clinical associations, unified reg/noreg |
| `Auto_survival_tcga_state_volcano_reg_noreg.pdf` | `ref_outs/` | Cox volcano — EAC + ESCC, reg + noreg (4 pages) |
| `Auto_clinical_assoc_topmp_v2B_reg_noreg_combined.pdf` | `ref_outs/` | Clinical associations combined PDF (10 pages) |
| `Auto_clinical_variable_plots_topmp_v2B_reg_noreg.R` | `analysis/plotting/` | Clinical variable plots script |
| `Auto_test_mp_weight_threshold.R` | `analysis/clinical/` | Sensitivity test — gene-weight threshold strategies |
| `Auto_test_mp_weight_threshold_EAC_volcanoes.pdf` | `ref_outs/` | 3-panel volcano: min-length vs weight-threshold vs all genes |
| `Auto_test_sample_type_filter.R` | `analysis/clinical/` | Sensitivity test — TCGA sample type code filter |
| `Auto_test_sample_type_filter_EAC_volcanoes.pdf` | `ref_outs/` | 2-panel volcano: all codes vs primary only |
| `enrichment_annotation.R` (modified) | `analysis/enrichment/` | Added Barretts\_Oesophagus reference; reordered Adult Epithelium terms; wider PNG output |
| `enrich_Adult_Epithelium.png` | `ref_outs/` | Reordered Adult Epithelium enrichment heatmap |
| `PDOs_enrichment_annotation.pdf` | `PDOs_outs/` | PDO MP annotation across 10 databases |
| `pdo_jaccard.pdf` / `pdo_meta_complexheatmap.pdf` | `PDOs_outs/` | PDO Jaccard self-similarity + Spearman correlation heatmaps |

---

## Decisions Made

1. **reg vs noreg:** Both modes retained for all downstream outputs. State identities are sufficiently stable that either mode is defensible; reg is the primary mode, noreg shown for completeness.
2. **Approach B threshold:** Z-score > 0.5 for dominant MP; cells below threshold in all non-CC MPs → Unresolved.
3. **Hybrid plot:** Pairwise only (2-way hybrids); >2-class hybrids excluded from edge rendering to keep the plot interpretable.
4. **Survival sensitivity:** Gene-weight threshold strategy does not materially alter Cox results; "all genes" strategy (option c) adopted as the standard. Primary-tumour-only filter (sample_type_code = "01") confirmed as appropriate for TCGA analysis.
5. **Barretts reference:** Integrated into `enrichment_annotation.R`; pending `enrich_dev.rds` rerun to propagate scores.
6. **PDO NMF:** Optimal nMP selected; full annotation completed. Global pheatmap not shown in update — only Jaccard and meta ComplexHeatmap retained as informative cross-dataset comparisons.

---

## Next Steps

- [ ] Rerun `enrichment_annotation.R` after `enrich_dev.rds` is updated to include Barretts\_Oesophagus scores
- [ ] Review Barretts enrichment heatmap once available
- [ ] Compare PDO MPs to scRef nMP=19 MPs directly (Jaccard + cosine similarity matrix)
- [ ] Run PDO state assignment using scRef state definitions (project scRef state scores onto PDOs)
- [ ] Refine unresolved cell subclassification labels based on 3CA MP identity
- [ ] Decide on final state labels for manuscript (reg mode primary)
- [ ] TCGA survival: consider multivariate Cox (adjust for stage/treatment) for key states
- [ ] Convert slides to PPTX once PDF is approved
