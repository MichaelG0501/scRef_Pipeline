# Progress Update — 9 March 2026

**Project:** OAC scATLAS — Single-cell Reference Pipeline  
**Cycle:** Since previous update cycle (post 2 March outputs)

---

## Context

This cycle focused on finalising state-definition outputs from `Auto_states_topmp_v2.R` (Approach A and B), extending hybrid-state decomposition with `Auto_states_topmp_v2_hybridB.R`, and integrating survival/clinical associations from `Auto_survival_clinical_mps_v2.R` with cohort-separated reporting. For presentation, MP15 investigation outputs were intentionally excluded (already presented), and all clinical/survival visuals were constrained to the requested EAC/ESCC split with an EAC-focused clinical panel.

## Summary of Work

- Added an overview slide introducing the two state-definition approaches (A vs B) with bullet-point explanations of the hybrid classification rules.
- Completed side-by-side state-definition reporting for two assignment approaches (A vs B), with per-cell heatmaps shown one approach per page and paired proportion comparison.
- Generated and summarised hybrid-only follow-up outputs (pairwise and multi-hybrid subclasses), including per-cell heatmap and combined mean/proportion view.
- Produced cohort-separated TCGA Cox volcano outputs (`EAC`, `ESCC`) and excluded the combined cohort panel from slides.
- Selected KM display as an EAC example with explicit note that only one state was significant.
- Restricted clinical-correlation visualisation to EAC panel only.
- Added one UMAP example by extracting only slide 3 from `Auto_topmp_v2_realstates_plushybrid_umap_top12.pdf`.

## Key Outputs

### New Scripts

| File | Purpose |
|------|---------|
| `analysis/cell_states/Auto_states_topmp_v2.R` | Real-state assignment in two approaches (A/B) and corresponding heatmap/proportion outputs |
| `analysis/cell_states/Auto_states_topmp_v2_hybridB.R` | Hybrid-cell subtype decomposition and hybrid-focused visual summaries |
| `analysis/clinical/Auto_survival_clinical_mps_v2.R` | Cohort-split TCGA Cox volcano, EAC/ESCC KM outputs, and clinical association analysis |

### Machine-Readable Summaries (this cycle)

| File | Description |
|------|-------------|
| `updates/09Mar/summaries/Auto_topmp_v2_summary.csv` | State counts/proportions for Approach A and B |
| `updates/09Mar/summaries/Auto_topmp_v2_hybridB_summary.csv` | Hybrid subtype counts/proportions (10 pairwise + 1 multi-hybrid class) |
| `updates/09Mar/summaries/Auto_survival_clinical_mps_v2_summary.csv` | Cohort sizes and counts of significant survival/clinical associations |

### Figures Used in 09Mar Slide Deck (11 slides)

| Slide | Figure | Source |
|-------|--------|--------|
| 1 | Overview: Two Approaches to State Definition | (text slide) |
| 2 | `stateA_percell_heatmap.pdf` | `ref_outs/Auto_topmp_v2_heatmap_A.pdf` |
| 3 | `stateB_percell_heatmap.pdf` | `ref_outs/Auto_topmp_v2_heatmap_B.pdf` |
| 4 | Mean heatmap A and B (side-by-side) | pages 1-2 of `ref_outs/Auto_topmp_v2_mean_heatmap.pdf` |
| 5 | Proportions A and B (side-by-side) | `ref_outs/Auto_topmp_v2_proportion_A.pdf`, `Auto_topmp_v2_proportion_B.pdf` |
| 6 | UMAP example | page 3 from `ref_outs/Auto_topmp_v2_realstates_plushybrid_umap_top12.pdf` |
| 7 | `hybrid_percell_heatmap.pdf` | `ref_outs/Auto_topmp_v2_hybridB_heatmap.pdf` |
| 8 | Hybrid mean heatmap + proportions (side-by-side) | `ref_outs/Auto_topmp_v2_hybridB_mean_heatmap.pdf`, `Auto_topmp_v2_hybridB_proportion.pdf` |
| 9 | Volcano EAC + ESCC (side-by-side) | `ref_outs/Auto_survival_tcga_volcano_EAC.pdf`, `Auto_survival_tcga_volcano_ESCC.pdf` |
| 10 | KM EAC Secretory-related (page 5) | `ref_outs/Auto_survival_tcga_state_km_EAC.pdf` (page 5) |
| 11 | Clinical correlation EAC | `ref_outs/Auto_clinical_assoc_EAC.pdf` |

## Decisions Made

1. **Overview slide added:** Brief bullet-point introduction of both approaches (A vs B) for context, matching style from previous updates.
2. **State-definition ordering:** per-cell heatmap for Approach A presented first, then Approach B on separate pages before aggregate comparisons.
3. **Approach comparison shown where most informative:** proportion plots are presented as A vs B side-by-side; mean heatmap now extracted as separate pages and shown side-by-side.
4. **Hybrid section:** hybrid per-cell heatmap kept standalone; hybrid mean heatmap and proportion combined in one slide.
5. **Survival visual policy:** volcano shown as EAC and ESCC side-by-side only; combined cohort omitted.
6. **KM emphasis:** only EAC KM panel is shown as an example (page 5 of the PDF), with explicit note that one state is significant.
7. **Clinical panel scope:** only EAC clinical-correlation figure is included in slides.
8. **MP15 policy for this update:** MP15 investigation results excluded from deck and markdown interpretation details (already presented previously).

## Next Steps

- [ ] Finalise interpretation text for the selected final state-definition approach (A/B/hybrid integration) in manuscript notes.
- [ ] Validate whether the significant KM state remains robust under alternative split thresholds.
- [ ] Prioritise EAC-specific clinical covariates with strongest MP/state associations for focused follow-up analyses.

---

*Generated: 9 March 2026*
