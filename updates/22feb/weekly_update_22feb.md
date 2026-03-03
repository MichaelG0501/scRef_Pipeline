# Weekly Progress Update — 22–28 February 2026

**Project:** Single-cell Multiomic Profiling of Intratumoural Heterogeneity in OAC  
**Repository:** `scRef_Pipeline` (Public scRNA-seq Atlas)  
**Author:** Shengbo Gong

---

## Context

This week focused on downstream analysis of the OAC scATLAS: exploring NMF metaprogramme resolutions, defining malignant cell states, extending NMF to non-malignant TME cell types, running enrichment and correlation analyses, and performing clinical validation. Additionally, infrastructure for LaTeX slide compilation on HPC was set up and comprehensive methods slides were created.

**Key working resolution:** nMP = 19 (confirmed this week after systematic evaluation of nMP 8–30).

---

## Day-by-Day Summary

### Sunday 23 Feb — NMF Resolution Exploration & Cell State Definition

**NMF resolution sweep:**
- Generated ComplexHeatmap visualisations for all resolutions nMP = 8 to 30
- Two variants per resolution: standard (top genes) and full (all MP genes)
- 46 high-resolution PNG heatmaps stored in `ref_outs/Metaprogrammes_Results/`
- Files: `geneNMF_metaprograms_nMP_{8-30}_complexheatmap.png` and `*_full.png`

**Cell state definition:**
- Explored multiple filtering strategies for assigning cells to states based on MP activation:
  - With/without cell cycle MPs
  - Different cosine similarity thresholds (0.1, 0.2)
  - Single-cell heatmap visualisation
- Iterative PDFs generated:
  - `MP_state_singlecell_heatmap.pdf` / `*_adjusted.pdf`
  - `MP_state_cc_all.pdf`, `MP_state_cc_filter.pdf`
  - `MP_state_remove_filter.pdf`, `MP_state_nofilter_cc.pdf`
  - `MP_state_02filter_nocc.pdf`
  - **`MP_state_selected.pdf`** — final chosen states
- 3CA MP comparison: `3CA_MPs_heatmap_by_state.pdf`
- **Result: 4 putative malignant cell states confirmed** (Classic Proliferative, Squamous Transition, Intestinal Metaplasia, Plastic/Tolerant)

### Monday 24 Feb — Non-Malignant TME NMF & Annotation

**TME GeneNMF pipeline** (`analysis/non_mali_nmf/`):
- Created complete GeneNMF + annotation pipeline for 7 non-malignant cell types:
  - CD4+ T cells, CD8+ T cells, NK cells, Macrophages, Fibroblasts, Endothelial cells, Plasma cells
- Scripts created:
  - `Auto_geneNMF_celltype.R` — GeneNMF analysis per cell type
  - `Auto_anno_celltype.R` — Enrichment annotation per cell type
  - `Auto_geneNMF_master.sh` — PBS master orchestrator
  - Per-cell-type PBS wrappers (`Auto_geneNMF_cd4.sh`, `Auto_geneNMF_cd8.sh`, etc.)
  - Per-cell-type annotation PBS wrappers (`Auto_anno_cd4.sh`, `Auto_anno_cd8.sh`, etc.)
  - `Auto_anno_master.sh` — annotation master orchestrator

**Per cell type outputs** (in `ref_outs/nmf_{celltype}/`):
- `tmdata_all.rds` — merged Seurat object
- `MP_outs_default.rds` — GeneNMF metaprogramme results
- `UCell_default.rds` — UCell scores
- `GO_outs.rds` — GO enrichment
- `metaprograms_heatmap.png` — ComplexHeatmap
- `vln_origident.png` — violin plots by sample
- Enrichment PNGs: `enrich_Hallmark_{ct}.png`, `enrich_GO_{ct}.png`, `enrich_MPs_3CA_{ct}.png`

**Summary outputs:**
- `ref_outs/nmf_summaries/` — combined enrichment PDFs and heatmaps for all 7 cell types
- `ref_outs/MP_metrics_summary.pdf` — MP metrics across cell types

**Developmental enrichment** (`analysis/developmental/developmental.R`):
- New script for enrichment against 6 developmental databases
- Output: `enrich_dev.rds`, per-stage RDS files in `analysis/developmental/per_stage/`

**Enrichment annotation** (`analysis/example_anno.R`, `analysis/terms_overlap.R`):
- Updated `cluster_enrich.rds` (74.5 MB) with enrichments across all 8 databases
- Annotation PDFs: `New_nMP15_anno.pdf`, `New_nMP19_anno.pdf`

**Documentation:**
- Created `AI_RULES.md` — HPC execution and file safety rules for AI agents
- Updated `AGENTS.md` knowledge base

### Tuesday 25 Feb — MP Correlation v1 & Comparative Analysis

**MP correlation v1** (`analysis/Auto_MP_correlation.R`):
- UCell scoring of all enrichment reference terms → `UCell_ref_terms.rds` (36.6 MB)
- Per-sample Spearman rank correlation between MP UCell scores and term UCell scores
- Fisher Z meta-analysis for aggregate significance
- Output: `Auto_MP_correlation_heatmaps.pdf`, `Auto_MP_correlation_results.rds`
- PBS wrapper: `analysis/Auto_MP_correlation.sh`

**nMP=19 UCell scoring** (`analysis/Auto_MP19_analysis.R`):
- UCell scoring for all 19 metaprogrammes with silhouette filtering
- Output: `Metaprogrammes_Results/UCell_nMP19_filtered.rds` (5.7 MB)
- Heatmap: `nMP19_UCellHeatmap_20k_subset.pdf` (20k cell subset)

**Comparative analyses:**
- `analysis/MP_analysis_sc.R` — single-cell MP correlation (per-sample Spearman, Fisher Z)
- `analysis/compare_pdos_sc.R` — PDO vs. primary tumour MP comparison (Jaccard + UCell correlation)
- `analysis/TCGA_data.R` — TCGA-ESCA bulk RNA-seq validation of cell state signatures

### Wednesday 26 Feb — MP Correlation v2 & Enrichment Heatmaps

**MP correlation v2** (`analysis/Auto_MP_correlation_v2.R`):
- Refined term selection for improved correlation analysis
- Two runs: default MPs and nMP=19 resolution
- Outputs:
  - `UCell_ref_terms_v2.rds` (83.6 MB), `UCell_ref_terms_v2_MP19.rds` (99.7 MB)
  - `Auto_MP_correlation_heatmaps_v2.pdf`, `Auto_MP_correlation_results_v2.rds`
  - `Auto_MP_correlation_heatmaps_v2_MP19.pdf`, `Auto_MP_correlation_results_v2_MP19.rds`

**Enrichment heatmap PNGs** (all in `ref_outs/`):
- `enrich_Hallmark.png`, `enrich_GO.png`, `enrich_MPs_3CA.png`
- `enrich_Early_Embryogenesis.png`, `enrich_Organogenesis_major.png`, `enrich_Organogenesis_sub.png`
- `enrich_Normal_Development_long.png`, `enrich_Normal_Development_short.png`
- `enrich_Adult_Epithelium.png`

### Thursday 27 Feb — Survival Analysis, Methods Slides, LaTeX Infrastructure

**Survival analysis** (`analysis/cibersort_result.R`):
- Kaplan-Meier survival across 3 feature types × multiple split methods:
  - DEG GSVA, MP gene GSVA, CIBERSORTx deconvolution
  - Median, tertile, optimal cutpoint splits
- Applied to 4 cell states and per-MP scores
- DEG filter pattern: `group_by(cluster) %>% slice_max(n=100)`
- Output: 6 survival PDFs in `ref_outs/`

**Methods LaTeX slides:**
- Created `ref_outs/Presentation/methods_slides.tex` — 454-line Beamer deck
- Covers: dataset collection, preprocessing, normalisation, QC, annotation, doublet filtering, malignant classification, NMF analysis, MP annotation, cell state definition, correlation methods
- 22 slides with figures, tables, and TikZ diagrams
- PBS compile script: `ref_outs/Presentation/compile_slides.sh`
- Compiled PDF: `ref_outs/Presentation/methods_slides.pdf`

**LaTeX infrastructure:**
- Custom `lualatex` wrapper at `~/bin/lualatex` for HPC
- Compiles via /tmp for 12× speedup (82s NFS → 7s local)
- VS Code LaTeX Workshop integration via `.vscode/settings.json`
- Full documentation: `analysis/Auto_instruction_latex.md`

### Friday 28 Feb — Documentation & Weekly Update

- Updated `AGENTS.md` with all new analysis scripts and data objects
- Created this weekly progress document

---

## Key Output Files Created This Week

### Data Objects (`.rds`)

| File | Size | Description |
|------|------|-------------|
| `ref_outs/cluster_enrich.rds` | 74.5 MB | MP enrichment across 8 databases |
| `ref_outs/UCell_ref_terms.rds` | 36.6 MB | UCell scores for reference terms (v1) |
| `ref_outs/UCell_ref_terms_v2.rds` | 83.6 MB | UCell scores for reference terms (v2) |
| `ref_outs/UCell_ref_terms_v2_MP19.rds` | 99.7 MB | UCell scores for nMP=19 terms |
| `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds` | 5.7 MB | Silhouette-filtered nMP=19 UCell scores |
| `ref_outs/Auto_MP_correlation_results.rds` | 15.6 KB | Correlation results v1 |
| `ref_outs/Auto_MP_correlation_results_v2.rds` | 31 KB | Correlation results v2 |
| `ref_outs/Auto_MP_correlation_results_v2_MP19.rds` | 72.6 KB | Correlation results v2 nMP=19 |
| `ref_outs/nmf_{celltype}/tmdata_all.rds` | 278 MB–1.1 GB | TME merged Seurat objects (×7) |
| `ref_outs/nmf_{celltype}/MP_outs_default.rds` | varies | TME GeneNMF results (×7) |
| `ref_outs/nmf_{celltype}/UCell_default.rds` | varies | TME UCell scores (×7) |
| `ref_outs/nmf_{celltype}/GO_outs.rds` | varies | TME GO enrichment (×7) |
| `analysis/developmental/enrich_dev.rds` | 175 KB | Developmental enrichment |
| `analysis/developmental/per_stage/*.rds` | varies | Per-stage enrichment (×6) |

### Visualisations

| File | Description |
|------|-------------|
| `ref_outs/Metaprogrammes_Results/*_complexheatmap*.png` (×46) | NMF resolution heatmaps |
| `ref_outs/enrich_*.png` (×9) | Enrichment heatmaps (Hallmark, GO, 3CA, developmental) |
| `ref_outs/nmf_summaries/metaprograms_heatmap_*.png` (×7) | TME cell type heatmaps |
| `ref_outs/nmf_summaries/enrich_combined_*.pdf` (×7) | TME combined enrichment |
| `ref_outs/nmf_{celltype}/vln_origident.png` (×7) | TME violin plots |
| `ref_outs/nmf_{celltype}/enrich_*.png` (~21 total) | TME per-type enrichment |
| `ref_outs/MP_state_*.pdf` (×7) | Cell state exploration PDFs |
| `ref_outs/3CA_MPs_heatmap_by_state.pdf` | 3CA comparison by state |
| `ref_outs/New_nMP15_anno.pdf`, `New_nMP19_anno.pdf` | MP annotation PDFs |
| `ref_outs/nMP19_UCellHeatmap_20k_subset.pdf` | nMP=19 UCell heatmap |
| `ref_outs/Auto_MP_correlation_heatmaps*.pdf` (×4) | Correlation heatmaps |
| `ref_outs/MP_metrics_summary.pdf` | MP metrics summary |
| `ref_outs/Presentation/methods_slides.pdf` | Methods presentation |

### Scripts Created/Modified

**Note:** "Modified" means only specific chunks or functionality were added this week — the full script is not new work from this week. "New" means the entire file was created this week.

| Script | New/Modified | What changed this week |
|--------|-------------|----------------------|
| `analysis/non_mali_nmf/Auto_geneNMF_celltype.R` | New | GeneNMF for TME cell types |
| `analysis/non_mali_nmf/Auto_anno_celltype.R` | New | Annotation for TME cell types |
| `analysis/non_mali_nmf/Auto_geneNMF_*.sh` (×8) | New | PBS wrappers |
| `analysis/non_mali_nmf/Auto_anno_*.sh` (×8) | New | PBS annotation wrappers |
| `analysis/Auto_MP19_analysis.R` | New | nMP=19 focused analysis |
| `analysis/Auto_MP_correlation.R` | Modified | Added per-sample Spearman + Fisher Z meta-analysis |
| `analysis/Auto_MP_correlation_v2.R` | New | MP correlation v2 |
| `analysis/Auto_MP_correlation.sh` | New | PBS wrapper for correlation |
| `analysis/cibersort_result.R` | Modified | Added KM survival for 3 feature types × split methods |
| `analysis/compare_pdos_sc.R` | Modified | Updated with latest MP definitions |
| `analysis/MP_analysis_sc.R` | Modified | Updated with nMP=19 definitions |
| `analysis/TCGA_data.R` | Modified | Added GSVA scoring of state-specific gene sets |
| `analysis/example_anno.R` | Modified | Added developmental databases to enrichment |
| `analysis/terms_overlap.R` | Modified | Extended to 8 databases (from original 3) |
| `analysis/developmental/developmental.R` | New | Developmental enrichment |
| `ref_outs/Presentation/methods_slides.tex` | New | Methods slides |
| `ref_outs/Presentation/compile_slides.sh` | New | PBS compile script |
| `analysis/Auto_instruction_latex.md` | New | LaTeX setup guide |
| `AI_RULES.md` | New | HPC rules for AI agents |
| `AGENTS.md` | Modified | Added analysis scripts docs, data objects, patterns |

---

## Decisions Made This Week

1. **nMP = 19 confirmed** as the working resolution after evaluating all resolutions 8–30 via ComplexHeatmap and silhouette scores.
2. **4 malignant cell states** defined: Classic Proliferative, Squamous Transition, Intestinal Metaplasia, Plastic/Tolerant — using cosine similarity to binary MP templates with cell cycle regression.
3. **Cell cycle exclusion** from state definition (regressed from correlated MPs instead).
4. **MP correlation v2** preferred over v1 due to refined term selection and coverage of nMP=19 resolution.
5. **TME NMF via GeneNMF** applied to all 7 major non-malignant cell types to characterise the tumour microenvironment.

---

## Open Questions / Next Steps

- [ ] Finalise cell state clinical correlation (pre vs. post treatment proportions)
- [ ] Integrate CosMx spatial transcriptomics data with MP scores
- [ ] In-house snRNA-seq/snATAC-seq pre/post comparison
- [ ] Expand TCGA validation (multivariable Cox models)
- [ ] TME NMF biological interpretation — which TME programs correlate with malignant states?
- [ ] Add results slides to methods deck
- [ ] Draft manuscript results figures
- [ ] Prepare for supervisor meeting

---

*This document is intended to be referenced in future sessions via `weekly_updates/` to maintain continuity across chat sessions.*
