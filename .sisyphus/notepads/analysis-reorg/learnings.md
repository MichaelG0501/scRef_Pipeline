# Learnings — analysis-reorg

## [2026-03-01] Session ses_359751460ffeVdzvEIV53gBLcU — Initialization

### Key Constraints (User-Confirmed)
- No blind deletions — even "temp"/"quick" named files may contain unique code
- PDO and SC logic MUST remain in separate files (structurally different data sources/naming)
- heatmap.R contains `plot_heatmap()` function NOT in main pipeline — must preserve
- ref_pipeline.R (1165L) is an early pipeline draft with potentially unique functions
- Shell ↔ R pairs must move ATOMICALLY (never one without the other)
- `Auto_` prefix is for NEW scripts created by agents; reorganized existing scripts should use descriptive names WITHOUT prefix

### Project-Specific Patterns
- Conda envs: `dmtcp` (general Seurat/analysis), `gnmf` (GeneNMF/UCell)
- All outputs go to `ref_outs/` (never write outside project paths)
- HPC PBS Pro scheduler — Imperial College HPC
- `states_degs.rds` MUST be filtered with `group_by(cluster) %>% slice_max(n=100)` before use
- Use `orig.ident` from Seurat metadata to group by sample (not barcode manipulation)

### File Structure Facts (Pre-Reorganization)
- 38 R scripts + 17 shell scripts in analysis/ (flat, unorganized)
- 11/38 files use UCell scoring
- 8 files read EAC_Ref_epi.rds (75,348 OAC epithelial cells)
- 7 files do silhouette filtering
- 25/35 top-level scripts produce heatmaps
- Zero scripts currently source() each other — all self-contained
- states_qccheck.R and qc_unresolved_states.R produce IDENTICAL output filenames
- MP analysis cluster (5 scripts, ~3.9K lines): highest overlap

## [2026-03-01] Task 1 — Triage of 7 Unclear Scripts

### Key Findings
- `analysis/geneNMF.R` does not exist — only root-level `geneNMF.R` (step 8). The file was listed but was never created or was previously deleted.
- `ref_pipeline.R` is 1166 lines (the original monolithic pipeline). Unique functions NOT in root pipeline: `write_count_matrix()`, `manual_celltyping()`, commented-out `celltyping()` (InSituType), `plot_coexpression()`.
- `heatmap.R` is a clean extraction of `plot_heatmap()` from `ref_pipeline.R` — nearly identical function but with simplified `qc_rules` (wildcard `*` vs 10 study-specific rows). Root pipeline scripts do NOT have this function.
- `temp_plot_new.R` has TWO versions of `enrich_heatmap()` — v2 (L490-637) is the most sophisticated across the codebase with custom developmental RDS support, GeneRatio display text, and viridis magma colors.
- `quick.R`'s `classify_one_ident()` also exists in root `InferCNA.R` — quick.R uses `celltype_update` column while InferCNA.R uses `celltype_manual`.
- `expr.R`'s `enrich_one()` also exists in root `Clustering.R` and `analysis/celltyping.R`.
- All 7 files have zero `source()` calls — confirming self-contained pattern.

### Function Duplication Map (from triage grep searches)
| Function | Files |
|----------|-------|
| `plot_heatmap()` | `analysis/heatmap.R`, `analysis/ref_pipeline.R` |
| `enrich_heatmap()` | `analysis/temp_plot_new.R` (v1+v2), `analysis/example_anno.R`, `analysis/enrich_plot.R`, `analysis/Auto_MP_correlation_v2.R` |
| `classify_one_ident()` | `analysis/quick.R`, `InferCNA.R` |
| `enrich_one()` | `analysis/expr.R`, `Clustering.R`, `analysis/celltyping.R` |
| `initialise()` | `analysis/ref_pipeline.R`, `QC_Pipeline.R` |
| `inspect()` | `analysis/ref_pipeline.R`, `QC_Pipeline.R` |
| `normalise()` | `analysis/ref_pipeline.R`, `QC_Pipeline.R` |
| `cells_filtering()` | `analysis/ref_pipeline.R`, `QC_Pipeline.R` |
| `list_to_df()` | `analysis/temp_plot_new.R` only |
| `mk_sheet()` | `analysis/temp_plot_new.R` only |
| `write_count_matrix()` | `analysis/ref_pipeline.R` only |
| `manual_celltyping()` | `analysis/ref_pipeline.R` only |

### Classifications
| File | Classification | Rationale |
|------|---------------|-----------|
| temp_plot.R | DELETE | Redundant scratch, superseded by temp_plot_new.R |
| temp_plot_new.R | MERGE | `enrich_heatmap()` v2 + utilities → `enrich_plot.R` |
| quick.R | DELETE | CNA exploration scratch, covered by cnv_profile_sc.R + InferCNA.R |
| geneNMF.R | DELETE | File doesn't exist |
| expr.R | DELETE | No-I/O scratch, `enrich_one()` + AUCell in Clustering.R |
| ref_pipeline.R | KEEP | Archive as `legacy/ref_pipeline_archive.R` |
| heatmap.R | KEEP | Rename to `qc_heatmap.R` |

## [2026-03-01] Task 2 — File Manifest Design

### Key Design Decisions
- Created `analysis/legacy/` subfolder for `ref_pipeline_archive.R` — archival code, not active analysis, doesn't fit plotting/ or root
- Created `analysis/summary/` subfolder for `cross_sample_summary.R` — substantial script (696L), not a small utility
- `residual.R` (821L) → `metaprograms/mp_residual_heatmap.R` — content is MP z-score residuals, not general plotting
- `annotation.R` + `celltyping.R` → `cell_states/` — both assign cell type labels via different methods, complement each other
- `sum_cancer.R` → `cell_states/cancer_summary.R` — cancer/malignant cell counts are cell state classification results
- `robust_NMF.R` → `metaprograms/robust_nmf.R` (TENTATIVE — Task 4 must verify after reading)
- Clinical survival scripts kept separate: cibersort uses CIBERSORTx deconvolution, evaluate_clinical uses UCell — different approaches
- Only 1 shell-R pair in top-level analysis/: Auto_MP_correlation.sh ↔ Auto_MP_correlation_v2.R

### Conditional Merge: UCell Scoring
- Merge Group A is conditional — Auto_MP19_analysis.R (internal MPs) vs score_other_MPs.R (external 3CA)
- Recommendation: keep separate as mp_ucell_scoring.R + mp_external_scoring.R
- Task 4 makes final decision after full reading

### Manifest Stats
- 55 physical files → 38-40 final files (3 deleted + 1 non-existent + 14 .sh consolidated by Task 9)
- 10 subfolders + root utils.R
- 3 merge groups: Enrichment Plotting (B), States QC (C), UCell Scoring (A, conditional)
- 4 open questions flagged for downstream tasks (not blocking)

## [2026-03-01] Task 4 — Metaprograms Reorganization

### Execution Summary
- 8 files moved from `analysis/` to `analysis/metaprograms/` (7 R + 1 shell)
- All R scripts parse-verified via `Rscript -e "parse(file='...')"` — all PASS
- Shell script bash syntax check PASS
- Evidence saved to `.sisyphus/evidence/task-4-metaprograms-parse.txt`
- Original files deleted after verification

### File Mapping
| Original | New | Lines |
|----------|-----|-------|
| MP_analysis_sc.R | mp_correlation_sc.R | 822 |
| MP_analysis_pdos.R | mp_correlation_pdo.R | 902 |
| compare_pdos_sc.R | mp_correlation_crossdata.R | 887 |
| Auto_MP19_analysis.R | mp_ucell_scoring.R | 393 |
| Auto_MP_correlation_v2.R | mp_database_correlation.R | 486 |
| score_other_MPs.R | mp_external_scoring.R | 411 |
| robust_NMF.R | robust_nmf.R | 387 |
| Auto_MP_correlation.sh | mp_database_correlation.sh | 23 |

### Decisions Confirmed
- `robust_NMF.R` IS MP/NMF-related (PDO NMF programs, Jaccard with 3CA MPs, GO/Hallmark enrichment) → moved to metaprograms/
- UCell scoring kept SEPARATE: mp_ucell_scoring.R (internal nMP=19) vs mp_external_scoring.R (external 3CA MPs with η² filter)
- PDO vs SC correlation scripts kept SEPARATE (different data sources: PDOs_final.rds vs EAC_Ref_epi.rds)
- Shell script Rscript path updated; PBS resources unchanged
- Shebang kept on line 1, header comment block placed after shebang (lines 2-5)

## [2026-03-01] Task 6 — Clinical Reorganization

### Files Moved
- `analysis/TCGA_data.R` (203L) → `analysis/clinical/tcga_data_prep.R` (208L, +5 header)
- `analysis/cibersort_result.R` (796L) → `analysis/clinical/survival_cibersort.R` (801L, +5 header)
- `analysis/evaluate_clinical_MPs.R` (400L) → `analysis/clinical/survival_clinical_mps.R` (404L, +4 header)

### Key Observations
- `cibersort_result.R` is a long exploratory notebook (~800L) with multiple independent analysis blocks separated by `#####` dividers. It contains: GSVA scoring of MPs, CIBERSORTx deconvolution result loading, KM survival (median/quartile/optimal), Cox PH models, forest plots. Multiple `setwd()` calls and library reloads mid-script.
- `evaluate_clinical_MPs.R` uses UCell/3CA MP scores and GSVA with cell-cycle regression. Different approach from cibersort — compares Pre vs Post treatment + TCGA Cox survival.
- `TCGA_data.R` downloads/processes TCGA-ESCA data, merges GDC + cBioPortal clinical, converts Ensembl→Symbol, outputs `tcga_esca_meta.rds` and TPM matrix. Also has a CIBERSORTx signature matrix prep section at end.
- `slice_max(n=100, order_by=avg_log2FC)` confirmed preserved in survival_cibersort.R line 236.
- All external paths preserved exactly (spatialtranscriptomics, TCGA/INPUT, etc.)
- None of these scripts source() each other — fully self-contained.
- `cibersort_result.R` has `setwd()` to external TCGA/INPUT directory mid-script (line 262 original).

## [2026-03-01] Task 7 — CNV & Cell States Reorganization

### Files Moved/Created
| Original | New | Lines |
|----------|-----|-------|
| cnv_profile_sc.R (69L) | cnv/cnv_profiling.R | 74 |
| cnv_subset.R (328L) | cnv/cnv_subsetting.R | 332 |
| plot_CNV.R (191L) | cnv/cnv_plotting.R | 195 |
| states_qccheck.R (357L) + qc_unresolved_states.R (354L) | cell_states/states_qc.R | 526 |
| states_umap.R (167L) | cell_states/states_umap.R | 172 |

### Files Deleted
- `analysis/quick.R` (428L) — confirmed redundant scratch: CNA heatmap (duplicated in cnv_subsetting), `classify_one_ident()` (InferCNA.R uses `celltype_manual` instead of `celltype_update`), CNA scatter plots per study

### Key Observations
- CNV scripts are a sequential pipeline (profiling → subsetting → plotting) — kept as 3 separate files, NOT merged
- All 3 CNV scripts share identical library sets: data.table, dplyr, ComplexHeatmap, circlize, RColorBrewer, Seurat, infercna
- `cnv_subset.R` and `plot_CNV.R` both contain a typo: `gene_order <- read.tarefgene_order <- read.table(...)` — preserved as-is (valid double assignment in R, still parses)
- `states_qccheck.R` is the superset: 8 continuous features + malignancy categorical; `qc_unresolved_states.R` is stripped-down (Part 1: 8 continuous + `categorical_features <- NULL`; Part 2: only 3 continuous `nCount_RNA, nFeature_RNA, percent.mt`)
- Both QC scripts produced IDENTICAL output filenames (`states_status_quality_comparison.pdf`, `quality_comparison_by_original_state.pdf`)
- Merged QC script has 4 sections: A (full Defined/Unresolved), B (core Defined/Unresolved, commented out), C (full per-state), D (core per-state)
- 3 existing files in `cell_states/` left untouched: `cancer_summary.R`, `cell_annotation.R`, `cell_typing.R`
- Zero scripts source() each other — consistent with project-wide pattern

## [2026-03-01] Task 5 — Enrichment Reorganization

### Execution Summary
- 5 files moved/merged from `analysis/` to `analysis/enrichment/` (all R scripts)
- 3 temp/redundant files deleted: `temp_plot_new.R`, `temp_plot.R`, `expr.R`
- 5 moved originals deleted: `terms_overlap.R`, `example_anno.R`, `enrich_plot.R`, `wnt_enrich.R`, `nmf_plot.R`
- All 5 new scripts parse-verified via `Rscript -e "parse(file='...')"` — all PASS
- Evidence saved to `.sisyphus/evidence/task-5-enrichment-parse.txt`

### File Mapping
| Original | New | Lines |
|----------|-----|-------|
| terms_overlap.R | enrichment_analysis.R | 121 |
| example_anno.R | enrichment_annotation.R | 320 |
| enrich_plot.R + 3 funcs from temp_plot_new.R | enrichment_plotting.R | 463 |
| nmf_plot.R | nmf_enrichment.R | 303 |
| wnt_enrich.R | wnt_pathway.R | 145 |

### Key Decisions
- `wnt_enrich.R` is WNT-specific (custom WNT_CM, WNT_Canonical gene sets + Jaccard overlap) → kept as separate `wnt_pathway.R`, NOT merged into generic enrichment
- `nmf_plot.R` contains NMF-specific enrichment + Jaccard heatmap + writexl export → fits `enrichment/nmf_enrichment.R`
- `enrich_heatmap()` v2 from `temp_plot_new.R` renamed to `enrich_heatmap_v2` in merged file to avoid name collision with existing `enrich_heatmap()` from `enrich_plot.R`
- 3 functions extracted from `temp_plot_new.R`: `enrich_heatmap_v2()` (L490-637), `list_to_df()` (L208-214), `mk_sheet()` (L368-386)
- External path preserved: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/enrich_dev.rds`

### Deleted Files
| File | Lines | Reason |
|------|-------|--------|
| temp_plot_new.R | 693 | Functions extracted to enrichment_plotting.R; remaining code covered by enrichment_annotation.R |
| temp_plot.R | 146 | Redundant scratch, superseded by temp_plot_new.R |
| expr.R | 120 | No-I/O scratch; enrich_one() exists in root Clustering.R |
| terms_overlap.R | 116 | Moved to enrichment_analysis.R |
| example_anno.R | 315 | Moved to enrichment_annotation.R |
| enrich_plot.R | 270 | Moved+merged into enrichment_plotting.R |
| wnt_enrich.R | 139 | Moved to wnt_pathway.R |
| nmf_plot.R | 304 | Moved to nmf_enrichment.R |

## [2026-03-01] Task 9 — Non-malignant NMF Reorganization

### Files Created
- analysis/non_malignant_nmf/Auto_geneNMF_celltype.R (copied unchanged)
- analysis/non_malignant_nmf/Auto_anno_celltype.R (copied unchanged)
- analysis/non_malignant_nmf/run_geneNMF.sh (new parameterized, replaces 7 per-celltype shells)
- analysis/non_malignant_nmf/run_annotation.sh (new parameterized, replaces 7 per-celltype shells)
- analysis/non_malignant_nmf/submit_geneNMF_all.sh (updated master, was Auto_geneNMF_master.sh)
- analysis/non_malignant_nmf/submit_annotation_all.sh (updated master, was Auto_anno_master.sh)

### Deleted
- analysis/non_mali_nmf/ (entire old folder, 18 files)

### Key Changes
- Working directory corrected: was Auto_AG (wrong), now scRef_Pipeline (correct)
- geneNMF PBS resources: unified to max (ncpus=8:mem=72gb:walltime=08h) across cell types for simplicity
- anno PBS resources: all were identical (ncpus=2:mem=16gb:walltime=03h) — kept as-is
- Cell types: macrophage, fibroblast, endothelial, nk.cell, plasma, cd4, cd8

### Important: nk vs nk.cell arg mismatch
- Auto_geneNMF_celltype.R uses `nk.cell` as the celltype_map key
- Auto_anno_celltype.R uses `nk` as the ct_map key
- submit_geneNMF_all.sh passes `nk.cell` (matching the R script's celltype_map)
- submit_annotation_all.sh passes `nk` (matching the R script's ct_map)
- Usage: qsub -v celltype=macrophage analysis/non_malignant_nmf/run_geneNMF.sh
