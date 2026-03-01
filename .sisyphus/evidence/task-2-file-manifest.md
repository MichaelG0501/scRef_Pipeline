# Task 2 — File Manifest

Generated: 2026-03-01
Files accounted for: 55 (38 R + 17 sh)

---

## Subfolder Structure

```
analysis/
├── utils.R                              # [NEW] Shared utilities (Task 3 skeleton, Task 11 populates)
├── metaprograms/                        # 7 files (6 R + 1 sh)
│   ├── mp_correlation_sc.R              #   SC MP correlation (Spearman + Fisher Z)
│   ├── mp_correlation_pdo.R             #   PDO MP correlation
│   ├── mp_correlation_crossdata.R       #   SC vs PDO Jaccard comparison
│   ├── mp_ucell_scoring.R              #   UCell scoring + silhouette filtering
│   ├── mp_database_correlation.R        #   MP-to-database-term correlation
│   ├── mp_database_correlation.sh       #   PBS wrapper for above
│   └── mp_residual_heatmap.R            #   MP residual/z-score heatmaps
├── enrichment/                          # 5 files (5 R)
│   ├── enrichment_analysis.R            #   Gene overlap enrichment computation
│   ├── enrichment_plotting.R            #   Enrichment heatmap visualization (merged)
│   ├── enrichment_annotation.R          #   MP annotation vs databases
│   ├── wnt_pathway.R                    #   WNT-specific enrichment
│   └── nmf_enrichment.R                 #   NMF result + enrichment plotting
├── clinical/                            # 3 files (3 R)
│   ├── tcga_data_prep.R                 #   TCGA download/processing
│   ├── survival_cibersort.R             #   KM survival (CIBERSORTx, GSVA, MP gene)
│   └── survival_clinical_mps.R          #   Clinical MP evaluation + survival
├── cnv/                                 # 3 files (3 R)
│   ├── cnv_profiling.R                  #   CNV profile computation
│   ├── cnv_subsetting.R                 #   CNV data subsetting
│   └── cnv_plotting.R                   #   CNV visualization (PDOs)
├── cell_states/                         # 5 files (5 R)
│   ├── states_qc.R                      #   State QC checks (merged from 2 scripts)
│   ├── states_umap.R                    #   State UMAP visualization
│   ├── cell_annotation.R               #   Cell type annotation (Seurat-based)
│   ├── cell_typing.R                    #   Cell type assignment (marker-based)
│   └── cancer_summary.R                 #   Cancer cell summary statistics
├── plotting/                            # 4 files (4 R)
│   ├── publication_umap.R               #   Publication-quality UMAPs
│   ├── gene_expression_heatmap.R        #   Gene expression comparison heatmaps
│   ├── clinical_variable_plots.R        #   Clinical variable by study plots
│   └── qc_heatmap.R                     #   Per-sample marker expression QC heatmap
├── summary/                             # 1 file (1 R)
│   └── cross_sample_summary.R           #   Cross-sample summary statistics
├── legacy/                              # 1 file (1 R)
│   └── ref_pipeline_archive.R           #   Original monolithic pipeline (archival)
├── non_malignant_nmf/                   # 18 files → 6 files after Task 9 (2 R + 4 sh)
│   ├── [Task 9 handles internals — folder rename from non_mali_nmf/]
│   └── [14 per-cell-type .sh → 2 parameterized .sh + 2 updated masters]
└── developmental/                       # 1 file (1 R) — already organized
    └── developmental.R                  #   Embryogenesis marker mappings
```

**File counts by subfolder (post-reorganization, excluding non_mali_nmf Task 9 changes):**
- metaprograms/: 7 (6R + 1sh)
- enrichment/: 5 (5R)
- clinical/: 3 (3R)
- cnv/: 3 (3R)
- cell_states/: 5 (5R)
- plotting/: 4 (4R)
- summary/: 1 (1R)
- legacy/: 1 (1R)
- non_malignant_nmf/: 18 → 6 (2R + 4sh, Task 9)
- developmental/: 1 (1R)
- root (utils.R): 1 (1R, NEW)
- **Total destination files**: 37 moved/renamed + 1 new = 38 (from 55 originals; 4 deleted, 2 merged, 18 in non_mali_nmf handled by Task 9 but all accounted)

---

## File-by-File Manifest

### DELETE (4 files)

| # | Original Path | Lines | Reason | Source |
|---|---|---|---|---|
| 1 | `analysis/temp_plot.R` | 146 | Redundant scratch; no file I/O, all code superseded by temp_plot_new.R | Task 1 |
| 2 | `analysis/quick.R` | 428 | Exploratory CNA scratch; `classify_one_ident()` in root InferCNA.R, CNV viz in cnv_profile_sc.R + plot_CNV.R | Task 1 |
| 3 | `analysis/expr.R` | 120 | No-I/O scratch; `enrich_one()` exists in root Clustering.R, AUCell workflow in Clustering.R | Task 1 |
| 4 | `analysis/geneNMF.R` | N/A | File does not exist (only root-level geneNMF.R exists) | Task 1 |

### MOVE — Content unchanged, renamed/relocated (28 files)

| # | Original Path | Final Path | Notes |
|---|---|---|---|
| 5 | `analysis/MP_analysis_sc.R` | `analysis/metaprograms/mp_correlation_sc.R` | SC MP Spearman + Fisher Z + ComplexHeatmap (827L) |
| 6 | `analysis/MP_analysis_pdos.R` | `analysis/metaprograms/mp_correlation_pdo.R` | PDO version of MP correlation (898L) |
| 7 | `analysis/compare_pdos_sc.R` | `analysis/metaprograms/mp_correlation_crossdata.R` | SC vs PDO Jaccard + pheatmap (911L) |
| 8 | `analysis/Auto_MP_correlation_v2.R` | `analysis/metaprograms/mp_database_correlation.R` | MP vs database terms correlation (482L) |
| 9 | `analysis/Auto_MP_correlation.sh` | `analysis/metaprograms/mp_database_correlation.sh` | PBS wrapper: 8 CPU, 128GB, 4h. **Must update Rscript path inside to `analysis/metaprograms/mp_database_correlation.R`** |
| 10 | `analysis/Auto_MP19_analysis.R` | `analysis/metaprograms/mp_ucell_scoring.R` | UCell scoring nMP=19 + silhouette + Jaccard (389L) |
| 11 | `analysis/score_other_MPs.R` | `analysis/metaprograms/mp_ucell_scoring.R` | **SEE MERGE GROUP A** — scores with external 3CA MPs. If Task 4 determines content is too different from Auto_MP19_analysis.R, keep as separate `mp_external_scoring.R` |
| 12 | `analysis/residual.R` | `analysis/metaprograms/mp_residual_heatmap.R` | MP heatmaps with z-score residuals (821L). Libraries: ComplexHeatmap, proxy, Seurat. Produces MP_heatmap_states_subset_z_residual.pdf. Content is MP-centric, not general plotting. |
| 13 | `analysis/terms_overlap.R` | `analysis/enrichment/enrichment_analysis.R` | Gene overlap enrichment → cluster_enrich.rds (116L) |
| 14 | `analysis/enrich_plot.R` | `analysis/enrichment/enrichment_plotting.R` | Enrichment heatmap + developmental data (270L). **Reads external path EAC_Ref_all/00_merged/developmental/ — must preserve.** Also receives merged content from temp_plot_new.R (see MERGE GROUP B) |
| 15 | `analysis/example_anno.R` | `analysis/enrichment/enrichment_annotation.R` | MP annotation vs Hallmark/GO/3CA/Pan-Cancer (315L). 9 enrichment heatmap PNGs |
| 16 | `analysis/wnt_enrich.R` | `analysis/enrichment/wnt_pathway.R` | WNT pathway enrichment (139L) |
| 17 | `analysis/nmf_plot.R` | `analysis/enrichment/nmf_enrichment.R` | NMF results + enrichment plotting via clusterProfiler (304L) |
| 18 | `analysis/TCGA_data.R` | `analysis/clinical/tcga_data_prep.R` | TCGA data download/processing (208L) |
| 19 | `analysis/cibersort_result.R` | `analysis/clinical/survival_cibersort.R` | KM survival: CIBERSORTx + GSVA + MP gene (963L) |
| 20 | `analysis/evaluate_clinical_MPs.R` | `analysis/clinical/survival_clinical_mps.R` | Clinical MP evaluation + survival (400L) |
| 21 | `analysis/cnv_profile_sc.R` | `analysis/cnv/cnv_profiling.R` | CNV profile computation (68L) |
| 22 | `analysis/cnv_subset.R` | `analysis/cnv/cnv_subsetting.R` | CNV data subsetting (328L) |
| 23 | `analysis/plot_CNV.R` | `analysis/cnv/cnv_plotting.R` | PDOs CNV visualization (191L) |
| 24 | `analysis/states_umap.R` | `analysis/cell_states/states_umap.R` | State UMAP → columnar_states_UMAP.pdf (166L) |
| 25 | `analysis/annotation.R` | `analysis/cell_states/cell_annotation.R` | Cell annotation, Seurat-based (508L) |
| 26 | `analysis/celltyping.R` | `analysis/cell_states/cell_typing.R` | Cell type assignment via markers (198L). Note: contains `enrich_one()` also in Clustering.R — inline copy OK for standalone use |
| 27 | `analysis/sum_cancer.R` | `analysis/cell_states/cancer_summary.R` | Cancer cell summary (183L). Thematically fits cell_states — it summarizes malignant/non-malignant counts |
| 28 | `analysis/beaut_umap.R` | `analysis/plotting/publication_umap.R` | Publication-quality UMAPs (428L) |
| 29 | `analysis/gene_expr_compare.R` | `analysis/plotting/gene_expression_heatmap.R` | Gene expression comparison heatmap (344L). Produces Heatmap_Matched_Colors_Dediff.pdf |
| 30 | `analysis/plot_clinical.R` | `analysis/plotting/clinical_variable_plots.R` | Clinical variable by study (342L). Produces CellState_PerVariable_ByStudy.pdf |
| 31 | `analysis/heatmap.R` | `analysis/plotting/qc_heatmap.R` | Per-sample marker expression QC heatmap (407L). Contains unique `plot_heatmap()` function |
| 32 | `analysis/summary.R` | `analysis/summary/cross_sample_summary.R` | Cross-sample summary statistics (696L). Reads EAC_Ref_filtered.rds |
| 33 | `analysis/ref_pipeline.R` | `analysis/legacy/ref_pipeline_archive.R` | Original monolithic pipeline (1166L). Archival — unique functions: `write_count_matrix()`, `manual_celltyping()`, commented InSituType/coexpression code |
| 34 | `analysis/robust_NMF.R` | `analysis/metaprograms/robust_nmf.R` | NMF-related analysis (383L). Tentative assignment to metaprograms — **Task 4 must verify after reading** |
| 35 | `analysis/developmental/developmental.R` | `analysis/developmental/developmental.R` | No change — already correctly organized (403L). Task 10 may rename if warranted |

### MERGE (2 merge targets, consuming 3 original files → 2 targets)

| # | Original Path | Target Path | Content to Preserve |
|---|---|---|---|
| 36 | `analysis/temp_plot_new.R` | → `analysis/enrichment/enrichment_plotting.R` | **MERGE GROUP B**: `enrich_heatmap()` v2 (L490-637) + `list_to_df()` (L208-214) + `mk_sheet()` (L368-386). See detailed merge group below. |
| 37 | `analysis/states_qccheck.R` | → `analysis/cell_states/states_qc.R` | **MERGE GROUP C**: Full content (357L) — base script for merged QC |
| 38 | `analysis/qc_unresolved_states.R` | → `analysis/cell_states/states_qc.R` | **MERGE GROUP C**: Full content (354L) — merge with states_qccheck.R. Both produce identical output filenames. |

### Non-Malignant NMF — Folder Rename + Task 9 Internals (18 files)

| # | Original Path | Status |
|---|---|---|
| 39 | `analysis/non_mali_nmf/Auto_geneNMF_celltype.R` | Task 9 handles — folder rename to `analysis/non_malignant_nmf/` |
| 40 | `analysis/non_mali_nmf/Auto_anno_celltype.R` | Task 9 handles — folder rename to `analysis/non_malignant_nmf/` |
| 41 | `analysis/non_mali_nmf/Auto_geneNMF_master.sh` | Task 9 handles — rename + update to use parameterized runner |
| 42 | `analysis/non_mali_nmf/Auto_anno_master.sh` | Task 9 handles — rename + update to use parameterized runner |
| 43 | `analysis/non_mali_nmf/Auto_geneNMF_cd4.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 44 | `analysis/non_mali_nmf/Auto_geneNMF_cd8.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 45 | `analysis/non_mali_nmf/Auto_geneNMF_endothelial.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 46 | `analysis/non_mali_nmf/Auto_geneNMF_fibroblast.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 47 | `analysis/non_mali_nmf/Auto_geneNMF_macrophage.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 48 | `analysis/non_mali_nmf/Auto_geneNMF_nk.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 49 | `analysis/non_mali_nmf/Auto_geneNMF_plasma.sh` | Task 9 handles — replaced by parameterized `run_geneNMF.sh` |
| 50 | `analysis/non_mali_nmf/Auto_anno_cd4.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 51 | `analysis/non_mali_nmf/Auto_anno_cd8.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 52 | `analysis/non_mali_nmf/Auto_anno_endothelial.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 53 | `analysis/non_mali_nmf/Auto_anno_fibroblast.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 54 | `analysis/non_mali_nmf/Auto_anno_macrophage.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 55 | `analysis/non_mali_nmf/Auto_anno_nk.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |
| 56 | `analysis/non_mali_nmf/Auto_anno_plasma.sh` | Task 9 handles — replaced by parameterized `run_annotation.sh` |

**Note**: File #4 (geneNMF.R) does not physically exist but was in the original tracking list. Including it brings the logical file count to 55 original entries (38 R + 17 sh), matching the `find` count of 55 actual files (since geneNMF.R is not on disk). The 55 physical files from `find` are accounted for by rows 1-3 and 5-56 (55 rows for real files). Row 4 is the non-existent geneNMF.R — a no-op deletion.

---

## Merge Groups (Detailed)

### Merge Group A: MP UCell Scoring (CONDITIONAL)

**Decision point for Task 4**: `Auto_MP19_analysis.R` (389L) and `score_other_MPs.R` (407L) both do UCell scoring but for different gene sets:
- `Auto_MP19_analysis.R` scores with the nMP=19 internal metaprograms + silhouette filtering
- `score_other_MPs.R` scores with external 3CA MP gene sets

**Option 1 (Merge)**: Combine into single `mp_ucell_scoring.R` with clearly separated sections for internal vs external scoring.
**Option 2 (Keep separate)**: `mp_ucell_scoring.R` (internal nMP=19) + `mp_external_scoring.R` (3CA external).

**Recommendation**: Keep separate (Option 2) — internal and external scoring serve different analytical purposes and produce different output files. Task 4 makes the final call after full reading.

If Option 2: row #11 changes to `analysis/metaprograms/mp_external_scoring.R` and row #10 stays as `mp_ucell_scoring.R`. Total metaprograms/ becomes 8 files (7R + 1sh).

### Merge Group B: Enrichment Plotting (temp_plot_new.R → enrichment_plotting.R)

**Source**: `analysis/temp_plot_new.R` (693L)
**Target**: `analysis/enrichment/enrichment_plotting.R` (base: `analysis/enrich_plot.R`, 270L)

**Content to extract from temp_plot_new.R and merge into target:**
1. `enrich_heatmap()` v2 function (L490-637) — the most sophisticated version with:
   - Custom developmental RDS support
   - GeneRatio display text
   - Significance filtering
   - Viridis magma color palette
   - **This REPLACES the simpler `enrich_heatmap()` in enrich_plot.R**
2. `list_to_df(gene_list)` (L208-214) — converts named gene list to padded data frame
3. `mk_sheet(df, sheet_name)` (L368-386) — formats marker genes for Excel export

**Content to DISCARD from temp_plot_new.R:**
- PDO boxplots (L1-40) — ad-hoc, covered by compare_pdos_sc.R
- PDO MP3 heatmap (L44-69) — ad-hoc, sample-specific
- ComplexHeatmap UCell block (L71-118) — simple heatmap, not unique workflow
- enrich_heatmap() v1 (L272-345) — superseded by v2
- Jaccard PDO×SC overlap (L121-227) — covered by compare_pdos_sc.R
- clusterProfiler enrichment code (L230-480) — covered by example_anno.R/enrichment_annotation.R

**Important**: enrich_plot.R reads from external path `EAC_Ref_all/00_merged/developmental/`. This path must be preserved in the merged version.

**Note**: `enrich_heatmap()` also exists in `example_anno.R` and `Auto_MP_correlation_v2.R`. After merge, these other files may reference a local copy or eventually source from utils.R (Task 11 decides if used ≥3 times).

### Merge Group C: States QC (states_qccheck.R + qc_unresolved_states.R → states_qc.R)

**Sources**: 
- `analysis/states_qccheck.R` (357L) — produces states_status_quality_comparison.pdf + quality_comparison_by_original_state.pdf
- `analysis/qc_unresolved_states.R` (354L) — produces SAME output filenames

**Target**: `analysis/cell_states/states_qc.R`

**Merge strategy**: Task 7 must read both scripts fully. One is likely a refined version of the other (given identical output filenames). The merge should:
1. Keep the more complete/refined version as the base
2. Incorporate any unique analyses from the other version as clearly marked sections
3. Resolve the identical-output-filename collision (either parameterize output names or keep only the better version's output logic)

**Known risk**: If one script handles only "unresolved" states while the other handles all states, they may be complementary rather than duplicative. Task 7 determines this.

---

## Shell ↔ R Pair Tracking

| Shell Script | R Script | Current | Final |
|---|---|---|---|
| `analysis/Auto_MP_correlation.sh` | `analysis/Auto_MP_correlation_v2.R` | Calls `Rscript analysis/Auto_MP_correlation_v2.R` | Must update to `Rscript analysis/metaprograms/mp_database_correlation.R` |

**No other shell↔R pairs exist in top-level analysis/** (the 16 other .sh files are all in non_mali_nmf/, handled by Task 9).

---

## Design Decisions & Rationale

### 1. `legacy/` subfolder for ref_pipeline.R

**Decision**: Create `analysis/legacy/` and place `ref_pipeline_archive.R` there.

**Rationale**: 
- `ref_pipeline.R` is archival code — the original monolithic pipeline now split across root-level QC_Pipeline.R, Clustering.R, Annotation.R, Expr_filtering.R
- It contains unique functions (`write_count_matrix()`, `manual_celltyping()`, commented InSituType/coexpression code) worth preserving as reference
- Does NOT belong in `plotting/` — it's not a plotting script, even though it contains `plot_heatmap()` (the clean extraction of that function is heatmap.R → `plotting/qc_heatmap.R`)
- Does NOT belong in analysis root — it would clutter the root and confuse users about what's active vs archival
- `legacy/` clearly communicates "preserved for reference, not active analysis"

### 2. `residual.R` → `metaprograms/mp_residual_heatmap.R`

**Decision**: Place in metaprograms/, not plotting/.

**Rationale**: Despite producing heatmaps, the content is MP-centric — it computes z-score residuals of MP expression across cell states. The analysis logic is inseparable from MP interpretation. Libraries include `proxy` (distance computation) which is analytical, not just visualization.

### 3. `annotation.R` and `celltyping.R` → `cell_states/`

**Decision**: Both go to `cell_states/` as `cell_annotation.R` and `cell_typing.R`.

**Rationale**: Both scripts assign cell type labels. `annotation.R` (508L) uses Seurat + dplyr/purrr for annotation. `celltyping.R` (198L) uses readxl + marker-based assignment. They complement rather than overlap — different approaches to the same goal. Both belong with cell state analysis. They are NOT pure plotting scripts.

### 4. `summary.R` → `summary/cross_sample_summary.R`

**Decision**: Create `analysis/summary/` subfolder rather than leaving at root.

**Rationale**: 
- `summary.R` (696L) is a substantial script, not a small utility that belongs at root
- It reads per-sample data across all of `by_samples/` — it's a cross-cutting analysis
- `sum_cancer.R` could also go here, but its focus on malignant cell counts is more specific to cell states
- A `summary/` subfolder keeps the analysis root clean (only `utils.R` there)

### 5. `sum_cancer.R` → `cell_states/cancer_summary.R`

**Decision**: Place in cell_states/ rather than summary/.

**Rationale**: `sum_cancer.R` (183L) summarizes cancer/malignant cell counts — this is thematically about cell state classification results, not general cross-sample statistics. If it turns out to be more of a general summary after Task 8 reads it, it can be moved to summary/.

### 6. `robust_NMF.R` → `metaprograms/robust_nmf.R` (tentative)

**Decision**: Tentatively assign to metaprograms/. Task 4 verifies.

**Rationale**: Name suggests NMF-related analysis. The pre-existing research context says "383L, Unknown — may be NMF-related." If Task 4 finds it's unrelated to metaprograms (e.g., it's about the root NMF pipeline step), this assignment should be revised.

### 7. Keeping clinical survival scripts separate

**Decision**: `cibersort_result.R` and `evaluate_clinical_MPs.R` remain separate files.

**Rationale**: 
- `cibersort_result.R` (963L) uses CIBERSORTx deconvolution + GSVA — specific computational approach
- `evaluate_clinical_MPs.R` (400L) evaluates MPs with UCell against clinical variables — different analytical framing
- The plan says "OR keep them separate if they do fundamentally different things" — they do
- Task 6 makes the final call after full reading

---

## Accounting Verification

### By action:
- **DELETE**: 4 files (#1-4; one is non-existent)
- **MOVE**: 28 files (#5-35; content unchanged, renamed/relocated)
- **MERGE**: 3 source files → 2 target files (#36-38)
- **Non-mali NMF (Task 9)**: 18 files (#39-56)
- **NEW**: 1 file (utils.R — created by Task 3)

### Check: 4 + 28 + 3 + 18 = 53 original file entries
- Plus 2 merge targets already counted in the MOVE section (enrichment_plotting.R at #14 and states_qc.R implicit in #37-38)
- Row #4 (geneNMF.R) is non-existent — a no-op
- **55 physical files from `find` = rows 1-3 + 5-56 = 55 ✓**

### By subfolder (destination file count):
| Subfolder | R files | sh files | Total |
|---|---|---|---|
| metaprograms/ | 7 (or 8 if Merge Group A → Option 2) | 1 | 8-9 |
| enrichment/ | 5 | 0 | 5 |
| clinical/ | 3 | 0 | 3 |
| cnv/ | 3 | 0 | 3 |
| cell_states/ | 5 | 0 | 5 |
| plotting/ | 4 | 0 | 4 |
| summary/ | 1 | 0 | 1 |
| legacy/ | 1 | 0 | 1 |
| developmental/ | 1 | 0 | 1 |
| non_malignant_nmf/ | 2 | 4 (post Task 9) | 6 |
| root (utils.R) | 1 (NEW) | 0 | 1 |
| **Total** | **33-34 + 1 new** | **5** | **38-40** |

Files removed: 3 deleted + 1 non-existent + 14 per-cell-type .sh (Task 9) = 18 removed
Net: 55 original → 38-40 final files (including 1 new utils.R)

---

## Open Questions [ASK USER]

### Q1: Merge Group A — mp_ucell_scoring.R scope
`Auto_MP19_analysis.R` (internal nMP=19 UCell scoring) and `score_other_MPs.R` (external 3CA scoring) — should they merge into one `mp_ucell_scoring.R` or stay as two files (`mp_ucell_scoring.R` + `mp_external_scoring.R`)? **Recommendation**: Keep separate (Option 2). Task 4 makes final call after full reading.

### Q2: robust_NMF.R placement
`robust_NMF.R` (383L) is marked as "Unknown" in pre-existing research. Tentatively assigned to `metaprograms/robust_nmf.R`. If it turns out to be about the root-level NMF pipeline step (steps 7-8), it may belong elsewhere or in legacy/. **Task 4 must verify after reading.**

### Q3: nmf_plot.R subfolder
`nmf_plot.R` (304L) does NMF result plotting WITH enrichment (clusterProfiler + msigdbr). Placed in `enrichment/nmf_enrichment.R` because the enrichment analysis is the primary content. If Task 5 finds it's primarily NMF visualization with minimal enrichment, it could move to `metaprograms/`. **Task 5 decides after reading.**

### Q4: sum_cancer.R vs summary.R grouping
`sum_cancer.R` is placed in `cell_states/` (cancer cell focus) while `summary.R` is in `summary/` (cross-sample focus). If Task 8 finds `sum_cancer.R` is actually a general summary that just happens to focus on cancer fractions, it should join `summary/`. **Task 8 decides after reading.**
