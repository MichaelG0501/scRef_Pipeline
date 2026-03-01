# Task 1 — Triage Manifest

Generated: 2026-03-01
Files assessed: 7 (6 existing + 1 non-existent)

---

## File: analysis/temp_plot.R
**Purpose**: Ad-hoc scratch script containing three disconnected code blocks for plotting PDO metaprogram scores: (1) Violin plots of MP scores by Batch with Wilcoxon tests, (2) Boxplots of MP scores by Clinical Response and Mandard TRG score, (3) Faceted boxplots of MP3 scores by clinical response per patient. None of the code saves output to files — all uses `print()` to display interactively.
**Unique Functions**: None defined. All inline code, no function definitions.
**Inputs**:
- `mp.genes` — assumed pre-loaded in environment (from geneNMF metaprogram object)
- `pdos` — assumed pre-loaded Seurat object (PDO data)
- `batch_A` — assumed pre-loaded Seurat object (subset of PDO data)
- No file reads (`readRDS`, `read.csv`, etc.) in this script
**Outputs**:
- No files saved. All output is `print()` to console/viewer.
**Library Dependencies**: Seurat, ggpubr, patchwork, ggplot2, dplyr
**Overlaps With**: `analysis/temp_plot_new.R` (same theme — PDO MP score plotting with statistical tests, but temp_plot_new.R is a much more developed version), `analysis/compare_pdos_sc.R` (PDO analysis)
**Classification**: DELETE
**Rationale**: This is a rough scratch pad with three disconnected code blocks that never save output. All plotting patterns here (violin plots, boxplots with Wilcoxon tests, faceted patient plots) appear in more complete form in `temp_plot_new.R`. The code relies entirely on pre-loaded objects with no file I/O, making it non-reproducible. No unique functions or analyses.

---

## File: analysis/temp_plot_new.R
**Purpose**: Large multi-section exploratory script (~693 lines) covering several distinct analyses: (1) PDO MP score boxplots by Treatment×Response with Wilcoxon tests (L1-40), (2) PDO MP3 gene heatmap by Batch with top-changing genes (L44-69), (3) UCell score ComplexHeatmap across all epithelial cells sorted by study (L71-118), (4) Jaccard overlap matrix between PDO and SC metaprograms with Fisher test + pheatmap (L121-227), (5) Gene enrichment analysis using clusterProfiler + custom developmental references with `enrichGO`/`enricher` (L230-480), (6) Two versions of `enrich_heatmap()` function (v1 at L272-345, v2 at L490-637) for plotting enrichment results as pheatmaps, (7) Export of enrichment PNGs and Excel gene lists.
**Unique Functions**:
- `list_to_df(gene_list)` (L208-214): Converts named list of gene vectors to a data frame padded with NAs. Unique to this file.
- `enrich_heatmap(cluster_enrich, element, ...)` (v1: L272-345; v2: L490-637): Extracts enrichment results, builds -log10(padj) matrix, renders pheatmap with gap logic. Two versions — v2 is more sophisticated with custom RDS handling, GeneRatio display, significance filtering, and viridis magma colors. This function also exists in `analysis/example_anno.R` and `analysis/enrich_plot.R` (similar but not identical).
- `mk_sheet(df, sheet_name)` (L368-386): Formats marker gene data for Excel export. Unique to this file.
**Inputs**:
- `tmdata_all` — assumed pre-loaded (or `readRDS("EAC_Ref_epi.rds")` at L71)
- `readRDS("EAC_Ref_epi.rds")` (L71) — main epithelial Seurat object
- `readRDS("UCell_default.rds")` (L72) — MP UCell scores
- `geneNMF.metaprograms` — assumed pre-loaded (used at L75)
- `MP_pdo`, `MP_sc` — assumed pre-loaded geneNMF objects
- `readRDS("/rds/general/ephemeral/.../ref_outs/MP_outs_default.rds")` (L420)
- `read.csv("/rds/general/project/.../PDOs/Count_Matrix/New_NMFs.csv")` (L407)
- `marker_term2name`, `marker_term2gene` — assumed pre-loaded
- Custom developmental RDS files from `/rds/general/project/.../EAC_Ref_all/00_merged/developmental/per_stage/` (L424-430)
- `pdos` — assumed pre-loaded Seurat object
**Outputs**:
- `"pdo_sc_metaprograms.xlsx"` (L226) — Excel with PDO and SC MP gene lists
- `"marker_gene_list.xlsx"` (L392) — Excel with marker gene lists
- Multiple PNG files: `enrich_Hallmark.png`, `enrich_GO.png`, `enrich_MPs_3CA.png`, `enrich_Early_Embryogenesis.png`, `enrich_Normal_Development_long.png`, `enrich_Normal_Development_short.png`, `enrich_Organogenesis_major.png`, `enrich_Organogenesis_sub.png`, `enrich_Adult_Epithelium.png` (L657-693)
**Library Dependencies**: Seurat, dplyr, ggpubr, ggplot2, viridis, ComplexHeatmap, circlize, pheatmap, writexl, clusterProfiler, org.Hs.eg.db, msigdbr, enrichplot, tidyr
**Overlaps With**:
- `analysis/example_anno.R` — enrichment + heatmap plotting (has its own `enrich_heatmap`)
- `analysis/enrich_plot.R` — enrichment heatmap plotting
- `analysis/terms_overlap.R` — gene overlap enrichment → `cluster_enrich.rds`
- `analysis/compare_pdos_sc.R` — PDO vs SC Jaccard overlap
- `analysis/temp_plot.R` — simpler version of the PDO plotting
**Classification**: MERGE
**Rationale**: This is an exploratory development script that evolved into several distinct analyses. The unique content worth preserving:
- **`enrich_heatmap()` v2 function** (L490-637): The most sophisticated version with custom RDS support, GeneRatio display, significance filtering — should be merged into `analysis/enrich_plot.R` or `analysis/example_anno.R` as the canonical enrichment heatmap function.
- **`list_to_df()` and `mk_sheet()`**: Utility functions for Excel export — merge into whichever enrichment script retains the Excel export workflow.
- **Jaccard PDO×SC overlap block** (L121-227): Already covered more cleanly by `analysis/compare_pdos_sc.R`.
- **PDO boxplots** (L1-69): Ad-hoc clinical exploration, not unique enough to preserve separately.
- **ComplexHeatmap UCell block** (L71-118): Simple heatmap that could be a standalone analysis but is better merged into an MP analysis script.
- **Target**: Merge `enrich_heatmap()` v2, `list_to_df()`, `mk_sheet()` into `analysis/enrich_plot.R`. Remainder is redundant with existing scripts.

---

## File: analysis/quick.R
**Purpose**: Multi-section exploratory scratch script (~428 lines) covering four distinct analyses: (1) CNV profile heatmap visualization with ComplexHeatmap — reads CNA output, bins genomic data by 100-gene windows, creates a binned CNV heatmap with chromosome annotations and cluster boundaries (L1-252), (2) Per-sample malignancy classification using `classify_one_ident()` function with SD-based thresholds on CNA correlation and signal (L255-310), (3) Summary bar plot of malignant cell fractions per sample (L312-329), (4) Per-sample CNA scatter plots with threshold lines, faceted by study in 2×5 PDF pages (L331-428).
**Unique Functions**:
- `classify_one_ident(df_ident)` (L263-287): Classifies epithelial cells as malignant/cna_unresolved/non_malignant based on SD thresholds relative to non-epithelial reference cells. **Also exists in `InferCNA.R` (root)** — so NOT unique to this file, but the implementation here uses `celltype_update` instead of `celltype_manual`.
- `make_summary(d)` (L357-368): Computes per-sample summary of cells passing CNA thresholds. Unique to this file.
- `plot_sample(d_samp, smp)` (L371-390): Single-sample CNA scatter plot with threshold lines. Unique to this file.
**Inputs**:
- `readRDS("all_outs.rds")` (L17) — CNA output matrix
- `load("meta.RData")` (L23) — cell metadata
- `read_excel("/rds/general/project/.../Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)` (L32)
- `read.table("/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt")` (L41-44)
- `readRDS("cell_summary.rds")` (L46)
- `readRDS("all_filtered_outs.rds")` (L47)
**Outputs**:
- `saveRDS(cell_clust, "malignant_T_clust.RDS")` (L107) — hierarchical clustering result
- PDF: `"<sample>_T_CNV_profile_infercnv_style.pdf"` (L243-252)
- PNG: `"malignant_fraction_plot.png"` (L329)
- PDF: `"all.pdf"` (L397-428) — CNA scatter plots
**Library Dependencies**: data.table, dplyr, ComplexHeatmap, circlize, RColorBrewer, grid, Seurat, infercna, ggplot2, readxl, purrr
**Overlaps With**:
- `analysis/cnv_profile_sc.R` — CNV profile heatmap visualization
- `analysis/cnv_subset.R` — CNV subsetting
- `analysis/plot_CNV.R` — CNV plotting
- `InferCNA.R` (root) — `classify_one_ident()` function
- `Malignancy.R` (root) — malignancy classification
**Classification**: DELETE
**Rationale**: This is an exploratory scratch script. The CNV heatmap visualization overlaps with `analysis/cnv_profile_sc.R` and `analysis/plot_CNV.R`. The `classify_one_ident()` function already exists in the root `InferCNA.R` pipeline script (and uses a slightly different column name `celltype_update` vs `celltype_manual`). The malignant fraction bar plot and CNA scatter plots are ad-hoc visualizations that aren't part of any downstream workflow. The `make_summary()` and `plot_sample()` helper functions are trivial and tightly coupled to this exploratory analysis. No unique analysis logic that isn't covered by existing pipeline and analysis scripts.

---

## File: analysis/geneNMF.R
**Purpose**: FILE DOES NOT EXIST. No file found at `analysis/geneNMF.R`. Only the root-level `geneNMF.R` (step 8 of the pipeline) exists.
**Unique Functions**: N/A
**Inputs**: N/A
**Outputs**: N/A
**Library Dependencies**: N/A
**Overlaps With**: N/A
**Classification**: DELETE (non-existent)
**Rationale**: The file does not exist. It may have been previously deleted or was listed in error. The root-level `geneNMF.R` is the canonical step-8 NMF script. No action needed.

---

## File: analysis/expr.R
**Purpose**: Exploratory per-sample cell-typing script (~120 lines) covering two distinct approaches: (1) UMAP visualization of module scores, cluster assignments, and marker expression using ggplot2 with patchwork layout (L1-56), (2) Fisher's exact test enrichment of DE genes against cell type marker sets via `enrich_one()` function (L62-80), (3) AUCell-based cell typing — builds rankings from expression matrix, computes AUC scores per cell type marker set, assigns cell type labels to clusters by median AUC (L86-118). The script is not self-contained: it relies on pre-loaded objects (`tmdata`, `mod_cols`, `de`, `setsN`, `ct_names`, etc.) with no file I/O.
**Unique Functions**:
- `enrich_one(cl)` (L62-78): Fisher's exact test enrichment of DE genes per cluster against cell type marker gene sets. **Also exists in `Clustering.R` (root)** — so not unique.
**Inputs**:
- No file reads. All objects assumed pre-loaded: `tmdata` (Seurat), `mod_cols`, `de`, `setsN`, `universe`, `ct_names`
**Outputs**:
- No files saved. All output is interactive (`combined_plot` display).
**Library Dependencies**: patchwork, ggplot2, scales (implicit via `scales::squish`), dplyr, purrr, tidyr, AUCell (implicit via `AUCell_buildRankings`, `AUCell_calcAUC`, `getAUC`), Seurat
**Overlaps With**:
- `Clustering.R` (root) — `enrich_one()` function and AUCell-based cell typing workflow
- `analysis/celltyping.R` — cell type annotation approaches
- `analysis/beaut_umap.R` — UMAP visualization
**Classification**: DELETE
**Rationale**: This is a snippet/scratch file with no file I/O (entirely dependent on pre-loaded workspace objects). The `enrich_one()` function exists in the root `Clustering.R`. The AUCell workflow and UMAP visualization are standard Seurat operations also present in `Clustering.R`. No unique analysis or functions that aren't already in the pipeline.

---

## File: analysis/ref_pipeline.R
**Purpose**: Early monolithic draft of the entire QC + annotation pipeline (~1166 lines). Contains the full workflow in a single file: CSV loading → Seurat creation → QC inspection plots → mito filtering → CPM normalization → housekeeping/gene-count cell filtering → DoubletFinder filtering (commented out) → manual cell typing by marker expression → `plot_heatmap()` for per-sample marker expression visualization → InSituType celltyping (commented out) → co-expression filtering (commented out) → object merging and UMAP (commented out). This was the precursor to the modular pipeline split across `QC_Pipeline.R`, `Clustering.R`, `Annotation.R`, `Expr_filtering.R`.
**Unique Functions**:
- `initialise(names)` (L86-103): Reads CSV count matrices, creates Seurat objects. **Also in `QC_Pipeline.R`**.
- `inspect(tmdata_list)` (L105-192): Creates violin plots of nFeature, nCount, percent.mt per sample. **Also in `QC_Pipeline.R`**.
- `normalise(tmdata_list)` (L194-207): CPM normalization + log2 transform. **Also in `QC_Pipeline.R`**.
- `doublets_filtering(tmdata, sampleid)` (L209-254): DoubletFinder workflow per sample. **Also in `QC_Pipeline.R`** (commented out in both).
- `cells_filtering(tmdata_list, rules)` (L292-366): Filters cells by gene count and HK expression. **Also in `QC_Pipeline.R`**.
- `write_count_matrix(filtered)` (L368-384): Exports count matrices to CSV. **Unique to this file** — not in `QC_Pipeline.R`.
- `manual_celltyping(tmdata_list)` (L485-523): Assigns cell types by mean marker expression. **Unique to this file** — the main pipeline uses a different approach in `Annotation.R`.
- `plot_heatmap(tmdata, sampleid, identity, reorder)` (L597-897): Large marker expression heatmap per sample with legends and grob layout. **Also in `analysis/heatmap.R`** (nearly identical copy).
- (Commented out) `celltyping()` — InSituType-based annotation. Unique concept but commented out.
- (Commented out) `plot_coexpression()` — Co-expression heatmap. Unique concept but commented out.
**Inputs**:
- CSV count matrices from `/rds/general/project/.../EAC_Ref_all/00_counts_matrix_all/<name>.csv`
- `readLines(paste0("names_tmdata", batch, ".txt"))` (L33)
- Reads data from `read_excel(Summary_EAC_Ref.xlsx)` (commented out, L24)
- `qc_rules` data frame defined inline (L65-82)
**Outputs**:
- PDF: `paste0("Inspections", batch, ".pdf")` (L182)
- PDF: `paste0("cells_filtering", batch, ".pdf")` (L357)
- CSV: `paste0("filtering_summary", batch, ".csv")` (L446)
- RDS: `paste0("/rds/general/project/.../EAC_Ref_list", batch, ".rds")` (L527)
- PNG: `paste0(identity, "_", sampleid, ".png")` — from plot_heatmap (L892)
- (Commented out) CSV count matrices, merged RDS, clustering PNG
**Library Dependencies**: Seurat, dplyr, DoubletFinder, patchwork, ggplot2, foreach, doParallel, ComplexHeatmap, circlize, gridExtra, grid, tidyr, tibble, purrr, harmony, reticulate, readxl (implicit), data.table (via initialise), scales (implicit)
**Overlaps With**:
- `QC_Pipeline.R` (root) — `initialise()`, `inspect()`, `normalise()`, `cells_filtering()`, `doublets_filtering()`
- `Annotation.R` (root) — cell typing approach (different method)
- `analysis/heatmap.R` — `plot_heatmap()` function (nearly identical)
**Classification**: KEEP
**Rationale**: While most functions here are duplicated in the root pipeline scripts (`QC_Pipeline.R`), this file contains several unique elements worth preserving as historical reference:
1. `write_count_matrix()` — not in `QC_Pipeline.R`
2. `manual_celltyping()` — different cell typing approach from `Annotation.R`
3. `plot_heatmap()` — shared with `heatmap.R` but this is the original location
4. Commented-out `celltyping()` (InSituType) and `plot_coexpression()` — unique approaches not elsewhere
5. The `qc_rules` data frame with per-study thresholds is defined here with 10 study-specific rows (vs the wildcard in heatmap.R)
6. Full commented-out merge+UMAP workflow (L979-1166) showing the original integration approach

**Suggested filename**: `ref_pipeline_archive.R` in a `legacy/` or `archive/` subfolder. This is archival code, not active analysis.

---

## File: analysis/heatmap.R
**Purpose**: Standalone script for generating per-sample marker expression QC heatmaps using `plot_heatmap()`. Contains: (1) setup of QC parameters and marker gene lists, (2) the full `plot_heatmap()` function definition (~300 lines), (3) a driver loop that loads annotated Seurat objects from `by_samples/` and calls `plot_heatmap()` with both "strict" and "loose" co-expression filtering. The heatmap shows expression of cell-type marker genes (grouped by cell type), housekeeping expression, and gene count per cell — all as ComplexHeatmap grobs assembled into a final composite PNG.
**Unique Functions**:
- `plot_heatmap(tmdata, sampleid, identity, reorder)` (L95-395): Large function (~300 lines) that creates a composite marker expression heatmap per sample. Takes a Seurat object, generates ComplexHeatmap grobs for each marker gene group × cell type, adds housekeeping and gene count tracks, builds legends, and assembles into a final PNG. **This function is NOT in the main pipeline** (`QC_Pipeline.R`, `Clustering.R`, etc.). It also exists in `ref_pipeline.R` (nearly identical) but this file is the cleaner, standalone version with a simpler wildcard `qc_rules`.
**Inputs**:
- `readRDS(paste0("by_samples/", sample, "/", sample, "_anno.rds"))` (L401) — per-sample annotated Seurat objects
**Outputs**:
- PNG: `paste0(identity, "_", sampleid, ".png")` (L390) — large composite heatmap (80×50 inches at 400 dpi)
**Library Dependencies**: Seurat, dplyr, DoubletFinder, patchwork, ggplot2, foreach, doParallel, ComplexHeatmap, circlize, gridExtra, grid, tidyr, tibble, purrr
**Overlaps With**:
- `analysis/ref_pipeline.R` — contains nearly identical `plot_heatmap()` (this is the cleaner extraction)
**Classification**: KEEP
**Rationale**: Contains `plot_heatmap()`, a ~300-line function **confirmed NOT to exist in the main pipeline**. This is the clean, standalone version extracted from the earlier `ref_pipeline.R` monolith. Key differences from `ref_pipeline.R` version: uses `qc_rules` with wildcard `*` pattern (simpler), color scale at 0.6×max vs 0.8×max. The driver section (L398-406) shows the intended usage pattern. **Suggested filename**: `qc_heatmap.R` in a `plotting/` subfolder, since it's a QC visualization tool.

---

# Summary Table

| # | File | Lines | Classification | Target/Action |
|---|------|-------|---------------|---------------|
| 1 | analysis/temp_plot.R | 146 | DELETE | Redundant scratch; superseded by temp_plot_new.R |
| 2 | analysis/temp_plot_new.R | 693 | MERGE | Merge `enrich_heatmap()` v2 + `list_to_df()` + `mk_sheet()` → `analysis/enrich_plot.R` |
| 3 | analysis/quick.R | 428 | DELETE | Exploratory CNA scratch; covered by cnv_profile_sc.R, plot_CNV.R, InferCNA.R |
| 4 | analysis/geneNMF.R | N/A | DELETE | File does not exist |
| 5 | analysis/expr.R | 120 | DELETE | No-I/O scratch; `enrich_one()` in Clustering.R, AUCell in Clustering.R |
| 6 | analysis/ref_pipeline.R | 1166 | KEEP | Archive as `legacy/ref_pipeline_archive.R` — unique: `write_count_matrix()`, `manual_celltyping()`, commented InSituType/merge code |
| 7 | analysis/heatmap.R | 407 | KEEP | Rename to `qc_heatmap.R` — unique `plot_heatmap()` not in main pipeline |
