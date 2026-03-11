# AGENTS.md — scRef_Pipeline

Single-cell RNA-seq QC and analysis pipeline for OAC (Oesophageal Adenocarcinoma) public reference datasets. Runs on Imperial College HPC (PBS Pro scheduler). All computation is in **R** with **bash** PBS wrappers.

## Repository Structure

```
*.R              — Core pipeline scripts (steps 1-8), executed via Rscript
*_master.sh      — PBS orchestrators that loop over samples and submit per-sample jobs
N_<Step>.sh      — PBS job scripts for each pipeline step
analysis/        — Downstream analysis R scripts (ad-hoc, not part of main pipeline)
ref_outs/        — All pipeline outputs (gitignored)
temp/            — PBS stdout/stderr logs
AI_RULES.md      — HPC execution & file safety rules (READ THIS FIRST)
```

## Pipeline Execution Order

| Step | Shell Script        | R Script           | Scope       | Conda Env |
|------|--------------------|--------------------|-------------|-----------|
| 1    | `1_QC_Pipeline.sh` | `QC_Pipeline.R`    | All samples | dmtcp     |
| 2    | `2_master.sh` → `2_Clustering.sh` | `Clustering.R` | Per-sample | dmtcp |
| 3    | `3_Annotation.sh`  | `Annotation.R`     | All samples | dmtcp     |
| 4    | `4_Expr_filtering.sh` | `Expr_filtering.R` | All samples | dmtcp |
| 5    | `5_master.sh` → `5_InferCNA.sh` | `InferCNA.R` | Per-sample | dmtcp |
| 6    | `6_master.sh` → `6_Malignancy.sh` | `Malignancy.R` | Per-sample | dmtcp |
| 7    | `7_master.sh` → `7_NMF.sh` | `NMF.R`        | Per-sample  | dmtcp     |
| 8    | `8_geneNMF.sh`     | `geneNMF.R`        | All samples | gnmf      |

## Build / Run / Test Commands

There is no build system, linter, or test suite. All execution is via PBS `qsub`.

```bash
# Submit a single-step job
qsub 1_QC_Pipeline.sh

# Submit per-sample jobs (master scripts handle throttling to max 46 concurrent)
qsub 2_master.sh

# Submit a per-sample job manually (pass sample name as PBS variable)
qsub -v sample="Walker_2025_OAC26T1" -N Walker_2025_OAC26T1 2_Clustering.sh

# Run R interactively (for analysis scripts or debugging)
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/summary/cross_sample_summary.R

# For GeneNMF / UCell scripts, use the gnmf environment instead
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf
```

## HPC & File Safety Rules (from AI_RULES.md)

These rules are **mandatory** for any agent operating in this repo:

1. **Working directory**: All outputs go to `ref_outs/`. Never write outside project paths.
2. **Conda init**: Always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating envs.
3. **Interactive first**: Tasks under 8 cores / 64 GB → write only the `.R` script, no `.sh` wrapper. User runs interactively.
4. **PBS required**: Heavy tasks → must create PBS `.sh` script with `#PBS` resource headers.
5. **File naming**: New persistent files MUST be prefixed with `Auto_` (e.g., `Auto_analysis.R`).
6. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
7. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
8. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
9. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).

### PBS Job Template
```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD
Rscript <script>.R
echo $(date +%T)
```

## Code Style Guidelines

### R Scripts

**Imports**: `library()` calls at top of file, one per line. Use double quotes for package names in `library("pkg")` or no quotes `library(pkg)` — both are used; match the file you're editing. No `require()`.

**Working directory**: Each script calls `setwd()` near the top to set the output directory (typically `ref_outs/`). All file paths are then relative to that.

**CLI arguments**: Per-sample scripts receive the sample name via `commandArgs(trailingOnly = TRUE)`:
```r
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
```

**Variable naming**: `snake_case` for variables and functions. Short descriptive names. Examples:
- `tmdata`, `tmdata_list`, `tmdata_annotated` — Seurat objects
- `sample_dirs`, `names_tmdata` — sample identifiers
- `markers_list`, `markers_ranked` — marker gene collections
- `step1_calls`, `step2_calls`, `final_calls` — intermediate results

**Function naming**: `snake_case`, descriptive verbs: `initialise()`, `inspect()`, `normalise()`, `cells_filtering()`, `score_markers()`, `classify_one_ident()`, `combine_marker_scores()`.

**Data structures**: Heavy use of Seurat objects (`CreateSeuratObject`, metadata in `@meta.data`, assay layers `counts`, `data`, `CPM`). Lists of Seurat objects keyed by sample name.

**Tidyverse style**: Pipe-heavy with `%>%` and `|>` (both used). `dplyr` verbs (`mutate`, `filter`, `group_by`, `summarize`), `tidyr::pivot_longer/wider`, `purrr::map/imap`.

**Plotting**: `ggplot2` with `theme_minimal()` or `theme_classic()`. Plots saved via `ggsave()` or `pdf()`/`png()` + `dev.off()`. Composite layouts with `patchwork` (`|`, `/`) or `gridExtra::grid.arrange()`.

**Parallelism in R**: `mclapply()` from `parallel` package, or `foreach`/`doParallel`. NMF uses `nmf.options(parallel = 6)`.

**Output format**: `.rds` files for data, `.png` / `.pdf` for plots, `.csv` for summary tables. Saved to `ref_outs/` or `ref_outs/by_samples/<sample>/`.

### File Organization Pattern

Per-sample outputs follow: `ref_outs/by_samples/<Author>_<Year>_<SampleID>/`
- `<sample>.rds` — post-QC Seurat object
- `<sample>_anno.rds` — annotated
- `<sample>_epi.rds` — epithelial subset with CNA
- `<sample>_epi_f.rds` — with malignancy labels
- `<sample>_rank4_9_nrun10.RDS` — NMF results
- Sentinel files: `no_ref`, `no_epi`, `no_cell`, `no_cancer` — skip markers

### Shell Scripts

**Shebang**: `#!/bin/bash` — always first line.

**PBS directives**: `#PBS -l select=...`, `#PBS -l walltime=...`, `#PBS -N <name>`. Resource sizing varies by step.

**Timestamps**: Every script prints `echo $(date +%T)` at start and end.

**Module loading**: `module purge` then `module load tools/dev` before conda.

**Master script pattern** (for per-sample parallelism):
```bash
for sample_folder in ref_outs/by_samples/*_*_*/; do
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do sleep 180; done
  sample=$(basename "$sample_folder")
  qsub -v sample=$sample -N $sample <step>.sh
done
```

### Error Handling

- Guard clauses with `stop()` for fatal conditions (e.g., no cells, no epithelial).
- Sentinel `.rds` files written before `stop()` to mark samples to skip in subsequent steps.
- `tryCatch()` for subsetting operations that might produce empty results.
- NULL assignment + skip pattern: `if (ncol(obj) <= 1) { obj <- NULL }`.

### Key R Packages

**Core**: Seurat, dplyr, tidyr, purrr, ggplot2, readxl, Matrix, parallel
**Specialised**: infercna, NMF, GeneNMF, UCell, ComplexHeatmap, fgsea, msigdbr, reticulate
**Plotting**: patchwork, gridExtra, ggrepel, RColorBrewer, circlize, scales, colorspace

### Naming Conventions for Cell Types

Lowercase with dots: `t.cell`, `b.cell`, `nk.cell`, `macrophage`, `fibroblast`, `endothelial`, `epithelial`, `plasma`, `dendritic`, `mast`. Multi-labels joined with `|` (e.g., `t.cell|nk.cell`). Unknown/ambiguous cells labelled `unresolved`.

### Git Conventions

- `.gitignore`: ignores `ref_outs/` and PBS job output files (`*.[oe]*`)
- Commit messages: short, lowercase descriptions
- No tests, no CI/CD

## Downstream Analysis

Analysis scripts are organised into themed subfolders under `analysis/`. Outputs go to `ref_outs/`. All scripts are self-contained (zero `source()` dependencies between scripts).

### Key Shared Data Objects
All paths are relative to the project root unless absolute paths are specified.

| Object | Absolute Path | Description |
| :--- | :--- | :--- |
| Main epithelial Seurat | `ref_outs/EAC_Ref_epi.rds` | 75,348 OAC epithelial cells: merged studies |
| MP UCell scores (default) | `ref_outs/UCell_default.rds` | 75348 x 9 MP scores |
| MP UCell scores (filtered) | `ref_outs/UCell_default_filtered.rds` | Silhouette-filtered MPs |
| NMF metaprograms nMP=19 | `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` | Primary working MP object |
| Enrichment results | `ref_outs/cluster_enrich.rds` | MP x database enrichment results |
| Sample metadata (epi) | `ref_outs/meta_full_epi.rds` | Per-cell metadata with clinical variables |
| State assignments | `ref_outs/state_temp.rds` | Named vector: barcode to manual_states |
| 3CA MP UCell scores | `ref_outs/UCell_3CA_MPs.rds` | 3CA reference MP scores |
| PDO UCell scores | `ref_outs/UCell_pdos.rds` | PDO metaprogram scores: external project |
| Ref term UCell scores | `ref_outs/UCell_ref_terms.rds` | UCell scores for reference database terms |
| All sample metadata | `ref_outs/all_meta.rds` | Merged metadata for all cell types |
| NMF geneNMF outs | `ref_outs/geneNMF_outs.rds` | Raw geneNMF output object |
| NMF MP outs default | `ref_outs/MP_outs_default.rds` | Processed MP outputs |
| DEG cache | `ref_outs/states_degs.rds` | FindAllMarkers result for 4 cell states |

### `analysis/utils.R` — Reference Patterns (NOT sourced)
Documents 4 canonical patterns used across 3+ scripts each. Not sourced by any script — reference only.
- `filter_silhouette(mp_genes, sil_scores)` — remove MPs with silhouette < 0
- `score_ucell(seurat_obj, gene_lists)` — UCell scoring wrapper
- `load_epi_data(path)` — load EAC_Ref_epi.rds
- `plot_correlation_heatmap(cor_matrix)` — ComplexHeatmap correlation wrapper

### `analysis/metaprograms/` — Metaprogram Analysis
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `mp_correlation_sc.R` | SC MP UCell scores, k-means states, per-study Spearman correlation, ComplexHeatmap | `EAC_Ref_epi.rds`, `UCell_default.rds` | heatmap PDFs, `state_temp.rds` |
| `mp_correlation_pdo.R` | PDO MP UCell scores, k-means states, per-PDO Spearman correlation | PDO-specific RDS | heatmap PDFs |
| `mp_correlation_crossdata.R` | Compare PDO vs SC MPs via Jaccard similarity and UCell correlation | `UCell_pdos.rds`, `UCell_default.rds` | cross-dataset comparison PDFs |
| `mp_ucell_scoring.R` | Score nMP=19 MPs on all epi cells, silhouette filtering, Jaccard self-similarity | `geneNMF_metaprograms_nMP_19.rds`, `EAC_Ref_epi.rds` | `UCell_nMP19_filtered.rds`, heatmap PDFs |
| `mp_database_correlation.R` | UCell-score database reference terms, per-sample Spearman correlation vs MPs | `cluster_enrich.rds`, `EAC_Ref_epi.rds` | `UCell_ref_terms_v2_MP19.rds`, correlation PDFs |
| `mp_external_scoring.R` | Score cells against external 3CA MP gene sets with η² signal filter | `EAC_Ref_epi.rds`, 3CA gene lists | `UCell_3CA_MPs.rds` |
| `robust_nmf.R` | PDO NMF programs, Jaccard similarity vs 3CA MPs, GO/Hallmark enrichment | PDO NMF RDS | enrichment heatmaps |
| `mp_database_correlation.sh` | PBS wrapper for `mp_database_correlation.R` (ncpus=8, mem=128gb, walltime=4h, dmtcp env) | — | submits R job |

**Run:** `qsub analysis/metaprograms/mp_database_correlation.sh`

### `analysis/enrichment/` — Gene Overlap Enrichment
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `enrichment_analysis.R` | Compute gene overlap enrichment (TERM2GENE) across all databases | geneNMF MP gene lists | `cluster_enrich.rds` |
| `enrichment_annotation.R` | MP annotation via enrichment against Hallmarks, GO, 3CA, Pan-Cancer, developmental databases | `geneNMF_metaprograms_nMP_19.rds`, `cluster_enrich.rds` | enrichment heatmap PNGs |
| `enrichment_plotting.R` | `enrich_heatmap()` v1 and v2 functions, `list_to_df()`, `mk_sheet()` — shared plotting utilities for enrichment results | enrichment data frames | heatmap plots, Excel sheets |
| `nmf_enrichment.R` | NMF-specific enrichment, Jaccard heatmap, writexl export | NMF result RDS, gene databases | enrichment heatmaps, `.xlsx` |
| `wnt_pathway.R` | WNT pathway-specific analysis (WNT_CM, WNT_Canonical custom gene sets, Jaccard overlap) | `cluster_enrich.rds` | WNT enrichment plots |

### `analysis/clinical/` — Survival & Clinical Analysis
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `tcga_data_prep.R` | Download/process TCGA-ESCA data, merge GDC + cBioPortal clinical, Ensembl→Symbol conversion, CIBERSORTx signature prep | GDC/cBioPortal APIs | `tcga_esca_meta.rds`, TPM matrix |
| `survival_cibersort.R` | KM survival analysis (median/quartile/optimal splits), Cox PH, forest plots using DEG GSVA, MP GSVA, CIBERSORTx deconvolution | `meta_full_epi.rds`, `state_temp.rds`, `states_degs.rds` | 6 survival PDFs in `ref_outs/` |
| `survival_clinical_mps.R` | Compare Pre vs Post treatment MPs, TCGA Cox survival using UCell/3CA MP scores + cell-cycle regression | `UCell_3CA_MPs.rds`, TCGA data | survival comparison PDFs |

**Critical:** Always filter `states_degs.rds` with `group_by(cluster) %>% slice_max(n=100, order_by=avg_log2FC)` before use.

### `analysis/cnv/` — Copy Number Variation
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `cnv_profiling.R` | Run infercna to infer CNV profiles on epithelial subsets | per-sample `_epi.rds` files | CNV profiles |
| `cnv_subsetting.R` | Subset cells by CNA score, ComplexHeatmap CNA visualisation | CNV profile outputs | CNA heatmaps, filtered RDS |
| `cnv_plotting.R` | Publication-ready CNV heatmaps with circlize colour scales | CNV matrices | CNV plot PDFs |

Note: all 3 CNV scripts share identical library set: `data.table, dplyr, ComplexHeatmap, circlize, RColorBrewer, Seurat, infercna`.

### `analysis/cell_states/` — Cell Type & State Assignment
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `cell_annotation.R` | Automated cell type annotation using marker gene scoring | per-sample annotated RDS | updated `celltype` metadata |
| `cell_typing.R` | Manual cell type assignment via `classify_one_ident()` and `combine_marker_scores()` | annotated RDS, marker lists | updated cell type calls |
| `state_assignment.R` | MP z-score residual computation (regress out cell-cycle MPs), cosine-similarity state assignment (SC: 5 MPs → 4 states; PDO: 4 MPs → 6 states), ComplexHeatmap visualisations | `EAC_Ref_epi.rds`, MP UCell scores | state assignment vectors, heatmaps |
| `states_qc.R` | QC comparison of Defined vs Unresolved cells across 8 continuous + categorical features (merged from 2 original scripts) | `EAC_Ref_epi.rds`, state assignments | `states_status_quality_comparison.pdf` |
| `states_umap.R` | UMAP embedding coloured by manual states | `EAC_Ref_epi.rds`, state assignments | UMAP PDFs |
| `cancer_summary.R` | Summary statistics of cancer/malignant cell counts per sample | per-sample `_epi_f.rds` files | summary CSV/plots |
| `Auto_states_cluster.R` | Louvain clustering on CC-regressed, Z-normed MP scores (PCA → FindNeighbors → FindClusters at res 0.5/0.8/1.0), UMAP, 5 visualisations | `EAC_Ref_epi.rds`, `geneNMF_metaprograms_nMP_19.rds`, `UCell_nMP19_filtered.rds` | `Auto_cluster_states.rds`, `Auto_cluster_umap_embeddings.rds`, `Auto_cluster_mp_adj.rds`, cluster PDFs |
| `Auto_states_topmp.R` | Assign cells to dominant non-CC MP (max Z-score, threshold 0.5 → Unresolved), 4 visualisations | `EAC_Ref_epi.rds`, `geneNMF_metaprograms_nMP_19.rds`, `UCell_nMP19_filtered.rds` | `Auto_topmp_states.rds`, `Auto_topmp_mp_adj.rds`, topmp PDFs |
| `Auto_states_comparison.R` | Compare cluster vs top-MP states: confusion matrix, ARI/NMI, side-by-side UMAP, bootstrap stability (100×80%), Cramér's V study-bias, DEG + fgsea Hallmark coherence | `Auto_cluster_states.rds`, `Auto_topmp_states.rds`, `Auto_cluster_umap_embeddings.rds`, `Auto_cluster_mp_adj.rds`, `Auto_topmp_mp_adj.rds` | `Auto_comparison_summary.csv`, comparison PDFs |
| `Auto_states_cluster.sh` | PBS wrapper for `Auto_states_cluster.R` (ncpus=8, mem=96gb, walltime=4h, dmtcp env) | — | submits R job |
| `Auto_states_topmp.sh` | PBS wrapper for `Auto_states_topmp.R` (ncpus=4, mem=64gb, walltime=2h, dmtcp env) | — | submits R job |
| `Auto_states_comparison.sh` | PBS wrapper for `Auto_states_comparison.R` (ncpus=8, mem=96gb, walltime=4h, dmtcp env) | — | submits R job |

### `analysis/plotting/` — Publication Figures
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `publication_umap.R` | Publication-ready UMAP figures (DimPlot, FeaturePlot) with polished themes | `EAC_Ref_epi.rds`, state metadata | UMAP PDFs |
| `gene_expression_heatmap.R` | Gene expression heatmap using ComplexHeatmap (averaged by state/cluster) | `EAC_Ref_epi.rds`, gene lists | expression heatmap PDFs |
| `clinical_variable_plots.R` | Boxplots and violin plots of clinical variables (treatment, stage) by cell state | `meta_full_epi.rds`, state assignments | clinical variable PDFs |
| `qc_heatmap.R` | QC metric heatmap via `plot_heatmap()` — unique function not in main pipeline | per-sample RDS, QC metadata | QC heatmap PDFs |

### `analysis/summary/` — Cross-Sample Summary
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `cross_sample_summary.R` | Summary statistics (cell counts, QC metrics, composition) across all samples | `ref_outs/by_samples/` per-sample RDS files | summary tables and plots |

**Run interactively:** `Rscript analysis/summary/cross_sample_summary.R`

### `analysis/developmental/` — Developmental Gene Reference
| File | Purpose | Key Inputs | Key Outputs |
| :--- | :--- | :--- | :--- |
| `developmental.R` | Build developmental gene reference from 4 external xlsx files (Early Embryogenesis, Organogenesis, Normal Development, Adult Epithelium) | xlsx files in spatialtranscriptomics project | `enrich_dev.rds`, per-stage RDS files |

**Output path:** `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/enrich_dev.rds`
This path is also used by `enrichment/enrichment_plotting.R`.

### `analysis/non_malignant_nmf/` — Non-Malignant Cell NMF
Runs GeneNMF on non-malignant cell types (macrophage, fibroblast, endothelial, nk.cell, plasma, cd4, cd8).

| File | Purpose |
| :--- | :--- |
| `Auto_geneNMF_celltype.R` | GeneNMF analysis per cell type (receives `celltype` arg) |
| `Auto_anno_celltype.R` | Annotate NMF programs per cell type (receives `celltype` arg) |
| `run_geneNMF.sh` | Parameterised PBS job runner for geneNMF (ncpus=8, mem=72gb, walltime=8h, gnmf env) |
| `run_annotation.sh` | Parameterised PBS job runner for annotation (ncpus=2, mem=16gb, walltime=3h, dmtcp env) |
| `submit_geneNMF_all.sh` | Master: submits geneNMF jobs for all 7 cell types (throttled ≤46 concurrent) |
| `submit_annotation_all.sh` | Master: submits annotation jobs for all 7 cell types (checks RDS existence first) |

**Run:** `qsub analysis/non_malignant_nmf/submit_geneNMF_all.sh`  
Or per cell type: `qsub -v celltype=macrophage analysis/non_malignant_nmf/run_geneNMF.sh`

**Note:** geneNMF R uses `nk.cell` as celltype key; annotation R uses `nk` as ct_map key. Submit scripts pass the correct key for each.

### `analysis/legacy/` — Archived Code
| File | Purpose |
| :--- | :--- |
| `ref_pipeline_archive.R` | Archive of original monolithic pipeline draft (1166L). Contains unique functions: `write_count_matrix()`, `manual_celltyping()`, `plot_coexpression()`. Not intended for execution. |

### Critical Recurring Patterns

**MP Silhouette Filtering**
Filter out MPs with a silhouette score below 0 before any downstream analysis. This is a strict requirement.
```r
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
```

**Sample Identity**
Use `orig.ident` from the Seurat metadata to group by sample. Avoid barcode manipulation like `sub("_[^_]+$", "", ...)`.

**Environment Usage**
- `gnmf` conda env: Use for UCell scoring and GeneNMF scripts.
- `dmtcp` conda env: Use for general Seurat and analysis tasks.

**Metaprogram Resolution**
The pipeline explores nMP range 8 to 30. **nMP=19** is the current selected working resolution.

## NotebookLM Skill (HPC Prerequisites)

The NotebookLM skill uses a Chromium browser via Patchright. On the Imperial HPC, the system `libgbm` library is missing — **you must load the Mesa module before any NotebookLM operation** (queries, auth checks, etc.):

```bash
module load Mesa/23.1.4-GCCcore-12.3.0
```

This applies to both interactive use and PBS jobs that invoke NotebookLM.

- **Skill location**: `/rds/general/user/sg3723/home/.claude/skills/notebooklm`
- **Auth state**: `~/.claude/skills/notebooklm/data/browser_state/state.json`
- **Auth setup requires X11**: Use MobaXterm (has built-in X11 forwarding) to run `python scripts/run.py auth_manager.py setup` — a visible browser window is needed for Google login.
- **Re-auth**: If auth expires, re-run setup via MobaXterm with the Mesa module loaded.

## AGENTS.md Living Document Rule

Future agents must update this file when they:
- Find a new analysis script and define its purpose.
- Locate new input or output file paths.
- Spot a recurring pattern or technical hurdle.
- Create a new `Auto_` script.

Append new findings to the appropriate section. Don't rewrite existing documentation unless fixing an error.

## Subagent Model Tier Policy (MANDATORY)

When delegating work to subagents or background tasks, you **MUST** restrict model choices based on the tier of the primary model you are currently running as. Check your own model identity and apply the corresponding rule:

### Tier Classification

**Free Tier** (zero cost):
- `opencode/big-pickle`
- `opencode/minimax-m2.5-free`
- `opencode/trinity-large-preview-free`
- `github-copilot/gpt-4.1`
- `github-copilot/gpt-4o`
- `github-copilot/gpt-5-mini`

**0.33X Tier** (reduced cost):
- `github-copilot/gemini-3-flash-preview`
- `github-copilot/claude-haiku-4.5`
- `github-copilot/gpt-5.1-codex-mini`
- `github-copilot/grok-code-fast-1`

**All Other Models** = Paid Tier (full cost)

### Delegation Rules

| Your Primary Model Tier | Allowed Subagent Models |
|------------------------|------------------------|
| **Free** | Free tier models ONLY |
| **0.33X** | Free + 0.33X tier models |
| **Paid** (any other) | Any model — no restrictions |

### How to Apply

1. **Identify your tier**: Check which model you are running as. If it appears in the Free list → you are Free tier. If in the 0.33X list → you are 0.33X tier. Otherwise → Paid tier.
2. **Constrain task() calls**: When using the `task()` tool or delegating to subagents, prefer models from your allowed tier. If a category default model exceeds your tier, override it with a model from your allowed list.
3. **Recommended free-tier subagent assignments**:
   - Quick/trivial tasks → `github-copilot/gpt-5-mini`
   - Code implementation → `github-copilot/gpt-4.1`
   - Reasoning/analysis → `github-copilot/gpt-4o`
   - Research/exploration → `github-copilot/gpt-4.1`
   - Fallback/general → `opencode/big-pickle`
4. **Do not escalate**: Never use a paid model as a subagent when your primary is free or 0.33X tier. This is a hard rule.

## Uploading Files to Google Drive (rclone)

rclone is configured with a remote named `gdrive`. Always upload into the `IMPERIAL/` folder:

```bash
module load rclone
rclone copy <local_file> gdrive:IMPERIAL/ --progress
```

To upload all `.pdf`, `.png`, and `.csv` files generated within the last 8 hours from `ref_outs/`:

```bash
module load rclone
find ref_outs -name '*.pdf' -o -name '*.png' -o -name '*.csv' | xargs -I{} find {} -mmin -480 2>/dev/null | while read f; do rclone copy "$f" gdrive:IMPERIAL/ --progress; done
```

Or more simply with `find`'s `-mmin` flag:

```bash
module load rclone
find ref_outs \( -name '*.pdf' -o -name '*.png' -o -name '*.csv' \) -mmin -480 -exec rclone copy {} gdrive:IMPERIAL/ --progress \;
```

- Remote name: `gdrive`
- Target folder: `gdrive:IMPERIAL/` (always use this, not the root)
- Must `module load rclone` before any rclone command on HPC

## Progress Updates (Bi-Weekly)

Updates generated every **Monday and Thursday**. Workflow: `updates/weekly_update_workflow.md`.

**Discovery:** `git status --short` — untracked files = new work since last commit.
Associated outputs identified from script header comments (`# Output:` lines).

**Slide style:** Result-focused. Large figures, 1–2 sentence key points per slide.
The `.md` companion doc retains full technical detail.

**Directory structure:** Each update cycle gets its own subfolder:
`updates/DDmon/` containing `figures/`, `summaries/`, `.tex`, `.md`, and `_plan.md`.

**Machine-readable summaries (convention for all future `Auto_` scripts):**
Every script that produces plots must also save a small (< 100 KB) `.csv` or `.txt`
summary of key metrics directly into `updates/new_updates/summaries/` so AI agents
can read results on the login node without loading heavy `.rds` files, create folder if not exist.
