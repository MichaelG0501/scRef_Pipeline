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
```

## HPC & File Safety Rules (from AI_RULES.md)

These rules are **mandatory** for any agent operating in this repo:

1. **Working directory**: All outputs go to `ref_outs/`. Never write outside project paths.
2. **Conda init**: Always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating envs.
3. **Interactive first**: Tasks under 8 cores / 64 GB → write only the `.R` script, no `.sh` wrapper. User runs interactively.
4. **PBS required**: Heavy tasks → must create PBS `.sh` script with `#PBS` resource headers.
5. **Live Logging**: Always use live streaming log file mode by adding `#PBS -koed` to the submission script. This ensures standard out and standard error are written to their final destination as the job is running, allowing for real-time monitoring from login nodes.
6. **Output file naming**: New persistent generated outputs should keep the historical `Auto_` prefix when useful for provenance. New analysis script filenames should instead be descriptive and follow the `analysis/ANALYSIS_MAP.md` naming rules.
7. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
8. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
9. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
10. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).

### PBS Job Template
```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
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
**Specialised**: infercna, NMF, GeneNMF, UCell, ComplexHeatmap, fgsea, msigdbr, reticulate, monocle3, SeuratWrappers
**Plotting**: patchwork, gridExtra, ggrepel, RColorBrewer, circlize, scales, colorspace

### Naming Conventions for Cell Types

Lowercase with dots: `t.cell`, `b.cell`, `nk.cell`, `macrophage`, `fibroblast`, `endothelial`, `epithelial`, `plasma`, `dendritic`, `mast`. Multi-labels joined with `|` (e.g., `t.cell|nk.cell`). Unknown/ambiguous cells labelled `unresolved`.

### Git Conventions

- `.gitignore`: ignores `ref_outs/` and PBS job output files (`*.[oe]*`)
- Commit messages: short, lowercase descriptions
- No tests, no CI/CD

## Downstream Analysis

Analysis scripts are organised into themed subfolders under `analysis/`. The authoritative script inventory, dependency map, run order, terminal-figure list, legacy/delete-candidate list, and untracked-file notes live in `analysis/ANALYSIS_MAP.md`. Future agents must update that map whenever scripts are added, renamed, moved, superseded, or given new downstream outputs.

### Required Analysis Script Standards

Every active analysis script must start with a registry/header block that states:
- script status: `active`, `terminal`, `legacy`, or `delete-candidate`
- short description
- methodology file path under `analysis/methodology/`
- exact inputs, including absolute external files or download requirements
- exact outputs, grouped as `intermediate/`, `tables/`, `figures/`, `logs/`, or `reports/` for new long-running scripts
- cache/replot behavior for expensive workflows
- run command and conda environment

New script filenames should be descriptive and should not use `Auto_` as the script prefix. Use `legacy_` for retained comparison/historical scripts and `delete_` for scripts recommended for manual deletion. Do not delete scripts directly unless the user explicitly asks. Existing output files may keep their historical `Auto_` prefixes when downstream scripts already depend on them.

### Shared Configuration And Helpers

Shared constants live in `analysis/shared/scRef_config.R`. Shared helper functions live in `analysis/shared/scRef_helpers.R`. New scripts should reuse or copy from these files for:
- preferred state definition: Approach B, noreg
- common input paths and output directories
- state/MP order and colors
- plot dimensions and presentation-readable typography defaults
- metadata column names
- common thresholds
- output-tier helpers
- lightweight run-summary logging

Scripts may source these helper files when the dependency is documented in the script header. For fully self-contained HPC scripts, copying a small stable helper block is acceptable, but avoid copy-paste drift for shared constants.

### Current Key Data Objects

| Object | Path | Notes |
| :--- | :--- | :--- |
| Main epithelial Seurat | `ref_outs/EAC_Ref_epi.rds` | merged OAC epithelial reference |
| MP object nMP=19 | `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` | primary MP gene set object |
| Filtered MP UCell | `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds` | silhouette-filtered MP scores |
| 3CA MP UCell | `ref_outs/UCell_3CA_MPs.rds` | pan-cancer MP scores |
| Epithelial metadata | `ref_outs/meta_full_epi.rds` | per-cell clinical/sample metadata |
| Approach B noreg states | `ref_outs/Auto_topmp_v2_noreg_states_B.rds` | preferred upstream state definition |
| Approach B noreg MP matrix | `ref_outs/Auto_topmp_v2_noreg_mp_adj.rds` | z-normalised state-definition MP matrix |
| Final state vector | `ref_outs/Auto_final_states.rds` | preferred downstream state object |
| Enrichment results | `ref_outs/cluster_enrich.rds` | MP x database enrichment object |
| Final MP SCENIC directory | `ref_outs/final_mp_scenic/` | selected cells, regulons, networks |

Do not use `ref_outs/state_temp.rds`, `ref_outs/Auto_topmp_states.rds`, `ref_outs/Auto_cluster_states.rds`, or `ref_outs/Auto_topmp_v2_states_B.rds` for new downstream work.

### Current Run Order

1. Score MPs: `analysis/metaprograms/mp_ucell_scoring.R` and `analysis/metaprograms/mp_3ca_ucell_scoring.R`
2. Define states: `analysis/cell_states/state_definition_approach_b_reg_noreg.R`
3. Create final states: `analysis/cell_states/final_state_unresolved_relabel.R`
4. Generate state/MP figures: `state_mp_sample_abundance.R`, `final_state_overall_proportions.R`, `top_diverse_sample_state_umap.R`, `basal_smg_mp_signature_heatmap.R`, `final_state_marker_discovery.R`, `final_mp_scenic.R`
5. Run pseudotime/hybrid workflows as needed: `pseudotime_top_diverse_samples.R`, `pseudotime_state_distance_matrix.R`, `hybrid_pairwise_percell_heatmap.R`, `hybrid_pairwise_distance_nodeplot.R`, `pseudotime_trajectory_report_partA.R`
6. Run clinical/bulk workflows: `tcga_mp_state_survival_reg_noreg.R`, `geo_survival_data_prep.R`, `geo_survival_mp_state_survival.R`, `bulk_tcga_geo_qc.R`, `bulk_tcga_geo_integrated_survival.R`, `tcga_mp_state_survival_qc_filtered.R`, `clinical_association_final_boxplots.R`
7. Run CNV/subclone workflows: `cnv_profiling.R`, `cnv_subsetting.R`, `cnv_plotting.R`, `cnv_malignant_subclone_mp_heatmap.R`, `cna_subclone_expression_correlation.R`
8. Optional external/reference workflows: developmental, enrichment, non-malignant NMF, and spatial scripts documented in `analysis/ANALYSIS_MAP.md`

### Output Tiers, Caches, And Logs

New long-running scripts should write outputs into clear tiers beneath their analysis output directory:
- `intermediate/`: cacheable RDS/model/matrix objects used for replotting
- `tables/`: final CSV/TSV/XLSX outputs
- `figures/`: final PDF/PNG/SVG figures
- `logs/`: run summaries and session information
- `reports/`: multi-page narrative PDFs or report bundles

Plot-only changes should be reproducible from `intermediate/` without rerunning heavy computation. Support `SCREF_FORCE_REBUILD=TRUE` to ignore caches and `SCREF_REPLOT_ONLY=TRUE` where practical. Long-running scripts should save a lightweight run summary with start/end time, input files, output files, parameters, cache reuse status, and session/package versions when relevant.

### Plot Readability

All final figures must be readable on PowerPoint slides. Increase figure size or split pages rather than shrinking labels to unreadable sizes. Explicitly tune font size, legend text, legend title, row/column names, point size, line width, and heatmap legend sizes relative to the output dimensions.

### Methodology Documentation

Methodology files live under `analysis/methodology/` with subfolders matching the corresponding `analysis/` script folder. Folder-level methodology is acceptable for related scripts; complex workflows should also get script-specific methodology files. The methodology should document the implemented logic in operational detail, including external downloads or absolute reference files.

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
- `dmtcp` conda env: Use for general Seurat and analysis tasks.
- `gnmf` conda env: Use for GeneNMF package.

**Adult oesophagus external reference**
- `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/` is a very large Matrix Market dataset (`EoE.mtx` plus metadata), so interactive MP scoring should subset epithelial barcodes before UCell scoring and cache the sampled subset under `ref_outs/Auto_external_epi_mp_ucell/`.

**Metaprogram Resolution**
The pipeline explores nMP range 8 to 30. **nMP=19** is the current selected working resolution.

**reg vs noreg Sensitivity Analysis**
Evaluate the impact of cell-cycle (CC) regression on state assignments and downstream associations (survival, clinical). 
- **reg**: Z-score MP scores *after* regressing out CC MPs (CC_G1S, CC_G2M).
- **noreg**: Z-score MP scores directly.
Required for establishing the robustness of state-linked clinical findings.

**Final decision: noreg + Approach B.** All new scripts (state_mp_sample_abundance, pseudotime_top_diverse_samples, final_state_unresolved_relabel, hybrid_pairwise_percell_heatmap, cibersortx_sc_reference_export) use **noreg Approach B only** — no reg/noreg parameterisation.

**Final MP SCENIC panel.** The curated final-MP regulatory workflow uses 14 scATLAS MPs (`MP1/7/9/2/17/14/5/10/8/13/12/18/16/15`) plus the 3 retained 3CA MPs from unresolved relabeling (coverage thresholds: `n_samples >= 50`, `n_studies >= 6`, `pct_cells >= 1`). All exported plots should label columns/states as full `MP# + description` strings rather than short aliases.

**SCENIC databases.** `final_mp_scenic.R` expects human cisTarget `.feather` databases via `SCENIC_DB_DIR` (or the default `ref_outs/final_mp_scenic/cistarget_databases`). The script supports `prepare_only=true` to validate cell selection and final-MP gene sets when SCENIC packages or databases are not yet installed.

**SCENIC compatibility note (25 Mar 2026).** The files in `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases/` with names like `hg38_*full_tx_v10_clust*.feather` are not compatible with the R `SCENIC`/`RcisTarget` workflow here because they lack the required `features` index column. A verified RcisTarget-compatible replacement directory is `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`, containing:
- `hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather`
- `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather`
On this cluster, `SCENIC` also needs two small runtime compatibility patches:
- fallback from `motifAnnotations_hgnc` to `motifAnnotations_hgnc_v9`
- sparse-aware `geneFiltering()` because the package-level base `rowSums()` call fails on `dgCMatrix`
`final_mp_scenic.R` now patches both automatically before calling the SCENIC workflow.

**Final-state downstream plotting**
- Prefer `Auto_final_states.rds` for downstream clinical/state association plots when it exists. This object includes the unresolved relabel updates and the current 8-state naming (`Immune Infiltrating`, `Basal to Intestinal Metaplasia`, `3CA_EMT_and_Protein_maturation`).

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
- Create a new analysis script.

Append new findings to the appropriate section. Don't rewrite existing documentation unless fixing an error.

## 18 Mar 2026 Naming + RDS Compatibility Update

- **State naming update applied**: use `Basal to Intestinal Metaplasia` (previously `Basal to Intest. Meta` / `Barretts Metaplasia`) and `Immune Infiltrating` (previously `Immune Infiltrated`) across new analysis outputs.
- **MP/state definitions unchanged**: MP memberships and state compositions remain identical; only labels changed.
- **Compatibility risk**: existing `.rds` outputs generated before this rename may still carry old labels; mixed usage can break joins, plotting order, pseudotime roots, and survival feature mapping.
- **Progressive relabeling update (23 Mar 2026)**:
  - `3CA_mp_30 Respiration 1` relabeled cells merged into `Classic Proliferative`.
  - `3CA_mp_12 Protein maturation` and `3CA_mp_17 EMT III` merged into `3CA_EMT_and_Protein_maturation`.
  - Final results stored in `ref_outs/Auto_final_states.rds`. This is the recommended object for all downstream clinical/abundance analysis.
- **Regeneration guidance after rename**:
  1. Regenerate `Auto_topmp_v2_{reg,noreg}_states_B.rds` and dependent outputs.
  2. Regenerate unresolved relabel outputs (`unresolved_states/*`) so newly relabeled states are consistent.
  3. Regenerate pseudotime outputs (`ref_outs/pseudotime/*`) with updated root definitions.
  4. Regenerate abundance and hybrid outputs so state/MP ordering and labels stay synchronized.

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

## Uploading Files to Google Drive (rclone) (ONLY do this when required explicitly)

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

**Machine-readable summaries (convention for all future analysis scripts):**
Every script that produces plots must also save a small (< 100 KB) `.csv` or `.txt`
summary of key metrics directly into `updates/new_updates/summaries/` so AI agents
can read results on the login node without loading heavy `.rds` files, create folder if not exist.

## 25 Mar 2026 Pseudotime State-Distance Scripts

- `analysis/cell_states/pseudotime_state_distance_matrix.R`
  - Purpose: Part A only, all-sample Monocle3 state-distance workflow for the 5 primary `noreg` Approach B states.
  - Inputs: `ref_outs/EAC_Ref_epi.rds`, `ref_outs/Auto_topmp_v2_noreg_states_B.rds`
  - Sample inclusion defaults: root state `>= 20` cells, at least one second state `>= 20` cells, total primary-state cells `>= 80`, and at least 2 represented states over threshold.
  - Methods compared:
    - directed pseudotime distance after iterating each state as Monocle root
    - principal-graph geodesic distance using state medoid-cell projections
    - principal-graph geodesic distance using state centroid projections
    - UMAP centroid Euclidean baseline
  - Additional plotting output: also orders each valid sample with `Basal to Intestinal Metaplasia` as the root and saves a state-only trajectory page into `ref_outs/state_distance_pseudotime/Auto_state_pseudotime_all_valid_samples.pdf`, plus per-sample pseudotime vectors in `ref_outs/state_distance_pseudotime/sample_state_trajectories/`.
  - Outputs: `ref_outs/state_distance_pseudotime/` containing per-sample CSVs, final matrix CSV/RDS files, the all-sample basal-root trajectory PDF, per-sample pseudotime `.rds` files, and a comparison heatmap PDF.

- `analysis/cell_states/hybrid_pairwise_distance_nodeplot.R`
  - Purpose: hybrid pairwise node plot with state positions derived from a chosen biological distance matrix instead of fixed circular spacing.
  - Inputs: `ref_outs/Auto_topmp_v2_noreg_states_B.rds`, `ref_outs/Auto_topmp_v2_noreg_group_max.rds`, `ref_outs/state_distance_pseudotime/Auto_state_distance_matrices.rds`
  - Layout: compares classical MDS (`cmdscale`) and non-metric MDS (`isoMDS`), then keeps the 2D fit with the fewest nearest-neighbor mismatches and lowest residual distance error after a single global scale factor.
  - Outputs: `ref_outs/task6_hybrid_pairwise_distance/` with per-method nodeplot/heatmap PDFs, edge CSVs, layout CSVs, fit-diagnostic CSVs (`Auto_task6_hybrid_pairwise_fit_pairs_<method>.csv`, `Auto_task6_hybrid_pairwise_fit_nearest_<method>.csv`), plus combined all-method PDFs: `Auto_task6_hybrid_pairwise_nodeplot_all_methods.pdf` and `Auto_task6_hybrid_pairwise_distance_heatmap_all_methods.pdf`

- `analysis/cell_states/pseudotime_state_distance_matrix.sh`
  - PBS wrapper for `pseudotime_state_distance_matrix.R` (`ncpus=8`, `mem=128gb`, `walltime=12h`, `dmtcp` env)

- `analysis/cell_states/hybrid_pairwise_distance_nodeplot.sh`
  - PBS wrapper for `hybrid_pairwise_distance_nodeplot.R` (`ncpus=4`, `mem=32gb`, `walltime=2h`, `dmtcp` env)

- `analysis/cell_states/Auto_submit_state_distance_and_nodeplot.sh`
  - Convenience submitter: submits the distance job first, then submits the nodeplot job with `afterok` dependency and optional `distance_method` argument.

## PBS Command Locations

- On this HPC, PBS binaries are available under `/opt/pbs/bin/`.
- Common commands:
  - `qsub`: `/opt/pbs/bin/qsub`
  - `qstat`: `/opt/pbs/bin/qstat`
  - `qdel`: `/opt/pbs/bin/qdel`
  - `pbsnodes`: `/opt/pbs/bin/pbsnodes`
- In non-interactive shells, these commands may not be on `PATH` by default, so use the absolute path if needed.

####################
## 22 Jun 2026 Refined MP Split Correlation Ordered Heatmap

- `analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R`
  - Purpose: terminal correlation heatmap aligned to `refined_mp_nmf_ordered_heatmap.R`; finalized MP blocks keep the same program-resolution height/width as the ordered NMF heatmap.
  - Critical plotting detail: rows/columns are expanded to original NMF programs and filled by split/full MP UCell correlations, so merged final MPs such as `MP2+` show all constituent pre-merge sub-MPs internally while full/single MPs such as `MP9` and `MP7j` remain single internal blocks.
  - Styling: final MP dotted boxes and label spacing match the ordered NMF heatmap; internal split-MP boxes use thin grey borders; extra label gaps are before `MP17`, `MP7r`, `MP8c`, `MP8b`, `MP12b`, `MP15a`, and `MP15b`.
  - Outputs: `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.{pdf,png}`, final/sub-block CSVs, `intermediate/refined_mp_split_correlation_ordered_matrix.rds`, and `updates/new_updates/summaries/refined_mp_split_correlation_ordered_heatmap_summary.csv`.
  - Methodology: `analysis/methodology/metaprograms/refined_mp_split_correlation_ordered_heatmap_methodology.md`.
####################

####################
## 21 Jun 2026 Refined MP NMF Ordered Heatmap

- `analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R`
  - Purpose: terminal plot-only heatmap that reorders the original GeneNMF programme similarity matrix by finalized merged refined MPs after `mp_refinement_submp.R` and `mp_refinement_merge_correlated_submps.R`.
  - Strict order: `MP7j, MP9, MP1, MP2+, MP17, MP8+, MP10+, MP14, MP5+, MP7r, MP7v, MP10e, MP16+, MP18, MP8c, MP15c, MP12c, MP2v, MP8e, MP12a, MP13, MP7+, MP7h, MP8b, MP12b, MP15a, MP15b`.
  - Critical plotting detail: dotted diagonal borders are generated from contiguous `merged_refined_mp` runs in `merged_refined_mp_assignments.rds`, so the borders reflect finalized refined MPs rather than original nMP19 clusters.
  - Outputs: `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.{pdf,png}`, `tables/refined_mp_nmf_ordered_blocks.csv`, `intermediate/refined_mp_nmf_ordered_similarity.rds`, and `updates/new_updates/summaries/refined_mp_nmf_ordered_heatmap_summary.csv`.
  - Methodology: `analysis/methodology/metaprograms/refined_mp_nmf_ordered_heatmap_methodology.md`.
####################

## 3 Jun 2026 CNA Subclone Functional Exclusivity Update

- **Exclusivity Logic**: Updated `analysis/cnv/cnv_malignant_subclone_mp_heatmap.R` with refined "functional exclusivity" definitions.
    - **Cell States**: Focus on 5 main biological states. A state is "exclusive" if it is present in exactly one subclone.
    - **Metaprogrammes (MPs)**: An MP is "exclusive" if its mean UCell score is > 0.10 in exactly one subclone.
- **Visualization**: Redesigned Page 1 of the cohort summary PDF (`Auto_malignant_subclone_mp_cohort_summary.pdf`) with neutral color palettes (gray/black) to avoid collisions with state-label colors. Added a "State exclusivity" percentage bar plot alongside subclone counts and MP exclusivity.
- **Maintenance**: Ported temporary regeneration logic into the canonical pipeline script to ensure consistency in future runs.
- **Resource Discovery**: Confirmed that `BiocParallel` workers must be restricted to 2 even if more cores are allocated, due to cluster-level environment limits (`_R_CHECK_LIMIT_CORES_`).

####################
## 5 Jun 2026 scATLAS Numbat Raw Expression Concordance

- `analysis/cnv/scatlas_numbat_raw_expression_concordance.R`
  - Purpose: compare Numbat raw expression-roll CNA values from `gexp_roll_wide.tsv.gz` against unfiltered per-sample InferCNA `_outs.rds` matrices for matched cells.
  - Important distinction: visible structure in `exp_roll_clust.png` is raw expression-CNA structure, not necessarily an accepted Numbat final subclone. Samples with Numbat terminal statuses such as `No clones remain after filtering by size` can still have two raw-expression clusters in this validation plot.
  - Outputs: `ref_outs/Auto_scatlas_numbat/raw_expression_concordance/` containing combined and per-sample Numbat-vs-InferCNA heatmaps, raw-expression cluster labels, and arm-level concordance summaries.
  - PBS wrapper: `analysis/cnv/scatlas_numbat_raw_expression_concordance.sh` (`dmtcp`, 4 cores, 192 GB, 12 h, live logging).
- Methodology: `analysis/methodology/cnv/scatlas_numbat_raw_expression_concordance_methodology.md`.
####################

####################
## 2 Jun 2026 scATLAS Numbat Terminal No-Subclone Handling

- `analysis/cnv/Auto_scatlas_numbat_run_sample.R` treats Numbat terminal biological statuses such as `No clones remain after filtering by size` and `No CNV remains after filtering by LLR in pseudobulks` as valid no-subclone outcomes for this validation workflow. It writes a done file, empty clone/joint summaries, and an RDS summary with `terminal_no_subclone = TRUE` instead of failing the PBS dependency chain.
- `analysis/cnv/Auto_scatlas_numbat_conservative_recut.R` records these samples as `terminal_no_subclone` in `ref_outs/Auto_scatlas_numbat/conservative_clones/Auto_scatlas_numbat_conservative_clone_summary.csv` and does not attempt to re-cut a missing Numbat tree.
- Do not loosen Numbat thresholds solely to force clone calls for these samples; the requested validation layer is conservative and should avoid over-fragmenting weak/no-CNV samples.
####################
## 29 May 2026 scATLAS Raw Redownload And Numbat Validation

- Raw-data rebuild scripts now live in `analysis/raw_data/`:
  - `Auto_download_alcindor_srr.sh`: downloads Alcindor `SRR27335925`-`SRR27335944` with `fasterq-dump --split-files --include-technical`.
  - `Auto_download_carroll_ega.sh`: downloads Carroll EGA dataset `EGAD00001009401` with `pyega3`; use `EGA_CREDENTIAL_JSON` for credentials and do not copy credentials into scripts.
  - `Auto_cellranger_alcindor_bam.sh` and `Auto_cellranger_carroll_bam.sh`: rerun Cell Ranger with `--create-bam=true` for Numbat pileup while otherwise matching the historical `cellranger count` logic.
  - `Auto_stage_validate_scatlas_cellranger_outputs.R/.sh`: stage new `filtered_feature_bc_matrix` outputs into historical `matrix_all/<sample>_filtered` structure, optionally export dense count CSVs using the original `write.sh` logic, and require exact sparse-matrix identity to the live historical matrices.
  - `Auto_submit_scatlas_raw_rebuild.sh`: submits download jobs first, then dependent Cell Ranger jobs.
- Raw FASTQs and BAM-producing Cell Ranger outputs are staged under `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/` by explicit user request for this workflow.
- scATLAS Numbat scripts now live in `analysis/cnv/`:
  - `Auto_scatlas_numbat_export_inputs.R`: writes `ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv`.
  - `Auto_prepare_scatlas_numbat_container.sh`, `Auto_run_scatlas_numbat_pileup.sh`, `Auto_scatlas_numbat_run_sample.R`, and `Auto_run_scatlas_numbat_sample.sh`: mirror the PDO Numbat settings (`max_iter=2`, `gamma=20`, `init_k=3`, `min_cells=50`).
  - `Auto_scatlas_numbat_conservative_recut.R`: default direct validation cut is `SCATLAS_NUMBAT_CONSERVATIVE_N_CUT=3`, with minor clones below `max(20 cells, 3%)` merged into major clones.
  - `Auto_00_submit_scatlas_numbat.sh`: submits the full Numbat dependency chain, but now refuses to run until `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/validation/Auto_scatlas_cellranger_matrix_validation.csv` exists and contains only `ok` rows.
- Methodology files:
  - `analysis/methodology/raw_data/scatlas_raw_redownload_numbat_methodology.md`
  - `analysis/methodology/cnv/scatlas_numbat_methodology.md`
####################

####################
####################
## 22 May 2026 Conference Poster Requested Revisions

- `analysis/publication/poster_requested_revisions.R`
  - Purpose: poster-specific replot layer for the requested Single Cell Biology 2026 A0 poster revisions. It creates minimal SVG workflow schematics, revised atlas/NMF/enrichment/TME/PDO/spatial/FLOT/survival/targeting panels, and copies final assets into `ref_outs/Auto_conference_poster_plan/assets/publication/`.
  - Inputs: final cached scRef outputs, Visium mapping outputs, TME interaction CSVs, TCGA TPM/meta and Cox results, and PDO crossdata/FLOT/drug-reversal outputs.
  - Outputs: `ref_outs/publication/requested_revisions/` plus poster-ready assets in `ref_outs/publication/assets/` and `ref_outs/Auto_conference_poster_plan/assets/publication/`.
  - Redundant assets: no files are deleted; removed poster panels are copied with `delete_` prefixes and listed in `ref_outs/publication/requested_revisions/tables/delete_redundant_asset_manifest.csv`.
  - Methodology: `analysis/methodology/publication/poster_requested_revisions_methodology.md`.
- Current live poster HTML: `ref_outs/Auto_conference_poster_plan/oac_poster_final.html`; older poster drafts are listed in `ref_outs/Auto_conference_poster_plan/legacy_html_manifest.md`.
- 23 May 2026 update: requested-revision plotting behaviour was moved into the canonical `analysis/publication/poster_section*.R` scripts. `analysis/publication/run_poster_publication_figures.sh` no longer runs `poster_requested_revisions.R`, so manually curated poster assets such as `assets/Schematic_overall.svg`, `assets/Schematic_Anno.svg`, and the merged atlas UMAP/barplot are not overwritten by the wrapper.
####################
## 1 Jun 2026 scATLAS Raw Redownload And Numbat Validation Notes

- `analysis/raw_data/Auto_cellranger_carroll_bam.sh`
  - Carroll EGA FASTQs include multiple flowcell-specific sample prefixes for some biological samples, e.g. `<sample>_HGHWKBGXH` and `<sample>_HGHYKBGXH`. Cell Ranger 8 requires `--sample` when multiple prefixes are present in one FASTQ folder; the wrapper derives all prefixes from filenames and passes them as a comma-separated `--sample` value.
  - The Carroll rebuild now processes all historical `EAC-` and `BARR-` samples by default so outputs can be validated exactly against the 54 live `Carroll_2023/matrix_all/*_filtered` directories. Numbat input export is the later step that restricts to tumour samples.
- `analysis/raw_data/Auto_cellranger_carroll_bam_single.sh`
  - Single-sample recovery wrapper for cancelled/failed Carroll array elements. It uses the same Cell Ranger binary, GRCh38-2024-A transcriptome, FASTQ symlink layout, comma-separated multi-flowcell `--sample` logic, and `--create-bam=true` as the array wrapper.
- Validation gate:
  - `analysis/cnv/Auto_00_submit_scatlas_numbat.sh` must not be run until `analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.sh` has written `Auto_scatlas_cellranger_matrix_validation.csv` with all rows `status == "ok"` and no validation-failures CSV.
- Alcindor SRA recovery:
  - `analysis/raw_data/Auto_download_alcindor_srr_array.sh` is a per-accession recovery wrapper for unfinished SRR downloads. It keeps `fasterq-dump --split-files --include-technical` identical to the sequential downloader and uses `pigz` only to speed compression.
  - `analysis/raw_data/Auto_cellranger_alcindor_bam.sh` stages non-destructive 10x-style symlinks under `Alcindor_2025/fastq_cellranger/<SRR>/` because Cell Ranger 8 does not accept raw `fasterq-dump` names like `<SRR>_1.fastq.gz` with `--sample=<SRR>`.
####################

####################
## 25 Mar 2026 Non-Malignant MP Cross-Celltype Correlation Script

- `analysis/non_malignant_nmf/mp_cross_celltype_correlations.R`
  - Purpose: build a full cross-celltype MP co-occurrence network across the available compartments (`cancer`, `macrophage`, `fibroblast`, `endothelial`, `nk`, `plasma`, `cd4`, `cd8`) using the complete `EAC_Ref_merged_strict.rds` atlas, then visualise positive and negative associations in a Fig. 5a-style network layout.
  - MP filtering: applies the standard silhouette filter to every compartment, then keeps MPs with positive-sample coverage `> 5` at the active cutoff before pairwise correlation testing.
  - Adjusted-score rule: a cell is MP-positive when `UCell > cutoff`; default cutoff is `0.25` (not `0.5`), because `0.5` is overly sparse for the current UCell score ranges. The script writes a cutoff-sensitivity CSV/PDF so this can be re-tuned.
  - Cancer compartment: uses `ref_outs/EAC_Ref_merged_strict.rds` directly, with cancer cells defined as the malignant epithelial subset identified from `ref_outs/meta_full_epi.rds` (`malignancy %in% c("malignant_level_1", "malignant_level_2")`).
  - Cancer MP labels: user-facing outputs replace raw cancer MP IDs with the curated cancer MP descriptions, while non-malignant compartments remain labeled as `MP# + celltype`.
  - Correlation universe: computes Pearson and Spearman correlations between all MP pairs from different cell types, not only cancer-versus-TME pairs. For each cell-type pair, only studies with at least 10 shared samples are retained.
  - Ligand-receptor rule: uses the `All.Pairs` sheet from `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`; retains only `Pair.Evidence %in% c("literature supported", "putative")` and removes all `EXCLUDED*` rows before matching. Positive edges are annotated when a ligand or receptor is in the top 4,000 ranked genes of one node and the paired receptor or ligand is in the connected node MP genes.
  - Inputs: `ref_outs/EAC_Ref_merged_strict.rds`, `ref_outs/meta_full_epi.rds`, `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`, `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`, per-celltype non-malignant outputs `ref_outs/nmf_*/MP_outs_default.rds`, `UCell_default.rds`, and ligand-receptor pairs from `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`
  - Caching: large intermediate objects are cached in `ref_outs/non_malignant_mp_correlations/cache/` as `Auto_cross_celltype_step1_compartment_cache.rds`, `Auto_cross_celltype_step2_correlation_cache.rds`, and `Auto_cross_celltype_step3_lr_cache.rds`; set environment variable `AUTO_MPXCELL_FORCE_REBUILD=TRUE` to ignore existing caches and rebuild.
  - Outputs: `ref_outs/non_malignant_mp_correlations/` containing per-compartment adjusted-score CSVs, `Auto_cross_celltype_node_summary.csv`, `Auto_cross_celltype_cutoff_sensitivity.csv/pdf`, shared-sample summaries, all/positive/negative correlation CSVs, a focal-celltype-paged bubble PDF, network PDFs, LR pair tables, `Auto_cross_celltype_positive_edges_lr_annotated.csv`, an LR-annotated positive network PDF, and one formatted focal-celltype workbook (`Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx`).
  - Summary output: `updates/new_updates/summaries/Auto_mp_cross_celltype_correlations_summary.csv`
  - Ligand-receptor note: prioritises `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`; if unavailable, the script falls back to local candidate files or `RAMILOWSKI_LR_PATH`, and otherwise writes `Auto_cross_celltype_ligand_receptor_status.csv` with `missing`.
  - Methodology: `analysis/methodology/non_malignant_nmf/mp_cross_celltype_correlations_methodology.md`

- `analysis/non_malignant_nmf/mp_cross_celltype_correlations.sh`
  - PBS wrapper for `mp_cross_celltype_correlations.R` (`ncpus=8`, `mem=128gb`, `walltime=8h`, `dmtcp` env) for the direct `EAC_Ref_merged_strict.rds` run; optional PBS variable `cutoff` is forwarded as the first R argument (for example `qsub -v cutoff=0.2 analysis/non_malignant_nmf/mp_cross_celltype_correlations.sh`).
####################
####################
## 25 Mar 2026 GEO Bulk Survival Update

- `analysis/clinical/geo_survival_data_prep.R`
  - Purpose: download and prepare external GEO bulk-expression cohorts for MP/state survival analysis.
  - Current datasets:
    - `GSE19417` (public GEO survival metadata available)
    - `GSE13898` (public clinicopathology only; no public GEO survival metadata)
  - GEO downloads are cached under `ref_outs/geo_survival/raw/`.
  - Probe-level series matrices are collapsed to gene-symbol matrices using GEO platform annotation (`Gene symbol`) and a highest-variance probe-per-gene rule.
  - Outputs: `ref_outs/geo_survival/Auto_geo_survival_dataset_manifest.csv/.rds`, per-dataset metadata RDS/CSV, gene-level expression RDS, and probe-map CSVs.

- `analysis/clinical/geo_survival_mp_state_survival.R`
  - Purpose: run GEO bulk MP/state survival analysis using the `tcga_mp_state_survival_reg_noreg.R` structure adapted for GEO bulk input.
  - Shared loops retained:
    - DGE vs reference overlap plots
    - `split_method` loop: `continuous`, `median`, `q1q4`
    - `mode` argument parsing (currently forced to `noreg`, matching the current default downstream choice)
  - GEO-adapted 4-method design:
    - `full_cohort_reference`
    - `eac_only_reference`
    - `full_cohort_dge`
    - `eac_only_dge`
  - Output directory: `ref_outs/geo_survival/geo_task2_survival/`
  - Volcano layout:
    - one combined EAC-only PDF with MP pages first and State pages second
    - one page per split method
    - each page is side-by-side: left = `eac_only_reference`, right = `eac_only_dge`
  - Summary CSV: `updates/new_updates/summaries/Auto_geo_survival_clinical_mps_v2_reg_noreg_summary.csv`
- `analysis/clinical/geo_survival_mp_state_survival.sh`
  - PBS wrapper for `geo_survival_mp_state_survival.R`
  - Resources: `ncpus=8`, `mem=96gb`, `walltime=06:00:00`, env `dmtcp`

- GEO dataset-specific notes:
  - `GSE19417`
    - Platform: `GPL4372`
    - Public supplement includes `SurvivalDays` and `SurvivalCensoring_Kaplan-Meier`
    - Survival-ready subset = `Oesophageal adenocarcinoma` rows with non-missing survival
- `GSE13898`
  - Platform: `GPL6102`
  - GEO public files include clinicopathology and expression only
  - Prepared for reference but not marked `analysis_ready_for_survival`

####################
####################
## 30 Mar 2026 Cross-platform Bulk QC + Integrated Survival

- `analysis/clinical/bulk_tcga_geo_qc.R`
  - Purpose: harmonize TCGA whole-bulk RNA-seq and GEO `GSE19417` bulk microarray for cross-platform QC before pooled survival analysis.
  - Harmonization rules:
    - shared genes only (`TCGA ∩ GEO`)
    - TCGA transform = `log2(TPM + 1)`
    - GEO transform = keep supplied processed log-scale matrix
    - per-dataset standardization = row-wise z-score
  - QC strategy:
    - PCA on the top variable shared genes after harmonization
    - per-dataset QC pages: one page for `TCGA`, one page for `GEO_GSE19417`
    - each page combines PCA, histology-consistency metrics, and expression-strength metrics (`median`, `IQR`, `breadth`)
    - QC is binary only (`Keep` / `Remove`); there is no separate review tier
    - automatic removal is moderately stringent: low-expression outliers are removed directly, and histology mismatches are removed when also inconsistent in local/global PCA structure
  - Outputs: `ref_outs/bulk_crossplatform/` with `Auto_bulk_crossplatform_qc_review.pdf`, before/after PCA PDFs, full sample QC table, removed/retained CSVs, preprocessing summary CSV, QC objects RDS, and summary CSV `updates/new_updates/summaries/Auto_bulk_tcga_geo_qc_summary.csv`

- `analysis/clinical/bulk_tcga_geo_integrated_survival.R`
  - Purpose: rerun MP/state survival after QC on a harmonized TCGA + GEO bulk matrix, rather than reusing platform-specific scores.
  - Scoring rules:
    - recompute scores directly from harmonized expression using gene-set mean expression
    - do not reuse previous GSVA outputs
    - include both reference gene sets and DGE-derived gene sets
  - Retained 4-method structure:
    - `full_cohort_reference`
    - `eac_only_reference`
    - `full_cohort_dge`
    - `eac_only_dge`
  - Survival outputs:
    - per-dataset continuous/median/q1q4 Cox
    - pooled Cox with `dataset` covariate
    - continuous score-by-dataset interaction tests for direction heterogeneity
    - one combined EAC-only volcano PDF with MP pages first and State pages second
    - one page per split method
    - each page is side-by-side: left = `eac_only_reference`, right = `eac_only_dge`
  - Output directory: `ref_outs/bulk_crossplatform/survival/`
  - Key outputs: `Auto_bulk_crossplatform_survival_results.csv`, `Auto_bulk_crossplatform_direction_summary.csv`, `Auto_bulk_crossplatform_interaction_results.csv`, `Auto_bulk_crossplatform_survival_volcano_eac_only.pdf`, MP score distribution PDF, and summary CSV `updates/new_updates/summaries/Auto_bulk_tcga_geo_integrated_survival_summary.csv`
- `analysis/clinical/bulk_tcga_geo_integrated_survival.sh`
  - PBS wrapper for `bulk_tcga_geo_integrated_survival.R`
  - Resources: `ncpus=8`, `mem=96gb`, `walltime=06:00:00`, env `dmtcp`

- `analysis/clinical/tcga_mp_state_survival_qc_filtered.R`
  - Purpose: rerun the TCGA whole-bulk survival volcano workflow on the QC-retained TCGA samples only.
  - Input filter: keeps `dataset == "TCGA"` and `integration_keep == TRUE` from `ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_sample_table.csv`
  - Methods retained:
    - `malignant_filtered_reference`
    - `malignant_filtered_dge`
    - `whole_filtered_reference`
    - `whole_filtered_dge`
  - Volcano layout:
    - one combined PDF with MP pages first and State pages second
    - one page per split method
    - 2 x 2 layout: malignant reference, malignant DGE, whole-bulk reference, whole-bulk DGE
  - Outputs: `ref_outs/task2_filtered_survival/Auto_task2_filtered_survival_volcano_methods_reg_noreg.pdf`, filtered Cox CSV, and summary CSV `updates/new_updates/summaries/Auto_survival_clinical_mps_v2_reg_noreg_filtered_summary.csv`
- `analysis/clinical/tcga_mp_state_survival_qc_filtered.sh`
  - PBS wrapper for `tcga_mp_state_survival_qc_filtered.R`
  - Resources: `ncpus=8`, `mem=128gb`, `walltime=08:00:00`, env `dmtcp`
- `analysis/clinical/tcga_mp_state_survival_reg_noreg.sh`
  - PBS wrapper for `tcga_mp_state_survival_reg_noreg.R`
  - Resources: `ncpus=8`, `mem=128gb`, `walltime=08:00:00`, env `dmtcp`

- `analysis/metaprograms/mp_3ca_ucell_scoring.sh`
  - PBS wrapper for `mp_3ca_ucell_scoring.R`
  - Resources: `ncpus=8` (reduced to `ncores=2` in R), `mem=128gb`, `walltime=08:00:00`, env `dmtcp`

####################
####################
## 25 Mar 2026 External Epithelial MP UCell Heatmap

- `analysis/developmental/external_epi_mp_ucell_heatmap.R`
  - Purpose: score the filtered scATLAS nMP=19 metaprograms with UCell in three external epithelial datasets, then summarise mean MP activity per matched epithelial cell type using the same adult-epithelium / Barretts row structure as the developmental enrichment workflow.
  - Environment: `dmtcp` (uses `UCell`)
  - Inputs:
    - `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Stomach/data_9_9_annotated_seurat_all_ut.rds`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Barretts/alldatahighquality.rds`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/metadata/EoE_meta.txt`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE.mtx`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE_cell.tsv`
    - `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/expression/63f53992d91a88956d36dc4f/EoE_gene.tsv`
  - Outputs:
    - `ref_outs/Auto_external_epi_mp_ucell/Auto_external_epi_mp_ucell_summary.csv`
    - `ref_outs/Auto_external_epi_mp_ucell/Auto_external_epi_mp_ucell_summary.rds`
    - `ref_outs/Auto_external_epi_mp_ucell/Auto_external_epi_mp_ucell_heatmap_adult_epithelium.pdf`
    - `ref_outs/Auto_external_epi_mp_ucell/Auto_external_epi_mp_ucell_heatmap_barretts_oesophagus.pdf`
    - `ref_outs/Auto_external_epi_mp_ucell/Auto_external_epi_mp_ucell_heatmap_combined.pdf`
    - `updates/new_updates/summaries/Auto_external_epi_mp_ucell_summary.csv`
  - Adult oesophagus note: the source is a very large Matrix Market file (~393k cells), so the script samples up to `max_cells_per_type` cells per epithelial label before UCell scoring and caches the filtered subset under `ref_outs/Auto_external_epi_mp_ucell/Auto_cache/`.
  - Barretts note: row mapping uses `cell_type_secondary` so the output can retain the original fine categories (`Suprabasal_Dividing`, `Undifferentiated_Dividing`, `Endocrine_NEUROG3`, submucosal gland terms, etc.).
####################

####################
## 28 Apr 2026 State And MP Ordering Standard

All new plots, heatmaps, summaries, and statistical tables that display finalized epithelial states or scATLAS MPs must use the biologically defined order below rather than alphabetical ordering.

- Defined state order:
  1. `Classic Proliferative`
  2. `Basal to Intestinal Metaplasia`
  3. `SMG-like Metaplasia`
  4. `Stress-adaptive`
  5. `Immune Infiltrating`
  6. 3CA EMT / pan-cancer relabelled states, with `3CA_EMT_and_Protein_maturation` first when present
- `Unresolved` and `Hybrid` are technical labels, not defined biological states. Only include them when the analysis explicitly needs them, and append them after the defined-state order.
- MP order follows state order with cell-cycle MPs first. Within each block, use the `geneNMF.metaprograms$programs.tree$order` / `programs.clusters` tree order, after filtering out negative-silhouette MPs.
- Cell-cycle MPs: `MP1`, `MP7`, `MP9`.
- For top-MP assignment used as a state-like transcriptomic label, exclude cell-cycle MPs and assign the maximum from the z-normalised Approach-B noreg MP score matrix (`Auto_topmp_v2_noreg_mp_adj.rds`), matching the state-assignment score scale.
- Reuse the established colours from `analysis/cell_states/state_mp_sample_abundance.R` for MP and state labels whenever those labels are plotted across samples.
####################

####################
## Nature-Figure Publication Skill (Selective Application)

A `nature-figure` skill is installed at `/rds/general/user/sg3723/home/nature-skills/nature-figure/`. It enforces Nature-journal visual standards. **Agents must exercise judgment to apply this skill primarily to scripts producing final, sharable results.**

### When to Apply (Clever Selection)

Do **not** apply this to every R script. Focus on scripts that synthesize data across samples or produce "Final" visualizations. Prioritize:
- Final cohort-level summary plots (abundance, survival, clinical associations)
- Cross-dataset comparison figures
- Any active terminal script producing figures explicitly intended for manuscript inclusion or presentation slides

### How to Apply (R Backend — PDF Priority)

1. **Figure contract**: Define the claim and evidence hierarchy first.
2. **Typography**: Use 6.5pt Arial (Nature standard) via `theme_nature_contract()`.
3. **Export policy (PDF Priority)**: **PDF is the preferred format.** SVG is not required unless requested. Use `grDevices::cairo_pdf()` to ensure font embedding.
   ```r
   save_pub_pdf <- function(plot, filename, width_mm = 183, height_mm = 120) {
     w <- width_mm / 25.4; h <- height_mm / 25.4
     grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
     if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
       ComplexHeatmap::draw(plot, merge_legend = TRUE)
     } else {
       print(plot)
     }
     dev.off()
   }
   ```
4. **Color & IA**: Use restrained palettes and follow the **overview → deviation → relationship** information architecture.

### Reference Files
- `~/nature-skills/nature-figure/SKILL.md` — full skill specification
- `~/nature-skills/nature-figure/references/r-workflow.md` — R-specific patterns

### Exceptions
- **Diagnostic/QC scripts**: Step 1-6 pipeline outputs, internal QC heatmaps, and debugging plots should use standard Seurat/ggplot2 defaults to save time.
- **Development/Test scripts**: `delete_*.R` scripts.
####################

####################
## 15 May 2026 TCGA Reconstruction And Gender Validation Workflow

- `analysis/TCGA/tcga_esca_reconstruct_data.R`
  - Purpose: reconstruct the deleted TCGA-ESCA bulk RNA-seq input set without using the old spatial TCGA directory. The script queries the current public GDC API for open `TCGA-ESCA` RNA-seq `STAR - Counts` files, downloads each file under `ref_outs/TCGA/esca_gdc_reconstruction/raw/gdc_files/`, verifies file size and MD5 checksum, pulls cBioPortal `esca_tcga_gdc` patient/sample clinical attributes, and builds a gene-symbol TPM matrix.
  - Outputs: `ref_outs/TCGA/esca_gdc_reconstruction/` with `raw/`, `intermediate/`, `tables/`, and `logs/`; compatibility copies at `ref_outs/tcga_esca_meta.rds` and `ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`.
  - PBS wrapper: `analysis/TCGA/tcga_esca_reconstruct_data.sh` (`ncpus=2`, `mem=32gb`, `walltime=12h`, `#PBS -koed`, `dmtcp` env).
  - Reuse controls: verified files are reused; `SCREF_TCGA_SKIP_DOWNLOAD=TRUE` processes existing downloads only; `SCREF_TCGA_OVERWRITE_BAD=TRUE` is required to replace a checksum-failed generated file.

- `analysis/TCGA/tcga_gender_state_mp_validation.R`
  - Purpose: validate scRef Female-vs-Male MP and final-state associations in TCGA EAC primary bulk RNA-seq. It computes scRef sample-level MP mean UCell/state proportions, TCGA log2(TPM+1) GSVA MP/state scores, Wilcoxon tests, Cliff's delta, and direction concordance.
  - Outputs: `ref_outs/TCGA/gender_validation/` with cached GSVA scores, feature statistics, concordance tables, and `Auto_tcga_gender_scRef_concordance.pdf/png`; compact AI-readable summary at `updates/new_updates/summaries/Auto_tcga_gender_scRef_concordance_summary.csv`.
  - PBS wrapper: `analysis/TCGA/tcga_gender_state_mp_validation.sh` (`ncpus=4`, `mem=64gb`, `walltime=4h`, `#PBS -koed`, `dmtcp` env).

- Methodology: `analysis/methodology/TCGA/tcga_reconstruction_and_gender_validation_methodology.md`
####################

####################
## 19 May 2026 CNA Subclone Expression Correlation (Merged)

- `analysis/cnv/cna_subclone_expression_correlation.R`
  - Purpose: single self-contained script for arm-level CNA profiling, dominant clone analysis, pairwise CNA-expression distances, consensus CNA clustering, OAC/OCCAMS/TCGA annotation, recurrent event testing, and presentation-quality visualisations. Merged from the previous v1 computation + v2 visualisation split.
  - Inputs: `ref_outs/Auto_malignant_subclone_mp/` CSVs, per-sample `_outs.rds` inferCNA matrices, `OAC_CNV.xlsx`, `41588_2018_331_MOESM3_ESM.xlsx`.
  - Outputs: `ref_outs/Auto_cna_subclone_expression/` with `figures/`, `tables/`, and `rds/Auto_cna_subclone_expression_results.rds`.
  - Compact summary: `updates/new_updates/summaries/Auto_cna_subclone_expression_summary.csv`.
  - Cache: saves intermediate RDS after computation; supports `SCREF_REPLOT_ONLY=TRUE` to skip computation and replot from cache.
  - Key plotting decisions: consensus heatmap hides row names and cluster row labels; recurrent CNA event summaries use recurrent + next ranked events split into MPs, six states (excluding Hybrid/Unresolved), and QC/CNA metrics; significance with group-level BH FDR point size and stars; largest-subclone standardized boxplots; chr8q/MYC 3-group comparison; pairwise CNA-distance sample-centered Spearman.
  - Supersedes: `analysis/cnv/legacy_cna_subclone_expression_visuals_v2.R` (retained as legacy reference).
  - Methodology: `analysis/methodology/cnv/cna_subclone_expression_correlation_methodology.md`.
####################

####################
## 21 May 2026 Conference Poster Publication Replot Layer

- New publication-only scripts live under `analysis/publication/` and should be used for A0 poster/manuscript-style figure replots rather than editing upstream biological analysis scripts.
- Shared helper: `analysis/publication/publication_helpers.R`
  - Defines canonical poster state order, state colours from `analysis/shared/scRef_config.R`, MP order, MP-to-state grouping, PDF/PNG export helpers, and placeholder/status helpers.
- Section scripts:
  - `poster_section1_atlas_metaprograms.R`
  - `poster_section2_genetics_regulons.R`
  - `poster_section3_tme_interactions.R`
  - `poster_section4_pdo_concordance.R`
  - `poster_section5_lineage_validation.R`
  - `poster_section6_flot_remodelling.R`
  - `poster_section7_survival_targeting.R`
- Wrapper: `analysis/publication/run_poster_publication_figures.sh` (`ncpus=4`, `mem=32gb`, `walltime=4h`, `#PBS -koed`, `dmtcp` env). It reruns all section scripts and copies generated assets into `ref_outs/Auto_conference_poster_plan/assets/publication/`.
- Dependency wrapper: `analysis/publication/replot_after_scenic.sh` (`ncpus=2`, `mem=16gb`, `walltime=1h`, `#PBS -koed`, `dmtcp` env). Submit with `qsub -W depend=afterok:<scenic_jobid>` to refresh the SCENIC publication panel after `final_mp_scenic.sh` completes.
- Outputs:
  - Sectioned figures/tables/logs under `ref_outs/publication/<section>/`.
  - Bulk poster-ready copies under `ref_outs/publication/assets/` and `ref_outs/Auto_conference_poster_plan/assets/publication/`.
- Known current placeholders/status files:
  - SCENIC regulon network: upstream `ref_outs/final_mp_scenic/` was absent.
  - scRef pseudotime/state-distance node plot: upstream `ref_outs/state_distance_pseudotime/` and `task6_hybrid_pairwise_distance/` were absent.
  - Visium validation: no final Visium state mapping PNG/PDF was found under `ref_outs`.
  - TCGA/GEO survival: no final survival/Cox/volcano output was found under `ref_outs`.
- Methodology: `analysis/methodology/publication/poster_figures_methodology.md`.
####################

####################
## 15 May 2026 MP Database Correlation Barretts Update

- `analysis/metaprograms/mp_database_correlation.R`
  - Purpose update: regenerate MP-vs-reference UCell correlation heatmaps using the same retained-MP tree order as `analysis/enrichment/enrichment_annotation.R`.
  - New reference included: `Barretts_Oesophagus` from `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/enrich_dev_Barretts_Oesophagus.rds`.
  - Self-contained enrichment behavior: if `ref_outs/cluster_enrich.rds` is missing or does not contain all requested databases, the script rebuilds it using the same GO/Hallmark/3CA/custom enrichment setup as `enrichment_annotation.R`.
  - UCell behavior: large Barretts signatures require `maxRank >= 5000`; the script now sets this automatically and saves `ref_outs/UCell_ref_terms_v2_MP19.rds` after each scoring batch.
  - Outputs: `ref_outs/cluster_enrich.rds`, `ref_outs/UCell_ref_terms_v2_MP19.rds`, `ref_outs/Auto_MP_correlation_heatmaps_v2_MP19.pdf`, `ref_outs/Auto_MP_correlation_results_v2_MP19.rds`, and `ref_outs/Auto_MP_database_reference_gene_lists_Adult_Epithelium_Barretts_Oesophagus.xlsx`.

- `analysis/metaprograms/mp_database_correlation.sh`
  - PBS wrapper for the regeneration (`ncpus=4`, `mem=128gb`, `walltime=8h`, `#PBS -koed`, `dmtcp` env).
####################

####################
## 27 May 2026 TCGA CNA Recurrent Event Validation

- `analysis/TCGA/tcga_cna_recurrent_event_validation.R`
  - Purpose: validate the top 8 scRef recurrent chromosome-arm CNA events in TCGA EAC bulk RNA-seq using the GDC segment file `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg`.
  - Arm calling: computes length-weighted segment means per hg38 chromosome arm and matches segment IDs like `TCGA-2H-A9GF-01` to RNA-seq barcodes like `TCGA-2H-A9GF-01A` using the first 15 TCGA barcode characters.
  - Thresholding: scans absolute arm mean thresholds 0.05-0.30 against curated `OAC_CNV.xlsx` arm events; current selected threshold is `0.12`.
  - Annotation workbook interpretation:
    - `OAC_CNV.xlsx` is a curated gain/loss arm-driver table with rank, cytoband, genes, approximate frequency, pathway, and clinical relevance.
    - `41588_2018_331_MOESM3_ESM.xlsx` ST1/ST2 are GISTIC peak sheets; ST3/ST4 are deletion driver sheets; ST5/ST6 are amplification driver sheets.
    - The parser locates the real `Gene`/`hgnc_symbol` header row in driver sheets instead of relying on a fixed skip count. ST5 row 25 `ENSG00000136997 / MYC*` is retained as high-confidence `MYC` on `gain_chr8q`.
  - Plotting: mirrors `analysis/cnv/cna_subclone_expression_correlation.R` presentation style: same MP/state order, separate MP/state pages, top 8 scRef recurrent events in dotplots, top 4 in boxplots, standardized event deltas, and large slide-readable labels.
  - Outputs: `ref_outs/TCGA/cna_recurrent_event_validation/` with `intermediate/`, `tables/`, `figures/`, and `logs/`; compact summary at `updates/new_updates/summaries/Auto_tcga_cna_recurrent_event_validation_summary.csv`.
  - Current result: 87 matched TCGA EAC primary samples; 0 top-scRef recurrent event TCGA MP/state associations at feature-type FDR < 0.10; 5 TCGA-discovered associations at feature-type FDR < 0.10, none recurrent in scRef. Conclusion: TCGA supports the scRef trend that recurrent arm CNAs are not robustly coupled to MP/state programmes.
  - Methodology: `analysis/methodology/TCGA/tcga_cna_recurrent_event_validation_methodology.md`.
####################

####################
## 28 May 2026 State UMAP Dispersion And Co-Localisation

- `analysis/cell_states/state_umap_dispersion_colocalisation.R`
  - Purpose: calculates per-sample UMAP dispersion and same-label nearest-neighbour co-localisation for the five primary noreg Approach B states, excluding `Unresolved` and `Hybrid`, then repeats the analysis within `Basal to Intestinal Metaplasia` using basal-only UMAPs labelled by top basal MP.
  - Basal MP labels: use the top z-normalised noreg MP among `MP17`, `MP14`, `MP5`, `MP10`, and `MP8` for each basal cell.
  - Inputs: `ref_outs/EAC_Ref_epi.rds`, `ref_outs/Auto_topmp_v2_noreg_states_B.rds`, `ref_outs/Auto_topmp_v2_noreg_mp_adj.rds`, and `ref_outs/state_distance_pseudotime/sample_state_trajectories/*_pseudotime_states.rds`.
  - Regeneration: if trajectory RDS files are missing, the script reruns `analysis/cell_states/pseudotime_state_distance_matrix.R` first.
  - Outputs: `ref_outs/state_umap_dispersion_colocalisation/` with `intermediate/`, `tables/`, `figures/`, and `logs/`; compact summary at `updates/new_updates/summaries/Auto_state_umap_dispersion_colocalisation_summary.csv`.
  - PBS wrapper: `analysis/cell_states/state_umap_dispersion_colocalisation.sh` (`ncpus=8`, `mem=128gb`, `walltime=12h`, `#PBS -koed`, `dmtcp` env).
  - Methodology: `analysis/methodology/cell_states/state_umap_dispersion_colocalisation_methodology.md`.
####################

####################
## 28 May 2026 Unified Developmental MP Validation

- `analysis/developmental/developmental_mp_enrichment_unified.R`
  - Purpose: unified developmental validation workflow for scATLAS metaprogrammes using three evidence layers: original `clusterProfiler::enricher()` overlap enrichment, original epithelial scATLAS expression-correlation enrichment, and original-style external annotated reference-celltype MP UCell scoring when annotated expression matrices are available.
  - Runtime core: The script now natively contains the original-script-aligned method logic directly alongside the ranked-reference construction.
  - Ranked gene handling: rebuilds ranked developmental reference genes from the source workbooks and writes `ref_outs/Auto_developmental_mp_enrichment_unified/tables/Auto_developmental_reference_rank_audit.csv`. Top-50 mode uses the first 50 ranked genes per term, then runs the same original enrichment/correlation logic. All-gene mode uses the original per-stage `TERM2GENE` references and validates against the original outputs.
  - Rank sources: Early embryogenesis/adult oesophagus/Barretts use source workbook marker order with available DE statistics; Organogenesis `S1D` uses marker columns labelled as DEGs ordered by z-score; Normal development long uses `fold.change`/`qval`; Normal development short is literature-marker order and is flagged as not a differential ranked list.
  - Validation outputs: `Auto_developmental_validation_overlap_all_vs_original.csv` and `Auto_developmental_validation_correlation_all_vs_original.csv` should be checked after reruns; the all-gene values are expected to have zero difference from `cluster_enrich.rds` and `Auto_MP_correlation_results_v2_MP19.rds`.
  - Method 3 availability: early embryogenesis processed Seurat data were downloaded from the paper dataset portal; organogenesis raw counts/cell/gene annotations were downloaded from GEO `GSE157329`; adult stomach uses the local developmental `Stomach.rds`. Normal fetal development, adult oesophagus, and Barretts were not scored in the 28 May 2026 run because directly usable annotated expression matrices were not available locally/directly.
  - Outputs: `ref_outs/Auto_developmental_mp_enrichment_unified/` with `intermediate/`, `tables/`, `figures/`, and `logs/`. Main PDFs are `Auto_developmental_mp_unified_top50.pdf` and `Auto_developmental_mp_top50_vs_all_overlap_correlation.pdf`.
  - Cache controls: `SCREF_FORCE_REBUILD=TRUE` rebuilds UCell caches; `SCREF_REPLOT_ONLY=TRUE` replots from existing tables; `SCREF_UCELL_CORES` and `SCREF_MAX_CELLS_PER_TYPE` control runtime.
  - PBS wrapper: `analysis/developmental/developmental_mp_enrichment_unified.sh` (`ncpus=8`, `mem=128gb`, `walltime=12h`, `#PBS -koed`, `dmtcp` env).
  - Methodology: `analysis/methodology/developmental/developmental_mp_enrichment_unified_methodology.md`.
####################

####################
## 4 Jun 2026 Developmental External Reference Downloads

- `analysis/developmental/developmental_mp_enrichment_unified.R` now scores additional processed annotated method-3 references:
  - Adult oesophagus: Broad SCP1242 downloaded under `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/` using the user-supplied temporary `generate_curl_config` URL and `curl -K cfg.txt`; `cfg.txt` is retained for provenance. The workflow streams a sampled subset from `EoE.mtx` using `EoE_meta.txt$cell_type_anno`.
  - Normal development stomach: Descartes/Fred Hutch direct public downloads under `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/`: `Stomach_gene_count.RDS`, `df_cell.RDS`, and `df_gene.RDS`. Only stomach cells are scored, with `Main_cluster_name` mapped to the 16 stomach terms in both long and short normal-development references.
  - Barretts oesophagus: Esophagus Cancer Cell Atlas high-quality combined RDS downloaded to `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/barretts/alldatahighquality.rds`, scored by `cell_type_secondary`.
- New audit output: `ref_outs/Auto_developmental_mp_enrichment_unified/tables/Auto_developmental_reference_celltype_coverage.csv` confirms expected external annotation terms are present and scored for each dataset source.
####################

####################
## 7 Jun 2026 MP Refinement Sub-MP Workflow

- `analysis/metaprograms/mp_refinement_submp.R`
  - Purpose: optional refinement of nMP19 metaprogrammes by keeping MPs with silhouette `>= 0.2`, removing MPs with silhouette `< 0`, and splitting MPs with silhouette `> 0` and `< 0.2`.
  - Split selection: evaluates mean silhouette for `k = 2..10`; if the raw optimum is `> 6`, the selected split is capped to `k = 5`.
  - Critical implementation detail: GeneNMF's `normVector()` normalizes by vector sum, not L2 norm. Sub-MP gene-list derivation must keep this behavior to avoid near-empty sub-MP signatures and weak/incorrect UCell correlations.
  - Outputs: `ref_outs/Metaprogrammes_Results/mp_refinement/` with `intermediate/`, `tables/`, `figures/`, and `logs/`. Key figures are `mp_splitting_diagnostics.pdf`, `refined_mp_correlation_heatmap.pdf`, and `refined_mp_jaccard_heatmap.pdf`.
  - Current run: 39 final state-associated refined MPs in the terminal heatmaps; all refined gene sets are non-empty. The refined UCell cache is `75348 x 39`.
  - Methodology: `analysis/methodology/metaprograms/mp_refinement_methodology.md`.
####################

####################
## 25 Jun 2026 scATLAS RNA Velocity Workflow

- `analysis/trajectory/scatlas_velocity_submit.sh` submits the full raw-BAM velocity chain for scATLAS samples with BAMs under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files`.
- Workflow structure mirrors the PDO velocity scripts: metadata/barcode export, GRCh38/RepeatMasker reference prep, epithelial barcode BAM filtering/sorting, velocyto loom generation, scVelo visualisation, and R nodeplots.
- Direction derivation uses the five primary noreg Approach B states only: `Classic Proliferative`, `Basal to Intestinal Metaplasia`, `SMG-like Metaplasia`, `Stress-adaptive`, and `Immune Infiltrating`. Final states are retained for UMAP coloring.
- Outputs live under `ref_outs/Auto_velocity_scATLAS/` with `tables/`, `figures/`, `logs/`, `h5ad/`, `looms/`, `coord/`, `barcodes/`, and `ref/` folders. Lightweight update summary is `updates/new_updates/summaries/Auto_scatlas_velocity_direction_summary.csv`.
- Methodology: `analysis/methodology/trajectory/scatlas_velocity_methodology.md`.
####################
