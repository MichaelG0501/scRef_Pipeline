# AGENTS.md — scRef_Pipeline

Single-cell RNA-seq QC and analysis pipeline for OAC public reference datasets on Imperial College PBS Pro. Core computation is R with bash PBS wrappers. Read this file and `AI_RULES.md` before acting.

## Repository and pipeline

```text
*.R / N_<Step>.sh   core steps 1–8 and PBS jobs
*_master.sh         per-sample PBS orchestration
analysis/           downstream analyses and tracked essential setup scripts
job_scripts/        local submission-only wrappers mirroring analysis/ (gitignored)
analysis/ANALYSIS_MAP.md  canonical analysis inventory/dependency map
analysis/methodology/     operational methods for complex workflows
ref_outs/           pipeline and analysis outputs (gitignored)
temp/               historical PBS logs; new scheduler logs go in the project root
```

| Step | Entry point | R script | Scope | Env |
|---|---|---|---|---|
| 1 | `1_QC_Pipeline.sh` | `QC_Pipeline.R` | all samples | dmtcp |
| 2 | `2_master.sh` | `Clustering.R` | per sample | dmtcp |
| 3 | `3_Annotation.sh` | `Annotation.R` | all samples | dmtcp |
| 4 | `4_Expr_filtering.sh` | `Expr_filtering.R` | all samples | dmtcp |
| 5 | `5_master.sh` | `InferCNA.R` | per sample | dmtcp |
| 6 | `6_master.sh` | `Malignancy.R` | per sample | dmtcp |
| 7 | `7_master.sh` | `NMF.R` | per sample | dmtcp |
| 8 | `8_geneNMF.sh` | `geneNMF.R` | all samples | gnmf |

## Mandatory HPC rules

1. Never run analytical, memory-heavy, CPU-heavy, or substantial I/O work on a login node. Submit it through PBS. Login-node checks must be lightweight (headers, file existence/size, CSV/TXT summaries, syntax parsing, job status).
2. Every PBS job uses `#PBS -koed` for live logging and has explicit CPU, memory, walltime, and name directives.
3. Initialize conda exactly before activation:

   ```bash
   eval "$(~/miniforge3/bin/conda shell.bash hook)"
   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
   ```

4. Use `dmtcp` for general Seurat/analysis and `gnmf` for GeneNMF.
5. PBS binaries may not be on `PATH`; use `/opt/pbs/bin/qsub`, `/opt/pbs/bin/qstat`, and `/opt/pbs/bin/qdel`.
6. Do not exceed 46 concurrent jobs. Per-sample submitters must throttle with the established `qstat | grep sg3723` loop.
7. If a task requires submission, submit it yourself and monitor until exit status and expected outputs are verified. A zero exit code alone is insufficient.
8. Strictly forbid fallbacks: must use only the first and best option for any analysis. If not available or if it fails, the script should stop immediately.

####################
## Commit-ready analysis code and figure organization

- Keep `analysis/` limited to important, maintained scientific workflows and indispensable shared helpers or setup dependencies intended for Git. Consolidate small plotting extensions into their main analysis when practical. Do not keep a separate persistent script for each repair or verification step.
- One-off preflight checks, audits, reference replays, validation probes, PDF render checks, and audit-only package installers are temporary session tools. Run them inline or place them outside `analysis/` in a Git-ignored project-local scratch directory such as `Auto_session_checks/`. They must not become required production dependencies or appear as active pipeline steps in the analysis inventory. Keep useful validation results in the output `reports/` tier.
- This rule takes precedence over older placement or naming guidance for temporary checks. A package installer is an essential tracked setup dependency only when the scientific workflow itself needs it; an audit-only renderer does not qualify. Never delete files without explicit user permission. If deletion is prohibited, retain temporary tools in the ignored session archive outside `analysis/`.
- Organize figures by purpose: a small, explicitly selected `Auto_final/` set; `Auto_supporting/` for secondary evidence; `Auto_diagnostics/` for QC; and, when needed, `Auto_archive/` for retained duplicates or obsolete exports. Use clear descriptive filenames, order final figures numerically, and provide a compact figure index explaining each file's role. Do not present all generated PDFs as necessary final figures.
- Avoid exporting both combined panels and duplicate standalone panels by default. Keep scientific source data and tests reproducible independently of temporary auditing tools. Before handoff, check that production scripts do not source temporary files and that the proposed Git changes contain no session-only tools.
####################

####################
## Visualization readability for presentations and papers

- Wherever possible and space allows, use generous, clearly readable font sizes for figures intended for PowerPoint presentations and inclusion in papers. Apply this to titles, axis titles, tick labels, legends, annotations, panel labels, and heatmap row/column labels.
- Judge readability at the intended final display or printed size, including after a figure is placed on a slide or reduced to a manuscript column. A large export canvas alone does not guarantee readable text after scaling.
- Use available space to enlarge text and improve spacing without overlapping or clipped labels and legends. Increasing font size must never cause a figure to be split across additional pages. Preserve the intended comparison on one page; for clinical associations, keep all MPs on one page per clinical variable, as in OCCAMS. Adjust layout, label rotation/wrapping, margins, canvas dimensions and font sizes within that page. Only split pages when the user explicitly requests it, not to accommodate larger text. Keep point sizes and line widths clear at the intended final size.
- Choose font sizes appropriate to the output dimensions and intended use; follow journal requirements when supplied. Check the exported figure at its intended presentation or publication size before handoff.
####################

####################
## Stable current-result directories

- Each maintained analysis has one stable, descriptive output directory under the project's canonical output root. Put current files directly in `figures/`, `tables/`, `reports/`, `logs/`, and, when needed, `intermediate/` and `reference/`. Do not create empty tiers without a purpose.
- Do not create timestamped `Auto_run_*` directories, numbered attempts, `v2`/`v3` result directories, or competing final/latest copies for routine reruns or plot repairs. Timestamps and job IDs belong in logs, not in the current result path.
- Authorized reruns update the analysis's documented current outputs. Keep a single selected final figure set and one compact current-results index with working relative links. Do not overwrite unrelated files or bypass an explicit user restriction on modification.
- Current reports describe the implemented method, results, validation and unresolved scientific limitations. Exclude obsolete run inventories, failed-attempt narratives, superseded job IDs and debugging chronology from the current results index and scientific reports.
- When previous attempts must be retained, move them to an ignored project-local session archive outside the active analysis output directory, such as `Auto_session_checks/<analysis>/Auto_archive/`. Do not leave old run directories mixed with current results. Never delete files without explicit permission, and do not move files outside the project without authorization.
- Keep indispensable source data and references in the current live output tiers. Current scripts, reports and downstream analyses must not depend on the session archive or temporary audit scripts.
- When reorganizing outputs, update the producing script, registry, methodology, analysis map and links together. Verify current file paths and figure/table integrity before handoff. This rule governs new work and scoped cleanup; it does not authorize unrelated migrations.
####################

####################
## Analysis job scripts and PBS logs

- Store every new submission-only shell wrapper for an R script under `analysis/` in the project-root `job_scripts/` directory. Create the directory and required subdirectories when needed. Mirror the R script's path relative to `analysis/` and keep its exact basename, replacing `.R` with `.sh`: `analysis/cell_states/example.R` uses `job_scripts/cell_states/example.sh`. Preserve deeper nesting. Do not create submission-only wrappers inside `analysis/`.
- Before running or rerunning an analysis, first check that matching `job_scripts/` path. Reuse the existing launcher after checking its R command, arguments, environment, resources, and log destinations; update it there when needed. If it is absent, create it at that exact path rather than inventing another name. This basename rule takes precedence over generic shell-file naming prefixes.
- Keep these launchers minimal: PBS directives, timestamps, standard module/conda activation, a change to the project root, and the R invocation. Record the resource requirements, environment, arguments, and run command in the R script's registry so the launcher can be recreated. Keep analytical logic in R; it must run directly on an allocated compute node without reading or depending on the launcher.
- Submit from the project root, using the wrapper's `job_scripts/...` path, and change to that same root inside the job. Use `#PBS -koed` for live logging and explicit absolute `#PBS -o <project_root>/` and `#PBS -e <project_root>/` destinations. PBS stdout/stderr belong directly in the project root, alongside `analysis/` and `job_scripts/`; never place them in `analysis/`, `job_scripts/`, or an analysis output/log subdirectory. This rule concerns scheduler stdout/stderr; R-generated scientific logs and reports may still use the documented output tiers.
- Ignore `/job_scripts/` in the project's `.gitignore`. These launchers are local, reconstructable execution aids and should not be staged or committed. Keep project-root PBS stdout/stderr untracked as well.
- Essential shell scripts that implement nontrivial environment/tool installation, required input preparation, Python/helper orchestration, or other indispensable workflow logic must remain in the appropriate `analysis/` subdirectory, tracked by Git and maintained with the analysis. For example, scRef's `analysis/cell_states/Auto_drug_reversal/Auto_setup_drug_reversal_env.sh` is an essential setup script. Standard module/conda activation alone is not an exception. Do not hide indispensable logic or helper scripts solely in ignored `job_scripts/`; a disposable launcher may call a tracked essential script.
- This policy governs future launcher creation and updates. Existing pipeline entry points and historical wrappers are not automatically moved, deleted, or renamed; migrate them only as part of an explicitly scoped change.
####################

Minimal wrapper:

```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/
#PBS -e /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd "$WD"
Rscript "analysis/<subdirectory>/<script>.R"
echo $(date +%T)
```

## Storage and file safety

- Project scripts and every persistent input/output belong under `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/`.
- Final RDS/CSV/PDF/PNG/SVG/XLSX files, downstream inputs, per-cell scores, state vectors, gene lists, enrichment results, normalized matrices, source data, run summaries, and anything needed to replot must be in live.
- Ephemeral storage is allowed only for exceptionally large, fully reconstructable same-script caches. If another script reads it, or it is needed to reproduce a figure without rerunning its producer, it must also be live.
- When uncertain, save to both locations and document which copy is canonical.
- Scripts creating an ephemeral path must call `dir.create(..., recursive=TRUE, showWarnings=FALSE)`.
- Never delete files unless the user explicitly requests deletion. Superseded scripts are retained with `legacy_` filenames; manual deletion candidates use `delete_`.
- Preserve user work in a dirty worktree. Do not overwrite unrelated edits.
- Test scripts are named `delete_<description>.R` and are removed only with explicit permission.

Before every `saveRDS`, `write.csv`, workbook save, or figure save, ask: is it consumed downstream, needed to replot, or critical source data? If yes, the canonical copy is live.

## Current MP and state contract

The authoritative details and run order are in `analysis/ANALYSIS_MAP.md`. New downstream work uses only the centred refined noreg workflow:

- MP genes: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
- MP UCell: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`
- MP grouping/descriptions: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`
- States: `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds`
- Adjusted MP matrix: `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds`
- Group maxima: `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_group_max.rds`
- Ranked current state markers: `ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv`, produced by `analysis/cell_states/final_state_marker_discovery.R`

The final panel has 17 MPs. Cell-cycle rows are MP1, MP5, and MP13+. State groups are `Classic proliferation`, `Squamous-to-intestinal`, `Glandular-to-intestinal`, `Stress-adaptive`, and `Cancer-cell immune mimicry`; technical states are `Unresolved` and `Hybrid`. Approach B thresholds are maximum group score `<0.5` for unresolved and top-minus-second group gap `<0.3` for hybrid.

MP2x and MP11c are excluded for coverage below three samples. MP18a is the explicit documented quality/coverage exclusion. Do not reintroduce them downstream.

Do not use `state_temp.rds`, `Auto_topmp_states.rds`, `Auto_cluster_states.rds`, `Auto_topmp_v2_*`, `Auto_final_states.rds`, uncentred `geneNMF_metaprograms_nMP_19.rds`, or `UCell_nMP19_filtered.rds` in new/current analyses. They are legacy comparison objects.

Sample identity is always `orig.ident`; never infer samples by trimming cell barcodes.

## Analysis script registry standard

Every active or terminal analytical R/Python script starts with an accurate registry containing:

- status (`active`, `terminal`, `legacy`, or `delete-candidate`);
- exact script path;
- concise operational description;
- detailed methodology path for complex workflows, or `Methodology: not required` with a reason for a simple deterministic script;
- exact live inputs and unavoidable external inputs/downloads;
- exact outputs grouped into `intermediate/`, `tables/`, `figures/`, `logs/`, or `reports/` for long workflows;
- cache/replot behavior and relevant environment variables;
- PBS/run command and conda environment.

A methodology must document implemented thresholds, choice rationale, assignment/statistical unit, filtering, cache semantics, limitations, and validation. Do not create thin methodology files that merely repeat a simple header.

New long-running scripts use clear output tiers and support `SCREF_FORCE_REBUILD=TRUE`; support `SCREF_REPLOT_ONLY=TRUE` when persistent intermediates permit plot-only reruns. Plotting scripts also write a compact (<100 KB) CSV/TXT result summary to `updates/new_updates/summaries/`.

Whenever a script is added, renamed, superseded, or changes inputs/outputs, update `analysis/ANALYSIS_MAP.md` in the same change. Add only genuinely recurring technical findings to this file; workflow-specific detail belongs in the script header, map, or methodology.

## R and shell conventions

- `library()` calls are at the top; do not use `require()`.
- Match the surrounding quote and pipe style. Use snake_case for variables/functions.
- Per-sample R scripts read `commandArgs(trailingOnly=TRUE)` and use the first argument as the sample.
- Filter MP collections before downstream use according to the current centred refinement outputs; historical uncentred analyses use silhouette `<0` filtering only for provenance.
- Use Seurat accessors where practical; metadata grouping uses `orig.ident`.
- Plot text, legends, rows/columns, points, and line widths must be readable on PowerPoint slides. Adjust dimensions, layout and font sizes without adding pages solely to enlarge text; preserve one page per intended comparison unless the user explicitly requests splitting.
- New code inserted into an existing file is enclosed in 20-hash comment blocks where the language permits it. Existing lines may be changed when the user explicitly authorizes the requested repair; avoid unrelated rewrites.
- Shell scripts start with `#!/bin/bash`, print start/end times, purge/load modules, initialize conda, and quote paths/variables.
- Sentinel files (`no_ref`, `no_epi`, `no_cell`, `no_cancer`) and existing guard/skip behavior must be preserved.

## Current complex-workflow references

- Centred GeneNMF/refinement: `analysis/methodology/metaprograms/centred_refinement_methodology.md`
- Centred state definition: `analysis/methodology/metaprograms/centred_refined_state_definition_noreg_methodology.md`
- Centred state markers: `analysis/methodology/cell_states/final_state_marker_discovery_methodology.md`
- Basal/SMG MP distances: `analysis/methodology/cell_states/basal_smg_mp_distance_methodology.md`
- Final-MP SCENIC: `analysis/methodology/cell_states/final_mp_scenic_methodology.md`
- Clinical/bulk: `analysis/methodology/clinical/clinical_bulk_and_association_methodology.md`
####################
- OCCAMS bulk survival/clinical associations: `analysis/methodology/clinical/Auto_occams_bulk_clinical_methodology.md`; current analyses use exactly the 281 subjects with `Annotation_Compatible_Below_50_Flag == no` and persistent scores/metadata under `ref_outs/OCCAMS/clinical/`.
####################
- CNV/Numbat: `analysis/methodology/cnv/`
- Visium HD final annotation: `analysis/methodology/spatial/visium_hd_final_annotation_methodology.md`
- RNA velocity: `analysis/methodology/trajectory/scatlas_velocity_methodology.md`
- Manuscript replot: `analysis/methodology/publication/oac_scatlas_paper_replot_methodology.md`

The Visium HD annotation method is threshold-sensitive; preserve its documented production settings and do not infer replacements from older `legacy_visiumhd` scripts.

The general Yates Visium/Xenium mappers require external `processed/*.h5ad` and
`sample_info.csv` inputs. The formerly configured spatial-transcriptomics
`Zenodo_upload` root was absent at the August 2026 audit; do not invent another
root—supply it explicitly when the data are restored.

For R SCENIC, use the RcisTarget-compatible mc9nr databases under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`. `final_mp_scenic.R` contains the required motif-annotation and sparse-filtering compatibility patches.

<!-- #################### -->
## OCCAMS bulk RNA-seq contract

- Controlled OCCAMS BAMs and the reconstructable raw featureCounts table remain under `/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/`; persistent metadata, mappings, QC, references, and the final compressed matrix are under `ref_outs/OCCAMS/`.
- The EGA credential JSON is external to Git at `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json` with mode `600`. OCCAMS download scripts reference this path; never embed or print its values.
- All 300 downloaded BAMs are GRCh37. There are 156 `chr`-named and 144 non-`chr`-named BAMs, handled with `analysis/OCCAMS/grch37_contig_aliases.csv` and pinned GENCODE v19 GRCh37.p13 annotation.
- The current downstream count matrix is `ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz`. Never use the historical ephemeral `counts/OCCAMS_RNAseq_raw_counts.txt`, which was quantified with a GRCh38 GTF against GRCh37 BAMs.
- Subject `90687ad9d09f7677ab6fc83f75b78bd3fed8a9934e071f071306007c79277dab` is retained but explicitly QC-flagged (30.11% annotation-compatible; 69.34% no-feature). Use `ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv` for sensitivity exclusions.
- Detailed run order and validation are in `analysis/methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md`.
<!-- #################### -->

## External operations

- Upload to Google Drive only when the user explicitly requests it. Load `rclone` and upload only to `gdrive:IMPERIAL/`.
- NotebookLM operations require `module load Mesa/23.1.4-GCCcore-12.3.0`; authentication setup requires an X11 session. Use it only when explicitly relevant.

## Git and status handoff

- `ref_outs/` and PBS stdout/stderr are ignored.
- Commit messages are short and lowercase when the user requests a commit.
- Report renamed legacy files, submitted job IDs/status, validation performed, and exact persistent outputs. Never claim a workflow is current when its executable inputs still point to legacy objects.
