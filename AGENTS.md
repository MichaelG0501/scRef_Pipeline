# AGENTS.md — scRef_Pipeline

Single-cell RNA-seq QC and analysis pipeline for OAC public reference datasets on Imperial College PBS Pro. Core computation is R with bash PBS wrappers. Read this file and `AI_RULES.md` before acting.

## Repository and pipeline

```text
*.R / N_<Step>.sh   core steps 1–8 and PBS jobs
*_master.sh         per-sample PBS orchestration
analysis/           downstream analyses
analysis/ANALYSIS_MAP.md  canonical analysis inventory/dependency map
analysis/methodology/     operational methods for complex workflows
ref_outs/           pipeline and analysis outputs (gitignored)
temp/               PBS logs
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

Minimal wrapper:

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
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd "$WD"
Rscript <script>.R
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
- Plot text, legends, rows/columns, points, and line widths must be readable on PowerPoint slides. Increase dimensions or paginate instead of shrinking labels.
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
