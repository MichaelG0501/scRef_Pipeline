####################
# State UMAP Dispersion And Co-Localisation Methodology
####################

## Purpose

`analysis/cell_states/state_umap_dispersion_colocalisation.R` quantifies whether noreg Approach B epithelial states occupy compact or dispersed regions of per-sample UMAP space, and whether each state is locally self-colocalised with cells of the same label.

The same workflow is then repeated within the `Basal to Intestinal Metaplasia` state after rebuilding basal-only per-sample UMAPs and labelling basal cells by the top state-defining MP among `MP17`, `MP14`, `MP5`, `MP10`, and `MP8`.

## Inputs

- `ref_outs/EAC_Ref_epi.rds`
- `ref_outs/Auto_topmp_v2_noreg_states_B.rds`
- `ref_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- `ref_outs/state_distance_pseudotime/sample_state_trajectories/*_pseudotime_states.rds`

If the trajectory RDS files are missing, the script regenerates them by running `analysis/cell_states/pseudotime_state_distance_matrix.R` before continuing.

## State And MP Labels

The state-level analysis uses the five primary noreg Approach B states and excludes `Unresolved` and `Hybrid`.

For every primary-state cell, the script also records its top state-defining MP by selecting the maximum z-normalised noreg MP score within that state's MP group. The basal-focused analysis keeps only basal-to-intestinal-metaplasia cells and labels them by the top basal MP.

## UMAP Scope

State-level UMAPs are rebuilt separately for each sample using only cells present in the saved pseudotime trajectory output and assigned to one of the five primary states.

Basal MP UMAPs are rebuilt separately for each sample using only basal-to-intestinal-metaplasia cells. This creates a basal-only subclustering view rather than reusing the broader state UMAP.

Both UMAP scopes use the same Seurat preprocessing pattern as the pseudotime scripts: `NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`, and `RunUMAP`.

## Dispersion Metric

For each sample and label:

1. The label centroid is the mean UMAP coordinate of cells with that label.
2. Each cell's raw dispersion is its Euclidean distance from the label centroid.
3. The primary reported dispersion is this distance divided by the sample's median all-cell distance from the sample centroid.

This normalisation makes dispersion more comparable across samples whose UMAP coordinate scales differ.

The script also reports per-label convex hull area and normalises it by the sample-level UMAP hull area.

## Co-Localisation Metric

For each cell, the script finds its `k` nearest neighbours in the relevant UMAP, excluding the cell itself. The default is `k = 15` and can be changed with `SCREF_COLOCAL_K`.

The same-label neighbour score is:

`number of nearest neighbours with the same label / k`

Scores near `1` indicate strong local co-localisation. Scores near `0` indicate local mixing or scattering among other labels.

## Outputs

All analysis outputs are written under:

`ref_outs/state_umap_dispersion_colocalisation/`

Output tiers:

- `intermediate/`: cached UMAP coordinates and per-cell metric RDS objects
- `tables/`: per-cell, per-sample, and overall CSV summaries
- `figures/`: boxplots, dispersion-vs-colocalisation scatter plots, and per-sample UMAP pages
- `logs/`: run summary RDS/TXT files

A compact summary is also written to:

`updates/new_updates/summaries/Auto_state_umap_dispersion_colocalisation_summary.csv`

## Cache Controls

- `SCREF_FORCE_REBUILD=TRUE`: rebuild UMAP coordinates and metrics even when caches exist.
- `SCREF_REPLOT_ONLY=TRUE`: reuse cached metrics and regenerate final tables/figures only.
- `SCREF_COLOCAL_K=<integer>`: set nearest-neighbour count for co-localisation.
- `SCREF_STATE_UMAP_MIN_CELLS=<integer>`: minimum primary-state cells per sample for state-level UMAPs.
- `SCREF_BASAL_UMAP_MIN_CELLS=<integer>`: minimum basal cells per sample for basal-only top-MP UMAPs.

## Run Command

Use the `dmtcp` environment:

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cell_states/state_umap_dispersion_colocalisation.R
```

For the full cohort, use the PBS wrapper:

```bash
/opt/pbs/bin/qsub analysis/cell_states/state_umap_dispersion_colocalisation.sh
```
