####################
# scATLAS RNA Velocity Methodology
####################

## Purpose

This workflow runs RNA velocity for the subset of scATLAS epithelial samples that have raw Cell Ranger BAM files available under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files`. It mirrors the PDO velocity workflow structure: metadata export, raw BAM filtering, velocyto loom generation, scVelo modelling, and state-direction node plots.

## Inputs

- `ref_outs/EAC_Ref_epi.rds`: merged epithelial Seurat object with UMAP coordinates and `orig.ident`.
- `ref_outs/Auto_final_states.rds`: finalized state labels for plotting.
- `ref_outs/Auto_topmp_v2_noreg_states_B.rds`: noreg Approach B state calls used for five-primary-state direction scoring.
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files/*_possorted_genome_bam.bam`: raw Cell Ranger BAMs for the raw-data subset.
- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz`: GRCh38 gene annotation for velocyto.
- UCSC hg38 RepeatMasker table or a pre-existing `repeatmasker.hg38.gtf`, used as the velocyto repeat mask.

## Workflow

1. `scatlas_velocity_metadata.R` loads the epithelial object and state vectors, keeps only samples whose raw BAM exists, exports per-cell metadata, and writes one QC barcode file per sample.
2. `scatlas_velocity_prepare_refs.py` expands the Cell Ranger gene GTF and prepares `repeatmasker.hg38.gtf`.
3. `scatlas_velocity_filter_sort.sh` filters each raw BAM to epithelial QC barcodes and coordinate-sorts the result.
4. `scatlas_velocity_run_velocyto.sh` runs velocyto on each filtered BAM to produce one loom per sample.
5. `scatlas_velocity_scvelo_visualise.py` joins loom barcodes back to scRef cell IDs, runs stochastic scVelo per sample, saves h5ad caches, and computes source-state velocity alignment toward target-state centroids.
6. `scatlas_velocity_nodeplots.R` redraws grouped node plots from the scVelo CSV tables.

## State Direction Logic

Direction scoring is restricted to the five primary noreg Approach B states:

- Classic Proliferative
- Basal to Intestinal Metaplasia
- SMG-like Metaplasia
- Stress-adaptive
- Immune Infiltrating

For each sample and source state, the workflow computes the mean velocity vector in UMAP space. For every target state present in the sample, it computes the cosine alignment between the source-state mean velocity vector and the vector from source-state centroid to target-state centroid. Positive values indicate velocity alignment from the source toward that target.

## Outputs

- `ref_outs/Auto_velocity_scATLAS/barcodes/`: per-sample QC barcode lists.
- `ref_outs/Auto_velocity_scATLAS/coord/`: filtered coordinate-sorted BAMs and BAI files.
- `ref_outs/Auto_velocity_scATLAS/looms/`: velocyto loom outputs.
- `ref_outs/Auto_velocity_scATLAS/h5ad/`: cached per-sample scVelo objects.
- `ref_outs/Auto_velocity_scATLAS/tables/`: metadata, sample manifest, state nodes, direction edges, top target per source, grouped summaries, and failure logs.
- `ref_outs/Auto_velocity_scATLAS/figures/`: per-sample velocity visualisations and state-direction node plots.
- `updates/new_updates/summaries/Auto_scatlas_velocity_direction_summary.csv`: lightweight summary for progress updates.

## Run Command

Submit the full dependency chain:

```bash
/opt/pbs/bin/qsub analysis/trajectory/scatlas_velocity_submit.sh
```

The submitter uses the `dmtcp` conda environment for R metadata/node plotting and the `velocity` environment for velocyto/scVelo.
