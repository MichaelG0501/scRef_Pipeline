####################
# scATLAS Numbat Subclone Methodology
####################

Status: active workflow methodology.

Purpose: perform haplotype-aware Numbat subclone analysis for Carroll 2023 tumour scRNA-seq samples and Alcindor 2025 scRNA-seq samples using the same core approach as the PDO Numbat workflow.

Inputs:
- `ref_outs/EAC_Ref_epi.rds` merged epithelial scATLAS Seurat object; sample-specific matrices are exported by splitting the raw count layer by `orig.ident`.
- BAM-producing Cell Ranger outputs staged under `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/<dataset>/cellranger/<raw_sample>/outs/`.
- Numbat container `docker://pkharchenkolab/numbat-rbase:latest`.

Output structure:
- `ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv`
- `ref_outs/Auto_scatlas_numbat/by_samples/<sample>/input/`
- `ref_outs/Auto_scatlas_numbat/by_samples/<sample>/<sample>_allele_counts.tsv.gz`
- `ref_outs/Auto_scatlas_numbat/by_samples/<sample>/numbat/`
- `ref_outs/Auto_scatlas_numbat/conservative_clones/`

Numbat settings:
- `lambdas_ref = ref_hca`
- `genome = "hg38"`
- `max_iter = 2`
- `t = 1e-5`
- `gamma = 20`
- `init_k = 3`
- `min_cells = 50`
- `call_clonal_loh = TRUE`
- `ncores = 12`

Terminal no-subclone handling:
- Samples for which Numbat returns terminal biological filter statuses such as `No clones remain after filtering by size` or `No CNV remains after filtering by LLR in pseudobulks` are retained as valid no-subclone outcomes.
- For these samples, `Auto_scatlas_numbat_run_sample.R` writes a done file, an empty clone/joint summary, and an RDS summary with `terminal_no_subclone = TRUE` rather than failing the dependency chain.
- `Auto_scatlas_numbat_conservative_recut.R` records these samples as `terminal_no_subclone` in the final conservative-clone summary and does not attempt to cut a missing tree.

Conservative clone layer:
- Raw Numbat outputs are retained untouched.
- `Auto_scatlas_numbat_conservative_recut.R` reads the final Numbat tree and re-cuts with `SCATLAS_NUMBAT_CONSERVATIVE_N_CUT`, default `3`.
- Clones below `max(20 cells, 3% of cells)` are merged into the best major clone by posterior probability where possible.
- This is the default validation layer to avoid over-fragmented sample-specific subclone calls.

Run command:
- After raw FASTQ and BAM Cell Ranger rebuild outputs exist, submit with `bash analysis/cnv/Auto_00_submit_scatlas_numbat.sh`.
