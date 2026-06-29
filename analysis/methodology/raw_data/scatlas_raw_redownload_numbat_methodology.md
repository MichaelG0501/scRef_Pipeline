####################
# scATLAS Raw Redownload For Numbat Methodology
####################

Status: active upstream methodology.

Purpose: rebuild raw FASTQ and BAM-producing Cell Ranger outputs for the Carroll 2023 and Alcindor 2025 scATLAS datasets so Numbat can run SNP pileup and haplotype-aware subclone calling.

Inputs:
- Alcindor 2025 SRA runs: `SRR27335925` through `SRR27335944`.
- Carroll 2023 EGA dataset: `EGAD00001009401`.
- Carroll EGA credential JSON supplied at runtime by `EGA_CREDENTIAL_JSON`; the scripts default to the existing local credential path but never copy credentials into the repository.
- Cell Ranger executable supplied by `CELLRANGER_BIN`.
- Cell Ranger GRCh38 transcriptome supplied by `TRANSCRIPTOME`.

Outputs:
- `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Alcindor_2025/fastq/`
- `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Alcindor_2025/cellranger/`
- `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Carroll_2023/ega_download/`
- `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Carroll_2023/fastq_by_sample/`
- `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Carroll_2023/cellranger/`

Operational logic:
- Alcindor FASTQs are downloaded with `fasterq-dump --split-files --include-technical`, matching the original SRR fetch approach. The sequential wrapper uses the original gzip-style compression; `Auto_download_alcindor_srr_array.sh` can be used for recovery/acceleration of unfinished accessions and only changes compression implementation to `pigz` when available, not FASTQ content.
- Alcindor Cell Ranger uses non-destructive symlink staging under `Alcindor_2025/fastq_cellranger/<SRR>/` because `fasterq-dump` emits `<SRR>_1.fastq.gz` and `<SRR>_2.fastq.gz`, while Cell Ranger 8 requires 10x-style FASTQ names such as `<SRR>_S1_L001_R1_001.fastq.gz` when `--sample=<SRR>` is supplied.
- Carroll FASTQs are downloaded with `pyega3 fetch EGAD00001009401`.
- Carroll FASTQs are organised by symlink into per-sample directories. Existing downloaded FASTQs are not moved. Samples sequenced across multiple flowcells retain their flowcell-specific FASTQ prefixes, and Cell Ranger is called with a comma-separated `--sample` list of all prefixes for that sample.
- Cell Ranger is rerun with `--create-bam=true`, `--chemistry=auto`, and GRCh38 2024-A by default, matching the original Cell Ranger version and reference assumptions while producing `possorted_genome_bam.bam`.
- Apart from BAM creation and the explicit multi-flowcell `--sample` selection required by Cell Ranger 8, Cell Ranger commands follow the historical processing logic from `00_scripts`: sample directories are named by the original accession/sample prefix, `cellranger count` is run with `--id=<sample>`, `--fastqs=<sample_dir>`, `--transcriptome=<GRCh38-2024-A>`, and `--chemistry=auto`.
- The default Carroll processing regex keeps all historical `EAC-` and `BARR-` samples so the staged outputs can be validated against all 54 live Carroll `matrix_all` directories. The downstream Numbat manifest then restricts subclone analysis to the tumour-relevant samples.
- `Auto_stage_validate_scatlas_cellranger_outputs.R` mirrors the old post-processing by staging `filtered_feature_bc_matrix` outputs as `matrix_all/<sample>_filtered`. It then reads both the new staged matrix and the live historical matrix with `Read10X()` and requires exact sparse-matrix identity before downstream Numbat submission.
- Dense count CSVs equivalent to `00_counts_matrix_all` can be exported under the staging root by setting `SCATLAS_EXPORT_COUNT_CSV=TRUE` when running the validation wrapper.

Run command:
- Submit the full dependency chain with `qsub analysis/raw_data/Auto_submit_scatlas_raw_rebuild.sh`.
