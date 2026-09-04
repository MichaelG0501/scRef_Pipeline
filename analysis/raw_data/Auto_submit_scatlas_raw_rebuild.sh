#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/raw_data/Auto_submit_scatlas_raw_rebuild.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=01:00:00
#PBS -N Auto_raw_submit
#PBS -koed

set -euo pipefail

####################
# Convenience submitter for the scATLAS raw-data rebuild. It submits downloads
# first, then BAM-producing Cell Ranger jobs with afterok dependencies.
####################

echo $(date +%T)
module purge
module load tools/dev

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
RAW_ROOT="${RAW_ROOT:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files}"
EGA_CREDENTIAL_JSON="${EGA_CREDENTIAL_JSON:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Carroll_2023/ega.json}"
CELLRANGER_BIN="${CELLRANGER_BIN:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-8.0.1/bin/cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A}"
LOG_DIR="${RAW_ROOT}/logs/Auto_raw_rebuild_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"
cd "$WD"

qsub_bin="/opt/pbs/bin/qsub"

jid_alc=$("$qsub_bin" -v RAW_ROOT="$RAW_ROOT" -o "${LOG_DIR}/Auto_download_alcindor_srr.log" -e "${LOG_DIR}/Auto_download_alcindor_srr.err" analysis/raw_data/Auto_download_alcindor_srr.sh)
echo "Submitted Alcindor SRA download: ${jid_alc}"

jid_car=$("$qsub_bin" -v RAW_ROOT="$RAW_ROOT",EGA_CREDENTIAL_JSON="$EGA_CREDENTIAL_JSON" -o "${LOG_DIR}/Auto_download_carroll_ega.log" -e "${LOG_DIR}/Auto_download_carroll_ega.err" analysis/raw_data/Auto_download_carroll_ega.sh)
echo "Submitted Carroll EGA download: ${jid_car}"

jid_alc_cr=$("$qsub_bin" -W depend=afterok:${jid_alc} -v RAW_ROOT="$RAW_ROOT",CELLRANGER_BIN="$CELLRANGER_BIN",TRANSCRIPTOME="$TRANSCRIPTOME" -o "${LOG_DIR}/Auto_cellranger_alcindor_bam.log" -e "${LOG_DIR}/Auto_cellranger_alcindor_bam.err" analysis/raw_data/Auto_cellranger_alcindor_bam.sh)
echo "Submitted Alcindor Cell Ranger array: ${jid_alc_cr}"

jid_car_cr=$("$qsub_bin" -W depend=afterok:${jid_car} -v RAW_ROOT="$RAW_ROOT",CELLRANGER_BIN="$CELLRANGER_BIN",TRANSCRIPTOME="$TRANSCRIPTOME" -o "${LOG_DIR}/Auto_cellranger_carroll_bam.log" -e "${LOG_DIR}/Auto_cellranger_carroll_bam.err" analysis/raw_data/Auto_cellranger_carroll_bam.sh)
echo "Submitted Carroll Cell Ranger array: ${jid_car_cr}"

echo $(date +%T)
