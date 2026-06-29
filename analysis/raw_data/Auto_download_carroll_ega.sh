#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=72:00:00
#PBS -N Auto_Car_EGA
#PBS -koed

set -euo pipefail

####################
# Download Carroll 2023 raw FASTQs from EGA dataset EGAD00001009401. The
# credential JSON is referenced by path and is never copied into this repo.
####################

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/egaenv

RAW_ROOT="${RAW_ROOT:-/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat}"
STUDY_DIR="${RAW_ROOT}/Carroll_2023"
DOWNLOAD_DIR="${STUDY_DIR}/ega_download"
LOG_DIR="${STUDY_DIR}/logs"
EGA_CREDENTIAL_JSON="${EGA_CREDENTIAL_JSON:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Carroll_2023/ega.json}"
EGA_DATASET="${EGA_DATASET:-EGAD00001009401}"

mkdir -p "$DOWNLOAD_DIR" "$LOG_DIR"

if [[ ! -s "$EGA_CREDENTIAL_JSON" ]]; then
  echo "ERROR: missing EGA credential JSON: $EGA_CREDENTIAL_JSON"
  echo "Set EGA_CREDENTIAL_JSON to a valid pyega3 credential file."
  exit 1
fi

cd "$DOWNLOAD_DIR"
pyega3 -cf "$EGA_CREDENTIAL_JSON" fetch "$EGA_DATASET"

printf "dataset\t%s\nfinished\t%s\n" "$EGA_DATASET" "$(date -Is)" > "${STUDY_DIR}/Auto_Carroll_2023_EGA_download_done.txt"
echo $(date +%T)
