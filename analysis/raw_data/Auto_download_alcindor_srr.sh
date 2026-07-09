#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=72:00:00
#PBS -N Auto_Alc_SRR
#PBS -koed

set -euo pipefail

####################
# Download Alcindor 2025 raw FASTQs from SRA/ENA accessions for the scATLAS
# Numbat validation rebuild. Accessions are the 20 SRR runs already present in
# the original processed Alcindor dataset.
####################

echo $(date +%T)
module purge
module load tools/prod
module load SRA-Toolkit/3.0.3-gompi-2022a

RAW_ROOT="${RAW_ROOT:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files}"
STUDY_DIR="${RAW_ROOT}/Alcindor_2025"
FASTQ_DIR="${STUDY_DIR}/fastq"
LOG_DIR="${STUDY_DIR}/logs"
SRA_BIN="${SRA_BIN:-$(dirname "$(command -v fasterq-dump)")}"

mkdir -p "$FASTQ_DIR" "$LOG_DIR"

ACCESSION_FILE="${STUDY_DIR}/Auto_Alcindor_2025_SRR_accessions.txt"
if [[ ! -s "$ACCESSION_FILE" ]]; then
  for n in $(seq 27335925 27335944); do
    echo "SRR${n}"
  done > "$ACCESSION_FILE"
fi

if [[ ! -x "${SRA_BIN}/fasterq-dump" ]]; then
  echo "ERROR: fasterq-dump not found or not executable: ${SRA_BIN}/fasterq-dump"
  echo "Set SRA_BIN to the directory containing fasterq-dump."
  exit 1
fi

while read -r srr; do
  [[ -z "$srr" ]] && continue
  sample_dir="${FASTQ_DIR}/${srr}"
  done_file="${sample_dir}/Auto_${srr}_download_done.txt"
  mkdir -p "$sample_dir"
  if [[ -s "$done_file" ]] && compgen -G "${sample_dir}/${srr}*.fastq*" >/dev/null; then
    echo "Existing FASTQs found for ${srr}; skipping."
    continue
  fi
  echo "Downloading ${srr}"
  cd "$sample_dir"
  "${SRA_BIN}/fasterq-dump" "$srr" --split-files --include-technical -p -e 4
  gzip -f ./*.fastq
  printf "sample\t%s\nfinished\t%s\n" "$srr" "$(date -Is)" > "$done_file"
done < "$ACCESSION_FILE"

echo $(date +%T)
