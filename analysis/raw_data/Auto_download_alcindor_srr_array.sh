#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=24:00:00
#PBS -N Auto_Alc_SRR_array
#PBS -J 1-20
#PBS -koed

set -euo pipefail

####################
# Per-accession Alcindor 2025 SRA recovery downloader. It uses the same
# fasterq-dump arguments as the sequential downloader but compresses with pigz
# when available so recovery reruns can finish within scheduler walltime.
####################

echo $(date +%T)
module purge
module load tools/prod
module load SRA-Toolkit/3.0.3-gompi-2022a

RAW_ROOT="${RAW_ROOT:-/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat}"
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

srr=$(sed -n "${PBS_ARRAY_INDEX}p" "$ACCESSION_FILE")
if [[ -z "$srr" ]]; then
  echo "No SRR accession for PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX}; exiting cleanly."
  echo $(date +%T)
  exit 0
fi

sample_dir="${FASTQ_DIR}/${srr}"
done_file="${sample_dir}/Auto_${srr}_download_done.txt"
mkdir -p "$sample_dir"

if [[ -s "$done_file" ]] && compgen -G "${sample_dir}/${srr}*.fastq.gz" >/dev/null; then
  echo "Existing completed FASTQs found for ${srr}; skipping."
  echo $(date +%T)
  exit 0
fi

if compgen -G "${sample_dir}/${srr}*.fastq" >/dev/null || compgen -G "${sample_dir}/${srr}*.fastq.gz" >/dev/null; then
  echo "ERROR: partial FASTQs already exist for ${srr} without a done marker: ${sample_dir}"
  echo "Inspect this directory before recovery to avoid overwriting an in-progress download."
  exit 1
fi

if [[ ! -x "${SRA_BIN}/fasterq-dump" ]]; then
  echo "ERROR: fasterq-dump not found or not executable: ${SRA_BIN}/fasterq-dump"
  exit 1
fi

echo "Downloading ${srr}"
cd "$sample_dir"
"${SRA_BIN}/fasterq-dump" "$srr" --split-files --include-technical -p -e 4

if command -v pigz >/dev/null 2>&1; then
  pigz -p 4 -f ./*.fastq
else
  gzip -f ./*.fastq
fi

printf "sample\t%s\nfinished\t%s\n" "$srr" "$(date -Is)" > "$done_file"
echo $(date +%T)
