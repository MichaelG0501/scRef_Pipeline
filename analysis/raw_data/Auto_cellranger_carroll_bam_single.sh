#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/raw_data/Auto_cellranger_carroll_bam_single.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N Auto_Car_CR_single
#PBS -koed

set -euo pipefail

####################
# Single-sample Carroll 2023 Cell Ranger rerun wrapper. This is used for
# targeted recovery of failed array elements while keeping the same FASTQ
# grouping, Cell Ranger version, transcriptome, and BAM-enabled count command.
####################

echo $(date +%T)
module purge
module load tools/dev

RAW_ROOT="${RAW_ROOT:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files}"
STUDY_DIR="${RAW_ROOT}/Carroll_2023"
FASTQ_BY_SAMPLE="${STUDY_DIR}/fastq_by_sample"
CELLRANGER_ROOT="${STUDY_DIR}/cellranger"
LOG_DIR="${STUDY_DIR}/logs"
CELLRANGER_BIN="${CELLRANGER_BIN:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-8.0.1/bin/cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A}"

if [[ -z "${sample:-}" ]]; then
  echo "ERROR: set sample=<Carroll sample name> with qsub -v"
  exit 1
fi
if [[ ! -x "$CELLRANGER_BIN" ]]; then
  echo "ERROR: missing Cell Ranger executable: $CELLRANGER_BIN"
  exit 1
fi
if [[ ! -d "$TRANSCRIPTOME" ]]; then
  echo "ERROR: missing Cell Ranger transcriptome: $TRANSCRIPTOME"
  exit 1
fi
if [[ ! -d "${FASTQ_BY_SAMPLE}/${sample}" ]]; then
  echo "ERROR: missing sample FASTQ directory: ${FASTQ_BY_SAMPLE}/${sample}"
  exit 1
fi

mkdir -p "$CELLRANGER_ROOT" "$LOG_DIR"

mapfile -t sample_prefixes < <(
  find "${FASTQ_BY_SAMPLE}/${sample}" -maxdepth 1 -type l -name '*.fastq.gz' -printf '%f\n' |
    sed -E 's/^(.+)_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz$/\1/' |
    sort -u
)
if [[ "${#sample_prefixes[@]}" -eq 0 ]]; then
  echo "ERROR: no Cell Ranger sample prefixes found in ${FASTQ_BY_SAMPLE}/${sample}"
  exit 1
fi
sample_arg=$(IFS=,; echo "${sample_prefixes[*]}")
echo "Using Cell Ranger --sample=${sample_arg}"

cd "$CELLRANGER_ROOT"
if [[ -s "${CELLRANGER_ROOT}/${sample}/outs/possorted_genome_bam.bam" ]]; then
  echo "Existing BAM found for ${sample}; skipping Cell Ranger."
  echo $(date +%T)
  exit 0
fi

"$CELLRANGER_BIN" count \
  --id="$sample" \
  --fastqs="${FASTQ_BY_SAMPLE}/${sample}" \
  --sample="$sample_arg" \
  --transcriptome="$TRANSCRIPTOME" \
  --create-bam=true \
  --chemistry=auto \
  --localcores=8 \
  --localmem=120

if [[ ! -s "${CELLRANGER_ROOT}/${sample}/outs/possorted_genome_bam.bam" ]]; then
  echo "ERROR: Cell Ranger BAM was not created for ${sample}"
  exit 1
fi

echo $(date +%T)
