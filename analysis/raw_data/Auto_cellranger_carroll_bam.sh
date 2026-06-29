#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N Auto_Car_CR
#PBS -J 1-54
#PBS -koed

set -euo pipefail

####################
# Organise downloaded Carroll 2023 EGA FASTQs by 10x sample name using
# symlinks, then re-run Cell Ranger with BAM creation enabled. By default this
# processes all historical EAC/BARR samples that feed Carroll_2023 matrix_all.
####################

echo $(date +%T)
module purge
module load tools/dev

RAW_ROOT="${RAW_ROOT:-/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat}"
STUDY_DIR="${RAW_ROOT}/Carroll_2023"
DOWNLOAD_DIR="${STUDY_DIR}/ega_download"
FASTQ_BY_SAMPLE="${STUDY_DIR}/fastq_by_sample"
CELLRANGER_ROOT="${STUDY_DIR}/cellranger"
LOG_DIR="${STUDY_DIR}/logs"
MANIFEST="${STUDY_DIR}/Auto_Carroll_2023_sample_manifest.txt"
CARROLL_SAMPLE_REGEX="${CARROLL_SAMPLE_REGEX:-^(EAC-|BARR-)}"
CELLRANGER_BIN="${CELLRANGER_BIN:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-8.0.1/bin/cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A}"

mkdir -p "$FASTQ_BY_SAMPLE" "$CELLRANGER_ROOT" "$LOG_DIR"

if [[ ! -x "$CELLRANGER_BIN" ]]; then
  echo "ERROR: missing Cell Ranger executable: $CELLRANGER_BIN"
  echo "Set CELLRANGER_BIN before qsub."
  exit 1
fi
if [[ ! -d "$TRANSCRIPTOME" ]]; then
  echo "ERROR: missing Cell Ranger transcriptome: $TRANSCRIPTOME"
  echo "Set TRANSCRIPTOME before qsub."
  exit 1
fi
if [[ ! -d "$DOWNLOAD_DIR" ]]; then
  echo "ERROR: missing EGA download directory: $DOWNLOAD_DIR"
  exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
  tmp_manifest="${MANIFEST}.tmp"
  : > "$tmp_manifest"
  while IFS= read -r fq; do
    base=$(basename "$fq")
    sample=$(printf "%s\n" "$base" | sed -E 's/^(.+)_H[^_]+_S[0-9]+_.+\.fastq\.gz$/\1/')
    if [[ "$sample" == "$base" ]]; then
      continue
    fi
    if [[ ! "$sample" =~ $CARROLL_SAMPLE_REGEX ]]; then
      continue
    fi
    mkdir -p "${FASTQ_BY_SAMPLE}/${sample}"
    ln -sf "$fq" "${FASTQ_BY_SAMPLE}/${sample}/${base}"
    echo "$sample" >> "$tmp_manifest"
  done < <(find "$DOWNLOAD_DIR" -type f -name '*.fastq.gz' | sort)
  sort -u "$tmp_manifest" > "$MANIFEST"
fi

sample=$(sed -n "${PBS_ARRAY_INDEX}p" "$MANIFEST")
if [[ -z "$sample" ]]; then
  echo "No Carroll sample for PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX}; exiting cleanly."
  echo $(date +%T)
  exit 0
fi

if [[ ! -d "${FASTQ_BY_SAMPLE}/${sample}" ]]; then
  echo "ERROR: missing sample FASTQ directory: ${FASTQ_BY_SAMPLE}/${sample}"
  exit 1
fi

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
