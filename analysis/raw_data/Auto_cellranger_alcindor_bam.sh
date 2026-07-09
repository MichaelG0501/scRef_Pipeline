#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N Auto_Alc_CR
#PBS -J 1-20
#PBS -koed

set -euo pipefail

####################
# Re-run Cell Ranger for Alcindor 2025 with BAM creation enabled so Numbat can
# perform SNP pileup/phasing from the resulting possorted BAM.
####################

echo $(date +%T)
module purge
module load tools/dev

RAW_ROOT="${RAW_ROOT:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files}"
STUDY_DIR="${RAW_ROOT}/Alcindor_2025"
FASTQ_DIR="${STUDY_DIR}/fastq"
CR_FASTQ_ROOT="${STUDY_DIR}/fastq_cellranger"
CELLRANGER_ROOT="${STUDY_DIR}/cellranger"
LOG_DIR="${STUDY_DIR}/logs"
CELLRANGER_BIN="${CELLRANGER_BIN:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-8.0.1/bin/cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A}"

mkdir -p "$CELLRANGER_ROOT" "$LOG_DIR" "$CR_FASTQ_ROOT"

ACCESSION_FILE="${STUDY_DIR}/Auto_Alcindor_2025_SRR_accessions.txt"
if [[ ! -s "$ACCESSION_FILE" ]]; then
  for n in $(seq 27335925 27335944); do
    echo "SRR${n}"
  done > "$ACCESSION_FILE"
fi

sample=$(sed -n "${PBS_ARRAY_INDEX}p" "$ACCESSION_FILE")
if [[ -z "$sample" ]]; then
  echo "ERROR: no sample for PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX}"
  exit 1
fi

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
if [[ ! -d "${FASTQ_DIR}/${sample}" ]]; then
  echo "ERROR: missing FASTQ directory: ${FASTQ_DIR}/${sample}"
  exit 1
fi
for read_no in 1 2; do
  src="${FASTQ_DIR}/${sample}/${sample}_${read_no}.fastq.gz"
  if [[ ! -s "$src" ]]; then
    echo "ERROR: missing FASTQ file: $src"
    exit 1
  fi
done

####################
# fasterq-dump creates SRR-style names (<SRR>_1.fastq.gz and
# <SRR>_2.fastq.gz), but Cell Ranger expects 10x-style FASTQ names. Stage
# symlinks without modifying the downloaded FASTQs.
####################
cr_fastq_dir="${CR_FASTQ_ROOT}/${sample}"
mkdir -p "$cr_fastq_dir"
for read_no in 1 2; do
  src="${FASTQ_DIR}/${sample}/${sample}_${read_no}.fastq.gz"
  if [[ "$read_no" == "1" ]]; then
    dst="${cr_fastq_dir}/${sample}_S1_L001_R1_001.fastq.gz"
  else
    dst="${cr_fastq_dir}/${sample}_S1_L001_R2_001.fastq.gz"
  fi
  if [[ -L "$dst" ]]; then
    current_target=$(readlink "$dst")
    if [[ "$current_target" != "$src" ]]; then
      echo "ERROR: existing symlink points elsewhere: $dst -> $current_target"
      exit 1
    fi
  elif [[ -e "$dst" ]]; then
    echo "ERROR: existing non-symlink staged FASTQ would be overwritten: $dst"
    exit 1
  else
    ln -s "$src" "$dst"
  fi
done

cd "$CELLRANGER_ROOT"
if [[ -s "${CELLRANGER_ROOT}/${sample}/outs/possorted_genome_bam.bam" ]]; then
  echo "Existing BAM found for ${sample}; skipping Cell Ranger."
  echo $(date +%T)
  exit 0
fi

"$CELLRANGER_BIN" count \
  --id="$sample" \
  --fastqs="$cr_fastq_dir" \
  --sample="$sample" \
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
