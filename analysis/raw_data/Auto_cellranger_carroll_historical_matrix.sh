#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N Auto_Car_CR_hist
#PBS -J 1-54
#PBS -koed

set -euo pipefail

####################
# Re-run Carroll 2023 Cell Ranger in a separate output root using the visible
# historical processing logic from Carroll_2023/submit.sh:
# FASTQs grouped by sample directory, no --sample filter, chemistry auto, and
# --create-bam=false. Existing BAM-producing outputs are not modified.
####################

echo $(date +%T)
module purge
module load tools/dev

RAW_ROOT="${RAW_ROOT:-/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files}"
STUDY_DIR="${RAW_ROOT}/Carroll_2023"
DOWNLOAD_DIR="${STUDY_DIR}/ega_download"
FASTQ_BY_SAMPLE="${STUDY_DIR}/fastq_by_sample"
CELLRANGER_ROOT="${CARROLL_HIST_CELLRANGER_ROOT:-${STUDY_DIR}/cellranger_historical_matrix}"
LOG_DIR="${STUDY_DIR}/logs"
MANIFEST="${STUDY_DIR}/Auto_Carroll_2023_sample_manifest.txt"
CARROLL_SAMPLE_REGEX="${CARROLL_SAMPLE_REGEX:-^(EAC-|BARR-)}"
CELLRANGER_BIN="${CELLRANGER_BIN:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-9.0.1/bin/cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A}"
CELLRANGER_CHEMISTRY="${CARROLL_CELLRANGER_CHEMISTRY:-auto}"

mkdir -p "$FASTQ_BY_SAMPLE" "$CELLRANGER_ROOT" "$LOG_DIR"

if [[ ! -x "$CELLRANGER_BIN" ]]; then
  echo "ERROR: missing Cell Ranger executable: $CELLRANGER_BIN"
  exit 1
fi
if [[ ! -d "$TRANSCRIPTOME" ]]; then
  echo "ERROR: missing Cell Ranger transcriptome: $TRANSCRIPTOME"
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
    sample=$(printf "%s\n" "$base" | grep -oP '^.*(?=_H.*_S[0-9])' || true)
    if [[ -z "$sample" ]]; then
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

if [[ -n "${sample:-}" ]]; then
  carroll_sample="$sample"
else
  carroll_sample=$(sed -n "${PBS_ARRAY_INDEX}p" "$MANIFEST")
fi

if [[ -z "$carroll_sample" ]]; then
  echo "No Carroll sample for PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX:-unset}; exiting cleanly."
  echo $(date +%T)
  exit 0
fi
if [[ ! -d "${FASTQ_BY_SAMPLE}/${carroll_sample}" ]]; then
  echo "ERROR: missing sample FASTQ directory: ${FASTQ_BY_SAMPLE}/${carroll_sample}"
  exit 1
fi
cellranger_sample_prefixes=$(
  find "${FASTQ_BY_SAMPLE}/${carroll_sample}" -maxdepth 1 \( -type l -o -type f \) \
    | xargs -r -n 1 basename \
    | grep -oP '^.*(?=_S[0-9])' \
    | sort -u \
    | paste -sd, - || true
)
if [[ -z "$cellranger_sample_prefixes" ]]; then
  echo "ERROR: could not infer Cell Ranger --sample prefixes for ${carroll_sample}"
  exit 1
fi

cd "$CELLRANGER_ROOT"
if [[ -d "${CELLRANGER_ROOT}/${carroll_sample}/outs/filtered_feature_bc_matrix" ]]; then
  echo "Existing historical-matrix Cell Ranger output found for ${carroll_sample}; skipping."
  echo $(date +%T)
  exit 0
fi

"$CELLRANGER_BIN" count \
  --id="$carroll_sample" \
  --fastqs="${FASTQ_BY_SAMPLE}/${carroll_sample}" \
  --transcriptome="$TRANSCRIPTOME" \
  --sample="$cellranger_sample_prefixes" \
  --create-bam=false \
  --chemistry="$CELLRANGER_CHEMISTRY"

if [[ ! -d "${CELLRANGER_ROOT}/${carroll_sample}/outs/filtered_feature_bc_matrix" ]]; then
  echo "ERROR: Cell Ranger filtered_feature_bc_matrix was not created for ${carroll_sample}"
  exit 1
fi

echo $(date +%T)
