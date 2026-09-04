#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/trajectory/scatlas_velocity_run_velocyto.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=36:00:00
#PBS -N scAtlas_Velocyto
#PBS -koed

set -euo pipefail

####################
# Run velocyto on the filtered, coordinate-sorted scATLAS sample BAM.
####################

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OUT="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_velocity_scATLAS"
OUT_LIVE="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_velocity_scATLAS"
export PYTHONPYCACHEPREFIX="${OUT}/tmp_pycache"

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>"
  exit 1
fi

BC="${OUT_LIVE}/barcodes/${sample}_qc_barcodes.tsv"
BAM="${OUT}/coord/${sample}.qc.coord.bam"
LOOM_DIR="${OUT}/looms/${sample}"
GENE_GTF="${OUT}/ref/genes.GRCh38-2024-A.gtf"
MASK_GTF="${OUT}/ref/repeatmasker.hg38.gtf"

if [[ ! -s "$BC" ]]; then
  echo "ERROR: missing barcode file: $BC"
  exit 1
fi
if [[ ! -f "$BAM" ]]; then
  echo "ERROR: missing filtered BAM: $BAM"
  exit 1
fi
if [[ ! -s "$GENE_GTF" || ! -s "$MASK_GTF" ]]; then
  echo "ERROR: missing velocity reference files under ${OUT}/ref"
  exit 1
fi
LOOM_DIR_LIVE="${OUT_LIVE}/looms/${sample}"
if find "$LOOM_DIR" -maxdepth 1 -type f -name "*.loom" -size +0c 2>/dev/null | grep -q .; then
  echo "Existing loom found for $sample in ephemeral; skipping velocyto."
  mkdir -p "$LOOM_DIR_LIVE"
  cp -n "$LOOM_DIR"/*.loom "$LOOM_DIR_LIVE/" 2>/dev/null || true
  echo $(date +%T)
  exit 0
fi
if find "$LOOM_DIR_LIVE" -maxdepth 1 -type f -name "*.loom" -size +0c 2>/dev/null | grep -q .; then
  echo "Existing loom found for $sample in live; skipping velocyto."
  echo $(date +%T)
  exit 0
fi

mkdir -p "$LOOM_DIR" "${OUT}/logs"
cd "$WD"

python analysis/trajectory/scatlas_velocyto_run.py \
  --bcfile "$BC" \
  --outputfolder "$LOOM_DIR" \
  --sampleid "$sample" \
  --mask "$MASK_GTF" \
  --samtools-threads 8 \
  --samtools-memory 7000 \
  -vv \
  "$BAM" \
  "$GENE_GTF"

####################
# Copy loom to live/ for persistent storage (critical downstream input)
####################
mkdir -p "$LOOM_DIR_LIVE"
cp "$LOOM_DIR"/*.loom "$LOOM_DIR_LIVE/"
echo "Copied loom to $LOOM_DIR_LIVE"

echo $(date +%T)
