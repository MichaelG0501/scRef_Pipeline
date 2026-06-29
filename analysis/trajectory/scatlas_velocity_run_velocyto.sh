#!/bin/bash
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

WD="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_velocity_scATLAS"
export PYTHONPYCACHEPREFIX="${OUT}/tmp_pycache"

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>"
  exit 1
fi

BC="${OUT}/barcodes/${sample}_qc_barcodes.tsv"
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
if find "$LOOM_DIR" -maxdepth 1 -type f -name "*.loom" -size +0c 2>/dev/null | grep -q .; then
  echo "Existing loom found for $sample; skipping velocyto."
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

echo $(date +%T)
