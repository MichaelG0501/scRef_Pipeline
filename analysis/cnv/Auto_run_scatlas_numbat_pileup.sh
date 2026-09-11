#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/Auto_run_scatlas_numbat_pileup.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=12:mem=96gb
#PBS -l walltime=36:00:00
#PBS -N Auto_scNBPile
#PBS -koed

set -euo pipefail

####################
# Run Numbat pileup_and_phase.R for one scATLAS sample.
####################

echo $(date +%T)
module purge

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>"
  exit 1
fi

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_scatlas_numbat"
MANIFEST="${OUT}/Auto_scatlas_numbat_manifest.csv"
SIF="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_scatlas_numbat/Auto_numbat-rbase_latest.sif"
NCORES="${NCORES:-12}"
RLIB="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_scatlas_numbat/Rlib"

mkdir -p "${OUT}/logs" "${OUT}/singularity_cache" "${OUT}/tmp"
export SINGULARITY_CACHEDIR="${OUT}/singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export TMPDIR="${OUT}/tmp"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: missing manifest: $MANIFEST"
  exit 1
fi
if [[ ! -f "$SIF" ]]; then
  echo "ERROR: missing Numbat container: $SIF"
  exit 1
fi

line=$(awk -F, -v s="$sample" 'NR > 1 && $1 == s {print; exit}' "$MANIFEST")
if [[ -z "$line" ]]; then
  echo "ERROR: sample not found in manifest: $sample"
  exit 1
fi

BAM=$(awk -F, '{print $4}' <<< "$line")
BC=$(awk -F, '{print $5}' <<< "$line")
SAMPLE_OUT=$(awk -F, '{print $8}' <<< "$line")
ALLELE=$(awk -F, '{print $10}' <<< "$line")

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: missing BAM: $BAM"
  exit 1
fi
if [[ ! -s "$BC" ]]; then
  echo "ERROR: missing barcode file: $BC"
  exit 1
fi
if [[ -s "$ALLELE" ]]; then
  echo "Existing allele counts found for $sample; skipping pileup/phasing."
  echo $(date +%T)
  exit 0
fi

mkdir -p "$SAMPLE_OUT"

apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" \
  Rscript /numbat/inst/bin/pileup_and_phase.R \
    --label "$sample" \
    --samples "$sample" \
    --bams "$BAM" \
    --barcodes "$BC" \
    --outdir "$SAMPLE_OUT" \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores "$NCORES" \
    --UMItag Auto \
    --cellTAG CB

if [[ ! -s "$ALLELE" ]]; then
  echo "ERROR: allele output was not created: $ALLELE"
  exit 1
fi

echo $(date +%T)
