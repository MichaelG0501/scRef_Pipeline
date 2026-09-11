#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/Auto_run_scatlas_numbat_sample.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=12:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N Auto_scNBRun
#PBS -koed

set -euo pipefail

####################
# Run Numbat clone inference for one scATLAS sample inside the prepared
# container.
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

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: missing Numbat container: $SIF"
  exit 1
fi

cd "$WD"
apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" \
  Rscript analysis/cnv/Auto_scatlas_numbat_run_sample.R "$sample" "$MANIFEST" "$NCORES"

echo $(date +%T)
