#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/Auto_run_scatlas_numbat_conservative_recut.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N Auto_scNBRecut
#PBS -koed

set -euo pipefail

####################
# Run conservative Numbat tree re-cut after all sample-level Numbat jobs finish.
####################

echo $(date +%T)
module purge

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_scatlas_numbat"
SIF="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_scatlas_numbat/Auto_numbat-rbase_latest.sif"
RLIB="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_scatlas_numbat/Rlib"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: missing Numbat container: $SIF"
  exit 1
fi

cd "$WD"
apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" \
  Rscript analysis/cnv/Auto_scatlas_numbat_conservative_recut.R

echo $(date +%T)
