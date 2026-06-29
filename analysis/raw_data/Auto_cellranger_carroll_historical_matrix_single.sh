#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N Auto_Car_CR_hist_one
#PBS -koed

set -euo pipefail

####################
# Single-sample launcher for the Carroll historical Cell Ranger matrix rerun.
# Pass sample=<Carroll sample> with qsub; outputs go to the same separate
# historical matrix root as Auto_cellranger_carroll_historical_matrix.sh.
####################

if [[ -z "${sample:-}" ]]; then
  echo "ERROR: submit with -v sample=<Carroll sample>"
  exit 1
fi

PBS_ARRAY_INDEX=1 sample="$sample" bash /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/raw_data/Auto_cellranger_carroll_historical_matrix.sh
