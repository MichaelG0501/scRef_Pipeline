#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N Auto_scRawVal
#PBS -koed

set -euo pipefail

####################
# Validate new BAM-producing Cell Ranger matrices against the historical
# matrix_all outputs and optionally export dense CSVs using the original
# write.sh logic.
####################

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
cd "$WD"

Rscript analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.R

echo $(date +%T)
