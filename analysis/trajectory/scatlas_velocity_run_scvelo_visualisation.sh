#!/bin/bash
#PBS -l select=1:ncpus=8:mem=160gb
#PBS -l walltime=24:00:00
#PBS -N scAtlas_scVelo
#PBS -koed

set -euo pipefail

####################
# Run scVelo visualisation and state-direction summary tables.
####################

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity

WD="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_velocity_scATLAS"
export PYTHONPYCACHEPREFIX="${OUT}/tmp_pycache"
cd "$WD"

python analysis/trajectory/scatlas_velocity_scvelo_visualise.py

echo $(date +%T)
