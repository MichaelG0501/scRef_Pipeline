#!/bin/bash
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -l walltime=02:00:00
#PBS -N scAtlas_VelNodes
#PBS -koed

set -euo pipefail

####################
# Draw R nodeplots from scVelo transition CSVs.
####################

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
cd "$WD"

Rscript analysis/trajectory/scatlas_velocity_nodeplots.R

echo $(date +%T)
