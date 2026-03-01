#!/bin/bash
####################
# Moved from: analysis/Auto_MP_correlation.sh
# Reorganized as part of analysis/ restructuring
####################
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=4:00:00
#PBS -N MP_corr

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd ${WD}

Rscript analysis/metaprograms/mp_database_correlation.R

echo $(date +%T)
