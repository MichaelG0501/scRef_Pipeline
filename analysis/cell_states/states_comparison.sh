#!/bin/bash
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=4:00:00
#PBS -N Auto_states_compare

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd ${WD}

Rscript analysis/cell_states/Auto_states_comparison.R

echo $(date +%T)
