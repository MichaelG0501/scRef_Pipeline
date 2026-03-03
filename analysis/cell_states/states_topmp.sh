#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=2:00:00
#PBS -N Auto_states_topmp

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd ${WD}

Rscript analysis/cell_states/Auto_states_topmp.R

echo $(date +%T)
