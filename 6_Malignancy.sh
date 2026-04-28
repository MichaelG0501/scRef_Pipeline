#!/bin/bash
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=01:00:00
#PBS -N Auto_malignancy
#PBS -koed

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD

Rscript Malignancy.R "$sample"

echo $(date +%T)
