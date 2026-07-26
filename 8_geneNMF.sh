#!/bin/bash
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=8:00:00
#PBS -N genenmf
#PBS -koed

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline

cd ${WD}
Rscript geneNMF.R

echo $(date +%T)
