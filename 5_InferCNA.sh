#!/bin/bash
#PBS -l select=1:ncpus=2:mem=25gb
#PBS -l walltime=01:00:00
#PBS -N infercna

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript InferCNA.R $sample

echo $(date +%T)
