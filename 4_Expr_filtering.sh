#!/bin/bash
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=24:00:00
#PBS -N expr_filtering

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript Expr_filtering.R

echo $(date +%T)
