#!/bin/bash
#PBS -l select=1:ncpus=4:mem=128gb
#PBS -l walltime=8:00:00
#PBS -N annotation
#PBS -koed

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline

cd $WD

Rscript Annotation.R

echo $(date +%T)
