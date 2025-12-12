#!/bin/bash
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -l walltime=00:30:00
#PBS -N malignancy

echo $(date +%T)

module purge
module load tools/dev
module load anaconda3/personal
source activate dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript Malignancy.R $sample

echo $(date +%T)
