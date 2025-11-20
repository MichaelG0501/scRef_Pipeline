#!/bin/bash
#PBS -l select=1:ncpus=8:mem=921gb
#PBS -l walltime=8:00:00
#PBS -N qc_pipeline

echo $(date +%T)

module purge
module load tools/dev
module load anaconda3/personal
source activate dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd ${WD}
Rscript QC_Pipeline.R

echo $(date +%T)
