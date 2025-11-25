#!/bin/bash
#PBS -l select=1:ncpus=4:mem=128gb
#PBS -l walltime=8:00:00
#PBS -N annotation

echo $(date +%T)

module purge
module load tools/dev
module load anaconda3/personal
source activate dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript Annotation.R

echo $(date +%T)
