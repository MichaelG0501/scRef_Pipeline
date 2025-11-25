#!/bin/bash
#PBS -l select=1:ncpus=2:mem=45gb
#PBS -l walltime=1:00:00
#PBS -N clustering

echo $(date +%T)

module purge
module load tools/dev
module load anaconda3/personal
source activate dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript Clustering.R $sample

echo $(date +%T)
