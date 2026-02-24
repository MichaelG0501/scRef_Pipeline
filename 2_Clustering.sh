#!/bin/bash
#PBS -l select=1:ncpus=2:mem=45gb
#PBS -l walltime=1:00:00
#PBS -N clustering

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript Clustering.R $sample

echo $(date +%T)
