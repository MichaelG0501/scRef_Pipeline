#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=03:00:00
#PBS -N anno_nk

echo "Start: $(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG
cd $WD

Rscript Auto_anno_celltype.R nk

echo "End: $(date +%T)"
