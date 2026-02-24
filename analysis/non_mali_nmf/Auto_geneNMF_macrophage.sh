#!/bin/bash
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=06:00:00
#PBS -N nmf_macro

echo "Start: $(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf

WD=/rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG
cd $WD

Rscript Auto_geneNMF_celltype.R macrophage

echo "End: $(date +%T)"
