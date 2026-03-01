#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=03:00:00
#PBS -N anno_${celltype}

echo "Start: $(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD

Rscript analysis/non_malignant_nmf/Auto_anno_celltype.R ${celltype}

echo "End: $(date +%T)"
