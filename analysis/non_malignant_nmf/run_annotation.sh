#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/run_annotation.sh
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=03:00:00
#PBS -N anno_${celltype}
#PBS -koed

echo "Start: $(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD

Rscript analysis/non_malignant_nmf/nmf_celltype_annotation.R ${celltype}

echo "End: $(date +%T)"
