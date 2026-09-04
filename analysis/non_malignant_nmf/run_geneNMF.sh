#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/run_geneNMF.sh
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

#PBS -l select=1:ncpus=8:mem=72gb
#PBS -l walltime=08:00:00
#PBS -N nmf_${celltype}
#PBS -koed
set -eo pipefail

echo "Start: $(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD

Rscript analysis/non_malignant_nmf/nmf_celltype_geneNMF.R ${celltype}

echo "End: $(date +%T)"
