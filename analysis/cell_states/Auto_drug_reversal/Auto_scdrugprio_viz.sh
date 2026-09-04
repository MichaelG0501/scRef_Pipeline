#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_scdrugprio_viz.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=01:00:00
#PBS -N Auto_scdrugprio_viz
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD

Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_scdrugprio_visuals.R

echo $(date +%T)
