#!/bin/bash
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_run_visium_hd_summary.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=00:30:00
#PBS -N visium_summary
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD

Rscript analysis/spatial/visium_hd_summary.R
echo $(date +%T)
