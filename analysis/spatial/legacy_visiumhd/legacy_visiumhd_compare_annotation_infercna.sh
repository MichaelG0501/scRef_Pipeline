#!/bin/bash
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_compare_annotation_infercna.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -l walltime=01:00:00
#PBS -N visiumhd_compare_infercna
#PBS -koed
echo $(date +%T)
set -e
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD
Rscript analysis/spatial/visiumhd_compare_annotation_infercna.R
echo $(date +%T)
