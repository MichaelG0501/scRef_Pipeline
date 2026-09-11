#!/bin/bash
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_manual_annotation_threshold_calibration.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=01:00:00
#PBS -N visiumhd_manual_calibration
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd $WD
python analysis/spatial/visiumhd_manual_annotation_threshold_calibration.py \
    --annotation-dir "$WD/ref_outs/visium_hd_outs/tables" \
    --output-dir "$WD/ref_outs/visium_hd_outs/tables" \
    --sample-names SUR1231 FFPEA1 FFPED1
echo $(date +%T)
