#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/map_scatlas_states_visium.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=06:00:00
#PBS -N visium_scatlas_states
#PBS -koed
set -eo pipefail
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
DATASET_ROOT=$WD/ref_outs/Yates_spatial/source_data/Zenodo_upload
OUT=$WD/ref_outs/visium_scatlas_states
cd "$WD"
conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/export_scatlas_visium_signatures.R "$OUT"
conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp
python analysis/spatial/map_scatlas_states_visium.py --dataset-root "$DATASET_ROOT" --output-dir "$OUT"
echo $(date +%T)
