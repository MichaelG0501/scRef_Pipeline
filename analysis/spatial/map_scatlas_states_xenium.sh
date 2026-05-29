#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/map_scatlas_states_xenium.sh
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=08:00:00
#PBS -N Auto_xenium_states
#PBS -koed
#PBS -o temp/Auto_map_scATLAS_states_xenium.out
#PBS -e temp/Auto_map_scATLAS_states_xenium.err

set -eo pipefail

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
DATASET_ROOT=/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Zenodo_upload
OUTPUT_DIR=/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_xenium_states

cd $WD
mkdir -p temp

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/export_scatlas_visium_signatures.R "$OUTPUT_DIR"

conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp
python analysis/spatial/map_scatlas_states_xenium.py \
  --dataset-root "$DATASET_ROOT" \
  --output-dir "$OUTPUT_DIR" \
  --signature-dir "$OUTPUT_DIR" \
  --top-n 100 \
  --min-genes 5 \
  --threshold 0.5 \
  --threshold-quantile 0.95 \
  --min-threshold 0.5 \
  --hybrid-gap 0.3

echo $(date +%T)
