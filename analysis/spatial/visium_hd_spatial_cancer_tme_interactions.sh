#!/bin/bash
####################
# Retained in place because file moves and deletion are prohibited.
# Current launcher: job_scripts/spatial/visium_hd_spatial_cancer_tme_interactions.sh
# This historical launcher is intentionally disabled.
printf '%s\n' 'Use job_scripts/spatial/visium_hd_spatial_cancer_tme_interactions.sh' >&2
exit 2
####################
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/visium_hd_spatial_cancer_tme_interactions.sh
#   Methodology: not required (PBS wrapper; the invoked R workflow documents the method)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Run the current Visium HD cancer-TME spatial validation workflow.
####################
#PBS -l select=1:ncpus=8:mem=160gb
#PBS -l walltime=24:00:00
#PBS -N vh_cancer_tme
#PBS -koed
#PBS -o temp/visium_hd_spatial_cancer_tme_interactions.out
#PBS -e temp/visium_hd_spatial_cancer_tme_interactions.err

set -euo pipefail

echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd "$WD"
mkdir -p temp

export SCREF_NCORES=8
export SCREF_FORCE_REBUILD="${SCREF_FORCE_REBUILD:-TRUE}"
export SCREF_SPATIAL_PERMUTATIONS="${SCREF_SPATIAL_PERMUTATIONS:-4999}"
Rscript analysis/spatial/visium_hd_spatial_cancer_tme_interactions.R

echo "$(date +%T)"
