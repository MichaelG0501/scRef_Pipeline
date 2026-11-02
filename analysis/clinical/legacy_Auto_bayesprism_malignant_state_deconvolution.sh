#!/bin/bash
####################
# Historical launcher retained for provenance; not a current submission path.
# Use job_scripts/clinical/Auto_bayesprism_bulk_deconvolution.sh.
####################
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N bayesprism_malignant
#PBS -koed

####################
set -euo pipefail
####################

echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd "$WD"
Rscript analysis/clinical/legacy_Auto_bayesprism_malignant_state_deconvolution.R
echo "$(date +%T)"
