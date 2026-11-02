#!/bin/bash
####################
# Historical launcher retained for provenance; not a current submission path.
# Use job_scripts/clinical/Auto_bayesprism_bulk_deconvolution.sh.
####################
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=00:30:00
#PBS -N occams_bp_qc
#PBS -koed
set -euo pipefail
echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
Rscript analysis/clinical/legacy_Auto_occams_bayesprism_qc_plots.R
echo "$(date +%T)"
