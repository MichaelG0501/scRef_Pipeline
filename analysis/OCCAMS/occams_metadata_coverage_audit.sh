#!/bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=01:00:00
#PBS -N occams_meta_audit
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_metadata_coverage_audit.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS metadata coverage audit"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
cd "$WD"
python "$WD/analysis/OCCAMS/occams_metadata_coverage_audit.py"
echo "$(date +%T) OCCAMS metadata coverage audit complete"
