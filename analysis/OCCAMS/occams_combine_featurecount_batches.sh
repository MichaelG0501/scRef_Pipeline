#!/bin/bash
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=04:00:00
#PBS -N occams_fc_combine
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_combine_featurecount_batches.log

set -euo pipefail
echo "$(date +%T) Starting OCCAMS featureCounts batch combination"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
cd "$WD"
python "$WD/analysis/OCCAMS/occams_combine_featurecount_batches.py"
SCREF_FORCE_REBUILD=TRUE python "$WD/analysis/OCCAMS/occams_build_count_matrix.py"
echo "$(date +%T) OCCAMS featureCounts batch combination complete"
