#!/bin/bash
set -euo pipefail

echo "$(date +%T) Starting OCCAMS featureCounts batch submission"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
cd "$WD"
python "$WD/analysis/OCCAMS/occams_prepare_featurecount_batches.py"
while [[ $(/opt/pbs/bin/qstat | grep -c sg3723 || true) -ge 46 ]]; do
    sleep 55
done
/opt/pbs/bin/qsub "$WD/analysis/OCCAMS/occams_featurecounts_grch37_batch.sh"
echo "$(date +%T) OCCAMS featureCounts batch submission complete"
