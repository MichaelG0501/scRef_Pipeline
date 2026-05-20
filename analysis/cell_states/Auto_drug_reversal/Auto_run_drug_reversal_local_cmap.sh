#!/bin/bash
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=04:00:00
#PBS -N Auto_Drug_CMap
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
ENV_PREFIX=$WD/ref_outs/Auto_drug_reversal/conda/Auto_drug_reversal
if [[ -d "$ENV_PREFIX" ]]; then
  source activate "$ENV_PREFIX"
else
  source activate /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_drug_reversal/conda/Auto_drug_reversal
fi
cd $WD
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_local_cmap.R
echo $(date +%T)
