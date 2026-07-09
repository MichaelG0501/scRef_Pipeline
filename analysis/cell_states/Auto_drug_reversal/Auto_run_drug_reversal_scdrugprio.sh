#!/bin/bash
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=72:00:00
#PBS -N Auto_scDrugPrio
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
export AUTO_SCDRUGPRIO_PPI=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/ppi.txt
export AUTO_SCDRUGPRIO_DRUG_TARGETS=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/all_drug_targets_drug_bank.txt
export AUTO_SCDRUGPRIO_PHARMA_EFFECT=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/all_drug_targets_drug_bank.txt
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_scdrugprio.R
echo $(date +%T)
