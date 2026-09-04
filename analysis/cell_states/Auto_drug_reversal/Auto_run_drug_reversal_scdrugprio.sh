#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_scdrugprio.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=72:00:00
#PBS -N Auto_scDrugPrio
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
ENV_PREFIX=$WD/ref_outs/Auto_drug_reversal/conda/Auto_drug_reversal
source activate "$ENV_PREFIX"
cd $WD
export AUTO_SCDRUGPRIO_PPI=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/ppi.txt
export AUTO_SCDRUGPRIO_DRUG_TARGETS=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/all_drug_targets_drug_bank.txt
export AUTO_SCDRUGPRIO_PHARMA_EFFECT=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/all_drug_targets_drug_bank.txt
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_scdrugprio.R
echo $(date +%T)
