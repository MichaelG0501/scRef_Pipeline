#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_asgard.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N Auto_ASGARD
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
ENV_PREFIX=$WD/ref_outs/Auto_drug_reversal/conda/Auto_drug_reversal
source activate "$ENV_PREFIX"
cd $WD
export AUTO_ASGARD_DRUG_RESPONSE=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference/stomach_rankMatrix.txt
export AUTO_ASGARD_GENE_INFO=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference/stomach_gene_info.txt
export AUTO_ASGARD_DRUG_INFO=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference/stomach_drug_info.txt
if [[ -f ref_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.sh ]]; then
  source ref_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.sh
fi
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_asgard.R
echo $(date +%T)
