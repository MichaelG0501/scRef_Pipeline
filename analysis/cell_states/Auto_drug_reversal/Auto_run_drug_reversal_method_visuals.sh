#!/bin/bash
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -l walltime=02:00:00
#PBS -N Auto_Reversal_Viz
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_method_visuals.R
echo $(date +%T)
