#!/bin/bash
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=04:00:00
#PBS -N Auto_xenium_cmp
#PBS -koed
#PBS -o temp/Auto_compare_xenium_state_methods.out
#PBS -e temp/Auto_compare_xenium_state_methods.err

set -eo pipefail

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD
mkdir -p temp ref_outs/Auto_xenium_state_method_comparison

python analysis/spatial/Auto_compare_xenium_state_methods.py

echo $(date +%T)
