#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=02:00:00
#PBS -N Auto_xenium_replot
#PBS -koed
#PBS -o temp/Auto_replot_xenium_states_highres.out
#PBS -e temp/Auto_replot_xenium_states_highres.err

set -eo pipefail

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD
mkdir -p temp ref_outs/Auto_xenium_state_replot_highres

python analysis/spatial/Auto_replot_xenium_states_highres.py

echo $(date +%T)
