#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=06:00:00
#PBS -N visium_scatlas_states
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
OUT=$WD/ref_outs/visium_scatlas_states
cd $WD
Rscript analysis/spatial/export_scatlas_visium_signatures.R $OUT
python analysis/spatial/map_scatlas_states_visium.py --output-dir $OUT
echo $(date +%T)
