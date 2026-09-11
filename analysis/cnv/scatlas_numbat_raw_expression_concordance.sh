#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/scatlas_numbat_raw_expression_concordance.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=192gb
#PBS -l walltime=12:00:00
#PBS -N scatlas_raw_numbat_concordance
#PBS -koed
set -eo pipefail

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd ${WD}

sample_arg=${sample_arg:-all}
max_plot_cells=${max_plot_cells:-1200}
gene_bin_size=${gene_bin_size:-100}

Rscript analysis/cnv/scatlas_numbat_raw_expression_concordance.R "${sample_arg}" "${max_plot_cells}" "${gene_bin_size}"

echo $(date +%T)
