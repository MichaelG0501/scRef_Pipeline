#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=8:00:00
#PBS -N master

echo $(date +%T)

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

for sample_folder in ref_outs/by_samples/*_*_*/; do
  while [[ $(qstat | wc -l) -gt 46 ]]; do
    sleep 180
  done
  sample=$(basename "$sample_folder")
  qsub -v sample=$sample -N $sample 2_Clustering.sh
done

echo $(date +%T)
