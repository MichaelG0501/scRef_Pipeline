#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/submit_geneNMF_all.sh
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

# Submit NMF jobs for all non-malignant cell types
# Usage: bash analysis/non_malignant_nmf/submit_geneNMF_all.sh

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD

# Note: R script expects "nk.cell" (not "nk") as the celltype key
for ct in macrophage fibroblast endothelial nk.cell plasma cd4 cd8; do
  qsub -v celltype=$ct -N nmf_${ct} analysis/non_malignant_nmf/run_geneNMF.sh
  echo "Submitted run_geneNMF.sh for $ct"
done
