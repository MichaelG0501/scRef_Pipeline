#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/non_malignant_nmf/submit_annotation_all.sh
#   Methodology: analysis/methodology/non_malignant_nmf/non_malignant_nmf_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

# Submit annotation jobs for non-malignant cell types (only if NMF output exists)
# Usage: bash analysis/non_malignant_nmf/submit_annotation_all.sh

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
BASE=${WD}/ref_outs
cd $WD

declare -A FOLDER_MAP
FOLDER_MAP[macrophage]="nmf_macrophage"
FOLDER_MAP[fibroblast]="nmf_fibroblast"
FOLDER_MAP[endothelial]="nmf_endothelial"
FOLDER_MAP[plasma]="nmf_plasma"
FOLDER_MAP[cd4]="nmf_cd4"
FOLDER_MAP[cd8]="nmf_cd8"
FOLDER_MAP[nk]="nmf_nk"

for ct in macrophage fibroblast endothelial plasma cd4 cd8 nk; do
  rds_file="${BASE}/${FOLDER_MAP[$ct]}/MP_outs_default.rds"
  if [[ -f "$rds_file" ]]; then
    qsub -v celltype=$ct -N anno_${ct} analysis/non_malignant_nmf/run_annotation.sh
    echo "Submitted run_annotation.sh for $ct"
  else
    echo "SKIPPED ${ct}: ${rds_file} not found"
  fi
done
