#!/bin/bash
# Master script to submit annotation jobs for all cell types
# Only submits if MP_outs_default.rds exists for that cell type

WD=/rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG
BASE=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs

declare -A FOLDER_MAP
FOLDER_MAP[macrophage]="nmf_macrophage"
FOLDER_MAP[fibroblast]="nmf_fibroblast"
FOLDER_MAP[endothelial]="nmf_endothelial"
FOLDER_MAP[plasma]="nmf_plasma"
FOLDER_MAP[cd4]="nmf_cd4"
FOLDER_MAP[cd8]="nmf_cd8"
FOLDER_MAP[nk]="nmf_nk"

cd $WD

for ct in macrophage fibroblast endothelial plasma cd4 cd8 nk; do
  rds_file="${BASE}/${FOLDER_MAP[$ct]}/MP_outs_default.rds"
  if [[ -f "$rds_file" ]]; then
    qsub Auto_anno_${ct}.sh
    echo "Submitted Auto_anno_${ct}.sh"
  else
    echo "SKIPPED ${ct}: ${rds_file} not found"
  fi
done
