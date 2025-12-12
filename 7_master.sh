#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=8:00:00
#PBS -N master

echo $(date +%T)

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

missing_samples=()
done_samples=()
no_cancer_samples=()
no_epi_samples=()
no_cell_samples=()
no_ref_samples=()

for sample_folder in ref_outs/by_samples/*_*_*/; do
  while [[ $(qstat | wc -l) -gt 46 ]]; do
    sleep 180
  done

  sample=$(basename "$sample_folder")

  rds_file="ref_outs/by_samples/$sample/${sample}_rank4_9_nrun10.RDS"
  no_cancer="ref_outs/by_samples/$sample/no_cancer"
  no_epi="ref_outs/by_samples/$sample/no_epi"
  no_cell="ref_outs/by_samples/$sample/no_cell"
  no_ref="ref_outs/by_samples/$sample/no_ref"

  if [[ ! -f "$rds_file" && ! -f "$no_cancer" && ! -f "$no_epi" && ! -f "$no_cell" && ! -f "$no_ref" ]]; then
    qsub -v sample="$sample" -N "$sample" 7_NMF.sh
    missing_samples+=("$sample")
  else
    [[ -f "$rds_file" ]] && done_samples+=("$sample")
    [[ -f "$no_cancer" ]] && no_cancer_samples+=("$sample")
    [[ -f "$no_epi" ]]  && no_epi_samples+=("$sample")
    [[ -f "$no_cell" ]] && no_cell_samples+=("$sample")
    [[ -f "$no_ref" ]]  && no_ref_samples+=("$sample")
  fi
done

echo
echo "Jobs submitted (with epithelial cells): ${#missing_samples[@]}"
((${#missing_samples[@]})) && printf '  %s\n' "${missing_samples[@]}"

echo
echo "Completed (has RDS): ${#done_samples[@]}"
((${#done_samples[@]})) && printf '  %s\n' "${done_samples[@]}"

echo
echo "Skip marker no_cancer: ${#no_cancer_samples[@]}"
((${#no_cancer_samples[@]})) && printf '  %s\n' "${no_cancer_samples[@]}"

echo
echo "Skip marker no_epi: ${#no_epi_samples[@]}"
((${#no_epi_samples[@]})) && printf '  %s\n' "${no_epi_samples[@]}"

echo
echo "Skip marker no_cell: ${#no_cell_samples[@]}"
((${#no_cell_samples[@]})) && printf '  %s\n' "${no_cell_samples[@]}"

echo
echo "Skip marker no_ref: ${#no_ref_samples[@]}"
((${#no_ref_samples[@]})) && printf '  %s\n' "${no_ref_samples[@]}"

echo
echo "$(date +%T)"
