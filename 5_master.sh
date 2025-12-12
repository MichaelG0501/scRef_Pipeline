#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=8:00:00
#PBS -N master

echo $(date +%T)

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

missing_samples=()
done_samples=()
no_ref_samples=()
no_epi_samples=()
no_cell_samples=()

for sample_folder in ref_outs/by_samples/*_*_*/; do
  while [[ $(qstat | wc -l) -gt 46 ]]; do
    sleep 180
  done

  sample=$(basename "$sample_folder")

  rds_file="ref_outs/by_samples/$sample/${sample}_epi.rds"
  no_ref="ref_outs/by_samples/$sample/no_ref"
  no_epi="ref_outs/by_samples/$sample/no_epi"
  no_cell="ref_outs/by_samples/$sample/no_cell"

  if [[ ! -f "$rds_file" && ! -f "$no_ref" && ! -f "$no_epi" && ! -f "$no_cell" ]]; then
    echo "Submitting job for $sample"
    qsub -v sample="$sample" -N "$sample" 5_InferCNA.sh
    missing_samples+=("$sample")
  else
    [[ -f "$rds_file" ]] && done_samples+=("$sample")
    [[ -f "$no_ref" ]]  && no_ref_samples+=("$sample")
    [[ -f "$no_epi" ]]  && no_epi_samples+=("$sample")
    [[ -f "$no_cell" ]] && no_cell_samples+=("$sample")
  fi
done

echo
echo "Jobs submitted (no outputs/markers): ${#missing_samples[@]}"
((${#missing_samples[@]})) && printf '  %s\n' "${missing_samples[@]}"

echo
echo "Completed (has RDS): ${#done_samples[@]}"
((${#done_samples[@]})) && printf '  %s\n' "${done_samples[@]}"

echo
echo "Skip marker no_ref: ${#no_ref_samples[@]}"
((${#no_ref_samples[@]})) && printf '  %s\n' "${no_ref_samples[@]}"

echo
echo "Skip marker no_epi: ${#no_epi_samples[@]}"
((${#no_epi_samples[@]})) && printf '  %s\n' "${no_epi_samples[@]}"

echo
echo "Skip marker no_cell: ${#no_cell_samples[@]}"
((${#no_cell_samples[@]})) && printf '  %s\n' "${no_cell_samples[@]}"

echo
echo "$(date +%T)"