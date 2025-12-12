#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=8:00:00
#PBS -N master

echo $(date +%T)

module purge
module load tools/dev
module load anaconda3/personal
source activate dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline

cd $WD

Rscript -e ' 

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
signatures <- c()
sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]  # match *_*_*
for (sample in sample_dirs) {
  if (!file.exists(paste0("by_samples/", sample, "/", sample, "_signatures.rds"))) {
    next
  }
  sig <- readRDS(paste0("by_samples/", sample, "/", sample, "_signatures.rds"))
  signatures <- unique(c(signatures, sig))
}
write.table(signatures, file = "cancer_signatures.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

'

missing_samples=()
no_epi_samples=()
no_cell_samples=()

for sample_folder in ref_outs/by_samples/*_*_*/; do
  while [[ $(qstat | wc -l) -gt 46 ]]; do
    sleep 180
  done

  sample=$(basename "$sample_folder")

  no_epi="ref_outs/by_samples/$sample/no_epi"
  no_cell="ref_outs/by_samples/$sample/no_cell"

  if [[ ! -f "$no_epi" && ! -f "$no_cell" ]]; then
    qsub -v sample="$sample" -N "$sample" 6_Malignancy.sh
    missing_samples+=("$sample")
  else
    [[ -f "$no_epi" ]]  && no_epi_samples+=("$sample")
    [[ -f "$no_cell" ]] && no_cell_samples+=("$sample")
  fi
done

echo
echo "Jobs submitted (with epithelial cells): ${#missing_samples[@]}"
((${#missing_samples[@]})) && printf '  %s\n' "${missing_samples[@]}"

echo
echo "Skip marker no_epi: ${#no_epi_samples[@]}"
((${#no_epi_samples[@]})) && printf '  %s\n' "${no_epi_samples[@]}"

echo
echo "Skip marker no_cell: ${#no_cell_samples[@]}"
((${#no_cell_samples[@]})) && printf '  %s\n' "${no_cell_samples[@]}"

echo
echo "$(date +%T)"
