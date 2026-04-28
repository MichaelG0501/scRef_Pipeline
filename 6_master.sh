#!/bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -N Auto_malignancy_master
#PBS -koed

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd $WD

Rscript -e '
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
signatures <- c()
sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]
for (sample in sample_dirs) {
  sig_path <- paste0("by_samples/", sample, "/", sample, "_signatures.rds")
  if (!file.exists(sig_path)) next
  sig <- readRDS(sig_path)
  signatures <- unique(c(signatures, sig))
}
write.table(signatures, file = "cancer_signatures.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
'

submitted=0
done_count=0
skipped=0

for sample_folder in ref_outs/by_samples/*_*_*/; do
  sample=$(basename "$sample_folder")
  epi_file="ref_outs/by_samples/$sample/${sample}_epi.rds"
  epi_f_file="ref_outs/by_samples/$sample/${sample}_epi_f.rds"
  no_ref="ref_outs/by_samples/$sample/no_ref"
  no_epi="ref_outs/by_samples/$sample/no_epi"
  no_cell="ref_outs/by_samples/$sample/no_cell"

  if [[ -f "$epi_f_file" ]]; then
    done_count=$((done_count + 1))
    continue
  fi

  if [[ ! -f "$epi_file" || -f "$no_ref" || -f "$no_epi" || -f "$no_cell" ]]; then
    skipped=$((skipped + 1))
    continue
  fi

  while [[ $(/opt/pbs/bin/qstat -u sg3723 | awk 'NR > 5 {print $1}' | wc -l) -gt 46 ]]; do
    sleep 180
  done

  echo "Submitting malignancy job for $sample"
  /opt/pbs/bin/qsub -koed -v sample="$sample" -N "$sample" Auto_6_Malignancy.sh
  submitted=$((submitted + 1))
done

echo "Submitted malignancy jobs: $submitted"
echo "Already complete: $done_count"
echo "Skipped: $skipped"
echo $(date +%T)
