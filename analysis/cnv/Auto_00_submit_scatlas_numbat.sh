#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/Auto_00_submit_scatlas_numbat.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
set -euo pipefail

####################
# Submit the full scATLAS Numbat workflow: export inputs, prepare container,
# per-sample SNP pileup, per-sample Numbat, then conservative tree re-cut.
####################

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_scatlas_numbat"
cd "$WD"

mkdir -p "${OUT}/logs"
run_tag=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${OUT}/logs/Auto_scatlas_numbat_run_${run_tag}"
mkdir -p "$LOG_DIR"

echo $(date +%T)
module purge
module load tools/dev

manifest="${OUT}/Auto_scatlas_numbat_manifest.csv"
raw_validation_prefix="${SCATLAS_VALIDATION_PREFIX:-Auto_scatlas_cellranger_matrix_validation}"
raw_validation_csv="/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/validation/${raw_validation_prefix}.csv"
raw_validation_fail="/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/validation/${raw_validation_prefix}_failures.csv"

if [[ ! -f "$raw_validation_csv" ]]; then
  echo "ERROR: missing raw Cell Ranger validation CSV: $raw_validation_csv"
  echo "Run qsub analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.sh after Cell Ranger completes."
  exit 1
fi
if awk -F, 'NR > 1 && $3 != "ok" {bad=1} END {exit bad ? 0 : 1}' "$raw_validation_csv"; then
  echo "ERROR: raw Cell Ranger validation contains non-ok rows: $raw_validation_csv"
  if [[ -s "$raw_validation_fail" ]]; then
    echo "Failure details: $raw_validation_fail"
  fi
  exit 1
fi

if [[ ! -f "$manifest" ]]; then
  eval "$(~/miniforge3/bin/conda shell.bash hook)"
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
  Rscript analysis/cnv/Auto_scatlas_numbat_export_inputs.R
fi

if [[ ! -f "$manifest" ]]; then
  echo "ERROR: missing manifest after export: $manifest"
  exit 1
fi

throttle() {
  while [[ $(/opt/pbs/bin/qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
}

sanitize_job_name() {
  local x="$1"
  x="${x//[^A-Za-z0-9_]/_}"
  echo "${x:0:14}"
}

throttle
jid_img=$(/opt/pbs/bin/qsub \
  -o "${LOG_DIR}/Auto_prepare_scatlas_numbat_container.log" \
  -e "${LOG_DIR}/Auto_prepare_scatlas_numbat_container.err" \
  analysis/cnv/Auto_prepare_scatlas_numbat_container.sh)
echo "Submitted Numbat container preparation: ${jid_img}"

run_jobs=()
while IFS=, read -r sample dataset raw_sample bam barcodes_file count_rds metadata_csv sample_out_dir numbat_dir allele_file clone_post_file joint_post_file n_cells has_bam has_barcode_file; do
  [[ "$sample" == "sample" ]] && continue
  short_name=$(sanitize_job_name "$sample")

  throttle
  jid_pile=$(/opt/pbs/bin/qsub \
    -W depend=afterok:${jid_img} \
    -v sample="$sample" \
    -N "scNBp_${short_name}" \
    -o "${LOG_DIR}/Auto_scatlas_numbat_pileup_${sample}.log" \
    -e "${LOG_DIR}/Auto_scatlas_numbat_pileup_${sample}.err" \
    analysis/cnv/Auto_run_scatlas_numbat_pileup.sh)
  echo "Submitted Numbat pileup ${sample}: ${jid_pile}"

  throttle
  jid_run=$(/opt/pbs/bin/qsub \
    -W depend=afterok:${jid_pile} \
    -v sample="$sample" \
    -N "scNBr_${short_name}" \
    -o "${LOG_DIR}/Auto_scatlas_numbat_run_${sample}.log" \
    -e "${LOG_DIR}/Auto_scatlas_numbat_run_${sample}.err" \
    analysis/cnv/Auto_run_scatlas_numbat_sample.sh)
  echo "Submitted Numbat run ${sample}: ${jid_run}"
  run_jobs+=("$jid_run")
done < "$manifest"

if [[ ${#run_jobs[@]} -eq 0 ]]; then
  echo "ERROR: manifest did not contain samples."
  exit 1
fi

dep=$(IFS=:; echo "${run_jobs[*]}")
throttle
jid_recut=$(/opt/pbs/bin/qsub \
  -W depend=afterok:${dep} \
  -o "${LOG_DIR}/Auto_scatlas_numbat_conservative_recut.log" \
  -e "${LOG_DIR}/Auto_scatlas_numbat_conservative_recut.err" \
  analysis/cnv/Auto_run_scatlas_numbat_conservative_recut.sh)
echo "Submitted dependent conservative re-cut: ${jid_recut}"
echo $(date +%T)
