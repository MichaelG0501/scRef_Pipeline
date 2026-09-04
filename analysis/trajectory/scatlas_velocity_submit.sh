#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/trajectory/scatlas_velocity_submit.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N scAtlas_VelSub
#PBS -koed

set -euo pipefail

####################
# Convenience submitter for scATLAS RNA velocity. It exports metadata, prepares
# velocyto references, submits per-sample filter/sort and velocyto jobs, then
# submits scVelo visualisation and R nodeplots with afterok dependencies.
####################

echo $(date +%T)
module purge
module load tools/dev

WD="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
OUT_EPH="${WD}/ref_outs/Auto_velocity_scATLAS"
OUT_LIVE="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_velocity_scATLAS"
mkdir -p "${OUT_LIVE}/logs" "${OUT_LIVE}/looms" "${OUT_LIVE}/h5ad"
mkdir -p "${OUT_EPH}/looms" "${OUT_EPH}/h5ad" "${OUT_EPH}/coord" "${OUT_EPH}/ref" "${OUT_EPH}/tmp_pycache" "${OUT_EPH}/logs"
cd "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
manifest="${OUT_LIVE}/tables/Auto_scatlas_velocity_sample_manifest.csv"
metadata="${OUT_LIVE}/tables/Auto_scatlas_velocity_cell_metadata.csv"
if [[ "${SCATLAS_FORCE_VELOCITY_METADATA:-FALSE}" == "TRUE" || ! -s "$manifest" || ! -s "$metadata" ]]; then
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
  Rscript analysis/trajectory/scatlas_velocity_metadata.R
else
  echo "Reusing existing velocity metadata: $manifest"
fi

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity
export PYTHONPYCACHEPREFIX="${OUT_EPH}/tmp_pycache"
gene_gtf="${OUT_EPH}/ref/genes.GRCh38-2024-A.gtf"
mask_gtf="${OUT_EPH}/ref/repeatmasker.hg38.gtf"
if [[ "${SCATLAS_FORCE_VELOCITY_REFS:-FALSE}" == "TRUE" || ! -s "$gene_gtf" || ! -s "$mask_gtf" ]]; then
  python analysis/trajectory/scatlas_velocity_prepare_refs.py
else
  echo "Reusing existing velocity references under ${OUT_EPH}/ref"
fi

if [[ ! -f "$manifest" ]]; then
  echo "ERROR: missing manifest after metadata export: $manifest"
  exit 1
fi

qsub_bin="/opt/pbs/bin/qsub"
qstat_bin="/opt/pbs/bin/qstat"

throttle() {
  while [[ $("$qstat_bin" | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
}

sanitize_job_name() {
  local x="$1"
  x="${x//[^A-Za-z0-9_]/_}"
  echo "${x:0:14}"
}

velocyto_jobs=()

while IFS=, read -r sample dataset study bam bai barcodes_file n_cells has_bam has_bai; do
  [[ "$sample" == "sample" ]] && continue
  if [[ "$has_bam" != "TRUE" && "$has_bam" != "true" ]]; then
    echo "Skipping ${sample}: manifest has_bam=${has_bam}"
    continue
  fi
  short_name=$(sanitize_job_name "$sample")

  throttle
  jid_sort=$("$qsub_bin" \
    -v sample="$sample" \
    -N "Srt_${short_name}" \
    -o "${OUT_LIVE}/logs/Auto_filter_sort_${sample}.log" \
    -e "${OUT_LIVE}/logs/Auto_filter_sort_${sample}.err" \
    analysis/trajectory/scatlas_velocity_filter_sort.sh)
  echo "Submitted filter/sort ${sample}: ${jid_sort}"

  throttle
  jid_vel=$("$qsub_bin" \
    -W depend=afterok:${jid_sort} \
    -v sample="$sample" \
    -N "Vel_${short_name}" \
    -o "${OUT_LIVE}/logs/Auto_velocyto_${sample}.log" \
    -e "${OUT_LIVE}/logs/Auto_velocyto_${sample}.err" \
    analysis/trajectory/scatlas_velocity_run_velocyto.sh)
  echo "Submitted velocyto ${sample}: ${jid_vel}"
  velocyto_jobs+=("$jid_vel")
done < "$manifest"

if [[ ${#velocyto_jobs[@]} -eq 0 ]]; then
  echo "ERROR: no velocyto jobs were submitted."
  exit 1
fi

dep=$(IFS=:; echo "${velocyto_jobs[*]}")
throttle
jid_vis=$("$qsub_bin" \
  -W depend=afterok:${dep} \
  -o "${OUT_LIVE}/logs/Auto_scvelo_visualisation.log" \
  -e "${OUT_LIVE}/logs/Auto_scvelo_visualisation.err" \
  analysis/trajectory/scatlas_velocity_run_scvelo_visualisation.sh)
echo "Submitted dependent scVelo visualisation: ${jid_vis}"

throttle
jid_nodes=$("$qsub_bin" \
  -W depend=afterok:${jid_vis} \
  -o "${OUT_LIVE}/logs/Auto_velocity_nodeplots.log" \
  -e "${OUT_LIVE}/logs/Auto_velocity_nodeplots.err" \
  analysis/trajectory/scatlas_velocity_run_nodeplots.sh)
echo "Submitted dependent R nodeplots: ${jid_nodes}"

echo $(date +%T)
