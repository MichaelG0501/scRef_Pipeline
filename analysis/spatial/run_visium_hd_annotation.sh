#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/run_visium_hd_annotation.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
####################
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=36:00:00
#PBS -N visiumhd_final_annotation
#PBS -koed

# Analysis registry:
#   Status: active
#   Description: Run the final two-method Visium HD annotation workflow.
#   Methodology: analysis/methodology/spatial/visium_hd_final_annotation_methodology.md
#   Inputs: analysis/spatial/visium_hd_samples.tsv and paths listed within it.
#   Outputs: ref_outs/visium_hd_outs/{rctd,intermediate,tables,figures,logs}/.
#   Run: /opt/pbs/bin/qsub analysis/spatial/run_visium_hd_annotation.sh
####################

####################
set -euo pipefail

echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline

OUT=$WD/ref_outs/visium_hd_outs
RCTD_OUT=$OUT/rctd
cd "$WD"

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visium_hd_rctd_doublet_detection.R \
    --manifest "$MANIFEST" \
    --output-dir "$RCTD_OUT" \
    --reference-path "$WD/ref_outs/EAC_Ref_merged.rds" \
    --legacy-cache-dir "$OUT/legacy_visiumhd/rctd" \
    --min-umis 100 \
    --force-rebuild "${SCREF_FORCE_RCTD:-FALSE}" \
    --max-cores 8

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
annotation_args=()
if [[ "${SCREF_REUSE_COMPLETE:-FALSE}" == "TRUE" ]]; then
    annotation_args+=(--reuse-complete)
fi
if [[ -n "${SCREF_FORCE_SAMPLES:-}" ]]; then
    read -r -a force_samples <<< "${SCREF_FORCE_SAMPLES//,/ }"
    annotation_args+=(--force-samples "${force_samples[@]}")
fi
python analysis/spatial/visium_hd_celltype_annotation.py \
    --manifest "$MANIFEST" \
    --output-dir "$OUT" \
    --rctd-dir "$RCTD_OUT" \
    --methods binned segmented \
    "${annotation_args[@]}"

echo "$(date +%T)"
####################
