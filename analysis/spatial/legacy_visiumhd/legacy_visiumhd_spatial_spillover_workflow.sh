#!/bin/bash
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_spatial_spillover_workflow.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=36:00:00
#PBS -N visiumhd_spatial
#PBS -koed
echo $(date +%T)
set -e
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
OUT=$WD/ref_outs/visium_hd_outs
MALIGNANCY_OUT=$OUT/malignancy
cd $WD

SEGMENTED_INPUTS=(
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs"
)
NAMES=("SUR1231" "FFPEA1" "FFPED1")

####################
# Annotation is independently calibrated per sample. Set the reuse variable
# only after checking the persisted calibration and FFPED1 lineage summary.
####################
if [[ "${SCREF_REUSE_SPATIAL_ANNOTATION:-FALSE}" != "TRUE" ]]; then
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/visiumhd_spatial_spillover_annotation.py \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --output-dir "$OUT" \
        --min-counts 200 \
        --max-mt 15 \
        --leiden-resolution 2
fi
####################

####################
# Annotation-only mode makes the calibration/recovery audit available before
# committing resources to InferCNA.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "annotation_only" ]]; then
    echo $(date +%T)
    exit 0
fi
####################

####################
# Downstream-only recovery reuses completed spatial InferCNA cell tables after
# a mapping or plotting failure.
if [[ "${SCREF_RUN_MODE:-all}" != "downstream_only" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_infercna_malignancy.R \
        --mode spatial \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --output-dir "$MALIGNANCY_OUT" \
        --min-reference-cells 20 \
        --min-epithelial-cells 30 \
        --cancer-signature-path "$WD/ref_outs/cancer_signatures.txt" \
        --cancer-signature-threshold 1 \
        --cna-sd-k 1
fi
####################

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
python analysis/spatial/process_visium_hd.py \
    --inputs "${SEGMENTED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --mode spatial \
    --stage map \
    --signature-dir "$OUT" \
    --output-dir "$OUT" \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --threshold 0.5 \
    --hybrid-gap 0.3

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
    --mode spatial \
    --inputs "${SEGMENTED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --annotation-dir "$OUT/tables" \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --output-dir "$OUT"

Rscript analysis/spatial/visiumhd_compare_annotation_infercna.R
echo $(date +%T)
