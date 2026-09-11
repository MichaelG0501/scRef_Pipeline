#!/bin/bash
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_run_visium_hd_states.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=36:00:00
#PBS -N visium_hd_scatlas_states
#PBS -koed
echo $(date +%T)
set -e
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
OUT=$WD/ref_outs/visium_hd_outs
RCTD_OUT=$OUT/rctd
MALIGNANCY_OUT=$OUT/malignancy
cd $WD

####################
# Avoid expensive signature export and high-memory stages for the standalone
# evidence audit; it only reads final cached binned malignancy tables.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "keratinocyte_evidence_audit" ]]; then
    Rscript analysis/spatial/visiumhd_keratinocyte_evidence_audit.R \
        --samples SUR1231 FFPEA1 FFPED1 \
        --malignancy-dir "$MALIGNANCY_OUT"
    echo $(date +%T)
    exit 0
fi
####################

Rscript analysis/spatial/export_scatlas_visiumhd_signatures.R $OUT

BINNED_INPUTS=(
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um"
)
SEGMENTED_INPUTS=(
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs"
    "/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs"
)
NAMES=(
    "SUR1231"
    "FFPEA1"
    "FFPED1"
)

####################
# Validate the custom RCTD-doublet/manual-singlet annotation before a full
# InferCNA run. Existing RCTD tables must already be present and are reused.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "custom_annotation" ]]; then
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/process_visium_hd.py \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode custom \
        --stage annotate \
        --signature-dir $OUT \
        --output-dir $OUT \
        --rctd-dir "$RCTD_OUT" \
        --min-counts 1 \
        --max-mt 100 \
        --leiden-resolution 6
    echo $(date +%T)
    exit 0
fi
####################

####################
# A fast, PBS-safe annotation-only mode supports checking corrected segmented
# reference compartments before expensive InferCNA reruns.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "segmented_annotation" ]]; then
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    echo "Refreshing manually annotated segmented cells only..."
    python analysis/spatial/process_visium_hd.py \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode segmented \
        --stage annotate \
        --signature-dir $OUT \
        --output-dir $OUT \
        --min-counts 200 \
        --max-mt 15 \
        --leiden-resolution 1
    echo $(date +%T)
    exit 0
fi

####################
# Reuse completed segmented InferCNA tables when only state mapping and
# annotation diagnostics need recovery after a downstream plotting failure.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "segmented_map_diagnostics" ]]; then
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/process_visium_hd.py \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode segmented \
        --stage map \
        --signature-dir $OUT \
        --output-dir $OUT \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --threshold 0.5 \
        --hybrid-gap 0.3
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
        --mode segmented \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --output-dir "$OUT"
    echo $(date +%T)
    exit 0
fi
####################

####################
# Reuse completed binned outputs and a corrected segmented annotation table
# when only segmented InferCNA/state/diagnostic outputs need regeneration.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "segmented_downstream" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_infercna_malignancy.R \
        --mode segmented \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --output-dir "$MALIGNANCY_OUT" \
        --min-reference-cells 20 \
        --min-epithelial-cells 30 \
        --cancer-signature-path "$WD/ref_outs/cancer_signatures.txt" \
        --cancer-signature-threshold 1 \
        --cna-sd-k 1
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/process_visium_hd.py \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode segmented \
        --stage map \
        --signature-dir $OUT \
        --output-dir $OUT \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --threshold 0.5 \
        --hybrid-gap 0.3
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
        --mode segmented \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --output-dir "$OUT"
    echo $(date +%T)
    exit 0
fi

####################
# Correct cached binned CNA calls for segmented keratinocyte-dominant bins,
# then regenerate only binned state tables/maps and diagnostics.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "binned_keratinocyte_reclassify" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_reclassify_binned_keratinocyte_normals.R \
        --samples "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT"
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/process_visium_hd.py \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode binned \
        --stage map \
        --signature-dir $OUT \
        --output-dir $OUT \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --threshold 0.5 \
        --hybrid-gap 0.3
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
        --mode binned \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --output-dir "$OUT"
    echo $(date +%T)
    exit 0
fi

####################
# Reclassify keratinocyte-dominant bins from cached full InferCNA profiles,
# then regenerate binned state maps/diagnostics without CNA reinference.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "binned_keratinocyte_profile" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_profile_classify_keratinocyte_bins.R \
        --samples "${NAMES[@]}" \
        --malignancy-dir "$MALIGNANCY_OUT"
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    python analysis/spatial/process_visium_hd.py \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode binned \
        --stage map \
        --signature-dir $OUT \
        --output-dir $OUT \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --threshold 0.5 \
        --hybrid-gap 0.3
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
        --mode binned \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --output-dir "$OUT"
    echo $(date +%T)
    exit 0
fi

####################
# Plot-only keratinocyte versus epithelial evidence audit from final cached
# binned malignancy tables.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "keratinocyte_evidence_audit" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_keratinocyte_evidence_audit.R \
        --samples "${NAMES[@]}" \
        --malignancy-dir "$MALIGNANCY_OUT"
    echo $(date +%T)
    exit 0
fi


####################
# RCTD distinguishes singlet, doublet, and reject bins before epithelial
# restriction. Reuse the completed RCTD tables by default; set
# SCREF_REUSE_RCTD=FALSE only when a new RCTD run is explicitly required.
####################
if [[ "${SCREF_REUSE_RCTD:-TRUE}" != "TRUE" ]] || \
   [[ ! -f "$RCTD_OUT/tables/Auto_SUR1231_binned_rctd_annotations.csv.gz" ]] || \
   [[ ! -f "$RCTD_OUT/tables/Auto_FFPEA1_binned_rctd_annotations.csv.gz" ]] || \
   [[ ! -f "$RCTD_OUT/tables/Auto_FFPED1_binned_rctd_annotations.csv.gz" ]]; then
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_rctd_annotation.R \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --output-dir "$RCTD_OUT" \
        --reference-path "$WD/ref_outs/EAC_Ref_merged.rds" \
        --min-umis 200 \
        --max-cores 8
else
    echo "Reusing completed RCTD annotation tables."
fi

####################
# Custom 16 um annotation: RCTD is retained only for non-singlet bins, while
# RCTD singlets use the segmented manual marker/cluster/guard workflow.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "custom" ]]; then
    ####################
    # A calibrated annotation-only job can be reused while regenerating the
    # expensive custom CNA, mapping, and diagnostic outputs.
    ####################
    if [[ "${SCREF_REUSE_CUSTOM_ANNOTATION:-FALSE}" != "TRUE" ]]; then
        conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
        python analysis/spatial/process_visium_hd.py \
            --inputs "${BINNED_INPUTS[@]}" \
            --sample-names "${NAMES[@]}" \
            --mode custom \
            --stage annotate \
            --signature-dir $OUT \
            --output-dir $OUT \
            --rctd-dir "$RCTD_OUT" \
            --min-counts 1 \
            --max-mt 100 \
            --leiden-resolution 6
    else
        echo "Reusing calibrated complete custom annotation tables."
    fi
    ####################
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    ####################
    # InferCNA deliberately returns non-zero when any sample lacks two normal
    # reference types. Preserve that per-sample guard while continuing maps
    # for completed samples and annotation diagnostics for every sample.
    ####################
    set +e
    Rscript analysis/spatial/visiumhd_infercna_malignancy.R \
        --mode custom \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --output-dir "$MALIGNANCY_OUT" \
        --min-reference-cells 20 \
        --min-epithelial-cells 30 \
        --cancer-signature-path "$WD/ref_outs/cancer_signatures.txt" \
        --cancer-signature-threshold 1 \
        --cna-sd-k 1
    custom_cna_status=$?
    set -e
    custom_summary="$MALIGNANCY_OUT/tables/Auto_visiumhd_custom_infercna_malignancy_summary.csv"
    if [[ ! -f "$custom_summary" ]]; then
        echo "Custom InferCNA did not write its summary; aborting downstream steps."
        exit "$custom_cna_status"
    fi
    mapfile -t CUSTOM_COMPLETE_NAMES < <(awk -F',' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                header = $i
                gsub(/\"/, "", header)
                if (header == "status") status_column = i
            }
            next
        }
        {
            status_value = $status_column
            gsub(/\"/, "", status_value)
            if (status_value == "complete") {
                sample_value = $1
                gsub(/\"/, "", sample_value)
                print sample_value
            }
        }
    ' "$custom_summary")
    CUSTOM_COMPLETE_INPUTS=()
    for index in "${!NAMES[@]}"; do
        for complete_name in "${CUSTOM_COMPLETE_NAMES[@]}"; do
            if [[ "${NAMES[$index]}" == "$complete_name" ]]; then
                CUSTOM_COMPLETE_INPUTS+=("${BINNED_INPUTS[$index]}")
            fi
        done
    done
    ####################
    conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
    if [[ "${#CUSTOM_COMPLETE_NAMES[@]}" -gt 0 ]]; then
        python analysis/spatial/process_visium_hd.py \
            --inputs "${CUSTOM_COMPLETE_INPUTS[@]}" \
            --sample-names "${CUSTOM_COMPLETE_NAMES[@]}" \
            --mode custom \
            --stage map \
            --signature-dir $OUT \
            --output-dir $OUT \
            --malignancy-dir "$MALIGNANCY_OUT" \
            --threshold 0.5 \
            --hybrid-gap 0.3
    else
        echo "No custom samples completed InferCNA; state maps were not generated."
    fi
    conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
    Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
        --mode custom \
        --inputs "${BINNED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --annotation-dir "$OUT/tables" \
        --malignancy-dir "$MALIGNANCY_OUT" \
        --output-dir "$OUT"
    if [[ "$custom_cna_status" -ne 0 ]]; then
        echo "Custom InferCNA was incomplete; see $custom_summary for per-sample reference statistics."
    fi
    echo $(date +%T)
    exit 0
fi
####################

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
echo "Preparing RCTD-gated 16um binned annotations..."
python analysis/spatial/process_visium_hd.py \
    --inputs "${BINNED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --mode binned \
    --stage annotate \
    --signature-dir $OUT \
    --output-dir $OUT \
    --rctd-dir "$RCTD_OUT"

####################
# InferCNA selects exactly two distinct abundant normal reference types. State
# mapping includes malignant levels 1 and 2 after the signature rescue.
####################
conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visiumhd_infercna_malignancy.R \
    --mode binned \
    --inputs "${BINNED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --annotation-dir "$OUT/tables" \
    --output-dir "$MALIGNANCY_OUT" \
    --min-reference-cells 20 \
    --min-epithelial-cells 30 \
    --cancer-signature-path "$WD/ref_outs/cancer_signatures.txt" \
    --cancer-signature-threshold 1 \
    --cna-sd-k 1

####################
# Finalize keratinocyte-dominant RCTD epithelial bins from their complete CNA
# profiles before state mapping. This replaces the provisional spatial screen.
####################
Rscript analysis/spatial/visiumhd_profile_classify_keratinocyte_bins.R \
    --samples "${NAMES[@]}" \
    --malignancy-dir "$MALIGNANCY_OUT"
####################

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
echo "Mapping scATLAS states in malignant epithelial bins..."
python analysis/spatial/process_visium_hd.py \
    --inputs "${BINNED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --mode binned \
    --stage map \
    --signature-dir $OUT \
    --output-dir $OUT \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --threshold 0.5 \
    --hybrid-gap 0.3

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
    --mode binned \
    --inputs "${BINNED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --annotation-dir "$OUT/tables" \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --output-dir "$OUT"

####################
# Rebuild only the common-reference binned CNA, maps, and diagnostics when the
# segmented branch is already complete or is being regenerated independently.
####################
if [[ "${SCREF_RUN_MODE:-all}" == "binned_downstream" ]]; then
    echo $(date +%T)
    exit 0
fi
####################

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter
if [[ "${SCREF_REUSE_SEGMENTED_ANNOTATION:-FALSE}" != "TRUE" ]]; then
    echo "Preparing manually annotated segmented cells..."
    python analysis/spatial/process_visium_hd.py \
        --inputs "${SEGMENTED_INPUTS[@]}" \
        --sample-names "${NAMES[@]}" \
        --mode segmented \
        --stage annotate \
        --signature-dir $OUT \
        --output-dir $OUT \
        --min-counts 200 \
        --max-mt 15 \
        --leiden-resolution 1
else
    echo "Reusing corrected segmented annotation tables."
fi

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visiumhd_infercna_malignancy.R \
    --mode segmented \
    --inputs "${SEGMENTED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --annotation-dir "$OUT/tables" \
    --output-dir "$MALIGNANCY_OUT" \
    --min-reference-cells 20 \
    --min-epithelial-cells 30 \
    --cancer-signature-path "$WD/ref_outs/cancer_signatures.txt" \
    --cancer-signature-threshold 1 \
    --cna-sd-k 1

conda activate /rds/general/user/sg3723/home/miniforge3/envs/jupyter

echo "Mapping scATLAS states in malignant epithelial segmented cells..."
python analysis/spatial/process_visium_hd.py \
    --inputs "${SEGMENTED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --mode segmented \
    --stage map \
    --signature-dir $OUT \
    --output-dir $OUT \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --threshold 0.5 \
    --hybrid-gap 0.3

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/spatial/visiumhd_annotation_diagnostics.R \
    --mode segmented \
    --inputs "${SEGMENTED_INPUTS[@]}" \
    --sample-names "${NAMES[@]}" \
    --annotation-dir "$OUT/tables" \
    --malignancy-dir "$MALIGNANCY_OUT" \
    --output-dir "$OUT"

echo $(date +%T)
