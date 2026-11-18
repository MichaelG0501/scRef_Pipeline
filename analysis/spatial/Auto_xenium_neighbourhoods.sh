#!/bin/bash
####################
# Analysis registry:
# Status: active
# Script: analysis/spatial/Auto_xenium_neighbourhoods.sh
# Description: Essential compute-node orchestration of Python input preparation,
#   R raw-count UCell scoring, Python spatial permutations and donor meta-analysis.
# Methodology: analysis/methodology/spatial/xenium_visiumhd_meta_methodology.md
# Inputs/outputs: exact live paths in spatial_platform_neighbourhoods.py and
#   xenium_hd_ucell_scores.R registries; no downloads or temporary dependencies.
# Run: /opt/pbs/bin/qsub -v SCREF_FORCE_REBUILD=TRUE
#   job_scripts/spatial/spatial_platform_neighbourhoods.sh
# Direct allocated-node run: bash analysis/spatial/Auto_xenium_neighbourhoods.sh
# Resources: 8 CPUs, 128 GB, 24 hours. Environments: dmtcp and bidcell_temp.
# Cache: exactly one of SCREF_FORCE_REBUILD, SCREF_CONTACT_ONLY, SCREF_META_ONLY or
#   SCREF_REPLOT_ONLY must be TRUE. Full rebuild recalculates raw-count scores;
#   contact-only reuses audited scores and repeats contact permutations only;
#   meta-only pools completed section tests and exports figures; replot reads
#   current meta tables only. No silent cache or analytical fallback.
#   SCREF_FIGURE_SCOPE=cross-platform-final restricts replotting to the final
#   cross-platform PDF and its index; all is the default.
#   SCREF_FIGURE_SCOPE=final-comparisons updates the two final bubble PDFs,
#   figure-source/summary tables and indexes without supporting figures.
####################
set -eo pipefail
echo "Start: $(date '+%Y-%m-%d %H:%M:%S')"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline
cd "$WD"
####################
export SCREF_NCORES="${SCREF_NCORES:-8}"
export SCREF_SPATIAL_PERMUTATIONS="${SCREF_SPATIAL_PERMUTATIONS:-4999}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
run_modes=0
for flag_name in SCREF_FORCE_REBUILD SCREF_CONTACT_ONLY SCREF_META_ONLY SCREF_REPLOT_ONLY; do
    flag_value="${!flag_name:-FALSE}"
    case "$flag_value" in
        TRUE) run_modes=$((run_modes + 1)) ;;
        FALSE) ;;
        *) echo "Invalid $flag_name: expected TRUE or FALSE" >&2; exit 1 ;;
    esac
done
if [[ "$run_modes" -ne 1 ]]; then
    echo 'Set exactly one of SCREF_FORCE_REBUILD, SCREF_CONTACT_ONLY, SCREF_META_ONLY or SCREF_REPLOT_ONLY to TRUE.' >&2
    exit 1
fi
conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp
if [[ "${SCREF_REPLOT_ONLY:-FALSE}" == TRUE ]]; then
    python analysis/spatial/spatial_platform_neighbourhoods.py --stage replot --figure-scope "${SCREF_FIGURE_SCOPE:-all}"
elif [[ "${SCREF_META_ONLY:-FALSE}" == TRUE ]]; then
    python analysis/spatial/spatial_platform_neighbourhoods.py --stage meta
else
    python analysis/spatial/spatial_platform_neighbourhoods.py --stage validate
    python analysis/spatial/spatial_platform_neighbourhoods.py --stage prepare
    if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == TRUE ]]; then
        conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
        Rscript analysis/spatial/xenium_hd_ucell_scores.R
        conda activate /rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp
    fi
    python analysis/spatial/spatial_platform_neighbourhoods.py --stage contacts
    if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == TRUE ]]; then
        python analysis/spatial/spatial_platform_neighbourhoods.py --stage meta
    fi
fi
####################
echo "End: $(date '+%Y-%m-%d %H:%M:%S')"
