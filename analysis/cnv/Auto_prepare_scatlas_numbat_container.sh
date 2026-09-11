#!/bin/bash
####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/Auto_prepare_scatlas_numbat_container.sh
#   Methodology: not required (PBS/submit wrapper; method is documented by the invoked analysis script)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Execution wrapper; resources, dependencies, and arguments are defined below.
####################
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N Auto_scNBImg
#PBS -koed

set -euo pipefail

####################
# Prepare the Numbat container and install the package into a local library if
# the container image does not expose it directly.
####################

echo $(date +%T)
module purge

WD_LIVE="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
WD_EPH="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"

OUT_LIVE="${WD_LIVE}/ref_outs/Auto_scatlas_numbat"
OUT_EPH="${WD_EPH}/ref_outs/Auto_scatlas_numbat"

SIF="${OUT_LIVE}/Auto_numbat-rbase_latest.sif"
RLIB="${OUT_LIVE}/Rlib"
NUMBAT_SRC="${OUT_LIVE}/Auto_numbat_source"

mkdir -p "${OUT_EPH}/singularity_cache" "${OUT_EPH}/tmp" "${OUT_LIVE}/logs" "$RLIB"
export SINGULARITY_CACHEDIR="${OUT_EPH}/singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export TMPDIR="${OUT_EPH}/tmp"
export APPTAINER_TMPDIR="${TMPDIR}"

if [[ ! -f "$SIF" ]]; then
  apptainer build --mksquashfs-args "-no-xattrs" "$SIF" docker://pkharchenkolab/numbat-rbase:latest
fi

if [[ ! -f "${NUMBAT_SRC}/DESCRIPTION" ]]; then
  mkdir -p "$NUMBAT_SRC"
  apptainer exec --cleanenv -B /rds:/rds "$SIF" bash -lc "cd /numbat && tar cf - ." | tar xf - -C "$NUMBAT_SRC"
fi

if ! apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" Rscript -e 'quit(status = ifelse(requireNamespace("numbat", quietly=TRUE), 0, 1))'; then
  apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" R CMD INSTALL -l "$RLIB" "$NUMBAT_SRC"
fi

apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" Rscript -e 'stopifnot(requireNamespace("numbat", quietly=TRUE)); library(numbat); cat(as.character(packageVersion("numbat")), "\n"); stopifnot(file.exists("/numbat/inst/bin/pileup_and_phase.R")); stopifnot(exists("run_numbat")); stopifnot(exists("ref_hca"))'

echo $(date +%T)
