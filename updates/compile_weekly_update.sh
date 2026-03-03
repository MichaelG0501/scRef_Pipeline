#!/bin/bash
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=00:20:00
#PBS -N compile_weekly_update
echo $(date +%T)

module purge
module load tools/dev
module load texlive/20230313-GCC-13.2.0

export FONTCONFIG_PATH=/etc/fonts
fc-cache -f ~/.fonts/ 2>/dev/null
luaotfload-tool --update --force

# Usage: qsub -v TEXFILE=updates/02mar/update_02mar.tex compile_weekly_update.sh
# Or set TEXFILE below for a default.
TEXFILE=${TEXFILE:-updates/02mar/update_02mar.tex}

WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd "$WD"

# Two passes for TOC and cross-references
lualatex --interaction=nonstopmode "$TEXFILE"
lualatex --interaction=nonstopmode "$TEXFILE"

echo "Done: $(date +%T)"
echo "Output: ${TEXFILE%.tex}.pdf"
