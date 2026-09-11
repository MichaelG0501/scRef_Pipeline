#!/bin/bash
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=24:00:00
#PBS -N occams_merge_bams
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_merge_subject_bams.log

set -euo pipefail
echo "$(date +%T) Starting OCCAMS subject BAM preparation"
module purge
module load tools/dev
module load SAMtools/1.18-GCC-12.3.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OCCAMS_EPHEMERAL="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS"
COMMANDS="$OCCAMS_EPHEMERAL/merge_commands.sh"
MERGED="$OCCAMS_EPHEMERAL/merged_bams"
LOG_DIR="$WD/ref_outs/OCCAMS/logs"
SUMMARY_DIR="$WD/updates/new_updates/summaries"
cd "$WD"
mkdir -p "$LOG_DIR" "$SUMMARY_DIR"
python "$WD/analysis/OCCAMS/occams_prepare_subject_bams.py"
[[ -s "$COMMANDS" ]] || { echo "ERROR: merge command file missing" >&2; exit 1; }
bash "$COMMANDS"
n_bams="$(find "$MERGED" -maxdepth 1 \( -type f -o -type l \) -name '*.bam' -printf '.' | wc -c)"
[[ "$n_bams" -eq 282 ]] || { echo "ERROR: expected 282 subject BAMs, observed $n_bams" >&2; exit 1; }
{
    echo "OCCAMS subject BAM merge"
    echo "completed_at=$(date --iso-8601=seconds)"
    echo "subject_bams=282"
    echo "result=PASS"
} > "$LOG_DIR/occams_merge_subject_bams_summary.txt"
cp "$LOG_DIR/occams_merge_subject_bams_summary.txt" "$SUMMARY_DIR/occams_merge_subject_bams_summary.txt"
cat "$LOG_DIR/occams_merge_subject_bams_summary.txt"
echo "$(date +%T) OCCAMS subject BAM preparation complete"
