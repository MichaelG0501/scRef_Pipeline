#!/bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=01:00:00
#PBS -N occams_archive_meta
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_archive_source_metadata.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS source-metadata archive"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
SRC="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS"
SOURCE_DIR="$WD/ref_outs/OCCAMS/source_data"
EGA_DIR="$SOURCE_DIR/ega_metadata_files"
TABLE_DIR="$WD/ref_outs/OCCAMS/tables"
HISTORICAL_LOG_DIR="$WD/ref_outs/OCCAMS/logs/historical"
SUMMARY_DIR="$WD/updates/new_updates/summaries"

cd "$WD"
mkdir -p "$SOURCE_DIR" "$EGA_DIR" "$TABLE_DIR" "$HISTORICAL_LOG_DIR" "$SUMMARY_DIR"

cp -p "$SRC/OCCAMS_metadata.csv" "$SOURCE_DIR/OCCAMS_metadata.csv"
cp -p "$SRC/OCCAMS_datasets.txt" "$SOURCE_DIR/OCCAMS_datasets.txt"
cp -p "$SRC"/ega_metadata_files/*.csv "$EGA_DIR/"
cp -p "$SRC/matched_samples_summary.csv" "$TABLE_DIR/matched_samples_summary.csv"
cp -p "$SRC/subject_bam_mapping.csv" "$TABLE_DIR/subject_bam_mapping.csv"
cp -p "$SRC/merge_commands.sh" "$HISTORICAL_LOG_DIR/merge_commands.sh"
cp -p "$SRC/pyega3_output.log" "$HISTORICAL_LOG_DIR/pyega3_output.log"
for log_file in "$SRC"/occams_merge_bams.[eo]* "$SRC"/occams_featurecounts.[eo]*; do
    [[ -f "$log_file" ]] || continue
    cp -p "$log_file" "$HISTORICAL_LOG_DIR/"
done

n_ega_csv="$(find "$EGA_DIR" -maxdepth 1 -type f -name '*.csv' -printf '.' | wc -c)"
{
    echo "OCCAMS source-metadata archive"
    echo "completed_at=$(date --iso-8601=seconds)"
    echo "ega_metadata_csv_files=$n_ega_csv"
    echo "master_metadata_bytes=$(stat -c '%s' "$SOURCE_DIR/OCCAMS_metadata.csv")"
    echo "matched_summary_rows=$(wc -l < "$TABLE_DIR/matched_samples_summary.csv")"
    echo "subject_mapping_rows=$(wc -l < "$TABLE_DIR/subject_bam_mapping.csv")"
    echo "credential_archived_in_repository=no"
    echo "credential_external_path=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json"
    echo "result=PASS"
} > "$WD/ref_outs/OCCAMS/logs/occams_archive_source_metadata_summary.txt"
cp "$WD/ref_outs/OCCAMS/logs/occams_archive_source_metadata_summary.txt" "$SUMMARY_DIR/occams_archive_source_metadata_summary.txt"
cat "$WD/ref_outs/OCCAMS/logs/occams_archive_source_metadata_summary.txt"
echo "$(date +%T) OCCAMS source-metadata archive complete"
