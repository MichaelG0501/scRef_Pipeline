#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=72:00:00
#PBS -N occams_fc_grch37
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_featurecounts_grch37.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS GRCh37 featureCounts quantification"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OCCAMS_EPHEMERAL="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS"
FEATURECOUNTS="/sw-eb/software/Subread/2.0.6-GCC-12.3.0/bin/featureCounts"
GTF="$WD/ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz"
GTF_MD5="bd83e28270e595d3bde6bfcb21c9748f"
ALIASES="$WD/analysis/OCCAMS/grch37_contig_aliases.csv"
MAPPING="$OCCAMS_EPHEMERAL/subject_bam_mapping.csv"
OUT_DIR="$OCCAMS_EPHEMERAL/counts"
RAW_COUNTS="$OUT_DIR/OCCAMS_RNAseq_GRCh37_raw_counts.txt"
TABLE_DIR="$WD/ref_outs/OCCAMS/tables"
LOG_DIR="$WD/ref_outs/OCCAMS/logs"
SUMMARY_DIR="$WD/updates/new_updates/summaries"

cd "$WD"
mkdir -p "$OUT_DIR" "$TABLE_DIR" "$LOG_DIR" "$SUMMARY_DIR"

[[ -x "$FEATURECOUNTS" ]] || { echo "ERROR: featureCounts executable missing: $FEATURECOUNTS" >&2; exit 1; }
[[ -s "$GTF" ]] || { echo "ERROR: GRCh37 GTF missing: $GTF" >&2; exit 1; }
[[ "$(md5sum "$GTF" | awk '{print $1}')" == "$GTF_MD5" ]] || { echo "ERROR: GRCh37 GTF checksum mismatch" >&2; exit 1; }
[[ -s "$ALIASES" ]] || { echo "ERROR: chromosome alias file missing: $ALIASES" >&2; exit 1; }
[[ -s "$MAPPING" ]] || { echo "ERROR: subject-BAM mapping missing: $MAPPING" >&2; exit 1; }

mapfile -t bam_files < <(awk -F',' 'NR > 1 {gsub(/\r/, "", $3); print $3}' "$MAPPING")
if [[ ${#bam_files[@]} -ne 282 ]]; then
    echo "ERROR: expected 282 mapped BAMs, observed ${#bam_files[@]}" >&2
    exit 1
fi
for bam in "${bam_files[@]}"; do
    [[ -e "$bam" ]] || { echo "ERROR: mapped BAM is missing: $bam" >&2; exit 1; }
done

if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == "TRUE" || ! -s "$RAW_COUNTS" || ! -s "$RAW_COUNTS.summary" ]]; then
    "$FEATURECOUNTS" \
        -T 8 \
        -p --countReadPairs \
        -s 0 \
        -t exon \
        -g gene_id \
        -A "$ALIASES" \
        -a "$GTF" \
        -o "$RAW_COUNTS" \
        "${bam_files[@]}"
fi

n_genes="$(awk '!/^#/ {if (++n > 1) genes++} END {print genes+0}' "$RAW_COUNTS")"
n_columns="$(awk '!/^#/ {print NF; exit}' "$RAW_COUNTS")"
n_samples=$((n_columns - 6))
n_qc_samples="$(awk -F'\t' 'NR == 1 {print NF-1; exit}' "$RAW_COUNTS.summary")"

if [[ "$n_samples" -ne 282 || "$n_qc_samples" -ne 282 || "$n_genes" -le 0 ]]; then
    echo "ERROR: corrected featureCounts validation failed: genes=$n_genes counts_samples=$n_samples qc_samples=$n_qc_samples" >&2
    exit 1
fi

cp "$RAW_COUNTS.summary" "$TABLE_DIR/OCCAMS_RNAseq_GRCh37_raw_counts.txt.summary"

{
    echo "OCCAMS GRCh37 featureCounts run"
    echo "completed_at=$(date --iso-8601=seconds)"
    echo "annotation=GENCODE_v19_GRCh37.p13"
    echo "gtf_md5=$GTF_MD5"
    echo "genes=$n_genes"
    echo "samples=$n_samples"
    echo "raw_counts_ephemeral=$RAW_COUNTS"
    echo "result=PASS"
} > "$LOG_DIR/occams_featurecounts_grch37_summary.txt"
cp "$LOG_DIR/occams_featurecounts_grch37_summary.txt" "$SUMMARY_DIR/occams_featurecounts_grch37_summary.txt"
cat "$LOG_DIR/occams_featurecounts_grch37_summary.txt"

python "$WD/analysis/OCCAMS/occams_build_count_matrix.py"

echo "$(date +%T) OCCAMS GRCh37 featureCounts quantification complete"
