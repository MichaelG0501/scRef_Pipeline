#!/bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=02:00:00
#PBS -N occams_ref_audit
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_bam_reference_audit.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS BAM reference audit"

module purge
module load tools/dev
module load SAMtools/1.18-GCC-12.3.0
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OCCAMS_EPHEMERAL="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS"
OUT_DIR="$WD/ref_outs/OCCAMS/reference_audit"
TABLE_DIR="$OUT_DIR/tables"
LOG_DIR="$OUT_DIR/logs"
SUMMARY_DIR="$WD/updates/new_updates/summaries"
AUDIT_CSV="$TABLE_DIR/occams_bam_reference_audit.csv"
SUMMARY_TXT="$LOG_DIR/occams_bam_reference_audit_summary.txt"
UPDATE_SUMMARY="$SUMMARY_DIR/occams_bam_reference_audit_summary.txt"

cd "$WD"
mkdir -p "$TABLE_DIR" "$LOG_DIR" "$SUMMARY_DIR"

printf '%s\n' "file_accession,bam_name,bam_path,header_status,sort_order,contig_naming,chr1_length,inferred_assembly,n_contigs" > "$AUDIT_CSV"

n_accessions=0
n_bams=0
n_invalid=0
n_grch37=0
n_grch38=0
n_other=0
n_chr=0
n_nochr=0
n_mixed=0

for accession_dir in "$OCCAMS_EPHEMERAL"/EGAF*; do
    [[ -d "$accession_dir" ]] || continue
    file_accession="$(basename "$accession_dir")"
    n_accessions=$((n_accessions + 1))

    mapfile -t bam_files < <(find "$accession_dir" -maxdepth 1 -type f -name '*.bam' -print | sort)
    if [[ ${#bam_files[@]} -ne 1 ]]; then
        printf '%s\n' "$file_accession,,,$([[ ${#bam_files[@]} -eq 0 ]] && echo missing || echo multiple),,,,,${#bam_files[@]}" >> "$AUDIT_CSV"
        n_invalid=$((n_invalid + 1))
        continue
    fi

    bam_path="${bam_files[0]}"
    bam_name="$(basename "$bam_path")"
    n_bams=$((n_bams + 1))

    if ! header="$(samtools view -H "$bam_path")"; then
        printf '%s\n' "$file_accession,$bam_name,$bam_path,unreadable,,,,," >> "$AUDIT_CSV"
        n_invalid=$((n_invalid + 1))
        continue
    fi

    sort_order="$(printf '%s\n' "$header" | awk -F'\t' '$1 == "@HD" {for (i=2; i<=NF; i++) if ($i ~ /^SO:/) {sub(/^SO:/, "", $i); print $i; exit}}')"
    [[ -n "$sort_order" ]] || sort_order="unknown"

    chr1_length="$(printf '%s\n' "$header" | awk -F'\t' '$1 == "@SQ" {sn=""; ln=""; for (i=2; i<=NF; i++) {if ($i ~ /^SN:/) {sn=$i; sub(/^SN:/, "", sn)}; if ($i ~ /^LN:/) {ln=$i; sub(/^LN:/, "", ln)}}; if (sn == "chr1" || sn == "1") {print ln; exit}}')"
    n_contigs="$(printf '%s\n' "$header" | awk -F'\t' '$1 == "@SQ" {n++} END {print n+0}')"
    has_chr="$(printf '%s\n' "$header" | awk -F'\t' '$1 == "@SQ" {for (i=2; i<=NF; i++) if ($i ~ /^SN:chr[0-9XYM]+$/) found=1} END {print found+0}')"
    has_nochr="$(printf '%s\n' "$header" | awk -F'\t' '$1 == "@SQ" {for (i=2; i<=NF; i++) if ($i ~ /^SN:([0-9]+|X|Y|MT)$/) found=1} END {print found+0}')"

    if [[ "$has_chr" -eq 1 && "$has_nochr" -eq 1 ]]; then
        contig_naming="mixed"
        n_mixed=$((n_mixed + 1))
    elif [[ "$has_chr" -eq 1 ]]; then
        contig_naming="chr"
        n_chr=$((n_chr + 1))
    elif [[ "$has_nochr" -eq 1 ]]; then
        contig_naming="nochr"
        n_nochr=$((n_nochr + 1))
    else
        contig_naming="other"
    fi

    case "$chr1_length" in
        249250621)
            inferred_assembly="GRCh37"
            n_grch37=$((n_grch37 + 1))
            ;;
        248956422)
            inferred_assembly="GRCh38"
            n_grch38=$((n_grch38 + 1))
            ;;
        *)
            inferred_assembly="other_or_missing"
            n_other=$((n_other + 1))
            ;;
    esac

    printf '%s\n' "$file_accession,$bam_name,$bam_path,ok,$sort_order,$contig_naming,$chr1_length,$inferred_assembly,$n_contigs" >> "$AUDIT_CSV"
done

{
    echo "OCCAMS BAM reference audit"
    echo "completed_at=$(date --iso-8601=seconds)"
    echo "accession_directories=$n_accessions"
    echo "bam_files=$n_bams"
    echo "invalid_or_non_singleton_accessions=$n_invalid"
    echo "grch37_bams=$n_grch37"
    echo "grch38_bams=$n_grch38"
    echo "other_or_missing_assembly_bams=$n_other"
    echo "chr_named_bams=$n_chr"
    echo "nochr_named_bams=$n_nochr"
    echo "mixed_named_bams=$n_mixed"
    if [[ "$n_accessions" -eq 300 && "$n_bams" -eq 300 && "$n_invalid" -eq 0 && "$n_grch37" -eq 300 && "$n_grch38" -eq 0 && "$n_other" -eq 0 ]]; then
        echo "result=PASS_ALL_300_GRCH37"
    else
        echo "result=FAIL_OR_MIXED_BUILD"
        exit_code=1
    fi
} | tee "$SUMMARY_TXT" "$UPDATE_SUMMARY"

if [[ "${exit_code:-0}" -ne 0 ]]; then
    exit "$exit_code"
fi

echo "$(date +%T) OCCAMS BAM reference audit complete"
