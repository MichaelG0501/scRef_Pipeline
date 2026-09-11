#!/bin/bash
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=04:00:00
#PBS -N occams_grch37_test
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_prepare_grch37_reference.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS GRCh37 reference preparation and three-BAM validation"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
SOURCE_DIR="$WD/ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37"
TABLE_DIR="$WD/ref_outs/OCCAMS/reference_audit/tables"
LOG_DIR="$WD/ref_outs/OCCAMS/reference_audit/logs"
SUMMARY_DIR="$WD/updates/new_updates/summaries"
TEST_DIR="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/counts/grch37_three_bam_validation"
GTF="$SOURCE_DIR/gencode.v19.annotation.gtf.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
EXPECTED_MD5="bd83e28270e595d3bde6bfcb21c9748f"
ALIASES="$WD/analysis/OCCAMS/grch37_contig_aliases.csv"
FEATURECOUNTS="/sw-eb/software/Subread/2.0.6-GCC-12.3.0/bin/featureCounts"
CHR_BAM="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/EGAF00001902639/SLX-15341.D703_D507.bam"
NOCHR_BAM="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/EGAF00002559344/SLX-16453.D702-D503.bam"
MIXED_BAM="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/merged_bams/1791449d2f01848f.bam"
TEST_COUNTS="$TEST_DIR/occams_grch37_three_bam_counts.txt"
TEST_TABLE="$TABLE_DIR/occams_grch37_three_bam_featurecounts_summary.tsv"
SUMMARY_TXT="$LOG_DIR/occams_grch37_three_bam_validation_summary.txt"
UPDATE_SUMMARY="$SUMMARY_DIR/occams_grch37_three_bam_validation_summary.txt"

cd "$WD"
mkdir -p "$SOURCE_DIR" "$TABLE_DIR" "$LOG_DIR" "$SUMMARY_DIR" "$TEST_DIR"

if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == "TRUE" || ! -s "$GTF" ]]; then
    curl --fail --location --retry 3 --output "$GTF" "$GTF_URL"
fi

actual_md5="$(md5sum "$GTF" | awk '{print $1}')"
if [[ "$actual_md5" != "$EXPECTED_MD5" ]]; then
    echo "ERROR: GENCODE v19 GTF MD5 mismatch: expected $EXPECTED_MD5, observed $actual_md5" >&2
    exit 1
fi
gzip -t "$GTF"

{
    echo "source_url=$GTF_URL"
    echo "expected_md5=$EXPECTED_MD5"
    echo "observed_md5=$actual_md5"
    echo "sha256=$(sha256sum "$GTF" | awk '{print $1}')"
    echo "bytes=$(stat -c '%s' "$GTF")"
    echo "assembly=GRCh37.p13"
    echo "annotation=GENCODE_v19_comprehensive_CHR"
} > "$SOURCE_DIR/gencode.v19.annotation.provenance.txt"

"$FEATURECOUNTS" \
    -T 4 \
    -p --countReadPairs \
    -s 0 \
    -t exon \
    -g gene_id \
    -A "$ALIASES" \
    -a "$GTF" \
    -o "$TEST_COUNTS" \
    "$CHR_BAM" "$NOCHR_BAM" "$MIXED_BAM"

awk -F'\t' 'BEGIN {OFS="\t"}
    NR == 1 {
        for (i=2; i<=NF; i++) sample[i]=$i
        next
    }
    {
        for (i=2; i<=NF; i++) total[i]+=$i
        if ($1 == "Assigned") for (i=2; i<=NF; i++) assigned[i]=$i
        if ($1 == "Unassigned_NoFeatures") for (i=2; i<=NF; i++) nofeatures[i]=$i
        if ($1 == "Unassigned_MultiMapping") for (i=2; i<=NF; i++) multimapping[i]=$i
    }
    END {
        print "sample", "assigned_pct", "no_features_pct", "multi_mapping_pct", "annotation_compatible_pct"
        for (i=2; i<=length(sample)+1; i++) {
            print sample[i], 100*assigned[i]/total[i], 100*nofeatures[i]/total[i], 100*multimapping[i]/total[i], 100*(assigned[i]+multimapping[i])/total[i]
        }
    }' "$TEST_COUNTS.summary" > "$TEST_TABLE"

min_assigned="$(awk -F'\t' 'NR==2 {m=$2} NR>1 && $2<m {m=$2} END {print m+0}' "$TEST_TABLE")"
max_no_features="$(awk -F'\t' 'NR==2 {m=$3} NR>1 && $3>m {m=$3} END {print m+0}' "$TEST_TABLE")"
min_annotation_compatible="$(awk -F'\t' 'NR==2 {m=$5} NR>1 && $5<m {m=$5} END {print m+0}' "$TEST_TABLE")"

exit_code=0
{
    echo "OCCAMS GRCh37 three-BAM featureCounts validation"
    echo "completed_at=$(date --iso-8601=seconds)"
    echo "gtf_md5=$actual_md5"
    echo "tested_chr_named_bam=$(basename "$CHR_BAM")"
    echo "tested_nochr_named_bam=$(basename "$NOCHR_BAM")"
    echo "tested_mixed_named_bam=$(basename "$MIXED_BAM")"
    echo "minimum_assigned_pct=$min_assigned"
    echo "maximum_no_features_pct=$max_no_features"
    echo "minimum_annotation_compatible_pct=$min_annotation_compatible"
    if awk -v compatible="$min_annotation_compatible" -v nofeatures="$max_no_features" 'BEGIN {exit !(compatible >= 50 && nofeatures <= 40)}'; then
        echo "result=PASS"
    else
        echo "result=FAIL"
        exit_code=1
    fi
} > "$SUMMARY_TXT"
cp "$SUMMARY_TXT" "$UPDATE_SUMMARY"
cat "$SUMMARY_TXT"

if [[ "${exit_code:-0}" -ne 0 ]]; then
    exit "$exit_code"
fi

echo "$(date +%T) OCCAMS GRCh37 reference preparation and three-BAM validation complete"
