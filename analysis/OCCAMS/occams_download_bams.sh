#!/bin/bash
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -l walltime=48:00:00
#PBS -N occams_ega_bams
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_download_bams.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS EGA BAM download"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/egaenv

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
OCCAMS_EPHEMERAL="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS"
MATCHED="$WD/ref_outs/OCCAMS/tables/matched_samples_summary.csv"
EGA_CREDENTIALS="/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json"
FAILURES="$OCCAMS_EPHEMERAL/occams_download_failures.tsv"
TABLE_DIR="$WD/ref_outs/OCCAMS/tables"
cd "$WD"
mkdir -p "$OCCAMS_EPHEMERAL" "$TABLE_DIR"

[[ -r "$EGA_CREDENTIALS" ]] || { echo "ERROR: EGA credential JSON is unreadable: $EGA_CREDENTIALS" >&2; exit 1; }
[[ -s "$MATCHED" ]] || { echo "ERROR: matched sample summary is missing: $MATCHED" >&2; exit 1; }
mapfile -t file_accessions < <(python - "$MATCHED" <<'PY'
import csv
import sys
with open(sys.argv[1], newline="") as handle:
    rows = csv.DictReader(handle)
    accessions = {
        row["File_Accession"].strip()
        for row in rows
        if row["Library_Strategy"].strip() == "RNA-Seq"
        and row["Phenotype"].strip() == "tumor"
    }
for accession in sorted(accessions):
    print(accession)
PY
)
[[ ${#file_accessions[@]} -eq 302 ]] || { echo "ERROR: expected 302 target EGAFs, observed ${#file_accessions[@]}" >&2; exit 1; }

printf '%s\n' $'File_Accession\tExit_Status' > "$FAILURES"
for file_accession in "${file_accessions[@]}"; do
    existing=("$OCCAMS_EPHEMERAL/$file_accession"/*.bam)
    if [[ ${#existing[@]} -eq 1 && -s "${existing[0]}" ]]; then
        continue
    fi
    if pyega3 -cf "$EGA_CREDENTIALS" fetch --output-dir "$OCCAMS_EPHEMERAL" "$file_accession"; then
        downloaded=("$OCCAMS_EPHEMERAL/$file_accession"/*.bam)
        if [[ ${#downloaded[@]} -ne 1 || ! -s "${downloaded[0]}" ]]; then
            printf '%s\t%s\n' "$file_accession" "missing_bam_after_success" >> "$FAILURES"
        fi
    else
        printf '%s\t%s\n' "$file_accession" "pyega3_failure" >> "$FAILURES"
    fi
done

cp "$FAILURES" "$TABLE_DIR/occams_download_failures.tsv"
python "$WD/analysis/OCCAMS/occams_metadata_coverage_audit.py"
echo "$(date +%T) OCCAMS EGA BAM download complete"
