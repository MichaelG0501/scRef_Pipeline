#!/bin/bash
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
#PBS -N occams_ega_meta
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp/occams_download_ega_metadata.log

set -euo pipefail

echo "$(date +%T) Starting OCCAMS EGA metadata download"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
DATASET_LIST="$WD/analysis/OCCAMS/OCCAMS_datasets.txt"
EGA_CREDENTIALS="/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json"
OUT_DIR="$WD/ref_outs/OCCAMS/source_data/ega_metadata_files"
cd "$WD"
mkdir -p "$OUT_DIR"

[[ -r "$EGA_CREDENTIALS" ]] || { echo "ERROR: EGA credential JSON is unreadable: $EGA_CREDENTIALS" >&2; exit 1; }
mapfile -t credential_values < <(python - "$EGA_CREDENTIALS" <<'PY'
import json
import sys
with open(sys.argv[1]) as handle:
    credentials = json.load(handle)
print(credentials["username"])
print(credentials["password"])
PY
)
[[ ${#credential_values[@]} -eq 2 ]] || { echo "ERROR: malformed EGA credential JSON" >&2; exit 1; }

token_response="$(curl --fail --silent --show-error --request POST \
    'https://idp.ega-archive.org/realms/EGA/protocol/openid-connect/token' \
    --data-urlencode 'client_id=metadata-api' \
    --data-urlencode 'grant_type=password' \
    --data-urlencode "username=${credential_values[0]}" \
    --data-urlencode "password=${credential_values[1]}")"
access_token="$(python -c 'import json,sys; print(json.load(sys.stdin)["access_token"])' <<< "$token_response")"
unset credential_values token_response
[[ -n "$access_token" ]] || { echo "ERROR: EGA access token was empty" >&2; exit 1; }

n_datasets=0
while IFS= read -r dataset_accession; do
    dataset_accession="${dataset_accession//$'\r'/}"
    [[ -n "$dataset_accession" ]] || continue
    n_datasets=$((n_datasets + 1))
    for endpoint in samples mappings/study_experiment_run_sample mappings/sample_file; do
        case "$endpoint" in
            samples) suffix="samples" ;;
            mappings/study_experiment_run_sample) suffix="study_experiment_run_sample" ;;
            mappings/sample_file) suffix="sample_file" ;;
        esac
        output="$OUT_DIR/${dataset_accession}_${suffix}.csv"
        if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == "TRUE" || ! -s "$output" ]]; then
            curl --fail --silent --show-error \
                --header "Authorization: Bearer $access_token" \
                --header 'Accept: text/csv' \
                "https://metadata.ega-archive.org/datasets/${dataset_accession}/${endpoint}" \
                --output "$output"
        fi
        [[ -s "$output" ]] || { echo "ERROR: empty EGA metadata output: $output" >&2; exit 1; }
    done
done < "$DATASET_LIST"
unset access_token

n_csv="$(find "$OUT_DIR" -maxdepth 1 -type f -name '*.csv' -printf '.' | wc -c)"
[[ "$n_datasets" -eq 52 && "$n_csv" -eq 156 ]] || { echo "ERROR: expected 52 datasets/156 CSVs, observed $n_datasets/$n_csv" >&2; exit 1; }
echo "datasets=52 metadata_csv_files=156 result=PASS"
echo "$(date +%T) OCCAMS EGA metadata download complete"
