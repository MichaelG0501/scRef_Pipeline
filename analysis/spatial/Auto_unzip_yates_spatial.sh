#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=01:00:00
#PBS -N Auto_unzip_yates
#PBS -koed
#PBS -o temp/Auto_unzip_yates.out
#PBS -e temp/Auto_unzip_yates.err

set -e
echo "Start unzip: $(date '+%Y-%m-%d %H:%M:%S')"

TARGET_DIR=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Yates_spatial/source_data
cd "$TARGET_DIR"

export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE

echo "Extracting Zenodo_upload.zip with zipbomb detection disabled..."
unzip -q -o -x "__MACOSX/*" Zenodo_upload.zip

echo "Verifying unzipped contents..."
ls -la
if [ -d "Zenodo_upload" ]; then
    ls -la Zenodo_upload
    if [ -d "Zenodo_upload/processed" ]; then
        ls -lh Zenodo_upload/processed
    fi
fi

echo "Done unzip: $(date '+%Y-%m-%d %H:%M:%S')"
