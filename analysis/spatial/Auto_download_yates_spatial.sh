#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=02:00:00
#PBS -N Auto_dl_yates
#PBS -koed
#PBS -o temp/Auto_dl_yates.out
#PBS -e temp/Auto_dl_yates.err

set -euo pipefail
echo "Start: $(date '+%Y-%m-%d %H:%M:%S')"

TARGET_DIR=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Yates_spatial/source_data
mkdir -p "$TARGET_DIR"
mkdir -p /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/temp
cd "$TARGET_DIR"

echo "Downloading README.md..."
curl -sSL -o README.md "https://zenodo.org/records/15341263/files/README.md?download=1"

echo "Downloading Zenodo_upload.zip (~4.4 GB)..."
curl -sSL --retry 3 --retry-delay 5 -o Zenodo_upload.zip "https://zenodo.org/records/15341263/files/Zenodo_upload.zip?download=1"

echo "Verifying downloaded archive..."
ls -lh Zenodo_upload.zip

echo "Inspecting archive contents..."
unzip -l Zenodo_upload.zip | head -n 30

echo "Extracting Zenodo_upload.zip..."
unzip -q -o Zenodo_upload.zip

# Ensure Zenodo_upload directory layout is standardised
if [ -d "Zenodo_upload" ]; then
    echo "Zenodo_upload directory exists after unzipping."
else
    echo "Creating Zenodo_upload directory from extracted files..."
    mkdir -p Zenodo_upload
    if [ -d "processed" ]; then
        cp -r processed Zenodo_upload/
    fi
    if [ -f "sample_info.csv" ]; then
        cp sample_info.csv Zenodo_upload/
    fi
    if [ -d "raw" ]; then
        cp -r raw Zenodo_upload/
    fi
fi

echo "Summary of files in $TARGET_DIR:"
ls -lh "$TARGET_DIR"
if [ -d "$TARGET_DIR/Zenodo_upload" ]; then
    ls -lh "$TARGET_DIR/Zenodo_upload"
    if [ -d "$TARGET_DIR/Zenodo_upload/processed" ]; then
        ls -lh "$TARGET_DIR/Zenodo_upload/processed"
    fi
fi

echo "Done: $(date '+%Y-%m-%d %H:%M:%S')"
