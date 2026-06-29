#!/bin/bash
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=24:00:00
#PBS -N scAtlas_VelSort
#PBS -koed

set -euo pipefail

####################
# Filter a raw Cell Ranger BAM to epithelial QC barcodes and coordinate-sort it
# for velocyto.
####################

echo $(date +%T)
module purge
module load tools/prod
module load SAMtools/1.22.1-GCC-14.2.0

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>"
  exit 1
fi

WD="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
OUT="${WD}/ref_outs/Auto_velocity_scATLAS"
manifest="${OUT}/tables/Auto_scatlas_velocity_sample_manifest.csv"
line=$(awk -F, -v s="$sample" 'NR > 1 && $1 == s {print; exit}' "$manifest")
if [[ -z "$line" ]]; then
  echo "ERROR: sample not found in manifest: $sample"
  exit 1
fi

BAM=$(awk -F, '{print $4}' <<< "$line")
BC=$(awk -F, '{print $6}' <<< "$line")
OUT_BAM="${OUT}/coord/${sample}.qc.coord.bam"
TMP_PREFIX="${OUT}/tmp_sort/${sample}.qc.coord"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: missing input BAM: $BAM"
  exit 1
fi
if [[ ! -s "$BC" ]]; then
  echo "ERROR: missing barcode file: $BC"
  exit 1
fi
if [[ -f "$OUT_BAM" ]]; then
  samtools quickcheck -v "$OUT_BAM"
  echo "Existing filtered BAM passed quickcheck for $sample; skipping."
  echo $(date +%T)
  exit 0
fi

mkdir -p "${OUT}/coord" "${OUT}/tmp_sort" "${OUT}/logs"
cd "$WD"

samtools view -h "$BAM" | \
awk -v bcfile="$BC" 'BEGIN { while ((getline line < bcfile) > 0) keep[line]=1 }
  /^@/ { print; next }
  {
    cb="";
    for (i=12; i<=NF; i++) {
      if ($i ~ /^CB:Z:/) cb=substr($i,6);
    }
    if (cb in keep) print;
  }' | \
samtools sort -@ 8 -m 7000M -T "$TMP_PREFIX" -O BAM -o "$OUT_BAM" -

samtools index -@ 8 "$OUT_BAM"
samtools quickcheck -v "$OUT_BAM"
samtools view -H "$OUT_BAM" | head -5

echo $(date +%T)
