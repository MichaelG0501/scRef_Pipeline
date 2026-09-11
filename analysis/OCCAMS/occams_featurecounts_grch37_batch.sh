#!/bin/bash
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=36:00:00
#PBS -N occams_fc_batch
#PBS -J 1-8
#PBS -koed

set -euo pipefail
echo "$(date +%T) Starting OCCAMS GRCh37 featureCounts batch ${PBS_ARRAY_INDEX}"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
FEATURECOUNTS="/sw-eb/software/Subread/2.0.6-GCC-12.3.0/bin/featureCounts"
GTF="$WD/ref_outs/OCCAMS/source_data/GENCODE_v19_GRCh37/gencode.v19.annotation.gtf.gz"
GTF_MD5="bd83e28270e595d3bde6bfcb21c9748f"
ALIASES="$WD/analysis/OCCAMS/grch37_contig_aliases.csv"
MANIFEST="$WD/ref_outs/OCCAMS/tables/occams_featurecounts_batch_manifest.csv"
OUT_DIR="/rds/general/project/tumourheterogeneity1/ephemeral/OCCAMS/counts/grch37_batches"
batch_number="$(printf '%02d' "$PBS_ARRAY_INDEX")"
output="$OUT_DIR/OCCAMS_RNAseq_GRCh37_batch_${batch_number}.txt"
cd "$WD"
mkdir -p "$OUT_DIR"

[[ -x "$FEATURECOUNTS" && -s "$GTF" && -s "$ALIASES" && -s "$MANIFEST" ]] || { echo "ERROR: required executable/reference/manifest missing" >&2; exit 1; }
[[ "$(md5sum "$GTF" | awk '{print $1}')" == "$GTF_MD5" ]] || { echo "ERROR: GTF checksum mismatch" >&2; exit 1; }
mapfile -t bam_files < <(awk -F',' -v batch="$PBS_ARRAY_INDEX" 'NR>1 && $1 == batch {gsub(/\r/, "", $4); print $4}' "$MANIFEST")
[[ ${#bam_files[@]} -gt 0 ]] || { echo "ERROR: batch $PBS_ARRAY_INDEX has no BAMs" >&2; exit 1; }
for bam in "${bam_files[@]}"; do [[ -e "$bam" ]] || { echo "ERROR: missing BAM: $bam" >&2; exit 1; }; done

if [[ "${SCREF_FORCE_REBUILD:-FALSE}" == "TRUE" || ! -s "$output" || ! -s "$output.summary" ]]; then
    "$FEATURECOUNTS" -T 4 -p --countReadPairs -s 0 -t exon -g gene_id \
        -A "$ALIASES" -a "$GTF" -o "$output" "${bam_files[@]}"
fi
n_samples="$(awk '!/^#/ {print NF-6; exit}' "$output")"
n_qc="$(awk -F'\t' 'NR==1 {print NF-1; exit}' "$output.summary")"
[[ "$n_samples" -eq ${#bam_files[@]} && "$n_qc" -eq ${#bam_files[@]} ]] || { echo "ERROR: batch output column mismatch" >&2; exit 1; }
echo "batch=$PBS_ARRAY_INDEX samples=$n_samples result=PASS"
echo "$(date +%T) OCCAMS GRCh37 featureCounts batch ${PBS_ARRAY_INDEX} complete"
