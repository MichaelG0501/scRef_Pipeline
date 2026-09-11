#!/bin/bash
#PBS -l select=1:ncpus=2:mem=24gb
#PBS -l walltime=24:00:00
#PBS -N Auto_Reconstruct_Drug_Refs
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev

# 1. Download ASGARD L1000 references
REF_ROOT=/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000
RAW_DIR=$REF_ROOT/raw
PLAIN_DIR=$REF_ROOT/plain
mkdir -p "$RAW_DIR" "$PLAIN_DIR"

download_one() {
  url=$1
  file=$2
  if [[ ! -s "$RAW_DIR/$file" ]]; then
    curl -L --fail --show-error --continue-at - --output "$RAW_DIR/$file" "$url"
  fi
  plain=${file%.gz}
  if [[ ! -s "$PLAIN_DIR/$plain" ]]; then
    gzip -dc "$RAW_DIR/$file" > "$PLAIN_DIR/$plain"
  fi
}

echo "Downloading ASGARD LINCS reference files..."
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz" "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz" "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz" "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz" "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz" "GSE92742_Broad_LINCS_cell_info.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz" "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz" "GSE92742_Broad_LINCS_sig_info.txt.gz"

# 2. Extract scDrugPrio reference datasets from the scDrugPrio package
echo "Extracting scDrugPrio PPI and drug-target tables..."

if [[ ! -s "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/ppi.txt" ]] || [[ ! -s "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/all_drug_targets_drug_bank.txt" ]]; then
  module load tools/dev
  eval "$(~/miniforge3/bin/conda shell.bash hook)"
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

  cat << 'EOF' > /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/extract_scdrugprio_refs.R
out_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/"

message("Downloading lit_ppi...")
load(url("https://raw.githubusercontent.com/SDTC-CPMed/scDrugPrio/main/data/lit_ppi.rda"))
write.table(lit_ppi, file=file.path(out_dir, "ppi.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

message("Downloading drug_bank_example_data...")
load(url("https://raw.githubusercontent.com/SDTC-CPMed/scDrugPrio/main/data/drug_bank_example_data.rda"))
write.table(drug_bank_example_data, file=file.path(out_dir, "all_drug_targets_drug_bank.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
EOF

  Rscript /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/extract_scdrugprio_refs.R
  rm /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/extract_scdrugprio_refs.R
else
  echo "scDrugPrio reference datasets already extracted."
fi

echo $(date +%T)
