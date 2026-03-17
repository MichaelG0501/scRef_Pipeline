####################
# Auto_cibersortx_reference.R
# Generate CIBERSORTx single-cell reference file from the full scATLAS
# (EAC_Ref_merged_strict.rds) using all cell types for unbiased deconvolution.
#
# Cell type labelling strategy:
#   - Epithelial cells present in meta_full_epi.rds → "Malignant"
#   - Epithelial cells NOT in meta_full_epi.rds → "Non_malignant_epithelial"
#   - Other cell types → celltype_update (dots replaced with underscores)
#   - "unresolved_inconsistent" cells → EXCLUDED
#
# Downsampling: min(100, all) per type, proportional fill to 3000 total.
#
# Input:
#   ref_outs/EAC_Ref_merged_strict.rds
#   ref_outs/meta_full_epi.rds
#
# Output:
#   ref_outs/cibersortx/CIBERSORTx_sc_reference.txt
#   ref_outs/cibersortx/CIBERSORTx_cell_labels.csv
#   updates/new_updates/summaries/Auto_cibersortx_reference_summary.csv
####################

library(Seurat)
library(data.table)
library(dplyr)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# Create output directory
####################
dir.create("cibersortx/", recursive = TRUE, showWarnings = FALSE)

####################
# Load data
####################
message("Loading EAC_Ref_merged_strict.rds ...")
merged <- readRDS("EAC_Ref_merged_strict.rds")

message("Loading meta_full_epi.rds ...")
meta_full_epi <- readRDS("meta_full_epi.rds")

####################
# Assign cell type labels for CIBERSORTx
####################
ct_raw <- as.character(merged$celltype_update)
names(ct_raw) <- Cells(merged)

# Identify epithelial cells present in meta_full_epi (these are the malignant/analysed ones)
epi_cells <- names(ct_raw)[ct_raw == "epithelial"]
epi_in_meta <- intersect(epi_cells, rownames(meta_full_epi))
epi_not_in_meta <- setdiff(epi_cells, rownames(meta_full_epi))

# Build label vector
cell_labels <- ct_raw
cell_labels[epi_in_meta] <- "Malignant"
cell_labels[epi_not_in_meta] <- "Non_malignant_epithelial"

# Replace dots with underscores (CIBERSORTx splits on dots)
cell_labels <- gsub("\\.", "_", cell_labels)

# Exclude unresolved_inconsistent cells
exclude_mask <- cell_labels == "unresolved_inconsistent"
message(sprintf("Excluding %d unresolved_inconsistent cells", sum(exclude_mask)))
keep_cells <- names(cell_labels)[!exclude_mask]
cell_labels <- cell_labels[keep_cells]
merged <- merged[, keep_cells]

# Report cell type distribution
label_table <- sort(table(cell_labels), decreasing = TRUE)
message("Cell type distribution:")
print(label_table)

####################
# Downsample: min(100, all) per type, proportional fill to 3000 total
####################
set.seed(42)
TARGET_TOTAL <- 3000
MIN_PER_TYPE <- 100

label_counts <- table(cell_labels)
all_types <- names(label_counts)

# Phase 1: guarantee min(100, all) per type
phase1 <- sapply(all_types, function(ct) {
  min(as.integer(label_counts[ct]), MIN_PER_TYPE)
})
names(phase1) <- all_types
phase1_total <- sum(phase1)

# Phase 2: distribute remaining slots proportionally
remaining <- TARGET_TOTAL - phase1_total
if (remaining > 0) {
  # Proportional to count, excluding already-saturated types
  eligible <- label_counts - phase1
  eligible[eligible < 0] <- 0
  if (sum(eligible) > 0) {
    extra_frac <- as.numeric(eligible) / sum(as.numeric(eligible))
    extra <- round(extra_frac * remaining)
    # Cap at available cells
    extra <- pmin(extra, as.integer(eligible))
    names(extra) <- all_types
    final_n <- phase1 + extra
  } else {
    final_n <- phase1
  }
} else {
  final_n <- phase1
}

message("Downsampling targets per type:")
print(final_n)
message(sprintf("Total target: %d", sum(final_n)))

# Sample cells
sampled_cells <- unlist(lapply(all_types, function(ct) {
  ct_cells <- names(cell_labels)[cell_labels == ct]
  n <- min(final_n[ct], length(ct_cells))
  if (n >= length(ct_cells)) return(ct_cells)
  sample(ct_cells, n)
}))

message(sprintf("Sampled %d cells total", length(sampled_cells)))

# Subset
merged_sub <- merged[, sampled_cells]
labels_sub <- cell_labels[sampled_cells]

####################
# Extract raw counts (NON-log, NON-normalized)
####################
mat <- GetAssayData(merged_sub, assay = "RNA", layer = "counts")
mat <- as.matrix(mat)

# Remove genes with zero expression across all sampled cells
gene_sums <- rowSums(mat)
mat <- mat[gene_sums > 0, , drop = FALSE]
message(sprintf("Genes retained: %d (removed %d zero-sum genes)",
                nrow(mat), sum(gene_sums == 0)))

####################
# Write CIBERSORTx sc_reference.txt
# Format: Row 1 = cell type labels, Col 1 = gene symbols
# Tab-delimited, raw counts
####################
message("Writing CIBERSORTx reference file ...")

# Header line: GeneSymbol followed by cell type labels (one per cell)
header_line <- paste(c("GeneSymbol", as.character(labels_sub)), collapse = "\t")

dt <- data.table(GeneSymbol = rownames(mat), mat)

outfile <- "cibersortx/CIBERSORTx_sc_reference.txt"
cat(header_line, file = outfile, sep = "\n")
fwrite(dt, file = outfile, sep = "\t", append = TRUE, col.names = FALSE)

message(sprintf("Reference file written: %s (%d genes x %d cells)",
                outfile, nrow(mat), ncol(mat)))

####################
# Save cell labels for reference
####################
labels_df <- data.frame(
  cell = sampled_cells,
  label = labels_sub,
  stringsAsFactors = FALSE
)
write.csv(labels_df, "cibersortx/CIBERSORTx_cell_labels.csv", row.names = FALSE)

####################
# Verify mixture file exists
####################
mixture_path <- "/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
if (file.exists(mixture_path)) {
  message(sprintf("Mixture file confirmed: %s", mixture_path))
  # Also copy to cibersortx/ for convenience
  file.copy(mixture_path, "cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt", overwrite = FALSE)
  message("Copied mixture file to cibersortx/ directory")
} else {
  warning("Mixture file NOT found at: ", mixture_path)
}

####################
# Summary CSV
####################
summary_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/updates/new_updates/summaries/"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- data.frame(
  cell_type = names(table(labels_sub)),
  n_cells_original = as.integer(label_counts[gsub("_", ".", names(table(labels_sub)))]),
  n_cells_sampled = as.integer(table(labels_sub)),
  stringsAsFactors = FALSE
)
# Fix: the original count lookup needs exact matching
summary_df$n_cells_original <- sapply(summary_df$cell_type, function(ct) {
  as.integer(label_counts[ct])
})
summary_df <- summary_df %>% arrange(desc(n_cells_sampled))

write.csv(summary_df, file.path(summary_dir, "Auto_cibersortx_reference_summary.csv"), row.names = FALSE)

message("CIBERSORTx reference generation complete.")
message("Upload the following files to CIBERSORTx web portal:")
message("  1. Single-cell reference: ref_outs/cibersortx/CIBERSORTx_sc_reference.txt")
message("  2. Mixture file: ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
message("  3. Use S-mode (batch correction) for high-resolution deconvolution")


##################

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
genes <- unique(as.vector(unlist(mp.genes)))
writeLines(genes, "cibersortx/CIBERSORTx_gene_subset.txt")
