####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/cibersortx_sc_reference_export.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

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
####################
# org.Hs.eg.db for Ensembl → gene symbol conversion
####################
library(org.Hs.eg.db)

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
message(sprintf("Genes retained (before symbol conversion): %d (removed %d zero-sum genes)",
                nrow(mat), sum(gene_sums == 0)))

####################
# Convert Ensembl gene IDs to gene symbols
# The Seurat object has a mix of symbols and ENSG IDs as feature names.
# The mixture file uses gene symbols only (from tcga_data_prep.R via org.Hs.eg.db).
# Mismatched IDs cause S-mode batch correction to stall (near-zero gene overlap).
####################
gene_names <- rownames(mat)
is_ensg <- grepl("^ENSG[0-9]", gene_names)
message(sprintf("Gene IDs: %d Ensembl, %d symbols", sum(is_ensg), sum(!is_ensg)))

if (sum(is_ensg) > 0) {
  ensg_ids <- gene_names[is_ensg]
  mapped_symbols <- mapIds(org.Hs.eg.db,
                           keys = ensg_ids,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
  # Replace ENSG names with mapped symbols; drop unmapped (NA)
  new_names <- gene_names
  new_names[is_ensg] <- mapped_symbols[ensg_ids]
  
  # Remove genes that failed to map
  unmapped <- is.na(new_names)
  message(sprintf("Mapped %d Ensembl IDs to symbols; %d unmapped (dropped)",
                  sum(is_ensg & !unmapped), sum(unmapped)))
  mat <- mat[!unmapped, , drop = FALSE]
  new_names <- new_names[!unmapped]
  
  # Handle duplicates: keep the gene with highest mean expression
  if (any(duplicated(new_names))) {
    dup_count <- sum(duplicated(new_names))
    message(sprintf("Resolving %d duplicate gene symbols by max mean expression", dup_count))
    mean_expr <- rowMeans(mat)
    keep_idx <- tapply(seq_along(new_names), new_names, function(idx) {
      idx[which.max(mean_expr[idx])]
    })
    keep_idx <- unlist(keep_idx)
    mat <- mat[keep_idx, , drop = FALSE]
    new_names <- new_names[keep_idx]
  }
  
  rownames(mat) <- new_names
}
message(sprintf("Final gene count after symbol conversion: %d", nrow(mat)))

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
# Generate Merged Class File
# Groups 11 cell types into 5 broader classes to improve
# the samples:cell-types ratio for High-Resolution mode.
# Format: single tab-delimited line, same order as signature matrix phenotypes.
# The signature matrix is built from sc_reference, so phenotype order matches
# the unique cell type labels in order of first appearance in the header.
####################
message("Generating merged class file ...")

# Define mapping from fine cell types to merged classes
merge_map <- c(
  b_cell                  = "Lymphocyte",
  t_cell                  = "Lymphocyte",
  nk_cell                 = "Lymphocyte",
  plasma                  = "Lymphocyte",
  macrophage              = "Myeloid",
  dendritic               = "Myeloid",
  mast                    = "Myeloid",
  fibroblast              = "Stromal",
  endothelial             = "Stromal",
  Malignant               = "Malignant",
  Non_malignant_epithelial = "Epithelial"
)

# Get cell type order from the reference header (order of first appearance)
ref_type_order <- unique(as.character(labels_sub))
message("Signature matrix phenotype order: ", paste(ref_type_order, collapse = ", "))

# Map to merged classes
merged_classes <- merge_map[ref_type_order]
if (any(is.na(merged_classes))) {
  warning("Unmapped cell types: ", paste(ref_type_order[is.na(merged_classes)], collapse = ", "))
}

# Write single-line tab-delimited merged class file
merged_line <- paste(merged_classes, collapse = "\t")
merged_outfile <- "cibersortx/CIBERSORTx_merged_classes.txt"
cat(merged_line, file = merged_outfile, sep = "\n")
message(sprintf("Merged class file written: %s (%d types -> %d classes)",
                merged_outfile, length(ref_type_order), length(unique(merged_classes))))
message("  Mapping: ", paste(ref_type_order, "->", merged_classes, collapse = ", "))

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
message("  1. Single-cell reference (S-mode): ref_outs/cibersortx/CIBERSORTx_sc_reference.txt")
message("  2. Mixture file: ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
message("  3. Merged class file: ref_outs/cibersortx/CIBERSORTx_merged_classes.txt")
message("  4. Gene subset file: ref_outs/cibersortx/CIBERSORTx_gene_subset.txt")
message("  5. Use S-mode batch correction, disable QN (RNA-seq data)")


##################

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

unresolved_relabeled_path <- "unresolved_states/Auto_unresolved_relabel_states.rds"
new_state_gene_sets <- list()

if (file.exists(unresolved_relabeled_path)) {
  state_rel <- readRDS(unresolved_relabeled_path)
  
  candidate_new_states <- setdiff(
    unique(as.character(state_rel)),
    c(names(state_groups), "Unresolved", "Hybrid", NA)
  )
  
  nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
  
  if (file.exists(nmf3ca_path) && length(candidate_new_states) > 0) {
    MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
    MP_list <- as.list(MP_df)
    MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
    
    names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
    clean_map <- setNames(clean_3ca_name(names(MP_list)), names(MP_list))
    
    keep_cols <- names(clean_map)[clean_map %in% candidate_new_states]
    if (length(keep_cols) > 0) {
      new_state_gene_sets <- MP_list[keep_cols]
      names(new_state_gene_sets) <- clean_map[keep_cols]
    }
  }
}

genes <- unique(unlist(c(mp.genes, new_state_gene_sets), use.names = FALSE))
writeLines(genes, "cibersortx/CIBERSORTx_gene_subset.txt")

# ####################
# # OPTIONAL: Filter Reference & Mixture for High-Resolution Speed (Tutorial 5)
# # Highly recommended for web-portal runs to avoid timeouts.
# ####################
# message("Generating Filtered versions of inputs for High-Resolution fast run...")
# # 1. Load the mixture (produced by tcga_data_prep.R)
# mix_full <- fread("cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
# mix_f <- mix_full[GeneSymbol %in% genes]
# fwrite(mix_f, "cibersortx/CIBERSORTx_Mixture_Filtered.txt", sep = "\t")
# 
# # 2. Filter the sc_reference
# # Header matches 'labels_sub', so we reconstruct the filtered table
# ref_f <- dt[GeneSymbol %in% genes]
# out_ref_f <- "cibersortx/CIBERSORTx_sc_reference_Filtered.txt"
# cat(header_line, file = out_ref_f, sep = "\n")
# fwrite(ref_f, out_ref_f, sep = "\t", append = TRUE, col.names = FALSE)
# 
# message("  Fast/Filtered files written to ref_outs/cibersortx/ directory.")
# ####################
