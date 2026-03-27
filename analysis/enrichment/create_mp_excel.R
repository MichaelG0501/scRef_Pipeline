####################
# Extract MP genes and create Excel summary
# For PDO metaprograms
# Output: MP_genes_summary.xlsx
####################
library(openxlsx)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# Load MP outputs
MP_sc  <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

# Get mp_tree_order (same logic as Auto_mp_correlation_pdo.R)
tree_order <- MP_sc$programs.tree$order
ordered_clusters <- MP_sc$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order)]
mp_tree_order <- rev(mp_tree_order)
mp_tree_order <- paste0("MP", mp_tree_order)

# Get metrics for filtering
metrics <- MP_sc$metaprograms.metrics
rownames(metrics) <- paste0("MP", seq_len(nrow(metrics)))

# Identify filtered MPs
sil_bad <- rownames(metrics)[metrics$silhouette < 0]
drop_mps <- unique(c(sil_bad))

# All MPs
all_mps <- paste0("MP", seq_len(nrow(metrics)))
keep_mps <- setdiff(all_mps, drop_mps)

# Filter keep_mps to follow mp_tree_order
keep_ordered <- mp_tree_order[mp_tree_order %in% keep_mps]

# MP description mapping from the R script
mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

# Get all MP genes
mp_genes <- MP_sc$metaprograms.genes
names(mp_genes) <- paste0("MP", seq_len(length(mp_genes)))

# Get descriptions for all MPs
get_desc <- function(mp_name) {
  if (mp_name %in% names(mp_descriptions)) {
    return(mp_descriptions[mp_name])
  } else {
    return(paste0(mp_name, "_unknown"))
  }
}

# Split into removed and retained (both following mp_tree_order)
removed_mps <- mp_tree_order[mp_tree_order %in% drop_mps]
retained_mps <- keep_ordered

# Build MP matrix function
build_mp_matrix <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)
  
  max_g <- max(sapply(mp_names_vec, function(x) length(mp_genes[[x]])))
  n_mp <- length(mp_names_vec)
  
  # Row 1: MP names
  # Row 2: Description
  # Rows 3+: Genes
  n_rows <- max_g + 2
  
  mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
  for (i in seq_along(mp_names_vec)) {
    mp <- mp_names_vec[i]
    mat[1, i] <- mp  # MP name
    mat[2, i] <- get_desc(mp)  # Description
    genes <- mp_genes[[mp]]
    if (length(genes) > 0) {
      mat[3:(length(genes)+2), i] <- genes
    }
  }
  
  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

# Create dataframes (retained on left, removed on right)
df_ret <- build_mp_matrix(retained_mps)
df_rem <- build_mp_matrix(removed_mps)

# Calculate column positions
col_retained_start <- 1
col_removed_start <- if (!is.null(df_ret)) ncol(df_ret) + 2 else 1

# Create workbook
wb <- createWorkbook()
addWorksheet(wb, "MP_Genes")

# Write section labels in row 1
if (!is.null(df_ret) && nrow(df_ret) > 0) {
  writeData(wb, sheet = 1, x = "RETAINED MPs", startCol = 1, startRow = 1)
}
if (!is.null(df_rem) && nrow(df_rem) > 0) {
  writeData(wb, sheet = 1, x = "REMOVED MPs", startCol = col_removed_start, startRow = 1)
}

# Write data (starting row 2)
if (!is.null(df_ret)) {
  writeData(wb, sheet = 1, x = df_ret, startCol = 1, startRow = 2, colNames = FALSE)
}
if (!is.null(df_rem)) {
  writeData(wb, sheet = 1, x = df_rem, startCol = col_removed_start, startRow = 2, colNames = FALSE)
}

# Style: Section labels bold large
section_style <- createStyle(fontSize = 14, textDecoration = "bold", fgFill = "#FFC000")
if (!is.null(df_ret)) {
  addStyle(wb, sheet = 1, section_style, rows = 1, cols = 1, gridExpand = TRUE)
}
if (!is.null(df_rem)) {
  addStyle(wb, sheet = 1, section_style, rows = 1, cols = col_removed_start, gridExpand = TRUE)
}

# Style: MP names (row 2) bold
mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
if (!is.null(df_ret)) {
  addStyle(wb, sheet = 1, mp_name_style, rows = 2, cols = 1:ncol(df_ret), gridExpand = TRUE)
}
if (!is.null(df_rem)) {
  addStyle(wb, sheet = 1, mp_name_style, rows = 2, cols = col_removed_start:(col_removed_start + ncol(df_rem) - 1), gridExpand = TRUE)
}

# Style: Description row (row 3) with light gray background
desc_style <- createStyle(fgFill = "#F2F2F2")
if (!is.null(df_ret)) {
  addStyle(wb, sheet = 1, desc_style, rows = 3, cols = 1:ncol(df_ret), gridExpand = TRUE)
}
if (!is.null(df_rem)) {
  addStyle(wb, sheet = 1, desc_style, rows = 3, cols = col_removed_start:(col_removed_start + ncol(df_rem) - 1), gridExpand = TRUE)
}

# Set column widths
total_cols <- max(c(if(is.null(df_ret)) 0 else ncol(df_ret), if(is.null(df_rem)) 0 else ncol(df_rem)))
for (i in 1:total_cols) {
  setColWidths(wb, 1, cols = i, widths = 25)
}

# Save
output_path <- "MP_genes_summary.xlsx"
saveWorkbook(wb, output_path, overwrite = TRUE)

message("Saved: ", output_path)
message("")
message("Summary (following mp_tree_order):")
message("  Retained MPs (left): ", paste(retained_mps, collapse = ", "))
message("  Removed MPs (right): ", paste(removed_mps, collapse = ", "))
