####################
# Analysis registry:
#   Status: active
#   Script: analysis/enrichment/merged_refined_mp_annotation_excel_export.R
#   Methodology: not required (direct table export of existing final MP objects)
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Export final centred refined MP genes, parent memberships, and
#     current full MP descriptions to a two-sheet workbook.
#   Inputs:
#     - ref_outs/Metaprogrammes_Results/centred/optimal_nMP.rds
#     - ref_outs/Metaprogrammes_Results/centred/geneNMF_metaprograms_nMP_<optimal>.rds
#     - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   Outputs: ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/merged_refined_MP_genes_summary.xlsx.
#   Cache/replot: deterministic workbook export; no hidden cache.
#   Run: qsub analysis/enrichment/merged_refined_mp_annotation_excel_export.sh
#   Conda env: dmtcp
####################

####################
# Extract merged refined MP genes and create Excel summary with 2 pages
# Output: Metaprogrammes_Results/centred/mp_refinement/tables/merged_refined_MP_genes_summary.xlsx
####################
library(openxlsx)
library(dplyr)

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

# Load inputs
optimal_nMP <- readRDS("Metaprogrammes_Results/centred/optimal_nMP.rds")
merged_mp_genes <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")
geneNMF.metaprograms <- readRDS(paste0("Metaprogrammes_Results/centred/geneNMF_metaprograms_nMP_", optimal_nMP, ".rds"))

# MP description mapping
mp_desc_map <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

get_desc <- function(mp_name) {
  if (mp_name %in% names(mp_desc_map)) {
    return(mp_desc_map[[mp_name]])
  }
  parent <- parent_id(mp_name)
  if (parent %in% names(mp_desc_map)) {
    return(mp_desc_map[[parent]])
  } else {
    return(paste0(parent, "_unknown"))
  }
}

merged_mps_ordered <- names(merged_mp_genes)

# Get baseline tree order
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order)]

# Hard-coded user-specified ordering for exact MP order
user_mp_order <- c(
  "MP1", "MP5", "MP13+",
  "MP2+",
  "MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+",
  "MP8+", "MP8b", "MP16", "MP18b", "MP17",
  "MP12",
  "MP15"
)
page1_order <- user_mp_order[user_mp_order %in% merged_mps_ordered]

# Logic for page 2 (Grouped_By_Parent)
# For by parent page, there should be small gap between distinct mps
all_parents <- parent_id(merged_mps_ordered)
unique_parents <- unique(all_parents)
unique_parent_ints <- as.integer(gsub("\\D", "", unique_parents))
ordered_unique_parents <- unique_parents[order(match(unique_parent_ints, mp_tree_order))]

page2_order <- character(0)
for (i in seq_along(ordered_unique_parents)) {
  p <- ordered_unique_parents[i]
  feats <- merged_mps_ordered[parent_id(merged_mps_ordered) == p]
  is_main <- !grepl("[a-z]$", feats)
  main_feat <- feats[is_main]
  if (length(main_feat) == 0) {
    main_feat <- p
  }
  sub_feats <- sort(feats[!is_main])
  page2_order <- c(page2_order, main_feat, sub_feats)
  if (i < length(ordered_unique_parents)) {
    page2_order <- c(page2_order, "GAP")
  }
}

build_mp_matrix <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)
  
  get_genes <- function(mp) {
    if (mp == "GAP") return(character(0))
    res <- merged_mp_genes[[mp]]
    if (is.null(res)) {
      nm <- gsub("\\+", "", mp)
      if (!is.null(geneNMF.metaprograms$metaprograms.genes[[nm]])) {
        res <- geneNMF.metaprograms$metaprograms.genes[[nm]]
      }
    }
    return(res)
  }
  
  max_g <- max(sapply(mp_names_vec, function(x) length(get_genes(x))))
  n_mp <- length(mp_names_vec)
  
  n_rows <- max_g + 2
  
  mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
  for (i in seq_along(mp_names_vec)) {
    mp <- mp_names_vec[i]
    if (mp == "GAP") {
      mat[1, i] <- ""
      mat[2, i] <- ""
    } else {
      mat[1, i] <- mp
      mat[2, i] <- get_desc(mp)
      genes <- get_genes(mp)
      if (length(genes) > 0) {
        mat[3:(length(genes)+2), i] <- genes
      }
    }
  }
  
  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

# Create DataFrames
df_p1 <- build_mp_matrix(page1_order)
df_p2 <- build_mp_matrix(page2_order)

wb <- createWorkbook()

mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
desc_style <- createStyle(fgFill = "#F2F2F2")

add_sheet <- function(wb, sheet_name, df, order_vec) {
  addWorksheet(wb, sheet_name)
  sheet_idx <- length(names(wb))
  
  if (!is.null(df)) {
    # Write data starting from row 1 (no RETAINED/REMOVED header)
    writeData(wb, sheet = sheet_idx, x = df, startCol = 1, startRow = 1, colNames = FALSE)
    
    # Set styles and widths per column
    for (i in seq_along(order_vec)) {
      if (order_vec[i] == "GAP") {
        setColWidths(wb, sheet_idx, cols = i, widths = 3)
      } else {
        setColWidths(wb, sheet_idx, cols = i, widths = 25)
        addStyle(wb, sheet = sheet_idx, mp_name_style, rows = 1, cols = i, gridExpand = TRUE)
        addStyle(wb, sheet = sheet_idx, desc_style, rows = 2, cols = i, gridExpand = TRUE)
      }
    }
  }
}

add_sheet(wb, "Split_Separated", df_p1, page1_order)
add_sheet(wb, "Grouped_By_Parent", df_p2, page2_order)

outdir <- "Metaprogrammes_Results/centred/mp_refinement/tables"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(outdir, "merged_refined_MP_genes_summary.xlsx")
saveWorkbook(wb, output_path, overwrite = TRUE)

message("Saved: ", output_path)
