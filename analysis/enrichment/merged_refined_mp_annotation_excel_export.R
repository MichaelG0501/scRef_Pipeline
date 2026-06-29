####################
# Analysis registry:
#   Status: active
#   Script: analysis/enrichment/merged_refined_mp_annotation_excel_export.R
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below.
####################

####################
# Extract merged refined MP genes and create Excel summary with 2 pages
# Output: Metaprogrammes_Results/mp_refinement/tables/merged_refined_MP_genes_summary.xlsx
####################
library(openxlsx)
library(dplyr)

project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

# Load inputs
merged_mp_genes <- readRDS("Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_genes.rds")
cor_matrices <- readRDS("Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_correlation_matrices.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mean_rho <- cor_matrices$mean_rho

# MP description mapping
mp_desc_map <- c(
  "MP1"  = "G2M Cell Cycle", "MP9" = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition", "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.", "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.", "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair", "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)", "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi."
)

submp_desc_map <- c(
  "MP7+"  = "Fanconi/HR repair progenitor",
  "MP7h"  = "Replication-stress dormant epithelial",
  "MP7j"  = "DNA damage response",
  "MP7r"  = "Stem-like glandular duct progenitor",
  "MP7v"  = "Mucous secretory progenitor",
  "MP10+" = "Metabolic columnar epithelium",
  "MP10e" = "Inflammatory mucous-secretory columnar epithelium",
  "MP8+"  = "Intestinal metaplasia",
  "MP8b"  = "Quiescent glandular-metabolic progenitor",
  "MP8c"  = "NF-κB inflammatory cycling glandular progenitor",
  "MP8e"  = "Cycling intestinal–columnar progenitor",
  "MP12a" = "Enteroendocrine-primed progenitor",
  "MP12b" = "Enteroendocrine differentiation",
  "MP12c" = "Cycling glandular–intestinal progenitor",
  "MP2+"  = "MYC proliferation",
  "MP2v"  = "EMT-V cycling invasive progenitor",
  "MP15a" = "T/NK-like epithelial immune mimicry",
  "MP15b" = "Type I IFN–activated EMT-primed epithelial",
  "MP15c" = "Type II IFN / NF-κB peak inflammatory epithelial",
  "MP5+"  = "Epithelial IFN Resp.",
  "MP16+" = "Secretory Diff. (Gastric)"
)

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

get_desc <- function(mp_name) {
  if (mp_name %in% names(submp_desc_map)) {
    return(submp_desc_map[[mp_name]])
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
  "MP7j", "MP9", "MP1", "MP2+", "MP17", "MP8+", "MP10+", "MP14", "MP5+",
  "MP7r", "MP7v", "MP10e", "MP16+", "MP18",
  "MP8c", "MP15c", "MP12c", "MP2v", "MP8e", "MP12a", "MP13", 
  "MP7+", "MP7h", "MP8b", "MP12b", "MP15a", "MP15b"
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

output_path <- "Metaprogrammes_Results/mp_refinement/tables/merged_refined_MP_genes_summary.xlsx"
saveWorkbook(wb, output_path, overwrite = TRUE)

message("Saved: ", output_path)
