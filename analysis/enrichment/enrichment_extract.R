#!/usr/bin/env Rscript
# Temp script to extract top N significant enrichment genes
# Early Embryogenesis (3), Stomach in Normal_Development_long (3), Organogenesis_sub (1), Organogenesis_major (4)

library(openxlsx)
library(dplyr)

# Load cluster enrichment results
cluster_enrich <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cluster_enrich.rds")

# Function to extract enrichment results
extract_enrich <- function(enrich_obj, top_n = 100) {
  if (is.null(enrich_obj)) return(NULL)
  res <- tryCatch(enrich_obj@result, error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # Filter significant
  res_sig <- res[res$p.adjust < 0.05, ]
  if (nrow(res_sig) == 0) return(NULL)
  
  # Get top n by p.adjust
  res_top <- res_sig %>% arrange(p.adjust) %>% head(top_n)
  
  # Extract overlap genes
  res_top$OverlapGenes <- sapply(strsplit(res_top$geneID, "/"), function(x) paste(x, collapse = ", "))
  
  return(res_top)
}

# Function: get TOP N significant across ALL MPs combined
process_category_topN <- function(cluster_enrich, category, top_n = 3) {
  all_results <- list()
  
  for (mp_name in names(cluster_enrich)) {
    enrich_obj <- cluster_enrich[[mp_name]][[category]]
    if (is.null(enrich_obj)) next
    
    res <- extract_enrich(enrich_obj, top_n = 100)
    if (is.null(res)) next
    
    res$MP <- mp_name
    all_results[[length(all_results) + 1]] <- res
  }
  
  if (length(all_results) == 0) return(NULL)
  
  combined <- bind_rows(all_results)
  
  # Sort by p.adjust and take TOP N overall
  combined <- combined %>% arrange(p.adjust) %>% head(top_n)
  
  return(combined)
}

# Function: filter for Stomach-related terms only
process_stomach_topN <- function(cluster_enrich, category, top_n = 3) {
  all_results <- list()
  
  for (mp_name in names(cluster_enrich)) {
    enrich_obj <- cluster_enrich[[mp_name]][[category]]
    if (is.null(enrich_obj)) next
    
    res <- extract_enrich(enrich_obj, top_n = 100)
    if (is.null(res)) next
    
    # Filter for Stomach-related terms
    res_stomach <- res[grep("Stomach", res$ID, ignore.case = TRUE), ]
    if (nrow(res_stomach) == 0) next
    
    res_stomach$MP <- mp_name
    all_results[[length(all_results) + 1]] <- res_stomach
  }
  
  if (length(all_results) == 0) return(NULL)
  
  combined <- bind_rows(all_results)
  
  # Sort by p.adjust and take TOP N overall
  combined <- combined %>% arrange(p.adjust) %>% head(top_n)
  
  return(combined)
}

# Extract requested categories
message("\n=== Processing categories ===")

message("1. Early_Embryogenesis (top 3)...")
early_embryo <- process_category_topN(cluster_enrich, "Early_Embryogenesis", top_n = 3)

message("2. Stomach in Normal_Development_long (top 3)...")
stomach_norm_dev <- process_stomach_topN(cluster_enrich, "Normal_Development_long", top_n = 3)

message("3. Organogenesis_sub (top 1)...")
organo_sub <- process_category_topN(cluster_enrich, "Organogenesis_sub", top_n = 1)

message("4. Organogenesis_major (top 4)...")
organo_major <- process_category_topN(cluster_enrich, "Organogenesis_major", top_n = 4)

# Create Excel workbook - only required columns
wb <- createWorkbook()

# Add Early Embryogenesis sheet
if (!is.null(early_embryo) && nrow(early_embryo) > 0) {
  early_embryo_clean <- early_embryo %>% 
    select(MP, ID, p.adjust, GeneRatio, OverlapGenes) %>%
    rename(
      TermID = ID,
      Pvalue = p.adjust,
      OverlappingGenes = OverlapGenes
    )
  addWorksheet(wb, "Early_Embryogenesis")
  writeData(wb, "Early_Embryogenesis", early_embryo_clean)
}

# Add Stomach Normal Development sheet
if (!is.null(stomach_norm_dev) && nrow(stomach_norm_dev) > 0) {
  stomach_clean <- stomach_norm_dev %>%
    select(MP, ID, p.adjust, GeneRatio, OverlapGenes) %>%
    rename(
      TermID = ID,
      Pvalue = p.adjust,
      OverlappingGenes = OverlapGenes
    )
  addWorksheet(wb, "Stomach_Normal_Dev_long")
  writeData(wb, "Stomach_Normal_Dev_long", stomach_clean)
}

# Add Organogenesis_sub sheet
if (!is.null(organo_sub) && nrow(organo_sub) > 0) {
  organo_sub_clean <- organo_sub %>%
    select(MP, ID, p.adjust, GeneRatio, OverlapGenes) %>%
    rename(
      TermID = ID,
      Pvalue = p.adjust,
      OverlappingGenes = OverlapGenes
    )
  addWorksheet(wb, "Organogenesis_sub")
  writeData(wb, "Organogenesis_sub", organo_sub_clean)
}

# Add Organogenesis_major sheet
if (!is.null(organo_major) && nrow(organo_major) > 0) {
  organo_major_clean <- organo_major %>%
    select(MP, ID, p.adjust, GeneRatio, OverlapGenes) %>%
    rename(
      TermID = ID,
      Pvalue = p.adjust,
      OverlappingGenes = OverlapGenes
    )
  addWorksheet(wb, "Organogenesis_major")
  writeData(wb, "Organogenesis_major", organo_major_clean)
}

# Save workbook
output_file <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/enrichment_topgenes_summary.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

message(paste("\nSaved Excel file to:", output_file))

# Print summary
message("\n=== TOP N SIGNIFICANT Results ===")
if (!is.null(early_embryo)) {
  message("\n--- Early Embryogenesis (top 3) ---")
  print(early_embryo %>% select(MP, ID, p.adjust, GeneRatio))
}
if (!is.null(stomach_norm_dev)) {
  message("\n--- Stomach Normal Development_long (top 3) ---")
  print(stomach_norm_dev %>% select(MP, ID, p.adjust, GeneRatio))
}
if (!is.null(organo_sub)) {
  message("\n--- Organogenesis_sub (top 1) ---")
  print(organo_sub %>% select(MP, ID, p.adjust, GeneRatio))
}
if (!is.null(organo_major)) {
  message("\n--- Organogenesis_major (top 4) ---")
  print(organo_major %>% select(MP, ID, p.adjust, GeneRatio))
}
