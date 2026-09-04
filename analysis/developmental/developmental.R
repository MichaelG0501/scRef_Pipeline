####################
# Analysis registry:
#   Status: active
#   Script: analysis/developmental/developmental.R
#   Methodology: analysis/methodology/developmental/developmental_reference_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Build ordered TERM2GENE/TERM2NAME developmental references for
#     embryogenesis, organogenesis, normal development, adult epithelium, and Barrett's oesophagus.
#   Inputs:
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/{Early embryogenesis.xls,Organogenesis.xlsx,Normal development.xlsx,Oesophagus.xlsx,Stomach.rds}
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Barretts/science.abd1449_Table_S{2,4,5,7}.xlsx
#   Outputs:
#     - ref_outs/developmental_reference/enrich_dev.rds
#     - ref_outs/developmental_reference/per_stage/enrich_dev_<reference>.rds
#   Cache/replot: deterministic reference rebuild; outputs are persistent inputs to centred step 04.
#   Run: qsub analysis/developmental/developmental.sh
#   Conda env: dmtcp
####################

####################
# developmental.R
# Builds developmental gene reference from 4 external xlsx files:
# Early Embryogenesis, Organogenesis, Normal Development, Adult Epithelium.
# Outputs: ref_outs/enrich_dev.rds — used by enrichment_plotting.R
# Environment: dmtcp
####################

ALL_TERM2GENE <- list()
ALL_TERM2NAME <- list()

####################
# Persistent current output paths
####################
developmental_out_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/developmental_reference"
individual_dir <- file.path(developmental_out_dir, "per_stage")
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

################################################################
################## Early embryogenesis #########################

library(readxl)
library(dplyr)

xlsx_path  <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Early embryogenesis.xls"
sheet_name <- " temp.human.mk.tsv"
marker_term2gene <- read_excel(xlsx_path, sheet = sheet_name) %>%
  transmute(
    term = as.character(cluster),
    gene = as.character(gene),
    padj = as.numeric(p_val_adj)
  ) %>%
  filter(!is.na(padj), padj < 0.05, term != "", gene != "") %>%
  distinct(term, gene)
lineage_order <- c(
  "Z4cell", "8 cell", "Morula", "Prelineage", 
  "ICM", "Epiblast", "PriS", "Mesoderm", "Axial Mes", 
  "AdvMes", "Erythroblasts", "DE", "HEP", "ExE_Mes", 
  "Amnion", "Hypoblast", "YSE", 
  "TE", "CTB", "EVT", "STB"
)
marker_term2gene$term <- factor(marker_term2gene$term, levels = lineage_order)
marker_term2gene <- marker_term2gene %>% arrange(term)

marker_term2gene$term <- paste0(gsub(" ", "_", marker_term2gene$term), "..", "Early_Embryogenesis")
marker_term2name <- marker_term2gene %>%
  distinct(term) %>%
  transmute(term, name = term)

ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_term2gene 
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_term2name

################################################################
################## Organogenesis ###############################

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

xlsx_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Organogenesis.xlsx"
sheet_name <- "S1D"
marker_term2gene <- read_excel(xlsx_path, sheet = sheet_name) %>%
  # 1. Rename and select specific columns, then everything from column 4 onwards
  dplyr::select(
    ID = Type_id, 
    Major = Developmental_system, 
    term = Final_annotation, 
    4:ncol(.)
  ) %>%
  dplyr::mutate(across(everything(), as.character)) %>%
  # 2. Pivot ONLY the gene columns (exclude the 3 metadata columns)
  tidyr::pivot_longer(
    cols = -c(ID, Major, term), 
    values_to = "gene",
    values_drop_na = TRUE
  ) %>%
  # 3. Clean up strings
  dplyr::mutate(
    term = stringr::str_trim(term),
    ID = stringr::str_trim(ID),
    Major = stringr::str_trim(Major),
    gene = stringr::str_remove_all(gene, "^'+|'+$") %>% stringr::str_trim()
  ) %>%
  # 4. Filter out empty rows and Excel junk
  dplyr::filter(
    !is.na(term), term != "", 
    !is.na(gene), gene != "", 
    !stringr::str_detect(gene, "^\\.\\.\\.") 
  ) %>%
  # 5. Keep unique combinations
  dplyr::distinct(ID, Major, term, gene)

# 1. Define the specific order for Major systems
major_order <- c(
  "neural progenitor", "neuron", "epidermis", "sensory neuron", 
  "schwann", "craniofacial", "head mesoderm", "somite", "IM", 
  "somatic LPM", "limb", "splanchnic LPM", "endothelium", "blood", 
  "endoderm", "PGC", "epithelium", "fibroblast"
)

marker_major2gene <- marker_term2gene %>%
  # Filter to keep only the majors in your list
  dplyr::filter(Major %in% major_order) %>%
  # Map all sub-cluster genes to the Major term name
  dplyr::transmute(
    term = paste0(gsub(" ", "_", Major), "..Organogenesis_major"),
    gene = gene
  ) %>%
  # THE KEY STEP: Unique genes for the final term
  dplyr::distinct(term, gene) %>%
  # Apply the factor levels for the requested order
  dplyr::mutate(
    term = factor(term, levels = paste0(gsub(" ", "_", major_order), "..Organogenesis_major"))
  ) %>%
  dplyr::arrange(term)

marker_major2name <- marker_major2gene %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = as.character(term), name = as.character(term))


specific_order <- c(
  "foregut/esophagus", "thymus", "lung proximal epithelium and trachea",
  "lung distal epithelium", "hepatocyte-1", "hepatocyte-2", "stomach",
  "pancreas", "duodenum", "undefined", "epithelium-1", "epithelium-2",
  "epithelium-3", "epithelium-4", "epithelium-5"
)

# 2. Create the exact level names for sorting (OUTSIDE of mutate)
# This logic matches your naming rules
final_levels <- ifelse(
  grepl("epithelium", specific_order, ignore.case = TRUE),
  paste0(gsub(" ", "_", specific_order), "..Organogenesis_sub"),
  paste0(gsub(" ", "_", specific_order), "_Endoderm..Organogenesis_sub")
)

# 3. Process the dataframe
marker_specific2gene <- marker_term2gene %>%
  # Filter for the terms in your list
  dplyr::filter(term %in% specific_order) %>%
  # Ensure genes are unique within each specific term
  dplyr::distinct(term, gene) %>%
  # Apply conditional naming logic to create term_new
  dplyr::mutate(
    term_new = dplyr::if_else(
      grepl("epithelium", term, ignore.case = TRUE),
      paste0(gsub(" ", "_", term), "..Organogenesis_sub"),
      paste0(gsub(" ", "_", term), "_Endoderm..Organogenesis_sub")
    )
  ) %>%
  # Convert the new names into a factor using the levels we pre-calculated
  dplyr::mutate(
    term = factor(term_new, levels = final_levels)
  ) %>%
  # Arrange by the factor level and clean up
  dplyr::arrange(term) %>%
  dplyr::select(term, gene)

# 4. Prepare the name mapping
marker_specific2name <- marker_specific2gene %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = as.character(term), name = as.character(term))

ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_major2gene
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_major2name

ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_specific2gene
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_specific2name

################################################################
################## Normal development ##########################

library(readxl)
library(dplyr)
library(stringr)
library(tibble)

xlsx_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Normal development.xlsx"

# 1. Define the ground truth mapping and order
# Counts: Stomach (16), Intestine (12), Pancreas (14), Lung (13), Liver (9) = 64 Total
mapping_list <- data.frame(
  organ = c(
    rep("Stomach", 16), 
    rep("Intestine", 12), 
    rep("Pancreas", 14), 
    rep("Lung", 13), 
    rep("Liver", 9)
  ),
  original_term = c(
    # Stomach (16)
    "Goblet cells", "Parietal and chief cells", "Squamous epithelial cells", "Stromal cells", 
    "MUC13_DMBT1 positive cells", "Lymphoid cells", "Vascular endothelial cells", 
    "PDE1C_ACSM3 positive cells", "Myeloid cells", "Erythroblasts", "ENS glia", 
    "ENS neurons", "Ciliated epithelial cells", "Neuroendocrine cells", 
    "Lymphatic endothelial cells", "Mesothelial cells",
    # Intestine (12)
    "Intestinal epithelial cells", "Stromal cells", "ENS glia", "ENS neurons", 
    "Myeloid cells", "Lymphoid cells", "Vascular endothelial cells", "Smooth muscle cells", 
    "Chromaffin cells", "Erythroblasts", "Lymphatic endothelial cells", "Mesothelial cells",
    # Pancreas (14)
    "Acinar cells", "Stromal cells", "Ductal cells", "Lymphoid cells", 
    "Vascular endothelial cells", "Islet endocrine cells", "ENS glia", "ENS neurons", 
    "Erythroblasts", "Myeloid cells", "CCL19_CCL21 positive cells", "Mesothelial cells", 
    "Lymphatic endothelial cells", "Smooth muscle cells",
    # Lung (13)
    "Bronchiolar and alveolar epithelial cells", "Stromal cells", "Ciliated epithelial cells", 
    "Neuroendocrine cells", "Squamous epithelial cells", "Visceral neurons", 
    "Myeloid cells", "Lymphoid cells", "Megakaryocytes", "Vascular endothelial cells", 
    "Lymphatic endothelial cells", "Mesothelial cells", "CSH1_CSH2 positive cells",
    # Liver (9)
    "Vascular endothelial cells", "Lymphoid cells", "Myeloid cells", "Megakaryocytes", 
    "Erythroblasts", "Stellate cells", "Hepatoblasts", "Mesothelial cells", "Hematopoietic stem cells"
  )
)

# 2. Pre-calculate the Factor Levels to maintain the sequence
ordered_terms <- mapping_list %>%
  dplyr::mutate(
    formatted_name = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development")
  ) %>%
  dplyr::pull(formatted_name)

# 3. Process Table S4 (Long format gene lists)
marker_term2gene_long <- read_excel(xlsx_path, sheet = "Table_S4", skip = 1) %>%
  dplyr::select(
    original_term = `max.cluster`,
    gene_short_name,
    padj = qval
  ) %>%
  dplyr::inner_join(mapping_list, by = "original_term", relationship = "many-to-many") %>%
  dplyr::filter(as.numeric(padj) < 0.05) %>%
  dplyr::mutate(
    gene = stringr::str_remove_all(gene_short_name, "^'+|'+$"),
    term_string = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development"),
    term = factor(term_string, levels = unique(ordered_terms))
  ) %>%
  dplyr::distinct(term, gene) %>%
  dplyr::arrange(term) %>%
  dplyr::mutate(
    term = factor(paste0(as.character(term), "_long"), 
                  levels = paste0(levels(term), "_long"))
  )

# 4. Process Table S3 (Comma-separated genes in Column D)
marker_term2gene_S3 <- read_excel(xlsx_path, sheet = "Table_S3", skip = 1) %>%
  dplyr::select(
    organ = Organ,
    original_term = `Main cell type annotation`,
    gene_string = `Gene markers supporting annotation`
  ) %>%
  dplyr::inner_join(mapping_list, by = c("organ", "original_term")) %>%
  tidyr::separate_rows(gene_string, sep = ",\\s*") %>%
  dplyr::mutate(
    gene = stringr::str_trim(gene_string),
    term_string = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development"),
    term = factor(term_string, levels = unique(ordered_terms))
  ) %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::distinct(term, gene) %>%
  dplyr::arrange(term) %>%
  dplyr::mutate(
    term = factor(paste0(as.character(term), "_short"), 
                  levels = paste0(levels(term), "_short"))
  )

# 5. Prepare term2name mappings (Preserving Factor Levels)
marker_term2name_long <- marker_term2gene_long %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = term, name = as.character(term))

marker_term2name_S3 <- marker_term2gene_S3 %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = term, name = as.character(term))

# Append to your global collections
ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_term2gene_long
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_term2name_long

ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_term2gene_S3
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_term2name_S3

################################################################
################## Adult Epithelium ############################

library(readxl)
library(dplyr)
library(purrr)
xlsx_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Oesophagus.xlsx"
sheets <- excel_sheets(xlsx_path)

library(Seurat)
tmdata <- readRDS("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds")
epi <- subset(tmdata, major_clusters == "epi")
Idents(epi) <- tmdata@meta.data$subcluster.v2
markers <- FindAllMarkers(
  epi,
  only.pos = TRUE,       # only positive markers (upregulated genes)
  min.pct = 0.25,        # expressed in at least 25% of cells in cluster
  logfc.threshold = 0.25 # minimum log fold-change
)

####################
# Preserve existing oesophagus order (requested unchanged)
####################
oesophagus_order <- c("Quiescent basal cell", "Basal cell (cycling)", "Suprabasal", "Apical cell") %>% 
  gsub(" ", "_", .) %>% 
  paste0("_Oesophagus..Adult_Epithelium")

####################
# Enforce strict stomach naming and ordering for Adult_Epithelium
####################
stomach_rename_map <- c(
  "GKN+F" = "GKN+F_PylG_Stomach..Adult_Epithelium",
  "ADH1+GKN1-F" = "ADH1+GKN1-F_PylG_Stomach..Adult_Epithelium",
  "PG/Neck1" = "PG/Neck1_PylG_Stomach..Adult_Epithelium",
  "PG/Neck2" = "PG/Neck2_PylG_Stomach..Adult_Epithelium",
  "NE1" = "NE1_PylG_Stomach..Adult_Epithelium",
  "PC" = "PC_PylG_Stomach..Adult_Epithelium",
  "Chief" = "Chief_FundG_Stomach..Adult_Epithelium",
  "Pr_epi" = "Pr_epi_FundG_Stomach..Adult_Epithelium",
  "NE2" = "NE2_PylG/IntestMeta_Stomach..Adult_Epithelium",
  "Ent" = "Ent_IntestMeta_Stomach..Adult_Epithelium",
  "Gob" = "Gob_IntestMeta_Stomach..Adult_Epithelium"
)

stomach_order <- c(
  "GKN+F_PylG_Stomach..Adult_Epithelium",
  "ADH1+GKN1-F_PylG_Stomach..Adult_Epithelium",
  "PG/Neck1_PylG_Stomach..Adult_Epithelium",
  "PG/Neck2_PylG_Stomach..Adult_Epithelium",
  "NE1_PylG_Stomach..Adult_Epithelium",
  "PC_PylG_Stomach..Adult_Epithelium",
  "Chief_FundG_Stomach..Adult_Epithelium",
  "Pr_epi_FundG_Stomach..Adult_Epithelium",
  "NE2_PylG/IntestMeta_Stomach..Adult_Epithelium",
  "Ent_IntestMeta_Stomach..Adult_Epithelium",
  "Gob_IntestMeta_Stomach..Adult_Epithelium"
)

# Combine for the final factor levels
final_combined_order <- c(oesophagus_order, stomach_order)

# 2. Process Oesophagus Data
term2gene_oesophagus <- map_dfr(sheets, \(sh) {
  read_excel(xlsx_path, sheet = sh) |>
    transmute(term = sh, gene = .data[["gene"]], padj = .data[["p_val_adj"]]) |>
    filter(!is.na(gene), gene != "", !is.na(padj), padj < 0.05) |>
    distinct(term, gene)
}) %>%
  mutate(
    term = paste0(gsub(" ", "_", term), "_Oesophagus..Adult_Epithelium")
  )

# 3. Process Stomach Data
marker_term2gene_stomach <- markers %>%
  as_tibble() %>%
  dplyr::transmute(
    term = as.character(cluster),
    gene = as.character(gene),
    padj = as.numeric(p_val_adj)
  ) %>%
  dplyr::filter(!is.na(padj), padj < 0.05, term != "", gene != "") %>%
  dplyr::distinct(term, gene) %>%
  mutate(
    term = dplyr::recode(term, !!!stomach_rename_map, .default = NA_character_)
  ) %>%
  dplyr::filter(!is.na(term))

# 4. Combine and Apply Factor Levels for Ordering
marker_term2gene <- bind_rows(term2gene_oesophagus, marker_term2gene_stomach) %>%
  mutate(
    # Set the order based on the custom biological sequence defined above
    term = factor(term, levels = final_combined_order)
  ) %>%
  # Filter out any terms that weren't in your defined order (optional)
  filter(!is.na(term)) %>%
  # Physically sort the dataframe
  arrange(term)

# 5. Create Name Mapping (Preserving Factor Levels)
marker_term2name <- marker_term2gene %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = term, name = as.character(term))

# 6. Append to Global Lists
ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_term2gene 
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_term2name

################################################################
################## Barretts Oesophagus #########################

barrett_base_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Barretts/"

barretts_groups <- list(
  list(
    group_name = "Normal_Esophagus",
    file = "science.abd1449_Table_S4.xlsx",
    sheets = c("Basal" = "Basal",
               "Suprabasal" = "Suprabasal",
               "Suprabasal_Dividing" = "Suprabasal_Dividing",
               "Intermediate" = "Intermediate",
               "Superficial" = "Superficial")
  ),
  list(
    group_name = "Normal_Gastric",
    file = "science.abd1449_Table_S5.xlsx",
    sheets = c("Undifferentiated" = "Undifferentiated",
               "Undifferentiated_Dividing" = "Undifferentiated_Dividing",
               "Foveolar_Intermediate" = "Foveolar_Intermediate",
               "Foveolar_differentiated" = "Foveolar_differentiated",
               "Chief" = "Chief",
               "Parietal" = "Parietal",
               "Endocrine_GHRL" = "Endocrine_GHRL",
               "Endocrine_CHGA" = "Endocrine_CHGA",
               "Endocrine_NEUROD1" = "Endocrine_NEUROD1")
  ),
  list(
    group_name = "Barretts_Esophagus",
    file = "science.abd1449_Table_S7.xlsx",
    sheets = c("Columnar_Undifferentiated" = "Columnar_Undifferentiated",
               "Columnar_Dividing" = "Columnar_Undifferentiated_Divid",
               "Endocrine_NEUROG3" = "Endocrine_NEUROG3",
               "Columnar_Intermediate" = "Columnar_Intermediate",
               "Columnar_differentiated" = "Columnar_differentiated",
               "Goblet" = "Goblet")
  ),
  list(
    group_name = "Submucosal_Glands",
    file = "science.abd1449_Table_S2.xlsx",
    sheets = c("Duct_Intercalating" = "Duct_Intercalating",
               "Oncocytes" = "Oncocytes",
               "Mucous" = "Mucous")
  )
)

barrett_order <- c()
marker_term2gene_barrett_list <- list()

for (grp in barretts_groups) {
  xlsx_path <- file.path(barrett_base_dir, grp$file)
  available_sheets <- readxl::excel_sheets(xlsx_path)
  
  for (canonical in names(grp$sheets)) {
    sheet_name <- grp$sheets[[canonical]]
    
    actual_sheet <- sheet_name
    if (!actual_sheet %in% available_sheets) {
      match_idx <- grep(sheet_name, available_sheets, ignore.case=TRUE)
      if(length(match_idx) > 0){
         actual_sheet <- available_sheets[match_idx[1]]
      } else {
         warning(paste("Sheet not found:", sheet_name, "in", grp$file))
         next
      }
    }
    
    final_term_name <- paste0(canonical, "_", grp$group_name, "..Barretts_Oesophagus")
    barrett_order <- c(barrett_order, final_term_name)
    
    df <- readxl::read_excel(xlsx_path, sheet = actual_sheet)
    
    gene_col <- c("gene", "Symbol", "Genename", "Gene")
    gene_col <- gene_col[gene_col %in% names(df)][1]
    
    padj_col <- c("p_val_adj", "FDR", "qval", "adj.P.Val", "padj")
    padj_col <- padj_col[padj_col %in% names(df)][1]
    
    if (is.na(gene_col) || is.na(padj_col)) {
      warning(paste0("Missing gene/padj columns in sheet: ", actual_sheet, " of file ", grp$file))
      next
    }
    
    tmp_df <- df %>%
      dplyr::transmute(
        term = final_term_name,
        gene = as.character(.data[[gene_col]]),
        padj = suppressWarnings(as.numeric(.data[[padj_col]]))
      ) %>%
      dplyr::filter(!is.na(gene), gene != "", !is.na(padj), padj < 0.05) %>%
      dplyr::distinct(term, gene)
      
    marker_term2gene_barrett_list[[length(marker_term2gene_barrett_list) + 1]] <- tmp_df
  }
}

marker_term2gene_barrett <- dplyr::bind_rows(marker_term2gene_barrett_list) %>%
  dplyr::mutate(term = factor(term, levels = barrett_order)) %>%
  dplyr::arrange(term) %>%
  dplyr::filter(!is.na(term))

marker_term2name_barrett <- marker_term2gene_barrett %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = term, name = as.character(term))

ALL_TERM2GENE[[length(ALL_TERM2GENE) + 1]] <- marker_term2gene_barrett
ALL_TERM2NAME[[length(ALL_TERM2NAME) + 1]] <- marker_term2name_barrett

################################################################
################## Combined ####################################

master_levels <- purrr::map(ALL_TERM2GENE, ~ as.character(unique(.x$term))) %>% 
  purrr::flatten_chr()

# Combine all rows
combined_TERM2GENE <- dplyr::bind_rows(ALL_TERM2GENE) %>%
  dplyr::mutate(
    # Set the factor levels based on the global addition order
    term = factor(term, levels = master_levels)
  ) %>%
  # Arrange to physically sort the rows by the factor levels
  dplyr::arrange(term) %>%
  dplyr::distinct(term, gene)

# Create the name mapping preserving the factor
combined_TERM2NAME <- combined_TERM2GENE %>%
  dplyr::distinct(term) %>%
  dplyr::transmute(term = term, name = as.character(term))

# Final Reference List
marker_ref_all <- list(
  TERM2GENE = combined_TERM2GENE,
  TERM2NAME = combined_TERM2NAME
)

# Save the master file
saveRDS(marker_ref_all, file = file.path(developmental_out_dir, "enrich_dev.rds"))

# Define the output directory
if(!dir.exists(individual_dir)) dir.create(individual_dir, recursive = TRUE)

# Loop through and save each index
for (i in seq_along(ALL_TERM2GENE)) {
  
  # 1. Get the first term as a character string
  term_example <- as.character(ALL_TERM2GENE[[i]]$term[1])
  
  # 2. Extract the suffix (everything after the last '..')
  # If '..' is not found, it defaults to the full name
  suffix <- sub(".*\\.\\.", "", term_example)
  
  # 3. Construct the specific list for this element
  current_ref <- list(
    TERM2GENE = ALL_TERM2GENE[[i]],
    TERM2NAME = ALL_TERM2NAME[[i]]
  )
  
  # 4. Save with the requested naming format: enrich_dev_Suffix.rds
  # Using the index 'i' is still recommended to prevent overwriting if 
  # two list elements share the same suffix (like multiple 'Normal_Development' parts)
  file_name <- file.path(individual_dir, paste0("enrich_dev_", suffix, ".rds"))
  
  saveRDS(current_ref, file = file_name)
}
