ALL_TERM2GENE <- list()
ALL_TERM2NAME <- list()

################################################################
################## Early embryogenesis #########################

library(readxl)
library(dplyr)

xlsx_path  <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/developmental/Early embryogenesis.xls"
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

xlsx_path <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/developmental/Organogenesis.xlsx"
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

xlsx_path <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/developmental/Normal development.xlsx"

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
xlsx_path <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/developmental/Oesophagus.xlsx"
sheets <- excel_sheets(xlsx_path)

library(Seurat)
tmdata <- readRDS("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/developmental/Stomach.rds")
epi <- subset(tmdata, major_clusters == "epi")
Idents(epi) <- tmdata@meta.data$subcluster.v2
markers <- FindAllMarkers(
  epi,
  only.pos = TRUE,       # only positive markers (upregulated genes)
  min.pct = 0.25,        # expressed in at least 25% of cells in cluster
  logfc.threshold = 0.25 # minimum log fold-change
)

oesophagus_order <- c("Quiescent basal cell", "Basal cell (cycling)", "Suprabasal", "Apical cell") %>% 
  gsub(" ", "_", .) %>% 
  paste0("_Oesophagus..Adult_Epithelium")

stomach_order <- c(
  "GKN+F", "ADH1+GKN1-F", "PG/Neck1", "PG/Neck2", "Chief", 
  "PC", "Ent", "Gob", "NE1", "NE2", "Pr_epi") %>%
  gsub(" ", "_", .) %>% 
    paste0("_Stomach..Adult_Epithelium")

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
    term = paste0(gsub(" ", "_", term), "_Stomach..Adult_Epithelium")
  )

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
saveRDS(marker_ref_all, file = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/enrich_dev.rds")

# Define the output directory
individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
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
  file_name <- paste0(individual_dir, "enrich_dev_", suffix, ".rds")
  
  saveRDS(current_ref, file = file_name)
}
