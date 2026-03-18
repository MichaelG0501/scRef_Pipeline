####################
# Moved from: analysis/TCGA_data.R
# Reorganized as part of analysis/ restructuring
####################
library(readr)
library(dplyr)

setwd("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT")
ss <- read_tsv("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/TCGA/gdc_samples.tsv", show_col_types = FALSE)
names(ss) <- make.names(names(ss))

ss2 <- ss %>%
  transmute(
    file_id   = File.ID,
    file_name = File.Name,
    project   = Project.ID,
    case_barcode = Case.ID,         # e.g., TCGA-LN-A5U5
    sample_barcode = Sample.ID,     # e.g., TCGA-LN-A5U5-01A
    sample_type_code = substr(Sample.ID, 14, 15)
  ) #%>% filter(sample_type_code == "01")  # keep Primary Tumor only


library(stringr)

clinical <- read_tsv("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/TCGA/clinical.tsv")
grep("submitter_id$", colnames(clinical), value = TRUE)
grep("days_to_death|days_to_last_follow|vital_status", colnames(clinical), value = TRUE, ignore.case = TRUE)
clin <- clinical %>%
  transmute(
    case_barcode = cases.submitter_id, 
    type = cases.disease_type, 
    vital_status = demographic.vital_status,
    days_to_death = as.numeric(demographic.days_to_death),
    days_to_last_follow_up = as.numeric(diagnoses.days_to_last_follow_up), 
    Stage = diagnoses.ajcc_pathologic_stage,
    Stage_Simple = case_when(
      str_detect(Stage, "Stage IV") ~ "Stage IV",
      str_detect(Stage, "Stage III") ~ "Stage III",
      str_detect(Stage, "Stage II") ~ "Stage II",
      str_detect(Stage, "Stage I") ~ "Stage I",
      TRUE ~ NA_character_ # Set "Not Reported" to NA
    ),
    Grade = diagnoses.tumor_grade, # Usually G1, G2, G3
    Gender = demographic.gender,
    Location = diagnoses.site_of_resection_or_biopsy
  ) %>%
  mutate(
    OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS_event = ifelse(tolower(vital_status) == "dead", 1, 0), 
    Grade = factor(Grade, levels = c("G1", "G2", "G3")), 
    Stage_Simple = factor(Stage_Simple, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")), 
  ) %>%
  arrange(case_barcode, desc(!is.na(days_to_death))) %>%
  distinct(case_barcode, .keep_all = TRUE) %>%   # one row per patient
  filter(!is.na(OS_time))

# ---------------------------------------------------------------------------
# Merge additional clinical data from cBioPortal
# Key: cases.submitter_id (clinical.tsv) == Patient ID (cbioportal)
# ---------------------------------------------------------------------------
cbio_raw <- read_tsv(
  "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/TCGA/clinical_cbioportal.tsv",
  show_col_types = FALSE
)

# Confirm ID overlap before joining
clin_ids  <- unique(clin$case_barcode)
cbio_ids  <- unique(cbio_raw$`Patient ID`)
missing_in_cbio <- setdiff(clin_ids, cbio_ids)
missing_in_clin <- setdiff(cbio_ids, clin_ids)
message("cBioPortal ID check — in clin but not cbioportal (", length(missing_in_cbio), "): ",
        if (length(missing_in_cbio) == 0) "None" else paste(missing_in_cbio, collapse = ", "))
message("cBioPortal ID check — in cbioportal but not clin (", length(missing_in_clin), "): ",
        if (length(missing_in_clin) == 0) "None" else paste(missing_in_clin, collapse = ", "))

# Select and rename cBioPortal columns; one row per patient
cbio <- cbio_raw %>%
  transmute(
    case_barcode                  = `Patient ID`,
    # Overlapping fields — cBioPortal values replace clinical.tsv where they differ
    Stage                         = `AJCC Pathologic Stage`,
    Stage_Simple = case_when(
      str_detect(Stage, "Stage IV")  ~ "Stage IV",
      str_detect(Stage, "Stage III") ~ "Stage III",
      str_detect(Stage, "Stage II")  ~ "Stage II",
      str_detect(Stage, "Stage I")   ~ "Stage I",
      TRUE ~ NA_character_
    ),
    vital_status                  = `Patient's Vital Status`,
    OS_months                     = `Overall Survival (Months)`,
    DFS_months                    = `Disease Free (Months)`,
    DFS_status                    = `Disease Free Status`,
    Gender                        = Sex,
    Age_at_diagnosis              = `Diagnosis Age`,
    # Additional novel fields
    Alcohol_history               = `Alcohol History Documented`,
    Ethnicity                     = `Ethnicity Category`,
    Race                          = `Race Category`,
    Fraction_genome_altered       = `Fraction Genome Altered`,
    Mutation_count                = `Mutation Count`,
    TMB_nonsynonymous             = `TMB (nonsynonymous)`,
    Smoking_pack_years            = `Person Cigarette Smoking History Pack Year Value`,
    Prior_malignancy              = `Prior Malignancy`,
    Prior_treatment               = `Prior Treatment`,
    AJCC_pathologic_T             = `AJCC Pathologic T-Stage`,
    AJCC_pathologic_N             = `AJCC Pathologic N-Stage`,
    AJCC_pathologic_M             = `AJCC Pathologic M-Stage`,
    Year_of_diagnosis             = `Year of Diagnosis`
  ) %>%
  distinct(case_barcode, .keep_all = TRUE)  # one row per patient

# Left-join: all clin patients kept; cBioPortal columns appended / overlapping replaced
clin <- clin %>%
  # Drop columns that will be superseded by cBioPortal values
  dplyr::select(-Stage, -Stage_Simple, -vital_status, -Gender) %>%
  left_join(cbio, by = "case_barcode") %>%
  mutate(
    OS_event    = ifelse(tolower(vital_status) == "dead", 1, 0),
    Stage_Simple = factor(Stage_Simple, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
  )

meta <- ss2 %>%
  left_join(clin, by = "case_barcode")

####################
# Keep both major ESCA histologies (EAC + ESCC) instead of EAC-only
####################
meta <- meta[meta$type %in% c("Adenomas and Adenocarcinomas", "Squamous Cell Neoplasms") & !is.na(meta$type), ]
saveRDS(meta, "tcga_esca_meta.rds")

files <- list.files(path = "../", pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", 
                    recursive = TRUE, full.names = TRUE)
file_map <- tibble(path = files, file_name = basename(files))
meta2 <- meta %>%
  left_join(file_map, by = "file_name")
sum(is.na(meta2$path))  # should be 0 if everything matched

######################################
library(purrr)
read_one <- function(p) {
  x <- read_tsv(p, comment = "#", show_col_types = FALSE)
  x <- x %>% filter(!gene_id %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous"))
  x$gene_id <- sub("\\..*$", "", x$gene_id)  # drop version
  x %>% dplyr::select(gene_id, tpm_unstranded)
}

expr_list <- lapply(meta2$path, read_one)

gene_ids <- expr_list[[1]]$gene_id
expr_columns <- map(expr_list, ~ .x$tpm_unstranded) %>% bind_cols()
expr <- bind_cols(gene_id = gene_ids, expr_columns)
colnames(expr) <- c("gene_id", meta2$sample_barcode) 

library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
library(data.table)

# Assuming 'expr' is your existing dataframe with column 1 as "gene_id" (ENSG...)
# and columns 2:N as samples.

# 1. Map Ensembl IDs to Gene Symbols
# This removes version numbers if present (e.g., ENSG000001.1 -> ENSG000001)
expr$gene_id_clean <- sub("\\..*$", "", expr$gene_id)

# Map to Symbols
symbols <- mapIds(org.Hs.eg.db, keys = expr$gene_id_clean, 
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# 2. Assign Symbols & Filter
expr$GeneSymbol <- symbols
# Remove genes that didn't map to a symbol (NA)
expr_clean <- expr[!is.na(expr$GeneSymbol), ]

# 3. Aggregate Duplicates (Sum TPMs for same Symbol) using data.table
dt_expr <- as.data.table(expr_clean)
# Drop the old ID columns before aggregating
cols_to_remove <- c("gene_id", "gene_id_clean")
dt_expr[, (cols_to_remove) := NULL]

# Sum TPMs by GeneSymbol (Fastest method)
final_bulk <- dt_expr[, lapply(.SD, sum), by = GeneSymbol, .SDcols = names(dt_expr)[-ncol(dt_expr)]]

# 4. Write to File
fwrite(final_bulk, "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt", sep = "\t")

saveRDS(expr_clean, "tcga_esca_tpm_matrix.rds")

#############################

library(Seurat)
library(data.table)

# --- 1. Downsample (CRITICAL for CIBERSORTx limits) ---
# CIBERSORTx has a ~2GB upload limit and works best with ~5,000 cells max total.
Idents(tmdata_all) <- "manual_state"
# Define your states of interest
my_states <- c("Classic Proliferative", "Basal to Intest. Meta", "Intestinal Metaplasia", "Stress-adaptive")
clean_seurat <- subset(tmdata_all, idents = my_states)
clean_seurat <- subset(clean_seurat, downsample = 500)

mat <- GetAssayData(clean_seurat, assay = "RNA", slot = "counts")
mat <- as.matrix(mat) 

headers <- c("GeneSymbol", as.character(Idents(clean_seurat)))
dt <- data.table(GeneSymbol = rownames(mat), mat)

outfile <- "test.txt"
cat(paste(headers, collapse = "\t"), file = outfile, sep = "\n")
fwrite(dt, file = outfile, sep = "\t", append = TRUE, col.names = FALSE)

message("Done! Upload '", outfile, "' to CIBERSORTx.")
