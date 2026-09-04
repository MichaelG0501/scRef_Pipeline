####################
# Analysis registry:
#   Status: active
#   Script: analysis/TCGA/tcga_esca_reconstruct_data.R
#   Methodology: analysis/methodology/TCGA/tcga_reconstruction_and_gender_validation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     - GDC API: TCGA-ESCA open RNA-Seq STAR-count gene expression files
#     - cBioPortal API: esca_tcga_gdc patient and sample clinical attributes
#   Outputs:
#     - ref_outs/TCGA/esca_gdc_reconstruction/raw/
#         Auto_gdc_star_counts_file_metadata.csv
#         Auto_gdc_manifest.txt
#         Auto_cbioportal_*_clinical_long.csv
#         Auto_cbioportal_*_clinical_wide.csv
#         gdc_files/<file_id>/*.rna_seq.augmented_star_gene_counts.tsv
#     - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/
#         Auto_tcga_esca_meta.rds/csv
#         Auto_tcga_esca_tpm_matrix.rds
#         Auto_tcga_esca_gene_key.csv
#     - ref_outs/TCGA/esca_gdc_reconstruction/tables/
#         TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#         Auto_tcga_esca_reconstruction_sample_summary.csv
#     - ref_outs/TCGA/esca_gdc_reconstruction/logs/
#         Auto_tcga_esca_reconstruction_run_summary.rds/txt
#     - compatibility copies:
#         ref_outs/tcga_esca_meta.rds
#         ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   Cache/replot behavior:
#     - Existing verified GDC files are reused.
#     - Set SCREF_TCGA_SKIP_DOWNLOAD=TRUE to process already-downloaded files only.
#     - Set SCREF_TCGA_OVERWRITE_BAD=TRUE only to replace a failed checksum file.
#   Run:
#     Rscript analysis/TCGA/tcga_esca_reconstruct_data.R
#   Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(httr2)
library(jsonlite)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_helpers.R")

setwd(SCREF_PROJECT_DIR)

####################
# 1) Paths and run settings
####################
script_path <- file.path(SCREF_ANALYSIS_DIR, "TCGA", "tcga_esca_reconstruct_data.R")
out_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "esca_gdc_reconstruction")
tiers <- ensure_output_tiers(out_dir)
raw_dir <- file.path(out_dir, "raw")
gdc_file_dir <- file.path(raw_dir, "gdc_files")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gdc_file_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SCREF_REF_OUTS_DIR, "cibersortx"), recursive = TRUE, showWarnings = FALSE)

gdc_files_endpoint <- "https://api.gdc.cancer.gov/files"
gdc_data_endpoint <- "https://api.gdc.cancer.gov/data/"
cbio_base <- "https://www.cbioportal.org/api"
study_id <- "esca_tcga_gdc"

options(timeout = as.numeric(Sys.getenv("SCREF_TCGA_DOWNLOAD_TIMEOUT", unset = "1800")))
skip_download <- tolower(Sys.getenv("SCREF_TCGA_SKIP_DOWNLOAD", unset = "FALSE")) %in% c("true", "1", "yes", "y")
overwrite_bad <- tolower(Sys.getenv("SCREF_TCGA_OVERWRITE_BAD", unset = "FALSE")) %in% c("true", "1", "yes", "y")
max_files_env <- Sys.getenv("SCREF_TCGA_MAX_FILES", unset = "")
max_files <- if (nzchar(max_files_env)) as.integer(max_files_env) else NA_integer_

run_summary <- start_run_summary(
  script = script_path,
  input_files = c(gdc_files_endpoint, cbio_base),
  output_files = c(
    file.path(tiers[["intermediate"]], "Auto_tcga_esca_meta.rds"),
    file.path(tiers[["tables"]], "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"),
    file.path(SCREF_REF_OUTS_DIR, "tcga_esca_meta.rds"),
    file.path(SCREF_REF_OUTS_DIR, "cibersortx", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
  ),
  parameters = list(
    study_id = study_id,
    skip_download = skip_download,
    overwrite_bad = overwrite_bad,
    max_files = ifelse(is.na(max_files), "all", max_files)
  )
)

####################
# 2) Small utilities
####################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

clean_missing <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "N/A", "not reported", "Not Reported", "'--", "--", "[Not Available]", "[Not Applicable]")] <- NA_character_
  x
}

clean_numeric <- function(x) {
  suppressWarnings(as.numeric(clean_missing(x)))
}

first_present <- function(...) {
  vals <- list(...)
  for (val in vals) {
    val <- clean_missing(val)
    if (length(val) > 0 && !is.na(val[1])) return(val[1])
  }
  NA_character_
}

coalesce_clean <- function(...) {
  vals <- lapply(list(...), clean_missing)
  if (length(vals) == 0) {
    return(NA_character_)
  }
  out <- vals[[1]]
  for (val in vals[-1]) {
    out <- ifelse(is.na(out), val, out)
  }
  out
}

add_missing_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA_character_
    }
  }
  df
}

normalise_sex <- function(x) {
  x <- tolower(clean_missing(x))
  dplyr::case_when(
    x %in% c("female", "f") ~ "Female",
    x %in% c("male", "m") ~ "Male",
    TRUE ~ NA_character_
  )
}

normalise_stage <- function(stage) {
  stage <- clean_missing(stage)
  dplyr::case_when(
    str_detect(stage, regex("Stage IV", ignore_case = TRUE)) ~ "Stage IV",
    str_detect(stage, regex("Stage III", ignore_case = TRUE)) ~ "Stage III",
    str_detect(stage, regex("Stage II", ignore_case = TRUE)) ~ "Stage II",
    str_detect(stage, regex("Stage I", ignore_case = TRUE)) ~ "Stage I",
    TRUE ~ NA_character_
  )
}

infer_histology <- function(type_vec, detailed_vec = NA_character_) {
  joined <- tolower(paste(clean_missing(type_vec), clean_missing(detailed_vec)))
  out <- rep("Other", length(joined))
  out[str_detect(joined, "adeno")] <- "EAC"
  out[str_detect(joined, "squamous")] <- "ESCC"
  out
}

verify_download <- function(path, expected_size = NA_real_, expected_md5 = NA_character_) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  if (!is.na(expected_size)) {
    observed_size <- file.info(path)$size
    if (is.na(observed_size) || observed_size != expected_size) {
      return(FALSE)
    }
  }
  if (!is.na(expected_md5) && nzchar(expected_md5)) {
    observed_md5 <- unname(tools::md5sum(path))
    if (is.na(observed_md5) || observed_md5 != expected_md5) {
      return(FALSE)
    }
  }
  TRUE
}

####################
# 3) Query current GDC STAR-count file metadata
####################
flatten_gdc_hit <- function(hit) {
  case <- hit$cases[[1]] %||% list()
  sample <- case$samples[[1]] %||% list()
  project <- case$project %||% list()
  sample_barcode <- sample$submitter_id %||% NA_character_
  tibble(
    file_id = hit$file_id %||% hit$id %||% NA_character_,
    file_name = hit$file_name %||% NA_character_,
    md5sum = hit$md5sum %||% NA_character_,
    file_size = as.numeric(hit$file_size %||% NA_real_),
    file_state = hit$state %||% NA_character_,
    data_format = hit$data_format %||% NA_character_,
    access = hit$access %||% NA_character_,
    project = project$project_id %||% "TCGA-ESCA",
    case_id = case$case_id %||% NA_character_,
    case_barcode = case$submitter_id %||% NA_character_,
    gdc_disease_type = case$disease_type %||% NA_character_,
    sample_barcode = sample_barcode,
    cbio_sample_id = substr(sample_barcode, 1, 16),
    sample_type_code = substr(sample_barcode, 14, 15),
    sample_type = sample$sample_type %||% NA_character_,
    tissue_type = sample$tissue_type %||% NA_character_,
    tumor_descriptor = sample$tumor_descriptor %||% NA_character_,
    specimen_type = sample$specimen_type %||% NA_character_,
    preservation_method = sample$preservation_method %||% NA_character_
  )
}

query_gdc_star_counts <- function() {
  filters <- list(
    op = "and",
    content = list(
      list(op = "=", content = list(field = "cases.project.project_id", value = list("TCGA-ESCA"))),
      list(op = "=", content = list(field = "data_category", value = list("Transcriptome Profiling"))),
      list(op = "=", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
      list(op = "=", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
      list(op = "=", content = list(field = "experimental_strategy", value = list("RNA-Seq"))),
      list(op = "=", content = list(field = "access", value = list("open")))
    )
  )

  fields <- paste(
    c(
      "file_id", "file_name", "md5sum", "file_size", "state", "data_format", "access",
      "cases.case_id", "cases.submitter_id", "cases.disease_type", "cases.project.project_id",
      "cases.samples.submitter_id", "cases.samples.sample_type", "cases.samples.tissue_type",
      "cases.samples.tumor_descriptor", "cases.samples.specimen_type", "cases.samples.preservation_method"
    ),
    collapse = ","
  )

  body <- list(
    filters = filters,
    fields = fields,
    format = "JSON",
    size = 5000,
    expand = "cases,cases.samples,cases.project"
  )

  response <- request(gdc_files_endpoint) |>
    req_body_json(body, auto_unbox = TRUE) |>
    req_retry(max_tries = 4) |>
    req_perform()

  parsed <- resp_body_json(response, simplifyVector = FALSE)
  hits <- parsed$data$hits
  if (length(hits) == 0) {
    stop("GDC query returned no TCGA-ESCA STAR-count files.")
  }
  bind_rows(lapply(hits, flatten_gdc_hit)) |>
    arrange(case_barcode, sample_barcode, file_name)
}

gdc_meta <- query_gdc_star_counts()
if (!is.na(max_files)) {
  gdc_meta <- head(gdc_meta, max_files)
}

metadata_path <- file.path(raw_dir, "Auto_gdc_star_counts_file_metadata.csv")
manifest_path <- file.path(raw_dir, "Auto_gdc_manifest.txt")
write.csv(gdc_meta, metadata_path, row.names = FALSE)
write_tsv(
  gdc_meta |>
    transmute(id = file_id, filename = file_name, md5 = md5sum, size = file_size, state = file_state),
  manifest_path
)
message("GDC metadata rows: ", nrow(gdc_meta))

####################
# 4) Pull cBioPortal patient and sample clinical attributes
####################
# bypassing API due to 503 errors; using local TSV
cbio_df <- fread("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/clinical_cbioportal.tsv", sep="\t", header=TRUE, data.table=FALSE)
colnames(cbio_df) <- gsub("[^A-Za-z0-9_]", "_", colnames(cbio_df))


####################
# 5) Download missing GDC files
####################
download_one_gdc_file <- function(file_id, file_name, md5sum, file_size) {
  dest_dir <- file.path(gdc_file_dir, file_id)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dest_dir, file_name)

  if (verify_download(dest, file_size, md5sum)) {
    return(tibble(file_id = file_id, path = dest, status = "exists_verified"))
  }

  if (file.exists(dest) && !overwrite_bad) {
    stop(
      "Existing GDC file failed size/checksum validation and SCREF_TCGA_OVERWRITE_BAD is not TRUE: ",
      dest
    )
  }

  if (skip_download) {
    stop("Missing or invalid GDC file while SCREF_TCGA_SKIP_DOWNLOAD is TRUE: ", dest)
  }

  url <- paste0(gdc_data_endpoint, file_id)
  last_error <- NULL
  for (attempt in seq_len(4)) {
    message("Downloading GDC file ", file_id, " (attempt ", attempt, "/4)")
    ok <- tryCatch(
      {
        utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE, method = "libcurl")
        TRUE
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        FALSE
      }
    )
    if (ok && verify_download(dest, file_size, md5sum)) {
      return(tibble(file_id = file_id, path = dest, status = "downloaded_verified"))
    }
    Sys.sleep(10 * attempt)
  }
  stop("Failed to download verified GDC file ", file_id, ". Last error: ", last_error %||% "checksum/size mismatch")
}

download_status <- vector("list", nrow(gdc_meta))
for (i in seq_len(nrow(gdc_meta))) {
  download_status[[i]] <- download_one_gdc_file(
    file_id = gdc_meta$file_id[i],
    file_name = gdc_meta$file_name[i],
    md5sum = gdc_meta$md5sum[i],
    file_size = gdc_meta$file_size[i]
  )
  if (i %% 10 == 0 || i == nrow(gdc_meta)) {
    message("Verified ", i, " / ", nrow(gdc_meta), " GDC files.")
  }
}
download_status <- bind_rows(download_status)
write.csv(download_status, file.path(raw_dir, "Auto_gdc_download_status.csv"), row.names = FALSE)

gdc_meta <- gdc_meta |>
  left_join(download_status |> select(file_id, path, download_status = status), by = "file_id")

####################
# 6) Build harmonized cBioPortal/GDC metadata
####################
patient_clin <- cbio_df %>%
  mutate(across(everything(), clean_missing)) %>%
  transmute(
    case_barcode = Patient_ID,
    Disease_Type = Disease_Type,
    vital_status = coalesce_clean(Patient_s_Vital_Status, ifelse(str_detect(Overall_Survival_Status, "^1:"), "Dead", "Alive")),
    OS_months = clean_numeric(Overall_Survival__Months_),
    OS_status = Overall_Survival_Status,
    OS_time = OS_months * 30.4375,
    OS_event = ifelse(str_detect(tolower(OS_status), "deceased|dead") | tolower(vital_status) == "dead", 1L, 0L),
    DFS_months = clean_numeric(Disease_Free__Months_),
    DFS_status = Disease_Free_Status,
    Gender = normalise_sex(Sex),
    Age_at_diagnosis = clean_numeric(Diagnosis_Age),
    Stage = AJCC_Pathologic_Stage,
    Stage_Simple = normalise_stage(Stage),
    Grade = NA_character_,
    AJCC_pathologic_T = AJCC_Pathologic_T_Stage,
    AJCC_pathologic_N = AJCC_Pathologic_N_Stage,
    AJCC_pathologic_M = AJCC_Pathologic_M_Stage,
    Race = Race_Category,
    Ethnicity = Ethnicity_Category,
    Prior_malignancy = Prior_Malignancy,
    Prior_treatment = Prior_Treatment,
    Year_of_diagnosis = clean_numeric(Year_of_Diagnosis)
  ) %>% distinct(case_barcode, .keep_all = TRUE)

sample_clin <- cbio_df %>%
  mutate(across(everything(), clean_missing)) %>%
  transmute(
    case_barcode = Patient_ID,
    cbio_sample_id = Sample_ID,
    Cancer_Type = Cancer_Type,
    Cancer_Type_Detailed = Cancer_Type_Detailed,
    Oncotree_Code = Oncotree_Code,
    Sample_Type_cBioPortal = Sample_Type,
    Mutation_count = clean_numeric(Mutation_Count),
    Fraction_genome_altered = clean_numeric(Fraction_Genome_Altered),
    TMB_nonsynonymous = clean_numeric(TMB__nonsynonymous_),
    Alcohol_history = Alcohol_History_Documented,
    Smoking_pack_years = clean_numeric(Person_Cigarette_Smoking_History_Pack_Year_Value),
    Primary_Diagnosis = Primary_Diagnosis,
    Patient_Primary_Tumor_Site = Patient_Primary_Tumor_Site
  ) %>% distinct(case_barcode, cbio_sample_id, .keep_all = TRUE)

meta <- gdc_meta |>
  left_join(patient_clin, by = "case_barcode") |>
  left_join(sample_clin, by = c("case_barcode", "cbio_sample_id")) |>
  mutate(
    type = if_else(!is.na(gdc_disease_type), gdc_disease_type, Disease_Type),
    HistologyGroup = infer_histology(type, Cancer_Type_Detailed),
    Stage_Simple = factor(Stage_Simple, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    Grade = factor(Grade, levels = c("G1", "G2", "G3", "G4")),
    Gender = factor(Gender, levels = c("Female", "Male"))
  ) |>
  arrange(case_barcode, sample_barcode, file_name)

meta_rds <- file.path(tiers[["intermediate"]], "Auto_tcga_esca_meta.rds")
meta_csv <- file.path(tiers[["intermediate"]], "Auto_tcga_esca_meta.csv")
saveRDS(meta, meta_rds)
write.csv(meta, meta_csv, row.names = FALSE)
saveRDS(meta, file.path(SCREF_REF_OUTS_DIR, "tcga_esca_meta.rds"))

####################
# 7) Build gene-symbol TPM mixture matrix
####################
read_one_star_file <- function(path) {
  x <- read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE)
  if (!"gene_id" %in% colnames(x)) {
    stop("Missing gene_id column in STAR-count file: ", path)
  }
  tpm_col <- intersect(c("tpm_unstranded", "tpm"), colnames(x))[1]
  if (is.na(tpm_col)) {
    stop("No TPM column found in STAR-count file: ", path)
  }
  gene_symbol_col <- intersect(c("gene_name", "GeneSymbol", "gene"), colnames(x))[1]
  if (is.na(gene_symbol_col)) {
    gene_symbol_col <- "gene_id"
  }
  gene_type_col <- intersect(c("gene_type", "gene_biotype"), colnames(x))[1]
  if (is.na(gene_type_col)) {
    x$gene_type_tmp <- NA_character_
    gene_type_col <- "gene_type_tmp"
  }
  x |>
    filter(!gene_id %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")) |>
    transmute(
      gene_id = sub("\\..*$", "", gene_id),
      gene_symbol = clean_missing(.data[[gene_symbol_col]]),
      gene_type = clean_missing(.data[[gene_type_col]]),
      tpm = as.numeric(.data[[tpm_col]])
    )
}

if (anyDuplicated(meta$sample_barcode)) {
  duplicated_samples <- unique(meta$sample_barcode[duplicated(meta$sample_barcode)])
  stop("Duplicate GDC sample barcodes found; refusing to build ambiguous matrix: ", paste(duplicated_samples, collapse = ", "))
}

if (any(is.na(meta$path)) || any(!file.exists(meta$path))) {
  missing_paths <- meta$path[is.na(meta$path) | !file.exists(meta$path)]
  stop("Some downloaded GDC files are missing; first missing path: ", missing_paths[1])
}

message("Reading first STAR-count file to establish gene key.")
first_expr <- read_one_star_file(meta$path[1])
gene_key <- first_expr |> select(gene_id, gene_symbol, gene_type)
expr_mat <- matrix(NA_real_, nrow = nrow(gene_key), ncol = nrow(meta))
colnames(expr_mat) <- meta$sample_barcode
rownames(expr_mat) <- gene_key$gene_id
expr_mat[, 1] <- first_expr$tpm

if (nrow(meta) > 1) {
  for (i in 2:nrow(meta)) {
    message("Reading STAR-count file ", i, " / ", nrow(meta), ": ", meta$sample_barcode[i])
    this_expr <- read_one_star_file(meta$path[i])
    if (!identical(this_expr$gene_id, gene_key$gene_id)) {
      idx <- match(gene_key$gene_id, this_expr$gene_id)
      if (any(is.na(idx))) {
        stop("Gene IDs in ", meta$path[i], " do not match the reference file.")
      }
      expr_mat[, i] <- this_expr$tpm[idx]
    } else {
      expr_mat[, i] <- this_expr$tpm
    }
  }
}

valid_gene <- !is.na(gene_key$gene_symbol) & nzchar(gene_key$gene_symbol)
dt_expr <- as.data.table(expr_mat[valid_gene, , drop = FALSE])
dt_expr[, GeneSymbol := gene_key$gene_symbol[valid_gene]]
setcolorder(dt_expr, c("GeneSymbol", setdiff(colnames(dt_expr), "GeneSymbol")))
final_bulk <- dt_expr[, lapply(.SD, sum, na.rm = TRUE), by = GeneSymbol, .SDcols = setdiff(colnames(dt_expr), "GeneSymbol")]
setorder(final_bulk, GeneSymbol)

tpm_matrix <- as.matrix(final_bulk[, -1, with = FALSE])
rownames(tpm_matrix) <- final_bulk$GeneSymbol

gene_key_path <- file.path(tiers[["intermediate"]], "Auto_tcga_esca_gene_key.csv")
tpm_matrix_rds <- file.path(tiers[["intermediate"]], "Auto_tcga_esca_tpm_matrix.rds")
mixture_path <- file.path(tiers[["tables"]], "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
compat_mixture_path <- file.path(SCREF_REF_OUTS_DIR, "cibersortx", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")

write.csv(gene_key, gene_key_path, row.names = FALSE)
saveRDS(tpm_matrix, tpm_matrix_rds)
fwrite(final_bulk, mixture_path, sep = "\t")
fwrite(final_bulk, compat_mixture_path, sep = "\t")

####################
# 8) Compact sample summary and run summary
####################
sample_summary <- meta |>
  count(HistologyGroup, sample_type_code, sample_type, Gender, name = "n_files") |>
  arrange(HistologyGroup, sample_type_code, Gender)
write.csv(
  sample_summary,
  file.path(tiers[["tables"]], "Auto_tcga_esca_reconstruction_sample_summary.csv"),
  row.names = FALSE
)

write_run_summary(
  run_summary,
  file.path(tiers[["logs"]], "Auto_tcga_esca_reconstruction_run_summary.rds")
)

message("TCGA ESCA reconstruction complete.")
message("Metadata: ", meta_rds)
message("TPM mixture: ", mixture_path)
message("Compatibility metadata: ", file.path(SCREF_REF_OUTS_DIR, "tcga_esca_meta.rds"))
message("Compatibility mixture: ", compat_mixture_path)
