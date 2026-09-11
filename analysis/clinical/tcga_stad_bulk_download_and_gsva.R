####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/tcga_stad_bulk_download_and_gsva.R
#   Methodology: analysis/methodology/clinical/centred_stomach_bulk_location_gsva_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Downloads TCGA-STAD bulk RNA-seq STAR-count files and per-sample metadata.
#   Scores the final centred MPs and their two requested state-level gene sets
#   in TCGA-STAD bulk RNA-seq, then compares score distributions by primary
#   tumour location within the stomach. TCGA EAC samples are also scored and
#   included as a distinct location group for comparison.
#
# Inputs:
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#   - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds
#
# Outputs:
#   - EAC_Ref_all/00_merged/stomach_bulk/matrices/
#   - EAC_Ref_all/00_merged/stomach_bulk/metadata/
#   - ref_outs/Metaprogrammes_Results/centred/bulk_tcga_stad_location_gsva/
#
# Run command:
#   qsub analysis/clinical/tcga_stad_bulk_download_and_gsva.sh
#
# Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(GSVA)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(httr2)
  library(jsonlite)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
base_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/stomach_bulk"
raw_dir <- file.path(base_dir, "raw")
gdc_file_dir <- file.path(raw_dir, "gdc_files")
metadata_dir <- file.path(base_dir, "metadata")
matrix_dir <- file.path(base_dir, "matrices")
table_dir <- file.path(base_dir, "tables")
log_dir <- file.path(base_dir, "logs")

for (dir_path in c(raw_dir, gdc_file_dir, metadata_dir, matrix_dir, table_dir, log_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

project_id <- "TCGA-STAD"
cbio_study_id <- "stad_tcga_gdc"
gdc_files_endpoint <- "https://api.gdc.cancer.gov/files"
gdc_cases_endpoint <- "https://api.gdc.cancer.gov/cases"
gdc_data_endpoint <- "https://api.gdc.cancer.gov/data/"
cbio_base <- "https://www.cbioportal.org/api"

options(timeout = as.numeric(Sys.getenv("SCREF_STAD_DOWNLOAD_TIMEOUT", unset = "1800")))
skip_download <- tolower(Sys.getenv("SCREF_STAD_SKIP_DOWNLOAD", unset = "FALSE")) %in% c("true", "1", "yes", "y")
overwrite_bad <- tolower(Sys.getenv("SCREF_STAD_OVERWRITE_BAD", unset = "FALSE")) %in% c("true", "1", "yes", "y")
max_files_env <- Sys.getenv("SCREF_STAD_MAX_FILES", unset = "")
max_files <- if (nzchar(max_files_env)) as.integer(max_files_env) else NA_integer_

run_start <- Sys.time()
messages <- character()
log_msg <- function(...) {
  msg <- paste0(...)
  messages <<- c(messages, paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  message(msg)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  x
}

clean_missing <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "N/A", "not reported", "Not Reported", "'--", "--", "[Not Available]", "[Not Applicable]", "Unknown")] <- NA_character_
  x
}

clean_numeric <- function(x) suppressWarnings(as.numeric(clean_missing(x)))

coalesce_clean <- function(...) {
  vals <- lapply(list(...), clean_missing)
  if (length(vals) == 0) return(NA_character_)
  out <- vals[[1]]
  for (val in vals[-1]) {
    out <- ifelse(is.na(out), val, out)
  }
  out
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

clean_names <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep = "_")
}

ensure_columns <- function(df, cols) {
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- rep(NA_character_, nrow(df))
    }
  }
  df
}

flatten_list_columns <- function(df) {
  list_cols <- vapply(df, is.list, logical(1))
  for (col in names(df)[list_cols]) {
    df[[col]] <- vapply(
      df[[col]],
      function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        if (length(x) == 1 && !is.list(x)) return(as.character(x))
        jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")
      },
      character(1)
    )
  }
  df
}

verify_download <- function(path, expected_size = NA_real_, expected_md5 = NA_character_) {
  if (!file.exists(path)) return(FALSE)
  if (!is.na(expected_size)) {
    observed_size <- file.info(path)$size
    if (is.na(observed_size) || observed_size != expected_size) return(FALSE)
  }
  if (!is.na(expected_md5) && nzchar(expected_md5)) {
    observed_md5 <- unname(tools::md5sum(path))
    if (is.na(observed_md5) || observed_md5 != expected_md5) return(FALSE)
  }
  TRUE
}

api_post_json <- function(url, body, max_tries = 4) {
  request(url) |>
    req_body_json(body, auto_unbox = TRUE) |>
    req_retry(max_tries = max_tries) |>
    req_perform() |>
    resp_body_json(simplifyVector = FALSE)
}

api_get_json <- function(url, query = list(), max_tries = 4) {
  request(url) |>
    req_url_query(!!!query) |>
    req_retry(max_tries = max_tries) |>
    req_perform() |>
    resp_body_json(simplifyVector = FALSE)
}

flatten_gdc_file_hit <- function(hit) {
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
    project = project$project_id %||% project_id,
    case_id = case$case_id %||% NA_character_,
    case_barcode = case$submitter_id %||% NA_character_,
    gdc_disease_type = case$disease_type %||% NA_character_,
    gdc_primary_site = case$primary_site %||% NA_character_,
    sample_id = sample$sample_id %||% NA_character_,
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
      list(op = "=", content = list(field = "cases.project.project_id", value = list(project_id))),
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
      "cases.case_id", "cases.submitter_id", "cases.disease_type", "cases.primary_site",
      "cases.project.project_id", "cases.samples.sample_id", "cases.samples.submitter_id",
      "cases.samples.sample_type", "cases.samples.tissue_type", "cases.samples.tumor_descriptor",
      "cases.samples.specimen_type", "cases.samples.preservation_method"
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
  parsed <- api_post_json(gdc_files_endpoint, body)
  hits <- parsed$data$hits
  if (length(hits) == 0) stop("GDC query returned no TCGA-STAD STAR-count files.")
  bind_rows(lapply(hits, flatten_gdc_file_hit)) |>
    arrange(case_barcode, sample_barcode, file_name)
}

primary_diagnosis <- function(diagnoses) {
  if (length(diagnoses) == 0) return(list())
  is_primary <- vapply(diagnoses, function(x) isTRUE(x$diagnosis_is_primary_disease), logical(1))
  if (any(is_primary)) return(diagnoses[[which(is_primary)[1]]])
  class_primary <- vapply(diagnoses, function(x) identical(tolower(x$classification_of_tumor %||% ""), "primary"), logical(1))
  if (any(class_primary)) return(diagnoses[[which(class_primary)[1]]])
  diagnoses[[1]]
}

flatten_gdc_case_hit <- function(hit) {
  demo <- hit$demographic %||% list()
  diag <- primary_diagnosis(hit$diagnoses %||% list())
  tibble(
    case_id = hit$case_id %||% hit$id %||% NA_character_,
    case_barcode = hit$submitter_id %||% NA_character_,
    gdc_primary_site_case = hit$primary_site %||% NA_character_,
    gdc_disease_type_case = hit$disease_type %||% NA_character_,
    Gender_gdc = normalise_sex(demo$gender %||% NA_character_),
    Race_gdc = clean_missing(demo$race %||% NA_character_),
    Ethnicity_gdc = clean_missing(demo$ethnicity %||% NA_character_),
    vital_status_gdc = clean_missing(demo$vital_status %||% NA_character_),
    days_to_death_gdc = clean_numeric(demo$days_to_death %||% NA_character_),
    age_at_diagnosis_days_gdc = clean_numeric(diag$age_at_diagnosis %||% NA_character_),
    age_at_diagnosis_years_gdc = age_at_diagnosis_days_gdc / 365.25,
    year_of_diagnosis_gdc = clean_numeric(diag$year_of_diagnosis %||% NA_character_),
    primary_diagnosis_gdc = clean_missing(diag$primary_diagnosis %||% NA_character_),
    tumor_location_gdc = clean_missing(diag$tissue_or_organ_of_origin %||% NA_character_),
    tissue_or_organ_of_origin_gdc = clean_missing(diag$tissue_or_organ_of_origin %||% NA_character_),
    site_of_resection_or_biopsy_gdc = clean_missing(diag$site_of_resection_or_biopsy %||% NA_character_),
    icd_10_code_gdc = clean_missing(diag$icd_10_code %||% NA_character_),
    tumor_grade_gdc = clean_missing(diag$tumor_grade %||% NA_character_),
    ajcc_pathologic_stage_gdc = clean_missing(diag$ajcc_pathologic_stage %||% NA_character_),
    Stage_Simple_gdc = normalise_stage(diag$ajcc_pathologic_stage %||% NA_character_),
    ajcc_pathologic_t_gdc = clean_missing(diag$ajcc_pathologic_t %||% NA_character_),
    ajcc_pathologic_n_gdc = clean_missing(diag$ajcc_pathologic_n %||% NA_character_),
    ajcc_pathologic_m_gdc = clean_missing(diag$ajcc_pathologic_m %||% NA_character_),
    residual_disease_gdc = clean_missing(diag$residual_disease %||% NA_character_),
    prior_malignancy_gdc = clean_missing(diag$prior_malignancy %||% NA_character_),
    prior_treatment_gdc = clean_missing(diag$prior_treatment %||% NA_character_),
    days_to_last_follow_up_gdc = clean_numeric(diag$days_to_last_follow_up %||% NA_character_)
  )
}

query_gdc_cases <- function() {
  filters <- list(op = "=", content = list(field = "project.project_id", value = list(project_id)))
  fields <- paste(
    c(
      "case_id", "submitter_id", "primary_site", "disease_type",
      "demographic.gender", "demographic.race", "demographic.ethnicity",
      "demographic.vital_status", "demographic.days_to_death",
      "diagnoses.diagnosis_is_primary_disease", "diagnoses.classification_of_tumor",
      "diagnoses.tissue_or_organ_of_origin", "diagnoses.site_of_resection_or_biopsy",
      "diagnoses.primary_diagnosis", "diagnoses.icd_10_code", "diagnoses.tumor_grade",
      "diagnoses.ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_t",
      "diagnoses.ajcc_pathologic_n", "diagnoses.ajcc_pathologic_m",
      "diagnoses.residual_disease", "diagnoses.prior_malignancy", "diagnoses.prior_treatment",
      "diagnoses.age_at_diagnosis", "diagnoses.year_of_diagnosis",
      "diagnoses.days_to_last_follow_up"
    ),
    collapse = ","
  )
  body <- list(
    filters = filters,
    fields = fields,
    format = "JSON",
    size = 2000,
    expand = "demographic,diagnoses"
  )
  parsed <- api_post_json(gdc_cases_endpoint, body)
  hits <- parsed$data$hits
  if (length(hits) == 0) stop("GDC case query returned no TCGA-STAD cases.")
  bind_rows(lapply(hits, flatten_gdc_case_hit)) |>
    distinct(case_barcode, .keep_all = TRUE) |>
    arrange(case_barcode)
}

clinical_long_to_wide <- function(x, id_col) {
  if (nrow(x) == 0) return(tibble())
  x |>
    transmute(
      id = .data[[id_col]],
      clinical_attribute = clinicalAttributeId,
      value = value
    ) |>
    filter(!is.na(id), !is.na(clinical_attribute)) |>
    distinct(id, clinical_attribute, .keep_all = TRUE) |>
    pivot_wider(names_from = clinical_attribute, values_from = value) |>
    rename(!!id_col := id)
}

fetch_cbio_clinical <- function(clinical_type) {
  url <- paste0(cbio_base, "/studies/", cbio_study_id, "/clinical-data")
  parsed <- api_get_json(
    url,
    query = list(clinicalDataType = clinical_type, projection = "DETAILED")
  )
  if (length(parsed) == 0) return(tibble())
  flatten_list_columns(as_tibble(bind_rows(parsed)))
}

final_tpm_rds <- file.path(matrix_dir, "Auto_tcga_stad_tpm_matrix_gene_symbol.rds")
final_count_rds <- file.path(matrix_dir, "Auto_tcga_stad_unstranded_counts_matrix_gene_symbol.rds")

if (file.exists(final_tpm_rds) && file.size(final_tpm_rds) > 1000 &&
    file.exists(final_count_rds) && file.size(final_count_rds) > 1000) {
  log_msg("Final TCGA-STAD matrices already exist. Skipping download phase.")
} else {
  log_msg("Querying GDC STAR-count file metadata.")
  log_msg("Querying GDC STAR-count file metadata.")
  gdc_meta <- query_gdc_star_counts()
  if (!is.na(max_files)) {
    gdc_meta <- head(gdc_meta, max_files)
  }
  write.csv(gdc_meta, file.path(raw_dir, "Auto_gdc_star_counts_file_metadata.csv"), row.names = FALSE)
  write_tsv(
    gdc_meta |>
      transmute(id = file_id, filename = file_name, md5 = md5sum, size = file_size, state = file_state),
    file.path(raw_dir, "Auto_gdc_manifest.txt")
  )
  log_msg("GDC STAR-count files in manifest: ", nrow(gdc_meta))
  
  log_msg("Querying GDC case clinical metadata.")
  gdc_cases <- query_gdc_cases()
  write.csv(gdc_cases, file.path(raw_dir, "Auto_gdc_cases_clinical.csv"), row.names = FALSE)
  
  log_msg("Querying cBioPortal clinical metadata for ", cbio_study_id, ".")
  cbio_patient_long <- fetch_cbio_clinical("PATIENT")
  cbio_sample_long <- fetch_cbio_clinical("SAMPLE")
  
  if (nrow(cbio_patient_long) > 0) {
    write.csv(cbio_patient_long, file.path(raw_dir, "Auto_cbioportal_stad_tcga_gdc_patient_clinical_long.csv"), row.names = FALSE)
  }
  if (nrow(cbio_sample_long) > 0) {
    write.csv(cbio_sample_long, file.path(raw_dir, "Auto_cbioportal_stad_tcga_gdc_sample_clinical_long.csv"), row.names = FALSE)
  }
  
  cbio_patient_wide <- clinical_long_to_wide(cbio_patient_long, "patientId")
  cbio_sample_wide <- clinical_long_to_wide(cbio_sample_long, "sampleId")
  
  if (nrow(cbio_patient_wide) > 0) {
    colnames(cbio_patient_wide) <- c("case_barcode", paste0("cbio_patient_", clean_names(colnames(cbio_patient_wide)[-1])))
    write.csv(cbio_patient_wide, file.path(raw_dir, "Auto_cbioportal_stad_tcga_gdc_patient_clinical_wide.csv"), row.names = FALSE)
  }
  if (nrow(cbio_sample_wide) > 0) {
    colnames(cbio_sample_wide) <- c("cbio_sample_id", paste0("cbio_sample_", clean_names(colnames(cbio_sample_wide)[-1])))
    write.csv(cbio_sample_wide, file.path(raw_dir, "Auto_cbioportal_stad_tcga_gdc_sample_clinical_wide.csv"), row.names = FALSE)
  }
  
  cbio_patient_wide <- ensure_columns(
    cbio_patient_wide,
    c("case_barcode", "cbio_patient_PRIMARY_SITE_PATIENT", "cbio_patient_BIOPSY_SITE", "cbio_patient_SEX", "cbio_patient_PATH_STAGE")
  )
  cbio_sample_wide <- ensure_columns(
    cbio_sample_wide,
    c("cbio_sample_id", "cbio_sample_SAMPLE_TYPE")
  )
  
  download_one_gdc_file <- function(file_id, file_name, md5sum, file_size) {
    dest_dir <- file.path(gdc_file_dir, file_id)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(dest_dir, file_name)
  
    if (verify_download(dest, file_size, md5sum)) {
      return(tibble(file_id = file_id, path = dest, status = "exists_verified"))
    }
    if (file.exists(dest) && !overwrite_bad) {
      stop("Existing GDC file failed size/checksum validation and SCREF_STAD_OVERWRITE_BAD is not TRUE: ", dest)
    }
    if (skip_download) {
      stop("Missing or invalid GDC file while SCREF_STAD_SKIP_DOWNLOAD is TRUE: ", dest)
    }
  
    url <- paste0(gdc_data_endpoint, file_id)
    last_error <- NULL
    for (attempt in seq_len(4)) {
      log_msg("Downloading GDC file ", file_id, " (attempt ", attempt, "/4).")
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
  
  log_msg("Downloading/verifying GDC STAR-count files.")
  download_status <- vector("list", nrow(gdc_meta))
  for (i in seq_len(nrow(gdc_meta))) {
    download_status[[i]] <- download_one_gdc_file(
      file_id = gdc_meta$file_id[i],
      file_name = gdc_meta$file_name[i],
      md5sum = gdc_meta$md5sum[i],
      file_size = gdc_meta$file_size[i]
    )
    if (i %% 10 == 0 || i == nrow(gdc_meta)) {
      log_msg("Verified ", i, " / ", nrow(gdc_meta), " GDC files.")
    }
  }
  download_status <- bind_rows(download_status)
  write.csv(download_status, file.path(raw_dir, "Auto_gdc_download_status.csv"), row.names = FALSE)
  
  gdc_meta <- gdc_meta |>
    left_join(download_status |> select(file_id, path, download_status = status), by = "file_id")
  
  sample_metadata <- gdc_meta |>
    left_join(gdc_cases, by = c("case_id", "case_barcode")) |>
    left_join(cbio_patient_wide, by = "case_barcode") |>
    left_join(cbio_sample_wide, by = "cbio_sample_id") |>
    mutate(
      tumor_location = coalesce_clean(
        tissue_or_organ_of_origin_gdc,
        site_of_resection_or_biopsy_gdc,
        cbio_patient_PRIMARY_SITE_PATIENT,
        cbio_patient_BIOPSY_SITE
      ),
      Gender = coalesce_clean(Gender_gdc, cbio_patient_SEX),
      Stage = coalesce_clean(ajcc_pathologic_stage_gdc, cbio_patient_PATH_STAGE),
      Stage_Simple = normalise_stage(Stage),
      Grade = tumor_grade_gdc,
      primary_diagnosis = primary_diagnosis_gdc
    ) |>
    arrange(case_barcode, sample_barcode, file_name)
  
  metadata_csv <- file.path(metadata_dir, "Auto_tcga_stad_sample_metadata.csv")
  metadata_rds <- file.path(metadata_dir, "Auto_tcga_stad_sample_metadata.rds")
  write.csv(sample_metadata, metadata_csv, row.names = FALSE)
  saveRDS(sample_metadata, metadata_rds)
  
  read_one_star_file <- function(path) {
    x <- read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE)
    if (!"gene_id" %in% colnames(x)) {
      stop("Missing gene_id column in STAR-count file: ", path)
    }
    count_col <- intersect(c("unstranded", "count", "counts"), colnames(x))[1]
    tpm_col <- intersect(c("tpm_unstranded", "tpm"), colnames(x))[1]
    gene_symbol_col <- intersect(c("gene_name", "GeneSymbol", "gene"), colnames(x))[1]
    gene_type_col <- intersect(c("gene_type", "gene_biotype"), colnames(x))[1]
    if (is.na(count_col)) stop("No unstranded count column found in STAR-count file: ", path)
    if (is.na(tpm_col)) stop("No TPM column found in STAR-count file: ", path)
    if (is.na(gene_symbol_col)) gene_symbol_col <- "gene_id"
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
        unstranded_count = as.numeric(.data[[count_col]]),
        tpm = as.numeric(.data[[tpm_col]])
      )
  }
  
  if (anyDuplicated(sample_metadata$sample_barcode)) {
    duplicated_samples <- unique(sample_metadata$sample_barcode[duplicated(sample_metadata$sample_barcode)])
    stop("Duplicate GDC sample barcodes found; refusing to build ambiguous matrices: ", paste(duplicated_samples, collapse = ", "))
  }
  if (any(is.na(sample_metadata$path)) || any(!file.exists(sample_metadata$path))) {
    missing_paths <- sample_metadata$path[is.na(sample_metadata$path) | !file.exists(sample_metadata$path)]
    stop("Some downloaded GDC files are missing; first missing path: ", missing_paths[1])
  }
  
  log_msg("Reading STAR-count files and building expression matrices.")
  first_expr <- read_one_star_file(sample_metadata$path[1])
  gene_key <- first_expr |> select(gene_id, gene_symbol, gene_type)
  count_mat <- matrix(NA_real_, nrow = nrow(gene_key), ncol = nrow(sample_metadata))
  tpm_mat <- matrix(NA_real_, nrow = nrow(gene_key), ncol = nrow(sample_metadata))
  colnames(count_mat) <- sample_metadata$sample_barcode
  colnames(tpm_mat) <- sample_metadata$sample_barcode
  rownames(count_mat) <- gene_key$gene_id
  rownames(tpm_mat) <- gene_key$gene_id
  count_mat[, 1] <- first_expr$unstranded_count
  tpm_mat[, 1] <- first_expr$tpm
  
  if (nrow(sample_metadata) > 1) {
    for (i in 2:nrow(sample_metadata)) {
      if (i %% 25 == 0 || i == nrow(sample_metadata)) {
        log_msg("Reading STAR-count file ", i, " / ", nrow(sample_metadata), ".")
      }
      this_expr <- read_one_star_file(sample_metadata$path[i])
      if (!identical(this_expr$gene_id, gene_key$gene_id)) {
        idx <- match(gene_key$gene_id, this_expr$gene_id)
        if (any(is.na(idx))) {
          stop("Gene IDs in ", sample_metadata$path[i], " do not match the reference file.")
        }
        count_mat[, i] <- this_expr$unstranded_count[idx]
        tpm_mat[, i] <- this_expr$tpm[idx]
      } else {
        count_mat[, i] <- this_expr$unstranded_count
        tpm_mat[, i] <- this_expr$tpm
      }
    }
  }
  
  collapse_to_symbol <- function(mat, gene_key, value_name) {
    valid_gene <- !is.na(gene_key$gene_symbol) & nzchar(gene_key$gene_symbol)
    dt <- as.data.table(mat[valid_gene, , drop = FALSE])
    dt[, GeneSymbol := gene_key$gene_symbol[valid_gene]]
    setcolorder(dt, c("GeneSymbol", setdiff(colnames(dt), "GeneSymbol")))
    out <- dt[, lapply(.SD, sum, na.rm = TRUE), by = GeneSymbol, .SDcols = setdiff(colnames(dt), "GeneSymbol")]
    setorder(out, GeneSymbol)
    attr(out, "value_name") <- value_name
    out
  }
  
  count_symbol <- collapse_to_symbol(count_mat, gene_key, "unstranded_count")
  tpm_symbol <- collapse_to_symbol(tpm_mat, gene_key, "tpm")
  
  count_symbol_matrix <- as.matrix(count_symbol[, -1, with = FALSE])
  rownames(count_symbol_matrix) <- count_symbol$GeneSymbol
  tpm_symbol_matrix <- as.matrix(tpm_symbol[, -1, with = FALSE])
  rownames(tpm_symbol_matrix) <- tpm_symbol$GeneSymbol
  
  write.csv(gene_key, file.path(matrix_dir, "Auto_tcga_stad_gene_key.csv"), row.names = FALSE)
  saveRDS(count_symbol_matrix, file.path(matrix_dir, "Auto_tcga_stad_unstranded_counts_matrix_gene_symbol.rds"))
  saveRDS(tpm_symbol_matrix, file.path(matrix_dir, "Auto_tcga_stad_tpm_matrix_gene_symbol.rds"))
  fwrite(count_symbol, file.path(matrix_dir, "Auto_tcga_stad_unstranded_counts_matrix_gene_symbol.tsv"), sep = "\t")
  fwrite(tpm_symbol, file.path(matrix_dir, "Auto_tcga_stad_tpm_matrix_gene_symbol.tsv"), sep = "\t")
  
  sample_summary <- sample_metadata |>
    count(sample_type_code, sample_type, tissue_type, tumor_location, Gender, Stage_Simple, name = "n_files") |>
    arrange(sample_type_code, sample_type, tumor_location, Gender, Stage_Simple)
  write.csv(sample_summary, file.path(table_dir, "Auto_tcga_stad_download_sample_summary.csv"), row.names = FALSE)
  
  check_summary <- tibble(
    metric = c(
      "gdc_star_count_files",
      "downloaded_or_verified_files",
      "metadata_rows",
      "metadata_rows_with_tumor_location",
      "gene_symbol_rows",
      "samples_in_matrices",
      "run_start",
      "run_end"
    ),
    value = c(
      nrow(gdc_meta),
      sum(download_status$status %in% c("exists_verified", "downloaded_verified")),
      nrow(sample_metadata),
      sum(!is.na(sample_metadata$tumor_location)),
      nrow(tpm_symbol),
      ncol(tpm_symbol_matrix),
      format(run_start, "%Y-%m-%d %H:%M:%S"),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
  write.csv(check_summary, file.path(table_dir, "Auto_tcga_stad_download_check_summary.csv"), row.names = FALSE)
  
  writeLines(messages, file.path(log_dir, "Auto_tcga_stad_bulk_download_messages.log"))
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "Auto_tcga_stad_bulk_download_session_info.txt"))
  
  log_msg("TCGA-STAD bulk download complete.")
  log_msg("Metadata: ", metadata_csv)
  log_msg("TPM matrix: ", file.path(matrix_dir, "Auto_tcga_stad_tpm_matrix_gene_symbol.tsv"))
  log_msg("Count matrix: ", file.path(matrix_dir, "Auto_tcga_stad_unstranded_counts_matrix_gene_symbol.tsv"))
}

####################
# Phase 2: GSVA Analysis
####################

out_dir <- file.path(
  project_dir, "ref_outs", "Metaprogrammes_Results", "centred",
  "bulk_tcga_stad_location_gsva"
)

for (tier in c("intermediate", "tables", "figures", "logs")) {
  dir.create(file.path(out_dir, tier), recursive = TRUE, showWarnings = FALSE)
}

refined_mp_path <- file.path(project_dir, "ref_outs", "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_mp_genes.rds")
tpm_path <- file.path(matrix_dir, "Auto_tcga_stad_tpm_matrix_gene_symbol.rds")
metadata_path <- file.path(metadata_dir, "Auto_tcga_stad_sample_metadata.rds")

required_inputs <- c(refined_mp_path, tpm_path, metadata_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))
}

####################
# Final centred MP/state definitions copied from
# analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R.
####################
state_groups <- list(
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17")
)

mp_desc <- c(
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
  "MP17" = "Immune-interactive glandular progenitor"
)

state_plot_labels <- c(
  "Basal to intestinal metaplasia" = "Basal-to-intestinal MPs",
  "SMG to intestinal metaplasia" = "SMG-to-intestinal MPs"
)

min_samples_per_location <- 10L
cardia_location <- "Cardia, NOS"
distal_noncardia_locations <- c("Gastric antrum", "Body of stomach", "Fundus of stomach")
####################

make_location_palette <- function(locations) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999", "#E41A1C")
  setNames(base[seq_along(locations)], locations)
}

significance_label <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
}

compute_global_stats <- function(data, score_type) {
  stats <- lapply(split(data, data$feature), function(d) {
    location_summary <- d %>%
      group_by(tumor_location) %>%
      summarise(
        n_samples = n_distinct(sample_barcode),
        median_score = median(score, na.rm = TRUE),
        mean_score = mean(score, na.rm = TRUE),
        .groups = "drop"
      )
    test_result <- tryCatch(kruskal.test(score ~ tumor_location, data = d), error = function(e) NULL)
    data.frame(
      score_type = score_type,
      state = as.character(d$state[1]),
      feature = as.character(d$feature[1]),
      test = "Kruskal-Wallis",
      n_samples = n_distinct(d$sample_barcode),
      n_locations = n_distinct(d$tumor_location),
      p_value = if (is.null(test_result)) NA_real_ else test_result$p.value,
      location_summary = paste0(
        location_summary$tumor_location, " (n=", location_summary$n_samples,
        ", median=", formatC(location_summary$median_score, format = "f", digits = 3), ")",
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(stats) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"), significance = significance_label(p_adj))
}

compute_binary_stats <- function(data, score_type) {
  stats <- lapply(split(data, data$feature), function(d) {
    group_summary <- d %>%
      group_by(tumor_location) %>%
      summarise(
        n_samples = n_distinct(sample_barcode),
        median_score = median(score, na.rm = TRUE),
        mean_score = mean(score, na.rm = TRUE),
        .groups = "drop"
      )
    test_result <- tryCatch(
      suppressWarnings(wilcox.test(score ~ tumor_location, data = d, exact = FALSE)),
      error = function(e) NULL
    )
    proximal_median <- group_summary$median_score[group_summary$tumor_location == "Proximal cardia"]
    distal_median <- group_summary$median_score[group_summary$tumor_location == "Distal non-cardia"]
    data.frame(
      score_type = score_type,
      state = as.character(d$state[1]),
      feature = as.character(d$feature[1]),
      test = "Wilcoxon rank-sum",
      n_samples = n_distinct(d$sample_barcode),
      n_locations = n_distinct(d$tumor_location),
      p_value = if (is.null(test_result)) NA_real_ else test_result$p.value,
      median_difference_distal_minus_cardia = distal_median - proximal_median,
      location_summary = paste0(
        group_summary$tumor_location, " (n=", group_summary$n_samples,
        ", median=", formatC(group_summary$median_score, format = "f", digits = 3), ")",
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(stats) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"), significance = significance_label(p_adj))
}

compute_pairwise_stats <- function(data, score_type) {
  stats <- lapply(split(data, data$feature), function(d) {
    locations <- sort(unique(as.character(d$tumor_location)))
    comparisons <- combn(locations, 2, simplify = FALSE)
    bind_rows(lapply(comparisons, function(pair) {
      d_pair <- d %>% filter(tumor_location %in% pair)
      test_result <- tryCatch(
        suppressWarnings(wilcox.test(score ~ tumor_location, data = d_pair, exact = FALSE)),
        error = function(e) NULL
      )
      data.frame(
        score_type = score_type,
        state = as.character(d_pair$state[1]),
        feature = as.character(d_pair$feature[1]),
        location_1 = pair[1],
        location_2 = pair[2],
        n_location_1 = sum(d_pair$tumor_location == pair[1]),
        n_location_2 = sum(d_pair$tumor_location == pair[2]),
        p_value = if (is.null(test_result)) NA_real_ else test_result$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })
  bind_rows(stats) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"), significance = significance_label(p_adj))
}

plot_location_boxplot <- function(data, stats, title_text, y_label, feature_labels, palette, subtitle_text) {
  annotation <- data %>%
    group_by(feature) %>%
    summarise(y_pos = max(score, na.rm = TRUE), .groups = "drop") %>%
    left_join(stats %>% select(feature, significance), by = "feature")
  y_span <- diff(range(data$score, na.rm = TRUE))
  annotation$y_pos <- annotation$y_pos + max(0.03, 0.08 * y_span)

  ggplot(data, aes(x = feature, y = score, fill = tumor_location)) +
    geom_boxplot(
      position = position_dodge(width = 0.8), width = 0.68,
      outlier.shape = NA, alpha = 0.82, linewidth = 0.4, colour = "black"
    ) +
    geom_point(
      aes(colour = tumor_location),
      position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.8),
      size = 0.9, alpha = 0.6, show.legend = FALSE
    ) +
    geom_text(
      data = annotation %>% filter(significance != ""),
      aes(x = feature, y = y_pos, label = significance),
      inherit.aes = FALSE, fontface = "bold", size = 4
    ) +
    scale_fill_manual(values = palette) +
    scale_colour_manual(values = palette) +
    scale_x_discrete(labels = feature_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL, y = y_label, fill = "Tumour location"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "grey35"),
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
      axis.line.x = element_blank(),
      legend.position = "top",
      plot.margin = margin(10, 18, 10, 10)
    )
}

####################
# Read and harmonise TCGA-STAD primary-tumour bulk RNA-seq with metadata.
####################
expr_mat <- readRDS(tpm_path)
metadata <- readRDS(metadata_path)
mp_genes <- readRDS(refined_mp_path)

required_metadata <- c("sample_barcode", "sample_type_code", "tumor_location")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0) {
  stop("Bulk metadata missing required column(s): ", paste(missing_metadata, collapse = ", "))
}

primary_meta <- metadata %>%
  filter(sample_type_code == "01", !is.na(tumor_location), tumor_location != "") %>%
  distinct(sample_barcode, .keep_all = TRUE)
common_samples <- intersect(colnames(expr_mat), primary_meta$sample_barcode)
if (length(common_samples) < 20) {
  stop("Fewer than 20 primary-tumour samples overlap the TPM matrix and metadata.")
}
primary_meta <- primary_meta %>%
  filter(sample_barcode %in% common_samples) %>%
  arrange(match(sample_barcode, common_samples))
expr_mat <- expr_mat[, primary_meta$sample_barcode, drop = FALSE]

rownames(expr_mat) <- toupper(trimws(rownames(expr_mat)))
valid_genes <- !is.na(rownames(expr_mat)) & nzchar(rownames(expr_mat))
expr_mat <- expr_mat[valid_genes, , drop = FALSE]
expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
expr_mat <- log2(expr_mat + 1)

####################
# Score TCGA EAC samples for comparison (prepare matrix)
####################
tcga_recon_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction"
tcga_meta_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_meta.rds")
tcga_matrix_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_tpm_matrix.rds")

if (!file.exists(tcga_meta_path) || !file.exists(tcga_matrix_path)) {
  stop("TCGA ESCA reconstruction missing in live storage.")
}
meta_tcga <- readRDS(tcga_meta_path)
tpm_tcga <- readRDS(tcga_matrix_path)

infer_histology <- function(type_vec, detailed_vec = NA_character_) {
  t <- tolower(paste(as.character(type_vec), as.character(detailed_vec)))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}
detailed_vec <- if ("Cancer_Type_Detailed" %in% colnames(meta_tcga)) meta_tcga$Cancer_Type_Detailed else NA_character_
type_vec <- if ("type" %in% colnames(meta_tcga)) meta_tcga$type else detailed_vec
meta_tcga$HistologyGroup <- infer_histology(type_vec, detailed_vec)

eac_samples <- meta_tcga %>%
  filter(sample_type_code == "01", HistologyGroup == "EAC") %>%
  pull(sample_barcode)
common_eac <- intersect(eac_samples, colnames(tpm_tcga))
if (length(common_eac) < 5) stop("Too few EAC samples found.")
tpm_eac <- tpm_tcga[, common_eac, drop = FALSE]

expr_eac <- log2(tpm_eac + 1)
expr_eac[!is.finite(expr_eac)] <- 0
rownames(expr_eac) <- toupper(trimws(rownames(expr_eac)))
valid_genes <- !is.na(rownames(expr_eac)) & nzchar(rownames(expr_eac))
expr_eac <- expr_eac[valid_genes, , drop = FALSE]
expr_eac <- rowsum(expr_eac, group = rownames(expr_eac), reorder = FALSE)

####################
# Combine STAD and EAC matrices for pooled GSVA
####################
# Intersect genes to ensure the background ranking environment is identical
common_genes <- intersect(rownames(expr_mat), rownames(expr_eac))
expr_combined <- cbind(expr_mat[common_genes, , drop=FALSE], expr_eac[common_genes, , drop=FALSE])

requested_mps <- unlist(state_groups, use.names = FALSE)
missing_mps <- setdiff(requested_mps, names(mp_genes))
if (length(missing_mps) > 0) {
  stop("Final centred MP gene list is missing requested MP(s): ", paste(missing_mps, collapse = ", "))
}

mp_gene_sets <- lapply(mp_genes[requested_mps], function(genes) {
  intersect(rownames(expr_combined), toupper(unique(as.character(genes))))
})
state_gene_sets <- lapply(state_groups, function(mps) unique(unlist(mp_gene_sets[mps], use.names = FALSE)))
all_gene_sets <- c(mp_gene_sets, state_gene_sets)
gene_set_sizes <- data.frame(
  score_type = c(rep("MP", length(mp_gene_sets)), rep("State", length(state_gene_sets))),
  feature = names(all_gene_sets),
  genes_in_refined_list = c(
    lengths(mp_genes[requested_mps]),
    lengths(lapply(state_groups, function(mps) unique(unlist(mp_genes[mps], use.names = FALSE))))
  ),
  genes_matched_in_bulk = lengths(all_gene_sets),
  stringsAsFactors = FALSE
)
if (any(gene_set_sizes$genes_matched_in_bulk < 5)) {
  bad_sets <- gene_set_sizes$feature[gene_set_sizes$genes_matched_in_bulk < 5]
  stop("Too few bulk-matched genes for: ", paste(bad_sets, collapse = ", "))
}
write.csv(gene_set_sizes, file.path(out_dir, "tables", "Auto_tcga_stad_centred_gene_set_sizes.csv"), row.names = FALSE)
####################

gsva_scores_combined <- GSVA::gsva(expr_combined, all_gene_sets, method = "gsva", kcdf = "Gaussian")
saveRDS(gsva_scores_combined, file.path(out_dir, "intermediate", "Auto_tcga_stad_centred_mp_state_gsva_scores.rds"))

score_df_all <- as.data.frame(t(gsva_scores_combined)) %>%
  mutate(sample_barcode = rownames(.))

# Split out the EAC scores
score_df_eac <- score_df_all %>%
  filter(sample_barcode %in% common_eac) %>%
  mutate(tumor_location = "TCGA EAC")

gsva_scores <- gsva_scores_combined[, colnames(expr_mat), drop=FALSE]

####################

location_counts_all <- primary_meta %>% count(tumor_location, name = "n_primary_tumours") %>% arrange(desc(n_primary_tumours), tumor_location)
eligible_locations <- location_counts_all %>%
  filter(n_primary_tumours >= min_samples_per_location) %>%
  pull(tumor_location)
excluded_locations <- location_counts_all %>% filter(!tumor_location %in% eligible_locations)
write.csv(location_counts_all, file.path(out_dir, "tables", "Auto_tcga_stad_primary_tumour_location_counts.csv"), row.names = FALSE)
write.csv(excluded_locations, file.path(out_dir, "tables", "Auto_tcga_stad_excluded_small_location_groups.csv"), row.names = FALSE)

score_df <- as.data.frame(t(gsva_scores)) %>%
  mutate(sample_barcode = rownames(.)) %>%
  left_join(primary_meta %>% select(sample_barcode, tumor_location), by = "sample_barcode") %>%
  filter(tumor_location %in% eligible_locations) %>%
  bind_rows(score_df_eac)

eligible_locations_combined <- c(eligible_locations, "TCGA EAC")

binary_score_df <- as.data.frame(t(gsva_scores)) %>%
  mutate(sample_barcode = rownames(.)) %>%
  left_join(primary_meta %>% select(sample_barcode, tumor_location), by = "sample_barcode") %>%
  mutate(
    tumor_location = case_when(
      tumor_location == cardia_location ~ "Proximal cardia",
      tumor_location %in% distal_noncardia_locations ~ "Distal non-cardia",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(tumor_location)) %>%
  bind_rows(score_df_eac) %>%
  mutate(tumor_location = factor(tumor_location, levels = c("Proximal cardia", "Distal non-cardia", "TCGA EAC")))

binary_excluded_locations <- primary_meta %>%
  filter(!tumor_location %in% c(cardia_location, distal_noncardia_locations)) %>%
  count(tumor_location, name = "n_primary_tumours") %>%
  arrange(desc(n_primary_tumours), tumor_location)
write.csv(binary_excluded_locations, file.path(out_dir, "tables", "Auto_tcga_stad_cardia_vs_distal_excluded_locations.csv"), row.names = FALSE)

mp_long <- score_df %>%
  pivot_longer(cols = all_of(requested_mps), names_to = "feature", values_to = "score") %>%
  mutate(
    state = vapply(feature, function(mp) names(state_groups)[vapply(state_groups, function(x) mp %in% x, logical(1))][1], character(1)),
    feature = factor(feature, levels = requested_mps),
    tumor_location = factor(tumor_location, levels = eligible_locations_combined)
  )
state_long <- score_df %>%
  pivot_longer(cols = all_of(names(state_groups)), names_to = "feature", values_to = "score") %>%
  mutate(
    state = feature,
    feature = factor(feature, levels = names(state_groups)),
    tumor_location = factor(tumor_location, levels = eligible_locations_combined)
  )

mp_long_cardia_distal <- binary_score_df %>%
  pivot_longer(cols = all_of(requested_mps), names_to = "feature", values_to = "score") %>%
  mutate(
    state = vapply(feature, function(mp) names(state_groups)[vapply(state_groups, function(x) mp %in% x, logical(1))][1], character(1)),
    feature = factor(feature, levels = requested_mps)
  )
state_long_cardia_distal <- binary_score_df %>%
  pivot_longer(cols = all_of(names(state_groups)), names_to = "feature", values_to = "score") %>%
  mutate(
    state = feature,
    feature = factor(feature, levels = names(state_groups))
  )

mp_global_stats <- compute_global_stats(mp_long, "MP GSVA")
mp_pairwise_stats <- compute_pairwise_stats(mp_long, "MP GSVA")
state_global_stats <- compute_global_stats(state_long, "State-union GSVA")
state_pairwise_stats <- compute_pairwise_stats(state_long, "State-union GSVA")
mp_cardia_distal_stats <- compute_binary_stats(mp_long_cardia_distal, "MP GSVA")
state_cardia_distal_stats <- compute_binary_stats(state_long_cardia_distal, "State-union GSVA")

write.csv(mp_long, file.path(out_dir, "tables", "Auto_tcga_stad_centred_mp_gsva_sample_scores.csv"), row.names = FALSE)
write.csv(state_long, file.path(out_dir, "tables", "Auto_tcga_stad_centred_state_union_gsva_sample_scores.csv"), row.names = FALSE)
write.csv(mp_global_stats, file.path(out_dir, "tables", "Auto_tcga_stad_centred_mp_gsva_global_stats.csv"), row.names = FALSE)
write.csv(mp_pairwise_stats, file.path(out_dir, "tables", "Auto_tcga_stad_centred_mp_gsva_pairwise_stats.csv"), row.names = FALSE)
write.csv(state_global_stats, file.path(out_dir, "tables", "Auto_tcga_stad_centred_state_union_gsva_global_stats.csv"), row.names = FALSE)
write.csv(state_pairwise_stats, file.path(out_dir, "tables", "Auto_tcga_stad_centred_state_union_gsva_pairwise_stats.csv"), row.names = FALSE)
write.csv(mp_long_cardia_distal, file.path(out_dir, "tables", "Auto_tcga_stad_cardia_vs_distal_mp_gsva_sample_scores.csv"), row.names = FALSE)
write.csv(state_long_cardia_distal, file.path(out_dir, "tables", "Auto_tcga_stad_cardia_vs_distal_state_union_gsva_sample_scores.csv"), row.names = FALSE)
write.csv(mp_cardia_distal_stats, file.path(out_dir, "tables", "Auto_tcga_stad_cardia_vs_distal_mp_gsva_stats.csv"), row.names = FALSE)
write.csv(state_cardia_distal_stats, file.path(out_dir, "tables", "Auto_tcga_stad_cardia_vs_distal_state_union_gsva_stats.csv"), row.names = FALSE)

location_palette <- make_location_palette(eligible_locations_combined)
mp_pdf <- file.path(out_dir, "figures", "Auto_tcga_stad_centred_mp_gsva_by_location_boxplots.pdf")
pdf(mp_pdf, width = 17, height = 9, useDingbats = FALSE)
for (state_name in names(state_groups)) {
  state_mps <- state_groups[[state_name]]
  labels <- setNames(paste0(state_mps, "\n", mp_desc[state_mps]), state_mps)
  print(plot_location_boxplot(
    mp_long %>% filter(state == state_name),
    mp_global_stats %>% filter(state == state_name),
    paste0(state_plot_labels[[state_name]], " by stomach site"),
    "GSVA enrichment score", labels, location_palette,
    "Primary tumours; BH-adjusted Kruskal-Wallis p-values."
  ))
}
dev.off()

state_labels <- setNames(
  c("Basal to intestinal\nmetaplasia", "SMG to intestinal\nmetaplasia"),
  names(state_groups)
)
state_pdf <- file.path(out_dir, "figures", "Auto_tcga_stad_centred_state_union_gsva_by_location_boxplots.pdf")
pdf(state_pdf, width = 10, height = 8, useDingbats = FALSE)
print(plot_location_boxplot(
  state_long, state_global_stats,
  "State score by stomach site",
  "GSVA enrichment score", state_labels, location_palette,
  "Primary tumours; BH-adjusted Kruskal-Wallis p-values."
))
dev.off()

binary_palette <- c("Proximal cardia" = "#D55E00", "Distal non-cardia" = "#0072B2", "TCGA EAC" = "#E41A1C")
binary_mp_pdf <- file.path(out_dir, "figures", "Auto_tcga_stad_cardia_vs_distal_mp_gsva_boxplots.pdf")
pdf(binary_mp_pdf, width = 15, height = 9, useDingbats = FALSE)
for (state_name in names(state_groups)) {
  state_mps <- state_groups[[state_name]]
  labels <- setNames(paste0(state_mps, "\n", mp_desc[state_mps]), state_mps)
  print(plot_location_boxplot(
    mp_long_cardia_distal %>% filter(state == state_name),
    mp_cardia_distal_stats %>% filter(state == state_name),
    paste0(state_plot_labels[[state_name]], ": cardia vs distal"),
    "GSVA enrichment score", labels, binary_palette,
    "Primary tumours; BH-adjusted Wilcoxon p-values."
  ))
}
dev.off()

binary_state_pdf <- file.path(out_dir, "figures", "Auto_tcga_stad_cardia_vs_distal_state_union_gsva_boxplots.pdf")
pdf(binary_state_pdf, width = 10, height = 8, useDingbats = FALSE)
print(plot_location_boxplot(
  state_long_cardia_distal, state_cardia_distal_stats,
  "State score: cardia vs distal",
  "GSVA enrichment score", state_labels, binary_palette,
  "Primary tumours; BH-adjusted Wilcoxon p-values."
))
dev.off()

summary_table <- bind_rows(
  mp_global_stats %>% select(score_type, state, feature, n_samples, n_locations, p_value, p_adj, significance, location_summary),
  state_global_stats %>% select(score_type, state, feature, n_samples, n_locations, p_value, p_adj, significance, location_summary)
)
write.csv(summary_table, file.path(out_dir, "tables", "Auto_tcga_stad_centred_location_gsva_summary.csv"), row.names = FALSE)

run_summary <- data.frame(
  metric = c("primary_tumours_in_metadata", "primary_tumours_scored", "eligible_locations", "excluded_location_groups", "cardia_vs_distal_samples", "cardia_vs_distal_excluded_samples", "MP_features", "state_features"),
  value = c(nrow(primary_meta), nrow(score_df), length(eligible_locations), nrow(excluded_locations), nrow(binary_score_df), sum(binary_excluded_locations$n_primary_tumours), length(requested_mps), length(state_groups)),
  stringsAsFactors = FALSE
)
write.csv(run_summary, file.path(out_dir, "logs", "Auto_tcga_stad_centred_location_gsva_run_summary.csv"), row.names = FALSE)
message("Completed TCGA-STAD centred MP/state GSVA location analysis.")
