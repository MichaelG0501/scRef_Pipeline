####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/geo_survival_data_prep.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_geo_survival_data_prep.R
#
# Download and prepare GEO bulk-expression datasets for external survival analysis.
# Current datasets:
#   - GSE19417 (public survival metadata available in GEO)
#   - GSE13898 (public clinicopathology only; no public survival metadata in GEO)
#
# Inputs:
#   GEO series matrix, clinical supplement, and platform annotation files
#
# Outputs:
#   ref_outs/geo_survival/raw/*
#   ref_outs/geo_survival/Auto_geo_survival_dataset_manifest.csv
#   ref_outs/geo_survival/Auto_geo_survival_dataset_manifest.rds
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_meta.csv
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#   ref_outs/geo_survival/Auto_GSE19417_probe_map.csv
#   ref_outs/geo_survival/Auto_GSE13898_meta.rds
#   ref_outs/geo_survival/Auto_GSE13898_meta.csv
#   ref_outs/geo_survival/Auto_GSE13898_expr_gene.rds
#   ref_outs/geo_survival/Auto_GSE13898_probe_map.csv
####################

library(data.table)
library(dplyr)
library(matrixStats)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- "geo_survival"
raw_dir <- file.path(out_dir, "raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) {
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }
  dest
}

clean_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\([^)]*\\)", "", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

parse_geo_line <- function(line) {
  vals <- strsplit(line, "\t", fixed = TRUE)[[1]]
  vals <- gsub('^"|"$', "", vals)
  vals
}

derive_series_key <- function(raw_key, values) {
  base_key <- sub("^!Sample_", "", raw_key)
  if (identical(base_key, "title")) {
    return("sample_title")
  }
  if (grepl("^characteristics_", base_key)) {
    first_value <- values[which(!is.na(values) & nzchar(values))[1]]
    if (length(first_value) == 1 && !is.na(first_value) && grepl(":", first_value, fixed = TRUE)) {
      base_key <- sub(":.*$", "", first_value)
    }
  }
  clean_name(base_key)
}

read_series_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  sample_lines <- character(0)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || identical(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_")) sample_lines <- c(sample_lines, line)
  }

  parsed <- lapply(sample_lines, parse_geo_line)
  geo_idx <- which(vapply(parsed, function(x) identical(x[1], "!Sample_geo_accession"), logical(1)))
  if (length(geo_idx) == 0) stop("Could not find !Sample_geo_accession in: ", path)

  sample_ids <- parsed[[geo_idx[1]]][-1]
  meta <- data.frame(sample_geo_accession = sample_ids, stringsAsFactors = FALSE, check.names = FALSE)
  key_counts <- integer(0)

  for (vals in parsed) {
    raw_key <- vals[1]
    if (identical(raw_key, "!Sample_geo_accession")) next
    key <- derive_series_key(raw_key, vals[-1])
    if (!key %in% names(key_counts)) key_counts[key] <- 0L
    key_counts[key] <- key_counts[key] + 1L
    if (key_counts[key] > 1L) key <- paste0(key, "_", key_counts[key])
    values <- vals[-1]
    length(values) <- length(sample_ids)
    meta[[key]] <- values
  }

  meta
}

read_geo_table <- function(path, start_marker, end_marker) {
  cmd <- paste0(
    "gzip -dc ", shQuote(path),
    " | awk '/^", start_marker, "$/{flag=1;next}/^", end_marker, "$/{flag=0} flag'"
  )
  fread(cmd = cmd, data.table = FALSE, quote = "\"", check.names = FALSE)
}

clean_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "---")] <- NA_character_
  x <- sub("///.*$", "", x)
  x <- sub("//.*$", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "---")] <- NA_character_
  x
}

strip_characteristic_prefix <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[^:]+:\\s*", "", x)
  x[x %in% c("", "NA", "---")] <- NA_character_
  x
}

collapse_expression_by_gene <- function(expr_df, annot_df) {
  probe_ids <- expr_df[[1]]
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- probe_ids

  annot_map <- annot_df %>%
    transmute(
      probe_id = ID,
      gene_symbol = clean_gene_symbol(`Gene symbol`)
    )

  gene_symbol <- annot_map$gene_symbol[match(rownames(expr_mat), annot_map$probe_id)]
  keep <- !is.na(gene_symbol) & nzchar(gene_symbol)
  expr_mat <- expr_mat[keep, , drop = FALSE]
  gene_symbol <- gene_symbol[keep]

  probe_vars <- matrixStats::rowVars(expr_mat, na.rm = TRUE)
  probe_means <- rowMeans(expr_mat, na.rm = TRUE)
  probe_vars[is.na(probe_vars)] <- -Inf
  probe_means[is.na(probe_means)] <- -Inf

  ord <- order(gene_symbol, -probe_vars, -probe_means, rownames(expr_mat))
  expr_mat <- expr_mat[ord, , drop = FALSE]
  gene_symbol <- gene_symbol[ord]
  probe_vars <- probe_vars[ord]
  probe_means <- probe_means[ord]

  keep_idx <- !duplicated(gene_symbol)
  collapsed <- expr_mat[keep_idx, , drop = FALSE]
  rownames(collapsed) <- gene_symbol[keep_idx]

  probe_counts <- table(gene_symbol)
  probe_map <- data.frame(
    gene_symbol = gene_symbol[keep_idx],
    probe_id = rownames(expr_mat)[keep_idx],
    row_variance = probe_vars[keep_idx],
    row_mean = probe_means[keep_idx],
    n_probes = as.integer(probe_counts[gene_symbol[keep_idx]]),
    stringsAsFactors = FALSE
  )

  list(expr_mat = collapsed, probe_map = probe_map)
}

infer_histology <- function(x) {
  y <- tolower(as.character(x))
  out <- rep("Other", length(y))
  out[grepl("oesophageal adenocarcinoma|esophageal adenocarcinoma", y)] <- "EAC"
  out[grepl("squamous", y)] <- "ESCC"
  out[grepl("gastric", y)] <- "Gastric"
  out
}

save_dataset_outputs <- function(dataset_id, meta_df, expr_obj) {
  meta_rds <- file.path(out_dir, paste0("Auto_", dataset_id, "_meta.rds"))
  meta_csv <- file.path(out_dir, paste0("Auto_", dataset_id, "_meta.csv"))
  expr_rds <- file.path(out_dir, paste0("Auto_", dataset_id, "_expr_gene.rds"))
  probe_csv <- file.path(out_dir, paste0("Auto_", dataset_id, "_probe_map.csv"))

  saveRDS(meta_df, meta_rds)
  write.csv(meta_df, meta_csv, row.names = FALSE)
  saveRDS(expr_obj$expr_mat, expr_rds)
  write.csv(expr_obj$probe_map, probe_csv, row.names = FALSE)

  list(
    meta_rds = meta_rds,
    meta_csv = meta_csv,
    expr_rds = expr_rds,
    probe_csv = probe_csv
  )
}

prepare_gse19417 <- function() {
  matrix_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19417/matrix/GSE19417_series_matrix.txt.gz",
    file.path(raw_dir, "GSE19417_series_matrix.txt.gz")
  )
  annot_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/annot/GPL4372.annot.gz",
    file.path(raw_dir, "GPL4372.annot.gz")
  )
  clinical_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19417/suppl/GSE19417_De-identified_Data_76-patients_13Oct2010.txt.gz",
    file.path(raw_dir, "GSE19417_De-identified_Data_76-patients_13Oct2010.txt.gz")
  )

  series_meta <- read_series_header(matrix_path)
  if (!"patient_id" %in% colnames(series_meta)) {
    series_meta$patient_id <- sub("^Esophageal Adenocarcinoma\\s+", "", series_meta$sample_title)
  }
  if ("patient_id" %in% colnames(series_meta)) {
    series_meta$patient_id <- strip_characteristic_prefix(series_meta$patient_id)
  }
  if ("gender" %in% colnames(series_meta)) {
    series_meta$gender <- strip_characteristic_prefix(series_meta$gender)
  }
  if ("survival_days" %in% colnames(series_meta)) {
    series_meta$survival_days <- suppressWarnings(as.numeric(strip_characteristic_prefix(series_meta$survival_days)))
  }
  if ("tumor_histology" %in% colnames(series_meta)) {
    series_meta$tumor_histology <- strip_characteristic_prefix(series_meta$tumor_histology)
  }
  if ("tumor_differentiation" %in% colnames(series_meta)) {
    series_meta$tumor_differentiation <- strip_characteristic_prefix(series_meta$tumor_differentiation)
  }
  if ("positive_nodes" %in% colnames(series_meta)) {
    series_meta$positive_nodes <- suppressWarnings(as.numeric(strip_characteristic_prefix(series_meta$positive_nodes)))
  }
  series_meta$dataset <- "GSE19417"

  expr_df <- read_geo_table(matrix_path, "!series_matrix_table_begin", "!series_matrix_table_end")
  annot_df <- read_geo_table(annot_path, "!platform_table_begin", "!platform_table_end")
  expr_obj <- collapse_expression_by_gene(expr_df, annot_df)

  clinical_raw <- read.delim(
    gzfile(clinical_path),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "n/a", "N/A")
  )

  clinical_df <- data.frame(
    patient_id = clinical_raw[["Patient_ID"]],
    gender = clinical_raw[["Gender"]],
    survival_days = suppressWarnings(as.numeric(clinical_raw[["SurvivalDays"]])),
    survival_censoring = clinical_raw[["SurvivalCensoring_Kaplan-Meier"]],
    positive_nodes = suppressWarnings(as.numeric(clinical_raw[["PositiveNodes"]])),
    tumor_histology = clinical_raw[["TumorHistology"]],
    tumor_differentiation = clinical_raw[["TumorDifferentiation"]],
    tumour_length_cm = suppressWarnings(as.numeric(clinical_raw[["TumourLength_cm"]])),
    tumour_width_cm = suppressWarnings(as.numeric(clinical_raw[["TumourWidth_cm"]])),
    tumor_volume_mm3 = suppressWarnings(as.numeric(clinical_raw[["TumorVolume_mm3"]])),
    stringsAsFactors = FALSE
  )

  meta_df <- series_meta %>%
    left_join(clinical_df, by = "patient_id", suffix = c("_series", "")) %>%
    mutate(
      gender = coalesce(gender, gender_series),
      survival_days_series = suppressWarnings(as.numeric(survival_days_series)),
      positive_nodes = coalesce(positive_nodes, positive_nodes_series),
      tumor_histology = coalesce(tumor_histology, tumor_histology_series),
      tumor_differentiation = coalesce(tumor_differentiation, tumor_differentiation_series),
      OS_time = coalesce(survival_days, survival_days_series),
      OS_event = ifelse(
        is.na(OS_time),
        NA_real_,
        ifelse(
          !is.na(survival_censoring) & tolower(trimws(survival_censoring)) == "rightcensored",
          0,
          1
        )
      ),
      HistologyGroup = infer_histology(tumor_histology),
      survival_public = TRUE,
      analysis_ready_for_survival = HistologyGroup == "EAC" & !is.na(OS_time) & !is.na(OS_event)
    )

  out_paths <- save_dataset_outputs("GSE19417", meta_df, expr_obj)

  data.frame(
    dataset = "GSE19417",
    pubmed_id = "20621683",
    platform = "GPL4372",
    title = "Human esophageal adenocarcinomas",
    survival_public = TRUE,
    analysis_ready_for_survival = TRUE,
    n_profiled_samples = ncol(expr_obj$expr_mat),
    n_gene_symbols = nrow(expr_obj$expr_mat),
    n_public_metadata_rows = nrow(meta_df),
    n_survival_rows = sum(!is.na(meta_df$OS_time)),
    n_eac_survival_rows = sum(meta_df$analysis_ready_for_survival),
    meta_rds = out_paths$meta_rds,
    expr_rds = out_paths$expr_rds,
    notes = "Use EAC-only rows with public GEO survival metadata; GSVA can be scored on full cohort and EAC-only subsets.",
    stringsAsFactors = FALSE
  )
}

prepare_gse13898 <- function() {
  matrix_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13898/matrix/GSE13898_series_matrix.txt.gz",
    file.path(raw_dir, "GSE13898_series_matrix.txt.gz")
  )
  annot_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6102/annot/GPL6102.annot.gz",
    file.path(raw_dir, "GPL6102.annot.gz")
  )
  clinical_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13898/suppl/GSE13898_EAC_MDACC_DepSysB_ClinicalInformation_update.txt.gz",
    file.path(raw_dir, "GSE13898_EAC_MDACC_DepSysB_ClinicalInformation_update.txt.gz")
  )

  series_meta <- read_series_header(matrix_path)
  if (!"pathology" %in% colnames(series_meta)) {
    series_meta$pathology <- NA_character_
  }
  series_meta$pathology <- strip_characteristic_prefix(series_meta$pathology)
  series_meta$array_id <- sub("-[0-9]+$", "", series_meta$sample_title)
  series_meta$dataset <- "GSE13898"
  series_meta$HistologyGroup <- ifelse(grepl("adenocarcinoma", tolower(series_meta$pathology)), "EAC", "Other")
  series_meta$survival_public <- FALSE
  series_meta$analysis_ready_for_survival <- FALSE
  series_meta$OS_time <- NA_real_
  series_meta$OS_event <- NA_real_

  expr_df <- read_geo_table(matrix_path, "!series_matrix_table_begin", "!series_matrix_table_end")
  annot_df <- read_geo_table(annot_path, "!platform_table_begin", "!platform_table_end")
  expr_obj <- collapse_expression_by_gene(expr_df, annot_df)

  clinical_raw <- read.delim(
    gzfile(clinical_path),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "n/a", "N/A")
  )

  clinical_df <- data.frame(
    array_id = clinical_raw[["Array ID"]],
    age = suppressWarnings(as.numeric(clinical_raw[["Age"]])),
    race = clinical_raw[["Race"]],
    gender_public = clinical_raw[["Gender"]],
    t_stage = clinical_raw[["T"]],
    n_stage = clinical_raw[["N"]],
    m_stage = clinical_raw[["M"]],
    tumor_grade = trimws(clinical_raw[["Tumor Grade"]]),
    group_public = trimws(clinical_raw[["Group"]]),
    tp63 = clinical_raw[["TP63"]],
    cluster_public = clinical_raw[["Cluster"]],
    surgical_stage = clinical_raw[["surgical STAGE"]],
    stringsAsFactors = FALSE
  )

  meta_df <- series_meta %>%
    left_join(clinical_df, by = "array_id")

  out_paths <- save_dataset_outputs("GSE13898", meta_df, expr_obj)

  data.frame(
    dataset = "GSE13898",
    pubmed_id = "21152079",
    platform = "GPL6102",
    title = "Robust prognostic biomarkers for EAC identified by systems-level characterization of tumor transcriptome",
    survival_public = FALSE,
    analysis_ready_for_survival = FALSE,
    n_profiled_samples = ncol(expr_obj$expr_mat),
    n_gene_symbols = nrow(expr_obj$expr_mat),
    n_public_metadata_rows = nrow(meta_df),
    n_survival_rows = 0,
    n_eac_survival_rows = 0,
    meta_rds = out_paths$meta_rds,
    expr_rds = out_paths$expr_rds,
    notes = "GEO public supplement contains clinicopathology only; no public per-sample survival metadata in GEO.",
    stringsAsFactors = FALSE
  )
}

manifest_df <- bind_rows(
  prepare_gse19417(),
  prepare_gse13898()
)

write.csv(
  manifest_df,
  file.path(out_dir, "Auto_geo_survival_dataset_manifest.csv"),
  row.names = FALSE
)
saveRDS(
  manifest_df,
  file.path(out_dir, "Auto_geo_survival_dataset_manifest.rds")
)

message("Prepared GEO survival inputs:")
print(manifest_df[, c("dataset", "survival_public", "analysis_ready_for_survival", "n_profiled_samples", "n_gene_symbols", "n_eac_survival_rows")])
