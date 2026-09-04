####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/Auto_drug_reversal/legacy_Auto_drug_reversal_clue_fallback.R
#   Methodology: not required (superseded exploratory drug-reversal branch)
#   Map: analysis/ANALYSIS_MAP.md
####################

####################
# Auto_drug_reversal_clue_fallback.R
#
# Direct CLUE/CMap L1000 query fallback for top 150 up/down scRef signatures.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(tibble)
})

####################
# setup
####################

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

base_dir <- "Auto_drug_reversal"
input_dir <- file.path(base_dir, "clue_inputs")
run_label <- Sys.getenv("AUTO_CLUE_RUN_LABEL", "")
out_dir <- file.path(base_dir, "clue_fallback")
if (nzchar(run_label)) out_dir <- file.path(out_dir, run_label)
raw_dir <- file.path(out_dir, "raw_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

query_name <- Sys.getenv("AUTO_CLUE_QUERY_NAME", paste0("Auto_scRef_drug_reversal_", format(Sys.Date(), "%Y%m%d")))
tool_id <- Sys.getenv("AUTO_CLUE_TOOL_ID", "sig_fastgutc_tool")
poll_seconds <- as.integer(Sys.getenv("AUTO_CLUE_POLL_SECONDS", "60"))
max_polls <- as.integer(Sys.getenv("AUTO_CLUE_MAX_POLLS", "60"))
max_genes <- as.integer(Sys.getenv("AUTO_CLUE_MAX_GENES", "75"))
keep_states <- Sys.getenv("AUTO_CLUE_KEEP_STATES", "")
keep_states <- if (nzchar(keep_states)) strsplit(keep_states, ",", fixed = TRUE)[[1]] else character()
existing_job_id <- Sys.getenv("AUTO_CLUE_JOB_ID", "")

####################
# helpers
####################

write_status <- function(status, detail, job_id = NA_character_, download_url = NA_character_) {
  fwrite(
    data.frame(
      step = "clue_fallback",
      status = status,
      detail = detail,
      job_id = job_id,
      download_url = download_url,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_clue_fallback_status.csv")
  )
}

read_body_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

redact_response <- function(txt) {
  txt <- gsub('"api_key":"[^"]+"', '"api_key":"REDACTED"', txt)
  api_key <- Sys.getenv("CLUE_API_KEY", Sys.getenv("CLUE_KEY", ""))
  if (nzchar(api_key)) txt <- gsub(api_key, "REDACTED", txt, fixed = TRUE)
  txt
}

limit_gmt <- function(infile, outfile, max_genes, keep_states = character()) {
  lines <- readLines(infile, warn = FALSE)
  if (length(keep_states) > 0) {
    lines <- lines[vapply(strsplit(lines, "\t", fixed = TRUE), function(x) x[1] %in% keep_states, logical(1))]
  }
  limited <- vapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) <= 2) return(line)
    genes <- fields[-seq_len(2)]
    genes <- genes[seq_len(min(length(genes), max_genes))]
    paste(c(fields[1:2], genes), collapse = "\t")
  }, character(1))
  writeLines(limited, outfile)
  outfile
}

extract_job_id <- function(x) {
  candidates <- c(
    x$job_id,
    x$result$job_id,
    x$result$params$job_id
  )
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (length(candidates) == 0) NA_character_ else as.character(candidates[[1]])
}

extract_download_url <- function(x) {
  candidates <- c(
    x$download_url,
    x$result$download_url
  )
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (length(candidates) == 0) return(NA_character_)
  out <- as.character(candidates[[1]])
  if (startsWith(out, "//")) out <- paste0("https:", out)
  out
}

parse_gct <- function(path) {
  raw_lines <- readLines(path, n = 3, warn = FALSE)
  skip <- ifelse(length(raw_lines) >= 2 && grepl("^#", raw_lines[1]), 2, 0)
  df <- fread(path, skip = skip, data.table = FALSE)
  if (nrow(df) == 0) return(NULL)
  id_col <- colnames(df)[1]
  desc_col <- intersect(c("Description", "desc", "pert_iname", "pert_desc"), colnames(df))[1]
  value_cols <- setdiff(colnames(df), c(id_col, desc_col))
  bind_rows(lapply(value_cols, function(state_col) {
    tibble(
      state_query = state_col,
      drug = if (!is.na(desc_col)) as.character(df[[desc_col]]) else as.character(df[[id_col]]),
      clue_id = as.character(df[[id_col]]),
      connectivity_score = suppressWarnings(as.numeric(df[[state_col]]))
    )
  }))
}

annotate_targets <- function(drugs, api_key) {
  drugs <- unique(drugs[!is.na(drugs) & nzchar(drugs)])
  if (length(drugs) == 0) return(tibble())

  bind_rows(lapply(split(drugs, ceiling(seq_along(drugs) / 50)), function(chunk) {
    filter <- toJSON(
      list(
        where = list(pert_iname = list(inq = as.list(chunk))),
        fields = list(pert_iname = TRUE, name = TRUE)
      ),
      auto_unbox = TRUE
    )
    res <- GET(
      "https://api.clue.io/api/rep_drug_target/",
      query = list(filter = filter),
      add_headers(user_key = api_key, Accept = "application/json")
    )
    if (status_code(res) >= 300) return(tibble())
    txt <- content(res, as = "text", encoding = "UTF-8")
    x <- fromJSON(txt, flatten = TRUE)
    if (length(x) == 0) return(tibble())
    as_tibble(x) %>%
      transmute(
        drug = as.character(pert_iname),
        target_gene = as.character(name)
      )
  })) %>%
    filter(!is.na(drug), !is.na(target_gene)) %>%
    distinct()
}

####################
# input checks
####################

up_gmt <- file.path(input_dir, "Auto_clue_up_entrez.gmt")
down_gmt <- file.path(input_dir, "Auto_clue_down_entrez.gmt")
if (!file.exists(up_gmt) || !file.exists(down_gmt)) {
  write_status("missing_inputs", "Missing CLUE GMT files. Run Auto_drug_reversal_inputs.R first.")
  quit(save = "no", status = 0)
}

api_key <- Sys.getenv("CLUE_API_KEY", Sys.getenv("CLUE_KEY", ""))
if (!nzchar(api_key)) {
  write_status(
    "missing_api_key",
    "No CLUE_API_KEY or CLUE_KEY environment variable is set. GMT files are ready under Auto_drug_reversal/clue_inputs/."
  )
  quit(save = "no", status = 0)
}

submit_up_gmt <- file.path(out_dir, paste0("Auto_clue_submitted_up_top", max_genes, "_entrez.gmt"))
submit_down_gmt <- file.path(out_dir, paste0("Auto_clue_submitted_down_top", max_genes, "_entrez.gmt"))
up_gmt <- limit_gmt(up_gmt, submit_up_gmt, max_genes, keep_states)
down_gmt <- limit_gmt(down_gmt, submit_down_gmt, max_genes, keep_states)

####################
# submit and poll
####################

if (nzchar(existing_job_id)) {
  job_id <- existing_job_id
  submit_json <- list(job_id = job_id)
} else {
  message("Submitting CLUE L1000 Touchstone query.")

  payload <- list(
    tool_id = tool_id,
    data_type = "L1000",
    name = query_name,
    dataset = "Touchstone",
    ignoreWarnings = TRUE
  )
  payload[["uptag-cmapfile"]] <- read_body_text(up_gmt)
  payload[["dntag-cmapfile"]] <- read_body_text(down_gmt)

  submit_res <- POST(
    "https://api.clue.io/api/jobs",
    add_headers(user_key = api_key, `Content-Type` = "application/json", Accept = "application/json"),
    body = payload,
    encode = "json"
  )

  submit_txt <- content(submit_res, as = "text", encoding = "UTF-8")
  writeLines(redact_response(submit_txt), file.path(out_dir, "Auto_clue_submit_response.json"))

  if (status_code(submit_res) >= 300) {
    write_status("submit_failed", paste("HTTP", status_code(submit_res), redact_response(submit_txt)))
    quit(save = "no", status = 0)
  }

  submit_json <- fromJSON(submit_txt, simplifyVector = FALSE)
  job_id <- extract_job_id(submit_json)
  if (is.na(job_id) || !nzchar(job_id)) {
    write_status("submit_failed", "CLUE response did not contain a job_id.")
    quit(save = "no", status = 0)
  }
}

message("Polling CLUE job: ", job_id)
download_url <- extract_download_url(submit_json)
status_json <- submit_json

for (i in seq_len(max_polls)) {
  if (!is.na(download_url) && nzchar(download_url)) break
  Sys.sleep(poll_seconds)
  poll_res <- GET(
    paste0("https://api.clue.io/api/jobs/findByJobId/", job_id),
    add_headers(user_key = api_key, Accept = "application/json")
  )
  poll_txt <- content(poll_res, as = "text", encoding = "UTF-8")
  writeLines(redact_response(poll_txt), file.path(out_dir, "Auto_clue_poll_response.json"))
  if (status_code(poll_res) >= 300) next
  status_json <- fromJSON(poll_txt, simplifyVector = FALSE)
  if (!is.null(status_json$status) && identical(status_json$status, "error")) {
    msg <- status_json$errorMessage
    if (is.null(msg)) msg <- "CLUE job returned status=error."
    write_status("clue_error", msg, job_id = job_id)
    quit(save = "no", status = 0)
  }
  download_url <- extract_download_url(status_json)
}

if (is.na(download_url) || !nzchar(download_url)) {
  write_status("poll_incomplete", "CLUE job did not provide a download URL within the polling window.", job_id = job_id)
  quit(save = "no", status = 0)
}

tar_path <- file.path(out_dir, "Auto_clue_results.tar.gz")
download.file(download_url, tar_path, mode = "wb", quiet = FALSE)
untar(tar_path, exdir = raw_dir)

####################
# parse results
####################

gct_files <- list.files(raw_dir, pattern = "\\.gct$", recursive = TRUE, full.names = TRUE)
if (length(gct_files) == 0) {
  write_status("downloaded_unparsed", "CLUE results downloaded, but no .gct files were found for lightweight parsing.", job_id, download_url)
  quit(save = "no", status = 0)
}

parsed <- bind_rows(lapply(gct_files, parse_gct))
if (nrow(parsed) == 0) {
  write_status("downloaded_unparsed", "CLUE .gct files were found but could not be parsed.", job_id, download_url)
  quit(save = "no", status = 0)
}

targets <- annotate_targets(parsed$drug, api_key)
target_summary <- targets %>%
  group_by(drug) %>%
  summarise(target_genes = paste(sort(unique(target_gene)), collapse = ";"), .groups = "drop")

ranked <- parsed %>%
  mutate(
    state = gsub("^Auto_|_up$|_down$", "", state_query),
    pipeline = "CLUE_FALLBACK",
    score = connectivity_score
  ) %>%
  filter(is.finite(score)) %>%
  group_by(state, drug) %>%
  summarise(
    score = mean(score, na.rm = TRUE),
    clue_id = clue_id[which.min(connectivity_score)][1],
    .groups = "drop"
  ) %>%
  left_join(target_summary, by = "drug") %>%
  arrange(state, score, drug) %>%
  group_by(state) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  transmute(
    state,
    pipeline,
    drug,
    rank,
    score,
    p_value = NA_real_,
    fdr = NA_real_,
    moa = NA_character_,
    target_genes = target_genes,
    clue_id
  )

fwrite(ranked, file.path(out_dir, "Auto_clue_ranked_drugs.csv"))
write_status("complete", paste("Parsed", nrow(ranked), "CLUE fallback drug-state rows."), job_id, download_url)

message("CLUE fallback complete.")
