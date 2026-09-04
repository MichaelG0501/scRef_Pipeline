####################
# Analysis registry:
#   Status: active, terminal
#   Script: analysis/TCGA/tcga_cna_recurrent_event_validation.R
#   Methodology: analysis/methodology/TCGA/tcga_cna_recurrent_event_validation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg
#     - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#     - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_tpm_matrix.rds
#     - ref_outs/TCGA/gender_validation/intermediate/Auto_tcga_gender_gsva_scores_centred17.rds, if present
#     - ref_outs/Auto_cna_subclone_expression/threshold_1/tables/Auto_v2_recomputed_recurrent_cna_event_summary.csv
#     - ref_outs/OAC_CNV.xlsx
#     - ref_outs/41588_2018_331_MOESM3_ESM.xlsx
#     - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds, if GSVA cache absent
#   Outputs:
#     - ref_outs/TCGA/cna_recurrent_event_validation/intermediate/
#     - ref_outs/TCGA/cna_recurrent_event_validation/tables/
#     - ref_outs/TCGA/cna_recurrent_event_validation/figures/
#     - ref_outs/TCGA/cna_recurrent_event_validation/logs/
#     - updates/new_updates/summaries/Auto_tcga_cna_recurrent_event_validation_summary.csv
#   Cache/replot behavior:
#     - TCGA arm calls and TCGA feature scores are cached and reused unless SCREF_FORCE_REBUILD=TRUE.
#     - SCREF_REPLOT_ONLY=TRUE reuses caches and regenerates tables/figures.
#   Run:
#     Rscript analysis/TCGA/tcga_cna_recurrent_event_validation.R
#   Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GSVA)
  library(readxl)
  library(scales)
  library(stringr)
  library(tidyr)
})

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_helpers.R")

setwd(SCREF_PROJECT_DIR)
set.seed(42)

####################
# 1) Paths and constants
####################
script_path <- file.path(SCREF_ANALYSIS_DIR, "TCGA", "tcga_cna_recurrent_event_validation.R")
out_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "cna_recurrent_event_validation")
tiers <- ensure_output_tiers(out_dir)
summary_dir <- SCREF_SUMMARY_DIR
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

segment_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg"
tcga_recon_dir <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "esca_gdc_reconstruction")
tcga_meta_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_meta.rds")
tcga_matrix_path <- file.path(tcga_recon_dir, "intermediate", "Auto_tcga_esca_tpm_matrix.rds")
tcga_mixture_path <- file.path(tcga_recon_dir, "tables", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
gender_gsva_cache_path <- file.path(SCREF_REF_OUTS_DIR, "TCGA", "gender_validation", "intermediate", "Auto_tcga_gender_gsva_scores_centred17.rds")
local_gsva_cache_path <- file.path(tiers[["intermediate"]], "Auto_tcga_cna_validation_gsva_scores_centred17.rds")
arm_cache_path <- file.path(tiers[["intermediate"]], "Auto_tcga_weighted_arm_cna_calls.rds")
sc_cna_dir <- file.path(SCREF_REF_OUTS_DIR, "Auto_cna_subclone_expression", "threshold_1")
sc_recurrent_path <- file.path(sc_cna_dir, "tables", "Auto_v2_recomputed_recurrent_cna_event_summary.csv")
sc_event_tests_path <- file.path(sc_cna_dir, "tables", "Auto_v2_recurrent_cna_event_feature_tests.csv")
oac_cnv_path <- file.path(SCREF_REF_OUTS_DIR, "OAC_CNV.xlsx")
occams_path <- file.path(SCREF_REF_OUTS_DIR, "41588_2018_331_MOESM3_ESM.xlsx")

candidate_thresholds <- seq(0.05, 0.30, by = 0.01)
recurrent_min_sample_fraction <- 0.15
recurrent_min_samples <- 3L
max_discovered_events <- 30L

chrom_levels <- paste0("chr", c(1:22, "X"))
arm_levels <- as.vector(rbind(paste0(chrom_levels, "p"), paste0(chrom_levels, "q")))

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

chrom_sizes <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
  chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
  chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
  chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468, chrX = 156040895
)

required_paths <- c(segment_path, tcga_meta_path, sc_recurrent_path, oac_cnv_path, occams_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = ", "))
}

run_summary <- start_run_summary(
  script = script_path,
  input_files = c(segment_path, tcga_meta_path, tcga_matrix_path, gender_gsva_cache_path,
                  sc_recurrent_path, sc_event_tests_path, oac_cnv_path, occams_path),
  output_files = c(
    file.path(tiers[["tables"]], "Auto_tcga_arm_cna_long.csv"),
    file.path(tiers[["tables"]], "Auto_tcga_cna_threshold_optimization.csv"),
    file.path(tiers[["tables"]], "Auto_scRef_recurrent_events_tcga_feature_tests.csv"),
    file.path(tiers[["tables"]], "Auto_tcga_discovered_recurrent_event_feature_tests.csv"),
    file.path(tiers[["figures"]], "Auto_tcga_cna_event_association_dotplots.pdf"),
    file.path(summary_dir, "Auto_tcga_cna_recurrent_event_validation_summary.csv")
  ),
  parameters = list(
    cohort = "TCGA EAC primary tumors with RNA-seq and segment CNA",
    threshold_grid = paste(range(candidate_thresholds), collapse = "-"),
    recurrent_min_sample_fraction = recurrent_min_sample_fraction,
    recurrent_min_samples = recurrent_min_samples
  )
)

mp_descriptions <- SCREF_MP_DESCRIPTIONS
mp_order <- SCREF_MP_ORDER
state_order <- SCREF_PRIMARY_STATE_ORDER

p_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

####################
# 2) Annotation workbook parsers
####################
normalise_gene_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\*$", "", x)
  x <- gsub("[^A-Z0-9.-]", "", x)
  x[nchar(x) == 0] <- NA_character_
  x
}

split_gene_string <- function(x) {
  x <- as.character(x)
  x <- gsub("\\(([^)]*)\\)", ",\\1", x)
  x <- gsub("\\[|\\]", "", x)
  pieces <- unlist(strsplit(x, "[,;]+"))
  out <- character()
  for (piece in pieces) {
    piece <- trimws(piece)
    if (!nzchar(piece) || is.na(piece)) next
    parts <- unlist(strsplit(piece, "/"))
    parts <- trimws(parts)
    if (length(parts) > 1) {
      prefix <- sub("^([A-Za-z]+).*", "\\1", parts[1])
      expanded <- vapply(seq_along(parts), function(i) {
        p <- parts[i]
        if (i > 1 && grepl("^[0-9]+[A-Za-z]*$", p) && grepl("^[A-Za-z]+", parts[1])) {
          paste0(prefix, p)
        } else {
          p
        }
      }, character(1))
      out <- c(out, expanded)
    } else {
      out <- c(out, parts)
    }
  }
  out <- gsub("\\s+", "", out)
  out <- normalise_gene_symbol(out)
  out <- out[!is.na(out)]
  out <- out[!grepl("^HSA-", out)]
  out <- out[!grepl("^[0-9]+$", out)]
  unique(out)
}

cytoband_to_arm <- function(cytoband) {
  cytoband <- as.character(cytoband)
  cytoband <- gsub("^chr", "", cytoband, ignore.case = TRUE)
  m <- regexec("^([0-9]+|X|Y)([pq])", cytoband)
  parts <- regmatches(cytoband, m)
  vapply(parts, function(x) {
    if (length(x) < 3) return(NA_character_)
    chr <- paste0("chr", x[2])
    if (!chr %in% chrom_levels) return(NA_character_)
    paste0(chr, x[3])
  }, character(1))
}

parse_peak_coord <- function(x) {
  x <- as.character(x)
  m <- regexec("(chr[0-9XY]+):([0-9]+)-([0-9]+)", x)
  parts <- regmatches(x, m)
  if (length(parts[[1]]) < 4) {
    return(data.frame(chr = NA_character_, start = NA_real_, end = NA_real_, arm = NA_character_))
  }
  chr <- parts[[1]][2]
  start <- as.numeric(parts[[1]][3])
  end <- as.numeric(parts[[1]][4])
  mid <- mean(c(start, end), na.rm = TRUE)
  arm <- if (chr %in% names(centromere_pos)) {
    paste0(chr, ifelse(mid <= centromere_pos[[chr]], "p", "q"))
  } else {
    NA_character_
  }
  data.frame(chr = chr, start = start, end = end, arm = arm)
}

parse_oac_cnv <- function(path) {
  raw <- readxl::read_excel(path, sheet = 1, col_names = FALSE)
  rows <- list()
  current <- NA_character_
  row_i <- 1L
  for (i in seq_len(nrow(raw))) {
    first <- as.character(raw[[1]][i])
    if (is.na(first)) next
    if (grepl("^CNV gain", first, ignore.case = TRUE)) {
      current <- "gain"
      next
    }
    if (grepl("^CNV loss", first, ignore.case = TRUE)) {
      current <- "loss"
      next
    }
    rank_num <- suppressWarnings(as.integer(first))
    if (!is.na(rank_num) && !is.na(current)) {
      genes_raw <- as.character(raw[[4]][i])
      genes <- split_gene_string(genes_raw)
      rows[[row_i]] <- data.frame(
        source = "OAC_CNV_curated",
        workbook = basename(path),
        sheet = "Sheet1",
        excel_row = i,
        direction = current,
        rank = rank_num,
        cytoband = as.character(raw[[2]][i]),
        region = as.character(raw[[3]][i]),
        genes_raw = genes_raw,
        gene = if (length(genes) > 0) genes else NA_character_,
        frequency_label = as.character(raw[[5]][i]),
        pathway = as.character(raw[[6]][i]),
        clinical_relevance = as.character(raw[[7]][i]),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  bind_rows(rows) |>
    mutate(
      arm = cytoband_to_arm(.data$cytoband),
      event_id = paste(.data$direction, .data$arm, sep = "_")
    )
}

parse_occams_gistic_peaks <- function(path, sheet, direction) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE)
  first_col <- as.character(raw[[1]])
  cytoband_row <- which(tolower(first_col) == "cytoband")[1]
  q_row <- which(tolower(first_col) == "q value")[1]
  residual_q_row <- which(tolower(first_col) == "residual q value")[1]
  boundary_row <- which(tolower(first_col) == "wide peak boundaries")[1]
  gene_row <- which(tolower(first_col) == "genes in wide peak")[1]
  if (is.na(cytoband_row) || is.na(boundary_row) || is.na(gene_row)) return(data.frame())
  rows <- list()
  row_i <- 1L
  for (j in 2:ncol(raw)) {
    cytoband <- as.character(raw[[j]][cytoband_row])
    if (is.na(cytoband) || !nzchar(cytoband)) next
    genes <- as.character(raw[[j]][(gene_row + 1):nrow(raw)])
    genes <- normalise_gene_symbol(genes[!is.na(genes) & nzchar(genes)])
    genes <- genes[!is.na(genes)]
    genes <- genes[!grepl("^HSA-", genes)]
    coord <- parse_peak_coord(raw[[j]][boundary_row])
    rows[[row_i]] <- data.frame(
      source = "OCCAMS_GISTIC_peak",
      workbook = basename(path),
      sheet = sheet,
      excel_column = j,
      direction = direction,
      cytoband = cytoband,
      q_value = suppressWarnings(as.numeric(raw[[j]][q_row])),
      residual_q_value = suppressWarnings(as.numeric(raw[[j]][residual_q_row])),
      wide_peak_boundaries = as.character(raw[[j]][boundary_row]),
      chr = coord$chr,
      start = coord$start,
      end = coord$end,
      arm = ifelse(!is.na(coord$arm), coord$arm, cytoband_to_arm(cytoband)),
      gene = if (length(genes) > 0) genes else NA_character_,
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
  bind_rows(rows) |>
    mutate(event_id = paste(.data$direction, .data$arm, sep = "_"))
}

parse_occams_driver_sheet <- function(path, sheet, direction, confidence) {
  raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE)
  header_row <- which(as.character(raw[[1]]) == "Gene" & as.character(raw[[2]]) == "hgnc_symbol")[1]
  if (is.na(header_row)) return(data.frame())
  x <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  names(x) <- as.character(unlist(raw[header_row, ]))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  peak_col <- grep("GISTIC Peak", colnames(x), value = TRUE)[1]
  if (is.na(peak_col)) return(data.frame())
  x$excel_row <- seq(header_row + 1, nrow(raw))
  out <- x |>
    filter(!is.na(.data$hgnc_symbol), nzchar(as.character(.data$hgnc_symbol))) |>
    mutate(
      source = paste0("OCCAMS_", confidence, "_driver"),
      workbook = basename(path),
      sheet = sheet,
      direction = direction,
      ensembl_id = as.character(.data$Gene),
      gene_raw = as.character(.data$hgnc_symbol),
      gene = normalise_gene_symbol(.data$hgnc_symbol),
      peak = as.character(.data[[peak_col]]),
      cytoband = sub(" .*", "", .data$peak)
    ) |>
    filter(!is.na(.data$gene))
  coords <- bind_rows(lapply(out$peak, parse_peak_coord))
  bind_cols(out, coords) |>
    mutate(
      arm = ifelse(!is.na(.data$arm), .data$arm, cytoband_to_arm(.data$cytoband)),
      event_id = paste(.data$direction, .data$arm, sep = "_")
    ) |>
    transmute(
      source = .data$source,
      workbook = .data$workbook,
      sheet = .data$sheet,
      excel_row = .data$excel_row,
      direction = .data$direction,
      ensembl_id = .data$ensembl_id,
      gene_raw = .data$gene_raw,
      gene = .data$gene,
      cytoband = .data$cytoband,
      arm = .data$arm,
      peak = .data$peak,
      event_id = .data$event_id
    )
}

message("Parsing OAC_CNV and OCCAMS annotation workbooks")
oac_cnv <- parse_oac_cnv(oac_cnv_path)
occams_gistic <- bind_rows(
  parse_occams_gistic_peaks(occams_path, "ST1 GISTIC amplification peaks", "gain"),
  parse_occams_gistic_peaks(occams_path, "ST2 GISTIC Deletion peaks", "loss")
)
occams_drivers <- bind_rows(
  parse_occams_driver_sheet(occams_path, "ST3 High Confidence Del Drivers", "loss", "high_confidence_deletion"),
  parse_occams_driver_sheet(occams_path, "ST4 Candidate Del Drivers", "loss", "candidate_deletion"),
  parse_occams_driver_sheet(occams_path, "ST5 High confidence Amp Drivers", "gain", "high_confidence_amplification"),
  parse_occams_driver_sheet(occams_path, "ST6 Candidate Amp Drivers", "gain", "candidate_amplification")
)

write.csv(oac_cnv, file.path(tiers[["tables"]], "Auto_oac_cnv_curated_reparsed.csv"), row.names = FALSE)
write.csv(occams_gistic, file.path(tiers[["tables"]], "Auto_occams_gistic_peaks_reparsed.csv"), row.names = FALSE)
write.csv(occams_drivers, file.path(tiers[["tables"]], "Auto_occams_cna_driver_genes_reparsed.csv"), row.names = FALSE)

cnv_annotation_long <- bind_rows(
  oac_cnv |>
    transmute(source, workbook, sheet, excel_row, direction, arm, event_id, cytoband,
              gene, gene_raw = gene, frequency_label, pathway, clinical_relevance),
  occams_drivers |>
    transmute(source, workbook, sheet, excel_row, direction, arm, event_id, cytoband,
              gene, gene_raw, frequency_label = NA_character_, pathway = NA_character_,
              clinical_relevance = NA_character_),
  occams_gistic |>
    transmute(source, workbook, sheet, excel_row = NA_integer_, direction, arm, event_id, cytoband,
              gene, gene_raw = gene, frequency_label = NA_character_, pathway = NA_character_,
              clinical_relevance = NA_character_)
) |>
  filter(!is.na(.data$arm), .data$arm %in% arm_levels)

summarise_event_genes <- function(gene, source, n = 24L) {
  gene <- as.character(gene)
  source <- as.character(source)
  ok <- !is.na(gene) & nzchar(gene)
  gene <- gene[ok]
  source <- source[ok]
  priority <- sort(unique(gene[!grepl("GISTIC_peak", source)]))
  background <- sort(setdiff(unique(gene[grepl("GISTIC_peak", source)]), priority))
  paste(head(c(priority, background), n), collapse = ", ")
}

event_annotation <- cnv_annotation_long |>
  group_by(.data$event_id, .data$direction, .data$arm) |>
  summarise(
    annotation_sources = paste(sort(unique(.data$source)), collapse = "; "),
    known_genes = summarise_event_genes(.data$gene, .data$source, 24L),
    raw_gene_notes = paste(sort(unique(.data$gene_raw[!is.na(.data$gene_raw) & .data$gene_raw != .data$gene])), collapse = ", "),
    cytobands = paste(sort(unique(.data$cytoband[!is.na(.data$cytoband)])), collapse = ", "),
    pathways = paste(sort(unique(.data$pathway[!is.na(.data$pathway) & nzchar(.data$pathway)])), collapse = "; "),
    clinical_relevance = paste(sort(unique(.data$clinical_relevance[!is.na(.data$clinical_relevance) & nzchar(.data$clinical_relevance)])), collapse = "; "),
    .groups = "drop"
  )

write.csv(cnv_annotation_long, file.path(tiers[["tables"]], "Auto_cnv_annotation_long_reparsed.csv"), row.names = FALSE)
write.csv(event_annotation, file.path(tiers[["tables"]], "Auto_cnv_event_annotation_summary_reparsed.csv"), row.names = FALSE)

####################
# 3) Weighted TCGA segment means by chromosome arm
####################
message("Computing weighted TCGA arm CNA means")
build_arm_calls <- function() {
  seg <- fread(segment_path)
  colnames(seg) <- make.names(colnames(seg))
  if (!all(c("ID", "chrom", "loc.start", "loc.end", "seg.mean") %in% colnames(seg))) {
    stop("Segment file does not have expected columns ID/chrom/loc.start/loc.end/seg.mean")
  }
  seg <- seg |>
    mutate(
      seg_sample_id = as.character(.data$ID),
      sample_key = substr(.data$seg_sample_id, 1, 15),
      chr = paste0("chr", gsub("^chr", "", as.character(.data$chrom), ignore.case = TRUE)),
      start = as.numeric(.data$loc.start),
      end = as.numeric(.data$loc.end),
      seg_mean = as.numeric(.data$seg.mean)
    ) |>
    filter(.data$chr %in% chrom_levels, is.finite(.data$start), is.finite(.data$end), is.finite(.data$seg_mean)) |>
    mutate(
      centromere = centromere_pos[.data$chr],
      chr_size = chrom_sizes[.data$chr],
      p_overlap = pmax(0, pmin(.data$end, .data$centromere) - .data$start + 1),
      q_overlap = pmax(0, pmin(.data$end, .data$chr_size) - pmax(.data$start, .data$centromere + 1) + 1)
    )

  arm_piece <- bind_rows(
    seg |>
      filter(.data$p_overlap > 0) |>
      transmute(seg_sample_id, sample_key, arm = paste0(.data$chr, "p"), segment_overlap_bp = .data$p_overlap, seg_mean),
    seg |>
      filter(.data$q_overlap > 0) |>
      transmute(seg_sample_id, sample_key, arm = paste0(.data$chr, "q"), segment_overlap_bp = .data$q_overlap, seg_mean)
  )

  arm_long <- arm_piece |>
    group_by(.data$seg_sample_id, .data$sample_key, .data$arm) |>
    summarise(
      arm_mean = weighted.mean(.data$seg_mean, .data$segment_overlap_bp, na.rm = TRUE),
      covered_bp = sum(.data$segment_overlap_bp, na.rm = TRUE),
      n_segments = n(),
      .groups = "drop"
    ) |>
    mutate(arm = factor(.data$arm, levels = arm_levels)) |>
    arrange(.data$sample_key, .data$arm)

  all_sample_arm <- tidyr::crossing(
    distinct(arm_long, .data$seg_sample_id, .data$sample_key),
    arm = factor(arm_levels, levels = arm_levels)
  )
  all_sample_arm |>
    left_join(arm_long, by = c("seg_sample_id", "sample_key", "arm")) |>
    mutate(
      arm = as.character(.data$arm),
      arm_mean = ifelse(is.na(.data$arm_mean), 0, .data$arm_mean),
      covered_bp = ifelse(is.na(.data$covered_bp), 0, .data$covered_bp),
      n_segments = ifelse(is.na(.data$n_segments), 0L, .data$n_segments)
    )
}

arm_cache <- load_or_build_cache(
  arm_cache_path,
  build_fun = build_arm_calls,
  label = "tcga_weighted_arm_cna_calls"
)
run_summary <- add_cache_status(run_summary, arm_cache$label, arm_cache$path, arm_cache$reused)
tcga_arm_long_all <- arm_cache$value

meta_tcga <- readRDS(tcga_meta_path) |>
  mutate(
    sample_key = substr(.data$sample_barcode, 1, 15),
    HistologyGroup = ifelse(is.na(.data$HistologyGroup), "Other", .data$HistologyGroup)
  )

tcga_eac_meta <- meta_tcga |>
  filter(.data$sample_type_code == "01", .data$HistologyGroup == "EAC") |>
  arrange(.data$sample_key, .data$sample_barcode) |>
  distinct(.data$sample_key, .keep_all = TRUE)

tcga_arm_long <- tcga_arm_long_all |>
  inner_join(
    tcga_eac_meta |>
      select(.data$sample_key, .data$sample_barcode, .data$case_barcode, .data$HistologyGroup,
             .data$sample_type_code, .data$sample_type),
    by = "sample_key"
  )

if (n_distinct(tcga_arm_long$sample_barcode) < 20) {
  stop("Too few TCGA EAC primary samples have segment arm calls: ", n_distinct(tcga_arm_long$sample_barcode))
}

write.csv(tcga_arm_long, file.path(tiers[["tables"]], "Auto_tcga_arm_cna_long.csv"), row.names = FALSE)

####################
# 4) Threshold optimization and event calls
####################
message("Selecting TCGA arm-call threshold")
reference_events <- oac_cnv |>
  filter(!is.na(.data$event_id), .data$arm %in% arm_levels) |>
  distinct(.data$event_id, .data$direction, .data$arm) |>
  mutate(reference_source = "OAC_CNV_curated")

summarise_events_at_threshold <- function(arm_long, threshold) {
  n_total <- n_distinct(arm_long$sample_barcode)
  min_n <- max(recurrent_min_samples, ceiling(recurrent_min_sample_fraction * n_total))
  bind_rows(
    arm_long |>
      filter(.data$arm_mean >= threshold) |>
      group_by(.data$arm) |>
      summarise(direction = "gain", n_samples_event = n_distinct(.data$sample_barcode), .groups = "drop"),
    arm_long |>
      filter(.data$arm_mean <= -threshold) |>
      group_by(.data$arm) |>
      summarise(direction = "loss", n_samples_event = n_distinct(.data$sample_barcode), .groups = "drop")
  ) |>
    mutate(
      event_id = paste(.data$direction, .data$arm, sep = "_"),
      n_samples_total = n_total,
      frac_samples_event = .data$n_samples_event / .data$n_samples_total,
      threshold = threshold,
      is_recurrent = .data$n_samples_event >= min_n
    ) |>
    select(.data$threshold, .data$event_id, .data$direction, .data$arm,
           .data$n_samples_event, .data$n_samples_total, .data$frac_samples_event, .data$is_recurrent)
}

threshold_event_summaries <- bind_rows(lapply(candidate_thresholds, function(thr) {
  summarise_events_at_threshold(tcga_arm_long, thr)
}))

threshold_optimization <- bind_rows(lapply(candidate_thresholds, function(thr) {
  pred <- threshold_event_summaries |>
    filter(.data$threshold == thr, .data$is_recurrent) |>
    pull(.data$event_id)
  ref <- reference_events$event_id
  tp <- length(intersect(pred, ref))
  fp <- length(setdiff(pred, ref))
  fn <- length(setdiff(ref, pred))
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  f1 <- ifelse(is.finite(precision + recall) && precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  data.frame(
    threshold = thr,
    n_tcga_recurrent_events = length(pred),
    n_reference_events = length(ref),
    true_positive = tp,
    false_positive = fp,
    false_negative = fn,
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = ifelse(length(union(pred, ref)) > 0, length(intersect(pred, ref)) / length(union(pred, ref)), NA_real_),
    stringsAsFactors = FALSE
  )
})) |>
  arrange(desc(.data$f1), desc(.data$jaccard), abs(.data$threshold - 0.10), .data$threshold) |>
  mutate(selected = row_number() == 1L)

selected_threshold <- threshold_optimization$threshold[threshold_optimization$selected][1]
message("Selected threshold: ", selected_threshold)

write.csv(threshold_event_summaries, file.path(tiers[["tables"]], "Auto_tcga_cna_event_threshold_sensitivity.csv"), row.names = FALSE)
write.csv(threshold_optimization, file.path(tiers[["tables"]], "Auto_tcga_cna_threshold_optimization.csv"), row.names = FALSE)

tcga_arm_calls <- tcga_arm_long |>
  mutate(
    threshold = selected_threshold,
    cna_call = case_when(
      .data$arm_mean >= selected_threshold ~ 1L,
      .data$arm_mean <= -selected_threshold ~ -1L,
      TRUE ~ 0L
    )
  )

tcga_discovered_events <- summarise_events_at_threshold(tcga_arm_long, selected_threshold) |>
  filter(.data$is_recurrent) |>
  left_join(event_annotation, by = c("event_id", "direction", "arm")) |>
  mutate(
    known_genes = ifelse(is.na(.data$known_genes), "", .data$known_genes),
    annotation_sources = ifelse(is.na(.data$annotation_sources), "", .data$annotation_sources)
  ) |>
  arrange(desc(.data$frac_samples_event), .data$event_id)

write.csv(tcga_arm_calls, file.path(tiers[["tables"]], "Auto_tcga_arm_cna_calls_selected_threshold.csv"), row.names = FALSE)
write.csv(tcga_discovered_events, file.path(tiers[["tables"]], "Auto_tcga_recurrent_cna_events.csv"), row.names = FALSE)

####################
# 5) TCGA MP/state scores
####################
message("Loading or computing TCGA MP/state feature scores")
make_feature_labels <- function(mp_features, state_features) {
  mp_labels <- paste0(mp_features, " ", mp_descriptions[mp_features])
  mp_labels[is.na(mp_labels)] <- label_mps(mp_features[is.na(mp_labels)])
  names(mp_labels) <- mp_features
  state_labels <- setNames(state_features, state_features)
  c(mp_labels, state_labels)
}

run_gsva_scores <- function(expr_mat, gene_sets) {
  gene_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(expr_mat)))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]
  if (length(gene_sets) == 0) stop("No gene sets retained at >=5 genes for GSVA.")
  suppressMessages(gsva(expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian"))
}

build_gene_sets <- function() {
  mp_genes <- readRDS(SCREF_MP_GENES_RDS)
  mp_genes <- mp_genes[intersect(mp_order, names(mp_genes))]
  mp_sets <- mp_genes
  state_groups <- SCREF_STATE_GROUPS
  state_sets <- lapply(state_groups, function(mps) {
    unique(unlist(mp_genes[intersect(mps, names(mp_genes))], use.names = FALSE))
  })
  state_sets <- state_sets[sapply(state_sets, length) >= 5]
  state_sets <- state_sets[intersect(state_order, names(state_sets))]
  list(mp_sets = mp_sets, state_sets = state_sets)
}

load_or_build_tcga_scores <- function() {
  if (file.exists(gender_gsva_cache_path) && reuse_cache_enabled()) {
    obj <- readRDS(gender_gsva_cache_path)
    if (all(c("mp_scores", "state_scores", "mp_sets", "state_sets") %in% names(obj))) {
      message("Reusing gender validation GSVA cache: ", gender_gsva_cache_path)
      return(list(value = obj, reused = TRUE, label = "tcga_gender_gsva_scores", path = gender_gsva_cache_path))
    }
  }
  load_or_build_cache(
    local_gsva_cache_path,
    build_fun = function() {
      sets <- build_gene_sets()
      if (file.exists(tcga_matrix_path)) {
        tpm_mat <- readRDS(tcga_matrix_path)
      } else if (file.exists(tcga_mixture_path)) {
        tpm_df <- fread(tcga_mixture_path)
        tpm_mat <- as.matrix(tpm_df[, -1, with = FALSE])
        rownames(tpm_mat) <- tpm_df[[1]]
      } else {
        stop("Missing TCGA TPM matrix. Run analysis/TCGA/tcga_esca_reconstruct_data.R first.")
      }
      common_samples <- intersect(colnames(tpm_mat), meta_tcga$sample_barcode)
      tpm_mat <- tpm_mat[, common_samples, drop = FALSE]
      expr_log <- log2(tpm_mat + 1)
      expr_log[!is.finite(expr_log)] <- 0
      list(
        mp_scores = t(run_gsva_scores(expr_log, sets$mp_sets)),
        state_scores = t(run_gsva_scores(expr_log, sets$state_sets)),
        mp_sets = sets$mp_sets,
        state_sets = sets$state_sets,
        transform = "log2(TPM + 1)"
      )
    },
    label = "tcga_cna_validation_gsva_scores"
  )
}

gsva_cache <- load_or_build_tcga_scores()
run_summary <- add_cache_status(run_summary, gsva_cache$label, gsva_cache$path, gsva_cache$reused)
gsva_obj <- gsva_cache$value

mp_sets <- gsva_obj$mp_sets
state_sets <- gsva_obj$state_sets
mp_score_mat <- as.matrix(gsva_obj$mp_scores)
state_score_mat <- as.matrix(gsva_obj$state_scores)
feature_labels <- make_feature_labels(colnames(mp_score_mat), colnames(state_score_mat))

score_samples <- Reduce(intersect, list(tcga_eac_meta$sample_barcode, rownames(mp_score_mat), rownames(state_score_mat)))
tcga_score_meta <- tcga_eac_meta |>
  filter(.data$sample_barcode %in% score_samples) |>
  arrange(match(.data$sample_barcode, score_samples)) |>
  select(.data$sample_barcode, .data$sample_key, .data$case_barcode, .data$HistologyGroup,
         .data$sample_type_code, .data$sample_type)
tcga_score_df <- bind_cols(
  tcga_score_meta,
  as.data.frame(mp_score_mat[tcga_score_meta$sample_barcode, , drop = FALSE], check.names = FALSE),
  as.data.frame(state_score_mat[tcga_score_meta$sample_barcode, , drop = FALSE], check.names = FALSE)
)

feature_info <- bind_rows(
  data.frame(feature = colnames(mp_score_mat), feature_type = "MP", stringsAsFactors = FALSE),
  data.frame(feature = colnames(state_score_mat), feature_type = "State", stringsAsFactors = FALSE)
) |>
  mutate(
    feature_label = feature_labels[.data$feature],
    feature_group = case_when(
      .data$feature_type == "MP" ~ "Metaprogrammes",
      .data$feature_type == "State" ~ "Centred states",
      TRUE ~ .data$feature_type
    ),
    feature_order = case_when(
      .data$feature_type == "MP" ~ match(.data$feature, mp_order),
      .data$feature_type == "State" ~ match(.data$feature, state_order),
      TRUE ~ NA_integer_
    )
  ) |>
  filter(!is.na(.data$feature_order)) |>
  arrange(factor(.data$feature_group, levels = c("Metaprogrammes", "Centred states")), .data$feature_order)

mp_features <- feature_info |> filter(.data$feature_group == "Metaprogrammes") |> pull(.data$feature)
state_features <- feature_info |> filter(.data$feature_group == "Centred states") |> pull(.data$feature)
plot_features <- c(mp_features, state_features)
feature_group <- setNames(feature_info$feature_group, feature_info$feature)
feature_label_map <- setNames(feature_info$feature_label, feature_info$feature)

score_long <- bind_rows(
  tcga_score_df |>
    select(.data$sample_barcode, .data$sample_key, .data$case_barcode, all_of(colnames(mp_score_mat))) |>
    pivot_longer(cols = all_of(colnames(mp_score_mat)), names_to = "feature", values_to = "feature_value") |>
    mutate(feature_type = "MP"),
  tcga_score_df |>
    select(.data$sample_barcode, .data$sample_key, .data$case_barcode, all_of(colnames(state_score_mat))) |>
    pivot_longer(cols = all_of(colnames(state_score_mat)), names_to = "feature", values_to = "feature_value") |>
    mutate(feature_type = "State")
) |>
  filter(.data$feature %in% plot_features) |>
  left_join(feature_info, by = c("feature", "feature_type")) |>
  group_by(.data$feature) |>
  mutate(
    feature_mean = mean(.data$feature_value, na.rm = TRUE),
    feature_sd = sd(.data$feature_value, na.rm = TRUE),
    feature_z = (.data$feature_value - .data$feature_mean) / .data$feature_sd,
    feature_z = ifelse(is.finite(.data$feature_z), .data$feature_z, 0)
  ) |>
  ungroup()

write.csv(feature_info, file.path(tiers[["tables"]], "Auto_tcga_mp_state_features_tested.csv"), row.names = FALSE)

####################
# 6) Event association tests
####################
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  diffs <- outer(x, y, "-")
  (sum(diffs > 0) - sum(diffs < 0)) / (length(x) * length(y))
}

####################
# Prioritise familiar CNA genes in compact event labels so key drivers are not hidden.
####################
cna_gene_priority <- c(
  "MYC", "ERBB2", "CCNE1", "CDKN2A", "SMAD4", "TP53", "KRAS", "GATA6",
  "EGFR", "MET", "VEGFA", "PIK3CA", "APC", "ARID1A", "AURKA", "ZNF217",
  "MDM2", "CDK6", "FGFR1", "SOX2", "TERT", "RNF43", "SMAD2", "CDH1"
)

prioritise_known_genes <- function(known_genes, max_genes = 3L) {
  vapply(as.character(known_genes), function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    genes <- trimws(unlist(strsplit(x, ",")))
    genes <- gsub("\\*$", "", genes)
    genes <- genes[nzchar(genes)]
    genes <- unique(genes)
    priority <- cna_gene_priority[cna_gene_priority %in% genes]
    ordered <- unique(c(priority, genes))
    shown <- head(ordered, max_genes)
    suffix <- if (length(ordered) > max_genes) paste0(" +", length(ordered) - max_genes) else ""
    paste0(paste(shown, collapse = ", "), suffix)
  }, character(1))
}
####################

event_label <- function(event_id, known_genes = NULL) {
  direction <- sub("_.*$", "", event_id)
  arm <- sub("^[^_]+_", "", event_id)
  base <- paste0(ifelse(direction == "gain", "Gain ", "Loss "), arm)
  known_genes <- as.character(known_genes)
  known_genes[is.na(known_genes)] <- ""
  known_genes <- prioritise_known_genes(known_genes)
  paste0(base, ifelse(nzchar(known_genes), paste0(" (", known_genes, ")"), ""))
}

collapse_annotation_columns <- function(df) {
  if ("known_genes.y" %in% colnames(df)) {
    left <- if ("known_genes.x" %in% colnames(df)) df$known_genes.x else ""
    right <- df$known_genes.y
    df$known_genes <- ifelse(!is.na(right) & nzchar(right), right, left)
  } else if (!"known_genes" %in% colnames(df)) {
    df$known_genes <- ""
  }
  if ("annotation_sources.y" %in% colnames(df)) {
    left <- if ("annotation_sources.x" %in% colnames(df)) df$annotation_sources.x else ""
    right <- df$annotation_sources.y
    df$annotation_sources <- ifelse(!is.na(right) & nzchar(right), right, left)
  } else if (!"annotation_sources" %in% colnames(df)) {
    df$annotation_sources <- ""
  }
  df$known_genes[is.na(df$known_genes)] <- ""
  df$annotation_sources[is.na(df$annotation_sources)] <- ""
  df[, setdiff(colnames(df), c("known_genes.x", "known_genes.y", "annotation_sources.x", "annotation_sources.y")), drop = FALSE]
}

event_presence_for <- function(events_df) {
  events_df <- events_df |>
    select(.data$event_id, .data$direction, .data$arm, everything()) |>
    distinct(.data$event_id, .keep_all = TRUE)
  bind_rows(lapply(seq_len(nrow(events_df)), function(i) {
    ev <- events_df[i, ]
    tcga_arm_calls |>
      filter(.data$arm == ev$arm) |>
      transmute(
        event_id = ev$event_id,
        direction = ev$direction,
        arm = ev$arm,
        sample_barcode = .data$sample_barcode,
        sample_key = .data$sample_key,
        arm_mean = .data$arm_mean,
        cna_call = .data$cna_call,
        event_present = if (ev$direction == "gain") .data$cna_call == 1L else .data$cna_call == -1L
      )
  }))
}

test_event_features <- function(events_df, analysis_label) {
  events_df <- events_df |>
    select(.data$event_id, .data$direction, .data$arm, everything()) |>
    distinct(.data$event_id, .keep_all = TRUE)
  presence <- event_presence_for(events_df) |>
    inner_join(score_long, by = c("sample_barcode", "sample_key")) |>
    left_join(events_df |> select(.data$event_id, .data$known_genes, .data$annotation_sources), by = "event_id") |>
    mutate(
      analysis = analysis_label,
      event_label = event_label(.data$event_id, .data$known_genes)
    )

  tests <- presence |>
    filter(is.finite(.data$feature_value), !is.na(.data$event_present)) |>
    group_by(.data$analysis, .data$event_id, .data$event_label, .data$direction, .data$arm,
             .data$known_genes, .data$annotation_sources, .data$feature_type, .data$feature_group,
             .data$feature, .data$feature_label) |>
    group_modify(function(.x, .y) {
      event_vals <- .x$feature_value[.x$event_present]
      no_event_vals <- .x$feature_value[!.x$event_present]
      event_z <- .x$feature_z[.x$event_present]
      no_event_z <- .x$feature_z[!.x$event_present]
      p_value <- if (length(event_vals) >= 3 && length(no_event_vals) >= 3) {
        suppressWarnings(tryCatch(wilcox.test(event_vals, no_event_vals, exact = FALSE)$p.value, error = function(e) NA_real_))
      } else {
        NA_real_
      }
      tibble(
        n_event = length(event_vals),
        n_no_event = length(no_event_vals),
        median_event = median(event_vals, na.rm = TRUE),
        median_no_event = median(no_event_vals, na.rm = TRUE),
        mean_event = mean(event_vals, na.rm = TRUE),
        mean_no_event = mean(no_event_vals, na.rm = TRUE),
        delta_median = median_event - median_no_event,
        delta_mean = mean_event - mean_no_event,
        median_event_z = median(event_z, na.rm = TRUE),
        median_no_event_z = median(no_event_z, na.rm = TRUE),
        unpaired_median_std_delta = median_event_z - median_no_event_z,
        cliffs_delta = cliffs_delta(event_vals, no_event_vals),
        p_value = p_value
      )
    }) |>
    ungroup() |>
    group_by(.data$analysis, .data$feature_type) |>
    mutate(p_adj_by_feature_type = p.adjust(.data$p_value, method = "BH")) |>
    ungroup() |>
    group_by(.data$analysis) |>
    mutate(p_adj_global = p.adjust(.data$p_value, method = "BH")) |>
    ungroup() |>
    mutate(
      sig_label = p_to_stars(.data$p_adj_by_feature_type),
      neglog10_fdr = pmin(-log10(pmax(.data$p_adj_by_feature_type, 1e-12)), 12),
      primary_delta = .data$unpaired_median_std_delta,
      primary_p_adj_group = .data$p_adj_by_feature_type
    ) |>
    arrange(.data$p_adj_global, .data$p_value)

  list(values = presence, tests = tests)
}

sc_recurrent <- read.csv(sc_recurrent_path, check.names = FALSE) |>
  filter(.data$is_recurrent) |>
  mutate(
    arm = as.character(.data$arm),
    direction = as.character(.data$direction),
    event_id = as.character(.data$event_id)
  ) |>
  left_join(event_annotation, by = c("event_id", "direction", "arm")) |>
  collapse_annotation_columns() |>
  arrange(desc(.data$frac_samples_event), desc(.data$frac_subclones_event), .data$event_id) |>
  head(8)

sc_boxplot_events <- sc_recurrent |> head(4) |> pull(.data$event_id)

sc_validation <- test_event_features(sc_recurrent, "scRef_recurrent_events")
write.csv(sc_validation$values, file.path(tiers[["tables"]], "Auto_scRef_recurrent_events_tcga_feature_values.csv"), row.names = FALSE)
write.csv(sc_validation$tests, file.path(tiers[["tables"]], "Auto_scRef_recurrent_events_tcga_feature_tests.csv"), row.names = FALSE)

tcga_discovery_events_for_test <- tcga_discovered_events |>
  head(max_discovered_events)
tcga_discovery <- test_event_features(tcga_discovery_events_for_test, "tcga_discovered_recurrent_events")

sc_recurrent_lookup <- sc_recurrent |>
  select(sc_event_id = .data$event_id, sc_is_recurrent = .data$is_recurrent,
         sc_frac_samples_event = .data$frac_samples_event, sc_n_samples_event = .data$n_samples_event)

tcga_discovery_tests <- tcga_discovery$tests |>
  left_join(sc_recurrent_lookup, by = c("event_id" = "sc_event_id")) |>
  mutate(
    sc_is_recurrent = ifelse(is.na(.data$sc_is_recurrent), FALSE, .data$sc_is_recurrent),
    sc_frac_samples_event = ifelse(is.na(.data$sc_frac_samples_event), 0, .data$sc_frac_samples_event)
  )

write.csv(tcga_discovery$values, file.path(tiers[["tables"]], "Auto_tcga_discovered_recurrent_event_feature_values.csv"), row.names = FALSE)
write.csv(tcga_discovery_tests, file.path(tiers[["tables"]], "Auto_tcga_discovered_recurrent_event_feature_tests.csv"), row.names = FALSE)

significant_overlap <- tcga_discovery_tests |>
  filter(!is.na(.data$p_adj_by_feature_type), .data$p_adj_by_feature_type < 0.10) |>
  arrange(.data$p_adj_by_feature_type, desc(abs(.data$cliffs_delta))) |>
  select(.data$event_id, .data$event_label, .data$direction, .data$arm, .data$feature_type,
         .data$feature, .data$feature_label, .data$n_event, .data$n_no_event,
         .data$cliffs_delta, .data$delta_median, .data$p_value, .data$p_adj_by_feature_type,
         .data$sc_is_recurrent, .data$sc_frac_samples_event, .data$known_genes, .data$annotation_sources)
write.csv(significant_overlap, file.path(tiers[["tables"]], "Auto_tcga_discovered_significant_event_overlap_scRef.csv"), row.names = FALSE)

if (file.exists(sc_event_tests_path)) {
  sc_event_tests <- read.csv(sc_event_tests_path, check.names = FALSE)
  state_feature_lookup <- c(
    "state__Classic_proliferation" = "Classic proliferation",
    "state__Basal_to_intestinal_metaplasia" = "Basal to intestinal metaplasia",
    "state__SMG_to_intestinal_metaplasia" = "SMG to intestinal metaplasia",
    "state__Stress_adaptive" = "Stress adaptive",
    "state__Cancer_cell_immune_mimicry" = "Cancer-cell immune mimicry"
  )
  sc_event_tests <- sc_event_tests |>
    mutate(
      tcga_feature = case_when(
        grepl("^mp__", .data$feature) ~ sub("^mp__", "", .data$feature),
        .data$feature %in% names(state_feature_lookup) ~ unname(state_feature_lookup[.data$feature]),
        TRUE ~ .data$feature
      )
    )
  sc_tcga_concordance <- sc_validation$tests |>
    select(.data$event_id, .data$feature, .data$feature_type, .data$feature_label,
           tcga_cliffs_delta = .data$cliffs_delta, tcga_delta_median = .data$delta_median,
           tcga_primary_delta = .data$primary_delta,
           tcga_p_value = .data$p_value, tcga_p_adj_by_feature_type = .data$p_adj_by_feature_type) |>
    inner_join(
      sc_event_tests |>
        select(.data$event_id, feature = .data$tcga_feature,
               sc_primary_delta = .data$primary_delta,
               sc_primary_p_adj_group = .data$primary_p_adj_group,
               sc_unpaired_p_adj_group = .data$unpaired_p_adj_group),
      by = c("event_id", "feature")
    ) |>
    mutate(
      direction_concordant = case_when(
        !is.finite(.data$tcga_cliffs_delta) | !is.finite(.data$sc_primary_delta) ~ NA,
        sign(.data$tcga_cliffs_delta) == sign(.data$sc_primary_delta) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  write.csv(sc_tcga_concordance, file.path(tiers[["tables"]], "Auto_scRef_tcga_event_feature_concordance.csv"), row.names = FALSE)
} else {
  sc_tcga_concordance <- data.frame()
}

####################
# 6b) Rectangular MP-by-CNA validation heatmaps
####################
all_arm_fdr_cutoff <- 0.05
recurrent_fdr_cutoff <- 0.05

sc_results_path_for_arm_mp <- file.path(sc_cna_dir, "rds", "Auto_cna_subclone_expression_results.rds")

safe_spearman <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(list(rho = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  cr <- suppressWarnings(tryCatch(cor.test(x[ok], y[ok], method = "spearman"), error = function(e) NULL))
  if (is.null(cr)) return(list(rho = NA_real_, p_value = NA_real_, n = sum(ok)))
  list(rho = unname(cr$estimate), p_value = cr$p.value, n = sum(ok))
}

compute_sc_arm_mp_spearman <- function(sc_results_obj, arm_order, mp_order_for_tests) {
  sc_arm_matrix <- as.matrix(sc_results_obj$arm_matrix)
  sc_arm_matrix[!is.finite(sc_arm_matrix)] <- 0
  sc_features_for_tests <- as.data.frame(sc_results_obj$features)
  if (!"subclone_id" %in% colnames(sc_features_for_tests) &&
      all(c("sample", "subclone") %in% colnames(sc_features_for_tests))) {
    sc_features_for_tests <- sc_features_for_tests |>
      mutate(subclone_id = paste(.data$sample, .data$subclone, sep = "::"))
  }
  ids <- intersect(rownames(sc_arm_matrix), sc_features_for_tests$subclone_id)
  bind_rows(lapply(arm_order, function(arm_name) {
    if (!arm_name %in% colnames(sc_arm_matrix)) return(NULL)
    bind_rows(lapply(mp_order_for_tests, function(mp_name) {
      mp_col <- paste0("mp__", mp_name)
      if (!mp_col %in% colnames(sc_features_for_tests)) return(NULL)
      mp_val <- sc_features_for_tests[[mp_col]][match(ids, sc_features_for_tests$subclone_id)]
      res <- safe_spearman(sc_arm_matrix[ids, arm_name], mp_val)
      data.frame(
        dataset = "scRef subclones",
        arm = arm_name,
        mp = mp_name,
        rho = res$rho,
        p_value = res$p_value,
        n = res$n,
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    mutate(
      p_adj = p.adjust(.data$p_value, method = "BH"),
      sig_fdr = !is.na(.data$p_adj) & .data$p_adj < all_arm_fdr_cutoff,
      sig_label = p_to_stars(.data$p_adj)
    )
}

compute_tcga_arm_mp_spearman <- function(arm_long_df, score_df, arm_order, mp_order_for_tests) {
  bind_rows(lapply(arm_order, function(arm_name) {
    arm_df <- arm_long_df |>
      filter(.data$arm == arm_name) |>
      select(.data$sample_barcode, .data$sample_key, .data$arm_mean)
    bind_rows(lapply(mp_order_for_tests, function(mp_name) {
      if (!mp_name %in% colnames(score_df)) return(NULL)
      d <- arm_df |>
        inner_join(score_df |> select(.data$sample_barcode, .data$sample_key, all_of(mp_name)),
                   by = c("sample_barcode", "sample_key"))
      res <- safe_spearman(d$arm_mean, d[[mp_name]])
      data.frame(
        dataset = "TCGA bulk",
        arm = arm_name,
        mp = mp_name,
        rho = res$rho,
        p_value = res$p_value,
        n = res$n,
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    mutate(
      p_adj = p.adjust(.data$p_value, method = "BH"),
      sig_fdr = !is.na(.data$p_adj) & .data$p_adj < all_arm_fdr_cutoff,
      sig_label = p_to_stars(.data$p_adj)
    )
}

make_concordance_flags <- function(df, sc_prefix, tcga_prefix, fdr_cutoff) {
  df |>
    mutate(
      sc_sig = !is.na(.data[[paste0(sc_prefix, "_p_adj")]]) & .data[[paste0(sc_prefix, "_p_adj")]] < fdr_cutoff,
      tcga_sig = !is.na(.data[[paste0(tcga_prefix, "_p_adj")]]) & .data[[paste0(tcga_prefix, "_p_adj")]] < fdr_cutoff,
      same_trend = is.finite(.data[[paste0(sc_prefix, "_effect")]]) &
        is.finite(.data[[paste0(tcga_prefix, "_effect")]]) &
        .data[[paste0(sc_prefix, "_effect")]] != 0 &
        .data[[paste0(tcga_prefix, "_effect")]] != 0 &
        sign(.data[[paste0(sc_prefix, "_effect")]]) == sign(.data[[paste0(tcga_prefix, "_effect")]]),
      same_trend_sig_any = .data$same_trend & (.data$sc_sig | .data$tcga_sig),
      same_trend_sig_both = .data$same_trend & .data$sc_sig & .data$tcga_sig
    )
}

build_all_arm_tile_plot <- function(plot_df, title_text, fill_title,
                                    x_text_size = 15, y_text_size = 14,
                                    base_size = 20) {
  max_abs <- max(abs(plot_df$effect_for_fill), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  if (!"star_label" %in% colnames(plot_df)) {
      plot_df <- plot_df |>
        mutate(
          star_label = ifelse(!is.na(.data$p_adj) & .data$p_adj < 0.05, "\u2605", "")
        )
  }

  ggplot(plot_df, aes(x = .data$arm, y = .data$mp_label)) +
    geom_tile(aes(fill = .data$effect_for_fill), colour = "grey85", linewidth = 0.4, width = 1, height = 1) +
    geom_text(
      data = plot_df |> filter(.data$star_label != ""),
      aes(label = .data$star_label), colour = "black", size = 8, vjust = 0.4, hjust = 0.5
    ) +
    facet_grid(. ~ dataset) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-0.25, 0.25), oob = scales::squish, na.value = "white", 
      name = fill_title
    ) +
    labs(title = title_text, x = "Chromosome arm level CNA events", y = "MP expression") +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = x_text_size, colour = "black"),
      axis.text.y = element_text(size = y_text_size, colour = "black"),
      axis.title.x = element_text(size = base_size + 1, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = base_size + 1, face = "bold", margin = margin(r = 8)),
      strip.text = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      panel.spacing.x = grid::unit(1.5, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank()
    )
}

if (file.exists(sc_results_path_for_arm_mp)) {
  message("Computing all-arm continuous CNA vs MP Spearman validation")
  sc_results_for_arm_mp <- readRDS(sc_results_path_for_arm_mp)
  arm_mp_features <- intersect(mp_order, mp_features)
  sc_all_arm_mp_tests <- compute_sc_arm_mp_spearman(sc_results_for_arm_mp, arm_levels, arm_mp_features)
  tcga_all_arm_mp_tests <- compute_tcga_arm_mp_spearman(tcga_arm_long, tcga_score_df, arm_levels, arm_mp_features)

  all_arm_concordance <- sc_all_arm_mp_tests |>
    select(.data$arm, .data$mp, sc_effect = .data$rho, sc_p_value = .data$p_value,
           sc_p_adj = .data$p_adj, sc_n = .data$n, sc_sig_fdr = .data$sig_fdr) |>
    inner_join(
      tcga_all_arm_mp_tests |>
        select(.data$arm, .data$mp, tcga_effect = .data$rho, tcga_p_value = .data$p_value,
               tcga_p_adj = .data$p_adj, tcga_n = .data$n, tcga_sig_fdr = .data$sig_fdr),
      by = c("arm", "mp")
    ) |>
    make_concordance_flags("sc", "tcga", all_arm_fdr_cutoff) |>
    mutate(
      mp_label = feature_label_map[.data$mp],
      mp_label = ifelse(is.na(.data$mp_label), .data$mp, .data$mp_label)
    )

  write.csv(sc_all_arm_mp_tests, file.path(tiers[["tables"]], "Auto_scRef_all_arm_mp_spearman_tests.csv"), row.names = FALSE)
  write.csv(tcga_all_arm_mp_tests, file.path(tiers[["tables"]], "Auto_tcga_all_arm_mp_spearman_tests.csv"), row.names = FALSE)
  write.csv(all_arm_concordance, file.path(tiers[["tables"]], "Auto_scRef_tcga_all_arm_mp_spearman_concordance.csv"), row.names = FALSE)

  all_arm_plot_df <- bind_rows(
    sc_all_arm_mp_tests |> transmute(dataset = "scATLAS subclones", arm, mp, effect = rho, p_adj, sig_fdr),
    tcga_all_arm_mp_tests |> transmute(dataset = "OAC TCGA dataset", arm, mp, effect = rho, p_adj, sig_fdr)
  ) |>
    left_join(
      all_arm_concordance |>
        select(.data$arm, .data$mp, sc_sig = .data$sc_sig_fdr, tcga_sig = .data$tcga_sig_fdr),
      by = c("arm", "mp")
    ) |>
    mutate(
      dataset = factor(.data$dataset, levels = c("scATLAS subclones", "OAC TCGA dataset")),
      arm = factor(.data$arm, levels = arm_levels),
      mp_label = feature_label_map[.data$mp],
      mp_label = ifelse(is.na(.data$mp_label), .data$mp, .data$mp_label),
      mp_label = factor(.data$mp_label, levels = feature_label_map[arm_mp_features]),
      effect_for_fill = case_when(
        dataset == "scATLAS subclones" & sc_sig ~ effect,
        dataset == "OAC TCGA dataset" & sc_sig ~ effect,
        TRUE ~ NA_real_
      ),
      star_label = case_when(
        dataset == "scATLAS subclones" & sc_sig ~ "\u2605",
        dataset == "OAC TCGA dataset" & sc_sig & tcga_sig ~ "\u2605",
        TRUE ~ ""
      )
    )

  p_all_arm_rectangles <- build_all_arm_tile_plot(
    all_arm_plot_df,
    "Chromosome arm level CNA events vs MP expression",
    "Spearman rho",
    x_text_size = 15,
    y_text_size = 14,
    base_size = 20
  )
  ggsave(
    file.path(tiers[["figures"]], "Auto_scRef_tcga_all_arm_mp_spearman_rectangles.pdf"),
    p_all_arm_rectangles, width = 25, height = 10, device = cairo_pdf
  )

  all_arm_all_trends_plot_df <- bind_rows(
    sc_all_arm_mp_tests |> transmute(dataset = "scATLAS subclones", arm, mp, effect = rho, p_adj, sig_fdr),
    tcga_all_arm_mp_tests |> transmute(dataset = "OAC TCGA dataset", arm, mp, effect = rho, p_adj, sig_fdr)
  ) |>
    left_join(
      all_arm_concordance |>
        select(.data$arm, .data$mp, sc_sig = .data$sc_sig_fdr, tcga_sig = .data$tcga_sig_fdr),
      by = c("arm", "mp")
    ) |>
    mutate(
      dataset = factor(.data$dataset, levels = c("scATLAS subclones", "OAC TCGA dataset")),
      arm = factor(.data$arm, levels = arm_levels),
      mp_label = feature_label_map[.data$mp],
      mp_label = ifelse(is.na(.data$mp_label), .data$mp, .data$mp_label),
      mp_label = factor(.data$mp_label, levels = feature_label_map[arm_mp_features]),
      effect_for_fill = .data$effect,
      star_label = ifelse(.data$sig_fdr, "\u2605", "")
    )

  p_all_arm_all_trends <- build_all_arm_tile_plot(
    all_arm_all_trends_plot_df,
    "Chromosome arm level CNA events vs MP expression: all trends",
    "Spearman rho",
    x_text_size = 15,
    y_text_size = 14,
    base_size = 20
  )
  ggsave(
    file.path(tiers[["figures"]], "Auto_scRef_tcga_all_arm_mp_spearman_rectangles_all_trends.pdf"),
    p_all_arm_all_trends, width = 25, height = 10, device = cairo_pdf
  )

  if (file.exists(sc_event_tests_path)) {
    message("Creating rectangular recurrent-event MP validation heatmap")
    sc_event_tests_for_rectangles <- read.csv(sc_event_tests_path, check.names = FALSE) |>
      mutate(
        tcga_feature = case_when(
          grepl("^mp__", .data$feature) ~ sub("^mp__", "", .data$feature),
          TRUE ~ .data$feature
        )
      )

    sc_recurrent_mp_rect <- sc_event_tests_for_rectangles |>
      filter(.data$event_id %in% sc_recurrent$event_id, .data$feature_group == "Metaprogrammes") |>
      transmute(
        event_id = .data$event_id,
        event_label = .data$event_label,
        arm = .data$arm,
        direction = .data$direction,
        mp = .data$tcga_feature,
        sc_effect = .data$primary_delta,
        sc_p_value = .data$unpaired_p_value,
        sc_p_adj = .data$primary_p_adj_group
      )

    tcga_recurrent_mp_rect <- sc_validation$tests |>
      filter(.data$feature_group == "Metaprogrammes") |>
      transmute(
        event_id = .data$event_id,
        event_label = as.character(.data$event_label),
        arm = .data$arm,
        direction = .data$direction,
        mp = .data$feature,
        tcga_effect = .data$primary_delta,
        tcga_p_value = .data$p_value,
        tcga_p_adj = .data$primary_p_adj_group
      )

    recurrent_mp_concordance <- sc_recurrent_mp_rect |>
      inner_join(tcga_recurrent_mp_rect, by = c("event_id", "arm", "direction", "mp"), suffix = c("_sc", "_tcga")) |>
      mutate(event_label = ifelse(!is.na(.data$event_label_tcga), .data$event_label_tcga, .data$event_label_sc)) |>
      select(.data$event_id, .data$event_label, .data$arm, .data$direction, .data$mp,
             .data$sc_effect, .data$sc_p_value, .data$sc_p_adj,
             .data$tcga_effect, .data$tcga_p_value, .data$tcga_p_adj) |>
      make_concordance_flags("sc", "tcga", recurrent_fdr_cutoff)

    write.csv(recurrent_mp_concordance, file.path(tiers[["tables"]], "Auto_scRef_tcga_recurrent_event_mp_concordance_rectangles.csv"), row.names = FALSE)

    recurrent_event_levels <- sc_recurrent |>
      mutate(event_label_for_rectangles = event_label(.data$event_id, .data$known_genes)) |>
      pull(.data$event_label_for_rectangles)

    recurrent_plot_df <- bind_rows(
      recurrent_mp_concordance |>
        transmute(dataset = "scATLAS subclones", event_id, event_label, mp, effect = sc_effect,
                  p_adj = sc_p_adj, sig_fdr = sc_sig, sc_sig = sc_sig, tcga_sig = tcga_sig),
      recurrent_mp_concordance |>
        transmute(dataset = "OAC TCGA dataset", event_id, event_label, mp, effect = tcga_effect,
                  p_adj = tcga_p_adj, sig_fdr = tcga_sig, sc_sig = sc_sig, tcga_sig = tcga_sig)
    ) |>
      mutate(
        dataset = factor(.data$dataset, levels = c("scATLAS subclones", "OAC TCGA dataset")),
        arm = factor(.data$event_label, levels = recurrent_event_levels),
        mp_label = feature_label_map[.data$mp],
        mp_label = ifelse(is.na(.data$mp_label), .data$mp, .data$mp_label),
        mp_label = factor(.data$mp_label, levels = feature_label_map[arm_mp_features]),
        effect_for_fill = case_when(
          dataset == "scATLAS subclones" & sc_sig ~ effect,
          dataset == "OAC TCGA dataset" & sc_sig ~ effect,
          TRUE ~ NA_real_
        ),
        star_label = case_when(
          dataset == "scATLAS subclones" & sc_sig ~ "\u2605",
          dataset == "OAC TCGA dataset" & sc_sig & tcga_sig ~ "\u2605",
          TRUE ~ ""
        )
      )

    p_recurrent_rectangles <- build_all_arm_tile_plot(
      recurrent_plot_df,
      "Top recurrent scRef CNA event MP associations and TCGA validation",
      "Standardized\nevent delta",
      x_text_size = 15,
      y_text_size = 14,
      base_size = 20
    ) + labs(x = "Recurrent chromosome-arm event")

    ggsave(
      file.path(tiers[["figures"]], "Auto_scRef_tcga_recurrent_event_mp_rectangles.pdf"),
      p_recurrent_rectangles, width = 25, height = 10, device = cairo_pdf
    )
  } else {
    recurrent_mp_concordance <- data.frame()
    message("Skipping recurrent rectangular heatmap because scRef recurrent event tests are missing: ", sc_event_tests_path)
  }

  rectangle_summary <- bind_rows(
    data.frame(
      summary_type = "all_arm_continuous",
      item = "scRef_fdr_lt_0.05",
      value = sum(sc_all_arm_mp_tests$sig_fdr, na.rm = TRUE),
      detail = paste0(nrow(sc_all_arm_mp_tests), " scRef arm-MP Spearman tests"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      summary_type = "all_arm_continuous",
      item = "tcga_fdr_lt_0.05",
      value = sum(tcga_all_arm_mp_tests$sig_fdr, na.rm = TRUE),
      detail = paste0(nrow(tcga_all_arm_mp_tests), " TCGA arm-MP Spearman tests"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      summary_type = "all_arm_continuous",
      item = "same_trend_significant_either_dataset",
      value = sum(all_arm_concordance$same_trend_sig_any, na.rm = TRUE),
      detail = paste0(sum(all_arm_concordance$same_trend_sig_both, na.rm = TRUE), " significant in both datasets"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      summary_type = "recurrent_rectangles",
      item = "same_trend_significant_either_dataset",
      value = if (exists("recurrent_mp_concordance")) sum(recurrent_mp_concordance$same_trend_sig_any, na.rm = TRUE) else NA_real_,
      detail = if (exists("recurrent_mp_concordance")) paste0(sum(recurrent_mp_concordance$same_trend_sig_both, na.rm = TRUE), " significant in both datasets") else "not computed",
      stringsAsFactors = FALSE
    )
  )
  write.csv(rectangle_summary, file.path(tiers[["tables"]], "Auto_tcga_cna_rectangle_heatmap_summary.csv"), row.names = FALSE)
} else {
  message("Skipping all-arm MP association heatmap because scRef result cache is missing: ", sc_results_path_for_arm_mp)
}
####################

validation_conclusion <- data.frame(
  conclusion_item = c(
    "event_universe",
    "threshold_definition",
    "validated_trend",
    "tcga_discovery_overlap"
  ),
  conclusion = c(
    paste0("Validated the same top ", n_distinct(sc_validation$tests$event_id),
           " scRef recurrent CNA events in the TCGA dotplot and the top ",
           length(sc_boxplot_events), " in the boxplot."),
    paste0("Selected threshold ", selected_threshold,
           " by maximizing F1 = 2*precision*recall/(precision+recall), where precision = curated OAC events recovered / TCGA recurrent events and recall = curated OAC events recovered / curated OAC events."),
    if (sum(sc_validation$tests$p_adj_by_feature_type < 0.10, na.rm = TRUE) == 0) {
      "No TCGA MP/state association reached feature-group FDR < 0.10 for the top scRef recurrent CNA events; this validates the scRef conclusion that recurrent arm events are not robustly coupled to MP/state programmes."
    } else {
      "At least one top scRef recurrent CNA event showed a TCGA MP/state association at feature-group FDR < 0.10; inspect Auto_scRef_recurrent_events_tcga_feature_tests.csv for direction and effect size."
    },
    if (nrow(significant_overlap) == 0) {
      "No TCGA-discovered recurrent event-feature associations reached feature-group FDR < 0.10."
    } else {
      paste0(nrow(significant_overlap), " TCGA-discovered event-feature association(s) reached feature-group FDR < 0.10; ",
             sum(significant_overlap$sc_is_recurrent, na.rm = TRUE), " were recurrent in scRef.")
    }
  ),
  stringsAsFactors = FALSE
)
write.csv(validation_conclusion, file.path(tiers[["tables"]], "Auto_tcga_cna_validation_conclusion.csv"), row.names = FALSE)

####################
# 7) Figures and summaries
####################
event_meta <- sc_recurrent |>
  mutate(
    event_label = event_label(.data$event_id, .data$known_genes),
    event_label = factor(.data$event_label, levels = event_label(.data$event_id, .data$known_genes))
  )
event_order <- event_meta$event_label

sc_validation$tests <- sc_validation$tests |>
  mutate(
    event_label = factor(event_label(.data$event_id, .data$known_genes), levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label_map[plot_features])))
  )
sc_validation$values <- sc_validation$values |>
  mutate(
    event_label = factor(event_label(.data$event_id, .data$known_genes), levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label_map[plot_features]))),
    event_group = factor(ifelse(.data$event_present, "Event-positive", "Event-negative"),
                         levels = c("Event-negative", "Event-positive"))
  )

event_bar <- event_meta |>
  mutate(event_label = factor(.data$event_label, levels = rev(as.character(.data$event_label)))) |>
  ggplot(aes(.data$frac_samples_event, .data$event_label, fill = .data$direction)) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(aes(label = percent(.data$frac_samples_event, accuracy = 1)), hjust = -0.08, size = 8) +
  scale_x_continuous(labels = percent, limits = c(0, min(1, max(event_meta$frac_samples_event) * 1.18))) +
  scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
  labs(title = "Top scRef recurrent arm-level CNA events validated in TCGA",
       x = "Fraction of samples with event in at least one subclone", y = NULL, fill = NULL) +
  theme_classic(base_size = 26) +
  theme(plot.title = element_text(face = "bold", size = 30),
        axis.text = element_text(size = 22),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12)),
        legend.position = "top",
        legend.text = element_text(size = 20))

plot_event_assoc <- function(group_name, title_suffix) {
  group_features <- plot_features[feature_group[plot_features] == group_name]
  df <- sc_validation$tests |>
    filter(.data$feature_group == group_name) |>
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label_map[group_features])))
  ggplot(df, aes(.data$event_label, .data$feature_label)) +
    geom_point(aes(size = .data$neglog10_fdr, fill = .data$primary_delta),
               shape = 21, color = "black", stroke = 0.45, alpha = 0.95) +
    geom_text(aes(label = .data$sig_label), size = 6.2, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         na.value = "grey90") +
    scale_size_continuous(range = c(3.5, 11), limits = c(0, 12)) +
    labs(
      title = paste0("TCGA validation of scRef recurrent CNA events: ", title_suffix),
      x = NULL,
      y = NULL,
      fill = "Median standardized\nevent delta",
      size = "-log10(FDR)"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15)
    )
}

event_page_chunks <- function(event_ids, page_size = 4L) {
  if (length(event_ids) == 0) return(list(character(0)))
  split(event_ids, ceiling(seq_along(event_ids) / page_size))
}

plot_event_boxplots <- function(group_name, title_suffix, event_ids) {
  group_features <- plot_features[feature_group[plot_features] == group_name]
  df <- sc_validation$values |>
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) |>
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label_map[group_features]),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])))
  sig_df <- sc_validation$tests |>
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) |>
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label_map[group_features]),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)]))) |>
    left_join(
      df |>
        group_by(.data$event_id, .data$feature) |>
        summarise(star_y = max(.data$feature_z, na.rm = TRUE) + 0.25, .groups = "drop"),
      by = c("event_id", "feature")
    )
  ggplot(df, aes(x = .data$feature_label, y = .data$feature_z, fill = .data$event_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.88,
                 position = position_dodge(width = 0.72), linewidth = 0.55) +
    geom_point(aes(color = .data$event_group),
               position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0.02,
                                               dodge.width = 0.72),
               alpha = 0.36, size = 1.3, show.legend = FALSE) +
    geom_text(data = sig_df, aes(x = .data$feature_label, y = .data$star_y, label = .data$sig_label),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~event_label, ncol = 2) +
    scale_fill_manual(values = c("Event-negative" = "grey72", "Event-positive" = "#B2182B")) +
    scale_color_manual(values = c("Event-negative" = "grey45", "Event-positive" = "#B2182B")) +
    labs(
      title = paste0("TCGA feature distributions by scRef recurrent CNA event: ", title_suffix),
      x = NULL,
      y = "Standardized TCGA GSVA score",
      fill = NULL
    ) +
    theme_classic(base_size = 20) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(face = "bold", size = 24),
      legend.position = "top",
      legend.text = element_text(size = 16)
    )
}

pdf(file.path(tiers[["figures"]], "Auto_tcga_cna_event_association_dotplots.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(event_bar)
print(plot_event_assoc("Metaprogrammes", "all metaprogrammes"))
print(plot_event_assoc("Centred states", "five centred states excluding Hybrid and Unresolved"))
dev.off()

pdf(file.path(tiers[["figures"]], "Auto_tcga_cna_event_boxplots.pdf"),
    width = 22, height = 14, useDingbats = FALSE)
for (event_page in event_page_chunks(sc_boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Metaprogrammes", "all metaprogrammes", event_page))
}
for (event_page in event_page_chunks(sc_boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Centred states", "five centred states excluding Hybrid and Unresolved", event_page))
}
dev.off()

p_threshold <- threshold_optimization |>
  arrange(.data$threshold) |>
  ggplot(aes(.data$threshold, .data$f1)) +
  geom_line(linewidth = 0.7, colour = "#2C7BB6") +
  geom_point(aes(fill = .data$selected), shape = 21, size = 2.6, colour = "black") +
  geom_vline(xintercept = selected_threshold, linetype = "dashed", colour = "#B2182B") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#B2182B"), guide = "none") +
  labs(title = "TCGA arm-call threshold selection", x = "Absolute arm mean threshold", y = "F1 vs curated OAC arm events") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(colour = "black"))
ggsave(file.path(tiers[["figures"]], "Auto_tcga_cna_threshold_optimization.pdf"),
       p_threshold, width = 8, height = 5, device = cairo_pdf)

summary_df <- bind_rows(
  data.frame(
    summary_type = "threshold",
    item = "selected_abs_arm_mean_threshold",
    value = selected_threshold,
    detail = paste0("F1=", signif(threshold_optimization$f1[threshold_optimization$selected][1], 3),
                    "; recurrent events=", threshold_optimization$n_tcga_recurrent_events[threshold_optimization$selected][1]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "cohort",
    item = "tcga_eac_primary_samples_with_segments_and_scores",
    value = length(intersect(unique(tcga_arm_calls$sample_barcode), unique(score_long$sample_barcode))),
    detail = "Matched by first 15 TCGA barcode characters between segment ID and RNA-seq sample_barcode",
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "scRef_validation",
    item = "scRef_recurrent_event_feature_tests_fdr_lt_0.10",
    value = sum(sc_validation$tests$p_adj_by_feature_type < 0.10, na.rm = TRUE),
    detail = paste0(n_distinct(sc_validation$tests$event_id), " scRef recurrent events tested"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "tcga_discovery",
    item = "tcga_discovered_event_feature_tests_fdr_lt_0.10",
    value = sum(tcga_discovery_tests$p_adj_by_feature_type < 0.10, na.rm = TRUE),
    detail = paste0(n_distinct(tcga_discovery_tests$event_id), " TCGA recurrent events tested"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "tcga_discovery",
    item = "significant_tcga_pairs_with_scRef_recurrent_event",
    value = sum(significant_overlap$sc_is_recurrent, na.rm = TRUE),
    detail = paste0(nrow(significant_overlap), " TCGA event-feature pairs at feature-type FDR < 0.10"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "annotation_check",
    item = "occams_st5_myc_8q_driver_present",
    value = as.numeric(any(occams_drivers$sheet == "ST5 High confidence Amp Drivers" &
                            occams_drivers$gene == "MYC" &
                            occams_drivers$arm == "chr8q")),
    detail = paste(occams_drivers |>
                     filter(.data$sheet == "ST5 High confidence Amp Drivers", .data$gene == "MYC") |>
                     transmute(detail = paste0(.data$ensembl_id, "/", .data$gene_raw, " at ", .data$cytoband, " ", .data$peak)) |>
                     pull(.data$detail),
                   collapse = "; "),
    stringsAsFactors = FALSE
  )
)

write.csv(summary_df, file.path(tiers[["tables"]], "Auto_tcga_cna_recurrent_event_validation_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(summary_dir, "Auto_tcga_cna_recurrent_event_validation_summary.csv"), row.names = FALSE)

####################
# 8) Consensus Heatmap Comparison & 8q/MYC Boxplots
####################
message("Generating consensus heatmap comparison and 8q/MYC boxplots")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(gridExtra)
  library(cluster)
})

# Load SC data
sc_results_path <- file.path(sc_cna_dir, "rds", "Auto_cna_subclone_expression_results.rds")
if (file.exists(sc_results_path)) {
  sc_results <- readRDS(sc_results_path)
  
  sc_arm_matrix_valid <- as.matrix(sc_results$arm_matrix)
  sc_arm_matrix_valid[!is.finite(sc_arm_matrix_valid)] <- 0
  valid_sc_rows <- rownames(sc_arm_matrix_valid)[rowSums(is.finite(sc_arm_matrix_valid)) == ncol(sc_arm_matrix_valid)]
  sc_arm_matrix_valid <- sc_arm_matrix_valid[valid_sc_rows, , drop = FALSE]
  
  sc_cna_cluster <- sc_results$cna_cluster
  sc_features <- as.data.frame(sc_results$features)
  
  sc_hc <- hclust(dist(sc_arm_matrix_valid), method = "ward.D2")
  sc_cluster_df <- data.frame(
    subclone_id = names(sc_cna_cluster),
    cna_cluster = unname(sc_cna_cluster),
    stringsAsFactors = FALSE
  ) |>
    left_join(sc_features |> select(.data$subclone_id, .data$sample, .data$subclone,
                                  .data$subclone_fraction, .data$dominance_class),
              by = "subclone_id")
  sc_cluster_df <- sc_cluster_df[match(rownames(sc_arm_matrix_valid), sc_cluster_df$subclone_id), , drop = FALSE]
  
  sc_cluster_cols <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(unique(sc_cluster_df$cna_cluster))),
    sort(unique(sc_cluster_df$cna_cluster))
  )
  sc_dominance_cols <- c(
    "single_subclone" = "grey70",
    "largest_not_significant" = "#FDB863",
    "significant_dominant" = "#B2182B"
  )
  
  chr_from_arm <- sub("[pq]$", "", arm_levels)
  
  plot_sc_consensus_heatmap <- function() {
    row_meta <- sc_cluster_df[sc_hc$order, , drop = FALSE]
    mat <- sc_arm_matrix_valid[sc_hc$order, arm_levels, drop = FALSE]
    
    row_ha <- rowAnnotation(
      Cluster = row_meta$cna_cluster,
      col = list(
        Cluster = sc_cluster_cols
      ),
      annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        Cluster = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12))
      ),
      simple_anno_size = unit(5, "mm")
    )
    Heatmap(
      mat,
      name = "Mean SC CNA",
      col = colorRamp2(c(-0.18, 0, 0.18), c("#2166AC", "white", "#B2182B")),
      height = unit(3.5, "inches"),
      left_annotation = row_ha,
      row_split = factor(row_meta$cna_cluster, levels = sort(unique(row_meta$cna_cluster))),
      row_title = NULL,
      row_title_gp = gpar(fontsize = 0, col = NA),
      row_gap = unit(1.2, "mm"),
      column_split = factor(chr_from_arm, levels = chrom_levels),
      column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 12),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      heatmap_legend_param = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      use_raster = TRUE,
      raster_quality = 4,
      border = FALSE,
      rect_gp = gpar(col = NA),
      column_title = "scATLAS Subclones"
    )
  }
  
  # TCGA matrix
  tcga_arm_wide <- tcga_arm_long |>
    select(.data$sample_barcode, .data$arm, .data$arm_mean) |>
    pivot_wider(names_from = .data$arm, values_from = .data$arm_mean, values_fill = 0) |>
    as.data.frame()
  rownames(tcga_arm_wide) <- tcga_arm_wide$sample_barcode
  tcga_arm_matrix <- as.matrix(tcga_arm_wide[, arm_levels, drop = FALSE])
  tcga_arm_matrix[!is.finite(tcga_arm_matrix)] <- 0
  
  tcga_hc <- hclust(dist(tcga_arm_matrix), method = "ward.D2")
  
  choose_consensus_k <- function(mat, max_k = 6L) {
    if (nrow(mat) < 4) return(1L)
    max_k <- min(max_k, nrow(mat) - 1L)
    d <- dist(mat)
    scores <- lapply(2:max_k, function(k) {
      cl <- cutree(hclust(d, method = "ward.D2"), k = k)
      if (any(table(cl) < 2)) return(data.frame(k = k, silhouette = NA_real_))
      sil <- tryCatch(mean(cluster::silhouette(cl, d)[, "sil_width"]), error = function(e) NA_real_)
      data.frame(k = k, silhouette = sil)
    }) |> bind_rows()
    scores <- scores |> filter(is.finite(.data$silhouette))
    if (nrow(scores) == 0) return(min(3L, nrow(mat)))
    scores$k[which.max(scores$silhouette)]
  }
  
  tcga_best_k <- choose_consensus_k(tcga_arm_matrix, max_k = 6L)
  tcga_cna_cluster <- if (tcga_best_k <= 1L) {
    setNames(rep("CNA cluster 1", nrow(tcga_arm_matrix)), rownames(tcga_arm_matrix))
  } else {
    setNames(paste0("CNA cluster ", cutree(tcga_hc, k = tcga_best_k)), rownames(tcga_arm_matrix))
  }
  
  tcga_cluster_df <- data.frame(
    sample_barcode = names(tcga_cna_cluster),
    cna_cluster = unname(tcga_cna_cluster),
    stringsAsFactors = FALSE
  )
  
  tcga_cluster_cols <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(unique(tcga_cluster_df$cna_cluster))),
    sort(unique(tcga_cluster_df$cna_cluster))
  )
  
  plot_tcga_consensus_heatmap <- function() {
    row_meta <- tcga_cluster_df[tcga_hc$order, , drop = FALSE]
    mat <- tcga_arm_matrix[tcga_hc$order, arm_levels, drop = FALSE]
    
    row_ha <- rowAnnotation(
      Cluster = row_meta$cna_cluster,
      col = list(
        Cluster = tcga_cluster_cols
      ),
      annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
      show_annotation_name = TRUE,
      annotation_name_side = "bottom",
      annotation_legend_param = list(
        Cluster = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12))
      ),
      simple_anno_size = unit(5, "mm")
    )
    Heatmap(
      mat,
      name = "Mean TCGA CNA",
      col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
      height = unit(3.5, "inches"),
      left_annotation = row_ha,
      row_split = factor(row_meta$cna_cluster, levels = sort(unique(row_meta$cna_cluster))),
      row_title = NULL,
      row_title_gp = gpar(fontsize = 0, col = NA),
      row_gap = unit(1.2, "mm"),
      column_split = factor(chr_from_arm, levels = chrom_levels),
      column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 12),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      heatmap_legend_param = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      use_raster = TRUE,
      raster_quality = 4,
      border = FALSE,
      rect_gp = gpar(col = NA),
      column_title = "TCGA Patients"
    )
  }
  
  pdf(file.path(tiers[["figures"]], "Auto_tcga_and_scRef_consensus_cna_heatmap.pdf"),
      width = 16, height = 9, useDingbats = FALSE)
  draw(plot_sc_consensus_heatmap() %v% plot_tcga_consensus_heatmap(), ht_gap = unit(0.8, "cm"), heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  # 8q / MYC Boxplots
  sc_arm_long <- as.data.frame(sc_results$arm_long) |>
    filter(.data$subclone != "Unresolved") |>
    mutate(
      cna_call = case_when(
        .data$arm_mean >= 0.05 ~ 1L,
        .data$arm_mean <= -0.05 ~ -1L,
        TRUE ~ 0L
      )
    )
  sc_chr8q_myc <- sc_arm_long |>
    filter(.data$arm == "chr8q") |>
    select(.data$sample, .data$subclone, .data$subclone_id, .data$cna_call) |>
    left_join(sc_features |> select(.data$sample, .data$subclone, .data$subclone_id, .data[["mp__MP2+"]], .data$n_cells),
              by = c("sample", "subclone", "subclone_id")) |>
    mutate(
      chr8q_group = case_when(
        .data$cna_call == 1L ~ "8q gain",
        .data$cna_call == -1L ~ "8q loss",
        TRUE ~ "No 8q CNA"
      ),
      chr8q_group = factor(.data$chr8q_group, levels = c("8q loss", "No 8q CNA", "8q gain"))
    )
  
  if (nrow(sc_chr8q_myc) > 0) {
    sc_group_n <- sc_chr8q_myc |>
      group_by(.data$chr8q_group) |>
      summarise(n = n(), .groups = "drop")
    
    sc_chr8q_myc <- sc_chr8q_myc |>
      left_join(sc_group_n, by = "chr8q_group") |>
      mutate(chr8q_group_label = paste0(.data$chr8q_group, "\n(n=", .data$n, ")"))
    
    sc_chr8q_myc$chr8q_group_label <- factor(sc_chr8q_myc$chr8q_group_label, 
        levels = unique(sc_chr8q_myc$chr8q_group_label[order(sc_chr8q_myc$chr8q_group)]))
        
    levels_sc <- levels(sc_chr8q_myc$chr8q_group_label)
    comps_sc <- list(c(levels_sc[1], levels_sc[2]), c(levels_sc[2], levels_sc[3]), c(levels_sc[1], levels_sc[3]))
    
    p_gain8q_myc_sc <- ggplot(sc_chr8q_myc, aes(x = .data$chr8q_group_label, y = .data[["mp__MP2+"]], fill = .data$chr8q_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.88, linewidth = 0.8, width = 0.62) +
      geom_point(aes(size = .data$n_cells), position = position_jitter(width = 0.14, height = 0),
                 alpha = 0.42, color = "black") +
      ggpubr::stat_compare_means(comparisons = comps_sc, label = "p.signif", 
                                 symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "NS")),
                                 size = 7, vjust = -0.2, step.increase = 0.12, tip.length = 0.02) +
      scale_fill_manual(values = c("8q loss" = "#2166AC", "No 8q CNA" = "grey72", "8q gain" = "#B2182B")) +
      scale_size_continuous(range = c(1.8, 6)) +
      labs(
        title = "scATLAS: chr8q CNA vs MYC MP",
        x = NULL,
        y = "Subclone mean MP2+ score",
        fill = "chr8q CNA",
        size = "Cells"
      ) +
      theme_classic(base_size = 22) +
      theme(
        plot.title = element_text(face = "bold", size = 26),
        axis.text = element_text(size = 22, color = "black"),
        axis.title = element_text(size = 24, face = "bold", color = "black"),
        legend.position = "none"
      )
  } else {
    p_gain8q_myc_sc <- NULL
  }
  
  tcga_chr8q_myc <- tcga_arm_calls |>
    filter(.data$arm == "chr8q") |>
    select(.data$sample_barcode, .data$sample_key, .data$cna_call) |>
    left_join(tcga_score_df |> select(.data$sample_barcode, .data$sample_key, .data[["MP2+"]]), by = c("sample_barcode", "sample_key")) |>
    mutate(
      chr8q_group = case_when(
        .data$cna_call == 1L ~ "8q gain",
        .data$cna_call == -1L ~ "8q loss",
        TRUE ~ "No 8q CNA"
      ),
      chr8q_group = factor(.data$chr8q_group, levels = c("8q loss", "No 8q CNA", "8q gain"))
    )
  
  if (nrow(tcga_chr8q_myc) > 0) {
    tcga_group_n <- tcga_chr8q_myc |>
      group_by(.data$chr8q_group) |>
      summarise(n = n(), .groups = "drop")
      
    tcga_chr8q_myc <- tcga_chr8q_myc |>
      left_join(tcga_group_n, by = "chr8q_group") |>
      mutate(chr8q_group_label = paste0(.data$chr8q_group, "\n(n=", .data$n, ")"))
      
    tcga_chr8q_myc$chr8q_group_label <- factor(tcga_chr8q_myc$chr8q_group_label, 
        levels = unique(tcga_chr8q_myc$chr8q_group_label[order(tcga_chr8q_myc$chr8q_group)]))
        
    levels_tcga <- levels(tcga_chr8q_myc$chr8q_group_label)
    comps_tcga <- list(c(levels_tcga[1], levels_tcga[2]), c(levels_tcga[2], levels_tcga[3]), c(levels_tcga[1], levels_tcga[3]))
    
    p_gain8q_myc_tcga <- ggplot(tcga_chr8q_myc, aes(x = .data$chr8q_group_label, y = .data[["MP2+"]], fill = .data$chr8q_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.88, linewidth = 0.8, width = 0.62) +
      geom_point(position = position_jitter(width = 0.14, height = 0),
                 alpha = 0.42, color = "black", size = 3) +
      ggpubr::stat_compare_means(comparisons = comps_tcga, label = "p.signif", 
                                 symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "NS")),
                                 size = 7, vjust = -0.2, step.increase = 0.12, tip.length = 0.02) +
      scale_fill_manual(values = c("8q loss" = "#2166AC", "No 8q CNA" = "grey72", "8q gain" = "#B2182B")) +
      labs(
        title = "TCGA: chr8q CNA vs MYC MP",
        x = NULL,
        y = "Patient mean MP2+ score",
        fill = "chr8q CNA"
      ) +
      theme_classic(base_size = 22) +
      theme(
        plot.title = element_text(face = "bold", size = 26),
        axis.text = element_text(size = 22, color = "black"),
        axis.title = element_text(size = 24, face = "bold", color = "black"),
        legend.position = "none"
      )
  } else {
    p_gain8q_myc_tcga <- NULL
  }
  
  if (!is.null(p_gain8q_myc_sc) && !is.null(p_gain8q_myc_tcga)) {
    pdf(file.path(tiers[["figures"]], "Auto_tcga_and_scRef_gain_chr8q_myc_mp.pdf"), width = 16, height = 9, useDingbats = FALSE)
    grid.arrange(p_gain8q_myc_sc, p_gain8q_myc_tcga, ncol = 2)
    dev.off()
  }
} else {
  message("Warning: Auto_cna_subclone_expression_results.rds not found. Skipping comparative plots.")
}

run_summary$parameters$selected_threshold <- selected_threshold
write_run_summary(run_summary, file.path(tiers[["logs"]], "Auto_tcga_cna_recurrent_event_validation_run_summary.rds"))

message("Saved TCGA CNA recurrent event validation outputs to: ", out_dir)
