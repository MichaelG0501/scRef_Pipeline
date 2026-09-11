####################
# Auto_scatlas_numbat_conservative_recut.R
#
# Analysis registry
# Status: active terminal/audit workflow.
# Script: analysis/cnv/Auto_scatlas_numbat_conservative_recut.R
# Short description: cut Numbat trees to a conservative clone layer so scATLAS
# validation avoids over-fragmented subclone calls.
# Methodology: analysis/methodology/cnv/scatlas_numbat_methodology.md
# Inputs:
# - ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv
# - ref_outs/Auto_scatlas_numbat/by_samples/<sample>/numbat/treeML_<iter>.rds
# - ref_outs/Auto_scatlas_numbat/by_samples/<sample>/numbat/clone_post_<iter>.tsv
# Outputs:
# - ref_outs/Auto_scatlas_numbat/conservative_clones/Auto_scatlas_numbat_conservative_clone_summary.csv
# - ref_outs/Auto_scatlas_numbat/conservative_clones/by_samples/<sample>/Auto_<sample>_numbat_conservative_clone_post.csv
# Cache/replot behavior: set SCATLAS_FORCE_REBUILD=TRUE to refresh outputs.
# Run command: Rscript analysis/cnv/Auto_scatlas_numbat_conservative_recut.R
# Conda/container environment: Numbat container or R environment with numbat.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(numbat)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
out_root <- file.path(root_dir, "ref_outs")
setwd(out_root)

preferred_n_cut <- as.integer(Sys.getenv("SCATLAS_NUMBAT_CONSERVATIVE_N_CUT", "3"))
min_clone_frac <- as.numeric(Sys.getenv("SCATLAS_NUMBAT_CONSERVATIVE_MIN_FRAC", "0.03"))
min_clone_cells_floor <- as.integer(Sys.getenv("SCATLAS_NUMBAT_CONSERVATIVE_MIN_CELLS", "20"))

manifest_path <- "Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv"
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)
manifest <- fread(manifest_path)

out_dir <- "Auto_scatlas_numbat/conservative_clones"
by_sample_dir <- file.path(out_dir, "by_samples")
dir.create(by_sample_dir, recursive = TRUE, showWarnings = FALSE)

final_iter_from <- function(numbat_dir, prefix = "treeML") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.rds")))
  if (length(files) == 0) return(NA_integer_)
  iter <- suppressWarnings(as.integer(sub(paste0("^", prefix, "_([0-9]+)\\.rds$"), "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

####################
terminal_no_subclone_status <- function(sample_id, numbat_dir) {
  done_file <- file.path(numbat_dir, paste0("Auto_", sample_id, "_numbat_done.txt"))
  if (!file.exists(done_file)) return(NULL)
  done_lines <- readLines(done_file, warn = FALSE)
  if (!any(grepl("^terminal_no_subclone\\tTRUE$", done_lines))) return(NULL)
  status_line <- done_lines[grepl("^status\\t", done_lines)]
  n_cells_line <- done_lines[grepl("^n_cells\\t", done_lines)]
  data.frame(
    sample = sample_id,
    status = "terminal_no_subclone",
    numbat_status = sub("^status\\t", "", if (length(status_line) > 0) status_line[1] else ""),
    n_cells = suppressWarnings(as.integer(sub("^n_cells\\t", "", if (length(n_cells_line) > 0) n_cells_line[1] else NA_character_))),
    stringsAsFactors = FALSE
  )
}
####################

as_vertex_df <- function(g) {
  attrs <- vertex_attr(g)
  as.data.frame(lapply(attrs, function(x) {
    if (is.list(x)) vapply(x, function(y) paste(y, collapse = ","), character(1)) else x
  }), stringsAsFactors = FALSE)
}

read_genotype_matrix <- function(path) {
  geno <- fread(path)
  if (!("cell" %in% colnames(geno))) stop("Missing cell column in genotype matrix: ", path)
  cells <- geno$cell
  geno$cell <- NULL
  mat <- as.matrix(geno)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- cells
  mat
}

tree_clone_counts <- function(gtree) {
  vertex_df <- as_vertex_df(gtree)
  leaf <- as.logical(vertex_df$leaf)
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  sort(table(clone[leaf]), decreasing = TRUE)
}

merge_minor_clone_post <- function(clone_post, gtree, min_clone_cells) {
  clone_post <- as.data.frame(clone_post)
  clone_post$clone_opt_raw_recut <- clone_post$clone_opt
  clone_post$GT_opt_raw_recut <- clone_post$GT_opt
  clone_post$p_opt_raw_recut <- clone_post$p_opt

  clone_raw <- as.character(clone_post$clone_opt_raw_recut)
  clone_counts <- sort(table(clone_raw), decreasing = TRUE)
  major <- names(clone_counts)[as.integer(clone_counts) >= min_clone_cells]
  if (length(major) == 0) major <- names(clone_counts)[1]

  p_cols <- paste0("p_", major)
  available_p_cols <- p_cols[p_cols %in% colnames(clone_post)]
  major_assign <- clone_raw
  minor_idx <- which(!clone_raw %in% major)
  if (length(minor_idx) > 0) {
    if (length(available_p_cols) > 0) {
      p_mat <- as.matrix(clone_post[minor_idx, available_p_cols, drop = FALSE])
      storage.mode(p_mat) <- "numeric"
      best_col <- max.col(p_mat, ties.method = "first")
      major_assign[minor_idx] <- sub("^p_", "", available_p_cols[best_col])
    } else {
      major_assign[minor_idx] <- major[1]
    }
  }

  vertex_df <- as_vertex_df(gtree)
  clone_gt <- vertex_df %>%
    filter(!is.na(.data$clone), nzchar(as.character(.data$clone))) %>%
    group_by(clone = as.character(.data$clone)) %>%
    summarise(GT = dplyr::first(as.character(.data$GT)), .groups = "drop")
  clone_gt_map <- setNames(clone_gt$GT, clone_gt$clone)

  clone_post$clone_opt <- major_assign
  clone_post$minor_clone_merged <- clone_raw != major_assign
  clone_post$GT_opt <- unname(clone_gt_map[as.character(clone_post$clone_opt)])
  clone_post$GT_opt[is.na(clone_post$GT_opt)] <- clone_post$GT_opt_raw_recut[is.na(clone_post$GT_opt)]
  clone_post
}

summary_rows <- list()
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample[i]
  numbat_dir <- manifest$numbat_dir[i]
  iter <- final_iter_from(numbat_dir, "treeML")
  message("Conservative re-cut: ", sample_id)

  if (!is.finite(iter)) {
    terminal_status <- terminal_no_subclone_status(sample_id, numbat_dir)
    if (!is.null(terminal_status)) {
      summary_rows[[sample_id]] <- terminal_status
    } else {
      summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_treeML", stringsAsFactors = FALSE)
    }
    next
  }

  tree_file <- file.path(numbat_dir, paste0("treeML_", iter, ".rds"))
  geno_file <- file.path(numbat_dir, paste0("geno_", iter, ".tsv"))
  clone_file <- file.path(numbat_dir, paste0("clone_post_", iter, ".tsv"))
  missing <- c(tree_file, geno_file, clone_file)[!file.exists(c(tree_file, geno_file, clone_file))]
  if (length(missing) > 0) {
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_input", missing = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
    next
  }

  tree_ml <- readRDS(tree_file)
  P <- read_genotype_matrix(geno_file)
  n_cells <- nrow(P)
  min_clone_cells <- max(min_clone_cells_floor, ceiling(n_cells * min_clone_frac))
  g_cut <- get_gtree(tree_ml, P, n_cut = preferred_n_cut, max_cost = 0)
  clone_counts <- tree_clone_counts(g_cut)

  clone_post <- fread(clone_file)
  clone_post_recut <- merge_minor_clone_post(clone_post, g_cut, min_clone_cells)

  sample_out <- file.path(by_sample_dir, sample_id)
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)
  out_clone_file <- file.path(sample_out, paste0("Auto_", sample_id, "_numbat_conservative_clone_post.csv"))
  out_tree_file <- file.path(sample_out, paste0("Auto_", sample_id, "_tree_final_conservative.rds"))
  fwrite(clone_post_recut, out_clone_file)
  saveRDS(g_cut, out_tree_file)

  recut_counts <- sort(table(as.character(clone_post_recut$clone_opt)), decreasing = TRUE)
  summary_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "ok",
    numbat_iteration = iter,
    n_cut = preferred_n_cut,
    n_cells = n_cells,
    min_clone_cells_threshold = min_clone_cells,
    raw_tree_n_clones = length(clone_counts),
    conservative_n_clones = length(recut_counts),
    largest_clone_cells = as.integer(max(recut_counts)),
    smallest_clone_cells = as.integer(min(recut_counts)),
    output_clone_file = normalizePath(out_clone_file, mustWork = TRUE),
    output_tree_file = normalizePath(out_tree_file, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_tbl <- bind_rows(summary_rows)
fwrite(summary_tbl, file.path(out_dir, "Auto_scatlas_numbat_conservative_clone_summary.csv"))
message("Wrote conservative clone summary.")
