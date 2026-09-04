####################
# Analysis registry:
#   Status: active
#   Script: analysis/shared/scRef_helpers.R
#   Methodology: not required (utility functions only)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# scRef_helpers.R
#
# Shared helper functions for downstream scRef analysis scripts.
#
# Purpose:
#   Provide small reusable functions for output-tier directories, cache-aware
#   rebuilds, run summaries, MP filtering, state ordering, and PPT-readable
#   plotting defaults.
#
# Input:
#   Optional source: analysis/shared/scRef_config.R
#
# Output:
#   None unless a workflow calls the logging/output helper functions.
####################

if (!exists("SCREF_PROJECT_DIR")) {
  source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(as.character(x[1]))) {
    return(y)
  }
  x
}

ensure_output_tiers <- function(base_dir, tiers = SCREF_OUTPUT_TIERS) {
  paths <- file.path(base_dir, tiers)
  for (path in paths) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  stats::setNames(paths, tiers)
}

output_tier_path <- function(base_dir, tier, filename = NULL) {
  if (!tier %in% SCREF_OUTPUT_TIERS) {
    stop("Unknown output tier: ", tier, ". Expected one of: ", paste(SCREF_OUTPUT_TIERS, collapse = ", "))
  }
  dir.create(file.path(base_dir, tier), recursive = TRUE, showWarnings = FALSE)
  if (is.null(filename)) {
    return(file.path(base_dir, tier))
  }
  file.path(base_dir, tier, filename)
}

reuse_cache_enabled <- function(force_env = "SCREF_FORCE_REBUILD") {
  value <- Sys.getenv(force_env, unset = "FALSE")
  !tolower(value) %in% c("true", "1", "yes", "y")
}

load_or_build_cache <- function(cache_path, build_fun, label = basename(cache_path), force = !reuse_cache_enabled()) {
  if (file.exists(cache_path) && !force) {
    message("Reusing cache: ", cache_path)
    return(list(value = readRDS(cache_path), reused = TRUE, label = label, path = cache_path))
  }
  message("Building cache: ", cache_path)
  value <- build_fun()
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, cache_path)
  list(value = value, reused = FALSE, label = label, path = cache_path)
}

filter_silhouette_mps <- function(mp_genes, sil_scores, threshold = SCREF_STATE_THRESHOLDS$mp_silhouette_min) {
  bad_mps <- which(sil_scores < threshold)
  if (length(bad_mps) == 0) {
    return(mp_genes)
  }
  bad_mp_names <- paste0("MP", bad_mps)
  mp_genes[!names(mp_genes) %in% bad_mp_names]
}

state_levels_current <- function(states, include_technical = TRUE) {
  states <- unique(as.character(states))
  levels <- c(
    SCREF_PRIMARY_STATE_ORDER[SCREF_PRIMARY_STATE_ORDER %in% states],
    SCREF_FINAL_EXTRA_STATE_ORDER[SCREF_FINAL_EXTRA_STATE_ORDER %in% states]
  )
  other <- setdiff(states, c(levels, SCREF_TECHNICAL_STATE_ORDER))
  levels <- c(levels, sort(other))
  if (include_technical) {
    levels <- c(levels, SCREF_TECHNICAL_STATE_ORDER[SCREF_TECHNICAL_STATE_ORDER %in% states])
  }
  levels
}

state_colours_current <- function(states) {
  levels <- state_levels_current(states)
  missing <- setdiff(levels, names(SCREF_STATE_COLOURS))
  colours <- SCREF_STATE_COLOURS
  if (length(missing) > 0) {
    if (requireNamespace("scales", quietly = TRUE)) {
      colours <- c(colours, stats::setNames(scales::hue_pal()(length(missing)), missing))
    } else {
      colours <- c(colours, stats::setNames(grDevices::rainbow(length(missing)), missing))
    }
  }
  colours[levels]
}

label_mps <- function(mp_vec, descriptions = SCREF_MP_DESCRIPTIONS) {
  mp_vec <- as.character(mp_vec)
  desc <- descriptions[mp_vec]
  ifelse(is.na(desc), mp_vec, paste0(mp_vec, " - ", desc))
}

z_normalise_by_sample_study <- function(mat, sample_var, study_var) {
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("tibble", quietly = TRUE)) {
    stop("z_normalise_by_sample_study requires dplyr and tibble.")
  }
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df |>
    dplyr::group_by(.data$.study) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(colnames(mat)), ~ stats::sd(.x, na.rm = TRUE)), .groups = "drop") |>
    tibble::column_to_rownames(".study") |>
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df |>
    dplyr::group_by(.data$.sample) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) |>
    dplyr::ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

ppt_theme <- function(base_size = SCREF_DEFAULT_PLOT$base_size) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ppt_theme requires ggplot2.")
  }
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = SCREF_DEFAULT_PLOT$axis_text_size, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text = ggplot2::element_text(size = SCREF_DEFAULT_PLOT$legend_text_size),
      legend.title = ggplot2::element_text(size = SCREF_DEFAULT_PLOT$legend_title_size),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size)
    )
}

save_ppt_pdf <- function(plot, path, width = SCREF_DEFAULT_PLOT$slide_width, height = SCREF_DEFAULT_PLOT$slide_height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

start_run_summary <- function(script, input_files = character(), output_files = character(), parameters = list()) {
  list(
    script = script,
    start_time = as.character(Sys.time()),
    end_time = NA_character_,
    input_files = input_files,
    output_files = output_files,
    parameters = parameters,
    cache_reused = list(),
    session = NULL
  )
}

add_cache_status <- function(run_summary, label, path, reused) {
  run_summary$cache_reused[[label]] <- list(path = path, reused = reused)
  run_summary
}

write_run_summary <- function(run_summary, path, include_session = TRUE) {
  run_summary$end_time <- as.character(Sys.time())
  if (include_session) {
    run_summary$session <- utils::capture.output(utils::sessionInfo())
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(run_summary, path)
  text_path <- sub("\\.rds$", ".txt", path)
  lines <- c(
    paste0("script: ", run_summary$script),
    paste0("start_time: ", run_summary$start_time),
    paste0("end_time: ", run_summary$end_time),
    "",
    "input_files:",
    paste0("  - ", run_summary$input_files),
    "",
    "output_files:",
    paste0("  - ", run_summary$output_files),
    "",
    "parameters:",
    paste0("  - ", names(run_summary$parameters), ": ", unlist(run_summary$parameters, use.names = FALSE))
  )
  writeLines(lines, text_path)
  invisible(path)
}
