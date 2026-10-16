####################
# Analysis registry:
#   Status: active
#   Script: analysis/TCGA/Auto_tcga_workflow_helpers.R
#   Description: shared storage, execution, cache and figure-export guards for TCGA workflows.
#   Methodology: not required (deterministic I/O and validation helpers only).
#   Inputs: explicit paths supplied by the registered TCGA callers.
#   Outputs: caller-selected live intermediate RDS files, figures and source snapshots.
#   Cache/replot: SCREF_FORCE_REBUILD rebuilds caches; SCREF_REPLOT_ONLY renders saved figure grobs only.
#   Run: sourced by registered analysis/TCGA/*.R callers through their Auto_*.sh PBS wrappers; conda: dmtcp.
####################

####################
# Register ggplot/gtable drawing methods before reading cached grobs in standalone callers.
library(ggplot2)
####################

tcga_flag <- function(name) {
  tolower(Sys.getenv(name, unset = "FALSE")) %in% c("true", "1", "yes", "y")
}

tcga_execution_guard <- function() {
  if (!nzchar(Sys.getenv("PBS_JOBID"))) {
    stop("TCGA workflows require a PBS allocation. Submit the registered PBS wrapper.")
  }
  if (tcga_flag("SCREF_FORCE_REBUILD") && tcga_flag("SCREF_REPLOT_ONLY")) {
    stop("SCREF_FORCE_REBUILD and SCREF_REPLOT_ONLY cannot both be TRUE.")
  }
}

tcga_require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "))
}

tcga_live_path <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  if (!startsWith(paste0(parent, "/"), paste0(normalizePath(SCREF_PROJECT_DIR), "/"))) {
    stop("Persistent TCGA paths must resolve inside the live project: ", path)
  }
  if (file.exists(path) && !startsWith(normalizePath(path), paste0(normalizePath(SCREF_PROJECT_DIR), "/"))) {
    stop("Persistent TCGA file resolves outside the live project: ", path)
  }
  path
}

tcga_stage_input <- function(source_path, live_path) {
  live_path <- tcga_live_path(live_path)
  if (!file.exists(live_path)) {
    tcga_require_files(source_path)
    if (!file.copy(source_path, live_path, overwrite = FALSE)) stop("Input import failed: ", source_path)
  }
  live_path
}

tcga_fingerprint <- function(paths) {
  tcga_require_files(paths)
  tools::md5sum(paths)
}

tcga_cache <- function(path, inputs, build_fun, label = basename(path)) {
  path <- tcga_live_path(path)
  fingerprint <- tcga_fingerprint(inputs)
  if (file.exists(path) && !tcga_flag("SCREF_FORCE_REBUILD")) {
    value <- readRDS(path)
    if (!identical(attr(value, "tcga_input_fingerprint"), fingerprint)) {
      stop("Stale or unversioned cache: ", path, ". Set SCREF_FORCE_REBUILD=TRUE.")
    }
    return(list(value = value, reused = TRUE, label = label, path = path))
  }
  if (tcga_flag("SCREF_REPLOT_ONLY")) stop("Replot requires an existing valid cache: ", path)
  value <- build_fun()
  attr(value, "tcga_input_fingerprint") <- fingerprint
  saveRDS(value, path)
  list(value = value, reused = FALSE, label = label, path = path)
}

tcga_validate_gene_sets <- function(mp_genes) {
  expected <- names(SCREF_MP_DESCRIPTIONS)
  if (anyDuplicated(names(mp_genes)) || !setequal(names(mp_genes), expected) || length(mp_genes) != 17L) {
    stop("Expected exactly the current centred refined 17 MPs; excluded MPs must not enter TCGA scoring.")
  }
  grouping <- read.csv(SCREF_MP_GROUPING_CSV, check.names = FALSE)
  expected_states <- c(setNames(rep("Cell cycle", length(SCREF_CC_MPS)), SCREF_CC_MPS),
                       setNames(rep(names(SCREF_STATE_GROUPS), lengths(SCREF_STATE_GROUPS)),
                                unlist(SCREF_STATE_GROUPS, use.names = FALSE)))
  if (!all(c("mp", "state") %in% names(grouping)) || anyDuplicated(grouping$mp) ||
      !setequal(grouping$mp, expected) ||
      !identical(unname(expected_states[grouping$mp]), as.character(grouping$state))) {
    stop("Current MP grouping table and configured five-state contract disagree.")
  }
  mp_genes[expected]
}

tcga_validate_scores <- function(obj, mp_sets, state_sets, samples) {
  if (!identical(obj$mp_sets, mp_sets) || !identical(obj$state_sets, state_sets) ||
      !identical(obj$transform, "log2(TPM + 1)") ||
      !setequal(colnames(obj$mp_scores), names(mp_sets)) ||
      !setequal(colnames(obj$state_scores), names(state_sets)) ||
      !all(samples %in% rownames(obj$mp_scores)) || !all(samples %in% rownames(obj$state_scores)) ||
      any(!is.finite(obj$mp_scores)) || any(!is.finite(obj$state_scores))) {
    stop("GSVA cache does not match current MP/state sets or samples. Rebuild its producer with SCREF_FORCE_REBUILD=TRUE.")
  }
}

tcga_figure_bundle <- list()

tcga_save_figure <- function(filename, plots, width, height) {
  path <- tcga_live_path(file.path(tiers[["figures"]], filename))
  grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  grobs <- NULL
  tryCatch({
    ####################
    # Measure text on the final Cairo device before freezing the plot layout.
    grobs <- lapply(plots, function(plot) {
      if (is.function(plot)) plot <- plot()
      if (inherits(plot, "ggplot")) ggplot2::ggplotGrob(plot) else plot
    })
    ####################
    ####################
    # Layout measurement can initialize an otherwise empty first device page.
    # Draw into that page; advance only between actual figure pages.
    for (page_index in seq_along(grobs)) {
      grob <- grobs[[page_index]]
      if (page_index > 1L) grid::grid.newpage()
      ####################
      ####################
      if (inherits(grob, "Heatmap") || inherits(grob, "HeatmapList")) {
        ComplexHeatmap::draw(grob, newpage = FALSE, ht_gap = grid::unit(0.8, "cm"),
                             heatmap_legend_side = "right", annotation_legend_side = "right",
                             padding = grid::unit(rep(12, 4), "mm"))
      } else {
        grid::grid.draw(grob)
      }
      ####################
    }
  }, finally = grDevices::dev.off())
  tcga_figure_bundle[[filename]] <<- list(grobs = grobs, width = width, height = height)
  invisible(path)
}

tcga_replot <- function(bundle_path, summary_path, log_path) {
  tcga_require_files(c(bundle_path, summary_path))
  bundle <- readRDS(bundle_path)
  if (!is.list(bundle$figures) || !length(bundle$figures)) stop("Invalid figure bundle: ", bundle_path)
  for (filename in names(bundle$figures)) {
    figure <- bundle$figures[[filename]]
    tcga_save_figure(filename, figure$grobs, figure$width, figure$height)
  }
  if (!is.null(bundle$png)) {
    ggplot2::ggsave(tcga_live_path(file.path(tiers[["figures"]], bundle$png$filename)),
                   bundle$png$plot, width = bundle$png$width, height = bundle$png$height, dpi = 300)
  }
  summary <- read.csv(summary_path, check.names = FALSE)
  summary_copy <- tcga_live_path(file.path(SCREF_SUMMARY_DIR, basename(summary_path)))
  write.csv(summary, summary_copy, row.names = FALSE)
  if (file.info(summary_copy)$size >= 100000) stop("Compact TCGA summary must be <100 KB.")
  run <- start_run_summary(script_path, c(bundle_path, summary_path),
                           file.path(tiers[["figures"]], names(bundle$figures)),
                           list(replot_only = TRUE, analysis_recomputed = FALSE))
  write_run_summary(run, log_path)
  message("Replotted saved TCGA figures without downloads, GSVA, tests or clustering.")
}
