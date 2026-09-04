####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/publication/legacy_publication_helpers.R
#   Methodology: analysis/methodology/README.md
#   Inputs: analysis/shared/scRef_config.R, analysis/shared/scRef_helpers.R
#   Outputs: helper functions only
#   Run command: source from publication plotting scripts
#   Conda environment: dmtcp
####################

source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/shared/scRef_helpers.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(grid)
})

# ── Output directories ──
PUB_BASE_DIR <- file.path(SCREF_REF_OUTS_DIR, "publication")
PUB_ASSET_DIR <- file.path(PUB_BASE_DIR, "assets")
POSTER_ASSET_DIR <- file.path(SCREF_REF_OUTS_DIR, "Auto_conference_poster_plan", "assets", "publication")
dir.create(PUB_BASE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PUB_ASSET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(POSTER_ASSET_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Poster state order: exclude 3CA_EMT_and_Protein_maturation from all poster figures ──
PUB_STATE_ORDER <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating"
)
PUB_STATE_COLOURS <- SCREF_STATE_COLOURS[PUB_STATE_ORDER]

# ── MP ordering: CC first, then MYC, then Basal…, SMG, Stress, Immune ──
PUB_MP_ORDER <- c("MP1", "MP7", "MP9",    # Cell-cycle / DNA repair
                  "MP2",                    # Classic Proliferative (MYC)
                  "MP17", "MP14", "MP5", "MP10", "MP8",  # Basal to Intestinal Metaplasia
                  "MP18", "MP16",           # SMG-like Metaplasia
                  "MP13", "MP12",           # Stress-adaptive
                  "MP15")                   # Immune Infiltrating

# ── MP → state group mapping ──
PUB_MP_TO_STATE <- unlist(lapply(names(SCREF_STATE_GROUPS), function(state) {
  stats::setNames(rep(state, length(SCREF_STATE_GROUPS[[state]])), SCREF_STATE_GROUPS[[state]])
}))
PUB_MP_TO_STATE[c("MP1", "MP9", "MP7")] <- "Cell-cycle / DNA repair"

PUB_MP_STATE_ORDER <- c("Cell-cycle / DNA repair", PUB_STATE_ORDER)
PUB_MP_STATE_COLOURS <- c(
  "Cell-cycle / DNA repair" = "#7A7A7A",
  PUB_STATE_COLOURS
)

# ── A0 poster-optimised ggplot theme ──
# At A0 portrait (84.1 × 118.9 cm), individual panels are typically ~35–40cm wide
# displayed at ~50%, so base_size 14–15 pt ensures ≥10–12pt at final display
pub_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "#111827", family = ""),
      axis.text = element_text(size = rel(0.85), colour = "#111827"),
      axis.title = element_text(size = rel(1.0), face = "bold"),
      legend.title = element_text(size = rel(0.88), face = "bold"),
      legend.text = element_text(size = rel(0.82)),
      strip.text = element_text(size = rel(0.95), face = "bold"),
      plot.title = element_text(size = rel(1.15), face = "bold", colour = "#111827",
                                margin = margin(b = 6)),
      plot.subtitle = element_text(size = rel(0.9), colour = "#374151"),
      plot.caption = element_text(size = rel(0.72), colour = "#4B5563"),
      panel.border = element_rect(colour = "#CBD5E1", fill = NA, linewidth = 0.4),
      axis.line = element_line(colour = "#111827", linewidth = 0.4),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ── Section output directory helper ──
pub_section_dir <- function(section) {
  out <- file.path(PUB_BASE_DIR, section)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
  out
}

# ── Save publication-quality ggplot figure ──
save_pub_gg <- function(plot, section, name, width, height, dpi = 450, copy_asset = TRUE) {
  out_dir <- pub_section_dir(section)
  pdf_path <- file.path(out_dir, "figures", paste0(name, ".pdf"))
  png_path <- file.path(out_dir, "figures", paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  if (copy_asset) {
    file.copy(png_path, file.path(PUB_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(PUB_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
    # Also copy to poster asset folder
    file.copy(png_path, file.path(POSTER_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(POSTER_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
  }
  invisible(c(pdf = pdf_path, png = png_path))
}

# ── Save ComplexHeatmap/grid-based figures ──
save_pub_grid <- function(draw_expr, section, name, width, height, dpi = 450, copy_asset = TRUE) {
  out_dir <- pub_section_dir(section)
  pdf_path <- file.path(out_dir, "figures", paste0(name, ".pdf"))
  png_path <- file.path(out_dir, "figures", paste0(name, ".png"))

  cairo_pdf(pdf_path, width = width, height = height)
  eval(draw_expr)
  dev.off()

  png(png_path, width = width, height = height, units = "in", res = dpi, bg = "white")
  eval(draw_expr)
  dev.off()

  if (copy_asset) {
    file.copy(png_path, file.path(PUB_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(PUB_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
    file.copy(png_path, file.path(POSTER_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(POSTER_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
  }
  invisible(c(pdf = pdf_path, png = png_path))
}

# ── Status logging ──
write_pub_status <- function(section, panel, status, detail = "") {
  out_dir <- pub_section_dir(section)
  path <- file.path(out_dir, "logs", paste0(panel, "_status.csv"))
  readr::write_csv(
    tibble::tibble(panel = panel, status = status, detail = detail, timestamp = as.character(Sys.time())),
    path
  )
  invisible(path)
}

# ── MP label helpers ──
mp_label <- function(mp) {
  desc <- SCREF_MP_DESCRIPTIONS[mp]
  ifelse(is.na(desc), mp, paste0(mp, "\n", desc))
}

short_mp_label <- function(mp) {
  desc <- SCREF_MP_DESCRIPTIONS[mp]
  ifelse(is.na(desc), mp, paste0(mp, " ", desc))
}

####################
# Poster-specific plotting helpers.
####################
pub_blank_title <- function(plot) {
  plot + labs(title = NULL, subtitle = NULL, caption = NULL)
}

pub_clean_term <- function(x) {
  x <- gsub("\\.\\.[^.]+$", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("\\+", "+", x, fixed = TRUE)
  x
}

pub_mp_state <- function(mp) {
  out <- PUB_MP_TO_STATE[as.character(mp)]
  out[is.na(out)] <- "Other"
  out
}

pub_copy_to_assets <- function(src, name = basename(src), section = NULL) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, file.path(PUB_ASSET_DIR, name), overwrite = TRUE)
  file.copy(src, file.path(POSTER_ASSET_DIR, name), overwrite = TRUE)
  if (!is.null(section)) {
    out_dir <- pub_section_dir(section)
    file.copy(src, file.path(out_dir, "figures", name), overwrite = TRUE)
  }
  TRUE
}

pub_save_current_plot <- function(section, name, width, height, dpi = 450, copy_asset = TRUE) {
  out_dir <- pub_section_dir(section)
  pdf_path <- file.path(out_dir, "figures", paste0(name, ".pdf"))
  png_path <- file.path(out_dir, "figures", paste0(name, ".png"))
  dev.copy(cairo_pdf, filename = pdf_path, width = width, height = height)
  dev.off()
  dev.copy(png, filename = png_path, width = width, height = height, units = "in", res = dpi, bg = "white")
  dev.off()
  if (copy_asset) {
    file.copy(png_path, file.path(PUB_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(PUB_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
    file.copy(png_path, file.path(POSTER_ASSET_DIR, paste0(name, ".png")), overwrite = TRUE)
    file.copy(pdf_path, file.path(POSTER_ASSET_DIR, paste0(name, ".pdf")), overwrite = TRUE)
  }
  invisible(c(pdf = pdf_path, png = png_path))
}

# ── State name cleaner (old → new) ──
clean_state <- function(x) {
  dplyr::case_when(
    x == "Basal to Intest. Meta" ~ "Basal to Intestinal Metaplasia",
    x == "Basal_to_Intest_Meta" ~ "Basal to Intestinal Metaplasia",
    x == "Barretts Metaplasia" ~ "Basal to Intestinal Metaplasia",
    x == "SMG_like_Metaplasia" ~ "SMG-like Metaplasia",
    x == "Stress_adaptive" ~ "Stress-adaptive",
    x == "Immune Infiltrated" ~ "Immune Infiltrating",
    TRUE ~ as.character(x)
  )
}

# ── Missing required figure input ──
abort_missing_figure <- function(section, name, title, message, width = 8, height = 5) {
  stop(paste(title, message), call. = FALSE)
}

# ── Significance labels ──
pub_sig_label <- function(p_adj) {
  case_when(
    is.na(p_adj) ~ "",
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}
