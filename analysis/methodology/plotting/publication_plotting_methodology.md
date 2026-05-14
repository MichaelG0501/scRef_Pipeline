# Publication Plotting Methodology

This document covers general plotting scripts under `analysis/plotting/`.

Plotting scripts are terminal figure-generation scripts. They should not create state, MP, or clinical objects that downstream analysis consumes.

Presentation-readability rules:

1. Plots must be legible on PowerPoint slides.
2. Use sufficiently large axis text, legends, labels, and point sizes relative to figure dimensions.
3. Prefer `cairo_pdf()` or `ggsave(..., device = cairo_pdf)` for PDF output when fonts matter.
4. Do not shrink row names or legends to unreadable sizes to fit a figure; increase figure size or split pages.
5. For heatmaps, tune cell size, row label size, column label size, and legend size explicitly.

For common defaults, use or copy from:

- `analysis/shared/scRef_config.R`
- `analysis/shared/scRef_helpers.R`
