library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(ggrepel)
library(patchwork)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

## ------------------------------------------------------------
## Colour palettes
## ------------------------------------------------------------
# For studies / technologies (viridis-based)
default_viridis_colors <- viridis_pal()(10)
shuffled_order <- c(1, 6, 2, 3, 9, 5, 8, 4, 7, 10)
shuffled_colors <- default_viridis_colors[shuffled_order]

before_df <- read.csv("filtering_summary_EAC_Ref_t.csv", header = TRUE, check.names = FALSE)[, -1]
merged_obj <- readRDS("EAC_Ref_merged_loose.rds")
after_df  <- count(merged_obj@meta.data, orig.ident, name = "Final")

study_colors <- c(
  "Alcindor_2025" = "steelblue",
  "Baek_2025"     = "darkorange",
  "Carroll_2023"  = "peru",
  "Croft_2022"    = "forestgreen",
  "Ju_2025"       = "royalblue",
  "Lambroia_2024" = "chocolate",
  "Strasser_2025" = "darkgreen",
  "Walker_2025"   = "dodgerblue",
  "Wu_2018"       = "sienna",
  "Yates_2025"    = "seagreen"
)

before_df <- before_df %>%
  rename(
    orig.ident           = sample,
    Raw                  = raw,
    Mitochondrial_Filter = `mito_DNA\npercentage < `,
    Gene_Count_Filter    = `number of\ngenes`,
    Housekeeping_Filter  = `housekeeping\nexpression > `
  ) %>%
  mutate(
    study = sub("^(.*?_\\d{4}).*", "\\1", orig.ident)
  )

after_df <- merged_obj@meta.data %>%
  count(orig.ident, name = "Singlets")

before_df <- before_df %>%
  left_join(after_df, by = "orig.ident") %>%
  mutate(Singlets = replace_na(Singlets, 0L))

## ------------------------------------------------------------
## 2) GLOBAL FILTERING BAR PLOT (All studies combined)
## ------------------------------------------------------------

global_data <- before_df %>%
  summarise(
    Raw       = sum(Raw,       na.rm = TRUE),
    Mito_orig = sum(Mitochondrial_Filter, na.rm = TRUE),
    Gene_orig = sum(Gene_Count_Filter, na.rm = TRUE),
    HK_orig   = sum(Housekeeping_Filter, na.rm = TRUE),
    Singlets  = sum(Singlets,  na.rm = TRUE)
  ) %>%
  pivot_longer(
    everything(),
    names_to  = "step",
    values_to = "cell_count"
  ) %>%
  mutate(
    step = factor(
      recode(step,
             "Raw"       = "Raw",
             "Mito_orig" = "Mito_DNA < ",
             "Gene_orig" = "Min_genes > ",
             "HK_orig"   = "HK_expr > ",
             "Singlets"  = "Final singlets"
      ),
      levels = c("Raw", "Mito_DNA < ", "Min_genes > ", "HK_expr > ", "Final singlets")
    )
  )

global_bar <- ggplot(global_data, aes(x = step, y = cell_count, fill = step)) +
  geom_col(width = 0.7, color = "gray20", linewidth = 0.25) +
  geom_text(
    aes(label = comma(cell_count)),
    vjust = -0.5,
    size = 3.2
  ) +
  scale_fill_brewer(palette = "Blues", direction = -1, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Total Cells Remaining After Each Filter (all studies)",
    x = NULL,
    y = "Number of Cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("global_filter_bar.png", global_bar, width = 8, height = 4, dpi = 300)

## ------------------------------------------------------------
## 3) PER-SAMPLE BEFORE vs AFTER (log10), coloured by study
## ------------------------------------------------------------

counts <- before_df %>%
  transmute(
    orig.ident,
    before = Raw,
    after  = Singlets,
    study
  )

counts_long <- counts %>%
  pivot_longer(
    cols = c("before", "after"),
    names_to = "filter_status",
    values_to = "cell_count"
  ) %>%
  mutate(
    filter_status = factor(filter_status, levels = c("before", "after")),
    orig.ident    = factor(orig.ident)
  )

# Sort samples by study
counts_long <- counts_long %>%
  mutate(orig.ident = factor(
    orig.ident,
    levels = counts %>% arrange(study) %>% pull(orig.ident)
  ))

axis_label_colors <- ifelse(
  levels(counts_long$orig.ident) %in% dplyr::filter(counts, after == 0)$orig.ident,
  "red",
  "black"
)

per_sample_bar <- ggplot(
  counts_long,
  aes(x = orig.ident, y = cell_count + 1, fill = study, alpha = filter_status)
) +
  geom_col(
    width    = 0.75,
    position = position_dodge(width = 0.85),
    color    = "gray20",
    linewidth = 0.15
  ) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001, 10001, 100001),
    labels = c("0", "10", "100", "1k", "10k", "100k")
  ) +
  scale_fill_manual(values = study_colors, drop = FALSE) +
  scale_alpha_manual(
    name   = "",
    values = c("before" = 1.0, "after" = 0.6),
    labels = c("Raw", "Final")
  ) +
  labs(
    title = "Cells Before vs After Filtering (per sample)",
    x = NULL,
    y = "Number of Cells (log10 scale)",
    fill = "Study"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = axis_label_colors
    ),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, face = "italic"),
    legend.position  = "top",
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  )

ggsave("per_sample_before_after.png", per_sample_bar, width = 12, height = 6, dpi = 300)

## ------------------------------------------------------------
## 4) PER-STUDY FILTERING SUMMARY with study-specific thresholds
## ------------------------------------------------------------

qc_rules <- data.frame(
  pattern = c(
    "Alcindor_2025",
    "Baek_2025",
    "Carroll_2023",
    "Croft_2022",
    "Ju_2025",
    "Lambroia_2024",
    "Strasser_2025",
    "Walker_2025",
    "Wu_2018",
    "Yates_2025"
  ),
  mito   = c(25, 15, 25, 10, 40, 10, 10, 25, 1, 20),
  nGenes = c(300, 300, 300, 500, 300, 500, 500, 300, 5000, 300),
  hk     = c(1, 3, 2, 3, 2, 4, 3, 2, 5, 0.5),
  stringsAsFactors = FALSE
)

study_summary <- before_df %>%
  group_by(study) %>%
  summarise(
    Raw       = sum(Raw,       na.rm = TRUE),
    Mito_orig = sum(Mitochondrial_Filter, na.rm = TRUE),
    Gene_orig = sum(Gene_Count_Filter, na.rm = TRUE),
    HK_orig   = sum(Housekeeping_Filter, na.rm = TRUE),
    Singlets  = sum(Singlets,  na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  left_join(qc_rules, by = c("study" = "pattern"))

study_long <- study_summary %>%
  pivot_longer(
    cols = c(Raw, Mito_orig, Gene_orig, HK_orig, Singlets),
    names_to  = "step",
    values_to = "cell_count"
  ) %>%
  mutate(
    step_order = case_when(
      step == "Raw"       ~ 1L,
      step == "Mito_orig" ~ 2L,
      step == "Gene_orig" ~ 3L,
      step == "HK_orig"   ~ 4L,
      step == "Singlets"  ~ 5L
    ),
    step_label = case_when(
      step == "Raw"       ~ "Raw",
      step == "Mito_orig" ~ paste0("Mito_DNA < ", mito),
      step == "Gene_orig" ~ paste0("Min_genes > ", nGenes),
      step == "HK_orig"   ~ paste0("HK_expr > ", hk),
      step == "Singlets"  ~ "Final singlets",
      TRUE ~ step
    ),
    step_f = factor(
      step_order,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Raw", "Mito", "Gene", "HK", "Singlets")
    )
  )

per_study_bar <- ggplot(
  study_long,
  aes(x = step_f, y = cell_count, fill = step_f)
) +
  geom_col(width = 0.65, color = "gray20", linewidth = 0.25) +
  geom_text(
    aes(label = comma(cell_count), y = cell_count),
    vjust = -0.4,
    size = 3
  ) +
  facet_wrap(~ study, scales = "free_y") +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.12, 0.14))
  ) +
  scale_fill_brewer(
    palette = "Blues",
    direction = -1,
    guide = "none"
  ) +
  labs(
    title = "Cells Remaining After Each Filter (by study, with study-specific thresholds)",
    x = NULL,
    y = "Number of Cells"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    strip.text    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 25, hjust = 1),
    plot.margin   = margin(10, 10, 20, 10),
    panel.spacing = unit(1, "lines")
  ) +
  coord_cartesian(clip = "off")

per_study_bar <- per_study_bar +
  geom_text(
    data = study_long,
    aes(x = step_f, y = 0, label = step_label),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 2.7
  )

ggsave("per_study_filtering.png", per_study_bar, width = 16, height = 10, dpi = 300)


tech_pie_data <- summary_df %>%
  group_by(study) %>%
  summarise(Final = sum(Final), .groups = "drop") %>%
  mutate(
    fraction = Final / sum(Final),
    percent  = percent(fraction, accuracy = 0.1),
    label    = paste0(study, "\n", comma(Final), " (", percent, ")")
  ) %>%
  arrange(desc(study)) %>%
  mutate(ypos = cumsum(fraction) - 0.5 * fraction)

tech_pie <- ggplot(tech_pie_data, aes(x = 1, y = fraction, fill = study)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.95) +  # Extended radial space
  geom_segment(aes(x = 1.02, xend = 1.15, y = ypos, yend = ypos),
               inherit.aes = FALSE, color = "grey50", linewidth = 0.3) +
  geom_label_repel(
    data = tech_pie_data,
    aes(x = 1.55, y = ypos, label = label),  # Pushed further out
    inherit.aes = FALSE,
    direction = "y", nudge_x = 0.1, box.padding = 0.45, point.padding = 0,
    segment.color = "grey50", min.segment.length = 0, size = 3.3,
    max.overlaps = Inf, label.r = unit(0.1, "lines")
  ) +
  scale_fill_manual(values = shuffled_colors) +
  labs(title = "Final Distribution by Technology") +
  theme_void(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(5, 12, 5, 5)
  )

ggsave("study_contribution.png", per_study_bar, width = 16, height = 10, dpi = 300)

## ------------------------------------------------------------
## 5) CELL TYPE COMPOSITION (stacked by study using merged_obj)
## ------------------------------------------------------------

meta_df <- merged_obj@meta.data

celltype_study_df <- meta_df %>%
  count(celltype_update, study, name = "n") %>%
  mutate(
    study_short = str_extract(study, "^([^_]+_[^_]+)")
  )

celltype_df <- celltype_study_df %>%
  group_by(celltype_update) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    prop          = n / sum(n),
    percent_label = paste0("N=", comma(n), "\n(", sprintf("%.1f%%", 100 * prop), ")")
  ) %>%
  arrange(desc(n)) %>%
  mutate(celltype_update = factor(celltype_update, levels = celltype_update))

celltype_study_df <- celltype_study_df %>%
  mutate(
    celltype_update = factor(celltype_update, levels = levels(celltype_df$celltype_update))
  )

study_levels <- sort(unique(celltype_study_df$study_short))
palette_for_plot <- setNames(
  shuffled_colors[seq_along(study_levels)],
  study_levels
)

celltype_bar <- ggplot(celltype_study_df,
                       aes(x = celltype_update, y = n, fill = study_short)) +
  geom_col(color = "gray20", linewidth = 0.25, width = 0.8) +
  geom_text(
    data = celltype_df,
    aes(x = celltype_update, y = n, label = percent_label),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.5,
    lineheight = 1.1
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = palette_for_plot, name = "Study", drop = FALSE) +
  labs(
    title = "Cell Type Counts (stacked by study)",
    x = NULL,
    y = "Number of Cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, size = 10, face = "bold"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  )

ggsave("celltype_stacked_by_study.png", celltype_bar, width = 10, height = 5, dpi = 300)

library(colorspace)
group_palette <- rainbow_hcl(10, c = 120, l = 60)

ct_pie_data <- celltype_df %>%
  mutate(fraction = n / sum(n)) %>%
  arrange(desc(celltype_update)) %>%
  mutate(ypos = cumsum(fraction) - 0.5 * fraction)

celltype_pie <- ggplot(ct_pie_data, aes(x = 1, y = fraction, fill = celltype_update)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.6) +
  geom_segment(aes(x = 1.02, xend = 1.12, y = ypos, yend = ypos),
               inherit.aes = FALSE, color = "grey70", linewidth = 0.25) +
  scale_fill_manual(values = group_palette, name = "") +
  labs(title = "Cell Type Composition") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    plot.margin = margin(5, 12, 5, 5)
  )

ggsave("celltype_composition_pie.png", celltype_pie, width = 8, height = 6, dpi = 300)