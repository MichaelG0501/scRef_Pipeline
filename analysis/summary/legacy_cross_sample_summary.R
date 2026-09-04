####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/summary/legacy_cross_sample_summary.R
#   Methodology: analysis/methodology/summary/summary_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Historical interactive QC-summary script with external ephemeral
#     input and caller-directory outputs; superseded by current pipeline summaries.
#   Inputs: historical filtering_summary_EAC_Ref_t.csv and external ephemeral EAC_Ref_filtered.rds.
#   Outputs: historical caller-directory PNG figures; no current dependency.
#   Run: provenance only; not for submission.
#   Conda env: dmtcp
####################

####################
# Moved from: analysis/summary.R
# Reorganized as part of analysis/ restructuring
####################
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)   # For combining plots
library(ggrepel)     # For non-overlapping labels
library(stringr)

library(scales) 
default_viridis_colors <- viridis_pal()(10)
shuffled_order <- c(1, 6, 2, 3, 9, 5, 8, 4, 7, 10)
shuffled_colors <- default_viridis_colors[shuffled_order]

# -----------------------------
# Step 1: Load and Clean Data
# -----------------------------
before_df <- read.csv("filtering_summary_EAC_Ref_t.csv", header = TRUE, check.names = FALSE)[, -1]
merged_obj <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/EAC_Ref_filtered.rds")
after_df  <- count(merged_obj@meta.data, orig.ident, name = "Final")

summary_df <- full_join(before_df, after_df, by = c("sample" = "orig.ident")) %>%
  rename(
    orig.ident           = sample,
    Raw                  = raw,
    Mitochondrial_Filter = `mito_DNA\npercentage < `,
    Gene_Count_Filter    = `number of\ngenes`,
    Housekeeping_Filter  = `housekeeping\nexpression > `
  ) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  select(orig.ident, Raw, Mitochondrial_Filter, Gene_Count_Filter, Housekeeping_Filter, Final) %>%
  mutate(study = str_extract(orig.ident, "^([^_]+_[^_]+)"))

# ------------------------------------------
# Pair 1: Filtering Bar + Technology Pie
# ------------------------------------------

bar_data <- summary_df %>%
  summarise(across(where(is.numeric), sum)) %>%
  pivot_longer(cols = everything(), names_to = "filter_step", values_to = "cell_count") %>%
  mutate(filter_step = recode(filter_step,
                              "Raw"                  = "Raw",
                              "Mitochondrial_Filter" = "Mito_DNA < 15",
                              "Gene_Count_Filter"    = "Min_genes > 500",
                              "Housekeeping_Filter"  = "HK_expr > 3",
                              "Final"                = "Final (Singlets)"
  )) %>%
  mutate(filter_step = factor(filter_step, levels = c(
    "Raw", "Mito_DNA < 15", "Min_genes > 500", "HK_expr > 3", "Final (Singlets)"
  )))
bar_plot <- ggplot(bar_data, aes(x = filter_step, y = cell_count, fill = filter_step)) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.25) +
  geom_text_repel(
    aes(label = comma(cell_count), y = cell_count),
    nudge_y = 0.04 * max(bar_data$cell_count, na.rm = TRUE),
    direction = "y", box.padding = 0.25, point.padding = 0.25,
    segment.color = NA, min.segment.length = 0, size = 3.2, max.overlaps = Inf
  ) +
  scale_fill_brewer(palette = "Blues", direction = -1, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)), labels = comma) +
  labs(title = "Total Cells Remaining After Each Filter", x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  )

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

pair1_combo <- (bar_plot + tech_pie) + plot_layout(ncol = 2, widths = c(1.6, 1))
ggsave("pair1_combo.png", pair1_combo, width = 12, height = 6, dpi = 300)

# ------------------------------------------
# Pair 2: Cell Type Bar + Cell Type Pie
# ------------------------------------------

# celltype_df <- lapply(tmdata_annotated, function(obj) {
#   data.frame(celltype_update = obj$celltype_update, stringsAsFactors = FALSE)
# }) %>%
#   bind_rows() %>%
#   count(celltype_update, name = "n") %>%
#   mutate(
#     celltype_update = factor(celltype_update, levels = sort(unique(celltype_update))),
#     percent_label = paste0("N=", comma(n), "\n(", sprintf("%.1f%%", 100 * n / sum(n)), ")")
#   )

library(dplyr)
library(scales)

celltype_df <- as.data.frame(table(EAC_Ref_merged$celltype_update), stringsAsFactors = FALSE) %>%
  rename(celltype_update = Var1, n = Freq) %>%
  filter(celltype_update != "unresolved_inconsistent") %>%
  arrange(desc(n)) %>%   # 🔹 sort by counts
  mutate(
    celltype_update = factor(celltype_update, levels = celltype_update),
    percent_label = paste0(
      "N=", comma(n),
      "\n(", sprintf("%.1f%%", 100 * n / sum(n)), ")"
    )
  )
celltype_df <- celltype_df[celltype_df$celltype_update != "unresolved_inconsistent", ]

library(colorspace)
group_palette <- rainbow_hcl(9, c = 120, l = 60)

# celltype_study_df <- lapply(tmdata_annotated, function(obj) {
#   data.frame(
#     celltype_update = obj$celltype_update,
#     study = sapply(strsplit(obj$orig.ident, "_"), function(x) paste(x[1:2], collapse = "_")),
#     stringsAsFactors = FALSE
#   )
# }) %>%
#   bind_rows() %>%
#   count(celltype_update, study, name = "n") %>%
#   mutate(
#     celltype_update = factor(celltype_update, levels = levels(celltype_df$celltype_update))
#   )

celltype_study_df <- as.data.frame(table(EAC_Ref_merged$celltype_update, EAC_Ref_merged$study), stringsAsFactors = FALSE) %>%
  rename(celltype_update = Var1, study = Var2, n = Freq) %>%
  # Ensure the factor order matches the total_df so the x-axis is correct
  mutate(
    celltype_update = factor(celltype_update, levels = levels(celltype_df$celltype_update))
  )
celltype_study_df <- celltype_study_df[!is.na(celltype_study_df$celltype_update), ]

# --- MODIFIED BAR PLOT ---

library(ggplot2)
library(dplyr)
library(scales)
library(stringr)

# --- OPTIONAL: wrap long celltype labels nicely (tweak width as needed) ---
celltype_study_df2 <- celltype_study_df %>%
  mutate(celltype_label = str_wrap(as.character(celltype_update), width = 18))

celltype_df2 <- celltype_df %>%
  mutate(celltype_label = str_wrap(as.character(celltype_update), width = 18))

# Ensure ordering is preserved after making wrapped labels
celltype_study_df2 <- celltype_study_df2 %>%
  mutate(celltype_label = factor(celltype_label, levels = celltype_df2$celltype_label))

celltype_df2 <- celltype_df2 %>%
  mutate(celltype_label = factor(celltype_label, levels = celltype_df2$celltype_label))

# --- PLOT ---
celltype_bar <- ggplot(celltype_study_df2, aes(x = celltype_label, y = n, fill = study)) +
  # stacked bars with subtle border
  geom_col(width = 0.78, linewidth = 0.35) +
  
  # add an outer outline to the full stack (gives a crisp “publication” look)
  geom_col(
    data = celltype_df2,
    aes(x = celltype_label, y = n),
    inherit.aes = FALSE,
    width = 0.78,
    fill = NA,
    color = "gray15",
    linewidth = 0.6
  ) +
  
  # total labels above each bar
  geom_text(
    data = celltype_df2,
    aes(x = celltype_label, y = n, label = percent_label),
    inherit.aes = FALSE,
    vjust = -0.25,
    size = 4.6,           # publish-friendly (adjust for your output size)
    fontface = "bold",
    lineheight = 1.0,
    color = "gray10"
  ) +
  
  # color palette (yours)
  scale_fill_manual(values = shuffled_colors, name = "Study") +
  
  # axes
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.18))  # extra headroom for labels
  ) +
  
  labs(
    title = "Cell Type Composition by Study",
    x = NULL,
    y = "Number of Cells"
  ) +
  
  # classic + publication tweaks
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    
    # subtle y gridlines for readability (classic theme has none)
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # legend formatting
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width  = unit(0.45, "cm"),
    
    # margins so labels never clip in saved figures
    plot.margin = margin(10, 18, 10, 10)
  )

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

pair2_combo <- (celltype_bar + celltype_pie) + plot_layout(ncol = 2, widths = c(1.5, 1))
pair2_combo


########################################################################

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Read the file and drop the first unnamed column
before_df <- read.table("filtering_summary_EAC_Ref.csv", 
                        sep = ",", 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE)

# Drop the first column (unnamed index)
before_df <- before_df[, -1]

# Rename columns safely
before_df <- before_df %>%
  rename(orig.ident = `sample`, before = `raw`)

# Step 2: Count cells after filtering from Seurat object
after_df <- count(merged_obj@meta.data, orig.ident, name = "after")

# Step 3: Merge before and after counts
counts <- full_join(before_df, after_df, by = "orig.ident") %>%
  replace_na(list(before = 0, after = 0)) %>%
  mutate(study = str_extract(orig.ident, "^([^_]+_[^_]+)"))

# Step 6: Reshape for plotting
counts_long <- counts %>%
  pivot_longer(
    cols = c("before", "after"),
    names_to = "filter_status",
    values_to = "cell_count"
  ) %>%
  mutate(
    filter_status = factor(filter_status, levels = c("before", "after")),
    orig.ident = factor(orig.ident)
  )

# Step 7: Sort samples by technology
counts_long <- counts_long %>%
  mutate(orig.ident = factor(orig.ident, levels = counts %>%
                               arrange(study) %>%
                               pull(orig.ident)))

study_colors <- c(
  "Alcindor_2025" = "steelblue",      # matches snSeq
  "Baek_2025"     = "darkorange",     # matches GEMX
  "Croft_2022"    = "forestgreen",    # matches Multiome
  "Ju_2025"       = "royalblue",      # variant of steelblue
  "Lambroia_2024" = "chocolate",      # variant of darkorange
  "Strasser_2025" = "darkgreen",      # variant of forestgreen
  "Walker_2025"   = "dodgerblue",     # another blue tone
  "Wu_2018"       = "sienna",         # earthy orange-brown
  "Yates_2025"    = "seagreen"        # green variant
)

axis_label_colors <- ifelse(
  levels(counts_long$orig.ident) %in% filter(counts, after == 0)$orig.ident,
  "red",
  "black"
)

# Step 9: Plot
ggplot(
  counts_long,
  aes(x = orig.ident, y = cell_count + 1, fill = study, alpha = filter_status)
) +
  geom_col(width=0.75, position = position_dodge(width = 0.85), color = "gray20", linewidth = 0.15) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001, 10001, 100001),
    labels = c("0", "10", "100", "1k", "10k", "100k")
  ) +
  scale_fill_manual(values = shuffled_colors) +
  scale_alpha_manual(
    name = "",
    values = c("before" = 1.0, "after" = 0.6),
    labels = c("Raw", "Final")
  ) +
  labs(
    x = NULL,
    y = "Number of Cells (log10 scale)",
    fill = "Technology"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = axis_label_colors
    ),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  )





##############################################################

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

# --------------------------------
# 0) Load data
# --------------------------------
before_df <- read.csv(
  "filtering_summary_EAC_Ref_t.csv",
  header = TRUE,
  check.names = FALSE
)[, -1]

cell_counts <- sapply(tmdata_annotated, ncol)
cell_counts_df <- data.frame(
  sample = names(cell_counts),
  final_singlets = as.integer(cell_counts),
  stringsAsFactors = FALSE
)

before_df <- merge(
  before_df,
  cell_counts_df,
  by = "sample",
  all.x = TRUE  # keep all rows from before_df
)

# qc rules you showed
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

# --------------------------------
# 1) Add study prefix
# --------------------------------
# sample: "Wu_2018_EA01-Cancer" -> study: "Wu_2018"
before_df <- before_df %>%
  mutate(study = sub("^(.*?_\\d{4}).*", "\\1", sample)) %>%
  rename(
    Raw       = raw,
    Mito_orig = `mito_DNA\npercentage < 20`,
    Gene_orig = `number of\ngenes > 500`,
    HK_orig   = `housekeeping\nexpression > 3`
  )

# =========================================================
# A) GLOBAL PLOT (all samples together)
# =========================================================

global_data <- before_df %>%
  summarise(
    Raw       = sum(Raw, na.rm = TRUE),
    Mito_orig = sum(Mito_orig, na.rm = TRUE),
    Gene_orig = sum(Gene_orig, na.rm = TRUE),
    HK_orig   = sum(HK_orig, na.rm = TRUE), 
    Singlets   = sum(final_singlets, na.rm = TRUE), 
  ) %>%
  pivot_longer(everything(),
               names_to = "step",
               values_to = "cell_count") %>%
  mutate(
    step = recode(step,
                  "Raw"       = "Raw",
                  "Mito_orig" = "Mito_DNA < 20",
                  "Gene_orig" = "Min_genes > 500",
                  "HK_orig"   = "HK_expr > 3", 
                  "Singlets"  = "Singlets"
    ),
    step = factor(step, levels = c("Raw", "Mito_DNA < 20", "Min_genes > 500", "HK_expr > 3", "Singlets"))
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

print(global_bar)

# =========================================================
# B) PER-STUDY PLOT (dynamic labels from qc_rules)
# =========================================================

# 1. aggregate per study
study_summary <- before_df %>%
  group_by(study) %>%
  summarise(
    Raw       = sum(Raw, na.rm = TRUE),
    Mito_orig = sum(Mito_orig, na.rm = TRUE),
    Gene_orig = sum(Gene_orig, na.rm = TRUE),
    HK_orig   = sum(HK_orig, na.rm = TRUE), 
    Singlets   = sum(final_singlets, na.rm = TRUE), 
    .groups   = "drop"
  ) %>%
  # bring in thresholds
  left_join(qc_rules, by = c("study" = "pattern"))

# 2. long format
study_long <- study_summary %>%
  pivot_longer(
    cols = c(Raw, Mito_orig, Gene_orig, HK_orig, Singlets), 
    names_to = "step",
    values_to = "cell_count"
  ) %>%
  mutate(
    # order for plotting
    step_order = case_when(
      step == "Raw"       ~ 1L,
      step == "Mito_orig" ~ 2L,
      step == "Gene_orig" ~ 3L,
      step == "HK_orig"   ~ 4L, 
      step == "Singlets"  ~ 5L, 
    ),
    # per-study label using that study's thresholds
    step_label = case_when(
      step == "Raw"       ~ "Raw",
      step == "Mito_orig" ~ paste0("Mito_DNA < ", mito),
      step == "Gene_orig" ~ paste0("Min_genes > ", nGenes),
      step == "HK_orig"   ~ paste0("HK_expr > ", hk), 
      step == "Singlets"   ~ paste0("Final singlets"), 
      TRUE ~ step
    )
  )

# we need step_label to be a factor *within study*, but ggplot wants global levels.
# So we create a stable order: Raw, Mito (any), Gene (any), HK (any)
# The actual printed text will come from aes(label=...) anyway.
study_long <- study_long %>%
  mutate(
    step_f = factor(
      step_order,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Raw", "Mito", "Gene", "HK", "Singlets")
    )
  )


library(ggplot2)
library(scales)

# Base plot
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
  facet_wrap(
    ~ study,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.12, 0.14))  # Add space below and above
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
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(10, 10, 20, 10),  # Extra bottom margin
    panel.spacing = unit(1, "lines")
  ) +
  coord_cartesian(clip = "off")  # Allow text outside panel

# Add x-axis labels showing actual thresholds
per_study_bar <- per_study_bar +
  geom_text(
    data = study_long,
    aes(x = step_f, y = 0, label = step_label),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 2.7
  )

# Print the plot
print(per_study_bar)





library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# --- 0) Derive the short study from the sample name ---
# e.g., "Alcindor_2025_SRR27335925" -> "Alcindor_2025"
count_mat2 <- count_mat %>%
  mutate(
    study_short = sub("^([^_]+)_([^_]+).*", "\\1_\\2", study)
  )

# --- 1) Recode combined cell-type columns into single bins ---
# nk.cell|t.cell -> nk.cell ; t.cell|nk.cell -> t.cell
# (silently skip if a column is missing)
add_if_present <- function(df, to, from) {
  if (from %in% names(df)) {
    df[[to]] <- df[[to]] + df[[from]]
    df[[from]] <- NULL
  }
  df
}

count_mat2 <- count_mat2 %>%
  mutate(
    nk.cell = ifelse(is.na(nk.cell), 0L, nk.cell),
    t.cell  = ifelse(is.na(t.cell),  0L, t.cell)
  ) %>%
  add_if_present("nk.cell", "nk.cell|t.cell") %>%
  add_if_present("t.cell",  "t.cell|nk.cell")

# --- 2) Long format and aggregate (colSums per study x celltype) ---
long_df <- count_mat2 %>%
  pivot_longer(
    cols = -c(study, study_short),
    names_to = "celltype",
    values_to = "n"
  )

celltype_study_df <- long_df %>%
  group_by(celltype, study = study_short) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop")

# --- 3) Totals across studies for labels and consistent celltype order ---
celltype_df <- celltype_study_df %>%
  group_by(celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    prop = n / sum(n),
    percent_label = paste0("N=", comma(n), "\n(", sprintf("%.1f%%", 100 * prop), ")")
  ) %>%
  arrange(desc(n)) %>%
  mutate(celltype = factor(celltype, levels = celltype))

# apply same order to stacked data
celltype_study_df <- celltype_study_df %>%
  mutate(celltype = factor(celltype, levels = levels(celltype_df$celltype)))

library(stringr)

# --- build a named palette that matches the pie's study levels ---
study_levels <- unique(count_mat2$study_short)             # short study, same as pie
palette_for_plot <- setNames(shuffled_colors[seq_along(study_levels)], study_levels)

# --- add a short study to the stacked-bar data (keep default x order) ---
celltype_study_df <- as.data.frame(table(merged_obj$celltype_manual, merged_obj$study),
                                   stringsAsFactors = FALSE) %>%
  dplyr::rename(celltype_manual = Var1, study = Var2, n = Freq) %>%
  dplyr::mutate(
    celltype_manual = factor(celltype_manual, levels = levels(celltype_df$celltype_manual)),
    study_short = str_extract(study, "^([^_]+_[^_]+)"),             # e.g. Alcindor_2025
    study_short = factor(study_short, levels = study_levels)        # lock to pie order
  )

celltype_bar <- ggplot(celltype_study_df, aes(x = celltype_manual, y = n, fill = study_short)) + 
  geom_col(color = "gray20", linewidth = 0.25, width = 0.8) +
  geom_text(
    data = celltype_df,
    aes(x = celltype_manual, y = n, label = percent_label),
    inherit.aes = FALSE,
    vjust = -0.5, size = 6.2, lineheight = 1.1
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = palette_for_plot, name = "Study", drop = FALSE) +
  labs(title = "Cell Types Counts", y = "Counts") + 
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, size = 16, face = "bold"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  )
