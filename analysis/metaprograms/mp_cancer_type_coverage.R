#!/usr/bin/env Rscript
####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/mp_cancer_type_coverage.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Summarize cancer-type breadth of the external 3CA metaprogram
#     catalogue and plot general/shared/specific coverage categories.
#   Inputs: ref_outs/43018_2025_957_MOESM3_ESM.csv.
#   Outputs: ref_outs/Auto_mp_cancer_type_coverage_heatmap_v3.{pdf,png} and
#     ref_outs/Auto_mp_cancer_type_coverage_summary_v3.csv.
#   Cache/replot: deterministic from the source CSV; no hidden cache.
#   Run: Rscript analysis/metaprograms/mp_cancer_type_coverage.R (PBS required)
#   Conda env: dmtcp
####################


library(tidyverse)
library(patchwork)

setwd("ref_outs")

input_file <- "43018_2025_957_MOESM3_ESM.csv"

cat(strrep("=", 60), "\n")
cat("Metaprogram Cancer-Type Coverage (Gavish-Style Layout)\n")
cat(strrep("=", 60), "\n")

####################
# Load and prepare presence matrix
####################
raw_df <- read.csv(input_file)

coverage_matrix <- raw_df %>%
  select(cancer_type, meta_program) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = cancer_type,
    values_from = present,
    values_fill = 0L
  ) %>%
  column_to_rownames("meta_program") %>%
  as.matrix()

mp_breadth <- rowSums(coverage_matrix)
cancer_type_coverage <- colSums(coverage_matrix)

####################
# Apply categories based on breadth
# general: >12, shared: 6-12, specific: <=5
####################
mp_category <- case_when(
  mp_breadth > 12 ~ "General",
  mp_breadth > 5 ~ "Shared",
  TRUE ~ "Specific"
)
names(mp_category) <- names(mp_breadth)

####################
# Order columns and rows
####################
# Order cancer types by coverage (descending)
cancer_types_ordered <- names(sort(cancer_type_coverage, decreasing = TRUE))

# Order MPs by category and then by breadth
mp_order_df <- tibble(
  mp = names(mp_breadth),
  breadth = as.numeric(mp_breadth),
  category = as.character(mp_category)
) %>%
  mutate(
    category_rank = case_when(
      category == "Specific" ~ 3L,  # Specific at the top
      category == "Shared" ~ 2L,    # Shared in middle
      TRUE ~ 1L                     # General at the bottom (reversed from before to match Gavish panel b where General is at bottom)
    )
  ) %>%
  arrange(desc(category_rank), breadth, mp)

mps_ordered <- mp_order_df$mp
coverage_ordered <- coverage_matrix[mps_ordered, cancer_types_ordered, drop = FALSE]

# Define separator positions (dashed lines)
n_specific <- sum(mp_category == "Specific")
n_shared <- sum(mp_category == "Shared")
n_general <- sum(mp_category == "General")

sep1 <- n_specific + 0.5
sep2 <- n_specific + n_shared + 0.5

# Greyscale colors to match Gavish panel b
category_colors <- c(
  "Specific" = "#D9D9D9",  # Light grey
  "Shared" = "#737373",    # Medium grey
  "General" = "#000000"    # Black
)

####################
# Build plotting data frames
####################
heat_df <- as.data.frame(coverage_ordered) %>%
  rownames_to_column("MP") %>%
  pivot_longer(cols = -MP, names_to = "Cancer_Type", values_to = "Present") %>%
  mutate(
    MP = factor(MP, levels = rev(mps_ordered)),
    Cancer_Type = factor(Cancer_Type, levels = cancer_types_ordered),
    Present = factor(Present)
  )

right_bar_df <- tibble(
  MP = factor(mps_ordered, levels = rev(mps_ordered)),
  Count = as.numeric(mp_breadth[mps_ordered]),
  Category = factor(mp_category[mps_ordered], levels = c("Specific", "Shared", "General"))
)

####################
# Publication-style panel components
####################
# Use a single professional color for presence (dark purple/blue)
presence_color <- "#440154" 

heatmap_panel <- ggplot(heat_df, aes(x = Cancer_Type, y = MP, fill = Present)) +
  geom_tile(color = NA) + # NO GRID
  scale_fill_manual(values = c("0" = "white", "1" = presence_color), guide = "none") +
  # Adding horizontal dashed separators between categories
  geom_hline(yintercept = c(n_general + 0.5, n_general + n_shared + 0.5), 
             linetype = "dashed", color = "black", linewidth = 0.4) +
  labs(x = "Cancer type", y = "MPs") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(5, 5, 5, 5)
  )

right_bar <- ggplot(right_bar_df, aes(x = Count, y = MP, fill = Category)) +
  geom_col(width = 0.8, color = NA) +
  scale_fill_manual(values = category_colors, guide = guide_legend(ncol = 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 25, 5)) +
  # Dashed separators matching the heatmap
  geom_hline(yintercept = c(n_general + 0.5, n_general + n_shared + 0.5), 
             linetype = "dashed", color = "black", linewidth = 0.4) +
  labs(x = "Number of cancer types", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 0)
  )

####################
# Compose layout
####################
# Move legend settings to individual plots and combine.
heatmap_panel <- heatmap_panel + theme(legend.position = "bottom")
right_bar <- right_bar + theme(legend.position = "bottom")

# Remove panel labels and collect guides to the bottom to avoid legend clipping.
final_plot <- (heatmap_panel + right_bar) +
  plot_layout(widths = c(4, 1), guides = "collect")

# Save outputs
pdf_file <- "Auto_mp_cancer_type_coverage_heatmap_v3.pdf"
ggsave(pdf_file, final_plot, width = 12, height = 8.5, units = "in", dpi = 320, bg = "white")

png_file <- "Auto_mp_cancer_type_coverage_heatmap_v3.png"
ggsave(png_file, final_plot, width = 12, height = 8.5, units = "in", dpi = 320, bg = "white")

# Summary data
summary_df <- tibble(
  Metaprogram = names(mp_breadth),
  Cancer_Type_Count = as.numeric(mp_breadth),
  Category = as.character(mp_category)
) %>%
  arrange(Category, desc(Cancer_Type_Count), Metaprogram)

write.csv(summary_df, "Auto_mp_cancer_type_coverage_summary_v3.csv", row.names = FALSE)

cat("Figure written:", pdf_file, "\n")
cat("Figure written:", png_file, "\n")
cat("Summary written: Auto_mp_cancer_type_coverage_summary_v3.csv\n")
cat("Done.\n")
