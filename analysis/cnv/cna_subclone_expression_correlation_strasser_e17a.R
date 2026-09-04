####################
# Analysis registry:
#   Status: active
#   Script: analysis/cnv/cna_subclone_expression_correlation_strasser_e17a.R
#   Methodology: analysis/methodology/cnv/cna_subclone_expression_correlation_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Generates a single-sample MP score boxplot by CNA subclone for Strasser_2025_E17A, styled exactly after the Lee et al. 2020 reference image.
#   Inputs:
#     - ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv
#     - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   Outputs:
#     - figures/Auto_Strasser_2025_E17A_subclone_mp_scores.pdf
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(stringr)
})

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_helpers.R")

setwd(SCREF_PROJECT_DIR)
set.seed(42)

out_dir <- file.path(SCREF_REF_OUTS_DIR, "Auto_cna_subclone_expression", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading cell metadata...")
cells_path <- file.path(SCREF_REF_OUTS_DIR, "Auto_malignant_subclone_mp", "Auto_malignant_subclone_cells.csv")
cells <- read.csv(cells_path, check.names = FALSE)
cells_e17a <- cells %>% filter(sample == "Strasser_2025_E17A")

if(nrow(cells_e17a) == 0) stop("No cells found for Strasser_2025_E17A")

message("Loading UCell scores...")
ucell_path <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "centred", "mp_refinement", "intermediate", "merged_refined_ucell_scores.rds")
ucell <- readRDS(ucell_path)

common_cells <- intersect(rownames(ucell), cells_e17a$cell)
cells_e17a <- cells_e17a %>% filter(cell %in% common_cells)

message("Z-scoring UCell scores across the cohort...")
ucell_scaled <- as.data.frame(scale(ucell))

requested_mps <- c("MP1", "MP2+", "MP8+", "MP16", "MP12", "MP15")
missing_mps <- setdiff(requested_mps, colnames(ucell_scaled))
if(length(missing_mps) > 0) stop("Missing MPs in UCell object: ", paste(missing_mps, collapse = ", "))

ucell_sub <- ucell_scaled[cells_e17a$cell, requested_mps, drop = FALSE]
ucell_sub$cell <- rownames(ucell_sub)

df <- cells_e17a %>%
  select(cell, subclone) %>%
  left_join(ucell_sub, by = "cell")

df_long <- df %>%
  pivot_longer(cols = starts_with("MP"), names_to = "MP", values_to = "Score")

message("Loading pre-calculated tests...")
tests_path <- file.path(SCREF_REF_OUTS_DIR, "Auto_malignant_subclone_mp", "Auto_malignant_subclone_mp_subclone_tests.csv")
tests_df <- read.csv(tests_path, check.names = FALSE)
e17a_tests <- tests_df %>% 
  filter(sample == "Strasser_2025_E17A", mp %in% requested_mps) %>%
  distinct(mp, p_adj) %>%
  mutate(
    p.signif = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

mp_descriptions <- SCREF_MP_DESCRIPTIONS

# Replace MP names with descriptions, formatted for plotting
df_long <- df_long %>%
  mutate(
    MP_Label = mp_descriptions[MP],
    MP_Label = factor(MP_Label, levels = mp_descriptions[requested_mps]),
    subclone = factor(subclone)
  )

message("Formatting statistical tests...")
stat.test <- df_long %>%
  group_by(MP_Label) %>%
  wilcox_test(Score ~ subclone) %>%
  add_xy_position(x = "MP_Label", dodge = 0.72) %>%
  left_join(
    e17a_tests %>% 
      mutate(
        MP_Label = mp_descriptions[mp],
        MP_Label = factor(MP_Label, levels = mp_descriptions[requested_mps])
      ) %>%
      select(MP_Label, p.signif_adj = p.signif),
    by = "MP_Label"
  ) %>%
  mutate(p.signif = p.signif_adj)

message("Generating plot...")
p <- ggplot(df_long, aes(x = MP_Label, y = Score, fill = subclone)) +
  geom_boxplot(
    outlier.shape = 16, outlier.size = 1.2, outlier.color = "black",
    linewidth = 0.8, width = 0.62, position = position_dodge(width = 0.72), color = "black"
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 16, size = 1.8, color = "red", 
    position = position_dodge(width = 0.72), show.legend = FALSE
  ) +
  stat_pvalue_manual(
    stat.test, label = "p.signif", tip.length = 0.02, size = 6, 
    bracket.size = 0.6, step.increase = 0, hide.ns = FALSE, vjust = 0.5
  ) +
  scale_fill_manual(
    values = c("Subclone A" = "#1A9850", "Subclone B" = "#984EA3", "Subclone C" = "grey80"),
    labels = c("Subclone A" = "Subclone A", "Subclone B" = "Subclone B", "Subclone C" = "Subclone C")
  ) +
  labs(
    title = "Strasser_2025_E17A",
    x = "MPs",
    y = "Scores"
  ) +
  scale_y_continuous(limits = c(min(df_long$Score) - 0.2, max(stat.test$y.position) + 0.5)) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 20),
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, hjust = 0.5),
    legend.key.size = unit(1.2, "cm"),
    legend.background = element_rect(fill = "transparent")
  ) +
  guides(fill = guide_legend(label.position = "left", label.hjust = 0.5))

pdf_path <- file.path(out_dir, "Auto_Strasser_2025_E17A_subclone_mp_scores.pdf")
ggsave(pdf_path, plot = p, width = 10.5, height = 6.5, useDingbats = FALSE)

message("Saved figure to ", pdf_path)
