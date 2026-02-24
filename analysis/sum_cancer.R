sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]  # match *_*_*
library(dplyr)
library(stringr)

epi_list <- list()
meta_all <- data.frame()
for (sample in sample_dirs) {
  sample_dir <- file.path("by_samples", sample)
  file_path <- file.path(sample_dir, paste0(sample, "_epi_f.rds"))
  if (!file.exists(file_path)) {
    print(paste("Epithelial does not exist for:", sample))
    next
  }
  epi <- readRDS(file_path)
  epi_list[[sample]] <- epi
  meta_epi <- epi@meta.data
  meta_all <- bind_rows(meta_all, meta_epi)
  print(paste("Loaded for:", sample))
}

#meta_all$malignancy <- ifelse(meta_all$malignancy %in% c("malignant_level_1", "malignant_level_2"), meta_all$malignancy, "unresolved")

save <- meta_all

bad_origidents <- meta_all %>%
  group_by(orig.ident) %>%
  tally() %>%
  filter(n < 30)

bad_origidents <- bad_origidents$orig.ident

meta_all <- meta_all %>%
  mutate(
    study = str_replace(orig.ident, "^([^_]+_[^_]+)_.*", "\\1"),
    malignancy_simple = case_when(
      grepl("^malignant", malignancy) ~ "malignant",
      TRUE ~ "non_malignant"
    )
  ) %>%
  filter(!orig.ident %in% bad_origidents)


plot_df <- meta_all %>%
  group_by(orig.ident) %>%
  summarise(
    malignant_prop = mean(malignancy_simple == "malignant")
  )

ggplot(plot_df, aes(x = malignant_prop, y = reorder(orig.ident, malignant_prop))) +
  geom_col() +
  theme_bw() +
  labs(title = "Malignant proportion per orig.ident",
       x = "Proportion malignant", y = "")


study_df <- meta_all %>%
  group_by(study, malignancy) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 0) %>%
  group_by(study) %>%
  mutate(prop = n / sum(n))

ggplot(study_df, aes(x = study, y = prop, fill = malignancy)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Malignancy composition per study",
       y = "Proportion", x = "Study", fill = "Class")

ggplot(study_df, aes(x = study, y = prop, fill = malignancy)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Malignancy composition per study",
    y = "Proportion",
    x = "Study",
    fill = "Class"
  ) +
  scale_fill_manual(
    values = c(
      malignant_level_1 = "#E41A1C",  # red
      malignant_level_2 = "#FB6A6A",  # lighter red (similar hue)
      unresolved        = "grey30"    # dark grey
    )
  )





library(dplyr)
library(stringr)
library(ggplot2)

bad_origidents <- meta_all %>%
  group_by(orig.ident) %>%
  tally() %>%
  filter(n < 30)

bad_origidents <- bad_origidents$orig.ident

meta_all <- meta_all %>%
  mutate(
    study = str_replace(orig.ident, "^([^_]+_[^_]+)_.*", "\\1"),
    malignancy_simple = case_when(
      grepl("^malignant", malignancy) ~ "malignant",
      TRUE ~ "non_malignant"
    )
  ) %>%
  filter(!orig.ident %in% bad_origidents)

summary_tbl <- meta_all %>%
  group_by(study, orig.ident) %>%
  summarise(
    n_total = n(),
    n_malignant = sum(malignancy_simple == "malignant", na.rm = TRUE),
    pct_malignant = 100 * n_malignant / n_total,
    .groups = "drop"
  )


ggplot(summary_tbl, aes(x = orig.ident, y = pct_malignant, fill = study)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = round(pct_malignant * n_total / 100, 1)),
            vjust = -0.2, size = 3) +
  facet_wrap(~ study, scales = "free_x") +
  labs(
    x = NULL,
    y = "% malignant cells",
    title = "Fraction of malignant cells per sample by study"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",    # color equals facet; hide redundant legend
    plot.title = element_text(face = "bold"), 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )






library(readxl)
library(dplyr)
library(purrr)
library(stringr)

# Load marker data
markers <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx", sheet = 1)
markers <- markers[markers$specificity > 0.2 & markers$cell_type == "Malignant", ]
combine_marker_scores <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }
  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)
  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}
markers_ranked <- combine_marker_scores(markers, w_specificity = 0.2, w_sensitivity = 0.8)


DoHeatmap(
  object = epi,
  features = markers_ranked$gene[1:200],
  slot = "data",
  group.by = "malignant_clus",
  size = 3
) + scale_fill_viridis_c()


DoHeatmap(
  object = epi,
  features = cs_genes,
  slot = "data",
  group.by = "malignant_clus",
  size = 3
) + scale_fill_viridis_c()

