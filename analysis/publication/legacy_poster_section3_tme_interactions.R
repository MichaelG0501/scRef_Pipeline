####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/publication/legacy_poster_section3_tme_interactions.R
#   Methodology: not required (legacy poster assembly)
#   Inputs:
#     ref_outs/non_malignant_mp_correlations/01_cancer_mps_cross_only/Auto_celltype_interaction_dotmap_data.csv
#     ref_outs/non_malignant_mp_correlations/01_cancer_mps_cross_only/Auto_celltype_positive_edges_lr_annotated.csv
#   Outputs:
#     ref_outs/publication/section3/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section3_tme_interactions.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(readxl)
  library(forcats)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section3"
out_dir <- pub_section_dir(section)
workbook_file <- file.path(SCREF_REF_OUTS_DIR, "non_malignant_mp_correlations", "02_cancer_mps_cross_and_within", "Auto_cancer_tme_interactions_TME_centric.xlsx")

cat("Generating TME-cancer MP dotmap from TME-centric workbook...\n")

if (!file.exists(workbook_file)) {
  make_placeholder(section, "s3_tme_dotmap", "TME-cancer MP dotmap",
                   "Missing Auto_cancer_tme_interactions_TME_centric.xlsx.\nRun mp_cross_celltype_correlations.R first.")
  quit(save = "no")
}

wide <- read_excel(workbook_file, sheet = "TME-Cancer Interactions", col_names = TRUE)

# Remove "within" part (i.e. cancer vs cancer columns)
cols_to_keep <- names(wide)[!grepl("^cancer ", names(wide))]

# Reverse map labels back to MPs using the description map
label_to_mp <- setNames(names(SCREF_MP_DESCRIPTIONS), SCREF_MP_DESCRIPTIONS)

dot <- lapply(cols_to_keep, function(nm) {
  vals <- as.character(wide[[nm]])
  vals <- na.omit(vals)
  
  stat_idx <- grep("^-log10\\(p\\):", vals)
  if (length(stat_idx) == 0) return(NULL)
  
  labels <- vals[stat_idx - 1]
  stats <- vals[stat_idx]
  
  tibble(
    tme_celltype = sub(" MP.*$", "", nm),
    tme_mp = nm,
    cancer_label = labels,
    cancer_mp = unname(label_to_mp[labels]),
    neglog10_p = suppressWarnings(as.numeric(sub(".*-log10\\(p\\): ([0-9.]+).*", "\\1", stats))),
    support = suppressWarnings(as.numeric(sub(".*Support: ([0-9.]+)%.*", "\\1", stats)))
  )
}) |> bind_rows() |>
  filter(!is.na(cancer_mp))

mp_group_df <- tibble(mp = PUB_MP_ORDER, group = pub_mp_state(PUB_MP_ORDER)) |>
  mutate(group = factor(group, levels = PUB_MP_STATE_ORDER)) |>
  arrange(group, match(mp, PUB_MP_ORDER)) |>
  mutate(group_id = as.integer(group),
         x = row_number() + cumsum(c(0, diff(group_id) != 0)) * 0.7)

dot <- dot |>
  filter(cancer_mp %in% mp_group_df$mp) |>
  group_by(tme_celltype, cancer_mp) |>
  slice_max(neglog10_p, n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(mp_group_df |> select(cancer_mp = mp, x, group), by = "cancer_mp")

tme_levels <- c("fibroblast", "macrophage", "endothelial", "cd8", "cd4", "nk", "plasma")
dot <- dot |>
  mutate(tme_celltype = factor(tme_celltype, levels = rev(tme_levels[tme_levels %in% unique(tme_celltype)])),
         y = as.numeric(tme_celltype))

ann <- mp_group_df |> filter(mp %in% unique(dot$cancer_mp)) |> mutate(y = length(levels(dot$tme_celltype)) + 0.65)
dot_plot <- ggplot(dot, aes(x, y)) +
  geom_tile(data = ann, aes(x, y), inherit.aes = FALSE,
            width = 0.92, height = 0.18,
            fill = PUB_MP_STATE_COLOURS[as.character(ann$group)], colour = NA) +
  geom_point(aes(size = support, fill = neglog10_p),
             shape = 21, colour = "#111111", stroke = 0.28, alpha = 0.96) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = expression(-log[10](p))) +
  scale_size_continuous(range = c(2.2, 7.5), name = "Sample\nSupport") +
  scale_x_continuous(breaks = mp_group_df$x[mp_group_df$mp %in% unique(dot$cancer_mp)],
                     labels = short_mp_label(mp_group_df$mp[mp_group_df$mp %in% unique(dot$cancer_mp)]),
                     expand = expansion(mult = c(0.02, 0.04)),
                     position = "top") +
  scale_y_continuous(breaks = seq_along(levels(dot$tme_celltype)),
                     labels = levels(dot$tme_celltype),
                     expand = expansion(mult = c(0.02, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "TME compartment") +
  pub_theme(12) +
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "right",
        plot.margin = margin(8, 10, 8, 8))
save_pub_gg(dot_plot, section, "s3_tme_dotmap", width = 9.3, height = 4.7)
write_csv(dot, file.path(out_dir, "tables", "s3_tme_dotmap_data.csv"))

cat("Section 3 complete.\n")
