####################
# Auto_task7_basal_smg_relationship.R
# Per-sample relationship analysis between:
# - Basal to Intestinal Metaplasia
# - SMG-like Metaplasia
# using both state proportions and max MP-group expression.
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "task7"
out_dir <- paste0(task_prefix, "_basal_smg_relationship")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")
# Note: all_meta.rds (covering all cell types) is missing, 
# so we use the total epithelial population per sample for now.
# total_counts_all <- readRDS("all_meta.rds") 

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)

common_cells <- intersect(Cells(tmdata_all), names(state_B))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))
tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

meta_df <- data.frame(
  cell = common_cells,
  sample = as.character(tmdata_all$orig.ident[common_cells]),
  study = as.character(tmdata_all$study[common_cells]),
  state = as.character(state_B[common_cells]),
  stringsAsFactors = FALSE
)

total_counts_epi <- meta_df %>% 
  count(sample, name = "n_total_sample")

target_states <- c("Basal to Intestinal Metaplasia", "SMG-like Metaplasia")

prop_df <- meta_df %>%
  count(sample, study, state, name = "n_state") %>%
  left_join(total_counts_epi, by = "sample") %>%
  mutate(prop = n_state / n_total_sample) %>%
  select(sample, study, state, prop, n_state, n_total_sample)

# Prepare wide version for scatter plots (Basal and SMG specific)
target_wide <- prop_df %>%
  filter(state %in% target_states) %>%
  pivot_wider(names_from = state, values_from = c(prop, n_state), values_fill = 0)

if (!all(c("prop_Basal to Intestinal Metaplasia", "prop_SMG-like Metaplasia") %in% colnames(target_wide))) {
  stop("Missing required proportion columns for Basal/SMG")
}

cor_prop <- suppressWarnings(cor.test(
  target_wide$`prop_Basal to Intestinal Metaplasia`,
  target_wide$`prop_SMG-like Metaplasia`,
  method = "spearman"
))

p_scatter_prop <- ggplot(target_wide, aes(`prop_Basal to Intestinal Metaplasia`, `prop_SMG-like Metaplasia`, color = study)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Basal to Intestial Metaplasia vs SMG-like Metaplasia (state proportion)",
    subtitle = paste0("Spearman rho = ", round(unname(cor_prop$estimate), 3), ", p = ", signif(cor_prop$p.value, 3)),
    x = "Basal to Intestinal Metaplasia proportion",
    y = "SMG-like Metaplasia proportion",
    color = "Study"
  ) +
  theme_classic(base_size = 12)

expr_df <- data.frame(
  cell = common_cells,
  sample = as.character(tmdata_all$orig.ident[common_cells]),
  study = as.character(tmdata_all$study[common_cells]),
  basal_max = apply(mp_adj_noncc[common_cells, intersect(state_groups[["Basal to Intestinal Metaplasia"]], colnames(mp_adj_noncc)), drop = FALSE], 1, max),
  smg_max = apply(mp_adj_noncc[common_cells, intersect(state_groups[["SMG-like Metaplasia"]], colnames(mp_adj_noncc)), drop = FALSE], 1, max),
  stringsAsFactors = FALSE
)

expr_sample <- expr_df %>%
  group_by(sample, study) %>%
  summarise(
    basal_max_mean = mean(basal_max, na.rm = TRUE),
    smg_max_mean = mean(smg_max, na.rm = TRUE),
    basal_max_median = median(basal_max, na.rm = TRUE),
    smg_max_median = median(smg_max, na.rm = TRUE),
    .groups = "drop"
  )

cor_expr <- suppressWarnings(cor.test(expr_sample$basal_max_mean, expr_sample$smg_max_mean, method = "spearman"))

p_scatter_expr <- ggplot(expr_sample, aes(basal_max_mean, smg_max_mean, color = study)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  labs(
    title = "Basal to Intestial Metaplasia vs SMG-like Metaplasia (max MP expression)",
    subtitle = paste0("Spearman rho = ", round(unname(cor_expr$estimate), 3), ", p = ", signif(cor_expr$p.value, 3)),
    x = "Basal-group max MP expression (sample mean)",
    y = "SMG-group max MP expression (sample mean)",
    color = "Study"
  ) +
  theme_classic(base_size = 12)

# Prepare matrix for full abundance Heatmap
all_present_states <- intersect(
  c("Classic Proliferative", "Basal to Intestinal Metaplasia", "Stress-adaptive", 
    "SMG-like Metaplasia", "Immune Infiltrating", "Unresolved", "Hybrid"),
  unique(prop_df$state)
)

heat_wide <- prop_df %>%
  select(sample, study, state, prop) %>%
  pivot_wider(names_from = state, values_from = prop, values_fill = 0) %>%
  arrange(study, desc(`Basal to Intestinal Metaplasia`))

mat_prop_all <- as.matrix(heat_wide[, all_present_states])
rownames(mat_prop_all) <- heat_wide$sample

# Annotation
anno_df <- data.frame(Study = heat_wide$study)
rownames(anno_df) <- heat_wide$sample

col_study <- setNames(
  DiscretePalette(length(unique(anno_df$Study)), palette = "polychrome"),
  unique(anno_df$Study)
)

ha_row <- rowAnnotation(
  Study = anno_df$Study,
  col = list(Study = col_study),
  show_annotation_name = FALSE
)

# Colors for proportion (normalized to all cells)
# Max proportion might be small, let's check max
max_p <- max(mat_prop_all, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_p/2, max_p), c("white", "orange", "firebrick"))

ht <- Heatmap(
  mat_prop_all,
  name = "Proportion of\nall cells",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  right_annotation = ha_row,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (mat_prop_all[i, j] > 0.001) {
      grid.text(sprintf("%.1f%%", mat_prop_all[i, j] * 100), x, y, gp = gpar(fontsize = 5))
    }
  },
  column_title = "Epithelial state distribution (as % of all epithelial cells)",
  border = TRUE
)

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_basal_smg_relationship_plots.pdf")), width = 10, height = 12, onefile = TRUE)
print(p_scatter_prop)
print(p_scatter_expr)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

write.csv(target_wide, file.path(out_dir, paste0("Auto_", task_prefix, "_sample_proportions_basal_smg.csv")), row.names = FALSE)
write.csv(expr_sample, file.path(out_dir, paste0("Auto_", task_prefix, "_sample_expression_basal_smg.csv")), row.names = FALSE)

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  task = task_prefix,
  n_samples = nrow(target_wide),
  spearman_prop_rho = unname(cor_prop$estimate),
  spearman_prop_p = cor_prop$p.value,
  spearman_expr_rho = unname(cor_expr$estimate),
  spearman_expr_p = cor_expr$p.value,
  stringsAsFactors = FALSE
), file.path(summary_dir, paste0("Auto_", task_prefix, "_basal_smg_relationship_summary.csv")), row.names = FALSE)

cat("Saved task7 outputs in:", out_dir, "\n")
