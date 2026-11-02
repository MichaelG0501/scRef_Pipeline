####################
# Retained for provenance. Current entry points are
# analysis/clinical/Auto_bayesprism_bulk_deconvolution.R and
# analysis/clinical/Auto_malignant_bulk_survival.R.
# This historical workflow is not part of the current execution order.
# Its reference misclassified state-unassigned malignant cells as normal.
####################
####################
# Analysis registry:
# Status: legacy
# Script: analysis/clinical/legacy_Auto_occams_bayesprism_qc_plots.R
# Description: Replot-only OCCAMS BayesPrism numerical QC dashboard.
# Methodology: analysis/methodology/clinical/bayesprism_malignant_state_deconvolution_methodology.md
# Inputs: ref_outs/OCCAMS/bayesprism/tables/Auto_OCCAMS_bayesprism_{state_fractions,mp_scores,gene_coverage}.csv;
#   ref_outs/OCCAMS/tables/OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv;
#   current MP grouping and state colours from analysis/shared/scRef_config.R.
# Outputs under ref_outs/OCCAMS/bayesprism/:
#   figures/Auto_OCCAMS_bayesprism_qc_overview.{pdf,png}
#   tables/Auto_OCCAMS_bayesprism_qc_{sample_order,mp_heatmap,state_summary}.csv
#   reports/Auto_OCCAMS_bayesprism_visual_qc.txt
#   updates/new_updates/summaries/Auto_OCCAMS_bayesprism_visual_qc.csv
# Cache: complete outputs are skipped unless SCREF_FORCE_REBUILD=TRUE or
#   SCREF_REPLOT_ONLY=TRUE. All calculations use persistent small CSV inputs.
# Run: /opt/pbs/bin/qsub analysis/clinical/legacy_Auto_occams_bayesprism_qc_plots.sh
# Environment: dmtcp. No fitting, survival tests, or external downloads.
####################
library(ggplot2)
library(patchwork)
library(data.table)
library(scales)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
root <- file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "bayesprism")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
report_dir <- file.path(root, "reports")
for (path in c(table_dir, figure_dir, report_dir, SCREF_SUMMARY_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
prefix <- "Auto_OCCAMS_bayesprism"
pdf_path <- file.path(figure_dir, paste0(prefix, "_qc_overview.pdf"))
png_path <- file.path(figure_dir, paste0(prefix, "_qc_overview.png"))
order_path <- file.path(table_dir, paste0(prefix, "_qc_sample_order.csv"))
heatmap_path <- file.path(table_dir, paste0(prefix, "_qc_mp_heatmap.csv"))
state_path <- file.path(table_dir, paste0(prefix, "_qc_state_summary.csv"))
report_path <- file.path(report_dir, paste0(prefix, "_visual_qc.txt"))
summary_path <- file.path(SCREF_SUMMARY_DIR, paste0(prefix, "_visual_qc.csv"))
outputs <- c(pdf_path, png_path, order_path, heatmap_path, state_path, report_path, summary_path)
force_rebuild <- toupper(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE")) == "TRUE"
replot_only <- toupper(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE")) == "TRUE"
if (all(file.exists(outputs)) && !force_rebuild && !replot_only) quit(save = "no")

fractions <- fread(file.path(table_dir, paste0(prefix, "_state_fractions.csv")))
scores <- fread(file.path(table_dir, paste0(prefix, "_mp_scores.csv")))
coverage <- fread(file.path(table_dir, paste0(prefix, "_gene_coverage.csv")))
grouping <- fread(SCREF_MP_GROUPING_CSV)
qc <- fread(file.path(SCREF_REF_OUTS_DIR, "OCCAMS", "tables", "OCCAMS_RNAseq_GRCh37_sample_qc_summary.csv"))
states <- c(SCREF_PRIMARY_STATE_ORDER, SCREF_TECHNICAL_STATE_ORDER)
fraction_columns <- paste0("fraction_tumour__", states)
bulk_columns <- paste0("fraction_bulk__", states)
expected <- qc[Annotation_Compatible_Below_50_Flag == "no", Subject_ID]
stopifnot(nrow(fractions) == 281L, nrow(scores) == 281L, nrow(grouping) == 17L,
          !anyDuplicated(fractions$sample_id), !anyDuplicated(scores$sample_id),
          setequal(fractions$sample_id, expected), setequal(scores$sample_id, expected),
          setequal(names(scores)[-1L], grouping$mp), setequal(coverage$mp, grouping$mp))
mp_order <- grouping$mp
score_matrix <- as.matrix(scores[match(fractions$sample_id, scores$sample_id), ..mp_order])
tumour_matrix <- as.matrix(fractions[, ..fraction_columns])
bulk_matrix <- as.matrix(fractions[, ..bulk_columns])
max_closure_error <- max(abs(rowSums(tumour_matrix) - 1))
stopifnot(all(is.finite(score_matrix)), all(is.finite(tumour_matrix)),
          all(tumour_matrix >= 0 & tumour_matrix <= 1), max_closure_error < 1e-6,
          max(abs(rowSums(bulk_matrix) - fractions$malignant_fraction_initial)) < 1e-6,
          all(is.finite(fractions$malignant_fraction_final)),
          all(fractions$malignant_fraction_final >= 0 & fractions$malignant_fraction_final <= 1))
score_sd <- apply(score_matrix, 2, sd)
stopifnot(all(score_sd > 0))
z_scores <- sweep(sweep(score_matrix, 2, colMeans(score_matrix), "-"), 2, score_sd, "/")
fractions[, sample_rank := frank(malignant_fraction_final, ties.method = "first")]
order_table <- fractions[order(sample_rank), .(sample_id, sample_rank,
  malignant_fraction_initial, malignant_fraction_final)]
state_long <- melt(fractions[, c("sample_id", "sample_rank", fraction_columns), with = FALSE],
  id.vars = c("sample_id", "sample_rank"), variable.name = "state", value.name = "fraction")
state_long[, state := factor(sub("fraction_tumour__", "", state, fixed = TRUE), levels = states)]
state_summary <- state_long[, .(mean_fraction = mean(fraction), median_fraction = median(fraction)), by = state]
heatmap <- as.data.table(z_scores)
heatmap[, sample_id := fractions$sample_id]
heatmap <- melt(heatmap, id.vars = "sample_id", variable.name = "mp", value.name = "z_score")
heatmap[, sample_rank := fractions$sample_rank[match(sample_id, fractions$sample_id)]]
heatmap[, mp := factor(mp, levels = rev(mp_order))]
coverage[, coverage_fraction := genes_present / genes_defined]
coverage[, mp := factor(mp, levels = rev(mp_order))]
fwrite(order_table, order_path)
fwrite(heatmap, heatmap_path)
fwrite(state_summary, state_path)

theme_set(theme_minimal(base_size = 15) + theme(
  panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 12, colour = "#526170"),
  axis.title = element_text(size = 14), plot.title.position = "plot",
  legend.title = element_text(size = 12), legend.text = element_text(size = 11),
  plot.margin = margin(12, 14, 12, 14)))
ink <- "#17324D"
teal <- "#168A8A"
median_final <- median(fractions$malignant_fraction_final)
correlation <- cor(fractions$malignant_fraction_initial, fractions$malignant_fraction_final, method = "spearman")
p_a <- ggplot(fractions, aes(malignant_fraction_initial, malignant_fraction_final)) +
  geom_abline(slope = 1, intercept = 0, colour = "#AAB6C1", linetype = 2, linewidth = 0.7) +
  geom_point(colour = teal, alpha = 0.65, size = 2.2) +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1)) +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
  coord_fixed() + labs(title = "A  Malignant RNA contribution",
    subtitle = sprintf("281 samples | final median %.1f%% | Spearman rho = %.2f", median_final * 100, correlation),
    x = "Initial posterior", y = "Updated posterior")
p_b <- ggplot(coverage, aes(coverage_fraction, mp)) +
  geom_segment(aes(x = 0.9, xend = coverage_fraction, yend = mp), colour = "#D7E4E8", linewidth = 1.1) +
  geom_point(size = 3, colour = teal) +
  geom_text(aes(label = paste0(genes_present, "/", genes_defined)), hjust = -0.28, size = 3.6, colour = ink) +
  scale_x_continuous(labels = label_percent(), limits = c(0.9, 1.05), breaks = c(0.9, 0.95, 1)) +
  labs(title = "B  MP genes captured", subtitle = "Present / defined genes | all 17 MPs retained",
    x = "Coverage in inferred tumour expression", y = NULL) +
  theme(panel.grid.major.y = element_blank())
p_c <- ggplot(state_long, aes(sample_rank, fraction, fill = state)) +
  geom_col(width = 1) + scale_fill_manual(values = SCREF_STATE_COLOURS, drop = FALSE,
    labels = function(x) vapply(x, function(s) paste(strwrap(s, width = 24), collapse = "\n"), character(1))) +
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 70, 140, 210, 281)) +
  labs(title = "C  Composition within the malignant compartment",
    subtitle = "Five biological states + Hybrid / Unresolved | initial state posterior",
    x = "Samples ordered by updated malignant RNA fraction", y = "Within-tumour RNA fraction", fill = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = "bottom", panel.grid = element_blank())
p_d <- ggplot(heatmap, aes(sample_rank, mp, fill = z_score)) +
  geom_raster() + scale_fill_gradient2(low = "#3B6FB6", mid = "#FAFAF7", high = "#C8504B",
    midpoint = 0, limits = c(-2.5, 2.5), oob = squish, name = "Within-MP\nz-score") +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 70, 140, 210, 281)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "D  Tumour-only MP activity", subtitle = "Gaussian GSVA on log2(CPM + 1) | standardized within each MP",
    x = "Same sample order as C", y = NULL) + theme(panel.grid = element_blank())
dashboard <- ((p_a | p_b) / (p_c | p_d)) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(title = "OCCAMS | BayesPrism output quality control",
    subtitle = "281 subjects     43,122 inferred genes     17 metaprograms     5 malignant states",
    caption = paste("Numerical QC passed; biological accuracy is not established.",
      "21,171 epithelial cells without current state labels were retained as nonmalignant in this trial.",
      "State fractions are RNA contributions; they are not validated cell counts. Heatmap colours are clipped at +/-2.5 SD.", sep = "\n"),
    theme = theme(plot.title = element_text(size = 25, face = "bold", colour = ink),
      plot.subtitle = element_text(size = 16, colour = "#526170"),
      plot.caption = element_text(size = 12, hjust = 0, colour = "#526170")))
ggsave(pdf_path, dashboard, width = 20, height = 14, bg = "white")
ggsave(png_path, dashboard, width = 20, height = 14, dpi = 200, bg = "white")
summary <- data.table(subjects = nrow(fractions), mps = ncol(score_matrix),
  median_final_malignant_fraction = median_final, initial_final_spearman = correlation,
  min_mp_gene_coverage = min(coverage$coverage_fraction), max_tumour_fraction_closure_error = max_closure_error,
  mean_stress_adaptive_fraction = state_summary[state == "Stress-adaptive", mean_fraction])
fwrite(summary, summary_path)
writeLines(c("OCCAMS BayesPrism visual QC: numerical checks PASS", capture.output(print(summary)),
  "The initial/final scatter measures update magnitude; it is not ground-truth validation.",
  "Sample order is identical in composition and MP panels; no clustering or significance testing is applied.",
  "Mean Stress-adaptive contribution is reported for follow-up, not treated as evidence of accuracy.",
  "The fitted reference assigns state-less epithelial cells to Non_malignant_epithelial; this can bias tumour estimates.",
  "Source CSVs, sample ordering, state summary and untruncated heatmap z-scores are saved in tables/.",
  paste("Figures:", pdf_path, png_path)), report_path)
stopifnot(all(file.exists(outputs)), all(file.info(outputs)$size > 0))
message("Visual QC complete: ", pdf_path)
