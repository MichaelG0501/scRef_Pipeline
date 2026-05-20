####################
# Auto_cna_subclone_expression_visuals_v2.R
#
# Script status: active, terminal.
# Description: presentation-focused v2 visualisation and statistics for
# CNA-subclone expression/state associations. This reads the completed v1
# results and regenerates corrected recurrent-event tests plus slide-readable
# plots.
# Methodology: analysis/methodology/cnv/Auto_cna_subclone_expression_visuals_v2_methodology.md
#
# Inputs:
#   - ref_outs/Auto_cna_subclone_expression/rds/Auto_cna_subclone_expression_results.rds
#
# Outputs:
#   - tables/: recurrent CNA event tests, per-sample deltas, largest-subclone
#     tests, cohort-level event feature values, secondary per-sample deltas,
#     largest-subclone tests, pairwise CNA-distance tests, run summary.
#   - figures/: consensus CNA heatmap, recurrent-event association dot plots,
#     recurrent-event boxplots, chr8q/MYC CNA-group plot, largest-subclone and pairwise
#     distance summaries.
#   - rds/: Auto_v2_visualisation_results.rds.
#   - logs/: none; this is a short replot/statistics workflow run
#     interactively per AGENTS.md.
#   - updates/new_updates/summaries/: compact CSV summary for agent-readable
#     update generation.
# Cache/replot behavior: always reads the v1 RDS cache and regenerates all v2
#   tables and figures; no heavy upstream CNV or expression objects are rebuilt.
# Run command: Rscript analysis/cnv/Auto_cna_subclone_expression_visuals_v2.R
# Conda environment: dmtcp.
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(scales)
  library(grid)
  library(gridExtra)
})

options(stringsAsFactors = FALSE)
set.seed(42)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

event_call_threshold <- 0.05
legacy_event_call_threshold <- 0.10
recurrent_min_sample_fraction <- 0.15
recurrent_min_samples <- 3L

in_rds <- "Auto_cna_subclone_expression/rds/Auto_cna_subclone_expression_results.rds"
if (!file.exists(in_rds)) stop("Missing v1 result RDS: ", in_rds)

out_dir <- "Auto_cna_subclone_expression_v2"
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
rds_dir <- file.path(out_dir, "rds")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(in_rds)
features <- as.data.frame(res$features)
arm_long <- as.data.frame(res$arm_long)
arm_matrix <- as.matrix(res$arm_matrix)
arm_call_matrix <- as.matrix(res$arm_call_matrix)
dominance_df <- as.data.frame(res$dominance)
pairwise_df <- as.data.frame(res$pairwise)
event_summary <- as.data.frame(res$event_summary)
cna_cluster <- res$cna_cluster

cell_path <- "Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv"
if (!file.exists(cell_path)) stop("Missing malignant subclone cell table: ", cell_path)
cell_df <- read.csv(cell_path, check.names = FALSE)

gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)
chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- read.table(
  gene_order_path,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "end"),
  stringsAsFactors = FALSE
) %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start)

max_plot_cells <- 1200L
cna_colour_limit <- 0.15

subclone_palette <- c(
  "Subclone A" = "#D73027",
  "Subclone B" = "#4575B4",
  "Subclone C" = "#1A9850",
  "Subclone D" = "#984EA3",
  "Subclone E" = "#FF7F00",
  "Subclone F" = "#A65628"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

mp_cols <- c(
  "MP1_G2M Cell Cycle" = "#B0B0B0",
  "MP7_DNA Damage Repair" = "#999999",
  "MP9_G1S Cell Cycle" = "#C0C0C0",
  "MP2_MYC-related Proliferation" = "#E41A1C",
  "MP17_Basal-like Transition" = "#4DAF4A",
  "MP14_Hypoxia Adapted Epi." = "#8DA0CB",
  "MP5_Epithelial IFN Resp." = "#66C2A5",
  "MP10_Columnar Diff." = "#A6D854",
  "MP8_Intestinal Diff." = "#FC8D62",
  "MP18_Secretory Diff. (Intest.)" = "#FF7F00",
  "MP16_Secretory Diff. (Gastric)" = "#FFD92F",
  "MP13_Hypoxic Inflam. Epi." = "#984EA3",
  "MP12_Neuro-responsive Epi" = "#E78AC3",
  "MP15_Immune Infiltration" = "#377EB8"
)

write_table <- function(x, filename) {
  write.csv(x, file.path(table_dir, filename), row.names = FALSE)
}

save_rds <- function(x, filename) {
  saveRDS(x, file.path(rds_dir, filename))
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(!is.na(x)) < 2) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_weighted_mean <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

p_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_feature_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values)]
  missing <- setdiff(values, names(cols))
  if (length(missing) > 0) {
    extra <- setNames(colorRampPalette(brewer.pal(8, palette))(length(missing)), missing)
    cols <- c(cols, extra)
  }
  cols[values]
}

assert_named_palette <- function(cols, label) {
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  if (length(cols) == 0) stop("No colours available for ", label)
  if (any(is.na(cols)) || any(!nzchar(cols))) stop("Invalid colours for ", label)
  cols
}

subclone_colours <- function(values) {
  assert_named_palette(complete_palette(subclone_palette, values, "Set2"), "Subclone")
}

mp_descriptions <- c(
  "MP1" = "G2M Cell Cycle",
  "MP7" = "DNA Damage Repair",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8" = "Intestinal Diff.",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP12" = "Neuro-responsive Epi",
  "MP15" = "Immune Infiltration"
)

mp_order <- c("MP1", "MP7", "MP9", "MP2", "MP17", "MP14", "MP5", "MP10",
              "MP8", "MP18", "MP16", "MP13", "MP12", "MP15")
mp_features <- paste0("mp__", mp_order)
mp_features <- mp_features[mp_features %in% colnames(features)]

state_display <- c(
  "Classic Proliferative" = "state__Classic_Proliferative",
  "Basal to Intestinal Metaplasia" = "state__Basal_to_Intestinal_Metaplasia",
  "SMG-like Metaplasia" = "state__SMG_like_Metaplasia",
  "Stress-adaptive" = "state__Stress_adaptive",
  "Immune Infiltrating" = "state__Immune_Infiltrating",
  "3CA EMT/Protein maturation" = "state__3CA_EMT_and_Protein_maturation"
)
state_features <- unname(state_display[state_display %in% colnames(features)])

qc_features <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "cc_score", "cs_score",
    "cna_signal", "cna_cor", "cna_burden_mean_abs", "n_gain_arms", "n_loss_arms"),
  colnames(features)
)

feature_group <- c(
  setNames(rep("Metaprogrammes", length(mp_features)), mp_features),
  setNames(rep("Six states", length(state_features)), state_features),
  setNames(rep("QC / CNA metrics", length(qc_features)), qc_features)
)
plot_features <- names(feature_group)

feature_label <- function(feature) {
  out <- as.character(feature)
  is_mp <- grepl("^mp__", out)
  mp <- sub("^mp__", "", out[is_mp])
  out[is_mp] <- paste0(mp, " ", mp_descriptions[mp])
  state_lookup <- setNames(names(state_display), unname(state_display))
  is_state <- out %in% names(state_lookup)
  out[is_state] <- state_lookup[out[is_state]]
  out <- gsub("^nCount_RNA$", "nCount_RNA", out)
  out <- gsub("^nFeature_RNA$", "nFeature_RNA", out)
  out <- gsub("^percent.mt$", "percent.mt", out)
  out <- gsub("^cc_score$", "CC score", out)
  out <- gsub("^cs_score$", "CS score", out)
  out <- gsub("^cna_signal$", "CNA signal", out)
  out <- gsub("^cna_cor$", "CNA correlation", out)
  out <- gsub("^cna_burden_mean_abs$", "CNA burden", out)
  out <- gsub("^n_gain_arms$", "No. gained arms", out)
  out <- gsub("^n_loss_arms$", "No. lost arms", out)
  out
}

event_label <- function(event_id, known_genes = NULL) {
  direction <- sub("_.*$", "", event_id)
  arm <- sub("^[^_]+_", "", event_id)
  base <- paste0(ifelse(direction == "gain", "Gain ", "Loss "), arm)
  if (!is.null(known_genes)) {
    kg <- as.character(known_genes)
    kg[is.na(kg)] <- ""
    kg <- ifelse(nzchar(kg), paste0(" (", substr(kg, 1, 36), ")"), "")
    return(paste0(base, kg))
  }
  base
}

sample_order_from_id <- function(x) {
  factor(x, levels = sort(unique(as.character(x))))
}

features <- features %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = as.character(.data$subclone_id)
  )

arm_long <- arm_long %>%
  mutate(
    sample = as.character(.data$sample),
    subclone = as.character(.data$subclone),
    subclone_id = as.character(.data$subclone_id),
    arm = as.character(.data$arm),
    cna_call_legacy = .data$cna_call,
    cna_call = case_when(
      .data$arm_mean >= event_call_threshold ~ 1L,
      .data$arm_mean <= -event_call_threshold ~ -1L,
      TRUE ~ 0L
    )
  )

event_summary_legacy <- event_summary
event_known_genes <- event_summary_legacy %>%
  select(.data$event_id, .data$known_genes) %>%
  distinct()

summarise_events <- function(direction_label, call_value) {
  arm_long %>%
    filter(.data$cna_call == call_value) %>%
    group_by(.data$arm) %>%
    summarise(
      direction = direction_label,
      n_samples_event = n_distinct(.data$sample),
      n_subclones_event = n_distinct(.data$subclone_id),
      samples_event = paste(sort(unique(.data$sample)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(event_id = paste0(.data$direction, "_", .data$arm))
}

event_summary <- bind_rows(
  summarise_events("gain", 1L),
  summarise_events("loss", -1L)
) %>%
  mutate(
    n_samples_total = n_distinct(features$sample),
    n_subclones_total = nrow(features),
    frac_samples_event = .data$n_samples_event / .data$n_samples_total,
    frac_subclones_event = .data$n_subclones_event / .data$n_subclones_total,
    is_recurrent = .data$n_samples_event >= recurrent_min_samples &
      .data$frac_samples_event >= recurrent_min_sample_fraction
  ) %>%
  left_join(event_known_genes, by = "event_id") %>%
  mutate(known_genes = ifelse(is.na(.data$known_genes), "", .data$known_genes)) %>%
  arrange(desc(.data$is_recurrent), desc(.data$frac_samples_event), desc(.data$frac_subclones_event))

threshold_sensitivity <- bind_rows(lapply(c(0, 0.02, 0.03, 0.05, 0.07, 0.08, 0.10, 0.12), function(thr) {
  bind_rows(
    arm_long %>%
      group_by(.data$arm) %>%
      summarise(
        threshold = thr,
        direction = "gain",
        n_subclones_event = sum(.data$arm_mean >= .env$thr, na.rm = TRUE),
        n_samples_event = n_distinct(.data$sample[.data$arm_mean >= .env$thr]),
        .groups = "drop"
      ),
    arm_long %>%
      group_by(.data$arm) %>%
      summarise(
        threshold = thr,
        direction = "loss",
        n_subclones_event = sum(.data$arm_mean <= -.env$thr, na.rm = TRUE),
        n_samples_event = n_distinct(.data$sample[.data$arm_mean <= -.env$thr]),
        .groups = "drop"
      )
  ) %>%
    mutate(event_id = paste0(.data$direction, "_", .data$arm))
})) %>%
  mutate(
    n_samples_total = n_distinct(features$sample),
    n_subclones_total = nrow(features),
    frac_samples_event = .data$n_samples_event / .data$n_samples_total,
    frac_subclones_event = .data$n_subclones_event / .data$n_subclones_total
  )
write_table(threshold_sensitivity, "Auto_v2_cna_event_threshold_sensitivity.csv")
write_table(event_summary, "Auto_v2_recomputed_recurrent_cna_event_summary.csv")

ranked_events <- event_summary %>%
  arrange(desc(.data$is_recurrent), desc(.data$frac_samples_event), desc(.data$frac_subclones_event))
base_recurrent_n <- sum(ranked_events$is_recurrent, na.rm = TRUE)
n_events_to_plot <- min(nrow(ranked_events), 8L)
recurrent_events <- ranked_events %>%
  head(n_events_to_plot) %>%
  pull(.data$event_id)
boxplot_events <- head(recurrent_events, min(4L, length(recurrent_events)))

event_meta <- event_summary %>%
  filter(.data$event_id %in% recurrent_events) %>%
  mutate(
    event_label = event_label(.data$event_id, .data$known_genes),
    event_label = factor(.data$event_label, levels = event_label(.data$event_id, .data$known_genes))
  )

event_presence_for <- function(event_id) {
  event_direction <- sub("_.*$", "", event_id)
  event_arm <- sub("^[^_]+_", "", event_id)
  arm_long %>%
    filter(.data$arm == .env$event_arm) %>%
    transmute(
      event_id = .env$event_id,
      sample = .data$sample,
      subclone = .data$subclone,
      subclone_id = .data$subclone_id,
      event_present = if (.env$event_direction == "gain") .data$cna_call == 1L else .data$cna_call == -1L
    ) %>%
    distinct(.data$event_id, .data$sample, .data$subclone, .data$subclone_id, .keep_all = TRUE)
}

event_presence <- bind_rows(lapply(recurrent_events, event_presence_for)) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

write_table(event_presence, "Auto_v2_recurrent_cna_event_subclone_presence.csv")

feature_sd <- vapply(plot_features, function(feat) safe_sd(features[[feat]]), numeric(1))
feature_sd[!is.finite(feature_sd) | feature_sd == 0] <- NA_real_
feature_mean <- vapply(plot_features, function(feat) safe_mean(features[[feat]]), numeric(1))

message("Computing corrected recurrent CNA event feature tests")

event_feature_values <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_df <- event_presence %>%
    filter(.data$event_id == ev) %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    value <- suppressWarnings(as.numeric(ev_df[[feat]]))
    data.frame(
      event_id = ev,
      sample = ev_df$sample,
      subclone = ev_df$subclone,
      subclone_id = ev_df$subclone_id,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      event_present = ev_df$event_present,
      feature_value = value,
      feature_z = (value - feature_mean[[feat]]) / feature_sd[[feat]],
      n_cells = ev_df$n_cells,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

event_sample_deltas <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_df <- event_presence %>%
    filter(.data$event_id == ev) %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    bind_rows(lapply(split(ev_df, ev_df$sample), function(sdf) {
      if (!all(c(TRUE, FALSE) %in% sdf$event_present)) return(NULL)
      event_value <- safe_weighted_mean(sdf[[feat]][sdf$event_present], sdf$n_cells[sdf$event_present])
      no_event_value <- safe_weighted_mean(sdf[[feat]][!sdf$event_present], sdf$n_cells[!sdf$event_present])
      data.frame(
        event_id = ev,
        sample = sdf$sample[1],
        feature = feat,
        feature_label = feature_label(feat),
        feature_group = feature_group[[feat]],
        event_value = event_value,
        no_event_value = no_event_value,
        delta = event_value - no_event_value,
        std_delta = (event_value - no_event_value) / feature_sd[[feat]],
        n_event_subclones = sum(sdf$event_present, na.rm = TRUE),
        n_no_event_subclones = sum(!sdf$event_present, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id")

event_feature_tests_v2 <- bind_rows(lapply(recurrent_events, function(ev) {
  ev_presence <- event_presence %>% filter(.data$event_id == ev)
  ev_df <- ev_presence %>%
    left_join(features, by = c("sample", "subclone", "subclone_id"))
  bind_rows(lapply(plot_features, function(feat) {
    d <- event_sample_deltas %>% filter(.data$event_id == ev, .data$feature == feat)
    event_values <- suppressWarnings(as.numeric(ev_df[[feat]][ev_df$event_present]))
    no_event_values <- suppressWarnings(as.numeric(ev_df[[feat]][!ev_df$event_present]))
    unpaired_p <- if (n_distinct(ev_df$event_present) == 2) {
      tryCatch(wilcox.test(event_values, no_event_values)$p.value,
               error = function(e) NA_real_)
    } else {
      NA_real_
    }
    paired_p <- if (nrow(d) >= 3 && sum(is.finite(d$delta)) >= 3) {
      tryCatch(wilcox.test(d$delta, mu = 0)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    data.frame(
      event_id = ev,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      n_subclones_event = sum(ev_df$event_present, na.rm = TRUE),
      n_subclones_no_event = sum(!ev_df$event_present, na.rm = TRUE),
      n_paired_samples = nrow(d),
      mean_event = safe_mean(event_values),
      mean_no_event = safe_mean(no_event_values),
      median_event = safe_median(event_values),
      median_no_event = safe_median(no_event_values),
      unpaired_delta = safe_mean(event_values) - safe_mean(no_event_values),
      unpaired_median_delta = safe_median(event_values) - safe_median(no_event_values),
      unpaired_median_std_delta = (safe_median(event_values) - safe_median(no_event_values)) / feature_sd[[feat]],
      paired_median_delta = if (nrow(d) > 0) median(d$delta, na.rm = TRUE) else NA_real_,
      paired_mean_delta = if (nrow(d) > 0) mean(d$delta, na.rm = TRUE) else NA_real_,
      paired_median_std_delta = if (nrow(d) > 0) median(d$std_delta, na.rm = TRUE) else NA_real_,
      unpaired_p_value = unpaired_p,
      paired_p_value = paired_p,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  left_join(event_meta %>% select(.data$event_id, .data$event_label, .data$direction, .data$arm,
                                  .data$known_genes, .data$frac_samples_event),
            by = "event_id") %>%
  group_by(.data$feature_group) %>%
  mutate(
    paired_p_adj_group = p.adjust(.data$paired_p_value, method = "BH"),
    unpaired_p_adj_group = p.adjust(.data$unpaired_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    paired_p_adj_global = p.adjust(.data$paired_p_value, method = "BH"),
    unpaired_p_adj_global = p.adjust(.data$unpaired_p_value, method = "BH"),
    sig_label = p_to_stars(.data$unpaired_p_adj_group),
    neglog10_fdr = pmin(-log10(pmax(.data$unpaired_p_adj_group, 1e-12)), 12),
    primary_delta = .data$unpaired_median_std_delta,
    primary_p_adj_group = .data$unpaired_p_adj_group
  )

write_table(event_feature_values, "Auto_v2_recurrent_cna_event_feature_values.csv")
write_table(event_sample_deltas, "Auto_v2_recurrent_cna_event_per_sample_feature_deltas.csv")
write_table(event_feature_tests_v2, "Auto_v2_recurrent_cna_event_feature_tests.csv")

message("Computing largest-subclone effects with standardized x-axis")

dominant_deltas_v2 <- bind_rows(lapply(plot_features, function(feat) {
  bind_rows(lapply(split(features, features$sample), function(df) {
    if (nrow(df) < 2) return(NULL)
    top <- df %>% filter(.data$is_largest_subclone)
    rest <- df %>% filter(!.data$is_largest_subclone)
    if (nrow(top) != 1 || nrow(rest) == 0) return(NULL)
    rest_value <- safe_weighted_mean(rest[[feat]], rest$n_cells)
    data.frame(
      sample = top$sample,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      dominant_value = top[[feat]],
      rest_weighted_value = rest_value,
      delta = top[[feat]] - rest_value,
      std_delta = (top[[feat]] - rest_value) / feature_sd[[feat]],
      dominance_class = top$dominance_class,
      stringsAsFactors = FALSE
    )
  }))
}))

dominant_tests_v2 <- dominant_deltas_v2 %>%
  group_by(.data$feature, .data$feature_label, .data$feature_group) %>%
  summarise(
    n_samples = n(),
    median_delta = median(.data$delta, na.rm = TRUE),
    median_std_delta = median(.data$std_delta, na.rm = TRUE),
    pct_positive_delta = mean(.data$delta > 0, na.rm = TRUE),
    wilcox_p_value = if (n() >= 3 && sum(is.finite(.data$delta)) >= 3) {
      tryCatch(wilcox.test(.data$delta, mu = 0)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(.data$feature_group) %>%
  mutate(wilcox_p_adj_group = p.adjust(.data$wilcox_p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    wilcox_p_adj_global = p.adjust(.data$wilcox_p_value, method = "BH"),
    sig_label = p_to_stars(.data$wilcox_p_adj_group),
    neglog10_fdr = pmin(-log10(pmax(.data$wilcox_p_adj_group, 1e-12)), 12)
  )

write_table(dominant_deltas_v2, "Auto_v2_largest_subclone_per_sample_feature_deltas.csv")
write_table(dominant_tests_v2, "Auto_v2_largest_subclone_feature_tests.csv")

message("Computing pairwise CNA-distance tests for all feature groups")

cna_distance_metrics <- c("cna_abs_mean", "cna_euclidean", "cna_call_discordance")
cna_metric_labels <- c(
  cna_abs_mean = "Mean absolute CNA distance",
  cna_euclidean = "Euclidean CNA distance",
  cna_call_discordance = "Arm-call discordance"
)

pairwise_feature_tests_v2 <- bind_rows(lapply(cna_distance_metrics, function(cna_metric) {
  bind_rows(lapply(plot_features, function(feat) {
    endpoint <- paste0("abs_delta__", feat)
    if (!endpoint %in% colnames(pairwise_df)) return(NULL)
    df <- pairwise_df %>%
      filter(is.finite(.data[[cna_metric]]), is.finite(.data[[endpoint]]))
    cor_res <- if (nrow(df) >= 5 && safe_sd(df[[cna_metric]]) > 0 && safe_sd(df[[endpoint]]) > 0) {
      tryCatch(cor.test(df[[cna_metric]], df[[endpoint]], method = "spearman"), error = function(e) NULL)
    } else {
      NULL
    }
    centered <- df %>%
      group_by(.data$sample) %>%
      mutate(
        cna_centered = .data[[cna_metric]] - mean(.data[[cna_metric]], na.rm = TRUE),
        endpoint_centered = .data[[endpoint]] - mean(.data[[endpoint]], na.rm = TRUE)
      ) %>%
      ungroup()
    centered_cor <- if (nrow(centered) >= 5 &&
                        safe_sd(centered$cna_centered) > 0 &&
                        safe_sd(centered$endpoint_centered) > 0) {
      tryCatch(cor.test(centered$cna_centered, centered$endpoint_centered, method = "spearman"),
               error = function(e) NULL)
    } else {
      NULL
    }
    data.frame(
      cna_metric = cna_metric,
      cna_metric_label = cna_metric_labels[[cna_metric]],
      endpoint = endpoint,
      feature = feat,
      feature_label = feature_label(feat),
      feature_group = feature_group[[feat]],
      n_pairs = nrow(df),
      spearman_rho = if (!is.null(cor_res)) unname(cor_res$estimate) else NA_real_,
      spearman_p_value = if (!is.null(cor_res)) cor_res$p.value else NA_real_,
      sample_centered_rho = if (!is.null(centered_cor)) unname(centered_cor$estimate) else NA_real_,
      sample_centered_p_value = if (!is.null(centered_cor)) centered_cor$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  group_by(.data$feature_group, .data$cna_metric) %>%
  mutate(
    sample_centered_p_adj_group = p.adjust(.data$sample_centered_p_value, method = "BH"),
    spearman_p_adj_group = p.adjust(.data$spearman_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    sample_centered_p_adj_global = p.adjust(.data$sample_centered_p_value, method = "BH"),
    spearman_p_adj_global = p.adjust(.data$spearman_p_value, method = "BH"),
    sig_label = p_to_stars(.data$sample_centered_p_adj_group),
    neglog10_fdr = pmin(-log10(pmax(.data$sample_centered_p_adj_group, 1e-12)), 12)
  )

write_table(pairwise_feature_tests_v2, "Auto_v2_pairwise_cna_distance_all_feature_tests.csv")

message("Creating v2 consensus heatmap")

chrom_levels <- c(paste0("chr", 1:22), "chrX")
arm_levels <- colnames(arm_matrix)
chr_from_arm <- sub("[pq]$", "", arm_levels)
arm_matrix[!is.finite(arm_matrix)] <- 0
valid_rows <- rownames(arm_matrix)[rowSums(is.finite(arm_matrix)) == ncol(arm_matrix)]
arm_matrix_valid <- arm_matrix[valid_rows, , drop = FALSE]
arm_call_wide_v2 <- arm_long %>%
  select(.data$subclone_id, .data$arm, .data$cna_call) %>%
  pivot_wider(names_from = .data$arm, values_from = .data$cna_call, values_fill = 0) %>%
  as.data.frame()
rownames(arm_call_wide_v2) <- arm_call_wide_v2$subclone_id
arm_call_matrix_v2 <- as.matrix(arm_call_wide_v2[, setdiff(colnames(arm_call_wide_v2), "subclone_id"), drop = FALSE])
missing_call_arms <- setdiff(arm_levels, colnames(arm_call_matrix_v2))
if (length(missing_call_arms) > 0) {
  arm_call_matrix_v2 <- cbind(
    arm_call_matrix_v2,
    matrix(0L, nrow = nrow(arm_call_matrix_v2), ncol = length(missing_call_arms),
           dimnames = list(rownames(arm_call_matrix_v2), missing_call_arms))
  )
}
arm_call_matrix_v2 <- arm_call_matrix_v2[rownames(arm_matrix_valid), arm_levels, drop = FALSE]
hc <- hclust(dist(arm_matrix_valid), method = "ward.D2")
cluster_df <- data.frame(
  subclone_id = names(cna_cluster),
  cna_cluster = unname(cna_cluster),
  stringsAsFactors = FALSE
) %>%
  left_join(features %>% select(.data$subclone_id, .data$sample, .data$subclone,
                                .data$subclone_fraction, .data$dominance_class),
            by = "subclone_id")
cluster_df <- cluster_df[match(rownames(arm_matrix_valid), cluster_df$subclone_id), , drop = FALSE]

cluster_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(cluster_df$cna_cluster))),
  sort(unique(cluster_df$cna_cluster))
)
dominance_cols <- c(
  "single_subclone" = "grey70",
  "largest_not_significant" = "#FDB863",
  "significant_dominant" = "#B2182B"
)

plot_consensus_heatmap_v2 <- function() {
  row_meta <- cluster_df[hc$order, , drop = FALSE]
  mat <- arm_matrix_valid[hc$order, arm_levels, drop = FALSE]
  row_ha <- rowAnnotation(
    Cluster = row_meta$cna_cluster,
    Dominance = row_meta$dominance_class,
    `Clone fraction` = row_meta$subclone_fraction,
    col = list(
      Cluster = cluster_cols,
      Dominance = dominance_cols,
      `Clone fraction` = colorRamp2(c(0, 0.5, 1), c("white", "#FDB863", "#B2182B"))
    ),
    annotation_name_gp = gpar(fontsize = 13, fontface = "bold"),
    annotation_legend_param = list(
      Cluster = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      Dominance = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
      `Clone fraction` = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12))
    ),
    simple_anno_size = unit(5, "mm")
  )
  Heatmap(
    mat,
    name = "Mean CNA",
    col = colorRamp2(c(-0.18, 0, 0.18), c("#2166AC", "white", "#B2182B")),
    left_annotation = row_ha,
    row_split = factor(row_meta$cna_cluster, levels = sort(unique(row_meta$cna_cluster))),
    row_title = NULL,
    row_title_gp = gpar(fontsize = 0, col = NA),
    row_gap = unit(1.2, "mm"),
    column_split = factor(chr_from_arm, levels = chrom_levels),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA)
  )
}

pdf(file.path(figure_dir, "Auto_v2_cna_consensus_heatmap_no_row_labels.pdf"),
    width = 17, height = 10, useDingbats = FALSE)
draw(plot_consensus_heatmap_v2(), heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

prepare_cna_matrix_v2 <- function(outs, cells) {
  cells <- intersect(cells, colnames(outs))
  common_genes <- intersect(rownames(outs), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(outs[go$gene, cells, drop = FALSE])
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  signal <- rowMeans(abs(mat), na.rm = TRUE)
  keep_signal <- signal >= as.numeric(quantile(signal, probs = 1 / 3, na.rm = TRUE))
  list(mat = mat[keep_signal, , drop = FALSE], gene_order = go[keep_signal, , drop = FALSE])
}

make_binned_cna_v2 <- function(cna_mat, go, bin_size = 100L) {
  go2 <- go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bin_levels <- unique(go2$bin)
  bins <- split(go2$.row, factor(go2$bin, levels = bin_levels))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  row_chr <- sub("_.*$", "", rownames(binned))
  list(mat = binned, chr = row_chr)
}

order_cells_by_subclone_v2 <- function(cna_mat, subclone) {
  split_cells <- split(names(subclone), factor(subclone, levels = unique(subclone)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    d <- dist(t(cna_mat[, cells, drop = FALSE]))
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

sample_plot_cells_v2 <- function(cells, subclone, max_cells = 1200L) {
  cells <- intersect(cells, names(subclone))
  if (length(cells) <= max_cells) return(cells)
  split_cells <- split(cells, factor(subclone[cells], levels = unique(subclone[cells])))
  target <- pmax(20L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

event_matrix_for_cells <- function(meta_plot, event_ids) {
  subclone_ids <- paste(meta_plot$sample, meta_plot$subclone, sep = "::")
  out <- do.call(rbind, lapply(event_ids, function(ev) {
    event_direction <- sub("_.*$", "", ev)
    event_arm <- sub("^[^_]+_", "", ev)
    calls <- arm_call_matrix_v2[subclone_ids, event_arm]
    present <- if (event_direction == "gain") calls == 1L else calls == -1L
    ifelse(present, ifelse(event_direction == "gain", 1L, -1L), 0L)
  }))
  rownames(out) <- as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])
  colnames(out) <- rownames(meta_plot)
  out
}

plot_sample_cell_heatmap <- function(sample_id) {
  sample_cells <- cell_df %>%
    filter(.data$sample == sample_id, .data$subclone != "Unresolved")
  if (nrow(sample_cells) == 0) return(NULL)
  outs_path <- file.path("by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  if (!file.exists(outs_path)) return(NULL)
  outs <- readRDS(outs_path)
  cells <- intersect(sample_cells$cell, colnames(outs))
  if (length(cells) == 0) return(NULL)
  sample_cells <- sample_cells[match(cells, sample_cells$cell), , drop = FALSE]
  subclone <- sample_cells$subclone
  names(subclone) <- sample_cells$cell
  cna_prepped <- prepare_cna_matrix_v2(outs, cells)
  cna_order <- order_cells_by_subclone_v2(cna_prepped$mat, subclone)
  plot_cells <- sample_plot_cells_v2(cna_order, subclone, max_plot_cells)
  binned <- make_binned_cna_v2(cna_prepped$mat[, cna_order, drop = FALSE], cna_prepped$gene_order)
  meta_plot <- sample_cells[match(plot_cells, sample_cells$cell), , drop = FALSE]
  rownames(meta_plot) <- meta_plot$cell
  meta_plot$subclone <- factor(meta_plot$subclone, levels = unique(subclone[cna_order]))
  meta_plot$top_mp_label <- factor(meta_plot$top_mp_label, levels = unique(meta_plot$top_mp_label))
  meta_plot$state_label <- factor(meta_plot$state_label, levels = unique(meta_plot$state_label))
  mat <- binned$mat[, rownames(meta_plot), drop = FALSE]
  row_chr <- factor(binned$chr, levels = unique(binned$chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr))), levels(row_chr))
  subclone_cols <- subclone_colours(meta_plot$subclone)
  topmp_cols <- mp_cols[names(mp_cols) %in% unique(as.character(meta_plot$top_mp_label))]
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(meta_plot$state_label))]
  missing_topmp <- setdiff(unique(as.character(meta_plot$top_mp_label)), names(topmp_cols))
  if (length(missing_topmp) > 0) topmp_cols <- c(topmp_cols, setNames(hue_pal()(length(missing_topmp)), missing_topmp))
  missing_states <- setdiff(unique(as.character(meta_plot$state_label)), names(local_state_cols))
  if (length(missing_states) > 0) local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  top_ha <- HeatmapAnnotation(
    Subclone = meta_plot$subclone,
    TopMP = meta_plot$top_mp_label,
    State = meta_plot$state_label,
    col = list(Subclone = subclone_cols, TopMP = topmp_cols, State = local_state_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(3, "mm"),
    na_col = "grey90"
  )
  left_ha <- rowAnnotation(
    Chr = row_chr,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )
  cna_ht <- Heatmap(
    mat,
    name = "CNA",
    col = colorRamp2(c(-cna_colour_limit, 0, cna_colour_limit), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_ha,
    left_annotation = left_ha,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    column_split = factor(meta_plot$subclone, levels = unique(meta_plot$subclone)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_title_gp = gpar(fontsize = 7),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA),
    column_title = sample_id
  )
  event_mat <- event_matrix_for_cells(meta_plot, recurrent_events)
  event_ht <- Heatmap(
    event_mat,
    name = "Event",
    col = c("-1" = "#2166AC", "0" = "white", "1" = "#B2182B"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      labels = c("Loss present", "Absent", "Gain present"),
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    height = unit(max(2.2, 0.32 * nrow(event_mat)), "cm"),
    border = TRUE,
    rect_gp = gpar(col = "grey70", lwd = 0.4)
  )
  cna_ht %v% event_ht
}

pdf(file.path(figure_dir, "Auto_v2_per_sample_heatmap_recurrent_events.pdf"),
    width = 14, height = 9, useDingbats = FALSE)
for (sample_id in sort(unique(features$sample))) {
  ht <- plot_sample_cell_heatmap(sample_id)
  if (!is.null(ht)) draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
}
dev.off()

message("Creating v2 recurrent CNA event plots")

event_order <- event_meta$event_label
event_feature_tests_v2 <- event_feature_tests_v2 %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features))))
  )
event_sample_deltas <- event_sample_deltas %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features))))
  )
event_feature_values <- event_feature_values %>%
  mutate(
    event_label = factor(.data$event_label, levels = event_order),
    feature_label = factor(.data$feature_label, levels = rev(unique(feature_label(plot_features)))),
    event_group = factor(ifelse(.data$event_present, "Event-positive", "Event-negative"),
                         levels = c("Event-negative", "Event-positive"))
  )

plot_event_assoc <- function(group_name, title_suffix) {
  df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  ggplot(df, aes(.data$event_label, .data$feature_label)) +
    geom_point(aes(size = .data$neglog10_fdr, fill = .data$primary_delta),
               shape = 21, color = "black", stroke = 0.45, alpha = 0.95) +
    geom_text(aes(label = .data$sig_label), size = 6.2, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         na.value = "grey90") +
    scale_size_continuous(range = c(3.5, 11), limits = c(0, 12)) +
    labs(
      title = paste0("Recurrent CNA event associations: ", title_suffix),
      x = NULL,
      y = NULL,
      fill = "Median standardized\nevent delta",
      size = "-log10(FDR)"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15)
    )
}

event_page_chunks <- function(event_ids, page_size = 4L) {
  if (length(event_ids) == 0) return(list(character(0)))
  split(event_ids, ceiling(seq_along(event_ids) / page_size))
}

plot_event_boxplots <- function(group_name, title_suffix, event_ids) {
  df <- event_feature_values %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])))
  sig_df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)]))) %>%
    left_join(
      df %>%
        group_by(.data$event_id, .data$feature) %>%
        summarise(star_y = max(.data$feature_z, na.rm = TRUE) + 0.25, .groups = "drop"),
      by = c("event_id", "feature")
    )
  ggplot(df, aes(x = .data$feature_label, y = .data$feature_z, fill = .data$event_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.88,
                 position = position_dodge(width = 0.72), linewidth = 0.55) +
    geom_point(aes(color = .data$event_group),
               position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0.02,
                                               dodge.width = 0.72),
               alpha = 0.36, size = 1.3, show.legend = FALSE) +
    geom_text(data = sig_df, aes(x = .data$feature_label, y = .data$star_y, label = .data$sig_label),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~event_label, ncol = 2) +
    scale_fill_manual(values = c("Event-negative" = "grey72", "Event-positive" = "#B2182B")) +
    scale_color_manual(values = c("Event-negative" = "grey45", "Event-positive" = "#B2182B")) +
    labs(
      title = paste0("Recurrent CNA event feature distributions: ", title_suffix),
      x = NULL,
      y = "Standardized subclone feature value",
      fill = NULL
    ) +
    theme_classic(base_size = 20) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(face = "bold", size = 24),
      legend.position = "top",
      legend.text = element_text(size = 16)
    )
}

plot_event_sample_deltas <- function(group_name, title_suffix, event_ids) {
  df <- event_sample_deltas %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)])))
  sig_df <- event_feature_tests_v2 %>%
    filter(.data$feature_group == group_name, .data$event_id %in% event_ids) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = feature_label(plot_features[feature_group == group_name])),
           event_label = factor(as.character(.data$event_label),
                                levels = as.character(event_meta$event_label[match(event_ids, event_meta$event_id)]))) %>%
    left_join(
      df %>%
        group_by(.data$event_id, .data$feature) %>%
        summarise(star_y = max(.data$std_delta, na.rm = TRUE) + 0.20, .groups = "drop"),
      by = c("event_id", "feature")
    )
  ggplot(df, aes(x = .data$feature_label, y = .data$std_delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.55) +
    geom_boxplot(outlier.shape = NA, width = 0.62, fill = "grey84", color = "black", linewidth = 0.55) +
    geom_point(position = position_jitter(width = 0.12, height = 0), alpha = 0.48,
               size = 1.8, color = "#B2182B") +
    geom_text(data = sig_df, aes(x = .data$feature_label, y = .data$star_y, label = .data$sig_label),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~event_label, ncol = 2) +
    labs(
      title = paste0("Per-sample recurrent CNA event deltas: ", title_suffix),
      x = NULL,
      y = "Standardized paired delta"
    ) +
    theme_classic(base_size = 20) +
    theme(
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(face = "bold", size = 24)
    )
}

event_bar <- event_meta %>%
  mutate(event_label = factor(.data$event_label, levels = rev(as.character(.data$event_label)))) %>%
  ggplot(aes(.data$frac_samples_event, .data$event_label, fill = .data$direction)) +
  geom_col(color = "black", linewidth = 0.35, width = 0.72) +
  geom_text(aes(label = percent(.data$frac_samples_event, accuracy = 1)), hjust = -0.08, size = 5) +
  scale_x_continuous(labels = percent, limits = c(0, min(1, max(event_meta$frac_samples_event) * 1.18))) +
  scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC")) +
  labs(title = "Recurrent arm-level CNA events used for association testing",
       x = "Fraction of samples with event in at least one subclone", y = NULL, fill = NULL) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 24),
        axis.text = element_text(size = 16),
        legend.position = "top")

pdf(file.path(figure_dir, "Auto_v2_recurrent_cna_event_associations_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(event_bar)
print(plot_event_assoc("Metaprogrammes", "all metaprogrammes"))
print(plot_event_assoc("Six states", "six states excluding Hybrid and Unresolved"))
print(plot_event_assoc("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

pdf(file.path(figure_dir, "Auto_v2_recurrent_cna_event_boxplots_all_features.pdf"),
    width = 22, height = 14, useDingbats = FALSE)
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Metaprogrammes", "all metaprogrammes", event_page))
}
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("Six states", "six states excluding Hybrid and Unresolved", event_page))
}
for (event_page in event_page_chunks(boxplot_events, page_size = 4L)) {
  print(plot_event_boxplots("QC / CNA metrics", "QC and CNA metrics", event_page))
}
dev.off()

pdf(file.path(figure_dir, "Auto_v2_recurrent_cna_event_per_sample_deltas.pdf"),
    width = 22, height = 14, useDingbats = FALSE)
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("Metaprogrammes", "all metaprogrammes", event_page))
}
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("Six states", "six states excluding Hybrid and Unresolved", event_page))
}
for (event_page in event_page_chunks(recurrent_events, page_size = 4L)) {
  print(plot_event_sample_deltas("QC / CNA metrics", "QC and CNA metrics", event_page))
}
dev.off()

chr8q_myc <- arm_long %>%
  filter(.data$arm == "chr8q") %>%
  select(.data$sample, .data$subclone, .data$subclone_id, .data$cna_call) %>%
  left_join(features %>% select(.data$sample, .data$subclone, .data$subclone_id, .data$mp__MP2, .data$n_cells),
            by = c("sample", "subclone", "subclone_id")) %>%
  mutate(
    chr8q_group = case_when(
      .data$cna_call == 1L ~ "8q gain",
      .data$cna_call == -1L ~ "8q loss",
      TRUE ~ "No 8q CNA"
    ),
    chr8q_group = factor(.data$chr8q_group, levels = c("8q loss", "No 8q CNA", "8q gain"))
  )
if (nrow(chr8q_myc) > 0) {
  group_n <- chr8q_myc %>%
    group_by(.data$chr8q_group) %>%
    summarise(n = n(), y = max(.data$mp__MP2, na.rm = TRUE) + 0.02, .groups = "drop")
  p_gain8q_myc <- ggplot(chr8q_myc, aes(.data$chr8q_group, .data$mp__MP2, fill = .data$chr8q_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.88, linewidth = 0.7, width = 0.62) +
    geom_point(aes(size = .data$n_cells), position = position_jitter(width = 0.14, height = 0),
               alpha = 0.42, color = "black") +
    geom_text(data = group_n, aes(x = .data$chr8q_group, y = .data$y, label = paste0("n=", .data$n)),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    scale_fill_manual(values = c("8q loss" = "#2166AC", "No 8q CNA" = "grey72", "8q gain" = "#B2182B")) +
    scale_size_continuous(range = c(1.8, 6)) +
    labs(
      title = "chr8q CNA status versus MYC-related proliferation MP",
      x = NULL,
      y = "Subclone mean MP2 score",
      fill = "chr8q CNA",
      size = "Cells"
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 24),
      axis.text = element_text(size = 18),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15),
      legend.position = "right"
    )
  ggsave(file.path(figure_dir, "Auto_v2_gain_chr8q_myc_mp_per_sample.pdf"),
         p_gain8q_myc, width = 12, height = 8, useDingbats = FALSE)
}

message("Creating v2 largest-subclone and pairwise-distance plots")

dominant_tests_v2 <- dominant_tests_v2 %>%
  mutate(
    feature_label = factor(.data$feature_label,
                           levels = rev(unique(feature_label(plot_features)))),
    feature_group = factor(.data$feature_group, levels = c("Metaprogrammes", "Six states", "QC / CNA metrics"))
  )

plot_dominant_group <- function(group_name, title_suffix) {
  df <- dominant_deltas_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  sig_df <- dominant_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))),
           star_x = max(df$std_delta, na.rm = TRUE) + 0.15)
  ggplot(df, aes(.data$std_delta, .data$feature_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey86", color = "black", linewidth = 0.65) +
    geom_point(aes(color = .data$dominance_class),
               position = position_jitter(height = 0.12, width = 0), alpha = 0.58, size = 2.2) +
    geom_text(data = sig_df, aes(x = .data$star_x, y = .data$feature_label, label = .data$sig_label),
              inherit.aes = FALSE, size = 6.5, fontface = "bold") +
    scale_color_manual(
      values = dominance_cols,
      breaks = c("single_subclone", "largest_not_significant", "significant_dominant"),
      labels = c("Single subclone", "Largest not dominant", "Dominant largest")
    ) +
    labs(
      title = paste0("Largest subclone versus other subclones: ", title_suffix),
      x = "Standardized largest-minus-rest delta",
      y = NULL,
      color = "Largest-subclone class"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15),
      legend.position = "top"
    )
}

pdf(file.path(figure_dir, "Auto_v2_largest_subclone_effects_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(plot_dominant_group("Metaprogrammes", "all metaprogrammes"))
print(plot_dominant_group("Six states", "six states excluding Hybrid and Unresolved"))
print(plot_dominant_group("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

pairwise_feature_tests_v2 <- pairwise_feature_tests_v2 %>%
  mutate(
    cna_metric_label = factor(.data$cna_metric_label, levels = cna_metric_labels[cna_distance_metrics]),
    feature_group = factor(.data$feature_group, levels = c("Metaprogrammes", "Six states", "QC / CNA metrics"))
  )

plot_pairwise_group <- function(group_name, title_suffix) {
  df <- pairwise_feature_tests_v2 %>%
    filter(.data$feature_group == group_name) %>%
    mutate(feature_label = factor(as.character(.data$feature_label),
                                  levels = rev(feature_label(plot_features[feature_group == group_name]))))
  ggplot(df, aes(.data$cna_metric_label, .data$feature_label)) +
    geom_point(aes(size = .data$neglog10_fdr, fill = .data$sample_centered_rho),
               shape = 21, color = "black", stroke = 0.45) +
    geom_text(aes(label = .data$sig_label), size = 6.1, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), na.value = "grey90") +
    scale_size_continuous(range = c(3.5, 11), limits = c(0, 12)) +
    labs(
      title = paste0("Pairwise CNA distance versus subclone divergence: ", title_suffix),
      x = NULL,
      y = NULL,
      fill = "Sample-centered rho",
      size = "-log10(FDR)"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 24),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 15)
    )
}

pdf(file.path(figure_dir, "Auto_v2_pairwise_cna_distance_all_features.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
print(plot_pairwise_group("Metaprogrammes", "all metaprogrammes"))
print(plot_pairwise_group("Six states", "six states excluding Hybrid and Unresolved"))
print(plot_pairwise_group("QC / CNA metrics", "QC and CNA metrics"))
dev.off()

run_summary_v2 <- data.frame(
  metric = c(
    "analysed_samples",
    "analysed_subclones",
    "recurrent_events_plotted",
    "features_plotted",
    "mp_features",
    "state_features_excluding_hybrid_unresolved",
    "qc_cna_features",
    "pairwise_subclone_comparisons"
  ),
  value = c(
    dplyr::n_distinct(features$sample),
    nrow(features),
    length(recurrent_events),
    length(plot_features),
    length(mp_features),
    length(state_features),
    length(qc_features),
    nrow(pairwise_df)
  )
)

write_table(run_summary_v2, "Auto_v2_run_summary.csv")

summary_dir <- file.path("..", "updates", "new_updates", "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
event_count_summary <- event_feature_tests_v2 %>%
  distinct(.data$event_id, .data$event_label, .data$n_subclones_event,
           .data$n_subclones_no_event, .data$n_paired_samples) %>%
  mutate(
    summary_type = "recurrent_event_count",
    item = .data$event_id,
    metric = "event_positive_event_negative_paired_samples",
    value = paste(.data$n_subclones_event, .data$n_subclones_no_event,
                  .data$n_paired_samples, sep = "/"),
    detail = as.character(.data$event_label)
  ) %>%
  select(.data$summary_type, .data$item, .data$metric, .data$value, .data$detail)

top_assoc_summary <- event_feature_tests_v2 %>%
  filter(is.finite(.data$unpaired_p_adj_group)) %>%
  arrange(.data$unpaired_p_adj_group) %>%
  head(10) %>%
  mutate(
    summary_type = "top_recurrent_event_association",
    item = paste(.data$event_id, .data$feature, sep = "__"),
    metric = "unpaired_median_standardized_delta_group_fdr",
    value = paste(signif(.data$unpaired_median_std_delta, 4),
                  signif(.data$unpaired_p_adj_group, 4), sep = "/"),
    detail = paste(.data$event_label, .data$feature_label, sep = " | ")
  ) %>%
  select(.data$summary_type, .data$item, .data$metric, .data$value, .data$detail)

compact_summary <- bind_rows(
  run_summary_v2 %>%
    mutate(
      summary_type = "run",
      item = .data$metric,
      detail = ""
    ) %>%
    transmute(.data$summary_type, .data$item, metric = "value",
              value = as.character(.data$value), .data$detail),
  event_count_summary,
  top_assoc_summary
)
write.csv(compact_summary,
          file.path(summary_dir, "Auto_cna_subclone_expression_v2_summary.csv"),
          row.names = FALSE)

save_rds(
  list(
    event_presence = event_presence,
    event_feature_values = event_feature_values,
    event_sample_deltas = event_sample_deltas,
    event_feature_tests = event_feature_tests_v2,
    dominant_deltas = dominant_deltas_v2,
    dominant_tests = dominant_tests_v2,
    pairwise_feature_tests = pairwise_feature_tests_v2
  ),
  "Auto_v2_visualisation_results.rds"
)

message("Done. V2 outputs saved under ", out_dir)
