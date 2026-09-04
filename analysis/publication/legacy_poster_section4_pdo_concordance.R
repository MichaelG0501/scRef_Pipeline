####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/publication/legacy_poster_section4_pdo_concordance.R
#   Methodology: not required (legacy poster assembly)
#   Inputs:
#     PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#     PDOs_Pipeline/PDOs_outs/Auto_PDO_final_states.rds
#     PDOs_Pipeline/PDOs_outs/UCell_scores_filtered.rds
#     PDOs_Pipeline/PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/*.csv
#     PDOs_Pipeline/PDOs_outs/Auto_pdo_sn_matched_pair_comparison/*.csv
#     ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   Outputs:
#     ref_outs/publication/section4/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section4_pdo_concordance.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(ggpubr)
  library(scales)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section4"
out_dir <- pub_section_dir(section)
pdo_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

# PDO MP descriptions and state mapping
pdo_mp_desc <- c(
  "MP6" = "G2M Cell Cycle", "MP7" = "DNA repair", "MP5" = "MYC-related Proliferation",
  "MP1" = "G2M checkpoint", "MP3" = "G1S Cell Cycle", "MP8" = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.", "MP9" = "ECM Remodeling Epi.", "MP4" = "Intestinal Metaplasia"
)
# PDO MP order: CC first, MYC, then lineage/stress
pdo_mp_order <- c("MP6", "MP3", "MP1", "MP7",   # Cell cycle
                   "MP5",                          # MYC / Classic Proliferative
                   "MP4",                          # Basal to Intestinal Metaplasia
                   "MP10", "MP9",                  # Stress-adaptive
                   "MP8")                          # SMG-like Metaplasia (Columnar Progenitor)
pdo_mp_to_state <- c(
  "MP6" = "Cell-cycle / DNA repair", "MP3" = "Cell-cycle / DNA repair",
  "MP1" = "Cell-cycle / DNA repair", "MP7" = "Cell-cycle / DNA repair",
  "MP5" = "Classic Proliferative",
  "MP4" = "Basal to Intestinal Metaplasia",
  "MP10" = "Stress-adaptive", "MP9" = "Stress-adaptive",
  "MP8" = "SMG-like Metaplasia"
)
pdo_state_order_ext <- c(PUB_STATE_ORDER, "Unresolved", "Hybrid")
pdo_state_colours_ext <- c(PUB_STATE_COLOURS,
                           "Unresolved" = "grey80",
                           "Hybrid" = "#111827")

####################
# FIGURE 1: Cross-system pan-cancer MP mean-score scatter.
# PDO NMF and PDO-MP correlation heatmaps are no longer active poster outputs.
####################
cat("Generating PDO cross-system scatter...\n")
ucell_3ca_pdo_file <- file.path(pdo_dir, "UCell_3CA_MPs.rds")
ucell_3ca_sc_file <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds"

if (file.exists(ucell_3ca_pdo_file) && file.exists(ucell_3ca_sc_file)) {
  ucell_3ca_pdo <- readRDS(ucell_3ca_pdo_file)
  ucell_3ca_sc <- readRDS(ucell_3ca_sc_file)
  
  pdo_cols <- colnames(ucell_3ca_pdo)
  sc_cols <- colnames(ucell_3ca_sc)
  common_cols <- intersect(pdo_cols, sc_cols)

  pdo_mean_scores <- colMeans(ucell_3ca_pdo[, common_cols, drop = FALSE], na.rm = TRUE)
  sc_mean_scores <- colMeans(ucell_3ca_sc[, common_cols, drop = FALSE], na.rm = TRUE)

  common_mps <- intersect(names(pdo_mean_scores), names(sc_mean_scores))
  comp_df <- data.frame(MP = common_mps, PDO_score = pdo_mean_scores[common_mps], scRef_score = sc_mean_scores[common_mps])

  # Add status and labeling threshold
  comp_df$Label <- ifelse(comp_df$PDO_score >= 0.1 | comp_df$scRef_score >= 0.1, sub("^X3CA_mp_", "", comp_df$MP), NA)
  comp_df$Status <- ifelse(comp_df$PDO_score < 0.1 & comp_df$scRef_score < 0.1, "Low", "Significant")

  # Determine max limit for synced axes
  max_limit <- max(c(comp_df$scRef_score, comp_df$PDO_score), na.rm = TRUE) * 1.05

  p <- ggplot(comp_df, aes(x = scRef_score, y = PDO_score)) +
    geom_vline(xintercept = 0.1, linetype = "dotted", color = "black", linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "black", linewidth = 0.4, alpha = 0.5) +
    geom_point(aes(color = Status), size = 2.0, alpha = 0.8) +
    scale_color_manual(values = c("Low" = "grey60", "Significant" = "black")) +
    geom_text_repel(aes(label = Label), size = 3, box.padding = 0.4, point.padding = 0.3, force = 2, max.overlaps = 15, min.segment.length = 0, na.rm = TRUE) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", fill = "red", alpha = 0.1) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    xlim(0, max_limit) + 
    ylim(0, max_limit) +
    labs(title = NULL,
         subtitle = NULL,
         x = "scATLAS mean Score",
         y = "PDO mean Score") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  save_pub_gg(p, section, "s4_crossdata_scatter", width = 5.8, height = 5.8)
  write_csv(comp_df, file.path(out_dir, "tables", "s4_crossdata_scatter.csv"))
} else {
  abort_missing_figure(section, "s4_crossdata_scatter", "PDO cross-system scatter",
                   "Missing PDO or scRef UCell_3CA_MPs.rds inputs.")
}

# ==============================================================================
# FIGURE 3: Side-by-side state proportion: scATLAS (left) vs PDO (right)
# ==============================================================================
cat("Generating state occurrence comparison...\n")
sc_states <- readRDS(SCREF_STATE_NOREG_RDS)
pdo_states_file <- file.path(pdo_dir, "Auto_PDO_states_noreg.rds")


if (file.exists(pdo_states_file)) {
  pdo_states <- readRDS(pdo_states_file)
  prop_df <- bind_rows(
    tibble(system = "scATLAS\ntumours", state = clean_state(sc_states)),
    tibble(system = "PDO\nmodel", state = clean_state(pdo_states))
  ) |>
    mutate(state = ifelse(is.na(state) | state == "", "Unresolved", state),
           state = clean_state(state)) |>
    count(system, state, name = "n") |>
    complete(system, state = pdo_state_order_ext, fill = list(n = 0)) |>
    group_by(system) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(system = factor(system, levels = c("scATLAS\ntumours", "PDO\nmodel")),
           state = factor(state, levels = rev(pdo_state_order_ext)))

  p <- ggplot(prop_df, aes(system, pct, fill = state)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = pdo_state_colours_ext, breaks = pdo_state_order_ext, name = "State") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Proportion (%)") +
    pub_theme(13) +
    theme(
      axis.title.y = element_text(face = "plain"),
      axis.text.x = element_text(face = "plain"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.5, "lines")
    )
  save_pub_gg(p, section, "s4_state_occurrence", width = 5.5, height = 5)
  write_csv(prop_df, file.path(out_dir, "tables", "s4_state_occurrence_pre_unresolved_relabel.csv"))
} else {
  abort_missing_figure(section, "s4_state_occurrence", "PDO state occurrence",
                   "Missing PDO state vector.")
}

# ==============================================================================
# FIGURE 5: Matched PDO vs snRNA-seq tumour concordance (before unresolved relabel)
# ==============================================================================
cat("Generating matched PDO-snseq concordance...\n")
pdo_states_raw_file <- file.path(pdo_dir, "Auto_PDO_states_noreg.rds")
sn_states_raw_file <- "/rds/general/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs/Auto_topmp_v2_noreg_states_B.rds"

if (file.exists(pdo_states_raw_file) && file.exists(sn_states_raw_file)) {
  pdo_states_raw <- readRDS(pdo_states_raw_file)
  sn_states_raw <- readRDS(sn_states_raw_file)

  pdo_sur680 <- pdo_states_raw[grepl("^SUR680T3_PDO", names(pdo_states_raw))]
  pdo_sur791 <- pdo_states_raw[grepl("^SUR791T3_PDO", names(pdo_states_raw))]
  sn_sur680 <- sn_states_raw[grepl("^H_post_T1_biopsy", names(sn_states_raw))]
  sn_sur791 <- sn_states_raw[grepl("^L_post_T1_biopsy", names(sn_states_raw))]

  make_matched_df <- function(state_vec, patient, mod, samp) {
    df <- data.frame(label = as.character(state_vec), stringsAsFactors = FALSE)
    df$label[df$label == "Basal to Intest. Meta"] <- "Basal to Intestinal Metaplasia"
    df$label[grepl("__", df$label)] <- "Hybrid"
    df |> 
      count(label) |> 
      mutate(
        patient_id = patient,
        modality = mod,
        sample = samp,
        total_n = sum(n),
        pct = 100 * n / total_n
      )
  }

  matched <- bind_rows(
    make_matched_df(pdo_sur680, "SUR680", "PDO\nmodel", "SUR680T3_PDO"),
    make_matched_df(sn_sur680, "SUR680", "snSeq matched\ntumour", "H_post_T1_biopsy"),
    make_matched_df(pdo_sur791, "SUR791", "PDO\nmodel", "SUR791T3_PDO"),
    make_matched_df(sn_sur791, "SUR791", "snSeq matched\ntumour", "L_post_T1_biopsy")
  ) |>
    mutate(label = clean_state(label),
           label = ifelse(is.na(label) | label == "", "Unresolved", label),
           label = factor(label, levels = rev(pdo_state_order_ext))) |>
    mutate(modality = factor(modality, levels = c("snSeq matched\ntumour", "PDO\nmodel")))

  p <- ggplot(matched, aes(x = modality, y = pct, fill = label)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    facet_wrap(~ patient_id, scales = "free_x") +
    scale_fill_manual(values = pdo_state_colours_ext, breaks = pdo_state_order_ext, drop = FALSE, name = "State") +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.03)),
      labels = function(x) paste0(x, "%")
    ) +
    labs(x = NULL, y = "Proportion (%)", fill = NULL) +
    pub_theme(13) +
    theme(
      axis.title.y = element_text(face = "plain"),
      axis.text.x = element_text(face = "plain"),
      strip.text = element_text(face = "plain", size = 13),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "none"
    )
  save_pub_gg(p, section, "s4_matched_concordance", width = 5.5, height = 5)
  write_csv(matched, file.path(out_dir, "tables", "s4_matched_concordance.csv"))
} else {
  abort_missing_figure(section, "s4_matched_concordance", "PDO state concordance",
                   "Missing PDO concordance inputs.")
}

cat("Section 4 complete.\n")
