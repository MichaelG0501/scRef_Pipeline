####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_tcga_survival_clinical_mps_v2_old.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_survival_clinical_mps_v2.R
# Structured survival and clinical association analysis for updated MPs/states.
#
# Goals:
# 1) 14 MP survival volcano plots (continuous Cox) in TCGA, split by EAC vs ESCC.
# 2) Kaplan-Meier plots for 5 states (state-level) in TCGA, split by EAC vs ESCC.
# 3) Clinical association tests (location, age, stage and other variables) for MPs/states.
#
# Input:  ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#         ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#         ref_outs/Auto_topmp_v2_states_B.rds
# Output: ref_outs/Auto_survival_tcga_volcano_EAC.pdf
#         ref_outs/Auto_survival_tcga_volcano_ESCC.pdf
#         ref_outs/Auto_survival_tcga_cox_results.csv
#         ref_outs/Auto_survival_tcga_state_km_EAC.pdf
#         ref_outs/Auto_survival_tcga_state_km_ESCC.pdf
#         ref_outs/Auto_clinical_assoc_mps_states.csv
#         ref_outs/Auto_clinical_assoc_EAC.pdf
#         ref_outs/Auto_clinical_assoc_ESCC.pdf
#         updates/DDmon/summaries/Auto_survival_clinical_mps_v2_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(GSVA)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

cat("=== Auto survival clinical MPs v2 ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load MP definitions (14 retained MPs after silhouette filter)
####################
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)
cat("Retained MPs:", paste(retained_mps, collapse = ", "), "\n")

####################
# Canonical MP display names (for clear reporting/plots)
####################
mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC Proliferation",
  "MP17" = "Basal-like Trans.",
  "MP14" = "Hypoxia Adapted",
  "MP5"  = "Epithelial IFN",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Intest.",
  "MP16" = "Secretory Gastric",
  "MP15" = "Immune Infilt.",
  "MP12" = "Neuro-responsive"
)

####################
# MP/state ordering (tree-based) for consistent plotting/reporting
####################
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

####################
# 2) Load TCGA metadata + TPM
####################
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

####################
# 3) Build GSVA input for retained 14 MPs
####################
gsva_sets <- lapply(mp.genes, unique)
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]

if (length(gsva_sets) == 0) stop("No valid GSVA sets after filtering")

cat(sprintf("Running GSVA for %d MP sets...\n", length(gsva_sets)))
gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

surv_data <- meta_tcga %>%
  inner_join(gsva_df, by = "sample_barcode") %>%
  filter(sample_type_code == "01")

####################
# 4) Harmonize histology labels: EAC vs ESCC vs Other
####################
infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}

surv_data$HistologyGroup <- infer_histology(surv_data$type)

cat("TCGA histology counts after join:\n")
print(table(surv_data$HistologyGroup, useNA = "ifany"))

####################
# 5) Continuous Cox for each MP, by cohort (EAC/ESCC only)
####################
run_cox_for_group <- function(df, group_name, mp_cols) {
  out <- list()
  for (mp in mp_cols) {
    if (!mp %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[mp]]))
    if (nrow(d) < 20) next
    if (var(d[[mp]], na.rm = TRUE) == 0) next

    fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", mp, "`")), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next

    ss <- summary(fit)
    out[[mp]] <- data.frame(
      cohort = group_name,
      MP = mp,
      beta = ss$coefficients[1, "coef"],
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  res <- bind_rows(out)
  res$padj <- p.adjust(res$P_value, method = "BH")
  res
}

mp_cols <- intersect(names(gsva_sets), colnames(surv_data))
cox_eac <- run_cox_for_group(surv_data %>% filter(HistologyGroup == "EAC"), "EAC", mp_cols)
cox_escc <- run_cox_for_group(surv_data %>% filter(HistologyGroup == "ESCC"), "ESCC", mp_cols)
cox_res <- bind_rows(cox_eac, cox_escc)
if (nrow(cox_res) > 0) {
  cox_res$MP_label <- ifelse(
    !is.na(mp_descriptions[cox_res$MP]),
    unname(mp_descriptions[cox_res$MP]),
    cox_res$MP
  )
}

if (nrow(cox_eac) == 0) cat("Warning: EAC cohort produced no valid Cox rows (check sample size/variance).\n")
if (nrow(cox_escc) == 0) cat("Warning: ESCC cohort produced no valid Cox rows (check tcga_data_prep.R and sample counts).\n")

write.csv(cox_res, "Auto_survival_tcga_cox_results.csv", row.names = FALSE)

plot_volcano <- function(df, title_text, out_pdf) {
  if (nrow(df) == 0) return(invisible(NULL))
  # use raw P-values for plotting/labeling (keep padj calculation in results)
  df <- df %>% mutate(sig = P_value < 0.05, neglog10 = -log10(P_value), log2HR = log2(HR))
  labs_df <- df %>% arrange(P_value) %>% filter(sig) %>% slice_head(n = 12)
  labs_df$MP_plot <- ifelse(!is.na(mp_descriptions[labs_df$MP]), unname(mp_descriptions[labs_df$MP]), labs_df$MP)
  df$MP_plot <- ifelse(!is.na(mp_descriptions[df$MP]), unname(mp_descriptions[df$MP]), df$MP)
  
  p <- ggplot(df, aes(x = log2HR, y = neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(data = df, aes(label = MP_plot), size = 3.2, max.overlaps = 20) +
    theme_minimal(base_size = 13) +
    labs(title = title_text, x = "log2(HR)", y = "-log10(p)")

  ggsave(out_pdf, p, width = 9, height = 7)
}

plot_volcano(cox_eac, "TCGA survival volcano (EAC)", "Auto_survival_tcga_volcano_EAC.pdf")
plot_volcano(cox_escc, "TCGA survival volcano (ESCC)", "Auto_survival_tcga_volcano_ESCC.pdf")

####################
# 6) State-level KM for 5 states (TCGA projection via state MPs)
####################
state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

for (nm in names(state_groups)) {
  mps <- intersect(state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}

state_cols <- ordered_group_names[ordered_group_names %in% colnames(surv_data)]

km_per_state <- function(df, cohort_name, out_pdf) {
  pdf(out_pdf, width = 8, height = 7, useDingbats = FALSE)
  for (st in state_cols) {
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[st]]))
    if (nrow(d) < 20) next
    med <- median(d[[st]], na.rm = TRUE)
    d$Group <- ifelse(d[[st]] >= med, "High", "Low")
    d$Group <- factor(d$Group, levels = c("Low", "High"))
    if (length(unique(d$Group)) < 2) next

    fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = d)
    # continuous Cox model p-value (state as continuous variable)
    p_cont <- tryCatch({
      fitc <- coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", st, "`")), data = d)
      ss <- summary(fitc)
      ss$coefficients[1, "Pr(>|z|)"]
    }, error = function(e) NA_real_)

    # build subtitle with continuous Cox p-value and sample sizes
    n_total <- nrow(d)
    n_low <- sum(d$Group == "Low")
    n_high <- sum(d$Group == "High")
    subtitle_text <- sprintf("Cox (continuous) p = %s | n = %d (Low=%d, High=%d)",
                            ifelse(is.na(p_cont), "NA", format.pval(p_cont, digits = 3)),
                            n_total, n_low, n_high)

    p <- ggsurvplot(
      fit,
      data = d,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = FALSE,
      title = paste0(cohort_name, " | ", st, " (median split)"),
      xlab = "Days",
      ylab = "Overall survival"
    )

    # attach subtitle with Cox continuous p and sample sizes
    p$plot <- p$plot + ggtitle(p$plot$labels$title, subtitle = subtitle_text)
    print(p)
  }
  dev.off()
}

km_per_state(surv_data %>% filter(HistologyGroup == "EAC"), "EAC", "Auto_survival_tcga_state_km_EAC.pdf")
km_per_state(surv_data %>% filter(HistologyGroup == "ESCC"), "ESCC", "Auto_survival_tcga_state_km_ESCC.pdf")

####################
# 7) Clinical associations for MPs and states (split by EAC/ESCC)
####################
candidate_numeric <- c("Age_at_diagnosis", "age_at_diagnosis", "OS_time", "Mutation_count", "TMB_nonsynonymous")
candidate_categorical <- c(
  "Location", "Stage_Simple", "Grade", "Gender", "Alcohol_history", "Ethnicity", "Race",
  "AJCC_pathologic_T", "AJCC_pathologic_N", "AJCC_pathologic_M", "Prior_treatment", "Prior_malignancy"
)

num_vars <- intersect(candidate_numeric, colnames(surv_data))
cat_vars <- intersect(candidate_categorical, colnames(surv_data))

# Explicit, non-ambiguous feature naming for downstream matching
feature_df <- dplyr::bind_rows(
  data.frame(
    feature = mp_tree_order_names[mp_tree_order_names %in% mp_cols],
    feature_type = "MP",
    stringsAsFactors = FALSE
  ),
  data.frame(
    feature = state_cols,
    feature_type = "State",
    stringsAsFactors = FALSE
  )
)
feature_df$feature_label <- ifelse(
  feature_df$feature_type == "MP",
  ifelse(!is.na(mp_descriptions[feature_df$feature]), unname(mp_descriptions[feature_df$feature]), feature_df$feature),
  feature_df$feature
)

assoc_rows <- list()
for (coh in c("EAC", "ESCC")) {
  surv_sub <- surv_data %>% filter(HistologyGroup == coh)
  for (i in seq_len(nrow(feature_df))) {
    feat <- feature_df$feature[i]
    ftype <- feature_df$feature_type[i]
    feature_label <- feature_df$feature_label[i]

    # numeric clinical variables -> Spearman
    for (v in num_vars) {
      d <- surv_sub %>% filter(!is.na(.data[[feat]]), !is.na(.data[[v]]))
      if (nrow(d) < 20) next
      tt <- suppressWarnings(cor.test(d[[feat]], d[[v]], method = "spearman"))
      assoc_rows[[paste(coh, feat, v, "num", sep = "|")]] <- data.frame(
        cohort = coh,
        feature = feature_label,
        feature_raw = feat,
        feature_type = ftype,
        feature_key = paste0(ftype, "::", feat),
        clinical_var = v,
        test_type = "spearman",
        estimate = unname(tt$estimate),
        statistic = unname(tt$statistic),
        p_value = tt$p.value,
        n = nrow(d),
        stringsAsFactors = FALSE
      )
    }

    # categorical clinical variables -> Kruskal-Wallis
    for (v in cat_vars) {
      d <- surv_sub %>% filter(!is.na(.data[[feat]]), !is.na(.data[[v]]))
      if (nrow(d) < 20) next
      if (length(unique(d[[v]])) < 2) next
      kw <- kruskal.test(as.formula(paste0("`", feat, "` ~ `", v, "`")), data = d)
      assoc_rows[[paste(coh, feat, v, "cat", sep = "|")]] <- data.frame(
        cohort = coh,
        feature = feature_label,
        feature_raw = feat,
        feature_type = ftype,
        feature_key = paste0(ftype, "::", feat),
        clinical_var = v,
        test_type = "kruskal",
        estimate = NA_real_,
        statistic = unname(kw$statistic),
        p_value = kw$p.value,
        n = nrow(d),
        stringsAsFactors = FALSE
      )
    }
  }
}

assoc_df <- if (length(assoc_rows) > 0) bind_rows(assoc_rows) else data.frame()
if (nrow(assoc_df) > 0) {
  assoc_df <- assoc_df %>%
    group_by(cohort) %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}
write.csv(assoc_df, "Auto_clinical_assoc_mps_states.csv", row.names = FALSE)

####################
# 7b) Clinical association visualisations (separate EAC/ESCC)
####################
plot_assoc_by_cohort <- function(df, cohort_name, out_pdf) {
  d <- df %>% filter(cohort == cohort_name)
  if (nrow(d) == 0) return(invisible(NULL))

  sig_label_fun <- function(p_value) {
    dplyr::case_when(
      is.na(p_value) ~ "NS",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  }

  d <- d %>%
    mutate(
      sig = sig_label_fun(p_value),
      neglog10p = -log10(p_value),
      neglog10padj = -log10(p_value),
      neglog10padj = ifelse(is.infinite(neglog10padj), NA_real_, neglog10padj),
      effect_plot = case_when(
        test_type == "spearman" ~ estimate,
        test_type == "kruskal" ~ statistic / sqrt(pmax(n, 1)),
        TRUE ~ NA_real_
      ),
      test_label = case_when(
        test_type == "spearman" ~ "Spearman",
        test_type == "kruskal" ~ "Kruskal-Wallis",
        TRUE ~ test_type
      ),
      feature_type = factor(feature_type, levels = c("MP", "State"))
    )

  # Only keep associations that are significant after adjustment
  d <- d %>% filter(!is.na(p_value) & p_value < 0.05)
  if (nrow(d) == 0) return(invisible(NULL))

  d <- d %>% filter(!is.na(neglog10padj))
  if (nrow(d) == 0) return(invisible(NULL))

  mp_label_order <- feature_df %>%
    filter(feature_type == "MP") %>%
    pull(feature_label)
  mp_label_order <- rev(mp_label_order)
  state_label_order <- feature_df %>%
    filter(feature_type == "State") %>%
    pull(feature_label)
  state_label_order <- rev(state_label_order)

  # Order clinical variables alphabetically
  clin_order <- sort(unique(d$clinical_var))

  d$feature <- factor(d$feature, levels = c(mp_label_order, state_label_order))
  d$clinical_var <- factor(d$clinical_var, levels = clin_order)

  p <- ggplot(d, aes(x = clinical_var, y = feature)) +
    geom_point(aes(size = neglog10padj, color = effect_plot), alpha = 0.9) +
    geom_text(
      data = subset(d, sig != "NS"),
      aes(label = sig),
      color = "black",
      size = 2.8,
      vjust = -0.9
    ) +
    scale_size_continuous(name = expression(-log[10]("adj. p")), range = c(2.5, 9)) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Effect") +
    facet_grid(feature_type ~ test_label, scales = "free", space = "free") +
    labs(
      title = paste0("Clinical associations (", cohort_name, ")"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 15),
      legend.position = "right"
    )

  ggsave(out_pdf, p, width = 14, height = 10)
}

plot_assoc_by_cohort(assoc_df, "EAC", "Auto_clinical_assoc_EAC.pdf")
plot_assoc_by_cohort(assoc_df, "ESCC", "Auto_clinical_assoc_ESCC.pdf")

####################
# 8) Compact summary for updates
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates",
  format(Sys.Date(), "%d%b"),
  "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "Auto_survival_clinical_mps_v2_summary.csv")

summary_df <- data.frame(
  metric = c(
    "n_tcga_primary",
    "n_eac",
    "n_escc",
    "n_mp_sets_tested",
    "n_signif_eac_pval",
    "n_signif_escc_pval",
    "n_clin_assoc_eac_pval",
    "n_clin_assoc_escc_pval"
  ),
  value = c(
    nrow(surv_data),
    sum(surv_data$HistologyGroup == "EAC", na.rm = TRUE),
    sum(surv_data$HistologyGroup == "ESCC", na.rm = TRUE),
    length(mp_cols),
    sum(cox_eac$P_value < 0.05, na.rm = TRUE),
    sum(cox_escc$P_value < 0.05, na.rm = TRUE),
    ifelse(nrow(assoc_df) > 0, sum(assoc_df$cohort == "EAC" & assoc_df$p_value < 0.05, na.rm = TRUE), 0),
    ifelse(nrow(assoc_df) > 0, sum(assoc_df$cohort == "ESCC" & assoc_df$p_value < 0.05, na.rm = TRUE), 0)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Saved: Auto_survival_tcga_volcano_EAC.pdf\n")
cat("Saved: Auto_survival_tcga_volcano_ESCC.pdf\n")
cat("Saved: Auto_survival_tcga_cox_results.csv\n")
cat("Saved: Auto_survival_tcga_state_km_EAC.pdf\n")
cat("Saved: Auto_survival_tcga_state_km_ESCC.pdf\n")
cat("Saved: Auto_clinical_assoc_mps_states.csv\n")
cat("Saved: Auto_clinical_assoc_EAC.pdf\n")
cat("Saved: Auto_clinical_assoc_ESCC.pdf\n")
cat(sprintf("Saved: %s\n", summary_path))
cat("=== Auto survival clinical MPs v2 complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
