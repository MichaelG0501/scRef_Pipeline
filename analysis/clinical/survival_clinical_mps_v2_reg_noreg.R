####################
# Auto_survival_clinical_mps_v2_reg_noreg.R
# Unified TCGA survival workflow for reg + noreg Approach B state contexts.
#
# Produces shared PDFs where each mode is rendered in separate pages:
#   - state volcano (EAC, ESCC)
#   - state KM (EAC, ESCC)
####################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(GSVA)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
state_reg <- readRDS("Auto_topmp_v2_reg_states_B.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

gsva_sets <- lapply(mp.genes, unique)
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]
gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

surv_data <- meta_tcga %>% inner_join(gsva_df, by = "sample_barcode") %>% filter(sample_type_code == "01")

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}
surv_data$HistologyGroup <- infer_histology(surv_data$type)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intest. Meta" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrated"    = c("MP15")
)
for (nm in names(state_groups)) {
  mps <- intersect(state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}
state_cols <- intersect(names(state_groups), colnames(surv_data))

run_cox_for_group <- function(df, features, cohort_name, mode_name) {
  out <- list()
  for (feat in features) {
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`")), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      mode = mode_name,
      cohort = cohort_name,
      feature = feat,
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, title_text) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(sig = P_value < 0.05, neglog10 = -log10(P_value), log2HR = log2(HR))
  ggplot(df, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.85) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(aes(label = feature), size = 3.1, max.overlaps = 100) +
    theme_minimal(base_size = 13) +
    labs(title = title_text, x = "log2(HR)", y = "-log10(p)")
}

all_cox <- list()

pdf("Auto_survival_tcga_state_volcano_reg_noreg.pdf", width = 9, height = 7)
for (mode_name in requested_modes) {
  for (coh in c("EAC", "ESCC")) {
    cox_df <- run_cox_for_group(
      surv_data %>% filter(HistologyGroup == coh),
      state_cols,
      cohort_name = coh,
      mode_name = mode_name
    )
    all_cox[[paste(mode_name, coh, sep = "_")]] <- cox_df
    p <- plot_volcano(cox_df, paste0("[", mode_name, "] TCGA state survival volcano (", coh, ")"))
    if (!is.null(p)) print(p)
  }
}
dev.off()

pdf("Auto_survival_tcga_state_km_reg_noreg.pdf", width = 8, height = 7)
for (mode_name in requested_modes) {
  for (coh in c("EAC", "ESCC")) {
    dcoh <- surv_data %>% filter(HistologyGroup == coh)
    for (st in state_cols) {
      d <- dcoh %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[st]]))
      if (nrow(d) < 20) next
      med <- median(d[[st]], na.rm = TRUE)
      d$Group <- factor(ifelse(d[[st]] >= med, "High", "Low"), levels = c("Low", "High"))
      if (length(unique(d$Group)) < 2) next
      fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = d)
      p <- ggsurvplot(
        fit,
        data = d,
        risk.table = TRUE,
        pval = TRUE,
        conf.int = FALSE,
        title = paste0("[", mode_name, "] ", coh, " | ", st, " (median split)"),
        xlab = "Days",
        ylab = "Overall survival"
      )
      print(p)
    }
  }
}
dev.off()

cox_res <- bind_rows(all_cox)
write.csv(cox_res, "Auto_survival_tcga_state_cox_results_reg_noreg.csv", row.names = FALSE)

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_df <- cox_res %>%
  group_by(mode, cohort) %>%
  summarise(n_sig_p005 = sum(P_value < 0.05, na.rm = TRUE), .groups = "drop")
write.csv(summary_df, file.path(summary_dir, "Auto_survival_clinical_mps_v2_reg_noreg_summary.csv"), row.names = FALSE)

message("Saved unified reg+noreg state survival volcano/KM outputs.")
