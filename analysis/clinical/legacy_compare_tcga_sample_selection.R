####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_compare_tcga_sample_selection.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_test_sample_type_filter.R
# Test whether restricting to primary tumour (sample_type_code == "01") vs
# including all sample types affects continuous Cox survival results (EAC only).
#
# Goals:
#   Compare two TCGA sample-type filtering strategies for EAC volcano plots:
#     all)  No sample_type_code filter (include all codes: 01, 06, 11, etc.)
#     01)   Only primary solid tumour (sample_type_code == "01") — current default
#
#   Uses all MP genes per MP (strategy c from weight-threshold test).
#   Produces two separate volcano plots (EAC) in a single side-by-side PDF.
#
# Input:  ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
# Output: ref_outs/Auto_test_sample_type_filter_EAC_volcanoes.pdf
#         updates/new_updates/summaries/Auto_test_sample_type_filter_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(GSVA)
library(patchwork)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

cat("=== Auto_test_sample_type_filter ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load MP definitions and apply silhouette filter
####################
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps  <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)
cat("Retained MPs:", paste(retained_mps, collapse = ", "), "\n")

####################
# Canonical MP display names
####################
mp_descriptions <- c(
  "MP1"  = "G2M_cycle",
  "MP2"  = "MYC_prolif",
  "MP5"  = "IFN_response",
  "MP7"  = "S_cycle",
  "MP8"  = "Intestinal_diff",
  "MP9"  = "G1S_cycle",
  "MP10" = "Columnar_diff",
  "MP12" = "Neuro_epithelial",
  "MP13" = "Partial_EMT",
  "MP14" = "Hypoxia_epithelial",
  "MP15" = "T_NK_infiltration",
  "MP16" = "Secretory_diff",
  "MP17" = "Squamous_transition",
  "MP18" = "Adaptive_secretory"
)

####################
# 2) Load TCGA data
####################
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df    <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat   <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

####################
# 3) Build GSVA gene sets (all genes — strategy c)
####################
gsva_sets <- lapply(mp.genes, unique)
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]

if (length(gsva_sets) == 0) stop("No valid GSVA sets after filtering")
cat(sprintf("Running GSVA on %d MP sets...\n", length(gsva_sets)))

gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df     <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

####################
# 4) Histology helper
####################
infer_histology <- function(type_vec) {
  t   <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno",    t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}

####################
# 5) Merge TCGA meta with GSVA scores (no sample_type_code filter yet)
####################
merged_all <- meta_tcga %>%
  inner_join(gsva_df, by = "sample_barcode")

merged_all$HistologyGroup <- infer_histology(merged_all$type)

cat("Sample type code distribution in merged data:\n")
print(table(merged_all$sample_type_code, useNA = "ifany"))

cat("\nHistology group counts (all sample types):\n")
print(table(merged_all$HistologyGroup, useNA = "ifany"))

####################
# 6) Subset strategies
####################
# all sample types
surv_all <- merged_all
# primary tumour only (current default)
surv_01  <- merged_all %>% filter(sample_type_code == "01")

cat(sprintf("\nEAC n (all types): %d\n", sum(surv_all$HistologyGroup == "EAC", na.rm = TRUE)))
cat(sprintf("EAC n (01 only):   %d\n", sum(surv_01$HistologyGroup  == "EAC", na.rm = TRUE)))

####################
# 7) Cox helper (returns data.frame of results for EAC subset)
####################
run_cox_eac <- function(surv_data, filter_label) {
  eac_data <- surv_data %>% filter(HistologyGroup == "EAC")
  mp_cols  <- intersect(names(gsva_sets), colnames(eac_data))

  cat(sprintf("  [%s] EAC n = %d, MPs available = %d\n",
              filter_label, nrow(eac_data), length(mp_cols)))

  out <- list()
  for (mp in mp_cols) {
    d <- eac_data %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[mp]]))
    if (nrow(d) < 20) next
    if (var(d[[mp]], na.rm = TRUE) == 0) next

    fit <- try(
      coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", mp, "`")), data = d),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next

    ss <- summary(fit)
    out[[mp]] <- data.frame(
      filter_label = filter_label,
      MP           = mp,
      beta         = ss$coefficients[1, "coef"],
      HR           = ss$coefficients[1, "exp(coef)"],
      P_value      = ss$coefficients[1, "Pr(>|z|)"],
      n            = fit$n,
      events       = fit$nevent,
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) return(data.frame())
  res      <- bind_rows(out)
  res$padj <- p.adjust(res$P_value, method = "BH")
  res$MP_label <- ifelse(
    !is.na(mp_descriptions[res$MP]),
    unname(mp_descriptions[res$MP]),
    res$MP
  )
  res
}

####################
# 8) Run both strategies
####################
cox_all <- run_cox_eac(surv_all, "All sample types")
cox_01  <- run_cox_eac(surv_01,  "Primary only (01)")

####################
# 9) Volcano plot function (returns ggplot object)
####################
make_volcano <- function(df, title_text) {
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No data", size = 6) +
        theme_void() +
        ggtitle(title_text)
    )
  }

  df <- df %>%
    mutate(
      sig      = P_value < 0.05,
      neglog10 = -log10(P_value),
      log2HR   = log2(HR)
    )

  p <- ggplot(df, aes(x = log2HR, y = neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "firebrick3"),
      guide  = "none"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype   = "dashed",
      color      = "grey45",
      linewidth  = 0.4
    ) +
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      color      = "grey45",
      linewidth  = 0.4
    ) +
    geom_text_repel(
      aes(label = MP_label),
      size         = 3.2,
      max.overlaps = 20
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title    = title_text,
      subtitle = sprintf(
        "n_sig (p<0.05): %d / %d MPs  |  n samples = %d",
        sum(df$sig, na.rm = TRUE),
        nrow(df),
        df$n[1]
      ),
      x = "log2(HR)",
      y = "-log10(p)"
    )

  p
}

####################
# 10) Compose and save two-panel PDF
####################
p_all <- make_volcano(cox_all, "EAC — All sample types")
p_01  <- make_volcano(cox_01,  "EAC — Primary tumour only (01)")

out_pdf <- "Auto_test_sample_type_filter_EAC_volcanoes.pdf"

pdf(out_pdf, width = 16, height = 7, useDingbats = FALSE)
print(p_all | p_01)
dev.off()

cat(sprintf("Saved: %s\n", out_pdf))

####################
# 11) Machine-readable summary (updates/new_updates/summaries/ convention)
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "Auto_test_sample_type_filter_summary.csv")

make_summary_row <- function(df, label) {
  data.frame(
    filter_label = label,
    n_mps_tested = nrow(df),
    n_sig_pval   = if (nrow(df) > 0) sum(df$P_value < 0.05, na.rm = TRUE) else 0L,
    n_sig_padj   = if (nrow(df) > 0) sum(df$padj   < 0.05, na.rm = TRUE) else 0L,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(
  make_summary_row(cox_all, "All sample types"),
  make_summary_row(cox_01,  "Primary only (01)")
)
write.csv(summary_df, summary_path, row.names = FALSE)
cat(sprintf("Saved summary: %s\n", summary_path))

cat("=== Auto_test_sample_type_filter complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
