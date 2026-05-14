####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_compare_mp_selection.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_test_mp_weight_threshold.R
# Test whether MP gene-weight threshold affects continuous Cox survival results.
#
# Goals:
#   Compare three gene-selection strategies for GSVA scoring (EAC only):
#     a) Top-N genes per MP, where N = min(lengths of all retained MP gene lists)
#        using ranked weights from geneNMF.metaprograms$metaprograms.genes.weights
#     b) Genes with weight >= min 10th-highest weight across all retained MPs
#        i.e. threshold = min(sapply(weights, `[`, 10))
#     c) All genes per MP (current default — no weight filtering)
#
#   Produces three separate volcano plots (EAC, one per strategy) side by side
#   in a single PDF for visual comparison.
#
# Input:  ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#         /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
# Output: ref_outs/Auto_test_mp_weight_threshold_EAC_volcanoes.pdf
#         updates/new_updates/summaries/Auto_test_mp_weight_threshold_summary.csv
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

cat("=== Auto_test_mp_weight_threshold ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load MP definitions and apply silhouette filter
####################
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
mp.weights <- geneNMF.metaprograms$metaprograms.genes.weights

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_names <- paste0("MP", bad_mps)
  mp.genes   <- mp.genes[!names(mp.genes) %in% bad_names]
  mp.weights <- mp.weights[!names(mp.weights) %in% bad_names]
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
# 2) Derive gene thresholds from weights
#    Weights lists are assumed to be already sorted descending (highest weight first).
#    If not, sort them here.
####################
mp.weights <- lapply(mp.weights, function(w) sort(w, decreasing = TRUE))

# a) Min-length threshold: N = minimum number of genes across all retained MPs
n_min <- min(sapply(mp.weights, length))
cat(sprintf("Strategy a: top-%d genes per MP (min MP length)\n", n_min))

# b) Weight-value threshold: minimum of the 10th-highest weight across all retained MPs
#    (AGENTS.md note: the minimum MP has 20 genes, so index 10 is always valid)
w10_per_mp <- sapply(mp.weights, function(w) w[10])
w_threshold <- min(w10_per_mp, na.rm = TRUE)
cat(sprintf("Strategy b: weight threshold = %.6f (min 10th weight across MPs)\n", w_threshold))

# c) All genes — no filtering (current default)
cat("Strategy c: all genes per MP (no weight filter)\n")

####################
# Helper: build gene set from weights using a selection strategy
####################
build_gene_sets <- function(strategy) {
  lapply(mp.weights, function(w) {
    if (strategy == "a") {
      names(w)[seq_len(min(n_min, length(w)))]
    } else if (strategy == "b") {
      names(w)[w >= w_threshold]
    } else {
      # c: all genes (use original mp.genes list to stay consistent)
      names(w)
    }
  })
}

strategies <- c("a", "b", "c")
strategy_labels <- c(
  "a" = sprintf("a) Top-%d genes (min-length)", n_min),
  "b" = sprintf("b) Weight >= %.4f (min 10th wt)", w_threshold),
  "c" = "c) All genes (no filter)"
)

####################
# 3) Load TCGA data
####################
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df    <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat   <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

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
# 5) Cox helper (EAC only)
####################
run_cox_eac <- function(gsva_sets, meta_tcga, tpm_mat, strategy_label) {
  # Filter gene sets to those with >= 5 overlapping genes in TPM
  gsva_sets <- lapply(gsva_sets, unique)
  gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
  gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]

  if (length(gsva_sets) == 0) {
    cat(sprintf("  [%s] No valid GSVA sets — skipping.\n", strategy_label))
    return(data.frame())
  }

  cat(sprintf("  [%s] Running GSVA on %d sets...\n", strategy_label, length(gsva_sets)))
  gene_set_sizes <- sapply(gsva_sets, length)
  cat(sprintf("  Gene set sizes (min/median/max): %d / %.0f / %d\n",
              min(gene_set_sizes), median(gene_set_sizes), max(gene_set_sizes)))

  gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
  gsva_df <- as.data.frame(t(gsva_scores))
  gsva_df$sample_barcode <- rownames(gsva_df)

  surv_data <- meta_tcga %>%
    inner_join(gsva_df, by = "sample_barcode") %>%
    filter(sample_type_code == "01")

  surv_data$HistologyGroup <- infer_histology(surv_data$type)
  eac_data <- surv_data %>% filter(HistologyGroup == "EAC")

  cat(sprintf("  [%s] EAC n = %d\n", strategy_label, nrow(eac_data)))

  mp_cols <- intersect(names(gsva_sets), colnames(eac_data))
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
      strategy    = strategy_label,
      MP          = mp,
      beta        = ss$coefficients[1, "coef"],
      HR          = ss$coefficients[1, "exp(coef)"],
      P_value     = ss$coefficients[1, "Pr(>|z|)"],
      n           = fit$n,
      events      = fit$nevent,
      n_genes_gsva = gene_set_sizes[mp],
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) return(data.frame())
  res      <- bind_rows(out)
  res$padj <- p.adjust(res$P_value, method = "BH")
  res
}

####################
# 6) Run all three strategies
####################
all_results <- list()
for (strat in strategies) {
  gene_sets <- build_gene_sets(strat)
  all_results[[strat]] <- run_cox_eac(gene_sets, meta_tcga, tpm_mat, strategy_labels[strat])
}

# Add MP display labels
for (strat in strategies) {
  df <- all_results[[strat]]
  if (nrow(df) == 0) next
  df$MP_label <- ifelse(
    !is.na(mp_descriptions[df$MP]),
    unname(mp_descriptions[df$MP]),
    df$MP
  )
  all_results[[strat]] <- df
}

####################
# 7) Volcano plot function (returns ggplot object)
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

  # Label all points (same style as source script: geom_text_repel on all df)
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
        "n_sig (p<0.05): %d / %d MPs",
        sum(df$sig, na.rm = TRUE),
        nrow(df)
      ),
      x = "log2(HR)",
      y = "-log10(p)"
    )

  p
}

####################
# 8) Compose and save three-panel PDF
####################
plots <- lapply(strategies, function(strat) {
  make_volcano(all_results[[strat]], strategy_labels[strat])
})

out_pdf <- "Auto_test_mp_weight_threshold_EAC_volcanoes.pdf"

pdf(out_pdf, width = 21, height = 7, useDingbats = FALSE)
print(plots[[1]] | plots[[2]] | plots[[3]])
dev.off()

cat(sprintf("Saved: %s\n", out_pdf))

####################
# 9) Machine-readable summary (updates/new_updates/summaries/ convention)
####################
summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "Auto_test_mp_weight_threshold_summary.csv")

summary_rows <- lapply(strategies, function(strat) {
  df <- all_results[[strat]]
  data.frame(
    strategy      = strategy_labels[strat],
    n_mps_tested  = nrow(df),
    n_sig_pval    = if (nrow(df) > 0) sum(df$P_value < 0.05, na.rm = TRUE) else 0L,
    n_sig_padj    = if (nrow(df) > 0) sum(df$padj   < 0.05, na.rm = TRUE) else 0L,
    min_genes     = if (nrow(df) > 0) min(df$n_genes_gsva,  na.rm = TRUE) else NA_integer_,
    max_genes     = if (nrow(df) > 0) max(df$n_genes_gsva,  na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
})
summary_df <- bind_rows(summary_rows)
write.csv(summary_df, summary_path, row.names = FALSE)
cat(sprintf("Saved summary: %s\n", summary_path))

cat("=== Auto_test_mp_weight_threshold complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
