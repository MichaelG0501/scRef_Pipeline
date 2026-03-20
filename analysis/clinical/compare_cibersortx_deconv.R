################################################################################
# Auto_cibersortx_mp_comparison.R
# Compare MPs and retained pan-cancer MPs (3CA) in deconvoluted vs bulk TCGA ESCA
# Order by MP tree/state definition order.
################################################################################
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)

# Setup paths
PWD <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(PWD)

DECON_PATH <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt"
MIXTURE_PATH <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
MP19_PATH <- "ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
PAN_CANCER_MP_PATH <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
RETAINED_3CA_CSV <- "ref_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv"

OUT_DIR <- "ref_outs/cibersortx/analysis/mp_comparison_refined"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
cat("Loading CIBERSORTx results...\n")
decon <- read.delim(DECON_PATH, header = TRUE, row.names = 1, check.names = FALSE)
decon[is.na(decon)] <- 0

cat("Loading mixture...\n")
mixture <- read.delim(MIXTURE_PATH, header = TRUE, row.names = 1, check.names = FALSE)

# Sync samples
common_samples <- intersect(colnames(decon), colnames(mixture))
decon <- decon[, common_samples]
mixture <- mixture[, common_samples]
cat(sprintf("Samples matched: %d\n", length(common_samples)))

# ------------------------------------------------------------------------------
# 2. Get Gene Sets (Retained 3CA and nMP=19)
# ------------------------------------------------------------------------------
cat("Loading MP gene sets & descriptions...\n")

# A. Pipeline MPs (nMP=19)
mp19_obj <- readRDS(MP19_PATH)
mp19_genes_raw <- mp19_obj$metaprograms.genes
# Filter silhouette
bad_mps <- which(mp19_obj$metaprograms.metrics$silhouette < 0)
if(length(bad_mps) > 0) {
  mp19_genes_raw <- mp19_genes_raw[!names(mp19_genes_raw) %in% paste0("MP", bad_mps)]
}

# Standard descriptions
pipe_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

# Tree order from survival script
ordered_mp_list <- c(
  "MP2",
  "MP17", "MP14", "MP5", "MP10", "MP8",
  "MP13", "MP12",
  "MP18", "MP16",
  "MP15",
  "MP1", "MP7", "MP9"
)
extra_mps <- setdiff(names(mp19_genes_raw), ordered_mp_list)
retained_mps_order <- c(ordered_mp_list, extra_mps)
retained_mps_order <- retained_mps_order[retained_mps_order %in% names(mp19_genes_raw)]

# Map to descriptions and store in requested order
mp19_genes <- list()
for(mp in retained_mps_order) {
  desc <- pipe_descriptions[mp]
  label <- if(is.na(desc)) mp else paste(mp, desc)
  mp19_genes[[label]] <- mp19_genes_raw[[mp]]
}

# B. Pan-Cancer MPs (3CA)
ca3_df <- read.csv(PAN_CANCER_MP_PATH, check.names = FALSE)
ca3_full_names <- colnames(ca3_df)
ca3_list <- as.list(ca3_df)
ca3_list <- lapply(ca3_list, function(x) x[x != "" & !is.na(x)])
names(ca3_list) <- ca3_full_names

# Identify retained 3CA ones from Task 4 CSV
cov_df <- read.csv(RETAINED_3CA_CSV)
retained_3ca_ids <- cov_df %>%
  filter(n_samples >= 50, n_studies >= 6, pct_cells >= 1) %>%
  pull(mp_label)

match_mp_num <- function(id) {
  num <- regmatches(id, regexpr("mp_[0-9]+", id))
  if (length(num) == 0) return(NA)
  as.numeric(gsub("mp_", "", num))
}
retained_nums <- sapply(retained_3ca_ids, match_mp_num)
retained_nums <- retained_nums[!is.na(retained_nums)]

retained_full_names <- ca3_full_names[sapply(ca3_full_names, function(h) {
  num_header <- regmatches(h, regexpr("MP[0-9]+", h))
  if (length(num_header) > 0) {
    n <- as.numeric(gsub("MP", "", num_header))
    return(n %in% retained_nums)
  }
  FALSE
})]

retained_sets_3ca <- ca3_list[retained_full_names]
cat(sprintf("Retained 3CA MPs: %d\n", length(retained_sets_3ca)))

# Combine all sets
all_gene_sets <- c(mp19_genes, retained_sets_3ca)
# Intersect with all available genes
available_genes <- rownames(decon)
all_gene_sets_filtered <- lapply(all_gene_sets, function(gs) intersect(gs, available_genes))
retained_sets <- all_gene_sets_filtered[sapply(all_gene_sets_filtered, length) >= 5]
cat(sprintf("Final scoreable MPs: %d\n", length(retained_sets)))

# ------------------------------------------------------------------------------
# 3. GSVA Scoring
# ------------------------------------------------------------------------------
cat("Computing GSVA scores (Bulk and Deconvoluted)...\n")
mixture_log <- as.matrix(log2(mixture + 1))
decon_log <- as.matrix(log2(decon + 1))

gsva_bulk <- gsva(mixture_log, retained_sets, method="gsva", kcdf="Gaussian", verbose=FALSE)
gsva_decon <- gsva(decon_log, retained_sets, method="gsva", kcdf="Gaussian", verbose=FALSE)

# Intersection to ensure common rows
common_mps <- intersect(rownames(gsva_bulk), rownames(gsva_decon))
# Enforce the defined order (Pipeline in its tree order, then 3CA)
ordered_common <- intersect(names(all_gene_sets), common_mps)

gsva_bulk <- gsva_bulk[ordered_common, ]
gsva_decon <- gsva_decon[ordered_common, ]

# ------------------------------------------------------------------------------
# 4. Heatmap Comparison
# ------------------------------------------------------------------------------
cat("Creating combined heatmap...\n")
bulk_z <- t(apply(gsva_bulk, 1, scale))
decon_z <- t(apply(gsva_decon, 1, scale))
colnames(bulk_z) <- colnames(gsva_bulk)
colnames(decon_z) <- colnames(gsva_decon)

# Common Sample Order (cluster samples by decon)
dist_decon <- dist(t(decon_z))
hc_samples <- hclust(dist_decon, method = "ward.D2")
sample_order <- hc_samples$order

col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Row split by source
row_split <- factor(ifelse(ordered_common %in% names(mp19_genes), "Our", "3CA"), levels=c("Our", "3CA"))

ht_bulk <- Heatmap(bulk_z[, sample_order, drop=FALSE], name = "Bulk-Z", col = col_fun,
                   column_order = seq_along(sample_order),
                   cluster_rows = FALSE, # Manual Order
                   row_split = row_split,
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 14),
                   column_title = "Raw TCGA Bulk Mixture Activity")

ht_decon <- Heatmap(decon_z[, sample_order, drop=FALSE], name = "Decon-Z", col = col_fun,
                    column_order = seq_along(sample_order),
                    cluster_rows = FALSE,
                    row_split = row_split,
                    show_column_names = FALSE,
                    row_names_gp = gpar(fontsize = 14),
                    column_title = "Deconvoluted Malignant-only Activity")

pdf(file.path(OUT_DIR, "mp_activity_comparison_heatmap_ordered.pdf"), width = 16, height = 10)
draw(
  ht_bulk + ht_decon,
  heatmap_legend_side = "right",
  merge_legend = TRUE,
  ht_gap = unit(12, "mm")
)
dev.off()

# ------------------------------------------------------------------------------
# 5. Correlation & Summary
# ------------------------------------------------------------------------------
cor_results <- data.frame(MP = ordered_common, stringsAsFactors = FALSE)
cor_results$cor <- sapply(seq_along(ordered_common), function(i) {
  cor(gsva_bulk[i, ], gsva_decon[i, ], method="pearson", use="complete.obs")
})
write.csv(cor_results, file=file.path(OUT_DIR, "mp_correlations_ordered.csv"), row.names = FALSE)

cat("Analysis complete. Results in:", OUT_DIR, "\n")
