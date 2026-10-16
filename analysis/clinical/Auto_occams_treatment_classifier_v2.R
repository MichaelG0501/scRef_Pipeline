####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/Auto_occams_treatment_classifier_v2.R
#   Description: Extracts Pre vs Post treatment signatures from a balanced subset of 
#     scATLAS studies (avoiding batch confounding), scores OCCAMS bulk via UCell, 
#     and validates against strict clinical metadata.
#   Map: analysis/ANALYSIS_MAP.md
# Conda env: dmtcp
####################

library(Seurat)
library(edgeR)
library(UCell)
library(ggplot2)
library(dplyr)
library(readxl)

WD <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(WD)
source("analysis/shared/scRef_config.R")

out_root <- "ref_outs/OCCAMS/treatment_classifier_v2"
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "intermediate"), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# STEP 1: Balanced DE
# ============================================================================
sig_path <- file.path(out_root, "intermediate", "Auto_balanced_signatures.rds")

if (!file.exists(sig_path)) {
  message("Loading scATLAS epithelial object...")
  epi <- readRDS(SCREF_EPI_RDS)
  
  xl <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
  xl$orig.ident <- paste(xl$Author, xl$Year, xl[["Sample Name"]], sep = "_")
  tx_lookup <- setNames(xl$Treatment, xl$orig.ident)
  author_lookup <- setNames(xl$Author, xl$orig.ident)
  
  epi$treatment_raw <- unname(tx_lookup[epi$orig.ident])
  epi$author <- unname(author_lookup[epi$orig.ident])
  
  # Bin into Pre/Post
  epi$treatment_bin <- ifelse(epi$treatment_raw == "Post", "Post", "Pre")
  
  # STRICT FILTERING: Only keep studies that have BOTH Pre and Post samples!
  # Based on our previous check: Alcindor, Carroll, Walker, Yates
  valid_studies <- c("Alcindor", "Carroll", "Walker", "Yates")
  valid_cells <- rownames(epi@meta.data)[epi$author %in% valid_studies & !is.na(epi$treatment_bin)]
  epi_sub <- subset(epi, cells = valid_cells)
  
  message("Pseudobulking balanced subset...")
  pb <- AggregateExpression(epi_sub, group.by = "orig.ident", assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA
  colnames(pb) <- gsub("-", "_", colnames(pb))
  
  sample_tx <- setNames(ifelse(tx_lookup[colnames(pb)] == "Post", "Post", "Pre"), colnames(pb))
  sample_study <- setNames(author_lookup[colnames(pb)], colnames(pb))
  
  valid_samps <- !is.na(sample_tx) & !is.na(sample_study)
  pb <- pb[, valid_samps, drop=FALSE]
  
  group <- factor(sample_tx[valid_samps], levels = c("Pre", "Post"))
  study <- factor(sample_study[valid_samps])
  design <- model.matrix(~ study + group)
  
  message("Running edgeR...")
  dge <- DGEList(counts = pb, group = group)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = "groupPost")
  de_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  de_results$gene <- rownames(de_results)
  
  # Filter out Ribo/Mito just in case they dominate the bulk ranks
  de_results <- de_results[!grepl("^RP[SL]|^MT-", de_results$gene), ]
  
  # Take the top 100 genes by P-value for each direction to guarantee a robust signature size
  sig_pre  <- de_results %>% filter(logFC < 0) %>% arrange(PValue) %>% head(100) %>% pull(gene)
  sig_post <- de_results %>% filter(logFC > 0) %>% arrange(PValue) %>% head(100) %>% pull(gene)
  
  signatures <- list(PreTreatment = sig_pre, PostTreatment = sig_post)
  saveRDS(signatures, sig_path)
  
  sig_df <- rbind(
    data.frame(signature="PreTreatment", gene=sig_pre, rank=seq_along(sig_pre)),
    data.frame(signature="PostTreatment", gene=sig_post, rank=seq_along(sig_post))
  )
  write.csv(sig_df, file.path(out_root, "tables", "Auto_balanced_signature_genes.csv"), row.names = FALSE)
  
  rm(epi, epi_sub, pb); gc()
} else {
  message("Loading cached signatures...")
  signatures <- readRDS(sig_path)
}

# ============================================================================
# STEP 2: Score OCCAMS Bulk
# ============================================================================
message("Loading OCCAMS bulk matrix...")
occams_logcpm <- readRDS("ref_outs/OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds")
occams_mat <- as.matrix(occams_logcpm)

# Filter signatures to present genes
signatures$PreTreatment <- intersect(signatures$PreTreatment, rownames(occams_mat))
signatures$PostTreatment <- intersect(signatures$PostTreatment, rownames(occams_mat))

message("Running UCell...")
ucell_scores <- ScoreSignatures_UCell(matrix = occams_mat, features = signatures, ncores = 4)
colnames(ucell_scores) <- sub("_UCell$", "", colnames(ucell_scores))

df <- data.frame(
  subject = rownames(ucell_scores),
  pre_score = ucell_scores[, "PreTreatment"],
  post_score = ucell_scores[, "PostTreatment"],
  stringsAsFactors = FALSE
)

# Because raw UCell scores have different baselines depending on the exact genes,
# comparing them directly via delta (Pre - Post) can be skewed if one set of genes 
# is universally higher expressed. To fix this, we Z-score the signatures across the cohort
# before taking the delta. This centers the cohort.
df$pre_z <- scale(df$pre_score)[,1]
df$post_z <- scale(df$post_score)[,1]
df$delta_z <- df$pre_z - df$post_z
df$predicted_source <- ifelse(df$delta_z > 0, "Likely_Biopsy", "Likely_Surgical")

# ============================================================================
# STEP 3: Validate
# ============================================================================
meta <- read.csv("ref_outs/OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv", stringsAsFactors = FALSE)
meta$strict_naive <- ifelse(
  meta$pathology_neoadjuvant_therapy == "no", 
  "Definite_Naive", "Ambiguous"
)
meta$strict_naive[is.na(meta$strict_naive)] <- "Ambiguous"

df <- merge(df, meta[, c("subject_id", "strict_naive")], by.x = "subject", by.y = "subject_id")

cm <- table(Predicted = df$predicted_source, Ground_Truth = df$strict_naive)

sink(file.path(out_root, "reports", "Auto_v2_treatment_summary.txt"))
cat("=== Balanced Treatment Classifier V2 ===\n\n")
cat("Signatures derived only from studies with both Pre and Post samples (Alcindor, Carroll, Walker, Yates).\n")
cat("Ribo/Mito genes excluded.\n")
cat("UCell scores are Z-scaled before delta calculation to remove baseline bulk expression bias.\n\n")
cat("Confusion Matrix:\n")
print(cm)
sink()

write.csv(df, file.path(out_root, "tables", "Auto_v2_classification.csv"), row.names=FALSE)

# Plots
p1 <- ggplot(df, aes(x = delta_z, fill = strict_naive, color = strict_naive)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  geom_rug(data = subset(df, strict_naive == "Definite_Naive"), sides = "b", linewidth = 1) +
  scale_fill_manual(values = c("Definite_Naive" = "#4CAF50", "Ambiguous" = "#9E9E9E")) +
  scale_color_manual(values = c("Definite_Naive" = "#2E7D32", "Ambiguous" = "#616161")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Scaled Delta (Pre - Post) Distribution", x = "Z-scaled Delta") +
  theme_minimal()

ggsave(file.path(out_root, "figures", "Auto_v2_density.pdf"), p1, width=8, height=5)

p2 <- ggplot(df, aes(x = pre_score, y = post_score, color = strict_naive)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Definite_Naive" = "#4CAF50", "Ambiguous" = "#9E9E9E")) +
  labs(title = "Raw UCell Scores", x = "Pre Score", y = "Post Score") +
  theme_minimal()

ggsave(file.path(out_root, "figures", "Auto_v2_scatter_raw.pdf"), p2, width=6, height=5)
