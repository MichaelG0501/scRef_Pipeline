####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_tcga_survival_clinical_mps_old.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/evaluate_clinical_MPs.R
# Reorganized as part of analysis/ restructuring
####################
library(dplyr)
library(ggplot2)
library(tidyr)
library(survival)
library(survminer)
library(GSVA)
library(Seurat)
library(purrr) # For map_dfr
library(ggrepel)

# Configuration
REMOVE_CELL_CYCLE <- TRUE # Set to FALSE to skip CC regression
# ── 1. Pre vs Post Comparison ────────────────────────────────────────────────

cat("Loading Seurat Object...\n")
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
# tmdata_all <- readRDS("EAC_Ref_epi.rds")
# ucell_scores <- readRDS("UCell_3CA_MPs.rds")
# tmdata_all@meta.data <- cbind(tmdata_all@meta.data, ucell_scores)
# score_cols <- grep("^X3CA_mp|^3CA_mp", colnames(tmdata_all@meta.data), value = TRUE)

cat("Running Pre vs Post correlation...\n")

if(!"Treatment" %in% colnames(tmdata_all@meta.data)) {
  stop("Column 'Treatment' not found in tmdata_all metadata")
}

tmdata_all$Treatment_B <- ifelse(tmdata_all$Treatment == "Post", "Post", "Pre")
tmdata_all$Treatment_B <- factor(tmdata_all$Treatment_B, levels = c("Pre", "Post"))

score_cols <- grep("^X3CA_mp|^3CA_mp", colnames(tmdata_all@meta.data), value = TRUE)
if(length(score_cols) == 0) {
  cat("Warning: No ^X3CA_mp|^3CA_mp score columns found.\n")
}

# Identify studies with BOTH Pre and Post samples
valid_studies <- tmdata_all@meta.data %>%
  filter(!is.na(Treatment_B), !is.na(study)) %>%
  group_by(study) %>%
  summarise(n_treat = n_distinct(Treatment_B), .groups = 'drop') %>%
  filter(n_treat == 2) %>%
  pull(study)

cat(sprintf("Found %d studies with both Pre and Post samples.\n", length(valid_studies)))

# 1.1 Filter data and calculate Sample-level means
mp_sample_means <- tmdata_all@meta.data %>%
  filter(study %in% valid_studies, !is.na(Treatment_B), !is.na(orig.ident)) %>%
  group_by(orig.ident, Treatment_B) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE), .groups = 'drop')

# 1.2 Correlate each MP with Pre/Post using Wilcoxon testing (Global across valid studies)
prepost_res <- lapply(score_cols, function(mp) {
  df_sub <- mp_sample_means %>% filter(!is.na(.data[[mp]]))
  
  if (length(unique(df_sub$Treatment_B)) == 2) {
    test_res <- wilcox.test(df_sub[[mp]] ~ df_sub$Treatment_B)
    
    mean_pre <- mean(df_sub[[mp]][df_sub$Treatment_B == "Pre"], na.rm = TRUE)
    mean_post <- mean(df_sub[[mp]][df_sub$Treatment_B == "Post"], na.rm = TRUE)
    lfc <- log2((mean_post + 1e-6) / (mean_pre + 1e-6))
    
    data.frame(MP = mp, P_value_PrePost = test_res$p.value, LFC = lfc)
  } else {
    data.frame(MP = mp, P_value_PrePost = NA, LFC = NA)
  }
})

prepost_res <- bind_rows(prepost_res)
prepost_res$padj_PrePost <- p.adjust(prepost_res$P_value_PrePost, method = "BH")


# --- Prepare volcano data ---
volc_df <- prepost_res %>%
  filter(!is.na(LFC), !is.na(P_value_PrePost)) %>%
  mutate(
    neglog10FDR = -log10(P_value_PrePost),
    sig = P_value_PrePost < 0.05
  )

# Join with consistency data
volc2 <- volc_df %>%
  left_join(consistency, by = "MP") %>%
  mutate(n_same_dir = ifelse(is.na(n_same_dir), 0, n_same_dir))

# Select top significant points for labeling
top_lab <- volc2 %>% 
  filter(sig == TRUE) %>%
  arrange(P_value_PrePost) %>% 
  slice_head(n = 15)

# --- Create the publication-quality volcano plot ---
ggplot(volc2, aes(x = LFC, y = neglog10FDR)) +
  # Non-significant points (grey)
  geom_point(
    data = filter(volc2, sig == FALSE),
    color = "grey70",
    size = 3,
    alpha = 0.7
  ) +
  # Significant points (red)
  geom_point(
    data = filter(volc2, sig == TRUE),
    color = "red",
    size = 3.5,
    alpha = 0.85
  ) +
  # Labels with repel to avoid overlap
  geom_text_repel(
    data = top_lab,
    aes(label = MP),
    size = 3,
    color = "grey30",
    segment.color = "grey50",
    segment.size = 0.3,
    segment.alpha = 0.6,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    min.segment.length = 0,
    force = 2
  ) +
  # Reference lines (subtle)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # Labels and title
  
  labs(
    title = "Pre (baseline) vs Post Treatment",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # Clean minimal theme matching reference
  theme_minimal(base_size = 14) +
  theme(
    # Clean white background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Subtle grid
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    # Axis styling
    axis.line = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks = element_line(color = "grey40", linewidth = 0.3),
    axis.text = element_text(color = "grey20", size = 11),
    axis.title = element_text(color = "grey10", size = 12, face = "plain"),
    # Title styling
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey10"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    # Remove legend (color is self-explanatory)
    legend.position = "none",
    # Add margins
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # Set axis limits for balance (adjust as needed)
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))

library(dplyr)
library(ggplot2)

mp_plot <- "X3CA_mp_44.Alveolar"

df_bar <- tmdata_all@meta.data %>%
  filter(study %in% valid_studies,
         !is.na(Treatment_B),
         !is.na(.data[[mp_plot]])) %>%
  group_by(study, Treatment_B) %>%
  summarise(
    val = mean(.data[[mp_plot]], na.rm = TRUE),   # mean over CELLS
    n_cells = dplyr::n(),                         # CELLS per bar
    n_samples = n_distinct(orig.ident),           # optional
    .groups = "drop"
  )

ggplot(df_bar, aes(Treatment_B, val, fill = Treatment_B)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("cells=", n_cells)),
            vjust = -0.5, size = 3.3) +
  facet_wrap(~study, scales = "free_y") +
  labs(title = paste0(mp_plot, ": Pre vs Post per study"),
       x = NULL, y = "Mean UCell score") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))

# ── 2. TCGA Survival Enrichment ──────────────────────────────────────────────

cat("\nRunning TCGA Cox Model enrichment...\n")

# 2.1 Load TCGA data
meta2   <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df  <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

# 2.2 Build gene-set list
cat("Reading New_NMFs.csv...\n")
MP_df <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", check.names = FALSE)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))

gsva_sets <- lapply(MP_list, unique)
gsva_sets_filtered <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets_filtered <- gsva_sets_filtered[sapply(gsva_sets_filtered, length) >= 5]

# 2.3 Run GSVA
cat(sprintf("Running GSVA on %d gene sets...\n", length(gsva_sets_filtered)))
gsva_scores <- gsva(
  tpm_mat,
  gsva_sets_filtered,
  method = "gsva",
  kcdf   = "Gaussian"
)
gsva_scores_t <- t(gsva_scores)

# 2.4 Optional: Cell Cycle Regression on GSVA Scores
CC_FIXED <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)

# Standardize names to match the cleaned MP_list names
# The `make.names` above converts starting numbers, so we handle potential prefix mismatch
CC_FIXED <- make.names(sub("^X3CA_mp_", "3CA_mp_", CC_FIXED))
CC_FIXED <- make.names(sub("^3CA_mp_", "X3CA_mp_", CC_FIXED)) # Just in case make.names added X

current_mps <- colnames(gsva_scores_t)

if (REMOVE_CELL_CYCLE) {
  cc_mps <- intersect(CC_FIXED, current_mps)
  non_cc_mps <- setdiff(current_mps, cc_mps)
  
  cat(sprintf("\nCell-cycle MPs found for regression: %d\n", length(cc_mps)))
  
  if (length(cc_mps) == 0) {
    cat("None of the fixed CC MPs present in GSVA output — skipping regression.\n")
    final_gsva <- gsva_scores_t
    final_mps <- current_mps
  } else {
    cat(sprintf("Regressing out CC signal using %d CC MPs; keeping %d MPs.\n",
                length(cc_mps), length(non_cc_mps)))
    
    X_cc    <- as.matrix(gsva_scores_t[, cc_mps, drop = FALSE])
    Y_other <- as.matrix(gsva_scores_t[, non_cc_mps, drop = FALSE])
    
    # Efficient OLS residualization
    X <- cbind(Intercept = 1, X_cc)
    XtX_inv <- solve(crossprod(X))
    B <- XtX_inv %*% crossprod(X, Y_other)
    Y_hat <- X %*% B
    Y_resid <- Y_other - Y_hat
    
    final_gsva <- Y_resid
    final_mps <- non_cc_mps
  }
} else {
  cat("\nSkipping cell-cycle removal (REMOVE_CELL_CYCLE = FALSE).\n")
  final_gsva <- gsva_scores_t
  final_mps <- current_mps
}

####################

# 2.5 Prepare Survival Data
gsva_df <- as.data.frame(final_gsva)
gsva_df$sample_barcode <- rownames(gsva_df)

surv_data_full <- meta2 %>%
  inner_join(gsva_df, by = "sample_barcode")

surv_data_full <- surv_data_full[surv_data_full$sample_type_code == "01", ]

# 2.6 Run Cox model for each MP (Continuous)
cat("\nFitting continuous Cox models...\n")
cox_results <- list()

for (mp_name in final_mps) {
  # Ensure column name exists and is valid
  if(!mp_name %in% colnames(surv_data_full)) next
  
  df <- surv_data_full %>% 
    filter(!is.na(OS_time) & !is.na(OS_event) & !is.na(.data[[mp_name]]))
  
  if(nrow(df) < 20) next 
  
  score_val <- df[[mp_name]]
  
  # Check for continuous variance instead of categorical proportions
  if(var(score_val, na.rm = TRUE) == 0) {
    cat(sprintf("Skipping %s: No variance in continuous score\n", mp_name))
    next
  }
  
  # Fit Cox model using continuous form
  fit_cont <- as.formula(paste("Surv(OS_time, OS_event) ~ `", mp_name, "`", sep=""))
  
  # Try-catch to handle convergence issues
  cox_fit <- try(coxph(fit_cont, data = df), silent = TRUE)
  
  if(inherits(cox_fit, "try-error")) next
  
  p_val <- summary(cox_fit)$coefficients[1, 5]
  hr <- summary(cox_fit)$coefficients[1, 2]
  
  cox_results[[mp_name]] <- data.frame(
    MP = mp_name,
    HR = hr,
    P_value_Survival = p_val
  )
}

cox_res_df <- bind_rows(cox_results)
if(nrow(cox_res_df) > 0) {
  cox_res_df$padj_Survival <- p.adjust(cox_res_df$P_value_Survival, method = "BH")
}

# ── 3. Combine and Visualise Clinically Significant MPs ──────────────────────

cat("\nGenerating final outputs...\n")

# Plotting the most significant TCGA survival results
if (nrow(cox_res_df) > 0) {
  # Require ggrepel for non-overlapping labels
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
  
  sig_cox <- cox_res_df %>% filter(P_value_Survival < 0.05) %>% arrange(P_value_Survival)
  
  p_volc <- ggplot(cox_res_df, aes(x = log2(HR), y = -log10(P_value_Survival))) +
    geom_point(aes(color = P_value_Survival < 0.05), size = 3, alpha = 0.7) +
    theme_minimal(base_size = 14) +
    labs(title = "Overall survival association (Continous Cox)", 
         x = "Log2 Hazard Ratio", y = "-Log10(p-value)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = c("gray", "red"), guide = "none") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  if(nrow(sig_cox) > 0) {
    # geom_label_repel prevents overlap and keeps labels inside panel
    p_volc <- p_volc + geom_label_repel(
      data = head(sig_cox, 10),
      aes(label = MP),
      size = 3.5,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = Inf,
      fill = alpha("white", 0.8),
      label.size = NA
    )
  }
  p_volc
}

cat("Script complete. Results saved to Auto_AG.\n")

library(survival)
library(survminer)
library(dplyr)

mp_km <- "X3CA_mp_44.Alveolar"  # change if needed

df_km <- surv_data_full %>%
  filter(sample_type_code == "01",
         !is.na(OS_time), !is.na(OS_event), !is.na(.data[[mp_km]])) %>%
  mutate(group = ifelse(.data[[mp_km]] >= median(.data[[mp_km]], na.rm = TRUE),
                        "High", "Low"))

fit_km <- survfit(Surv(OS_time, OS_event) ~ group, data = df_km)

p_km <- ggsurvplot(
  fit_km, data = df_km,
  pval = TRUE, risk.table = TRUE,
  conf.int = FALSE,
  legend.title = mp_km, legend.labs = c("High", "Low"),
  title = paste0("Kaplan–Meier: ", mp_km, " (High vs Low, median split)"),
  xlab = "Days", ylab = "Overall survival probability"
)

print(p_km)

####################

selected_tcga <- cox_res_df %>% filter(P_value_Survival < 0.05) %>% pull(MP)
selected_prepost <- prepost_res %>% filter(P_value_PrePost < 0.05) %>% pull(MP)
selected <- union(selected_tcga, selected_prepost)
