####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/clinical/legacy_survival_cibersort.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/cibersort_result.R
# Reorganized as part of analysis/ restructuring
####################
library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(survival)
library(survminer)

# ── 1. Load data ──────────────────────────────────────────────────────────────
meta2   <- readRDS("tcga_esca_meta.rds")
tpm_df  <- data.table::fread("TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

# ── 2. Build gene-set list: one entry per MP ─────────────────────────────────
# mp.genes is assumed to be a named list, e.g.
#   mp.genes$MP1, mp.genes$MP2, ..., mp.genes$MP9
# Keep every MP as its own separate gene set
gsva_sets <- lapply(mp.genes, function(genes) unique(genes))
cat("Gene sets to score:\n")
print(names(gsva_sets))
cat("\nGenes per set:\n")
print(sapply(gsva_sets, length))

# ── 3. Filter gene sets to genes actually present in the TPM matrix ──────────
gsva_sets_filtered <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
cat("\nGenes matched in TCGA TPM matrix per set:\n")
print(sapply(gsva_sets_filtered, length))

# Drop any MP with too few genes (< 5) to score reliably
gsva_sets_filtered <- gsva_sets_filtered[sapply(gsva_sets_filtered, length) >= 5]

# ── 4. Run GSVA ──────────────────────────────────────────────────────────────
gsva_scores <- gsva(
  tpm_mat,
  gsva_sets_filtered,
  method = "gsva",
  kcdf   = "Gaussian"
)

# ── 5. Build analysis data frame (primary tumours only) ──────────────────────
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

surv_data <- meta2 %>%
  inner_join(gsva_df, by = "sample_barcode") %>%
  filter(sample_type_code == "01")        # keep primary tumours only

mp_names <- names(gsva_sets_filtered)

# ── 6. Per-sample "dominant MP" assignment ───────────────────────────────────
#    For every sample, which MP has the highest GSVA enrichment score?
mp_score_mat <- surv_data[, mp_names]
surv_data$dominant_MP <- mp_names[apply(mp_score_mat, 1, which.max)]
surv_data$max_gsva    <- apply(mp_score_mat, 1, max)

# ── 7. Per-sample "validated" MPs (score > 0 = enriched above background) ───
#    A more lenient view: which MPs are "active" (GSVA > 0) in each sample?
validated_long <- surv_data %>%
  select(sample_barcode, all_of(mp_names)) %>%
  pivot_longer(cols = all_of(mp_names),
               names_to  = "MP",
               values_to = "gsva_score") %>%
  mutate(is_enriched = gsva_score > 0)     # GSVA > 0 ⟹ enriched

# ── 8. PLOT A — Bar plot: how many samples have each MP as DOMINANT ──────────
dominant_counts <- surv_data %>%
  count(dominant_MP) %>%
  mutate(dominant_MP = factor(dominant_MP, levels = dominant_MP[order(-n)]))

p1 <- ggplot(dominant_counts, aes(x = dominant_MP, y = n, fill = dominant_MP)) +
  geom_col(colour = "black", show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
  labs(
    title    = "Dominant MP per TCGA-ESCA primary tumour",
    subtitle = "Each sample assigned to its highest-scoring MP",
    x = "Meta-Program", y = "Number of samples"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# ── 9. PLOT B — Bar plot: how often is each MP enriched (GSVA > 0)? ─────────
enrichment_summary <- validated_long %>%
  group_by(MP) %>%
  summarise(
    n_enriched   = sum(is_enriched),
    n_total      = n(),
    pct_enriched = 100 * n_enriched / n_total,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_enriched)) %>%
  mutate(MP = factor(MP, levels = MP))
enrichment_summary$MP <- factor(enrichment_summary$MP, levels = mp_tree_order)

library(ggplot2)

# ---- Harmonized MP palette (matched to your existing annotation style) ----
mp_cols <- c(
  "MP1" = "#377EB8",   # distinct blue
  "MP2" = "#66C2A5",   # teal
  "MP3" = "#8E63B0",   # metaplasia-like purple
  "MP4" = "#A6CEE3",   # soft blue
  "MP5" = "#D73027",   # prolif-like red
  "MP6" = "#B07CC6",   # second metaplasia purple
  "MP7" = "#1F78B4",   # strong blue (distinct)
  "MP8" = "#F28E2B",   # plastic-like orange
  "MP9" = "#4DAF4A"    # squamous-like green
)

# ---- Plot ----
p2 <- ggplot(enrichment_summary,
             aes(x = MP, y = pct_enriched, fill = MP)) +
  
  geom_col(width = 0.72,
           colour = "black",
           linewidth = 0.4,
           show.legend = FALSE) +
  
  geom_text(aes(label = paste0(round(pct_enriched, 1),
                               "%\n(n=", n_enriched, ")")),
            vjust = -0.6,
            size = 5) +
  
  scale_fill_manual(values = mp_cols) +
  
  labs(
    title = "MP Validation in TCGA-ESCA",
    y = "Samples Enriched (%)",
    x = NULL
  ) +
  
  coord_cartesian(ylim = c(0, 105)) +
  
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold",
                              size = 20,
                              hjust = 0.5),
    
    axis.title.y = element_text(face = "bold",
                                size = 16),
    
    axis.text.x = element_text(face = "bold",
                               size = 12),
    
    axis.text.y = element_text(size = 14),
    
    axis.line = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.7),
    
    panel.grid = element_blank()
  )

print(p2)



# ── 10. PLOT C — Heatmap of GSVA scores across all samples ──────────────────
library(pheatmap)

annotation_col <- data.frame(
  Dominant_MP = surv_data$dominant_MP,
  row.names   = surv_data$sample_barcode
)

score_matrix <- as.matrix(surv_data[, mp_names])
rownames(score_matrix) <- surv_data$sample_barcode

# Order samples by dominant MP for cleaner visualisation
sample_order <- order(surv_data$dominant_MP)

pheatmap(
  t(score_matrix[sample_order, ]),
  cluster_cols         = FALSE,
  cluster_rows         = TRUE,
  annotation_col       = annotation_col,
  show_colnames        = FALSE,
  color                = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main                 = "GSVA enrichment scores — all MPs across TCGA-ESCA tumours",
  fontsize_row         = 10
)

# ── 11. PLOT D — Violin + jitter of GSVA score distributions per MP ─────────
p4 <- ggplot(validated_long,
             aes(x = MP, y = gsva_score, fill = MP)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    title = "GSVA score distribution per MP in TCGA-ESCA",
    x = "Meta-Program", y = "GSVA enrichment score"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)

# ── 12. Summary table ───────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("SUMMARY: MP validation in TCGA-ESCA\n")
cat("══════════════════════════════════════════════════════\n\n")

summary_table <- enrichment_summary %>%
  left_join(dominant_counts, by = c("MP" = "dominant_MP")) %>%
  rename(n_dominant = n) %>%
  mutate(n_dominant = ifelse(is.na(n_dominant), 0, n_dominant)) %>%
  select(MP, n_enriched, pct_enriched, n_dominant, n_total)

print(as.data.frame(summary_table))

#############################

library(GSVA)
meta2 <- readRDS("tcga_esca_meta.rds")

tpm_df <- data.table::fread("TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt") # Or your specific file path
tpm_mat <- as.matrix(tpm_df[, -1]) # Remove GeneSymbol col to make matrix
rownames(tpm_mat) <- tpm_df$GeneSymbol

Idents(clean_seurat) <- "manual_state"
all_markers <- FindAllMarkers(
  clean_seurat, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.5,
  test.use = "wilcox"
)

dge_gsva_list <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  split(f = .$cluster) %>%
  lapply(function(x) x$gene)

gsva_scores <- gsva(tpm_mat, dge_gsva_list, method = "gsva", kcdf = "Gaussian")

gsva_sets <- list(
  "Intestinal Metaplasia" = unique(c(mp.genes$MP3, mp.genes$MP6)),
  "Basal to Intest. Meta" = mp.genes$MP9,
  "Classic Proliferative" = mp.genes$MP5,
  "Stress-adaptive"       = mp.genes$MP8
)

gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")

gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

surv_data <- meta2 %>%
  inner_join(gsva_df, by = "sample_barcode")
surv_data <- surv_data[surv_data$sample_type_code == "01", ]

##############################
# 1. Load Essential Libraries (Order matters)
library(survival)
library(dplyr)   # Fixes 'could not find function "rename"'
library(readr)   # Fixes 'could not find function "read_csv"'
library(ggplot2) # Used for plotting instead of survminer

setwd("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/TCGA/INPUT")
meta2 <- readRDS("tcga_esca_meta.rds")
results <- read_csv("CIBERSORTx_Job5_Results.csv", show_col_types = FALSE) %>%
  dplyr::rename(sample_barcode = Mixture)

surv_data <- meta2 %>%
  inner_join(results, by = "sample_barcode") %>%
  filter(`P-value` <= 0.05) 

surv_data <- surv_data[surv_data$sample_type_code == "01", ]

# ==============================================================================
# 1. SETUP & FRESH DATA LOAD (Prevents Double-Conversion Error)
# ==============================================================================
library(survival)
library(survminer)
library(dplyr)

cell_states <- c("Stress-adaptive", "Intestinal Metaplasia", "Basal to Intest. Meta", "Classic Proliferative")

pdf("CIBERSORT.pdf", width = 8, height = 10)

methods <- c("median", "quartile", "optimal")

for (m in methods) {
  
  for (state in cell_states) {
    
    # drop missing survival/state for safety
    df <- surv_data %>%
      filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[state]]))
    
    if (nrow(df) < 10) next
    
    if (m == "median") {
      
      median_val <- median(df[[state]], na.rm = TRUE)
      df$Group <- ifelse(df[[state]] >= median_val, "High", "Low")
      df$Group <- factor(df$Group, levels = c("Low", "High"))
      
      fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = df)
      
      p <- ggsurvplot(
        fit, data = df,
        palette = c("#E7B800", "#2E9FDF"),
        title = paste("Median split:", state),
        subtitle = paste("Median cutoff =", round(median_val, 2)),
        xlab = "Time (Months)", xscale = "d_m",
        break.time.by = 365.25, xlim = c(0, 1826),
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        legend.labs = c("Low (< Median)", "High (≥ Median)"),
        ggtheme = theme_bw()
      )
      
      print(p)
      
    } else if (m == "optimal") {
      
      # survminer needs complete rows
      df2 <- df %>% select(OS_time, OS_event, all_of(state))
      res.cut <- surv_cutpoint(df2, time = "OS_time", event = "OS_event", variables = state)
      df_cat <- surv_categorize(res.cut)
      df_cat$Group <- factor(df_cat[[state]], levels = c("low", "high"))
      
      fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = df_cat)
      
      p <- ggsurvplot(
        fit, data = df_cat,
        palette = c("#E7B800", "#2E9FDF"),
        title = paste("Optimal cutpoint:", state),
        subtitle = paste("Cut =", round(res.cut$cutpoint$cutpoint, 2)),
        xlab = "Time (Months)", xscale = "d_m",
        break.time.by = 365.25, xlim = c(0, 1826),
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        legend.labs = c("Low", "High"),
        ggtheme = theme_bw()
      )
      
      print(p)
      
    } else if (m == "quartile") {
      
      q_low  <- quantile(df[[state]], 0.25, na.rm = TRUE)
      q_high <- quantile(df[[state]], 0.75, na.rm = TRUE)
      
      df$Group <- NA_character_
      df$Group[df[[state]] <= q_low]  <- "Low Q1"
      df$Group[df[[state]] >= q_high] <- "High Q4"
      df <- df[!is.na(df$Group), ]
      df$Group <- factor(df$Group, levels = c("Low Q1", "High Q4"))
      
      # skip if too small after dropping middle
      if (min(table(df$Group)) < 5) next
      
      fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = df)
      
      p <- ggsurvplot(
        fit, data = df,
        palette = c("#E7B800", "#CC0000"),
        title = paste("Quartile extremes:", state),
        subtitle = paste0("Q1 ≤ ", round(q_low, 2), " vs Q4 ≥ ", round(q_high, 2)),
        xlab = "Time (Months)", xscale = "d_m",
        break.time.by = 365.25, xlim = c(0, 1826),
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        legend.labs = c("Low Q1", "High Q4"),
        ggtheme = theme_bw()
      )
      
      print(p)
    }
  }
  
  # separator page BETWEEN methods only (prevents white first page)
  if (m != tail(methods, 1)) {
    grid::grid.newpage()
  }
}

dev.off()




state = "Squamous_Transition"
q <- quantile(surv_data[[state]], probs = c(0.25, 0.75), na.rm = TRUE)
surv_data$Group <- NA
surv_data$Group[surv_data[[state]] <= q[1]] <- "Low"
surv_data$Group[surv_data[[state]] >= q[2]] <- "High"
surv_data$Group <- factor(surv_data$Group, levels = c("Low", "High"))
fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = surv_data)
median_summary <- summary(fit)$table
median_df <- data.frame(
  Group = rownames(median_summary),
  N = median_summary[, "records"],
  Events = median_summary[, "events"],
  Median_Survival_Days = median_summary[, "median"],
  Lower_95CI = median_summary[, "0.95LCL"],
  Upper_95CI = median_summary[, "0.95UCL"]
)
print(median_df)



######################################################
library(survival)
library(ggplot2)

# ---- settings ----
gender_col <- "Stage_Simple"  # adjust if your column name differs

valid_data <- surv_data
valid_data[[gender_col]] <- as.character(valid_data[[gender_col]])
valid_data[[gender_col]][is.na(valid_data[[gender_col]])] <- "Unknown"
valid_data[[gender_col]] <- factor(
  valid_data[[gender_col]],
  levels = c(levels(surv_data[[gender_col]]), "Unknown")
)
time_col  <- "OS_time"
event_col <- "OS_event"
form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", gender_col))
fit <- survfit(form, data = valid_data)

# Log-rank test
sd <- survdiff(form, data = valid_data)
p_val <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

# Build plot df (no survminer)
plot_df <- data.frame(
  time  = fit$time,
  surv  = fit$surv,
  strata = rep(names(fit$strata), fit$strata)
)
zeros <- data.frame(time = 0, surv = 1, strata = names(fit$strata))
plot_df <- rbind(zeros, plot_df)

# Nice labels: turn "Gender=Male" -> "Male"
plot_df$strata_label <- sub(paste0("^", gender_col, "="), "", plot_df$strata)

# N per gender for annotation
n_by_gender <- table(valid_data[[gender_col]])
n_text <- paste(names(n_by_gender), as.integer(n_by_gender), sep = "=", collapse = "  |  ")

p <- ggplot(plot_df, aes(x = time, y = surv, color = strata_label)) +
  geom_step(linewidth = 1) +
  theme_bw() +
  labs(
    title = paste(gender_col),
    subtitle = paste("Log-rank P =", format(p_val, digits = 4)),
    x = "Time (Days)",
    y = "Survival Probability",
    color = gender_col
  ) +
  ylim(0, 1) +
  annotate("text",
           x = max(plot_df$time, na.rm = TRUE) * 0.05,
           y = 0.1,
           label = paste0("N=", nrow(valid_data), "  (", n_text, ")"),
           hjust = 0)

print(p)


##########################

library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Prepare data with a safer separator
summary_df <- surv_data %>%
  select(all_of(c(gender_col, cell_states))) %>%
  group_by(.data[[gender_col]]) %>%
  summarise(
    across(
      all_of(cell_states),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        se   = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}.{.fn}"
    ),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(Gender = .data[[gender_col]]) %>%   # optional: standardize the group column name
  pivot_longer(
    cols = matches("\\.(mean|se)$"),
    names_to = c("State", ".value"),
    names_sep = "\\."
  )
x_labs <- summary_df %>%
  distinct(Gender, n) %>%
  arrange(Gender) %>%
  { setNames(paste0(.$Gender, "\n(n=", .$n, ")"), .$Gender) }

ggplot(summary_df, aes(x = Gender, y = mean, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  scale_x_discrete(labels = x_labs) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  labs(
    x = gender_col,      # keep your original label text if you want
    y = "Proportion",
    fill = "Cell State"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 11, face = "bold")
  )

################################

library(dplyr)
library(survival)
library(ggplot2)

# ----------------------------
# USER INPUTS
# ----------------------------
cell_states <- c("Stress-adaptive", "Intestinal Metaplasia", "Basal to Intest. Meta", "Classic Proliferative")
state_of_interest <- "Classic Proliferative"   # <- change to any column in surv_data

time_col  <- "OS_time"
event_col <- "OS_event"
stage_col <- "Stage_Simple"
grade_col <- "Grade"

# KM extremes settings (for visualization)
km_probs <- c(0.25, 0.75)   # extremes split; change to c(0.25, 0.75) if desired
km_min_per_group <- 0
km_min_events <- 0

# ----------------------------
# HELPERS
# ----------------------------
zscore_cols <- function(df, cols) {
  df %>% mutate(across(all_of(cols), ~ as.numeric(scale(as.numeric(.))), .names = "{.col}_z"))
}

tidy_cox <- function(fit, model_name) {
  s <- summary(fit)
  out <- data.frame(
    model = model_name,
    term = rownames(s$coef),
    HR = exp(s$coef[, "coef"]),
    CI_low = s$conf.int[, "lower .95"],
    CI_high = s$conf.int[, "upper .95"],
    p = s$coef[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  out$n <- fit$n
  out$events <- fit$nevent
  out
}

plot_forest <- function(tab, title_text) {
  tab2 <- tab %>%
    mutate(term = gsub("`", "", term)) %>%
    mutate(term = factor(term, levels = rev(unique(term))))
  
  ggplot(tab2, aes(x = term, y = HR)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_y_log10() +
    coord_flip() +
    theme_bw() +
    labs(title = title_text, x = NULL, y = "Hazard Ratio (log scale)")
}

plot_km_extremes_survminer <- function(df, state,
                                       probs = c(0.25, 0.75),
                                       min_per_group = 5,
                                       palette = c("#E7B800", "#CC0000"), # Low=yellow, High=red
                                       xlim_days = c(0, 1826),
                                       break_days = 365.25) {
  
  # keep complete cases for this state + survival
  df <- df %>%
    dplyr::filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[state]]))
  
  # compute cutoffs within this df
  q <- quantile(df[[state]], probs = probs, na.rm = TRUE, names = FALSE)
  if (!is.finite(q[1]) || !is.finite(q[2]) || q[1] >= q[2]) return(NULL)
  
  # Low = bottom tail, High = top tail (consistent with your quartile code)
  df$Group <- NA_character_
  df$Group[df[[state]] <= q[1]] <- "Low"
  df$Group[df[[state]] >= q[2]] <- "High"
  df <- df[!is.na(df$Group), , drop = FALSE]
  df$Group <- factor(df$Group, levels = c("Low", "High"))
  
  # skip if too small
  if (min(table(df$Group)) < min_per_group) return(NULL)
  
  fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = df)
  
  ggsurvplot(
    fit, data = df,
    palette = palette,                          # Low then High
    title = paste("Extremes split:", state),
    subtitle = paste0(
      round(probs[1] * 100), "% ≤ ", round(q[1], 2),
      " vs ",
      round(probs[2] * 100), "% ≥ ", round(q[2], 2)
    ),
    xlab = "Time (Months)", xscale = "d_m",      # months
    break.time.by = break_days, xlim = xlim_days,
    risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    legend.labs = c("Low", "High"),
    ggtheme = theme_bw()
  )
}


# ----------------------------
# MAIN ANALYSIS
# ----------------------------

# Basic checks
stopifnot(state_of_interest %in% colnames(surv_data))
stopifnot(all(cell_states %in% colnames(surv_data)))

# Ensure Stage/Grade are factors if present
surv_data <- surv_data %>%
  mutate(
    !!stage_col := factor(.data[[stage_col]]),
    !!grade_col := factor(.data[[grade_col]])
  )

# Z-score all states (adds *_z columns)
surv_data <- zscore_cols(surv_data, cell_states)

# Model 1: all states together
all_z <- paste0(cell_states, "_z")
df1 <- surv_data %>%
  dplyr::select(all_of(c(time_col, event_col, all_z))) %>%
  filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]])) %>%
  na.omit()

fit1 <- coxph(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ ",
                                paste(all_z, collapse = " + "))), data = df1)

tab1 <- tidy_cox(fit1, "Cox: all states (z)")

# Model 2: chosen state + stage
state_z <- paste0(state_of_interest, "_z")
df2 <- surv_data %>%
  select(all_of(c(time_col, event_col, state_z, stage_col))) %>%
  filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]])) %>%
  na.omit()

fit2 <- coxph(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ ",
                                state_z, " + ", stage_col)), data = df2)

tab2 <- tidy_cox(fit2, paste0("Cox: ", state_of_interest, " + Stage"))

# Model 3: chosen state + grade
df3 <- surv_data %>%
  select(all_of(c(time_col, event_col, state_z, grade_col))) %>%
  filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]])) %>%
  na.omit()

fit3 <- coxph(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ ",
                                state_z, " + ", grade_col)), data = df3)

tab3 <- tidy_cox(fit3, paste0("Cox: ", state_of_interest, " + Grade"))

# Combined significance table
sig_table <- bind_rows(tab1, tab2, tab3) %>%
  mutate(p_adj_BH = ave(p, model, FUN = function(x) p.adjust(x, method = "BH"))) %>%
  arrange(model, p)

print(sig_table)

# ----------------------------
# VISUALIZE (forest + KM)
# ----------------------------
print(plot_forest(tab1, "Cox model: all states (z-scored)"))
print(plot_forest(tab2, paste0("Cox model: ", state_of_interest, " (z) + Stage")))
print(plot_forest(tab3, paste0("Cox model: ", state_of_interest, " (z) + Grade")))

# KM overall for chosen state (extremes)
p <- plot_km_extremes_survminer(df1, state, probs = c(0.25, 0.75))
if (!is.null(p)) print(p)

###########################

median_val <- median(surv_data[[state]], na.rm = TRUE)
surv_data$Group <- ifelse(surv_data[[state]] >= median_val, "High", "Low")
surv_data$Group <- factor(surv_data$Group, levels = c("Low", "High"))

library(survival)
library(survminer)
library(dplyr)

surv_data_early <- surv_data %>%
  filter(sample_type_code == "01") %>%
  #  filter(is.na(Stage_Simple) | Stage_Simple != "Stage IV") %>%
  filter(!is.na(OS_time), !is.na(OS_event), !is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("Low","High")))

fit_early <- survfit(Surv(OS_time, OS_event) ~ Group, data = surv_data_early)

p <- ggsurvplot(
  fit_early,
  data = surv_data_early,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  title = "Overall Survival by Group (Stage IV excluded; NA stage included)",
  
  # ---- ADD THESE ----
  xlab = "Time (Months)",
  xscale = "d_m",            # convert days -> months
  break.time.by = 365.25,    # tick every year
  xlim = c(0, 1826),         # 5 years
  
  ylab = "Survival Probability",
  legend.title = "Group",
  legend.labs = c("Low", "High"),
  ggtheme = theme_bw()
)

print(p)


###########################

library(ggplot2)
library(rlang)
library(survminer)

# ---- Compute cut values ----
median_cut  <- median(surv_data[[state]], na.rm = TRUE)

# Optimal cut
res.cut <- surv_cutpoint(
  surv_data,
  time = "OS_time",
  event = "OS_event",
  variables = state
)
optimal_cut <- res.cut$cutpoint$cutpoint

# Quartiles
q1_cut <- quantile(surv_data[[state]], 0.25, na.rm = TRUE)
q4_cut <- quantile(surv_data[[state]], 0.75, na.rm = TRUE)

# ---- Build dataframe of cut lines ----
cut_df <- data.frame(
  cut_value = c(median_cut, optimal_cut, q1_cut, q4_cut),
  Cut_Type = factor(
    c("Median", "Optimal", "Q1 (25%)", "Q4 (75%)"),
    levels = c("Median", "Optimal", "Q1 (25%)", "Q4 (75%)")
  )
)

# ---- Density Plot ----
p_dist <- ggplot(surv_data, aes(x = !!sym(state))) +
  geom_density(fill = "grey85", color = "black", linewidth = 1) +
  
  # Add vertical cut lines
  geom_vline(
    data = cut_df,
    aes(xintercept = cut_value, color = Cut_Type),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "Median" = "red",
      "Optimal" = "orange",
      "Q1 (25%)" = "green",
      "Q4 (75%)" = "green"
    )
  ) +
  
  theme_classic() +
  labs(
    title = paste("Distribution of", state),
    subtitle = "Median, Optimal, and Quartile Cutpoints",
    x = "Score / Fraction",
    y = "Density",
    color = "Cut Type"
  )

print(p_dist)
