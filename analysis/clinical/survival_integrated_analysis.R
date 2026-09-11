####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/survival_integrated_analysis.R
#   Description: Integrated volcano plot of unadjusted survival association 
#     across TCGA, OCCAMS, and GEO datasets.
#   Methodology: adaptation of analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R
#   Map: analysis/ANALYSIS_MAP.md
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   Rscript analysis/clinical/survival_integrated_analysis.R
#
# Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(survival)
library(GSVA)
library(survminer)
library(grid)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")
setwd(SCREF_REF_OUTS_DIR)

out_dir <- "clinical_integrated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out
}

run_gsva <- function(expr_mat, gene_sets) {
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) < 10) return(NULL)
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  GSVA::gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
}

run_cox <- function(df, feature_cols, dataset_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    if (split_method == "median") {
      med_val <- median(d[[feat]], na.rm = TRUE)
      d$split_val <- factor(ifelse(d[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(d[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      d <- d %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(d) < 20) next
      d$split_val <- factor(ifelse(d[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      d$split_val <- as.numeric(scale(d[[feat]])) # Scaled continuous for comparability
    }
    
    if (nrow(d) < 20 || var(as.numeric(d$split_val), na.rm = TRUE) == 0) next
    form <- if (split_method == "continuous") {
      as.formula("Surv(OS_time, OS_event) ~ split_val")
    } else {
      as.formula("Surv(OS_time, OS_event) ~ split_val")
    }
    
    fit <- try(coxph(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      dataset = dataset_name,
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
      logHR = ss$coefficients[1, "coef"],
      SE = ss$coefficients[1, "se(coef)"],
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

run_optimal_cut_tcga <- function(df, feature_cols, dataset_name, feature_type) {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    probs <- seq(0.20, 0.80, by = 0.05)
    best_p <- 1
    best_cut <- NA
    best_prob <- NA
    best_ss <- NULL
    
    for (prob in probs) {
      cutpoint <- quantile(d[[feat]], probs = prob, na.rm = TRUE, names = FALSE)
      d$split_val <- factor(ifelse(d[[feat]] >= cutpoint, "High", "Low"), levels = c("Low", "High"))
      if (any(table(d$split_val) == 0)) next
      fit <- try(coxph(Surv(OS_time, OS_event) ~ split_val, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) next
      ss <- summary(fit)
      p_val <- ss$coefficients[1, "Pr(>|z|)"]
      if (p_val < best_p) {
        best_p <- p_val
        best_cut <- cutpoint
        best_prob <- prob
        best_ss <- ss
      }
    }
    
    if (is.na(best_cut)) next
    
    out[[feat]] <- data.frame(
      dataset = dataset_name,
      feature_type = feature_type,
      feature = feat,
      split_method = "optimal_cut",
      cutpoint = best_cut,
      quantile_prob = best_prob,
      logHR = best_ss$coefficients[1, "coef"],
      SE = best_ss$coefficients[1, "se(coef)"],
      HR = best_ss$coefficients[1, "exp(coef)"],
      P_value = best_p,
      n = best_ss$n,
      events = best_ss$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

project_optimal_cut_occams <- function(df, tcga_cuts_df, dataset_name) {
  out <- list()
  for (i in seq_len(nrow(tcga_cuts_df))) {
    feat <- tcga_cuts_df$feature[i]
    quantile_prob <- tcga_cuts_df$quantile_prob[i]
    feature_type <- tcga_cuts_df$feature_type[i]
    
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    cutpoint <- quantile(d[[feat]], probs = quantile_prob, na.rm = TRUE, names = FALSE)
    
    d$split_val <- factor(ifelse(d[[feat]] >= cutpoint, "High", "Low"), levels = c("Low", "High"))
    if (any(table(d$split_val) == 0)) next
    
    fit <- try(coxph(Surv(OS_time, OS_event) ~ split_val, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    
    out[[feat]] <- data.frame(
      dataset = dataset_name,
      feature_type = feature_type,
      feature = feat,
      split_method = "optimal_cut_projected",
      cutpoint = cutpoint,
      quantile_prob = quantile_prob,
      logHR = ss$coefficients[1, "coef"],
      SE = ss$coefficients[1, "se(coef)"],
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

run_multivariate_cox <- function(df, feature_cols, dataset_name, feature_type, split_method, covar_formula) {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    if (split_method == "median") {
      med_val <- median(d[[feat]], na.rm = TRUE)
      d$split_val <- factor(ifelse(d[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(d[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      d <- d %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(d) < 20) next
      d$split_val <- factor(ifelse(d[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      d$split_val <- as.numeric(scale(d[[feat]]))
    }
    
    if (nrow(d) < 20 || var(as.numeric(d$split_val), na.rm = TRUE) == 0) next
    
    form_str <- paste("Surv(OS_time, OS_event) ~ split_val", covar_formula)
    form <- as.formula(form_str)
    
    fit <- try(coxph(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    
    coef_name <- ifelse(split_method == "continuous", "split_val", "split_valHigh")
    if (!coef_name %in% rownames(ss$coefficients)) next
    
    out[[feat]] <- data.frame(
      dataset = dataset_name,
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
      logHR = ss$coefficients[coef_name, "coef"],
      SE = ss$coefficients[coef_name, "se(coef)"],
      HR = ss$coefficients[coef_name, "exp(coef)"],
      P_value = ss$coefficients[coef_name, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

run_random_effects_meta <- function(df, target_datasets) {
  df <- df %>% filter(dataset %in% target_datasets, !is.na(logHR), !is.na(SE), SE > 0)
  if (nrow(df) == 0 || !all(target_datasets %in% unique(df$dataset))) return(data.frame())
  yi <- df$logHR
  vi <- df$SE ^ 2
  k <- length(yi)
  wi <- 1 / vi
  mu_fe <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - mu_fe) ^ 2)
  c_val <- sum(wi) - (sum(wi ^ 2) / sum(wi))
  tau2 <- if (k > 1 && is.finite(c_val) && c_val > 0) max((Q - (k - 1)) / c_val, 0) else 0
  wi_re <- 1 / (vi + tau2)
  mu_re <- sum(wi_re * yi) / sum(wi_re)
  se_re <- sqrt(1 / sum(wi_re))
  z_re <- mu_re / se_re
  p_re <- 2 * pnorm(-abs(z_re))
  ci_low <- mu_re - 1.96 * se_re
  ci_high <- mu_re + 1.96 * se_re
  i2 <- if (k > 1 && Q > 0) max((Q - (k - 1)) / Q, 0) * 100 else 0
  q_p <- if (k > 1) pchisq(Q, df = k - 1, lower.tail = FALSE) else NA_real_

  data.frame(
    feature_type = df$feature_type[1],
    feature = df$feature[1],
    split_method = df$split_method[1],
    k = k,
    logHR = mu_re,
    SE = se_re,
    HR = exp(mu_re),
    lower_CI = exp(ci_low),
    upper_CI = exp(ci_high),
    P_value = p_re,
    tau2 = tau2,
    Q = Q,
    Q_p_value = q_p,
    I2 = i2,
    n_total = sum(df$n, na.rm = TRUE),
    events_total = sum(df$events, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_volcano_panel <- function(df, title_label) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + theme_void() + annotate("text", x=0, y=0, label="No models") + labs(title = title_label))
  }
  pdat <- df %>%
    mutate(sig = P_value < 0.05, log2HR = log2(HR), neglog10 = -log10(P_value))
  
  # Clean labels for plot
  pdat$label <- pdat$feature
  pdat$label[pdat$feature_type == "MP" & pdat$feature %in% names(SCREF_MP_DESCRIPTIONS)] <- 
    paste0(pdat$feature[pdat$feature_type == "MP" & pdat$feature %in% names(SCREF_MP_DESCRIPTIONS)], " ", 
           SCREF_MP_DESCRIPTIONS[pdat$feature[pdat$feature_type == "MP" & pdat$feature %in% names(SCREF_MP_DESCRIPTIONS)]])
  
  ggplot(pdat, aes(log2HR, neglog10)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_point(aes(shape = feature_type, colour = sig), size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 100) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = title_label, x = "log2(HR)", y = "-log10(p)", shape = "Feature Type") +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
}

message("Loading common genesets...")
mp_genes_path <- "Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds"
mp.genes <- readRDS(mp_genes_path)

state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Squamous-to-intestinal" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "Glandular-to-intestinal" = c("MP18b", "MP16", "MP17", "MP8b", "MP8+"),
  "Stress-adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

state.genes <- lapply(state_groups, function(mps) {
  unique(unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE))
})
state.genes <- state.genes[sapply(state.genes, length) >= 5]

ranked_markers_path <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv")
ranked_markers <- read.csv(ranked_markers_path, stringsAsFactors = FALSE)
dge_sets <- lapply(split(ranked_markers, ranked_markers$state), function(x) head(x$gene, 20))
names(dge_sets) <- paste0(names(dge_sets), " (DGE)")


all_cox_results <- list()
all_multivariate_results <- list()
all_meta_dfs <- list()
split_methods <- c("continuous", "median", "q1q4")
datasets <- c("TCGA", "OCCAMS", "GEO")

process_dataset <- function(dataset_name, expr_mat, meta_df, covar_formula = NULL) {
  message("Running GSVA for ", dataset_name, "...")
  mp_gs <- run_gsva(expr_mat, mp.genes)
  st_gs <- run_gsva(expr_mat, state.genes)
  dge_gs <- run_gsva(expr_mat, dge_sets)
  
  if (!is.null(mp_gs)) {
    mp_df <- as.data.frame(t(mp_gs)); mp_df$sample_id <- rownames(mp_df)
    meta_df <- meta_df %>% left_join(mp_df, by = "sample_id")
  }
  if (!is.null(st_gs)) {
    st_df <- as.data.frame(t(st_gs)); st_df$sample_id <- rownames(st_df)
    meta_df <- meta_df %>% left_join(st_df, by = "sample_id")
  }
  if (!is.null(dge_gs)) {
    dge_df <- as.data.frame(t(dge_gs)); dge_df$sample_id <- rownames(dge_df)
    meta_df <- meta_df %>% left_join(dge_df, by = "sample_id")
  }
  
  for (sm in split_methods) {
    if (!is.null(mp_gs)) {
      all_cox_results[[paste(dataset_name, "MP", sm)]] <<- run_cox(meta_df, rownames(mp_gs), dataset_name, "MP", sm)
      if (!is.null(covar_formula)) {
        all_multivariate_results[[paste(dataset_name, "MP", sm)]] <<- run_multivariate_cox(meta_df, rownames(mp_gs), dataset_name, "MP", sm, covar_formula)
      }
    }
    if (!is.null(st_gs)) {
      all_cox_results[[paste(dataset_name, "State", sm)]] <<- run_cox(meta_df, rownames(st_gs), dataset_name, "State", sm)
      if (!is.null(covar_formula)) {
        all_multivariate_results[[paste(dataset_name, "State", sm)]] <<- run_multivariate_cox(meta_df, rownames(st_gs), dataset_name, "State", sm, covar_formula)
      }
    }
    if (!is.null(dge_gs)) {
      all_cox_results[[paste(dataset_name, "DGE", sm)]] <<- run_cox(meta_df, rownames(dge_gs), dataset_name, "State DGE", sm)
      if (!is.null(covar_formula)) {
        all_multivariate_results[[paste(dataset_name, "DGE", sm)]] <<- run_multivariate_cox(meta_df, rownames(dge_gs), dataset_name, "State DGE", sm, covar_formula)
      }
    }
  }
  
  all_meta_dfs[[dataset_name]] <<- meta_df
}

# TCGA
message("Processing TCGA...")
tcga_meta_path <- "TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds"
meta_tcga <- readRDS(tcga_meta_path)
if (!"HistologyGroup" %in% colnames(meta_tcga)) meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
meta_tcga <- meta_tcga %>% 
  filter(sample_type_code == "01", HistologyGroup == "EAC") %>%
  mutate(
    sample_id = sample_barcode,
    Stage_Simple = ifelse(is.na(Stage_Simple) | Stage_Simple == "", "Unknown", Stage_Simple)
  )

tcga_whole_path <- "TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
tcga_df <- data.table::fread(tcga_whole_path)
tpm_tcga <- as.matrix(tcga_df[, -1]); rownames(tpm_tcga) <- tcga_df$GeneSymbol
tpm_tcga <- tpm_tcga[, colnames(tpm_tcga) %in% meta_tcga$sample_id, drop = FALSE]

process_dataset("TCGA", tpm_tcga, meta_tcga, "+ factor(Stage_Simple)")

# OCCAMS
message("Processing OCCAMS...")
occams_meta_path <- "OCCAMS/clinical/tables/Auto_OCCAMS_qc_pass_subject_metadata.csv"
meta_occams <- read.csv(occams_meta_path, stringsAsFactors = FALSE)
meta_occams <- meta_occams %>% 
  mutate(OS_time = os_days, OS_event = os_event, sample_id = subject_id) %>%
  mutate(
    T_Group = ifelse(pretreatment_t_stage %in% c("T0", "T1", "T2"), "Early",
              ifelse(pretreatment_t_stage %in% c("T3", "T4"), "Advanced", "Unknown")),
    N_Group = ifelse(pretreatment_n_stage_tnm7 == "N0", "N0",
              ifelse(grepl("^N", pretreatment_n_stage_tnm7), "N+", "Unknown")),
    M_Group = ifelse(pretreatment_m_stage == "M0", "M0",
              ifelse(grepl("^M", pretreatment_m_stage), "M1", "Unknown"))
  )

occams_expr_path <- "OCCAMS/clinical/intermediate/Auto_OCCAMS_qc_pass_logcpm_gene_symbols.rds"
expr_occams <- readRDS(occams_expr_path)
expr_occams <- expr_occams[, colnames(expr_occams) %in% meta_occams$sample_id, drop = FALSE]

process_dataset("OCCAMS", expr_occams, meta_occams, "+ factor(T_Group) + factor(N_Group) + factor(M_Group)")

# GEO
message("Processing GEO...")
geo_meta_path <- "geo_survival/Auto_GSE19417_meta.rds"
meta_geo <- readRDS(geo_meta_path)
meta_geo <- meta_geo %>% 
  filter(HistologyGroup == "EAC", analysis_ready_for_survival) %>%
  mutate(sample_id = sample_geo_accession)

geo_expr_path <- "geo_survival/Auto_GSE19417_expr_gene.rds"
expr_geo <- readRDS(geo_expr_path)
expr_geo <- expr_geo[, colnames(expr_geo) %in% meta_geo$sample_id, drop = FALSE]

process_dataset("GEO", expr_geo, meta_geo)

message("Generating integrated plot...")
final_results <- bind_rows(all_cox_results)
write.csv(final_results, file.path(out_dir, "survival_integrated_analysis_results.csv"), row.names = FALSE)

plot_list <- list()
for (ds in datasets) {
  for (sm in split_methods) {
    df <- final_results %>% filter(dataset == ds, split_method == sm)
    plot_list[[paste(ds, sm)]] <- make_volcano_panel(df, paste(ds, "-", sm))
  }
}

open_pdf_device(file.path(out_dir, "survival_integrated_analysis_volcano.pdf"), width = 20, height = 20)
grid.newpage()
grid.draw(arrangeGrob(grobs = plot_list, ncol = 3, nrow = 3))
dev.off()

message("Generating multivariate integrated plot...")
final_multi <- bind_rows(all_multivariate_results)
if (nrow(final_multi) > 0) {
  write.csv(final_multi, file.path(out_dir, "survival_multivariate_analysis_results.csv"), row.names = FALSE)
  
  plot_list_multi <- list()
  for (ds in c("TCGA", "OCCAMS")) {
    for (sm in split_methods) {
      df <- final_multi %>% filter(dataset == ds, split_method == sm)
      plot_list_multi[[paste(ds, sm)]] <- make_volcano_panel(df, paste(ds, "(Multivar) -", sm))
    }
  }
  
  open_pdf_device(file.path(out_dir, "survival_multivariate_analysis_volcano.pdf"), width = 20, height = 13.33)
  grid.newpage()
  grid.draw(arrangeGrob(grobs = plot_list_multi, ncol = 3, nrow = 2))
  dev.off()
}

message("Generating meta-analysis integrated plots...")
# Process meta analysis for unadjusted models
meta_tcga_occams <- final_results %>% 
  filter(dataset %in% c("TCGA", "OCCAMS")) %>%
  group_by(feature_type, feature, split_method) %>%
  group_split(.keep = TRUE) %>%
  lapply(run_random_effects_meta, target_datasets = c("TCGA", "OCCAMS")) %>%
  bind_rows() %>%
  mutate(dataset = "TCGA + OCCAMS")

meta_all_three <- final_results %>% 
  filter(dataset %in% c("TCGA", "OCCAMS", "GEO")) %>%
  group_by(feature_type, feature, split_method) %>%
  group_split(.keep = TRUE) %>%
  lapply(run_random_effects_meta, target_datasets = c("TCGA", "OCCAMS", "GEO")) %>%
  bind_rows() %>%
  mutate(dataset = "TCGA + OCCAMS + GEO")

final_meta <- bind_rows(meta_tcga_occams, meta_all_three)

if (nrow(final_meta) > 0) {
  write.csv(final_meta, file.path(out_dir, "survival_meta_analysis_results.csv"), row.names = FALSE)
  
  plot_list_meta1 <- list()
  plot_list_meta2 <- list()
  
  for (sm in split_methods) {
    df1 <- final_meta %>% filter(dataset == "TCGA + OCCAMS", split_method == sm)
    plot_list_meta1[[paste("TCGA+OCCAMS", sm)]] <- make_volcano_panel(df1, paste("TCGA+OCCAMS -", sm))
    
    df2 <- final_meta %>% filter(dataset == "TCGA + OCCAMS + GEO", split_method == sm)
    plot_list_meta2[[paste("All", sm)]] <- make_volcano_panel(df2, paste("TCGA+OCCAMS+GEO -", sm))
  }
  
  open_pdf_device(file.path(out_dir, "survival_meta_analysis_volcano.pdf"), width = 20, height = 6.66)
  grid.newpage()
  grid.draw(arrangeGrob(grobs = plot_list_meta1, ncol = 3, nrow = 1))
  grid.newpage()
  grid.draw(arrangeGrob(grobs = plot_list_meta2, ncol = 3, nrow = 1))
  dev.off()
}

message("Running optimal cutpoint analysis...")
tcga_meta <- all_meta_dfs[["TCGA"]]
occams_meta <- all_meta_dfs[["OCCAMS"]]

mp_features <- names(mp.genes)
state_features <- names(state.genes)
dge_features <- names(dge_sets)

tcga_opt_mp <- run_optimal_cut_tcga(tcga_meta, intersect(mp_features, colnames(tcga_meta)), "TCGA", "MP")
tcga_opt_st <- run_optimal_cut_tcga(tcga_meta, intersect(state_features, colnames(tcga_meta)), "TCGA", "State")
tcga_opt_dge <- run_optimal_cut_tcga(tcga_meta, intersect(dge_features, colnames(tcga_meta)), "TCGA", "State DGE")

tcga_opt_all <- bind_rows(tcga_opt_mp, tcga_opt_st, tcga_opt_dge)

if (nrow(tcga_opt_all) > 0) {
  occams_opt_all <- project_optimal_cut_occams(occams_meta, tcga_opt_all, "OCCAMS")
  
  opt_results <- bind_rows(tcga_opt_all, occams_opt_all)
  write.csv(opt_results, file.path(out_dir, "survival_optimal_cut_results.csv"), row.names = FALSE)
  
  p1 <- make_volcano_panel(tcga_opt_all, "TCGA - Optimal Cutpoint Discovery")
  p2 <- make_volcano_panel(occams_opt_all, "OCCAMS - TCGA Cutpoint Projection")
  
  open_pdf_device(file.path(out_dir, "survival_optimal_cut_projection_volcano.pdf"), width = 16, height = 8)
  
  # Page 1: Volcano plots
  grid.newpage()
  grid.draw(arrangeGrob(grobs = list(p1, p2), ncol = 2, nrow = 1))
  
  # Subsequent pages: KM plots for every feature
  for (feat in unique(tcga_opt_all$feature)) {
    tcga_res <- tcga_opt_all %>% filter(feature == feat)
    occams_res <- occams_opt_all %>% filter(feature == feat)
    if (nrow(tcga_res) == 0 || nrow(occams_res) == 0) next
    
    tcga_cut <- tcga_res$cutpoint[1]
    tcga_prob <- tcga_res$quantile_prob[1]
    tcga_p <- tcga_res$P_value[1]
    
    occams_cut <- occams_res$cutpoint[1]
    occams_prob <- occams_res$quantile_prob[1]
    occams_p <- occams_res$P_value[1]
    
    d_tcga <- tcga_meta %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    d_tcga$split_val <- factor(ifelse(d_tcga[[feat]] >= tcga_cut, "High", "Low"), levels = c("Low", "High"))
    fit_tcga <- survival::survfit(survival::Surv(OS_time, OS_event) ~ split_val, data = d_tcga)
    km_tcga <- survminer::ggsurvplot(fit_tcga, data = d_tcga, risk.table = TRUE, 
      pval = paste0("p = ", format.pval(tcga_p, digits = 3)), conf.int = FALSE,
      palette = c("#377EB8", "#E41A1C"), xlab = "Time (Days)", ylab = "Survival Probability",
      title = paste0("TCGA: ", feat, "\nOptimal Cut = ", signif(tcga_cut, 4), "\n(Top ", signif((1-tcga_prob)*100, 3), "%)"),
      legend.labs = c("Low", "High"), ggtheme = theme_minimal())
      
    d_occams <- occams_meta %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), OS_time > 0)
    d_occams$split_val <- factor(ifelse(d_occams[[feat]] >= occams_cut, "High", "Low"), levels = c("Low", "High"))
    fit_occams <- survival::survfit(survival::Surv(OS_time, OS_event) ~ split_val, data = d_occams)
    km_occams <- survminer::ggsurvplot(fit_occams, data = d_occams, risk.table = TRUE, 
      pval = paste0("p = ", format.pval(occams_p, digits = 3)), conf.int = FALSE,
      palette = c("#377EB8", "#E41A1C"), xlab = "Time (Days)", ylab = "Survival Probability",
      title = paste0("OCCAMS: ", feat, "\nProjected Cut = ", signif(occams_cut, 4), "\n(Top ", signif((1-occams_prob)*100, 3), "%)"),
      legend.labs = c("Low", "High"), ggtheme = theme_minimal())
      
    res <- survminer::arrange_ggsurvplots(list(km_tcga, km_occams), print = FALSE, ncol = 2, nrow = 1)
    grid.newpage()
    grid.draw(res)
  }
  
  dev.off()
}

message("Done.")
