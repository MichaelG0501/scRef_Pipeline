####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R
#   Methodology: adaptation of analysis/clinical/tcga_mp_state_survival_reg_noreg.R
#   Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Computes GSVA scores for centred metaprograms and state-union gene sets on
#   TCGA ESCA (EAC) primary tumours. Runs univariate Cox proportional hazards
#   models on the scores against Overall Survival (OS), generating volcano plots.
#
# Inputs:
#   - ref_outs/TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   - ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#   - ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   - ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv
#
# Outputs:
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_volcano_centred_mps.pdf
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_mp_cox_methods_splits.csv
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_optimal_cut_volcano_km.pdf
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_optimal_cut_km.pdf
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_optimal_cut_results.csv
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_model_data.rds
#   - ref_outs/Metaprogrammes_Results/centred/survival/Auto_task2_centred_mp_survival_run_summary.txt
#   - updates/new_updates/summaries/tcga_centred_mp_survival_summary.csv
#
# Cache/replot behavior:
#   Always rebuilds GSVA scores from log2(TPM+1) using the persistent TCGA matrix; no alternate-input fallback.
#
# Run command:
#   eval "$(~/miniforge3/bin/conda shell.bash hook)"
#   source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#   /opt/pbs/bin/qsub analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.sh
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
library(survminer)
library(GSVA)

source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

task_prefix <- "task2_centred_mp"
out_dir <- "Metaprogrammes_Results/centred/survival"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

write_grob_pdf <- function(path, grob_list, width, height) {
  grob_list <- grob_list[!vapply(grob_list, is.null, logical(1))]
  if (length(grob_list) == 0) {
    message("No grobs available to write: ", path)
    return()
  }
  open_pdf_device(path, width = width, height = height)
  for (g in grob_list) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  dev.off()
}

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out
}

run_gsva <- function(expr_mat, gene_sets) {
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) < 10) stop("Invalid GSVA expression matrix")
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) stop("No gene sets retain at least five TCGA genes")
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

run_cox <- function(df, feature_cols, mode_name, method_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) stop("Missing requested feature: ", feat)
    d <- df %>% filter(!is.na(OS_time), OS_time > 0, !is.na(OS_event), OS_event %in% c(0, 1), !is.na(.data[[feat]]), HistologyGroup == "EAC")
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) stop("Insufficient survival data for: ", feat)

    if (split_method == "median") {
      med_val <- median(d[[feat]], na.rm = TRUE)
      d$split_val <- factor(ifelse(d[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(d[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      d <- d %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(d) < 20) next
      d$split_val <- factor(ifelse(d[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      d[["split_val"]] <- as.numeric(scale(d[[feat]]))
    }

    if (nrow(d) < 20 || var(as.numeric(d[["split_val"]]), na.rm = TRUE) == 0) next
    form <- as.formula("Surv(OS_time, OS_event) ~ split_val")

    fit <- coxph(form, data = d)
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      mode = mode_name,
      method = method_name,
      cohort = "EAC",
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
      score_scaling = ifelse(split_method == "continuous", "per 1 SD increase", "high versus low"),
      model_formula = "Surv(OS_time, OS_event) ~ split_val",
      covariates = "none",
      HR = ss$coefficients[1, "exp(coef)"],
      CI_low = ss[["conf.int"]][1, "lower .95"],
      CI_high = ss[["conf.int"]][1, "upper .95"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, ttl) {
  if (nrow(df) == 0) return(NULL)
  pdat <- df %>%
    mutate(sig = P_value < 0.05, log2HR = log2(HR), neglog10 = -log10(P_value))
  ggplot(pdat, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_text_repel(aes(label = feature), size = 2.8, max.overlaps = 100) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = "log2(HR)", y = "-log10(p)")
}

make_tcga_page <- function(plot1, plot2, plot3, page_title) {
  gridExtra::arrangeGrob(
    plot1, plot2, plot3,
    ncol = 3,
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_feature_label <- function(x, feature_type, mp_desc) {
  if (feature_type == "MP") {
    d <- mp_desc[x]
    d[is.na(d) | d == ""] <- x[is.na(d) | d == ""]
    return(paste0(x, " ", d))
  }
  x
}

mp_desc <- SCREF_MP_DESCRIPTIONS

# Load Metadata
meta_path <- "TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds"
if (!file.exists(meta_path)) {
  stop("Metadata file not found: ", meta_path)
}
meta_tcga <- readRDS(meta_path)
required_meta_columns <- c("sample_barcode", "case_barcode", "sample_type_code", "HistologyGroup", "OS_time", "OS_event")
missing_meta_columns <- setdiff(required_meta_columns, colnames(meta_tcga))
if (length(missing_meta_columns) > 0) stop("Missing required TCGA metadata columns: ", paste(missing_meta_columns, collapse = ", "))


# Load Whole TCGA
whole_path <- "TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
if (!file.exists(whole_path)) {
  stop("Whole TPM file not found: ", whole_path)
}
tpm_df <- data.table::fread(whole_path)
tpm_whole <- as.matrix(tpm_df[, -1])
rownames(tpm_whole) <- tpm_df[["GeneSymbol"]]
storage.mode(tpm_whole) <- "numeric"
if (any(!is.finite(tpm_whole)) || any(tpm_whole < 0)) stop("TCGA TPM matrix contains invalid values")
if (!setequal(colnames(tpm_whole), meta_tcga[["sample_barcode"]])) stop("TCGA TPM and metadata sample identifiers differ")

mp_genes_path <- "Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds"
if(!file.exists(mp_genes_path)) {
    stop("MP genes file not found: ", mp_genes_path)
}
mp.genes <- readRDS(mp_genes_path)
retained_mps <- names(mp.genes)
if (length(mp.genes) != 17 || !setequal(retained_mps, names(SCREF_MP_DESCRIPTIONS))) stop("Expected exact current 17-MP panel")

state_groups <- SCREF_STATE_GROUPS

state.genes <- lapply(state_groups, function(mps) {
  unique(unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE))
})
state.genes <- state.genes[sapply(state.genes, length) >= 5]

all_res <- list()
split_methods <- c("continuous", "median", "q1q4")
panel_results_mp <- list()
panel_results_state <- list()

mode_name <- "noreg"
method_name <- "whole_tcga"
expr_mat <- log2(tpm_whole + 1)

mp_gs <- run_gsva(expr_mat, mp.genes)
st_gs <- run_gsva(expr_mat, state.genes)

if (!is.null(mp_gs) || !is.null(st_gs)) {
  merged_df <- meta_tcga %>% filter(sample_type_code == "01")
  if (anyDuplicated(merged_df[["sample_barcode"]]) || anyDuplicated(merged_df[["case_barcode"]])) stop("Primary-tumour metadata are not unique by sample and case")
  
  if (!is.null(mp_gs)) {
    mp_df <- as.data.frame(t(mp_gs))
    mp_df$sample_barcode <- rownames(mp_df)
    merged_df <- merged_df %>% left_join(mp_df, by = "sample_barcode")
    mp_cols <- intersect(colnames(as.data.frame(t(mp_gs))), colnames(merged_df))
  } else {
    mp_cols <- character(0)
  }
  
  if (!is.null(st_gs)) {
    st_df <- as.data.frame(t(st_gs))
    st_df$sample_barcode <- rownames(st_df)
    merged_df <- merged_df %>% left_join(st_df, by = "sample_barcode", suffix = c("", "_state"))
    st_cols <- intersect(colnames(as.data.frame(t(st_gs))), colnames(merged_df))
  } else {
    st_cols <- character(0)
  }
  
  for (sm in split_methods) {
    if (length(mp_cols) > 0) {
      cox_mp <- run_cox(merged_df, mp_cols, mode_name, method_name, "MP", split_method = sm)
      all_res[[paste("MP", sm, sep = "_")]] <- cox_mp
      
      this_mp <- cox_mp %>% filter(cohort == "EAC")
      if (nrow(this_mp) > 0) {
        mp_levels <- make_feature_label(retained_mps, "MP", mp_desc)
        this_mp$feature <- make_feature_label(this_mp$feature, "MP", mp_desc)
        all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
        this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
        panel_results_mp[[sm]] <- plot_volcano(
          this_mp,
          paste0("whole_tcga MP volcano (", sm, ")")
        )
      } else {
        panel_results_mp[[sm]] <- NULL
      }
    }
    
    if (length(st_cols) > 0) {
      cox_st <- run_cox(merged_df, st_cols, mode_name, method_name, "State", split_method = sm)
      all_res[[paste("State", sm, sep = "_")]] <- cox_st
      
      this_st <- cox_st %>% filter(cohort == "EAC")
      if (nrow(this_st) > 0) {
        this_st <- this_st %>% mutate(feature = factor(feature, levels = names(state_groups)))
        panel_results_state[[sm]] <- plot_volcano(
          this_st,
          paste0("whole_tcga State volcano (", sm, ")")
        )
      } else {
        panel_results_state[[sm]] <- NULL
      }
    }
  }
}

volcano_page_mp <- make_tcga_page(
  panel_results_mp[["continuous"]],
  panel_results_mp[["median"]],
  panel_results_mp[["q1q4"]],
  "Centred MP volcano: Whole TCGA"
)

volcano_page_state <- make_tcga_page(
  panel_results_state[["continuous"]],
  panel_results_state[["median"]],
  panel_results_state[["q1q4"]],
  "Centred State volcano: Whole TCGA"
)

write_grob_pdf(
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_centred_mps.pdf")),
  grob_list = list(volcano_page_mp, volcano_page_state),
  width = 18,
  height = 8
)

cox_res <- bind_rows(all_res)
if (nrow(cox_res) > 0) {
  cox_res$padj <- ave(
    cox_res$P_value,
    interaction(cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type, cox_res$split_method),
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(
    cox_res,
    file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_cox_methods_splits.csv")),
    row.names = FALSE
  )
  ####################
  summary_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/updates/new_updates/summaries"
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  compact_summary <- cox_res %>%
    group_by(feature_type, split_method) %>%
    summarise(
      n_features = n_distinct(feature),
      n_nominal_p_lt_0_05 = sum(P_value < 0.05, na.rm = TRUE),
      n_bh_fdr_lt_0_05 = sum(padj < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(
    compact_summary,
    file.path(summary_dir, "tcga_centred_mp_survival_summary.csv"),
    row.names = FALSE
  )
  ####################
}

message("Saved filtered TCGA whole-bulk MP volcano outputs for centred MPs.")


plot_optimal_cut_volcano_km <- function(model_data, feature_types, output_path, km_output_path) {
  specs <- data.frame(feature = unlist(feature_types, use.names = FALSE),
                      feature_type = rep(names(feature_types), lengths(feature_types)),
                      stringsAsFactors = FALSE)
  results <- dplyr::bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    feature <- specs$feature[i]
    d <- model_data[model_data$HistologyGroup == "EAC" & !is.na(model_data$OS_time) &
      !is.na(model_data$OS_event) & !is.na(model_data[[feature]]) & model_data$OS_time > 0, , drop = FALSE]
    if (nrow(d) < 50 || stats::sd(d[[feature]]) == 0) stop("Insufficient optimal-cut data for: ", feature)
    cuts <- unique(stats::quantile(d[[feature]], seq(0.20, 0.80, by = 0.05), na.rm = TRUE, names = FALSE))
    candidates <- dplyr::bind_rows(lapply(cuts, function(cutpoint) {
      d$score_group <- factor(ifelse(d[[feature]] >= cutpoint, "High", "Low"), levels = c("Low", "High"))
      if (any(table(d$score_group) == 0)) return(NULL)
      fit <- survival::coxph(survival::Surv(OS_time, OS_event) ~ score_group, data = d)
      ss <- summary(fit)
      data.frame(cutpoint = cutpoint, hazard_ratio = ss$conf.int[1, "exp(coef)"],
        ci_low = ss[["conf.int"]][1, "lower .95"], ci_high = ss[["conf.int"]][1, "upper .95"],
        p_value = ss$coefficients[1, "Pr(>|z|)"], n = fit$n, events = fit$nevent)
    }))
    if (nrow(candidates) == 0) return(NULL)
    cbind(specs[i, , drop = FALSE], candidates[order(candidates$p_value, candidates$cutpoint), ][1, ])
  }))
  if (nrow(results) == 0) stop("No valid TCGA optimal-cut survival models available")
  results <- results %>% mutate(covariates = "none", model_formula = "Surv(OS_time, OS_event) ~ score_group",
    p_value_interpretation = "exploratory unadjusted minimum Cox Wald p across searched cuts",
    nominal_p_lt_0_05 = p_value < 0.05, log2_hr = log2(hazard_ratio),
    neg_log10_p = -log10(pmax(p_value, .Machine[["double.xmin"]])), label = ifelse(feature_type == "MP", paste(feature, mp_desc[feature]), feature))
  make_panel <- function(x, title) ggplot(x, aes(log2_hr, neg_log10_p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_point(aes(colour = nominal_p_lt_0_05), size = 2.8, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 100) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) + labs(title = title, x = "log2(HR)", y = "-log10(minimum raw p)")
  grDevices::cairo_pdf(output_path, width = 18, height = 8, onefile = TRUE)
  panel_list <- lapply(names(feature_types), function(ft) {
    make_panel(results[results$feature_type == ft, , drop = FALSE], paste0(ft, " (Optimal Cut)"))
  })
  grid::grid.newpage()
  grid::grid.draw(gridExtra::arrangeGrob(grobs = panel_list, ncol = 3,
    top = grid::textGrob("TCGA EAC Survival (Exploratory Optimal Cutpoint; Univariable Cox)",
      gp = grid::gpar(fontsize = 14, fontface = "bold"))))
  dev.off()

  selected <- results %>% arrange(p_value) %>% group_by(hazard_ratio > 1) %>% slice_head(n = 5) %>% ungroup() %>% pull(feature) %>% unique()
  results[["selected_for_km"]] <- results[["feature"]] %in% selected
  grDevices::cairo_pdf(km_output_path, width = 10, height = 8, onefile = TRUE)
  first_km_page <- TRUE
  for (feature in selected) {
    d <- model_data[model_data$HistologyGroup == "EAC" & !is.na(model_data$OS_time) & !is.na(model_data$OS_event) &
      !is.na(model_data[[feature]]) & model_data$OS_time > 0, , drop = FALSE]
    result_index <- match(feature, results[["feature"]])
    cutpoint <- results[["cutpoint"]][result_index]
    selected_p <- results[["p_value"]][result_index]
    d$score_group <- factor(ifelse(d[[feature]] >= cutpoint, "High", "Low"), levels = c("Low", "High"))
    fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ score_group, data = d)
    title_label <- if (feature %in% names(mp_desc)) paste(feature, mp_desc[feature]) else feature
    km_plot <- survminer::ggsurvplot(fit, data = d, risk.table = TRUE, pval = paste0("Optimal-cut Cox Wald p = ", format.pval(selected_p, digits = 3)), conf.int = FALSE,
      palette = c("#377EB8", "#E41A1C"), xlab = "Time (Days)", ylab = "Overall Survival Probability",
      title = paste0(title_label, "\nOptimal cut = ", signif(cutpoint, 4)),
      legend.labs = c("Low", "High"), ggtheme = theme_minimal())
    print(km_plot, newpage = !first_km_page)
    first_km_page <- FALSE
  }
  dev.off()
  results
}

if (exists("merged_df") && exists("mp_cols") && exists("st_cols")) {
  ranked_markers_path <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv")
  if (!file.exists(ranked_markers_path)) stop("Missing ranked state markers: ", ranked_markers_path)
  ranked_markers <- read.csv(ranked_markers_path, stringsAsFactors = FALSE)
  dge_sets <- lapply(split(ranked_markers, ranked_markers$state), function(x) head(x$gene, 20))
  names(dge_sets) <- paste0(names(dge_sets), " (DGE)")
  dge_gs <- run_gsva(expr_mat, dge_sets)
  if (is.null(dge_gs)) stop("Unable to compute TCGA State DGE scores")
  dge_df <- as.data.frame(t(dge_gs))
  dge_df$sample_barcode <- rownames(dge_df)
  merged_df <- merged_df %>% left_join(dge_df, by = "sample_barcode")
  dge_cols <- setdiff(colnames(dge_df), "sample_barcode")
  saveRDS(merged_df, file.path(out_dir, paste0("Auto_", task_prefix, "_survival_model_data.rds")), compress = "xz")
  optimal_cut_results <- plot_optimal_cut_volcano_km(merged_df, list(MP = mp_cols, State = st_cols, `State DGE` = dge_cols),
    file.path(out_dir, paste0("Auto_", task_prefix, "_survival_optimal_cut_volcano_km.pdf")),
    file.path(out_dir, paste0("Auto_", task_prefix, "_survival_optimal_cut_km.pdf")))
  write.csv(optimal_cut_results, file.path(out_dir, paste0("Auto_", task_prefix, "_survival_optimal_cut_results.csv")), row.names = FALSE, na = "")
  eac_evaluable <- merged_df[merged_df[["HistologyGroup"]] == "EAC" & !is.na(merged_df[["OS_time"]]) & merged_df[["OS_time"]] > 0 & !is.na(merged_df[["OS_event"]]), , drop = FALSE]
  writeLines(c("TCGA centred survival run", paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "cohort=primary-tumour EAC", "expression_transform=log2(TPM+1)", "covariates=none", "optimal_cut_p_values=exploratory unadjusted minima", paste0("subjects=", nrow(eac_evaluable)), paste0("events=", sum(eac_evaluable[["OS_event"]])), "result=PASS"), file.path(out_dir, paste0("Auto_", task_prefix, "_survival_run_summary.txt")))
}
