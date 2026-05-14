####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/bulk_tcga_geo_meta_survival.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_bulk_tcga_geo_meta_survival.R
#
# Cross-platform bulk survival meta-analysis for TCGA-ESCA whole-bulk RNA-seq
# and GEO GSE19417 microarray. Recomputes harmonized MP/state scores directly
# from shared-gene expression after dataset-wise scaling, runs dataset-specific
# Cox models, and integrates TCGA + GEO using a random-effects meta-analysis.
#
# Inputs:
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_sample_table.csv
#   ref_outs/tcga_esca_meta.rds
#   ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_final_states.rds
#
# Outputs:
#   ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_dataset_specific_results.csv
#   ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_random_effects_results.csv
#   ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_direction_summary.csv
#   ref_outs/bulk_crossplatform/meta_survival/Auto_bulk_crossplatform_meta_volcano_eac_only.pdf
#   updates/new_updates/summaries/Auto_bulk_tcga_geo_meta_survival_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(survival)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

out_dir <- file.path("bulk_crossplatform", "meta_survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("../updates/new_updates/summaries", recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

write_grob_pdf <- function(path, grob_list, width, height) {
  grob_list <- grob_list[!vapply(grob_list, is.null, logical(1))]
  if (length(grob_list) == 0) stop("No grobs available to write: ", path)
  open_pdf_device(path, width = width, height = height)
  for (g in grob_list) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  dev.off()
}

infer_tcga_histology <- function(type_vec) {
  type_vec <- tolower(as.character(type_vec))
  out <- rep("Other", length(type_vec))
  out[grepl("adeno", type_vec)] <- "EAC"
  out[grepl("squamous", type_vec)] <- "ESCC"
  out
}

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1
  scaled <- sweep(mat, 1, row_means, "-")
  scaled <- sweep(scaled, 1, row_sds, "/")
  scaled[is.na(scaled)] <- 0
  scaled
}

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

score_gene_sets <- function(expr_mat, gene_sets) {
  valid_sets <- gene_sets[vapply(gene_sets, function(gs) {
    length(intersect(unique(gs), rownames(expr_mat))) >= 5
  }, logical(1))]
  if (length(valid_sets) == 0) return(NULL)
  score_mat <- vapply(valid_sets, function(gs) {
    keep_genes <- intersect(unique(gs), rownames(expr_mat))
    colMeans(expr_mat[keep_genes, , drop = FALSE])
  }, numeric(ncol(expr_mat)))
  if (is.null(dim(score_mat))) {
    score_mat <- matrix(score_mat, nrow = ncol(expr_mat), dimnames = list(colnames(expr_mat), names(valid_sets)))
  }
  score_mat <- t(score_mat)
  colnames(score_mat) <- colnames(expr_mat)
  rownames(score_mat) <- names(valid_sets)
  score_mat
}

run_single_dataset_cox <- function(df, feature_cols, dataset_name, method_name, feature_type, split_method) {
  out <- list()
  dat <- df %>%
    filter(
      dataset == dataset_name,
      histology_group == "EAC",
      analysis_ready_for_survival,
      !is.na(OS_time),
      !is.na(OS_event)
    )
  for (feat in feature_cols) {
    if (!feat %in% colnames(dat)) next
    dd <- dat %>% filter(!is.na(.data[[feat]]))
    if (nrow(dd) < 20 || var(dd[[feat]], na.rm = TRUE) == 0) next

    if (split_method == "median") {
      med_val <- median(dd[[feat]], na.rm = TRUE)
      dd$split_val <- factor(ifelse(dd[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(dd[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      dd <- dd %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(dd) < 20) next
      dd$split_val <- factor(ifelse(dd[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      dd$split_val <- dd[[feat]]
    }

    form <- if (split_method == "continuous") {
      as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`"))
    } else {
      as.formula("Surv(OS_time, OS_event) ~ split_val")
    }
    fit <- try(coxph(form, data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      dataset = dataset_name,
      method = method_name,
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
      logHR = ss$coefficients[1, "coef"],
      SE = ss$coefficients[1, "se(coef)"],
      HR = ss$coefficients[1, "exp(coef)"],
      lower_CI = ss$conf.int[1, "lower .95"],
      upper_CI = ss$conf.int[1, "upper .95"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

run_random_effects_meta <- function(df) {
  df <- df %>% filter(!is.na(logHR), !is.na(SE), SE > 0)
  required_datasets <- c("TCGA", "GEO_GSE19417")
  if (nrow(df) == 0 || !all(required_datasets %in% unique(df$dataset))) return(data.frame())
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
    method = df$method[1],
    feature_type = df$feature_type[1],
    feature = df$feature[1],
    split_method = df$split_method[1],
    k = k,
    meta_logHR = mu_re,
    meta_SE = se_re,
    meta_HR = exp(mu_re),
    meta_lower_CI = exp(ci_low),
    meta_upper_CI = exp(ci_high),
    meta_P_value = p_re,
    tau2 = tau2,
    Q = Q,
    Q_p_value = q_p,
    I2 = i2,
    n_total = sum(df$n, na.rm = TRUE),
    events_total = sum(df$events, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

plot_volcano <- function(df, title_text) {
  if (nrow(df) == 0) return(NULL)
  pdat <- df %>%
    mutate(
      log2HR = log2(meta_HR),
      neglog10 = -log10(meta_P_value),
      sig = meta_P_value < 0.05
    )
  ggplot(pdat, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_text_repel(aes(label = feature), size = 2.7, max.overlaps = 100) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = "log2(meta HR)", y = "-log10(meta p)")
}

make_placeholder_plot <- function(title_text) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0, label = "No model available", size = 5) +
    labs(title = title_text)
}

make_volcano_page <- function(left_plot, right_plot, page_title) {
  gridExtra::arrangeGrob(
    left_plot %||% make_placeholder_plot("Reference"),
    right_plot %||% make_placeholder_plot("DGE"),
    ncol = 2,
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_feature_label <- function(x, feature_type, mp_desc) {
  feature_type <- feature_type[1]
  if (feature_type == "MP") {
    d <- mp_desc[x]
    d[is.na(d)] <- x[is.na(d)]
    return(paste0(x, " ", d))
  }
  x
}

make_reference_state_sets <- function(mp.genes, state_groups, pan_mp_sets, candidate_new_states, new_state_gene_sets, state_order, retained_mps) {
  canonical_state_sets <- lapply(state_groups, function(mps) {
    genes_nmf <- unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE)
    genes_3ca <- unlist(pan_mp_sets[intersect(mps, names(pan_mp_sets))], use.names = FALSE)
    unique(c(genes_nmf, genes_3ca))
  })
  canonical_state_sets <- canonical_state_sets[vapply(canonical_state_sets, length, numeric(1)) >= 5]

  extra_state_mp_sets <- list()
  is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
  if (length(is_mp) > 0) extra_state_mp_sets <- mp.genes[is_mp]

  state_sets <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
  extra <- setdiff(names(state_sets), c(state_order, "Unresolved", "Hybrid"))
  extra_3ca <- names(new_state_gene_sets)[names(new_state_gene_sets) %in% extra]
  extra_retained <- retained_mps[retained_mps %in% setdiff(extra, extra_3ca)]
  other_extra <- setdiff(extra, c(extra_3ca, extra_retained))
  ord <- c(state_order, extra_3ca, extra_retained, other_extra)
  state_sets[ord[ord %in% names(state_sets)]]
}

make_dge_sets <- function(tmdata_all, state_noreg, state_rel, retained_mps, mp.genes, pan_mp_sets, state_groups, candidate_new_states, new_state_gene_sets, state_order) {
  common_cells <- intersect(Cells(tmdata_all), names(state_noreg))
  tmp <- subset(tmdata_all, cells = common_cells)
  tmp$state_for_dge <- state_noreg[Cells(tmp)]

  if (!is.null(state_rel)) {
    rel_cells <- intersect(Cells(tmp), names(state_rel))
    tmp$state_for_dge[rel_cells] <- state_rel[rel_cells]
  }

  tmp$state_for_dge <- as.character(tmp$state_for_dge)
  keep_state_cells <- unlist(lapply(split(Cells(tmp), tmp$state_for_dge), function(v) {
    sample(v, min(length(v), 600))
  }), use.names = FALSE)
  tmp <- subset(tmp, cells = keep_state_cells)
  DefaultAssay(tmp) <- "RNA"
  tmp <- FindVariableFeatures(tmp, nfeatures = 4000, verbose = FALSE)

  Idents(tmp) <- tmp$state_for_dge
  markers_state <- FindAllMarkers(
    tmp,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident = 2000,
    features = VariableFeatures(tmp),
    verbose = FALSE
  )

  state_sets <- markers_state %>%
    filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 100, with_ties = FALSE) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  state_list <- setNames(state_sets$genes, state_sets$cluster)

  ref_state <- make_reference_state_sets(
    mp.genes = mp.genes,
    state_groups = state_groups,
    pan_mp_sets = pan_mp_sets,
    candidate_new_states = candidate_new_states,
    new_state_gene_sets = new_state_gene_sets,
    state_order = state_order,
    retained_mps = retained_mps
  )
  missing_states <- setdiff(c(names(state_groups), candidate_new_states), names(state_list))
  for (st in missing_states) {
    if (!is.null(ref_state[[st]])) state_list[[st]] <- ref_state[[st]]
  }
  state_list <- state_list[names(ref_state)[names(ref_state) %in% names(state_list)]]

  ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
  keep_ucell_cols <- retained_mps[retained_mps %in% colnames(ucell_scores)]
  ucell_scores <- ucell_scores[intersect(rownames(ucell_scores), Cells(tmp)), keep_ucell_cols, drop = FALSE]
  topmp <- colnames(ucell_scores)[max.col(ucell_scores, ties.method = "first")]
  names(topmp) <- rownames(ucell_scores)
  state_by_mp <- tmp$state_for_dge[names(topmp)]
  topmp <- topmp[!state_by_mp %in% c("Unresolved", "Hybrid")]

  tmp_mp <- subset(tmp, cells = names(topmp))
  keep_mp_cells <- unlist(lapply(split(Cells(tmp_mp), topmp[Cells(tmp_mp)]), function(v) {
    sample(v, min(length(v), 600))
  }), use.names = FALSE)
  tmp_mp <- subset(tmp_mp, cells = keep_mp_cells)
  tmp_mp$mp_for_dge <- topmp[Cells(tmp_mp)]
  tmp_mp <- FindVariableFeatures(tmp_mp, nfeatures = 4000, verbose = FALSE)

  Idents(tmp_mp) <- tmp_mp$mp_for_dge
  markers_mp <- FindAllMarkers(
    tmp_mp,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident = 2000,
    features = VariableFeatures(tmp_mp),
    verbose = FALSE
  )

  mp_sets <- markers_mp %>%
    filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 100, with_ties = FALSE) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  mp_list <- setNames(mp_sets$genes, mp_sets$cluster)

  missing_mps <- setdiff(retained_mps, names(mp_list))
  if (length(missing_mps) > 0) {
    all_source_mps <- c(mp.genes, pan_mp_sets)
    avail_miss <- missing_mps[missing_mps %in% names(all_source_mps)]
    if (length(avail_miss) > 0) {
      mp_list <- c(mp_list, all_source_mps[avail_miss])
    }
  }
  mp_list <- mp_list[retained_mps[retained_mps %in% names(mp_list)]]

  list(state = state_list, mp = mp_list)
}

qc_table <- fread(file.path("bulk_crossplatform", "Auto_bulk_crossplatform_qc_sample_table.csv")) %>%
  as.data.frame()

tcga_meta <- readRDS("tcga_esca_meta.rds") %>%
  transmute(
    sample_id = sample_barcode,
    dataset = "TCGA",
    histology_group = infer_tcga_histology(type),
    source_histology = type,
    analysis_ready_for_survival = histology_group == "EAC" & !is.na(OS_time) & !is.na(OS_event),
    OS_time = OS_time,
    OS_event = OS_event
  )

geo_meta <- readRDS(file.path("geo_survival", "Auto_GSE19417_meta.rds")) %>%
  transmute(
    sample_id = sample_geo_accession,
    dataset = "GEO_GSE19417",
    histology_group = HistologyGroup,
    source_histology = HistologyGroup,
    analysis_ready_for_survival = analysis_ready_for_survival,
    OS_time = OS_time,
    OS_event = OS_event
  )

meta_all <- bind_rows(tcga_meta, geo_meta) %>%
  left_join(qc_table %>% select(sample_id, qc_remove, qc_reason), by = "sample_id") %>%
  mutate(
    qc_remove = ifelse(is.na(qc_remove), FALSE, qc_remove),
    qc_reason = ifelse(is.na(qc_reason), "", qc_reason),
    dataset = factor(dataset, levels = c("TCGA", "GEO_GSE19417"))
  )

keep_samples <- meta_all$sample_id[!meta_all$qc_remove]

tcga_df <- fread("cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tcga_expr <- as.matrix(tcga_df[, -1, with = FALSE])
rownames(tcga_expr) <- tcga_df[[1]]
storage.mode(tcga_expr) <- "numeric"
tcga_expr <- log2(tcga_expr + 1)

geo_expr <- readRDS(file.path("geo_survival", "Auto_GSE19417_expr_gene.rds"))
geo_expr <- as.matrix(geo_expr)
storage.mode(geo_expr) <- "numeric"

shared_genes <- intersect(rownames(tcga_expr), rownames(geo_expr))
shared_genes <- sort(shared_genes)
tcga_shared <- tcga_expr[shared_genes, intersect(colnames(tcga_expr), keep_samples), drop = FALSE]
geo_shared <- geo_expr[shared_genes, intersect(colnames(geo_expr), keep_samples), drop = FALSE]

tcga_keep_ids <- intersect(colnames(tcga_shared), meta_all$sample_id[meta_all$dataset == "TCGA" & !meta_all$qc_remove])
geo_keep_ids <- intersect(colnames(geo_shared), meta_all$sample_id[meta_all$dataset == "GEO_GSE19417" & !meta_all$qc_remove])
tcga_shared <- tcga_shared[, tcga_keep_ids, drop = FALSE]
geo_shared <- geo_shared[, geo_keep_ids, drop = FALSE]

eac_tcga_ids <- meta_all$sample_id[meta_all$dataset == "TCGA" & !meta_all$qc_remove & meta_all$histology_group == "EAC"]
eac_geo_ids <- meta_all$sample_id[meta_all$dataset == "GEO_GSE19417" & !meta_all$qc_remove & meta_all$histology_group == "EAC"]

expr_by_method <- list(
  full_cohort_reference = cbind(row_zscore(tcga_shared), row_zscore(geo_shared)),
  full_cohort_dge = cbind(row_zscore(tcga_shared), row_zscore(geo_shared)),
  eac_only_reference = cbind(
    row_zscore(tcga_shared[, intersect(colnames(tcga_shared), eac_tcga_ids), drop = FALSE]),
    row_zscore(geo_shared[, intersect(colnames(geo_shared), eac_geo_ids), drop = FALSE])
  ),
  eac_only_dge = cbind(
    row_zscore(tcga_shared[, intersect(colnames(tcga_shared), eac_tcga_ids), drop = FALSE]),
    row_zscore(geo_shared[, intersect(colnames(geo_shared), eac_geo_ids), drop = FALSE])
  )
)

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

ordered_mp_list <- c(
  "MP2",
  "MP17", "MP14", "MP5", "MP10", "MP8",
  "MP13", "MP12",
  "MP18", "MP16",
  "MP15",
  "MP1", "MP7", "MP9"
)
extra_mps <- setdiff(names(mp.genes), ordered_mp_list)
retained_mps <- c(ordered_mp_list, extra_mps)

state_groups <- list(
  "Classic Proliferative" = c("MP2", "X3CA_mp_30.Respiration.1"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating"
)

mp_desc <- c(
  "MP1" = "G2M Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP5" = "Epithelial IFN Resp.",
  "MP7" = "DNA Damage Repair",
  "MP8" = "Intestinal Diff.",
  "MP9" = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)",
  "X3CA_mp_12.Protein.maturation" = "Protein maturation",
  "X3CA_mp_17.EMT.III" = "EMT III",
  "X3CA_mp_30.Respiration.1" = "Respiration 1"
)

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
state_rel <- if (file.exists("Auto_final_states.rds")) readRDS("Auto_final_states.rds") else NULL

pan_mp_sets <- list()
new_state_gene_sets <- list()
candidate_new_states <- character(0)
nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
if (file.exists(nmf3ca_path)) {
  MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
  MP_list <- as.list(MP_df)
  MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
  names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))

  target_3ca_mps <- c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1")
  pan_mp_sets <- MP_list[intersect(target_3ca_mps, names(MP_list))]

  if (!is.null(state_rel)) {
    candidate_new_states <- setdiff(unique(as.character(state_rel)), c(names(state_groups), "Unresolved", "Hybrid", NA))
    if (length(candidate_new_states) > 0) {
      clean_map <- setNames(clean_3ca_name(names(MP_list)), names(MP_list))
      for (st in candidate_new_states) {
        if (st == "3CA_EMT_and_Protein_maturation") {
          keep_mps <- intersect(c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III"), names(MP_list))
          new_state_gene_sets[[st]] <- unique(unlist(MP_list[keep_mps], use.names = FALSE))
        } else {
          orig_name <- names(clean_map)[clean_map == st][1]
          if (!is.na(orig_name)) new_state_gene_sets[[st]] <- MP_list[[orig_name]]
        }
      }
    }
  }
}

reference_state_sets <- make_reference_state_sets(
  mp.genes = mp.genes,
  state_groups = state_groups,
  pan_mp_sets = pan_mp_sets,
  candidate_new_states = candidate_new_states,
  new_state_gene_sets = new_state_gene_sets,
  state_order = state_order,
  retained_mps = retained_mps
)

dge_sets <- make_dge_sets(
  tmdata_all = tmdata_all,
  state_noreg = state_noreg,
  state_rel = state_rel,
  retained_mps = retained_mps,
  mp.genes = mp.genes,
  pan_mp_sets = pan_mp_sets,
  state_groups = state_groups,
  candidate_new_states = candidate_new_states,
  new_state_gene_sets = new_state_gene_sets,
  state_order = state_order
)

all_dataset_results <- list()
split_methods <- c("continuous", "median", "q1q4")
method_order <- names(expr_by_method)
volcano_pages_mp <- list()
volcano_pages_state <- list()

for (split_method in split_methods) {
  panel_results <- list(MP = list(), State = list())

  for (method_name in method_order) {
    expr_mat <- expr_by_method[[method_name]]
    is_dge <- grepl("_dge$", method_name)
    mp_sets <- if (is_dge) c(dge_sets$mp, pan_mp_sets) else c(mp.genes, pan_mp_sets)
    state_sets <- if (is_dge) dge_sets$state else reference_state_sets

    mp_scores <- score_gene_sets(expr_mat, mp_sets)
    state_scores <- score_gene_sets(expr_mat, state_sets)
    if (!is.null(mp_scores)) {
      mp_df <- as.data.frame(t(mp_scores))
      mp_df$sample_id <- rownames(mp_df)
    } else {
      mp_df <- data.frame(sample_id = colnames(expr_mat))
    }
    if (!is.null(state_scores)) {
      state_df <- as.data.frame(t(state_scores))
      state_df$sample_id <- rownames(state_df)
    } else {
      state_df <- data.frame(sample_id = colnames(expr_mat))
    }

    analysis_df <- meta_all %>%
      filter(sample_id %in% colnames(expr_mat), !qc_remove) %>%
      left_join(mp_df, by = "sample_id") %>%
      left_join(state_df, by = "sample_id")

    mp_features <- names(mp_sets)[names(mp_sets) %in% colnames(analysis_df)]
    state_features <- names(state_sets)[names(state_sets) %in% colnames(analysis_df)]

    res_df <- bind_rows(
      run_single_dataset_cox(analysis_df, mp_features, "TCGA", method_name, "MP", split_method),
      run_single_dataset_cox(analysis_df, mp_features, "GEO_GSE19417", method_name, "MP", split_method),
      run_single_dataset_cox(analysis_df, state_features, "TCGA", method_name, "State", split_method),
      run_single_dataset_cox(analysis_df, state_features, "GEO_GSE19417", method_name, "State", split_method)
    )
    if (nrow(res_df) > 0) {
      res_df$feature_label <- ifelse(
        res_df$feature_type == "MP",
        make_feature_label(res_df$feature, "MP", mp_desc),
        res_df$feature
      )
      all_dataset_results[[paste(method_name, split_method, sep = "_")]] <- res_df
    }
  }
}

dataset_results_df <- bind_rows(all_dataset_results)
if (nrow(dataset_results_df) == 0) {
  stop("No dataset-specific survival results were generated.")
}

dataset_results_df <- dataset_results_df %>%
  group_by(dataset, method, feature_type, split_method) %>%
  mutate(padj = p.adjust(P_value, method = "BH")) %>%
  ungroup()

meta_results_df <- dataset_results_df %>%
  group_by(method, feature_type, feature, split_method) %>%
  group_split(.keep = TRUE) %>%
  lapply(run_random_effects_meta) %>%
  bind_rows()
if (nrow(meta_results_df) == 0) {
  stop("No random-effects meta-analysis results were generated.")
}

meta_results_df <- meta_results_df %>%
  mutate(feature_label = ifelse(
    feature_type == "MP",
    make_feature_label(feature, "MP", mp_desc),
    feature
  )) %>%
  group_by(method, feature_type, split_method) %>%
  mutate(meta_padj = p.adjust(meta_P_value, method = "BH")) %>%
  ungroup()

for (split_method in split_methods) {
  plot_df <- meta_results_df %>% filter(split_method == !!split_method)
  this_mp_ref <- plot_df %>% filter(feature_type == "MP", method == "eac_only_reference") %>% mutate(feature = feature_label)
  this_mp_dge <- plot_df %>% filter(feature_type == "MP", method == "eac_only_dge") %>% mutate(feature = feature_label)
  this_st_ref <- plot_df %>% filter(feature_type == "State", method == "eac_only_reference") %>% mutate(feature = feature_label)
  this_st_dge <- plot_df %>% filter(feature_type == "State", method == "eac_only_dge") %>% mutate(feature = feature_label)

  volcano_pages_mp[[split_method]] <- make_volcano_page(
    if (nrow(this_mp_ref) > 0) plot_volcano(this_mp_ref, paste0("Meta-analysis EAC MP: reference (", split_method, ")")) else NULL,
    if (nrow(this_mp_dge) > 0) plot_volcano(this_mp_dge, paste0("Meta-analysis EAC MP: DGE (", split_method, ")")) else NULL,
    paste0("Random-effects meta-analysis EAC MP volcano: ", split_method)
  )
  volcano_pages_state[[split_method]] <- make_volcano_page(
    if (nrow(this_st_ref) > 0) plot_volcano(this_st_ref, paste0("Meta-analysis EAC State: reference (", split_method, ")")) else NULL,
    if (nrow(this_st_dge) > 0) plot_volcano(this_st_dge, paste0("Meta-analysis EAC State: DGE (", split_method, ")")) else NULL,
    paste0("Random-effects meta-analysis EAC State volcano: ", split_method)
  )
}

ordered_volcano_pages <- c(
  unname(volcano_pages_mp[split_methods]),
  unname(volcano_pages_state[split_methods])
)
write_grob_pdf(
  file.path(out_dir, "Auto_bulk_crossplatform_meta_volcano_eac_only.pdf"),
  grob_list = ordered_volcano_pages,
  width = 14,
  height = 7.5
)

direction_summary <- dataset_results_df %>%
  filter(split_method == "continuous") %>%
  select(dataset, method, feature_type, feature, HR, P_value, padj) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(HR, P_value, padj),
    names_sep = "__"
  ) %>%
  left_join(
    meta_results_df %>%
      filter(split_method == "continuous") %>%
      select(method, feature_type, feature, meta_HR, meta_P_value, meta_padj, I2, tau2, Q, Q_p_value),
    by = c("method", "feature_type", "feature")
  )

fwrite(dataset_results_df, file.path(out_dir, "Auto_bulk_crossplatform_meta_dataset_specific_results.csv"))
fwrite(meta_results_df, file.path(out_dir, "Auto_bulk_crossplatform_meta_random_effects_results.csv"))
fwrite(direction_summary, file.path(out_dir, "Auto_bulk_crossplatform_meta_direction_summary.csv"))

summary_df <- data.frame(
  metric = c(
    "shared_genes",
    "qc_retained_eac_tcga",
    "qc_retained_eac_geo",
    "dataset_specific_results_total",
    "meta_results_total",
    "meta_nominal_hits",
    "meta_bh_hits"
  ),
  value = c(
    length(shared_genes),
    sum(meta_all$dataset == "TCGA" & !meta_all$qc_remove & meta_all$histology_group == "EAC" & meta_all$analysis_ready_for_survival),
    sum(meta_all$dataset == "GEO_GSE19417" & !meta_all$qc_remove & meta_all$histology_group == "EAC" & meta_all$analysis_ready_for_survival),
    nrow(dataset_results_df),
    nrow(meta_results_df),
    sum(meta_results_df$meta_P_value < 0.05, na.rm = TRUE),
    sum(meta_results_df$meta_padj < 0.05, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_df,
  "../updates/new_updates/summaries/Auto_bulk_tcga_geo_meta_survival_summary.csv"
)
