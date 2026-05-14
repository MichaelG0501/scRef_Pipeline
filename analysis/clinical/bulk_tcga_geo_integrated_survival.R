####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/bulk_tcga_geo_integrated_survival.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_bulk_tcga_geo_integrated_survival.R
#
# QC-aware cross-platform bulk survival analysis for TCGA-ESCA whole-bulk RNA-seq
# and GEO GSE19417 microarray. Recomputes harmonized MP/state scores directly
# from shared-gene expression after dataset-wise scaling and runs per-dataset and
# pooled Cox models.
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
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_survival_volcano_mp_<split>.pdf
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_survival_volcano_state_<split>.pdf
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_survival_results.csv
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_interaction_results.csv
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_direction_summary.csv
#   ref_outs/bulk_crossplatform/survival/Auto_bulk_crossplatform_mp_score_distributions.pdf
#   updates/new_updates/summaries/Auto_bulk_tcga_geo_integrated_survival_summary.csv
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

out_dir <- file.path("bulk_crossplatform", "survival")
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
      analysis = "dataset_specific",
      dataset = dataset_name,
      method = method_name,
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
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

run_combined_adjusted_cox <- function(df, feature_cols, method_name, feature_type, split_method) {
  out <- list()
  dat <- df %>%
    filter(
      histology_group == "EAC",
      analysis_ready_for_survival,
      !is.na(OS_time),
      !is.na(OS_event)
    ) %>%
    mutate(dataset = factor(dataset, levels = c("TCGA", "GEO_GSE19417")))
  for (feat in feature_cols) {
    if (!feat %in% colnames(dat)) next
    dd <- dat %>% filter(!is.na(.data[[feat]]))
    if (nrow(dd) < 30 || var(dd[[feat]], na.rm = TRUE) == 0) next

    if (split_method == "median") {
      med_val <- median(dd[[feat]], na.rm = TRUE)
      dd$split_val <- factor(ifelse(dd[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(dd[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      dd <- dd %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(dd) < 30) next
      dd$split_val <- factor(ifelse(dd[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      dd$split_val <- dd[[feat]]
    }

    form <- if (split_method == "continuous") {
      as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "` + dataset"))
    } else {
      as.formula("Surv(OS_time, OS_event) ~ split_val + dataset")
    }
    fit <- try(coxph(form, data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      analysis = "combined_adjusted",
      dataset = "Combined",
      method = method_name,
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
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

run_interaction_cox <- function(df, feature_cols, method_name, feature_type) {
  out <- list()
  dat <- df %>%
    filter(
      histology_group == "EAC",
      analysis_ready_for_survival,
      !is.na(OS_time),
      !is.na(OS_event)
    ) %>%
    mutate(dataset = factor(dataset, levels = c("TCGA", "GEO_GSE19417")))
  for (feat in feature_cols) {
    if (!feat %in% colnames(dat)) next
    dd <- dat %>% filter(!is.na(.data[[feat]]))
    if (nrow(dd) < 30 || var(dd[[feat]], na.rm = TRUE) == 0) next
    fit <- try(
      coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "` * dataset")), data = dd),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    coef_names <- rownames(ss$coefficients)
    int_idx <- grep(":", coef_names)
    if (length(int_idx) == 0) next
    out[[feat]] <- data.frame(
      method = method_name,
      feature_type = feature_type,
      feature = feat,
      interaction_HR = ss$coefficients[int_idx[1], "exp(coef)"],
      interaction_P_value = ss$coefficients[int_idx[1], "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, title_text) {
  if (nrow(df) == 0) return(NULL)
  pdat <- df %>%
    mutate(
      log2HR = log2(HR),
      neglog10 = -log10(P_value),
      sig = P_value < 0.05
    )
  ggplot(pdat, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_text_repel(aes(label = feature), size = 2.7, max.overlaps = 100) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = "log2(HR)", y = "-log10(p)")
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
qc_table$dataset <- factor(qc_table$dataset, levels = c("TCGA", "GEO_GSE19417"))

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

all_res <- list()
all_interactions <- list()
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

    res_parts <- list(
      run_single_dataset_cox(analysis_df, mp_features, "TCGA", method_name, "MP", split_method),
      run_single_dataset_cox(analysis_df, mp_features, "GEO_GSE19417", method_name, "MP", split_method),
      run_single_dataset_cox(analysis_df, state_features, "TCGA", method_name, "State", split_method),
      run_single_dataset_cox(analysis_df, state_features, "GEO_GSE19417", method_name, "State", split_method),
      run_combined_adjusted_cox(analysis_df, mp_features, method_name, "MP", split_method),
      run_combined_adjusted_cox(analysis_df, state_features, method_name, "State", split_method)
    )
    res_df <- bind_rows(res_parts)
    if (nrow(res_df) > 0) {
      res_df$feature_label <- ifelse(
        res_df$feature_type == "MP",
        make_feature_label(res_df$feature, "MP", mp_desc),
        res_df$feature
      )
      all_res[[paste(method_name, split_method, sep = "_")]] <- res_df
    }

    if (split_method == "continuous") {
      int_df <- bind_rows(
        run_interaction_cox(analysis_df, mp_features, method_name, "MP"),
        run_interaction_cox(analysis_df, state_features, method_name, "State")
      )
      if (nrow(int_df) > 0) {
        int_df$feature_label <- ifelse(
          int_df$feature_type == "MP",
          make_feature_label(int_df$feature, "MP", mp_desc),
          int_df$feature
        )
        all_interactions[[method_name]] <- int_df
      }
    }

    plot_df <- res_df %>%
      filter(analysis == "combined_adjusted") %>%
      mutate(feature = feature_label)
    this_mp <- plot_df %>% filter(feature_type == "MP")
    this_st <- plot_df %>% filter(feature_type == "State")
    if (nrow(this_mp) > 0) {
      panel_results$MP[[method_name]] <- plot_volcano(
        this_mp,
        paste0("Combined adjusted EAC: ", method_name, " MP volcano (", split_method, ")")
      )
    } else {
      panel_results$MP[[method_name]] <- NULL
    }
    if (nrow(this_st) > 0) {
      panel_results$State[[method_name]] <- plot_volcano(
        this_st,
        paste0("Combined adjusted EAC: ", method_name, " State volcano (", split_method, ")")
      )
    } else {
      panel_results$State[[method_name]] <- NULL
    }
  }

  volcano_pages_mp[[split_method]] <- make_volcano_page(
    panel_results$MP[["eac_only_reference"]],
    panel_results$MP[["eac_only_dge"]],
    paste0("Combined adjusted EAC MP volcano: ", split_method)
  )
  volcano_pages_state[[split_method]] <- make_volcano_page(
    panel_results$State[["eac_only_reference"]],
    panel_results$State[["eac_only_dge"]],
    paste0("Combined adjusted EAC State volcano: ", split_method)
  )
}

unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_continuous.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_median.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_q1q4.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_mp_continuous.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_mp_median.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_mp_q1q4.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_state_continuous.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_state_median.pdf"), force = TRUE)
unlink(file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_state_q1q4.pdf"), force = TRUE)

ordered_volcano_pages <- c(
  unname(volcano_pages_mp[split_methods]),
  unname(volcano_pages_state[split_methods])
)
write_grob_pdf(
  file.path(out_dir, "Auto_bulk_crossplatform_survival_volcano_eac_only.pdf"),
  grob_list = ordered_volcano_pages,
  width = 14,
  height = 7.5
)

results_df <- bind_rows(all_res)
if (nrow(results_df) == 0) {
  stop("No survival results were generated for the QC-retained combined dataset.")
}
results_df <- results_df %>%
  group_by(analysis, dataset, method, feature_type, split_method) %>%
  mutate(padj = p.adjust(P_value, method = "BH")) %>%
  ungroup()

interaction_df <- bind_rows(all_interactions)
if (nrow(interaction_df) > 0) {
  interaction_df <- interaction_df %>%
    group_by(method, feature_type) %>%
    mutate(interaction_padj = p.adjust(interaction_P_value, method = "BH")) %>%
    ungroup()
}

direction_summary <- results_df %>%
  filter(split_method == "continuous", analysis %in% c("dataset_specific", "combined_adjusted")) %>%
  select(analysis, dataset, method, feature_type, feature, HR, P_value, padj) %>%
  pivot_wider(
    names_from = c(analysis, dataset),
    values_from = c(HR, P_value, padj),
    names_sep = "__"
  )

if (nrow(interaction_df) > 0) {
  direction_summary <- direction_summary %>%
    left_join(
      interaction_df %>%
        select(method, feature_type, feature, interaction_HR, interaction_P_value, interaction_padj),
      by = c("method", "feature_type", "feature")
    )
}

fwrite(results_df, file.path(out_dir, "Auto_bulk_crossplatform_survival_results.csv"))
fwrite(direction_summary, file.path(out_dir, "Auto_bulk_crossplatform_direction_summary.csv"))
if (nrow(interaction_df) > 0) {
  fwrite(interaction_df, file.path(out_dir, "Auto_bulk_crossplatform_interaction_results.csv"))
}

score_plot_methods <- c("full_cohort_reference", "eac_only_reference")
score_plot_features <- c("MP1", "MP7", "MP9")
score_plot_df <- bind_rows(lapply(score_plot_methods, function(method_name) {
  expr_mat <- expr_by_method[[method_name]]
  mp_scores <- score_gene_sets(expr_mat, c(mp.genes, pan_mp_sets))
  if (is.null(mp_scores)) return(data.frame())
  tmp <- as.data.frame(t(mp_scores[intersect(score_plot_features, rownames(mp_scores)), , drop = FALSE]))
  tmp$sample_id <- rownames(tmp)
  tmp %>%
    left_join(meta_all, by = "sample_id") %>%
    filter(
      !qc_remove,
      histology_group == "EAC",
      analysis_ready_for_survival
    ) %>%
    mutate(method = method_name)
}))

if (nrow(score_plot_df) > 0) {
  long_score_plot_df <- score_plot_df %>%
    pivot_longer(cols = any_of(score_plot_features), names_to = "feature", values_to = "score")
  open_pdf_device(file.path(out_dir, "Auto_bulk_crossplatform_mp_score_distributions.pdf"), width = 10, height = 7)
  print(
    ggplot(long_score_plot_df, aes(dataset, score, fill = dataset)) +
      geom_violin(alpha = 0.6, linewidth = 0.3) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
      facet_grid(feature ~ method, scales = "free_y") +
      theme_classic(base_size = 12) +
      scale_fill_manual(values = c("TCGA" = "#4DAF4A", "GEO_GSE19417" = "#377EB8")) +
      labs(title = "QC-retained EAC MP score distributions", x = "Dataset", y = "Harmonized score")
  )
  dev.off()
}

summary_df <- data.frame(
  metric = c(
    "shared_genes",
    "qc_retained_eac_tcga",
    "qc_retained_eac_geo",
    "continuous_results_total",
    "continuous_nominal_hits",
    "continuous_bh_hits",
    "interaction_nominal_hits"
  ),
  value = c(
    length(shared_genes),
    sum(meta_all$dataset == "TCGA" & !meta_all$qc_remove & meta_all$histology_group == "EAC" & meta_all$analysis_ready_for_survival),
    sum(meta_all$dataset == "GEO_GSE19417" & !meta_all$qc_remove & meta_all$histology_group == "EAC" & meta_all$analysis_ready_for_survival),
    sum(results_df$split_method == "continuous"),
    sum(results_df$split_method == "continuous" & results_df$P_value < 0.05),
    sum(results_df$split_method == "continuous" & results_df$padj < 0.05),
    if (nrow(interaction_df) > 0) sum(interaction_df$interaction_P_value < 0.05) else 0
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_df,
  "../updates/new_updates/summaries/Auto_bulk_tcga_geo_integrated_survival_summary.csv"
)
