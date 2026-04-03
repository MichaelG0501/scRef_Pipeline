####################
# Auto_survival_clinical_mps_v2_reg_noreg_filtered.R
#
# TCGA whole-bulk MP/state survival workflow on QC-retained samples only.
# Mirrors the reference-vs-DGE structure of the main TCGA survival workflow,
# retains both malignant-only and whole-TCGA branches, and writes one combined
# volcano PDF with MP pages first and State pages second.
#
# Inputs:
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_sample_table.csv
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_reg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_final_states.rds
#   ref_outs/tcga_esca_meta.rds
#   ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#
# Outputs:
#   ref_outs/task2_filtered_survival/Auto_task2_filtered_survival_volcano_methods_reg_noreg.pdf
#   ref_outs/task2_filtered_survival/Auto_task2_filtered_survival_mp_state_cox_methods_splits.csv
#   updates/new_updates/summaries/Auto_survival_clinical_mps_v2_reg_noreg_filtered_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(survival)
library(GSVA)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "task2_filtered"
out_dir <- paste0(task_prefix, "_survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
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
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

run_cox <- function(df, feature_cols, mode_name, method_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), HistologyGroup == "EAC")
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
      d$split_val <- d[[feat]]
    }

    if (nrow(d) < 20 || var(as.numeric(d$split_val), na.rm = TRUE) == 0) next
    form <- if (split_method == "continuous") {
      as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`"))
    } else {
      as.formula("Surv(OS_time, OS_event) ~ split_val")
    }

    fit <- try(coxph(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      mode = mode_name,
      method = method_name,
      cohort = "EAC",
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

make_placeholder_plot <- function(title_text) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0, label = "No model available", size = 5) +
    labs(title = title_text)
}

make_tcga_page <- function(top_left, top_right, bottom_left, bottom_right, page_title) {
  gridExtra::arrangeGrob(
    top_left %||% make_placeholder_plot("Malignant reference"),
    top_right %||% make_placeholder_plot("Malignant DGE"),
    bottom_left %||% make_placeholder_plot("Whole TCGA reference"),
    bottom_right %||% make_placeholder_plot("Whole TCGA DGE"),
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

ordered_states_for_plot <- function(state_vec, state_order, retained_3ca_order, retained_mps) {
  extra <- setdiff(state_vec, c(state_order, "Unresolved", "Hybrid"))
  extra_3ca <- retained_3ca_order[retained_3ca_order %in% extra]
  extra_retained <- retained_mps[retained_mps %in% setdiff(extra, extra_3ca)]
  other_extra <- setdiff(extra, c(extra_3ca, extra_retained))
  c(state_order, extra_3ca, extra_retained, other_extra)
}

make_reference_state_sets <- function(mp.genes, pan_mp_sets, state_groups, candidate_new_states, new_state_gene_sets, state_order, retained_3ca_order, retained_mps) {
  canonical_state_sets <- lapply(state_groups, function(mps) {
    genes_nmf <- unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE)
    genes_3ca <- unlist(pan_mp_sets[intersect(mps, names(pan_mp_sets))], use.names = FALSE)
    unique(c(genes_nmf, genes_3ca))
  })
  canonical_state_sets <- canonical_state_sets[sapply(canonical_state_sets, length) >= 5]

  extra_state_mp_sets <- list()
  is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
  if (length(is_mp) > 0) extra_state_mp_sets <- mp.genes[is_mp]

  base_sets <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
  ord <- ordered_states_for_plot(names(base_sets), state_order, retained_3ca_order, retained_mps)
  base_sets[ord[ord %in% names(base_sets)]]
}

make_dge_sets <- function(mode_name, tmdata_all, state_reg, state_noreg, state_rel, retained_mps, mp.genes, pan_mp_sets, state_groups, candidate_new_states, new_state_gene_sets, state_order, retained_3ca_order) {
  state_vec <- if (mode_name == "reg") state_reg else state_noreg
  common_cells <- intersect(Cells(tmdata_all), names(state_vec))
  tmp <- subset(tmdata_all, cells = common_cells)
  tmp$state_for_dge <- state_vec[Cells(tmp)]

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

  reference_state_sets <- make_reference_state_sets(
    mp.genes = mp.genes,
    pan_mp_sets = pan_mp_sets,
    state_groups = state_groups,
    candidate_new_states = candidate_new_states,
    new_state_gene_sets = new_state_gene_sets,
    state_order = state_order,
    retained_3ca_order = retained_3ca_order,
    retained_mps = retained_mps
  )
  missing_states <- setdiff(c(names(state_groups), candidate_new_states), names(state_list))
  for (ms in missing_states) {
    if (!is.null(reference_state_sets[[ms]])) state_list[[ms]] <- reference_state_sets[[ms]]
  }
  state_list <- state_list[ordered_states_for_plot(names(state_list), state_order, retained_3ca_order, retained_mps)]
  state_list <- state_list[!is.na(names(state_list))]

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
  all_source_mps <- c(mp.genes, pan_mp_sets)
  avail_miss <- missing_mps[missing_mps %in% names(all_source_mps)]
  if (length(avail_miss) > 0) mp_list <- c(mp_list, all_source_mps[avail_miss])
  mp_list <- mp_list[retained_mps[retained_mps %in% names(mp_list)]]

  list(state = state_list, mp = mp_list)
}

requested_modes <- "noreg"

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

qc_table <- read.csv(file.path("bulk_crossplatform", "Auto_bulk_crossplatform_qc_sample_table.csv"), stringsAsFactors = FALSE)
filtered_ids <- qc_table$sample_id[qc_table$dataset == "TCGA" & qc_table$integration_keep]

meta_tcga <- readRDS("tcga_esca_meta.rds")
meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
meta_tcga <- meta_tcga %>% filter(sample_type_code == "01", sample_barcode %in% filtered_ids)

tpm_df <- fread("cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_whole <- as.matrix(tpm_df[, -1, with = FALSE])
rownames(tpm_whole) <- tpm_df[[1]]
tpm_whole <- tpm_whole[, intersect(colnames(tpm_whole), meta_tcga$sample_barcode), drop = FALSE]

mal_df <- fread("cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt")
tpm_mal <- as.matrix(mal_df[, -1, with = FALSE])
rownames(tpm_mal) <- mal_df[[1]]
tpm_mal[is.na(tpm_mal)] <- 0
tpm_mal <- tpm_mal[, intersect(colnames(tpm_mal), meta_tcga$sample_barcode), drop = FALSE]

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_reg <- readRDS("Auto_topmp_v2_reg_states_B.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
state_rel <- if (file.exists("Auto_final_states.rds")) readRDS("Auto_final_states.rds") else NULL

pan_mp_sets <- list()
new_state_gene_sets <- list()
candidate_new_states <- character(0)
retained_3ca_order <- character(0)
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
  retained_3ca_order <- names(new_state_gene_sets)
}

dge_by_mode <- list()
reference_state_by_mode <- list()
for (mode_name in requested_modes) {
  reference_state_by_mode[[mode_name]] <- make_reference_state_sets(
    mp.genes = mp.genes,
    pan_mp_sets = pan_mp_sets,
    state_groups = state_groups,
    candidate_new_states = candidate_new_states,
    new_state_gene_sets = new_state_gene_sets,
    state_order = state_order,
    retained_3ca_order = retained_3ca_order,
    retained_mps = retained_mps
  )
  dge_by_mode[[mode_name]] <- make_dge_sets(
    mode_name = mode_name,
    tmdata_all = tmdata_all,
    state_reg = state_reg,
    state_noreg = state_noreg,
    state_rel = state_rel,
    retained_mps = retained_mps,
    mp.genes = mp.genes,
    pan_mp_sets = pan_mp_sets,
    state_groups = state_groups,
    candidate_new_states = candidate_new_states,
    new_state_gene_sets = new_state_gene_sets,
    state_order = state_order,
    retained_3ca_order = retained_3ca_order
  )
}

all_res <- list()
split_methods <- c("continuous", "median", "q1q4")
method_order <- c(
  "malignant_filtered_reference",
  "malignant_filtered_dge",
  "whole_filtered_reference",
  "whole_filtered_dge"
)
volcano_pages_mp <- list()
volcano_pages_state <- list()

for (sm in split_methods) {
  panel_results_mp <- list()
  panel_results_state <- list()

  for (mode_name in requested_modes) {
    dge_sets <- dge_by_mode[[mode_name]]
    reference_state_sets <- reference_state_by_mode[[mode_name]]
    method_inputs <- list(
      malignant_filtered_reference = tpm_mal,
      malignant_filtered_dge = tpm_mal,
      whole_filtered_reference = tpm_whole,
      whole_filtered_dge = tpm_whole
    )

    for (method_name in method_order) {
      expr_mat <- method_inputs[[method_name]]
      use_dge <- grepl("_dge$", method_name)
      mp_sets <- if (use_dge) c(dge_sets$mp, pan_mp_sets) else c(mp.genes, pan_mp_sets)
      state_sets <- if (use_dge) dge_sets$state else reference_state_sets

      mp_gs <- run_gsva(expr_mat, mp_sets)
      st_gs <- run_gsva(expr_mat, state_sets)
      if (is.null(mp_gs) && is.null(st_gs)) next

      merged_df <- meta_tcga
      if (!is.null(mp_gs)) {
        mp_df <- as.data.frame(t(mp_gs))
        mp_df$sample_barcode <- rownames(mp_df)
        merged_df <- merged_df %>% left_join(mp_df, by = "sample_barcode")
      }
      if (!is.null(st_gs)) {
        st_df <- as.data.frame(t(st_gs))
        st_df$sample_barcode <- rownames(st_df)
        merged_df <- merged_df %>% left_join(st_df, by = "sample_barcode", suffix = c("", "_state"))
      }

      mp_cols <- if (!is.null(mp_gs)) intersect(colnames(as.data.frame(t(mp_gs))), colnames(merged_df)) else character(0)
      st_cols <- if (!is.null(st_gs)) intersect(colnames(as.data.frame(t(st_gs))), colnames(merged_df)) else character(0)

      cox_mp <- run_cox(merged_df, mp_cols, mode_name, method_name, "MP", split_method = sm)
      cox_st <- run_cox(merged_df, st_cols, mode_name, method_name, "State", split_method = sm)
      cox_all <- bind_rows(cox_mp, cox_st)
      all_res[[paste(sm, mode_name, method_name, sep = "|:")]] <- cox_all

      this_mp <- cox_mp %>% filter(cohort == "EAC")
      if (nrow(this_mp) > 0) {
        mp_levels <- c(make_feature_label(retained_mps, "MP", mp_desc), make_feature_label(names(pan_mp_sets), "MP", mp_desc))
        this_mp$feature <- make_feature_label(this_mp$feature, "MP", mp_desc)
        all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
        this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
        panel_results_mp[[paste(mode_name, method_name, sep = "|")]] <- plot_volcano(
          this_mp,
          paste0("[", sm, "] ", mode_name, " ", method_name, " MP volcano")
        )
      } else {
        panel_results_mp[[paste(mode_name, method_name, sep = "|")]] <- NULL
      }

      this_st <- cox_st %>%
        filter(cohort == "EAC") %>%
        filter(!feature %in% c("Unresolved", "Hybrid"))
      if (nrow(this_st) > 0) {
        this_st <- this_st %>%
          mutate(feature = factor(feature, levels = ordered_states_for_plot(unique(as.character(feature)), state_order, retained_3ca_order, retained_mps)))
        panel_results_state[[paste(mode_name, method_name, sep = "|")]] <- plot_volcano(
          this_st,
          paste0("[", sm, "] ", mode_name, " ", method_name, " State volcano")
        )
      } else {
        panel_results_state[[paste(mode_name, method_name, sep = "|")]] <- NULL
      }
    }
    volcano_pages_mp[[paste(mode_name, sm, sep = "|")]] <- make_tcga_page(
      panel_results_mp[[paste(mode_name, "malignant_filtered_reference", sep = "|")]],
      panel_results_mp[[paste(mode_name, "malignant_filtered_dge", sep = "|")]],
      panel_results_mp[[paste(mode_name, "whole_filtered_reference", sep = "|")]],
      panel_results_mp[[paste(mode_name, "whole_filtered_dge", sep = "|")]],
      paste0(mode_name, " MP volcano: ", sm)
    )
    volcano_pages_state[[paste(mode_name, sm, sep = "|")]] <- make_tcga_page(
      panel_results_state[[paste(mode_name, "malignant_filtered_reference", sep = "|")]],
      panel_results_state[[paste(mode_name, "malignant_filtered_dge", sep = "|")]],
      panel_results_state[[paste(mode_name, "whole_filtered_reference", sep = "|")]],
      panel_results_state[[paste(mode_name, "whole_filtered_dge", sep = "|")]],
      paste0(mode_name, " State volcano: ", sm)
    )
  }
}

unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_q1q4.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_q1q4.pdf")), force = TRUE)

ordered_volcano_pages <- c(
  unname(volcano_pages_mp[paste("noreg", split_methods, sep = "|")]),
  unname(volcano_pages_state[paste("noreg", split_methods, sep = "|")])
)
write_grob_pdf(
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_methods_reg_noreg.pdf")),
  grob_list = ordered_volcano_pages,
  width = 14,
  height = 10
)

cox_res <- bind_rows(all_res)
if (nrow(cox_res) > 0) {
  cox_res$padj <- ave(
    cox_res$P_value,
    interaction(cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type, cox_res$split_method),
    FUN = function(x) p.adjust(x, method = "BH")
  )
}
write.csv(
  cox_res,
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_cox_methods_splits.csv")),
  row.names = FALSE
)

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- if (nrow(cox_res) > 0) {
  cox_res %>%
    group_by(mode, split_method, method, cohort, feature_type) %>%
    summarise(n_tested = n(), n_sig_p005 = sum(P_value < 0.05, na.rm = TRUE), .groups = "drop")
} else {
  data.frame()
}
write.csv(
  summary_df,
  file.path(summary_dir, "Auto_survival_clinical_mps_v2_reg_noreg_filtered_summary.csv"),
  row.names = FALSE
)

message("Saved filtered TCGA whole-bulk MP/state volcano outputs.")
