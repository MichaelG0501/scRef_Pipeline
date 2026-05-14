####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/geo_survival_mp_state_survival.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_geo_survival_clinical_mps_v2_reg_noreg.R
#
# GEO bulk-expression survival comparison using the same MP/state framework as the
# TCGA reg/noreg script, adapted to datasets with public GEO survival metadata.
#
# Current survival-enabled dataset:
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#
# Inputs:
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_reg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_final_states.rds
#   ref_outs/geo_survival/Auto_geo_survival_dataset_manifest.csv
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#
# Outputs:
#   ref_outs/geo_survival/geo_task2_survival/Auto_geo_task2_dge_overlap_analysis.pdf
#   ref_outs/geo_survival/geo_task2_survival/Auto_geo_task2_survival_volcano_mp_<split>.pdf
#   ref_outs/geo_survival/geo_task2_survival/Auto_geo_task2_survival_volcano_state_<split>.pdf
#   ref_outs/geo_survival/geo_task2_survival/Auto_geo_task2_survival_mp_state_cox_methods_noreg_splits.csv
#   updates/new_updates/summaries/Auto_geo_survival_clinical_mps_v2_reg_noreg_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(survival)
library(GSVA)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

task_prefix <- "geo_task2"
out_dir <- file.path("geo_survival", paste0(task_prefix, "_survival"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Use Cairo-backed PDFs to avoid intermittent base-pdf colorspace failures on HPC.
open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}
####################

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

requested_modes <- "noreg"

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

run_gsva <- function(expr_mat, gene_sets) {
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) < 10) return(NULL)
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

run_cox <- function(df, feature_cols, dataset_name, mode_name, method_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>%
      filter(
        analysis_ready_for_survival,
        HistologyGroup == "EAC",
        !is.na(OS_time),
        !is.na(OS_event),
        !is.na(.data[[feat]])
      )
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
      dataset = dataset_name,
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
    annotate("text", x = 0, y = 0, label = "No significant model to plot", size = 5) +
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

pretty_method_name <- function(x) {
  recode(
    x,
    full_cohort_reference = "Full cohort + reference sets",
    eac_only_reference = "EAC-only + reference sets",
    full_cohort_dge = "Full cohort + DGE sets",
    eac_only_dge = "EAC-only + DGE sets",
    .default = x
  )
}

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

make_feature_label <- function(x, feature_type) {
  if (feature_type == "MP") {
    d <- mp_desc[x]
    d[is.na(d)] <- x[is.na(d)]
    return(paste0(x, " ", d))
  }
  x
}

manifest_df <- read.csv(file.path("geo_survival", "Auto_geo_survival_dataset_manifest.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    survival_public = survival_public %in% TRUE,
    analysis_ready_for_survival = analysis_ready_for_survival %in% TRUE
  ) %>%
  filter(analysis_ready_for_survival)

if (nrow(manifest_df) == 0) stop("No GEO datasets are marked as survival-ready in Auto_geo_survival_dataset_manifest.csv")

geo_datasets <- setNames(vector("list", length = nrow(manifest_df)), manifest_df$dataset)
for (i in seq_len(nrow(manifest_df))) {
  dataset_id <- manifest_df$dataset[i]
  meta_df <- readRDS(file.path("geo_survival", paste0("Auto_", dataset_id, "_meta.rds")))
  expr_mat <- readRDS(file.path("geo_survival", paste0("Auto_", dataset_id, "_expr_gene.rds")))
  keep_cols <- meta_df$sample_geo_accession[meta_df$sample_geo_accession %in% colnames(expr_mat)]
  expr_mat <- expr_mat[, keep_cols, drop = FALSE]
  geo_datasets[[dataset_id]] <- list(meta = meta_df, expr = expr_mat)
}

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_reg <- readRDS("Auto_topmp_v2_reg_states_B.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
final_states_path <- "Auto_final_states.rds"
state_rel <- if (file.exists(final_states_path)) readRDS(final_states_path) else NULL

nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
pan_mp_sets <- list()
new_state_gene_sets <- list()
candidate_new_states <- character(0)

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

retained_3ca_order <- names(new_state_gene_sets)

ordered_states_for_plot <- function(state_vec) {
  extra <- setdiff(state_vec, c(state_order, "Unresolved", "Hybrid"))
  extra_3ca <- retained_3ca_order[retained_3ca_order %in% extra]
  extra_retained <- retained_mps[retained_mps %in% setdiff(extra, extra_3ca)]
  other_extra <- setdiff(extra, c(extra_3ca, extra_retained))
  c(state_order, extra_3ca, extra_retained, other_extra)
}

make_reference_state_sets <- function() {
  canonical_state_sets <- lapply(state_groups, function(mps) {
    genes_NM <- unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE)
    genes_3CA <- unlist(pan_mp_sets[intersect(mps, names(pan_mp_sets))], use.names = FALSE)
    unique(c(genes_NM, genes_3CA))
  })
  canonical_state_sets <- canonical_state_sets[sapply(canonical_state_sets, length) >= 5]

  extra_state_mp_sets <- list()
  is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
  if (length(is_mp) > 0) extra_state_mp_sets <- mp.genes[is_mp]

  base_sets <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
  ord <- ordered_states_for_plot(names(base_sets))
  base_sets[ord[ord %in% names(base_sets)]]
}

make_dge_sets <- function(mode_name) {
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

  expected_states <- c(names(state_groups), candidate_new_states)
  missing_states <- setdiff(expected_states, names(state_list))
  if (length(missing_states) > 0) {
    reference_state_sets <- make_reference_state_sets()
    for (ms in missing_states) {
      if (!is.null(reference_state_sets[[ms]])) state_list[[ms]] <- reference_state_sets[[ms]]
    }
  }
  if (length(state_list) > 0) {
    ord <- ordered_states_for_plot(names(state_list))
    state_list <- state_list[ord[ord %in% names(state_list)]]
  }

  ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
  keep_ucell_cols <- retained_mps[retained_mps %in% colnames(ucell_scores)]
  ucell_scores <- ucell_scores[intersect(rownames(ucell_scores), Cells(tmp)), keep_ucell_cols, drop = FALSE]
  topmp <- colnames(ucell_scores)[max.col(ucell_scores, ties.method = "first")]
  names(topmp) <- rownames(ucell_scores)
  state_by_mp <- tmp$state_for_dge[names(topmp)]
  unresolved_or_hybrid <- state_by_mp %in% c("Unresolved", "Hybrid")
  topmp <- topmp[!unresolved_or_hybrid]

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

  mp_ord <- retained_mps[retained_mps %in% names(mp_list)]
  if (length(mp_ord) > 0) mp_list <- mp_list[mp_ord]

  list(state = state_list, mp = mp_list)
}

dge_by_mode <- list()
overlap_pdf <- file.path(out_dir, paste0("Auto_", task_prefix, "_dge_overlap_analysis.pdf"))
open_pdf_device(overlap_pdf, width = 12, height = 10)

for (mode_name in requested_modes) {
  message("Calculating DGE sets for mode: ", mode_name)
  dge_sets <- make_dge_sets(mode_name)
  dge_by_mode[[mode_name]] <- dge_sets

  ref_mp <- mp.genes
  all_dge_mps <- retained_mps[retained_mps %in% names(dge_sets$mp)]
  all_ref_mps <- retained_mps[retained_mps %in% names(ref_mp)]
  mp_labels_dge <- make_feature_label(all_dge_mps, "MP")
  mp_labels_ref <- make_feature_label(all_ref_mps, "MP")

  ref_state <- make_reference_state_sets()
  expected_state_order <- ordered_states_for_plot(unique(c(names(dge_sets$state), names(ref_state))))
  all_dge_states <- expected_state_order[expected_state_order %in% names(dge_sets$state)]
  all_ref_states <- expected_state_order[expected_state_order %in% names(ref_state)]

  mat_mp_inter <- matrix(0, nrow = length(all_dge_mps), ncol = length(all_ref_mps))
  mat_mp_jaccard <- mat_mp_inter
  for (i in seq_along(all_dge_mps)) {
    for (j in seq_along(all_ref_mps)) {
      inter <- intersect(dge_sets$mp[[all_dge_mps[i]]], ref_mp[[all_ref_mps[j]]])
      uni <- union(dge_sets$mp[[all_dge_mps[i]]], ref_mp[[all_ref_mps[j]]])
      mat_mp_inter[i, j] <- length(inter)
      mat_mp_jaccard[i, j] <- if (length(uni) > 0) length(inter) / length(uni) else 0
    }
  }
  rownames(mat_mp_jaccard) <- mp_labels_dge
  colnames(mat_mp_jaccard) <- mp_labels_ref

  right_anno_mp <- rowAnnotation(`DGE Size` = anno_barplot(sapply(all_dge_mps, function(x) length(dge_sets$mp[[x]])), gp = gpar(fill = "#E69F00")), annotation_name_rot = 90)
  top_anno_mp <- HeatmapAnnotation(`Ref Size` = anno_barplot(sapply(all_ref_mps, function(x) length(ref_mp[[x]])), gp = gpar(fill = "#56B4E9")))
  draw(
    Heatmap(
      mat_mp_jaccard,
      name = "Jaccard",
      column_title = paste0("MP DGE vs Original MP Overlap (", mode_name, ")"),
      col = colorRamp2(c(0, max(mat_mp_jaccard, na.rm = TRUE)), c("white", "#004488")),
      top_annotation = top_anno_mp,
      right_annotation = right_anno_mp,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      rect_gp = gpar(col = "gray80", lwd = 0.5),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (mat_mp_inter[i, j] > 0) {
          txt_col <- ifelse(mat_mp_jaccard[i, j] > 0.25, "white", "black")
          lbl <- if (all_dge_mps[i] == all_ref_mps[j]) {
            paste0(mat_mp_inter[i, j], " / ", length(ref_mp[[all_ref_mps[j]]]))
          } else {
            mat_mp_inter[i, j]
          }
          grid.text(lbl, x, y, gp = gpar(fontsize = 8, col = txt_col))
        }
      }
    )
  )

  mat_st_inter <- matrix(0, nrow = length(all_dge_states), ncol = length(all_ref_states), dimnames = list(all_dge_states, all_ref_states))
  mat_st_jaccard <- mat_st_inter
  for (i in all_dge_states) {
    for (j in all_ref_states) {
      inter <- intersect(dge_sets$state[[i]], ref_state[[j]])
      uni <- union(dge_sets$state[[i]], ref_state[[j]])
      mat_st_inter[i, j] <- length(inter)
      mat_st_jaccard[i, j] <- if (length(uni) > 0) length(inter) / length(uni) else 0
    }
  }

  right_anno_st <- rowAnnotation(`DGE Size` = anno_barplot(sapply(all_dge_states, function(x) length(dge_sets$state[[x]])), gp = gpar(fill = "#E69F00")), annotation_name_rot = 90)
  top_anno_st <- HeatmapAnnotation(`Ref Size` = anno_barplot(sapply(all_ref_states, function(x) length(ref_state[[x]])), gp = gpar(fill = "#56B4E9")))
  draw(
    Heatmap(
      mat_st_jaccard,
      name = "Jaccard",
      column_title = paste0("State DGE vs Reference State Overlap (", mode_name, ")"),
      col = colorRamp2(c(0, max(mat_st_jaccard, na.rm = TRUE)), c("white", "darkgreen")),
      top_annotation = top_anno_st,
      right_annotation = right_anno_st,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      rect_gp = gpar(col = "gray80", lwd = 0.5),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (mat_st_inter[i, j] > 0) {
          txt_col <- ifelse(mat_st_jaccard[i, j] > 0.25, "white", "black")
          lbl <- if (all_dge_states[i] == all_ref_states[j]) {
            paste0(mat_st_inter[i, j], " / ", length(ref_state[[all_ref_states[j]]]))
          } else {
            mat_st_inter[i, j]
          }
          grid.text(lbl, x, y, gp = gpar(fontsize = 8, col = txt_col))
        }
      }
    )
  )
}
dev.off()

all_res <- list()
split_methods <- c("continuous", "median", "q1q4")
volcano_method_order <- c("eac_only_reference", "eac_only_dge")
volcano_pages_mp <- list()
volcano_pages_state <- list()

for (sm in split_methods) {
  panel_results_mp <- list()
  panel_results_state <- list()

  for (dataset_name in names(geo_datasets)) {
    geo_meta <- geo_datasets[[dataset_name]]$meta
    expr_all <- geo_datasets[[dataset_name]]$expr
    expr_all <- expr_all[, geo_meta$sample_geo_accession[geo_meta$sample_geo_accession %in% colnames(expr_all)], drop = FALSE]
    eac_samples <- geo_meta$sample_geo_accession[geo_meta$HistologyGroup == "EAC" & geo_meta$sample_geo_accession %in% colnames(expr_all)]
    expr_eac <- expr_all[, eac_samples, drop = FALSE]

    for (mode_name in requested_modes) {
      dge_sets <- dge_by_mode[[mode_name]]
      reference_state_sets <- make_reference_state_sets()

      method_inputs <- list(
        full_cohort_reference = expr_all,
        eac_only_reference = expr_eac,
        full_cohort_dge = expr_all,
        eac_only_dge = expr_eac
      )
      method_order <- names(method_inputs)

      for (method_name in method_order) {
        expr_mat <- method_inputs[[method_name]]
        use_dge <- grepl("_dge$", method_name)
        mp_sets <- if (use_dge) c(dge_sets$mp, pan_mp_sets) else c(mp.genes, pan_mp_sets)
        state_sets <- if (use_dge) dge_sets$state else reference_state_sets

        mp_gs <- run_gsva(expr_mat, mp_sets)
        st_gs <- run_gsva(expr_mat, state_sets)
        if (is.null(mp_gs) && is.null(st_gs)) next

        merged_df <- geo_meta
        if (!is.null(mp_gs)) {
          mp_df <- as.data.frame(t(mp_gs))
          mp_df$sample_geo_accession <- rownames(mp_df)
          merged_df <- merged_df %>% left_join(mp_df, by = "sample_geo_accession")
        }
        if (!is.null(st_gs)) {
          st_df <- as.data.frame(t(st_gs))
          st_df$sample_geo_accession <- rownames(st_df)
          merged_df <- merged_df %>% left_join(st_df, by = "sample_geo_accession", suffix = c("", "_state"))
        }

        mp_cols <- if (!is.null(mp_gs)) intersect(colnames(as.data.frame(t(mp_gs))), colnames(merged_df)) else character(0)
        st_cols <- if (!is.null(st_gs)) intersect(colnames(as.data.frame(t(st_gs))), colnames(merged_df)) else character(0)

        cox_mp <- run_cox(merged_df, mp_cols, dataset_name, mode_name, method_name, "MP", split_method = sm)
        cox_st <- run_cox(merged_df, st_cols, dataset_name, mode_name, method_name, "State", split_method = sm)
        cox_all <- bind_rows(cox_mp, cox_st)
        all_res[[paste(dataset_name, sm, mode_name, method_name, sep = "|:")]] <- cox_all

        this_mp <- cox_mp %>% filter(cohort == "EAC")
        if (nrow(this_mp) > 0) {
          mp_levels <- c(make_feature_label(retained_mps, "MP"), make_feature_label(names(pan_mp_sets), "MP"))
          this_mp$feature <- make_feature_label(this_mp$feature, "MP")
          all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
          this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
          panel_results_mp[[paste(dataset_name, mode_name, method_name, sep = "|")]] <- plot_volcano(
            this_mp,
            paste0(dataset_name, " [", sm, "] ", pretty_method_name(method_name), " MP volcano")
          )
        } else {
          panel_results_mp[[paste(dataset_name, mode_name, method_name, sep = "|")]] <- NULL
        }

        this_st <- cox_st %>%
          filter(cohort == "EAC") %>%
          filter(!feature %in% c("Unresolved", "Hybrid"))
        if (nrow(this_st) > 0) {
          this_st <- this_st %>% mutate(feature = factor(feature, levels = ordered_states_for_plot(unique(as.character(feature)))))
          panel_results_state[[paste(dataset_name, mode_name, method_name, sep = "|")]] <- plot_volcano(
            this_st,
            paste0(dataset_name, " [", sm, "] ", pretty_method_name(method_name), " State volcano")
          )
        } else {
          panel_results_state[[paste(dataset_name, mode_name, method_name, sep = "|")]] <- NULL
        }
      }
    }
  }
  dataset_name <- names(geo_datasets)[1]
  mode_name <- requested_modes[1]
  left_key <- paste(dataset_name, mode_name, volcano_method_order[1], sep = "|")
  right_key <- paste(dataset_name, mode_name, volcano_method_order[2], sep = "|")
  volcano_pages_mp[[sm]] <- make_volcano_page(
    panel_results_mp[[left_key]],
    panel_results_mp[[right_key]],
    paste0(dataset_name, " MP volcano: ", sm)
  )
  volcano_pages_state[[sm]] <- make_volcano_page(
    panel_results_state[[left_key]],
    panel_results_state[[right_key]],
    paste0(dataset_name, " State volcano: ", sm)
  )
}

unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_q1q4.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_mp_q1q4.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_state_q1q4.pdf")), force = TRUE)

ordered_volcano_pages <- c(
  unname(volcano_pages_mp[split_methods]),
  unname(volcano_pages_state[split_methods])
)
write_grob_pdf(
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_eac_only.pdf")),
  grob_list = ordered_volcano_pages,
  width = 14,
  height = 7.5
)

cox_res <- bind_rows(all_res)
if (nrow(cox_res) > 0) {
  cox_res$padj <- ave(
    cox_res$P_value,
    interaction(cox_res$dataset, cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type, cox_res$split_method),
    FUN = function(x) p.adjust(x, method = "BH")
  )
}

write.csv(
  cox_res,
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_cox_methods_noreg_splits.csv")),
  row.names = FALSE
)

summary_dir <- file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- if (nrow(cox_res) > 0) {
  cox_res %>%
    group_by(dataset, mode, split_method, method, cohort, feature_type) %>%
    summarise(n_tested = n(), n_sig_p005 = sum(P_value < 0.05, na.rm = TRUE), .groups = "drop")
} else {
  data.frame()
}

write.csv(
  summary_df,
  file.path(summary_dir, "Auto_geo_survival_clinical_mps_v2_reg_noreg_summary.csv"),
  row.names = FALSE
)

message("Saved GEO bulk MP/state survival outputs (noreg).")
