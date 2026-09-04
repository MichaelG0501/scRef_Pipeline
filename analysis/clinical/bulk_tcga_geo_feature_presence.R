####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/bulk_tcga_geo_feature_presence.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_bulk_tcga_geo_feature_presence.R
#
# Cross-platform bulk feature-presence visualization for QC-retained EAC samples
# from TCGA-ESCA whole-bulk RNA-seq and GEO GSE19417 microarray. Recomputes
# harmonized MP/state scores directly from shared-gene expression after
# dataset-wise scaling, then visualizes whether each feature is represented by
# top-scoring samples in both datasets.
#
# Inputs:
#   ref_outs/bulk_crossplatform/Auto_bulk_crossplatform_qc_sample_table.csv
#   ref_outs/tcga_esca_meta.rds
#   ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   ref_outs/geo_survival/Auto_GSE19417_meta.rds
#   ref_outs/geo_survival/Auto_GSE19417_expr_gene.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#   ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#
# Outputs:
#   ref_outs/bulk_crossplatform/presence/Auto_bulk_crossplatform_feature_presence_eac_only.pdf
#   ref_outs/bulk_crossplatform/presence/Auto_bulk_crossplatform_feature_presence_summary.csv
#   ref_outs/bulk_crossplatform/presence/Auto_bulk_crossplatform_feature_presence_top_calls.csv
#   ref_outs/bulk_crossplatform/presence/Auto_bulk_crossplatform_feature_presence_top_samples.csv
#   updates/new_updates/summaries/Auto_bulk_tcga_geo_feature_presence_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

out_dir <- file.path("bulk_crossplatform", "presence")
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

  ucell_scores <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds")
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

summarize_presence <- function(score_mat, feature_type, method_name, meta_df, feature_order, mp_desc) {
  if (is.null(score_mat) || nrow(score_mat) == 0 || ncol(score_mat) == 0) {
    return(list(summary = data.frame(), top_calls = data.frame(), top_samples = data.frame()))
  }

  score_df <- as.data.frame(t(score_mat))
  score_df$sample_id <- rownames(score_df)
  feature_names <- rownames(score_mat)

  analysis_df <- meta_df %>%
    select(sample_id, dataset) %>%
    left_join(score_df, by = "sample_id")

  score_long <- analysis_df %>%
    pivot_longer(cols = all_of(feature_names), names_to = "feature", values_to = "score")

  score_wide <- as.matrix(score_df[, feature_names, drop = FALSE])
  top_idx <- max.col(score_wide, ties.method = "first")
  top_calls <- data.frame(
    sample_id = score_df$sample_id,
    top_feature = feature_names[top_idx],
    top_score = score_wide[cbind(seq_len(nrow(score_wide)), top_idx)],
    stringsAsFactors = FALSE
  ) %>%
    left_join(meta_df %>% select(sample_id, dataset), by = "sample_id")
  top_calls$method <- method_name
  top_calls$feature_type <- feature_type
  top_calls$feature_label <- make_feature_label(top_calls$top_feature, feature_type, mp_desc)

  summary_df <- score_long %>%
    left_join(top_calls %>% select(sample_id, top_feature), by = "sample_id") %>%
    group_by(dataset, feature) %>%
    summarise(
      n_samples = n(),
      n_top = sum(top_feature == feature, na.rm = TRUE),
      frac_top = ifelse(n() > 0, n_top / n(), NA_real_),
      n_positive = sum(score > 0, na.rm = TRUE),
      frac_positive = ifelse(n() > 0, n_positive / n(), NA_real_),
      median_score = median(score, na.rm = TRUE),
      q90_score = as.numeric(quantile(score, probs = 0.9, na.rm = TRUE)),
      max_score = max(score, na.rm = TRUE),
      .groups = "drop"
    )
  summary_df$method <- method_name
  summary_df$feature_type <- feature_type
  summary_df$feature_label <- make_feature_label(summary_df$feature, feature_type, mp_desc)
  summary_df$feature_label <- factor(
    summary_df$feature_label,
    levels = rev(make_feature_label(feature_order, feature_type, mp_desc))
  )
  summary_df <- summary_df %>%
    arrange(feature_label, dataset)

  top_samples <- score_long %>%
    group_by(dataset, feature) %>%
    arrange(desc(score), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    mutate(rank_within_feature = row_number()) %>%
    ungroup()
  top_samples$method <- method_name
  top_samples$feature_type <- feature_type
  top_samples$feature_label <- make_feature_label(top_samples$feature, feature_type, mp_desc)

  list(summary = summary_df, top_calls = top_calls, top_samples = top_samples)
}

plot_presence_dotplot <- function(summary_df, title_text) {
  ggplot(summary_df, aes(dataset, feature_label)) +
    geom_point(aes(size = n_top, color = frac_top), alpha = 0.95) +
    geom_text(aes(label = ifelse(n_top > 0, n_top, "")), size = 3.1, color = "black") +
    scale_size_continuous(range = c(2.5, 11), breaks = pretty(summary_df$n_top, n = 4)) +
    scale_color_gradient(low = "#D9D9D9", high = "#B2182B", limits = c(0, 1)) +
    theme_classic(base_size = 12) +
    labs(
      title = title_text,
      subtitle = "Point size and label: number of top-scoring EAC samples per dataset; color: fraction of dataset",
      x = NULL,
      y = NULL,
      size = "Top samples",
      color = "Fraction top"
    )
}

plot_presence_support <- function(summary_df, title_text) {
  ggplot(summary_df, aes(dataset, feature_label, fill = q90_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(n_positive, "/", n_samples)), size = 3.1, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_classic(base_size = 12) +
    labs(
      title = title_text,
      subtitle = "Tile fill: 90th percentile score; label: samples with score > 0 / total EAC samples",
      x = NULL,
      y = NULL,
      fill = "90th pct.\nscore"
    ) +
    theme(panel.border = element_blank(), axis.line = element_blank())
}

make_page <- function(left_plot, right_plot, page_title) {
  gridExtra::arrangeGrob(
    left_plot,
    right_plot,
    ncol = 2,
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

qc_table <- fread(file.path("bulk_crossplatform", "Auto_bulk_crossplatform_qc_sample_table.csv")) %>%
  as.data.frame()

tcga_meta <- readRDS("tcga_esca_meta.rds") %>%
  transmute(
    sample_id = sample_barcode,
    dataset = "TCGA",
    histology_group = infer_tcga_histology(type),
    analysis_ready_for_survival = histology_group == "EAC" & !is.na(OS_time) & !is.na(OS_event)
  )

geo_meta <- readRDS(file.path("geo_survival", "Auto_GSE19417_meta.rds")) %>%
  transmute(
    sample_id = sample_geo_accession,
    dataset = "GEO_GSE19417",
    histology_group = HistologyGroup,
    analysis_ready_for_survival = analysis_ready_for_survival
  )

meta_all <- bind_rows(tcga_meta, geo_meta) %>%
  left_join(qc_table %>% select(sample_id, qc_remove), by = "sample_id") %>%
  mutate(
    qc_remove = ifelse(is.na(qc_remove), FALSE, qc_remove),
    dataset = factor(dataset, levels = c("TCGA", "GEO_GSE19417"))
  ) %>%
  filter(!qc_remove, histology_group == "EAC", analysis_ready_for_survival)

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

tcga_ids <- intersect(colnames(tcga_expr), meta_all$sample_id[meta_all$dataset == "TCGA"])
geo_ids <- intersect(colnames(geo_expr), meta_all$sample_id[meta_all$dataset == "GEO_GSE19417"])

expr_eac <- cbind(
  row_zscore(tcga_expr[shared_genes, tcga_ids, drop = FALSE]),
  row_zscore(geo_expr[shared_genes, geo_ids, drop = FALSE])
)

mp.genes <- readRDS("Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")

ordered_mp_list <- c(
  "MP1", "MP5", "MP13+", "MP2+", "MP14", "MP3+", "MP6+", "MP11+",
  "MP9+", "MP10+", "MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP12", "MP15"
)
extra_mps <- setdiff(names(mp.genes), ordered_mp_list)
retained_mps <- c(ordered_mp_list, extra_mps)

state_groups <- list(
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

state_order <- c(
  "Classic proliferation", "Basal to intestinal metaplasia",
  "SMG to intestinal metaplasia", "Stress adaptive", "Cancer-cell immune mimicry"
)

####################
# Current descriptions are authoritative in the centred grouping table.
####################
grouping_current <- read.csv("Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv", check.names = FALSE)
mp_desc <- setNames(grouping_current$description, grouping_current$mp)

tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_noreg <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
state_rel <- NULL

pan_mp_sets <- list()
new_state_gene_sets <- list()
candidate_new_states <- character(0)
####################
# Current workflow excludes historical 3CA relabel features.
####################

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

method_sets <- list(
  reference = list(
    mp = c(mp.genes, pan_mp_sets),
    state = reference_state_sets
  ),
  dge = list(
    mp = c(dge_sets$mp, pan_mp_sets),
    state = dge_sets$state
  )
)

page_list <- list()
all_summary <- list()
all_top_calls <- list()
all_top_samples <- list()

for (method_name in names(method_sets)) {
  this_sets <- method_sets[[method_name]]
  mp_scores <- score_gene_sets(expr_eac, this_sets$mp)
  state_scores <- score_gene_sets(expr_eac, this_sets$state)

  mp_presence <- summarize_presence(
    score_mat = mp_scores,
    feature_type = "MP",
    method_name = method_name,
    meta_df = meta_all,
    feature_order = names(this_sets$mp),
    mp_desc = mp_desc
  )
  state_presence <- summarize_presence(
    score_mat = state_scores,
    feature_type = "State",
    method_name = method_name,
    meta_df = meta_all,
    feature_order = names(this_sets$state),
    mp_desc = mp_desc
  )

  all_summary[[paste0(method_name, "_mp")]] <- mp_presence$summary
  all_summary[[paste0(method_name, "_state")]] <- state_presence$summary
  all_top_calls[[paste0(method_name, "_mp")]] <- mp_presence$top_calls
  all_top_calls[[paste0(method_name, "_state")]] <- state_presence$top_calls
  all_top_samples[[paste0(method_name, "_mp")]] <- mp_presence$top_samples
  all_top_samples[[paste0(method_name, "_state")]] <- state_presence$top_samples

  page_list[[paste0("mp_", method_name)]] <- make_page(
    plot_presence_dotplot(mp_presence$summary, paste0("MP presence by top-scoring samples: ", method_name)),
    plot_presence_support(mp_presence$summary, paste0("MP score support across EAC samples: ", method_name)),
    paste0("Cross-platform EAC MP presence (", method_name, " gene sets)")
  )
  page_list[[paste0("state_", method_name)]] <- make_page(
    plot_presence_dotplot(state_presence$summary, paste0("State presence by top-scoring samples: ", method_name)),
    plot_presence_support(state_presence$summary, paste0("State score support across EAC samples: ", method_name)),
    paste0("Cross-platform EAC state presence (", method_name, " gene sets)")
  )
}

summary_df <- bind_rows(all_summary)
top_calls_df <- bind_rows(all_top_calls)
top_samples_df <- bind_rows(all_top_samples)

write_grob_pdf(
  file.path(out_dir, "Auto_bulk_crossplatform_feature_presence_eac_only.pdf"),
  grob_list = list(
    page_list[["mp_reference"]],
    page_list[["mp_dge"]],
    page_list[["state_reference"]],
    page_list[["state_dge"]]
  ),
  width = 14,
  height = 8
)

fwrite(summary_df, file.path(out_dir, "Auto_bulk_crossplatform_feature_presence_summary.csv"))
fwrite(top_calls_df, file.path(out_dir, "Auto_bulk_crossplatform_feature_presence_top_calls.csv"))
fwrite(top_samples_df, file.path(out_dir, "Auto_bulk_crossplatform_feature_presence_top_samples.csv"))

summary_metrics <- summary_df %>%
  group_by(method, feature_type, dataset) %>%
  summarise(
    features_with_top_support = sum(n_top > 0, na.rm = TRUE),
    total_features = n(),
    samples_in_dataset = unique(n_samples)[1],
    .groups = "drop"
  )

fwrite(
  summary_metrics,
  "../updates/new_updates/summaries/Auto_bulk_tcga_geo_feature_presence_summary.csv"
)
####################
# Manuscript status override (31 Jul 2026): legacy nMP19/old-state workflow.
# Do not redirect paths and treat this as current. The exact current-centred
# replacement contract is recorded in
# analysis/oac_scatlas_paper/figure04/figure04f_bulk_presence_PLACEHOLDER.md.
####################
