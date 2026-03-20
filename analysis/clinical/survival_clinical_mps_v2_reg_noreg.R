####################
# Auto_survival_clinical_mps_v2_reg_noreg.R
#
# Continuous-Cox survival comparison with volcano plots for 4 methods:
#   (a) GSVA on CIBERSORTx deconvoluted malignant-only TCGA profile
#   (b) GSVA on whole TCGA profile (legacy/reference)
#   (c) GSVA on DGE-derived gene sets on malignant-only TCGA profile
#   (d) GSVA on DGE-derived gene sets on whole TCGA profile
#
# For both MP and State features; includes newly relabeled unresolved-derived states.
#
# Inputs:
#   ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_reg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_reg_mp_adj.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#   ref_outs/unresolved_states/Auto_unresolved_relabel_states.rds (optional)
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt
#
# Outputs:
#   ref_outs/Auto_survival_tcga_mp_state_volcano_methods_reg_noreg.pdf
#   ref_outs/Auto_survival_tcga_mp_state_cox_methods_reg_noreg.csv
#   updates/new_updates/summaries/Auto_survival_clinical_mps_v2_reg_noreg_summary.csv
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(GSVA)
library(Seurat)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
task_prefix <- "task2"
out_dir <- paste0(task_prefix, "_survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
requested_modes <- if (length(args) >= 1 && nzchar(args[1])) unlist(strsplit(args[1], ",")) else c("reg", "noreg")
requested_modes <- intersect(c("reg", "noreg"), requested_modes)
if (length(requested_modes) == 0) stop("No valid modes requested. Use: reg,noreg or reg or noreg")
requested_modes <- intersect("noreg", requested_modes)
if (length(requested_modes) == 0) requested_modes <- "noreg"

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
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

run_cox <- function(df, feature_cols, mode_name, method_name, feature_type) {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), HistologyGroup %in% c("EAC"))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    for (coh in c("EAC")) {
      dd <- d %>% filter(HistologyGroup == coh)
      if (nrow(dd) < 20 || var(dd[[feat]], na.rm = TRUE) == 0) next
      fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`")), data = dd), silent = TRUE)
      if (inherits(fit, "try-error")) next
      ss <- summary(fit)
      out[[paste(feat, coh, sep = "|")]] <- data.frame(
        mode = mode_name,
        method = method_name,
        cohort = coh,
        feature_type = feature_type,
        feature = feat,
        HR = ss$coefficients[1, "exp(coef)"],
        P_value = ss$coefficients[1, "Pr(>|z|)"],
        n = fit$n,
        events = fit$nevent,
        stringsAsFactors = FALSE
      )
    }
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

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
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
  "MP1"  = "G2M Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP5"  = "Epithelial IFN Resp.",
  "MP7"  = "DNA Damage Repair",
  "MP8"  = "Intestinal Diff.",
  "MP9"  = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)"
)

# Reorder MPs to follow state definition order + MP number
# State order:
# 1. Classic Proliferative (MP2)
# 2. Basal to Intestinal (MP17, MP14, MP5, MP10, MP8)
# 3. Stress-adaptive (MP13, MP12)
# 4. SMG-like (MP18, MP16)
# 5. Immune Infiltrating (MP15)
# Remaining: MP1, MP7, MP9
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
retained_mps <- retained_mps[retained_mps %in% names(mp.genes)]

meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)

tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_whole <- as.matrix(tpm_df[, -1])
rownames(tpm_whole) <- tpm_df$GeneSymbol

mal_df <- data.table::fread("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt")
tpm_mal <- as.matrix(mal_df[, -1])
rownames(tpm_mal) <- mal_df$GeneSymbol
tpm_mal[is.na(tpm_mal)] <- 0

tmdata_all <- readRDS("EAC_Ref_epi.rds")
unresolved_relabeled_path <- "task4_unresolved_states/Auto_task4_unresolved_relabel_states.rds"

state_reg <- readRDS("Auto_topmp_v2_reg_states_B.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")

if (file.exists(unresolved_relabeled_path)) {
  state_rel <- readRDS(unresolved_relabeled_path)
} else {
  state_rel <- NULL
}

new_state_gene_sets <- list()
candidate_new_states <- character(0)
if (!is.null(state_rel)) {
  candidate_new_states <- setdiff(unique(as.character(state_rel)), c(names(state_groups), "Unresolved", "Hybrid", NA))
  nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
  if (file.exists(nmf3ca_path) && length(candidate_new_states) > 0) {
    MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
    MP_list <- as.list(MP_df)
    MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
    names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
    clean_map <- setNames(clean_3ca_name(names(MP_list)), names(MP_list))
    keep_cols <- names(clean_map)[clean_map %in% candidate_new_states]
    if (length(keep_cols) > 0) {
      new_state_gene_sets <- MP_list[keep_cols]
      names(new_state_gene_sets) <- clean_map[keep_cols]
    }
  }
}

pan_mp_sets <- new_state_gene_sets

retained_3ca_order <- names(new_state_gene_sets)

make_feature_label <- function(x, feature_type) {
  if (feature_type == "MP") {
    d <- mp_desc[x]
    d[is.na(d)] <- x[is.na(d)]
    return(paste0(x, " ", d))
  }
  x
}

ordered_states_for_plot <- function(state_vec) {
  extra <- setdiff(state_vec, c(state_order, "Unresolved", "Hybrid"))
  extra_3ca <- retained_3ca_order[retained_3ca_order %in% extra]
  extra_retained <- retained_mps[retained_mps %in% setdiff(extra, extra_3ca)]
  other_extra <- setdiff(extra, c(extra_3ca, extra_retained))
  c(state_order, extra_3ca, extra_retained, other_extra)
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

  # Fallback for missing states (not significant in DGE but we want them)
  expected_states <- c(names(state_groups), candidate_new_states)
  missing_states <- setdiff(expected_states, names(state_list))
  if (length(missing_states) > 0) {
    canonical_state_sets <- lapply(state_groups, function(mps) {
      mps_use <- mps[mps %in% names(mp.genes)]
      unique(unlist(mp.genes[mps_use], use.names = FALSE))
    })
    canonical_state_sets <- canonical_state_sets[sapply(canonical_state_sets, length) >= 5]
    
    extra_state_mp_sets <- list()
    is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
    if (length(is_mp) > 0) {
      extra_state_mp_sets <- mp.genes[is_mp]
    }
    
    base_sets <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
    for (ms in missing_states) {
      if (!is.null(base_sets[[ms]])) {
        state_list[[ms]] <- base_sets[[ms]]
      }
    }
  }

  if (length(state_list) > 0) {
    ord <- ordered_states_for_plot(names(state_list))
    state_list <- state_list[ord[ord %in% names(state_list)]]
  }

  mp_adj <- readRDS(paste0("Auto_topmp_v2_", mode_name, "_mp_adj.rds"))
  mp_adj <- mp_adj[intersect(rownames(mp_adj), Cells(tmp)), retained_mps[retained_mps %in% colnames(mp_adj)], drop = FALSE]
  topmp <- colnames(mp_adj)[max.col(mp_adj, ties.method = "first")]
  names(topmp) <- rownames(mp_adj)
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
  
  # Fallback for missing MPs (not significant in DGE)
  missing_mps <- setdiff(retained_mps, names(mp_list))
  if (length(missing_mps) > 0) {
    # only those actually in mp.genes
    avail_miss <- missing_mps[missing_mps %in% names(mp.genes)]
    if (length(avail_miss) > 0) {
      mp_list <- c(mp_list, mp.genes[avail_miss])
    }
  }

  mp_ord <- retained_mps[retained_mps %in% names(mp_list)]
  if (length(mp_ord) > 0) mp_list <- mp_list[mp_ord]

  list(state = state_list, mp = mp_list)
}

all_res <- list()

pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_volcano_methods_noreg.pdf")), width = 9, height = 7)
for (mode_name in requested_modes) {
  message("Building DGE gene sets for mode: ", mode_name)
  dge_sets <- make_dge_sets(mode_name)

  method_inputs <- list(
    malignant_only = tpm_mal,
    whole_tcga = tpm_whole,
    dge_based_malignant = tpm_mal,
    dge_based_whole = tpm_whole
  )

  method_order <- c("malignant_only", "whole_tcga", "dge_based_malignant", "dge_based_whole")
  panel_results <- list()

  for (method_name in method_order) {
    message("Running method: ", method_name, " | mode: ", mode_name)
    expr_mat <- method_inputs[[method_name]]

    use_dge <- grepl("^dge_based", method_name)
    mp_sets <- if (use_dge) c(dge_sets$mp, pan_mp_sets) else c(mp.genes, pan_mp_sets)
    state_sets <- if (use_dge) {
      dge_sets$state
    } else {
      # Build canonical state gene sets from MP gene sets (NOT MP IDs)
      canonical_state_sets <- lapply(state_groups, function(mps) {
        mps_use <- mps[mps %in% names(mp.genes)]
        unique(unlist(mp.genes[mps_use], use.names = FALSE))
      })
      canonical_state_sets <- canonical_state_sets[sapply(canonical_state_sets, length) >= 5]
      
      extra_state_mp_sets <- list()
      is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
      if (length(is_mp) > 0) {
        extra_state_mp_sets <- mp.genes[is_mp]
      }

      base_sets <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
      ord <- ordered_states_for_plot(names(base_sets))
      base_sets[ord[ord %in% names(base_sets)]]
    }

    mp_gs <- run_gsva(expr_mat, mp_sets)
    st_gs <- run_gsva(expr_mat, state_sets)
    if (is.null(mp_gs) && is.null(st_gs)) next

    merged_df <- meta_tcga %>% filter(sample_type_code == "01")
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

    cox_mp <- run_cox(merged_df, mp_cols, mode_name, method_name, "MP")
    cox_st <- run_cox(merged_df, st_cols, mode_name, method_name, "State")
    cox_all <- bind_rows(cox_mp, cox_st)
    all_res[[paste(mode_name, method_name, sep = "|:")]] <- cox_all

    # cache per-method panels for strict final order (4 MP + 4 State)
    this_mp <- cox_mp %>% filter(cohort == "EAC")
    if (nrow(this_mp) > 0) {
      mp_levels <- c(make_feature_label(retained_mps, "MP"), make_feature_label(names(pan_mp_sets), "MP"))
      this_mp$feature <- make_feature_label(this_mp$feature, "MP")
      # Robust levels to ensure no missing entries
      all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
      this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
      panel_results[[paste0("MP|", method_name)]] <- plot_volcano(this_mp, paste0("[", mode_name, "] ", method_name, " MP volcano (EAC)"))
    } else {
      panel_results[[paste0("MP|", method_name)]] <- NULL
    }

    this_st <- cox_st %>% filter(cohort == "EAC")
    if (nrow(this_st) > 0) {
      this_st <- this_st %>%
        filter(!feature %in% c("Unresolved", "Hybrid")) %>%
        mutate(feature = factor(feature, levels = ordered_states_for_plot(unique(as.character(feature)))))
      panel_results[[paste0("State|", method_name)]] <- plot_volcano(this_st, paste0("[", mode_name, "] ", method_name, " State volcano (EAC)"))
    } else {
      panel_results[[paste0("State|", method_name)]] <- NULL
    }
  }

  # print exactly 8 panels in requested order
  for (m in method_order) {
    p <- panel_results[[paste0("MP|", m)]]
    if (!is.null(p)) print(p)
  }
  for (m in method_order) {
    p <- panel_results[[paste0("State|", m)]]
    if (!is.null(p)) print(p)
  }
}
dev.off()

cox_res <- bind_rows(all_res)
if (nrow(cox_res) > 0) {
  cox_res$padj <- ave(cox_res$P_value, interaction(cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type),
                      FUN = function(x) p.adjust(x, method = "BH"))
}
write.csv(cox_res, file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_cox_methods_noreg.csv")), row.names = FALSE)

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  "updates", "new_updates", "summaries"
)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- if (nrow(cox_res) > 0) {
  cox_res %>%
    group_by(mode, method, cohort, feature_type) %>%
    summarise(n_tested = n(), n_sig_p005 = sum(P_value < 0.05, na.rm = TRUE), .groups = "drop")
} else {
  data.frame()
}
write.csv(summary_df, file.path(summary_dir, "Auto_survival_clinical_mps_v2_reg_noreg_summary.csv"), row.names = FALSE)

message("Saved 4-method MP/state continuous-Cox volcano outputs (noreg).")
