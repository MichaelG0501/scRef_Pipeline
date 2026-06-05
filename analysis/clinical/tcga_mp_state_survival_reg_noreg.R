####################
# Analysis registry:
#   Status: active
#   Script: analysis/clinical/tcga_mp_state_survival_reg_noreg.R
#   Methodology: analysis/methodology/clinical/clinical_bulk_and_association_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

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
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(survival)
library(GSVA)
library(Seurat)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")
task_prefix <- "task2"
out_dir <- paste0(task_prefix, "_survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

requested_modes <- "noreg"

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

run_cox <- function(df, feature_cols, mode_name, method_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), HistologyGroup %in% c("EAC"))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    # Apply split method
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
    
    for (coh in c("EAC")) {
      dd <- d %>% filter(HistologyGroup == coh)
      if (nrow(dd) < 20 || var(as.numeric(dd$split_val), na.rm = TRUE) == 0) next
      
      form <- if (split_method == "continuous") {
        as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`"))
      } else {
        as.formula("Surv(OS_time, OS_event) ~ split_val")
      }
      
      fit <- try(coxph(form, data = dd), silent = TRUE)
      if (inherits(fit, "try-error")) next
      ss <- summary(fit)
      out[[paste(feat, coh, sep = "|")]] <- data.frame(
        mode = mode_name,
        method = method_name,
        cohort = coh,
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
    geom_text_repel(aes(label = feature), size = 3.5, max.overlaps = 30, fontface = "bold") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 14) +
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
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
  )
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)

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
  "MP1"  = "G2M Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP5"  = "Epithelial IFN Resp.",
  "MP7"  = "DNA Damage Repair",
  "MP8"  = "Intestinal Diff.",
  "MP9"  = "G1S Cell Cycle",
  "MP10" = "Columnar Diff.",
  "MP12" = "Neuro-responsive Epi.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP15" = "Immune Infiltration",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP17" = "Basal-like Transition",
  "MP18" = "Secretory Diff. (Intest.)",
  "X3CA_mp_12.Protein.maturation" = "Protein maturation",
  "X3CA_mp_17.EMT.III"            = "EMT III",
  "X3CA_mp_30.Respiration.1"      = "Respiration 1"
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
retained_mps <- c(ordered_mp_list, extra_mps, "X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1")
retained_mps <- retained_mps[retained_mps %in% c(names(mp.genes), "X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1")]

meta_tcga <- readRDS("TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds")
if (!"HistologyGroup" %in% colnames(meta_tcga)) {
  meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
}

tpm_df <- data.table::fread("TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_whole <- as.matrix(tpm_df[, -1])
rownames(tpm_whole) <- tpm_df$GeneSymbol

if(file.exists("cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt")) { mal_df <- data.table::fread("cibersortx/CIBERSORTx_Job11_output/CIBERSORTxHiRes_Job11_Malignant_Window20.txt") } else { mal_df <- tpm_df }
tpm_mal <- as.matrix(mal_df[, -1])
rownames(tpm_mal) <- mal_df$GeneSymbol
tpm_mal[is.na(tpm_mal)] <- 0

tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states_path <- "Auto_final_states.rds"

state_reg <- if ("reg" %in% requested_modes) readRDS("Auto_topmp_v2_reg_states_B.rds") else NULL
state_noreg <- if ("noreg" %in% requested_modes) readRDS("Auto_topmp_v2_noreg_states_B.rds") else NULL

if (file.exists(final_states_path)) {
  state_rel <- readRDS(final_states_path)
} else {
  state_rel <- NULL
}

# Load 3CA MP gene sets for separate inclusion in volcano plots
nmf3ca_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
pan_mp_sets <- list() 
new_state_gene_sets <- list() 

if (file.exists(nmf3ca_path)) {
  MP_df <- read.csv(nmf3ca_path, check.names = FALSE)
  MP_list <- as.list(MP_df)
  MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
  names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
  
  # 1. Target 3CA MPs (separate for MP volcano)
  target_3ca_mps <- c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III", "X3CA_mp_30.Respiration.1")
  pan_mp_sets <- MP_list[intersect(target_3ca_mps, names(MP_list))]
  
  # 2. Respiration genes for state merge fallback
  respiration_genes <- if ("X3CA_mp_30.Respiration.1" %in% names(MP_list)) MP_list[["X3CA_mp_30.Respiration.1"]] else character(0)

  # 3. Handle additional candidate states from final relabeling
  if (!is.null(state_rel)) {
    candidate_new_states <- setdiff(unique(as.character(state_rel)), c(names(state_groups), "Unresolved", "Hybrid", NA))
    if (length(candidate_new_states) > 0) {
      clean_map <- setNames(clean_3ca_name(names(MP_list)), names(MP_list))
      for (st in candidate_new_states) {
        if (st == "3CA_EMT_and_Protein_maturation") {
          emt_prot_genes <- unique(unlist(MP_list[intersect(c("X3CA_mp_12.Protein.maturation", "X3CA_mp_17.EMT.III"), names(MP_list))]))
          new_state_gene_sets[[st]] <- emt_prot_genes
        } else {
          orig_name <- names(clean_map)[clean_map == st][1]
          if (!is.na(orig_name)) new_state_gene_sets[[st]] <- MP_list[[orig_name]]
        }
      }
    }
  }
}
retained_3ca_order <- names(new_state_gene_sets)

# pan_mp_sets was handles above

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
      genes_NM <- unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE)
      genes_3CA <- unlist(pan_mp_sets[intersect(mps, names(pan_mp_sets))], use.names = FALSE)
      unique(c(genes_NM, genes_3CA))
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

  ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
  ucell_scores <- ucell_scores[intersect(rownames(ucell_scores), Cells(tmp)), retained_mps[retained_mps %in% colnames(ucell_scores)], drop = FALSE]
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
  
  # Fallback for missing MPs (not significant in DGE)
  missing_mps <- setdiff(retained_mps, names(mp_list))
  if (length(missing_mps) > 0) {
    # Combine NM and 3CA source sets for fallback
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

# Pre-calculate DGE sets and Plot Overlap
dge_by_mode <- list()
message("Pre-calculating DGE sets and visualizing overlap...")
overlap_pdf <- file.path(out_dir, paste0("Auto_", task_prefix, "_dge_overlap_analysis.pdf"))
pdf(overlap_pdf, width = 12, height = 10)

for (mode_name in requested_modes) {
  message("  Calculating DGE sets for mode: ", mode_name)
  dge_sets <- make_dge_sets(mode_name)
  dge_by_mode[[mode_name]] <- dge_sets
  
  # ==========================================
  # Pre-processing: Enforce Strict Ordering & Labels
  # ==========================================
  # 1. MPs
  ref_mp <- mp.genes 
  all_dge_mps <- retained_mps[retained_mps %in% names(dge_sets$mp)]
  all_ref_mps <- retained_mps[retained_mps %in% names(ref_mp)]
  
  # Generate full descriptions for plotting labels
  mp_labels_dge <- make_feature_label(all_dge_mps, "MP")
  mp_labels_ref <- make_feature_label(all_ref_mps, "MP")
  
  # 2. States
  canonical_state_sets <- lapply(state_groups, function(mps) {
    # Match both original MPs and potential 3CA MPs in the groups
    mps_original <- intersect(mps, names(mp.genes))
    genes_original <- if(length(mps_original) > 0) unique(unlist(mp.genes[mps_original], use.names = FALSE)) else character(0)
    
    # Check for 3CA MPs in the state groups (e.g. Respiration 1)
    # We use new_state_gene_sets which was checked for these
    # but for canonical groups we should also check the already parsed respiration_genes
    mps_3ca <- intersect(mps, "3CA_mp_30 Respiration 1")
    genes_3ca <- if(length(mps_3ca) > 0) respiration_genes else character(0)
    
    unique(c(genes_original, genes_3ca))
  })
  canonical_state_sets <- canonical_state_sets[sapply(canonical_state_sets, length) >= 5]
  
  extra_state_mp_sets <- list()
  is_mp <- candidate_new_states[candidate_new_states %in% names(mp.genes)]
  if (length(is_mp) > 0) extra_state_mp_sets <- mp.genes[is_mp]
  
  ref_state <- c(canonical_state_sets, new_state_gene_sets, extra_state_mp_sets)
  
  expected_state_order <- ordered_states_for_plot(unique(c(names(dge_sets$state), names(ref_state))))
  all_dge_states <- expected_state_order[expected_state_order %in% names(dge_sets$state)]
  all_ref_states <- expected_state_order[expected_state_order %in% names(ref_state)]
  
  # ==========================================
  # 1. MP Overlap Heatmap (Strict Order + Diagonal Fractions)
  # ==========================================
  mat_mp_inter <- matrix(0, nrow = length(all_dge_mps), ncol = length(all_ref_mps))
  mat_mp_jaccard <- mat_mp_inter
  
  for (i in seq_along(all_dge_mps)) {
    for (j in seq_along(all_ref_mps)) {
      inter <- intersect(dge_sets$mp[[ all_dge_mps[i] ]], ref_mp[[ all_ref_mps[j] ]])
      uni <- union(dge_sets$mp[[ all_dge_mps[i] ]], ref_mp[[ all_ref_mps[j] ]])
      mat_mp_inter[i, j] <- length(inter)
      mat_mp_jaccard[i, j] <- if (length(uni) > 0) length(inter)/length(uni) else 0
    }
  }
  
  # Assign descriptive names for plot axes
  rownames(mat_mp_jaccard) <- mp_labels_dge
  colnames(mat_mp_jaccard) <- mp_labels_ref
  
  right_anno_mp <- rowAnnotation(`DGE Size` = anno_barplot(rowSums(mat_mp_inter > -1) * 0 + sapply(all_dge_mps, function(x) length(dge_sets$mp[[x]])), gp = gpar(fill = "#E69F00")), annotation_name_rot = 90)
  top_anno_mp <- HeatmapAnnotation(`Ref Size` = anno_barplot(colSums(mat_mp_inter > -1) * 0 + sapply(all_ref_mps, function(x) length(ref_mp[[x]])), gp = gpar(fill = "#56B4E9")))
  
  p1 <- Heatmap(mat_mp_jaccard, name = "Jaccard", 
                column_title = paste0("MP DGE vs Original MP Overlap"),
                col = colorRamp2(c(0, max(mat_mp_jaccard, na.rm = TRUE)), c("white", "#004488")),
                top_annotation = top_anno_mp,
                right_annotation = right_anno_mp,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_side = "left",
                rect_gp = gpar(col = "gray80", lwd = 0.5),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (mat_mp_inter[i, j] > 0) {
                    txt_col <- ifelse(mat_mp_jaccard[i, j] > 0.25, "white", "black")
                    # Check if it's the diagonal (matching MP to MP)
                    is_diag <- all_dge_mps[i] == all_ref_mps[j]
                    if (is_diag) {
                      # Format: Overlap / Ref Size
                      lbl <- paste0(mat_mp_inter[i, j], " / ", length(ref_mp[[ all_ref_mps[j] ]]))
                    } else {
                      lbl <- mat_mp_inter[i, j]
                    }
                    grid.text(lbl, x, y, gp = gpar(fontsize = 8, col = txt_col))
                  }
                })
  draw(p1)
  
  # ==========================================
  # 2. State Overlap Heatmap 
  # ==========================================
  mat_st_inter <- matrix(0, nrow = length(all_dge_states), ncol = length(all_ref_states),
                         dimnames = list(all_dge_states, all_ref_states))
  mat_st_jaccard <- mat_st_inter
  
  for (i in all_dge_states) {
    for (j in all_ref_states) {
      inter <- intersect(dge_sets$state[[i]], ref_state[[j]])
      uni <- union(dge_sets$state[[i]], ref_state[[j]])
      mat_st_inter[i, j] <- length(inter)
      mat_st_jaccard[i, j] <- if (length(uni) > 0) length(inter)/length(uni) else 0
    }
  }
  
  right_anno_st <- rowAnnotation(`DGE Size` = anno_barplot(rowSums(mat_st_inter > -1) * 0 + sapply(all_dge_states, function(x) length(dge_sets$state[[x]])), gp = gpar(fill = "#E69F00")), annotation_name_rot = 90)
  top_anno_st <- HeatmapAnnotation(`Ref Size` = anno_barplot(colSums(mat_st_inter > -1) * 0 + sapply(all_ref_states, function(x) length(ref_state[[x]])), gp = gpar(fill = "#56B4E9")))
  
  p2 <- Heatmap(mat_st_jaccard, name = "Jaccard", 
                column_title = paste0("State DGE vs Ref State Overlap"),
                col = colorRamp2(c(0, max(mat_st_jaccard, na.rm = TRUE)), c("white", "darkgreen")),
                top_annotation = top_anno_st,
                right_annotation = right_anno_st,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_side = "left",
                rect_gp = gpar(col = "gray80", lwd = 0.5),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (mat_st_inter[i, j] > 0) {
                    txt_col <- ifelse(mat_st_jaccard[i, j] > 0.25, "white", "black")
                    is_diag <- all_dge_states[i] == all_ref_states[j]
                    if (is_diag) {
                      lbl <- paste0(mat_st_inter[i, j], " / ", length(ref_state[[ all_ref_states[j] ]]))
                    } else {
                      lbl <- mat_st_inter[i, j]
                    }
                    grid.text(lbl, x, y, gp = gpar(fontsize = 8, col = txt_col))
                  }
                })
  draw(p2)
  
  # ==========================================
  # 3. Unified Stacked Bar Plots (Explicit Overlap Styling)
  # ==========================================
  
  build_overlap_df <- function(dge_list, ref_list, ordered_names, feature_type) {
    matched_names <- intersect(ordered_names, intersect(names(dge_list), names(ref_list)))
    summary_list <- list()
    for (m in matched_names) {
      dge_genes <- dge_list[[m]]; ref_genes <- ref_list[[m]]
      inter_len <- length(intersect(dge_genes, ref_genes))
      if ((length(dge_genes) + length(ref_genes)) > 0) {
        summary_list[[m]] <- data.frame(
          Identity = make_feature_label(m, feature_type), # Apply descriptions here
          Both = inter_len,
          DGE_Only = length(setdiff(dge_genes, ref_genes)),
          Ref_Only = length(setdiff(ref_genes, dge_genes))
        )
      }
    }
    if (length(summary_list) == 0) return(NULL)
    
    df <- bind_rows(summary_list) %>%
      pivot_longer(cols = c("Both", "DGE_Only", "Ref_Only"), names_to = "Category", values_to = "Count")
    
    # Map identity back to ordered factor based on descriptions
    ordered_labels <- make_feature_label(matched_names, feature_type)
    df$Identity <- factor(df$Identity, levels = rev(ordered_labels))
    
    # Rename categories to make the overlap obvious
    df$Category <- factor(df$Category, 
                          levels = c("Ref_Only", "Both", "DGE_Only"),
                          labels = c("Unique to Reference", "Shared (Overlap)", "Unique to DGE"))
    return(df)
  }
  
  # Plot MPs
  df_mp_overlap <- build_overlap_df(dge_sets$mp, ref_mp, retained_mps, "MP")
  if (!is.null(df_mp_overlap)) {
    p3 <- ggplot(df_mp_overlap, aes(x = Identity, y = Count, fill = Category)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.75) + # Added outlines to define the blocks clearly
      geom_text(aes(label = ifelse(Count > 0, Count, "")),       # Add numbers directly inside the bars
                position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
      scale_fill_manual(values = c("Unique to Reference" = "#D55E00", 
                                   "Shared (Overlap)" = "#009E73", 
                                   "Unique to DGE" = "#CC79A7")) +
      coord_flip() + theme_classic() +
      labs(title = paste0("MP Gene Overlap"), 
           x = "Metaprogram", y = "Number of Genes", fill = "") +
      theme(legend.position = "bottom", 
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, face = "italic", color = "gray30"))
    print(p3)
  }
  
  # Plot States
  df_st_overlap <- build_overlap_df(dge_sets$state, ref_state, expected_state_order, "State")
  if (!is.null(df_st_overlap)) {
    p4 <- ggplot(df_st_overlap, aes(x = Identity, y = Count, fill = Category)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.75) +
      geom_text(aes(label = ifelse(Count > 0, Count, "")), 
                position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
      scale_fill_manual(values = c("Unique to Reference" = "#D55E00", 
                                   "Shared (Overlap)" = "#009E73", 
                                   "Unique to DGE" = "#CC79A7")) +
      coord_flip() + theme_classic() +
      labs(title = paste0("State Gene Overlap"), 
           x = "State", y = "Number of Genes", fill = "") +
      theme(legend.position = "bottom", 
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, face = "italic", color = "gray30"))
    print(p4)
  }
}
dev.off()

all_res <- list()
split_methods <- c("continuous", "median", "q1q4")
volcano_pages_mp <- list()
volcano_pages_state <- list()

for (sm in split_methods) {
  message("Running survival analysis with split: ", sm)
  
  for (mode_name in requested_modes) {
    message("  Mode: ", mode_name)
    dge_sets <- dge_by_mode[[mode_name]]

    method_inputs <- list(
      malignant_only = tpm_mal,
      whole_tcga = tpm_whole,
      dge_based_malignant = tpm_mal,
      dge_based_whole = tpm_whole
    )

    method_order <- c("malignant_only", "whole_tcga", "dge_based_malignant", "dge_based_whole")
    panel_results <- list()

    for (method_name in method_order) {
      expr_mat <- method_inputs[[method_name]]
      use_dge <- grepl("^dge_based", method_name)
      mp_sets <- if (use_dge) dge_sets$mp else c(mp.genes, pan_mp_sets)
      state_sets <- if (use_dge) dge_sets$state else {
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

      cox_mp <- run_cox(merged_df, mp_cols, mode_name, method_name, "MP", split_method = sm)
      cox_st <- run_cox(merged_df, st_cols, mode_name, method_name, "State", split_method = sm)
      cox_all <- bind_rows(cox_mp, cox_st)
      all_res[[paste(sm, mode_name, method_name, sep = "|:")]] <- cox_all

      this_mp <- cox_mp %>% filter(cohort == "EAC")
      if (nrow(this_mp) > 0) {
        mp_levels <- c(make_feature_label(retained_mps, "MP"), make_feature_label(names(pan_mp_sets), "MP"))
        this_mp$feature <- make_feature_label(this_mp$feature, "MP")
        all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
        this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
        panel_results[[paste0("MP|", method_name)]] <- plot_volcano(this_mp, paste0("[", sm, "] ", mode_name, " ", method_name, " MP volcano"))
      } else {
        panel_results[[paste0("MP|", method_name)]] <- NULL
      }
      this_st <- cox_st %>% filter(cohort == "EAC")
      if (nrow(this_st) > 0) {
        this_st <- this_st %>% filter(!feature %in% c("Unresolved", "Hybrid")) %>%
          mutate(feature = factor(feature, levels = ordered_states_for_plot(unique(as.character(feature)))))
        panel_results[[paste0("State|", method_name)]] <- plot_volcano(this_st, paste0("[", sm, "] ", mode_name, " ", method_name, " State volcano"))
      } else {
        panel_results[[paste0("State|", method_name)]] <- NULL
      }
    }

    volcano_pages_mp[[paste(mode_name, sm, sep = "|")]] <- make_tcga_page(
      panel_results[["MP|malignant_only"]],
      panel_results[["MP|dge_based_malignant"]],
      panel_results[["MP|whole_tcga"]],
      panel_results[["MP|dge_based_whole"]],
      paste0(mode_name, " MP volcano: ", sm)
    )
    volcano_pages_state[[paste(mode_name, sm, sep = "|")]] <- make_tcga_page(
      panel_results[["State|malignant_only"]],
      panel_results[["State|dge_based_malignant"]],
      panel_results[["State|whole_tcga"]],
      panel_results[["State|dge_based_whole"]],
      paste0(mode_name, " State volcano: ", sm)
    )
  }
}

unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_continuous.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_median.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_q1q4.pdf")), force = TRUE)
unlink(file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_volcano_methods_noreg.pdf")), force = TRUE)

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
  cox_res$padj <- ave(cox_res$P_value, interaction(cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type, cox_res$split_method),
                      FUN = function(x) p.adjust(x, method = "BH"))
}
write.csv(cox_res, file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_state_cox_methods_splits.csv")), row.names = FALSE)

summary_dir <- file.path(
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
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
write.csv(summary_df, file.path(summary_dir, "Auto_survival_clinical_mps_v2_reg_noreg_summary.csv"), row.names = FALSE)

message("Saved 4-method MP/state continuous-Cox volcano outputs (noreg).")
