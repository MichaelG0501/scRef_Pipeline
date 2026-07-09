####################
# Analysis registry:
#   Status: active
#   Script: analysis/publication/poster_section7_survival_targeting.R
#   Inputs:
#     ref_outs/task2_filtered_survival/Auto_task2_filtered_survival_mp_state_cox_methods_splits.csv
#     ref_outs/Auto_drug_reversal/overlap_visuals/top_3_overlapped_drugs_final_methods.csv
#     ref_outs/Auto_drug_reversal/method_visuals/Auto_drug_reversal_final_asgard_cmap_top100_membership.csv
#   Outputs:
#     ref_outs/publication/section7/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section7_survival_targeting.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(data.table)
  library(survival)
  library(ggrepel)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section7"
out_dir <- pub_section_dir(section)
poor_states <- c("Classic Proliferative", "Stress-adaptive")
lineage_survival_states <- c("Basal to Intestinal Metaplasia", "SMG-like Metaplasia")
pdo_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

####################
# FIGURE 1: TCGA Kaplan-Meier curves for adverse states.
####################
cat("Generating TCGA KM survival curves...\n")
make_km_curve <- function(states, output_name) {
  expr_file <- file.path(SCREF_REF_OUTS_DIR, "cibersortx", "TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
  meta_file <- file.path(SCREF_REF_OUTS_DIR, "tcga_esca_meta.rds")
  cox_file <- file.path(SCREF_REF_OUTS_DIR, "task2_survival", "Auto_task2_survival_mp_state_cox_methods_splits.csv")
  if (!all(file.exists(expr_file, meta_file, SCREF_MP_OBJECT_RDS, cox_file))) return(NULL)
  expr_dt <- fread(expr_file)
  genes <- expr_dt[[1]]
  expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
  rownames(expr_mat) <- genes
  expr_mat <- log2(expr_mat + 1)
  meta_tcga <- readRDS(meta_file) |>
    filter(sample_type_code == "01", HistologyGroup == "EAC", !is.na(OS_time), !is.na(OS_event))
  geneNMF.metaprograms <- readRDS(SCREF_MP_OBJECT_RDS)
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  state_sets <- lapply(SCREF_STATE_GROUPS[states], function(mps) {
    unique(unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE))
  })
  state_sets <- state_sets[sapply(state_sets, length) >= 5]
  gs <- GSVA::gsva(expr_mat, state_sets, method = "gsva", kcdf = "Gaussian")
  gs_df <- as.data.frame(t(gs)) |> tibble::rownames_to_column("sample_barcode")
  dd <- meta_tcga |> left_join(gs_df, by = "sample_barcode")
  cox <- read_csv(cox_file, show_col_types = FALSE) |>
    filter(method == "whole_tcga", split_method == "continuous",
           feature_type == "State", feature %in% states)
           
  # Calculate N for High/Low using the first state to ensure consistent legend labels
  d_first <- dd |> filter(!is.na(.data[[states[1]]]))
  q_first <- quantile(d_first[[states[1]]], probs = c(0.25, 0.75), na.rm = TRUE)
  n_high <- sum(d_first[[states[1]]] >= q_first[2], na.rm = TRUE)
  n_low <- sum(d_first[[states[1]]] <= q_first[1], na.rm = TRUE)
  label_high <- sprintf("Q4 High (n=%d)", n_high)
  label_low <- sprintf("Q1 Low (n=%d)", n_low)

  km_long <- bind_rows(lapply(states, function(st) {
    d <- dd |> filter(!is.na(.data[[st]]))
    q_st <- quantile(d[[st]], probs = c(0.25, 0.75), na.rm = TRUE)
    d_q <- d |> filter(.data[[st]] <= q_st[1] | .data[[st]] >= q_st[2])
    d_q$group <- factor(ifelse(d_q[[st]] >= q_st[2], label_high, label_low), levels = c(label_low, label_high))
    fit <- survfit(Surv(OS_time / 30.44, OS_event) ~ group, data = d_q)
    ss <- summary(fit)
    tibble(time = ss$time, surv = ss$surv, strata = sub("^group=", "", ss$strata), state = st)
  }))
  cox_lab <- cox |>
    transmute(state = feature, label = sprintf("Cox HR %.2f\nP = %.3f", HR, P_value))
    
  strata_colors <- setNames(c("#4B5563", "#B2182B"), c(label_low, label_high))
  
  p <- ggplot(km_long, aes(time, surv, colour = strata)) +
    geom_step(linewidth = 0.82) +
    facet_wrap(~ state, nrow = 1) +
    geom_text(data = cox_lab, aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3.0, fontface = "bold") +
    scale_colour_manual(values = strata_colors, name = "State score") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = "Overall survival (months)", y = "Survival probability") +
    pub_theme(12) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 11))
  save_pub_gg(p, section, output_name, width = 7.5, height = 3.9)
}
if (is.null(tryCatch(make_km_curve(poor_states, "s7_km_survival"), error = function(e) NULL))) {
  abort_missing_figure(section, "s7_km_survival", "TCGA Kaplan-Meier survival",
                   "Required TCGA expression, metadata, MP object, or Cox table missing.", width = 7.5, height = 3.9)
}
if (is.null(tryCatch(make_km_curve(lineage_survival_states, "s5_lineage_km_survival"), error = function(e) NULL))) {
  abort_missing_figure(section, "s5_lineage_km_survival", "Lineage-state TCGA Kaplan-Meier survival",
                   "Required TCGA expression, metadata, MP object, or Cox table missing.", width = 7.5, height = 3.9)
}

####################
# FIGURE 2: Top state-specific markers for future ASO/CRISPR targeting.
####################
marker_file <- file.path(pdo_dir, "Auto_five_state_markers", "Auto_five_state_markers_final.csv")
global_file <- file.path(pdo_dir, "Auto_five_state_markers", "Auto_five_state_global_marker_screen.csv.gz")

if (file.exists(marker_file) && file.exists(global_file)) {
  top_markers <- read_csv(marker_file, show_col_types = FALSE) |>
    mutate(state = clean_state(state)) |>
    filter(state %in% poor_states, best_state_match) |>
    group_by(state) |>
    slice_max(ranking_score, n = 3, with_ties = FALSE) |>
    ungroup() |>
    pull(gene) |>
    unique()
    
  plot_data <- read_csv(global_file, show_col_types = FALSE) |>
    mutate(state = clean_state(state)) |>
    filter(gene %in% top_markers) |>
    group_by(gene) |>
    mutate(scaled_expr = as.numeric(scale(global_mean_state))) |>
    ungroup() |>
    filter(state %in% poor_states) |>
    mutate(state = factor(state, levels = poor_states),
           gene = factor(gene, levels = rev(top_markers)))

  p <- ggplot(plot_data, aes(state, gene)) +
    geom_point(aes(size = global_pct_state, fill = scaled_expr),
               shape = 21, colour = "#111111", stroke = 0.30) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = "Scaled\nexpression") +
    scale_size_continuous(range = c(3, 8), labels = scales::percent_format(accuracy = 1),
                          name = "State cells\nexpressing") +
    labs(x = NULL, y = NULL) +
    pub_theme(12) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1, face = "bold"))
    
  save_pub_gg(p, section, "s7_marker_specificity_targets", width = 4.8, height = 3.9)
  write_csv(plot_data, file.path(out_dir, "tables", "s7_marker_specificity_targets.csv"))
}

## ==============================================================================
# FIGURE 2: scAtlas and PDO L1000 signature reversal profiles
# ==============================================================================
cat("Generating L1000 reversal profiles dynamically for overlapped drugs...\n")

top_overlap_file <- file.path(SCREF_REF_OUTS_DIR, "Auto_drug_reversal", "overlap_visuals", "top_3_overlapped_drugs_per_state.csv")
if (file.exists(top_overlap_file)) {
  top_overlap <- read_csv(top_overlap_file, show_col_types = FALSE)
  
  load_l1000_reference <- function() {
    default_ref_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference"
    rank_dt <- fread(file.path(default_ref_dir, "stomach_rankMatrix.txt"))
    gene_dt <- fread(file.path(default_ref_dir, "stomach_gene_info.txt"))
    drug_dt <- fread(file.path(default_ref_dir, "stomach_drug_info.txt"))
    setnames(gene_dt, old = names(gene_dt)[1:2], new = c("probe_id", "gene_symbol"))
    instance_cols <- setdiff(colnames(rank_dt), "probe_id")
    rank_dt <- merge(rank_dt, gene_dt, by = "probe_id", all.x = FALSE, all.y = FALSE)
    rank_dt <- rank_dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
    rank_gene <- rank_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol, .SDcols = instance_cols]
    rank_mat <- as.matrix(rank_gene[, ..instance_cols])
    rownames(rank_mat) <- rank_gene$gene_symbol
    norm_mat <- (rank_mat - 1) / pmax(nrow(rank_mat) - 1, 1)

    instance_map <- data.table(instance_id = as.character(drug_dt$instance_id))
    if ("cmap_name" %in% names(drug_dt) && "catalog_name" %in% names(drug_dt)) {
      instance_map[, drug := mapply(
        FUN = function(cmap_name, catalog_name) sub(paste0("_", catalog_name, "$"), "", cmap_name),
        as.character(drug_dt$cmap_name),
        as.character(drug_dt$catalog_name),
        USE.NAMES = FALSE
      )]
    } else if ("cmap_name" %in% names(drug_dt)) {
      instance_map[, drug := sub("_[^_]+$", "", as.character(drug_dt$cmap_name))]
    } else {
      instance_map[, drug := as.character(instance_id)]
    }
    instance_map[, drug_key := tolower(trimws(drug))]
    list(norm_mat = norm_mat, instance_map = instance_map)
  }

  get_drug_gene_profile <- function(norm_mat, instance_map, drug_key, genes) {
    instance_ids <- instance_map[["instance_id"]][instance_map[["drug_key"]] == drug_key]
    instance_ids <- intersect(instance_ids, colnames(norm_mat))
    genes <- intersect(genes, rownames(norm_mat))
    if (length(instance_ids) == 0 || length(genes) == 0) return(NULL)
    out <- rowMeans(norm_mat[genes, instance_ids, drop = FALSE], na.rm = TRUE)
    tibble(gene = names(out), l1000_rank = as.numeric(out), n_l1000_instances = length(instance_ids))
  }

  l1000 <- tryCatch(load_l1000_reference(), error = function(e) NULL)
  
  if (!is.null(l1000)) {
    build_profile <- function(dataset_name, base_dir) {
      degs <- fread(file.path(base_dir, "Auto_drug_reversal", "Auto_drug_reversal_degs_all_states.csv.gz")) |>
        filter(state %in% poor_states) |> mutate(avg_logFC = avg_log2FC)
      sig_dt <- fread(file.path(base_dir, "Auto_drug_reversal", "Auto_drug_reversal_signature_top150.csv")) |>
        filter(state %in% poor_states)
        
      bind_rows(lapply(seq_len(nrow(top_overlap)), function(i) {
        state_name <- top_overlap$state[i]
        if (!state_name %in% poor_states) return(NULL)
        drug_key_i <- top_overlap$drug_key[i]
        genes <- sig_dt |> filter(state == state_name) |> pull(gene) |> unique()
        prof <- get_drug_gene_profile(l1000$norm_mat, l1000$instance_map, drug_key_i, genes)
        if (is.null(prof)) return(NULL)
        prof |>
          left_join(degs |> filter(state == state_name) |> dplyr::select(gene, avg_logFC), by = "gene") |>
          mutate(state = state_name, drug = top_overlap$drug[i], drug_key = drug_key_i, dataset = dataset_name,
                 direction = ifelse(avg_logFC >= 0, "State-up genes", "State-down genes"))
      }))
    }
    
    prof_sc <- build_profile("scAtlas", SCREF_REF_OUTS_DIR)
    
    if (nrow(prof_sc) > 0) {
      prof_plot <- prof_sc |>
        group_by(state, drug, drug_key, direction) |>
        summarise(mean_l1000_rank = mean(l1000_rank, na.rm = TRUE),
                  n_l1000_instances = first(n_l1000_instances), .groups = "drop") |>
        mutate(state = factor(state, levels = poor_states),
               direction = factor(direction, levels = c("State-up genes", "State-down genes")),
               drug_label = ifelse(nchar(drug) > 25, paste0(substr(drug, 1, 22), "..."), drug))

      p <- ggplot(prof_plot, aes(x = direction, y = mean_l1000_rank)) +
        geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40", linewidth = 0.6) +
        geom_line(aes(group = drug_key), color = "grey40", alpha = 0.4, linewidth = 0.5) +
        geom_point(aes(color = state, size = n_l1000_instances), alpha = 0.8) +
        geom_text_repel(
          aes(label = drug_label),
          size = 4,
          max.overlaps = 50,
          min.segment.length = 0,
          box.padding = 0.4,
          segment.alpha = 0.4
        ) +
        scale_color_manual(values = PUB_STATE_COLOURS) +
        scale_size_continuous(range = c(2, 6)) +
        scale_y_continuous(limits = c(0.1, 0.9), breaks = seq(0, 1, 0.2)) +
        facet_wrap(~ state, ncol = 2) +
        labs(x = NULL, y = "Mean normalized L1000 rank") +
        pub_theme(11) +
        theme(strip.text = element_text(face = "bold", size = 11),
              legend.position = "none", axis.title.y = element_text(face = "bold"))
              
      save_pub_gg(p, section, "s7_l1000_reversal_profiles", width = 8, height = 7)
      write_csv(prof_plot, file.path(out_dir, "tables", "s7_l1000_reversal_profiles.csv"))
    }
  }
}

cat("Section 7 complete.\n")
