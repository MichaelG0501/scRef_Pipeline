####################
# Auto_drug_reversal_method_visuals.R
#
# Unified visualization script for drug reversal analysis.
# Combines method-level screening, final consensus overlap, and predicted reversion.
# Supersedes: method_visuals, predicted_reversion, and consensus_visuals scripts.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
  library(grid)
  library(patchwork)
  library(futile.logger)
})

flog.threshold(ERROR) # silence VennDiagram logs

####################
# setup
####################

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

base_dir <- "Auto_drug_reversal"
out_dir <- file.path(base_dir, "method_visuals")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#A65628"
)

method_order <- c("ASGARD", "scDrugPrio", "CLUE_FALLBACK_LOCAL")
method_labels <- c(
  "ASGARD" = "ASGARD",
  "scDrugPrio" = "scDrugPrio",
  "CLUE_FALLBACK_LOCAL" = "Local CMap/L1000"
)

primary_method_order <- c("ASGARD", "CLUE_FALLBACK_LOCAL")

####################
# helpers
####################

safe_state_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

state_safe_map <- setNames(state_order, safe_state_name(state_order))

write_status <- function(status, detail) {
  fwrite(
    data.frame(
      step = "drug_reversal_unified_visuals",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_drug_reversal_unified_visual_status.csv")
  )
}

standardize_rank_file <- function(path, method_name) {
  if (!file.exists(path)) return(tibble())
  df <- fread(path)
  if (nrow(df) == 0) return(tibble())

  nms <- names(df)
  pick <- function(candidates) {
    hit <- intersect(candidates, nms)[1]
    if (is.na(hit)) {
      hit <- grep(paste(candidates, collapse = "|"), nms, ignore.case = TRUE, value = TRUE)[1]
    }
    hit
  }

  state_col <- pick(c("state", "State", "state_query"))
  drug_col <- pick(c("drug", "Drug", "pert_iname", "compound", "drug_name", "name"))
  rank_col <- pick(c("rank", "Rank"))
  score_col <- pick(c("score", "drug_score", "connectivity_score", "rank_metric"))
  p_col <- pick(c("p_value", "pvalue", "P.Value", "P"))
  fdr_col <- pick(c("fdr", "FDR", "adj.P.Val"))
  target_col <- pick(c("target_genes", "targets", "target", "counteracting_targets"))
  moa_col <- pick(c("moa", "MOA", "mechanism"))

  if (is.na(state_col) || is.na(drug_col)) return(tibble())

  out <- tibble(
    state = as.character(df[[state_col]]),
    method = method_name,
    drug = as.character(df[[drug_col]]),
    rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(df[[rank_col]])) else NA_real_,
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    fdr = if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]])) else NA_real_,
    target_genes = if (!is.na(target_col)) as.character(df[[target_col]]) else NA_character_,
    moa = if (!is.na(moa_col)) as.character(df[[moa_col]]) else NA_character_
  ) %>%
    mutate(
      state = ifelse(state %in% names(state_safe_map), state_safe_map[state], state),
      state = factor(state, levels = state_order),
      method = factor(method, levels = method_order),
      drug_key = tolower(trimws(drug))
    ) %>%
    filter(!is.na(state), !is.na(drug_key), nzchar(drug_key))

  if (all(is.na(out$rank))) {
    out <- out %>%
      arrange(state, desc(score), drug) %>%
      group_by(state) %>%
      mutate(rank = row_number()) %>%
      ungroup()
  }

  out %>%
    group_by(state, method, drug_key) %>%
    arrange(rank, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
}

parse_targets <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(character())
  y <- unlist(strsplit(x, "[;,|/]"))
  y <- gsub("\\s*\\([^)]*\\)", "", y)
  y <- gsub("[^A-Za-z0-9_.-]", "", y)
  unique(y[nzchar(y)])
}

load_drugbank_targets <- function(drug_keys) {
  drugbank_path <- Sys.getenv(
    "AUTO_DRUGBANK_TARGETS",
    "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/all_drug_targets_drug_bank.txt"
  )
  if (!file.exists(drugbank_path)) return(tibble())

  db <- fread(drugbank_path)
  db <- db %>%
    transmute(
      drug_key = tolower(trimws(drug_name)),
      target_gene = as.character(gene_symbol),
      action = as.character(target_organism),
      organism = as.character(drug_action)
    ) %>%
    filter(
      drug_key %in% drug_keys,
      organism == "Humans",
      !is.na(target_gene), nzchar(target_gene)
    ) %>%
    distinct()
  return(db)
}

shorten_drug <- function(x, n = 34) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

score_scale <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  out <- as.numeric(scale(x))
  out[is.na(out)] <- 0
  pmax(pmin(out, 3), -3)
}

load_l1000_reference <- function() {
  ref_cfg_path <- file.path(base_dir, "asgard_reference", "Auto_asgard_reference_paths.csv")
  default_ref_dir <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/DrugReference"
  if (file.exists(ref_cfg_path)) {
    ref_cfg <- fread(ref_cfg_path)
    rank_path <- ref_cfg$AUTO_ASGARD_DRUG_RESPONSE[1]
    gene_path <- ref_cfg$AUTO_ASGARD_GENE_INFO[1]
    drug_path <- ref_cfg$AUTO_ASGARD_DRUG_INFO[1]
  } else {
    rank_path <- file.path(default_ref_dir, "stomach_rankMatrix.txt")
    gene_path <- file.path(default_ref_dir, "stomach_gene_info.txt")
    drug_path <- file.path(default_ref_dir, "stomach_drug_info.txt")
  }
  
  rank_dt <- fread(rank_path)
  gene_dt <- fread(gene_path)
  drug_dt <- fread(drug_path)
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
  tibble(
    gene = names(out),
    l1000_rank = as.numeric(out),
    l1000_centered = (as.numeric(out) - 0.5) * 2,
    n_l1000_instances = length(instance_ids)
  )
}

get_assay_layer <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

load_state5 <- function() {
  cache <- file.path(base_dir, "cache", "Auto_drug_reversal_scRef_states.rds")
  if (file.exists(cache)) {
    obj <- readRDS(cache)
    obj <- subset(obj, subset = state %in% state_order)
    obj$state <- factor(as.character(obj$state), levels = state_order)
    return(obj)
  }

  sc_ref_all <- readRDS("EAC_Ref_epi.rds")
  state_labels <- readRDS("Auto_final_states.rds")
  DefaultAssay(sc_ref_all) <- "RNA"
  common_cells <- intersect(colnames(sc_ref_all), names(state_labels))
  keep_cells <- common_cells[state_labels[common_cells] %in% state_order]
  obj <- subset(sc_ref_all, cells = keep_cells)
  obj$state <- factor(as.character(state_labels[colnames(obj)]), levels = state_order)
  obj
}

####################
# load core data
####################

deg_path <- file.path(base_dir, "Auto_drug_reversal_degs_all_states.csv.gz")
sig_path <- file.path(base_dir, "Auto_drug_reversal_signature_top150.csv")

degs <- tibble()
signature_dt <- tibble()
if (file.exists(deg_path) && file.exists(sig_path)) {
  degs <- fread(deg_path) %>%
    mutate(state = factor(as.character(state), levels = state_order)) %>%
    filter(!is.na(state))
  logfc_col <- intersect(c("avg_logFC", "avg_log2FC", "logFC", "log2FC"), names(degs))[1]
  if (!is.na(logfc_col)) {
    degs <- degs %>% mutate(avg_logFC = .data[[logfc_col]])
  }
  signature_dt <- fread(sig_path) %>%
    mutate(state = factor(as.character(state), levels = state_order)) %>%
    filter(!is.na(state))
}

####################
# load rankings
####################

rankings <- bind_rows(
  standardize_rank_file(file.path(base_dir, "asgard", "Auto_asgard_ranked_drugs.csv"), "ASGARD"),
  standardize_rank_file(file.path(base_dir, "scdrugprio", "Auto_scdrugprio_ranked_drugs.csv"), "scDrugPrio"),
  standardize_rank_file(file.path(base_dir, "clue_fallback", "Auto_clue_ranked_drugs.csv"), "CLUE_FALLBACK_LOCAL")
)

if (nrow(rankings) == 0) {
  write_status("missing_rankings", "No method ranking files were available.")
  quit(save = "no", status = 0)
}

rankings <- rankings %>%
  group_by(method, state) %>%
  mutate(
    rank_percentile = percent_rank(rank),
    score_z = score_scale(score)
  ) %>%
  ungroup()

####################
# overlap logic
####################

top100 <- rankings %>%
  filter(rank <= 100) %>%
  mutate(method_chr = as.character(method)) %>%
  distinct(state, drug_key, drug, method_chr, rank)

membership <- top100 %>%
  group_by(state, drug_key) %>%
  summarise(
    drug = drug[which.min(rank)],
    methods = paste(sort(unique(method_chr)), collapse = " + "),
    n_methods = n_distinct(method_chr),
    best_rank = min(rank, na.rm = TRUE),
    mean_rank = mean(rank, na.rm = TRUE),
    .groups = "drop"
  )

# REFINED grouping: assign drug to its BEST consensus state (min mean rank)
drug_best_state <- membership %>%
  filter(n_methods >= 2) %>%
  group_by(drug_key) %>%
  arrange(mean_rank, match(as.character(state), state_order)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(drug_key, target_state = state)

membership <- membership %>%
  left_join(drug_best_state, by = "drug_key") %>%
  mutate(
    consensus_class = case_when(
      n_methods == 3 ~ "All three methods",
      n_methods == 2 ~ "Any two methods",
      TRUE ~ "Single method only"
    )
  )

fwrite(membership, file.path(out_dir, "Auto_drug_reversal_three_method_top100_membership.csv"))

primary_membership <- top100 %>%
  filter(method_chr %in% primary_method_order) %>%
  group_by(state, drug_key) %>%
  summarise(
    drug = drug[which.min(rank)],
    methods = paste(sort(unique(method_chr)), collapse = " + "),
    n_methods = n_distinct(method_chr),
    best_rank = min(rank, na.rm = TRUE),
    mean_rank = mean(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_methods == length(primary_method_order)) %>%
  mutate(
    target_state = state,
    consensus_class = "ASGARD + Local CMap/L1000",
    final_selection_rule = "ASGARD_and_Local_CMap_only_scDrugPrio_visual_only"
  )

fwrite(primary_membership, file.path(out_dir, "Auto_drug_reversal_final_asgard_cmap_top100_membership.csv"))

####################
# basic visuals (Method Heatmap, Overlap, Venns)
####################

# Per-method top-12
top_method <- rankings %>%
  group_by(method, state) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(
    state = factor(as.character(state), levels = state_order),
    method = factor(as.character(method), levels = method_order),
    drug_label = shorten_drug(drug),
    drug_axis = factor(paste(method, state, sprintf("%03d", rank), drug_label, sep = " | "),
                       levels = rev(paste(method, state, sprintf("%03d", rank), drug_label, sep = " | ")))
  )

p_method_heat <- ggplot(top_method, aes(x = state, y = drug_axis)) +
  geom_tile(aes(fill = pmin(rank, 100)), color = "white", linewidth = 0.25) +
  geom_text(aes(label = rank), size = 2.5) +
  scale_fill_gradient(low = "#CB181D", high = "#FEE0D2", name = "Rank") +
  scale_y_discrete(labels = function(x) sub("^([^|]+\\| ){3}", "", x)) +
  facet_wrap(~ method, scales = "free_y", ncol = 1, labeller = as_labeller(method_labels)) +
  labs(x = NULL, y = NULL, title = "Top predicted inhibitors by method and scRef state") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 35, hjust = 1, color = state_cols[levels(top_method$state)]),
    strip.text = element_text(face = "bold", size = 11)
  )
ggsave(file.path(out_dir, "Auto_drug_reversal_method_rank_heatmap.pdf"), p_method_heat, width = 13, height = 14)

# Overlap barplot
overlap_summary <- membership %>%
  count(state, methods, n_methods, consensus_class, name = "n_drugs") %>%
  arrange(state, desc(n_methods), desc(n_drugs))

p_overlap <- ggplot(overlap_summary, aes(x = methods, y = n_drugs, fill = consensus_class)) +
  geom_col() +
  geom_text(aes(label = n_drugs), vjust = -0.25, size = 2.8) +
  scale_fill_manual(values = c("All three methods" = "#CB181D", "Any two methods" = "#FD8D3C", "Single method only" = "#BDBDBD")) +
  facet_wrap(~ state, scales = "free_x", ncol = 1) +
  labs(x = NULL, y = "Top-100 drugs", fill = NULL, title = "Top-100 drug overlap across screening methods") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
ggsave(file.path(out_dir, "Auto_drug_reversal_overlap_summary_barplot.pdf"), p_overlap, width = 13, height = 11)

# Venns
venn_plots <- list()
for (st in state_order) {
  a <- top100 %>% filter(state == st, method_chr == "ASGARD") %>% pull(drug_key)
  s <- top100 %>% filter(state == st, method_chr == "scDrugPrio") %>% pull(drug_key)
  c_l <- top100 %>% filter(state == st, method_chr == "CLUE_FALLBACK_LOCAL") %>% pull(drug_key)

  v3 <- venn.diagram(x = list(ASGARD=a, scDrugPrio=s, CMap=c_l), category.names = c("ASGARD", "scDrugPrio", "CMap"),
                     fill = c("#6BAED6", "#FDAE6B", "#A1D99B"), alpha = 0.5, filename = NULL, disable.logging = TRUE)
  p3 <- wrap_elements(gTree(children = v3)) + ggtitle(paste(st, "- All 3 Methods")) + theme(plot.title=element_text(size=10, face="bold", hjust=0.5))

  venn_plots[[st]] <- p3
}
ggsave(file.path(out_dir, "Auto_drug_reversal_overlap_venns.pdf"), wrap_plots(venn_plots, ncol = 3), width = 16, height = 11)

####################
# consensus visuals (Rank Matrix, Targets, Rank-Rank, Dotplot)
####################

final_candidates <- primary_membership %>%
  group_by(state) %>%
  arrange(mean_rank, best_rank, drug, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

final_rank_matrix <- final_candidates %>%
  inner_join(rankings %>% select(state, drug_key, method, rank), by = c("state", "drug_key")) %>%
  mutate(
    drug_label = shorten_drug(drug),
    y_axis = factor(paste(state, sprintf("%02d", dense_rank(mean_rank)), drug_label, sep = " | "),
                    levels = rev(unique(paste(state, sprintf("%02d", dense_rank(mean_rank)), drug_label, sep = " | ")))),
    method = factor(as.character(method), levels = method_order)
  )

p_final_rank <- ggplot(final_rank_matrix, aes(x = method, y = y_axis)) +
  geom_tile(aes(fill = pmin(rank, 100)), color = "white", linewidth = 0.25) +
  geom_text(aes(label = rank), size = 2.5) +
  scale_fill_gradient(low = "#CB181D", high = "#FEE0D2", name = "Rank") +
  scale_x_discrete(labels = method_labels, drop = FALSE) +
  scale_y_discrete(labels = function(x) sub("^[^|]+ \\| [^|]+ \\| ", "", x)) +
  facet_wrap(~ state, scales = "free_y", ncol = 1) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0), axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "Auto_drug_reversal_final_overlap_rank_matrix.pdf"), p_final_rank, width = 9, height = 12)

# Rank-Rank
rank_rank <- rankings %>%
  filter(method %in% primary_method_order) %>%
  pivot_wider(id_cols = c(state, drug_key, drug), names_from = method, values_from = rank) %>%
  filter(!is.na(ASGARD), !is.na(CLUE_FALLBACK_LOCAL)) %>%
  mutate(consensus_top100 = ASGARD <= 100 & CLUE_FALLBACK_LOCAL <= 100)

if (nrow(rank_rank) == 0) {
  rank_rank <- tibble()
} else {
  rank_rank <- rank_rank %>% rename(CMap = CLUE_FALLBACK_LOCAL)
}

if (nrow(rank_rank) > 0) {
  top5_rank_keys <- primary_membership %>%
    group_by(state) %>%
    arrange(mean_rank, best_rank, drug, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    transmute(state_drug_key = paste(state, drug_key, sep = "||")) %>%
    pull(state_drug_key)

  p_rank <- ggplot(rank_rank, aes(x = ASGARD, y = CMap)) +
    geom_hline(yintercept = 100, linewidth = 0.25, linetype = "dashed", color = "grey45") +
    geom_vline(xintercept = 100, linewidth = 0.25, linetype = "dashed", color = "grey45") +
    geom_point(aes(color = consensus_top100), alpha = 0.55, size = 1.5) +
    geom_text_repel(
      data = rank_rank %>% filter(paste(state, drug_key, sep = "||") %in% top5_rank_keys),
      aes(label = shorten_drug(drug, 22)),
      size = 2.6, max.overlaps = 50, min.segment.length = 0, box.padding = 0.25
    ) +
    scale_x_reverse() + scale_y_reverse() +
    scale_color_manual(values = c("FALSE" = "grey72", "TRUE" = "#CB181D"), guide = "none") +
    facet_wrap(~ state, scales = "free", ncol = 2) +
    labs(x = "ASGARD rank (lower is better)", y = "Local CMap/L1000 rank (lower is better)",
         title = "Final consensus rank comparison: ASGARD versus Local CMap/L1000") +
    theme_classic(base_size = 11) + theme(strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"))
  ggsave(file.path(out_dir, "Auto_drug_reversal_rank_rank_scatter.pdf"), p_rank, width = 12, height = 10)
}

# Target Dotplot
consensus_top5 <- primary_membership %>%
  group_by(state) %>%
  arrange(mean_rank, best_rank, drug, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

# Extract targets from rankings and DrugBank
all_targets <- rankings %>%
  filter(drug_key %in% consensus_top5$drug_key) %>%
  group_by(drug_key) %>%
  summarise(target_genes = paste(unique(target_genes), collapse = ";"), .groups = "drop") %>%
  rowwise() %>%
  mutate(target_gene = list(parse_targets(target_genes))) %>%
  unnest(target_gene) %>%
  distinct(drug_key, target_gene)

db_targets <- load_drugbank_targets(unique(consensus_top5$drug_key))
all_targets <- bind_rows(all_targets, db_targets %>% select(drug_key, target_gene)) %>% distinct()

if (nrow(all_targets) > 0) {
  obj <- load_state5()
  expr <- get_assay_layer(obj)
  all_targets <- all_targets %>% filter(target_gene %in% rownames(expr))
  
  if (nrow(all_targets) > 0) {
    cells_by_state <- split(colnames(obj), as.character(obj$state))
    metrics <- bind_rows(lapply(names(cells_by_state), function(st) {
      m <- expr[unique(all_targets$target_gene), cells_by_state[[st]], drop=FALSE]
      tibble(target_gene = rownames(m), expression_state = st, mean_expr = rowMeans(m), pct_expr = rowMeans(m > 0))
    }))
    
    target_expr <- consensus_top5 %>%
      select(state, drug, drug_key, mean_rank, target_state) %>%
      left_join(all_targets, by = "drug_key", relationship = "many-to-many") %>%
      left_join(metrics, by = "target_gene", relationship = "many-to-many") %>%
      filter(!is.na(target_gene)) %>%
      mutate(target_state_flag = (as.character(expression_state) == as.character(state)),
             expression_state = factor(expression_state, levels = state_order))
    
    # Sort drugs by target state factor for the dotplot columns
    ordered_drugs_dot <- consensus_top5 %>%
      filter(drug_key %in% target_expr$drug_key) %>%
      distinct(drug, drug_key, target_state) %>%
      mutate(drug_label = shorten_drug(drug, 25)) %>%
      arrange(match(target_state, state_order), drug_label) %>%
      pull(drug_label)
    
    target_expr <- target_expr %>%
      mutate(drug_label = factor(shorten_drug(drug, 25), levels = ordered_drugs_dot)) %>%
      filter(!is.na(drug_label))
    
    # Sort genes by mean expression across all states to make the heatmap cleaner
    gene_order <- target_expr %>%
      group_by(target_gene) %>%
      summarise(max_expr = max(mean_expr), .groups = "drop") %>%
      arrange(max_expr) %>%
      pull(target_gene)
    
    target_expr$target_gene <- factor(target_expr$target_gene, levels = gene_order)
    
    p_dot <- ggplot(target_expr, aes(x = expression_state, y = target_gene)) +
      geom_point(aes(size = pct_expr, color = mean_expr), alpha = 0.88) +
      geom_point(data = target_expr %>% filter(target_state_flag), shape = 21, stroke = 0.8, size = 4.2, color = "black", fill = NA) +
      scale_color_viridis_c(option = "C", name = "Mean\nexpression") +
      scale_size(range = c(0.7, 5.2), labels = scales::percent_format(accuracy = 1), name = "Cells\nexpressing") +
      facet_grid(state ~ drug_label, scales = "free_y", space = "free_y") +
      labs(x = NULL, y = "Drug target gene", title = "Mechanism Target Expression: Top Consensus Inhibitors") +
      theme_classic(base_size = 10) +
      theme(strip.text.x = element_text(face = "bold", size = 6.5), strip.text.y = element_text(face = "bold", angle = 0),
            axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
    ggsave(file.path(out_dir, "Auto_drug_reversal_mechanism_target_dotplot.pdf"), p_dot, width = 22, height = 11, limitsize = FALSE)
  }
}

####################
# predicted reversion visuals (Scatter, Heatmap)
####################

l1000 <- tryCatch(load_l1000_reference(), error = function(e) NULL)
profile_all <- tibble()
if (!is.null(l1000) && nrow(consensus_top5) > 0) {
  profile_all <- bind_rows(lapply(seq_len(nrow(consensus_top5)), function(i) {
    state_name <- as.character(consensus_top5$state[i])
    drug_key_i <- consensus_top5$drug_key[i]
    genes <- signature_dt %>% filter(state == state_name) %>% pull(gene) %>% unique()
    drug_profile <- get_drug_gene_profile(l1000$norm_mat, l1000$instance_map, drug_key_i, genes)
    if (is.null(drug_profile)) return(NULL)
    disease <- degs %>% filter(state == state_name, gene %in% drug_profile$gene) %>% select(gene, avg_logFC)
    drug_profile %>% left_join(disease, by = "gene") %>%
      mutate(state = state_name, drug = consensus_top5$drug[i], drug_key = drug_key_i,
             direction = ifelse(avg_logFC >= 0, "State-up genes", "State-down genes"),
             l1000_centered = l1000_centered, abs_logfc_rank = rank(-abs(avg_logFC), ties.method = "first"))
  }))
}

if (nrow(profile_all) > 0) {
  ordered_drugs_pred <- consensus_top5 %>%
    distinct(drug, drug_key, target_state) %>%
    mutate(drug_label = shorten_drug(drug, 25)) %>%
    arrange(match(target_state, state_order), drug_label) %>%
    pull(drug_label) %>%
    unique()

  scatter_df <- profile_all %>%
    group_by(state, drug_key) %>% mutate(highlight_gene = abs_logfc_rank <= 20) %>% ungroup() %>%
    mutate(state = factor(state, levels = state_order), drug_label = factor(shorten_drug(drug, 25), levels = ordered_drugs_pred))

  p_scatter <- ggplot(scatter_df, aes(x = avg_logFC, y = -l1000_centered)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey65") + geom_vline(xintercept = 0, linewidth = 0.2, color = "grey65") +
    geom_point(aes(color = direction), alpha = 0.45, size = 1.1) +
    geom_point(data = scatter_df %>% filter(highlight_gene), color = "black", fill = "#FFD92F", shape = 21, size = 1.9, stroke = 0.25) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, color = "grey20") +
    scale_color_manual(values = c("State-up genes" = "#CB181D", "State-down genes" = "#2171B5")) +
    facet_grid(state ~ drug_label, scales = "free", space = "free") +
    labs(x = "scRef state-vs-rest logFC", y = "Predicted opposing LINCS coordinate", color = NULL,
         title = "Predicted anti-correlation: malignant genes versus reference drug signatures",
         subtitle = "Yellow points represent the top 20 state-defining genes by absolute logFC",
         caption = "Final drugs require ASGARD + Local CMap/L1000 top-100 support; scDrugPrio is visual-only.") +
    theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", hjust = 0), strip.text.x = element_text(face = "bold", size = 6.5), strip.text.y = element_text(face = "bold", angle = 0), legend.position = "top")
  ggsave(file.path(out_dir, "Auto_drug_reversal_predicted_anticorrelation_scatter.pdf"), p_scatter, width = 24, height = 10, limitsize = FALSE)

  heat_df <- profile_all %>%
    group_by(state, drug_key) %>% mutate(gene_rank = rank(-abs(avg_logFC), ties.method = "first")) %>% filter(gene_rank <= 40) %>% ungroup() %>%
    transmute(state = factor(state, levels = state_order), drug_label = factor(shorten_drug(drug, 25), levels = ordered_drugs_pred),
              gene, direction, `scRef malignant logFC` = pmax(pmin(avg_logFC, 2.5), -2.5),
              `Predicted treatment\n(opposing LINCS coordinate)` = pmax(pmin(-l1000_centered * 2.5, 2.5), -2.5)) %>%
    pivot_longer(cols = c(`scRef malignant logFC`, `Predicted treatment\n(opposing LINCS coordinate)`), names_to = "column", values_to = "value") %>%
    group_by(state, drug_label) %>% mutate(gene_axis = factor(gene, levels = rev(unique(gene[order(direction, -abs(value))])))) %>% ungroup()

  p_heat <- ggplot(heat_df, aes(x = column, y = gene_axis, fill = value)) +
    geom_tile(color = "white", linewidth = 0.12) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), name = "Scaled\nvalue") +
    facet_grid(state ~ drug_label, scales = "free", space = "free") +
    labs(x = NULL, y = "Top state-defining genes", title = "Predicted state-flipping heatmap",
         caption = "Final drugs require ASGARD + Local CMap/L1000 top-100 support; scDrugPrio is visual-only.") +
    theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 25, hjust = 1), strip.text.x = element_text(face = "bold", size = 6.5), strip.text.y = element_text(face = "bold", angle = 0))
  ggsave(file.path(out_dir, "Auto_drug_reversal_predicted_state_flipping_heatmap.pdf"), p_heat, width = 24, height = 12, limitsize = FALSE)

  # Signature Reversal Profile (Connected points style matching requested image)
  profile_summary <- profile_all %>%
    group_by(state, drug, drug_key, direction) %>%
    summarise(
      mean_l1000_rank = mean(l1000_rank, na.rm = TRUE),
      median_l1000_rank = median(l1000_rank, na.rm = TRUE),
      n_signature_genes = n(),
      n_l1000_instances = first(n_l1000_instances),
      .groups = "drop"
    ) %>%
    mutate(
      state = factor(state, levels = state_order),
      drug_label = shorten_drug(drug, 25),
      direction = factor(direction, levels = c("State-up genes", "State-down genes"))
    )

  # Canonical colors for the states
  state_colors <- c(
    "Classic Proliferative" = "#E41A1C",
    "Basal to Intestinal Metaplasia" = "#4DAF4A",
    "SMG-like Metaplasia" = "#FF7F00",
    "Stress-adaptive" = "#984EA3",
    "3CA_EMT_and_Protein_maturation" = "#377EB8"
  )

  p_profile <- ggplot(profile_summary, aes(x = direction, y = mean_l1000_rank)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_line(aes(group = drug_key), color = "grey40", alpha = 0.4, linewidth = 0.5) +
    geom_point(aes(color = state, size = n_l1000_instances), alpha = 0.8) +
    geom_text_repel(
      aes(label = drug_label),
      size = 2.2,
      max.overlaps = 50,
      min.segment.length = 0,
      box.padding = 0.2,
      segment.alpha = 0.3
    ) +
    scale_color_manual(values = state_colors) +
    scale_size_continuous(range = c(1.5, 4.5)) +
    scale_y_continuous(limits = c(0.1, 0.85), breaks = seq(0, 1, 0.1)) +
    facet_wrap(~ state, ncol = 5) +
    labs(
      x = NULL,
      y = "Mean normalized L1000 rank",
      title = "L1000 Signature Reversal Profiles: Larger rank = downregulated",
      subtitle = "Transition from state-up to state-down gene ranks; Dotted line (y=0.5) = random expectation",
      caption = "Downward slope indicates reversal: State-up genes are downregulated in drug, State-down genes are upregulated."
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "none",
      axis.title.y = element_text(face = "bold")
    )

  ggsave(file.path(out_dir, "Auto_drug_reversal_l1000_signature_reversal_profiles.pdf"), p_profile, width = 24, height = 8, limitsize = FALSE)
  fwrite(profile_summary, file.path(out_dir, "Auto_drug_reversal_l1000_signature_reversal_profiles.csv"))
}

write_status("complete", "Generated unified visualizations. Final consensus ignores scDrugPrio and requires ASGARD + Local CMap/L1000 top-100 support.")
message("Drug-reversal unified visualizations complete.")
