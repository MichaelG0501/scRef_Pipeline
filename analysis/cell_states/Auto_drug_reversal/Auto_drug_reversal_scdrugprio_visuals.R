####################
# Auto_drug_reversal_scdrugprio_visuals.R
#
# Focused visualization of scDrugPrio results with a full-hub PPI network overlay.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
})

####################
# setup
####################

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

base_dir <- "Auto_drug_reversal"
out_dir <- file.path(base_dir, "scdrugprio_visuals")
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
      step = "scdrugprio_visuals",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_scdrugprio_visuals_status.csv")
  )
}

shorten_drug <- function(x, n = 32) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

standardize_rank_file <- function(path, method_name) {
  if (!file.exists(path)) return(tibble())
  df <- fread(path)
  if (nrow(df) == 0) return(tibble())

  nms <- names(df)
  pick <- function(candidates) {
    hit <- intersect(candidates, nms)[1]
    if (is.na(hit)) hit <- grep(paste(candidates, collapse = "|"), nms, ignore.case = TRUE, value = TRUE)[1]
    hit
  }

  state_col <- pick(c("state", "State"))
  drug_col <- pick(c("drug", "Drug", "pert_iname", "compound", "drug_name"))
  rank_col <- pick(c("rank", "Rank"))
  score_col <- pick(c("score", "drug_score", "connectivity_score", "rank_metric"))
  target_col <- pick(c("target_genes", "targets", "target", "counteracting_targets"))
  moa_col <- pick(c("moa", "MOA", "mechanism"))

  if (is.na(state_col) || is.na(drug_col)) return(tibble())

  out <- tibble(
    state = as.character(df[[state_col]]),
    method = method_name,
    drug = as.character(df[[drug_col]]),
    rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(df[[rank_col]])) else NA_real_,
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    target_genes = if (!is.na(target_col)) as.character(df[[target_col]]) else NA_character_,
    moa = if (!is.na(moa_col)) as.character(df[[moa_col]]) else NA_character_
  ) %>%
    mutate(
      state = ifelse(state %in% names(state_safe_map), state_safe_map[state], state),
      state = factor(state, levels = state_order),
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
    group_by(state, drug_key) %>%
    arrange(rank, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
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
  if (!all(file.exists(c(rank_path, gene_path, drug_path)))) {
    stop("Rank matrix, gene info, or drug info path does not exist.")
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

load_drugbank_targets <- function(drug_keys) {
  drugbank_path <- Sys.getenv(
    "AUTO_DRUGBANK_TARGETS",
    "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/all_drug_targets_drug_bank.txt"
  )
  if (!file.exists(drugbank_path)) return(tibble())

  db <- fread(drugbank_path)
  required <- c("drugID", "drug_name", "gene_symbol", "target_organism", "drug_action")
  if (!all(required %in% names(db))) return(tibble())

  db %>%
    transmute(
      drug_key = tolower(trimws(drug_name)),
      target_gene = as.character(gene_symbol),
      action = as.character(target_organism),
      organism = as.character(drug_action)
    ) %>%
    filter(
      drug_key %in% drug_keys,
      organism == "Humans",
      !is.na(target_gene), nzchar(target_gene),
      !is.na(action), nzchar(action)
    ) %>%
    distinct()
}

map_symbols_to_entrez <- function(symbols) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(tibble(gene = character(), entrez = character()))
  }
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  mapped %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    transmute(gene = SYMBOL, entrez = ENTREZID)
}

####################
# load core data
####################

deg_path <- file.path(base_dir, "Auto_drug_reversal_degs_all_states.csv.gz")
sig_path <- file.path(base_dir, "Auto_drug_reversal_signature_top150.csv")
if (!file.exists(deg_path) || !file.exists(sig_path)) {
  write_status("missing_inputs", "Missing DEG or top-signature file.")
  quit(save = "no", status = 0)
}

degs <- fread(deg_path) %>%
  mutate(state = factor(as.character(state), levels = state_order)) %>%
  filter(!is.na(state))
logfc_col <- intersect(c("avg_logFC", "avg_log2FC", "logFC", "log2FC"), names(degs))[1]
if (is.na(logfc_col)) {
  write_status("missing_logfc", "Could not find a logFC column in DEG table.")
  quit(save = "no", status = 0)
}
degs <- degs %>% mutate(avg_logFC = .data[[logfc_col]])
signature_dt <- fread(sig_path) %>%
  mutate(state = factor(as.character(state), levels = state_order)) %>%
  filter(!is.na(state))

# FOCUS ONLY ON scDrugPrio
rank_file <- file.path(base_dir, "scdrugprio", "Auto_scdrugprio_direction_audit_all_drugs.csv")
if (!file.exists(rank_file)) {
  write_status("missing_scdrugprio", "Missing scDrugPrio ranking file.")
  quit(save = "no", status = 0)
}

rankings <- standardize_rank_file(rank_file, "scDrugPrio")
if (nrow(rankings) == 0) {
  write_status("empty_rankings", "scDrugPrio ranking file is empty.")
  quit(save = "no", status = 0)
}

audit_file <- file.path(base_dir, "scdrugprio", "Auto_scdrugprio_direction_audit_all_drugs.csv")
audit_rankings <- if (file.exists(audit_file)) {
  audit_df <- fread(audit_file)
  if ("network_rank" %in% names(audit_df) && !"rank" %in% names(audit_df)) {
    setnames(audit_df, "network_rank", "rank")
  }
  audit_tmp <- tempfile(fileext = ".csv")
  fwrite(audit_df, audit_tmp)
  standardize_rank_file(audit_tmp, "scDrugPrio")
} else {
  rankings
}

# Select scDrugPrio-only candidates for visualization highlighting.
final_candidates <- rankings %>%
  group_by(state) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

fwrite(final_candidates, file.path(out_dir, "Auto_scdrugprio_selected_drugs.csv"))

# Keep an unfiltered network-proximity view for biological audit plots.
audit_candidates <- audit_rankings %>%
  group_by(state) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()
fwrite(audit_candidates, file.path(out_dir, "Auto_scdrugprio_audit_network_candidates.csv"))

####################
# waterfall ranking plots
####################

highlight_keys <- final_candidates %>% transmute(state, drug_key, final_hit = TRUE)
waterfall <- audit_rankings %>%
  group_by(state) %>%
  arrange(rank, .by_group = TRUE) %>%
  mutate(
    rank_index = row_number(),
    rank_evidence = -log10(rank / max(rank, na.rm = TRUE)),
    rank_evidence = ifelse(is.finite(rank_evidence), rank_evidence, NA_real_)
  ) %>%
  ungroup() %>%
  left_join(highlight_keys, by = c("state", "drug_key")) %>%
  mutate(
    final_hit = ifelse(is.na(final_hit), FALSE, final_hit),
    state = factor(as.character(state), levels = state_order),
    drug_label = ifelse(final_hit, shorten_drug(drug, 22), NA_character_)
  )

p_waterfall <- ggplot(waterfall, aes(x = rank_index, y = rank_evidence)) +
  geom_col(aes(fill = final_hit), width = 0.9) +
  geom_text_repel(
    data = waterfall %>% filter(final_hit),
    aes(label = drug_label),
    size = 2.2,
    min.segment.length = 0,
    box.padding = 0.2,
    max.overlaps = 80
  ) +
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "#CB181D"), guide = "none") +
  facet_wrap(~ state, scales = "free_x", ncol = 3) +
  labs(
    x = "Drug rank (scDrugPrio)",
    y = "Rank evidence: -log10(rank / universe)",
    title = "scDrugPrio network-proximity candidates highlighted for visualization"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(file.path(out_dir, "Auto_scdrugprio_waterfall.pdf"), p_waterfall, width = 12, height = 8, useDingbats = FALSE)
ggsave(file.path(out_dir, "Auto_scdrugprio_waterfall.png"), p_waterfall, width = 12, height = 8, dpi = 300)

####################
# L1000-based predicted anti-correlation
####################

l1000 <- tryCatch(load_l1000_reference(), error = function(e) NULL)
profile_all <- tibble()

if (!is.null(l1000) && nrow(final_candidates) > 0) {
  profile_all <- bind_rows(lapply(seq_len(nrow(final_candidates)), function(i) {
    state_name <- as.character(final_candidates$state[i])
    drug_key_i <- final_candidates$drug_key[i]
    genes <- signature_dt %>%
      filter(state == state_name) %>%
      pull(gene) %>%
      unique()
    drug_profile <- get_drug_gene_profile(l1000$norm_mat, l1000$instance_map, drug_key_i, genes)
    if (is.null(drug_profile)) return(NULL)

    disease <- degs %>%
      filter(state == state_name, gene %in% drug_profile$gene) %>%
      select(gene, avg_logFC, p_val_adj)

    drug_profile %>%
      left_join(disease, by = "gene") %>%
      mutate(
        state = state_name,
        drug = final_candidates$drug[i],
        drug_key = drug_key_i,
        direction = ifelse(avg_logFC >= 0, "State-up genes", "State-down genes"),
        predicted_opposition = -sign(avg_logFC) * l1000_centered,
        abs_logfc_rank = rank(-abs(avg_logFC), ties.method = "first")
      )
  }))
}

if (nrow(profile_all) > 0) {
  scatter_df <- profile_all %>%
    group_by(state, drug_key) %>%
    mutate(highlight_gene = abs_logfc_rank <= 15) %>%
    ungroup() %>%
    mutate(
      state = factor(state, levels = state_order),
      drug_label = shorten_drug(drug, 25)
    )

  p_scatter <- ggplot(scatter_df, aes(x = avg_logFC, y = -l1000_centered)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey65") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey65") +
    geom_point(aes(color = direction), alpha = 0.5, size = 1.0) +
    geom_point(data = scatter_df %>% filter(highlight_gene), color = "black", fill = "#FFD92F", shape = 21, size = 1.6, stroke = 0.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "grey20") +
    scale_color_manual(values = c("State-up genes" = "#CB181D", "State-down genes" = "#2171B5")) +
    facet_grid(state ~ drug_label, scales = "free") +
    labs(
      x = "scRef state logFC",
      y = "Predicted treatment effect (inverted L1000)",
      color = NULL,
      title = "scDrugPrio: Predicted transcriptomic reversal (LINCS/L1000)"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.text.x = element_text(face = "bold", size = 6),
      strip.text.y = element_text(face = "bold", angle = 0),
      legend.position = "top"
    )

  ggsave(file.path(out_dir, "Auto_scdrugprio_reversion_scatter.pdf"), p_scatter, width = 16, height = 11, useDingbats = FALSE)
}

####################
# PPI targeted-hub overlay (REFINED: show full DEG web)
####################

ppi_path <- Sys.getenv("AUTO_SCDRUGPRIO_PPI", "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/ppi.txt")

if (file.exists(ppi_path) && nrow(audit_candidates) > 0) {
  target_tbl <- load_drugbank_targets(audit_candidates$drug_key)
  if (nrow(target_tbl) > 0) {
    # Match symbols to entrez (assuming PPI uses entrez based on scdrugprio script)
    all_symbols <- unique(c(degs$gene, target_tbl$target_gene))
    symbol_map <- map_symbols_to_entrez(all_symbols)
    
    ppi <- fread(ppi_path)
    setnames(ppi, names(ppi)[1:2], c("from", "to"))
    ppi[, `:=`(from = as.character(from), to = as.character(to))]

    network_rows <- lapply(state_order, function(state_name) {
      state_hits <- audit_candidates %>% filter(state == state_name)
      if (nrow(state_hits) == 0) return(NULL)

      # Top 250 DEGs
      state_deg <- degs %>%
        filter(state == state_name) %>%
        arrange(p_val_adj, desc(abs(avg_logFC))) %>%
        slice_head(n = 120) %>%
        left_join(symbol_map, by = "gene") %>%
        filter(!is.na(entrez))
      if (nrow(state_deg) == 0) return(NULL)

      # Target genes for top drugs
      state_targets <- target_tbl %>%
        filter(drug_key %in% state_hits$drug_key) %>%
        left_join(symbol_map, by = c("target_gene" = "gene")) %>%
        filter(!is.na(entrez))
      
      # Combined node set: DEGs + Targets
      keep_entrez <- unique(c(state_deg$entrez, state_targets$entrez))
      
      # REFINED: Show full interconnected web of these nodes
      edge_sub <- ppi[from %in% keep_entrez & to %in% keep_entrez]
      if (nrow(edge_sub) == 0) return(NULL)

      ppi_edges <- edge_sub %>% select(from, to) %>% mutate(edge_type = "PPI", action = "PPI")
      drug_edges <- state_targets %>% select(from = drug_key, to = entrez, action) %>% mutate(edge_type = "Drug-Target")

      all_edges <- bind_rows(ppi_edges, drug_edges) %>% distinct(from, to, edge_type, .keep_all = TRUE)

      graph <- igraph::graph_from_data_frame(all_edges, directed = FALSE)
      deg_vec <- igraph::degree(graph)
      
      coords <- as.data.frame(igraph::layout_with_fr(graph))
      colnames(coords) <- c("x", "y")
      coords$node_id <- igraph::V(graph)$name

      sym_vec <- setNames(symbol_map$gene, symbol_map$entrez)

      node_tbl <- coords %>%
        mutate(
          is_drug = node_id %in% state_hits$drug_key,
          is_target = node_id %in% state_targets$entrez,
          is_deg = node_id %in% state_deg$entrez
        ) %>%
        left_join(state_deg %>% select(entrez, avg_logFC), by = c("node_id" = "entrez")) %>%
        left_join(state_hits %>% select(drug_key, drug), by = c("node_id" = "drug_key")) %>%
        mutate(
          state = state_name,
          degree = unname(deg_vec[node_id]),
          gene_symbol = sym_vec[node_id],
          label = case_when(
            is_drug ~ drug,
            is_target ~ gene_symbol,
            degree >= quantile(degree, 0.96, na.rm = TRUE) ~ gene_symbol,
            TRUE ~ NA_character_
          )
        )

      edge_tbl <- igraph::as_data_frame(graph, what = "edges") %>%
        left_join(coords %>% select(from = node_id, x, y), by = "from") %>%
        rename(x_from = x, y_from = y) %>%
        left_join(coords %>% select(to = node_id, x, y), by = "to") %>%
        rename(x_to = x, y_to = y) %>%
        mutate(
          state = state_name,
          is_ppi = edge_type == "PPI",
          action_type = case_when(
            is_ppi ~ "PPI",
            grepl("inhibitor|antagonist|blocker|suppress", tolower(action)) ~ "Inhibitor",
            grepl("agonist|activator|stimulat|induc|posit", tolower(action)) ~ "Activator",
            TRUE ~ "Other/Unknown"
          )
        )

      list(nodes = node_tbl, edges = edge_tbl)
    })
    network_rows <- Filter(Negate(is.null), network_rows)

    if (length(network_rows) > 0) {
      node_all <- bind_rows(lapply(network_rows, `[[`, "nodes"))
      edge_all <- bind_rows(lapply(network_rows, `[[`, "edges"))
      
      p_net <- ggplot() +
        geom_segment(
          data = edge_all %>% filter(is_ppi),
          aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
          linewidth = 0.3,
          color = "grey60",
          alpha = 0.6
        ) +
        geom_segment(
          data = edge_all %>% filter(!is_ppi),
          aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = action_type),
          linewidth = 0.8,
          alpha = 0.9,
          linetype = "dashed"
        ) +
        scale_color_manual(values = c("Inhibitor" = "#E41A1C", "Activator" = "#4DAF4A", "Other/Unknown" = "grey30", "PPI" = "grey60"), name = "Drug Action") +
        geom_point(
          data = node_all %>% filter(!is_drug & !is_target),
          aes(x = x, y = y, fill = avg_logFC, size = degree),
          shape = 21, color = "black", stroke = 0.2, alpha = 0.9
        ) +
        geom_point(
          data = node_all %>% filter(is_target),
          aes(x = x, y = y, fill = avg_logFC, size = degree),
          shape = 23, color = "black", stroke = 0.8
        ) +
        geom_point(
          data = node_all %>% filter(is_drug),
          aes(x = x, y = y),
          shape = 22, fill = "#984EA3", color = "black", size = 5, stroke = 1.0
        ) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "scRef logFC") +
        scale_size(range = c(1.0, 5.0), name = "PPI degree") +
        geom_text_repel(
          data = node_all %>% filter(!is.na(label), nzchar(label)),
          aes(x = x, y = y, label = label, fontface = ifelse(is_drug, "bold", "plain")),
          size = 2.5,
          max.overlaps = 50,
          min.segment.length = 0,
          box.padding = 0.3
        ) +
        facet_wrap(~ state, scales = "free", ncol = 2) +
        labs(title = "scDrugPrio: Targeted-hub overlay (Full top-250 DEG web)") +
        theme_void(base_size = 9) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0),
          strip.text = element_text(face = "bold")
        )

      ggsave(file.path(out_dir, "Auto_scdrugprio_ppi_hub_refined.pdf"), p_net, width = 14, height = 12, useDingbats = FALSE)
      ggsave(file.path(out_dir, "Auto_scdrugprio_ppi_hub_refined.png"), p_net, width = 14, height = 12, dpi = 300)
    }
  }
}

write_status("complete", "Generated focused scDrugPrio visualizations with refined PPI hub overlay.")
message("Focused scDrugPrio visualizations complete.")
