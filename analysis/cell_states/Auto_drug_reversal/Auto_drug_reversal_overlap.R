####################
# Auto_drug_reversal_overlap.R
# Find and visualize overlapping chemical inhibitors between scRef and PDOs
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral"
scref_dir <- file.path(project_dir, "scRef_Pipeline", "ref_outs", "Auto_drug_reversal")
pdo_dir <- file.path(project_dir, "PDOs_Pipeline", "PDOs_outs", "Auto_drug_reversal")

out_dir <- file.path(scref_dir, "overlap_visuals")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_order <- c("Classic Proliferative", "Stress-adaptive")

# Read selected drugs from scDrugPrio
read_selected <- function(base_dir) {
  path <- file.path(base_dir, "scdrugprio_visuals", "Auto_scdrugprio_selected_drugs.csv")
  if(!file.exists(path)) return(NULL)
  df <- fread(path)
  df %>% filter(state %in% state_order)
}

scref_sel <- read_selected(scref_dir)
pdo_sel <- read_selected(pdo_dir)

# Read full rankings for waterfall
scref_rank <- fread(file.path(scref_dir, "scdrugprio", "Auto_scdrugprio_direction_audit_all_drugs.csv"))
if ("network_rank" %in% names(scref_rank) && !"rank" %in% names(scref_rank)) {
  setnames(scref_rank, "network_rank", "rank")
}
pdo_rank <- fread(file.path(pdo_dir, "scdrugprio", "Auto_scdrugprio_direction_audit_all_drugs.csv"))
if ("network_rank" %in% names(pdo_rank) && !"rank" %in% names(pdo_rank)) {
  setnames(pdo_rank, "network_rank", "rank")
}

scref_rank <- scref_rank %>% filter(state %in% state_order)
pdo_rank <- pdo_rank %>% filter(state %in% state_order)

scref_sel <- scref_sel %>% mutate(drug_key = tolower(trimws(drug)))
pdo_sel <- pdo_sel %>% mutate(drug_key = tolower(trimws(drug)))

overlap <- inner_join(
  scref_sel %>% select(state, drug_key, drug, scref_rank = rank, scref_score = score),
  pdo_sel %>% select(state, drug_key, pdo_rank = rank, pdo_score = score),
  by = c("state", "drug_key")
) %>% mutate(mean_rank = (scref_rank + pdo_rank)/2) %>%
  arrange(mean_rank)

print("Found overlaps:")
print(overlap)

top_overlap <- overlap %>% head(3)
fwrite(top_overlap, file.path(out_dir, "top_3_overlapped_drugs.csv"))

####################
# Waterfall plot
####################
# Standardize scref rankings
scref_rank_std <- scref_rank %>%
  mutate(
    drug_key = tolower(trimws(drug)),
    rank = as.numeric(rank)
  ) %>%
  group_by(state) %>%
  arrange(rank) %>%
  mutate(
    rank_index = row_number(),
    rank_evidence = -log10(rank / max(rank, na.rm = TRUE)),
    rank_evidence = ifelse(is.finite(rank_evidence), rank_evidence, NA_real_)
  ) %>%
  ungroup()

scref_rank_std <- scref_rank_std %>%
  mutate(
    final_hit = drug_key %in% top_overlap$drug_key,
    state = factor(as.character(state), levels = state_order),
    drug_label = ifelse(final_hit, drug, NA_character_)
  )

p_waterfall <- ggplot(scref_rank_std, aes(x = rank_index, y = rank_evidence)) +
  geom_col(aes(fill = final_hit), width = 0.9) +
  geom_text_repel(
    data = scref_rank_std %>% filter(final_hit),
    aes(label = drug_label),
    size = 3,
    min.segment.length = 0,
    box.padding = 0.5,
    max.overlaps = 80
  ) +
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "#CB181D"), guide = "none") +
  facet_wrap(~ state, scales = "free_x", ncol = 2) +
  labs(
    x = "Drug rank (scDrugPrio)",
    y = "Rank evidence: -log10(rank / universe)",
    title = "scDrugPrio Overlap: Rank Evidence (scRef perspective)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(file.path(out_dir, "overlap_waterfall.pdf"), p_waterfall, width = 8, height = 5, useDingbats = FALSE)
ggsave(file.path(out_dir, "overlap_waterfall.png"), p_waterfall, width = 8, height = 5, dpi = 300)


####################
# LINCS Scatter
####################
# Re-use helpers from scDrugPrio script
load_l1000_reference <- function() {
  default_ref_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000/DrugReference"
  rank_path <- file.path(default_ref_dir, "stomach_rankMatrix.txt")
  gene_path <- file.path(default_ref_dir, "stomach_gene_info.txt")
  drug_path <- file.path(default_ref_dir, "stomach_drug_info.txt")

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

l1000 <- tryCatch(load_l1000_reference(), error = function(e) NULL)
degs <- fread(file.path(scref_dir, "Auto_drug_reversal_degs_all_states.csv.gz")) %>% 
  filter(state %in% state_order) %>% mutate(avg_logFC = avg_log2FC)
signature_dt <- fread(file.path(scref_dir, "Auto_drug_reversal_signature_top150.csv")) %>%
  filter(state %in% state_order)

profile_all <- tibble()
if (!is.null(l1000) && nrow(top_overlap) > 0) {
  profile_all <- bind_rows(lapply(seq_len(nrow(top_overlap)), function(i) {
    state_name <- top_overlap$state[i]
    drug_key_i <- top_overlap$drug_key[i]
    
    genes <- signature_dt %>% filter(state == state_name) %>% pull(gene) %>% unique()
    drug_profile <- get_drug_gene_profile(l1000$norm_mat, l1000$instance_map, drug_key_i, genes)
    if (is.null(drug_profile)) return(NULL)

    disease <- degs %>%
      filter(state == state_name, gene %in% drug_profile$gene) %>%
      select(gene, avg_logFC, p_val_adj)

    drug_profile %>%
      left_join(disease, by = "gene") %>%
      mutate(
        state = state_name,
        drug = top_overlap$drug[i],
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
      state = factor(state, levels = state_order)
    )

  p_scatter <- ggplot(scatter_df, aes(x = avg_logFC, y = -l1000_centered)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey65") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey65") +
    geom_point(aes(color = direction), alpha = 0.5, size = 1.0) +
    geom_point(data = scatter_df %>% filter(highlight_gene), color = "black", fill = "#FFD92F", shape = 21, size = 1.6, stroke = 0.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "grey20") +
    scale_color_manual(values = c("State-up genes" = "#CB181D", "State-down genes" = "#2171B5")) +
    facet_grid(state ~ drug, scales = "free") +
    labs(
      x = "scRef state logFC",
      y = "Predicted treatment effect (inverted L1000)",
      color = NULL,
      title = "Overlapped Drugs: Predicted transcriptomic reversal (LINCS/L1000)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.text.x = element_text(face = "bold", size = 8),
      strip.text.y = element_text(face = "bold", angle = 0),
      legend.position = "top"
    )

  ggsave(file.path(out_dir, "overlap_reversion_scatter.pdf"), p_scatter, width = 12, height = 8, useDingbats = FALSE)
  ggsave(file.path(out_dir, "overlap_reversion_scatter.png"), p_scatter, width = 12, height = 8, dpi = 300)
}
print("Success! Outputs generated in overlap_visuals")

