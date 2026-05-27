####################
# Auto_drug_reversal_overlap_final.R
# Find and visualize overlapping chemical inhibitors between scRef and PDOs
# Using final ASGARD+CMap candidates
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
state_colors <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

shorten_drug <- function(x, n = 32) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

# Read final selection from method_visuals
read_selected <- function(base_dir) {
  path <- file.path(base_dir, "method_visuals", "Auto_drug_reversal_final_asgard_cmap_top100_membership.csv")
  if(!file.exists(path)) return(NULL)
  df <- fread(path)
  df %>% filter(state %in% state_order)
}

scref_sel <- read_selected(scref_dir)
pdo_sel <- read_selected(pdo_dir)

scref_sel <- scref_sel %>% mutate(drug_key = tolower(trimws(drug_key)))
pdo_sel <- pdo_sel %>% mutate(drug_key = tolower(trimws(drug_key)))

overlap <- inner_join(
  scref_sel %>% select(state, drug_key, drug, scref_mean_rank = mean_rank),
  pdo_sel %>% select(state, drug_key, pdo_mean_rank = mean_rank),
  by = c("state", "drug_key")
) %>% mutate(overall_mean_rank = (scref_mean_rank + pdo_mean_rank)/2)

print("Found overlaps in final selection:")
print(overlap)

# Select top 3 per state
top_overlap <- overlap %>% 
  group_by(state) %>%
  arrange(overall_mean_rank, .by_group=TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

fwrite(top_overlap, file.path(out_dir, "top_3_overlapped_drugs_per_state.csv"))

####################
# Rank Comparison Plot (scRef vs PDO)
####################
# We plot scRef mean rank on X and PDO mean rank on Y for the overlaps.

scatter_overlap <- overlap %>%
  mutate(
    is_top3 = drug_key %in% top_overlap$drug_key,
    state = factor(state, levels = state_order)
  )

p_rank <- ggplot(scatter_overlap, aes(x = scref_mean_rank, y = pdo_mean_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = is_top3), shape = 21, size = 3, color = "black", alpha = 0.8) +
  geom_text_repel(
    data = scatter_overlap %>% filter(is_top3),
    aes(label = drug),
    size = 3.5,
    fontface = "bold",
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#CB181D")) +
  scale_x_reverse() + 
  scale_y_reverse() +
  facet_wrap(~ state, scales = "free", ncol=2) +
  labs(
    x = "Mean Rank (scRef pipeline) [reversed]",
    y = "Mean Rank (PDOs pipeline) [reversed]",
    title = "Rank Comparison: Overlapping Chemical Inhibitors",
    subtitle = "Top-right corner = highly ranked in both pipelines. Red points are the 6 final selections (3 per state)."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "overlap_final_rank_comparison.pdf"), p_rank, width = 10, height = 5, useDingbats = FALSE)
ggsave(file.path(out_dir, "overlap_final_rank_comparison.png"), p_rank, width = 10, height = 5, dpi = 300)


####################
# LINCS Signature Reversal Profiles
####################
load_l1000_reference <- function() {
  default_ref_dir <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/DrugReference"
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
  # Signature Reversal Profile Plot
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

  p_profile <- ggplot(profile_summary, aes(x = direction, y = mean_l1000_rank)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_line(aes(group = drug_key), color = "grey40", alpha = 0.4, linewidth = 0.5) +
    geom_point(aes(color = state, size = n_l1000_instances), alpha = 0.8) +
    geom_text_repel(
      aes(label = drug_label),
      size = 3.5,
      max.overlaps = 50,
      min.segment.length = 0,
      box.padding = 0.3,
      segment.alpha = 0.4
    ) +
    scale_color_manual(values = state_colors) +
    scale_size_continuous(range = c(2, 6)) +
    scale_y_continuous(limits = c(0.1, 0.9), breaks = seq(0, 1, 0.1)) +
    facet_wrap(~ state, ncol = 2) +
    labs(
      x = NULL,
      y = "Mean normalized L1000 rank",
      title = "Top 3 Overlapped Drugs: L1000 Signature Reversal Profiles",
      subtitle = "Transition from state-up to state-down gene ranks; Dotted line (y=0.5) = random expectation",
      caption = "Downward slope indicates reversal: State-up genes are downregulated in drug, State-down genes are upregulated."
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none",
      axis.title.y = element_text(face = "bold")
    )

  ggsave(file.path(out_dir, "overlap_final_l1000_signature_reversal_profiles.pdf"), p_profile, width = 12, height = 7, useDingbats = FALSE)
  ggsave(file.path(out_dir, "overlap_final_l1000_signature_reversal_profiles.png"), p_profile, width = 12, height = 7, dpi = 300)
} else {
  print("Failed to map drugs to L1000 instance_map.")
}

print("Success! Outputs generated in overlap_visuals")

