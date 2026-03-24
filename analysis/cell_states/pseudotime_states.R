####################
# Auto_pseudotime_states.R
# Monocle3 pseudotime analysis for the noreg Approach B cell states.
#
# Part A: Per-sample pseudotime using the 5 defined states.
#   - Top 12 samples ranked by diversity (geometric mean across 5 states)
#   - Root = Basal to Intestinal Metaplasia state cells
#   - Labels: 5 defined states only (no Unresolved/Hybrid)
#
# Part B: Per-sample pseudotime within 3 state subsets:
#   - Basal to Intestinal Metaplasia: MP labels from topMP of {MP17, MP14, MP5, MP10, MP8}, root = MP17
#   - SMG-like Metaplasia: MP labels from topMP of {MP18, MP16}, root = MP18
#   - Stress-adaptive: MP labels from topMP of {MP13, MP12}, root = MP13
#   - Also top 12 samples per state subset (by cell count)
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#
# Output:
#   ref_outs/task1_pseudotime_updated_naming/partA/Auto_task1_partA_top12_pseudotime_states.pdf (12 pages)
#   ref_outs/task1_pseudotime_updated_naming/partA/<sample>_pseudotime_states.rds
#   ref_outs/task1_pseudotime_updated_naming/partB/<state>/Auto_task1_partB_<state>_top12_pseudotime.pdf (up to 12 pages)
#   ref_outs/task1_pseudotime_updated_naming/partB/<state>/<sample>_pseudotime_<state>.rds
#   updates/new_updates/summaries/Auto_pseudotime_states_summary.csv
####################

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# Create output directories
####################
task_prefix <- "task1"
out_root <- paste0(task_prefix, "_pseudotime_updated_naming")
dir.create(file.path(out_root, "partA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "partB", "Basal_to_Intestinal_Metaplasia"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "partB", "SMG_like_Metaplasia"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "partB", "Stress_adaptive"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "partB", "Basal_and_SMG_Metaplasia"), recursive = TRUE, showWarnings = FALSE)

summary_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/updates/new_updates/summaries/"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Constants
####################
state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8"
)

# Identify any extra states (e.g. 3CA relabeled)
extra_states <- setdiff(unique(as.character(state_B)), c(names(group_cols), "Unresolved", "Hybrid"))
if (length(extra_states) > 0) {
  new_cols <- setNames(scales::hue_pal()(length(extra_states)), extra_states)
  group_cols <- c(group_cols, new_cols)
}

mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

# State subset definitions: MPs and root MP
state_subsets <- list(
  "Basal to Intestinal Metaplasia" = list(
    mps = c("MP17", "MP14", "MP5", "MP10", "MP8"),
    root_mp = "MP17"
  ),
  "SMG-like Metaplasia" = list(
    mps = c("MP18", "MP16"),
    root_mp = "MP18"
  ),
  "Stress-adaptive" = list(
    mps = c("MP13", "MP12"),
    root_mp = "MP13"
  ),
  "Basal and SMG Metaplasia" = list(
    mps = c("MP17", "MP14", "MP5", "MP10", "MP8", "MP18", "MP16"),
    root_mp = "MP17",
    source_states = c("Basal to Intestinal Metaplasia", "SMG-like Metaplasia")
  )
)

state_dir_map <- c(
  "Basal to Intestinal Metaplasia" = "Basal_to_Intestinal_Metaplasia",
  "SMG-like Metaplasia" = "SMG_like_Metaplasia",
  "Stress-adaptive" = "Stress_adaptive",
  "Basal and SMG Metaplasia" = "Basal_and_SMG_Metaplasia"
)

####################
# Load data
message("Loading data ...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states_path <- "Auto_final_states.rds"
if (file.exists(final_states_path)) {
  state_B <- readRDS(final_states_path)
} else {
  state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
}
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")

# Align cells
common_cells <- intersect(names(state_B), Cells(tmdata_all))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))
tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

# Add state to metadata
tmdata_all$state_B <- state_B[Cells(tmdata_all)]

####################
# Helper: find root node closest to root state cells
####################
get_root_node <- function(cds, root_cells) {
  # Find principal graph node most represented among root cells
  cell_ids <- which(colnames(cds) %in% root_cells)
  if (length(cell_ids) == 0) {
    warning("No root cells found in CDS. Using first principal node.")
    return(igraph::V(principal_graph(cds)[["UMAP"]])$name[1])
  }

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  closest_vals <- closest_vertex[cell_ids, 1]
  root_vertex_table <- table(closest_vals)
  top_key <- names(which.max(root_vertex_table))

  graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  if (top_key %in% graph_nodes) {
    root_node <- top_key
  } else if (!is.na(suppressWarnings(as.numeric(top_key)))) {
    idx <- as.integer(top_key)
    idx <- max(1, min(length(graph_nodes), idx))
    root_node <- graph_nodes[idx]
  } else {
    root_node <- graph_nodes[1]
  }

  root_node
}

####################
# Helper: run Seurat preprocessing → monocle3 pseudotime
####################
run_pseudotime <- function(
  seurat_obj,
  group_col,
  root_cells,
  sample_id,
  left_title_prefix,
  colour_map,
  legend_labels,
  legend_title,
  n_cells
) {
  # Seurat preprocessing
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  n_pcs <- min(30, ncol(seurat_obj) - 1)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, verbose = FALSE)

  # Convert to CDS
  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  root_cells_use <- intersect(root_cells, colnames(cds))
  if (length(root_cells_use) == 0) {
    stop("No root cells overlap with CDS columns")
  }
  cds <- order_cells(cds, root_cells = root_cells_use)

  # Extract pseudotime
  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA

  # Plots
  p1 <- plot_cells(cds, color_cells_by = group_col, show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE, label_groups_by_cluster = FALSE, 
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
    labs(title = paste0(left_title_prefix, " - ", sample_id, " (n = ", n_cells, ")"), color = NULL) +
    theme_minimal()

  p1 <- p1 + scale_color_manual(
    values = colour_map,
    breaks = names(colour_map),
    labels = legend_labels[names(colour_map)],
    name = legend_title,
    na.value = "grey80",
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 4))
  )

  p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE, label_groups_by_cluster = FALSE, 
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
    scale_color_viridis_c() +
    labs(title = paste0("Pseudotime - ", sample_id), color = "Pseudotime") +
    theme_minimal()

  list(cds = cds, pseudotime = pt, plots = list(group = p1, pseudotime = p2))
}

####################
# PART A: Per-sample pseudotime with 5 defined states
####################
message("=== PART A: Per-sample state pseudotime ===")

# Compute diversity ranking
target_states <- setdiff(names(group_cols), c("Unresolved", "Hybrid"))
defined_cells <- names(state_B)[state_B %in% target_states]

state_df <- data.frame(
  cell = names(state_B),
  state = as.character(state_B),
  orig.ident = as.character(tmdata_all$orig.ident[names(state_B)]),
  stringsAsFactors = FALSE
)

sample_totals <- state_df %>% count(orig.ident, name = "total_n")

counts_long <- state_df %>%
  filter(state %in% target_states) %>%
  count(orig.ident, state, .drop = FALSE) %>%
  complete(orig.ident, state = target_states, fill = list(n = 0))

rank_df <- counts_long %>%
  group_by(orig.ident) %>%
  summarise(
    target_n = sum(n),
    geo_mean_score = exp(mean(log(n + 1))),
    .groups = "drop"
  ) %>%
  left_join(sample_totals, by = "orig.ident") %>%
  filter(target_n > 20) %>%
  arrange(desc(geo_mean_score)) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 12)

top12_samples <- rank_df$orig.ident
message(sprintf("Top 12 samples for Part A: %s", paste(top12_samples, collapse = ", ")))

partA_summary <- list()

pdf(file.path(out_root, "partA", paste0("Auto_", task_prefix, "_partA_top12_pseudotime_states.pdf")),
    width = 14, height = 6, onefile = TRUE)

for (i in seq_along(top12_samples)) {

  sample_id <- top12_samples[i]
  message(sprintf("Part A [%d/12]: %s", i, sample_id))

  # Subset to this sample's defined-state cells only
  sample_cells <- names(state_B)[state_B %in% target_states &
                                   tmdata_all$orig.ident[names(state_B)] == sample_id]

  if (length(sample_cells) < 50) {
    message(sprintf("  Skipping %s: only %d defined-state cells", sample_id, length(sample_cells)))
    next
  }

  sub_obj <- tmdata_all[, sample_cells]
  sub_obj$state_label <- factor(state_B[sample_cells], levels = target_states)

  state_counts_sample <- table(factor(as.character(state_B[sample_cells]), levels = target_states))
  state_legend_labels <- setNames(
    paste0(target_states, " (", as.integer(state_counts_sample[target_states]), ")"),
    target_states
  )

  # Root cells = Basal to Intestinal Metaplasia
  root_cells <- sample_cells[state_B[sample_cells] == "Basal to Intestinal Metaplasia"]
  if (length(root_cells) == 0) {
    message(sprintf("  Skipping %s: no Basal to Intestinal Metaplasia cells", sample_id))
    next
  }

  result <- tryCatch({
    run_pseudotime(
      seurat_obj = sub_obj,
      group_col = "state_label",
      root_cells = root_cells,
      sample_id = sample_id,
      left_title_prefix = "Cell states",
      colour_map = group_cols[target_states],
      legend_labels = state_legend_labels,
      legend_title = "State",
      n_cells = length(sample_cells)
    )
  }, error = function(e) {
    message(sprintf("  Error in %s: %s", sample_id, e$message))
    NULL
  })

  if (is.null(result)) next

  # add one page per sample in single Part-A PDF
  print(result$plots$group + result$plots$pseudotime)

  # Save pseudotime values
  saveRDS(result$pseudotime, file.path(out_root, "partA", paste0(sample_id, "_pseudotime_states.rds")))

  partA_summary[[sample_id]] <- data.frame(
    part = "A",
    sample = sample_id,
    rank = i,
    n_cells = length(sample_cells),
    n_root = length(root_cells),
    median_pseudotime = median(result$pseudotime, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  message(sprintf("  Done: %d cells, median PT = %.2f",
                  length(sample_cells),
                  median(result$pseudotime, na.rm = TRUE)))
}
dev.off()

####################
# PART B: Per-sample pseudotime within 3 state subsets
####################
message("\n=== PART B: Per-sample MP-level pseudotime within state subsets ===")

partB_summary <- list()
used_partB_samples <- top12_samples

for (state_name in names(state_subsets)) {
  subset_info <- state_subsets[[state_name]]
  mps_in_state <- subset_info$mps
  root_mp <- subset_info$root_mp

  message(sprintf("\n--- State: %s (root: %s) ---", state_name, root_mp))


  # Get cells belonging to this state
  source_states <- subset_info$source_states
  if (is.null(source_states)) source_states <- state_name
  
  state_cells <- names(state_B)[state_B %in% source_states]
  if (length(state_cells) == 0) {
    message(sprintf("  No cells for state %s, skipping", state_name))
    next
  }

  # Assign MP labels from topMP within relevant MPs only
  mp_cols_avail <- intersect(mps_in_state, colnames(mp_adj_noncc))
  if (length(mp_cols_avail) == 0) {
    message(sprintf("  No MP columns available for %s, skipping", state_name))
    next
  }

  mp_sub <- mp_adj_noncc[state_cells, mp_cols_avail, drop = FALSE]
  top_mp_idx <- max.col(mp_sub, ties.method = "first")
  top_mp_label <- mp_cols_avail[top_mp_idx]
  names(top_mp_label) <- state_cells

  # MP colour palette for this subset, strict state-defined order
  ordered_mps_state <- mps_in_state[mps_in_state %in% unique(top_mp_label)]
  if (length(ordered_mps_state) == 0) ordered_mps_state <- unique(top_mp_label)
  mp_subset_cols <- setNames(scales::hue_pal()(length(ordered_mps_state)), ordered_mps_state)

  # Find top samples by number of cells in this state, excluding Part-A samples
  sample_counts <- table(tmdata_all$orig.ident[state_cells])
  sample_counts <- sort(sample_counts, decreasing = TRUE)
  candidate_samples <- names(sample_counts)[sample_counts >= 30]
  candidate_samples <- setdiff(candidate_samples, used_partB_samples)
  top12_state_samples <- head(candidate_samples, 12)
  used_partB_samples <- c(used_partB_samples, top12_state_samples)

  message(sprintf("  %d total cells, top %d samples: %s",
                  length(state_cells), length(top12_state_samples),
                  paste(top12_state_samples, collapse = ", ")))

  out_dir <- file.path(out_root, "partB", state_dir_map[[state_name]])
  pdf(file.path(out_dir, paste0("Auto_", task_prefix, "_partB_", state_dir_map[[state_name]], "_top12_pseudotime.pdf")),
      width = 14, height = 6, onefile = TRUE)

  for (j in seq_along(top12_state_samples)) {
    sample_id <- top12_state_samples[j]
    message(sprintf("  Part B [%s] %d/%d: %s", state_name, j, length(top12_state_samples), sample_id))

    sample_state_cells <- state_cells[tmdata_all$orig.ident[state_cells] == sample_id]
    if (length(sample_state_cells) < 30) {
      message(sprintf("    Skipping: only %d cells", length(sample_state_cells)))
      next
    }

    sub_obj <- tmdata_all[, sample_state_cells]
    sample_mp_levels <- top_mp_label[sample_state_cells]
    sub_obj$mp_label <- factor(sample_mp_levels, levels = ordered_mps_state)

    mp_counts_sample <- table(factor(sample_mp_levels, levels = ordered_mps_state))
    mp_desc <- ifelse(
      is.na(mp_descriptions[ordered_mps_state]),
      ordered_mps_state,
      mp_descriptions[ordered_mps_state]
    )
    mp_legend_labels <- setNames(
      paste0(ordered_mps_state, " ", mp_desc, " (", as.integer(mp_counts_sample[ordered_mps_state]), ")"),
      ordered_mps_state
    )

    # Root cells = those assigned to root_mp
    root_cells <- sample_state_cells[top_mp_label[sample_state_cells] == root_mp]
    if (length(root_cells) == 0) {
      message(sprintf("    Skipping: no %s cells for root", root_mp))
      next
    }

    result <- tryCatch({
      run_pseudotime(
        seurat_obj = sub_obj,
        group_col = "mp_label",
        root_cells = root_cells,
        sample_id = sample_id,
        left_title_prefix = state_name,
        colour_map = mp_subset_cols,
        legend_labels = mp_legend_labels,
        legend_title = "MP",
        n_cells = length(sample_state_cells)
      )
    }, error = function(e) {
      message(sprintf("    Error: %s", e$message))
      NULL
    })

    if (is.null(result)) next

    # add one page per sample to single state PDF
    print(result$plots$group + result$plots$pseudotime)

    saveRDS(result$pseudotime, file.path(out_dir, paste0(sample_id, "_pseudotime_", state_name, ".rds")))

    partB_summary[[paste(state_name, sample_id, sep = "_")]] <- data.frame(
      part = "B",
      state_subset = state_name,
      sample = sample_id,
      rank = j,
      n_cells = length(sample_state_cells),
      n_root = length(root_cells),
      median_pseudotime = median(result$pseudotime, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    message(sprintf("    Done: %d cells, median PT = %.2f",
                    length(sample_state_cells),
                    median(result$pseudotime, na.rm = TRUE)))
  }

  dev.off()
}

####################
# Save summary CSV
####################
summary_A <- if (length(partA_summary) > 0) bind_rows(partA_summary) else data.frame()
summary_B <- if (length(partB_summary) > 0) bind_rows(partB_summary) else data.frame()
summary_all <- bind_rows(summary_A, summary_B)

write.csv(summary_all, file.path(summary_dir, paste0("Auto_", task_prefix, "_pseudotime_states_summary.csv")), row.names = FALSE)

message("\nPseudotime analysis complete.")
message(sprintf("Part A: %d samples processed", nrow(summary_A)))
message(sprintf("Part B: %d sample-state combinations processed", nrow(summary_B)))
