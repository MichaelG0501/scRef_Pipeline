####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/pseudotime_batch_correction.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_pseudotime_batch_correction.R
# Monocle3 pseudotime analysis for the noreg Approach B cell states,
# incorporating scVI and Harmony batch correction across ALL cells.
# Highly optimized for strict correct integration embedding carryover.
#
# Part A: Pseudotime using the 5 defined states (all cells)
#   - Root = Basal to Intestinal Metaplasia state cells
#
# Part B: Pseudotime within 3 state subsets (all cells in subset)
#   - Basal to Intestinal Metaplasia
#   - SMG-like Metaplasia
#   - Stress-adaptive
#   - Basal and SMG Metaplasia
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#   ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
####################

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(harmony)

reticulate::use_condaenv("dmtcp", conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda", required = TRUE)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

####################
# Constants
####################
task_prefix <- "task1_batch_correction"
out_root <- "pseudotime_batch_correction_states"

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
####################
message("Loading data ...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")

# Align cells
common_cells <- intersect(names(state_B), Cells(tmdata_all))
common_cells <- intersect(common_cells, rownames(mp_adj_noncc))
tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

tmdata_all$state_B <- state_B[Cells(tmdata_all)]

####################
# Pipeline functions
####################

run_pseudotime_batch_corrected <- function(
  seurat_obj, group_col, root_cells, batch_method, left_title_prefix,
  colour_map, legend_labels, legend_title, n_cells,
  out_dir_cds, out_dir_diag, save_prefix
) {
  message(sprintf("Running preprocess & batch correction: %s ...", batch_method))
  
  if (batch_method == "Harmony") {
    # Full preprocessing for Harmony
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    n_pcs <- min(30, ncol(seurat_obj) - 1)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE)
    dims_use <- min(15, n_pcs)
    
    if (length(Layers(seurat_obj[["RNA"]])) > 1) {
      seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
    }
    
    seurat_obj <- RunHarmony(
      seurat_obj, 
      group.by.vars = "study", 
      reduction = "pca",
      reduction.save = "harmony",
      assay.use = "RNA",
      theta = 4, # Tuned parameter for biological-preserving mixing
      plot_convergence = FALSE
    )
    # Strictly define downstream on integrated embedding
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:dims_use, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
    
    # Overwrite the default UMAP with the corrected one
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, reduction = "harmony", reduction.name = "UMAP", reduction.key = "UMAP_", verbose = FALSE)
    
  } else if (batch_method == "scVI") {
    
    # scVI prefers minimally-processed raw counts layer.
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    dims_use <- min(15, ncol(seurat_obj) - 1)
    
    if (length(Layers(seurat_obj[["RNA"]])) == 1) {
      seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$study)
    }
    
    # Using scVIIntegration via SeuratWrappers/Seurat v5
    seurat_obj <- IntegrateLayers(
      object = seurat_obj, 
      method = scVIIntegration,
      new.reduction = "scvi",
      batch_key = "study",
      conda_env = "dmtcp",
      verbose = FALSE
    )
    
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "scvi", dims = 1:dims_use, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
    
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, reduction = "scvi", reduction.name = "UMAP", reduction.key = "UMAP_", verbose = FALSE)
    
    seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  }
  
  # === CRITICAL STEP: UMAP OVERWRITE ===
  # To prevent as.cell_data_set() from picking up old / base UMAPs, 
  # strictly clear out anything not built specifically for this iteration.
  for (red in names(seurat_obj@reductions)) {
    if (red != "UMAP" && red != tolower(batch_method)) {
      seurat_obj[[red]] <- NULL
    }
  }
  
  # === DIAGNOSTICS GENERATION ===
  message("Generating and saving diagnostic QC embeddings ...")
  p_diag_study <- DimPlot(seurat_obj, group.by = "study", reduction = "UMAP", label = FALSE, alpha = 0.6) + 
    ggtitle(paste0(batch_method, " Integration by Study")) + 
    theme_minimal()
    
  p_diag_state <- DimPlot(seurat_obj, group.by = group_col, reduction = "UMAP", label = TRUE, repel = TRUE, alpha = 0.8) + 
    ggtitle(paste0(batch_method, " Integration by ", legend_title)) + 
    theme_minimal()
    
  pdf(file.path(out_dir_diag, paste0("Auto_diagnostics_", save_prefix, ".pdf")), width = 14, height = 6)
  print(p_diag_study | p_diag_state)
  dev.off()
  
  message("Saving intermediate integrated Seurat object ...")
  saveRDS(seurat_obj, file.path(out_dir_cds, paste0("Auto_", save_prefix, "_seurat_integrated.rds")))

  # === MONOCLE 3 PSEUDOTIME ===
  message("Converting to CDS and calculating Monocle3 pseudotime ...")
  
  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  
  root_cells_use <- intersect(root_cells, colnames(cds))
  if (length(root_cells_use) == 0) {
    stop("No root cells overlap with CDS columns")
  }
  
  cds <- order_cells(cds, root_cells = root_cells_use)
  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA
  
  p1 <- plot_cells(cds, color_cells_by = group_col, show_trajectory_graph = TRUE,
                   label_cell_groups = FALSE, label_groups_by_cluster = FALSE, 
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5) +
    labs(title = paste0(left_title_prefix, " via ", batch_method, " (n = ", n_cells, ")"), color = NULL) +
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
                   label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.5) +
    scale_color_viridis_c() +
    labs(title = paste0("Pseudotime (", batch_method, ")"), color = "Pseudotime") +
    theme_minimal()
    
  list(cds = cds, pseudotime = pt, plots = list(group = p1, pseudotime = p2))
}


####################
# Main Logic Loop over Batch Methods
####################
batch_methods <- c("Harmony", "scVI")

for (batch_method in batch_methods) {
  message(sprintf("\n===================================="))
  message(sprintf("=== Running Workflow for: %s ===", batch_method))
  message(sprintf("====================================\n"))
  
  out_root_method <- file.path(out_root, batch_method)
  
  # Create path layouts
  dir.create(file.path(out_root_method, "partA"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_root_method, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  
  for (d in state_dir_map) {
    dir.create(file.path(out_root_method, "partB", d), recursive = TRUE, showWarnings = FALSE)
  }
  
  # --- PART A ---
  message(sprintf("=== PART A: All-sample state pseudotime (%s) ===", batch_method))
  
  target_states <- names(state_groups)
  partA_cells <- names(state_B)[state_B %in% target_states]
  
  sub_obj_A <- tmdata_all[, partA_cells]
  sub_obj_A$state_label <- factor(state_B[partA_cells], levels = target_states)
  
  state_counts <- table(factor(as.character(state_B[partA_cells]), levels = target_states))
  state_legend_labels <- setNames(
    paste0(target_states, " (", as.integer(state_counts[target_states]), ")"),
    target_states
  )
  
  root_cells_A <- partA_cells[state_B[partA_cells] == "Basal to Intestinal Metaplasia"]
  
  result_A <- tryCatch({
    run_pseudotime_batch_corrected(
      seurat_obj = sub_obj_A,
      group_col = "state_label",
      root_cells = root_cells_A,
      batch_method = batch_method,
      left_title_prefix = "All Cells",
      colour_map = group_cols[target_states],
      legend_labels = state_legend_labels,
      legend_title = "State",
      n_cells = length(partA_cells),
      out_dir_cds = file.path(out_root_method, "partA"),
      out_dir_diag = file.path(out_root_method, "diagnostics"),
      save_prefix = paste0("partA_", batch_method)
    )
  }, error = function(e) {
    message(sprintf("  Error in Part A (%s): %s", batch_method, e$message))
    NULL
  })
  
  if (!is.null(result_A)) {
    pdf(file.path(out_root_method, "partA", sprintf("Auto_partA_%s_pseudotime.pdf", batch_method)), width = 14, height = 6)
    print(result_A$plots$group + result_A$plots$pseudotime)
    dev.off()
    
    # Save RDS output
    saveRDS(result_A$cds, file.path(out_root_method, "partA", sprintf("Auto_partA_%s_cds.rds", batch_method)))
  }
  
  # --- PART B ---
  message(sprintf("\n=== PART B: MP-level pseudotime within state subsets (%s) ===", batch_method))
  
  for (state_name in names(state_subsets)) {
    subset_info <- state_subsets[[state_name]]
    mps_in_state <- subset_info$mps
    root_mp <- subset_info$root_mp
    
    source_states <- subset_info$source_states
    if (is.null(source_states)) source_states <- state_name
    
    state_cells <- names(state_B)[state_B %in% source_states]
    if (length(state_cells) < 50) {
      message(sprintf("  Not enough cells for %s, skipping", state_name))
      next
    }
    
    mp_cols_avail <- intersect(mps_in_state, colnames(mp_adj_noncc))
    mp_sub <- mp_adj_noncc[state_cells, mp_cols_avail, drop = FALSE]
    top_mp_idx <- max.col(mp_sub, ties.method = "first")
    top_mp_label <- mp_cols_avail[top_mp_idx]
    names(top_mp_label) <- state_cells
    
    ordered_mps_state <- mps_in_state[mps_in_state %in% unique(top_mp_label)]
    if (length(ordered_mps_state) == 0) ordered_mps_state <- unique(top_mp_label)
    mp_subset_cols <- setNames(scales::hue_pal()(length(ordered_mps_state)), ordered_mps_state)
    
    sub_obj_B <- tmdata_all[, state_cells]
    sample_mp_levels <- top_mp_label[state_cells]
    sub_obj_B$mp_label <- factor(sample_mp_levels, levels = ordered_mps_state)
    
    mp_counts <- table(factor(sample_mp_levels, levels = ordered_mps_state))
    mp_desc <- ifelse(
      is.na(mp_descriptions[ordered_mps_state]),
      ordered_mps_state,
      mp_descriptions[ordered_mps_state]
    )
    mp_legend_labels <- setNames(
      paste0(ordered_mps_state, " ", mp_desc, " (", as.integer(mp_counts[ordered_mps_state]), ")"),
      ordered_mps_state
    )
    
    root_cells_B <- state_cells[top_mp_label[state_cells] == root_mp]
    
    out_dir_B <- file.path(out_root_method, "partB", state_dir_map[[state_name]])
    
    result_B <- tryCatch({
      run_pseudotime_batch_corrected(
        seurat_obj = sub_obj_B,
        group_col = "mp_label",
        root_cells = root_cells_B,
        batch_method = batch_method,
        left_title_prefix = state_name,
        colour_map = mp_subset_cols,
        legend_labels = mp_legend_labels,
        legend_title = "MP",
        n_cells = length(state_cells),
        out_dir_cds = out_dir_B,
        out_dir_diag = file.path(out_root_method, "diagnostics"),
        save_prefix = paste0("partB_", state_dir_map[[state_name]], "_", batch_method)
      )
    }, error = function(e) {
      message(sprintf("  Error in Part B %s (%s): %s", state_name, batch_method, e$message))
      NULL
    })
    
    if (!is.null(result_B)) {
      pdf(file.path(out_dir_B, sprintf("Auto_partB_%s_%s_pseudotime.pdf", state_dir_map[[state_name]], batch_method)), width = 14, height = 6)
      print(result_B$plots$group + result_B$plots$pseudotime)
      dev.off()
      
      saveRDS(result_B$cds, file.path(out_dir_B, sprintf("Auto_partB_%s_%s_cds.rds", state_dir_map[[state_name]], batch_method)))
    }
  }
}

message("\nBatch-corrected pseudotime analysis complete.")
