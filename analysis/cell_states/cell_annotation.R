####################
# Moved from: analysis/annotation.R
# Reorganized as part of analysis/ restructuring
####################
library(Seurat)
library(dplyr)
library(purrr)


sample_dirs <- list.dirs(path = "./by_samples", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]  # match *_*_*

tmdata_list <- list()

for (sample in sample_dirs) {
  rds_path <- file.path("/rds/general/project/tumourheterogeneity1/ephemeral/EAC_Ref_all/00_merged/by_samples", sample, paste0(sample, ".rds"))
  if (file.exists(rds_path)) {
    tmdata_list[[sample]] <- readRDS(rds_path)
    if (ncol(tmdata_list[[sample]]) < 50) {
      tmdata_list[[sample]]$celltype_final <- rep("unresolved", ncol(tmdata_list[[sample]]))
    }
  } else {
    warning(paste("Missing RDS file for sample:", sample))
  }
}


# for (sample in sample_dirs) {
#   obj <- tmdata_list[[sample]]
#   id_counts <- table(Idents(obj))
#   
#   # Find identities with fewer than 30 cells
#   low_ids <- names(id_counts[id_counts < 50])
#   
#   # Update celltype_final for those cells
#   cells_to_update <- WhichCells(obj, idents = low_ids)
#   obj$celltype_final[cells_to_update] <- paste0("unresolved_", obj$celltype_final[cells_to_update])
#   
#   # Save back to list
#   tmdata_list[[sample]] <- obj
# }


ct_names <- names(markers_list)
tmdata <- AddModuleScore(tmdata, features = unname(markers_list), name = "mod_")
mod_cols <- paste0("mod_", seq_along(markers_list)); names(mod_cols) <- ct_names




Idents(tmdata) <- tmdata$celltype_final

ct_names <- names(setsN)
all_cts  <- unique(as.character(tmdata$celltype_final))
ordered_cts <- c(intersect(ct_names, all_cts), setdiff(all_cts, ct_names))
tmdata$celltype_final <- factor(tmdata$celltype_final, levels = ordered_cts)

# 1) Keep only markers present in the object; drop empty sets
features_list <- setsN %>%
  imap(~ intersect(.x, rownames(tmdata))) %>%
  keep(~ length(.x) > 0)


K <- 10
features_list <- lapply(features_list, function(v) head(v, K))

p <- DotPlot(
  tmdata,
  features = features_list,
  group.by = "celltype_final",
  scale = TRUE,            # color: scaled average expression (per gene)
  dot.scale = 5
) +
  RotatedAxis() +
  ggplot2::labs(
    x = NULL, y = NULL,
    color = "Scaled avg expr", size = "% expressing"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.major = ggplot2::element_line(size = 0.2, color = "grey90")
  )

print(p)

# 3) Save (optional)
# ggsave("marker_bubbleplot.pdf", p, width = 10, height = 12)
# ggsave("marker_bubbleplot.png", p, width = 10, height = 12, dpi = 300)






library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# ---- params ----
gene <- "ELF3"            # change as needed
ct_col <- "celltype_final" # column in meta.data with your final cell types
min_cells_per_group <- 0   # drop study×celltype groups with < this many cells

# ---- summarize one Seurat object ----
summarize_one <- function(obj, study_id, gene, ct_col) {
  if (!gene %in% rownames(obj)) {
    warning(paste0("Skipping ", study_id, " (gene '", gene, "' not found)."))
    return(NULL)
  }
  meta <- obj@meta.data
  if (!ct_col %in% colnames(meta)) {
    stop(paste0("'", ct_col, "' not found in meta.data for study: ", study_id))
  }
  df <- FetchData(obj, vars = c(gene, ct_col))
  colnames(df) <- c("expr", "celltype")
  df %>%
    group_by(celltype) %>%
    summarize(
      n = n(),
      pct = 100 * mean(expr > 0, na.rm = TRUE),
      avg = mean(expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n >= min_cells_per_group) %>%
    mutate(study = study_id)
}

sumdf <- imap_dfr(tmdata_list, ~ summarize_one(.x, .y, gene, ct_col))
sumdf <- sumdf[sumdf$celltype %in% names(setsN), ]

ct_order <- sumdf %>%
  group_by(celltype) %>%
  summarize(mean_avg = mean(avg, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_avg)) %>%
  pull(celltype)

sumdf$celltype <- factor(sumdf$celltype, levels = names(setsN))
sumdf$study    <- factor(sumdf$study, levels = names(tmdata_list))  # preserves your list order

# ---- bubble plot ----
p <- ggplot(sumdf, aes(x = study, y = celltype)) +
  geom_point(aes(size = pct,
                 color = avg)) +
  scale_size(name = "% expressing", range = c(0.5, 8), limits = c(0, 100)) +
  labs(
    x = NULL, y = NULL,
    title = paste0("Bubble plot for ", gene, " across studies"),
    color = "Avg expr"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(size = 0.2, color = "grey90")
  )

print(p)

# ---- optional: save ----
# ggsave(paste0("bubble_by_study_", gene, ".pdf"), p, width = 10, height = 6)
# ggsave(paste0("bubble_by_study_", gene, ".png"), p, width = 10, height = 6, dpi = 300)



library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

method <- "celltype_update"

celltypes <- names(setsN)
all_genes <- unique(unlist(setsN))

# ---- 0. helper: exclusion rules -----------------------------------------
# This returns which OFF-celltypes to ignore, given the *home* celltype
get_excluded_off_cts <- function(home_ct) {
  # base blank
  excluded <- character(0)
  
  # 1) B-cell markers: ignore plasma
  # adjust the pattern to your actual plasma name, e.g. "Plasma", "Plasmablast", etc.
  if (grepl("^b.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "plasma")
  }

  # 2) T-cell markers: ignore NK
  if (grepl("^plasma", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "b.cell")
  }

  if (grepl("^t.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "nk.cell")
  }

  if (grepl("^nk.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "t.cell")
  }
  
  # 3) --- add new rules here ---
  # e.g.
  # if (home_ct == "Myeloid") excluded <- c(excluded, "DC", "Mono")
  
  excluded
}

# ---- 1. Compute % expression for each sample × celltype × gene ----
pct_tbl <- imap_dfr(tmdata_list, function(obj, sid) {
  meta <- obj@meta.data
  genes <- intersect(all_genes, rownames(obj))
  expr  <- GetAssayData(obj, slot = "data")[genes, , drop = FALSE]
  meta <- meta %>% pull(!!sym(method))
  split_cells <- split(colnames(obj), meta)
  
  map_dfr(names(split_cells), function(ct) {
    cells <- split_cells[[ct]]
    if (length(cells) == 0) return(NULL)
    pct <- Matrix::rowMeans(expr[, cells, drop = FALSE] > 1) * 100
    tibble(sample = sid, celltype = ct, gene = genes, pct_expr = pct)
  })
})

# keep only target celltypes
pct_tbl <- pct_tbl %>% filter(celltype %in% celltypes)
pct_tbl <- pct_tbl[pct_tbl$celltype != "dendritic", ]

expr_off_threshold   <- 15   # per-sample for OFF celltypes
expr_home_threshold  <- 30   # per-sample for HOME celltype
home_min_pct_samps   <- 15   # % of samples (home) that must pass
off_max_pct_samps    <- 30   # % of samples (off) that may pass


expr_off_threshold   <- 20   # per-sample for OFF celltypes
expr_home_threshold  <- 50   # per-sample for HOME celltype
home_min_pct_samps   <- 50   # % of samples (home) that must pass
off_max_pct_samps    <- 20   # % of samples (off) that may pass

# marker mapping: (home celltype, gene)
markers_long <- tibble(celltype = names(setsN), gene = setsN) %>% unnest(gene)

pct_with_homeflag <- pct_tbl %>%
  left_join(markers_long %>% mutate(is_home = TRUE),
            by = c("gene", "celltype")) %>%
  mutate(
    is_home = coalesce(is_home, FALSE),
    pass    = ifelse(
      is_home,
      pct_expr >= expr_home_threshold,   # stricter for home
      pct_expr >= expr_off_threshold     # 15% for off
    )
  )

gene_ct_summary <- pct_with_homeflag %>%
  group_by(gene, celltype, is_home) %>%
  summarise(
    n_total     = sum(!is.na(pct_expr)),
    n_pass      = sum(pass, na.rm = TRUE),
    pct_samples = ifelse(n_total > 0, 100 * n_pass / n_total, 0),
    .groups = "drop"
  )

gene_ct_long <- gene_ct_summary %>%
  select(gene, celltype, pct_samples, is_home)

# --- 2) Home % for each (gene, home celltype) from setsN ---
home_tbl <- markers_long %>%
  left_join(gene_ct_long, by = c("gene", "celltype")) %>%
  mutate(home_pct = coalesce(pct_samples, 0)) %>%
  select(gene, celltype, home_pct)

# --- 3) Max OFF-celltype % with exclusion rules ----------------------------
off_tbl_raw <- gene_ct_long %>%
  rename(off_celltype = celltype, off_pct = pct_samples) %>%
  inner_join(markers_long, by = "gene")  # brings in the *home* celltype as `celltype`

off_tbl_filtered <- off_tbl_raw %>%
  rowwise() %>%
  mutate(
    drop_these = list(get_excluded_off_cts(celltype)),
    keep_row   = !(off_celltype %in% drop_these)
  ) %>%
  ungroup() %>%
  filter(off_celltype != celltype, keep_row) %>%
  select(-drop_these, -keep_row)

off_tbl <- off_tbl_filtered %>%
  group_by(gene, celltype) %>%
  summarise(off_pct_max = max(off_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(off_pct_max = ifelse(is.finite(off_pct_max), off_pct_max, 0))

# --- 4) Exclusivity call ---------------------------------------------------
exclusive_tbl <- home_tbl %>%
  left_join(off_tbl, by = c("gene", "celltype")) %>%
  mutate(
    off_pct_max  = coalesce(off_pct_max, 0),
    is_exclusive = (home_pct >= home_min_pct_samps) & (off_pct_max <= off_max_pct_samps)
  ) %>%
  arrange(desc(is_exclusive), desc(home_pct), off_pct_max)

save <- exclusive_tbl %>%
  filter(is_exclusive) %>%
  select(gene, home_celltype = celltype, home_pct, off_pct_max)


save <- save %>%
  mutate(score = home_pct - off_pct_max) %>%
  group_by(home_celltype) %>%
  mutate(score_rank = rank(-score, ties.method = "first")) %>%
  slice_max(order_by = score, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(home_celltype, desc(score))

markers_list <- split(save$gene, save$home_celltype)

markers_list$housekeeping <- housekeeping
markers <- unlist(markers_list)


#z_cut <- 1
#z_min <- 0.4
gap_cut <- 0.8

# allowed "close" pairs
allowed_pairs <- list(
  c("t.cell", "nk.cell"),
  c("nk.cell", "t.cell"),
  c("b.cell", "plasma"),
  c("plasma", "b.cell")
)

# to collect per-sample % tables
pct_list   <- list()
count_list <- list()

for (sample in sample_dirs) {
  
  tmdata <- tmdata_list[[sample]]
  
  Idents(tmdata) <- tmdata$seurat_clusters
  clusters <- levels(Idents(tmdata))
  available_genes <- rownames(tmdata)
  
  markers_in_data <- lapply(markers_list, function(gene_set) intersect(gene_set, available_genes))
  markers_in_data <- Filter(function(v) length(v) > 0, markers_in_data)
  ct_names <- names(markers_in_data)
  
  score_mat <- matrix(NA_real_,
                      nrow = ncol(tmdata),
                      ncol = length(ct_names),
                      dimnames = list(colnames(tmdata), ct_names))
  
  for (cl in clusters) {
    cells_cl <- WhichCells(tmdata, idents = cl)
    if (length(cells_cl) == 0) next
    tm_sub <- subset(tmdata, cells = cells_cl)
    mtx <- GetAssayData(tm_sub, slot = "data")
    
    # for each cell type, pick top 8 expressed markers *within this cluster*
    cl_features <- lapply(markers_in_data, function(genes) {
      g <- intersect(genes, rownames(mtx))
      if (length(g) == 0) {
        character(0)
      } else {
        # mean expression across cells in this cluster
        m <- Matrix::rowMeans(mtx[g, , drop = FALSE])
        g[order(m, decreasing = TRUE)][seq_len(min(4, length(g)))]
      }
    })
    
    # keep only cell types with at least one expressed gene in this cluster
    keep <- vapply(cl_features, function(v) length(v) > 0, logical(1))
    if (!any(keep)) next
    cl_features_kept <- cl_features[keep]
    kept_ct <- names(cl_features_kept)
    
    tm_sub <- AddModuleScore(
      object   = tm_sub,
      features = cl_features_kept,
      name     = "mod_tmp_",
      assay    = DefaultAssay(tm_sub)
    )
    
    score_cols <- paste0("mod_tmp_", seq_along(cl_features_kept))
    scdf <- tm_sub@meta.data[, score_cols, drop = FALSE]
    colnames(scdf) <- kept_ct
    score_mat[cells_cl, kept_ct] <- as.matrix(scdf[rownames(scdf), kept_ct, drop = FALSE])
  }
  
  for (ct in ct_names) {
    tmdata@meta.data[[paste0("mod_", ct)]] <- score_mat[colnames(tmdata), ct]
  }
  
  mod_cols <- setNames(paste0("mod_", ct_names), ct_names)
  
  # 1) per-cluster median module scores and z
  scores_long <- tmdata@meta.data %>%
    mutate(cluster = tmdata$seurat_clusters) %>%
    group_by(cluster) %>%
    summarize(across(all_of(mod_cols), median, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
    mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
    group_by(cluster) %>%
    mutate(z = as.numeric(scale(score))) %>%   # z across celltypes for THIS cluster
    ungroup() %>%
    select(cluster, cell_type, score, z)
  
  # 2) step-2 call per cluster (new logic + pair exception)
  # final rule: step2 = { ... }
  step2_calls <- scores_long %>%
    group_by(cluster) %>%
    summarize(
      cell_types = list(cell_type),
      zs = list(z),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      all_z = list(unlist(zs)),
      all_ct = list(unlist(cell_types)),
      max_idx = which.max(all_z),
      top_ct  = all_ct[max_idx],
      top_z   = all_z[max_idx],
      # NEW: margin rule (top must be >= 0.8 higher than all others to call top_ct)
      step2 = {
        az  <- all_z
        act <- all_ct
        other_idx <- setdiff(seq_along(az), max_idx)
        
        if (length(other_idx) == 0) {
          # only one cell type scored in this cluster
          top_ct
        } else {
          margins <- top_z - az[other_idx]
          # "close" others are those within 0.8 of the top (i.e., margin < 0.8)
          close_idx <- other_idx[margins < 0.8]
          
          if (length(close_idx) == 0) {
            # top is >= 0.8 higher than ALL others -> call top_ct
            top_ct
          } else if (length(close_idx) == 1) {
            # exactly one close competitor -> allow pair if in allowed_pairs, else unresolved
            other_ct <- act[close_idx]
            pair_vec <- c(top_ct, other_ct)
            is_allowed <- any(vapply(allowed_pairs, function(p) identical(p, pair_vec), logical(1)))
            if (is_allowed) paste0(top_ct, "|", other_ct) else "unresolved"
          } else {
            # multiple close competitors -> unresolved
            "unresolved"
          }
        }
      }
    ) %>%
    ungroup() %>%
    select(cluster, step2)
  
  # 3) map back to cells
  cl_map <- step2_calls$step2
  names(cl_map) <- step2_calls$cluster
  
  tmdata$celltype_update <- as.character(cl_map[as.character(tmdata$seurat_clusters)])
  
  # write back
  tmdata_list[[sample]] <- tmdata
  
  # 4a) per-sample % table
  tab <- table(tmdata$celltype_update)
  pct <- 100 * tab / sum(tab)
  pct_df <- data.frame(
    study    = sample,
    celltype = names(pct),
    pct      = as.numeric(pct),
    stringsAsFactors = FALSE
  )
  pct_list[[sample]] <- pct_df
  
  # 4b) per-sample COUNT table (new)
  count_df <- data.frame(
    study    = sample,
    celltype = names(tab),
    count    = as.integer(tab),
    stringsAsFactors = FALSE
  )
  count_list[[sample]] <- count_df
  
  print(sample)
}

# ----------------- final matrices -----------------
library(dplyr)
library(tidyr)

# pct matrix (studies as rows, celltypes as cols)
pct_all <- bind_rows(pct_list)
pct_mat <- pct_all %>%
  tidyr::pivot_wider(
    names_from  = celltype,
    values_from = pct,
    values_fill = 0
  )

# count matrix (studies as rows, celltypes as cols)
count_all <- bind_rows(count_list)
count_mat <- count_all %>%
  tidyr::pivot_wider(
    names_from  = celltype,
    values_from = count,
    values_fill = 0
  )

# pct_mat and count_mat now both available