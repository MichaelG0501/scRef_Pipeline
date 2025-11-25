library(Seurat)
library(dplyr)
library(tidyr)
library(parallel)
library(Matrix)
library(readxl)
library(ggplot2)
library(patchwork)

reticulate::use_condaenv("dmtcp", conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda", required = TRUE)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

tmdata_annotated <- list()

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^[^/]+_[^/]+_[^/]+$", sample_dirs)]

marker_columns <- c("epithelial", "endothelial", "fibroblast", "b.cell", "t.cell", "macrophage",
                    "dendritic", "mast", "neutrophil", "erythrocyte", "keratinocyte", "lymph")

for (sample in sample_dirs) {
  tmdata <- readRDS(file.path("by_samples", sample, paste0(sample, "_anno.rds")))
  tmdata@meta.data[, marker_columns] <- NULL
  ct <- tmdata@meta.data$celltype_update
  tmdata@meta.data$celltype_update[grepl("plasma\\|", ct, ignore.case = TRUE)] <- "plasma"
  tmdata@meta.data$celltype_update[grepl("b\\.cell\\|", ct, ignore.case = TRUE)] <- "b.cell"
  tmdata@meta.data$celltype_update[grepl("nk\\.cell\\|", ct, ignore.case = TRUE)] <- "nk.cell"
  tmdata@meta.data$celltype_update[grepl("t\\.cell\\|", ct, ignore.case = TRUE)] <- "t.cell"
  tmdata_annotated[[sample]] <- tmdata
}

score_markers <- function(tmdata, top_n = 10) {
  expr <- as.matrix(tmdata@assays$RNA$data)
  tp_markers_list <- markers_list
  for (marker in names(tp_markers_list)) {
    tp_markers_list[[marker]] <- tp_markers_list[[marker]][tp_markers_list[[marker]] %in% rownames(expr)]
  }
  tp_markers_list$housekeeping <- NULL
  ct_mtx <- matrix(0, ncol = length(tp_markers_list), nrow = ncol(expr),
                   dimnames = list(colnames(expr), names(tp_markers_list)))
  for (celltype in names(tp_markers_list)) {
    marker_genes <- tp_markers_list[[celltype]]
    marker_expr <- expr[marker_genes, , drop = FALSE]  # keep matrix even if 1 row
    top_mean <- apply(marker_expr, 2, function(v) {
      k <- min(top_n, length(v))
      mean(sort(v, decreasing = TRUE)[seq_len(k)])
    })
    ct_mtx[, celltype] <- top_mean
  }
  return(ct_mtx)
}

fibroblast <- c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN")
macrophage <- c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68")
mast <- c("MS4A2", "CPA3", "TPSB2", "TPSAB1")
epithelial <- c("KRT7", "MUC1", "KRT19", "EPCAM")
t.cell <- c("CD3E", "CD3D", "CD2", "CD3G") #, "PTPRC")
b.cell <- c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1")
nk.cell <- c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1")
plasma <- c("MZB1", "JCHAIN", "DERL3")
dendritic <- c("CLEC10A", "CCR7", "CD86")
endothelial <- c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5")

lymph <- c("CCL21")
erythrocyte <- c("HBA1", "HBA2", "HBB")
keratinocyte <- c("FLG", "IVL")
neutrophil <- c("CTSG", "ELANE", "MPO", "AZU1")

markers_list <- list(
  erythrocyte = erythrocyte, 
  keratinocyte = keratinocyte,
  lymph = lymph, 
  neutrophil = neutrophil, 
  endothelial = endothelial,
  epithelial = epithelial,
  fibroblast = fibroblast, 
  b.cell = b.cell, 
  plasma = plasma, 
  dendritic = dendritic, 
  macrophage = macrophage,
  mast = mast,
  nk.cell = nk.cell, 
  t.cell = t.cell
)

for (name in names(tmdata_annotated)) {
  score_all <- score_markers(tmdata_annotated[[name]], top_n = 10)
  score_all[is.na(score_all)] <- 0
  high_marker_count <- rowSums(score_all > 1, na.rm = TRUE)
  high_markers <- score_all > 1
  
  non_exclusive_groups <- list(
    c("macrophage", "mast", "dendritic"),
    c("t.cell", "nk.cell", "dendritic"), 
    c("b.cell", "plasma")
  )
  
  special_case <- apply(high_markers, 1, function(cell_row) {
    high_ct <- names(which(cell_row))
    if (length(high_ct) > 0) {
      is_non_exclusive <- any(sapply(non_exclusive_groups, function(group) {
        all(high_ct %in% group) && length(high_ct) <= length(group)
      }))
      is_non_exclusive
    } else {
      FALSE
    }
  })
  
  tmdata_annotated[[name]]@meta.data$coexpression <- ifelse(high_marker_count > 1 & !special_case, "doublet", "singlet")
  
  top_two_diff <- apply(score_all, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    top_markers <- names(x)[ord]
    top1 <- top_markers[1]
    conflicting_groups <- non_exclusive_groups[sapply(non_exclusive_groups, function(g) top1 %in% g)]
    if (length(conflicting_groups) > 0) {
      conflicting_markers <- unique(unlist(conflicting_groups))
      other_idx <- which(!top_markers %in% conflicting_markers)[1]
      diff <- x[ord[1]] - x[top_markers[other_idx]]
    } else {
      diff <- x[ord[1]] - x[ord[2]]
    }
    diff
  })
  
  tmdata_annotated[[name]]@meta.data$coexpression_loose <- ifelse((high_marker_count > 1) & (top_two_diff < 1), "doublet", "singlet")
  
  resolve_marker_col <- function(i) {
    mc <- tmdata_annotated[[name]]@meta.data[i, "celltype_update"]
    sp <- list(c("t.cell","nk.cell"), c("nk","t"), c("plasma","b.cell"), c("b","plasma"))
    p <- sp[sapply(sp, function(x) mc == paste(x, collapse="|"))]
    if (length(p)) names(which.max(score_all[i, p[[1]], drop=TRUE])) else mc
  }
  
  tmdata_annotated[[name]]@meta.data$marker_expression <- mapply(
    function(celltype, row_index) {
      row_vals <- score_all[row_index, , drop = FALSE]
      max_all <- max(as.numeric(row_vals), na.rm = TRUE)
      
      if (celltype != "unresolved") {
        marker_col <- resolve_marker_col(row_index)
        if (!is.na(marker_col) &&
            marker_col %in% colnames(score_all) &&
            score_all[row_index, marker_col] > 2) {
          "good"
        } else if (max_all > 2) {
          tmdata_annotated[[name]]@meta.data[row_index, "celltype_update"] <<- "unresolved_inconsistent"
          "good_inconsistent"
        } else {
          "poor"
        }
      } else {
        if (max_all > 2) "good_unresolved" else "poor_unresolved"
      }
    },
    tmdata_annotated[[name]]@meta.data$celltype_update,
    seq_len(nrow(tmdata_annotated[[name]]@meta.data))
  )
  
  print(name)
}

allowed_pairs <- list(
  c("t.cell", "nk.cell"),
  c("nk.cell", "t.cell"),
  c("b.cell", "plasma"),
  c("plasma", "b.cell")
)
gap_cut <- 0.4

for (name in names(tmdata_annotated)) {
  
  print(paste0("checking good_xx from sample: ", name))
  mask <- tmdata_annotated[[name]]@meta.data$coexpression == "singlet" &
    !(tmdata_annotated[[name]]@meta.data$marker_expression %in% c("poor", "poor_unresolved"))
  
  marker_expression_table <- table(
    tmdata_annotated[[name]]@meta.data$seurat_clusters[mask],
    tmdata_annotated[[name]]@meta.data$marker_expression[mask]
  )
  cluster_tbc <- integer(0)
  if ("good_inconsistent" %in% colnames(marker_expression_table)) {
    frac_inconsistent <- marker_expression_table[, "good_inconsistent"] / rowSums(marker_expression_table)
    cluster_tbc <- which(frac_inconsistent > 0.1 & rowSums(marker_expression_table) >= 100)
  }
  
  if ("good_unresolved" %in% colnames(marker_expression_table)) {
    cluster_tbc <- c(cluster_tbc, which(marker_expression_table[, "good_unresolved"] > 50))
  }
  cluster_tbc <- unique(cluster_tbc)
  
  if (length(cluster_tbc) == 0) next
  for (cluster in cluster_tbc) {
    tmdata <- subset(tmdata_annotated[[name]], subset = seurat_clusters == cluster & coexpression == "singlet")
    ncells <- ncol(tmdata)
    tmdata <- FindVariableFeatures(tmdata, nfeatures = 3000)
    tmdata <- ScaleData(tmdata, verbose = FALSE)
    ncells <- ncol(tmdata)
    nfeatures <- nrow(tmdata)
    max_pcs <- min(ncells, nfeatures, 50)
    tmdata <- RunPCA(tmdata, npcs = max_pcs, verbose = FALSE)
    pc_use <- 1:min(30, max_pcs)
    k_use  <- min(30, max(10, round(sqrt(ncells) - 2)))
    tmdata <- FindNeighbors(tmdata, dims = pc_use, k.param = k_use, verbose = FALSE)
    tmdata <- FindClusters(tmdata, resolution = 0.8, algorithm = 4, verbose = FALSE)
    tmdata <- RunUMAP(tmdata, dims = pc_use, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
    
    available_genes <- rownames(tmdata)
    markers_in_data <- lapply(markers_list, function(gene_set) intersect(gene_set, available_genes))
    markers_in_data <- Filter(function(v) length(v) > 0, markers_in_data)
    ct_names <- names(markers_in_data)
    tmdata <- AddModuleScore(tmdata, features = unname(markers_in_data), name = "mod_")
    mod_cols <- paste0("mod_", seq_along(markers_in_data)); names(mod_cols) <- ct_names
    
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
    ct <- tmdata@meta.data$celltype_update
    tmdata@meta.data$celltype_update[grepl("plasma\\|", ct, ignore.case = TRUE)] <- "plasma"
    tmdata@meta.data$celltype_update[grepl("b\\.cell\\|", ct, ignore.case = TRUE)] <- "b.cell"
    tmdata@meta.data$celltype_update[grepl("nk\\.cell\\|", ct, ignore.case = TRUE)] <- "nk.cell"
    tmdata@meta.data$celltype_update[grepl("t\\.cell\\|", ct, ignore.case = TRUE)] <- "t.cell"
    score_all <- score_markers(tmdata, top_n = 10)
    
    resolve_marker_col_s <- function(i) {
      mc <- tmdata@meta.data[i, "celltype_update"]
      sp <- list(c("t.cell","nk.cell"), c("nk","t"), c("plasma","b.cell"), c("b","plasma"))
      p <- sp[sapply(sp, function(x) mc == paste(x, collapse="|"))]
      if (length(p)) names(which.max(score_all[i, p[[1]], drop=TRUE])) else mc
    }
    
    tmdata@meta.data$marker_expression <- mapply(
      function(celltype, row_index) {
        row_vals <- score_all[row_index, , drop = FALSE]
        max_all <- max(as.numeric(row_vals), na.rm = TRUE)
        
        if (celltype != "unresolved") {
          marker_col <- resolve_marker_col_s(row_index)
          if (!is.na(marker_col) &&
              marker_col %in% colnames(score_all) &&
              score_all[row_index, marker_col] > 2) {
            "good"
          } else if (max_all > 2) {
            tmdata_annotated[[name]]@meta.data[row_index, "celltype_update"] <<- "unresolved_inconsistent"
            "good_inconsistent"
          } else {
            "poor"
          }
        } else {
          if (max_all > 2) "good_unresolved" else "poor_unresolved"
        }
      },
      tmdata@meta.data$celltype_update,
      seq_len(nrow(tmdata@meta.data))
    )
    
    tmdata_annotated[[name]]@meta.data[colnames(tmdata), c("celltype_update","marker_expression")] <- tmdata@meta.data[, c("celltype_update","marker_expression")]
  }
}

n_clusters <- 8
mclapply(names(tmdata_annotated), function(nm) {
  sample_dir <- file.path("by_samples/", nm)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
  saveRDS(tmdata_annotated[[nm]], file.path(sample_dir, paste0(nm, "_anno.rds")), compress = FALSE)
  NULL
}, mc.cores = n_clusters, mc.preschedule = FALSE)

###########################################################################
###########################################################################

for (name in names(tmdata_annotated)) {
  
  tmdata_annotated[[name]] <- tryCatch({
    tmp <- subset(tmdata_annotated[[name]],
                  subset = (marker_expression == "good" | is.na(marker_expression)) & coexpression == "singlet")
    if (ncol(tmp) <= 1) {
      message(name, ": ≤1 cell after subset — skipping")
      NULL
    } else {
      tmp
    }
  }, error = function(e) {
    message(name, ": Error during subset — skipping")
    NULL
  })
}

gc()

all_genes <- Reduce(union, lapply(tmdata_annotated, function(obj) {
  rownames(GetAssayData(obj, layer = "counts"))
}))

pad_matrix <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                               dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, zero_mat)
  }
  # Ensure the same row order
  mat <- mat[all_genes, , drop = FALSE]
  return(mat)
}

counts_list <- lapply(names(tmdata_annotated), function(id) {
  mat <- GetAssayData(tmdata_annotated[[id]], layer = "counts")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

cpm_list <- lapply(names(tmdata_annotated), function(id) {
  mat <- GetAssayData(tmdata_annotated[[id]], layer = "CPM")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

lognorm_list <- lapply(names(tmdata_annotated), function(id) {
  mat <- GetAssayData(tmdata_annotated[[id]], layer = "data")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

meta_list <- lapply(names(tmdata_annotated), function(id) {
  meta <- tmdata_annotated[[id]]@meta.data[, c(names(tmdata_annotated[[id]]@meta.data)[1:4], "celltype_update")]
  rownames(meta) <- paste(id, rownames(meta), sep = "_")
  return(meta)
})

combined_counts <- do.call(cbind, counts_list)
combined_cpm <- do.call(cbind, cpm_list)
combined_lognorm <- do.call(cbind, lognorm_list)
combined_meta <- do.call(rbind, meta_list)

merged_obj <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)

merged_obj@assays$RNA$CPM <- combined_cpm
merged_obj@assays$RNA$data <- combined_lognorm
rm(tmdata_annotated)

merged_obj@meta.data$study <- sapply(strsplit(merged_obj@meta.data$orig.ident, "_"), function(x) paste(x[1:2], collapse = "_"))

merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- FindNeighbors(merged_obj, dims = 1:20)

merged_obj <- FindClusters(merged_obj, resolution = 0.8, algorithm = 1)
merged_obj$leiden_clusters <- Idents(merged_obj)

merged_obj <- RunUMAP(merged_obj, dims = 1:20)

library(readxl)
data <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
data <- data %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))
cell_names <- colnames(merged_obj)
merged_obj@meta.data <- merged_obj@meta.data %>%
  left_join(data, by = "orig.ident")
rownames(merged_obj@meta.data) <- cell_names

saveRDS(merged_obj, paste0("EAC_Ref_merged.rds"))

p1 <- DimPlot(merged_obj, group.by = "leiden_clusters", label = TRUE) + ggtitle("Louvain Clustering")
p2 <- DimPlot(merged_obj, group.by = "study", label = FALSE) + ggtitle("Study")
p3 <- DimPlot(merged_obj, group.by = "celltype_update", label = FALSE) + ggtitle("Celltype")
p4 <- DimPlot(merged_obj, group.by = "Treatment", label = FALSE) + ggtitle("Treatment Status")

combined_plot <- (p1 | p2) / (p3 | p4)

ggsave("EAC_Ref_filtered.png", plot = combined_plot, width = 12, height = 6, dpi = 300)
