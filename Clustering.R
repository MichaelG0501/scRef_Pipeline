library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(Seurat)
library(ggplot2)
library(purrr)
library(tidyr)
library(readxl)
library(stringr)

reticulate::use_condaenv("dmtcp", conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda", required = TRUE)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
print(sample)

tmdata <- readRDS(paste0("by_samples/", sample, "/", sample, ".rds"))

ncells <- ncol(tmdata)

if (ncells > 50) {
  tmdata <- FindVariableFeatures(tmdata, nfeatures = 3000)
  tmdata <- ScaleData(tmdata, verbose = FALSE)
  tmdata <- RunPCA(tmdata, npcs = 50, verbose = FALSE)
  
  k_use  <- min(30, max(10, round(sqrt(ncells) - 2)))
  pc_use <- 1:30
  tmdata <- FindNeighbors(tmdata, dims = pc_use, k.param = k_use, verbose = FALSE)
  tmdata <- FindClusters(tmdata, resolution = 0.8, algorithm = 4, verbose = FALSE)
  
  tmdata <- RunUMAP(tmdata, dims = pc_use, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
} else {
  tmdata <- FindVariableFeatures(tmdata, nfeatures = 3000, verbose = FALSE)
  tmdata <- ScaleData(tmdata, features = VariableFeatures(tmdata), verbose = FALSE)
  
  max_pcs <- min(length(VariableFeatures(tmdata)), ncells - 1)
  npcs    <- max(2, max_pcs)
  tmdata  <- RunPCA(tmdata, npcs = npcs, verbose = FALSE)
  pc_use <- 1:npcs
  k_use  <- max(5, min(15, ncells - 1))
  tmdata <- FindNeighbors(tmdata, dims = pc_use, k.param = k_use, verbose = FALSE)
  tmdata <- FindClusters(tmdata, resolution = 0.8, algorithm = 4, verbose = FALSE)
  
  umap_neighbors <- max(5, min(15, ncells - 1))
  tmdata <- RunUMAP(tmdata, dims = pc_use, n.neighbors = umap_neighbors, min.dist = 0.3, verbose = FALSE)
}

markers <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx", sheet = 1)
celltype_map <- c(
  "Fibroblast"  = "fibroblast",
  "Macrophage"  = "macrophage",
  "Mast"        = "mast",
  "B cell"      = "b.cell",
  "T cell"      = "t.cell",
  "Dendritic"   = "dendritic",
  "Endothelial" = "endothelial",
  "Epithelial"  = "epithelial",
  "NK cell"     = "nk.cell",
  "Plasma"      = "plasma"
)

combine_marker_scores <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }
  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)
  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

markers <- markers[markers$specificity > 0.2 & markers$cell_type != "Malignant", ]

markers_list <- markers %>%
  mutate(cell_type = recode(cell_type, !!!celltype_map)) %>%
  split(.$cell_type)

markers_ranked <- lapply(markers_list, function(df) {
  combine_marker_scores(df, w_specificity = 0.2, w_sensitivity = 0.8)
})


N <- 100
lfc_th   <- 1
z_cut <- 1.75
z_min <- 1
diffone <- 4

Idents(tmdata) <- tmdata$seurat_clusters

setsN <- markers_ranked |>
  imap(~ .x %>% arrange(desc(Combined)) %>% slice_head(n = N) %>% pull(gene) %>% intersect(rownames(tmdata)))

if (ncells > 50) {
  
  de <- FindAllMarkers(tmdata, only.pos = TRUE, min.pct = 0.8, logfc.threshold = lfc_th)
  universe <- rownames(tmdata)
  
  enrich_one <- function(cl){
    de_genes <- de %>% filter(cluster == cl) %>% pull(gene) %>% unique()
    
    res <- imap_dfr(setsN, ~{
      a <- length(intersect(.x, de_genes)); b <- length(setdiff(de_genes, .x))
      c <- length(setdiff(.x, de_genes));   d <- length(universe) - a - b - c
      p <- fisher.test(matrix(c(a,b,c,d), 2, 2), alternative = "greater")$p.value
      tibble(cell_type = .y, overlap = a, pval = p)
    }) %>%
      mutate(padj = p.adjust(pval, "BH"), cluster = cl)
    
    # keep only significant
    sig <- res %>%
      filter(padj <= 0.05) %>%
      arrange(desc(overlap), padj)
    
    # no significant → unknown
    if (nrow(sig) == 0) {
      return(tibble(cluster = cl, step1 = "unknown"))
    }
    
    # top overlap
    top_ov <- sig$overlap[1]
    
    # ---- CASE 1: only 1 significant ----
    if (nrow(sig) == 1) {
      if (top_ov > 5) {
        return(tibble(cluster = cl, step1 = sig$cell_type[1]))
      } else {
        return(tibble(cluster = cl, step1 = "unknown"))
      }
    }
    
    # ---- CASE 2: >= 2 significant ----
    second_ov <- sig$overlap[2]
    
    if ((top_ov - second_ov) > diffone) {
      # top is clearly better → keep ONLY top, but only if >5
      if (top_ov > 5) {
        keep <- sig %>% slice(1)
      } else {
        return(tibble(cluster = cl, step1 = "unknown"))
      }
    } else {
      # multiple similar ones → keep all close to top, but only those with overlap > 5
      keep <- sig %>%
        filter((top_ov - overlap) <= diffone, overlap > diffone)
      if (nrow(keep) == 0) {
        return(tibble(cluster = cl, step1 = "unknown"))
      }
    }
    
    calls <- keep %>%
      arrange(desc(overlap), padj) %>%
      pull(cell_type) %>%
      unique()
    
    tibble(cluster = cl, step1 = paste(calls, collapse = "|"))
  }
  
  step1_calls <- unique(de$cluster) %>% map_dfr(enrich_one)
  
  # ----------------- Step-2: AddModuleScore on top-M specificity -----------------
  ct_names <- names(setsN)
  
  tmdata <- AddModuleScore(tmdata, features = unname(setsN), name = "mod_")
  mod_cols <- paste0("mod_", seq_along(setsN)); names(mod_cols) <- ct_names
  
  scores_long <- tmdata@meta.data %>%
    mutate(cluster = tmdata$seurat_clusters) %>%
    group_by(cluster) %>%
    summarize(across(all_of(mod_cols), median, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
    mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
    group_by(cluster) %>%
    mutate(z = as.numeric(scale(score))) %>%  # z-score across celltypes in this cluster
    ungroup() %>%
    select(cluster, cell_type, score, z)
  
  step2_calls <- scores_long %>%
    group_by(cluster) %>%
    summarize(
      cell_types = list(cell_type),
      zs         = list(z),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      # make them explicit list-columns in the rowwise context
      all_z  = list(unlist(zs)),
      all_ct = list(unlist(cell_types)),
      max_idx = which.max(all_z),
      top_ct  = all_ct[max_idx],
      top_z   = all_z[max_idx],
      good_idx = list(which(all_z >= z_min)),
      n_good   = length(good_idx),
      step2 = {
        az <- all_z
        # 1) all < 0.5 → unknown
        if (all(az < 0.5, na.rm = TRUE)) {
          "unknown"
        } else if (n_good == 1 && top_z >= z_cut && all(az[-max_idx] < z_min, na.rm = TRUE)) {
          # 2) single strong and everyone else low → single
          top_ct
        } else if (n_good >= 2) {
          # 3) multiple above z_min → join them
          paste(all_ct[good_idx], collapse = "|")
        } else {
          "unknown"
        }
      }
    ) %>%
    ungroup() %>%
    select(cluster, step2)
  
  # ------------------- COMBINE STEP 1 + STEP 2 ----------------------
  # step1 is like: cluster | step1="Bcell|Tcell" OR "unknown"
  # step2 is like: cluster | step2="Bcell" or "Bcell|Tcell" or "unknown"
  
  final_calls <- step1_calls %>%
    left_join(step2_calls, by = "cluster") %>%
    rowwise() %>%
    mutate(
      final = {
        s1 <- step1
        s2 <- step2
        
        # helper: is single?
        is_single <- function(x) !is.na(x) && x != "unknown" && !grepl("\\|", x)
        
        if (!is.na(s1) && !is.na(s2) && s1 == s2) {
          # both same
          s1
        } else if (is_single(s1) && (is.na(s2) || s2 == "unknown")) {
          # step1 single, step2 unknown
          s1
        } else if (is_single(s2) && (is.na(s1) || s1 == "unknown")) {
          # step2 single, step1 unknown
          s2
        } else {
          # either multi, or different singles, or one is unknown but other is multi
          paste0(s1, " <> ", s2)
        }
      }
    ) %>%
    ungroup()
  
  # attach to object
  lab_map <- final_calls$final; names(lab_map) <- final_calls$cluster
  tmdata$celltype_initial <- as.character(lab_map[as.character(tmdata$seurat_clusters)])
  
} else {
  tmdata$celltype_initial <- rep("unknown", ncol(tmdata))
}

saveRDS(tmdata, paste0("by_samples/", sample, "/", sample, ".rds"))