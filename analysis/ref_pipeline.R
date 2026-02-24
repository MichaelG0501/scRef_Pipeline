###############Loding required packages########################

library("Seurat")
library("dplyr")
library("Seurat")
library("DoubletFinder")
library("patchwork")
library("ggplot2")
library("foreach")
library("doParallel")
library("ComplexHeatmap")
library("circlize")
library("gridExtra")
library("grid")
library("tidyr")
library("tibble")
library("purrr")
library("harmony")

reticulate::use_condaenv("dmtcp", conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda", required = TRUE)

######################Setting parameters########################

# data <- read_excel("/rds/general/project/tumourheterogeneity1/ephemeral/EAC_Ref_all/00_merged/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
# data$study <- paste(data$Author, data$Year, data$`Sample Name`, sep = "_")
# tsample <- data$study[data$T_Status == "Tumour"]
# writeLines(tsample, "names_tmdata_EAC_Ref_t.txt")
# nsample <- data$study[data$T_Status != "Tumour"]
# writeLines(nsample, "names_tmdata_EAC_Ref_n.txt")

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/EAC_Ref_all/00_merged")
batch <- "_Carroll"
names_tmdata <- readLines(paste0("names_tmdata", batch, ".txt"))
max_mt = 25
min_ngenes = 300
min_hk_expr = 1.0
n_clusters = 8
celltyping_method = "celltype_manual"
color_scale <- colorRamp2(c(0, 6.5, 10), c("#D0D0D0", "red4", "red4"))
least_ncells = 3
ct_reorder <- c(
  "B_GC", "B_mature", "B_progenitor", "Plasma",    # b.cell
  "Dendritic_classical", "Dendritic_plasmacytoid", # dendritic
  "Endothelial",                                   # endothelial
  "Ductal", "Squamous", "Schwann", "Ciliated",     # epithelial
  "Erythroid",                                     # erythrocyte
  "Mural", "Fibroblast",                           # fibroblast
  "Macrophage", "Monocyte",                        # macrophage
  "Mast&Basophil",                                 # mast
  "Neutrophil",                                    # neutrophil
  "T&NK",                                          # t.cell
  "Neuron_excitatory", "Neuron_inhibitory", "Neuron_bipolar", # neuron
  "Astrocyte", "Muller", "Oligodendrocyte_mature", # glial
  "Hematopoietic", "Hepatocyte", "Melanocyte", "Rod", "Spermatocyte" # others
)

ct_reorder <- c("b.cell", "dendritic", "endothelial", "epithelial", 
                "erythrocyte", "fibroblast", "macrophage", "mast", 
                "neutrophil", "t.cell", "neuron", "glial", "others", "_others_")

ct_reorder <- c("b.cell", "dendritic", "endothelial", "epithelial", 
                "fibroblast", "macrophage", "mast", "nk.cell", "nk.cell|t.cell", 
                "plasma", "t.cell", "t.cell|nk.cell", "unresolved", "unreolved_inconsistent")

qc_rules <- data.frame(
  pattern = c(
    "Alcindor_2025",
    "Baek_2025",
    "Carroll_2023", 
    "Croft_2022",
    "Ju_2025",
    "Lambroia_2024",
    "Strasser_2025",
    "Walker_2025",
    "Wu_2018",
    "Yates_2025"
  ),
  mito   = c(25, 15, 25, 10, 40, 10, 10, 25, 1, 20),
  nGenes = c(300, 300, 300, 500, 300, 500, 500, 300, 5000, 300),
  hk     = c(1, 3, 2, 3, 2, 4, 3, 2, 5, 0.5),
  stringsAsFactors = FALSE
)

################################################################

initialise <- function(names) {
  
  tmdata_list <- list()
  for (name in names) {
    filename <- paste0("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_counts_matrix_all/", name, ".csv")
    tmdata <- data.table::fread(filename)
    Genes <- tmdata[[1]]
    counts <- as.matrix(tmdata[, -1])
    rownames(counts) <- Genes
    colnames(counts) <- make.unique(colnames(counts))
    rownames(counts) <- make.unique(rownames(counts))
    tmdata_list[[name]] <- CreateSeuratObject(counts = counts)
    tmdata_list[[name]]$orig.ident <- rep(name, dim(tmdata_list[[name]])[2])
    print(paste0("finished reading ", name))
  }
  
  return(tmdata_list)
}

inspect <- function(tmdata_list) {
  
  x_features_plot <- list()
  x_count_plot <- list()
  x_mito_plot <- list()
  
  for (name in names(tmdata_list)) {
    
    tmdata_list[[name]][["percent.mt"]] <- PercentageFeatureSet(tmdata_list[[name]], pattern = "^MT-")
    
    mean_nFeature <- mean(tmdata_list[[name]]$nFeature_RNA, na.rm = TRUE)
    median_nFeature <- median(tmdata_list[[name]]$nFeature_RNA, na.rm = TRUE)
    
    mean_nCount <- mean(tmdata_list[[name]]$nCount_RNA, na.rm = TRUE)
    median_nCount <- median(tmdata_list[[name]]$nCount_RNA, na.rm = TRUE)
    
    mean_percent_mt <- mean(tmdata_list[[name]]$percent.mt, na.rm = TRUE)
    median_percent_mt <- median(tmdata_list[[name]]$percent.mt, na.rm = TRUE)
    
    base_theme <- theme(
      text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
    
    # nFeature_RNA plot
    x_features_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_nFeature, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_nFeature, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_nFeature, label = paste("Mean:", round(mean_nFeature, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_nFeature, label = paste("Median:", round(median_nFeature, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
    # nCount_RNA plot
    x_count_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_nCount, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_nCount, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_nCount, label = paste("Mean:", round(mean_nCount, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_nCount, label = paste("Median:", round(median_nCount, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
    # percent.mt plot
    x_mito_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_percent_mt, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_percent_mt, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_percent_mt, label = paste("Mean:", round(mean_percent_mt, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_percent_mt, label = paste("Median:", round(median_percent_mt, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  }
  
  plot_chunks <- function(plot_list) {
    split(plot_list, ceiling(seq_along(plot_list) / 6))
  }
  
  get_grid_dims <- function(n) {
    if (n == 1) return(c(1, 1))
    if (n == 2) return(c(1, 2))
    if (n == 3) return(c(2, 2))
    if (n == 4) return(c(2, 2))
    return(c(2, 3)) # Default for larger sets
  }
  
  pdf(paste0("Inspections", batch, ".pdf"), width = 8, height = 11)
  
  for (plot_list in list(x_features_plot, x_count_plot, x_mito_plot)) {
    for (chunk in plot_chunks(plot_list)) {
      dims <- get_grid_dims(length(chunk))
      grid.arrange(grobs = chunk, ncol = dims[1], nrow = dims[2])
    }
  }
  
  dev.off()
}

normalise <- function(tmdata_list) {
  
  for (name in names(tmdata_list)) {
    CPM <- apply(tmdata_list[[name]]@assays$RNA$counts, 2, function(x) (x / sum(x)) * 1e6)
    CPM <- as(CPM, "dgCMatrix")
    tmdata_list[[name]]@assays$RNA$CPM <- CPM
    expr <- log2((CPM / 10) + 1)
    expr <- as(expr, "CsparseMatrix")
    tmdata_list[[name]]@assays$RNA$data <- expr
    print(paste0("finished normalising ", name))
  }
  
  return(tmdata_list)
}

doublets_filtering <- function(tmdata, sampleid) {
  
  tmdata <- NormalizeData(tmdata)
  tmdata <- FindVariableFeatures(tmdata, selection.method = "vst", nfeatures = 2000)
  tmdata <- ScaleData(tmdata)
  tmdata <- RunPCA(tmdata)
  tmdata <- FindNeighbors(tmdata, dims = 1:50)
  tmdata <- FindClusters(tmdata, resolution = 0.5)
  tmdata <- RunUMAP(tmdata, dims = 1:50)
  
  sweep_tmdata <- paramSweep(tmdata, PCs = 1:50, sct = T)
  sweep_summary <- summarizeSweep(sweep_tmdata)
  bcmvn <- find.pK(sweep_summary)
  pk_tmdata <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  ncells <- ncol(tmdata)
  if (ncells <= 1000) {
    DR <- 0.008
  } else if (ncells <= 5000) {
    DR <- 0.04
  } else if (ncells <= 10000) {
    DR <- 0.08
  } else if (ncells <= 20000) {
    DR <- 0.16
  } else if (ncells <= 30000) {
    DR <- 0.24
  } else {
    DR <- 0.24
  }
  homotypic <- modelHomotypic(tmdata$seurat_clusters)
  nExp <- round(DR * ncol(tmdata))
  nExp_adj <- round(nExp * (1 - homotypic))
  
  tmdata <- doubletFinder(tmdata, PCs = 1:50, pN = 0.25, pK = pk_tmdata,
                          nExp = nExp_adj, reuse.pANN = F, sct = F)
  
  SorD <- grep("DF", colnames(tmdata@meta.data), value = T)
  p1 <- DimPlot(tmdata, reduction = "umap")
  p2 <- DimPlot(tmdata, reduction = "umap", group.by = SorD)
  combined_plot <- p1 + p2
  
  singlet <- colnames(tmdata)[which(tmdata@meta.data[[SorD]] == "Singlet")]
  tmdata <- subset(tmdata, cells = singlet)
  
  return(list(tmdata = tmdata, plot = combined_plot))
}

# doublets_parallel <- function(tmdata_list) {
#
#   cl <- makeCluster(n_clusters)
#   registerDoParallel(cl)
#
#   clusterExport(cl, c("doublets_filtering"))
#
#   tmdata_updated <- foreach(name = names(tmdata_list),
#                             .packages = c("Seurat", "DoubletFinder", "ggplot2"),
#                             .errorhandling = 'pass') %dopar% {
#                               tryCatch({
#                                 doublets_filtering(
#                                   tmdata = tmdata_list[[name]],
#                                   sampleid = name
#                                 )
#                               }, error = function(e) {
#                                 list(tmdata = tmdata_list[[name]], plot = NULL)
#                               })
#                             }
#   stopCluster(cl)
#   names(tmdata_updated) <- names(tmdata_list)
#
#   pdf(paste0("doublets_filtering", batch, ".pdf"), width = 12, height = 8)
#   for (name in names(tmdata_updated)) {
#     if (!is.null(tmdata_updated[[name]]$plot)) {
#       print(tmdata_updated[[name]]$plot + ggtitle(name))
#     }
#   }
#   dev.off()
#
#   tmdata_output <- lapply(tmdata_updated, function(x) x$tmdata)
#   names(tmdata_output) <- names(tmdata_updated)
#
#   return(tmdata_output)
# }

cells_filtering <- function(tmdata_list, rules = qc_rules) {
  
  plot <- list()
  g_filter <- vector()
  hk_filter <- vector()
  for (name in names(tmdata_list)) {
    ##############
    match_idx <- which(sapply(qc_rules$pattern, function(p) grepl(p, name)))
    matched_row <- qc_rules[match_idx[1], ]
    ngenes <- matched_row$nGenes
    hkmean <- matched_row$hk
    ##############
    expr <- as.matrix(tmdata_list[[name]]@assays$RNA$data)
    n_genes <- colSums(expr > 0)
    hk_list <- c("ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
                 "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18")
    
    hk_list <- hk_list[hk_list %in% rownames(expr)]
    hk_expression <- expr[hk_list, , drop = FALSE]
    hk_mean <- colMeans(hk_expression)
    sl_cells_g <- n_genes >= ngenes
    g_filter[[name]] <- sum(sl_cells_g)
    sl_cells_hk <- hk_mean >= hkmean & sl_cells_g
    hk_filter[[name]] <- sum(sl_cells_hk)
    if (sum(sl_cells_hk) != 0) {
      tmdata_list[[name]] <- subset(tmdata_list[[name]], cells = names(sl_cells_hk)[sl_cells_hk])
    } else {
      tmdata_list[[name]] <- NULL
    }
    
    
    plot_data <- data.frame(hk_mean = hk_mean, n_genes = n_genes, valid = sl_cells_hk)
    
    p <- ggplot(plot_data, aes(x = n_genes, y = hk_mean, color = valid)) +
      geom_point() +
      scale_x_continuous(trans = "log10", labels = scales::comma) +
      scale_y_continuous(trans = "log10", labels = scales::comma) +
      scale_color_manual(values = c("lightgrey", "black")) +
      labs(x = "Number of Genes", y = "HK Mean", color = "Valid Cells") +
      theme_minimal() +
      theme(legend.position = "none", plot.title = element_text(size = 8)) +
      annotate("text", x = Inf, y = Inf, label = paste0("NCells passed: ", sum(sl_cells_hk)),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
      ggtitle(name) +
      geom_vline(xintercept = ngenes, linetype = "dashed", color = "red") +
      geom_hline(yintercept = hkmean, linetype = "dashed", color = "red") +
      annotate("text", x = ngenes, y = min(plot_data$hk_mean), label = paste0("Number of genes > ", ngenes),
               hjust = -0.1, vjust = 0, size = 3, color = "red") +
      annotate("text", x = min(plot_data$n_genes), y = hkmean, label = paste0("HK Mean > ", hkmean),
               hjust = 0, vjust = -0.5, size = 3, color = "red")
    
    plot[[name]] <- p
    print(paste0("finished cell filtering for  ", name))
  }
  
  split_plots <- function(plot_list) {
    split(plot_list, ceiling(seq_along(plot_list) / 6))
  }
  get_layout_dims <- function(n) {
    if (n == 1) return(c(1, 1))
    if (n == 2) return(c(1, 2))
    if (n == 3 || n == 4) return(c(2, 2))
    return(c(2, 3))  # Default layout for 5-6 plots
  }
  
  pdf(paste0("cells_filtering", batch, ".pdf"), width = 8, height = 11)
  
  for (chunk in split_plots(plot)) {
    dims <- get_layout_dims(length(chunk))
    grid.arrange(grobs = chunk, ncol = dims[1], nrow = dims[2])  # ncol and nrow reversed
  }
  dev.off()
  
  return(list(tmdata_list, g_filter, hk_filter))
}

write_count_matrix <- function(filtered) {
  
  filtered_names <- names(filtered)
  
  cl <- makeCluster(3)
  registerDoParallel(cl)
  
  foreach(name = names(filtered), .packages = "data.table") %dopar% {
    data.table::fwrite(
      as.data.frame(t(as.matrix(filtered[[name]]@assays$RNA$counts))),
      paste0("count_matrix/count_matrix_", name, ".csv"),
      row.names = TRUE
    )
  }
  
  stopCluster(cl)
}

##################################################################

raw <- initialise(names_tmdata)
x_filter <- sapply(raw, ncol)
filtered <- raw
rm(raw)
for (name in names(filtered)) {
  if (!any(grepl("^counts$", Layers(filtered[[name]]@assays$RNA)))) {
    filtered[[name]]@assays$RNA$counts <- filtered[[name]]@assays$RNA$`counts.Gene Expression`
    filtered[[name]]@assays$RNA$`counts.Gene Expression` = NULL
    filtered[[name]]@assays$RNA$`counts.Peaks` = NULL
  }
}
inspect(filtered)
#filtered <- doublets_parallel(filtered)
#print("finished finding doublets")
#db_filter <- sapply(filtered, ncol)
for (name in names(filtered)) {
  filtered[[name]][["percent.mt"]] <- PercentageFeatureSet(filtered[[name]], pattern = "^MT-")
  match_idx <- which(sapply(qc_rules$pattern, function(p) grepl(p, name)))
  matched_row <- qc_rules[match_idx[1], ]
  max_mt <- matched_row$mito
  if (sum(filtered[[name]]$percent.mt < max_mt) != 0) {
    filtered[[name]] <- subset(filtered[[name]], percent.mt < max_mt)
  } else {
    filtered[[name]] <- NULL
  }
}

for (name in names(filtered)) {
  if (ncol(filtered[[name]]) <= 1) {
    filtered[[name]] <- NULL
  }
}

mt_filter <- sapply(filtered, ncol)
filtered <- normalise(filtered)
cells_ft_outs <- cells_filtering(filtered, rules = qc_rules)
filtered <- cells_ft_outs[[1]]
g_filter <- cells_ft_outs[[2]]
hk_filter <- cells_ft_outs[[3]]

for (name in names(filtered)) {
  if (ncol(filtered[[name]]) <= 1) {
    filtered[[name]] <- NULL
  }
}

sm_table <- list(
  x_filter,
  mt_filter,
  g_filter,
  hk_filter
) |>
  map(~ tibble::enframe(.x, name = "sample", value = "value")) |>
  reduce(full_join, by = "sample") |>
  mutate(across(-sample, ~ tidyr::replace_na(.x, 0)))

colnames(sm_table) <- c("sample", "raw", paste0("mito_DNA\npercentage < ", max_mt),
                        paste0("number of\ngenes > ", min_ngenes), paste0("housekeeping\nexpression > ", min_hk_expr))
write.csv(sm_table, paste0("filtering_summary", batch, ".csv"))
#write_count_matrix(filtered)

##############################Markers list#####################################

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


housekeeping <- c("ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
                  "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18")

markers <- c(fibroblast, macrophage, mast, epithelial, t.cell, nk.cell,  
             b.cell, plasma, dendritic, endothelial, housekeeping)

markers_list <- list(
  b.cell = b.cell,
  dendritic = dendritic,
  endothelial = endothelial,
  epithelial = epithelial,
  fibroblast = fibroblast,
  macrophage = macrophage,
  mast = mast,
  nk.cell = nk.cell, 
  plasma = plasma, 
  t.cell = t.cell,
  housekeeping = housekeeping
)

#############################################################

manual_celltyping <- function(tmdata_list) {
  
  for (name in names(tmdata_list)) {
    if (ncol(tmdata_list[[name]]) > 1){
      expr <- as.matrix(tmdata_list[[name]]@assays$RNA$data)
      tp_markers_list <- markers_list
      for (marker in names(tp_markers_list)) {
        tp_markers_list[[marker]] <- tp_markers_list[[marker]][tp_markers_list[[marker]] %in% rownames(expr)]
      }
      tp_markers_list$housekeeping <- NULL
      ct_mtx <- matrix(0, ncol = length(tp_markers_list), nrow = ncol(expr),
                       dimnames = list(colnames(expr), names(tp_markers_list)))
      for (celltype in names(tp_markers_list)) {
        marker_genes <- tp_markers_list[[celltype]]
        marker_expr <- expr[marker_genes, ]
        marker_mean <- ifelse(rep(!is.null(nrow(marker_expr)), ncol(expr)),
                              colMeans(marker_expr), marker_expr)
        ct_mtx[, celltype] <- marker_mean
      }
      tmdata_list[[name]]@meta.data <- cbind(tmdata_list[[name]]@meta.data, ct_mtx)
      celltypes <- apply(ct_mtx, 1, function(row) {
        valid_row <- row[!is.na(row)]
        max_mean <- max(valid_row, na.rm = TRUE)
        if (is.infinite(max_mean)) {
          return("unclassified")
        }
        max_name <- names(valid_row)[valid_row == max_mean]
        if (length(max_name) > 1 || max_mean < 1) {
          return("unclassified")
        } else {
          return(max_name)
        }
      })
      tmdata_list[[name]]@meta.data$celltype_manual <- celltypes
    }
  }
  
  return(tmdata_list)
}

tmdata_annotated <- manual_celltyping(filtered)

saveRDS(tmdata_annotated, paste0("/rds/general/project/spatialtranscriptomics/ephemeral/EAC_Ref_list", batch, ".rds"))

################################################################
# 
# celltyping <- function(tmdata, sampleid) {
# 
#   set.seed(0)
#   counts <- as.matrix(tmdata@assays$RNA$counts)
#   neg <- rep(0, length(colnames(tmdata)))
#   names(neg) <- colnames(tmdata)
# 
#   pm_sub <- profile_matrix[is.element(rownames(profile_matrix), rownames(tmdata)), ]
# 
#   annotation <- do.call(rbind, (lapply(names(cellGroups), function(x)
#   {data.frame(celltype = rep(x, length(cellGroups[[x]])), row.names = cellGroups[[x]])})))
#   annotation <- annotation[colnames(pm_sub), , drop = FALSE]
# 
#   celltypes <- InSituType::insitutypeML(x = t(counts),
#                                         neg = neg,
#                                         reference_profiles = as.matrix(pm_sub))
# 
#   cols <- InSituType::colorCellTypes(freqs = table(celltypes$clust), palette = "brewers")
# 
#   flightpath <- InSituType::flightpath_layout(logliks = celltypes$logliks, profiles = celltypes$profiles)
#   png(paste0("Flightpath_", sampleid, ".png"), width = 20, height = 14, units = "in", res = 600)
#   par(mar = c(0,0,0,0))
#   plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[celltypes$clust])
#   text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
#   print("Flightpath created")
#   dev.off()
# 
#   tmdata@meta.data$celltype <- celltypes$clust[rownames(tmdata@meta.data)]
#   tmdata@meta.data$celltype_meta <- annotation[match(celltypes$clust, rownames(annotation)), ]
#   tmdata@meta.data$celltype_meta[is.na(tmdata@meta.data$celltype_meta)] <-
#     tmdata@meta.data$celltype[is.na(tmdata@meta.data$celltype_meta)]
# 
#   return(tmdata)
# }
# 
# celltyping_parallel <- function(tmdata_list) {
# 
#   cl <- makeCluster(n_clusters)
#   registerDoParallel(cl)
# 
#   clusterExport(cl, c("celltyping", "profile_matrix", "cellGroups"))
# 
#   tmdata_updated <- foreach(name = names(tmdata_list),
#                             .packages = c("Seurat","dplyr","patchwork","InSituType","pheatmap","ggplot2"),
#                             .errorhandling = 'pass') %dopar% {
#                               tryCatch({
#                                 celltyping(
#                                   tmdata = tmdata_list[[name]],
#                                   sampleid = name
#                                 )
#                               }, error = function(e) {
#                                 tmdata_list[[name]]
#                               })
#                             }
#   stopCluster(cl)
#   names(tmdata_updated) <- names(tmdata_list)
# 
#   return(tmdata_updated)
# }
# 
# #load("/rds/general/project/tumourheterogeneity1/live/CosMx/Esophagus_HCA.RData")
# #tmdata_annotated <- celltyping_parallel(filtered)
# print("finished celltyping")
# 
###################################################################

plot_heatmap <- function(tmdata, sampleid, identity, reorder) {
  
  expr_data <- FetchData(tmdata, vars = markers, layer = "data")
  missing_markers <- setdiff(markers, colnames(expr_data))
  for (gene in missing_markers) {
    expr_data[[gene]] <- 0
  }
  expr_data <- expr_data[, markers]
  for (name in names(markers_list)) {
    markers_list[[name]] <- markers_list[[name]][markers_list[[name]] %in% colnames(expr_data)]
  }
  expr_data$celltype <- tmdata[[identity]]
  expr_data$celltype <- sapply(expr_data$celltype, function(x) gsub(" ", "\n", x))
  
  expr_data <- expr_data[order(expr_data$celltype), ]
  hk_avg <- rowMeans(expr_data[, colnames(expr_data) %in% markers_list$housekeeping, drop = FALSE])
  hk_avg <- matrix(hk_avg, nrow = 1, dimnames = list("avg_hk", names(hk_avg)))
  nCounts <- tmdata$nFeature_RNA[rownames(expr_data)]
  nCounts <- matrix(nCounts, nrow = 1, dimnames = list("nGenes", names(nCounts)))
  
  # Prepare heatmap data
  heatplot <- t(as.matrix(expr_data[, 1:ncol(expr_data)-1]))
  
  # Custom colors
  custom_colors <- colorRamp2(c(0, round(0.8 * max(heatplot), 1), ceiling(max(heatplot))),
                              c("#D0D0D0", "red4", "red4"))
  
  heatmap_grobs <- list()
  len <- length(unique(expr_data$celltype))
  
  if (!reorder) {
    desired_order <- unique(expr_data$celltype)
  } else {
    desired_order <- ct_reorder
  }
  present_celltypes <- intersect(desired_order, unique(expr_data$celltype))
  
  match_idx <- which(sapply(qc_rules$pattern, function(p) grepl(p, tmdata$orig.ident[1])))
  matched_row <- qc_rules[match_idx[1], ]
  min_ngenes <- matched_row$nGenes
  min_hk_expr <- matched_row$hk
  max_mt <- matched_row$mito
  
  # Loop through each cell marker gene to create and draw heatmaps
  for (i in seq_along(markers_list)) {
    
    marker <- markers_list[[i]]
    marker_names <- rownames(heatplot[marker, , drop = FALSE])
    temp <- list()
    for (j in seq_along(present_celltypes)) {
      
      celltype <- present_celltypes[j]
      cells <- rownames(expr_data)[which(expr_data$celltype == celltype)]
      
      ht <- Heatmap(
        heatplot[marker, cells, drop = FALSE],
        col = custom_colors,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 40),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        use_raster = FALSE
      )
      
      # Capture the heatmap as a grob
      ht_grob <- grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(2, 1, 2, 1), "mm")))
      
      # Add the heatmap grob to the list
      temp[[j]] <- ht_grob
      
    }
    
    gene_label_grobs <- lapply(marker_names, function(name) {
      textGrob(
        label = name,
        x = unit(1, "npc"), just = "right", gp = gpar(fontsize = 40)
      )
    })
    
    # Arrange them into a single column
    gene_label_col <- arrangeGrob(
      grobs = gene_label_grobs,
      ncol = 1
    )
    
    temp <- c(list(gene_label_col), temp)
    
    if (j == length(present_celltypes)) {
      text_grob <- textGrob(
        names(markers_list)[i],
        gp = gpar(fontsize = 50, fontface = "bold")
      )
      
      rect_grob <- rectGrob(
        gp = gpar(fill = "grey", col = NA)
      )
      
      merged_grob <- gTree(children = gList(rect_grob, text_grob))
      
      # Arrange the text and heatmap grobs side by side
      temp[[j + 2]] <- merged_grob
      combined_grob <- do.call(
        arrangeGrob, c(temp, list(ncol = len + 2, widths = c((len+1)/23, rep(1, len), (len+1)/12))))
      
      heatmap_grobs[[length(heatmap_grobs) + 1]] <- combined_grob
    }
    
  }
  
  stats_grobs <- list()
  
  # Loop through each cell marker gene to create and draw heatmaps
  for (i in seq_along(list(nCounts, hk_avg))) {
    
    data <- list(nCounts, hk_avg)[[i]]
    data_names <- rownames(data)
    temp <- list()
    for (j in seq_along(present_celltypes)) {
      
      celltype <- present_celltypes[j]
      cells <- rownames(expr_data)[which(expr_data$celltype == celltype)]
      
      if (i == 1){
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name  = "Number of\ngenes detected",
          col   = colorRamp2(c(min_ngenes, round(0.7 * max(nCounts), -3), ceiling(max(nCounts) / 1000) * 1000),
                             c("#D0D0D0", "blue3", "blue3")),
          cluster_rows    = FALSE,
          cluster_columns = FALSE,
          show_row_names  = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      } else {
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name  = "Average\nhousekeeping\nexpression",
          col   = colorRamp2(c(0, min_hk_expr, rounded_value <- ceiling(max(hk_avg) * 10) / 10),
                             c("#D0D0D0", "#A0B8E6", "blue3")),
          cluster_rows    = FALSE,
          cluster_columns = FALSE,
          show_row_names  = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      }
      
      # Capture the heatmap as a grob
      ht_grob <- grid.grabExpr(draw(ht_add, newpage = FALSE, padding = unit(c(6, 1.5, 6, 1.5), "mm")))
      
      # Add the heatmap grob to the list
      temp[[j]] <- ht_grob
      
    }
    
    gene_label_grobs <- lapply(data_names, function(name) {
      textGrob(
        label = name,
        just = "center", gp = gpar(fontsize = 40)
      )
    })
    
    # Arrange them into a single column
    gene_label_col <- arrangeGrob(
      grobs = gene_label_grobs,
      ncol = 1
    )
    
    temp <- c(list(gene_label_col), temp)
    
    if (j == length(present_celltypes)) {
      text_grob <- textGrob(
        c("Number\nof genes", "Average\nhousekeeping\nexpression")[i],
        gp = gpar(fontsize = 40, fontface = "bold")
      )
      
      rect_grob <- rectGrob(
        gp = gpar(fill = "white", col = NA)
      )
      
      merged_grob <- gTree(children = gList(rect_grob, text_grob))
      
      # Arrange the text and heatmap grobs side by side
      
      temp[[j + 2]] <- merged_grob
      combined_grob <- do.call(
        arrangeGrob, c(temp, list(ncol = len + 2, widths = c((len+1)/23, rep(1, len), (len+1)/12))))
      stats_grobs[[length(stats_grobs) + 1]] <- combined_grob
    }
    
  }
  
  # Create cell type labels
  celltype_labels <- c("", as.character(present_celltypes), "", "")  # Add an empty label to match the grid size
  text_grobs <- lapply(celltype_labels, function(celltype) {
    textGrob(celltype, gp = gpar(fontsize = 50, fontface = "bold"), just = "center", rot = 0)
  })
  
  # Combine the text grobs into a single row for alignment with heatmaps
  text_grob_row <- do.call(
    arrangeGrob, c(text_grobs, list(ncol = len + 2, widths = c((len+1)/23, rep(1, len), (len+1)/12))))
  
  celltype_num <- c("", as.numeric(table(expr_data$celltype)[present_celltypes]), "Markers", "")
  num_grobs <- lapply(celltype_num, function(celltype_num) {
    textGrob(celltype_num, gp = gpar(fontsize = 40, fontface = "bold"), just = "center")
  })
  
  num_grob_row <- do.call(
    arrangeGrob, c(num_grobs, list(ncol = len + 2, widths = c((len+1)/23, rep(1, len), (len+1)/12))))
  
  title_grob <- textGrob(
    paste0("Mitochondria DNA < ",
           max_mt, "     &&     Minimum number of genes > ",
           min_ngenes, "     &&     Minimum HK expression > ", min_hk_expr),
    gp = gpar(fontsize = 40),
    just = "center"
  )
  
  legend_obj1 <- Legend(title = "\nExpression\nlevels (E)\n",
                        at = pretty(c(0, ceiling(max(heatplot))), n = 5),
                        labels_gp = gpar(fontsize = 50),
                        title_gp = gpar(fontsize = 50, fontface = "bold"),
                        grid_height = unit(220, "mm"),
                        legend_width = unit(40, "mm"),
                        legend_height = unit(220, "mm"),
                        grid_width = unit(40, "mm"),   # make the bar thicker
                        title_position = "topcenter",
                        col_fun = colorRamp2(c(0, round(0.8 * max(heatplot), 1), ceiling(max(heatplot))),
                                             c("#D0D0D0", "red4", "red4")))
  
  legend_obj2 <- Legend(title = "\nNumber\nof genes\n",
                        at = round(seq(min_ngenes, ceiling(max(nCounts) / 1000) * 1000, length.out = 4), -3),
                        labels_gp = gpar(fontsize = 50),
                        title_gp = gpar(fontsize = 50, fontface = "bold"),
                        grid_height = unit(220, "mm"),
                        legend_width = unit(40, "mm"),
                        legend_height = unit(220, "mm"),
                        grid_width = unit(40, "mm"),   # make the bar thicker
                        title_position = "topcenter",
                        col_fun = colorRamp2(c(min_ngenes, round(0.7 * max(nCounts), -3), ceiling(max(nCounts) / 1000) * 1000),
                                             c("#D0D0D0", "blue3", "blue3")))
  
  legend_obj3 <- Legend(title = "\nAverage\nhousekeeping\nexpression\n",
                        at = seq(floor(min_hk_expr), ceiling(max(hk_avg) * 10) / 10, length.out = 4) %>% round(0),
                        labels_gp = gpar(fontsize = 50),
                        title_gp = gpar(fontsize = 50, fontface = "bold"),
                        grid_height = unit(220, "mm"),
                        legend_width = unit(40, "mm"),
                        legend_height = unit(220, "mm"),
                        grid_width = unit(40, "mm"),   # make the bar thicker
                        title_position = "topcenter",
                        col_fun = colorRamp2(c(0, min_hk_expr, rounded_value <- ceiling(max(hk_avg) * 10) / 10),
                                             c("#D0D0D0", "#A0B8E6", "blue3")))
  
  main_content <- arrangeGrob(
    grobs = c(list(title_grob),
              list(text_grob_row),
              list(num_grob_row),
              heatmap_grobs,    # your heatmap rows
              stats_grobs),     # your additional stats rows
    ncol = 1,
    heights = c(2.5, 2, 1.5, rep(4, length(heatmap_grobs)), 4, 4)
  )
  
  legend_grob <- textGrob(
    label = sampleid,
    gp = gpar(fontsize = 35, fontface = "bold")
  )
  legend_grob0 <- textGrob(
    label = paste0("Total cell count: ", sum(expr_data$celltype %in% present_celltypes)),
    gp = gpar(fontsize = 35)
  )
  legend_grob1 <- grid.grabExpr(draw(legend_obj1))
  legend_grob2 <- grid.grabExpr(draw(legend_obj2))
  legend_grob3 <- grid.grabExpr(draw(legend_obj3))
  # Combine them vertically (adjust ncol if you want them side-by-side instead)
  legend_column <- arrangeGrob(legend_grob, legend_grob0, legend_grob1, legend_grob2, legend_grob3,
                               ncol = 1, heights = c(0.3, 0.2, 2, 2, 2))
  
  final_layout <- grid.arrange(
    main_content, legend_column,
    ncol = 2,
    widths = c(9.3, 1)  # You can tweak these numbers as needed.
  )
  
  png(paste0(identity, "_", sampleid, ".png"), width = 80, height = 50, units = "in", res = 400)
  grid.draw(final_layout)
  dev.off()
  
  return("pass")
}
# 
# plot_parallel <- function(tmdata_list) {
#   
#   cl <- makeCluster(n_clusters)
#   registerDoParallel(cl)
#   
#   clusterExport(cl, c("plot_heatmap", "markers", "markers_list", "color_scale",
#                       "celltyping_method", "max_mt", "min_ngenes", "min_hk_expr", "ct_reorder"))
#   
#   error_message <- foreach(name = names(tmdata_list),
#                            .packages = c("Seurat", "dplyr", "patchwork", "InSituType",
#                                          "pheatmap", "ggplot2", "ComplexHeatmap",
#                                          "circlize", "gridExtra", "grid"),
#                            .errorhandling = 'pass') %dopar% {
#                              tryCatch({
#                                plot_heatmap(
#                                  tmdata = tmdata_list[[name]],
#                                  sampleid = name,
#                                  identity = celltyping_method
#                                )
#                                NULL
#                              }, error = function(e) {
#                                e$message
#                              })
#                            }
#   stopCluster(cl)
#   names(error_message) <- names(tmdata_list)
#   failed_sampleid <- !sapply(error_message, is.null)
#   message_out <- ""
#   for (failed in names(error_message)[failed_sampleid]) {
#     message_out <- paste0(message_out, "pipeline failed for sample: ",
#                           failed, "\nwith error message: ", error_message[[failed]], "\n\n")
#   }
#   
#   return(message_out)
# }
#
# final <- plot_parallel(tmdata_annotated)
# cat(final)
# 
# ##########################################################################
# 
# plot_coexpression <- function(tmdata_list, identity, nfilter = 5) {
# 
#   for (name in names(tmdata_list)) {
# 
#     ct_mtx <- tmdata_list[[name]]@meta.data[, colnames(tmdata_list[[name]]@meta.data) %in% names(markers_list)]
#     ct_mtx$celltype <- tmdata_list[[name]]@meta.data[[identity]]
#     heatmap_data <- ct_mtx %>%
#       filter(celltype != "unclassified") %>%
#       group_by(celltype) %>%
#       filter(n() > nfilter) %>%
#       summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
#       as.data.frame()
#     rownames(heatmap_data) <- heatmap_data$celltype
#     heatmap_data <- heatmap_data[, -1]
# 
#     color_scale <- colorRamp2(c(0, 2, 5), c("#F4F4F4", "red2", "red3"))
#     if (dim(heatmap_data)[1] != 0) {
#       png(paste0("coexpression_", name, ".png"), width = 6, height = 6, units = "in", res = 500)
#       hm <- Heatmap(
#         t(as.matrix(heatmap_data)),
#         name = "Expression",
#         col = color_scale,
#         show_column_names = TRUE,
#         show_row_names = TRUE,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         show_row_dend = FALSE,
#         show_column_dend = FALSE,
#         column_names_side = "top",
#         column_names_rot = 45
#       )
#       draw(hm)
#       dev.off()
#     }
#   }
# }
# 
# plot_coexpression(tmdata_annotated, identity = celltyping_method, nfilter = least_ncells)

####################################END########################################

# tmdata_annotated <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/all_samples_list.rds")
# 
# marker_columns <- c("epithelial", "endothelial", "fibroblast", "b.cell", "t.cell", "macrophage",
#                     "dendritic", "mast", "neutrophil", "erythrocyte", "keratinocyte", "lymph")
# 
# for (name in names(tmdata_annotated)) {
#   
#   print(name)
#   celltypist <- read.table(paste0("celltypist/", name, "/predicted_labels.csv"), sep = ",", row.names = 1, header = TRUE)
#   if (!"majority_voting" %in% colnames(celltypist)) {
#     message(name, ": no majority_voting column")
#     tmdata_annotated[[name]] <- NULL
#     next
#   }
#   tmdata_annotated[[name]]$celltypist <- celltypist$majority_voting
#   
#   high_marker_count <- rowSums(tmdata_annotated[[name]]@meta.data[, marker_columns] > 1)
#   
#   tmdata_annotated[[name]]@meta.data$coexpression <- ifelse(high_marker_count > 1, "doublet", "singlet")
#   
#   top_two_diff <- apply(tmdata_annotated[[name]]@meta.data[, marker_columns], 1, function(x) {
#     sorted <- sort(x, decreasing = TRUE)
#     sorted[1] - sorted[2]
#   })
#   
#   tmdata_annotated[[name]]@meta.data$coexpression_loose <- ifelse((high_marker_count > 1) & (top_two_diff < 1), "doublet", "singlet")
#   
#   annotation_to_marker <- c(
#     "B_GC" = "b.cell", "B_mature" = "b.cell", "B_progenitor" = "b.cell", "Plasma" = "b.cell",
#     "Dendritic_classical" = "dendritic", "Dendritic_plasmacytoid" = "dendritic",
#     "Endothelial" = "endothelial",
#     "Ductal" = "epithelial", "Squamous" = "epithelial", "Schwann" = "epithelial", "Ciliated" = "epithelial",
#     "Erythroid" = "erythrocyte",
#     "Mural" = "fibroblast", "Fibroblast" = "fibroblast",
#     "Macrophage" = "macrophage", "Monocyte" = "macrophage",
#     "Mast&Basophil" = "mast",
#     "Neutrophil" = "neutrophil",
#     "T&NK" = "t.cell",
#     "Neuron_excitatory" = "_others_", "Neuron_inhibitory" = "_others_", "Neuron_bipolar" = "_others_",
#     "Astrocyte" = "_others_", "Muller" = "_others_", "Oligodendrocyte_mature" = "_others_",
#     "Hematopoietic" = "_others_", "Hepatocyte" = "_others_", "Melanocyte" = "_others_", "Rod" = "_others_", "Spermatocyte" = "_others_"
#   )
#   
#   tmdata_annotated[[name]]@meta.data$marker_expression <- mapply(function(celltype, row_index) {
#     if (celltype %in% names(annotation_to_marker)) {
#       marker_col <- annotation_to_marker[[celltype]]
#       if (marker_col != "_others_") {
#         if (!is.na(tmdata_annotated[[name]]@meta.data[row_index, marker_col])) {
#           if (tmdata_annotated[[name]]@meta.data[row_index, marker_col] > 1.5) {
#             return("good")
#           } else {
#             return("poor")
#           }
#         } else {
#           return(NA)
#         }
#       } else {
#         return(NA)
#       }
#     } else {
#       return(NA)  # Celltype not in annotation mapping
#     }
#   }, tmdata_annotated[[name]]@meta.data$celltypist, seq_len(nrow(tmdata_annotated[[name]]@meta.data)))
# 
#   tmdata_annotated[[name]] <- tryCatch({
#     tmp <- subset(tmdata_annotated[[name]],
#                   subset = (marker_expression == "good" | is.na(marker_expression)) & coexpression == "singlet")
#     if (ncol(tmp) <= 1) {
#       message(name, ": ≤1 cell after subset — skipping")
#       NULL
#     } else {
#       tmp
#     }
#   }, error = function(e) {
#     message(name, ": Error during subset — skipping")
#     NULL
#   })
# 
#   tmdata_annotated[[name]]$celltype_group <- annotation_to_marker[as.character(tmdata_annotated[[name]]$celltypist)]
# 
#   unmatched <- unique(tmdata_annotated[[name]]$celltypist[is.na(tmdata_annotated[[name]]$celltype_group)])
#   if (length(unmatched) > 0) {
#     warning("The following cell types in 'celltypist' were not found in the mapping vector and resulted in NAs:")
#     print(unmatched)
#   }
# }
# 
# saveRDS(tmdata_annotated, "/rds/general/project/spatialtranscriptomics/ephemeral/all_samples_list_filtered.rds")
# print(sapply(tmdata_annotated, ncol))
# 
# celltyping_method = "celltypist"
# final <- plot_parallel(tmdata_annotated)
# cat(final)
# # 
# # # ###############################################################################
# # #
# # # tmdata_annotated <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/all_samples_list_filtered.rds")
# # #
# all_genes <- Reduce(union, lapply(tmdata_annotated, function(obj) {
#   rownames(GetAssayData(obj, layer = "counts"))
# }))
# 
# pad_matrix <- function(mat, all_genes) {
#   missing_genes <- setdiff(all_genes, rownames(mat))
#   if (length(missing_genes) > 0) {
#     zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
#                                dimnames = list(missing_genes, colnames(mat)))
#     mat <- rbind(mat, zero_mat)
#   }
#   # Ensure the same row order
#   mat <- mat[all_genes, , drop = FALSE]
#   return(mat)
# }
# 
# counts_list <- lapply(names(tmdata_annotated), function(id) {
#   mat <- GetAssayData(tmdata_annotated[[id]], layer = "counts")
#   colnames(mat) <- paste(id, colnames(mat), sep = "_")
#   pad_matrix(mat, all_genes)
# })
# 
# cpm_list <- lapply(names(tmdata_annotated), function(id) {
#   mat <- GetAssayData(tmdata_annotated[[id]], layer = "CPM")
#   colnames(mat) <- paste(id, colnames(mat), sep = "_")
#   pad_matrix(mat, all_genes)
# })
# 
# lognorm_list <- lapply(names(tmdata_annotated), function(id) {
#   mat <- GetAssayData(tmdata_annotated[[id]], layer = "data")
#   colnames(mat) <- paste(id, colnames(mat), sep = "_")
#   pad_matrix(mat, all_genes)
# })
# 
# meta_list <- lapply(names(tmdata_annotated), function(id) {
#   meta <- tmdata_annotated[[id]]@meta.data[, c(names(tmdata_annotated[[id]]@meta.data)[1:4], "celltype_manual")]#, "celltypist")]
#   rownames(meta) <- paste(id, rownames(meta), sep = "_")
#   return(meta)
# })
# 
# # Combine
# combined_counts <- do.call(cbind, counts_list)
# combined_cpm <- do.call(cbind, cpm_list)
# combined_lognorm <- do.call(cbind, lognorm_list)
# combined_meta <- do.call(rbind, meta_list)
# 
# # Create object with counts
# merged_obj <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
# 
# merged_obj@assays$RNA$CPM <- combined_cpm
# merged_obj@assays$RNA$data <- combined_lognorm
# rm(tmdata_annotated)
# 
# merged_obj@meta.data$study <- sapply(strsplit(merged_obj@meta.data$orig.ident, "_"), function(x) paste(x[1:2], collapse = "_"))
# 
# #merged_obj <- NormalizeData(merged_obj)
# merged_obj <- FindVariableFeatures(merged_obj)
# merged_obj <- ScaleData(merged_obj)
# merged_obj <- RunPCA(merged_obj)
# merged_obj <- FindNeighbors(merged_obj, dims = 1:20)
# 
# merged_obj <- FindClusters(merged_obj, resolution = 0.8, algorithm = 1)  # Leiden
# merged_obj$leiden_clusters <- Idents(merged_obj)
# 
# merged_obj <- RunUMAP(merged_obj, dims = 1:20)
# 
# saveRDS(merged_obj, "/rds/general/project/spatialtranscriptomics/ephemeral/EAC_Ref.rds")
# 
# library(readxl)
# data <- read_excel("/rds/general/project/tumourheterogeneity1/ephemeral/EAC_Ref_all/00_merged/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
# data <- data %>%
#   mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))
# cell_names <- colnames(merged_obj)
# merged_obj@meta.data <- merged_obj@meta.data %>%
#   left_join(data, by = "orig.ident")
# rownames(merged_obj@meta.data) <- cell_names
# 
# #p1 <- DimPlot(merged_obj, group.by = "leiden_clusters", label = TRUE) + ggtitle("Louvain Clustering")
# p4 <- DimPlot(merged_obj, group.by = "celltype_manual", label = FALSE) + ggtitle("Celltype")
# p2 <- DimPlot(merged_obj, group.by = "T_Status", label = FALSE) + ggtitle("Sample Origin")
# p1 <- DimPlot(merged_obj, group.by = "study", label = FALSE) + ggtitle("Study")
# p3 <- DimPlot(merged_obj, group.by = "Technology", label = FALSE) + ggtitle("Technology")
# #p5 <- DimPlot(merged_obj, group.by = "Treatment", label = FALSE) + ggtitle("Treatment Status")
# 
# #combined_plot <- p1 + p2 + p3
# combined_plot <- (p1 | p2) / (p3 | p4)
# 
# ggsave("clustering_comparison_filtered_tsamples.png", plot = combined_plot, width = 12, height = 6, dpi = 300)