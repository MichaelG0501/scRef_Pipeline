###############Loding required packages########################

library("Seurat")
library("dplyr")
library("Seurat")
library("patchwork")
library("ggplot2")
library("foreach")
library("doParallel")
library("circlize")
library("gridExtra")
library("grid")
library("tidyr")
library("tibble")
library("purrr")
library("readxl")
library('parallel')

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

######################Setting parameters########################

data <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Summary_EAC_Ref.xlsx", sheet = 2, skip = 1)
data$study <- paste(data$Author, data$Year, data$`Sample Name`, sep = "_")
tsample <- data$study[data$T_Status == "Tumour"]
writeLines(tsample, "names_tmdata_EAC_Ref_t.txt")
#nsample <- data$study[data$T_Status != "Tumour"]
#writeLines(nsample, "names_tmdata_EAC_Ref_n.txt")

batch <- "_EAC_Ref_t"
names_tmdata <- readLines(paste0("names_tmdata", batch, ".txt"))
n_clusters = 8

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

colnames(sm_table) <- c("sample", "raw", "mito_DNA\npercentage < ", "number of\ngenes", "housekeeping\nexpression > ")
write.csv(sm_table, paste0("filtering_summary", batch, ".csv"))

saveRDS(filtered, paste0("EAC_Ref_list", batch, ".rds"))

out_dir <- "by_samples"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
mclapply(names(filtered), function(nm) {
  sample_dir <- file.path(out_dir, nm)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
  saveRDS(filtered[[nm]], file.path(sample_dir, paste0(nm, ".rds")), compress = FALSE)
  NULL
}, mc.cores = n_clusters, mc.preschedule = FALSE)
