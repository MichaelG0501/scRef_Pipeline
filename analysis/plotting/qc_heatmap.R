####################
# Analysis registry:
#   Status: active
#   Script: analysis/plotting/qc_heatmap.R
#   Methodology: analysis/methodology/plotting/publication_plotting_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/heatmap.R
# Reorganized as part of analysis/ restructuring
####################
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


max_mt = 15
min_ngenes = 200
min_hk_expr = 0
n_clusters = 8
celltyping_method = "first_type"
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
    "*"
  ),
  mito   = c(15),
  nGenes = c(200),
  hk     = c(0),
  stringsAsFactors = FALSE
)


setwd("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")


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
  custom_colors <- colorRamp2(c(0, round(0.6 * max(heatplot), 1), ceiling(max(heatplot))),
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
                        col_fun = colorRamp2(c(0, round(0.6 * max(heatplot), 1), ceiling(max(heatplot))),
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


samples <- c("Carroll_2023_EAC-JCNP_PreTx_tumour_frozen")

for (sample in samples) {
  tmdata <- readRDS(paste0("by_samples/", sample, "/", sample, "_anno.rds"))
  toplot <- subset(tmdata, subset = coexpression == "singlet" & marker_expression == "good")
  plot_heatmap(toplot, paste0(sample, "strict"), "celltype_update", reorder = TRUE)
  toplot <- subset(tmdata, subset = coexpression_loose == "singlet" & marker_expression == "good")
  plot_heatmap(toplot, paste0(sample, "loose"), "celltype_update", reorder = TRUE)
}

