####################
# Analysis registry:
#   Status: active
#   Script: analysis/publication/poster_section1_atlas_metaprograms.R
#   Inputs:
#     ref_outs/all_meta.rds
#     ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#     ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#     ref_outs/cluster_enrich.rds
#     ref_outs/Auto_topmp_v2_noreg_states_B.rds
#     ref_outs/Auto_topmp_v2_noreg_mp_adj.rds
#     ref_outs/meta_full_epi.rds
#   Outputs:
#     ref_outs/publication/section1/figures/*.pdf/png
#   Run command: Rscript analysis/publication/poster_section1_atlas_metaprograms.R
#   Conda environment: dmtcp
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(scattermore)
  library(grid)
  library(Seurat)
  library(patchwork)
})
source("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/publication/publication_helpers.R")

section <- "section1"
out_dir <- pub_section_dir(section)
ucell_file <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "UCell_nMP19_filtered.rds")
mp_group_df <- tibble(mp = PUB_MP_ORDER, group = pub_mp_state(PUB_MP_ORDER)) |>
  mutate(group = factor(group, levels = PUB_MP_STATE_ORDER)) |>
  arrange(group, match(mp, PUB_MP_ORDER)) |>
  mutate(group_id = as.integer(group),
         x = row_number() + cumsum(c(0, diff(group_id) != 0)) * 0.7)
mp_x_lookup <- setNames(mp_group_df$x, mp_group_df$mp)

# ==============================================================================
# FIGURE 1: UMAP & Vertical Barplot (from publication_umap.R)
# ==============================================================================
cat("Generating cell type UMAP and barplot...\n")
meta_file <- file.path(SCREF_REF_OUTS_DIR, "all_meta.rds")
atlas_seurat_file <- file.path(SCREF_REF_OUTS_DIR, "EAC_Ref_merged.rds")
plot_data <- NULL
if (file.exists(meta_file)) {
  all_meta <- readRDS(meta_file)
  if (all(c("umap_1", "umap_2", "celltype_update") %in% colnames(all_meta))) {
    plot_data <- data.frame(UMAP_1 = all_meta$umap_1, UMAP_2 = all_meta$umap_2,
                            celltype = all_meta$celltype_update)
  } else if (all(c("UMAP_1", "UMAP_2", "celltype_final") %in% colnames(all_meta))) {
    plot_data <- data.frame(UMAP_1 = all_meta$UMAP_1, UMAP_2 = all_meta$UMAP_2,
                            celltype = all_meta$celltype_final)
  }
}
if (is.null(plot_data) && file.exists(atlas_seurat_file)) {
  atlas <- readRDS(atlas_seurat_file)
  emb <- Seurat::Embeddings(atlas, "umap")
  plot_data <- data.frame(
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    celltype = atlas@meta.data[rownames(emb), "celltype_update", drop = TRUE],
    stringsAsFactors = FALSE
  )
  rm(atlas)
  gc()
}

if (!is.null(plot_data)) {
  plot_data <- plot_data |>
    filter(!is.na(celltype), !celltype %in% c("unresolved", "unresolved_inconsistent"))
  my_colors <- c(
    "b.cell" = "#E41A1C", "dendritic" = "#377EB8", "endothelial" = "#4DAF4A",
    "epithelial" = "#984EA3", "fibroblast" = "#FF7F00", "macrophage" = "#A65628",
    "mast" = "#F781BF", "nk.cell" = "#FFD700", "plasma" = "#999999", "t.cell" = "#00CED1"
  )
  set.seed(42)
  plot_data <- plot_data[sample(nrow(plot_data)), ]
  celltype_levels <- names(sort(table(plot_data$celltype), decreasing = TRUE))
  plot_data$celltype <- factor(plot_data$celltype, levels = celltype_levels)

  umap_p <- ggplot(plot_data, aes(UMAP_1, UMAP_2, color = celltype)) +
    geom_scattermore(pointsize = 4.0, alpha = 1, pixels = c(2048, 2048)) +
    scale_color_manual(values = my_colors[celltype_levels], name = "Cell Type") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    theme_void()

  bar_data <- plot_data |> count(celltype) |> mutate(pct = n / sum(n) * 100)
  bar_data$celltype <- factor(bar_data$celltype, levels = celltype_levels)
  max_n <- max(bar_data$n)

  bar_p <- ggplot(bar_data, aes(x = celltype, y = n)) +
    geom_col(aes(fill = celltype), width = 0.7, colour = "black", linewidth = 0.3) +
    geom_text(aes(y = n + max_n * 0.04, label = sprintf("%s\n(%.1f%%)", format(n, big.mark=","), pct)), 
              vjust = 0, size = 4.0, lineheight = 0.8) +
    scale_fill_manual(values = my_colors, name = "Cell Type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_blank(),
          legend.position = "none")

  umap_p <- umap_p + 
    theme(legend.position = "right",
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))
    
  umap_no_leg <- umap_p + theme(legend.position = "none") +
    plot_annotation(
      title = "OAC single cell atlas",
      subtitle = paste0("Final cell counts: ", format(nrow(plot_data), big.mark = ",")),
      theme = theme(plot.title = element_text(size = 20, face = "bold"),
                    plot.subtitle = element_text(size = 16))
    )
    
  leg <- cowplot::get_legend(umap_p)
  leg_plot <- cowplot::ggdraw(leg)
  
  save_pub_gg(umap_no_leg, section, "s1_atlas_celltype_umap_only", width = 8, height = 8, dpi = 600)
  save_pub_gg(bar_p, section, "s1_atlas_celltype_barplot", width = 7.5, height = 8, dpi = 600)
  save_pub_gg(leg_plot, section, "s1_atlas_celltype_legend", width = 4, height = 8, dpi = 600)
} else {
  abort_missing_figure(section, "s1_atlas_celltype_umap", "Atlas UMAP", "required atlas metadata input is unavailable")
  abort_missing_figure(section, "s1_atlas_celltype_barplot", "Atlas Barplot", "required atlas metadata input is unavailable")
}

# ==============================================================================
# FIGURE 2: NMF Clustering Heatmap (from geneNMF.R)
# ==============================================================================
cat("Generating NMF clustering heatmap...\n")
nmf_file <- file.path(SCREF_REF_OUTS_DIR, "Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")
if (file.exists(nmf_file)) {
  geneNMF.metaprograms <- readRDS(nmf_file)

  sil <- geneNMF.metaprograms$metaprograms.metrics$silhouette
  names(sil) <- paste0("MP", seq_along(sil))
  retained <- names(sil)[sil >= 0]

  programme_order <- geneNMF.metaprograms$programs.tree$order
  sim <- geneNMF.metaprograms$programs.similarity[programme_order, programme_order]
  clusters <- geneNMF.metaprograms$programs.clusters[programme_order]
  cluster_mps <- paste0("MP", clusters)
  runs <- rle(cluster_mps)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1
  retained_runs <- data.frame(
    mp = runs$values,
    start = run_start,
    end = run_end,
    stringsAsFactors = FALSE
  ) |>
    filter(mp %in% retained)
  nmf_label_desc <- c(
    "G2M Cell Cycle" = "G2M cycle",
    "G1S Cell Cycle" = "G1S cycle",
    "MYC-related Proliferation" = "MYC prolif.",
    "Basal-like Transition" = "Basal-like",
    "Hypoxia Adapted Epi." = "Hypoxia",
    "Epithelial IFN Resp." = "IFN resp.",
    "Columnar Diff." = "Columnar diff.",
    "Intestinal Diff." = "Intestinal diff.",
    "Hypoxic Inflam. Epi." = "Hypoxic inflam.",
    "DNA Damage Repair" = "DNA repair",
    "Secretory Diff. (Intest.)" = "Secretory int.",
    "Secretory Diff. (Gastric)" = "Secretory gastric",
    "Immune Infiltration" = "Immune infil.",
    "Neuro-responsive Epi." = "Neuro-resp."
  )
  desc <- SCREF_MP_DESCRIPTIONS[retained_runs$mp]
  retained_runs$label <- paste0(retained_runs$mp, ": ", desc)
  n_prog <- nrow(sim)

  col_fun <- colorRamp2(c(0, 0.16, 0.42), c("#FFF7F3", "#FB6A4A", "#67000D"))

  ht <- Heatmap(
    sim,
    name = "Jaccard",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,
    column_title = paste0("NMF Programmes (n = ", n_prog, ")"),
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 18, fontface = "bold"),
                                labels_gp = gpar(fontsize = 16),
                                direction = "horizontal",
                                legend_width = unit(8, "cm")),
    width = unit(160, "mm"),
    height = unit(160, "mm")
  )

  draw_expr <- quote({
    draw(ht, heatmap_legend_side = "bottom", padding = unit(c(18, 2, 34, 180), "mm"))
    decorate_heatmap_body("Jaccard", {
      raw_y <- 1 - ((retained_runs$start + retained_runs$end - 1) / (2 * n_prog))
      label_y <- raw_y
      
      # Grouping to enforce specific gaps between labels
      label_group <- c(1, 1, 2, 3, 3, 3, 3, 3, 4, 5, 6, 6, 7, 8)
      
      for (iter in 1:200) {
        for (i in seq_len(length(label_y) - 1)) {
          required_dist <- ifelse(label_group[i] == label_group[i+1], 0.035, 0.08)
          if (label_y[i] - label_y[i+1] < required_dist) {
            overlap <- required_dist - (label_y[i] - label_y[i+1])
            label_y[i] <- min(0.98, label_y[i] + overlap * 0.3)
            label_y[i+1] <- max(0.02, label_y[i+1] - overlap * 0.7)
          }
        }
      }
      label_y <- pmax(0.015, pmin(0.985, label_y))
      
      for (rr in seq_len(nrow(retained_runs))) {
        s <- retained_runs$start[rr]
        e <- retained_runs$end[rr]
        centre <- (s + e - 1) / (2 * n_prog)
        extent <- (e - s + 1) / n_prog
        box_y <- 1 - centre
        ly <- label_y[rr]
        grid.rect(x = unit(centre, "npc"), y = unit(1 - centre, "npc"),
                  width = unit(extent, "npc"), height = unit(extent, "npc"),
                  gp = gpar(fill = NA, col = "#111111", lwd = 2.0, lty = "dotted"))
        grid.lines(x = unit(c(centre + extent / 2, 1.008), "npc"), y = unit(c(box_y, ly), "npc"),
                   gp = gpar(col = "#111111", lwd = 2.0, lty = "dotted"))
        grid.text(retained_runs$label[rr], x = unit(1.012, "npc"), y = unit(ly, "npc"),
                  just = "left", gp = gpar(fontsize = 14, col = "#111111", lineheight = 0.84))
      }
    })
  })

  save_pub_grid(draw_expr, section, "s1_nmf_clustering_heatmap", width = 14, height = 8.5)
} else {
  abort_missing_figure(section, "s1_nmf_clustering_heatmap", "NMF Clustering Heatmap", "NMF Results missing")
}

# ==============================================================================
# FIGURE 3: Enrichment Annotation (from enrichment_annotation.R)
# ==============================================================================
cat("Generating enrichment heatmaps...\n")
enrich_file <- file.path(SCREF_REF_OUTS_DIR, "cluster_enrich.rds")
if (file.exists(enrich_file)) {
  enrich <- readRDS(enrich_file)
  
  # Helper to extract a specific database across all MPs
  extract_db <- function(enrich_list, db_name) {
    do.call(bind_rows, lapply(names(enrich_list), function(mp) {
      res <- enrich_list[[mp]][[db_name]]
      if (!is.null(res)) {
        df <- as.data.frame(res)
        if (nrow(df) > 0) {
          df$Cluster <- mp
          return(df)
        }
      }
      return(NULL)
    }))
  }
  
  enrich_3ca <- extract_db(enrich, "MPs_3CA")
  enrich_early <- extract_db(enrich, "Early_Embryogenesis")
  enrich_dev_long <- extract_db(enrich, "Normal_Development_long")
  enrich_adult <- extract_db(enrich, "Adult_Epithelium")
  enrich_barretts <- extract_db(enrich, "Barretts_Oesophagus")
  
  build_enrich_plot <- function(df, filename, palette_high, width, height) {
    if (nrow(df) == 0) {
      abort_missing_figure(section, filename, filename, "No requested enrichment rows found.")
      return(invisible(NULL))
    }
    df <- df |>
      mutate(MP = factor(Cluster, levels = PUB_MP_ORDER),
             x = mp_x_lookup[Cluster],
             state_group = factor(pub_mp_state(Cluster), levels = PUB_MP_STATE_ORDER),
             term = factor(term, levels = rev(unique(term))),
             group = factor(group, levels = unique(group)),
             neglogP = pmin(neglogP, 10))
    term_levels <- levels(df$term)
    plot_term_levels <- c(term_levels, " ")
    full_df <- tidyr::expand_grid(term = term_levels, MP = mp_group_df$mp) |>
      left_join(df |> select(group, term, MP, neglogP), by = c("term", "MP")) |>
      left_join(mp_group_df |> select(MP = mp, x, state_group = group), by = "MP") |>
      group_by(term) |> fill(group, .direction = "downup") |> ungroup() |>
      mutate(neglogP = ifelse(is.na(neglogP), 0, neglogP),
             term = factor(as.character(term), levels = plot_term_levels))
    top_group <- levels(df$group)[1]
    ann <- tidyr::expand_grid(group = factor(top_group, levels = levels(df$group)),
                              mp_group_df |> rename(mp_state = group)) |>
      mutate(term = factor(" ", levels = plot_term_levels))
    p <- ggplot(full_df, aes(x, term, fill = neglogP)) +
      geom_tile(width = 0.92, height = 0.92, colour = "white", linewidth = 0.18) +
      geom_tile(data = ann, aes(x = x, y = term),
                inherit.aes = FALSE, width = 0.92, height = 0.18,
                fill = PUB_MP_STATE_COLOURS[as.character(ann$mp_state)], colour = NA) +
      facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_fill_gradient(low = "white", high = "#B2182B", name = "-log10(FDR)",
                          limits = c(0, max(full_df$neglogP, na.rm = TRUE))) +
      scale_x_continuous(breaks = mp_group_df$x, 
                         labels = paste0(mp_group_df$mp, " ", SCREF_MP_DESCRIPTIONS[mp_group_df$mp]),
                         expand = expansion(mult = c(0.01, 0.01)),
                         position = "top") +
      scale_y_discrete(labels = function(x) ifelse(x == " ", "", x)) +
      labs(x = NULL, y = NULL) +
      coord_cartesian(clip = "off") +
      pub_theme(14) +
      theme(axis.text.x.top = element_text(angle = 30, hjust = 0, vjust = 0, size = 16, face = "bold"),
            axis.text.y = element_text(size = 16),
            strip.placement = "outside",
            strip.background = element_rect(fill = "#F1F5F9", colour = "#CBD5E1"),
            strip.text.y.left = element_text(angle = 0, size = 14, face = "bold"),
            panel.spacing.y = unit(3, "mm"),
            legend.position = "right",
            legend.title = element_text(size = 18, face = "bold"),
            legend.text = element_text(size = 16),
            legend.key.size = unit(2, "lines"),
            plot.margin = margin(t = 10, r = 80, b = 10, l = 10))
    save_pub_gg(p, section, filename, width = width + 2, height = height + 2)
    write_csv(df, file.path(out_dir, "tables", paste0(filename, ".csv")))
  }

  df_4a <- list()
  if (!is.null(enrich_3ca)) {
    df_4a[["3CA"]] <- enrich_3ca |>
      filter(Cluster %in% PUB_MP_ORDER,
             Description %in% c("3CA_mp3 Cell Cylce HMG-rich",
                                "3CA_mp7 Hypoxia",
                                "3CA_mp31 Respiration 2")) |>
      group_by(Cluster, Description) |>
      summarise(neglogP = max(-log10(p.adjust), na.rm = TRUE), .groups = "drop") |>
      mutate(group = "3CA Pan-cancer", term = pub_clean_term(Description))
  }
  if (!is.null(enrich_early)) {
    df_4a[["early"]] <- enrich_early |>
      filter(Cluster %in% PUB_MP_ORDER, grepl("PriS|primitive", Description, ignore.case = TRUE)) |>
      group_by(Cluster, Description) |>
      summarise(neglogP = max(-log10(p.adjust), na.rm = TRUE), .groups = "drop") |>
      mutate(group = "Embryogenesis", term = "Primitive streak")
  }
  if (!is.null(enrich_dev_long)) {
    df_4a[["dev"]] <- enrich_dev_long |>
      filter(Cluster %in% PUB_MP_ORDER,
             Description %in% c("MUC13_DMBT1_positive_cells_Stomach..Normal_Development_long",
                                "Lymphoid_cells_Stomach..Normal_Development_long",
                                "Squamous_epithelial_cells_Stomach..Normal_Development_long")) |>
      group_by(Cluster, Description) |>
      summarise(neglogP = max(-log10(p.adjust), na.rm = TRUE), .groups = "drop") |>
      mutate(group = "Development\n(Fetal Stomach)", term = pub_clean_term(Description))
  }

  df_4b <- list()
  if (!is.null(enrich_adult)) {
    adult_base <- enrich_adult |>
      filter(Cluster %in% PUB_MP_ORDER) |>
      mutate(neglogP = -log10(p.adjust))
    df_4b[["stomach"]] <- adult_base |>
      filter(grepl("PylG|GKN", Description, ignore.case = TRUE),
             !grepl("PG/Neck", Description, ignore.case = TRUE)) |>
      group_by(Cluster, Description) |> summarise(neglogP = max(neglogP), .groups = "drop") |>
      mutate(group = "Adult stomach", term = pub_clean_term(Description))
    df_4b[["gastric_im"]] <- adult_base |>
      filter(grepl("IntestMeta", Description, ignore.case = TRUE)) |>
      group_by(Cluster, Description) |> summarise(neglogP = max(neglogP), .groups = "drop") |>
      mutate(group = "Gastric IM", term = pub_clean_term(Description))
  }
  if (!is.null(enrich_barretts)) {
    barretts_base <- enrich_barretts |>
      filter(Cluster %in% PUB_MP_ORDER) |>
      mutate(neglogP = -log10(p.adjust))
    df_4b[["barretts"]] <- barretts_base |>
      filter(grepl("Columnar_Undifferentiated|Undifferentiated_Normal|Goblet|Columnar_Intermediate",
                   Description, ignore.case = TRUE)) |>
      group_by(Cluster, Description) |> summarise(neglogP = max(neglogP), .groups = "drop") |>
      mutate(group = "Barrett's epithelium", term = pub_clean_term(Description))
    df_4b[["smg"]] <- barretts_base |>
      filter(grepl("Mucous_Submucosal_Glands|Duct_Intercalating_Submucosal_Glands|Oncocytes_Submucosal_Glands",
                   Description, ignore.case = TRUE)) |>
      group_by(Cluster, Description) |> summarise(neglogP = max(neglogP), .groups = "drop") |>
      mutate(group = "Submucosal gland", term = pub_clean_term(Description))
  }
  
  combined_df <- bind_rows(bind_rows(df_4a), bind_rows(df_4b)) |>
    group_by(group, term, Cluster) |>
    summarise(neglogP = max(neglogP), .groups = "drop") |>
    mutate(term = sub("[ _]*Stomach[ _]*", "", term, ignore.case = TRUE)) |>
    mutate(term = sub("[ _]*Barrett'?s?[ _]*[Oo]?esophagus[ _]*", "", term, ignore.case = TRUE)) |>
    mutate(term = sub("[ _]*Submucosal[ _]*Gland(s)?[ _]*", "", term, ignore.case = TRUE)) |>
    mutate(term = trimws(term)) |>
    mutate(group = factor(group, levels = c("3CA Pan-cancer", "Embryogenesis", "Development\n(Fetal Stomach)",
                                            "Adult stomach", "Gastric IM",
                                            "Barrett's epithelium", "Submucosal gland"))) |>
    arrange(group, term)
    
  build_enrich_plot(combined_df, "s1_combined_enrichment", "#B2182B", 14.5, 6.5)
} else {
  abort_missing_figure(section, "s1_combined_enrichment", "Combined Enrichment", "required enrichment input is unavailable")
}

state_file <- file.path(SCREF_REF_OUTS_DIR, "Auto_topmp_v2_noreg_states_B.rds")

# ==============================================================================
# FIGURE 4: State per sample abundance plot (from state_mp_sample_abundance.R)
# ==============================================================================
cat("Generating state per sample abundance plot...\n")
epi_meta_file <- file.path(SCREF_REF_OUTS_DIR, "meta_full_epi.rds")
if (file.exists(epi_meta_file) && file.exists(state_file)) {
  meta <- readRDS(epi_meta_file)
  states <- readRDS(state_file)
  
  common <- intersect(rownames(meta), names(states))
  meta <- meta[common, ]
  meta$state <- states[common]
  
  ab_states <- c(PUB_STATE_ORDER, "Unresolved", "Hybrid")
  df_ab <- meta |> 
    filter(state %in% ab_states) |>
    count(orig.ident, state) |>
    group_by(orig.ident) |>
    mutate(pct = n / sum(n) * 100) |>
    ungroup()
  
  # Order samples by dominant state
  order_df <- df_ab |> 
    mutate(state_num = as.integer(factor(state, levels = ab_states))) |>
    group_by(orig.ident) |>
    summarise(dom = sum(pct * state_num)) |>
    arrange(dom)
  
  df_ab$orig.ident <- factor(df_ab$orig.ident, levels = order_df$orig.ident)
  df_ab$state <- factor(df_ab$state, levels = rev(ab_states))
  
  p_ab_main <- ggplot(df_ab, aes(x = orig.ident, y = pct, fill = state)) +
    geom_col(width = 1) +
    scale_fill_manual(values = SCREF_STATE_COLOURS) +
    labs(x = "Tumour Samples", y = "Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "none"
    )

  # Create pie chart
  df_pie <- df_ab |>
    group_by(state) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(pct = n / sum(n) * 100) |>
    arrange(desc(state)) |>
    mutate(ypos = cumsum(pct) - 0.5 * pct)

  p_pie <- ggplot(df_pie, aes(x = 1, y = pct, fill = state)) +
    geom_col(width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    xlim(c(0.5, 1.5)) +
    geom_text(aes(x = 1.2, y = ypos, label = ifelse(pct > 2, sprintf("%.1f%%", pct), ""),
                  color = ifelse(state == "Hybrid", "white", "black")), 
              size = 4.5, fontface = "bold") +
    scale_fill_manual(values = SCREF_STATE_COLOURS, name = "State") +
    scale_color_identity() +
    theme_void() +
    theme(
      legend.position = "top", 
      legend.title = element_text(size = 18, face = "bold", margin = margin(b = 10)),
      legend.text = element_text(size = 17, margin = margin(r = 15)),
      legend.spacing.x = unit(0.3, "cm"),
      legend.justification = "left",
      legend.box.just = "left",
      legend.key.size = unit(1.5, "lines"),
      legend.margin = margin(t = 15, l = 10),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    guides(fill = guide_legend(ncol = 2, reverse = TRUE, title.position = "top"))

  tmp <- ggplot_gtable(ggplot_build(p_pie))
  leg_idx <- which(grepl("guide-box", sapply(tmp$grobs, function(x) x$name)))
  leg_grob <- tmp$grobs[[leg_idx[1]]]
  
  p_pie_noleg <- p_pie + theme(legend.position = "none")
  
  # Stack legend over the pie chart to create a left column
  left_col <- cowplot::plot_grid(cowplot::ggdraw(leg_grob), p_pie_noleg, ncol = 1, rel_heights = c(0.6, 1.3))
  
  # Combine left column with the sample abundance plot
  combined_ab <- cowplot::plot_grid(left_col, p_ab_main, nrow = 1, rel_widths = c(2.8, 5))

  save_pub_gg(combined_ab, section, "s1_state_abundance_per_sample", width = 19, height = 5)
} else {
  abort_missing_figure(section, "s1_state_abundance_per_sample", "State Abundance", "required final-state input is unavailable")
}

cat("Section 1 complete.\n")
