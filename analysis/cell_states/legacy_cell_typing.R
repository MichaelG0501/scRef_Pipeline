####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/cell_states/legacy_cell_typing.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/celltyping.R
# Reorganized as part of analysis/ restructuring
####################
library(readxl)
library(dplyr)
library(purrr)
library(stringr)

# Load marker data
markers <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx", sheet = 1)

# Cell type mapping
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

combine_marker_scores <- function(df, w_specificity = 0.3, w_sensitivity = 0.7) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }
  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) / 
    (w_specificity + w_sensitivity)
  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

markers <- markers[markers$specificity > 0.4 & markers$cell_type != "Malignant", ]

markers_list <- markers %>%
  mutate(cell_type = recode(cell_type, !!!celltype_map)) %>%
  split(.$cell_type)

markers_ranked <- lapply(markers_list, function(df) {
  combine_marker_scores(df, w_specificity = 0.3, w_sensitivity = 0.7)
})

########################################################

library(dplyr)
library(purrr)
library(tidyr)
library(Seurat)
library(ggplot2)

N <- 100
lfc_th   <- 1
z_cut <- 1

# ----------------- Step-1: DE overlap enrichment -----------------
Idents(tmdata) <- tmdata$seurat_clusters
de <- FindAllMarkers(tmdata, only.pos = TRUE, min.pct = 0.1, logfc.threshold = lfc_th)

setsN <- markers_ranked |>
  imap(~ .x %>% arrange(desc(Combined)) %>% slice_head(n = N) %>% pull(gene) %>% intersect(rownames(tmdata)))

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
  sig <- res %>%
    filter(padj <= 0.05) %>%
    arrange(desc(overlap), padj)   # overlap primary, then padj
  
  if (nrow(sig) == 0) {
    return(tibble(cluster = cl, step1 = "unknown"))
  }
  top_ov <- sig$overlap[1]
  if (nrow(sig) == 1) {
    keep <- sig %>% slice(1)
  } else {
    second_ov <- sig$overlap[2]
    if ((top_ov - second_ov) > 5) {
      keep <- sig %>% slice(1)
    } else {
      keep <- sig %>% filter((top_ov - overlap) <= 5)
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
  pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
  mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
  group_by(cluster) %>%
  mutate(
    z = as.numeric(scale(score)),                      # z-score across celltypes in this cluster
  ) %>%
  ungroup() %>%
  select(cluster, cell_type, score, z)

# ---- VALIDATION of Step-1 calls with module scores ----
step1_long <- step1_calls %>%
  mutate(cell_type = strsplit(step1, "\\|")) %>%
  unnest(cell_type)

# join z-scores for only the candidates
val <- step1_long %>%
  left_join(scores_long, by = c("cluster", "cell_type")) %>%
  mutate(high = z >= z_cut)

# final: keep Step-1 if "high", else unresolved
final_calls <- val %>%
  group_by(cluster) %>%
  summarize(
    step1 = first(step1),  # original multi-call string for reference
    final = {
      cand <- cell_type[high & !is.na(high)]
      if (length(cand) == 0) {
        "unresolved"
      } else if (length(cand) == 1) {
        cand
      } else {
        ord <- order(z[high & !is.na(high)], decreasing = TRUE)
        paste(cand[ord], collapse = "|")
      }
    },
    .groups = "drop"
  )

 # attach to object
lab_map <- final_calls$final; names(lab_map) <- final_calls$cluster
tmdata$celltype_final <- as.character(lab_map[as.character(tmdata$seurat_clusters)])


plot_df <- cbind(
  as.data.frame(Embeddings(tmdata, "umap")), # UMAP_1 and UMAP_2
  tmdata@meta.data                          # All metadata, including scores and clusters
)
plot_df$seurat_clusters <- as.factor(plot_df$seurat_clusters)

plot_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), # Centered, bold
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11)
  )

mod_plots <- imap(mod_cols, function(mod_col_name, plot_title) {
  
  ggplot(plot_df, aes(x = umap_1, y = umap_2, color = .data[[mod_col_name]])) +
    geom_point(size = 0.5) +
    scale_color_gradientn(colors = c("grey90", "blue"), 
                          limits = c(-0.2, 0.8), 
                          oob = scales::squish) +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Score") +
    ggtitle(plot_title) + # Use the nice name as the title
    plot_theme # Apply the shared theme
})

umap_cluster <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = seurat_clusters)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
  ggtitle("UMAP by cluster") +
  plot_theme +
  # Make legend dots bigger so they are visible
  guides(color = guide_legend(override.aes = list(size = 3))) 

umap_cluster2 <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = celltype_final)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
  ggtitle("UMAP by cluster") +
  plot_theme +
  # Make legend dots bigger so they are visible
  guides(color = guide_legend(override.aes = list(size = 3))) 

combined_plot <- wrap_plots(c(mod_plots, list(umap_cluster), list(umap_cluster2)), ncol = 4)

combined_plot
