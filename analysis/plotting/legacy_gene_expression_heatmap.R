####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/plotting/legacy_gene_expression_heatmap.R
#   Methodology: analysis/methodology/plotting/publication_plotting_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Moved from: analysis/gene_expr_compare.R
# Reorganized as part of analysis/ restructuring
####################
genes_mp1 <- geneNMF.metaprograms$metaprograms.genes[["MP1"]]
genes_mp5 <- geneNMF.metaprograms$metaprograms.genes[["MP5"]]
genes_mp7 <- geneNMF.metaprograms$metaprograms.genes[["MP7"]]
genes_to_exclude <- unique(c(genes_mp1, genes_mp7))
selected_genes <- setdiff(genes_mp5, genes_to_exclude)
discarded_genes <- intersect(genes_mp5, genes_to_exclude)
gene_sets <- list(
  "Cell cycle genes MP5" = discarded_genes, 
  "Other genes MP5" = selected_genes
)
gene2module <- setNames(
  rep(names(gene_sets), times = lengths(gene_sets)),
  unlist(gene_sets, use.names = FALSE)
)
features <- names(gene2module)
DefaultAssay(tmdata_all) <- "RNA"
features_present <- features[features %in% rownames(tmdata_all)]

library(ComplexHeatmap)
library(circlize)
library(Seurat)
# Define the states you want to compare and their order
state_order <- c("Classic Proliferative", "Basal to Intest. Meta", "Intestinal Metaplasia", "Stress-adaptive", "Unassigned/Quienscent", "Unassigned")

manual_names <- levels(tmdata_all$manual_state)
manual_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Intestinal Metaplasia" = "#984EA3",
  "Stress-adaptive"       = "#FF7F00",
  "Unassigned/Quienscent" = "grey80",
  "Unassigned"            = "grey50"
)


# --- 3. Prepare Matrix & Rows ---
# Extract expression matrix

features <- unlist(gene_sets)
features_present <- features[features %in% rownames(tmdata_all)]
mat <- as.matrix(GetAssayData(tmdata_all, slot = "data")[features_present, ])

# Create Row Split (Gene Modules)
# Assuming 'gene2module' exists and maps genes to module names
row_split <- factor(gene2module[features_present], levels = names(gene_sets))

# Define Color Function
col_fun <- colorRamp2(c(0, 1.5, 3, 6), c("#FCFDBF", "#FEB078", "#B73779", "#000004"))

# --- 4. Annotation (Top Bars) ---
top_anno <- HeatmapAnnotation(
  # Bar 1: The State (Main Grouping)
  State = tmdata_all$manual_state,
  
  # Bar 2: The Sample (To see patient composition within each state)
  Sample = tmdata_all$orig.ident,
  
  col = list(
    State = manual_cols
  ),
  
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm"),
  gap = unit(1, "mm")
)

# --- 5. Generate Heatmap ---
ht <- Heatmap(
  mat,
  name = "Log-Expr",
  col = col_fun,
  
  # --- KEY CHANGE: SPLIT BY STATE ---
  column_split = tmdata_all$manual_state, 
  column_gap = unit(2, "mm"), # Clean gap between states
  
  column_title_rot = 30, # Keep state names horizontal
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_side = "top",
  
  # Row Settings
  row_split = row_split,
  row_gap = unit(2.5, "mm"),
  cluster_rows = FALSE,      # Keep genes in your set order
  cluster_columns = TRUE,    # Cluster cells within each state
  show_column_dend = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  # Style
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_names_gp = gpar(fontsize = 9, fontitalic = TRUE),
  top_annotation = top_anno,
  border = TRUE, 
  use_raster = TRUE, 
  raster_quality = 5
)

# --- 6. Save ---
pdf("Heatmap_temp.pdf", width = 16, height = 8)
draw(ht, 
     merge_legends = TRUE, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm")) # Adjusted padding
dev.off()


#######################################

# ---- 1) Define modules ----
# We keep them separate first to handle the overlap
squamous_genes <- c("CSTA", "SPINK5", "ZNF185", "KRT80", "NDRG2", "KRT17", "KRT15")
mp9_genes      <- geneNMF.metaprograms$metaprograms.genes$MP9

# ---- 2) Remove squamous genes from MP9 list (Assignment Priority) ----
# This ensures a gene is ONLY in one category
mp9_genes_unique <- setdiff(mp9_genes, squamous_genes)

gene_sets <- list(
  "Squamous epithelial cells" = squamous_genes,
  "Other MP9 genes"           = mp9_genes_unique
)

# ---- 3) Create a named vector: gene -> module ----
gene2module <- setNames(
  rep(names(gene_sets), times = lengths(gene_sets)),
  unlist(gene_sets, use.names = FALSE)
)

# Verify no duplicates
if(any(duplicated(names(gene2module)))) stop("Duplicate genes found across modules!")

# Ordered feature vector
features <- names(gene2module)

# ---- 4) Keep only genes present in the object ----
DefaultAssay(tmdata_all) <- "RNA"
features_present <- features[features %in% rownames(tmdata_all)]

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 1. Aesthetic Color Setup (Keep full list here)
all_sample_colors <- c(
  "SUR1072_Untreated_PDO" = "#BDD7EE", "SUR1072_Treated_PDO" = "#2171B5",
  "SUR1090_Untreated_PDO" = "#FDBE85", "SUR1090_Treated_PDO" = "#D94701",
  "SUR1070_Untreated_PDO" = "#E1BEE7", "SUR1070_Treated_PDO" = "#8E24AA",
  "SUR1181_Untreated_PDO" = "#F8BBD0", "SUR1181_Treated_PDO" = "#C2185B"
)

# 2. Subset & DYNAMIC CLEANUP
#tmdata_all4 <- subset(tmdata_all, subset = orig.ident %in% names(all_sample_colors))
tmdata_all4 <- subset(tmdata_all, subset = orig.ident %in% names(all_sample_colors) & manual_state == "Basal to Intest. Meta")

# --- DROP-IN FIX: Only keep samples that exist in the data ---
existing_samples <- intersect(names(all_sample_colors), unique(tmdata_all4$orig.ident))
sample_colors <- all_sample_colors[existing_samples]

# Re-factor to drop empty levels
tmdata_all4$orig.ident <- factor(tmdata_all4$orig.ident, levels = existing_samples)

# Create Patient Metadata based on existing samples
tmdata_all4$Patient <- gsub("_.*", "", tmdata_all4$orig.ident)
existing_patients <- unique(tmdata_all4$Patient)
tmdata_all4$Patient <- factor(tmdata_all4$Patient, levels = intersect(c("SUR1072", "SUR1090", "SUR1070", "SUR1181"), existing_patients))

# 3. Handle Rows (Genes) - Only if present
features <- unlist(gene_sets)
features_present <- features[features %in% rownames(tmdata_all4)]
row_split_vector <- rep(names(gene_sets), lengths(gene_sets))
names(row_split_vector) <- features
row_split <- factor(gene2module[features_present], levels = names(gene_sets))

# 4. Matrix & Color
mat <- as.matrix(GetAssayData(tmdata_all4, slot = "data")[features_present, ])
col_fun <- colorRamp2(c(0, 1.5, 3, 6), c("#FCFDBF", "#FEB078", "#B73779", "#000004"))

sample_levels <- levels(tmdata_all4$orig.ident)

# 1. Map your existing manual_cols to the annotation
# We use the specific color vector you provided
manual_names <- levels(tmdata_all$manual_state)
manual_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Intestinal Metaplasia" = "#984EA3",
  "Stress-adaptive"       = "#FF7F00",
  "Unassigned/Quienscent" = "grey80",
  "Unassigned"            = "grey50"
)

# 2. Update Top Annotation with matching colors
top_anno <- HeatmapAnnotation(
  # First Bar: Sample Color (using existing sample_colors)
  Sample_Color = tmdata_all4$orig.ident,
  
  # Second Bar: State (using your manual_cols mapping)
  State = tmdata_all4$manual_state,
  
  col = list(
    Sample_Color = sample_colors,
    State = manual_cols # Now perfectly matches your "manual_cols" scheme
  ),
  
  show_legend = c(Sample_Color = FALSE, State = TRUE),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm"),
  gap = unit(1, "mm")
)

n_slices <- length(levels(tmdata_all4$orig.ident))
col_gaps <- rep(1.5, n_slices - 1)
for(i in 1:(n_slices-1)) {
  # If the patient ID changes between slice i and i+1, make the gap larger
  p1 <- gsub("_.*", "", levels(tmdata_all4$orig.ident)[i])
  p2 <- gsub("_.*", "", levels(tmdata_all4$orig.ident)[i+1])
  if(p1 != p2) col_gaps[i] <- 5
}

# 3. Final Heatmap Call
ht <- Heatmap(
  mat,
  name = "Log-Expr",
  col = col_fun,
  
  # Grouping & Titles
  column_split = tmdata_all4$orig.ident,
  column_gap = unit(col_gaps, "mm"),
  column_title_rot = 20,
  column_title_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title_side = "top",
  
  # Structure
  row_split = row_split,
  row_gap = unit(2.5, "mm"),
  cluster_rows = FALSE, 
  cluster_columns = TRUE, 
  show_column_dend = FALSE, 
  show_column_names = FALSE, 
  show_row_names = TRUE,
  
  # Style
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_names_gp = gpar(fontsize = 9, fontitalic = TRUE),
  top_annotation = top_anno,
  border = TRUE, 
  use_raster = TRUE, 
  raster_quality = 5
)

# 4. Save to PDF
pdf("Heatmap_Matched_Colors_Dediff.pdf", width = 16, height = 8)
draw(ht, 
     merge_legends = TRUE, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 2, 15, 2), "mm"))
dev.off()

########################################

selected <- rownames(
  tmdata_pdos@meta.data[
    tmdata_pdos$orig.ident %in% c(
      "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
      "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
      "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
      "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
    ),
  ]
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# --- 1) Compute averages (KEEP MP column names intact) ---
df <- tmdata_pdos@meta.data[selected, ]

mp_cols <- mp_tree_order  # c("MP1","MP7",...)
avg_by_sample <- df %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(mp_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(
    SUR_id = str_extract(orig.ident, "^SUR\\d+"),
    Condition = case_when(
      str_detect(orig.ident, "Untreated") ~ "Untreated",
      str_detect(orig.ident, "Treated")   ~ "Treated",
      TRUE ~ NA_character_
    ),
    Condition = factor(Condition, levels = c("Untreated", "Treated"))
  )

# --- 2) Define x-axis order: all Untreated first (by SUR), then all Treated (by SUR) ---
sur_order <- avg_by_sample %>%
  distinct(SUR_id) %>%
  arrange(SUR_id) %>%
  pull(SUR_id)

x_levels <- c(
  avg_by_sample %>% filter(Condition == "Untreated") %>% arrange(factor(SUR_id, levels = sur_order)) %>% pull(orig.ident),
  avg_by_sample %>% filter(Condition == "Treated")   %>% arrange(factor(SUR_id, levels = sur_order)) %>% pull(orig.ident)
)

# --- 3) Long format for ggplot ---
bar_df <- avg_by_sample %>%
  mutate(orig.ident = factor(orig.ident, levels = x_levels)) %>%
  pivot_longer(cols = all_of(mp_cols), names_to = "MP", values_to = "score") %>%
  mutate(MP = factor(MP, levels = mp_cols))

# --- 4) Make two related palettes (base R) and assign a unique color per sample ---
samples_unt <- avg_by_sample %>% filter(Condition == "Untreated") %>% arrange(factor(SUR_id, levels = sur_order)) %>% pull(orig.ident)
samples_trt <- avg_by_sample %>% filter(Condition == "Treated")   %>% arrange(factor(SUR_id, levels = sur_order)) %>% pull(orig.ident)

pal_unt <- grDevices::hcl.colors(length(samples_unt)+2, palette = "Blues 3")
pal_trt <- grDevices::hcl.colors(length(samples_trt)+2, palette = "Reds 3")

sample_colors <- c(setNames(pal_unt, samples_unt), setNames(pal_trt, samples_trt))

# --- 5) Optional: pretty facet labels using mp_descriptions named vector ---
# mp_descriptions should be like c(MP1="...", MP2="...", ...)
facet_labs <- setNames(
  mp_descriptions[mp_cols],
  mp_cols
)

p_bar <- ggplot(bar_df, aes(x = orig.ident, y = score, fill = orig.ident)) +
  geom_col(width = 0.82, color = "black", linewidth = 0.25) +
  facet_wrap(~ MP, ncol = 5, scales = "free_y", labeller = as_labeller(facet_labs)) +
  scale_fill_manual(values = sample_colors, guide = "none") +
  labs(x = NULL, y = "Mean MP score", title = "Mean MP scores per sample") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.spacing = unit(0.9, "lines")
  )

p_bar
