####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/basal_mp_distance_matrix.R
#   Methodology: analysis/methodology/cell_states/basal_smg_mp_distance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
####################

####################
# basal_mp_distance_matrix.R
# Monocle3-based MP-to-MP distance comparison for the six basal and five SMG
# MPs defining the current centred refined noreg states.
#
# Within-state MP hybrids use a deliberately stricter top-two gap <0.1.
#
# Methods (same as pseudotime_state_distance_matrix.R):
#   1. Directed pseudotime gaps (median)
#   2. Directed pseudotime gaps (mean)
#   3. Principal-graph geodesic distance (medoid)
#   4. Principal-graph geodesic distance (centroid)
#   5. UMAP centroid Euclidean distance
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#   ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds
#
# Output:
#   ref_outs/basal_mp_distance/intermediate/basal_mp_assignments.rds
#   ref_outs/basal_mp_distance/tables/basal_mp_distance_summary.csv
#   ref_outs/basal_mp_distance/tables/basal_mp_*_matrix.csv
#   ref_outs/basal_mp_distance/intermediate/basal_mp_distance_matrices.rds
#   ref_outs/basal_mp_distance/tables/basal_mp_proportions.csv
#   ref_outs/basal_mp_distance/figures/basal_mp_distance_heatmaps.pdf
#   ref_outs/basal_mp_distance/figures/basal_mp_pseudotime_samples.pdf
#   ref_outs/basal_mp_distance/figures/basal_mp_proportions.pdf
#   ref_outs/basal_mp_distance/intermediate/sample_trajectories/<sample>_basal_pseudotime.rds
####################

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(igraph)
library(patchwork)

setwd("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs")

# --- Output dirs ---
base_dir <- "basal_mp_distance"
int_dir  <- file.path(base_dir, "intermediate")
tab_dir  <- file.path(base_dir, "tables")
fig_dir  <- file.path(base_dir, "figures")
traj_dir <- file.path(int_dir, "sample_trajectories")
for (d in c(int_dir, tab_dir, fig_dir, traj_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Constants ---
basal_mps <- c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+")
mp_labels <- c(
  "MP14"  = "MP14 Squamoid/basal transition",
  "MP3+"  = "MP3+ Basal-columnar invasive epithelium",
  "MP6+"  = "MP6+ Stress-reactive columnar epithelium",
  "MP11+" = "MP11+ Epithelial antiviral interferon response",
  "MP9+"  = "MP9+ Metabolic columnar epithelium",
  "MP10+" = "MP10+ Intestinal metaplasia"
)
mp_cols <- c(
  "MP14" = "#4DAF4A", "MP3+" = "#8DA0CB", "MP6+" = "#66C2A5",
  "MP11+" = "#FC8D62", "MP9+" = "#A6D854", "MP10+" = "#E78AC3", Hybrid = "black"
)
HYBRID_GAP <- 0.1
ROOT_MIN <- 15
OTHER_MIN <- 15
TOTAL_MIN <- 50
MIN_GROUPS <- 2

# --- Helpers (from pseudotime_state_distance_matrix.R) ---
safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_median <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else median(x) }

make_empty_matrix <- function(labs) {
  m <- matrix(NA_real_, length(labs), length(labs), dimnames = list(labs, labs)); diag(m) <- 0; m
}

matrix_to_long <- function(mat, method_name) {

  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    as_tibble() %>% rename(mp_a = Var1, mp_b = Var2, distance = Freq) %>%
    mutate(method = method_name)
}

make_symmetric_matrix <- function(sdf, mname, labs) {
  out <- make_empty_matrix(labs)
  df_use <- sdf %>% filter(method == mname) %>%
    mutate(a = pmin(mp_a, mp_b), b = pmax(mp_a, mp_b)) %>%
    group_by(a, b) %>% summarise(distance = mean(mean_distance, na.rm = TRUE), .groups = "drop")
  for (i in seq_len(nrow(df_use))) { out[df_use$a[i], df_use$b[i]] <- df_use$distance[i]; out[df_use$b[i], df_use$a[i]] <- df_use$distance[i] }
  diag(out) <- 0; out
}

get_medoid_cell <- function(udf) {
  if (nrow(udf) == 1) return(udf$cell[1])
  dm <- as.matrix(dist(udf[, c("UMAP_1","UMAP_2")])); udf$cell[which.min(rowSums(dm))]
}

compute_graph_weights <- function(g, gc) {
  edf <- igraph::as_data_frame(g, "edges")
  if (!nrow(edf)) return(g)
  w <- purrr::map2_dbl(edf$from, edf$to, ~sqrt(sum((gc[,.x] - gc[,.y])^2)))
  igraph::E(g)$weight <- w; g
}

nearest_graph_vertex <- function(pt, gc) colnames(gc)[which.min(apply(gc, 2, function(n) sqrt(sum((pt-n)^2))))]

coerce_vertex <- function(rv, g) {
  gn <- igraph::V(g)$name; rv <- as.character(rv)[1]
  if (is.na(rv) || rv == "") return(NA_character_)
  if (rv %in% gn) return(rv)
  ri <- as.integer(suppressWarnings(as.numeric(rv)))
  if (!is.na(ri) && ri >= 1 && ri <= length(gn)) return(gn[ri])
  NA_character_
}

get_root_node <- function(cds, root_cells) {
  cv <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(cv)) return(NULL)
  if (is.matrix(cv)) cv <- cv[,1,drop=TRUE]
  root_cells <- intersect(root_cells, names(cv))
  if (!length(root_cells)) return(NULL)
  rv <- names(sort(table(as.character(cv[root_cells])), decreasing=TRUE))[1]
  res <- coerce_vertex(rv, principal_graph(cds)[["UMAP"]])
  if (!is.na(res)) res else NULL
}

prep_monocle <- function(obj) {
  obj <- NormalizeData(obj, verbose=FALSE)
  obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
  obj <- ScaleData(obj, verbose=FALSE)
  np <- min(30, max(2, ncol(obj)-1))
  obj <- RunPCA(obj, npcs=np, verbose=FALSE)
  obj <- RunUMAP(obj, dims=1:min(15,np), verbose=FALSE)
  cds <- as.cell_data_set(obj)
  cds <- cluster_cells(cds, verbose=FALSE)
  cds <- learn_graph(cds, verbose=FALSE)
  cds
}

extract_graph <- function(cds) {
  g <- principal_graph(cds)[["UMAP"]]
  gc <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  cv <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(g)||is.null(gc)||is.null(cv)) return(NULL)
  if (is.matrix(cv)) cv <- cv[,1,drop=TRUE]
  g <- compute_graph_weights(g, gc)
  list(graph=g, coords=gc, closest=cv)
}

# --- Load & subset to Basal cells ---
message("Loading data ...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
state_vec  <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")
mp_adj     <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_mp_adj.rds")

common <- intersect(Cells(tmdata_all), intersect(names(state_vec), rownames(mp_adj)))
tmdata_all <- tmdata_all[, common]
state_vec  <- state_vec[common]
mp_adj     <- mp_adj[common, , drop=FALSE]

basal_cells <- names(state_vec)[state_vec == "Basal to intestinal metaplasia"]
message("Basal cells: ", length(basal_cells))

# --- Assign dominant MP & hybrid ---
basal_scores <- mp_adj[basal_cells, intersect(basal_mps, colnames(mp_adj)), drop=FALSE]
best_idx <- max.col(basal_scores, ties.method = "first")
best_val <- apply(basal_scores, 1, max)
sorted   <- t(apply(basal_scores, 1, sort, decreasing = TRUE))
gap      <- sorted[,1] - sorted[,2]

dom_mp <- colnames(basal_scores)[best_idx]
dom_mp[gap < HYBRID_GAP] <- "Hybrid"
names(dom_mp) <- basal_cells

tmdata_all$state_B <- state_vec[Cells(tmdata_all)]
tmdata_basal <- tmdata_all[, basal_cells]
tmdata_basal$dominant_mp <- dom_mp[Cells(tmdata_basal)]

saveRDS(dom_mp, file.path(int_dir, "basal_mp_assignments.rds"))
message("MP assignments — ", paste(names(table(dom_mp)), table(dom_mp), sep=":", collapse=", "))

# --- Proportion table ---
prop_df <- data.frame(cell=basal_cells, mp=dom_mp, sample=as.character(tmdata_basal$orig.ident), stringsAsFactors=FALSE)
overall <- prop_df %>% count(mp) %>% mutate(pct = 100*n/sum(n), sample="Overall")
per_samp <- prop_df %>% count(sample, mp) %>% group_by(sample) %>% mutate(pct=100*n/sum(n)) %>% ungroup()
write_csv(bind_rows(overall, per_samp), file.path(tab_dir, "basal_mp_proportions.csv"))

# Proportion barplot
mp_levels <- c(basal_mps, "Hybrid")
p_prop <- ggplot(overall, aes(x=mp, y=pct, fill=mp)) +
  geom_col(color="black", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)", pct, n)), vjust=-0.2, size=3.5) +
  scale_fill_manual(values=mp_cols, drop=FALSE) +
  scale_x_discrete(labels=function(x) ifelse(x %in% names(mp_labels), mp_labels[x], x)) +
  labs(title="Basal state: dominant MP proportions", x=NULL, y="% of Basal cells") +
  theme_minimal(base_size=13) + theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none") +
  ylim(0, max(overall$pct)*1.15)
ggsave(file.path(fig_dir, "basal_mp_proportions.pdf"), p_prop, width=9, height=6)

# --- Per-sample distance computation ---
meta_df <- data.frame(
  cell=basal_cells, orig.ident=prop_df$sample, mp_group=dom_mp,
  stringsAsFactors=FALSE
) %>% filter(mp_group %in% basal_mps)  # exclude Hybrid from distance calc

samp_counts <- meta_df %>% count(orig.ident, mp_group) %>%
  complete(orig.ident, mp_group=basal_mps, fill=list(n=0))

samp_summary <- samp_counts %>% group_by(orig.ident) %>%
  summarise(total=sum(n), n_groups=sum(n>=OTHER_MIN),
            valid=total>=TOTAL_MIN & n_groups>=MIN_GROUPS, .groups="drop")
valid_samples <- samp_summary %>% filter(valid) %>% pull(orig.ident)
message("Valid samples for distance: ", length(valid_samples))

directed_rows <- list()
geo_med_rows  <- list()
geo_cen_rows  <- list()
umap_cen_rows <- list()

# Use MP14 as default root (Squamoid/basal transition)
ROOT_MP <- "MP14"

pdf(file.path(fig_dir, "basal_mp_pseudotime_samples.pdf"), width=14, height=6, onefile=TRUE)

for (sid in valid_samples) {
  message("Processing: ", sid)
  scells <- meta_df %>% filter(orig.ident==sid) %>% pull(cell)
  # Also include Hybrid cells in trajectory building for continuity
  hybrid_cells <- basal_cells[dom_mp[basal_cells]=="Hybrid" & as.character(tmdata_basal$orig.ident[basal_cells])==sid]
  all_scells <- unique(c(scells, hybrid_cells))
  if (length(all_scells) < TOTAL_MIN) next

  sobj <- tmdata_basal[, all_scells]
  smeta <- data.frame(cell=all_scells, mp_group=dom_mp[all_scells],
                      orig.ident=rep(sid, length(all_scells)), stringsAsFactors=FALSE)

  cds <- tryCatch(prep_monocle(sobj), error=function(e) { message("Skip monocle: ", sid); NULL })
  if (is.null(cds)) next

  gb <- tryCatch(extract_graph(cds), error=function(e) NULL)
  if (is.null(gb)) next

  umap_mat <- reducedDims(cds)$UMAP
  udf <- data.frame(cell=rownames(umap_mat), UMAP_1=umap_mat[,1], UMAP_2=umap_mat[,2], stringsAsFactors=FALSE) %>%
    left_join(smeta, by="cell")

  # Only use non-hybrid cells for representative points
  udf_nh <- udf %>% filter(mp_group %in% basal_mps)
  present_mps <- intersect(basal_mps, unique(udf_nh$mp_group))
  if (length(present_mps) < 2) next

  rep_df <- udf_nh %>% count(mp_group, name="n_cells") %>%
    mutate(medoid_cell = vapply(mp_group, function(m) get_medoid_cell(udf_nh %>% filter(mp_group==m)), character(1))) %>%
    left_join(udf %>% select(cell, UMAP_1, UMAP_2), by=c("medoid_cell"="cell")) %>%
    rowwise() %>%
    mutate(
      cx = mean(udf_nh$UMAP_1[udf_nh$mp_group==mp_group], na.rm=TRUE),
      cy = mean(udf_nh$UMAP_2[udf_nh$mp_group==mp_group], na.rm=TRUE),
      med_v = coerce_vertex(gb$closest[medoid_cell], gb$graph),
      cen_v = nearest_graph_vertex(c(cx,cy), gb$coords)
    ) %>% ungroup()

  # Pairwise root-free distances
  for (ma in present_mps) {
    for (mb in present_mps) {
      ra <- rep_df %>% filter(mp_group==ma)
      rb <- rep_df %>% filter(mp_group==mb)
      if (!nrow(ra)||!nrow(rb)) next

      md <- suppressWarnings(igraph::distances(gb$graph, v=ra$med_v[1], to=rb$med_v[1], weights=igraph::E(gb$graph)$weight)[1,1])
      md <- if (is.finite(md)) as.numeric(md) else NA_real_

      cd <- suppressWarnings(igraph::distances(gb$graph, v=ra$cen_v[1], to=rb$cen_v[1], weights=igraph::E(gb$graph)$weight)[1,1])
      cd <- if (is.finite(cd)) as.numeric(cd) else NA_real_

      ed <- sqrt((ra$cx[1]-rb$cx[1])^2 + (ra$cy[1]-rb$cy[1])^2)

      geo_med_rows[[length(geo_med_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=md, stringsAsFactors=FALSE)
      geo_cen_rows[[length(geo_cen_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=cd, stringsAsFactors=FALSE)
      umap_cen_rows[[length(umap_cen_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=ed, stringsAsFactors=FALSE)
    }
  }

  # Directed pseudotime for all present MPs
  for (cur_root in present_mps) {
    root_cells <- smeta$cell[smeta$mp_group == cur_root]
    if (length(root_cells) >= ROOT_MIN) {
      root_node <- get_root_node(cds, root_cells)
      if (!is.null(root_node)) {
        ord_cds <- tryCatch(order_cells(cds, root_pr_nodes=root_node), error=function(e) NULL)
        if (!is.null(ord_cds)) {
          pt_vals <- pseudotime(ord_cds)
          pt_vals[is.infinite(pt_vals)] <- NA_real_

          # Plot only for ROOT_MP
          if (cur_root == ROOT_MP) {
            mp_count_str <- paste(present_mps, vapply(present_mps, function(m) sum(smeta$mp_group==m), integer(1)), sep=":", collapse=", ")
            p1 <- plot_cells(ord_cds, color_cells_by="pseudotime", show_trajectory_graph=TRUE,
                             label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, cell_size=0.8) +
              scale_color_viridis_c() + labs(title=paste0(sid, " | Pseudotime (MP14 root)")) + theme_minimal(base_size=11)

            colData(ord_cds)$dom_mp <- dom_mp[rownames(colData(ord_cds))]
            p2 <- plot_cells(ord_cds, color_cells_by="dom_mp", show_trajectory_graph=TRUE,
                             label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, cell_size=0.8) +
              scale_color_manual(values=mp_cols, na.value="grey80") +
              labs(title=paste0(sid, " | ", mp_count_str)) + theme_minimal(base_size=11)
            print(p2 + p1)

            saveRDS(pt_vals, file.path(traj_dir, paste0(sid, "_basal_pseudotime.rds")))
          }

          # Per-MP pseudotime stats
          pt_df <- data.frame(cell=names(pt_vals), pseudotime=as.numeric(pt_vals), stringsAsFactors=FALSE) %>%
            left_join(smeta, by="cell") %>% filter(mp_group %in% basal_mps)
          mp_pt <- pt_df %>% group_by(mp_group) %>%
            summarise(med_pt=safe_median(pseudotime), mean_pt=safe_mean(pseudotime), n=sum(is.finite(pseudotime)), .groups="drop")
          root_med <- mp_pt$med_pt[mp_pt$mp_group==cur_root]
          root_mean <- mp_pt$mean_pt[mp_pt$mp_group==cur_root]
          if (length(root_med) && length(root_mean)) {
            for (target in present_mps) {
              tr <- mp_pt %>% filter(mp_group==target)
              if (!nrow(tr)) next
              directed_rows[[length(directed_rows)+1]] <- data.frame(
                sample=sid, root_mp=cur_root, target_mp=target,
                median_distance=tr$med_pt[1]-root_med[1], mean_distance=tr$mean_pt[1]-root_mean[1],
                stringsAsFactors=FALSE)
            }
          }
        }
      }
    }
  }
}

dev.off()

# --- Aggregate summaries ---
dir_df  <- bind_rows(directed_rows)
gm_df   <- bind_rows(geo_med_rows)
gc_df   <- bind_rows(geo_cen_rows)
uc_df   <- bind_rows(umap_cen_rows)

write_csv(dir_df, file.path(tab_dir, "basal_mp_directed_pseudotime.csv"))
write_csv(gm_df,  file.path(tab_dir, "basal_mp_geodesic_medoid.csv"))
write_csv(gc_df,  file.path(tab_dir, "basal_mp_geodesic_centroid.csv"))
write_csv(uc_df,  file.path(tab_dir, "basal_mp_umap_centroid.csv"))

# Symmetrise
sym_fn <- function(df, id_a="mp_a", id_b="mp_b") {
  df %>% mutate(a=pmin(!!sym(id_a),!!sym(id_b)), b=pmax(!!sym(id_a),!!sym(id_b))) %>%
    group_by(sample, a, b) %>% summarise(distance=mean(distance, na.rm=TRUE), .groups="drop") %>%
    rename(mp_a=a, mp_b=b)
}

dir_pw_med <- dir_df %>% filter(root_mp!=target_mp) %>%
  mutate(mp_a=pmin(root_mp,target_mp), mp_b=pmax(root_mp,target_mp), distance=abs(median_distance)) %>%
  group_by(sample, mp_a, mp_b) %>% summarise(distance=mean(distance, na.rm=TRUE), .groups="drop")
dir_pw_mean <- dir_df %>% filter(root_mp!=target_mp) %>%
  mutate(mp_a=pmin(root_mp,target_mp), mp_b=pmax(root_mp,target_mp), distance=abs(mean_distance)) %>%
  group_by(sample, mp_a, mp_b) %>% summarise(distance=mean(distance, na.rm=TRUE), .groups="drop")
gm_sym <- sym_fn(gm_df)
gc_sym <- sym_fn(gc_df)
uc_sym <- sym_fn(uc_df)

agg_fn <- function(df, mname) {
  df %>% group_by(mp_a, mp_b) %>%
    summarise(method=mname, mean_distance=safe_mean(distance), median_distance=safe_median(distance),
              n_samples=sum(is.finite(distance)), .groups="drop")
}

summary_df <- bind_rows(
  agg_fn(dir_pw_med,  "directed_pseudotime_median"),
  agg_fn(dir_pw_mean, "directed_pseudotime_mean"),
  agg_fn(gm_sym, "geodesic_medoid"),
  agg_fn(gc_sym, "geodesic_centroid"),
  agg_fn(uc_sym, "umap_centroid_euclidean")
) %>% arrange(method, mp_a, mp_b)

write_csv(summary_df, file.path(tab_dir, "basal_mp_distance_summary.csv"))

# Build matrices
methods_list <- unique(summary_df$method)
mat_list <- setNames(lapply(methods_list, function(m) make_symmetric_matrix(summary_df, m, basal_mps)), methods_list)
saveRDS(mat_list, file.path(int_dir, "basal_mp_distance_matrices.rds"))

for (mn in names(mat_list)) {
  write.csv(mat_list[[mn]], file.path(tab_dir, paste0("basal_mp_", mn, "_matrix.csv")), quote=FALSE)
}

# --- Heatmap PDF ---
heat_long <- bind_rows(lapply(names(mat_list), function(m) matrix_to_long(mat_list[[m]], m)))

if (nrow(heat_long) > 0) {
  heat_long <- heat_long %>%
    mutate(
      mp_a_lab = ifelse(mp_a %in% names(mp_labels), mp_labels[mp_a], mp_a),
      mp_b_lab = ifelse(mp_b %in% names(mp_labels), mp_labels[mp_b], mp_b),
      mp_a_lab = factor(mp_a_lab, levels=mp_labels[basal_mps]),
      mp_b_lab = factor(mp_b_lab, levels=mp_labels[basal_mps])
    )

  p_heat <- ggplot(heat_long, aes(x=mp_a_lab, y=mp_b_lab, fill=distance)) +
    geom_tile(color="white", linewidth=0.4) +
    geom_text(aes(label=ifelse(is.na(distance), "NA", sprintf("%.2f", distance))), size=3.2) +
    scale_fill_gradient(low="#f7fbff", high="#084594", na.value="grey90") +
    facet_wrap(~method, scales="free") +
    theme_minimal(base_size=12) +
    theme(axis.title=element_blank(), axis.text.x=element_text(angle=40, hjust=1),
          panel.grid=element_blank(), strip.text=element_text(face="bold", size=11))

  ggsave(file.path(fig_dir, "basal_mp_distance_heatmaps.pdf"), p_heat, width=16, height=10)
}

message("Done. Outputs in: ", file.path(getwd(), base_dir))

# --- Nodeplot ---
library(ggrepel)
library(patchwork)

# Functions from hybrid_pairwise_distance_nodeplot.R
add_label_positions <- function(layout_df) {
  center_x <- mean(layout_df$x)
  center_y <- mean(layout_df$y)

  dx <- layout_df$x - center_x
  dy <- layout_df$y - center_y
  norms <- sqrt(dx ^ 2 + dy ^ 2)
  norms[norms == 0] <- 1

  ux <- dx / norms
  uy <- dy / norms

  span <- max(diff(range(layout_df$x)), diff(range(layout_df$y)))
  if (!is.finite(span) || span <= 0) {
    span <- 1
  }
  label_offset <- 0.10 * span

  layout_df$nudge_x <- label_offset * ux
  layout_df$nudge_y <- label_offset * uy
  
  layout_df$label_x <- layout_df$x + layout_df$nudge_x
  layout_df$label_y <- layout_df$y + layout_df$nudge_y

  layout_df$hjust <- ifelse(ux >= 0.15, 0, ifelse(ux <= -0.15, 1, 0.5))
  layout_df$vjust <- ifelse(uy >= 0.15, 0, ifelse(uy <= -0.15, 1, 0.5))
  
  layout_df
}

wrap_state_label <- function(state_name, labels_map, width = 18) {
  vapply(state_name, function(x) {
    if (x %in% names(labels_map)) {
      x <- labels_map[[x]]
    }
    paste(strwrap(x, width = width), collapse = "\n")
  }, character(1))
}

gap_mask <- gap < HYBRID_GAP
hybrid_c <- basal_cells[gap_mask]

if (length(hybrid_c) > 0) {
  top_2 <- t(apply(basal_scores[hybrid_c, , drop=FALSE], 1, function(x) names(sort(x, decreasing=TRUE))[1:2]))
  edf_base <- data.frame(from=top_2[,1], to=top_2[,2], stringsAsFactors=FALSE) %>%
    mutate(from_mp = pmin(from, to), to_mp = pmax(from, to)) %>%
    count(from_mp, to_mp, name="n_hybrid") %>%
    mutate(pct_hybrid = 100 * n_hybrid / length(basal_cells))
  
method_labels <- c(
  "directed_pseudotime_median" = "Directed pseudotime median",
  "directed_pseudotime_mean" = "Directed pseudotime mean",
  "geodesic_medoid" = "Graph geodesic medoid",
  "geodesic_centroid" = "Graph geodesic centroid",
  "umap_centroid_euclidean" = "UMAP centroid Euclidean"
)

nodeplots <- list()
for (m_name in names(mat_list)) {
  cur_mat <- mat_list[[m_name]]
  diag(cur_mat) <- NA
  # Scale matrix for better cmdscale
  cur_mat <- cur_mat / median(cur_mat[upper.tri(cur_mat)], na.rm=TRUE)
  
  mds <- tryCatch(cmdscale(as.dist(cur_mat), k=2), error=function(e) NULL)
  if (is.null(mds)) {
      mds <- matrix(rnorm(10), ncol=2)
      rownames(mds) <- rownames(cur_mat)
  }
  
  ndf <- data.frame(mp=rownames(mds), x=mds[,1], y=mds[,2], stringsAsFactors=FALSE)
  ndf <- ndf %>% left_join(overall, by="mp")
  ndf <- add_label_positions(ndf)
  ndf <- ndf %>%
    mutate(
      label_text = paste0(
        wrap_state_label(mp, labels_map = mp_labels),
        "\n",
        sprintf("%.1f%%", pct)
      )
    )
    
    edf <- edf_base %>%
      left_join(ndf %>% dplyr::select(mp, x, y), by=c("from_mp"="mp")) %>%
      left_join(ndf %>% dplyr::select(mp, x, y), by=c("to_mp"="mp"), suffix = c("", "end"))
      
    p_node <- ggplot() +
      geom_segment(
        data = edf,
        aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct_hybrid),
        color = "grey35",
        alpha = 0.8
      ) +
      geom_point(
        data = ndf,
        aes(x = x, y = y, size = pct, color = mp)
      ) +
      geom_text_repel(
        data = ndf,
        aes(
          x = x,
          y = y,
          label = label_text
        ),
        nudge_x = ndf$nudge_x,
        nudge_y = ndf$nudge_y,
        size = 2.5,
        fontface = "bold",
        lineheight = 0.95,
        segment.color = "grey70",
        segment.size = 0.3,
        min.segment.length = 0,
        box.padding = 0.4,
        point.padding = 0.5,
        force = 3
      ) +
      geom_label(
        data = edf,
        aes(
          x = (x + xend) / 2,
          y = (y + yend) / 2,
          label = paste0(sprintf("%.1f", pct_hybrid), "%")
        ),
        size = 2.0,
        fill = "white",
        label.size = 0,
        linewidth = 0,
        fontface = "bold"
      ) +
      scale_color_manual(values = mp_cols) +
      scale_size(range = c(6, 18), guide = "none") +
      scale_linewidth(range = c(0.5, 5.5), guide = "none") +
      coord_equal(clip = "off") +
      expand_limits(
        x = {
          pad <- 0.10 * max(diff(range(ndf$label_x)), diff(range(ndf$x)))
          c(min(ndf$label_x) - pad, max(ndf$label_x) + pad)
        },
        y = {
          pad <- 0.10 * max(diff(range(ndf$label_y)), diff(range(ndf$y)))
          c(min(ndf$label_y) - pad, max(ndf$label_y) + pad)
        }
      ) +
      theme_void(base_size = 12) +
      labs(
        title = method_labels[[m_name]]
      ) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(t = 8, b = 4)),
        plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
      )
      
    nodeplots[[m_name]] <- p_node
  }
  
  while(length(nodeplots) < 6) {
    nodeplots[[length(nodeplots) + 1]] <- ggplot() + theme_void()
  }
  
  combined_nodeplot <- wrap_plots(nodeplots[1:6], ncol = 3, nrow = 2) +
    plot_annotation(
      title = "Basal MP Nodeplots by Distance Method",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 8)),
        plot.margin = margin(8, 8, 8, 8)
      )
    )
  
  ggsave(file.path(fig_dir, "basal_mp_nodeplot_all_methods.pdf"), combined_nodeplot, width=16, height=10)
  message("Saved combined basal MP nodeplots.")
}

####################
# === SMG TO INTESTINAL METAPLASIA STATE ===
# Parallel analysis for the SMG state MPs
####################

smg_mps <- c("MP8+", "MP8b", "MP16", "MP18b", "MP17")
smg_mp_labels <- c(
  "MP8+"  = "MP8+ Glandular intestinal metaplasia",
  "MP8b"  = "MP8b Metabolic intestinal metaplasia",
  "MP16"  = "MP16 Mucous-secretory glandular epithelium",
  "MP18b" = "MP18b Mucous-secretory differentiation",
  "MP17"  = "MP17 Immune-interactive glandular progenitor"
)
smg_mp_cols <- c(
  "MP8+" = "#FF7F00", "MP8b" = "#FDBF6F", "MP16" = "#FF9896",
  "MP18b" = "#C5B0D5", "MP17" = "#C49C94", Hybrid = "black"
)
SMG_ROOT_MP <- "MP8+"

# --- SMG output dirs ---
smg_base_dir <- "smg_mp_distance"
smg_int_dir  <- file.path(smg_base_dir, "intermediate")
smg_tab_dir  <- file.path(smg_base_dir, "tables")
smg_fig_dir  <- file.path(smg_base_dir, "figures")
smg_traj_dir <- file.path(smg_int_dir, "sample_trajectories")
for (d in c(smg_int_dir, smg_tab_dir, smg_fig_dir, smg_traj_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- SMG cells ---
smg_cells <- names(state_vec)[state_vec == "SMG to intestinal metaplasia"]
message("\n=== SMG to intestinal metaplasia MP distance ===")
message("SMG cells: ", length(smg_cells))

if (length(smg_cells) >= TOTAL_MIN) {
  # --- Assign dominant MP & hybrid ---
  smg_scores <- mp_adj[smg_cells, intersect(smg_mps, colnames(mp_adj)), drop=FALSE]
  smg_best_idx <- max.col(smg_scores, ties.method = "first")
  smg_sorted   <- t(apply(smg_scores, 1, sort, decreasing = TRUE))
  smg_gap      <- smg_sorted[,1] - smg_sorted[,2]

  smg_dom_mp <- colnames(smg_scores)[smg_best_idx]
  smg_dom_mp[smg_gap < HYBRID_GAP] <- "Hybrid"
  names(smg_dom_mp) <- smg_cells

  tmdata_smg <- tmdata_all[, smg_cells]
  tmdata_smg$dominant_mp <- smg_dom_mp[Cells(tmdata_smg)]

  saveRDS(smg_dom_mp, file.path(smg_int_dir, "smg_mp_assignments.rds"))
  message("SMG MP assignments — ", paste(names(table(smg_dom_mp)), table(smg_dom_mp), sep=":", collapse=", "))

  # --- Proportion table ---
  smg_prop_df <- data.frame(cell=smg_cells, mp=smg_dom_mp, sample=as.character(tmdata_smg$orig.ident), stringsAsFactors=FALSE)
  smg_overall <- smg_prop_df %>% count(mp) %>% mutate(pct = 100*n/sum(n), sample="Overall")
  smg_per_samp <- smg_prop_df %>% count(sample, mp) %>% group_by(sample) %>% mutate(pct=100*n/sum(n)) %>% ungroup()
  write_csv(bind_rows(smg_overall, smg_per_samp), file.path(smg_tab_dir, "smg_mp_proportions.csv"))

  # Proportion barplot
  smg_mp_levels <- c(smg_mps, "Hybrid")
  p_smg_prop <- ggplot(smg_overall, aes(x=mp, y=pct, fill=mp)) +
    geom_col(color="black", linewidth=0.3) +
    geom_text(aes(label=sprintf("%.1f%%\n(n=%d)", pct, n)), vjust=-0.2, size=3.5) +
    scale_fill_manual(values=smg_mp_cols, drop=FALSE) +
    scale_x_discrete(labels=function(x) ifelse(x %in% names(smg_mp_labels), smg_mp_labels[x], x)) +
    labs(title="SMG state: dominant MP proportions", x=NULL, y="% of SMG cells") +
    theme_minimal(base_size=13) + theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none") +
    ylim(0, max(smg_overall$pct)*1.15)
  ggsave(file.path(smg_fig_dir, "smg_mp_proportions.pdf"), p_smg_prop, width=9, height=6)

  # --- Per-sample distance computation ---
  smg_meta_df <- data.frame(
    cell=smg_cells, orig.ident=smg_prop_df$sample, mp_group=smg_dom_mp,
    stringsAsFactors=FALSE
  ) %>% filter(mp_group %in% smg_mps)

  smg_samp_counts <- smg_meta_df %>% count(orig.ident, mp_group) %>%
    complete(orig.ident, mp_group=smg_mps, fill=list(n=0))

  smg_samp_summary <- smg_samp_counts %>% group_by(orig.ident) %>%
    summarise(total=sum(n), n_groups=sum(n>=OTHER_MIN),
              valid=total>=TOTAL_MIN & n_groups>=MIN_GROUPS, .groups="drop")
  smg_valid_samples <- smg_samp_summary %>% filter(valid) %>% pull(orig.ident)
  message("Valid SMG samples for distance: ", length(smg_valid_samples))

  smg_directed_rows <- list()
  smg_geo_med_rows  <- list()
  smg_geo_cen_rows  <- list()
  smg_umap_cen_rows <- list()

  pdf(file.path(smg_fig_dir, "smg_mp_pseudotime_samples.pdf"), width=14, height=6, onefile=TRUE)

  for (sid in smg_valid_samples) {
    message("Processing SMG sample: ", sid)
    scells <- smg_meta_df %>% filter(orig.ident==sid) %>% pull(cell)
    hybrid_c_smg <- smg_cells[smg_dom_mp[smg_cells]=="Hybrid" & as.character(tmdata_smg$orig.ident[smg_cells])==sid]
    all_scells <- unique(c(scells, hybrid_c_smg))
    if (length(all_scells) < TOTAL_MIN) next

    sobj <- tmdata_smg[, all_scells]
    smeta <- data.frame(cell=all_scells, mp_group=smg_dom_mp[all_scells],
                        orig.ident=rep(sid, length(all_scells)), stringsAsFactors=FALSE)

    cds <- tryCatch(prep_monocle(sobj), error=function(e) { message("Skip monocle: ", sid); NULL })
    if (is.null(cds)) next

    gb <- tryCatch(extract_graph(cds), error=function(e) NULL)
    if (is.null(gb)) next

    umap_mat <- reducedDims(cds)$UMAP
    udf <- data.frame(cell=rownames(umap_mat), UMAP_1=umap_mat[,1], UMAP_2=umap_mat[,2], stringsAsFactors=FALSE) %>%
      left_join(smeta, by="cell")

    udf_nh <- udf %>% filter(mp_group %in% smg_mps)
    present_mps <- intersect(smg_mps, unique(udf_nh$mp_group))
    if (length(present_mps) < 2) next

    rep_df <- udf_nh %>% count(mp_group, name="n_cells") %>%
      mutate(medoid_cell = vapply(mp_group, function(m) get_medoid_cell(udf_nh %>% filter(mp_group==m)), character(1))) %>%
      left_join(udf %>% select(cell, UMAP_1, UMAP_2), by=c("medoid_cell"="cell")) %>%
      rowwise() %>%
      mutate(
        cx = mean(udf_nh$UMAP_1[udf_nh$mp_group==mp_group], na.rm=TRUE),
        cy = mean(udf_nh$UMAP_2[udf_nh$mp_group==mp_group], na.rm=TRUE),
        med_v = coerce_vertex(gb$closest[medoid_cell], gb$graph),
        cen_v = nearest_graph_vertex(c(cx,cy), gb$coords)
      ) %>% ungroup()

    # Pairwise root-free distances
    for (ma in present_mps) {
      for (mb in present_mps) {
        ra <- rep_df %>% filter(mp_group==ma)
        rb <- rep_df %>% filter(mp_group==mb)
        if (!nrow(ra)||!nrow(rb)) next

        md <- suppressWarnings(igraph::distances(gb$graph, v=ra$med_v[1], to=rb$med_v[1], weights=igraph::E(gb$graph)$weight)[1,1])
        md <- if (is.finite(md)) as.numeric(md) else NA_real_

        cd <- suppressWarnings(igraph::distances(gb$graph, v=ra$cen_v[1], to=rb$cen_v[1], weights=igraph::E(gb$graph)$weight)[1,1])
        cd <- if (is.finite(cd)) as.numeric(cd) else NA_real_

        ed <- sqrt((ra$cx[1]-rb$cx[1])^2 + (ra$cy[1]-rb$cy[1])^2)

        smg_geo_med_rows[[length(smg_geo_med_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=md, stringsAsFactors=FALSE)
        smg_geo_cen_rows[[length(smg_geo_cen_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=cd, stringsAsFactors=FALSE)
        smg_umap_cen_rows[[length(smg_umap_cen_rows)+1]] <- data.frame(sample=sid, mp_a=ma, mp_b=mb, distance=ed, stringsAsFactors=FALSE)
      }
    }

    # Directed pseudotime
    for (cur_root in present_mps) {
      root_cells <- smeta$cell[smeta$mp_group == cur_root]
      if (length(root_cells) >= ROOT_MIN) {
        root_node <- get_root_node(cds, root_cells)
        if (!is.null(root_node)) {
          ord_cds <- tryCatch(order_cells(cds, root_pr_nodes=root_node), error=function(e) NULL)
          if (!is.null(ord_cds)) {
            pt_vals <- pseudotime(ord_cds)
            pt_vals[is.infinite(pt_vals)] <- NA_real_

            if (cur_root == SMG_ROOT_MP) {
              mp_count_str <- paste(present_mps, vapply(present_mps, function(m) sum(smeta$mp_group==m), integer(1)), sep=":", collapse=", ")
              p1 <- plot_cells(ord_cds, color_cells_by="pseudotime", show_trajectory_graph=TRUE,
                               label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, cell_size=0.8) +
                scale_color_viridis_c() + labs(title=paste0(sid, " | Pseudotime (MP8+ root)")) + theme_minimal(base_size=11)

              colData(ord_cds)$dom_mp <- smg_dom_mp[rownames(colData(ord_cds))]
              p2 <- plot_cells(ord_cds, color_cells_by="dom_mp", show_trajectory_graph=TRUE,
                               label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, cell_size=0.8) +
                scale_color_manual(values=smg_mp_cols, na.value="grey80") +
                labs(title=paste0(sid, " | ", mp_count_str)) + theme_minimal(base_size=11)
              print(p2 + p1)

              saveRDS(pt_vals, file.path(smg_traj_dir, paste0(sid, "_smg_pseudotime.rds")))
            }

            pt_df <- data.frame(cell=names(pt_vals), pseudotime=as.numeric(pt_vals), stringsAsFactors=FALSE) %>%
              left_join(smeta, by="cell") %>% filter(mp_group %in% smg_mps)
            mp_pt <- pt_df %>% group_by(mp_group) %>%
              summarise(med_pt=safe_median(pseudotime), mean_pt=safe_mean(pseudotime), n=sum(is.finite(pseudotime)), .groups="drop")
            root_med <- mp_pt$med_pt[mp_pt$mp_group==cur_root]
            root_mean <- mp_pt$mean_pt[mp_pt$mp_group==cur_root]
            if (length(root_med) && length(root_mean)) {
              for (target in present_mps) {
                tr <- mp_pt %>% filter(mp_group==target)
                if (!nrow(tr)) next
                smg_directed_rows[[length(smg_directed_rows)+1]] <- data.frame(
                  sample=sid, root_mp=cur_root, target_mp=target,
                  median_distance=tr$med_pt[1]-root_med[1], mean_distance=tr$mean_pt[1]-root_mean[1],
                  stringsAsFactors=FALSE)
              }
            }
          }
        }
      }
    }
  }

  dev.off()

  # --- SMG Aggregate summaries ---
  smg_dir_df  <- bind_rows(smg_directed_rows)
  smg_gm_df   <- bind_rows(smg_geo_med_rows)
  smg_gc_df   <- bind_rows(smg_geo_cen_rows)
  smg_uc_df   <- bind_rows(smg_umap_cen_rows)

  write_csv(smg_dir_df, file.path(smg_tab_dir, "smg_mp_directed_pseudotime.csv"))
  write_csv(smg_gm_df,  file.path(smg_tab_dir, "smg_mp_geodesic_medoid.csv"))
  write_csv(smg_gc_df,  file.path(smg_tab_dir, "smg_mp_geodesic_centroid.csv"))
  write_csv(smg_uc_df,  file.path(smg_tab_dir, "smg_mp_umap_centroid.csv"))

  smg_dir_pw_med <- smg_dir_df %>% filter(root_mp!=target_mp) %>%
    mutate(mp_a=pmin(root_mp,target_mp), mp_b=pmax(root_mp,target_mp), distance=abs(median_distance)) %>%
    group_by(sample, mp_a, mp_b) %>% summarise(distance=mean(distance, na.rm=TRUE), .groups="drop")
  smg_dir_pw_mean <- smg_dir_df %>% filter(root_mp!=target_mp) %>%
    mutate(mp_a=pmin(root_mp,target_mp), mp_b=pmax(root_mp,target_mp), distance=abs(mean_distance)) %>%
    group_by(sample, mp_a, mp_b) %>% summarise(distance=mean(distance, na.rm=TRUE), .groups="drop")
  smg_gm_sym <- sym_fn(smg_gm_df)
  smg_gc_sym <- sym_fn(smg_gc_df)
  smg_uc_sym <- sym_fn(smg_uc_df)

  smg_summary_df <- bind_rows(
    agg_fn(smg_dir_pw_med,  "directed_pseudotime_median"),
    agg_fn(smg_dir_pw_mean, "directed_pseudotime_mean"),
    agg_fn(smg_gm_sym, "geodesic_medoid"),
    agg_fn(smg_gc_sym, "geodesic_centroid"),
    agg_fn(smg_uc_sym, "umap_centroid_euclidean")
  ) %>% arrange(method, mp_a, mp_b)

  write_csv(smg_summary_df, file.path(smg_tab_dir, "smg_mp_distance_summary.csv"))

  smg_methods_list <- unique(smg_summary_df$method)
  smg_mat_list <- setNames(lapply(smg_methods_list, function(m) make_symmetric_matrix(smg_summary_df, m, smg_mps)), smg_methods_list)
  saveRDS(smg_mat_list, file.path(smg_int_dir, "smg_mp_distance_matrices.rds"))

  for (mn in names(smg_mat_list)) {
    write.csv(smg_mat_list[[mn]], file.path(smg_tab_dir, paste0("smg_mp_", mn, "_matrix.csv")), quote=FALSE)
  }

  # --- SMG Heatmap ---
  smg_heat_long <- bind_rows(lapply(names(smg_mat_list), function(m) matrix_to_long(smg_mat_list[[m]], m)))

  if (nrow(smg_heat_long) > 0) {
    smg_heat_long <- smg_heat_long %>%
      mutate(
        mp_a_lab = ifelse(mp_a %in% names(smg_mp_labels), smg_mp_labels[mp_a], mp_a),
        mp_b_lab = ifelse(mp_b %in% names(smg_mp_labels), smg_mp_labels[mp_b], mp_b),
        mp_a_lab = factor(mp_a_lab, levels=smg_mp_labels[smg_mps]),
        mp_b_lab = factor(mp_b_lab, levels=smg_mp_labels[smg_mps])
      )

    p_smg_heat <- ggplot(smg_heat_long, aes(x=mp_a_lab, y=mp_b_lab, fill=distance)) +
      geom_tile(color="white", linewidth=0.4) +
      geom_text(aes(label=ifelse(is.na(distance), "NA", sprintf("%.2f", distance))), size=3.2) +
      scale_fill_gradient(low="#f7fbff", high="#084594", na.value="grey90") +
      facet_wrap(~method, scales="free") +
      theme_minimal(base_size=12) +
      theme(axis.title=element_blank(), axis.text.x=element_text(angle=40, hjust=1),
            panel.grid=element_blank(), strip.text=element_text(face="bold", size=11))

    ggsave(file.path(smg_fig_dir, "smg_mp_distance_heatmaps.pdf"), p_smg_heat, width=16, height=10)
  }

  # --- SMG Nodeplot ---
  smg_gap_vals <- smg_sorted[,1] - smg_sorted[,2]
  smg_hybrid_c <- smg_cells[smg_gap_vals < HYBRID_GAP]

  if (length(smg_hybrid_c) > 0) {
    smg_top_2 <- t(apply(smg_scores[smg_hybrid_c, , drop=FALSE], 1, function(x) names(sort(x, decreasing=TRUE))[1:2]))
    smg_edf_base <- data.frame(from=smg_top_2[,1], to=smg_top_2[,2], stringsAsFactors=FALSE) %>%
      mutate(from_mp = pmin(from, to), to_mp = pmax(from, to)) %>%
      count(from_mp, to_mp, name="n_hybrid") %>%
      mutate(pct_hybrid = 100 * n_hybrid / length(smg_cells))

    smg_nodeplots <- list()
    for (m_name in names(smg_mat_list)) {
      cur_mat <- smg_mat_list[[m_name]]
      diag(cur_mat) <- NA
      cur_mat <- cur_mat / median(cur_mat[upper.tri(cur_mat)], na.rm=TRUE)

      mds <- tryCatch(cmdscale(as.dist(cur_mat), k=2), error=function(e) NULL)
      if (is.null(mds)) {
        mds <- matrix(rnorm(2*length(smg_mps)), ncol=2)
        rownames(mds) <- rownames(cur_mat)
      }

      ndf <- data.frame(mp=rownames(mds), x=mds[,1], y=mds[,2], stringsAsFactors=FALSE)
      ndf <- ndf %>% left_join(smg_overall, by="mp")
      ndf <- add_label_positions(ndf)
      ndf <- ndf %>%
        mutate(
          label_text = paste0(
            wrap_state_label(mp, labels_map = smg_mp_labels),
            "\n",
            sprintf("%.1f%%", pct)
          )
        )

      edf <- smg_edf_base %>%
        left_join(ndf %>% dplyr::select(mp, x, y), by=c("from_mp"="mp")) %>%
        left_join(ndf %>% dplyr::select(mp, x, y), by=c("to_mp"="mp"), suffix = c("", "end"))

      p_node <- ggplot() +
        geom_segment(
          data = edf,
          aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct_hybrid),
          color = "grey35", alpha = 0.8
        ) +
        geom_point(
          data = ndf,
          aes(x = x, y = y, size = pct, color = mp)
        ) +
        geom_text_repel(
          data = ndf,
          aes(x = x, y = y, label = label_text),
          nudge_x = ndf$nudge_x, nudge_y = ndf$nudge_y,
          size = 2.5, fontface = "bold", lineheight = 0.95,
          segment.color = "grey70", segment.size = 0.3,
          min.segment.length = 0, box.padding = 0.4,
          point.padding = 0.5, force = 3
        ) +
        geom_label(
          data = edf,
          aes(x = (x + xend) / 2, y = (y + yend) / 2,
              label = paste0(sprintf("%.1f", pct_hybrid), "%")),
          size = 2.0, fill = "white", label.size = 0,
          linewidth = 0, fontface = "bold"
        ) +
        scale_color_manual(values = smg_mp_cols) +
        scale_size(range = c(6, 18), guide = "none") +
        scale_linewidth(range = c(0.5, 5.5), guide = "none") +
        coord_equal(clip = "off") +
        expand_limits(
          x = {
            pad <- 0.10 * max(diff(range(ndf$label_x)), diff(range(ndf$x)))
            c(min(ndf$label_x) - pad, max(ndf$label_x) + pad)
          },
          y = {
            pad <- 0.10 * max(diff(range(ndf$label_y)), diff(range(ndf$y)))
            c(min(ndf$label_y) - pad, max(ndf$label_y) + pad)
          }
        ) +
        theme_void(base_size = 12) +
        labs(title = method_labels[[m_name]]) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(t = 8, b = 4)),
          plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
        )

      smg_nodeplots[[m_name]] <- p_node
    }

    while(length(smg_nodeplots) < 6) {
      smg_nodeplots[[length(smg_nodeplots) + 1]] <- ggplot() + theme_void()
    }

    smg_combined_nodeplot <- wrap_plots(smg_nodeplots[1:6], ncol = 3, nrow = 2) +
      plot_annotation(
        title = "SMG MP Nodeplots by Distance Method",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 8)),
          plot.margin = margin(8, 8, 8, 8)
        )
      )

    ggsave(file.path(smg_fig_dir, "smg_mp_nodeplot_all_methods.pdf"), smg_combined_nodeplot, width=16, height=10)
    message("Saved combined SMG MP nodeplots.")
  } else {
    message("No hybrid SMG cells detected — skipping SMG nodeplots.")
  }

  message("SMG MP distance analysis complete. Outputs in: ", file.path(getwd(), smg_base_dir))
} else {
  message("Too few SMG cells (", length(smg_cells), ") — skipping SMG distance analysis.")
}
