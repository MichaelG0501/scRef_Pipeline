####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/metaprograms/legacy_nmf_rank_selection_diagnostics.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

library(cluster)
library(ggplot2)
library(patchwork)

results_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results"
k_vals <- 8:30
avg_sil_widths <- numeric(length(k_vals))
wss_vals <- numeric(length(k_vals))

for (i in seq_along(k_vals)) {
	  k <- k_vals[i]
  rds_path <- file.path(results_dir, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))
    
    if (file.exists(rds_path)) {
	        mp_res <- readRDS(rds_path)
      
      # Pull distance matrix directly used during clustering 
      dist_mat <- as.dist(1 - mp_res$programs.similarity)
          cluster_assignments <- cutree(mp_res$programs.tree, k = k)
          
          # 1. Silhouette score
          sil <- silhouette(cluster_assignments, dist = dist_mat)
	      avg_sil_widths[i] <- summary(sil)$avg.width
	      
	      # 2. Within-Cluster Sum of Squares (WSS) based on Cosine Distance
	      wss_k <- 0
	          dist_m <- as.matrix(dist_mat)
	          for (clust_id in unique(cluster_assignments)) {
			        idx <- which(cluster_assignments == clust_id)
		        if (length(idx) > 1) {
				        # Add sum of squared pairwise distances within cluster divided by 2|C|
				        cluster_dist <- dist_m[idx, idx]
			        wss_k <- wss_k + sum(cluster_dist^2) / (2 * length(idx))
				      }
			    }
		      wss_vals[i] <- wss_k
		    } else {
			        avg_sil_widths[i] <- NA
		        wss_vals[i] <- NA
			  }
}

df_metrics <- data.frame(nMP = k_vals, Silhouette = avg_sil_widths, WSS = wss_vals)

# Print metrics table for easy inspection
print(df_metrics)

####################
# Find optimal nMP using inflection point (maximum curvature / knee detection)
# Uses the "kneedle" approach: for each metric curve, find the point of maximum
# distance from the line connecting the first and last data points.
# For silhouette (increasing): first major slowdown in gains
# For WSS (decreasing): classical elbow point
####################

find_knee <- function(x, y) {
  # Normalize to [0,1] range
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  # Line from first to last point
  x1 <- x_norm[1]; y1 <- y_norm[1]
  x2 <- x_norm[length(x_norm)]; y2 <- y_norm[length(y_norm)]
  # Distance from each point to this line
  dists <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
           sqrt((y2 - y1)^2 + (x2 - x1)^2)
  # Return the x value with maximum distance (the knee/elbow)
  return(x[which.max(dists)])
}

# Silhouette knee (first slowdown in increasing curve)
sil_knee <- find_knee(df_metrics$nMP, df_metrics$Silhouette)
# WSS elbow (classic elbow in decreasing curve)
wss_knee <- find_knee(df_metrics$nMP, df_metrics$WSS)

optimal_nMP <- sil_knee  # primary criterion
message(paste0("Silhouette inflection point: nMP = ", sil_knee))
message(paste0("WSS elbow point: nMP = ", wss_knee))
message(paste0("Selected optimal nMP: ", optimal_nMP))

p1 <- ggplot(df_metrics, aes(x = nMP, y = Silhouette)) +
	  geom_line(color = "steelblue", linewidth = 1) +
	    geom_point(color = "steelblue", size = 3) +
	    geom_vline(xintercept = sil_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
	    annotate("text", x = sil_knee + 0.5, y = max(df_metrics$Silhouette, na.rm = TRUE),
	             label = paste0("Inflection: ", sil_knee), hjust = 0, color = "red", size = 3.5) +
	      theme_minimal() +
	        scale_x_continuous(breaks = k_vals) +
		  labs(title = "Silhouette Analysis",
		              x = "Number of MetaPrograms (nMP)",
			             y = "Average Silhouette Width") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(df_metrics, aes(x = nMP, y = WSS)) +
	  geom_line(color = "darkred", linewidth = 1) +
	    geom_point(color = "darkred", size = 3) +
	    geom_vline(xintercept = wss_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
	    annotate("text", x = wss_knee + 0.5, y = max(df_metrics$WSS, na.rm = TRUE) * 0.95,
	             label = paste0("Elbow: ", wss_knee), hjust = 0, color = "red", size = 3.5) +
	      theme_minimal() +
	        scale_x_continuous(breaks = k_vals) +
		  labs(title = "Elbow Method (WSS)",
		              x = "Number of MetaPrograms (nMP)",
			             y = "Total Within Sum of Squares") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Use patchwork to display them side-by-side
p1 + p2
