####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/final_state_overall_proportions.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_overall_state_proportions.R
# Overall proportion barplot for noreg Approach B cell states.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   ref_outs/Auto_topmp_v2_noreg_states_B.rds
#
# Output:
#   ref_outs/Auto_overall_state_proportions.pdf
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

# Load data
message("Loading data ...")
tmdata_all <- readRDS("EAC_Ref_epi.rds")
final_states_path <- "Auto_final_states.rds"
if (file.exists(final_states_path)) {
  state_B <- readRDS(final_states_path)
} else {
  state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
}

# Common cells
common_cells <- intersect(names(state_B), Cells(tmdata_all))
state_B <- state_B[common_cells]

# State order and colors
group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

# Identify any extra states (e.g. 3CA relabeled)
extra_states <- setdiff(unique(as.character(state_B)), names(group_cols))
if (length(extra_states) > 0) {
  new_cols <- setNames(scales::hue_pal()(length(extra_states)), extra_states)
  group_cols <- c(group_cols, new_cols)
}

state_order <- names(group_cols)

# Calculate proportions
prop_df <- data.frame(
  state = factor(as.character(state_B), levels = state_order),
  stringsAsFactors = FALSE
)

overall <- prop_df %>%
  count(state, .drop = FALSE) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  mutate(label = sprintf("%.1f%%", pct))

# Plotting overall barplot (single bar)
message("Generating overall barplot ...")
p <- ggplot(overall, aes(x = "Overall", y = pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 6, 
            fontface = "bold",
            color = ifelse(overall$state %in% c("Hybrid"), "white", "black")) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(title = "scATLAS State Proportions",
       subtitle = paste0("Total cells: ", sum(overall$n)),
       x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        legend.title = element_text(size = 16),  # title size
        legend.text  = element_text(size = 14)) +   # label size
  coord_cartesian(expand = FALSE)

# Save output
pdf("Auto_overall_state_proportions.pdf", width = 7, height = 9)
print(p)
dev.off()

message("Finished. Plot saved to ref_outs/Auto_overall_state_proportions.pdf")
