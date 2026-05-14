####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/mp_3ca_ucell_scoring.R
#   Methodology: analysis/methodology/metaprograms/metaprogram_scoring_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################

####################
# Auto_3ca_mp_ucell_scoring.R
#
# Regenerate pan-cancer 3CA MP UCell scores for the epithelial scATLAS object.
#
# Input:
#   ref_outs/EAC_Ref_epi.rds
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#
# Output:
#   ref_outs/UCell_3CA_MPs.rds
#   ref_outs/Auto_3ca_mp_ucell_summary.csv
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
})

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs")

args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 8L

tmdata_all <- readRDS("EAC_Ref_epi.rds")
mp_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
if (!file.exists(mp_path)) stop("Missing 3CA MP gene list: ", mp_path)

mp_df <- read.csv(mp_path, check.names = FALSE)
mp_list <- as.list(mp_df)
mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
mp_list <- mp_list[lengths(mp_list) > 0]
names(mp_list) <- make.names(sub("^MP", "3CA_mp_", names(mp_list)))

tmdata_all <- AddModuleScore_UCell(tmdata_all, features = mp_list, ncores = ncores, name = "")
score_cols <- grep("^X3CA_mp|^3CA_mp", colnames(tmdata_all@meta.data), value = TRUE)
if (length(score_cols) == 0) stop("No 3CA UCell score columns were generated.")

ucell_scores <- tmdata_all@meta.data[, score_cols, drop = FALSE]
saveRDS(ucell_scores, file = "UCell_3CA_MPs.rds")

summary_df <- data.frame(
  n_cells = nrow(ucell_scores),
  n_mps = ncol(ucell_scores),
  min_score = min(as.matrix(ucell_scores), na.rm = TRUE),
  max_score = max(as.matrix(ucell_scores), na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(summary_df, "Auto_3ca_mp_ucell_summary.csv", row.names = FALSE)

message("Saved UCell_3CA_MPs.rds with ", nrow(ucell_scores), " cells and ", ncol(ucell_scores), " MPs.")
