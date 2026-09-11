####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_opnme_prioritisation.R
#   Description: In-silico prioritisation of 99 opnMe small molecules against
#     malignant cell states in both scAtlas and PDO (untreated only) datasets.
#     Two-layer data-driven approach per dataset: (1) per-state target gene
#     expression validation, (2) PROGENy pathway activity scoring per state.
#     Mechanism-based priors from the supervisor workbook are used as the
#     knowledge layer. Outputs ranked compound tables, publication figures,
#     and formatted Excel workbooks for both datasets.
#   Methodology: not required — deterministic scoring pipeline
#   Inputs:
#     - ref_outs/opnMe_prioritisation/opnMe_99_drug_OAC_state_predictions.xlsx
#     - ref_outs/EAC_Ref_epi.rds (scAtlas epithelial Seurat)
#     - ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#     - /rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/PDOs_merged.rds
#     - /rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#   Outputs:
#     tables/:
#       - opnme_scatlas_target_expression_by_state.csv
#       - opnme_pdo_target_expression_by_state.csv
#       - opnme_scatlas_pathway_activity_by_state.csv
#       - opnme_pdo_pathway_activity_by_state.csv
#       - opnme_scatlas_final_ranked.csv
#       - opnme_pdo_final_ranked.csv
#       - opnme_combined_ranking.csv
#     figures/:
#       - opnme_scatlas_fc_vs_expression.pdf/png
#       - opnme_pdo_fc_vs_expression.pdf/png
#       - opnme_scatlas_prioritisation_heatmap.pdf/png
#       - opnme_pdo_prioritisation_heatmap.pdf/png
#       - opnme_scatlas_pathway_state_heatmap.pdf/png
#       - opnme_pdo_pathway_state_heatmap.pdf/png
#       - opnme_rank_scatlas_vs_pdo.pdf/png
#     excel/:
#       - opnme_scatlas_prioritisation.xlsx
#       - opnme_pdo_prioritisation.xlsx
#   Cache: supports SCREF_FORCE_REBUILD=TRUE / SCREF_REPLOT_ONLY=TRUE
#   PBS: Auto_opnme_prioritisation.sh (dmtcp, 1 node, 8 cpus, 96 gb, 3h)
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(openxlsx)
})

####################
# setup
####################

scref_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
pdo_dir   <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(file.path(scref_dir, "ref_outs"))

out_dir   <- "opnMe_prioritisation"
table_dir <- file.path(out_dir, "tables")
fig_dir   <- file.path(out_dir, "figures")
cache_dir <- file.path(out_dir, "cache")
excel_dir <- file.path(out_dir, "excel")
for (d in c(out_dir, table_dir, fig_dir, cache_dir, excel_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

force_rebuild <- identical(Sys.getenv("SCREF_FORCE_REBUILD", "FALSE"), "TRUE")
replot_only   <- identical(Sys.getenv("SCREF_REPLOT_ONLY", "FALSE"), "TRUE")

set.seed(1471)

# scAtlas states
scatlas_state_order <- c(
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)
scatlas_state_colors <- c(
  "Classic proliferation"      = "#E41A1C",
  "Squamous-to-intestinal"     = "#377EB8",
  "Glandular-to-intestinal"    = "#4DAF4A",
  "Stress-adaptive"            = "#984EA3",
  "Cancer-cell immune mimicry" = "#FF7F00"
)

# PDO states
pdo_state_order <- c(
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Glandular differentiation",
  "Stress-adaptive",
  "ECM-remodelling",
  "Motile-cilia differentiation"
)
pdo_state_colors <- c(
  "Classic proliferation"        = "#E41A1C",
  "Columnar-to-intestinal"       = "#4DAF4A",
  "Glandular differentiation"    = "#FF7F00",
  "Stress-adaptive"              = "#984EA3",
  "ECM-remodelling"              = "#A65628",
  "Motile-cilia differentiation" = "#F781BF"
)

# Focus states in both datasets
primary_states <- c("Classic proliferation", "Stress-adaptive")

####################
# helpers
####################

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

####################
# Layer 1: Parse supervisor annotation workbook
####################

message("=== Layer 1: Loading supervisor annotation workbook ===")

workbook_path <- file.path(out_dir, "opnMe_99_drug_OAC_state_predictions.xlsx")
if (!file.exists(workbook_path)) {
  stop("Supervisor annotation workbook not found: ", workbook_path)
}

annot <- read_excel(workbook_path, sheet = "99-drug annotations")

state1_col <- grep("State 1", colnames(annot), value = TRUE)
state2_col <- grep("State 2", colnames(annot), value = TRUE)
state3_col <- grep("State 3", colnames(annot), value = TRUE)
state4_col <- grep("State 4", colnames(annot), value = TRUE)

drug_table <- annot %>%
  transmute(
    compound       = Compound,
    plate_well     = `Plate well`,
    target_class   = `Library target class`,
    primary_target = `Primary target(s)`,
    mechanism      = Mechanism,
    pathway        = `Pathway / biological process`,
    neg_control    = `Matched negative control`,
    oac_applicability = `OAC applicability`,
    oac_prior_relevance = `OAC prior relevance`,
    state1_score   = .data[[state1_col]],
    state2_score   = .data[[state2_col]],
    state3_score   = .data[[state3_col]],
    state4_score   = .data[[state4_col]],
    priority_tier  = `Priority tier`,
    composite_prior = `Composite prior score`,
    predicted_effect = `Dominant predicted state/effect`,
    expected_response = .data[[grep("Expected", colnames(annot), value = TRUE)[1]]],
    confidence     = `Prediction confidence`,
    evidence_note  = `Evidence/source note`,
    source_url     = `Source URL`
  )

message("  Loaded ", nrow(drug_table), " compounds")
message("  Priority tier distribution:")
print(table(drug_table$priority_tier))

####################
# Parse target gene symbols from primary_target column
####################

message("\n=== Parsing target gene symbols ===")

target_gene_map <- list(
  "AURKB (Aurora B)"                  = c("AURKB"),
  "GPR142"                            = c("GPR142"),
  "BDKRB1 (bradykinin B1 receptor)"   = c("BDKRB1"),
  "ADRB2 (\u03b22-adrenergic receptor)"   = c("ADRB2"),
  "PTK2/FAK"                          = c("PTK2"),
  "PHGDH"                             = c("PHGDH"),
  "CCR10"                             = c("CCR10"),
  "SLC9A1 / NHE1"                     = c("SLC9A1"),
  "CTSS (cathepsin S)"                = c("CTSS"),
  "KRAS/HRAS/NRAS switch I/II pocket" = c("KRAS", "HRAS", "NRAS"),
  "KRAS G12C"                         = c("KRAS"),
  "KRAS (pan-KRAS)"                   = c("KRAS"),
  "PLK1"                              = c("PLK1"),
  "SMARCA2/SMARCA4/PBRM1"             = c("SMARCA2", "SMARCA4", "PBRM1"),
  "SMARCA2 >> SMARCA4/PBRM1"          = c("SMARCA2", "SMARCA4", "PBRM1"),
  "FASN"                              = c("FASN"),
  "SOS1::KRAS"                        = c("SOS1", "KRAS"),
  "SOS1"                              = c("SOS1"),
  "IRAK4"                             = c("IRAK4"),
  "EGFR"                              = c("EGFR"),
  "RPS6KA1/2/3 (RSK1/2/3)"            = c("RPS6KA1", "RPS6KA2", "RPS6KA3"),
  "IKBKB / IKK\u03b2"                      = c("IKBKB"),
  "CDK8/cyclin C"                     = c("CDK8", "CCNC"),
  "ENPP2 / autotaxin"                 = c("ENPP2"),
  "BRD9/BRD7"                         = c("BRD9", "BRD7"),
  "BRD9 > BRD7"                       = c("BRD9", "BRD7"),
  "BRD9"                              = c("BRD9"),
  "BRD4 > BRD2/BRD3"                  = c("BRD4", "BRD2", "BRD3"),
  "BRD2/BRD3/BRD4 BET bromodomains"   = c("BRD2", "BRD3", "BRD4"),
  "BRD4"                              = c("BRD4"),
  "TGFBR1 / ALK5 (PDGFRA secondary at higher dose)" = c("TGFBR1", "PDGFRA"),
  "ABCC4 / MRP4"                      = c("ABCC4"),
  "MAPK14/p38\u03b1 (p38 family)"         = c("MAPK14"),
  "MLKL"                              = c("MLKL"),
  "ADRB3 (\u03b23-adrenergic receptor)"   = c("ADRB3"),
  "GSK3A/GSK3B"                       = c("GSK3A", "GSK3B"),
  "SYK"                               = c("SYK"),
  "BCL6"                              = c("BCL6"),
  "GCH1 (GTP cyclohydrolase I)"       = c("GCH1"),
  "NR3C1 / glucocorticoid receptor"   = c("NR3C1"),
  "TERT / telomerase"                 = c("TERT"),
  "BPTF bromodomain"                  = c("BPTF"),
  "NSD3/WHSC1L1 PWWP1"                = c("NSD3"),
  "VNN1/VNN2 (vanin-1/2)"             = c("VNN1", "VNN2"),
  "MMP13"                             = c("MMP13"),
  "ALOX5AP / FLAP"                    = c("ALOX5AP"),
  "FFAR1 / GPR40"                     = c("FFAR1"),
  "EPHX2 / soluble epoxide hydrolase" = c("EPHX2"),
  "OLR1 / LOX-1"                      = c("OLR1"),
  "PDE9A"                             = c("PDE9A"),
  "CCR1"                              = c("CCR1"),
  "LTB4R1 / BLT1"                     = c("LTB4R", "LTB4R1"),
  "NPY1R"                             = c("NPY1R"),
  "NPY2R"                             = c("NPY2R"),
  "HCRTR1 / OX1R"                     = c("HCRTR1"),
  "GRM1 / mGluR1"                     = c("GRM1"),
  "GRM2 / mGluR2"                     = c("GRM2"),
  "HTR2C (5-HT2C receptor)"           = c("HTR2C"),
  "GPR88"                             = c("GPR88"),
  "SCTR (secretin receptor)"          = c("SCTR"),
  "PTAFR (PAF receptor)"              = c("PTAFR"),
  "AMPK complexes"                    = c("PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2"),
  "CTSC / cathepsin C"                = c("CTSC"),
  "CTSC / cathepsin C assay substrate" = c("CTSC"),
  "ELANE / neutrophil elastase"       = c("ELANE"),
  "CMA1 (chymase)"                    = c("CMA1"),
  "PTPN5 / STEP"                      = c("PTPN5"),
  "SCN1A / Nav1.1"                    = c("SCN1A"),
  "SCN2A / Nav1.2"                    = c("SCN2A"),
  "SCNN1A/B/G (ENaC)"                 = c("SCNN1A", "SCNN1B", "SCNN1G"),
  "SLC13A5 / NaCT"                    = c("SLC13A5"),
  "CETP"                              = c("CETP"),
  "HSD17B13"                          = c("HSD17B13"),
  "KHK-A/KHK-C (ketohexokinase)"      = c("KHK"),
  "ITGAL/CD11a\u2013ICAM1 (LFA-1)"        = c("ITGAL", "ICAM1"),
  "GABRA5-containing GABA-A receptor" = c("GABRA5"),
  "NMDA receptor (GRIN complex)"      = c("GRIN1", "GRIN2A", "GRIN2B"),
  "SLO-1 / BK-type Ca2+-activated K+ channel (invertebrate tool)" = c("KCNMA1"),
  # Non-human / viral / parasitic targets
  "HCV NS3 protease"                  = character(0),
  "HCV NS3\u2013NS4A protease"            = character(0),
  "HCV NS3/4A protease (faldaprevir)" = character(0),
  "HCV NS5B polymerase"              = character(0),
  "HIV-1 capsid"                     = character(0),
  "HIV-1 integrase (allosteric/non-catalytic site)" = character(0),
  "HIV-1 reverse transcriptase"      = character(0),
  "Human cytomegalovirus DNA polymerase" = character(0),
  "Mycobacterial ClpC1 + BRDT ligand (BacPROTAC)" = c("BRDT"),
  "Plasmodium falciparum DPAP1"      = character(0),
  "TPP riboswitch RNA aptamer"       = character(0)
)

drug_table$target_genes <- vapply(drug_table$primary_target, function(tgt) {
  genes <- target_gene_map[[tgt]]
  if (is.null(genes)) return(NA_character_)
  if (length(genes) == 0) return("non_human")
  paste(genes, collapse = ";")
}, character(1))

unmapped <- drug_table %>% filter(is.na(target_genes))
if (nrow(unmapped) > 0) {
  message("WARNING: ", nrow(unmapped), " compounds have unmapped targets:")
  for (i in seq_len(nrow(unmapped))) {
    message("  ", unmapped$compound[i], " -> ", unmapped$primary_target[i])
  }
}

drug_genes <- drug_table %>%
  filter(!is.na(target_genes), target_genes != "non_human") %>%
  separate_rows(target_genes, sep = ";") %>%
  rename(gene = target_genes)

message("  ", length(unique(drug_genes$gene)), " unique target gene symbols across ",
        length(unique(drug_genes$compound)), " compounds with human targets")

####################
# Drug pathway to PROGENy mapping
####################

drug_pathway_to_progeny <- list(
  "EGFR"        = c("EGFR"),
  "RAS"         = c("MAPK"),
  "RAF"         = c("MAPK"),
  "MEK"         = c("MAPK"),
  "ERK"         = c("MAPK"),
  "MAPK"        = c("MAPK"),
  "PI3K"        = c("PI3K"),
  "mTOR"        = c("PI3K"),
  "FAK"         = c("VEGF"),
  "YAP"         = c("Hypoxia"),
  "EMT"         = c("TGFb"),
  "WNT"         = c("WNT"),
  "NF-kB"       = c("NFkB"),
  "NF-\u03baB"       = c("NFkB"),
  "IRAK"        = c("NFkB"),
  "TGF"         = c("TGFb"),
  "SMAD"        = c("TGFb"),
  "p38"         = c("MAPK"),
  "JAK"         = c("JAK-STAT"),
  "STAT"        = c("JAK-STAT"),
  "TNF"         = c("TNFa"),
  "necroptosis" = c("TNFa"),
  "Hypoxia"     = c("Hypoxia"),
  "VEGF"        = c("VEGF"),
  "p53"         = c("p53"),
  "estrogen"    = c("Estrogen"),
  "androgen"    = c("Androgen"),
  "trail"       = c("Trail")
)

get_progeny_for_pathway <- function(pw) {
  if (is.na(pw)) return(NULL)
  matched <- character(0)
  for (kw in names(drug_pathway_to_progeny)) {
    if (grepl(kw, pw, ignore.case = TRUE)) {
      matched <- c(matched, drug_pathway_to_progeny[[kw]])
    }
  }
  if (length(matched) == 0) return(NULL)
  unique(matched)
}

####################
# Generic functions for computing expression and pathway across any dataset
####################

compute_target_expression <- function(seurat_obj, state_vec, state_order_vec, dataset_label) {
  message("  Computing per-state target expression for ", dataset_label, "...")

  norm_data <- tryCatch(
    GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  )

  all_target_genes <- unique(drug_genes$gene)
  available_genes <- intersect(all_target_genes, rownames(norm_data))
  missing_genes <- setdiff(all_target_genes, rownames(norm_data))
  message("  Target genes found: ", length(available_genes), "/", length(all_target_genes))
  if (length(missing_genes) > 0) {
    message("  Missing: ", paste(missing_genes, collapse = ", "))
  }

  cell_state <- state_vec

  target_expr <- bind_rows(lapply(available_genes, function(g) {
    expr_vec <- as.numeric(norm_data[g, ])
    bind_rows(lapply(state_order_vec, function(s) {
      in_state <- which(cell_state == s)
      out_state <- which(cell_state != s)
      tibble(
        gene = g,
        state = s,
        mean_expr_in_state = mean(expr_vec[in_state]),
        mean_expr_out_state = mean(expr_vec[out_state]),
        pct_expressing_in = mean(expr_vec[in_state] > 0) * 100,
        pct_expressing_out = mean(expr_vec[out_state] > 0) * 100,
        log2fc = log2((mean(expr_vec[in_state]) + 0.01) /
                       (mean(expr_vec[out_state]) + 0.01)),
        n_cells_in = length(in_state),
        n_cells_out = length(out_state)
      )
    }))
  }))

  overall_expr <- bind_rows(lapply(available_genes, function(g) {
    expr_vec <- as.numeric(norm_data[g, ])
    tibble(
      gene = g,
      overall_mean_expr = mean(expr_vec),
      overall_pct_expressing = mean(expr_vec > 0) * 100
    )
  }))

  left_join(target_expr, overall_expr, by = "gene")
}

compute_progeny_scores <- function(seurat_obj, state_vec, state_order_vec, dataset_label) {
  if (!requireNamespace("progeny", quietly = TRUE)) {
    message("  progeny not available; skipping pathway scoring for ", dataset_label)
    return(NULL)
  }
  message("  Computing PROGENy pathway scores for ", dataset_label, "...")

  seurat_obj <- progeny::progeny(seurat_obj, scale = FALSE, organism = "Human",
                                  top = 500, perm = 1, return_assay = TRUE)
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = "progeny")

  progeny_mat <- tryCatch(
    GetAssayData(seurat_obj, assay = "progeny", layer = "scale.data"),
    error = function(e) GetAssayData(seurat_obj, assay = "progeny", slot = "scale.data")
  )

  cell_state <- state_vec

  pathway_scores <- bind_rows(lapply(rownames(progeny_mat), function(pw) {
    pw_vec <- as.numeric(progeny_mat[pw, ])
    bind_rows(lapply(state_order_vec, function(s) {
      in_state <- which(cell_state == s)
      out_state <- which(cell_state != s)
      tibble(
        pathway = pw,
        state = s,
        mean_score_in = mean(pw_vec[in_state]),
        mean_score_out = mean(pw_vec[out_state]),
        score_diff = mean(pw_vec[in_state]) - mean(pw_vec[out_state]),
        n_cells = length(in_state)
      )
    }))
  }))

  pathway_scores
}

####################
# Composite score function (no MOA; reweighted)
# Weights: mechanism prior (45%) + target expression (35%) + pathway activity (20%)
####################

compute_composite <- function(state_score, target_pct, target_fc, pw_score, conf) {
  prior_score <- state_score / 3

  expr_score <- ifelse(is.na(target_pct), 0, pmin(target_pct / 50, 1))
  fc_val <- ifelse(is.na(target_fc), 0, target_fc)
  fc_score <- pmin(pmax(fc_val / 1, 0), 1)
  target_score <- expr_score * 0.6 + fc_score * 0.4

  pathway_score <- ifelse(is.na(pw_score), 0.5, pmin(pmax((pw_score + 2) / 4, 0), 1))

  conf_weight <- ifelse(conf == "High", 1.0, ifelse(conf == "Medium", 0.7, 0.4))

  (prior_score * 0.45 + target_score * 0.35 + pathway_score * 0.20) * conf_weight
}

####################
# Build final scored table for a dataset
####################

build_scored_table <- function(target_expr, pathway_scores, state_order_vec, dataset_label) {
  message("\n=== Building scored table for ", dataset_label, " ===")

  final_table <- drug_table %>%
    mutate(
      has_human_target = !is.na(target_genes) & target_genes != "non_human",
      max_primary_state_score = pmax(state1_score, state4_score, na.rm = TRUE),
      dominant_primary_state = case_when(
        state1_score > state4_score ~ "State 1 (Proliferative)",
        state4_score > state1_score ~ "State 4 (Stress-adaptive)",
        state1_score == state4_score & state1_score >= 2 ~ "Dual (S1+S4)",
        TRUE ~ "Neither"
      ),
      state_specificity = case_when(
        state1_score >= 2 & state4_score >= 2 ~ "Both S1 and S4",
        state1_score >= 2 ~ "State 1-specific",
        state4_score >= 2 ~ "State 4-specific",
        state2_score >= 2 | state3_score >= 2 ~ "Other",
        TRUE ~ "Other"
      )
    )

  # Join target expression — map scAtlas "Classic proliferation" / PDO same
  if (!is.null(target_expr) && nrow(target_expr) > 0) {
    compound_target_expr <- drug_genes %>%
      left_join(target_expr, by = "gene", relationship = "many-to-many") %>%
      group_by(compound, state) %>%
      summarise(
        mean_target_expr_in = mean(mean_expr_in_state, na.rm = TRUE),
        mean_target_pct_in = mean(pct_expressing_in, na.rm = TRUE),
        mean_target_log2fc = mean(log2fc, na.rm = TRUE),
        n_targets_measured = sum(!is.na(mean_expr_in_state)),
        .groups = "drop"
      )

    s1_expr <- compound_target_expr %>%
      filter(state == "Classic proliferation") %>%
      select(compound,
             s1_target_expr = mean_target_expr_in,
             s1_target_pct = mean_target_pct_in,
             s1_target_log2fc = mean_target_log2fc)

    s4_expr <- compound_target_expr %>%
      filter(state == "Stress-adaptive") %>%
      select(compound,
             s4_target_expr = mean_target_expr_in,
             s4_target_pct = mean_target_pct_in,
             s4_target_log2fc = mean_target_log2fc)

    # Also get expression in ALL states for Excel
    all_state_expr <- compound_target_expr %>%
      select(compound, state, mean_target_pct_in, mean_target_expr_in, mean_target_log2fc)

    final_table <- final_table %>%
      left_join(s1_expr, by = "compound") %>%
      left_join(s4_expr, by = "compound")

    # Add per-state expression columns for Excel
    for (s in state_order_vec) {
      s_safe <- safe_name(s)
      s_data <- all_state_expr %>% filter(state == s)
      pct_col  <- paste0("pct_", s_safe)
      mean_col <- paste0("mean_", s_safe)
      fc_col   <- paste0("fc_", s_safe)
      final_table <- final_table %>%
        left_join(
          s_data %>% select(
            compound,
            !!pct_col  := mean_target_pct_in,
            !!mean_col := mean_target_expr_in,
            !!fc_col   := mean_target_log2fc
          ),
          by = "compound"
        )
    }
  }

  # Join PROGENy pathway scores
  if (!is.null(pathway_scores) && nrow(pathway_scores) > 0) {
    compound_pathway <- bind_rows(lapply(seq_len(nrow(drug_table)), function(i) {
      pw <- drug_table$pathway[i]
      progeny_pws <- get_progeny_for_pathway(pw)
      if (is.null(progeny_pws)) return(NULL)

      pw_data <- pathway_scores %>% filter(pathway %in% progeny_pws)
      if (nrow(pw_data) == 0) return(NULL)

      pw_data %>%
        group_by(state) %>%
        summarise(
          pathway_score_diff = mean(score_diff),
          matched_progeny = paste(unique(pathway), collapse = ";"),
          .groups = "drop"
        ) %>%
        mutate(compound = drug_table$compound[i])
    }))

    if (nrow(compound_pathway) > 0) {
      s1_pw <- compound_pathway %>%
        filter(state == "Classic proliferation") %>%
        select(compound, s1_pathway_score = pathway_score_diff)

      s4_pw <- compound_pathway %>%
        filter(state == "Stress-adaptive") %>%
        select(compound, s4_pathway_score = pathway_score_diff)

      final_table <- final_table %>%
        left_join(s1_pw, by = "compound") %>%
        left_join(s4_pw, by = "compound")
    }
  }

  # Ensure columns exist
  ensure_col <- function(df, col, default = NA_real_) {
    if (!col %in% colnames(df)) df[[col]] <- default
    df
  }
  for (col in c("s1_target_pct", "s1_target_log2fc", "s1_pathway_score",
                "s4_target_pct", "s4_target_log2fc", "s4_pathway_score",
                "s1_target_expr", "s4_target_expr")) {
    final_table <- ensure_col(final_table, col)
  }

  # Compute composite scores
  final_table <- final_table %>%
    mutate(
      s1_composite = compute_composite(state1_score, s1_target_pct, s1_target_log2fc,
                                        s1_pathway_score, confidence),
      s4_composite = compute_composite(state4_score, s4_target_pct, s4_target_log2fc,
                                        s4_pathway_score, confidence),
      max_composite = pmax(s1_composite, s4_composite, na.rm = TRUE)
    ) %>%
    arrange(desc(max_composite)) %>%
    mutate(overall_rank = row_number())

  # Add S1 and S4 ranks
  final_table <- final_table %>%
    mutate(
      s1_rank = rank(-s1_composite, ties.method = "min"),
      s4_rank = rank(-s4_composite, ties.method = "min")
    )

  final_table
}

####################
# Run computation for both datasets
####################

if (!replot_only) {

  ## ---- scAtlas ----
  message("\n========== scATLAS DATASET ==========\n")

  scatlas_expr_cache <- file.path(cache_dir, "scatlas_target_expression.rds")
  scatlas_pw_cache   <- file.path(cache_dir, "scatlas_progeny_scores.rds")

  if (file.exists(scatlas_expr_cache) && file.exists(scatlas_pw_cache) && !force_rebuild) {
    message("  Loading cached scAtlas results")
    scatlas_target_expr <- readRDS(scatlas_expr_cache)
    scatlas_pw_scores   <- readRDS(scatlas_pw_cache)
  } else {
    message("  Loading scAtlas epithelial Seurat...")
    epi <- readRDS("EAC_Ref_epi.rds")
    states <- readRDS("Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds")

    common_cells <- intersect(colnames(epi), names(states))
    keep_cells <- common_cells[states[common_cells] %in% scatlas_state_order]
    message("  Using ", length(keep_cells), " cells with defined states")

    epi <- subset(epi, cells = keep_cells)
    epi$state <- factor(states[colnames(epi)], levels = scatlas_state_order)
    DefaultAssay(epi) <- "RNA"
    Idents(epi) <- "state"

    scatlas_target_expr <- compute_target_expression(epi, epi$state, scatlas_state_order, "scAtlas")
    saveRDS(scatlas_target_expr, scatlas_expr_cache)

    if (file.exists(scatlas_pw_cache) && !force_rebuild) {
      message("  Loading cached scAtlas PROGENy scores")
      scatlas_pw_scores <- readRDS(scatlas_pw_cache)
    } else {
      scatlas_pw_scores <- compute_progeny_scores(epi, epi$state, scatlas_state_order, "scAtlas")
      if (!is.null(scatlas_pw_scores)) saveRDS(scatlas_pw_scores, scatlas_pw_cache)
    }

    rm(epi, states); invisible(gc())
  }

  fwrite(scatlas_target_expr, file.path(table_dir, "opnme_scatlas_target_expression_by_state.csv"))
  if (!is.null(scatlas_pw_scores)) fwrite(scatlas_pw_scores, file.path(table_dir, "opnme_scatlas_pathway_activity_by_state.csv"))

  scatlas_final <- build_scored_table(scatlas_target_expr, scatlas_pw_scores, scatlas_state_order, "scAtlas")
  fwrite(scatlas_final, file.path(table_dir, "opnme_scatlas_final_ranked.csv"))
  saveRDS(scatlas_final, file.path(cache_dir, "scatlas_final_ranked.rds"))

  ## ---- PDO (untreated only) ----
  message("\n========== PDO DATASET (UNTREATED) ==========\n")

  pdo_expr_cache <- file.path(cache_dir, "pdo_target_expression.rds")
  pdo_pw_cache   <- file.path(cache_dir, "pdo_progeny_scores.rds")

  if (file.exists(pdo_expr_cache) && file.exists(pdo_pw_cache) && !force_rebuild) {
    message("  Loading cached PDO results")
    pdo_target_expr <- readRDS(pdo_expr_cache)
    pdo_pw_scores   <- readRDS(pdo_pw_cache)
  } else {
    message("  Loading PDO merged Seurat...")
    pdo_obj <- readRDS(file.path(pdo_dir, "PDOs_outs", "PDOs_merged.rds"))
    pdo_states <- readRDS(file.path(pdo_dir, "PDOs_outs", "centred_mp_refinement",
                                     "centred_refined_noreg_states.rds"))

    # Filter to untreated only
    treatment_col <- NULL
    for (tc in c("Treatment", "treatment")) {
      if (tc %in% colnames(pdo_obj@meta.data)) {
        treatment_col <- tc
        break
      }
    }

    if (!is.null(treatment_col)) {
      untreated_cells <- colnames(pdo_obj)[pdo_obj@meta.data[[treatment_col]] == "Untreated"]
      message("  Untreated cells (from metadata): ", length(untreated_cells), " / ", ncol(pdo_obj))
    } else {
      # Fallback to orig.ident parsing (per PDO SCENIC script)
      message("  No treatment column found; deriving from orig.ident")
      pdo_obj$treatment <- ifelse(grepl("_Treated_", pdo_obj$orig.ident, ignore.case=TRUE) | 
                                    grepl("_FLOT_", pdo_obj$orig.ident, ignore.case=TRUE), 
                                  "Treated", "Untreated")
      untreated_cells <- colnames(pdo_obj)[pdo_obj$treatment == "Untreated"]
      message("  Untreated cells (derived): ", length(untreated_cells), " / ", ncol(pdo_obj))
    }

    common_cells <- Reduce(intersect, list(untreated_cells, names(pdo_states), colnames(pdo_obj)))
    keep_cells <- common_cells[pdo_states[common_cells] %in% pdo_state_order]
    message("  Using ", length(keep_cells), " untreated cells with defined states")

    pdo_obj <- subset(pdo_obj, cells = keep_cells)
    pdo_obj$state <- factor(pdo_states[colnames(pdo_obj)], levels = pdo_state_order)
    DefaultAssay(pdo_obj) <- "RNA"
    Idents(pdo_obj) <- "state"

    pdo_target_expr <- compute_target_expression(pdo_obj, pdo_obj$state, pdo_state_order, "PDO")
    saveRDS(pdo_target_expr, pdo_expr_cache)

    if (file.exists(pdo_pw_cache) && !force_rebuild) {
      message("  Loading cached PDO PROGENy scores")
      pdo_pw_scores <- readRDS(pdo_pw_cache)
    } else {
      pdo_pw_scores <- compute_progeny_scores(pdo_obj, pdo_obj$state, pdo_state_order, "PDO")
      if (!is.null(pdo_pw_scores)) saveRDS(pdo_pw_scores, pdo_pw_cache)
    }

    rm(pdo_obj, pdo_states); invisible(gc())
  }

  fwrite(pdo_target_expr, file.path(table_dir, "opnme_pdo_target_expression_by_state.csv"))
  if (!is.null(pdo_pw_scores)) fwrite(pdo_pw_scores, file.path(table_dir, "opnme_pdo_pathway_activity_by_state.csv"))

  pdo_final <- build_scored_table(pdo_target_expr, pdo_pw_scores, pdo_state_order, "PDO")
  fwrite(pdo_final, file.path(table_dir, "opnme_pdo_final_ranked.csv"))
  saveRDS(pdo_final, file.path(cache_dir, "pdo_final_ranked.rds"))

  ## ---- Combined ranking ----
  combined <- scatlas_final %>%
    select(compound, primary_target, mechanism, pathway, priority_tier, confidence,
           state_specificity, state1_score, state2_score, state3_score, state4_score,
           scatlas_s1_composite = s1_composite, scatlas_s4_composite = s4_composite,
           scatlas_s1_rank = s1_rank, scatlas_s4_rank = s4_rank,
           scatlas_overall_rank = overall_rank) %>%
    left_join(
      pdo_final %>%
        select(compound,
               pdo_s1_composite = s1_composite, pdo_s4_composite = s4_composite,
               pdo_s1_rank = s1_rank, pdo_s4_rank = s4_rank,
               pdo_overall_rank = overall_rank),
      by = "compound"
    ) %>%
    mutate(
      combined_s1_rank = scatlas_s1_rank + pdo_s1_rank,
      combined_s4_rank = scatlas_s4_rank + pdo_s4_rank,
      combined_overall = combined_s1_rank + combined_s4_rank
    ) %>%
    arrange(combined_overall)

  fwrite(combined, file.path(table_dir, "opnme_combined_ranking.csv"))
  saveRDS(combined, file.path(cache_dir, "combined_ranking.rds"))

} else {
  message("=== REPLOT_ONLY mode: loading cached results ===")
  scatlas_final <- readRDS(file.path(cache_dir, "scatlas_final_ranked.rds"))
  pdo_final     <- readRDS(file.path(cache_dir, "pdo_final_ranked.rds"))
  combined      <- readRDS(file.path(cache_dir, "combined_ranking.rds"))
  scatlas_pw_scores <- tryCatch(readRDS(file.path(cache_dir, "scatlas_progeny_scores.rds")), error = function(e) NULL)
  pdo_pw_scores     <- tryCatch(readRDS(file.path(cache_dir, "pdo_progeny_scores.rds")), error = function(e) NULL)
}

####################
# Figures
####################

message("\n=== Generating publication figures ===")

# --- Figure: Log2FC vs Expression scatter (one per dataset) ---

make_fc_vs_expr_plot <- function(final_tbl, dataset_label) {
  plot_data_s1 <- final_tbl %>%
    filter(has_human_target, state1_score >= 2) %>%
    mutate(
      target_pct = ifelse(is.na(s1_target_pct), 0, s1_target_pct),
      target_fc = ifelse(is.na(s1_target_log2fc), 0, s1_target_log2fc),
      composite = s1_composite,
      panel = "State 1: Classic proliferation"
    )

  plot_data_s4 <- final_tbl %>%
    filter(has_human_target, state4_score >= 2) %>%
    mutate(
      target_pct = ifelse(is.na(s4_target_pct), 0, s4_target_pct),
      target_fc = ifelse(is.na(s4_target_log2fc), 0, s4_target_log2fc),
      composite = s4_composite,
      panel = "State 4: Stress-adaptive"
    )

  plot_data <- bind_rows(plot_data_s1, plot_data_s4) %>%
    group_by(panel) %>%
    mutate(
      rank_fc = rank(-target_fc, ties.method = "min"),
      rank_pct = rank(-target_pct, ties.method = "min"),
      rank_comp = rank(-composite, ties.method = "min"),
      label = ifelse(rank_fc <= 3 | rank_pct <= 3 | rank_comp <= 5, compound, NA_character_)
    ) %>%
    ungroup() %>%
    mutate(panel = factor(panel, levels = c("State 1: Classic proliferation", "State 4: Stress-adaptive")))

  p <- ggplot(plot_data, aes(x = target_fc, y = target_pct)) +
    geom_point(aes(color = composite), size = 3, alpha = 0.8, stroke = 0.3) +
    scale_color_distiller(palette = "RdYlBu", name = "Composite\nscore") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    facet_wrap(~ panel, ncol = 2) +
    labs(
      x = "Target gene log2FC (state vs rest)",
      y = "Target gene expression (% cells in state)",
      title = paste0(dataset_label, ": Target Expression vs Log2FC")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = 3, max.overlaps = 20, min.segment.length = 0,
        box.padding = 0.4, segment.alpha = 0.4
      )
  }
  p
}

p_scatlas_fc <- make_fc_vs_expr_plot(scatlas_final, "scAtlas")
ggsave(file.path(fig_dir, "opnme_scatlas_fc_vs_expression.pdf"),
       p_scatlas_fc, width = 14, height = 7, useDingbats = FALSE)
ggsave(file.path(fig_dir, "opnme_scatlas_fc_vs_expression.png"),
       p_scatlas_fc, width = 14, height = 7, dpi = 300)

p_pdo_fc <- make_fc_vs_expr_plot(pdo_final, "PDO (untreated)")
ggsave(file.path(fig_dir, "opnme_pdo_fc_vs_expression.pdf"),
       p_pdo_fc, width = 14, height = 7, useDingbats = FALSE)
ggsave(file.path(fig_dir, "opnme_pdo_fc_vs_expression.png"),
       p_pdo_fc, width = 14, height = 7, dpi = 300)

# --- Figure: Prioritisation Heatmap (one per dataset) ---

make_prio_heatmap <- function(final_tbl, dataset_label) {
  hm_data <- final_tbl %>%
    filter(priority_tier %in% c("Tier 1", "Tier 2"), has_human_target) %>%
    arrange(desc(max_composite))

  if (nrow(hm_data) == 0) return(invisible(NULL))

  hm_mat <- hm_data %>%
    select(compound, state1_score, state2_score, state3_score, state4_score) %>%
    column_to_rownames("compound") %>%
    as.matrix()

  colnames(hm_mat) <- c("State 1\nProliferation", "State 2\nColumnar>Intest.",
                          "State 3\nGlandular>Intest.", "State 4\nStress-adaptive")

  conf_colors <- c("High" = "#2166AC", "Medium" = "#F4A582", "Low" = "#D6604D")
  tier_colors <- c("Tier 1" = "#1B7837", "Tier 2" = "#A6D96A")
  spec_colors <- c("Both S1 and S4" = "#7570B3",
                    "State 1-specific" = "#E41A1C",
                    "State 4-specific" = "#984EA3",
                    "Other" = "#999999")

  row_ha <- rowAnnotation(
    Tier = hm_data$priority_tier,
    Confidence = hm_data$confidence,
    Specificity = hm_data$state_specificity,
    `S1 composite` = anno_barplot(hm_data$s1_composite, bar_width = 0.7,
                                   gp = gpar(fill = "#E41A1C")),
    `S4 composite` = anno_barplot(hm_data$s4_composite, bar_width = 0.7,
                                   gp = gpar(fill = "#984EA3")),
    col = list(
      Tier = tier_colors,
      Confidence = conf_colors,
      Specificity = spec_colors
    ),
    annotation_name_gp = gpar(fontsize = 9),
    gap = unit(2, "mm")
  )

  score_col <- colorRamp2(c(0, 1, 2, 3), c("#F7F7F7", "#FEE08B", "#F46D43", "#A50026"))

  ht <- Heatmap(
    hm_mat,
    name = "Prior\nscore",
    col = score_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 35,
    right_annotation = row_ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.0f", hm_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    },
    row_title = "Compounds (ranked by composite score)",
    column_title = paste0(dataset_label, ": Mechanism Priors + Data Validation"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    width = unit(8, "cm"),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9))
  )

  prefix <- ifelse(grepl("PDO", dataset_label), "pdo", "scatlas")
  pdf(file.path(fig_dir, paste0("opnme_", prefix, "_prioritisation_heatmap.pdf")),
      width = 16, height = max(8, nrow(hm_data) * 0.28))
  draw(ht, padding = unit(c(5, 20, 5, 5), "mm"))
  dev.off()

  png(file.path(fig_dir, paste0("opnme_", prefix, "_prioritisation_heatmap.png")),
      width = 16, height = max(8, nrow(hm_data) * 0.28), units = "in", res = 300)
  draw(ht, padding = unit(c(5, 20, 5, 5), "mm"))
  dev.off()
}

make_prio_heatmap(scatlas_final, "scAtlas")
make_prio_heatmap(pdo_final, "PDO (untreated)")

# --- Figure: PROGENy Pathway Heatmap (one per dataset) ---

make_pathway_heatmap <- function(pw_scores, state_order_vec, dataset_label) {
  if (is.null(pw_scores) || nrow(pw_scores) == 0) return(invisible(NULL))

  pw_mat <- pw_scores %>%
    select(pathway, state, score_diff) %>%
    pivot_wider(names_from = state, values_from = score_diff) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  # Ensure column order matches state order
  pw_mat <- pw_mat[, intersect(state_order_vec, colnames(pw_mat)), drop = FALSE]

  pw_col <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

  ht_pw <- Heatmap(
    pw_mat,
    name = "Score\ndiff",
    col = pw_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 35,
    column_title = paste0(dataset_label, ": PROGENy Pathway Activity (State vs Rest)"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", pw_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    },
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9))
  )

  prefix <- ifelse(grepl("PDO", dataset_label), "pdo", "scatlas")
  pdf(file.path(fig_dir, paste0("opnme_", prefix, "_pathway_state_heatmap.pdf")),
      width = 10, height = 8)
  draw(ht_pw, padding = unit(c(5, 5, 5, 5), "mm"))
  dev.off()

  png(file.path(fig_dir, paste0("opnme_", prefix, "_pathway_state_heatmap.png")),
      width = 10, height = 8, units = "in", res = 300)
  draw(ht_pw, padding = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
}

make_pathway_heatmap(scatlas_pw_scores, scatlas_state_order, "scAtlas")
make_pathway_heatmap(pdo_pw_scores, pdo_state_order, "PDO (untreated)")

# --- Figure: Rank in scAtlas vs Rank in PDO scatter ---

message("  Generating rank comparison scatter...")

rank_plot_data <- combined %>%
  filter(!is.na(scatlas_overall_rank) & !is.na(pdo_overall_rank)) %>%
  mutate(
    combined_rank = rank(combined_overall, ties.method = "min"),
    label = ifelse(combined_rank <= 15, compound, NA_character_),
    specificity = state_specificity
  )

p_rank <- ggplot(rank_plot_data, aes(x = scatlas_overall_rank, y = pdo_overall_rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = specificity), size = 3, alpha = 0.8) +
  scale_x_reverse() +
  scale_y_reverse() +
  scale_color_manual(
    values = c("Both S1 and S4" = "#7570B3",
               "State 1-specific" = "#E41A1C",
               "State 4-specific" = "#984EA3",
               "Other" = "#999999"),
    name = "Specificity"
  ) +
  labs(
    x = "Rank in scAtlas",
    y = "Rank in PDO (untreated)",
    title = "Drug Rank Comparison: scAtlas vs PDO"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_rank <- p_rank +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 3, max.overlaps = 25, min.segment.length = 0,
      box.padding = 0.4, segment.alpha = 0.4
    )
}

ggsave(file.path(fig_dir, "opnme_rank_scatlas_vs_pdo.pdf"),
       p_rank, width = 10, height = 10, useDingbats = FALSE)
ggsave(file.path(fig_dir, "opnme_rank_scatlas_vs_pdo.png"),
       p_rank, width = 10, height = 10, dpi = 300)

####################
# Excel workbooks
####################

message("\n=== Generating Excel workbooks ===")

# Color functions
get_gradient_color <- function(val, min_val, max_val) {
  if (is.null(val) || is.na(val)) return(NULL)
  
  if (is.character(val)) {
    val <- as.numeric(gsub("%", "", val))
  }
  if (!is.numeric(val) || is.na(val)) return(NULL)
  
  val <- min(max(val, min_val), max_val)
  frac <- (val - min_val) / (max_val - min_val)
  cr <- colorRamp(c("#F7F7F7", "#FEE08B", "#F46D43", "#A50026"))
  rgb_vals <- cr(frac)
  list(
    fill = sprintf("#%02X%02X%02XFF", as.integer(rgb_vals[1]), as.integer(rgb_vals[2]), as.integer(rgb_vals[3])),
    font = if (frac > 0.72) "#FFFFFF" else "#000000"
  )
}
score_fill     <- function(val) get_gradient_color(val, 0, 3)
pct_fill       <- function(val) get_gradient_color(val, 0, 50)
mean_fill      <- function(val) get_gradient_color(val, 0, 2.5)
fc_fill        <- function(val) get_gradient_color(val, 0, 1.5)
composite_fill <- function(val) get_gradient_color(val, 0, 0.8)

pathway_fill <- function(val) {
  if (is.null(val) || is.na(val) || !is.numeric(val)) return(NULL)
  capped <- min(max(val, -2), 2)
  frac <- (capped + 2) / 4
  r <- as.integer(33 + frac * (178 - 33))
  g <- as.integer(102 + frac * (24 - 102))
  b <- as.integer(172 + frac * (43 - 172))
  list(
    fill = sprintf("#%02X%02X%02XFF", r, g, b),
    font = if (frac < 0.25 || frac > 0.75) "#FFFFFF" else "#000000"
  )
}

build_excel <- function(final_tbl, state_order_vec, dataset_label) {
  prefix <- ifelse(grepl("PDO", dataset_label, ignore.case = TRUE), "pdo", "scatlas")

  wb <- createWorkbook()

  # Prepare base columns
  base_cols <- c("compound", "primary_target", "mechanism", "pathway")

  # Build sheet data
  build_sheet_data <- function(tbl) {
    # Pathway score columns
    pw_s1 <- if ("s1_pathway_score" %in% colnames(tbl)) tbl$s1_pathway_score else rep(NA_real_, nrow(tbl))
    pw_s4 <- if ("s4_pathway_score" %in% colnames(tbl)) tbl$s4_pathway_score else rep(NA_real_, nrow(tbl))

    # State expression columns (for all states)
    pct_cols  <- list()
    mean_cols <- list()
    fc_cols   <- list()
    for (s in state_order_vec) {
      s_safe <- safe_name(s)
      pcol <- paste0("pct_", s_safe)
      mcol <- paste0("mean_", s_safe)
      fcol <- paste0("fc_", s_safe)
      pct_cols[[s]]  <- if (pcol %in% colnames(tbl)) round(tbl[[pcol]], 1) else rep(NA_real_, nrow(tbl))
      mean_cols[[s]] <- if (mcol %in% colnames(tbl)) round(tbl[[mcol]], 2) else rep(NA_real_, nrow(tbl))
      fc_cols[[s]]   <- if (fcol %in% colnames(tbl)) round(tbl[[fcol]], 2) else rep(NA_real_, nrow(tbl))
    }

    # Build data frame
    df <- data.frame(
      Compound = tbl$compound,
      `Primary Targets` = tbl$primary_target,
      `Proposed Target Genes` = tbl$target_genes,
      Mechanism = tbl$mechanism,
      Pathway = tbl$pathway,
      `PROGENy S1` = round(pw_s1, 2),
      `PROGENy S4` = round(pw_s4, 2),
      `Prior S1` = tbl$state1_score,
      `Prior S2` = tbl$state2_score,
      `Prior S3` = tbl$state3_score,
      `Prior S4` = tbl$state4_score,
      Tier = tbl$priority_tier,
      Confidence = tbl$confidence,
      Specificity = tbl$state_specificity,
      `Composite S1` = round(tbl$s1_composite, 3),
      `Composite S4` = round(tbl$s4_composite, 3),
      ` ` = rep("", nrow(tbl)),
      check.names = FALSE, stringsAsFactors = FALSE
    )

    # 1. Percentage Expression per state (% cells expressing target)
    for (s in state_order_vec) {
      short <- gsub(" ", "_", substr(s, 1, 15))
      df[[paste0("Pct_", short, " (%)")]] <- pct_cols[[s]]
    }

    df[["  "]] <- rep("", nrow(tbl))

    # 2. Mean Expression per state (log-normalized average expression)
    for (s in state_order_vec) {
      short <- gsub(" ", "_", substr(s, 1, 15))
      df[[paste0("Mean_", short)]] <- mean_cols[[s]]
    }

    df[["   "]] <- rep("", nrow(tbl))

    # 3. Log2 Fold Change per state (state vs rest)
    for (s in state_order_vec) {
      short <- gsub(" ", "_", substr(s, 1, 15))
      df[[paste0("FC_", short)]] <- fc_cols[[s]]
    }

    df
  }

  # --- Sheet 1: All 99 compounds ---
  addWorksheet(wb, "All 99 Compounds")
  sheet_data_all <- build_sheet_data(final_tbl)
  writeData(wb, "All 99 Compounds", sheet_data_all, headerStyle = createStyle(textDecoration = "bold"))

  # --- Sheet 2: Ranked by S1 composite ---
  s1_ranked <- final_tbl %>% arrange(desc(s1_composite))
  addWorksheet(wb, "Ranked S1 Proliferative")
  sheet_data_s1 <- build_sheet_data(s1_ranked)
  writeData(wb, "Ranked S1 Proliferative", sheet_data_s1, headerStyle = createStyle(textDecoration = "bold"))

  # --- Sheet 3: Ranked by S4 composite ---
  s4_ranked <- final_tbl %>% arrange(desc(s4_composite))
  addWorksheet(wb, "Ranked S4 Stress-adaptive")
  sheet_data_s4 <- build_sheet_data(s4_ranked)
  writeData(wb, "Ranked S4 Stress-adaptive", sheet_data_s4, headerStyle = createStyle(textDecoration = "bold"))

  # Apply styling to all sheets
  n_states <- length(state_order_vec)

  for (sheet_name in c("All 99 Compounds", "Ranked S1 Proliferative", "Ranked S4 Stress-adaptive")) {
    sheet_data <- switch(sheet_name,
                          "All 99 Compounds" = sheet_data_all,
                          "Ranked S1 Proliferative" = sheet_data_s1,
                          "Ranked S4 Stress-adaptive" = sheet_data_s4)

    nr <- nrow(sheet_data)

    apply_style <- function(wb, sheet, res, r, c) {
      if (!is.null(res)) {
        sty <- createStyle(fgFill = substr(res$fill, 1, 7), fontColour = res$font)
        addStyle(wb, sheet, sty, rows = r, cols = c)
      }
    }

    # Style PROGENy columns (6-7)
    for (row_i in seq_len(nr)) {
      for (col_i in 6:7) {
        apply_style(wb, sheet_name, pathway_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
      }
    }

    # Style prior score columns (8-11)
    for (row_i in seq_len(nr)) {
      for (col_i in 8:11) {
        apply_style(wb, sheet_name, score_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
      }
    }

    # Style composite score columns (15-16)
    for (row_i in seq_len(nr)) {
      for (col_i in 15:16) {
        apply_style(wb, sheet_name, composite_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
      }
    }

    # Block 1: Percentage Expression (cols 18 to 18 + n_states - 1)
    pct_start <- 18
    for (row_i in seq_len(nr)) {
      for (j in seq_len(n_states)) {
        col_i <- pct_start + j - 1
        if (col_i <= ncol(sheet_data)) {
          apply_style(wb, sheet_name, pct_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
        }
      }
    }

    # Block 2: Mean Expression (cols pct_start + n_states + 1 to pct_start + 2*n_states)
    mean_start <- pct_start + n_states + 1
    for (row_i in seq_len(nr)) {
      for (j in seq_len(n_states)) {
        col_i <- mean_start + j - 1
        if (col_i <= ncol(sheet_data)) {
          apply_style(wb, sheet_name, mean_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
        }
      }
    }

    # Block 3: FC columns (cols mean_start + n_states + 1 to mean_start + 2*n_states)
    fc_start <- mean_start + n_states + 1
    for (row_i in seq_len(nr)) {
      for (j in seq_len(n_states)) {
        col_i <- fc_start + j - 1
        if (col_i <= ncol(sheet_data)) {
          apply_style(wb, sheet_name, fc_fill(sheet_data[row_i, col_i]), row_i + 1, col_i)
        }
      }
    }

    # Auto-width
    setColWidths(wb, sheet_name, cols = 1:ncol(sheet_data), widths = "auto")
    # Make gap columns narrow
    gap1_col <- 17
    gap2_col <- pct_start + n_states
    gap3_col <- mean_start + n_states
    setColWidths(wb, sheet_name, cols = gap1_col, widths = 3)
    if (gap2_col <= ncol(sheet_data)) setColWidths(wb, sheet_name, cols = gap2_col, widths = 3)
    if (gap3_col <= ncol(sheet_data)) setColWidths(wb, sheet_name, cols = gap3_col, widths = 3)
  }

  out_path <- file.path(excel_dir, paste0("opnme_", prefix, "_prioritisation.xlsx"))
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message("  Saved: ", out_path)
}

build_excel(scatlas_final, scatlas_state_order, "scAtlas")
build_excel(pdo_final, pdo_state_order, "PDO")

####################
# Summary report
####################

message("\n=== Final Summary ===")

for (ds_name in c("scAtlas", "PDO")) {
  ft <- if (ds_name == "scAtlas") scatlas_final else pdo_final
  message("\n--- ", ds_name, " ---")
  message("  S1 Tier1 candidates (score >= 2): ", nrow(ft %>% filter(priority_tier == "Tier 1", state1_score >= 2)))
  message("  S4 Tier1 candidates (score >= 2): ", nrow(ft %>% filter(priority_tier == "Tier 1", state4_score >= 2)))
  message("  Top S1: ", paste(head((ft %>% arrange(desc(s1_composite)))$compound, 5), collapse = ", "))
  message("  Top S4: ", paste(head((ft %>% arrange(desc(s4_composite)))$compound, 5), collapse = ", "))
}

message("\n--- Cross-dataset top combined rank ---")
message("  Top 10: ", paste(head(combined$compound, 10), collapse = ", "))

summary_lines <- c(
  "opnMe 99-Compound Dual-Dataset Prioritisation Summary",
  paste("Date:", Sys.time()),
  "",
  "Datasets: scAtlas (epithelial) + PDO (untreated only)",
  paste("Total compounds:", nrow(drug_table)),
  "",
  "Layers used:",
  "  Layer 1: Mechanism-based priors (supervisor workbook 0-3 scores)",
  "  Layer 2: Target gene expression validation (per dataset)",
  "  Layer 3: PROGENy pathway activity scoring (per dataset)",
  "",
  paste("Top 10 combined rank:", paste(head(combined$compound, 10), collapse = ", "))
)
writeLines(summary_lines, file.path(table_dir, "opnme_prioritisation_summary.txt"))

updates_dir <- file.path("updates", "new_updates", "summaries")
dir.create(updates_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(summary_lines, file.path(updates_dir, "opnme_prioritisation_summary.txt"))

message("\n=== opnMe dual-dataset prioritisation complete ===")
message("Outputs in: ref_outs/opnMe_prioritisation/")
