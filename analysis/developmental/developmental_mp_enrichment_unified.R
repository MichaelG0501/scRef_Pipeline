####################
# Analysis registry:
#   Status: active
#   Script: analysis/developmental/developmental_mp_enrichment_unified.R
#   Description: Unified developmental reference validation for scATLAS metaprogrammes.
#   Methodology: analysis/methodology/developmental/developmental_mp_enrichment_unified_methodology.md
#   Inputs:
#     - ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#     - ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#     - ref_outs/EAC_Ref_epi.rds
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/*.xlsx
#     - /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach_epi_DGE_ordered_pretty.csv
#     - optional annotated external expression objects listed in reference_expression_sources
#   Outputs:
#     - intermediate/Auto_developmental_ranked_references.rds
#     - intermediate/Auto_developmental_ref_ucell_scores_<top50|all>.rds
#     - intermediate/Auto_developmental_external_mp_ucell_scores.rds
#     - tables/Auto_developmental_reference_rank_audit.csv
#     - tables/Auto_developmental_reference_expression_availability.csv
#     - tables/Auto_developmental_reference_celltype_coverage.csv
#     - tables/Auto_developmental_overlap_<top50|all>.csv
#     - tables/Auto_developmental_expression_correlation_<top50|all>.csv
#     - tables/Auto_developmental_external_mp_ucell_summary.csv
#     - figures/Auto_developmental_mp_unified_top50.pdf
#     - figures/Auto_developmental_mp_top50_vs_all_overlap_correlation.pdf
#     - logs/Auto_developmental_mp_enrichment_run_summary.txt
#   Cache/replot:
#     - Set SCREF_FORCE_REBUILD=TRUE to rebuild cached UCell score matrices.
#     - Set SCREF_REPLOT_ONLY=TRUE to reuse existing intermediate tables and only redraw PDFs.
#   Run command:
#     - source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
#     - Rscript analysis/developmental/developmental_mp_enrichment_unified.R
####################

library(Seurat)
library(UCell)
library(Matrix)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(clusterProfiler)
library(BiocParallel)
library(dplyr)

####################
# Paths and parameters
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
ref_out_dir <- file.path(project_dir, "ref_outs")
setwd(ref_out_dir)

out_dir <- file.path(ref_out_dir, "Auto_developmental_mp_enrichment_unified")
intermediate_dir <- file.path(out_dir, "intermediate")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
logs_dir <- file.path(out_dir, "logs")
summary_dir <- file.path(project_dir, "updates/new_updates/summaries")

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

start_time <- Sys.time()
force_rebuild <- identical(Sys.getenv("SCREF_FORCE_REBUILD"), "TRUE")
replot_only <- identical(Sys.getenv("SCREF_REPLOT_ONLY"), "TRUE")
max_cells_per_type <- as.integer(Sys.getenv("SCREF_MAX_CELLS_PER_TYPE", "5000"))
ucell_cores <- as.integer(Sys.getenv("SCREF_UCELL_CORES", "4"))
seed_base <- 20260528L

developmental_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental"
mp_path <- file.path(ref_out_dir, "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
mp_ucell_path <- file.path(ref_out_dir, "Metaprogrammes_Results/UCell_nMP19_filtered.rds")
epi_path <- file.path(ref_out_dir, "EAC_Ref_epi.rds")

reference_order <- c(
  "Early_Embryogenesis",
  "Organogenesis_major",
  "Organogenesis_sub",
  "Normal_Development_long",
  "Normal_Development_short",
  "Adult_Epithelium",
  "Barretts_Oesophagus"
)

reference_pretty <- c(
  Early_Embryogenesis = "Early embryogenesis",
  Organogenesis_major = "Organogenesis major",
  Organogenesis_sub = "Organogenesis endoderm/subtypes",
  Normal_Development_long = "Normal development long markers",
  Normal_Development_short = "Normal development literature markers",
  Adult_Epithelium = "Adult epithelium",
  Barretts_Oesophagus = "Barretts oesophagus"
)

mp_descriptions <- c(
  MP1 = "G2M Cell Cycle",
  MP9 = "G1S Cell Cycle",
  MP2 = "MYC-related Proliferation",
  MP17 = "Basal-like Transition",
  MP14 = "Hypoxia Adapted Epi.",
  MP5 = "Epithelial IFN Resp.",
  MP10 = "Columnar Diff.",
  MP8 = "Intestinal Diff.",
  MP13 = "Hypoxic Inflam. Epi.",
  MP7 = "DNA Damage Repair",
  MP18 = "Secretory Diff. (Intest.)",
  MP16 = "Secretory Diff. (Gastric)",
  MP15 = "Immune Infiltration",
  MP12 = "Neuro-responsive Epi"
)

heat_cols <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
cor_cols <- circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
ucell_cols <- heat_cols
na_col <- "#eeeeee"

####################
# Generic helpers
####################
clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove_all(x, "^'+|'+$")
  stringr::str_trim(x)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

make_term <- function(label, suffix) {
  paste0(gsub(" ", "_", label), "..", suffix)
}

make_safe_name <- function(reference, term, mode) {
  make.names(paste(reference, mode, term, sep = "__"), unique = FALSE)
}

write_csv_safe <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE, quote = TRUE)
}

stars_from_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}

cap_neglog10 <- function(p, cap = 8) {
  p <- ifelse(is.na(p) | p <= 0, NA_real_, p)
  pmin(-log10(p), cap)
}

dedupe_ranked_genes <- function(df) {
  df %>%
    filter(!is.na(reference), !is.na(term), !is.na(gene), gene != "") %>%
    arrange(reference, term, rank, desc(abs_logfc), p_adjust, gene) %>%
    group_by(reference, term, gene) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(reference, term) %>%
    arrange(rank, desc(abs_logfc), p_adjust, gene, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

split_reference_genes <- function(ranked_refs, mode = c("top50", "all")) {
  mode <- match.arg(mode)
  df <- ranked_refs
  if (mode == "top50") {
    df <- df %>% group_by(reference, term) %>% slice_head(n = 50) %>% ungroup()
  }
  split(df$gene, paste(df$reference, df$term, sep = "\t"))
}

reference_term_table <- function(ranked_refs) {
  ranked_refs %>%
    distinct(reference, term) %>%
    mutate(reference = factor(reference, levels = reference_order)) %>%
    arrange(reference, term) %>%
    mutate(reference = as.character(reference))
}

####################
# Ranked developmental reference construction
####################
build_early_embryogenesis <- function() {
  xlsx_path <- file.path(developmental_dir, "Early embryogenesis.xls")
  lineage_order <- c(
    "Z4cell", "8 cell", "Morula", "Prelineage", "ICM", "Epiblast", "PriS",
    "Mesoderm", "Axial Mes", "AdvMes", "Erythroblasts", "DE", "HEP",
    "ExE_Mes", "Amnion", "Hypoblast", "YSE", "TE", "CTB", "EVT", "STB"
  )

  readxl::read_excel(xlsx_path, sheet = " temp.human.mk.tsv") %>%
    mutate(source_row = row_number()) %>%
    transmute(
      reference = "Early_Embryogenesis",
      term = paste0(gsub(" ", "_", as.character(cluster)), "..Early_Embryogenesis"),
      term_order = match(as.character(cluster), lineage_order),
      gene = clean_gene(gene),
      p_adjust = safe_numeric(p_val_adj),
      logfc = safe_numeric(avg_log2FC),
      abs_logfc = abs(logfc),
      source_row = source_row,
      source_file = xlsx_path,
      source_sheet = " temp.human.mk.tsv",
      rank_basis = "source row order within lineage after p_val_adj < 0.05; sheet contains avg_log2FC and p_val_adj"
    ) %>%
    filter(!is.na(p_adjust), p_adjust < 0.05, !is.na(term_order)) %>%
    arrange(term_order, source_row) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

build_organogenesis <- function() {
  xlsx_path <- file.path(developmental_dir, "Organogenesis.xlsx")
  major_order <- c(
    "neural progenitor", "neuron", "epidermis", "sensory neuron",
    "schwann", "craniofacial", "head mesoderm", "somite", "IM",
    "somatic LPM", "limb", "splanchnic LPM", "endothelium", "blood",
    "endoderm", "PGC", "epithelium", "fibroblast"
  )
  specific_order <- c(
    "foregut/esophagus", "thymus", "lung proximal epithelium and trachea",
    "lung distal epithelium", "hepatocyte-1", "hepatocyte-2", "stomach",
    "pancreas", "duodenum", "undefined", "epithelium-1", "epithelium-2",
    "epithelium-3", "epithelium-4", "epithelium-5"
  )
  specific_levels <- ifelse(
    grepl("epithelium", specific_order, ignore.case = TRUE),
    paste0(gsub(" ", "_", specific_order), "..Organogenesis_sub"),
    paste0(gsub(" ", "_", specific_order), "_Endoderm..Organogenesis_sub")
  )

  raw <- readxl::read_excel(xlsx_path, sheet = "S1D") %>%
    select(
      ID = Type_id,
      Major = Developmental_system,
      term = Final_annotation,
      4:ncol(.)
    ) %>%
    mutate(source_row = row_number()) %>%
    pivot_longer(
      cols = -c(ID, Major, term, source_row),
      names_to = "marker_column",
      values_to = "gene",
      values_drop_na = TRUE
    ) %>%
    mutate(
      gene_rank = match(marker_column, unique(marker_column)),
      ID = str_trim(as.character(ID)),
      Major = str_trim(as.character(Major)),
      term = str_trim(as.character(term)),
      gene = clean_gene(gene)
    ) %>%
    filter(!is.na(gene), gene != "", !str_detect(gene, "^\\.\\.\\."))

  major <- raw %>%
    filter(Major %in% major_order) %>%
    transmute(
      reference = "Organogenesis_major",
      term = paste0(gsub(" ", "_", Major), "..Organogenesis_major"),
      term_order = match(Major, major_order),
      gene = gene,
      p_adjust = NA_real_,
      logfc = NA_real_,
      abs_logfc = NA_real_,
      source_row = source_row,
      marker_rank = gene_rank,
      source_file = xlsx_path,
      source_sheet = "S1D",
      rank_basis = "S1D marker columns are labelled as DEGs ordered by z-score; major terms aggregate subtypes by best marker-column rank"
    ) %>%
    group_by(reference, term, gene) %>%
    summarise(
      term_order = min(term_order, na.rm = TRUE),
      source_row = min(source_row, na.rm = TRUE),
      marker_rank = min(marker_rank, na.rm = TRUE),
      p_adjust = dplyr::first(p_adjust),
      logfc = dplyr::first(logfc),
      abs_logfc = dplyr::first(abs_logfc),
      source_file = dplyr::first(source_file),
      source_sheet = dplyr::first(source_sheet),
      rank_basis = dplyr::first(rank_basis),
      .groups = "drop"
    ) %>%
    arrange(term_order, marker_rank, source_row, gene) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  specific <- raw %>%
    filter(term %in% specific_order) %>%
    mutate(
      term_new = if_else(
        grepl("epithelium", term, ignore.case = TRUE),
        paste0(gsub(" ", "_", term), "..Organogenesis_sub"),
        paste0(gsub(" ", "_", term), "_Endoderm..Organogenesis_sub")
      )
    ) %>%
    transmute(
      reference = "Organogenesis_sub",
      term = term_new,
      term_order = match(term_new, specific_levels),
      gene = gene,
      p_adjust = NA_real_,
      logfc = NA_real_,
      abs_logfc = NA_real_,
      source_row = source_row,
      marker_rank = gene_rank,
      source_file = xlsx_path,
      source_sheet = "S1D",
      rank_basis = "S1D marker columns are labelled as DEGs ordered by z-score; subtype terms use marker-column rank"
    ) %>%
    group_by(reference, term, gene) %>%
    summarise(
      term_order = min(term_order, na.rm = TRUE),
      source_row = min(source_row, na.rm = TRUE),
      marker_rank = min(marker_rank, na.rm = TRUE),
      p_adjust = dplyr::first(p_adjust),
      logfc = dplyr::first(logfc),
      abs_logfc = dplyr::first(abs_logfc),
      source_file = dplyr::first(source_file),
      source_sheet = dplyr::first(source_sheet),
      rank_basis = dplyr::first(rank_basis),
      .groups = "drop"
    ) %>%
    arrange(term_order, marker_rank, source_row, gene) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  bind_rows(major, specific)
}

build_normal_development <- function() {
  xlsx_path <- file.path(developmental_dir, "Normal development.xlsx")
  mapping_list <- data.frame(
    organ = c(
      rep("Stomach", 16),
      rep("Intestine", 12),
      rep("Pancreas", 14),
      rep("Lung", 13),
      rep("Liver", 9)
    ),
    original_term = c(
      "Goblet cells", "Parietal and chief cells", "Squamous epithelial cells", "Stromal cells",
      "MUC13_DMBT1 positive cells", "Lymphoid cells", "Vascular endothelial cells",
      "PDE1C_ACSM3 positive cells", "Myeloid cells", "Erythroblasts", "ENS glia",
      "ENS neurons", "Ciliated epithelial cells", "Neuroendocrine cells",
      "Lymphatic endothelial cells", "Mesothelial cells",
      "Intestinal epithelial cells", "Stromal cells", "ENS glia", "ENS neurons",
      "Myeloid cells", "Lymphoid cells", "Vascular endothelial cells", "Smooth muscle cells",
      "Chromaffin cells", "Erythroblasts", "Lymphatic endothelial cells", "Mesothelial cells",
      "Acinar cells", "Stromal cells", "Ductal cells", "Lymphoid cells",
      "Vascular endothelial cells", "Islet endocrine cells", "ENS glia", "ENS neurons",
      "Erythroblasts", "Myeloid cells", "CCL19_CCL21 positive cells", "Mesothelial cells",
      "Lymphatic endothelial cells", "Smooth muscle cells",
      "Bronchiolar and alveolar epithelial cells", "Stromal cells", "Ciliated epithelial cells",
      "Neuroendocrine cells", "Squamous epithelial cells", "Visceral neurons",
      "Myeloid cells", "Lymphoid cells", "Megakaryocytes", "Vascular endothelial cells",
      "Lymphatic endothelial cells", "Mesothelial cells", "CSH1_CSH2 positive cells",
      "Vascular endothelial cells", "Lymphoid cells", "Myeloid cells", "Megakaryocytes",
      "Erythroblasts", "Stellate cells", "Hepatoblasts", "Mesothelial cells",
      "Hematopoietic stem cells"
    ),
    stringsAsFactors = FALSE
  )
  ordered_terms <- mapping_list %>%
    mutate(formatted_name = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development")) %>%
    pull(formatted_name)

  long <- readxl::read_excel(xlsx_path, sheet = "Table_S4", skip = 1) %>%
    mutate(source_row = row_number()) %>%
    transmute(
      original_term = as.character(max.cluster),
      gene = clean_gene(gene_short_name),
      p_adjust = safe_numeric(qval),
      logfc = safe_numeric(fold.change),
      source_row = source_row
    ) %>%
    inner_join(mapping_list, by = "original_term", relationship = "many-to-many") %>%
    filter(!is.na(p_adjust), p_adjust < 0.05, !is.na(gene), gene != "") %>%
    mutate(
      reference = "Normal_Development_long",
      term = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development_long"),
      term_order = match(sub("_long$", "", term), ordered_terms),
      abs_logfc = abs(logfc),
      source_file = xlsx_path,
      source_sheet = "Table_S4",
      rank_basis = "ranked within mapped cell type by descending fold.change, then qval"
    ) %>%
    arrange(term_order, desc(logfc), p_adjust, source_row, gene) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  short <- readxl::read_excel(xlsx_path, sheet = "Table_S3", skip = 1) %>%
    mutate(source_row = row_number()) %>%
    select(
      organ = Organ,
      original_term = `Main cell type annotation`,
      gene_string = `Gene markers supporting annotation`,
      source_row
    ) %>%
    inner_join(mapping_list, by = c("organ", "original_term")) %>%
    mutate(term_raw = paste0(gsub(" ", "_", original_term), "_", gsub(" ", "_", organ), "..Normal_Development")) %>%
    filter(term_raw %in% ordered_terms) %>%
    separate_rows(gene_string, sep = ",\\s*") %>%
    group_by(organ, original_term, term_raw) %>%
    mutate(marker_rank = row_number()) %>%
    ungroup() %>%
    transmute(
      reference = "Normal_Development_short",
      term = paste0(term_raw, "_short"),
      term_order = match(term_raw, ordered_terms),
      gene = clean_gene(gene_string),
      p_adjust = NA_real_,
      logfc = NA_real_,
      abs_logfc = NA_real_,
      source_row = source_row,
      marker_rank = marker_rank,
      source_file = xlsx_path,
      source_sheet = "Table_S3",
      rank_basis = "literature marker string order; this table is not a differential ranked gene list"
    ) %>%
    filter(!is.na(gene), gene != "") %>%
    arrange(term_order, marker_rank, gene) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  bind_rows(long, short)
}

build_adult_epithelium <- function() {
  oesophagus_path <- file.path(developmental_dir, "Oesophagus.xlsx")
  stomach_pretty_path <- file.path(developmental_dir, "Stomach_epi_DGE_ordered_pretty.csv")
  oesophagus_order <- c("Quiescent basal cell", "Basal cell (cycling)", "Suprabasal", "Apical cell")

  oesophagus <- map_dfr(readxl::excel_sheets(oesophagus_path), function(sh) {
    readxl::read_excel(oesophagus_path, sheet = sh) %>%
      mutate(source_row = row_number()) %>%
      transmute(
        reference = "Adult_Epithelium",
        term = paste0(gsub(" ", "_", sh), "_Oesophagus..Adult_Epithelium"),
        term_order = match(sh, oesophagus_order),
        gene = clean_gene(gene),
        p_adjust = safe_numeric(p_val_adj),
        logfc = safe_numeric(avg_log2FC),
        abs_logfc = abs(logfc),
        source_row = source_row,
        source_file = oesophagus_path,
        source_sheet = sh,
        rank_basis = "source row order within sheet after p_val_adj < 0.05; sheet contains avg_log2FC and p_val_adj"
      ) %>%
      filter(!is.na(p_adjust), p_adjust < 0.05, !is.na(gene), gene != "")
  }) %>%
    arrange(term_order, source_row) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  stomach <- read.csv(stomach_pretty_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
    group_by(mapped_cluster) %>%
    mutate(source_row = row_number()) %>%
    ungroup() %>%
    transmute(
      reference = "Adult_Epithelium",
      term = mapped_cluster,
      term_order = 100L + match(mapped_cluster, unique(mapped_cluster)),
      gene = clean_gene(gene),
      p_adjust = safe_numeric(p_val_adj),
      logfc = safe_numeric(avg_log2FC),
      abs_logfc = abs(logfc),
      source_row = source_row,
      source_file = stomach_pretty_path,
      source_sheet = "Stomach_epi_DGE_ordered_pretty.csv",
      rank_basis = "cached FindAllMarkers table ordered by avg_log2FC within epithelial stomach cluster"
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(p_adjust), p_adjust < 0.05) %>%
    arrange(term_order, source_row) %>%
    group_by(reference, term) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  bind_rows(oesophagus, stomach)
}

build_barretts <- function() {
  barrett_base_dir <- file.path(developmental_dir, "Barretts")
  barretts_groups <- list(
    list(
      group_name = "Normal_Esophagus",
      file = "science.abd1449_Table_S4.xlsx",
      sheets = c(
        Basal = "Basal",
        Suprabasal = "Suprabasal",
        Suprabasal_Dividing = "Suprabasal_Dividing",
        Intermediate = "Intermediate",
        Superficial = "Superficial"
      )
    ),
    list(
      group_name = "Normal_Gastric",
      file = "science.abd1449_Table_S5.xlsx",
      sheets = c(
        Undifferentiated = "Undifferentiated",
        Undifferentiated_Dividing = "Undifferentiated_Dividing",
        Foveolar_Intermediate = "Foveolar_Intermediate",
        Foveolar_differentiated = "Foveolar_differentiated",
        Chief = "Chief",
        Parietal = "Parietal",
        Endocrine_GHRL = "Endocrine_GHRL",
        Endocrine_CHGA = "Endocrine_CHGA",
        Endocrine_NEUROD1 = "Endocrine_NEUROD1"
      )
    ),
    list(
      group_name = "Barretts_Esophagus",
      file = "science.abd1449_Table_S7.xlsx",
      sheets = c(
        Columnar_Undifferentiated = "Columnar_Undifferentiated",
        Columnar_Dividing = "Columnar_Undifferentiated_Divid",
        Endocrine_NEUROG3 = "Endocrine_NEUROG3",
        Columnar_Intermediate = "Columnar_Intermediate",
        Columnar_differentiated = "Columnar_differentiated",
        Goblet = "Goblet"
      )
    ),
    list(
      group_name = "Submucosal_Glands",
      file = "science.abd1449_Table_S2.xlsx",
      sheets = c(
        Duct_Intercalating = "Duct_Intercalating",
        Oncocytes_MUC5B_Low = "Oncocytes",
        Mucous_MUC5B_High = "Mucous"
      )
    )
  )

  out <- list()
  term_counter <- 0L
  for (grp in barretts_groups) {
    xlsx_path <- file.path(barrett_base_dir, grp$file)
    available_sheets <- readxl::excel_sheets(xlsx_path)
    for (canonical in names(grp$sheets)) {
      sheet_name <- grp$sheets[[canonical]]
      actual_sheet <- sheet_name
      if (!actual_sheet %in% available_sheets) {
        match_idx <- grep(sheet_name, available_sheets, ignore.case = TRUE)
        if (length(match_idx) == 0) {
          warning("Barretts sheet not found: ", sheet_name, " in ", grp$file)
          next
        }
        actual_sheet <- available_sheets[match_idx[1]]
      }
      term_counter <- term_counter + 1L
      final_term_name <- paste0(canonical, "_", grp$group_name, "..Barretts_Oesophagus")
      df <- readxl::read_excel(xlsx_path, sheet = actual_sheet) %>% mutate(source_row = row_number())
      gene_col <- c("gene", "Symbol", "Genename", "Gene")
      gene_col <- gene_col[gene_col %in% names(df)][1]
      padj_col <- c("p_val_adj", "FDR", "qval", "adj.P.Val", "padj")
      padj_col <- padj_col[padj_col %in% names(df)][1]
      logfc_cols <- grep("_logFC$", names(df), value = TRUE)
      if (is.na(gene_col) || is.na(padj_col)) {
        warning("Missing Barretts gene/FDR columns in ", grp$file, " sheet ", actual_sheet)
        next
      }
      df$Auto_max_logfc <- if (length(logfc_cols) > 0) {
        do.call(pmax, c(lapply(df[logfc_cols], safe_numeric), list(na.rm = TRUE)))
      } else {
        NA_real_
      }
      out[[length(out) + 1L]] <- df %>%
        transmute(
          reference = "Barretts_Oesophagus",
          term = final_term_name,
          term_order = term_counter,
          gene = clean_gene(.data[[gene_col]]),
          p_adjust = safe_numeric(.data[[padj_col]]),
          logfc = Auto_max_logfc,
          abs_logfc = abs(logfc),
          source_row = source_row,
          source_file = xlsx_path,
          source_sheet = actual_sheet,
          rank_basis = "source row order within sheet after FDR < 0.05; sheets contain per-comparison logFC and FDR"
        ) %>%
        filter(!is.na(gene), gene != "", !is.na(p_adjust), p_adjust < 0.05) %>%
        arrange(source_row) %>%
        group_by(reference, term) %>%
        mutate(rank = row_number()) %>%
        ungroup()
    }
  }
  bind_rows(out)
}

ranked_refs_path <- file.path(intermediate_dir, "Auto_developmental_ranked_references.rds")
if (file.exists(ranked_refs_path) && !force_rebuild) {
  ranked_refs <- readRDS(ranked_refs_path)
} else {
  ranked_refs <- bind_rows(
    build_early_embryogenesis(),
    build_organogenesis(),
    build_normal_development(),
    build_adult_epithelium(),
    build_barretts()
  ) %>%
    mutate(
      reference = as.character(reference),
      term = as.character(term),
      gene = clean_gene(gene)
    ) %>%
    dedupe_ranked_genes() %>%
    filter(reference %in% reference_order)

  saveRDS(ranked_refs, ranked_refs_path)
}

rank_audit <- ranked_refs %>%
  group_by(reference, term, source_file, source_sheet, rank_basis) %>%
  summarise(
    n_genes_all = n_distinct(gene),
    n_genes_top50 = min(50L, n_distinct(gene)),
    min_rank = min(rank, na.rm = TRUE),
    max_rank = max(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(reference = factor(reference, levels = reference_order)) %>%
  arrange(reference, term) %>%
  mutate(reference = as.character(reference))

write_csv_safe(rank_audit, file.path(tables_dir, "Auto_developmental_reference_rank_audit.csv"))

####################
# MP loading and ordering
####################
geneNMF.metaprograms <- readRDS(mp_path)
mp_genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp_genes <- mp_genes[!names(mp_genes) %in% paste0("MP", bad_mps)]
}

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_genes)))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_order <- paste0("MP", mp_tree_order)
mp_order <- mp_order[mp_order %in% names(mp_genes)]
mp_genes <- mp_genes[mp_order]
mp_sizes <- vapply(mp_genes, length, integer(1))
mp_labels <- sapply(mp_order, function(mp) {
  desc <- if (mp %in% names(mp_descriptions)) mp_descriptions[[mp]] else ""
  paste0(mp, "\n", desc, "\nn=", mp_sizes[[mp]])
})
names(mp_labels) <- mp_order

####################

####################
# Original-aligned developmental MP enrichment core.
# This file is sourced by developmental_mp_enrichment_unified.R after inputs,
# ranked references, and MP ordering have been prepared.
####################

cat("Running original-aligned developmental MP enrichment core.\n")
library(dplyr)
filter <- dplyr::filter
arrange <- dplyr::arrange
group_by <- dplyr::group_by
slice_head <- dplyr::slice_head
ungroup <- dplyr::ungroup
transmute <- dplyr::transmute
mutate <- dplyr::mutate
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows
full_join <- dplyr::full_join
left_join <- dplyr::left_join
distinct <- dplyr::distinct
pull <- dplyr::pull
n_distinct <- dplyr::n_distinct
if_else <- dplyr::if_else

####################
# Original developmental TERM2GENE references
####################
developmental_ref_dir <- file.path(developmental_dir, "per_stage")
custom_files <- list.files(developmental_ref_dir, pattern = "\\.rds$", full.names = TRUE)
custom_refs_all <- lapply(custom_files, readRDS)
names(custom_refs_all) <- sub(".*enrich_dev_", "", basename(custom_files))
names(custom_refs_all) <- sub("\\.rds$", "", names(custom_refs_all))
custom_refs_all <- custom_refs_all[reference_order[reference_order %in% names(custom_refs_all)]]

if (length(custom_refs_all) == 0) {
  stop("No original per-stage developmental enrichment references were found.")
}

make_top50_custom_refs <- function(custom_refs) {
  ####################
  # Fixed: take top 50 genes per term directly from original TERM2GENE
  # (preserving original row order), not from re-parsed ranked_refs.
  # This ensures the top-50 gene sets exactly match what the original
  # per-stage enrichment files define.
  ####################
  out <- list()
  for (ref_name in names(custom_refs)) {
    term2name <- as.data.frame(custom_refs[[ref_name]]$TERM2NAME, stringsAsFactors = FALSE)
    colnames(term2name)[1:2] <- c("term", "name")
    term_order <- as.character(term2name$term)

    full_term2gene <- as.data.frame(custom_refs[[ref_name]]$TERM2GENE, stringsAsFactors = FALSE)
    colnames(full_term2gene)[1:2] <- c("term", "gene")

    term2gene <- full_term2gene %>%
      filter(term %in% term_order) %>%
      group_by(term) %>%
      slice_head(n = 50) %>%
      ungroup() %>%
      transmute(term = as.character(term), gene = as.character(gene))

    term2gene <- term2gene %>%
      mutate(term = factor(term, levels = term_order)) %>%
      arrange(term) %>%
      mutate(term = as.character(term))

    out[[ref_name]] <- list(
      TERM2GENE = term2gene,
      TERM2NAME = term2name
    )
  }
  out
}

custom_refs_top50 <- make_top50_custom_refs(custom_refs_all)

term_count_note <- bind_rows(lapply(names(custom_refs_all), function(ref_name) {
  all_tbl <- as.data.frame(custom_refs_all[[ref_name]]$TERM2GENE, stringsAsFactors = FALSE)
  top_tbl <- as.data.frame(custom_refs_top50[[ref_name]]$TERM2GENE, stringsAsFactors = FALSE)
  colnames(all_tbl)[1:2] <- c("term", "gene")
  colnames(top_tbl)[1:2] <- c("term", "gene")
  full_join(
    all_tbl %>% dplyr::count(term, name = "n_genes_all_original"),
    top_tbl %>% dplyr::count(term, name = "n_genes_top50"),
    by = "term"
  ) %>%
    mutate(
      reference = ref_name,
      top50_same_as_all = !is.na(n_genes_all_original) & n_genes_all_original <= 50,
      top50_note = ifelse(top50_same_as_all, "top50 equals all because term has fewer than 50 genes", "")
    ) %>%
    select(reference, term, n_genes_all_original, n_genes_top50, top50_same_as_all, top50_note)
})) %>%
  mutate(reference = factor(reference, levels = reference_order)) %>%
  arrange(reference, term) %>%
  mutate(reference = as.character(reference))

write_csv_safe(term_count_note, file.path(tables_dir, "Auto_developmental_top50_equals_all_terms.csv"))

####################
# Method 1: exact original clusterProfiler enrichment logic
####################
safe_enrich_custom <- function(genes, custom_ref) {
  ####################
  # Fixed: removed minGSSize/maxGSSize overrides so enricher() uses
  # its default values (minGSSize=10, maxGSSize=500), exactly matching
  # enrichment_annotation.R's enricher() calls.
  ####################
  tryCatch(
    clusterProfiler::enricher(
      gene = genes,
      TERM2GENE = custom_ref$TERM2GENE,
      TERM2NAME = custom_ref$TERM2NAME,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    ),
    error = function(e) {
      warning(conditionMessage(e))
      NULL
    }
  )
}

cluster_enrich_has_custom_refs <- function(cluster_enrich, custom_refs) {
  if (is.null(cluster_enrich) || length(cluster_enrich) == 0) {
    return(FALSE)
  }
  all(mp_order %in% names(cluster_enrich)) &&
    all(vapply(names(custom_refs), function(ref_name) {
      all(vapply(cluster_enrich[mp_order], function(x) ref_name %in% names(x), logical(1)))
    }, logical(1)))
}

build_custom_cluster_enrich <- function(mp_genes, custom_refs, cache_path = NULL) {
  mp_assignments <- geneNMF.metaprograms$programs.clusters
  valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_genes)))
  mp_assignments <- mp_assignments[mp_assignments %in% valid_cluster_ids & !is.na(mp_assignments)]

  out <- lapply(names(mp_genes), function(mp_name) {
    genes <- mp_genes[[mp_name]]
    mp_id <- as.numeric(gsub("\\D", "", mp_name))
    members <- names(mp_assignments)[mp_assignments == mp_id]
    message("Processing MP custom enrichment: ", mp_name)

    res_custom_list <- lapply(names(custom_refs), function(ref_name) {
      message("  -> Running custom enrichment: ", ref_name)
      safe_enrich_custom(genes, custom_refs[[ref_name]])
    })
    names(res_custom_list) <- names(custom_refs)

    c(list(rep_prog = mp_name, members = members, genes = genes), res_custom_list)
  })
  names(out) <- names(mp_genes)
  if (!is.null(cache_path)) {
    saveRDS(out, cache_path)
  }
  out
}

cluster_enrich_original_path <- file.path(ref_out_dir, "cluster_enrich.rds")
cache_path_all <- file.path(intermediate_dir, "Auto_developmental_cluster_enrich_all_original_logic.rds")

if (file.exists(cache_path_all) && !force_rebuild) {
  cluster_enrich_all <- readRDS(cache_path_all)
} else {
  cat("Local custom enrichment cache missing or invalid; rebuilding all-gene custom enrichment.\n")
  cluster_enrich_all <- build_custom_cluster_enrich(
    mp_genes,
    custom_refs_all,
    cache_path = cache_path_all
  )
}
cluster_enrich_all <- cluster_enrich_all[mp_order]

cluster_enrich_top50_path <- file.path(intermediate_dir, "Auto_developmental_cluster_enrich_top50_original_logic.rds")
if (file.exists(cluster_enrich_top50_path) && !force_rebuild) {
  cluster_enrich_top50 <- readRDS(cluster_enrich_top50_path)
} else {
  cluster_enrich_top50 <- build_custom_cluster_enrich(
    mp_genes,
    custom_refs_top50,
    cache_path = cluster_enrich_top50_path
  )
}
cluster_enrich_top50 <- cluster_enrich_top50[mp_order]

extract_overlap_table <- function(cluster_enrich, custom_refs, mode, cap = 7) {
  bind_rows(lapply(names(custom_refs), function(ref_name) {
    terms_use <- as.character(custom_refs[[ref_name]]$TERM2NAME$term)
    df_list <- lapply(names(cluster_enrich), function(prog) {
      er <- cluster_enrich[[prog]][[ref_name]]
      if (is.null(er)) {
        return(NULL)
      }
      r <- tryCatch(er@result, error = function(e) NULL)
      if (is.null(r) || nrow(r) == 0) {
        return(NULL)
      }
      term <- if ("Description" %in% colnames(r)) r$Description else r$ID
      data.frame(
        Program = prog,
        Term = term,
        padj = r$p.adjust,
        Overlap = r$GeneRatio,
        stringsAsFactors = FALSE
      )
    })
    df <- bind_rows(df_list)
    if (is.null(df) || nrow(df) == 0) {
      df <- data.frame(Program = character(), Term = character(), padj = numeric(), Overlap = character())
    }
    full_grid <- expand.grid(Term = terms_use, Program = mp_order, stringsAsFactors = FALSE)
    full_grid %>%
      left_join(df, by = c("Term", "Program")) %>%
      transmute(
        gene_mode = mode,
        reference = ref_name,
        term = Term,
        mp = Program,
        p_adjust = padj,
        neglog10_padj = replace_na(pmin(-log10(padj), cap), 0),
        display = replace_na(Overlap, "")
      )
  }))
}

overlap_top50 <- extract_overlap_table(cluster_enrich_top50, custom_refs_top50, "top50")
overlap_all <- extract_overlap_table(cluster_enrich_all, custom_refs_all, "all")
write_csv_safe(overlap_top50, file.path(tables_dir, "Auto_developmental_overlap_top50.csv"))
write_csv_safe(overlap_all, file.path(tables_dir, "Auto_developmental_overlap_all.csv"))

overlap_validation <- overlap_all %>%
  select(reference, term, mp, neglog10_padj, display) %>%
  left_join(
    extract_overlap_table(cluster_enrich_all, custom_refs_all, "original") %>%
      select(reference, term, mp, neglog10_padj_original = neglog10_padj, display_original = display),
    by = c("reference", "term", "mp")
  ) %>%
  group_by(reference) %>%
  summarise(
    n_cells = n(),
    max_abs_diff_neglog10 = max(abs(neglog10_padj - neglog10_padj_original), na.rm = TRUE),
    n_display_mismatch = sum(display != display_original, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(overlap_validation, file.path(tables_dir, "Auto_developmental_validation_overlap_all_vs_original.csv"))

####################
# Method 2: exact original term extraction, UCell scoring, and correlation logic
####################
extract_all_terms_custom <- function(cluster_enrich, db_name, custom_refs) {
  terms_use <- as.character(custom_refs[[db_name]]$TERM2NAME$term)
  gene_lists <- list()
  for (mp in names(cluster_enrich)) {
    er <- cluster_enrich[[mp]][[db_name]]
    if (is.null(er)) next
    gs <- tryCatch(er@geneSets, error = function(e) NULL)
    if (is.null(gs) || length(gs) == 0) next
    for (term in terms_use) {
      if (term %in% names(gene_lists)) next
      if (term %in% names(gs)) {
        gene_lists[[term]] <- gs[[term]]
      }
    }
    if (length(gene_lists) == length(terms_use)) break
  }
  gene_lists[intersect(terms_use, names(gene_lists))]
}

####################
# Fixed: extract gene lists from @geneSets of the enricher() result,
# exactly like mp_database_correlation.R's extract_all_terms() does.
# This ensures the gene lists are universe-filtered (matching
# the enricher's internal processing) rather than raw TERM2GENE genes.
####################
extract_gene_lists_from_enrichment <- function(cluster_enrich, custom_refs, ref_name) {
  terms_use <- as.character(custom_refs[[ref_name]]$TERM2NAME$term)
  gene_lists <- list()
  for (mp in names(cluster_enrich)) {
    er <- cluster_enrich[[mp]][[ref_name]]
    if (is.null(er)) next
    gs <- tryCatch(er@geneSets, error = function(e) NULL)
    if (is.null(gs) || length(gs) == 0) next
    for (term in terms_use) {
      if (term %in% names(gene_lists)) next
      if (term %in% names(gs)) {
        gene_lists[[term]] <- gs[[term]]
      }
    }
    if (length(gene_lists) == length(terms_use)) break
  }
  gene_lists[intersect(terms_use, names(gene_lists))]
}

gene_lists_from_custom_ref <- function(custom_ref) {
  term2gene <- as.data.frame(custom_ref$TERM2GENE, stringsAsFactors = FALSE)
  term2name <- as.data.frame(custom_ref$TERM2NAME, stringsAsFactors = FALSE)
  colnames(term2gene)[1:2] <- c("term", "gene")
  colnames(term2name)[1:2] <- c("term", "name")
  term_order <- as.character(term2name$term)
  split_genes <- split(term2gene$gene, term2gene$term)
  split_genes[intersect(term_order, names(split_genes))]
}

db_gene_lists_all <- lapply(names(custom_refs_all), function(ref_name) {
  gene_lists_from_custom_ref(custom_refs_all[[ref_name]])
})
names(db_gene_lists_all) <- names(custom_refs_all)

db_gene_lists_top50 <- lapply(names(custom_refs_top50), function(ref_name) {
  gene_lists_from_custom_ref(custom_refs_top50[[ref_name]])
})
names(db_gene_lists_top50) <- names(custom_refs_top50)

flatten_gene_lists_original <- function(db_gene_lists) {
  all_gene_lists <- list()
  term_to_db <- character(0)
  safe_name_to_term <- character(0)

  for (db in names(db_gene_lists)) {
    for (term in names(db_gene_lists[[db]])) {
      safe_name_base <- paste0(db, "__", make.names(term))
      safe_name <- tail(make.unique(c(names(all_gene_lists), safe_name_base), sep = "_dup"), 1)
      all_gene_lists[[safe_name]] <- unique(db_gene_lists[[db]][[term]])
      term_to_db[safe_name] <- db
      safe_name_to_term[safe_name] <- term
    }
  }
  list(all_gene_lists = all_gene_lists, term_to_db = term_to_db, safe_name_to_term = safe_name_to_term)
}

score_reference_terms_original <- function(tmdata, db_gene_lists, mode) {
  flat <- flatten_gene_lists_original(db_gene_lists)
  all_gene_lists <- flat$all_gene_lists
  safe_names <- names(all_gene_lists)
  if (length(safe_names) == 0) {
    stop("No gene lists available for UCell scoring in mode: ", mode)
  }

  all_gene_lists <- lapply(all_gene_lists, function(genes) intersect(genes, rownames(tmdata)))
  keep <- lengths(all_gene_lists) > 0
  all_gene_lists <- all_gene_lists[keep]
  safe_names <- names(all_gene_lists)
  flat$term_to_db <- flat$term_to_db[safe_names]
  flat$safe_name_to_term <- flat$safe_name_to_term[safe_names]

  cache_path <- file.path(intermediate_dir, paste0("Auto_developmental_ref_ucell_scores_original_", mode, ".rds"))
  original_cache_path <- file.path(ref_out_dir, "UCell_ref_terms_v2_MP19.rds")
  common_cells <- colnames(tmdata)

  ref_scores <- data.frame(row.names = common_cells)

  if (mode == "all" && file.exists(original_cache_path) && !force_rebuild) {
    original_scores <- readRDS(original_cache_path)
    common_cols <- intersect(colnames(original_scores), safe_names)
    if (length(common_cols) > 0) {
      ref_scores <- cbind(ref_scores, original_scores[common_cells, common_cols, drop = FALSE])
    }
  }

  if (file.exists(cache_path) && !force_rebuild) {
    cached <- readRDS(cache_path)
    common_cols <- intersect(colnames(cached), safe_names)
    new_cols <- setdiff(common_cols, colnames(ref_scores))
    if (length(new_cols) > 0) {
      ref_scores <- cbind(ref_scores, cached[common_cells, new_cols, drop = FALSE])
    }
  }

  batch_size <- 20L
  ucell_max_rank <- max(5000, max(lengths(all_gene_lists)) + 100)
  missing_safe_names <- setdiff(safe_names, colnames(ref_scores))
  n_batches <- ceiling(length(missing_safe_names) / batch_size)

  for (b in seq_len(n_batches)) {
    idx_start <- (b - 1L) * batch_size + 1L
    idx_end <- min(b * batch_size, length(missing_safe_names))
    batch_names <- missing_safe_names[idx_start:idx_end]
    batch_lists <- all_gene_lists[batch_names]
    cat(sprintf("  %s UCell batch %d/%d: scoring %d gene sets\n", mode, b, n_batches, length(batch_lists)))

    tmdata <- UCell::AddModuleScore_UCell(
      tmdata,
      features = batch_lists,
      ncores = ucell_cores,
      name = "",
      maxRank = ucell_max_rank
    )
    score_cols <- intersect(batch_names, colnames(tmdata@meta.data))
    if (length(score_cols) > 0) {
      ref_scores <- cbind(ref_scores, tmdata@meta.data[common_cells, score_cols, drop = FALSE])
      tmdata@meta.data <- tmdata@meta.data[, !(colnames(tmdata@meta.data) %in% score_cols), drop = FALSE]
      saveRDS(ref_scores, cache_path)
    }
  }
  ref_scores <- ref_scores[, safe_names[safe_names %in% colnames(ref_scores)], drop = FALSE]
  saveRDS(ref_scores, cache_path)
  c(flat, list(ref_scores = ref_scores))
}

compute_cross_cor_original <- function(ref_score_mat, mod_mat, sample_ids, samples, min_cells = 5) {
  n_ref <- nrow(ref_score_mat)
  n_mp <- nrow(mod_mat)
  cor_array <- array(
    NA_real_,
    dim = c(n_ref, n_mp, length(samples)),
    dimnames = list(rownames(ref_score_mat), rownames(mod_mat), samples)
  )

  for (smp in samples) {
    idx <- which(sample_ids == smp)
    if (length(idx) > min_cells) {
      cor_array[, , smp] <- suppressWarnings(cor(
        t(ref_score_mat[, idx, drop = FALSE]),
        t(mod_mat[, idx, drop = FALSE]),
        method = "spearman"
      ))
    }
  }

  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(0, n_ref, n_mp, dimnames = list(rownames(ref_score_mat), rownames(mod_mat)))
  p_vals <- matrix(1, n_ref, n_mp, dimnames = list(rownames(ref_score_mat), rownames(mod_mat)))

  for (i in seq_len(n_ref)) {
    for (j in seq_len(n_mp)) {
      z_scores <- z_array[i, j, ]
      z_scores <- z_scores[!is.na(z_scores)]
      if (length(z_scores) > 1) {
        test_res <- t.test(z_scores)
        mean_rho[i, j] <- tanh(mean(z_scores))
        p_vals[i, j] <- test_res$p.value
      }
    }
  }
  list(mean_rho = mean_rho, p_vals = p_vals)
}

run_correlation_original <- function(scored, mode) {
  mod_mat <- t(as.matrix(mp_ucell))
  mod_mat <- mod_mat[mp_order[mp_order %in% rownames(mod_mat)], , drop = FALSE]
  sample_ids <- as.character(tmdata_all$orig.ident[colnames(mod_mat)])
  samples <- unique(sample_ids)

  results <- list()
  tables <- list()
  for (db in names(custom_refs_all)) {
    db_prefix <- paste0(db, "__")
    db_cols <- grep(paste0("^", db_prefix), colnames(scored$ref_scores), value = TRUE)
    if (length(db_cols) == 0) {
      warning("No scored terms found for correlation database: ", db)
      next
    }

    ref_mat <- t(as.matrix(scored$ref_scores[, db_cols, drop = FALSE]))
    rownames(ref_mat) <- unname(scored$safe_name_to_term[rownames(ref_mat)])
    result <- compute_cross_cor_original(ref_mat, mod_mat, sample_ids, samples)
    padj <- matrix(
      p.adjust(as.vector(result$p_vals), method = "BH"),
      nrow = nrow(result$p_vals),
      dimnames = dimnames(result$p_vals)
    )

    term_order <- names(if (mode == "all") db_gene_lists_all[[db]] else db_gene_lists_top50[[db]])
    row_order <- c(intersect(term_order, rownames(result$mean_rho)), setdiff(rownames(result$mean_rho), term_order))
    result$mean_rho <- result$mean_rho[row_order, mp_order, drop = FALSE]
    result$p_vals <- result$p_vals[row_order, mp_order, drop = FALSE]
    padj <- padj[row_order, mp_order, drop = FALSE]

    results[[db]] <- list(mean_rho = result$mean_rho, p_vals = result$p_vals, padj = padj)

    tables[[db]] <- expand.grid(term = rownames(result$mean_rho), mp = colnames(result$mean_rho), stringsAsFactors = FALSE) %>%
      mutate(
        gene_mode = mode,
        reference = db,
        mean_rho = as.vector(result$mean_rho),
        p_value = as.vector(result$p_vals),
        p_adjust = as.vector(padj),
        stars = stars_from_p(p_adjust),
        n_samples = length(samples)
      ) %>%
      select(gene_mode, reference, term, mp, mean_rho, p_value, p_adjust, stars, n_samples)
  }
  list(results = results, table = bind_rows(tables))
}

if (!replot_only) {
  tmdata_all <- readRDS(epi_path)
  mp_ucell <- readRDS(mp_ucell_path)
  common_cells <- intersect(colnames(tmdata_all), rownames(mp_ucell))
  tmdata_all <- tmdata_all[, common_cells]
  mp_ucell <- mp_ucell[common_cells, mp_order, drop = FALSE]

  scored_top50 <- score_reference_terms_original(tmdata_all, db_gene_lists_top50, "top50")
  cor_top50_obj <- run_correlation_original(scored_top50, "top50")
  cor_top50 <- cor_top50_obj$table
  saveRDS(cor_top50_obj$results, file.path(intermediate_dir, "Auto_developmental_correlation_results_top50_original_logic.rds"))
  write_csv_safe(cor_top50, file.path(tables_dir, "Auto_developmental_expression_correlation_top50.csv"))

  scored_all <- score_reference_terms_original(tmdata_all, db_gene_lists_all, "all")
  cor_all_obj <- run_correlation_original(scored_all, "all")
  cor_all <- cor_all_obj$table
  saveRDS(cor_all_obj$results, file.path(intermediate_dir, "Auto_developmental_correlation_results_all_original_logic.rds"))
  write_csv_safe(cor_all, file.path(tables_dir, "Auto_developmental_expression_correlation_all.csv"))
} else {
  cor_top50 <- read.csv(file.path(tables_dir, "Auto_developmental_expression_correlation_top50.csv"), stringsAsFactors = FALSE)
  cor_all <- read.csv(file.path(tables_dir, "Auto_developmental_expression_correlation_all.csv"), stringsAsFactors = FALSE)
  cor_top50_obj <- list(results = readRDS(file.path(intermediate_dir, "Auto_developmental_correlation_results_top50_original_logic.rds")))
  cor_all_obj <- list(results = readRDS(file.path(intermediate_dir, "Auto_developmental_correlation_results_all_original_logic.rds")))
}

original_cor_path <- file.path(ref_out_dir, "Auto_MP_correlation_results_v2_MP19.rds")
if (file.exists(original_cor_path)) {
  original_cor <- readRDS(original_cor_path)
  cor_validation <- bind_rows(lapply(intersect(names(cor_all_obj$results), names(original_cor)), function(db) {
    new_mat <- cor_all_obj$results[[db]]$mean_rho
    old_mat <- original_cor[[db]]$mean_rho
    common_terms <- intersect(rownames(new_mat), rownames(old_mat))
    common_mps <- intersect(colnames(new_mat), colnames(old_mat))
    data.frame(
      reference = db,
      n_terms_compared = length(common_terms),
      n_mps_compared = length(common_mps),
      max_abs_diff_mean_rho = max(abs(new_mat[common_terms, common_mps, drop = FALSE] - old_mat[common_terms, common_mps, drop = FALSE]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
} else {
  cor_validation <- data.frame(reference = names(cor_all_obj$results), n_terms_compared = NA_integer_, n_mps_compared = NA_integer_, max_abs_diff_mean_rho = NA_real_)
}
write_csv_safe(cor_validation, file.path(tables_dir, "Auto_developmental_validation_correlation_all_vs_original.csv"))

####################
# Method 3: exact external MP UCell scoring helper logic, for annotated data available locally
####################
oesophagus_labels <- c("Quiescent basal cell", "Basal cell (cycling)", "Suprabasal", "Apical cell")
oesophagus_term_map <- setNames(paste0(gsub(" ", "_", oesophagus_labels), "_Oesophagus..Adult_Epithelium"), oesophagus_labels)
oesophagus_order <- unname(oesophagus_term_map)

stomach_rename_map <- c(
  "GKN+F" = "GKN+F_PylG_Stomach..Adult_Epithelium",
  "ADH1+GKN1-F" = "ADH1+GKN1-F_PylG_Stomach..Adult_Epithelium",
  "PG/Neck1" = "PG/Neck1_PylG_Stomach..Adult_Epithelium",
  "PG/Neck2" = "PG/Neck2_PylG_Stomach..Adult_Epithelium",
  "NE1" = "NE1_PylG_Stomach..Adult_Epithelium",
  "PC" = "PC_PylG_Stomach..Adult_Epithelium",
  "Chief" = "Chief_FundG_Stomach..Adult_Epithelium",
  "Pr_epi" = "Pr_epi_FundG_Stomach..Adult_Epithelium",
  "NE2" = "NE2_PylG/IntestMeta_Stomach..Adult_Epithelium",
  "Ent" = "Ent_IntestMeta_Stomach..Adult_Epithelium",
  "Gob" = "Gob_IntestMeta_Stomach..Adult_Epithelium"
)
stomach_order <- unname(stomach_rename_map)

barrett_term_map <- c(
  "Basal" = "Basal_Normal_Esophagus..Barretts_Oesophagus",
  "Suprabasal" = "Suprabasal_Normal_Esophagus..Barretts_Oesophagus",
  "Suprabasal_Dividing" = "Suprabasal_Dividing_Normal_Esophagus..Barretts_Oesophagus",
  "Intermediate" = "Intermediate_Normal_Esophagus..Barretts_Oesophagus",
  "Superficial" = "Superficial_Normal_Esophagus..Barretts_Oesophagus",
  "Undifferentiated" = "Undifferentiated_Normal_Gastric..Barretts_Oesophagus",
  "Undifferentiated_Dividing" = "Undifferentiated_Dividing_Normal_Gastric..Barretts_Oesophagus",
  "Foveolar_Intermediate" = "Foveolar_Intermediate_Normal_Gastric..Barretts_Oesophagus",
  "Foveolar_differentiated" = "Foveolar_differentiated_Normal_Gastric..Barretts_Oesophagus",
  "Chief" = "Chief_Normal_Gastric..Barretts_Oesophagus",
  "Parietal" = "Parietal_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_GHRL" = "Endocrine_GHRL_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_CHGA" = "Endocrine_CHGA_Normal_Gastric..Barretts_Oesophagus",
  "Endocrine_NEUROD1" = "Endocrine_NEUROD1_Normal_Gastric..Barretts_Oesophagus",
  "Columnar_Undifferentiated" = "Columnar_Undifferentiated_Barretts_Esophagus..Barretts_Oesophagus",
  "Columnar_Undifferentiated_Dividing" = "Columnar_Dividing_Barretts_Esophagus..Barretts_Oesophagus",
  "Endocrine_NEUROG3" = "Endocrine_NEUROG3_Barretts_Esophagus..Barretts_Oesophagus",
  "C1" = "Columnar_Intermediate_Barretts_Esophagus..Barretts_Oesophagus",
  "C2" = "Columnar_differentiated_Barretts_Esophagus..Barretts_Oesophagus",
  "Goblet" = "Goblet_Barretts_Esophagus..Barretts_Oesophagus",
  "Duct_Intercalating" = "Duct_Intercalating_Submucosal_Glands..Barretts_Oesophagus",
  "Oncocytes_MUC5B_Low" = "Oncocytes_Submucosal_Glands..Barretts_Oesophagus",
  "Mucous_MUC5B_High" = "Mucous_Submucosal_Glands..Barretts_Oesophagus"
)
barrett_order <- unname(barrett_term_map)

download_dir <- file.path(out_dir, "downloads")
organogenesis_required_files <- file.path(download_dir, "organogenesis", c(
  "GSE157329_cell_annotate.txt.gz",
  "GSE157329_gene_annotate.txt.gz",
  "GSE157329_raw_counts.mtx.gz"
))
early_required_files <- file.path(download_dir, "early_embryogenesis", "psd.R3.6.em.seurat.ob.rds")
normal_development_required_files <- file.path(download_dir, "normal_development", c(
  "Stomach_gene_count.RDS",
  "df_cell.RDS",
  "df_gene.RDS"
))
barretts_required_file <- file.path(download_dir, "barretts", "alldatahighquality.rds")
oesophagus_download_dir <- file.path(download_dir, "adult_oesophagus")

find_one_file <- function(root, pattern) {
  hits <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

oesophagus_meta_path <- find_one_file(oesophagus_download_dir, "^EoE_meta\\.txt$")
oesophagus_cell_path <- find_one_file(oesophagus_download_dir, "^EoE_cell\\.tsv$")
oesophagus_gene_path <- find_one_file(oesophagus_download_dir, "^EoE_gene\\.tsv$")
oesophagus_mtx_path <- find_one_file(oesophagus_download_dir, "^EoE\\.mtx(\\.gz)?$")
oesophagus_required_files <- c(oesophagus_meta_path, oesophagus_cell_path, oesophagus_gene_path, oesophagus_mtx_path)

####################
# Annotated expression download provenance:
# - adult oesophagus: Broad SCP1242 temporary authorized bulk-download config generated with:
#   curl "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1242&auth_code=vkAMi2rX&directory=all&context=study" -o cfg.txt
#   curl -K cfg.txt
#   The generated cfg.txt is retained under downloads/adult_oesophagus for provenance.
# - normal development stomach: direct Descartes/Fred Hutch RDS downloads from the atlas
#   page-data JSON, using only Stomach_gene_count.RDS plus df_cell.RDS/df_gene.RDS.
# - Barretts: direct high-quality combined RDS download from esophaguscancercellatlas.org.
####################
reference_expression_sources <- data.frame(
  reference = c("Early_Embryogenesis", "Organogenesis_major", "Organogenesis_sub", "Normal_Development_long", "Normal_Development_short", "Adult_Epithelium", "Adult_Epithelium", "Barretts_Oesophagus"),
  dataset_source = c("Early_Embryogenesis_integrated", "Organogenesis_GSE157329_major", "Organogenesis_GSE157329_sub", "Normal_Development_Stomach_long", "Normal_Development_Stomach_short", "Adult_Oesophagus", "Adult_Stomach", "Barretts_HighQuality"),
  dataset = c("early_embryogenesis", "organogenesis_4_6wk", "organogenesis_4_6wk", "normal_development_stomach_10_18wk", "normal_development_stomach_10_18wk", "adult_oesophagus", "adult_stomach", "barretts_high_quality"),
  paper_url = c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11725501/",
    "https://www.nature.com/articles/s41556-023-01108-w",
    "https://www.nature.com/articles/s41556-023-01108-w",
    "https://www.science.org/doi/10.1126/science.aba7721",
    "https://www.science.org/doi/10.1126/science.aba7721",
    "https://www.nature.com/articles/s41467-024-47647-0",
    "https://www.sciencedirect.com/science/article/pii/S2211124723012482",
    "https://www.science.org/doi/10.1126/science.abd1449"
  ),
  download_or_repository = c(
    "paper/data portal downloaded Seurat RDS to downloads/early_embryogenesis",
    "GEO GSE157329 and VisCello https://heoa.shinyapps.io/base/",
    "GEO GSE157329 and VisCello https://heoa.shinyapps.io/base/",
    "Descartes/Fred Hutch direct RDS: Stomach_gene_count.RDS, df_cell.RDS, df_gene.RDS",
    "Descartes/Fred Hutch direct RDS: Stomach_gene_count.RDS, df_cell.RDS, df_gene.RDS",
    "Broad Single Cell Portal SCP1242 temporary authorized bulk-download curl config",
    "local annotated Seurat object found under developmental reference directory",
    "Esophagus Cancer Cell Atlas direct high-quality combined RDS download"
  ),
  download_command = c(
    "downloaded from paper dataset portal as psd.R3.6.em.seurat.ob.rds",
    "downloaded GEO files GSE157329_cell_annotate.txt.gz, GSE157329_gene_annotate.txt.gz, GSE157329_raw_counts.mtx.gz",
    "downloaded GEO files GSE157329_cell_annotate.txt.gz, GSE157329_gene_annotate.txt.gz, GSE157329_raw_counts.mtx.gz",
    "wget -c https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/Stomach_gene_count.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS",
    "wget -c https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/Stomach_gene_count.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS",
    "curl Broad SCP1242 generate_curl_config URL with temporary auth_code to cfg.txt; curl -K cfg.txt",
    "local file /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds; no download in this script",
    "wget/curl https://cellgeni.cog.sanger.ac.uk/esophagus-cancer/alldatahighquality.rds"
  ),
  access_requirement = c(
    "direct/local cached RDS; original paper data portal",
    "direct public GEO/VisCello files",
    "direct public GEO/VisCello files",
    "direct public atlas download; no login",
    "direct public atlas download; no login",
    "temporary Broad auth_code supplied by user; no interactive login during this run",
    "local project file; no login",
    "direct public atlas download; no login"
  ),
  local_path = c(
    paste(early_required_files, collapse = ";"),
    paste(organogenesis_required_files, collapse = ";"),
    paste(organogenesis_required_files, collapse = ";"),
    paste(normal_development_required_files, collapse = ";"),
    paste(normal_development_required_files, collapse = ";"),
    paste(oesophagus_required_files, collapse = ";"),
    file.path(developmental_dir, "Stomach.rds"),
    barretts_required_file
  ),
  stringsAsFactors = FALSE
)
reference_expression_sources$available <- c(
  all(file.exists(early_required_files)),
  all(file.exists(organogenesis_required_files)),
  all(file.exists(organogenesis_required_files)),
  all(file.exists(normal_development_required_files)),
  all(file.exists(normal_development_required_files)),
  all(!is.na(oesophagus_required_files)) && all(file.exists(oesophagus_required_files)),
  file.exists(file.path(developmental_dir, "Stomach.rds")),
  file.exists(barretts_required_file)
)
write_csv_safe(reference_expression_sources, file.path(tables_dir, "Auto_developmental_reference_expression_availability.csv"))

sample_meta_by_term_original <- function(meta_df, term_levels, max_cells, seed) {
  meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)
  meta_df$term <- factor(meta_df$term, levels = term_levels)
  meta_df <- meta_df[!is.na(meta_df$term), , drop = FALSE]
  set.seed(seed)
  sampled_list <- lapply(split(meta_df, meta_df$term, drop = TRUE), function(df) {
    if (nrow(df) <= max_cells) return(df)
    df[sample(seq_len(nrow(df)), max_cells, replace = FALSE), , drop = FALSE]
  })
  sampled_meta <- bind_rows(sampled_list)
  sampled_meta$term <- factor(sampled_meta$term, levels = term_levels)
  sampled_meta <- sampled_meta[order(sampled_meta$term, sampled_meta$cell), , drop = FALSE]
  sampled_meta$term <- as.character(sampled_meta$term)
  sampled_meta
}

score_sampled_counts_original <- function(counts_mat, sampled_meta, full_meta, dataset_source, seed) {
  sampled_meta <- as.data.frame(sampled_meta, stringsAsFactors = FALSE)
  full_meta <- as.data.frame(full_meta, stringsAsFactors = FALSE)
  if (anyDuplicated(rownames(counts_mat))) {
    rownames(counts_mat) <- make.unique(rownames(counts_mat))
  }
  counts_mat <- counts_mat[, sampled_meta$cell, drop = FALSE]
  sampled_meta <- sampled_meta[match(colnames(counts_mat), sampled_meta$cell), , drop = FALSE]
  rownames(sampled_meta) <- sampled_meta$cell

  usable_mp_genes <- lapply(mp_genes, function(genes) intersect(genes, rownames(counts_mat)))
  empty_mps <- names(usable_mp_genes)[lengths(usable_mp_genes) == 0]
  if (length(empty_mps) > 0) {
    stop("No overlapping genes for ", dataset_source, " in: ", paste(empty_mps, collapse = ", "))
  }

  set.seed(seed)
  obj <- CreateSeuratObject(counts = counts_mat, meta.data = sampled_meta)
  obj <- AddModuleScore_UCell(
    obj,
    features = usable_mp_genes,
    slot = "counts",
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
    ncores = 1,
    name = "_UCell"
  )

  ucell_score_cols <- paste0(mp_order, "_UCell")
  ucell_scores <- obj@meta.data[, ucell_score_cols, drop = FALSE]
  colnames(ucell_scores) <- mp_order
  score_df <- cbind(sampled_meta[, c("cell", "term", "dataset_source"), drop = FALSE], as.data.frame(ucell_scores))
  mean_scores <- score_df %>% group_by(term) %>% summarise(across(all_of(mp_order), mean), .groups = "drop")
  total_counts <- full_meta %>% dplyr::count(term, name = "n_cells_total")
  scored_counts <- sampled_meta %>% dplyr::count(term, name = "n_cells_scored")

  full_join(total_counts, scored_counts, by = "term") %>%
    left_join(mean_scores, by = "term") %>%
    mutate(
      dataset_source = dataset_source,
      n_cells_total = ifelse(is.na(n_cells_total), 0L, n_cells_total),
      n_cells_scored = ifelse(is.na(n_cells_scored), 0L, n_cells_scored)
    ) %>%
    select(dataset_source, term, n_cells_total, n_cells_scored, all_of(mp_order))
}

normalise_annotation_label <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

term_base_label <- function(term) {
  gsub("_", " ", sub("_Stomach\\.\\..*$", "", as.character(term)))
}

make_normal_development_stomach_meta <- function(df_cell, count_cells, ref_name, dataset_source) {
  stomach_terms <- grep("_Stomach\\.\\.", as.character(custom_refs_all[[ref_name]]$TERM2NAME$term), value = TRUE)
  label_to_term <- setNames(stomach_terms, normalise_annotation_label(term_base_label(stomach_terms)))
  df_cell %>%
    filter(sample %in% count_cells, Organ == "Stomach", !is.na(Main_cluster_name)) %>%
    transmute(
      cell = as.character(sample),
      term = unname(label_to_term[normalise_annotation_label(Main_cluster_name)]),
      dataset_source = dataset_source
    ) %>%
    filter(!is.na(term))
}

subset_oesophagus_counts_original <- function(mtx_path, cell_path, gene_path, sampled_meta, cache_prefix) {
  cache_path <- file.path(intermediate_dir, paste0(cache_prefix, "_counts_subset.rds"))
  if (file.exists(cache_path) && !force_rebuild) {
    return(readRDS(cache_path))
  }

  all_cells <- data.table::fread(cell_path, header = FALSE, col.names = "cell")$cell
  all_genes <- data.table::fread(gene_path, header = FALSE, col.names = "gene")$gene
  keep_idx <- match(sampled_meta$cell, all_cells)
  if (any(is.na(keep_idx))) {
    stop("Adult oesophagus sampled cells missing from downloaded EoE_cell.tsv: ", sum(is.na(keep_idx)))
  }
  index_map <- data.frame(old_index = keep_idx, new_index = seq_along(keep_idx))
  map_path <- file.path(intermediate_dir, paste0(cache_prefix, "_cell_index_map.tsv"))
  triplet_path <- file.path(intermediate_dir, paste0(cache_prefix, "_filtered_triplets.tsv"))
  filtered_mtx_path <- file.path(intermediate_dir, paste0(cache_prefix, "_filtered.mtx"))
  write.table(index_map, map_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  awk_program <- "FNR==NR {map[$1]=$2; next} /^%/ {next} !seen_dims {seen_dims=1; next} NF==3 {if ($2 in map) print $1, map[$2], $3}"
  if (grepl("\\.gz$", mtx_path)) {
    cmd <- sprintf("zcat %s | awk '%s' %s - > %s", shQuote(mtx_path), awk_program, shQuote(map_path), shQuote(triplet_path))
  } else {
    cmd <- sprintf("awk '%s' %s %s > %s", awk_program, shQuote(map_path), shQuote(mtx_path), shQuote(triplet_path))
  }
  status <- system(cmd)
  if (!is.null(status) && status != 0) {
    stop("Failed to subset adult oesophagus MatrixMarket file with awk.")
  }
  n_entries <- as.integer(scan(text = system2("wc", c("-l", triplet_path), stdout = TRUE), what = character(), quiet = TRUE)[1])
  writeLines(c(
    "%%MatrixMarket matrix coordinate integer general",
    paste(length(all_genes), length(keep_idx), n_entries)
  ), filtered_mtx_path)
  file.append(filtered_mtx_path, triplet_path)
  counts <- Matrix::readMM(filtered_mtx_path)
  counts <- as(counts, "CsparseMatrix")
  rownames(counts) <- all_genes
  colnames(counts) <- sampled_meta$cell
  saveRDS(counts, cache_path)
  counts
}

score_external_available <- function() {
  external_cache <- file.path(intermediate_dir, "Auto_developmental_external_mp_ucell_scores_original_logic_v2.rds")
  if (file.exists(external_cache) && !force_rebuild) {
    return(readRDS(external_cache))
  }

  summaries <- list()
  if (all(file.exists(early_required_files))) {
    cat("Scoring downloaded early embryogenesis Seurat object with original external UCell logic.\n")
    early_obj <- readRDS(early_required_files)
    early_meta <- early_obj@meta.data
    early_meta$cell <- rownames(early_meta)
    early_label <- as.character(early_meta$rename_EML)
    early_label[early_label == "2-4 cell"] <- "Z4cell"
    early_terms <- as.character(custom_refs_all$Early_Embryogenesis$TERM2NAME$term)
    early_meta <- early_meta %>%
      mutate(term = paste0(gsub(" ", "_", early_label), "..Early_Embryogenesis")) %>%
      transmute(cell = cell, term = term, dataset_source = "Early_Embryogenesis_integrated") %>%
      filter(term %in% early_terms)
    if (nrow(early_meta) > 0) {
      early_sampled <- sample_meta_by_term_original(
        early_meta,
        term_levels = early_terms,
        max_cells = max_cells_per_type,
        seed = seed_base + 3L
      )
      early_counts <- early_obj@assays[[early_obj@active.assay]]@counts
      summaries[["Early_Embryogenesis_integrated"]] <- score_sampled_counts_original(
        counts_mat = early_counts,
        sampled_meta = early_sampled,
        full_meta = early_meta,
        dataset_source = "Early_Embryogenesis_integrated",
        seed = seed_base + 13L
      )
      rm(early_counts)
      gc()
    }
    rm(early_obj)
    gc()
  }

  organogenesis_counts_path <- organogenesis_required_files[3]
  if (all(file.exists(organogenesis_required_files))) {
    cat("Scoring downloaded organogenesis GSE157329 annotated raw counts with original external UCell logic.\n")
    organogenesis_meta <- data.table::fread(organogenesis_required_files[1])
    organogenesis_genes <- data.table::fread(organogenesis_required_files[2])

    major_terms <- as.character(custom_refs_all$Organogenesis_major$TERM2NAME$term)
    sub_terms <- as.character(custom_refs_all$Organogenesis_sub$TERM2NAME$term)
    organogenesis_meta <- organogenesis_meta %>%
      mutate(
        cell = as.character(cell_id),
        major_term = paste0(gsub(" ", "_", as.character(`developmental system`)), "..Organogenesis_major"),
        sub_base = as.character(final_annotation),
        sub_term = ifelse(
          grepl("epithelium", sub_base, ignore.case = TRUE),
          paste0(gsub(" ", "_", sub_base), "..Organogenesis_sub"),
          paste0(gsub(" ", "_", sub_base), "_Endoderm..Organogenesis_sub")
        )
      )

    organogenesis_major_meta <- organogenesis_meta %>%
      transmute(cell = cell, term = major_term, dataset_source = "Organogenesis_GSE157329_major") %>%
      filter(term %in% major_terms)
    organogenesis_sub_meta <- organogenesis_meta %>%
      transmute(cell = cell, term = sub_term, dataset_source = "Organogenesis_GSE157329_sub") %>%
      filter(term %in% sub_terms)
    if (nrow(organogenesis_major_meta) > 0 || nrow(organogenesis_sub_meta) > 0) {
      organogenesis_counts <- Matrix::readMM(gzfile(organogenesis_counts_path))
      organogenesis_counts <- as(organogenesis_counts, "CsparseMatrix")
      if (nrow(organogenesis_counts) == nrow(organogenesis_genes)) {
        rownames(organogenesis_counts) <- organogenesis_genes$gene_short_name
        colnames(organogenesis_counts) <- organogenesis_meta$cell
      } else if (ncol(organogenesis_counts) == nrow(organogenesis_genes)) {
        organogenesis_counts <- t(organogenesis_counts)
        rownames(organogenesis_counts) <- organogenesis_genes$gene_short_name
        colnames(organogenesis_counts) <- organogenesis_meta$cell
      } else {
        stop("Organogenesis counts dimensions do not match downloaded gene annotation.")
      }

      if (nrow(organogenesis_major_meta) > 0) {
        organogenesis_major_sampled <- sample_meta_by_term_original(
          organogenesis_major_meta,
          term_levels = major_terms,
          max_cells = max_cells_per_type,
          seed = seed_base + 4L
        )
        summaries[["Organogenesis_GSE157329_major"]] <- score_sampled_counts_original(
          counts_mat = organogenesis_counts,
          sampled_meta = organogenesis_major_sampled,
          full_meta = organogenesis_major_meta,
          dataset_source = "Organogenesis_GSE157329_major",
          seed = seed_base + 14L
        )
      }

      if (nrow(organogenesis_sub_meta) > 0) {
        organogenesis_sub_sampled <- sample_meta_by_term_original(
          organogenesis_sub_meta,
          term_levels = sub_terms,
          max_cells = max_cells_per_type,
          seed = seed_base + 5L
        )
        summaries[["Organogenesis_GSE157329_sub"]] <- score_sampled_counts_original(
          counts_mat = organogenesis_counts,
          sampled_meta = organogenesis_sub_sampled,
          full_meta = organogenesis_sub_meta,
          dataset_source = "Organogenesis_GSE157329_sub",
          seed = seed_base + 15L
        )
      }
      rm(organogenesis_counts)
      gc()
    }
  }

  if (all(file.exists(normal_development_required_files))) {
    cat("Scoring downloaded Descartes normal-development stomach counts with original external UCell logic.\n")
    normal_counts <- readRDS(normal_development_required_files[1])
    normal_counts <- as(normal_counts, "CsparseMatrix")
    normal_cell <- readRDS(normal_development_required_files[2])
    normal_gene <- readRDS(normal_development_required_files[3])
    normal_symbols <- normal_gene$gene_short_name[match(rownames(normal_counts), normal_gene$gene_id)]
    normal_symbols[is.na(normal_symbols) | normal_symbols == ""] <- rownames(normal_counts)[is.na(normal_symbols) | normal_symbols == ""]
    rownames(normal_counts) <- normal_symbols

    normal_long_terms <- grep("_Stomach\\.\\.", as.character(custom_refs_all$Normal_Development_long$TERM2NAME$term), value = TRUE)
    normal_short_terms <- grep("_Stomach\\.\\.", as.character(custom_refs_all$Normal_Development_short$TERM2NAME$term), value = TRUE)
    normal_long_meta <- make_normal_development_stomach_meta(
      df_cell = normal_cell,
      count_cells = colnames(normal_counts),
      ref_name = "Normal_Development_long",
      dataset_source = "Normal_Development_Stomach_long"
    )
    normal_short_meta <- make_normal_development_stomach_meta(
      df_cell = normal_cell,
      count_cells = colnames(normal_counts),
      ref_name = "Normal_Development_short",
      dataset_source = "Normal_Development_Stomach_short"
    )
    if (nrow(normal_long_meta) > 0) {
      normal_long_sampled <- sample_meta_by_term_original(
        normal_long_meta,
        term_levels = normal_long_terms,
        max_cells = max_cells_per_type,
        seed = seed_base + 6L
      )
      summaries[["Normal_Development_Stomach_long"]] <- score_sampled_counts_original(
        counts_mat = normal_counts,
        sampled_meta = normal_long_sampled,
        full_meta = normal_long_meta,
        dataset_source = "Normal_Development_Stomach_long",
        seed = seed_base + 16L
      )
    }
    if (nrow(normal_short_meta) > 0) {
      normal_short_sampled <- sample_meta_by_term_original(
        normal_short_meta,
        term_levels = normal_short_terms,
        max_cells = max_cells_per_type,
        seed = seed_base + 7L
      )
      summaries[["Normal_Development_Stomach_short"]] <- score_sampled_counts_original(
        counts_mat = normal_counts,
        sampled_meta = normal_short_sampled,
        full_meta = normal_short_meta,
        dataset_source = "Normal_Development_Stomach_short",
        seed = seed_base + 17L
      )
    }
    rm(normal_counts, normal_cell, normal_gene)
    gc()
  }

  if (all(!is.na(oesophagus_required_files)) && all(file.exists(oesophagus_required_files))) {
    cat("Scoring downloaded adult oesophagus SCP1242 MatrixMarket counts with original external UCell logic.\n")
    oesophagus_meta <- data.table::fread(oesophagus_meta_path)
    oesophagus_meta <- oesophagus_meta %>%
      filter(NAME != "TYPE") %>%
      transmute(
        cell = as.character(NAME),
        term = unname(oesophagus_term_map[as.character(cell_type_anno)]),
        dataset_source = "Adult_Oesophagus"
      ) %>%
      filter(!is.na(term))
    if (nrow(oesophagus_meta) > 0) {
      oesophagus_sampled <- sample_meta_by_term_original(
        oesophagus_meta,
        term_levels = oesophagus_order,
        max_cells = max_cells_per_type,
        seed = seed_base + 8L
      )
      oesophagus_counts <- subset_oesophagus_counts_original(
        mtx_path = oesophagus_mtx_path,
        cell_path = oesophagus_cell_path,
        gene_path = oesophagus_gene_path,
        sampled_meta = oesophagus_sampled,
        cache_prefix = paste0("Auto_adult_oesophagus_original_logic_max", max_cells_per_type)
      )
      summaries[["Adult_Oesophagus"]] <- score_sampled_counts_original(
        counts_mat = oesophagus_counts,
        sampled_meta = oesophagus_sampled,
        full_meta = oesophagus_meta,
        dataset_source = "Adult_Oesophagus",
        seed = seed_base + 18L
      )
      rm(oesophagus_counts)
      gc()
    }
  }

  stomach_path <- file.path(developmental_dir, "Stomach.rds")
  if (file.exists(stomach_path)) {
    cat("Scoring adult stomach annotated Seurat object with original external UCell logic.\n")
    adult_stomach <- readRDS(stomach_path)
    stomach_meta <- adult_stomach@meta.data
    stomach_meta$cell <- Cells(adult_stomach)
    stomach_meta <- stomach_meta %>%
      filter(major_clusters == "epi") %>%
      transmute(
        cell = cell,
        term = unname(stomach_rename_map[as.character(subcluster.v2)]),
        dataset_source = "Adult_Stomach"
      ) %>%
      filter(!is.na(term))
    stomach_sampled <- sample_meta_by_term_original(stomach_meta, stomach_order, max_cells_per_type, seed_base + 2L)
    stomach_counts <- GetAssayData(adult_stomach, layer = "counts")
    summaries[["Adult_Stomach"]] <- score_sampled_counts_original(
      counts_mat = stomach_counts,
      sampled_meta = stomach_sampled,
      full_meta = stomach_meta,
      dataset_source = "Adult_Stomach",
      seed = seed_base + 12L
    )
    rm(adult_stomach, stomach_counts)
    gc()
  }

  if (file.exists(barretts_required_file)) {
    cat("Scoring downloaded Barretts high-quality SingleCellExperiment with original external UCell logic.\n")
    barretts_sce <- readRDS(barretts_required_file)
    barretts_assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(barretts_sce)) "counts" else SummarizedExperiment::assayNames(barretts_sce)[1]
    barretts_counts <- SummarizedExperiment::assay(barretts_sce, barretts_assay_name)
    barretts_counts <- as(barretts_counts, "CsparseMatrix")
    barretts_meta <- as.data.frame(SummarizedExperiment::colData(barretts_sce), stringsAsFactors = FALSE)
    barretts_cell_ids <- colnames(barretts_sce)
    if (is.null(barretts_cell_ids) || length(barretts_cell_ids) != ncol(barretts_counts)) {
      barretts_cell_ids <- rownames(barretts_meta)
    }
    if (is.null(barretts_cell_ids) || length(barretts_cell_ids) != ncol(barretts_counts)) {
      barretts_cell_ids <- paste0("Barretts_cell_", seq_len(ncol(barretts_counts)))
    }
    colnames(barretts_counts) <- barretts_cell_ids
    barretts_meta$cell_id_for_scoring <- barretts_cell_ids
    barretts_meta <- barretts_meta %>%
      transmute(
        cell = cell_id_for_scoring,
        term = unname(barrett_term_map[as.character(cell_type_secondary)]),
        dataset_source = "Barretts_HighQuality"
      ) %>%
      filter(!is.na(term))
    if (nrow(barretts_meta) > 0) {
      barretts_sampled <- sample_meta_by_term_original(
        barretts_meta,
        term_levels = barrett_order,
        max_cells = max_cells_per_type,
        seed = seed_base + 9L
      )
      summaries[["Barretts_HighQuality"]] <- score_sampled_counts_original(
        counts_mat = barretts_counts,
        sampled_meta = barretts_sampled,
        full_meta = barretts_meta,
        dataset_source = "Barretts_HighQuality",
        seed = seed_base + 19L
      )
    }
    rm(barretts_sce, barretts_counts)
    gc()
  }

  out <- bind_rows(summaries)
  saveRDS(out, external_cache)
  out
}

external_summary <- score_external_available()
write_csv_safe(external_summary, file.path(tables_dir, "Auto_developmental_external_mp_ucell_summary.csv"))
write_csv_safe(external_summary, file.path(summary_dir, "Auto_developmental_external_mp_ucell_summary.csv"))

expected_external_terms <- bind_rows(
  data.frame(reference = "Early_Embryogenesis", dataset_source = "Early_Embryogenesis_integrated", term = as.character(custom_refs_all$Early_Embryogenesis$TERM2NAME$term), stringsAsFactors = FALSE),
  data.frame(reference = "Organogenesis_major", dataset_source = "Organogenesis_GSE157329_major", term = as.character(custom_refs_all$Organogenesis_major$TERM2NAME$term), stringsAsFactors = FALSE),
  data.frame(reference = "Organogenesis_sub", dataset_source = "Organogenesis_GSE157329_sub", term = as.character(custom_refs_all$Organogenesis_sub$TERM2NAME$term), stringsAsFactors = FALSE),
  data.frame(reference = "Normal_Development_long", dataset_source = "Normal_Development_Stomach_long", term = grep("_Stomach\\.\\.", as.character(custom_refs_all$Normal_Development_long$TERM2NAME$term), value = TRUE), stringsAsFactors = FALSE),
  data.frame(reference = "Normal_Development_short", dataset_source = "Normal_Development_Stomach_short", term = grep("_Stomach\\.\\.", as.character(custom_refs_all$Normal_Development_short$TERM2NAME$term), value = TRUE), stringsAsFactors = FALSE),
  data.frame(reference = "Adult_Epithelium", dataset_source = "Adult_Oesophagus", term = oesophagus_order, stringsAsFactors = FALSE),
  data.frame(reference = "Adult_Epithelium", dataset_source = "Adult_Stomach", term = stomach_order, stringsAsFactors = FALSE),
  data.frame(reference = "Barretts_Oesophagus", dataset_source = "Barretts_HighQuality", term = barrett_order, stringsAsFactors = FALSE)
)
external_counts_for_coverage <- external_summary %>%
  select(dataset_source, term, n_cells_total, n_cells_scored)
external_coverage <- expected_external_terms %>%
  left_join(reference_expression_sources %>% select(dataset_source, available, access_requirement, download_or_repository, local_path), by = "dataset_source") %>%
  left_join(external_counts_for_coverage, by = c("dataset_source", "term")) %>%
  mutate(
    present_in_annotation = !is.na(n_cells_total) & n_cells_total > 0,
    scored = !is.na(n_cells_scored) & n_cells_scored > 0,
    n_cells_total = ifelse(is.na(n_cells_total), 0L, n_cells_total),
    n_cells_scored = ifelse(is.na(n_cells_scored), 0L, n_cells_scored)
  ) %>%
  arrange(factor(reference, levels = reference_order), dataset_source, term)
write_csv_safe(external_coverage, file.path(tables_dir, "Auto_developmental_reference_celltype_coverage.csv"))
write_csv_safe(external_coverage, file.path(summary_dir, "Auto_developmental_reference_celltype_coverage.csv"))

####################
# Plotting: original heatmap styles with unified readable labels/layout
####################
mp_labels_single <- sapply(mp_order, function(mp) {
  desc <- if (mp %in% names(mp_descriptions)) mp_descriptions[[mp]] else ""
  paste0(mp, " ", desc, " (n=", mp_sizes[[mp]], ")")
})
names(mp_labels_single) <- mp_order

make_matrix_from_table_original <- function(df, value_col, reference, terms, fill = 0) {
  x <- df %>%
    filter(reference == !!reference, term %in% terms, mp %in% mp_order) %>%
    transmute(term, mp, value = .data[[value_col]]) %>%
    pivot_wider(names_from = mp, values_from = value) %>%
    as.data.frame()
  if (nrow(x) == 0) {
    return(matrix(fill, nrow = length(terms), ncol = length(mp_order), dimnames = list(terms, mp_order)))
  }
  rownames(x) <- x$term
  mat <- as.matrix(x[, intersect(mp_order, colnames(x)), drop = FALSE])
  missing_cols <- setdiff(mp_order, colnames(mat))
  if (length(missing_cols) > 0) {
    mat <- cbind(mat, matrix(fill, nrow = nrow(mat), ncol = length(missing_cols), dimnames = list(rownames(mat), missing_cols)))
  }
  missing_rows <- setdiff(terms, rownames(mat))
  if (length(missing_rows) > 0) {
    mat <- rbind(mat, matrix(fill, nrow = length(missing_rows), ncol = length(mp_order), dimnames = list(missing_rows, mp_order)))
  }
  mat <- mat[terms, mp_order, drop = FALSE]
  matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
}

make_display_matrix_original <- function(df, display_col, reference, terms) {
  x <- df %>%
    filter(reference == !!reference, term %in% terms, mp %in% mp_order) %>%
    transmute(term, mp, display = .data[[display_col]]) %>%
    pivot_wider(names_from = mp, values_from = display) %>%
    as.data.frame()
  if (nrow(x) == 0) {
    return(matrix("", nrow = length(terms), ncol = length(mp_order), dimnames = list(terms, mp_order)))
  }
  rownames(x) <- x$term
  mat <- as.matrix(x[, intersect(mp_order, colnames(x)), drop = FALSE])
  missing_cols <- setdiff(mp_order, colnames(mat))
  if (length(missing_cols) > 0) {
    mat <- cbind(mat, matrix("", nrow = nrow(mat), ncol = length(missing_cols), dimnames = list(rownames(mat), missing_cols)))
  }
  missing_rows <- setdiff(terms, rownames(mat))
  if (length(missing_rows) > 0) {
    mat <- rbind(mat, matrix("", nrow = length(missing_rows), ncol = length(mp_order), dimnames = list(missing_rows, mp_order)))
  }
  mat[terms, mp_order, drop = FALSE]
}

term_labels_for_mode <- function(terms, mode) {
  return(terms)
}

heatmap_with_text_original <- function(mat, text_mat, title, legend_title, col_fun, row_labels = rownames(mat),
                                       row_fontsize = NULL, cell_fontsize = NULL, value_format = NULL) {
  n_rows <- nrow(mat)
  if (is.null(row_fontsize)) row_fontsize <- max(10, min(14, 380 / max(n_rows, 1)))
  if (is.null(cell_fontsize)) cell_fontsize <- max(8, min(11, 260 / max(n_rows, 1)))
  Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    na_col = na_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_labels = row_labels,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_fontsize),
    row_names_max_width = unit(16, "cm"),
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_rot = 45,
    column_labels = mp_labels_single[colnames(mat)],
    column_names_gp = gpar(fontsize = 13, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
    cell_fun = function(j, i, x, y, width, height, fill) {
      label <- text_mat[i, j]
      if (!is.na(label) && label != "") {
        if (!is.null(value_format) && is.numeric(mat[i, j]) && is.finite(mat[i, j])) {
          label <- sprintf(value_format, mat[i, j])
        }
        grid.text(label, x, y, gp = gpar(fontsize = cell_fontsize, col = "black", fontface = "bold"))
      }
    }
  )
}

make_external_matrix <- function(terms) {
  if (is.null(external_summary) || nrow(external_summary) == 0) {
    return(NULL)
  }
  x <- external_summary %>%
    filter(term %in% terms) %>%
    select(term, all_of(intersect(mp_order, colnames(external_summary)))) %>%
    as.data.frame()
  if (nrow(x) == 0) {
    return(NULL)
  }
  rownames(x) <- x$term
  mat <- as.matrix(x[, intersect(mp_order, colnames(x)), drop = FALSE])
  missing_rows <- setdiff(terms, rownames(mat))
  if (length(missing_rows) > 0) {
    mat <- rbind(mat, matrix(NA_real_, nrow = length(missing_rows), ncol = length(mp_order), dimnames = list(missing_rows, mp_order)))
  }
  missing_cols <- setdiff(mp_order, colnames(mat))
  if (length(missing_cols) > 0) {
    mat <- cbind(mat, matrix(NA_real_, nrow = nrow(mat), ncol = length(missing_cols), dimnames = list(rownames(mat), missing_cols)))
  }
  mat[terms, mp_order, drop = FALSE]
}

plot_unified_pdf <- function() {
  pdf(file.path(figures_dir, "Auto_developmental_mp_unified_top50.pdf"), width = 30, height = 18, onefile = TRUE)
  for (reference in reference_order) {
    if (!reference %in% names(custom_refs_top50)) next
    terms <- as.character(custom_refs_top50[[reference]]$TERM2NAME$term)
    if (reference %in% c("Normal_Development_long", "Normal_Development_short")) {
      terms <- grep("_Stomach\\.\\.", terms, value = TRUE)
    }
    row_labels <- term_labels_for_mode(terms, "top50")
    ov_mat <- make_matrix_from_table_original(overlap_top50, "neglog10_padj", reference, terms, fill = 0)
    ov_txt <- make_display_matrix_original(overlap_top50, "display", reference, terms)
    cor_mat <- make_matrix_from_table_original(cor_top50, "mean_rho", reference, terms, fill = NA_real_)
    cor_txt <- make_display_matrix_original(cor_top50, "stars", reference, terms)
    cor_txt <- ifelse(is.na(cor_mat), "", paste0(sprintf("%.2f", cor_mat), ifelse(cor_txt == "", "", paste0("\n", cor_txt))))
    ext_mat <- make_external_matrix(terms)

    ov_ht <- heatmap_with_text_original(
      ov_mat,
      ov_txt,
      "Overlap, top 50 ranked genes",
      "Overlap\n-log10 FDR",
      circlize::colorRamp2(seq(0, 7, length.out = length(heat_cols)), heat_cols),
      row_labels = row_labels
    )
    cor_ht <- heatmap_with_text_original(
      cor_mat,
      cor_txt,
      "Expression correlation, top 50 ranked genes",
      "Mean Rho",
      cor_cols,
      row_labels = row_labels
    )
    if (!is.null(ext_mat)) {
      ext_txt <- ifelse(is.na(ext_mat), "", sprintf("%.2f", ext_mat))
      ext_max <- max(ext_mat, na.rm = TRUE)
      if (is.finite(ext_max) && ext_max > 0) {
        ext_cap <- quantile(ext_mat[ext_mat > 0], 0.98, na.rm = TRUE)
        if (is.na(ext_cap) || ext_cap <= 0) ext_cap <- ext_max
      } else {
        ext_cap <- 1
      }
      ext_ht <- heatmap_with_text_original(
        ext_mat,
        ext_txt,
        "Reference celltype MP UCell",
        "Mean UCell",
        circlize::colorRamp2(seq(0, ext_cap, length.out = length(ucell_cols)), ucell_cols),
        row_labels = row_labels,
        cell_fontsize = 8.5
      )
      draw(ov_ht + cor_ht + ext_ht, column_title = reference_pretty[[reference]], column_title_gp = gpar(fontsize = 20, fontface = "bold"))
    } else {
      draw(ov_ht + cor_ht, column_title = paste0(reference_pretty[[reference]], " (no annotated reference cells available)"), column_title_gp = gpar(fontsize = 20, fontface = "bold"))
    }
  }
  dev.off()
}

plot_top50_vs_all_pdf <- function() {
  pdf(file.path(figures_dir, "Auto_developmental_mp_top50_vs_all_overlap_correlation.pdf"), width = 30, height = 18, onefile = TRUE)
  for (reference in reference_order) {
    if (!reference %in% names(custom_refs_all)) next
    terms <- as.character(custom_refs_all[[reference]]$TERM2NAME$term)
    if (reference %in% c("Normal_Development_long", "Normal_Development_short")) {
      terms <- grep("_Stomach\\.\\.", terms, value = TRUE)
    }
    top_labels <- term_labels_for_mode(terms, "top50")

    ov_top_mat <- make_matrix_from_table_original(overlap_top50, "neglog10_padj", reference, terms, fill = 0)
    ov_all_mat <- make_matrix_from_table_original(overlap_all, "neglog10_padj", reference, terms, fill = 0)
    ov_top_txt <- make_display_matrix_original(overlap_top50, "display", reference, terms)
    ov_all_txt <- make_display_matrix_original(overlap_all, "display", reference, terms)
    reference_terms <- term_count_note %>% filter(reference == !!reference)
    n_terms <- nrow(reference_terms)
    n_under_50 <- sum(reference_terms$top50_same_as_all, na.rm = TRUE)
    is_under_50 <- (n_under_50 > (n_terms / 2))
    
    top50_ov_label <- if (is_under_50) "Overlap top 50\n(< 50 total genes)" else "Overlap top 50"
    top50_cor_label <- if (is_under_50) "Correlation top 50\n(< 50 total genes)" else "Correlation top 50"

    ht1 <- heatmap_with_text_original(
      ov_top_mat, ov_top_txt, top50_ov_label, "Overlap\n-log10 FDR",
      circlize::colorRamp2(seq(0, 7, length.out = length(heat_cols)), heat_cols),
      row_labels = top_labels
    )
    ht2 <- heatmap_with_text_original(
      ov_all_mat, ov_all_txt, "Overlap all genes", "Overlap\n-log10 FDR",
      circlize::colorRamp2(seq(0, 7, length.out = length(heat_cols)), heat_cols),
      row_labels = terms
    )
    draw(ht2 + ht1, column_title = paste0(reference_pretty[[reference]], ": overlap all genes vs top 50"), column_title_gp = gpar(fontsize = 20, fontface = "bold"), ht_gap = unit(3, "cm"))

    cor_top_mat <- make_matrix_from_table_original(cor_top50, "mean_rho", reference, terms, fill = NA_real_)
    cor_all_mat <- make_matrix_from_table_original(cor_all, "mean_rho", reference, terms, fill = NA_real_)
    cor_top_star <- make_display_matrix_original(cor_top50, "stars", reference, terms)
    cor_all_star <- make_display_matrix_original(cor_all, "stars", reference, terms)
    cor_top_txt <- ifelse(is.na(cor_top_mat), "", paste0(sprintf("%.2f", cor_top_mat), ifelse(cor_top_star == "", "", paste0("\n", cor_top_star))))
    cor_all_txt <- ifelse(is.na(cor_all_mat), "", paste0(sprintf("%.2f", cor_all_mat), ifelse(cor_all_star == "", "", paste0("\n", cor_all_star))))
    ht3 <- heatmap_with_text_original(
      cor_top_mat, cor_top_txt, top50_cor_label, "Mean Rho",
      cor_cols, row_labels = top_labels
    )
    ht4 <- heatmap_with_text_original(
      cor_all_mat, cor_all_txt, "Correlation all genes", "Mean Rho",
      cor_cols, row_labels = terms
    )
    draw(ht4 + ht3, column_title = paste0(reference_pretty[[reference]], ": correlation all genes vs top 50"), column_title_gp = gpar(fontsize = 20, fontface = "bold"), ht_gap = unit(3, "cm"))
  }
  dev.off()
}

plot_unified_pdf()
plot_top50_vs_all_pdf()

####################
# Run summary
####################
summary_lines <- c(
  paste("Started:", format(start_time)),
  paste("Finished:", format(Sys.time())),
  paste("SCREF_FORCE_REBUILD:", force_rebuild),
  paste("SCREF_REPLOT_ONLY:", replot_only),
  paste("MPs retained:", paste(mp_order, collapse = ", ")),
  paste("Developmental references:", paste(names(custom_refs_all), collapse = ", ")),
  paste("Developmental terms:", n_distinct(term_count_note$term)),
  paste("Ranked reference rows:", nrow(ranked_refs)),
  paste("Top50-equals-all terms:", sum(term_count_note$top50_same_as_all, na.rm = TRUE)),
  paste("External annotated expression datasets available:", sum(reference_expression_sources$available), "of", nrow(reference_expression_sources)),
  paste("Output directory:", out_dir)
)
writeLines(summary_lines, file.path(logs_dir, "Auto_developmental_mp_enrichment_run_summary.txt"))
writeLines(summary_lines, file.path(summary_dir, "Auto_developmental_mp_enrichment_run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
