####################
# Analysis registry (authoritative override):
#   Status: legacy; retained for provenance, no current downstream use
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_prepare_asgard_reference.R
#   Methodology: analysis/methodology/cell_states/Auto_drug_reversal_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: This workflow consumes the superseded pre-centred PDO state
#     or marker route. It must be redesigned against centred states before reuse.
####################

####################
# Auto_prepare_asgard_reference.R
#
# Build ASGARD tissue-specific L1000 rank-matrix references from downloaded GEO files.
####################

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
setwd(project_dir)

ref_root <- Sys.getenv(
  "AUTO_ASGARD_REF_ROOT",
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Auto_drug_reversal/asgard_l1000"
)
raw_dir <- file.path(ref_root, "raw")
plain_dir <- file.path(ref_root, "plain")
drug_ref_dir <- file.path(ref_root, "DrugReference")
status_dir <- file.path(project_dir, "ref_outs", "Auto_drug_reversal", "asgard_reference")

dir.create(drug_ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)

target_tissue <- Sys.getenv("AUTO_ASGARD_TISSUE", "stomach")
target_safe <- gsub(" ", "-", target_tissue)

write_status <- function(status, detail) {
  data.table(
    step = "asgard_reference",
    status = status,
    detail = detail,
    reference_root = ref_root,
    target_tissue = target_tissue
  ) |>
    fwrite(file.path(status_dir, "Auto_asgard_reference_status.csv"))
}

required_files <- c(
  cell_info = file.path(plain_dir, "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt"),
  gene_info = file.path(plain_dir, "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt"),
  sig_70138 = file.path(plain_dir, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"),
  sig_92742 = file.path(plain_dir, "GSE92742_Broad_LINCS_sig_info.txt"),
  gctx_70138 = file.path(plain_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"),
  gctx_92742 = file.path(plain_dir, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  write_status("missing_downloads", paste("Missing files:", paste(missing, collapse = "; ")))
  quit(save = "no", status = 0)
}

if (!requireNamespace("Asgard", quietly = TRUE)) {
  write_status("missing_package", "R package Asgard is not installed in the active environment.")
  quit(save = "no", status = 0)
}

message("Building ASGARD L1000 references under: ", drug_ref_dir)
Asgard::PrepareReference(
  cell.info = required_files[["cell_info"]],
  gene.info = required_files[["gene_info"]],
  GSE70138.sig.info = required_files[["sig_70138"]],
  GSE92742.sig.info = required_files[["sig_92742"]],
  GSE70138.gctx = required_files[["gctx_70138"]],
  GSE92742.gctx = required_files[["gctx_92742"]],
  Output.Dir = paste0(drug_ref_dir, "/")
)

rank_path <- file.path(drug_ref_dir, paste0(target_safe, "_rankMatrix.txt"))
gene_path <- file.path(drug_ref_dir, paste0(target_safe, "_gene_info.txt"))
drug_path <- file.path(drug_ref_dir, paste0(target_safe, "_drug_info.txt"))

if (!all(file.exists(c(rank_path, gene_path, drug_path)))) {
  available <- list.files(drug_ref_dir, pattern = "_rankMatrix\\.txt$", full.names = FALSE)
  write_status(
    "target_tissue_missing",
    paste(
      "Requested target tissue was not generated.",
      "Available rank matrices:",
      paste(available, collapse = ", ")
    )
  )
  quit(save = "no", status = 0)
}

path_table <- data.table(
  AUTO_ASGARD_TISSUE = target_tissue,
  AUTO_ASGARD_DRUG_RESPONSE = rank_path,
  AUTO_ASGARD_GENE_INFO = gene_path,
  AUTO_ASGARD_DRUG_INFO = drug_path
)
fwrite(path_table, file.path(status_dir, "Auto_asgard_reference_paths.csv"))

writeLines(
  c(
    paste0("export AUTO_ASGARD_TISSUE=", shQuote(target_tissue)),
    paste0("export AUTO_ASGARD_DRUG_RESPONSE=", shQuote(rank_path)),
    paste0("export AUTO_ASGARD_GENE_INFO=", shQuote(gene_path)),
    paste0("export AUTO_ASGARD_DRUG_INFO=", shQuote(drug_path))
  ),
  con = file.path(status_dir, "Auto_asgard_reference_paths.sh")
)

write_status("complete", paste("Prepared ASGARD reference for tissue:", target_tissue))
message("ASGARD reference preparation complete.")
