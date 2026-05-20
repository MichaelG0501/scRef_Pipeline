####################
# Auto_drug_reversal_inputs.R
#
# Prepare five-state scRef malignant-state reversal inputs for ASGARD,
# scDrugPrio, and direct CLUE/CMap fallback querying.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(edgeR)
})

####################
# setup
####################

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline"
setwd(file.path(project_dir, "ref_outs"))

out_dir <- "Auto_drug_reversal"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "asgard_inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "scdrugprio_inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "clue_inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "matrix"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "deg_checkpoints"), recursive = TRUE, showWarnings = FALSE)

set.seed(1471)

options(future.globals.maxSize = 50 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)

params <- list(
  min_pct = as.numeric(Sys.getenv("AUTO_DRUG_DEG_MIN_PCT", "0.01")),
  top_n = as.integer(Sys.getenv("AUTO_DRUG_SIGNATURE_TOP_N", "150")),
  force_degs = identical(Sys.getenv("AUTO_FORCE_DRUG_DEGS", "0"), "1"),
  deg_mode = Sys.getenv("AUTO_DRUG_DEG_MODE", "pseudobulk"),
  export_matrix = !identical(Sys.getenv("AUTO_EXPORT_DRUG_MATRIX", "1"), "0"),
  min_cells_per_pseudobulk = as.integer(Sys.getenv("AUTO_DRUG_PSEUDOBULK_MIN_CELLS", "20")),
  min_paired_samples = as.integer(Sys.getenv("AUTO_DRUG_PSEUDOBULK_MIN_SAMPLES", "3"))
)

####################
# helpers
####################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_state_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

pick_logfc_col <- function(df) {
  out <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), colnames(df))[1]
  if (is.na(out)) stop("Could not find a logFC column.")
  out
}

get_assay_layer <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

write_gmt <- function(gene_sets, file) {
  lines <- vapply(names(gene_sets), function(nm) {
    genes <- unique(as.character(gene_sets[[nm]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    paste(c(nm, "na", genes), collapse = "\t")
  }, character(1))
  writeLines(lines, con = file)
}

save_rds_atomic <- function(object, path) {
  tmp_path <- paste0(path, ".tmp")
  saveRDS(object, tmp_path)
  renamed <- file.rename(tmp_path, path)
  if (!renamed) {
    ok <- file.copy(tmp_path, path, overwrite = TRUE)
    unlink(tmp_path)
    if (!ok) stop("Failed to atomically write cache file: ", path)
  }
}

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    warning("AnnotationDbi/org.Hs.eg.db unavailable; writing CLUE GMT files with gene symbols.")
    return(setNames(symbols, symbols))
  }

  mapped <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  mapped <- mapped[!is.na(mapped$ENTREZID) & !duplicated(mapped$SYMBOL), , drop = FALSE]
  out <- mapped$ENTREZID
  names(out) <- mapped$SYMBOL
  out
}

write_status <- function(status, detail) {
  fwrite(
    data.frame(
      step = "inputs",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_drug_reversal_input_status.csv")
  )
}

run_pseudobulk_state_dge <- function(obj, counts_mat, state_name, all_states, params) {
  meta <- obj@meta.data
  meta$cell <- rownames(meta)
  meta$state <- as.character(meta$state)
  meta$orig.ident <- as.character(meta$orig.ident)
  meta <- meta %>%
    filter(state %in% all_states) %>%
    mutate(
      pb_group = ifelse(state == state_name, "target", "rest"),
      pb_id = paste(orig.ident, pb_group, sep = "___")
    )

  cells_use <- intersect(meta$cell, colnames(counts_mat))
  meta <- meta[match(cells_use, meta$cell), , drop = FALSE]
  counts_use <- counts_mat[, cells_use, drop = FALSE]

  pb_design <- Matrix::sparse.model.matrix(~ 0 + factor(meta$pb_id))
  colnames(pb_design) <- sub("^factor\\(meta\\$pb_id\\)", "", colnames(pb_design))
  pb_counts <- counts_use %*% pb_design

  pb_meta <- data.frame(
    pb_id = colnames(pb_counts),
    sample = sub("___(target|rest)$", "", colnames(pb_counts)),
    group = sub("^.*___", "", colnames(pb_counts)),
    n_cells = as.integer(table(meta$pb_id)[colnames(pb_counts)]),
    stringsAsFactors = FALSE
  )

  pb_meta <- pb_meta %>% filter(n_cells >= params$min_cells_per_pseudobulk)
  paired_samples <- pb_meta %>%
    count(sample, name = "n_groups") %>%
    filter(n_groups == 2) %>%
    pull(sample)
  pb_meta <- pb_meta %>% filter(sample %in% paired_samples) %>% arrange(sample, group)

  if (length(unique(pb_meta$sample)) < params$min_paired_samples) {
    warning("Skipping ", state_name, ": only ", length(unique(pb_meta$sample)), " paired pseudobulk samples.")
    return(NULL)
  }

  pb_counts <- as.matrix(pb_counts[, pb_meta$pb_id, drop = FALSE])
  colnames(pb_counts) <- pb_meta$pb_id

  pb_meta$group <- factor(pb_meta$group, levels = c("rest", "target"))
  pb_meta$sample <- factor(pb_meta$sample)
  design <- model.matrix(~ sample + group, data = pb_meta)
  coef_name <- "grouptarget"
  if (!coef_name %in% colnames(design) || qr(design)$rank < ncol(design)) {
    stop("Pseudobulk design is not estimable for ", state_name, ".")
  }

  dge <- DGEList(counts = pb_counts, group = pb_meta$group)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design, robust = TRUE)
  fit <- glmQLFit(dge, design = design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coef_name)

  tt <- topTags(qlf, n = Inf, sort.by = "none")$table %>%
    rownames_to_column("gene") %>%
    as_tibble()

  paired_cells <- meta %>% filter(orig.ident %in% paired_samples)
  target_cells <- paired_cells$cell[paired_cells$state == state_name]
  rest_cells <- paired_cells$cell[paired_cells$state %in% setdiff(all_states, state_name)]
  pct_1 <- Matrix::rowMeans(counts_mat[tt$gene, target_cells, drop = FALSE] > 0)
  pct_2 <- Matrix::rowMeans(counts_mat[tt$gene, rest_cells, drop = FALSE] > 0)

  tt %>%
    transmute(
      state = state_name,
      gene,
      p_val = PValue,
      avg_log2FC = logFC,
      pct.1 = as.numeric(pct_1[gene]),
      pct.2 = as.numeric(pct_2[gene]),
      p_val_adj = FDR,
      logCPM,
      F,
      n_paired_samples = length(unique(pb_meta$sample)),
      n_target_cells = length(target_cells),
      n_rest_cells = length(rest_cells),
      dge_method = "pseudobulk_edgeR_QLF_sample_blocked"
    ) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
}

####################
# data loading
####################

message("Loading scRef epithelial object and reconstructing defined malignant-state object.")

auto_state_cache <- file.path(out_dir, "cache", "Auto_drug_reversal_scRef_states.rds")

sc_ref_state <- NULL
if (file.exists(auto_state_cache) && !identical(Sys.getenv("AUTO_REBUILD_DRUG_STATE5", "0"), "1")) {
  sc_ref_state <- tryCatch(
    readRDS(auto_state_cache),
    error = function(e) {
      warning("Failed to read cached scRef state object ", auto_state_cache, ": ", conditionMessage(e))
      NULL
    }
  )
}

if (is.null(sc_ref_state)) {
  sc_ref_all <- readRDS("EAC_Ref_epi.rds")
  state_labels <- readRDS("Auto_final_states.rds")
  DefaultAssay(sc_ref_all) <- "RNA"

  common_cells <- intersect(colnames(sc_ref_all), names(state_labels))
  keep_cells <- common_cells[state_labels[common_cells] %in% state_order]
  if (length(keep_cells) < 100) {
    stop("Too few cells retained for drug reversal after applying final-state labels.")
  }

  meta_state <- sc_ref_all@meta.data[keep_cells, , drop = FALSE]
  meta_state$state <- factor(as.character(state_labels[keep_cells]), levels = state_order)
  if (!"batch" %in% colnames(meta_state)) {
    meta_state$batch <- if ("study" %in% colnames(meta_state)) {
      as.character(meta_state$study)
    } else {
      as.character(meta_state$orig.ident)
    }
  }

  counts_mat <- get_assay_layer(sc_ref_all, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
  keep_features <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) >= 10]
  counts_mat <- counts_mat[keep_features, , drop = FALSE]

  sc_ref_state <- CreateSeuratObject(counts = counts_mat, meta.data = meta_state, assay = "RNA")
  sc_ref_state <- NormalizeData(sc_ref_state, verbose = FALSE)
  sc_ref_state$state <- factor(as.character(sc_ref_state$state), levels = state_order)
  Idents(sc_ref_state) <- "state"

  rm(sc_ref_all, state_labels, counts_mat, meta_state)
  invisible(gc())

  save_rds_atomic(sc_ref_state, auto_state_cache)
}

if (!"batch" %in% colnames(sc_ref_state@meta.data)) {
  sc_ref_state$batch <- if ("study" %in% colnames(sc_ref_state@meta.data)) {
    as.character(sc_ref_state$study)
  } else {
    as.character(sc_ref_state$orig.ident)
  }
}

sc_ref_state$state <- factor(as.character(sc_ref_state$state), levels = state_order)
Idents(sc_ref_state) <- "state"

state_counts <- sc_ref_state@meta.data %>%
  count(state, orig.ident, batch, name = "n_cells") %>%
  arrange(state, orig.ident)
fwrite(state_counts, file.path(out_dir, "Auto_drug_reversal_state_sample_cell_counts.csv"))

metadata_out <- sc_ref_state@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, orig.ident, batch, state)
fwrite(metadata_out, file.path(out_dir, "Auto_drug_reversal_metadata.csv"))

####################
# matrix export
####################

if (params$export_matrix) {
  message("Exporting sparse five-state count matrix for external wrappers.")
  counts_out <- get_assay_layer(sc_ref_state, assay = "RNA", layer = "counts")
  Matrix::writeMM(counts_out, file.path(out_dir, "matrix", "Auto_drug_reversal_counts.mtx"))
  fwrite(
    data.frame(gene = rownames(counts_out), stringsAsFactors = FALSE),
    file.path(out_dir, "matrix", "Auto_drug_reversal_features.tsv"),
    sep = "\t",
    col.names = FALSE
  )
  fwrite(
    data.frame(cell = colnames(counts_out), stringsAsFactors = FALSE),
    file.path(out_dir, "matrix", "Auto_drug_reversal_barcodes.tsv"),
    sep = "\t",
    col.names = FALSE
  )
}

####################
# differential expression
####################

deg_file <- file.path(out_dir, "Auto_drug_reversal_degs_all_states.csv.gz")
deg_checkpoint_dir <- file.path(out_dir, "deg_checkpoints")

if (file.exists(deg_file) && !params$force_degs) {
  message("Loading cached drug-reversal DEG table.")
  all_degs <- fread(deg_file)
  if (!"dge_method" %in% colnames(all_degs) ||
      !all(all_degs$dge_method == "pseudobulk_edgeR_QLF_sample_blocked", na.rm = TRUE)) {
    message("Cached DEG table is not sample-blocked pseudobulk; rebuilding.")
    all_degs <- NULL
  }
}

if (exists("all_degs") && !is.null(all_degs)) {
  all_degs <- all_degs
} else if (params$deg_mode == "global") {
  message("Using existing descriptive global marker screen as a non-statistical fallback.")
  global_path <- file.path("Auto_five_state_markers", "Auto_five_state_global_marker_screen.csv.gz")
  if (!file.exists(global_path)) {
    stop("AUTO_DRUG_DEG_MODE=global requested, but global marker screen is missing.")
  }
  all_degs <- fread(global_path) %>%
    filter(state %in% state_order) %>%
    transmute(
      state,
      gene,
      p_val = NA_real_,
      avg_log2FC = global_mean_diff,
      pct.1 = global_pct_state,
      pct.2 = global_pct_other,
      p_val_adj = NA_real_,
      global_mean_state,
      global_mean_other,
      global_mean_diff,
      global_pct_state,
      global_pct_other,
      global_pct_delta
    )
  fwrite(all_degs, deg_file)
} else if (params$deg_mode == "pseudobulk") {
  message("Running sample-blocked pseudobulk edgeR QL state-vs-rest DGE for each finalized scRef state.")
  counts_for_dge <- get_assay_layer(sc_ref_state, assay = "RNA", layer = "counts")
  all_degs <- bind_rows(lapply(state_order, function(state_name) {
    checkpoint_file <- file.path(deg_checkpoint_dir, paste0("Auto_deg_", safe_state_name(state_name), ".csv.gz"))
    if (file.exists(checkpoint_file)) {
      checkpoint_dt <- tryCatch(fread(checkpoint_file), error = function(e) NULL)
      if (!is.null(checkpoint_dt) && nrow(checkpoint_dt) > 0) {
        if (!"dge_method" %in% colnames(checkpoint_dt) ||
            !all(checkpoint_dt$dge_method == "pseudobulk_edgeR_QLF_sample_blocked", na.rm = TRUE)) {
          message("  DEGs: ", state_name, " [checkpoint ignored: not pseudobulk]")
        } else {
        message("  DEGs: ", state_name, " [checkpoint]")
        write_status("running", paste("Reusing DEG checkpoint for", state_name))
        return(checkpoint_dt)
        }
      }
    }

    message("  Pseudobulk DEGs: ", state_name)
    write_status("running", paste("Running pseudobulk edgeR QL for", state_name))
    out <- run_pseudobulk_state_dge(
      obj = sc_ref_state,
      counts_mat = counts_for_dge,
      state_name = state_name,
      all_states = state_order,
      params = params
    )
    if (is.null(out) || nrow(out) == 0) return(NULL)
    fwrite(out, checkpoint_file)
    out
  }))
  fwrite(all_degs, deg_file)
} else {
  stop("Unsupported AUTO_DRUG_DEG_MODE='", params$deg_mode, "'. Use pseudobulk.")
}

if (nrow(all_degs) == 0) {
  stop("No DEG rows were produced.")
}

####################
# signatures and wrapper inputs
####################

message("Writing top up/down signatures and wrapper input files.")

ranked_degs <- all_degs %>%
  mutate(
    p_rank_value = ifelse(is.na(p_val_adj), 1, p_val_adj),
    abs_logfc = abs(avg_log2FC)
  )

signature_top <- bind_rows(lapply(state_order, function(state_name) {
  state_df <- ranked_degs %>% filter(state == state_name)
  up <- state_df %>%
    filter(avg_log2FC > 0) %>%
    arrange(p_rank_value, desc(avg_log2FC), desc(pct.1)) %>%
    slice_head(n = params$top_n) %>%
    mutate(direction = "up", signature_rank = row_number())
  down <- state_df %>%
    filter(avg_log2FC < 0) %>%
    arrange(p_rank_value, avg_log2FC, desc(pct.2)) %>%
    slice_head(n = params$top_n) %>%
    mutate(direction = "down", signature_rank = row_number())
  bind_rows(up, down)
})) %>%
  select(state, direction, signature_rank, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, everything())

fwrite(signature_top, file.path(out_dir, "Auto_drug_reversal_signature_top150.csv"))

asgard_gene_list <- lapply(state_order, function(state_name) {
  df <- ranked_degs %>%
    filter(state == state_name) %>%
    transmute(
      gene,
      score = avg_log2FC,
      adj.P.Val = ifelse(is.na(p_val_adj), 1, p_val_adj),
      P.Value = ifelse(is.na(p_val), adj.P.Val, p_val)
    ) %>%
    arrange(adj.P.Val, desc(abs(score)))
  out <- as.data.frame(df[, -1])
  rownames(out) <- df$gene
  out
})
names(asgard_gene_list) <- state_order
saveRDS(asgard_gene_list, file.path(out_dir, "asgard_inputs", "Auto_asgard_gene_list.rds"))

for (state_name in state_order) {
  state_safe <- safe_state_name(state_name)

  asgard_df <- asgard_gene_list[[state_name]] %>%
    rownames_to_column("gene")
  fwrite(
    asgard_df,
    file.path(out_dir, "asgard_inputs", paste0("Auto_asgard_deg_", state_safe, ".txt")),
    sep = "\t"
  )

  scdrug_df <- ranked_degs %>%
    filter(state == state_name) %>%
    transmute(
      gene,
      p_val = ifelse(is.na(p_val), 1, p_val),
      avg_logFC = avg_log2FC,
      pct.1 = pct.1,
      pct.2 = pct.2,
      p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj)
    ) %>%
    arrange(p_val_adj, desc(abs(avg_logFC)))
  fwrite(
    scdrug_df,
    file.path(out_dir, "scdrugprio_inputs", paste0("Auto_scdrugprio_deg_", state_safe, ".txt")),
    sep = "\t"
  )
}

up_symbol_sets <- lapply(state_order, function(state_name) {
  signature_top %>%
    filter(state == state_name, direction == "up") %>%
    arrange(signature_rank) %>%
    pull(gene)
})
names(up_symbol_sets) <- safe_state_name(state_order)

down_symbol_sets <- lapply(state_order, function(state_name) {
  signature_top %>%
    filter(state == state_name, direction == "down") %>%
    arrange(signature_rank) %>%
    pull(gene)
})
names(down_symbol_sets) <- safe_state_name(state_order)

all_symbols <- unique(c(unlist(up_symbol_sets), unlist(down_symbol_sets)))
symbol_to_entrez <- map_symbols_to_entrez(all_symbols)

up_entrez_sets <- lapply(up_symbol_sets, function(x) unname(symbol_to_entrez[x[!is.na(symbol_to_entrez[x])]]))
down_entrez_sets <- lapply(down_symbol_sets, function(x) unname(symbol_to_entrez[x[!is.na(symbol_to_entrez[x])]]))

write_gmt(up_symbol_sets, file.path(out_dir, "clue_inputs", "Auto_clue_up_symbols.gmt"))
write_gmt(down_symbol_sets, file.path(out_dir, "clue_inputs", "Auto_clue_down_symbols.gmt"))
write_gmt(up_entrez_sets, file.path(out_dir, "clue_inputs", "Auto_clue_up_entrez.gmt"))
write_gmt(down_entrez_sets, file.path(out_dir, "clue_inputs", "Auto_clue_down_entrez.gmt"))

fwrite(
  data.frame(
    symbol = names(symbol_to_entrez),
    entrez_id = unname(symbol_to_entrez),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "clue_inputs", "Auto_clue_symbol_to_entrez.csv")
)

####################
# anchor-gene diagnostic
####################

anchor_genes <- ranked_degs %>%
  mutate(
    p_val_adj_for_anchor = ifelse(is.na(p_val_adj), 1, p_val_adj),
    min_pct = pmin(pct.1, pct.2, na.rm = TRUE),
    anchor_score = min_pct - abs(avg_log2FC)
  ) %>%
  filter(abs(avg_log2FC) <= 0.10, pct.1 >= 0.20, pct.2 >= 0.20, p_val_adj_for_anchor >= 0.05) %>%
  group_by(state) %>%
  arrange(desc(anchor_score), desc(min_pct), .by_group = TRUE) %>%
  slice_head(n = 500) %>%
  ungroup() %>%
  select(state, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, anchor_score)

fwrite(anchor_genes, file.path(out_dir, "Auto_asgard_anchor_gene_diagnostic.csv"))

write_status(
  status = "complete",
  detail = paste0(
    "Prepared ", nrow(all_degs), " state-vs-rest rows using mode '",
    params$deg_mode, "' and top ", params$top_n, " up/down signatures for ",
    length(state_order), " states."
  )
)

message("Drug reversal input preparation complete.")
