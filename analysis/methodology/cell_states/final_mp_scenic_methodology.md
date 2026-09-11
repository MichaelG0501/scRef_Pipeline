# Final Centred-MP SCENIC Methodology

`analysis/cell_states/final_mp_scenic.R` identifies regulons associated with the final 17 centred refined MPs in malignant epithelial cells. It reads the persistent merged refined UCell matrix, final MP gene lists, current centred state vector, and raw RNA counts from `EAC_Ref_epi.rds`; it does not use the legacy nMP=19 or 3CA-relabeled state objects.

## MP assignment and cell sampling

UCell scores are centred within `orig.ident` and divided by the within-study standard deviation. Each cell is assigned to its highest adjusted MP. The script records the winning score and the difference from the second score, applies the configured minimum winner/gap filters, requires at least the configured minimum eligible cells per MP, and selects the strongest cells up to `cells_per_mp` per retained MP. The selection cache is versioned as `mp_cell_selection_final17.rds` so an earlier MP2x-containing cache cannot be reused.

The 17-MP panel is exactly the grouping table produced by centred refinement: MP1, MP5, MP13+, MP2+, MP14, MP3+, MP6+, MP11+, MP9+, MP10+, MP8+, MP8b, MP16, MP18b, MP17, MP12, and MP15. MP2x and MP11c are absent because coverage is below three samples; MP18a is the explicit upstream exclusion.

## SCENIC databases and compatibility

Human RcisTarget-compatible motif databases are supplied through `SCENIC_DB_DIR` or the script argument `db_dir`. The verified Imperial directory is `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`, containing the refseq-r80 500-bp TSS and 10-kb TSS mc9nr feather databases. Files without a `features` index are rejected.

Two cluster-specific compatibility patches run before SCENIC: motif annotations fall back from `motifAnnotations_hgnc` to `motifAnnotations_hgnc_v9`, and sparse-aware gene filtering replaces the package-level base `rowSums()` path that fails for `dgCMatrix`. `prepare_only=true` validates inputs, gene sets, selection, and database availability without starting network inference.

Raw counts are restricted to selected cells. Gene filtering requires the larger of 20 cells or 1% of selected cells and the larger of 20 counts or three times 1% of selected cells. GENIE3/SCENIC intermediates are restartable in the ephemeral SCENIC working directory because they are large and reconstructable; selected cells, gene sets, final regulon AUC, RSS matrices, targets, networks, tables, and figures are all copied/saved persistently in `ref_outs/final_mp_scenic/`.

## Interpretation and output

Regulon specificity scores and mean AUCell activity are summarized by MP and state. Network edges represent the script's configured RSS/activity/target filters and are descriptive regulatory associations, not causal validation. All plots and workbooks use full `MP + description` labels. A compact run summary records thresholds, selected cell/MP counts, databases, and produced files in `updates/new_updates/summaries/`.

SCENIC is computationally intensive and must be submitted through PBS. A completed run must be checked for non-empty AUC/RSS/target objects and expected PDF/CSV/XLSX outputs, not merely a zero shell exit code.
