# scATLAS Numbat raw-expression concordance methodology

This workflow compares the unfiltered expression-derived CNA profile from Numbat with the unfiltered InferCNA expression-CNA profile for the same scATLAS cells.

For each sample in `ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv`, the script reads Numbat `gexp_roll_wide.tsv.gz`. This is the matrix used for Numbat `exp_roll_clust.png`, before Numbat final posterior clone filtering or conservative tree re-cutting. Cells are mapped through `Auto_<sample>_numbat_cell_map.csv` so raw 10x barcodes align to the pipeline cell IDs.

The InferCNA comparison uses `ref_outs/by_samples/<sample>/<sample>_outs.rds` directly. No signal quantile filtering is applied. Genes are intersected between the Numbat raw expression-roll matrix, the InferCNA matrix, and the hg38 gene-order file, then ordered by chromosome and genomic coordinate. The plotted Numbat heatmap therefore uses the same underlying values as the raw Numbat expression-roll heatmap, while the InferCNA heatmap uses the corresponding unfiltered expression-CNA values for the same genes and cells.

The script clusters cells using Numbat raw expression-roll profiles binned by ordered genes only for stable cell ordering and a simple two-cluster raw-expression annotation. This raw-expression cluster is not treated as a Numbat final subclone call. Numbat final clone status is read separately from `Auto_<sample>_numbat_done.txt` and `Auto_<sample>_numbat_clone_post.csv`. A sample can therefore show a visible raw expression-CNA pattern while still being marked as `terminal_no_subclone` if no Numbat clone survived the final Numbat size/tree/CNV filters.

For concordance summaries, both matrices are averaged over ordered gene bins and then chromosome arms. The reported Spearman correlation is computed over matched cell-by-arm values. Where available, InferCNA subclone annotations are read from `ref_outs/Auto_malignant_subclone_mp/Auto_malignant_subclone_cells.csv` and used only as plot annotations and raw-cluster overlap summaries.

Outputs are written to `ref_outs/Auto_scatlas_numbat/raw_expression_concordance/`:

- `figures/Auto_scatlas_numbat_raw_expression_infercna_matched_heatmaps.pdf`: one page per sample.
- `figures/per_sample/Auto_<sample>_raw_expression_infercna_concordance.pdf/png`: per-sample pages.
- `tables/Auto_scatlas_numbat_raw_expression_infercna_summary.csv`: sample-level status, raw-cluster sizes, Numbat final clone status, and arm-level concordance.
- `tables/Auto_scatlas_numbat_raw_expression_cell_clusters.csv`: per-cell raw expression cluster and annotation labels.
- `logs/Auto_scatlas_numbat_raw_expression_concordance_run_summary.txt`: lightweight run summary.
