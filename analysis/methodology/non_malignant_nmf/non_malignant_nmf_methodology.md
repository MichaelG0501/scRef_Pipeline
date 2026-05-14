# Non-Malignant NMF Methodology

This document covers `analysis/non_malignant_nmf/`.

## Per-Celltype GeneNMF

`nmf_celltype_geneNMF.R` runs GeneNMF per non-malignant compartment. The PBS wrapper is `run_geneNMF.sh`.

Cell types:

- `macrophage`
- `fibroblast`
- `endothelial`
- `nk.cell`
- `plasma`
- `cd4`
- `cd8`

The submitter throttles jobs through PBS and keeps the established max-concurrency policy.

## Per-Celltype Annotation

`nmf_celltype_annotation.R` annotates per-celltype NMF programs. The PBS wrapper is `run_annotation.sh`.

Note: GeneNMF uses `nk.cell`, while annotation uses `nk` for the folder map.

## Cross-Celltype MP Correlations

`mp_cross_celltype_correlations.R` builds the cross-celltype MP network. It uses full scATLAS metadata, cancer MPs, non-malignant MP outputs, and ligand-receptor references.

Method details are in:

- `analysis/methodology/non_malignant_nmf/mp_cross_celltype_correlations_methodology.md`

Large intermediate objects are cached under `ref_outs/non_malignant_mp_correlations/cache/`. Use `AUTO_MPXCELL_FORCE_REBUILD=TRUE` to force rebuilding where supported.
