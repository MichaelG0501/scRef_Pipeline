# Enrichment And External Reference Methodology

This document covers `analysis/enrichment/` and the external epithelial reference scorer now located at `analysis/developmental/external_epi_mp_ucell_heatmap.R`.

## MP Enrichment

`enrichment_analysis.R` computes gene-overlap enrichment of MP gene lists against TERM2GENE-style reference databases.

`enrichment_annotation.R` visualizes MP annotation results against Hallmark, GO, 3CA, pan-cancer, developmental, adult epithelium, and Barretts references.

Critical rule:

- Filter out MPs with negative silhouette before interpreting enrichment or scoring outputs.

## Plotting Helpers

`enrichment_plotting_helpers.R` contains reusable enrichment heatmap helpers. New scripts may source this file if they explicitly document the dependency in the script header.

## External Epithelial MP Scoring

`developmental/external_epi_mp_ucell_heatmap.R` scores filtered scATLAS MPs in adult oesophagus, adult stomach, and Barretts epithelial references.

Important external inputs:

- `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Oesophagus/`
- `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Adult_Stomach/data_9_9_annotated_seurat_all_ut.rds`
- `/rds/general/project/spatialtranscriptomics/ephemeral/EAC_data/Barretts/alldatahighquality.rds`

The adult oesophagus Matrix Market dataset is large. The script must subset/sample epithelial barcodes before UCell scoring and cache the sampled subset under `ref_outs/Auto_external_epi_mp_ucell/Auto_cache/` so plot changes can be regenerated without repeating the heavy input read.
