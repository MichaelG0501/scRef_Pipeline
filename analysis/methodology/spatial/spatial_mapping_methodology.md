# Spatial Mapping Methodology

This document covers `analysis/spatial/`.

`export_scatlas_visium_signatures.R` exports scATLAS state/MP signatures for spatial mapping.

`map_scatlas_states_visium.py` maps scATLAS states to Visium data.

`map_scatlas_states_xenium.py` is currently untracked and should remain unstaged unless the user explicitly requests adding it.

Spatial scripts must document:

- required external spatial object paths
- signature files read from `ref_outs/`
- output locations
- whether they overwrite previous mapping results
- whether they can replot from cached mapped scores

## Visium HD Annotation And CNA Gate

`visiumhd_rctd_annotation.R`, `visiumhd_infercna_malignancy.R`, and the gated stages in `process_visium_hd.py` are documented in `analysis/methodology/spatial/visium_hd_annotation_cnv_methodology.md`. They annotate each modality before state mapping, retain only epithelial observations, define malignancy from InferCNA using two within-sample normal reference compartments, and map states only in malignant epithelial observations.
