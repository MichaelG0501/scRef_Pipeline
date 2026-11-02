# Spatial Mapping And Per-Platform Neighbourhood Methodology

This document covers the general Yates Visium/Xenium signature mapping and the superseded independent neighbourhood tests. Current Xenium/Visium HD contact inference and pooling are documented in `xenium_visiumhd_meta_methodology.md`; ordinary Visium is excluded from that workflow.

## Inputs And Signature Export

`export_scatlas_visium_signatures.R` exports the canonical centred refined 17-MP panel from `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds` and `centred_refined_mp_state_grouping.csv`. It writes ranked genes, state groups, descriptions, and plot order to the selected mapping output directory. Both platform wrappers call this exporter. The current state unions are Classic proliferation, Squamous-to-intestinal, Glandular-to-intestinal, Stress-adaptive, and Cancer-cell immune mimicry. MP1, MP5, and MP13+ are cell-cycle MPs and are excluded from state assignment and neighbourhood MP sources.

The restored Yates input root is `ref_outs/Yates_spatial/source_data/Zenodo_upload/`, containing `processed/Visium.h5ad`, `processed/Xenium.h5ad`, and `sample_info.csv`. The accompanying README describes five Visium samples and six Xenium patients. Mapping jobs write persistent annotations and plots under `ref_outs/visium_scatlas_states/` and `ref_outs/xenium_scatlas_states/`; reruns overwrite matching outputs. The mapping wrappers do not have plot-only mode, so rerunning them repeats scoring and plotting.

## Visium State Mapping

`map_scatlas_states_visium.py` reads the Cell2Location carcinoma count layer, normalises each spot to 10,000 counts, applies log1p, and computes gene-wise z-scores within each sample’s published `Smoothed label == Tumor` spots. Each MP score is the mean z-score of its top 100 ranked genes. The state-group score is the maximum adjusted MP score in that group. Approach B uses maximum group score 0.5 and top-minus-second gap 0.3; non-tumour spots are labelled `Normal/Mixed` and do not receive malignant state scores.

The Visium H5AD does not contain a `Cell_type` label. For neighbourhood tests, tumour source spots use `Smoothed label == Tumor`. Non-tumour target spots use the largest of the five Cell2Location estimates `Endothelial`, `Fibroblast`, `Lymphoid`, `Muscle`, and `Myeloid`; the winning estimate must be at least one expected cell. This yields one dominant target class per spot and does not resolve within-spot cell mixtures.

## Xenium State Mapping

`map_scatlas_states_xenium.py` reads the processed log1p(CP10K) expression matrix and the published `Cell_type` and spatial coordinates. It retains the top 100 genes per MP when at least five genes overlap the targeted panel. Within each patient, genes are centred/scaled using Carcinoma, Epithelial, and Epithelial (Pre-cancerous) cells. MP scores are centred within patient and scaled by the pooled score SD. Group scores are maxima over current state unions. Approach B assigns Unresolved when the maximum group score is below 0.5 and Hybrid when the top-two group gap is below 0.3. Only Carcinoma cells receive malignant states. The wrapper uses the restored live input root.

## Superseded Independent Neighbourhood Tests (Legacy)

`legacy_spatial_platform_neighbourhoods_pre_ucell.py` preserves the previous per-sample implementation for provenance only. Its 14 non-cell-cycle z-score sources and ordinary Visium scope were superseded on 14 October 2026. Its contact outputs were renamed with `legacy_` paths and must not be used downstream. The remaining description in this section records that historical method.

For Visium, nodes are spots and graph edges use the observed six-neighbour hex-lattice offsets in `array_row/array_col`: `(0,2)`, `(1,1)`, `(1,-1)`, and their reverses. Missing tissue-edge spots naturally have fewer neighbours. For Xenium, nodes are cells and graph edges are unique edges from per-patient Delaunay triangulation; duplicate coordinates or a failed triangulation stop the run. Source nodes are Visium Tumor spots or Xenium Carcinoma cells. Targets are the Visium dominant Cell2Location classes above or every non-Carcinoma Xenium `Cell_type`.

Each test counts observed source-to-target graph edges. Null counts come from 999 label permutations by default, with at least 999 required. Target labels move only within local tiles: 5-by-5 array-coordinate tiles for Visium and 200 µm coordinate tiles for Xenium. Four shifted tilings rotate across permutations. The run stops if fewer than half of targets can move in any tiling. The effect is `log2((observed contacts + 0.5)/(mean null contacts + 0.5))`; two-sided empirical p-values use the add-one correction. Benjamini–Hochberg FDR is calculated within sample and source kind. No samples or platforms are pooled in this workflow.

The interaction table and graph/mobility diagnostics are written under each platform’s `spatial_neighbourhoods/tables/`; per-sample effect heatmaps are under `figures/`. The compact summary is written to `updates/new_updates/summaries/`. `SCREF_REPLOT_ONLY=TRUE` regenerates figures and summary from the persistent interaction table. Otherwise each run recomputes the tests. `SCREF_FORCE_REBUILD=TRUE` is recorded and runs the full recomputation.

Visium spots contain mixed RNA and Cell2Location labels are estimated cell abundances, so Visium contacts indicate tissue-scale co-localisation. Xenium gives cell-level contacts but the targeted panel limits MP overlap; overlap below five genes excludes that MP before state mapping. The cell-centre Delaunay graph is a geometric neighbour definition, not a direct membrane-contact call.

## Visium HD Annotation And CNA Gate

`visiumhd_rctd_annotation.R`, `visiumhd_infercna_malignancy.R`, and the gated stages in `process_visium_hd.py` are documented in `analysis/methodology/spatial/visium_hd_annotation_cnv_methodology.md`. They annotate each modality before state mapping, retain only epithelial observations, define malignancy from InferCNA using two within-sample normal reference compartments, and map states only in malignant epithelial observations. The terminal cancer–TME interaction workflow and its queen-adjacency/local-permutation method are documented separately in `visium_hd_spatial_cancer_tme_interactions_methodology.md`.
