# Visium HD Annotation And CNA Methodology

## Scope

This workflow gates scATLAS MP/state mapping to malignant epithelial Visium HD observations. It supports RCTD-binned, custom-binned, and Space Ranger cell-segmented matrices. Canonical binned annotation is RCTD-based, whereas custom RCTD-singlet bins and segmented cells deliberately share the hierarchical manual-marker strategy described below.

## Inputs

- Binned: `square_016um/filtered_feature_bc_matrix.h5` and `spatial/tissue_positions.{parquet,csv}`.
- Segmented: `segmented_outputs/filtered_feature_cell_matrix.h5` and `cell_segmentations.geojson`.
- scRNA reference: `ref_outs/EAC_Ref_merged.rds` under the live project path.
- MP signatures: the centred refined MP gene-set RDS exported by `export_scatlas_visiumhd_signatures.R`.
- Gene order: `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`.

## Binned Annotation

`visiumhd_rctd_annotation.R` constructs an RCTD reference from `celltype_update` (falling back to compatible metadata columns), omits `unresolved_inconsistent`, combines `t.cell` and `nk.cell` into `t_nk.cell`, and samples at most 2,000 cells per reference type, except 6,000 epithelial cells. It runs RCTD in `doublet` mode on raw 16 um bin UMI counts after a 200-UMI filter.

Only `spot_class == "singlet"` and `first_type == "epithelial"` pass downstream. `doublet_certain`, `doublet_uncertain`, and `reject` observations remain in the saved RCTD table/object but are excluded from epithelial CNA/state mapping. This follows the RCTD documentation: doublet mode records `spot_class`, `first_type`, `second_type`, and doublet weights. At 16 um, this is an operational mixed-bin exclusion, not proof of a technical doublet.

## Segmented Annotation

####################
`process_visium_hd.py --stage annotate --mode segmented` applies the hierarchical manual-marker strategy from the Visium HD notebooks after QC (`total_counts >= 200`, mitochondrial fraction `<= 15%`). Raw marker means are retained only for audit. Marker module scores are calculated after CP10K/log1p normalisation, matching the documented scale of the Scanpy/Seurat score method, and are not thresholded for annotation. The former raw-UMI gates (`0.05`, `0.015`, `0.5`, `1.0`, `1.5`, and `2.5`) are retired because their scale changes with library size and is therefore not transferable between segmented cells and 16 um bins.

Cluster labels first require cluster-versus-rest enrichment of the same marker genes on CP10K expression. A marker is supported when its log2 fold change is at least the selected shared effect-size threshold, its within-cluster detection and detection increase pass the selected shared prevalence thresholds, and its one-sided hypergeometric detection-enrichment test has Benjamini-Hochberg adjusted `P <= 0.05`. Panels with at least three available genes require two supported markers. One- and two-gene panels require one supported marker, preventing obligatory dropout-driven false negatives for lymphatic (`CCL21`) and keratinocyte (`FLG`, `IVL`) while retaining their highly specific markers. Each output records the decision stage and supporting markers, and separate cluster/type and cluster/gene evidence CSVs retain all effect sizes, prevalences, and adjusted P values.

The decision is deliberately asymmetric after the protected first stage. Supported non-epithelial/non-fibroblast lineages are selected first. If none passes, fibroblast requires differential enrichment of at least two ECM-panel genes and therefore remains a specificity call despite stromal spillover. Only if fibroblast also fails can epithelial be assigned by positive control-adjusted CP10K module score plus co-detection of at least two epithelial markers. To reproduce the original notebook's complete annotation behavior without weakening the spillover guard, remaining low-evidence clusters use the highest protected-lineage raw marker mean when it exceeds `0.10`; the final residual is assigned to the larger epithelial or fibroblast raw marker mean. These two fallback stages occur only after all differential and residual evidence stages and are explicitly recorded as lower-confidence evidence.
####################

####################
It clusters observations using HVGs, depth/mitochondrial regression, PCA, cosine neighbours, and Leiden clustering. Cluster labels follow the marker-evidence hierarchy above; an epithelial/fibroblast fallback can never override a supported first-stage call.
####################

The marker panels include fibroblast, macrophage, mast, epithelial, T, B, NK, plasma, dendritic, endothelial, lymphatic, erythrocyte, keratinocyte, and neutrophil genes. Segmented cells receive no synthetic doublet call because this marker-based method is not a doublet detector; this status is explicitly recorded as `not_applicable_manual_segmented`.

## Custom Binned Annotation

####################
`process_visium_hd.py --stage annotate --mode custom` reuses the completed RCTD doublet-mode table only to separate singlet from non-singlet 16 um bins. Each RCTD non-singlet (`doublet_certain`, `doublet_uncertain`, or `reject`) remains excluded from downstream CNA/state mapping and is labelled from its RCTD `first_type|second_type` fields. RCTD singlets are not assigned by RCTD cell type: every singlet is passed through CP10K marker scoring, cluster-marker enrichment, and the same hierarchical lineage logic as segmented observations. The segmented-cell 200-UMI/15%-mitochondrial exclusion is not applied to custom bins because RCTD has already defined the analysed bin universe; this prevents excluded bins from being reintroduced as `unresolved` labels. No RCTD singlet type is used as an annotation fallback.
####################

####################
Custom singlets use the notebook's fine Leiden setting (`resolution=6`) to retain rare 16 um-bin compartments; segmented cells use resolution 1. Resolution controls only cluster granularity. Both representations use exactly the same biological decision rule and thresholds. For every cluster, all supported non-epithelial/non-fibroblast types are ranked first while epithelial and fibroblast evidence is ignored. Fibroblast-specific differential evidence is evaluated second, epithelial absolute residual evidence third, protected absolute evidence fourth, and the epithelial/fibroblast forced fallback last. This preserves the spillover guard while producing the complete labels expected from the original notebooks.

`visiumhd_manual_annotation_threshold_calibration.py` scans one shared grid of marker log2 fold-change, within-cluster prevalence, and prevalence-increase thresholds. A candidate is ineligible unless both segmented and custom annotations retain at least two of endothelial, macrophage, and fibroblast with at least 20 observations in every sample. Among eligible candidates, selection maximises `0.6 *` spatial macro-recall against the independently segmented representation plus `0.4 *` macro-recall against RCTD singlet labels. RCTD singlet labels are used only for this external audit and never become custom singlet annotations. Sample-level validation, reference counts, the complete grid, and the selected thresholds are persisted under `ref_outs/visium_hd_outs/tables/`.

The final 20 Jul 2026 calibration, rerun after retaining every RCTD singlet, selected log2 fold change `>= 1.0`, within-cluster detection `>= 0.10`, detection increase over the rest `>= 0.05`, and Benjamini-Hochberg adjusted `P <= 0.05`. Its pooled validation score was `0.5568` (`0.6 *` spatial macro-recall `0.5274` plus `0.4 *` RCTD macro-recall `0.6009`), mean spatial nearest-neighbour accuracy was `0.8179`, and mean RCTD-singlet accuracy was `0.7597`. The resolved fraction was `1.0`, and all six sample/representation combinations retained at least two eligible normal reference types. These are the production defaults for both manual modes; the full grid remains the sensitivity analysis rather than implying that the selected agreement is perfect.

Methodological anchors are the Scanpy `score_genes` documentation and cell-cycle example (normalise and log-transform before module scoring), the UCell paper's demonstration that rank-based signatures are more robust to dataset composition, SCINA/Sargent's positive marker-evidence framework, and the spacexr/RCTD doublet-mode definitions. The implementation does not copy a universal cutoff from those methods: the shared Visium HD thresholds above are selected from the persisted cross-resolution sensitivity grid because segmented cells and 16 um bins have different depth and capture characteristics.
####################

####################
The custom annotation table retains the RCTD source fields (`Auto_rctd_first_type`, `Auto_rctd_second_type`, `Auto_rctd_spot_class`), manual-QC status, manual cluster, CP10K module scores, raw marker means, and selected marker evidence. Therefore RCTD is still auditable for all bins, but only non-singlet labels are used from RCTD. Custom manually called keratinocytes are non-epithelial and are excluded before InferCNA rather than being corrected after CNA calling.
####################

## InferCNA Malignancy Classification

`visiumhd_infercna_malignancy.R` is run per sample and mode after annotation. Targets are all gated epithelial observations. For a fair binned-annotation comparison, binned RCTD and binned custom use the identical RCTD-singlet observation universe and the identical two most abundant eligible RCTD normal compartments; therefore a shared barcode has identical CNA signal/correlation coordinates and the same reference-derived thresholds in both binned outputs. After strict gene/barcode validation, custom reuses the finalized binned InferCNA matrix and recomputes its annotation-dependent epithelial, signature, and malignancy gates. The custom annotation changes which coordinates are treated as epithelial targets, not the CNA coordinate system. The segmented custom method uses its own segmented endothelial/macrophage/fibroblast references because cells and 16 um bins are distinct measurement units. Every reference type requires at least 20 observations; this retains the SUR1231 RCTD fibroblast/macrophage pair (112/21), while fewer than two eligible compartments still stops that sample before state mapping.

If fewer than two eligible reference groups exist, the script writes the per-sample annotation/reference statistics, skips CNA inference for that sample, and exits non-zero after all samples have been assessed. This prevents downstream state mapping with an invalid reference baseline.

Raw counts are converted to CPM and restricted to genome-ordered genes before InferCNA. `cnaScatterPlot()` is calculated with reference cells excluded from the tumour-profile average. A gated epithelial observation is `malignant_level_1` when both its CNA signal and CNA correlation exceed the reference mean plus one standard deviation. The two-metric requirement is retained, but one SD is used because single segmented cells have substantially broader normal-reference scatter than 16 um bins; `--cna-sd-k` records and permits explicit sensitivity changes. CNA-unresolved epithelial observations are then rescored using the `cancer_signatures.txt` procedure from `Malignancy.R`: among present signature genes, the 50 most expressed across the sample's epithelial targets are scored from `log1p(CPM / 100)`, with score `>= 1` called signature malignant. Only CNA-unresolved and signature-malignant observations become `malignant_level_2`; CNA-non-malignant observations are not rescued. The scatter thresholds, tier counts, proportions, per-observation calls, and saved CNA matrices are retained under `ref_outs/visium_hd_outs/malignancy/`.

The summary explicitly reports signature-positive epithelial cells stratified by CNA class. This avoids treating a small level-2 count as a failed signature calculation: signature-positive CNA-malignant cells are already level 1. For segmented cells only, a signature-positive CNA-non-malignant cell can become level 2 when the nearest epithelial 16 um bin from the same sample is CNA-malignant. This cross-resolution rescue requires two independent evidence layers and records the matched bin barcode/distance and `binned_cna_plus_signature` evidence label; signature-positive cells without that bin evidence remain unresolved.

For binned data, RCTD epithelial annotations remain unchanged even when keratinocytes are absent from the RCTD reference. The malignancy layer maps segmented cells to their nearest binned epithelial target (maximum distance 25 spatial pixels), and a bin with at least one mapped cell and keratinocyte fraction `>= 0.5` is assigned `normal_keratinocyte`. This overrides an otherwise level-1/2 malignancy call before state mapping, while preserving the original CNA values, RCTD annotation, and pre-override malignancy label for audit.

The initial spatial keratinocyte-normal label is subsequently replaced by full-CNA-profile classification from the cached InferCNA matrix. For keratinocyte-dominant bins, correlation is calculated against the centroid of non-keratinocyte malignant level-1 bins and against the normal-reference centroid. A keratinocyte bin is CNA-like malignant only when both its tumour-centroid correlation exceeds the 99th percentile of normal-reference correlations and its tumour-minus-normal correlation margin exceeds the 95th percentile of the normal references. Bins below both thresholds are normal-like; discordant bins are indeterminate and are excluded from state mapping. This avoids using scalar CNA scatter coordinates or cancer-signature score alone when those distributions overlap.

A sample with no segmented keratinocyte-dominant bins requires no profile correction: its malignancy calls are retained and the profile summary records `no_keratinocyte_targets`. This is a valid terminal status, provided malignant anchors and normal references remain sufficient. In the final run, FFPEA1 had zero keratinocyte targets, while SUR1231 and FFPED1 retained 1,551 and 1,481 targets for profile classification.

## State Mapping

The mapper scores MPs across all annotated epithelial observations per sample, preserving an epithelial-only centring baseline, then exports/plots state calls for both malignant tiers (`Auto_malignant == TRUE`). Output tables retain `Auto_malignancy`, CNA class, and signature score/status alongside the state calls; pre-state epithelial scores are cached under `intermediate/`; state maps are written under `figures/`.

`visiumhd_annotation_diagnostics.R` writes per-sample spatial annotation/malignancy plots and expression UMAP annotation/malignancy plots. It uses an existing Space Ranger projection where available and otherwise computes a Seurat UMAP once, caching the coordinates under `ref_outs/visium_hd_outs/intermediate/` for plot-only reruns.

## Execution

Submit `analysis/spatial/run_visium_hd_states.sh` with PBS. It executes signature export, reuses completed RCTD tables by default, regenerates binned/segmented annotations as applicable, runs both InferCNA modes, maps states, and writes diagnostics. `SCREF_REUSE_SEGMENTED_ANNOTATION=TRUE` reuses a deliberately refreshed segmented annotation table during a CNA-only rerun. It includes SUR1231, FFPEA1, and FFPED1 in every stage. Do not run the RCTD, InferCNA, or UMAP stages on a login node.

`SCREF_RUN_MODE=custom` runs the complete custom-binned branch: custom annotation, InferCNA, malignant state mapping, and spatial/UMAP diagnostics. It reuses the existing RCTD doublet tables unless `SCREF_REUSE_RCTD=FALSE` is explicitly set. Set `SCREF_REUSE_CUSTOM_ANNOTATION=TRUE` only after a successful calibrated annotation-only job to regenerate the expensive downstream tiers without reclustering.

`SCREF_RUN_MODE=binned_downstream` rebuilds only binned RCTD annotation, common-reference InferCNA, malignant state mapping, and binned diagnostics, then exits before segmented processing. `visiumhd_compare_annotation_infercna.R` reads the three completed cell tables and writes one sample per PDF page with common axes across `Binned RCTD`, `Binned custom`, and `Segmented custom`, plus paired binned-coordinate correlations proving whether the two binned modes are aligned.

The final comparison has paired binned CNA signal and CNA-correlation coefficients of `1.0` with maximum absolute difference `0` in SUR1231, FFPEA1, and FFPED1. Common binned InferCNA peaked near 129 GB and requires a 192 GB PBS request; segmented custom peaked near 105 GB and uses 128 GB. Exact-barcode per-sample cache validation supports recovery without repeating completed rolling-CNA calculations.

When one custom sample has fewer than two manual normal-reference types, its InferCNA guard remains active and its status/counts are retained in `Auto_visiumhd_custom_infercna_malignancy_summary.csv`. The launcher still produces custom spatial/UMAP diagnostics for all samples and maps only the custom samples whose CNA classification completed; it does not create an invalid SUR1231 state map without normal references.

For recovery after a segmented-only failure, `SCREF_RUN_MODE=segmented_downstream` reuses completed binned outputs and runs only segmented InferCNA, mapping, and diagnostics. `SCREF_RUN_MODE=segmented_annotation` refreshes only the manual segmented annotation table.

`SCREF_RUN_MODE=binned_keratinocyte_reclassify` applies the cached binned keratinocyte-normal correction and regenerates binned state maps/diagnostics without rerunning InferCNA.

`SCREF_RUN_MODE=binned_keratinocyte_profile` performs the final cached full-CNA-profile keratinocyte classification and regenerates binned state maps/diagnostics without rerunning InferCNA.

`visiumhd_keratinocyte_evidence_audit.R` creates `Auto_visiumhd_keratinocyte_vs_epithelial_evidence_audit.pdf` under `malignancy/figures/`. It presents two readable pages per sample: a thresholded four-panel InferCNA scatter and a CNA signal/correlation/cancer-signature distribution page. Each scalar density and boxplot panel has its active dotted threshold and numeric threshold label.

## External Guidance

- [spacexr RCTD documentation](https://github.com/dmcable/spacexr) describes doublet-mode `spot_class`, `first_type`, `second_type`, and doublet weights.
- [infercna documentation](https://jlaffy.github.io/infercna/articles/infercna.html) specifies that absolute reference correction requires two or more normal reference cell groups and documents `cnaScatterPlot`, CNA signal, and CNA correlation.
