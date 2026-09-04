# scRef Analysis Map

Last audited: 1 September 2026. This is the canonical inventory and dependency map for `analysis/`. Opening script registries provide exact file-level inputs, outputs, cache controls, run command, and environment. Complex methods live under `analysis/methodology/`; simple deterministic utilities may state that a separate methodology is not required.

Audit coverage: 170 R scripts, 28 Python scripts, and 101 shell/PBS scripts. Every script has an opening `Status`, exact `Script`, and `Methodology` field (or an explicit wrapper/simple-script exemption).

## Current analytical contract

- Project and every persistent input/output: `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/`.
- Exception: exceptionally large, fully reconstructable same-script caches may also use the matching ephemeral tree. Anything consumed by another script or needed to replot is persistent in live.
- Canonical MP genes: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`.
- Canonical MP scores: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`.
- Canonical MP grouping/descriptions: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`.
- Canonical state vector/matrices: `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_{states,mp_adj,group_max}.rds`.
- Canonical state method: centred refined noreg, Approach B (`max group <0.5` unresolved; top-two state gap `<0.3` hybrid).
- Final MP panel: 17 MPs. MP2x and MP11c fail coverage in at least three samples; MP18a is the explicit upstream exclusion.
- Legacy `Auto_final_states.rds`, `Auto_topmp_v2_*`, `UCell_nMP19_filtered.rds`, and uncentred `geneNMF_metaprograms_nMP_19.rds` are not valid inputs for current work.

## Canonical run order

1. `metaprograms/centred/01_centred_geneNMF.R` — centred per-sample GeneNMF programs and nMP 8:30 solutions.
2. `metaprograms/centred/02_nmf_rank_selection_diagnostics.R` — rank diagnostics and `optimal_nMP.rds` (currently 19); its colocated `.sh` is the PBS entry point.
3. `metaprograms/centred/03_mp_refinement_submp.R` — silhouette/coverage/gene-count triage and sub-MP splitting.
4. `metaprograms/centred/04_mp_refinement_merge_correlated_submps.R` — correlated sub-MP merging, final filtering, UCell, enrichment, and external scoring.
5. `metaprograms/centred/05_centred_refined_mp_ordered_heatmaps.R` — terminal final-panel heatmaps and grouping table.
6. `metaprograms/centred/06_centred_refined_state_definition_noreg.R` — terminal state vector, adjusted matrices, abundance tables, and state figures.
7. Current downstream scripts below consume only those centred outputs.

## Active and terminal scripts

### Shared

- `shared/scRef_config.R`: current live paths, final 17-MP/state definitions, colours, thresholds, and output conventions.
- `shared/scRef_helpers.R`: cache, output-tier, state-order, label, and run-summary helpers.
- `shared/legacy_utils_reference.R`: legacy helper reference only.

### Metaprograms

- Current chain: `centred/01_centred_geneNMF.R` through `centred/06_centred_refined_state_definition_noreg.R`.
- `centred/3ca_vs_refined_mp_correlation.R`: exploratory current-panel vs external 3CA score correlation.
- `centred/legacy_refined_mp_abundance_comparison.R`: historical centred-versus-uncentred abundance comparison; legacy because its retired uncentred UCell input is absent. The old unprefixed path is retained as a legacy compatibility copy because files must not be deleted.
- `centred/tcga_mp_survival_volcano_centred.R`: current centred TCGA-ESCA MP survival analysis.
- `centred_method_comparison_figures.R`: terminal centred-versus-historical comparison; historical inputs are comparison-only.
- `mp_3ca_ucell_scoring.R`: supporting scorer for the external 3CA panel, not a state-definition step.
- `mp_cancer_type_coverage.R`: external 3CA cancer-type coverage audit.

### Cell states, pseudotime, and regulons

- `state_mp_sample_abundance.R`: current sample-level state/top-MP abundance.
- `final_state_marker_discovery.R`: sample-aware recurrent within-sample Wilcoxon markers for the five current centred states; its ranked live table is the canonical state-marker signature source.
- `basal_smg_marker_mp_dissection_heatmap.R`: current basal/SMG marker expression by dominant within-state MP.
- `top_diverse_sample_state_umap.R`: current most-diverse-sample UMAP (rank 1, not rank 3).
- `pseudotime_top_diverse_samples.R`: current state and within-state MP pseudotime for diverse samples.
- `pseudotime_batch_correction.R`: Harmony/scVI sensitivity analysis using current states.
- `pseudotime_state_distance_matrix.R`: five-method sample-level state distance matrices.
- `pseudotime_trajectory_report_partA.R`: current per-sample trajectory reports and persistent trajectory assets.
- `basal_mp_distance_matrix.R`: six-basal/five-SMG within-state MP distance analysis; MP hybrid gap 0.1 and sample thresholds are documented separately.
- `hybrid_pairwise_distance_nodeplot.R`: hybrid pair network laid out from biological state-distance matrices.
- `state_umap_dispersion_colocalisation.R`: per-sample UMAP dispersion and neighbour colocalisation for current states/MPs.
- `final_mp_scenic.R`: current final-17-MP SCENIC; large restart caches may be ephemeral, all selected cells/results/plots are live.
- `final_mp_scenic_parse_overlap.R`: terminal overlap parsing for completed SCENIC outputs.
- PBS wrappers and `Auto_submit_*` files in this folder submit the corresponding R scripts; they do not define alternative methods.

### Clinical, bulk, and TCGA

- `TCGA/tcga_esca_reconstruct_data.R`: persistent TCGA-ESCA expression/metadata reconstruction.
- `TCGA/tcga_gender_state_mp_validation.R`: current centred gender association validation.
- `TCGA/tcga_cna_recurrent_event_validation.R`: current centred recurrent-CNA validation.
- `clinical/clinical_association_final_boxplots.R`: sample-level current MP/state clinical tests and boxplots.
- `clinical/clinical_association_final_stacked.R`: current state composition across clinical groups.
- `clinical/clinical_association_mp_ucell_plots.R`: sample-mean current MP UCell clinical summaries.
- `clinical/geo_survival_data_prep.R`: persistent GSE19417/GSE13898 downloads and gene-level matrices.
- `clinical/geo_survival_mp_state_survival.R`: current-centred GEO survival scoring.
- `clinical/bulk_tcga_geo_qc.R`: cross-platform QC and retained-sample manifest.
- `clinical/bulk_tcga_geo_integrated_survival.R`: current centred dataset-aware pooled Cox models.
- `clinical/bulk_tcga_geo_meta_survival.R`: current centred dataset-specific/random-effects meta-analysis.
- `clinical/bulk_tcga_geo_feature_presence.R`: current centred feature coverage visualization.
- `clinical/tcga_eac_escc_mp_state_compare.R`: current centred TCGA EAC-versus-ESCC comparison.
- `clinical/tcga_stad_bulk_download_and_gsva.R`: current centred TCGA-STAD location analysis.
- `clinical/cibersortx_sc_reference_export.R`: all-cell CIBERSORTx reference plus current centred gene subset.
<!-- #################### -->
- `clinical/Auto_occams_bulk_mp_survival.R` (PBS: `clinical/Auto_occams_bulk_mp_survival.sh`): current 281-subject QC-pass OCCAMS GRCh37 bulk normalization, 17 centred refined MP/five state-union GSVA scoring, subject-level metadata resolution, and OS Cox/volcano/KM analysis. Persistent normalized expression, scores, model tables, figures, and audit outputs are under `ref_outs/OCCAMS/clinical/`.
- `clinical/Auto_occams_bulk_clinical_associations.R` (PBS: `clinical/Auto_occams_bulk_clinical_associations.sh`): current subject-level OCCAMS MP/state boxplots and dominant-bulk-state stacked plots across core and OAC-specific clinical variables. It consumes the survival workflow's QC-pass scores and harmonized metadata and writes persistent statistics, variable inventory, plot data, and figures under `ref_outs/OCCAMS/clinical/`.
<!-- #################### -->

### CNV and subclones

- `legacy_cnv_profiling.R` and `legacy_cnv_plotting.R`: incomplete historical merged/PDO InferCNA fragments retained for provenance; current CNA inputs come from the core per-sample InferCNA pipeline.
- `legacy_cnv_subsetting.R`: historical incomplete Ju_2025 reference-spiking heatmap; not an active dependency.
- `cnv_malignant_subclone_mp_heatmap.R`: current centred MP/state profiles by malignant CNA subclone.
- `cna_subclone_expression_correlation.R`: sample-blocked robust subclone-expression effects and persistent tables.
- `cna_subclone_expression_correlation_strasser_e17a.R`: current single-sample validation display.
- `mp_chromosomal_mapping.R`: genomic-location mapping of current centred MP genes.
- `Auto_scatlas_numbat_export_inputs.R`, `Auto_scatlas_numbat_run_sample.R`, `Auto_scatlas_numbat_conservative_recut.R`: scATLAS Numbat chain.
- `scatlas_numbat_raw_expression_concordance.R`: terminal raw-expression/Numbat/InferCNA concordance.
- Numbat shell scripts are submit/container wrappers for those R stages.

### Spatial

- `visium_hd_rctd_doublet_detection.R`: binned RCTD metadata; calls do not gate the final custom annotation.
- `visium_hd_celltype_annotation.py`: canonical binned/segmented Visium HD annotation with one-pass targeted graph refinement.
- `visium_hd_binned_post_annotation_filter.R` and `visium_hd_postfilter_diagnostics.py`: final binned filtering and audit.
- `visium_hd_binned_malignancy.R`: binned epithelial malignancy classification.
- `export_visium_hd_binned_scatlas_signatures.R` and `visium_hd_binned_state_mapping.py`: current centred signature export/state mapping.
- `visium_hd_binned_state_colocalisation.R`: terminal spatial state abundance/colocalisation.
- `visium_hd_spatial_cancer_tme_interactions.R`: terminal spatial cancer–TME interaction tests.
- `visium_hd_binned_reference_selection_audit.R`: diagnostic reference-selection comparison.
- `export_scatlas_visium_signatures.R`, `map_scatlas_states_visium.py`, and `map_scatlas_states_xenium.py`: current centred general spatial exports/mapping.
- General Yates Visium/Xenium mapping additionally requires external `processed/Visium.h5ad`, `processed/Xenium.h5ad`, and `sample_info.csv`. The previously configured `spatialtranscriptomics/.../Zenodo_upload` root was absent during the 12 August 2026 audit; signatures are valid, but mapping must wait for those files or an explicit replacement root.
- `run_visium_hd_*.sh`, `map_scatlas_states_*.sh`, and audit wrappers submit the corresponding stages.
- `legacy_Auto_compare_xenium_state_methods.py` and `legacy_Auto_replot_xenium_states_highres.py` are legacy uncentred comparisons.

### Non-malignant NMF

- `non_malignant_nmf/nmf_celltype_geneNMF.R`: per-cell-type NMF.
- `non_malignant_nmf/nmf_celltype_annotation.R`: enrichment-based annotation.
- `non_malignant_nmf/mp_cross_celltype_correlations.R`: correlations against the current centred malignant MP panel.
- Shell scripts in this folder are PBS submit/orchestration wrappers.

### Trajectory / RNA velocity

- `trajectory/scatlas_velocity_prepare_refs.py`: persistent GTF/repeat reference preparation.
- `trajectory/scatlas_velocyto_run.py`: per-BAM velocyto runner.
- `trajectory/scatlas_velocity_metadata.R`: current centred-state/MP metadata export.
- `trajectory/scatlas_velocity_scvelo_visualise.py`: per-sample scVelo using current five states and 17-MP panel.
- `trajectory/scatlas_velocity_nodeplots.R`: terminal sample-aggregated velocity node plots.
- Shell scripts submit/filter/sort the corresponding stages. Large BAM/loom/scVelo working caches may be ephemeral; final matrices, tables, figures, and replot inputs are live.

### Publication layer

- `oac_scatlas_paper/00_verify_live_inputs.py`: lightweight 14-input live check.
- `oac_scatlas_paper/01_prepare_replot_source_data.R`: compact persistent figure source data; no analytical recomputation.
- `oac_scatlas_paper/02_validate_manuscript_outputs.py`: output container/inventory validation.
- `oac_scatlas_paper/figure01` through `figure05`: current terminal replot and conditional schematic scripts.
- `oac_scatlas_paper/shared/`: shared style/config only.
- Shell files submit the corresponding R replot scripts.

### Raw data and summary

- `raw_data/Auto_stage_validate_scatlas_cellranger_outputs.R`: stages and validates rebuilt raw matrices.
- Other `raw_data/*.sh`: download, Cell Ranger, staging, and submit wrappers; raw computation must run on PBS.
- `summary/legacy_cross_sample_summary.R` and `summary/legacy_summary_qc_plots.R`: historical interactive QC summaries; neither is a current dependency.

<!-- #################### -->
### OCCAMS bulk RNA-seq

- `OCCAMS/occams_download_ega_metadata.sh`, `occams_match_metadata.py`, and `occams_download_bams.sh`: credential-safe EGA metadata matching and controlled BAM acquisition. The credential is external at `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/ega.json`; it is never part of this repository.
- `OCCAMS/occams_metadata_coverage_audit.py` and `occams_bam_reference_audit.sh`: terminal audits confirming that all 300 downloaded EGAF BAMs map to OCCAMS metadata and all use GRCh37. The 282-subject cohort is complete despite two additional EGAF accessions returning EGA 403 responses for a represented subject.
- `OCCAMS/occams_prepare_subject_bams.py` and `occams_merge_subject_bams.sh`: deterministic 282-subject symlink/merge preparation in ephemeral storage with full SHA-ID mapping persisted live.
- `OCCAMS/occams_prepare_grch37_reference.sh`: pinned GENCODE v19 GRCh37.p13 download/checksum and `chr`/non-`chr`/mixed-name preflight.
- `OCCAMS/occams_prepare_featurecount_batches.py`, `occams_featurecounts_grch37_batch.sh`, and `occams_combine_featurecount_batches.py`: current eight-batch size-balanced fragment quantification and exact column recombination; `occams_featurecounts_grch37.sh` is the single-job fallback. `occams_build_count_matrix.py` then creates the final live matrix. Large reconstructable featureCounts tables remain ephemeral; the canonical downstream matrix is `ref_outs/OCCAMS/counts/OCCAMS_RNAseq_GRCh37_gene_counts.tsv.gz`.
- `OCCAMS/occams_validate_final_outputs.py`: terminal exhaustive validation of all count values, the 282-subject mapping, metadata/QC coverage, the pinned GTF checksum, and a persistent per-sample QC flag table.
- Never use the historical ephemeral `counts/OCCAMS_RNAseq_raw_counts.txt`; it used a GRCh38 GTF against GRCh37 BAMs.
<!-- #################### -->

### Enrichment and developmental support

- `enrichment/merged_refined_mp_annotation_excel_export.R`: current centred annotation workbook export.
- `enrichment/enrichment_result_extract.R`: exports the persistent current `cluster_enrich_centred.rds` object.
- `enrichment/legacy_nmf_enrichment.R` and `enrichment/legacy_wnt_pathway.R`: session-dependent historical fragments; centred step 04 is the current enrichment workflow.
- `developmental/developmental.R`: builds persistent combined/per-stage developmental references under `ref_outs/developmental_reference/`; centred step 04 consumes these live files.
- `developmental/legacy_developmental_mp_enrichment_unified.R`: legacy nMP=19 validation retained for provenance; centred step 04 is the current developmental/external annotation source.

## Legacy and delete-candidate policy

- Every `legacy_*` script is retained for provenance/comparison and must not feed a current workflow.
- The old uncentred refinement/scoring/correlation scripts under `metaprograms/`, old state/relabel/marker scripts under `cell_states/`, old poster scripts under `publication/`, interactive plotting fragments, old clinical survival scripts, and former Visium HD workflow were renamed with `legacy_` prefixes during the August 2026 audit.
- `spatial/legacy_visiumhd/` contains the former Visium HD pipeline; filenames also carry `legacy_` prefixes.
- `cell_states/Auto_drug_reversal/` contains the historical uncentred drug-reversal workflow; filenames carry `legacy_` prefixes.
- `delete_*` files remain manual deletion candidates. No file is deleted by agents without an explicit user request.
- A legacy script may read ephemeral/historical objects only for provenance. Its outputs are not current dependencies.

## Methodology index

Complex current workflows have detailed operational methods under:

- `methodology/metaprograms/centred_refinement_methodology.md`
- `methodology/metaprograms/centred_refined_mp_ordered_heatmaps_methodology.md`
- `methodology/metaprograms/centred_refined_state_definition_noreg_methodology.md`
- `methodology/cell_states/basal_smg_mp_distance_methodology.md`
- `methodology/cell_states/final_mp_scenic_methodology.md`
- `methodology/cell_states/final_state_marker_discovery_methodology.md`
- `methodology/clinical/clinical_bulk_and_association_methodology.md`
<!-- #################### -->
- `methodology/clinical/Auto_occams_bulk_clinical_methodology.md`
<!-- #################### -->
- `methodology/cnv/` (InferCNA, Numbat, subclone correlation)
- `methodology/spatial/visium_hd_final_annotation_methodology.md`
- `methodology/spatial/visium_hd_binned_*` and `visium_hd_spatial_cancer_tme_interactions_methodology.md`
- `methodology/trajectory/scatlas_velocity_methodology.md`
- `methodology/publication/oac_scatlas_paper_replot_methodology.md`
- `methodology/OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md`

The Visium HD final annotation methodology records every production threshold: segmented QC (`>=100` UMIs, `<15%` mitochondrial), gene/HVG/PCA/neighbour settings, size-aware Leiden 10/6 rule, raw cluster score `>0.25`, targeted trigger `>1.0` in at least 20 observations, local Leiden 1.0, one-pass/no-recursion behavior, and structural ambiguity gap `<=0.25`.

## Header and maintenance rules

An active/terminal analytical script registry must state status, exact script path, accurate description, methodology path or `not required`, exact inputs, tiered outputs, cache/replot behavior, run command, and environment. PBS wrappers must identify the script they submit and include `#PBS -koed`. Whenever a script is added, renamed, superseded, or changes dependencies/outputs, update this map in the same change.
