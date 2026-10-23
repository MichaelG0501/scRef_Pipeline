# scRef Analysis Map

Last audited: 14 September 2026. This is the canonical inventory and dependency map for `analysis/`. Opening script registries provide exact file-level inputs, outputs, cache controls, run command, and environment. Complex methods live under `analysis/methodology/`; simple deterministic utilities may state that a separate methodology is not required.

Audit coverage: 170 R scripts, 28 Python scripts, and 103 shell/PBS scripts. Every script has an opening `Status`, exact `Script`, and `Methodology` field (or an explicit wrapper/simple-script exemption).

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
- `centred/tcga_mp_survival_volcano_centred.R` (PBS: `centred/tcga_mp_survival_volcano_centred.sh`): current primary-tumour TCGA-ESCA EAC survival analysis using log2(TPM+1), exact final-17 MPs, current state unions and ranked state markers. Standard and exploratory optimal-cut Cox models are univariable with no purity or clinical covariates; persistent model tables and 18x8-inch volcano/10x8-inch KM outputs are under `ref_outs/Metaprogrammes_Results/centred/survival/`.
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
- `Auto_opnme_prioritisation.R`: in-silico prioritisation of 99 opnMe small molecules against five malignant cell states. Three-layer approach: (1) mechanism-based priors from the supervisor annotation workbook `opnMe_99_drug_OAC_state_predictions.xlsx`, (2) per-state target gene expression validation from `EAC_Ref_epi.rds` + centred states, (3) PROGENy pathway activity scoring per state, with optional MOA-class reversal evidence from the existing drug-reversal pipeline. Outputs ranked compound tables and publication figures to `ref_outs/opnMe_prioritisation/`. PBS wrapper: `Auto_opnme_prioritisation.sh`.

### Clinical, bulk, and TCGA

####################
- `clinical/Auto_bayesprism_bulk_deconvolution.R` (PBS: `job_scripts/clinical/Auto_bayesprism_bulk_deconvolution.sh`; dmtcp): corrected common scATLAS reference and raw-count BayesPrism for 281 QC-pass OCCAMS subjects and 88 workbook-selected TCGA EAC primary tumours. Upstream `malignant_level_1`/`malignant_level_2` calls define malignancy; the 21,171 malignant cells absent from the current state vector remain `Malignant state unassigned`, never normal epithelium. The common inference gene universe, initial state fractions, updated malignant-type fractions and inferred malignant expression are persistent under `ref_outs/bulk_malignant_deconvolution/{reference,intermediate,tables,reports,logs,figures/Auto_diagnostics}/`. Run `prepare`, then independent `OCCAMS` and `TCGA` tasks, then `qc`. Prepare uses 4 CPUs/96 GB/1 hour; fits use 32 CPUs/96 GB/12 hours; QC uses 1 CPU/8 GB/30 minutes. Prior OCCAMS-only BayesPrism fits are provenance and are not inputs. Methodology: `methodology/clinical/malignant_bulk_survival_methodology.md`.
- `clinical/Auto_malignant_bulk_survival.R` (PBS: `job_scripts/clinical/Auto_malignant_bulk_survival.sh`; dmtcp): current malignant-compartment OCCAMS/TCGA OS comparison using the exact existing 68/72 surgery-subgroup patient lists, current clinical workbooks and separately documented TCGA endpoint/demographic metadata. Uses identical current 17 MPs, five state MP unions, five state marker signatures and five inferred within-malignant fractions; paired whole-bulk scoring and distinct neoadjuvant analyses provide secondary context. Continuous Cox estimates, categorical stage adjustment, explicit eligibility/model diagnostics and FDR-controlled cross-cohort conjunction tests replace outcome-selected cutpoints for this workflow. Stable current outputs are under `ref_outs/clinical_malignant_survival/{intermediate,tables,reports,logs,figures/Auto_final,figures/Auto_supporting,figures/Auto_diagnostics}/`. The existing selected-survival source tables remain required by other maintained clinical analyses and are not migrated. Matched subgroup definitions do not establish identical treatment histories or survival origins. Methodology: `methodology/clinical/malignant_bulk_survival_methodology.md`.
####################

####################
- `clinical/survival_selected_cohorts.R` (optional minimal PBS launcher: `clinical/survival_selected_cohorts.sh`): standalone R entry point owning revision archiving, caching, analysis and plotting; no audit, Python, or PDF-validation rerun dependency. Selected-cohort OS/RFS stage Kaplan–Meier analysis and final-17-MP OS Cox volcanoes, including a distinct matched surgery-pathway page (68/72 RNA metadata matched; 67 OS-evaluable, 43 deaths). OCCAMS uses matched neoadjuvant 157 (77 biopsy/80 surgery specimens) and matched surgery-pathway 68, combined 225, after three exact exclusions. TCGA retains pathway 86 and surgery-Yes 72. Current figures/tables/cache/reports live directly under `ref_outs/clinical_selected_survival/`; superseded revisions are preserved in `archive/`. Full MP descriptions use fixed label columns; RFS assumptions and optimal-cut limitations remain in methods/tables. Methodology: `methodology/clinical/survival_selected_cohorts_methodology.md`.
- Selected-cohort extension (17 October 2026): `clinical/survival_selected_cohorts.R` also writes `figures/occams_neoadjuvant_biopsy_vs_surgery_boxplots.pdf` (two pages: MPs, then both five-state signature types) and `figures/selected_cohort_mp_and_state_os_volcanoes.pdf` (six cohort pages, all 27 features labelled in each of four OS panels). Additional input: `ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv`; state genes match the integrated script's MP unions and first 20 DGE rows per state. Independent-patient specimen tests compare 77 biopsy versus 80 surgery specimens using shared neoadjuvant-parent GSVA; BH correction is within signature family. Complete scores, gene sets, 648 Cox rows, optimal-cut grids, specimen source/statistics, and label layouts persist under `ref_outs/clinical_selected_survival/tables/`; v4 caches support replot-only. Current launcher: `job_scripts/clinical/survival_selected_cohorts.sh` (4 CPUs, 24 GB, 2 hours, dmtcp; PBS logs in project root), superseding the historical launcher path above. Original PDFs and archived revisions remain preserved.
- `clinical/legacy_validate_selected_survival_pdfs.py` and `clinical/legacy_survival_selected_cohorts.sh`: retired one-time PDF validator and complex wrapper, retained for provenance. Neither is called by current runs. The former independent direction-audit R script is absent; its historical results remain archived.
####################

- `TCGA/Auto_tcga_cna_gene_event_preflight.R` (PBS: `TCGA/Auto_tcga_cna_gene_event_preflight.sh`): gate for the targeted CNA/MP extension; checks manual InferCNA source lineage, reproduces cached arm means from original matrices, verifies `orig.ident` and inventories the official TCGA-ESCA GDC molecular profiles. Outputs are under `ref_outs/TCGA/cna_gene_event_validation/{intermediate,tables,reports,logs}/`; this audit performs no association fitting.
- `TCGA/Auto_tcga_cna_gene_event_validation.R` (PBS: `TCGA/Auto_tcga_cna_gene_event_validation.sh`): ten gain and ten loss gene/MP rectangles using TCGA gene GISTIC CNA and verified manual InferCNA subclone RNA. Four reference-style pages show gain/loss all trends and TCGA BH-FDR<0.05 masks; stars require TCGA significance and matching rho direction. Top-ten arm and gene pairs require matching nonzero rho directions and prioritize FDR support in both datasets then either dataset, then minimize the sum of ranks of maximum nominal P and negative minimum absolute rho across cohorts; no significance/rho cutoff. Rectangle footers are suppressed; definitions remain in source tables/methodology. Gene pairs use scATLAS RNA-tertile and TCGA GISTIC boxplots. Complete ranks, rectangle annotations, exact bin thresholds/observations, source/plot RDS, sensitivity, four PDFs, logs and reports are canonical under `ref_outs/TCGA/cna_gene_event_validation/`; compact summary in `updates/new_updates/summaries/Auto_targeted_cna_summary.csv`. Methodology: `methodology/TCGA/tcga_cna_gene_event_validation_methodology.md`. Independent PBS validator: `temp/Auto_tcga_gene_event_output_validation.sh`; pinned software and reports remain live. Manual source audit `4373597.pbs-7` passed; current rules supersede earlier nominal-P-only and direction-only-star policies. Current analysis `4374755.pbs-7` and validation `4374758.pbs-7` exited 0; four PDFs, 34 pages, all ranking/bin/mask/star and clipping checks passed.
- `TCGA/Auto_tcga_cnv_burden_mp_validation.R` (PBS: `TCGA/Auto_tcga_cnv_burden_mp_validation.sh`): per-cell manual InferCNA burden, exact current state-marker study mapping, within-sample median groups (>=20 cells each), current 17-MP UCell contrasts and TCGA EAC autosomal segment FGA/cohort median groups. Primary continuous-rho PDF has three pages: separated pooled/individual-study/TCGA heatmap, pooled scATLAS clinical-style all-MP boxplots, TCGA all-MP boxplots. Boxes use raw sample-group means and paired Wilcoxon in scATLAS, raw tumour GSVA and unpaired Wilcoxon in TCGA; BH correction across 17 MPs. A separate standardized-effect forest remains. Outer-tertile, middle-excluded sensitivity saves exact cutpoints/counts, sample means and pooled/per-study tests. Canonical sources, caches, full statistics, two PDFs and reports: `ref_outs/TCGA/cnv_burden_mp_validation/{intermediate,tables,figures,reports,logs}/`; compact summary: `updates/new_updates/summaries/Auto_burden_summary.csv`. User-requested paired/per-study boxplot PDFs and per-study forest PDF are no longer exported or present in the replot bundle. Methodology: `methodology/TCGA/tcga_cnv_burden_mp_validation_methodology.md`; independent PBS validator: `temp/Auto_tcga_burden_output_validation.sh`. No Numbat inputs or cell-level pseudoreplication. Current presentation validated by analysis `4375938.pbs-7` and final cached replot/PDF review `4376118.pbs-7`, both exit 0; prior independent R statistical checks passed and seven source/statistics table hashes remain identical after graphics-only corrections. Two PDFs/four content pages; three explicitly retired PDFs stay absent after replot.

- `clinical/Auto_occams_bayesprism_qc_plots.R` (PBS: `clinical/Auto_occams_bayesprism_qc_plots.sh`): four-panel numerical QC dashboard for OCCAMS BayesPrism; initial/updated tumour fraction, MP coverage, aligned state composition and MP heatmap. Persistent PDF/PNG, plot-data CSVs, report and compact summary are under `ref_outs/OCCAMS/bayesprism/` and `updates/new_updates/summaries/`. Methodology: `methodology/clinical/bayesprism_malignant_state_deconvolution_methodology.md`.

- `TCGA/tcga_esca_reconstruct_data.R` (PBS: `TCGA/Auto_tcga_esca_reconstruct_data.sh`): canonical TCGA-ESCA TPM/metadata reconstruction under `ref_outs/TCGA/esca_gdc_reconstruction/`; pinned downloads, clinical TSV snapshot and GDC query caches in `intermediate/source_data/`, matrix/source RDS in `intermediate/`, CSV/mixture/AJCC annotation template in `tables/`, run summaries in `logs/`. Retains live compatibility metadata/mixture copies. `SCREF_FORCE_REBUILD=TRUE` rebuilds query/expression caches; no replot mode. Existing historical `raw/` and intermediate CSVs are retained but no longer refreshed.
- `clinical/tcga_surgery_clinical_associations.R` (PBS: `job_scripts/clinical/tcga_surgery_clinical_associations.sh`; 2 CPUs, 8 GB, 1 hour, dmtcp): general clinical boxplots for exactly 72 TCGA EAC surgery-Yes patients from the corrected September-18 AJCC8 workbook, supplemented with cBioPortal demographics/exposures and exact reconstruction sample mapping. Reuses hash-checked selected-survival 86-parent GSVA scores (17 MPs and both five-state signatures). Prespecified important variables, >=10 patients/category, Wilcoxon/Kruskal–Wallis, BH per variable/signature family. Each clinical variable occupies one page containing all 17 MPs and one page containing both state definitions; fonts never cause extra pages. Stable root `ref_outs/TCGA/surgery_clinical_associations/`: two priority PDFs in `figures/Auto_final/`, two secondary PDFs in `figures/Auto_supporting/`, full source/statistical tables, metadata audit, cached replot data, current linked report, and canonical clinical TSV in `reference/`. Methods: `methodology/clinical/tcga_surgery_clinical_associations_methodology.md`. Existing sex-concordance and EAC/ESCC comparison workflows retain their distinct purposes.
- `TCGA/tcga_gender_state_mp_validation.R` (PBS: `TCGA/Auto_tcga_gender_state_mp_validation.sh`): depends on reconstruction; exact current centred 17-MP/five-state sex association validation. Canonical root `ref_outs/TCGA/gender_validation/`: fingerprinted GSVA scores, clinical workbook snapshot, sample-level plot source and rendered figure bundle in `intermediate/`; statistics in `tables/`, four-page PDF/PNG in `figures/`, run/replot summaries in `logs/`; compact CSV in `updates/new_updates/summaries/`. Force rebuild recomputes scores; replot-only renders saved grobs without analysis.
- `TCGA/tcga_cna_recurrent_event_validation.R` (PBS: `TCGA/Auto_tcga_cna_recurrent_event_validation.sh`): depends on reconstruction, the validated gender GSVA cache and `Auto_cna_subclone_expression/threshold_1/` tables/live RDS. Canonical root `ref_outs/TCGA/cna_recurrent_event_validation/`: segment snapshot, arm cache, complete plot source and figure bundle in `intermediate/`; associations/cluster assignments in `tables/`, eight PDF products in `figures/`, run/replot summaries in `logs/`; compact CSV in `updates/new_updates/summaries/`. No alternative score producer or missing-input panel skips. Force rebuild refreshes local arm calls; replot-only renders saved grobs. The user supplied the canonical annotation workbooks at `analysis/TCGA/OAC_CNV.xlsx` and `analysis/TCGA/41588_2018_331_MOESM3_ESM.xlsx` after the initial October 2026 audit. Production regeneration uses these exact live inputs and the existing reconstructed TPM/metadata; gender GSVA is rebuilt before CNA validation. Regeneration and final validation completed on 15 October 2026 (final CNA `4373346.pbs-7`, validation `4373347.pbs-7`, both exit 0); 87 EAC samples were tested at arm threshold 0.12. Exact provenance and output checks are in `ref_outs/TCGA/cna_recurrent_event_validation/reports/Auto_production_run_summary.txt`.
- `TCGA/Auto_tcga_workflow_helpers.R`: deterministic shared PBS/storage guards, one-time live input imports, fingerprinted cache validation, exact MP/group enforcement and rendered-grob and native fixed-order heatmap save/replot helpers. Standalone cache rendering explicitly initializes ggplot2/gtable drawing methods to prevent blank gtables when caller libraries are absent, and advances the Cairo device only between real pages so layout measurement cannot leave a blank first page. Methodology is not required for this utility; workflow-specific semantics are documented in `methodology/TCGA/`. Validation harness: `temp/Auto_tcga_contract_validation.sh`; isolated fixture products are retained under `ref_outs/TCGA/contract_validation/<job_id>/`. Production validation uses `temp/Auto_tcga_production_validation.sh` and writes inventories/page counts, PDF metadata and preview PNGs using the pinned `ref_outs/clinical_selected_survival/validation_dependencies/pymupdf_1.26.7/` package to each workflow's `reports/` tier.
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
- `clinical/Auto_bayesprism_malignant_state_deconvolution.R` (PBS: `clinical/Auto_bayesprism_malignant_state_deconvolution.sh`): OCCAMS-only BayesPrism trial using one malignant cell type with five current malignant states plus `Unresolved` and `Hybrid`; reports malignant fractions, state fractions, tumour-only expression, and tumour-only current MP scores. `SCREF_RESUME_EXTRACTION=TRUE` regenerates validated extraction/GSVA outputs and the final report from the fitted RDS without resampling. Methodology: `methodology/clinical/bayesprism_malignant_state_deconvolution_methodology.md`.
<!-- #################### -->
- `clinical/Auto_occams_bulk_mp_survival.R` (PBS: `clinical/Auto_occams_bulk_mp_survival.sh`): current 281-subject QC-pass OCCAMS GRCh37 bulk normalization, 17 centred refined MP/five state-union GSVA scoring, source-column-specific subject metadata, and OS Cox/volcano/KM analysis. Diagnostic, pretreatment, and resection-pathology fields and TNM editions are never coalesced; optimal-cut models are univariable. Persistent normalized expression, scores, model tables, figures, and audit outputs are under `ref_outs/OCCAMS/clinical/`.
- `clinical/Auto_occams_bulk_clinical_associations.R` (PBS: `job_scripts/clinical/Auto_occams_bulk_clinical_associations.sh`; 2 CPUs, 12 GB, 1 hour, dmtcp; logs in project root): corrected September-25 patient-master workbook is authoritative for all covered fields; only unavailable comorbidity remains from the older resolved metadata. Uses canonical 281-subject MP/MP-union/DGE GSVA matrices on a common scale, testing five selected OCCAMS cohorts (225/157/77/80/68) matching selected-survival's first five cohort pages. Ten PDFs under `ref_outs/OCCAMS/clinical/figures/selected_cohorts/` contain five cohort views per type, one eligible clinical variable per page; State pages show both signature families. Clinical/pathological AJCC stages, curated initial treatment pathway, and RNA specimen/treatment status are included. Cohort/level eligibility and BH correction are recalculated per cohort/variable/signature family. Persistent corrected metadata, source/change audit, cohort membership, input manifest, test tables, PDF inventory, and hash-validated plot cache remain in the clinical output tiers. Previous association-owned outputs are preserved in `archive/clinical_associations/`; upstream metadata/scores and unrelated figures are untouched. Methodology: `methodology/clinical/Auto_occams_bulk_clinical_methodology.md`.
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
- `visium_hd_spatial_cancer_tme_interactions.R` (PBS: `visium_hd_spatial_cancer_tme_interactions.sh`): terminal validation of scRNA-predicted cancer–TME interactions in three Visium HD sections. It uses the exact current 17-MP grouping/descriptions/order, UCell activity, immediate queen contact with 32/48 µm sensitivity rings, shifted 320 µm local permutations, empirical sample p-values, equal-weight pooling, and Spearman-direction candidate-family FDR. Persistent scores, memberships, eligibility/candidate tables, labelled validation, four-mode significant-only comparison, and three-section validated maps are under `ref_outs/visium_hd_outs/spatial_interactions/`.
- `visium_hd_binned_reference_selection_audit.R`: diagnostic reference-selection comparison.
- `export_scatlas_visium_signatures.R`, `map_scatlas_states_visium.py`, and `map_scatlas_states_xenium.py`: current centred 17-MP signature export and Yates state mapping. PBS wrappers are `map_scatlas_states_visium.sh` and `map_scatlas_states_xenium.sh`.
- Yates source files are staged under `ref_outs/Yates_spatial/source_data/Zenodo_upload/`: `processed/Visium.h5ad`, `processed/Xenium.h5ad`, and `sample_info.csv`. The README describes five Visium samples and six Xenium patients. Xenium wrapper and default input now use this restored live root.
- `spatial_platform_neighbourhoods.py` and `xenium_hd_ucell_scores.R` (PBS: `Auto_xenium_neighbourhoods.sh`): corrected Xenium/Visium HD workflow. Xenium uses raw-count R UCell for all 17 malignant MPs and the exact HD TME signatures, approved CAF/TAM/T-cell unions, CD4/CD8 programmes in the shared T-cell parent, calibrated pixel-to-micrometre coordinates (0.2125 µm/pixel), HD-equivalent physical Chebyshev queen windows (±16/32/48 µm), shifted 320 µm local nulls and 4,999 permutations. All four interaction modes feed per-sample signed-Stouffer meta-analysis with within-platform, balanced-platform and leave-one-sample sensitivities. Xenium persistent scores/memberships/contacts are under `ref_outs/xenium_scatlas_states/spatial_neighbourhoods/`; meta tables, figures and report are under `ref_outs/spatial_cross_platform/`. Methodology: `methodology/spatial/xenium_visiumhd_meta_methodology.md`.
- `legacy_spatial_platform_neighbourhoods_pre_ucell.py` and `legacy_Auto_visium_neighbourhoods.sh` preserve the superseded Visium-inclusive, Xenium z-score contact implementation. Its contact outputs and compact summaries have `legacy_` names and are excluded from the current workflow. Ordinary Visium contributes no interaction statistics to the meta-analysis.
- `run_visium_hd_*.sh`, `map_scatlas_states_*.sh`, and neighbourhood wrappers submit the corresponding stages.
- `legacy_Auto_compare_xenium_state_methods.py` and `legacy_Auto_replot_xenium_states_highres.py` are legacy uncentred comparisons.

### Non-malignant NMF

- `non_malignant_nmf/nmf_celltype_geneNMF.R`: per-cell-type NMF.
- `non_malignant_nmf/nmf_celltype_annotation.R`: enrichment-based annotation.
- `non_malignant_nmf/mp_cross_celltype_correlations.R` (PBS: `non_malignant_nmf/mp_cross_celltype_correlations.sh`): sample-blocked correlations against the current centred 17-MP panel and five noreg states; it sources current labels/descriptions from `shared/scRef_config.R`, requires current state markers, and writes the candidate tables consumed by spatial validation.
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
