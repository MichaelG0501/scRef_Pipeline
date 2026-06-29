# scRef Analysis Map

This document is the canonical map for `analysis/`. Update it whenever a script is added, renamed, superseded, or moved.

## Current Defaults

- Preferred state definition: `Approach B, noreg`
- Preferred final state object: `ref_outs/Auto_final_states.rds`
- Upstream noreg state object: `ref_outs/Auto_topmp_v2_noreg_states_B.rds`
- Upstream noreg MP score object: `ref_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- Shared constants: `analysis/shared/scRef_config.R`
- Shared helper functions: `analysis/shared/scRef_helpers.R`
- Output tiers for new long-running analyses: `intermediate/`, `tables/`, `figures/`, `logs/`, `reports/`

## Run Order

1. Metaprogram inputs
   - `analysis/metaprograms/mp_ucell_scoring.R`
   - `analysis/metaprograms/mp_3ca_ucell_scoring.R`
   - Optional MP refinement/splitting: `analysis/metaprograms/mp_refinement_submp.R`

2. Current epithelial state definition
   - `analysis/cell_states/state_definition_approach_b_reg_noreg.R`
   - Use the `noreg` outputs for all downstream state-like analyses.

3. Final state relabeling
   - `analysis/cell_states/final_state_unresolved_relabel.R`
   - Writes `ref_outs/Auto_final_states.rds`, the preferred final state vector.

4. Final-state terminal figures and tables
   - `analysis/cell_states/state_mp_sample_abundance.R`
   - `analysis/cell_states/final_state_overall_proportions.R`
   - `analysis/cell_states/top_diverse_sample_state_umap.R`
   - `analysis/cell_states/basal_smg_mp_signature_heatmap.R`
   - `analysis/cell_states/final_state_marker_discovery.R`
   - `analysis/cell_states/final_mp_scenic.R`

5. Pseudotime and hybrid analyses
   - `analysis/cell_states/pseudotime_top_diverse_samples.R`
   - `analysis/cell_states/pseudotime_state_distance_matrix.R`
   - `analysis/cell_states/hybrid_pairwise_percell_heatmap.R`
   - `analysis/cell_states/hybrid_pairwise_distance_nodeplot.R`
   - `analysis/cell_states/pseudotime_trajectory_report_partA.R`
   - `analysis/cell_states/pseudotime_batch_correction.R`

6. Clinical and bulk survival analyses
   - `analysis/clinical/tcga_data_prep.R`
   - `analysis/clinical/tcga_mp_state_survival_reg_noreg.R`
   - `analysis/clinical/geo_survival_data_prep.R`
   - `analysis/clinical/geo_survival_mp_state_survival.R`
   - `analysis/clinical/bulk_tcga_geo_qc.R`
   - `analysis/clinical/bulk_tcga_geo_integrated_survival.R`
   - `analysis/clinical/tcga_mp_state_survival_qc_filtered.R`
   - Terminal association plots: `analysis/clinical/clinical_association_final_boxplots.R`, `analysis/clinical/clinical_association_mp_ucell_plots.R`

7. CNV and subclone analyses
   - `analysis/cnv/cnv_profiling.R`
   - `analysis/cnv/cnv_subsetting.R`
   - `analysis/cnv/cnv_plotting.R`
   - `analysis/cnv/cnv_malignant_subclone_mp_heatmap.R`
   - `analysis/cnv/cna_subclone_expression_correlation.R`
   - `analysis/cnv/cna_subclone_expression_correlation_strasser_e17a.R`

8. Optional external references, enrichment, spatial, and non-malignant NMF
   - See folder-specific tables below.

## Dependency Map

| Script | Status | Key inputs | Key outputs | Downstream use |
| :--- | :--- | :--- | :--- | :--- |
| `shared/scRef_config.R` | shared | none | constants | all new scripts |
| `shared/scRef_helpers.R` | shared | `scRef_config.R` | helper functions | all new scripts |
| `shared/legacy_utils_reference.R` | legacy reference | none | none | do not source in new scripts |
| `metaprograms/mp_ucell_scoring.R` | active upstream | `EAC_Ref_epi.rds`, `geneNMF_metaprograms_nMP_19.rds` | `Metaprogrammes_Results/UCell_nMP19_filtered.rds` | state definition, MP plots, survival |
| `metaprograms/mp_3ca_ucell_scoring.R` | active upstream | `EAC_Ref_epi.rds`, 3CA `New_NMFs.csv` | `UCell_3CA_MPs.rds` | unresolved relabel, SCENIC |
| `metaprograms/mp_refinement_submp.R` | active optional refinement/terminal figure | `geneNMF_metaprograms_nMP_19.rds`, `geneNMF_outs.rds`, `EAC_Ref_epi.rds` | `Metaprogrammes_Results/mp_refinement/` refined sub-MP genes, UCell scores, split diagnostics, refined correlation/Jaccard PDFs and matrix tables | exploratory refined MP interpretation; does not replace canonical nMP19 scoring unless explicitly adopted |
| `metaprograms/delete_mp_external_scoring_superseded.R` | delete candidate | old 3CA workflow | old diagnostics | superseded by `mp_3ca_ucell_scoring.R` |
| `metaprograms/final_mp_correlation.R` | terminal figure | MP genes, UCell scores, final states if present | final MP correlation/Jaccard PDFs | terminal |
| `metaprograms/mp_database_correlation.R` | active upstream/terminal | `cluster_enrich.rds`, `EAC_Ref_epi.rds`, developmental per-stage references including `Adult_Epithelium` and `Barretts_Oesophagus` | ref-term UCell/correlation PDFs/RDS, Adult/Barretts reference gene-list XLSX | enrichment interpretation |
| `metaprograms/mp_database_correlation.sh` | PBS wrapper | `metaprograms/mp_database_correlation.R` | live PBS log/output for MP database correlation regeneration | run wrapper |
| `metaprograms/mp_pdo_sc_crossdataset_correlation.R` | terminal comparison | PDO UCell, scRef UCell | cross-dataset MP comparison PDFs | terminal |
| `metaprograms/mp_pancancer_correlation.R` | terminal comparison | scRef MPs, pan-cancer MP scores | pan-cancer correlation plots | terminal |
| `metaprograms/mp_cancer_type_coverage.R` | terminal figure | pan-cancer MP annotations | cancer-type coverage figures | terminal |
| `metaprograms/nmf_rank_selection_diagnostics.R` | diagnostic | NMF rank metrics | rank-selection plots | no downstream dependency |
| `metaprograms/pdo_robust_nmf_annotation.R` | optional upstream | PDO NMF outputs | PDO NMF enrichment plots | PDO-only comparisons |
| `metaprograms/legacy_mp_correlation_sc_kmeans_state_temp.R` | legacy | `UCell_default.rds` | `state_temp.rds`, old heatmaps | no new downstream use |
| `metaprograms/legacy_mp_correlation_pdo_kmeans_state_temp.R` | legacy | PDO UCell | old PDO k-means state heatmaps | no new downstream use |
| `cell_states/state_definition_approach_b_reg_noreg.R` | active upstream | `EAC_Ref_epi.rds`, filtered MP UCell, `meta_full_epi.rds`, cell-cycle genes | `Auto_topmp_v2_noreg_states_B.rds`, `Auto_topmp_v2_noreg_mp_adj.rds`, paired reg/noreg QC plots | final state relabel, abundance, pseudotime, clinical |
| `cell_states/final_state_unresolved_relabel.R` | active upstream | noreg Approach B states, 3CA UCell, MP UCell, TCGA inputs | `Auto_final_states.rds`, task4 relabel tables/figures | preferred state object for downstream scripts |
| `cell_states/state_mp_sample_abundance.R` | terminal figure/table | `EAC_Ref_epi.rds`, final or noreg states, noreg MP matrix | sample abundance PDF/summary | terminal |
| `cell_states/final_state_overall_proportions.R` | terminal figure | final or noreg states | overall state proportion PDF | terminal |
| `cell_states/top_diverse_sample_state_umap.R` | terminal figure | `EAC_Ref_epi.rds`, noreg states | top-diverse sample UMAP PDF | terminal |
| `cell_states/basal_smg_mp_signature_heatmap.R` | terminal figure/table | noreg states, noreg MP matrix, MP genes, cell-cycle genes | basal/SMG signature heatmap and bubble plots | terminal |
| `cell_states/final_state_marker_discovery.R` | terminal report/table | `EAC_Ref_epi.rds`, `Auto_final_states.rds` | six-state markers, heatmaps, tables | terminal and marker interpretation |
| `cell_states/final_mp_scenic.R` | heavy terminal workflow | final states, scRef/3CA UCell, MP genes, cisTarget DBs | SCENIC selected cells, regulons, networks | terminal |
| `cell_states/pseudotime_top_diverse_samples.R` | heavy terminal/intermediate | final or noreg states, noreg MP matrix | Monocle pseudotime RDS/PDF summaries | terminal; reusable pseudotime objects |
| `cell_states/pseudotime_state_distance_matrix.R` | heavy intermediate/terminal | `EAC_Ref_epi.rds`, noreg states | state-distance matrix RDS/CSV/PDF | input for distance nodeplot |
| `cell_states/hybrid_pairwise_distance_nodeplot.R` | terminal figure | noreg states/group max, state-distance matrices | distance-aware hybrid nodeplots | terminal |
| `cell_states/hybrid_pairwise_percell_heatmap.R` | terminal figure/table | noreg states, noreg MP matrix, UCell | hybrid subtype RDS/heatmaps/UMAPs | terminal |
| `cell_states/pseudotime_trajectory_report_partA.R` | terminal report | Part A pseudotime outputs | trajectory report PDF, CDS/projection cache | terminal; replot-friendly |
| `cell_states/pseudotime_batch_correction.R` | exploratory heavy | `EAC_Ref_epi.rds`, final or noreg states | Harmony/scVI pseudotime diagnostics | exploratory terminal |
| `cell_states/cell_annotation.R` | annotation utility | per-sample annotated RDS | updated cell type calls | historical/non-main downstream |
| `cell_states/cell_typing.R` | annotation utility | annotated RDS, marker lists | cell type calls | historical/non-main downstream |
| `cell_states/cancer_summary.R` | summary | per-sample epithelial/malignancy RDS | cancer count summaries | terminal summary |
| `cell_states/legacy_*` | legacy | varies | old comparison/diagnostic outputs | no new downstream use |
| `clinical/clinical_association_final_boxplots.R` | terminal figure/table | `meta_full_epi.rds`, final states, MP UCell, clinical xlsx | final association boxplot PDFs/stats | terminal |
| `clinical/clinical_association_final_stacked.R` | untracked terminal figure | `meta_full_epi.rds`, final states, clinical xlsx | final stacked association PDFs | terminal; currently unstaged |
| `clinical/clinical_association_mp_ucell_plots.R` | terminal figure/table | `meta_full_epi.rds`, MP UCell, final-state naming | MP clinical association PDFs | terminal |
| `clinical/tcga_mp_state_survival_reg_noreg.R` | active survival | TCGA meta/TPM, MP genes, reg/noreg states, final states | TCGA Cox/volcano outputs | superseded for QC-filtered TCGA by filtered script but still reference |
| `clinical/tcga_mp_state_survival_qc_filtered.R` | active survival | bulk QC retained TCGA samples, TCGA inputs, MP/state genes | QC-filtered TCGA survival outputs | terminal |
| `clinical/geo_survival_data_prep.R` | active upstream | GEO downloads/platform annotation | GEO meta/expression/cache files | GEO survival |
| `clinical/geo_survival_mp_state_survival.R` | active survival | GEO prepared data, MP/state gene sets | GEO Cox/volcano outputs | terminal |
| `clinical/bulk_tcga_geo_qc.R` | active upstream | TCGA TPM/meta, GEO expression/meta | cross-platform QC tables/RDS/PDF | integrated survival |
| `clinical/bulk_tcga_geo_integrated_survival.R` | active survival | cross-platform QC outputs, MP/state gene sets | pooled survival CSV/PDF | terminal |
| `clinical/bulk_tcga_geo_meta_survival.R` | active survival | harmonized TCGA/GEO metadata | metadata survival outputs | terminal |
| `clinical/bulk_tcga_geo_feature_presence.R` | terminal QC | TCGA/GEO expression, MP/state sets | feature-presence dotplots/tables | terminal |
| `clinical/cibersortx_sc_reference_export.R` | active export | `EAC_Ref_merged_strict.rds`, `meta_full_epi.rds` | CIBERSORTx reference matrix/labels | external CIBERSORTx |
| `clinical/tcga_data_prep.R` | external-data upstream | GDC/cBioPortal/TCGA files | TCGA meta/TPM matrix | TCGA survival |
| `clinical/survival_cibersort.R` | legacy/terminal | `state_temp.rds`, DEG cache, CIBERSORTx | old survival PDFs | no new state downstream use |
| `clinical/legacy_*` | legacy | varies | old comparison or superseded clinical outputs | no new downstream use |
| `cnv/cnv_malignant_subclone_mp_heatmap.R` | active terminal | per-sample InferCNA, epithelial RDS, UCell, state vectors | subclone MP/state PDFs and CSVs | terminal |
| `cnv/cna_subclone_expression_correlation.R` | active terminal | CNA/subclone outputs, inferCNA per-sample matrices, OAC/OCCAMS xlsx | `ref_outs/Auto_cna_subclone_expression/` with arm CNA, consensus heatmap, recurrent event associations, largest-subclone, pairwise distance figures/tables/RDS | terminal; merged from previous v1+v2 split |
| `cnv/cna_subclone_expression_correlation_strasser_e17a.R` | active terminal figure | `Auto_malignant_subclone_cells.csv`, `UCell_nMP19_filtered.rds` | single-sample MP/subclone boxplot PDF | terminal |
| `cnv/scatlas_numbat_raw_expression_concordance.R` | active terminal validation | scATLAS Numbat `gexp_roll_wide.tsv.gz`, per-sample InferCNA `_outs.rds`, Numbat cell maps, optional InferCNA subclone labels | `ref_outs/Auto_scatlas_numbat/raw_expression_concordance/` with matched raw Numbat-vs-InferCNA heatmaps, per-sample plots, raw-expression cluster tables, and concordance summary | terminal Numbat-vs-InferCNA validation; uses raw expression-roll values, not final Numbat joint-post filters |
| `cnv/mp_chromosomal_mapping.R` | untracked terminal | MP gene positions | MP chromosomal mapping figures | terminal; currently unstaged |
| `cnv/cnv_profiling.R` | core CNV | epithelial per-sample RDS | inferCNA profiles | CNV subsetting/plotting |
| `cnv/cnv_subsetting.R` | core CNV | inferCNA profiles | filtered CNA matrices | CNV plotting |
| `cnv/cnv_plotting.R` | core CNV terminal | CNA matrices | publication CNV plots | terminal |
| `developmental/developmental.R` | reference build | developmental xlsx files | `enrich_dev.rds`, per-stage RDS | enrichment/reference interpretation |
| `developmental/external_epi_mp_ucell_heatmap.R` | active terminal/upstream | MP genes, adult oesophagus/stomach, Barretts references | external epithelial MP UCell heatmaps/cache | terminal |
| `developmental/developmental_mp_enrichment_unified.R` | active terminal/upstream | MP genes/UCell, `EAC_Ref_epi.rds`, developmental ranked marker workbooks, available annotated external references | unified developmental MP overlap/correlation/reference-celltype UCell PDFs, top50-vs-all comparison, rank/external-data audit tables | terminal developmental validation |
| `developmental/developmental_mp_enrichment_unified.sh` | PBS wrapper | `developmental_mp_enrichment_unified.R` | live PBS log/output for unified developmental MP validation | run wrapper |
| `enrichment/enrichment_analysis.R` | upstream | MP gene lists, TERM2GENE DBs | `cluster_enrich.rds` | enrichment annotation |
| `enrichment/enrichment_annotation.R` | terminal | MP object, `cluster_enrich.rds`, reference DBs | enrichment heatmaps | terminal |
| `enrichment/enrichment_plotting_helpers.R` | helper/reference | enrichment data frames | helper functions/plots | reuse or reference |
| `enrichment/enrichment_result_extract.R` | utility | `cluster_enrich.rds` | extracted enrichment tables | terminal tables |
| `enrichment/mp_annotation_excel_export.R` | export | MP annotations | Excel-ready MP matrices | terminal export |
| `enrichment/nmf_enrichment.R` | terminal | NMF outputs, DBs | enrichment heatmaps/xlsx | terminal |
| `enrichment/wnt_pathway.R` | terminal | `cluster_enrich.rds` | WNT overlap plots | terminal |
| `non_malignant_nmf/nmf_celltype_geneNMF.R` | active upstream | non-malignant cell type, full atlas | per-celltype NMF outputs | annotation/cross-celltype |
| `non_malignant_nmf/nmf_celltype_annotation.R` | active upstream | per-celltype NMF outputs | per-celltype enrichment outputs | cross-celltype interpretation |
| `non_malignant_nmf/mp_cross_celltype_correlations.R` | active terminal | full atlas, cancer/non-malignant MP outputs, LR reference | cross-celltype MP networks and LR workbooks | terminal |
| `plotting/publication_umap.R` | terminal figure | `EAC_Ref_epi.rds`, state metadata | publication UMAP PDFs | terminal |
| `plotting/gene_expression_heatmap.R` | terminal figure | expression matrix/gene lists | expression heatmaps | terminal |
| `plotting/qc_heatmap.R` | terminal QC | per-sample RDS/QC metadata | QC heatmaps | terminal |
| `summary/cross_sample_summary.R` | terminal summary | per-sample RDS outputs | cross-sample summary tables/plots | terminal |
| `summary/legacy_summary_qc_plots.R` | legacy | merged loose object | old QC plots | no new downstream use |
| `spatial/export_scatlas_visium_signatures.R` | active export | MP object | Visium signature CSVs | spatial mapping |
| `spatial/map_scatlas_states_visium.py` | active spatial | Visium data/signatures | mapped scATLAS states | terminal spatial |
| `spatial/map_scatlas_states_xenium.py` | untracked spatial | Xenium data/signatures | mapped scATLAS states | terminal spatial; currently unstaged |

## Superseded Or No-Downstream Scripts

These are intentionally retained but prefixed so they are not confused with current workflows:

- `analysis/cell_states/legacy_state_definition_topmp_residual.R`
- `analysis/cell_states/legacy_state_definition_cluster_louvain.R`
- `analysis/cell_states/legacy_state_definition_cluster_vs_topmp_comparison.R`
- `analysis/cell_states/legacy_state_definition_topmp_hybrid_v2.R`
- `analysis/cell_states/legacy_manual_state_assignment.R`
- `analysis/cell_states/legacy_manual_state_umap.R`
- `analysis/cell_states/legacy_state_qc_manual_assignment.R`
- `analysis/metaprograms/legacy_mp_correlation_sc_kmeans_state_temp.R`
- `analysis/metaprograms/legacy_mp_correlation_pdo_kmeans_state_temp.R`
- `analysis/clinical/legacy_tcga_survival_clinical_mps_old.R`
- `analysis/clinical/legacy_tcga_survival_clinical_mps_v2_old.R`
- `analysis/clinical/legacy_clinical_association_topmp_v2B_reg_noreg.R`
- `analysis/clinical/legacy_clinical_variable_plots.R`

Delete-candidates are not removed by agents. They are prefixed for manual deletion:

- `analysis/metaprograms/delete_mp_external_scoring_superseded.R`

## Outdated Downstream Pointers To Avoid

- Do not use `ref_outs/state_temp.rds` for new state analyses.
- Do not use `ref_outs/Auto_topmp_states.rds` for downstream figures.
- Do not use `ref_outs/Auto_cluster_states.rds` for downstream figures.
- Do not use `ref_outs/Auto_topmp_v2_states_B.rds`; use `ref_outs/Auto_topmp_v2_noreg_states_B.rds` or, preferably, `ref_outs/Auto_final_states.rds`.
- Do not use old clinical state plots under `analysis/plotting/`; final clinical association scripts now live in `analysis/clinical/`.

## Cache And Replot Policy

Long-running scripts should write data-heavy or model-heavy intermediates to `intermediate/`, final CSV/TSV tables to `tables/`, plot files to `figures/`, run summaries to `logs/`, and multi-page narrative outputs to `reports/`. Replotting should be possible by reading `intermediate/` and regenerating only `figures/` or `reports/`.

Recommended environment toggles:

- `SCREF_FORCE_REBUILD=TRUE`: ignore cached intermediates.
- `SCREF_REPLOT_ONLY=TRUE`: read cached intermediates and regenerate plots only.

## Publication Poster Figures

| Script | Status | Inputs | Outputs | Notes |
|---|---|---|---|---|
| `publication/publication_helpers.R` | active helper | shared scRef config/helpers | helper functions only | canonical poster state/MP order, colours, saving, placeholders |
| `publication/poster_section1_atlas_metaprograms.R` | terminal publication figure | scRef nMP19 GeneNMF, UCell, final/noreg states, enrichment RDS | `ref_outs/publication/section1_atlas_metaprograms/` | NMF composition, MP correlation, state MP heatmap, state abundance, focused enrichment |
| `publication/poster_section2_genetics_regulons.R` | terminal publication figure | CNA subclone visualisation RDS, optional final SCENIC outputs | `ref_outs/publication/section2_genetics_regulons/` | CNA consensus/association replot; writes SCENIC placeholder when upstream result absent |
| `publication/poster_section3_tme_interactions.R` | terminal publication figure | cross-celltype MP correlation and LR CSVs | `ref_outs/publication/section3_tme_interactions/` | relaxed-threshold cancer-TME dotmap and simplified LR network |
| `publication/poster_section4_pdo_concordance.R` | terminal publication figure | PDO nMP13 GeneNMF, PDO states, 3CA crossdata, matched PDO/snRNA CSVs | `ref_outs/publication/section4_pdo_concordance/` | PDO NMF, scAtlas-vs-PDO state occurrence, pan-cancer concordance, matched-pair state bars |
| `publication/poster_section5_lineage_validation.R` | terminal publication figure | lineage marker CSVs, basal/SMG signature RDS, clinical summary CSV, optional pseudotime/spatial outputs | `ref_outs/publication/section5_lineage_validation/` | marker/signature bubbles and location association; writes placeholders for missing state-distance/Visium outputs |
| `publication/poster_section6_flot_remodelling.R` | terminal publication figure | PDO FLOT matched-response CSVs | `ref_outs/publication/section6_flot_remodelling/` | pathway response matrix and selected MP paired-expression plot |
| `publication/poster_section7_survival_targeting.R` | terminal publication figure | drug reversal overlap CSVs, optional survival outputs | `ref_outs/publication/section7_survival_targeting/` | inhibitor overlap plots; writes survival placeholder when upstream TCGA/GEO outputs absent |
| `publication/poster_requested_revisions.R` | legacy | final scRef/PDO/Visium/TME/TCGA poster outputs | `ref_outs/publication/requested_revisions/` and poster asset copies | superseded by in-place updates to `poster_section*.R`; retained for provenance only and not called by the wrapper |
| `publication/run_poster_publication_figures.sh` | PBS wrapper | all publication section scripts | `ref_outs/publication/` plus copied poster assets | re-runs canonical section replots and copies assets into the conference poster folder without overwriting manually curated schematics |
| `publication/replot_after_scenic.sh` | PBS dependency wrapper | completed `final_mp_scenic` outputs | section 2 SCENIC publication asset | submit with `afterok:<scenic_jobid>` to refresh the regulon panel automatically |

## Untracked Files

The following files or folders are intentionally not staged by agents unless the user explicitly requests it:

- `analysis/cell_states/Auto_drug_reversal/`
- `analysis/clinical/clinical_association_final_stacked.R`
- `analysis/clinical/run_clinical_association_final_boxplots.sh`
- `analysis/clinical/run_clinical_association_final_stacked.sh`
- `analysis/cnv/cna_subclone_expression_correlation.R`
- `analysis/cnv/cna_subclone_expression_correlation.sh`
- `analysis/cnv/legacy_cna_subclone_expression_visuals_v2.R`
- `analysis/cnv/cnv_malignant_subclone_mp_heatmap.sh`
- `analysis/cnv/mp_chromosomal_mapping.R`
- `analysis/spatial/map_scatlas_states_xenium.py`
- `analysis/spatial/map_scatlas_states_xenium.sh`

####################
## 22 Jun 2026 Refined MP Split Correlation Ordered Heatmap

- `analysis/metaprograms/refined_mp_split_correlation_ordered_heatmap.R`
  - Status: terminal plot-only/intermediate correlation figure.
  - Purpose: create a split-MP correlation heatmap aligned to `refined_mp_nmf_ordered_heatmap.R`, with finalized refined MP blocks retaining the same program-resolution height/width as the ordered NMF heatmap.
  - Block logic: rows/columns are expanded to original NMF-program resolution; finalized blocks follow the strict order `MP7j, MP9, MP1, MP2+, MP17, MP8+, MP10+, MP14, MP5+, MP7r, MP7v, MP10e, MP16+, MP18, MP8c, MP15c, MP12c, MP2v, MP8e, MP12a, MP13, MP7+, MP7h, MP8b, MP12b, MP15a, MP15b`; merged blocks such as `MP2+` are internally subdivided by their pre-merge `refined_submp` labels.
  - Correlation: colours are Fisher-Z averaged per-sample Spearman correlations from `refined_ucell_scores.rds`, using `EAC_Ref_epi.rds` sample labels, cached in `refined_mp_split_display_correlation_matrices.rds`.
  - Outputs: `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_split_correlation_ordered_heatmap.{pdf,png}`, final/sub-block tables, `intermediate/refined_mp_split_correlation_ordered_matrix.rds`, and `updates/new_updates/summaries/refined_mp_split_correlation_ordered_heatmap_summary.csv`.
  - Methodology: `analysis/methodology/metaprograms/refined_mp_split_correlation_ordered_heatmap_methodology.md`.
####################

####################
## 21 Jun 2026 Refined MP NMF Ordered Heatmap

- `analysis/metaprograms/refined_mp_nmf_ordered_heatmap.R`
  - Status: terminal plot-only figure.
  - Purpose: replot the original GeneNMF programme similarity matrix after MP refinement, enforcing the finalized merged refined MP order from `analysis/metaprograms/mp_refinement_merge_correlated_submps.R`.
  - Ordering: finalized refined MPs follow `MP7j, MP9, MP1, MP2+, MP17, MP8+, MP10+, MP14, MP5+, MP7r, MP7v, MP10e, MP16+, MP18, MP8c, MP15c, MP12c, MP2v, MP8e, MP12a, MP13, MP7+, MP7h, MP8b, MP12b, MP15a, MP15b`; programmes remain in original GeneNMF dendrogram order within each finalized MP block.
  - Dotted borders: diagonal dotted boxes are derived from contiguous runs of `merged_refined_mp`, so they mark finalized refined MPs rather than original nMP19 clusters.
  - Inputs: `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`, `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_assignments.rds`, and `ref_outs/Metaprogrammes_Results/mp_refinement/intermediate/merged_refined_mp_genes.rds`.
  - Outputs: `ref_outs/Metaprogrammes_Results/mp_refinement/figures/refined_mp_nmf_ordered_heatmap.{pdf,png}`, `tables/refined_mp_nmf_ordered_blocks.csv`, `intermediate/refined_mp_nmf_ordered_similarity.rds`, and `updates/new_updates/summaries/refined_mp_nmf_ordered_heatmap_summary.csv`.
  - Methodology: `analysis/methodology/metaprograms/refined_mp_nmf_ordered_heatmap_methodology.md`.
####################
## 15 May 2026 TCGA Reconstruction And Gender Validation

- `analysis/TCGA/tcga_esca_reconstruct_data.R`
  - Status: active upstream.
  - Purpose: rebuild the deleted TCGA-ESCA bulk RNA-seq input set by querying current GDC STAR-count file metadata, downloading verified open-access STAR-count TSVs, pulling cBioPortal `esca_tcga_gdc` patient/sample clinical attributes, and producing a gene-symbol TPM matrix plus harmonized metadata.
  - Inputs: GDC API TCGA-ESCA RNA-seq `STAR - Counts`; cBioPortal API study `esca_tcga_gdc`.
  - Outputs: `ref_outs/TCGA/esca_gdc_reconstruction/` with raw manifests/clinical tables/downloads, processed metadata and TPM matrix intermediates, final `TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`, run logs, and compatibility copies at `ref_outs/tcga_esca_meta.rds` and `ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`.
  - PBS wrapper: `analysis/TCGA/tcga_esca_reconstruct_data.sh` (`ncpus=2`, `mem=32gb`, `walltime=12h`, `dmtcp`, live logging).
  - Methodology: `analysis/methodology/TCGA/tcga_reconstruction_and_gender_validation_methodology.md`.

- `analysis/TCGA/tcga_gender_state_mp_validation.R`
  - Status: active terminal.
  - Purpose: validate scRef sex-associated MP and final-state signals in reconstructed TCGA EAC bulk RNA-seq by comparing Female-vs-Male effect directions between scRef sample-level UCell/state proportions and TCGA GSVA scores.
  - Inputs: reconstructed TCGA TPM/meta, `meta_full_epi.rds`, `Auto_final_states.rds`, filtered nMP19 UCell, optional 3CA UCell, nMP19 geneNMF object, and `Concise_Summary_EAC_Ref.xlsx`.
  - Outputs: `ref_outs/TCGA/gender_validation/` with cached GSVA scores, scRef/TCGA feature stats, concordance tables, and `Auto_tcga_gender_scRef_concordance.pdf/png`; compact summary at `updates/new_updates/summaries/Auto_tcga_gender_scRef_concordance_summary.csv`.
  - PBS wrapper: `analysis/TCGA/tcga_gender_state_mp_validation.sh` (`ncpus=4`, `mem=64gb`, `walltime=4h`, `dmtcp`, live logging).
  - Methodology: `analysis/methodology/TCGA/tcga_reconstruction_and_gender_validation_methodology.md`.
####################

####################
## 27 May 2026 TCGA CNA Recurrent Event Validation

- `analysis/TCGA/tcga_cna_recurrent_event_validation.R`
  - Status: active terminal.
  - Purpose: validate the top 8 scRef recurrent chromosome-arm CNA events in TCGA EAC bulk RNA-seq by deriving weighted arm-level gain/loss calls from GDC segment means, splitting patients by event presence, and testing TCGA MP/state GSVA score associations. It also discovers recurrent TCGA arm events and reports whether significant event-feature associations are also recurrent in scRef.
  - Inputs: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_scripts/TCGA/esca_tcga_gdc_segments.seg`, reconstructed TCGA TPM/meta and/or the cached gender-validation GSVA scores, scRef recurrent CNA event tables from `ref_outs/Auto_cna_subclone_expression/tables/`, `ref_outs/OAC_CNV.xlsx`, and `ref_outs/41588_2018_331_MOESM3_ESM.xlsx`.
  - Thresholding: scans absolute TCGA arm means from 0.05 to 0.30 and selects the threshold with best F1/Jaccard concordance to curated `OAC_CNV.xlsx` arm events; current run selected `0.12`.
  - Annotation correction: parses OCCAMS driver sheets by locating the real `Gene`/`hgnc_symbol` header row and normalizes display symbols while preserving raw symbols, so ST5 row 25 `ENSG00000136997 / MYC*` is retained as a high-confidence `MYC` driver on `gain_chr8q`.
  - Plotting: mirrors `analysis/cnv/cna_subclone_expression_correlation.R` presentation style: same MP/state order, separate MP/state pages, top 8 scRef recurrent events in dotplots, top 4 in boxplots, standardized event deltas, large slide-readable labels, and newer rectangular scRef-vs-TCGA MP heatmaps for recurrent events plus significant-only and all-trend continuous all-arm Spearman associations.
  - Outputs: `ref_outs/TCGA/cna_recurrent_event_validation/` with arm-call caches, reparsed CNA annotation tables, threshold sensitivity/optimization tables, scRef-validation event tests, TCGA-discovery event tests, dotplot/boxplot PDFs, rectangular MP/CNA validation heatmaps, run logs, conclusion CSV, and compact summary at `updates/new_updates/summaries/Auto_tcga_cna_recurrent_event_validation_summary.csv`.
  - Current result summary: 87 TCGA EAC primary samples matched; the top 8 scRef recurrent events had 0 TCGA MP/state associations at feature-type FDR < 0.10; TCGA-discovered events had 5 associations at feature-type FDR < 0.10, none recurrent in scRef. Conclusion: the TCGA validation supports the scRef trend that recurrent arm CNAs are not robustly coupled to MP/state programmes.
  - Methodology: `analysis/methodology/TCGA/tcga_cna_recurrent_event_validation_methodology.md`.
####################

####################
## 28 May 2026 State UMAP Dispersion And Co-Localisation

- `analysis/cell_states/state_umap_dispersion_colocalisation.R`
  - Status: active terminal/intermediate.
  - Purpose: quantify per-sample UMAP dispersion and same-label nearest-neighbour co-localisation for the five primary noreg Approach B states, excluding `Unresolved` and `Hybrid`, then repeat the analysis within `Basal to Intestinal Metaplasia` using basal-only UMAPs labelled by top basal state-defining MP (`MP17`, `MP14`, `MP5`, `MP10`, `MP8`).
  - Inputs: `ref_outs/EAC_Ref_epi.rds`, `ref_outs/Auto_topmp_v2_noreg_states_B.rds`, `ref_outs/Auto_topmp_v2_noreg_mp_adj.rds`, and `ref_outs/state_distance_pseudotime/sample_state_trajectories/*_pseudotime_states.rds`.
  - Regeneration: if the pseudotime trajectory RDS files are absent, the script reruns `analysis/cell_states/pseudotime_state_distance_matrix.R` before calculating metrics.
  - Outputs: `ref_outs/state_umap_dispersion_colocalisation/` with cached UMAP coordinates/metric RDS files in `intermediate/`, per-cell/per-sample/overall tables in `tables/`, colocalisation boxplots, dispersion boxplots, dispersion-vs-colocalisation scatter plots, and per-sample UMAP pages in `figures/`, plus compact summary `updates/new_updates/summaries/Auto_state_umap_dispersion_colocalisation_summary.csv`.
  - Cache controls: `SCREF_FORCE_REBUILD=TRUE`, `SCREF_REPLOT_ONLY=TRUE`, `SCREF_COLOCAL_K`, `SCREF_STATE_UMAP_MIN_CELLS`, and `SCREF_BASAL_UMAP_MIN_CELLS`.
  - PBS wrapper: `analysis/cell_states/state_umap_dispersion_colocalisation.sh` (`ncpus=8`, `mem=128gb`, `walltime=12h`, `dmtcp`, live logging).
  - Methodology: `analysis/methodology/cell_states/state_umap_dispersion_colocalisation_methodology.md`.
####################

####################
## 29 May 2026 scATLAS Raw Redownload And Numbat Subclone Validation

- `analysis/raw_data/Auto_download_alcindor_srr.sh` and `analysis/raw_data/Auto_download_alcindor_srr_array.sh`
  - Status: active upstream PBS downloader.
  - Purpose: redownload Alcindor 2025 raw FASTQs for `SRR27335925` through `SRR27335944` using the same `fasterq-dump --split-files --include-technical` approach as the original SRR fetch script. The array wrapper is for recovery/acceleration of unfinished accessions and uses `pigz` compression when available.
  - Outputs: `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Alcindor_2025/fastq/`.

- `analysis/raw_data/Auto_download_carroll_ega.sh`
  - Status: active upstream PBS downloader.
  - Purpose: redownload Carroll 2023 raw FASTQs from EGA dataset `EGAD00001009401` via `pyega3`; credentials are supplied by `EGA_CREDENTIAL_JSON` and are not copied into the repo.
  - Outputs: `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/Carroll_2023/ega_download/`.

- `analysis/raw_data/Auto_cellranger_alcindor_bam.sh`, `analysis/raw_data/Auto_cellranger_carroll_bam.sh`, and `analysis/raw_data/Auto_cellranger_carroll_bam_single.sh`
  - Status: active upstream PBS preprocessors.
  - Purpose: rerun Cell Ranger with `--create-bam=true` so Numbat can run SNP pileup/phasing, while preserving the original `cellranger count --id=<sample> --fastqs=<sample_dir> --transcriptome=<GRCh38-2024-A> --chemistry=auto` processing logic. Alcindor SRA FASTQs are staged as non-destructive 10x-style symlinks under `fastq_cellranger/` before Cell Ranger because raw `fasterq-dump` names are not accepted by Cell Ranger 8. Carroll samples sequenced across multiple flowcells pass a comma-separated `--sample` list of all flowcell-specific FASTQ prefixes; the single wrapper is for targeted recovery of individual failed/cancelled array elements.
  - Outputs: `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/<dataset>/cellranger/<sample>/outs/`.

- `analysis/raw_data/Auto_stage_validate_scatlas_cellranger_outputs.R` and `.sh`
  - Status: active upstream validation.
  - Purpose: copy new `filtered_feature_bc_matrix` outputs into the same `matrix_all/<sample>_filtered` structure as the historical workflow, optionally export dense count CSVs using the original `write.sh` `Read10X()`/`fwrite(as.matrix())` logic, and require exact sparse-matrix identity to the live historical matrices before Numbat.
  - Outputs: `/rds/general/project/spatialtranscriptomics/ephemeral/scRef_raw_numbat/{Alcindor_2025,Carroll_2023}/matrix_all/`, validation CSVs under `/validation/`, and optional `/00_counts_matrix_all/`.

- `analysis/raw_data/Auto_submit_scatlas_raw_rebuild.sh`
  - Status: active convenience submitter.
  - Purpose: submit EGA/SRA download jobs and dependent BAM-producing Cell Ranger arrays.
  - Methodology: `analysis/methodology/raw_data/scatlas_raw_redownload_numbat_methodology.md`.

- `analysis/cnv/Auto_scatlas_numbat_export_inputs.R`
  - Status: active upstream.
  - Purpose: export raw-barcode count matrices, cell maps, BAM paths, and barcode files for Carroll tumour and Alcindor samples.
  - Outputs: `ref_outs/Auto_scatlas_numbat/Auto_scatlas_numbat_manifest.csv` plus per-sample input folders.

- `analysis/cnv/Auto_prepare_scatlas_numbat_container.sh`, `Auto_run_scatlas_numbat_pileup.sh`, `Auto_scatlas_numbat_run_sample.R`, and `Auto_run_scatlas_numbat_sample.sh`
  - Status: active Numbat execution workflow.
  - Purpose: mirror the PDO Numbat approach using `pileup_and_phase.R` followed by `run_numbat(max_iter=2, gamma=20, init_k=3, min_cells=50)`.
  - Terminal biological statuses such as no surviving clone or no CNV after Numbat filtering are written as valid `terminal_no_subclone` sample summaries rather than rerun with looser clone-forming thresholds.
  - Outputs: per-sample allele counts and Numbat outputs under `ref_outs/Auto_scatlas_numbat/by_samples/<sample>/`.

- `analysis/cnv/Auto_scatlas_numbat_conservative_recut.R`
  - Status: active terminal/audit workflow.
  - Purpose: apply the direct conservative tree cut requested for validation, default `SCATLAS_NUMBAT_CONSERVATIVE_N_CUT=3`, with minor clone merging below `max(20 cells, 3%)`.
  - Terminal no-subclone samples are retained in the final summary as `terminal_no_subclone` and are not treated as missing tree failures.
  - Outputs: `ref_outs/Auto_scatlas_numbat/conservative_clones/`.

- `analysis/cnv/Auto_00_submit_scatlas_numbat.sh`
  - Status: active convenience submitter.
  - Purpose: submit manifest export, container preparation, per-sample pileup/Numbat jobs, and dependent conservative re-cut.
  - Methodology: `analysis/methodology/cnv/scatlas_numbat_methodology.md`.
####################

####################
## 25 Jun 2026 scATLAS RNA Velocity Workflow

- `analysis/trajectory/scatlas_velocity_metadata.R`
  - Purpose: export velocity metadata and per-sample QC barcode lists for only epithelial scATLAS samples with raw Cell Ranger BAMs under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/raw_bam_files`.
  - Inputs: `ref_outs/EAC_Ref_epi.rds`, `ref_outs/Auto_final_states.rds`, `ref_outs/Auto_topmp_v2_noreg_states_B.rds`, and raw BAM files.
  - Outputs: `ref_outs/Auto_velocity_scATLAS/tables/Auto_scatlas_velocity_cell_metadata.csv`, sample manifest, raw-BAM missing sample table, and `barcodes/*_qc_barcodes.tsv`.

- `analysis/trajectory/scatlas_velocity_prepare_refs.py`
  - Purpose: prepare velocyto GRCh38 gene GTF and hg38 RepeatMasker GTF, reusing a local RepeatMasker GTF when available before UCSC download fallback.

- `analysis/trajectory/scatlas_velocity_filter_sort.sh`, `analysis/trajectory/scatlas_velocyto_run.py`, and `analysis/trajectory/scatlas_velocity_run_velocyto.sh`
  - Purpose: filter raw Cell Ranger BAMs to epithelial QC barcodes, coordinate-sort them, and run velocyto per sample.
  - Outputs: `ref_outs/Auto_velocity_scATLAS/coord/` and `ref_outs/Auto_velocity_scATLAS/looms/`.

- `analysis/trajectory/scatlas_velocity_scvelo_visualise.py`
  - Purpose: run per-sample stochastic scVelo, cache `.h5ad` files, plot velocity streams/grids, and derive source-state to target-state directions from mean velocity alignment toward target centroids.
  - Direction states: five primary noreg Approach B states (`Classic Proliferative`, `Basal to Intestinal Metaplasia`, `SMG-like Metaplasia`, `Stress-adaptive`, `Immune Infiltrating`).
  - Outputs: `ref_outs/Auto_velocity_scATLAS/h5ad/`, scVelo cell metadata, state node CSVs, direction edge CSVs, top-direction tables, and velocity/state-direction PDFs.

- `analysis/trajectory/scatlas_velocity_nodeplots.R`
  - Purpose: terminal R nodeplots from scVelo direction CSVs, grouped across all raw-BAM scATLAS samples and by dataset.
  - Outputs: `ref_outs/Auto_velocity_scATLAS/figures/Auto_scatlas_velocity_nodeplot_by_dataset.pdf` and `updates/new_updates/summaries/Auto_scatlas_velocity_direction_summary.csv`.

- `analysis/trajectory/scatlas_velocity_submit.sh`
  - Purpose: submit the full PBS dependency chain with live logging: metadata/reference prep, per-sample filter/sort, per-sample velocyto, scVelo visualisation, and R nodeplots.
  - Methodology: `analysis/methodology/trajectory/scatlas_velocity_methodology.md`.
####################
