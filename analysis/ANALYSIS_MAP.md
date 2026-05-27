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
| `cnv/mp_chromosomal_mapping.R` | untracked terminal | MP gene positions | MP chromosomal mapping figures | terminal; currently unstaged |
| `cnv/cnv_profiling.R` | core CNV | epithelial per-sample RDS | inferCNA profiles | CNV subsetting/plotting |
| `cnv/cnv_subsetting.R` | core CNV | inferCNA profiles | filtered CNA matrices | CNV plotting |
| `cnv/cnv_plotting.R` | core CNV terminal | CNA matrices | publication CNV plots | terminal |
| `developmental/developmental.R` | reference build | developmental xlsx files | `enrich_dev.rds`, per-stage RDS | enrichment/reference interpretation |
| `developmental/external_epi_mp_ucell_heatmap.R` | active terminal/upstream | MP genes, adult oesophagus/stomach, Barretts references | external epithelial MP UCell heatmaps/cache | terminal |
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
