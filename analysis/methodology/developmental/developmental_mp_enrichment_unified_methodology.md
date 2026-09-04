# Unified Developmental MP Enrichment Methodology

`analysis/developmental/legacy_developmental_mp_enrichment_unified.R` records the historical nMP19 validation against developmental reference databases. It is retained for provenance; centred step 04 is the current workflow.
The runtime computational core is in `analysis/developmental/developmental_mp_enrichment_unified_original_aligned_core.R`; this was split out so the original-script-aligned logic is explicit and auditable.

## Inputs

- scATLAS MP definitions: `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- scATLAS MP UCell scores: `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- epithelial Seurat atlas: `ref_outs/EAC_Ref_epi.rds`
- developmental marker workbooks under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/`
- cached stomach epithelial marker table: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach_epi_DGE_ordered_pretty.csv`
- original enrichment references under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/`
- optional annotated external expression objects listed in `Auto_developmental_reference_expression_availability.csv`

## Ranked Reference Construction

The script rebuilds a ranked long table of developmental reference genes before scoring.

- Early embryogenesis uses the workbook row order within each lineage after filtering `p_val_adj < 0.05`; the sheet includes `avg_log2FC` and `p_val_adj`.
- Organogenesis uses `S1D`, where marker columns are explicitly described as DEGs ordered by z-score. Subtype terms use marker-column order. Major terms aggregate subtype markers by the best marker-column rank.
- Normal development long markers use `Table_S4`, filtered to `qval < 0.05` and ranked within mapped cell type by descending `fold.change`, then `qval`.
- Normal development short markers use `Table_S3` literature marker order. These are not differential ranked lists; the rank audit records this explicitly.
- Adult oesophagus uses source sheet order after `p_val_adj < 0.05`. Adult stomach uses `Stomach_epi_DGE_ordered_pretty.csv`, already ordered by marker strength from the cached epithelial differential-expression table.
- Barretts oesophagus uses source sheet order after FDR filtering, with per-sheet logFC values retained for audit.

The ranked reference table is cached at `intermediate/Auto_developmental_ranked_references.rds`; the audit table is written to `tables/Auto_developmental_reference_rank_audit.csv`.

## Method 1: Overlap Enrichment

Overlap enrichment now follows `analysis/enrichment/enrichment_annotation.R` exactly. For each MP and developmental reference, it calls `clusterProfiler::enricher()` with the same `TERM2GENE`, `TERM2NAME`, `pAdjustMethod = "BH"`, and `qvalueCutoff = 0.05` pattern as the original script.

- `all`: uses the original per-stage `TERM2GENE`/`TERM2NAME` RDS objects directly, so the all-gene overlap matrices validate exactly against `ref_outs/cluster_enrich.rds`.
- `top50`: constructs temporary per-stage `TERM2GENE` objects from the first 50 ranked genes per term, then runs the same `clusterProfiler::enricher()` logic.

The heatmap table construction mirrors `enrich_heatmap()` from the original script: all custom reference terms are retained in `TERM2NAME` order, values are `pmin(-log10(p.adjust), 7)` with missing results set to zero, and cell text is the original `GeneRatio` overlap string. `tables/Auto_developmental_validation_overlap_all_vs_original.csv` records the exact all-gene validation.

## Method 2: Expression Correlation

Developmental term gene sets are UCell-scored in the scATLAS epithelial Seurat object. For each sample (`orig.ident`), the script calculates Spearman correlations between developmental term UCell scores and MP UCell scores, then Fisher-z transforms the per-sample correlations and tests whether the mean differs from zero.

This follows `analysis/metaprograms/mp_database_correlation.R` exactly:

- custom term gene lists for all-gene mode are extracted from `er@geneSets` after the original enrichment object is loaded;
- reference terms are flattened with `db__make.names(term)` names, scored with `AddModuleScore_UCell(..., name = "", maxRank = max(5000, max(lengths(gene_sets)) + 100))`, and correlated against the filtered MP UCell matrix;
- sample grouping uses `orig.ident`;
- per-sample Spearman rho values are Fisher-z transformed, tested with `t.test()`, and BH-adjusted within each developmental database.

All-gene mode reuses `ref_outs/UCell_ref_terms_v2_MP19.rds` when available, matching the original cache. `tables/Auto_developmental_validation_correlation_all_vs_original.csv` records zero-difference validation against `ref_outs/Auto_MP_correlation_results_v2_MP19.rds`. Top50 mode uses the same correlation machinery with the top50 gene lists.

## Method 3: External Annotated Reference Scoring

The script checks all requested developmental references for annotated expression matrices and scores those that are available with the same helper logic as `analysis/developmental/external_epi_mp_ucell_heatmap.R`: create a Seurat object from counts, score filtered scATLAS MPs with `AddModuleScore_UCell(..., slot = "counts", name = "_UCell", BPPARAM = SerialParam(), ncores = 1)`, and summarise mean UCell per reference cell type after per-celltype downsampling.

Downloaded or used annotated references:

- Early embryogenesis: downloaded `psd.R3.6.em.seurat.ob.rds` from the paper dataset portal into `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/early_embryogenesis/` and scored `rename_EML` cell-type annotations after mapping `2-4 cell` to the reference `Z4cell` term.
- Organogenesis: downloaded GEO `GSE157329_cell_annotate.txt.gz`, `GSE157329_gene_annotate.txt.gz`, and `GSE157329_raw_counts.mtx.gz` into `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/`; scored major developmental-system and endoderm/subtype annotations separately.
- Normal development stomach: downloaded public Descartes/Fred Hutch `Stomach_gene_count.RDS`, `df_cell.RDS`, and `df_gene.RDS` into `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/`. Only stomach cells are scored. `df_cell$Main_cluster_name` supplies the 16 stomach annotations used by both long and short normal-development enrichment references, and `df_gene$gene_short_name` maps Ensembl row IDs to symbols.
- Adult oesophagus: downloaded Broad Single Cell Portal study `SCP1242` into `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/` using the temporary authorized curl config supplied by the user. The generated `cfg.txt` is retained for provenance. The workflow scores the raw `EoE.mtx` MatrixMarket file after streaming a per-celltype sampled subset, using `EoE_meta.txt$cell_type_anno`.
- Adult stomach: scored the local annotated Seurat object `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds`.
- Barretts oesophagus: downloaded direct public high-quality combined RDS `https://cellgeni.cog.sanger.ac.uk/esophagus-cancer/alldatahighquality.rds` into `ref_outs/Auto_developmental_mp_enrichment_unified/downloads/barretts/` and scored `cell_type_secondary`.

Below records download commands, access requirements, and local paths. The output is a mean MP UCell matrix per scored reference cell type.

| Reference | Dataset Source | Dataset | Paper URL | Download or Repository | Download Command | Access Requirement | Local Path | Available |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Early_Embryogenesis | Early_Embryogenesis_integrated | early_embryogenesis | https://pmc.ncbi.nlm.nih.gov/articles/PMC11725501/ | paper/data portal downloaded Seurat RDS to downloads/early_embryogenesis | downloaded from paper dataset portal as psd.R3.6.em.seurat.ob.rds | direct/local cached RDS; original paper data portal | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/early_embryogenesis/psd.R3.6.em.seurat.ob.rds | TRUE |
| Organogenesis_major | Organogenesis_GSE157329_major | organogenesis_4_6wk | https://www.nature.com/articles/s41556-023-01108-w | GEO GSE157329 and VisCello https://heoa.shinyapps.io/base/ | downloaded GEO files GSE157329_cell_annotate.txt.gz, GSE157329_gene_annotate.txt.gz, GSE157329_raw_counts.mtx.gz | direct public GEO/VisCello files | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_cell_annotate.txt.gz;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_gene_annotate.txt.gz;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_raw_counts.mtx.gz | TRUE |
| Organogenesis_sub | Organogenesis_GSE157329_sub | organogenesis_4_6wk | https://www.nature.com/articles/s41556-023-01108-w | GEO GSE157329 and VisCello https://heoa.shinyapps.io/base/ | downloaded GEO files GSE157329_cell_annotate.txt.gz, GSE157329_gene_annotate.txt.gz, GSE157329_raw_counts.mtx.gz | direct public GEO/VisCello files | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_cell_annotate.txt.gz;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_gene_annotate.txt.gz;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/organogenesis/GSE157329_raw_counts.mtx.gz | TRUE |
| Normal_Development_long | Normal_Development_Stomach_long | normal_development_stomach_10_18wk | https://www.science.org/doi/10.1126/science.aba7721 | Descartes/Fred Hutch direct RDS: Stomach_gene_count.RDS, df_cell.RDS, df_gene.RDS | wget -c https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/Stomach_gene_count.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS | direct public atlas download; no login | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/Stomach_gene_count.RDS;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_cell.RDS;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_gene.RDS | TRUE |
| Normal_Development_short | Normal_Development_Stomach_short | normal_development_stomach_10_18wk | https://www.science.org/doi/10.1126/science.aba7721 | Descartes/Fred Hutch direct RDS: Stomach_gene_count.RDS, df_cell.RDS, df_gene.RDS | wget -c https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/Stomach_gene_count.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_cell.RDS https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/downloads/data_summarize_fetus_data/df_gene.RDS | direct public atlas download; no login | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/Stomach_gene_count.RDS;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_cell.RDS;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/normal_development/df_gene.RDS | TRUE |
| Adult_Epithelium | Adult_Oesophagus | adult_oesophagus | https://www.nature.com/articles/s41467-024-47647-0 | Broad Single Cell Portal SCP1242 temporary authorized bulk-download curl config | curl Broad SCP1242 generate_curl_config URL with temporary auth_code to cfg.txt; curl -K cfg.txt | temporary Broad auth_code supplied by user; no interactive login during this run | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/SCP1242/metadata/EoE_meta.txt;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/SCP1242/expression/63f53992d91a88956d36dc4f/EoE_cell.tsv;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/SCP1242/expression/63f53992d91a88956d36dc4f/EoE_gene.tsv;<br>/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/adult_oesophagus/SCP1242/expression/63f53992d91a88956d36dc4f/EoE.mtx | TRUE |
| Adult_Epithelium | Adult_Stomach | adult_stomach | https://www.sciencedirect.com/science/article/pii/S2211124723012482 | local annotated Seurat object found under developmental reference directory | local file /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds; no download in this script | local project file; no login | /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/Stomach.rds | TRUE |
| Barretts_Oesophagus | Barretts_HighQuality | barretts_high_quality | https://www.science.org/doi/10.1126/science.abd1449 | Esophagus Cancer Cell Atlas direct high-quality combined RDS download | wget/curl https://cellgeni.cog.sanger.ac.uk/esophagus-cancer/alldatahighquality.rds | direct public atlas download; no login | /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_developmental_mp_enrichment_unified/downloads/barretts/alldatahighquality.rds | TRUE |

## Outputs

- `figures/Auto_developmental_mp_unified_top50.pdf`: one page per developmental reference. Pages include overlap and expression-correlation heatmaps, plus external reference-celltype UCell scoring when an annotated expression object is available.
- `figures/Auto_developmental_mp_top50_vs_all_overlap_correlation.pdf`: separate overlap and correlation comparison pages per developmental reference, comparing top-50 versus all-gene inputs for methods 1 and 2. Rows marked `*` have fewer than 50 all-gene markers, so top50 equals all.
- `tables/Auto_developmental_overlap_<top50|all>.csv`
- `tables/Auto_developmental_expression_correlation_<top50|all>.csv`
- `tables/Auto_developmental_external_mp_ucell_summary.csv`
- `tables/Auto_developmental_reference_celltype_coverage.csv`
- `tables/Auto_developmental_top50_equals_all_terms.csv`
- `tables/Auto_developmental_validation_overlap_all_vs_original.csv`
- `tables/Auto_developmental_validation_correlation_all_vs_original.csv`
- `logs/Auto_developmental_mp_enrichment_run_summary.txt`

## Rebuild Controls

- `SCREF_FORCE_REBUILD=TRUE`: rebuild cached reference UCell scores and external scores.
- `SCREF_REPLOT_ONLY=TRUE`: reuse existing tables and redraw PDFs.
- `SCREF_UCELL_CORES`: controls UCell cores for scATLAS reference scoring.
- `SCREF_MAX_CELLS_PER_TYPE`: controls downsampling for external annotated reference scoring.
