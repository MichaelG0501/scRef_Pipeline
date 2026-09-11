# TCGA-STAD Centred MP/State Location GSVA

## Purpose

Quantify the updated centred metaprogrammes and the Basal to intestinal metaplasia and SMG to intestinal metaplasia state gene-set expression across TCGA-STAD primary-tumour anatomical locations.

## Inputs

- TCGA-STAD gene-symbol TPM matrix: `EAC_Ref_all/00_merged/stomach_bulk/matrices/Auto_tcga_stad_tpm_matrix_gene_symbol.rds`.
- TCGA-STAD per-sample metadata: `EAC_Ref_all/00_merged/stomach_bulk/metadata/Auto_tcga_stad_sample_metadata.rds`.
- Current final centred MP genes: `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`.

## Gene sets and scoring

The final state grouping is copied from `analysis/metaprograms/centred/tcga_mp_survival_volcano_centred.R`.

- Basal to intestinal metaplasia: `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, and `MP10+`.
- SMG to intestinal metaplasia: `MP8+`, `MP8b`, `MP16`, `MP18b`, and `MP17`.

Bulk TPM values are transformed as `log2(TPM + 1)`. Each individual MP is scored with GSVA. Each state is also scored with GSVA using the union of its constituent final MP genes, avoiding duplicated genes receiving extra weight.

## Cohort and tests

Only TCGA sample-type code `01` (primary tumour) is analysed, so anatomical differences are not confounded by solid-tissue normal samples. The raw stomach location categories and their counts are exported. Groups with fewer than 10 primary tumours are retained in QC tables but excluded from plots and inferential tests; this prevents a single rare site from driving a nominal result.

For every MP and state score, a Kruskal-Wallis test compares eligible locations. Benjamini-Hochberg correction is applied across all 12 MP tests and, separately, across the two state tests. Pairwise Wilcoxon tests with BH correction are exported as supporting results.

## Cardia versus distal non-cardia comparison

The same GSVA scores are also compared between proximal cardia (`Cardia, NOS`) and resolved distal non-cardia tumours (gastric antrum, body of stomach, and fundus of stomach). `Stomach, NOS` and lesser-curvature NOS are not forced into either class and are exported as excluded from this binary comparison. Wilcoxon rank-sum tests with BH correction are performed across the 12 MPs and, separately, across the two state-union scores.

## Outputs

All outputs are under `ref_outs/Metaprogrammes_Results/centred/bulk_tcga_stad_location_gsva/`:

- `intermediate/`: GSVA score matrix.
- `tables/`: sample scores, location counts, gene-set match counts, global tests, pairwise tests, and summary.
- `figures/`: MP and state score boxplots by stomach location.
- `logs/`: compact run summary.
