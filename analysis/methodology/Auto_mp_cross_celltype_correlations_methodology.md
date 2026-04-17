# Auto Cross-Celltype MP Correlation Workflow

Generated: 2026-04-13 BST

## Goal

Reproduce the cross-celltype metaprogram (MP) co-occurrence analysis on the full single-cell atlas, then annotate significant positive cross-celltype associations with curated ligand-receptor (LR) support.

This workflow is designed to mirror the logic described for Fig. 5a in the linked Nature paper, with implementation adapted to this repository and to the available scATLAS objects.

## Analysis scope

- The analysis uses the full atlas `ref_outs/EAC_Ref_merged_strict.rds`, not the epithelial-only object.
- Cancer MPs are scored only in malignant epithelial cells.
- Non-malignant MPs are scored only within their matched TME compartment.
- Correlations are computed across samples using adjusted scores, where each adjusted score is the percentage of cells in a sample that are positive for a given MP.
- LR annotation is applied only to retained positive cross-celltype edges.

## Main inputs

- `ref_outs/EAC_Ref_merged_strict.rds`
  - Full atlas containing all cell types.
- `ref_outs/meta_full_epi.rds`
  - Epithelial metadata used to define malignant cells for the cancer compartment.
- `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
  - Cancer MP definitions.
- `ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
  - Cancer MP UCell scores.
- `ref_outs/nmf_macrophage/MP_outs_default.rds`
- `ref_outs/nmf_macrophage/UCell_default.rds`
- `ref_outs/nmf_fibroblast/MP_outs_default.rds`
- `ref_outs/nmf_fibroblast/UCell_default.rds`
- `ref_outs/nmf_endothelial/MP_outs_default.rds`
- `ref_outs/nmf_endothelial/UCell_default.rds`
- `ref_outs/nmf_nk/MP_outs_default.rds`
- `ref_outs/nmf_nk/UCell_default.rds`
- `ref_outs/nmf_plasma/MP_outs_default.rds`
- `ref_outs/nmf_plasma/UCell_default.rds`
- `ref_outs/nmf_cd4/MP_outs_default.rds`
- `ref_outs/nmf_cd4/UCell_default.rds`
- `ref_outs/nmf_cd8/MP_outs_default.rds`
- `ref_outs/nmf_cd8/UCell_default.rds`
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`
  - LR catalog, using sheet `All.Pairs`.

## Compartments included

- `cancer`
- `macrophage`
- `fibroblast`
- `endothelial`
- `nk`
- `plasma`
- `cd4`
- `cd8`

## Fixed parameters used by the script

- UCell positivity cutoff: default `0.25`
- Minimum positive-sample coverage for an MP to be retained: `> 5` samples
- Minimum shared samples per study for a celltype pair to be eligible: `10`
- Minimum number of eligible shared samples for a tested MP pair: `3`
- Positive edge threshold: Pearson `> 0`, Spearman `P < 0.05`, and `-log10(P) >= 4`
- Negative edge threshold: Pearson `< 0` and Spearman `P < 0.05`
- Top genes used for LR annotation: top `4,000` ranked genes per node
- LR evidence retained: only `literature supported` and `putative`
- LR evidence excluded: all rows with `Pair.Evidence` beginning `EXCLUDED`

## Why the UCell cutoff is 0.25

The original paper used a score threshold of `> 1` with `sigScores` from `scalop`, which is not directly comparable to the UCell score range used here. In this repository, the non-malignant UCell distributions are much narrower, and a `0.5` cutoff was too sparse for several compartments. The script therefore defaults to `0.25` and writes cutoff-sensitivity diagnostics so this choice can be inspected and retuned if needed.

## Workflow

### 1. Resolve project paths and create output directories

The script first resolves the active project root from the two expected HPC mount styles:

- `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline`
- `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline`

It then sets the working directory to `ref_outs/` and creates:

- `ref_outs/non_malignant_mp_correlations/`
- `ref_outs/non_malignant_mp_correlations/cache/`
- `updates/new_updates/summaries/`

### 2. Define display metadata and cancer MP descriptions

The script registers the eight compartments, their atlas celltype labels, input paths, and plotting colours.

For cancer MPs, user-facing plot and export labels are converted from raw MP IDs to descriptive labels using the internal lookup table:

- `MP1` -> `G2M Cell Cycle`
- `MP2` -> `MYC-related Proliferation`
- `MP5` -> `Epithelial IFN Resp.`
- `MP7` -> `DNA Damage Repair`
- `MP8` -> `Intestinal Diff.`
- `MP9` -> `G1S Cell Cycle`
- `MP10` -> `Columnar Diff.`
- `MP12` -> `Neuro-responsive Epi.`
- `MP13` -> `Hypoxic Inflam. Epi.`
- `MP14` -> `Hypoxia Adapted Epi.`
- `MP15` -> `Immune Infiltration`
- `MP16` -> `Secretory Diff. (Gastric)`
- `MP17` -> `Basal-like Transition`
- `MP18` -> `Secretory Diff. (Intest.)`

These descriptions are used wherever cancer MPs appear in node labels, network labels, and LR exports.

For cancer MP ordering in plots, the script follows the same logic used in `analysis/cell_states/states_unresolved_relabel.R`:

- all cancer cell-cycle MPs are grouped first
- the remaining cancer MPs then follow the state order
- this order is carried through the bubble plot pages and node-ordered outputs

### 3. Load the full atlas and derive per-cell sample metadata

The script loads:

- the full merged atlas from `EAC_Ref_merged_strict.rds`
- epithelial metadata from `meta_full_epi.rds`

It then creates a simple sample-level metadata table per cell containing:

- cell barcode
- `orig.ident` sample ID
- study label

If `study` is missing, the script derives it from the first two underscore-delimited tokens of the sample name.

### 4. Define cancer cells and matched non-malignant compartments

The cancer compartment is not taken as all epithelial cells. It is explicitly restricted to epithelial cells with:

- `malignancy == "malignant_level_1"` or
- `malignancy == "malignant_level_2"`

Non-malignant compartments are matched to the atlas using `celltype_update` in the full atlas:

- macrophage -> `macrophage`
- fibroblast -> `fibroblast`
- endothelial -> `endothelial`
- nk -> `nk.cell`
- plasma -> `plasma`
- cd4 -> `t.cell` with subtype restriction to `cd4`
- cd8 -> `t.cell` with subtype restriction to `cd8`

This ensures the analysis is built from the complete sample context, while still scoring each compartment only within its biologically relevant cells.

### 5. Load MP definitions and apply the silhouette filter

For each compartment, the script loads the MP object and corresponding UCell score matrix.

Before any downstream analysis, it applies the repository-standard silhouette filter:

- MPs with silhouette `< 0` are discarded.
- Only retained MPs are taken forward.

This is done for both cancer and all non-malignant compartments.

### 6. Match score matrices back to the full atlas

For each compartment, the script intersects:

- cells in the UCell score matrix
- cells available in the atlas
- cells that pass the compartment-specific filtering rule

This produces a compartment-specific cell set that is then used consistently for:

- adjusted score calculation
- sample coverage summaries
- positive-cell extraction for LR annotation

### 7. Calculate adjusted scores for each MP in each sample

For each compartment separately:

1. A cell is declared MP-positive when `UCell > cutoff`.
2. Cells are grouped by `orig.ident`.
3. For each sample and MP, the adjusted score is calculated as:

`100 x (number of MP-positive cells in that sample and compartment) / (total number of cells of that compartment in that sample)`

The script saves one adjusted-score CSV per compartment. These tables are the direct sample-level inputs to the correlation analysis.

### 8. Calculate MP coverage and cutoff sensitivity

For every retained MP, the script calculates:

- the number of samples with adjusted score `> 0`
- the percentage of samples with adjusted score `> 0`
- the percentage of cells positive at the current cutoff

The script repeats this across the cutoff grid:

- `0.10`
- `0.15`
- `0.20`
- `0.25`
- `0.30`
- `0.35`
- `0.40`
- `0.50`
- plus any user-provided cutoff

This produces:

- `Auto_cross_celltype_cutoff_sensitivity.csv`
- `Auto_cross_celltype_cutoff_sensitivity.pdf`

The purpose is to document how sensitive sample coverage is to the positivity threshold.

### 9. Retain only MPs with non-trivial sample coverage

After adjusted-score construction, an MP is retained for downstream pairwise testing only if it is positive in more than 5 samples.

This removes extremely sparse MPs that would otherwise create unstable or uninformative correlations.

### 10. Build the node table used across all downstream outputs

Each retained compartment-MP combination becomes one network node. The node summary stores:

- compartment
- celltype display name
- MP ID
- MP display name
- node label
- sample coverage counts and percentages
- mean positive-cell fraction
- coverage pass/fail

Cancer nodes use descriptive MP labels instead of raw MP IDs in the node label.

### 11. Define eligible sample overlap between every pair of different cell types

The script tests all pairwise combinations of different compartments.

For each celltype pair:

1. It finds samples present in both compartments.
2. It maps those samples to studies.
3. It counts the number of shared samples per study.
4. It retains only studies with at least 10 shared samples.
5. It keeps only samples from those eligible studies.

This reproduces the intended rule that cross-celltype associations should be evaluated only where there is sufficient matched tumour coverage.

### 12. Compute Pearson and Spearman correlations for every eligible MP pair

Within each eligible celltype pair, the script tests every retained MP from compartment A against every retained MP from compartment B.

For each MP pair, it computes across the eligible shared samples:

- Pearson correlation coefficient and P value
- Spearman correlation coefficient and P value
- `-log10(Spearman P value)`

The resulting master table is saved as:

- `Auto_cross_celltype_correlations_all.csv`

The script also writes:

- `Auto_cross_celltype_correlations_positive.csv`
- `Auto_cross_celltype_correlations_negative.csv`
- `Auto_cross_celltype_shared_sample_summary.csv`

### 13. Define positive and negative edge sets

Positive edges are defined as:

- Pearson `r > 0`
- Spearman `P < 0.05`
- `-log10(Spearman P) >= 4`

Negative edges are defined as:

- Pearson `r < 0`
- Spearman `P < 0.05`

This follows the logic described in the paper: stronger significance filtering for positive co-occurrence edges, while retaining all statistically significant negative edges.

### 14. Generate cross-celltype overview plots

The script writes:

- `Auto_cross_celltype_correlation_bubble.pdf`
  - bubble size = `-log10(Spearman P)`
  - bubble fill = Pearson `r`
  - one page per focal cell type
  - within each page, panels are split by celltype pair without bidirectional duplication
  - focal page order is `cancer`, `fibroblast`, `endothelial`, `cd8`, `cd4`, `macrophage`, `nk`
- `Auto_cross_celltype_positive_network.pdf`
- `Auto_cross_celltype_negative_network.pdf`

These plots use the node labels built earlier. Cancer nodes display the cancer MP description plus the `cancer` suffix. Celltype colours follow the publication UMAP palette used in `analysis/plotting/publication_umap.R`.

### 15. Load and clean the LR reference

The script loads the `All.Pairs` sheet from:

- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`

It standardizes the LR columns and derives:

- ligand symbol
- receptor symbol
- pair name
- pair source
- pair evidence
- support database flags
- support database count

Critically, it retains only:

- `literature supported`
- `putative`

It excludes all rows labelled:

- `EXCLUDED`
- `EXCLUDED not ligand`
- `EXCLUDED not receptor`
- any other `EXCLUDED*` variant

The LR status table records both total rows and retained rows after evidence filtering.

### 16. Extract positive cells for each retained positive edge

LR annotation is applied only to positive edges.

For each retained positive edge:

1. The script recovers the eligible studies for that edge.
2. It derives the corresponding eligible shared samples.
3. For each node, it identifies the cells in those shared samples with `UCell > cutoff` for that MP.
4. If either node has zero positive cells, the edge receives no LR annotation.

This ensures the LR step is restricted to the exact cells that contributed to the observed cross-sample association.

### 17. Rank expressed genes in the contributing cells

For each positive node:

- RNA expression is taken from the full atlas `RNA` assay.
- Mean expression is computed across the positive cells only.
- Genes are ranked by mean expression.
- The top 4,000 genes are retained.

This produces a ranked expressed-gene list representing the active transcriptional context of each node.

### 18. Match LR pairs in both biologically relevant directions

For each positive edge, the script checks both directions:

1. Driver top genes as ligands against target MP genes as receptors.
2. Driver top genes as receptors against target MP genes as ligands.

This is done for both node1 -> node2 and node2 -> node1, so each edge is interrogated symmetrically.

For every retained LR match, the output stores:

- ligand
- receptor
- pair name
- pair evidence
- pair source
- support database list
- support database count
- match mode
- edge identity
- node labels
- celltype pair
- Pearson and Spearman statistics
- positive cell counts for driver and target nodes

### 19. Summarise LR support at the edge level

The script aggregates LR rows into edge-level summaries that capture:

- total candidate LR pairs
- number of literature-supported pairs
- number of putative pairs
- maximum support database count
- top example retained pairs

These summaries are joined back onto the positive edge table and saved as:

- `Auto_cross_celltype_ligand_receptor_edge_summary.csv`
- `Auto_cross_celltype_positive_edges_lr_annotated.csv`

The raw pair-level table is saved as:

- `Auto_cross_celltype_ligand_receptor_pairs.csv`

### 20. Generate the LR-annotated positive network

When LR-supported edges are present, the script produces:

- `Auto_cross_celltype_positive_network_lr_annotated.pdf`

This is the positive network with top LR labels added to the strongest annotated edges.

It also writes:

- `Auto_cross_celltype_ligand_receptor_summary.pdf`

This summary PDF contains:

- a barplot of positive edges ranked by number of retained LR pairs
- a barplot of the most recurrent retained LR pairs across edges

### 21. Export the formatted focal-celltype Excel workbook

The script generates one formatted Excel workbook.

- `Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx`

Structure:

- one `Overview` sheet summarizing LR-supported edges
- one sheet per focal cell type without bidirectional duplication

Within each focal sheet:

- rows are grouped by partner cell type
- for example, in the `cancer` sheet, cancer-fibroblast rows appear together, then cancer-endothelial rows, then other partner compartments
- rows are then sorted by edge significance and LR support
- the `Direction` column is encoded as `1` when the focal cell type is the ligand side and `0` when the focal cell type is the receptor side
- database-source bookkeeping columns are intentionally omitted; only the retained evidence class (`literature supported` or `putative`) is shown

This workbook is the main review file when the reader wants to follow one focal compartment at a time.

The workbook is styled with:

- formatted headers
- frozen panes
- numeric formatting
- wrapped text in annotation fields
- colour scaling on `-log10 Spearman p`
- tab colours matched to the focal cell type where applicable

### 22. Write a final run summary

The script writes:

- `updates/new_updates/summaries/Auto_mp_cross_celltype_correlations_summary.csv`

This summary records:

- cache version
- active cutoff
- number of retained nodes
- number of pairwise tests
- number of positive and negative edges
- LR catalog retention statistics
- workbook paths

## Cache behaviour

Large intermediate steps are cached under:

- `ref_outs/non_malignant_mp_correlations/cache/`

The cache files are:

- `Auto_cross_celltype_step1_compartment_cache.rds`
- `Auto_cross_celltype_step2_correlation_cache.rds`
- `Auto_cross_celltype_step3_lr_cache.rds`

Each cache stores the current `cache_version`. On rerun:

- if the cache exists and the version matches, it is loaded
- otherwise, that step is rebuilt and the cache is overwritten

To force a full rebuild, set:

```bash
AUTO_MPXCELL_FORCE_REBUILD=TRUE
```

## Run commands

### Interactive run

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
AUTO_MPXCELL_FORCE_REBUILD=TRUE Rscript analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.R
```

### Interactive run with a custom UCell cutoff

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
AUTO_MPXCELL_FORCE_REBUILD=TRUE Rscript analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.R 0.20
```

### PBS submission

```bash
/opt/pbs/bin/qsub -v AUTO_MPXCELL_FORCE_REBUILD=TRUE analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.sh
```

### PBS submission with a custom cutoff

```bash
/opt/pbs/bin/qsub -v AUTO_MPXCELL_FORCE_REBUILD=TRUE,cutoff=0.20 analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.sh
```

## Expected outputs

All main outputs are written under:

- `ref_outs/non_malignant_mp_correlations/`

Key files:

- per-compartment adjusted-score CSVs
- `Auto_cross_celltype_node_summary.csv`
- `Auto_cross_celltype_cutoff_sensitivity.csv`
- `Auto_cross_celltype_cutoff_sensitivity.pdf`
- `Auto_cross_celltype_shared_sample_summary.csv`
- `Auto_cross_celltype_correlations_all.csv`
- `Auto_cross_celltype_correlations_positive.csv`
- `Auto_cross_celltype_correlations_negative.csv`
- `Auto_cross_celltype_correlation_bubble.pdf`
- `Auto_cross_celltype_positive_network.pdf`
- `Auto_cross_celltype_negative_network.pdf`
- `Auto_cross_celltype_ligand_receptor_status.csv`
- `Auto_cross_celltype_ligand_receptor_pairs.csv`
- `Auto_cross_celltype_ligand_receptor_edge_summary.csv`
- `Auto_cross_celltype_positive_edges_lr_annotated.csv`
- `Auto_cross_celltype_positive_network_lr_annotated.pdf`
- `Auto_cross_celltype_ligand_receptor_summary.pdf`
- `Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx`

## Reproducibility notes

- The analysis is sample-level, not cell-level, at the correlation stage.
- Cancer is scored only in malignant epithelial cells, even though the full atlas is used for sample matching and expression extraction.
- Non-malignant compartments are evaluated in their own matched celltype only.
- LR annotation is a support layer on top of the positive correlation graph; it does not alter edge significance.
- Old outputs should not be trusted after methodology changes unless the cache was rebuilt with the current `cache_version`.
