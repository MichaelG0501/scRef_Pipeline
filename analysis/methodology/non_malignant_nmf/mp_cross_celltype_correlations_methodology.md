# Auto MP and Cancer-State Cross-Celltype Correlation Methodology

Generated: 2026-04-17 BST

## 1. Goal and Scope

This workflow quantifies sample-level co-occurrence between cancer programs or finalized cancer states and non-malignant metaprograms (MPs), then annotates significant positive associations with curated ligand-receptor (LR) support.

The analysis is built on the full atlas, not the epithelial subset alone. Cancer is always restricted to malignant epithelial cells, while non-malignant MPs are evaluated only within their matched cell type.

The script now runs four related analysis modes in one pass:

1. `01_cancer_mps_cross_only`
   - Cancer represented by malignant-epithelial MP-positive cells.
   - Only cross-celltype correlations are tested.
2. `02_cancer_mps_cross_and_within`
   - Cancer represented by malignant-epithelial MP-positive cells.
   - Both cross-celltype and within-celltype correlations are tested.
3. `03_cancer_states_cross_only`
   - Cancer represented by finalized cancer-state labels.
   - Only cross-celltype correlations are tested.
4. `04_cancer_states_cross_and_within`
   - Cancer represented by finalized cancer-state labels.
   - Both cross-celltype and within-celltype correlations are tested.

Each mode writes to its own subfolder under `ref_outs/non_malignant_mp_correlations/`.

---

## 2. Core Inputs

- `ref_outs/EAC_Ref_merged_strict.rds`
  - Full scRNA-seq atlas containing all cell types and `celltype_update`.
- `ref_outs/meta_full_epi.rds`
  - Epithelial metadata used to define malignant epithelial cells.
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
- `ref_outs/Auto_final_states.rds`
  - Finalized per-cell cancer-state labels.
- `ref_outs/Auto_six_state_markers/Auto_six_state_markers_ranked.csv`
  - Ranked recurrent state-marker table used to define cancer-state gene sets for LR matching.
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`
  - Ligand-receptor reference catalog, using sheet `All.Pairs`.

---

## 3. Compartments Included

The workflow evaluates eight compartments:

- `cancer`
- `fibroblast`
- `endothelial`
- `cd8`
- `cd4`
- `macrophage`
- `nk`
- `plasma`

The display and page order is:

- `cancer`
- `fibroblast`
- `endothelial`
- `cd8`
- `cd4`
- `macrophage`
- `nk`
- `plasma`

---

## 4. Cancer and Non-Malignant Cell Definitions

### 4.1 Cancer cells

Cancer is never defined as all epithelial cells. It is restricted to epithelial cells with:

- `malignancy == "malignant_level_1"`, or
- `malignancy == "malignant_level_2"`

This malignant epithelial subset is the cancer denominator in both MP-based and state-based modes.

### 4.2 Non-malignant compartments

Non-malignant cells are taken from the full atlas using `celltype_update`, with T-cell subtype splitting where needed:

- fibroblast -> `fibroblast`
- endothelial -> `endothelial`
- macrophage -> `macrophage`
- nk -> `nk.cell`
- plasma -> `plasma`
- cd4 -> `t.cell` plus `cd4` subtype restriction
- cd8 -> `t.cell` plus `cd8` subtype restriction

This ensures that all sample matching uses the full atlas while each MP score is evaluated only in its matched cell class.

---

## 5. Cancer Representation in the Two Cancer Modes

### 5.1 Cancer MP modes

For the cancer MP modes, cancer cells are scored using the silhouette-filtered epithelial MP UCell matrix. User-facing cancer labels use descriptive cancer MP names, and the cancer MP plot order follows the same ordering logic used in `analysis/cell_states/final_state_unresolved_relabel.R`:

- all cell-cycle cancer MPs first
- then the remaining cancer MPs in the finalized state-linked order

### 5.2 Cancer state modes

For the cancer state modes, cancer cells are represented by finalized labels from `Auto_final_states.rds`.

The retained final states are:

- `Classic Proliferative`
- `Basal to Intestinal Metaplasia`
- `Stress-adaptive`
- `SMG-like Metaplasia`
- `Immune Infiltrating`
- `3CA_EMT_and_Protein_maturation`

This branch does not apply any UCell threshold to cancer-state positivity. A cancer cell is positive for a state if and only if it carries that finalized state label. The sample-level cancer-state adjusted score is therefore the percentage of malignant epithelial cells in a sample assigned to that state.

---

## 6. Common Fixed Parameters

These values are hard-coded in the current script:

- Non-cancer UCell positivity cutoff: `0.25`
- Minimum positive-sample coverage for a node to pass the coverage filter: `> 5` samples
- Minimum shared samples per study for a celltype pair to be eligible: `10`
- Minimum total eligible shared samples for a tested node pair: `3`
- Positive edge threshold: Pearson `> 0`, Spearman `P < 0.05`, and `-log10(P) >= 4`
- Negative edge threshold: Pearson `< 0` and Spearman `P < 0.05`
- Top ranked genes used for LR annotation: `4,000`
- Retained LR evidence classes: `literature supported`, `putative`
- Removed LR evidence classes: all `EXCLUDED*`
- Cancer-state marker genes retained for LR target matching: top `100` ranked genes per finalized state

---

## 7. Step-by-Step Workflow

### 7.1 Resolve paths and create output structure

The script resolves the project root from either HPC mount style:

- `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline`
- `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline`

It then sets the working directory to `ref_outs/` and creates:

- `ref_outs/non_malignant_mp_correlations/`
- one subfolder per analysis mode
- one cache directory per mode
- `updates/new_updates/summaries/`

### 7.2 Load the full atlas, malignant metadata, and state references

The script loads the full merged atlas once and extracts:

- per-cell metadata
- the RNA expression matrix used later for LR ranking
- the malignant epithelial cell set

For state modes it also loads:

- finalized per-cell state labels
- the ranked six-state marker table
- a state-specific gene set per finalized state, taking the top ranked recurrent markers after state-aware filtering

### 7.3 Load MP definitions and apply the silhouette filter

For every MP-scored compartment, the script loads the MP object and corresponding UCell score matrix. It applies the repository-standard silhouette filter before any downstream use:

- MPs with silhouette `< 0` are discarded
- only retained MPs proceed to adjusted-score construction, correlation, and LR annotation

### 7.4 Intersect each score matrix with the atlas cell set

For each compartment, the script intersects:

- cells present in the score matrix
- cells present in the atlas
- cells that satisfy the compartment-specific cell definition

This produces a single consistent cell set per compartment for all downstream steps.

### 7.5 Build per-sample adjusted scores

Adjusted scores are always percentages within a compartment and sample.

For MP-scored compartments:

1. A cell is considered positive when `UCell > 0.25`.
2. Cells are grouped by `orig.ident`.
3. For each MP in each sample, the adjusted score is:

`100 x (number of positive cells in sample) / (total cells of that compartment in sample)`

For finalized cancer states:

1. A one-hot matrix is constructed from the final labels.
2. State positivity is the assigned label itself, not a UCell threshold.
3. For each final state in each sample, the adjusted score is:

`100 x (number of malignant epithelial cells assigned to that state in sample) / (total malignant epithelial cells in sample)`

This distinction is important:

- cancer MP modes use UCell-thresholded cancer scores
- cancer state modes use assignment-based cancer-state fractions
- non-cancer compartments always remain MP-based and UCell-thresholded

### 7.6 Compute node coverage summaries

Each retained node stores:

- compartment
- display cell type
- node identifier
- node label
- underlying MP or finalized state
- sample coverage count
- denominator sample count
- coverage percentage
- whether coverage exceeds the minimum threshold

These summaries are saved as `Auto_celltype_node_summary.csv` within each analysis-mode folder.

### 7.7 Generate cutoff-sensitivity diagnostics

Cutoff sensitivity is meaningful only for UCell-scored compartments.

The script evaluates the cutoff grid:

- `0.10`
- `0.15`
- `0.20`
- `0.25`
- `0.30`
- `0.35`
- `0.40`
- `0.50`
- plus any user-provided cutoff

For each cutoff and each UCell-scored node, it records:

- positive sample coverage
- positive sample percentage
- positive cell fraction
- whether the node passes the sample-coverage filter

In the cancer-state modes, the cancer compartment is excluded from cutoff sensitivity because cancer states are assignment-based rather than threshold-based.

### 7.8 Define eligible sample overlap for each celltype pair

For each compartment pair:

1. The script identifies samples present in both compartments.
2. Those samples are mapped to studies.
3. Shared samples are counted per study.
4. Only studies with at least `10` shared samples are retained.
5. The corresponding eligible samples are kept for the actual correlation tests.

This rule is applied to:

- only cross-celltype pairs in the `cross_only` modes
- both cross-celltype and same-celltype pairs in the `cross_and_within` modes

### 7.9 Compute Pearson and Spearman correlations

For every eligible node pair, the script computes across the retained samples:

- Pearson correlation coefficient and P value
- Spearman correlation coefficient and P value
- `-log10(Spearman P)`

The full result table is written as `Auto_celltype_correlations_all.csv` inside each mode folder.

### 7.10 Define positive and negative edge sets

Positive edges are defined as:

- Pearson `r > 0`
- Spearman `P < 0.05`
- `-log10(Spearman P) >= 4`

Negative edges are defined as:

- Pearson `r < 0`
- Spearman `P < 0.05`

The stronger threshold for positive edges mirrors the published Fig. 5a logic.

### 7.11 Build bubble and network visualizations

Each mode writes:

- `Auto_celltype_correlation_bubble.pdf`
- `Auto_celltype_positive_network.pdf`
- `Auto_celltype_negative_network.pdf`
- `Auto_celltype_interaction_dotmap.pdf`
- `Auto_celltype_interaction_dotmap_data.csv`

Bubble plot design:

- one page per focal cell type
- focal cell type always on the x-axis
- all focal-to-partner comparisons shown on that focal page
- fixed panel size so sparse partner sets do not stretch the page
- cancer MP ordering follows the finalized state-linked order

Network plot design:

- node color follows the publication palette from `analysis/plotting/publication_umap.R`
- node size reflects positive sample coverage
- edge color and width reflect `-log10(Spearman P)`
- labels are laid out with force-directed coordinates and repel text

Interaction dotmap design:

- the PDF is paginated with one page per focal cell type
- each page places only the focal cell type's MPs/states on the rows
- columns contain all possible partner MPs/states, grouped into celltype facets
- same-celltype columns are included only in the cross-and-within modes
- cancer appears first, with cancer MPs labelled by description and ordered by the existing cancer MP tree/state-linked order
- non-cancer compartments follow `fibroblast`, `endothelial`, `cd8`, `cd4`, `macrophage`, `nk`, `plasma`
- only interactions meeting the final positive-network threshold are shown: Pearson `> 0`, Spearman `P < 0.05`, and `-log10(Spearman P) >= 4`
- dot fill represents Spearman rho
- dot size represents the percentage of eligible shared samples where both interacting nodes have adjusted score `> 0`
- interactions with retained LR support are marked by a black ring and cross, with no explanatory caption printed on the page
- the paired CSV stores the focal cell type, partner cell type, plotted row/column labels, support percentage, Spearman statistics, and LR-support flag

### 7.12 Load and filter the ligand-receptor catalog

The script reads the `All.Pairs` sheet from:

- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Ligand_Receptor_Pairs.xlsx`

The LR table is standardized to a common schema and then filtered so that only:

- `literature supported`
- `putative`

are retained.

All `EXCLUDED*` rows are removed before any LR matching.

### 7.13 Recover the positive cells contributing to each retained positive edge

LR annotation is applied only to retained positive edges.

For each positive edge:

1. The script recovers the eligible studies stored for that edge.
2. It derives the matching eligible samples.
3. It identifies the cells in those samples that contributed to the node definition.

Contribution logic depends on the node type:

- MP nodes: cells with `UCell > 0.25`
- finalized cancer-state nodes: cells carrying that final state label

This means the LR step uses exactly the cells that generated the sample-level adjusted score for that node.

### 7.14 Rank genes and test LR support

For each positive node:

- the RNA assay is used
- mean expression is calculated across the contributing cells
- genes are ranked by mean expression
- the top `4,000` genes are retained

For each positive edge, the script then checks both biologically relevant directions:

1. driver top genes as ligands against target node gene-set receptors
2. driver top genes as receptors against target node gene-set ligands

Cancer-state gene sets come from the ranked final-state marker table. MP gene sets come from the silhouette-filtered MP objects.

### 7.15 Summarise LR support

The raw LR table stores one row per retained LR match. The script then aggregates LR support to the edge level, including:

- number of matched LR rows
- number of distinct retained LR pairs
- number of literature-supported pairs
- number of putative pairs
- top example retained pairs

These are written as:

- `Auto_celltype_ligand_receptor_pairs.csv`
- `Auto_celltype_ligand_receptor_edge_summary.csv`
- `Auto_celltype_positive_edges_lr_annotated.csv`

### 7.16 Export the focal-celltype Excel workbook

Each mode writes:

- `Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx`

Workbook structure:

- `Overview` sheet
- one sheet per focal cell type

Within each focal sheet:

- all rows for that focal cell type are retained
- partner cell types are grouped together
- reverse-direction rows are not removed from the workbook
- `Direction` is written as `ligand` or `receptor`
- only the retained evidence class is kept from the LR catalog metadata
- database-support bookkeeping columns are omitted
- tab color matches the focal cell type

### 7.17 Write mode-level summaries and global summary

Each mode writes `Auto_celltype_mode_summary.csv`, including:

- analysis id
- cancer representation
- whether within-celltype tests were included
- cache version
- active UCell cutoff
- cancer positive-rule description
- number of nodes
- number of pairwise tests
- number of positive and negative edges
- cancer-specific diagnostic counts
- LR catalog retention counts
- workbook path

The script also combines all mode summaries into:

- `updates/new_updates/summaries/Auto_mp_cross_celltype_correlations_summary.csv`

The cancer-specific diagnostics are included to distinguish:

- lack of cancer-state coverage
- from lack of cancer-state edges passing the positive-network significance threshold

---

## 8. Interpreting the Cancer-State Modes

In the finalized cancer-state modes, there is no cancer-state positivity cutoff. If a cancer-state network shows few or no cancer-state edges in the positive graph, the correct interpretation is:

- first check that cancer-state coverage is high in `Auto_adjusted_scores_cancer.csv` and `Auto_celltype_node_summary.csv`
- then check `Auto_celltype_mode_summary.csv` for:
  - `n_cancer_positive_edges_p05`
  - `n_cancer_positive_edges_sig4`
  - `max_cancer_positive_spearman_significance`

If cancer-state nodes have strong sample coverage but `n_cancer_positive_edges_sig4 = 0`, then the absence of cancer-state edges in the positive network is due to the significance threshold, not to state sparsity or a missing positivity cutoff.

---

## 9. Cache Structure

Each analysis mode has its own cache directory:

- `01_cancer_mps_cross_only/cache/`
- `02_cancer_mps_cross_and_within/cache/`
- `03_cancer_states_cross_only/cache/`
- `04_cancer_states_cross_and_within/cache/`

The cached steps are:

- `Auto_celltype_step1_compartment_cache.rds`
- `Auto_celltype_step2_correlation_cache.rds`
- `Auto_celltype_step3_lr_cache.rds`

Each cache stores the current `cache_version`. If the version changes, the step is rebuilt automatically. A full rebuild can also be forced manually.

---

## 10. Run Commands

### Interactive full rebuild

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
AUTO_MPXCELL_FORCE_REBUILD=TRUE Rscript analysis/non_malignant_nmf/mp_cross_celltype_correlations.R
```

### Interactive run with a custom non-cancer UCell cutoff

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
AUTO_MPXCELL_FORCE_REBUILD=TRUE Rscript analysis/non_malignant_nmf/mp_cross_celltype_correlations.R 0.20
```

### PBS submission

```bash
/opt/pbs/bin/qsub -v AUTO_MPXCELL_FORCE_REBUILD=TRUE analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.sh
```

### PBS submission with a custom non-cancer UCell cutoff

```bash
/opt/pbs/bin/qsub -v AUTO_MPXCELL_FORCE_REBUILD=TRUE,cutoff=0.20 analysis/non_malignant_nmf/Auto_mp_cross_celltype_correlations.sh
```

---

## 11. Output Layout

All outputs are written under:

- `ref_outs/non_malignant_mp_correlations/`

Each of the four mode folders contains:

- `Auto_adjusted_scores_<compartment>.csv`
- `Auto_celltype_node_summary.csv`
- `Auto_celltype_cutoff_sensitivity.csv`
- `Auto_celltype_cutoff_sensitivity.pdf`
- `Auto_celltype_shared_sample_summary.csv`
- `Auto_celltype_correlations_all.csv`
- `Auto_celltype_correlations_positive.csv`
- `Auto_celltype_correlations_negative.csv`
- `Auto_celltype_correlation_bubble.pdf`
- `Auto_celltype_interaction_dotmap.pdf`
- `Auto_celltype_interaction_dotmap_data.csv`
- `Auto_celltype_positive_network.pdf`
- `Auto_celltype_negative_network.pdf`
- `Auto_celltype_ligand_receptor_status.csv`
- `Auto_celltype_ligand_receptor_pairs.csv`
- `Auto_celltype_ligand_receptor_edge_summary.csv`
- `Auto_celltype_positive_edges_lr_annotated.csv`
- `Auto_celltype_positive_network_lr_annotated.pdf`
- `Auto_celltype_ligand_receptor_summary.pdf`
- `Auto_cross_celltype_ligand_receptor_pairs_by_focal_celltype.xlsx`
- `Auto_celltype_mode_summary.csv`

---

## 12. Reproducibility Notes

- Correlations are computed at the sample level, not at the single-cell level.
- Cancer-state positivity is assignment-based, not threshold-based.
- Non-cancer compartments remain MP-based and use the active UCell cutoff.
- LR support annotates positive edges but does not change the correlation statistics themselves.
- If outputs look inconsistent with recent code changes, rebuild after updating `cache_version` or set `AUTO_MPXCELL_FORCE_REBUILD=TRUE`.

---

## 13. Whole-Celltype Direct-Annotation Abundance Modes (20 Jul 2026)

The workflow also writes two cancer-versus-whole-celltype modes:

1. `05_cancer_mps_vs_whole_celltypes`
2. `06_cancer_states_vs_whole_celltypes`

They reuse the cancer adjusted-score tables generated by modes `01` and `03`, respectively. Cancer is therefore still resolved into individual refined MPs or individual primary state labels. In contrast, each non-epithelial partner is represented once per sample by its final direct annotation, not by its NMF MPs.

The modes use `celltype_update` from the complete atlas directly. The retained non-epithelial annotations are fibroblast, endothelial, t.cell, macrophage, nk.cell, and plasma; `t.cell` remains one annotation-level compartment rather than creating separate CD4/CD8 whole-cell signatures.

The denominator is all confidently annotated atlas cells in a sample: cells with missing/empty `celltype_update` labels or labels starting with `unresolved` are excluded. No UCell scoring, gene-list activity, or cutoff is applied in these modes.

The per-sample score is `100 x cells with the target final celltype annotation / confidently annotated atlas cells`. Since `Expr_filtering.R` retains the confident annotations by marker-expression and singlet criteria, this is the direct final-annotation TME abundance requested for the aggregate modes.

Eligibility retains studies with at least 10 cancer-scored samples. Each cancer feature is correlated with each whole celltype abundance using Pearson and Spearman correlations across all eligible cancer-scored samples. Positive associations require Pearson `> 0` and Spearman `P < 0.05`; negative associations require Pearson `< 0` and Spearman `P < 0.05`. No additional score or `-log10(P)` cutoff is applied to these direct-abundance modes.

The aggregate modes use the same cross-compartment output contract and plotting helpers as modes `01`--`04`. Each folder therefore contains cancer and non-malignant adjusted-score exports, node coverage and shared-sample summaries, all/positive/negative correlation tables, a full bubble overview, positive/negative networks, the multi-page focal interaction dotmap, direct-abundance/cutoff documentation, Excel workbooks, LR-status files, and mode summary. The non-malignant inputs are one direct-abundance node each for fibroblast, endothelial, t.cell, macrophage, nk, and plasma; the previous CD4/CD8 split is intentionally replaced by the final `t.cell` annotation.

All eligible cross-compartment pairs are tested, not only cancer--TME pairs: this gives 123 tests for the 18 cancer-MP plus six whole-celltype mode, and 45 tests for the five cancer-state plus six whole-celltype mode. Direct abundance has no paired MP gene set, so LR outputs are retained with an explicit `not_applicable` status rather than fabricated ligand-receptor evidence.
