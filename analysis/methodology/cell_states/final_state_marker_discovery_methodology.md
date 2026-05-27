# Auto Six-State Marker Methodology

Generated: 2026-05-26 15:42:41 BST

## 1. Goal and Scope

The workflow aims to identify the top 5 most robust and specific markers for each finalized state. It prioritizes genes that are not only highly expressed but also consistently reproducible across the heterogeneous cohort of samples and studies, while maintaining clear specificity for a single state.

**Core Inputs:**
- `ref_outs/EAC_Ref_epi.rds`: The main epithelial Seurat object (75,348 cells).
- `ref_outs/Auto_final_states.rds`: Finalized six-state labels, with `Unresolved` and `Hybrid` excluded.

**Finalized States Retained:**
- `Classic Proliferative`
- `Basal to Intestinal Metaplasia`
- `Stress-adaptive`
- `SMG-like Metaplasia`
- `Immune Infiltrating`
- `3CA_EMT_and_Protein_maturation`

---

## 2. Six-State Subset and Re-embedding

To ensure marker analysis is focused on the finalized transcriptional landscape, a clean subset and re-embedding are performed.

### 2.1 Cell and Feature Selection
1. The global epithelial count matrix is subset to cells present in both inputs.
2. Cells with `Unresolved` or `Hybrid` labels are removed.
3. Genes detected in fewer than **20 cells** within this six-state subset are discarded.

### 2.2 Re-embedding Pipeline
The subsetted object is processed through a standard Seurat pipeline:
- **Normalization:** `NormalizeData` (standard log-normalization).
- **HVGs:** `FindVariableFeatures` (vst method, **3000 features**).
- **Scaling:** `ScaleData` on the HVGs.
- **PCA:** `RunPCA` (**30 PCs**).
- **Neighbors:** `FindNeighbors` on 30 PCs.
- **Clustering:** `FindClusters` (Leiden algorithm, **resolution 0.5**).
- **UMAP:** `RunUMAP` on 30 PCs.

---

## 3. Candidate Gene Screening (Pooled)

Because sample-wise DGE is computationally expensive, a global descriptive screen is used to define a tractable candidate set for detailed recurrence testing.

For each state:
1. Total cells are divided into **State** vs. **Rest**.
2. Mean expression and detection frequency (pct) are computed for both groups.
3. Genes are ranked by `global_mean_difference` (State - Rest).
4. The **top 1500** genes with positive mean difference are retained as candidates for that state.

---

## 4. Sample-Wise Recurrent DGE

The core of the reproducibility analysis is running DGE within individual samples to see how often a gene is a "hit."

### 4.1 Sample Eligibility
A sample is considered "eligible" to test a specific state only if:
- The sample contains **at least 20 cells** belonging to the target state.
- The sample contains **at least 20 cells** belonging to the other five states combined (the "Rest").

### 4.2 Differential Expression Testing
For each eligible sample and each of its qualified states:
- **Test:** Seurat `FindMarkers` (Wilcoxon rank-sum test).
- **Universe:** One state versus the other 5 states within that sample.
- **Genes:** Only the 1500 candidate genes identified in the pooled screen for that state.
- **Thresholds:** No hard expression or logFC gates are applied at the testing stage (`logfc.threshold = 0`, `min.pct = 0`) to capture all valid statistics.

### 4.3 Hit Definition
A gene is defined as a "hit" within a specific sample if:
- `p_val_adj < 0.05` (FDR corrected).
- `avg_log2FC > 0` (statistically higher in the target state).

---

## 5. Sample-Aware State Specificity

To ensure markers are globally specific across the atlas, a "specificity gap" is computed using sample-level medians.

1. For every gene in the marker summary, the mean expression is calculated within the target state for every sample eligible for that state.
2. The **median of these sample-level means** is computed, representing the "typical" expression of that gene in that state.
3. This is repeated for all six states.
4. **Specificity Gap:** The typical expression in the target state minus the maximum typical expression seen in any of the other five states.
5. A gene is considered a "best state match" only if the target state has the highest median expression.

---

## 6. Ranking and Final Selection

Genes are ranked within each state using a multi-component reproducibility and specificity score.

### 6.1 Ranking Metrics
Three metrics are computed per gene/state:
1. **Reproducibility Score:** 
   - `sample_recurrence`: Fraction of eligible samples that are DGE hits.
   - `study_recurrence`: Fraction of eligible studies that have at least one DGE hit sample.
   - `score = 0.5 * sample_recurrence + 0.5 * study_recurrence`.
2. **Effect Size:** Median `avg_log2FC` across all samples that were DGE hits.
3. **Specificity Gap:** The gap computed in Section 5.

### 6.2 The Ranking Score
Within each state, genes are assigned a `percent_rank` (0 to 1) for each of the three metrics above. The final **Ranking Score** is the sum of these three ranks.

### 6.3 Hard Selection Rules
To be considered for the final top 5, a gene MUST:
- Have at least one significant positive DGE hit in at least one sample (`hit_sample_n > 0`).
- Have the target state as its highest-expressing state globally (`best_state_match == TRUE`).
- Have a positive specificity gap (`specificity_gap > 0`).

The top 5 markers by **Ranking Score** per state are selected for the final panel.

---

## 7. Recurrence and Support Classification

- **Support Class:**
  - `multi-study`: Hit detected in samples from 2 or more studies.
  - `multi-sample_single-study`: Hit detected in 2 or more samples, but all within one study.
  - `single-sample_single-study`: Hit detected in only 1 sample.
- **Legacy Strict Check:** The markers are also checked against a "Legacy Strict" rule (`hit_sample_n >= 20% of eligible` and `hit_study_n >= 35% of eligible`). This highlights which public markers are extremely robust vs. those that are specific but lower sensitivity.

---

## 8. Heatmap Construction

- **Data:** Median of sample-level means per state (as computed in Section 5).
- **Z-scoring:** Values are Z-scored per row across the six states to highlight state-specific enrichment.

## 9. Current 3CA recurrence profile

```text
# A tibble: 5 × 7
  gene      hit_sample_n hit_study_n sample_recurrence study_recurrence
  <chr>            <int>       <int>             <dbl>            <dbl>
1 CD44                 2           1            0.0833            0.143
2 EXT1                 3           1            0.125             0.143
3 MYOF                 2           1            0.0833            0.143
4 MBNL2                2           1            0.0833            0.143
5 LINC03009            3           1            0.125             0.143
# ℹ 2 more variables: support_class <chr>,
#   passes_legacy_strict_recurrence <lgl>
```

## 10. State-level summary

```text
                          state cell_n eligible_sample_n eligible_study_n
          Classic Proliferative  12573                68                8
 Basal to Intestinal Metaplasia  16642                78                7
                Stress-adaptive  10169                59                7
            SMG-like Metaplasia   8173                53                7
            Immune Infiltrating   6056                51                7
 3CA_EMT_and_Protein_maturation   2903                24                7
 top_marker_n
            5
            5
            5
            5
            5
            5
```

## 11. Recurrence summary by state

```text
# A tibble: 6 × 7
  state                  top_marker_n multi_sample_marker_n multi_study_marker_n
  <chr>                         <int>                 <int>                <int>
1 3CA_EMT_and_Protein_m…            5                     5                    0
2 Basal to Intestinal M…            5                     5                    5
3 Classic Proliferative             5                     5                    5
4 Immune Infiltrating               5                     5                    5
5 SMG-like Metaplasia               5                     5                    5
6 Stress-adaptive                   5                     5                    5
# ℹ 3 more variables: median_sample_recurrence <dbl>,
#   median_study_recurrence <dbl>, n_passing_legacy_strict_recurrence <int>
```

## 12. Output Files

- `Auto_six_state_markers_final.csv`: The final top 5 markers per state with their ranking and recurrence stats.
- `Auto_six_state_markers_ranked.csv`: The full table of candidate genes ranked by the workflow.
- `Auto_six_state_markers_top5_recurrence_summary.csv`: Summary of hit counts and support classes.
- `Auto_six_state_markers_top5_sample_support.csv.gz`: per-sample support table for the final top-5 markers.
- `Auto_six_state_markers_top5_study_support.csv`: per-study support table for the final top-5 markers.
- `Auto_six_state_markers_top5_state_recurrence_summary.csv`: state-level sensitivity summary for the final top-5 markers.
- `Auto_six_state_marker_heatmap.pdf`: final publication-facing heatmap.
- `Auto_six_state_umap.pdf`: UMAP visualizations of the six-state subset.

