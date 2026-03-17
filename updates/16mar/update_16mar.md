# Progress Update — 16 March 2026

**Project:** scRef_Pipeline & PDOs_Pipeline

---

## Context
This cycle focused on comprehensively profiling MP and state abundances across diverse samples, utilizing a pan-cancer approach to subclassify "Unresolved" cells, and rigorously defining and comparing PDO states, especially across patient-matched treated vs. untreated pairs.

## Summary of Work
* Generated sample abundance plots for top MPs and states across the entire scRef dataset, sorted strategically by geometric mean diversity and alphabetical cluster.
* Refined classification of previously "Unresolved" epithelial cells by integrating pan-cancer (3CA) metaprograms, selecting the most robust MPs based on sample/study coverage.
* Derived defined PDO metaprogram-based states mapping analogous "Approach B" conventions.
* Conducted paired comparisons (untreated vs treated) for 4 specific PDO patients tracing matched shifts in global proportions, UMAP identities, and raw MP expression metrics. 
* Mapped TCGA EAC transcriptomic cohorts to evaluate comparative clinical metrics between overarching cell states and constituent MPs.

## Key Outputs

| File | Associated Script | Description |
|------|-------------------|-------------|
| `Auto_sample_abundance.pdf` | `Auto_sample_abundance.R` | Sorted multi-pane per sample abundance composition referencing MPs and Approach B states. | 
| `Auto_unresolved_relabel_heatmap.pdf` / `*.proportion.pdf` / `*.volcano.pdf` | `Auto_unresolved_relabel.R` | Granular characterization incorporating 3CA Pan-cancer classifications over original unresolved cohorts with corresponding updated heatmap and clinical significance charts. |
| `Auto_PDO_states_heatmap_B_noreg.pdf` | `Auto_PDO_states_analysis.R` | High-fidelity derived definitions and spatial mappings mapping unique PDO state characteristics. | 
| `Auto_PDO_survival_state_volcano_EAC_noreg.pdf` | `Auto_PDO_states_analysis.R` | Contrast establishing the predictive difference in TCGA-EAC cohort survivals between broad states vs discrete MPs. |
| `Auto_PDO_matched_state_proportions.pdf` / `*.UMAP_combined.pdf` / `*.MP_expression.pdf` | `Auto_PDO_matched_analysis.R` | Comprehensive paired evaluation evaluating untreated to matched treated pairs isolating spatial clusters and proportion drifts. | 

## Decisions Made
1. **Geometric Sort Sorting:** Embraced sorting composition plots iteratively via sample diversity logic to rapidly index highly heterogeneous conditions instead of relying on default alphabetically randomized approaches. 
2. **Expanding Unresolved Taxonomy:** Prioritized cross-cancer coverage limits over strict silhouette thresholds when pulling 3CA-based annotations to expand unresolved clusters. 
3. **Paired Isolate Examination:** Specifically subset the global PDO matrix to four complete sample pairings (e.g. `SUR1070_Treated_PDO` & `SUR1070_Untreated_PDO`) leveraging connected line graphs rather than dense generalized heatmaps for intuitive tracking. 

## Next Steps
- [ ] Incorporate pipeline iterations resolving deeper hybrid intersections spanning into subsequent cell groupings.
- [ ] Evaluate survival correlation models combining paired post-treatment shifts.
