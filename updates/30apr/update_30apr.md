# Weekly Progress Update — 30 April 2026

## Summary

This update covers five key areas of ongoing development:

1. **Parse Data Analysis** (18 slides)
   - Metaprogram optimization across nMP range
   - Enrichment annotation across multiple databases (3CA, Hallmarks, GO, developmental stages, adult tissues)
   - Correlation and Jaccard similarity analysis comparing Parse, PDO, and scATLAS datasets
   - MP activity boxplots stratified by dataset
   - Top metaprogram abundance distribution across samples

2. **CNA Subclones Analysis** (3 slides)
   - Methodology for identifying malignant subclones from CNA profiles
   - Per-sample metaprogram activity within subclone populations
   - Subclone state assignment via top-MP classification
   - Representative sample pages showing integration of CNA and metaprogram signals

3. **Pan-cancer MPs Expression** (1 slide)
   - 3CA pseudobulk correlation cross-data comparison

4. **Ongoing Tasks** (1 slide)
   - RNA velocity preprocessing for parse dataset (loom files generated)
   - Cell trajectory analysis in progress
   - CNV inference parse workflow (T0 vs PDO comparison)
   - Chemical inhibitor screening (ASGARD completed, manual mapping underway)
   - Spatial dataset download and preprocessing completed; implementing scATLAS state mapping
   - Signature heatmaps for basal-to-intestinal vs SMG metaplasia in development
   - PDO RNA velocity analysis planned

## Key Assets

### Parse Analysis
- Optimal nMP selection via multiple criteria
- Metaprograms heatmap and enrichment landscapes
- Cross-dataset validation (PDO, scATLAS, Parse)
- Per-sample composition and state distribution

### CNA Subclones
- Integrated CNA-metaprogram profiling
- Sample-level subclone characterization
- State assignment consistency checks

### Cross-data Integration
- Pan-cancer MP expression patterns
- 3CA reference comparison

## Next Steps

1. **Continue parse pipeline development**
   - Complete RNA velocity and trajectory analyses
   - Finalize CNV inference parsing

2. **Expand spatial mapping**
   - Map full scATLAS states to spatial datasets
   - Visualize spatial distribution of metaprogrammed states

3. **Chemical inhibitor integration**
   - Compare ASGARD and manual mapping results
   - Prepare consensus predictions

4. **PDO analysis acceleration**
   - Complete RNA velocity workflows
   - Integrate with treated/untreated state comparisons

---

**Generated:** 30 April 2026  
**Total Slides:** 26  
**Update Format:** PDF (Beamer)
