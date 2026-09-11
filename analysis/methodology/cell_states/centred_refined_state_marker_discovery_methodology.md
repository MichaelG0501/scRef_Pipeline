# Centred-refined state marker discovery methodology

This workflow discovers markers for the five current centred-refined noreg OAC epithelial state archetypes. It supersedes use of marker tables generated from `Auto_final_states.rds` for new manuscript work.

Cells assigned `Hybrid` or `Unresolved` are excluded from archetype-versus-archetype marker testing. They remain represented in state-abundance and uncertainty panels elsewhere.

For each target state, cells are aggregated into two pseudobulks per eligible sample: target-state cells and all other archetypal-state cells. A sample is eligible when both groups contain at least 20 cells. The contrast is fitted with edgeR quasi-likelihood models using sample as a blocking factor (`~ sample + target/rest`). The biological unit is therefore the sample, not the cell.

The complete result reports log fold change, FDR, the number of paired samples, and target-sample expression coverage. Main marker candidates require positive log fold change, FDR below 0.05, and CPM at least 1 in at least half of target pseudobulks. When a gene is eligible for more than one state, it is assigned to the state with the largest positive log fold change before selecting the top five display markers per state.

The terminal dot plot uses mean log-normalized single-cell expression only as a descriptive display layer. Dot size is the percentage of cells with non-zero normalized expression; colour is gene-wise z-scored mean expression. Inferential statistics come exclusively from the sample-blocked pseudobulk model.

All complete marker results, selected markers, display source data, contrast sample counts, and the replot cache are stored in the live project under `ref_outs/Metaprogrammes_Results/centred/state_markers/`. The full marker table is a downstream input for current state-level ligand-receptor support and manuscript figures.

