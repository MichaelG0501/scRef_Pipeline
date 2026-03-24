# Update Plan — 23 March 2026

## Objective
Create the progress update slides (PDF) and markdown summary for the supervisor meeting. Focus on finalized cell state definitions, TCGA survival comparison across methods, and batch-corrected pseudotime trajectories.

## Results to Include
1.  **scRef Heatmap**: Relabeled unresolved cells using 3CA MPs (8-state model).
2.  **State Proportions**: Overall and per-study stacking for scRef and PDO.
3.  **TCGA Survival**: Continuous vs Median vs Quantile splits comparisons.
4.  **DGE Overlap**: Side-by-side hallmark/pathway enrichment.
5.  **Cohort Comparison**: EAC/ESCC state activity.
6.  **Pseudotime**: Harmony and scVI batch integration partA trajectories.
7.  **3CA MP Correlation**: Cross-dataset (scRef vs PDO) bar and scatter plots.

## Steps
1.  Establish directory `updates/23Mar/` with `figures/` and `summaries/`.
2.  Discover latest results in `ref_outs/` and `PDOs_outs/` from this cycle's scripts (`states_unresolved_relabel.R`, `pseudotime_batch_correction_states.R`).
3.  Copy and rename figures to short, descriptive names suitable for LaTeX.
4.  Write `update_23Mar.tex` using the results-focused Beamer template.
5.  Compile LaTeX to PDF using `lualatex`.
6.  Write detailed `update_23Mar.md` progress document summarizing findings and decisions.
7.  Verify summary CSVs are present in `summaries/`.

## Note
User requested *not* to convert to PPTX until content is confirmed.
