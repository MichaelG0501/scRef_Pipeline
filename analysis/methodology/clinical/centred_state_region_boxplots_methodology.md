# Centred State Region Boxplots

## Purpose

Compare the centred refined metaprogrammes defining Basal to intestinal metaplasia and SMG to intestinal metaplasia across distal oesophagus, GOJ, and stomach, then compare the abundance of the two assigned states across those regions.

## Inputs

- `ref_outs/meta_full_epi.rds`
- `ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_ucell_scores.rds`
- `ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds`
- `Concise_Summary_EAC_Ref.xlsx`, sheet 3

## Definitions

The final centred grouping is copied from `tcga_mp_survival_volcano_centred.R`.

- Basal to intestinal metaplasia: `MP14`, `MP3+`, `MP6+`, `MP11+`, `MP9+`, `MP10+`.
- SMG to intestinal metaplasia: `MP8+`, `MP8b`, `MP16`, `MP18b`, `MP17`.

`Tumor Location` is retained verbatim and mapped to a three-level region variable only when it identifies distal oesophagus, GOJ/GEJ, or stomach/gastric tissue. Unmapped values are written to a QC table and excluded from region comparisons.

## Analysis

For each requested state, the MP activity analysis retains cells assigned to that state, averages each constituent MP's UCell score per sample, and compares sample-level values among regions. State abundance is the percentage of all epithelial cells per sample assigned to each requested state. Global regional tests use Kruskal-Wallis tests with Benjamini-Hochberg correction within each analysis output; pairwise Wilcoxon tests are also exported with BH correction within state and analysis type.

## Outputs

Outputs are under `ref_outs/Metaprogrammes_Results/centred/state_region_boxplots/`:

- `tables/`: raw-to-standardised location mapping, unmapped samples, sample-level values, global statistics, and pairwise statistics.
- `figures/`: one multi-page MP activity PDF and one state abundance PDF.
- `logs/`: compact run summary.

The combined inferential summary is also written to `updates/new_updates/summaries/centred_state_region_boxplots_summary.csv`.
