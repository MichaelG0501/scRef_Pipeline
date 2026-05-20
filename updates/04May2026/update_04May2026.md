# Progress Update — 04 May 2026

**Project:** scRef_Pipeline / PDOs_Pipeline

---

## Context
This cycle focused on copy number alteration (CNA) subclone analysis, exploring the differences between subclone states and QC metrics. We also evaluated Basal to Intestinal vs SMG-like metaplasia signatures, ran simulations for PDO marker selection, and finalized the candidates for the chemical inhibitors screening (drug reversal profiles).

## Summary of Work
* Completed the CNA subclone analysis and mapped them against metaprograms and states.
* Compared QC metrics between different CNA subclones to rule out technical confounders.
* Analyzed signature heatmaps for Basal to Intest vs SMG-like Metaplasia.
* Validated PDO marker expression specificity in SUR1090 and SUR1072.
* Conducted chemical inhibitors screening, selecting final candidates using rank-rank scatter and L1000 signature reversal profiles.
* Generated predicted state DEGs flipping heatmaps.

## Key Outputs
| File/Output | Description |
|---|---|
| `Auto_malignant_subclone_mp_cohort_summary.pdf` | CNA cohort summary and MP/States stats. |
| `Auto_basal_smg_mp_signature_heatmap.pdf` | Basal vs SMG metaplasia heatmap. |
| `Auto_drug_reversal_overlap_venns.pdf` | Final candidates overlap for chemical inhibitors. |
| `Auto_drug_reversal_l1000_signature_reversal_profiles.pdf` | L1000 signature reversal profiles. |

## Decisions Made
1. Prioritized specific chemical inhibitor candidates based on overlapping Venns and rank-rank scatter from drug reversal profiles.
2. Proceeded with L1000 signature reversal mappings for predicted state flipping.

## Next Steps
- [ ] Parse geneNMF high nMP setting and filter by trend/pattern
- [ ] RNA velocity parse (preprocessed raw bam files)
- [ ] Cell trajectory parse
- [ ] CNV inference parse - Compare T0 vs PDO
- [ ] States concordance analysis for different MP gene list selection
- [ ] Spatial datasets (implement code to map scATLAS states)
- [ ] PDOs RNA velocity analysis
- [ ] PDOs FLOT matched pairs geneNMF analysis
- [ ] PDOs subclone analysis
