# Progress Update — 18 May 2026

**Project:** scRef_Pipeline / PDOs_Pipeline

---

## Context
This cycle focused on deeper genomic mapping and CNA correlations, notably mapping MP genes to chromosomes and linking recurrent CNA events (like chr8q myc gain) with MP expression. We also finalized the chemical inhibitor screening comparisons between scRef and PDOs.

## Summary of Work
* Visualized Tumor Cell Stiffness Signatures and T2T4-high MP trends.
* Performed chromosomal mapping of MPs genes for both scRef and PDOs.
* Analyzed dominant subclones and their correlation with MP expression and pairwise CNA distances.
* Identified and mapped recurrent CNA events (e.g., chr8q gain) and associated them with specific MP expression and QC metrics.
* Validated TCGA gender association with MP expression (scRef concordance).
* Correlated MP expression with Enrichment Reference Markers (Adult Epithelium and Barrett's Dataset).
* Finalized Chemical Inhibitors Screening (L1000 signature reversal profiles) for both scRef and PDOs.

## Key Outputs
| File/Output | Description |
|---|---|
| `Auto_parse_tumour_stiffness_module_gene_heatmap_no_pdo_sur1090.pdf` | Tumor cell stiffness signatures heatmap. |
| `Auto_mp_chromosomal_mapping.pdf` | Chromosomal mapping of MP genes (scRef & PDOs). |
| `Auto_v2_largest_subclone_effects_all_features.pdf` | Dominant subclones and MP expression effects. |
| `Auto_v2_recurrent_cna_event_associations_all_features.pdf` | Recurrent CNA events and MP expression association. |
| `Auto_drug_reversal_l1000_signature_reversal_profiles.pdf` | Finalized chemical inhibitors screening. |

## Decisions Made
1. Chromosomal location impacts certain MP expressions, confirming the need to associate CNA patterns with states.
2. Validated the TCGA gender association to rule out bias in certain state classifications.
3. Finalized the chemical inhibitor candidates by comparing L1000 reversal profiles across both scRef and PDO pipelines.

## Next Steps
- [ ] Single cell conference poster
- [ ] Parse T2/T4 vs T0 DGE analysis and visualization
- [ ] scATLAS malignant cell states mapping in Xenium
- [ ] scATLAS chemical inhibitors screening
- [ ] Survival (meta) analysis in OCCAMS
- [ ] Clinical association validation/analysis in OCCAMS
- [ ] PDOs RNA velocity analysis
- [ ] PDOs Pseudotime analysis
- [ ] PDOs Numbat subclone analysis / validation
