# Progress Update — 28 May 2026

**Project:** scRef_Pipeline / scATLAS

---

## Context
This cycle focused on advancing several key analyses across the pipeline: validating recurrent CNA events using TCGA data, performing an in-depth investigation of MP18 CNA associations across subclones, quantifying state dispersion and co-localisation in UMAP space, mapping scATLAS states to Xenium spatial transcriptomics data, and unifying the developmental MP enrichment workflow.

## Summary of Work
- **MP18 CNA Deep Dive:** Investigated the significant association between MP18 (Secretory Diff. Intest.) and pairwise CNA distances across subclones to identify specific chromosome arm drivers and sample contributions.
- **TCGA Validation:** Executed recurrent event validation for CNA profiles using external TCGA data.
- **Cell State Dispersion:** Quantified how epithelial states occupy per-sample UMAP space by measuring distance to centroids (dispersion) and local mixing (co-localisation).
- **Spatial Mapping:** Developed python scripts to map scATLAS states onto Xenium and Visium spatial data.
- **Developmental Validation:** Created a unified script for developmental reference validation of scATLAS metaprogrammes, generating both "all genes" and "top 50" comparative enrichment heatmaps.

## Key Outputs
| File | Description |
| --- | --- |
| `analysis/cnv/Auto_mp18_cna_investigation.R` | Investigates MP18 association with pairwise CNA distances across subclones |
| `analysis/TCGA/tcga_cna_recurrent_event_validation.R` | TCGA recurrent event validation for CNA profiles |
| `analysis/cell_states/state_umap_dispersion_colocalisation.R` | Quantifies UMAP dispersion and state co-localisation |
| `analysis/cell_states/basal_mp_distance_matrix.R` | Generates Basal-to-IM metaprogram distance matrix nodeplots |
| `analysis/spatial/map_scatlas_states_xenium.py` | Maps scATLAS states to Xenium spatial transcriptomics data |
| `analysis/developmental/developmental_mp_enrichment_unified.R` | Unified developmental reference validation for scATLAS MPs |

## Decisions Made
1. **Unified Enrichment:** Adopted a unified approach for developmental MP enrichment to systematically compare "top 50" vs "all" marker genes against reference databases.
2. **Focus on MP18:** Prioritized MP18 for CNA deep dive due to its highly significant sample-centered correlation with absolute mean CNA.
3. **Spatial Metrics:** Leveraged existing scATLAS definitions to systematically map to high-resolution spatial contexts (Xenium/Visium).

## Next Steps
- [ ] Present MP18 CNA driver findings to the team.
- [ ] Review formatting of the unified developmental enrichment heatmaps.
- [ ] Refine visual parameters for the spatial state mappings.
- [ ] Finalize basal nodeplots and dispersion summary figures for the main manuscript.
