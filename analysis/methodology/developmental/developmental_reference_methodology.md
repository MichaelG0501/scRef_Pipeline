# Developmental Reference Methodology

`developmental.R` builds the developmental reference objects used by enrichment and annotation scripts.

Inputs are external xlsx files under the spatial transcriptomics project. The output reference path is:

- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/enrich_dev.rds`

Per-stage RDS files are stored under:

- `analysis/developmental/per_stage/`

The script should be treated as a reference-builder rather than a frequent plotting script. If developmental references are updated, downstream enrichment scripts should be rerun and the change should be recorded in `analysis/ANALYSIS_MAP.md`.
