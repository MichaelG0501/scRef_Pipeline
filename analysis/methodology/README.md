# Analysis Methodology Index

Complex active workflows must reference a methodology file in the opening registry. Methodology files are organized to mirror `analysis/` subfolders and describe implemented logic, thresholds, decisions, and limitations rather than an idealized method. Simple deterministic utilities or plot-only scripts may state `Methodology: not required` in the header; a thin methodology file that merely repeats the header is not useful.

## Required Script Header Fields

Every new or substantially edited script must start with:

- script name
- status: `active`, `terminal`, `legacy`, or `delete-candidate`
- short description
- methodology path, or a concise `not required` rationale for a simple script
- inputs, including external downloads or absolute reference files
- outputs, grouped as `intermediate/`, `tables/`, `figures/`, `logs/`, or `reports/` when applicable
- cache/replot behavior for long-running scripts
- run command and conda environment

## Folder Methodology Files

- `cell_states/state_workflows_methodology.md`
- `clinical/clinical_bulk_and_association_methodology.md`
- `cnv/cnv_workflows_methodology.md`
- `developmental/developmental_reference_methodology.md`
- `enrichment/enrichment_external_reference_methodology.md`
- `metaprograms/metaprogram_scoring_methodology.md`
- `non_malignant_nmf/non_malignant_nmf_methodology.md`
- `OCCAMS/occams_bulk_rnaseq_reconstruction_methodology.md`
- `plotting/publication_plotting_methodology.md`
- `spatial/spatial_mapping_methodology.md`
- `summary/summary_methodology.md`

Script-specific methodology files can live beside these folder-level files when a workflow is complex enough to need a dedicated document.
