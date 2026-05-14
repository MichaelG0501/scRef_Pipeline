# Analysis Methodology Index

Each active script must reference a methodology file in its opening header. Methodology files are organized to mirror `analysis/` subfolders and should describe implemented logic, not an idealized method.

## Required Script Header Fields

Every new or substantially edited script must start with:

- script name
- status: `active`, `terminal`, `legacy`, or `delete-candidate`
- short description
- methodology path
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
- `plotting/publication_plotting_methodology.md`
- `spatial/spatial_mapping_methodology.md`
- `summary/summary_methodology.md`

Script-specific methodology files can live beside these folder-level files when a workflow is complex enough to need a dedicated document.
