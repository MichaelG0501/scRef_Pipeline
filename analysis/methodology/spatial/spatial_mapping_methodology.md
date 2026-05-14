# Spatial Mapping Methodology

This document covers `analysis/spatial/`.

`export_scatlas_visium_signatures.R` exports scATLAS state/MP signatures for spatial mapping.

`map_scatlas_states_visium.py` maps scATLAS states to Visium data.

`map_scatlas_states_xenium.py` is currently untracked and should remain unstaged unless the user explicitly requests adding it.

Spatial scripts must document:

- required external spatial object paths
- signature files read from `ref_outs/`
- output locations
- whether they overwrite previous mapping results
- whether they can replot from cached mapped scores
