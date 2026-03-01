# Issues — analysis-reorg

## [2026-03-01] Session ses_359751460ffeVdzvEIV53gBLcU — Initialization

### Known Risks
1. **setwd() collisions**: 4+ scripts setwd() to external repos. Must wrap in setwd()/restore or use full paths when merging.
2. **Shell↔R cross-references**: Auto_MP_correlation.sh hardcodes `Rscript analysis/Auto_MP_correlation_v2.R` — must update when moving
3. **AGENTS.md paths**: 15+ hardcoded paths will break on move — AGENTS.md update is atomic with reorganization (Task 14)
4. **enrich_plot.R**: Reads from external path `EAC_Ref_all/00_merged/developmental/` — must preserve this path in merged version

### Unresolved Flags
(None yet — will be populated during Task 1 triage)
