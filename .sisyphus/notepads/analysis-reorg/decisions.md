# Decisions — analysis-reorg

## [2026-03-01] Session ses_359751460ffeVdzvEIV53gBLcU — Initialization

### Architecture Decisions
1. **Subfolder structure**: metaprograms/, enrichment/, clinical/, cnv/, cell_states/, plotting/, non_malignant_nmf/, developmental/
2. **utils.R**: Extract shared patterns used in 3+ scripts only (no premature DRY)
3. **Shell parameterization**: non_mali_nmf/ per-cell-type shells → 2 parameterized runners
4. **PDO vs SC**: Keep in SEPARATE files — different data sources, naming, logic
5. **Merge strategy**: Read every file fully before any merge decisions
6. **Naming**: Descriptive English names without `Auto_` prefix for reorganized files

### Pending Decisions (require reading files first)
- How to handle ref_pipeline.R (1165L) — unique functions unclear until read
- Whether cibersort_result.R and evaluate_clinical_MPs.R are mergeable
- What residual.R (821L) actually does — MP-adjacent or standalone?
- Whether qc_unresolved_states.R is truly redundant with states_qccheck.R
