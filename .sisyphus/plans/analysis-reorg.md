# Reorganize analysis/ Folder — Clean Modular Structure

## TL;DR

> **Quick Summary**: Reorganize 38 R scripts + 17 shell scripts in `analysis/` into themed subfolders, merging overlapping scripts, extracting shared utilities, and parameterizing duplicate shell wrappers. Every file is fully read before any decision.
> 
> **Deliverables**:
> - Themed subfolder structure under `analysis/` (metaprograms/, enrichment/, clinical/, cnv/, cell_states/, plotting/, non_malignant_nmf/, developmental/)
> - Merged scripts where overlap is heavy (MP analysis cluster, enrichment cluster, QC pair)
> - Shared `analysis/utils.R` with common patterns (silhouette filtering, heatmap helpers, UCell wrappers)
> - Parameterized shell scripts replacing 14 near-identical per-cell-type wrappers
> - Updated AGENTS.md documenting the new structure
> - Zero information/functionality loss — every unique analysis preserved
> 
> **Estimated Effort**: Large (38 files to read, assess, reorganize)
> **Parallel Execution**: YES — 6 waves (after sequential Wave 1 triage)
> **Critical Path**: Wave 1 (triage) → Wave 2-5 (parallel theme batches) → Wave 6 (integration + AGENTS.md)

---

## Context

### Original Request
Reorganize the analysis/ folder: read all unorganized scripts, understand what they do, then reorganize, merge redundant code, create clean modular structure with descriptive filenames, and update AGENTS.md.

### Interview Summary
**Key Discussions**:
- No blind deletions — even "temp" and "quick" named files may contain important unique code
- Legacy scripts (heatmap.R, ref_pipeline.R) contain unique functions not in main pipeline — must preserve
- PDO vs SC data sources/naming are totally different — must keep structurally distinct in merged scripts
- Shell scripts in non_mali_nmf/ are nearly identical — parameterize into one
- Subfolders by theme with descriptive English names
- Extract shared utilities into analysis/utils.R
- Full read of every file, batched by theme to manage context

**Research Findings**:
- 11/38 files use UCell scoring, 8 read EAC_Ref_epi.rds, 5 read UCell_default.rds
- Zero scripts source() each other — all currently self-contained
- states_qccheck.R and qc_unresolved_states.R produce identical output filenames
- 25/35 top-level scripts produce heatmaps, 7 do silhouette filtering
- Per-cell-type shell scripts differ only in PBS resource allocation
- MP analysis cluster (5 scripts, ~3.9K lines) has highest overlap — all do UCell + correlation + heatmap + silhouette

### Metis Review
**Identified Gaps** (addressed):
1. **AGENTS.md hardcoded paths**: 15+ paths like `analysis/example_anno.R` will break on move → AGENTS.md update is atomic with reorganization (Wave 6)
2. **Shell↔R cross-references**: `Auto_MP_correlation.sh` hardcodes `Rscript analysis/Auto_MP_correlation_v2.R`. All `.sh`+`.R` pairs must move atomically → enforced as guardrail
3. **Cross-project setwd() collisions**: 4+ scripts setwd() to external repos. During merge, must wrap in setwd()/restore pattern or use full paths → each merge task includes setwd() audit
4. **No behaviour verification post-merge**: No tests exist. → Smoke-run acceptance criterion (`Rscript --vanilla` exits 0, library loads succeed) for every merged script
5. **non_mali_nmf/ fate**: Existing subfolder with live PBS workflows → rename to `non_malignant_nmf/` to match theme naming, parameterize shells, update masters

---

## Work Objectives

### Core Objective
Transform the flat, ad-hoc `analysis/` folder into a clean, themed subfolder structure where every script has a descriptive name, redundant code is merged, shared patterns are extracted to a utility file, and AGENTS.md documents the complete final structure.

### Concrete Deliverables
- `analysis/utils.R` — shared utility functions
- `analysis/metaprograms/` — all MP analysis scripts (sc, pdo, correlation, scoring)
- `analysis/enrichment/` — gene overlap, enrichment, annotation scripts
- `analysis/clinical/` — TCGA, survival, CIBERSORTx, clinical variable scripts
- `analysis/cnv/` — CNV profiling, subsetting, plotting
- `analysis/cell_states/` — state assignment, QC, UMAP, annotation
- `analysis/plotting/` — general-purpose plotting scripts
- `analysis/non_malignant_nmf/` — cell-type NMF (parameterized shells)
- `analysis/developmental/` — developmental programme analysis
- Updated `AGENTS.md` with complete new structure

### Definition of Done
- [ ] Every original `.R` file accounted for (moved, merged, or deleted with justification)
- [ ] All merged scripts pass `Rscript --vanilla <script>.R` without error (syntax + library load check)
- [ ] All `.sh` scripts reference correct paths to their `.R` counterparts
- [ ] `analysis/utils.R` exists and is sourced by scripts that use shared patterns
- [ ] AGENTS.md "Downstream Analysis" section reflects new paths and structure
- [ ] Zero functionality loss — every unique analysis from original scripts is preserved

### Must Have
- Full read of every file before any reorganization decision
- PDO and SC logic clearly separated in merged metaprogram scripts
- Every `.sh` + `.R` pair moved atomically (never one without the other)
- setwd() collisions resolved in merged files (wrap or use full paths)
- Descriptive filenames that communicate purpose
- Each merged script has a clear header comment explaining what was merged and why

### Must NOT Have (Guardrails)
- **No blind deletions** — every file read before any delete decision
- **No functionality loss** — if a script does something unique, that code survives
- **No orphaned shell scripts** — every .sh must point to a valid .R after reorganization
- **No broken setwd() paths** — merged scripts must handle working directory correctly
- **No root pipeline changes** — steps 1-8 in project root are untouched
- **No ref_outs/ data changes** — output data files are never modified
- **No guessing** — if a file's purpose is unclear after reading, ASK the user before deciding
- **No over-abstraction** — utils.R contains only patterns used in 3+ scripts, not premature DRY
- **No Auto_ prefix on reorganized files** — the Auto_ prefix rule is for NEW scripts created by agents; reorganized existing scripts keep descriptive names without prefix

---

## Verification Strategy

> **ZERO HUMAN INTERVENTION** — ALL verification is agent-executed. No exceptions.

### Test Decision
- **Infrastructure exists**: NO (R analysis scripts, no test suite)
- **Automated tests**: None — no test framework
- **Framework**: N/A

### QA Policy
Every task includes agent-executed QA scenarios. Evidence saved to `.sisyphus/evidence/task-{N}-{scenario-slug}.{ext}`.

- **Syntax validation**: `Rscript --vanilla -e "source('<script>.R')"` or `Rscript -e "parse(file='<script>.R')"` — exits 0
- **File accounting**: Agent produces manifest mapping every original file → destination
- **Path integrity**: Grep all `.sh` files for `Rscript` calls, verify target `.R` exists
- **Structure verification**: `find analysis/ -type f | sort` matches expected layout

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Sequential — Triage + Foundation):
├── Task 1: Read + classify all "unclear" scripts [deep]
├── Task 2: Design final subfolder structure + file manifest [deep]
└── Task 3: Create analysis/utils.R skeleton [quick]

Wave 2 (Parallel — Theme Batches, after Wave 1):
├── Task 4: Metaprograms — read 6 scripts, merge + restructure [deep]
├── Task 5: Enrichment — read 5 scripts, merge + restructure [deep]
├── Task 6: Clinical/Survival — read 3 scripts, clean + restructure [deep]
├── Task 7: CNV + Cell States — read 6 scripts, merge QC pair + restructure [deep]
└── Task 8: Plotting + standalone — read remaining, rename + move [unspecified-high]

Wave 3 (Parallel — after Wave 2):
├── Task 9: Non-malignant NMF — parameterize shells, rename folder [unspecified-high]
├── Task 10: Developmental — assess, restructure if needed [quick]
└── Task 11: Populate utils.R with extracted shared functions [deep]

Wave 4 (Sequential — Integration):
├── Task 12: Wire utils.R into all scripts that need it [unspecified-high]
└── Task 13: Audit all shell script paths + fix cross-references [unspecified-high]

Wave 5 (Sequential — Documentation):
└── Task 14: Update AGENTS.md with complete new structure [writing]

Wave FINAL (Parallel — Verification):
├── Task F1: Plan compliance audit [oracle]
├── Task F2: Syntax validation of all reorganized scripts [unspecified-high]
├── Task F3: File accounting — every original file mapped [unspecified-high]
└── Task F4: Path integrity check — no orphaned references [deep]

Critical Path: Task 1 → Task 2 → Task 4 → Task 11 → Task 12 → Task 13 → Task 14 → F1-F4
Parallel Speedup: ~50% faster than sequential
Max Concurrent: 5 (Wave 2)
```

### Dependency Matrix

| Task | Depends On | Blocks |
|------|-----------|--------|
| 1    | —         | 2, 4-8 |
| 2    | 1         | 3, 4-8 |
| 3    | 2         | 11, 12 |
| 4    | 1, 2      | 11, 12, 13 |
| 5    | 1, 2      | 11, 12, 13 |
| 6    | 1, 2      | 11, 12, 13 |
| 7    | 1, 2      | 11, 12, 13 |
| 8    | 1, 2      | 12, 13 |
| 9    | 1, 2      | 13 |
| 10   | 1, 2      | 14 |
| 11   | 3, 4-7    | 12 |
| 12   | 11        | 13, 14 |
| 13   | 4-9, 12   | 14 |
| 14   | 12, 13    | F1-F4 |
| F1-F4| 14        | — |

### Agent Dispatch Summary

- **Wave 1**: 3 tasks — T1 → `deep`, T2 → `deep`, T3 → `quick`
- **Wave 2**: 5 tasks — T4-T7 → `deep`, T8 → `unspecified-high`
- **Wave 3**: 3 tasks — T9 → `unspecified-high`, T10 → `quick`, T11 → `deep`
- **Wave 4**: 2 tasks — T12 → `unspecified-high`, T13 → `unspecified-high`
- **Wave 5**: 1 task — T14 → `writing`
- **FINAL**: 4 tasks — F1 → `oracle`, F2-F3 → `unspecified-high`, F4 → `deep`

---

## TODOs


- [ ] 1. Read + Classify All Unclear Scripts (Triage)

  **What to do**:
  - Read EVERY LINE of these 7 scripts: `temp_plot.R` (146L), `temp_plot_new.R` (693L), `quick.R` (428L), `geneNMF.R` (?L), `expr.R` (120L), `ref_pipeline.R` (1165L), `heatmap.R` (407L)
  - For each script, document:
    - Primary purpose (1 sentence)
    - Unique functions/analyses not found elsewhere
    - Data inputs (readRDS) and outputs (saveRDS, plots)
    - Whether it overlaps with other known scripts
  - Classify each as: **KEEP** (unique, move to appropriate subfolder), **MERGE** (overlaps with X, combine into X), or **DELETE** (truly redundant, all functionality exists elsewhere)
  - If ANY script's purpose is unclear after full read, flag it with `[ASK USER]` — do NOT guess
  - Write classification to `.sisyphus/evidence/task-1-triage-manifest.md`

  **Must NOT do**:
  - Do not delete any files yet — classification only
  - Do not guess about unclear scripts — flag for user
  - Do not read files outside this list (other batches handle those)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires careful reading of ~3K lines with analytical judgment about code purpose
  - **Skills**: []
    - No special skills needed — pure R code reading and analysis

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 1 (sequential — must complete before all other waves)
  - **Blocks**: Tasks 2, 3, 4, 5, 6, 7, 8
  - **Blocked By**: None (can start immediately)

  **References**:

  **Pattern References**:
  - `analysis/temp_plot.R` — 146 lines, unknown purpose. Name suggests temporary but user says may contain important code
  - `analysis/temp_plot_new.R` — 693 lines, grep shows it produces enrichment PNGs identical to `example_anno.R`. May be a newer version.
  - `analysis/quick.R` — 428 lines, reads `all_outs.rds`, `meta.RData`, `cell_summary.rds`. Libraries: data.table, infercna, readxl.
  - `analysis/geneNMF.R` — Root-level copy. AGENTS.md mentions it "explores NMF results". May be empty or a stub.
  - `analysis/expr.R` — 120 lines, only loads `patchwork`. Likely a plotting fragment.
  - `analysis/ref_pipeline.R` — 1165 lines, loads Seurat+DoubletFinder+harmony. Appears to be an early pipeline draft. May contain unique functions.
  - `analysis/heatmap.R` — 407 lines, same library set as ref_pipeline.R. User confirmed it contains `plot_heatmap()` function NOT in main pipeline.

  **Acceptance Criteria**:
  - [ ] All 7 files read completely (every line)
  - [ ] Classification document produced at `.sisyphus/evidence/task-1-triage-manifest.md`
  - [ ] Each file has: purpose, unique code, inputs/outputs, overlap assessment, classification (KEEP/MERGE/DELETE)
  - [ ] Any unclear files flagged with [ASK USER]

  **QA Scenarios:**

  ```
  Scenario: Triage manifest completeness
    Tool: Bash
    Preconditions: Task 1 completed
    Steps:
      1. Read `.sisyphus/evidence/task-1-triage-manifest.md`
      2. Verify it contains sections for all 7 files
      3. Verify each section has: purpose, unique_code, inputs, outputs, overlaps, classification
      4. Verify no [ASK USER] flags remain unresolved (if any exist, task pauses for user input)
    Expected Result: 7/7 files classified, document well-structured
    Evidence: .sisyphus/evidence/task-1-triage-completeness.txt
  ```

  **Commit**: NO (triage only, no file changes)

---

- [ ] 2. Design Final Subfolder Structure + File Manifest

  **What to do**:
  - Using results from Task 1 (triage) + the pre-existing research on all other files, produce the COMPLETE file manifest
  - The manifest maps EVERY original file (38 R + 17 sh) to its final destination:
    - `analysis/metaprograms/` — MP analysis scripts (sc, pdo, correlation, scoring)
    - `analysis/enrichment/` — gene overlap, enrichment, annotation
    - `analysis/clinical/` — TCGA, survival, CIBERSORTx
    - `analysis/cnv/` — CNV profiling, subsetting, plotting
    - `analysis/cell_states/` — state assignment, QC, UMAP, annotation
    - `analysis/plotting/` — general-purpose plotting (beaut_umap, etc.)
    - `analysis/non_malignant_nmf/` — renamed from non_mali_nmf/, parameterized shells
    - `analysis/developmental/` — developmental programme analysis
    - `analysis/utils.R` — shared utility file at analysis root
  - For each file, specify: `original_path → action (move/merge/delete) → final_path`
  - For MERGE actions, specify which files merge into what, and the new descriptive filename
  - Propose descriptive final filenames. Examples:
    - `MP_analysis_sc.R` + `MP_analysis_pdos.R` → `analysis/metaprograms/mp_correlation_heatmap.R` (with clearly separated sc/pdo sections)
    - `states_qccheck.R` + `qc_unresolved_states.R` → `analysis/cell_states/states_qc.R`
  - Save manifest to `.sisyphus/evidence/task-2-file-manifest.md`

  **Must NOT do**:
  - Do not move or create any files yet — manifest only
  - Do not assign files to subfolders without considering their actual content (from Task 1 + research)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Requires synthesizing triage results + research findings into a coherent structure
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO (depends on Task 1)
  - **Parallel Group**: Wave 1 (sequential)
  - **Blocks**: Tasks 3, 4, 5, 6, 7, 8, 9, 10
  - **Blocked By**: Task 1

  **References**:

  **Pattern References**:
  - `.sisyphus/evidence/task-1-triage-manifest.md` — Task 1 output with classifications for 7 unclear files
  - `.sisyphus/drafts/analysis-reorg.md` — Interview draft with all research findings, line counts, overlap analysis
  - `AGENTS.md` "Confirmed Analysis Scripts" section — documented purposes for key scripts

  **Acceptance Criteria**:
  - [ ] Manifest accounts for all 55 original files (38 R + 17 sh)
  - [ ] Every file has: original_path, action, final_path
  - [ ] All merge groups specify the new combined filename
  - [ ] All subfolder assignments are consistent with theme groupings
  - [ ] Manifest saved to `.sisyphus/evidence/task-2-file-manifest.md`

  **QA Scenarios:**

  ```
  Scenario: Manifest covers all original files
    Tool: Bash
    Preconditions: Task 2 completed, manifest exists
    Steps:
      1. Count unique original_path entries in manifest
      2. Compare against `find analysis/ -type f \( -name '*.R' -o -name '*.sh' \) | wc -l`
      3. Verify counts match
    Expected Result: 55 files accounted for, 0 unaccounted
    Evidence: .sisyphus/evidence/task-2-manifest-audit.txt
  ```

  **Commit**: NO (manifest only, no file changes)

---

- [ ] 3. Create analysis/utils.R Skeleton

  **What to do**:
  - Create `analysis/utils.R` with:
    - Header comment explaining purpose ("Shared utility functions for analysis scripts")
    - Library imports section (common libraries used across 3+ scripts: dplyr, ggplot2, ComplexHeatmap, circlize, Seurat)
    - Placeholder function stubs for patterns identified in research:
      - `filter_silhouette(geneNMF_obj, threshold = 0)` — silhouette filtering (used in 7+ scripts)
      - `plot_correlation_heatmap(cor_matrix, ...)` — ComplexHeatmap wrapper for MP correlation
      - `score_ucell(seurat_obj, gene_lists, ...)` — UCell scoring wrapper
      - `load_epi_data()` — standardized loading of EAC_Ref_epi.rds + metadata
    - Each stub has a `# TODO: Extract from Task N` comment
  - This skeleton is populated in Task 11 after the actual scripts are read and merged

  **Must NOT do**:
  - Do not implement full functions yet — stubs only
  - Do not add functions used by fewer than 3 scripts (avoid premature DRY)

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Simple file creation with placeholder stubs
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (after Task 2)
  - **Parallel Group**: Wave 1 (can overlap with start of Wave 2)
  - **Blocks**: Task 11
  - **Blocked By**: Task 2 (needs manifest to know which utils are needed)

  **References**:

  **Pattern References**:
  - `.sisyphus/evidence/task-2-file-manifest.md` — manifest showing which scripts share patterns
  - Research finding: 7 scripts do silhouette filtering, 11 use UCell, 25 produce heatmaps

  **Acceptance Criteria**:
  - [ ] `analysis/utils.R` exists with header comment, library section, and 4+ stub functions
  - [ ] `Rscript -e "parse(file='analysis/utils.R')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: utils.R is valid R and contains expected stubs
    Tool: Bash
    Preconditions: Task 3 completed
    Steps:
      1. Run `Rscript -e "parse(file='analysis/utils.R')"` — verify exit code 0
      2. Grep for `filter_silhouette`, `plot_correlation_heatmap`, `score_ucell`, `load_epi_data`
      3. Verify all 4 function names present
    Expected Result: Parse succeeds, all 4 stubs found
    Evidence: .sisyphus/evidence/task-3-utils-skeleton.txt
  ```

  **Commit**: YES (group with Wave 1)
  - Message: `chore(analysis): triage scripts and design final structure`
  - Files: `analysis/utils.R`, `.sisyphus/evidence/task-1-*.md`, `.sisyphus/evidence/task-2-*.md`

---

- [ ] 4. Metaprograms — Read, Merge, Restructure (6 scripts, ~3.9K lines)

  **What to do**:
  - Read EVERY LINE of: `MP_analysis_sc.R` (827L), `MP_analysis_pdos.R` (898L), `compare_pdos_sc.R` (911L), `Auto_MP_correlation_v2.R` (482L), `Auto_MP19_analysis.R` (389L), `score_other_MPs.R` (407L)
  - Create `analysis/metaprograms/` directory
  - Design merged script structure. Proposed split (adjust based on actual code):
    - `mp_correlation_sc.R` — Single-cell MP correlation analysis (from MP_analysis_sc.R + relevant parts of compare_pdos_sc.R)
    - `mp_correlation_pdo.R` — PDO MP correlation analysis (from MP_analysis_pdos.R + relevant parts of compare_pdos_sc.R)
    - `mp_correlation_crossdata.R` — Cross-dataset comparison: SC vs PDO (unique parts of compare_pdos_sc.R that compare between datasets)
    - `mp_ucell_scoring.R` — UCell scoring + silhouette filtering (from Auto_MP19_analysis.R + score_other_MPs.R)
    - `mp_database_correlation.R` — MP-to-database-term correlation (from Auto_MP_correlation_v2.R)
  - **CRITICAL**: PDO data has different naming conventions and sources than SC. Keep PDO and SC logic in SEPARATE files, not mixed into one mega-script
  - Deduplicate shared logic: silhouette filtering, heatmap generation, Fisher Z meta-analysis → call from utils.R
  - Each merged script gets a header comment block listing which original scripts were merged
  - Handle `setwd()` conflicts: if scripts use different working directories, use full paths or wrap in local setwd/restore
  - Move `Auto_MP_correlation.sh` alongside its R script, update the `Rscript` path inside it

  **Must NOT do**:
  - Do not combine PDO and SC logic into a single function with a `type` argument — they are structurally different
  - Do not drop any unique analysis (e.g., Jaccard similarity in compare_pdos_sc.R)
  - Do not rename the shell script's job name or PBS resources

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Most complex merge batch — 3.9K lines, 6 files, subtle differences between SC/PDO logic
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 5, 6, 7, 8)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 11, 12, 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/MP_analysis_sc.R` — 827L. SC MP correlation using Spearman + Fisher Z meta-analysis + ComplexHeatmap. Reads EAC_Ref_epi.rds, UCell_default.rds, MP_outs_default.rds
  - `analysis/MP_analysis_pdos.R` — 898L. PDO version of above. Reads PDOs_final.rds and PDO-specific MP_outs_default.rds from PDOs_Pipeline/
  - `analysis/compare_pdos_sc.R` — 911L. Cross-comparison between PDO and SC MPs. Reads both UCell_default.rds AND UCell_pdos.rds. Uses Jaccard similarity + pheatmap.
  - `analysis/Auto_MP_correlation_v2.R` — 482L. Correlates MPs against database terms (Hallmarks, GO, 3CA). Reads cluster_enrich.rds, UCell_nMP19_filtered.rds, geneNMF_metaprograms_nMP_19.rds
  - `analysis/Auto_MP19_analysis.R` — 389L. UCell scoring for nMP=19 with silhouette filtering + Jaccard + heatmap. Writes UCell_nMP19_filtered.rds
  - `analysis/score_other_MPs.R` — 407L. Scores cells with external MP gene sets (3CA). Reads UCell_3CA_MPs.rds
  - `analysis/Auto_MP_correlation.sh` — PBS wrapper: 8 CPUs, 128GB, 4h. Calls `Rscript analysis/Auto_MP_correlation_v2.R`

  **API/Type References**:
  - `ref_outs/EAC_Ref_epi.rds` — 75,348 OAC epithelial cells Seurat object (central data)
  - `ref_outs/UCell_default.rds` — 75348 x 9 MP UCell scores
  - `ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` — Primary MP object
  - `ref_outs/MP_outs_default.rds` — Processed MP outputs

  **External References**:
  - `.sisyphus/evidence/task-2-file-manifest.md` — Manifest specifying final filenames for this group

  **Acceptance Criteria**:
  - [ ] All 6 original scripts fully read
  - [ ] `analysis/metaprograms/` directory created with reorganized scripts
  - [ ] PDO and SC analyses in separate files
  - [ ] All unique analyses preserved (Jaccard similarity, Fisher Z, database term correlation, 3CA scoring)
  - [ ] `Auto_MP_correlation.sh` moved to `analysis/metaprograms/` with updated Rscript path
  - [ ] Each merged script has header comment listing source files
  - [ ] `Rscript -e "parse(file='<each_new_script>')"` exits 0 for all new scripts

  **QA Scenarios:**

  ```
  Scenario: All metaprogram scripts parse correctly
    Tool: Bash
    Preconditions: Task 4 completed, metaprograms/ directory exists
    Steps:
      1. Run `find analysis/metaprograms/ -name '*.R' -exec Rscript -e "parse(file='{}')" \;`
      2. Verify all exit 0
      3. Run `bash -n analysis/metaprograms/Auto_MP_correlation.sh` if .sh was moved
    Expected Result: All scripts parse without error
    Evidence: .sisyphus/evidence/task-4-metaprograms-parse.txt

  Scenario: No functionality lost — key patterns preserved
    Tool: Bash (grep)
    Preconditions: Task 4 completed
    Steps:
      1. Grep all files in analysis/metaprograms/ for 'Jaccard|jaccard' — must find at least 1 match
      2. Grep for 'fisher.*z|Fisher.*Z|fisherz' — must find at least 1 match
      3. Grep for 'ScoreSignatures_UCell|UCell' — must find at least 2 matches
      4. Grep for 'Spearman|spearman' — must find at least 1 match
      5. Grep for 'UCell_3CA|3CA' — must find at least 1 match
    Expected Result: All 5 key patterns found in reorganized scripts
    Evidence: .sisyphus/evidence/task-4-metaprograms-functionality.txt

  Scenario: Shell script points to valid R file
    Tool: Bash
    Preconditions: Task 4 completed
    Steps:
      1. Extract Rscript path from the moved .sh file
      2. Verify that file exists at the extracted path
    Expected Result: Rscript target exists
    Evidence: .sisyphus/evidence/task-4-shell-integrity.txt
  ```

  **Commit**: YES (group with Wave 2)
  - Message: `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders`
  - Files: `analysis/metaprograms/*`
  - Pre-commit: `find analysis/metaprograms/ -name '*.R' -exec Rscript -e "parse(file='{}')" \;`

---

- [ ] 5. Enrichment — Read, Merge, Restructure (5 scripts, ~1.1K lines)

  **What to do**:
  - Read EVERY LINE of: `example_anno.R` (315L), `enrich_plot.R` (270L), `terms_overlap.R` (116L), `wnt_enrich.R` (139L), `nmf_plot.R` (304L)
  - Create `analysis/enrichment/` directory
  - Design merged structure:
    - `enrichment_analysis.R` — Gene overlap enrichment computation (from terms_overlap.R + core logic of example_anno.R)
    - `enrichment_plotting.R` — All enrichment visualization (from enrich_plot.R + plotting parts of example_anno.R + wnt_enrich.R enrichment plots + nmf_plot.R enrichment heatmaps)
  - If `wnt_enrich.R` does WNT-specific analysis beyond generic enrichment, keep it as a separate `wnt_pathway.R`
  - If `nmf_plot.R` does NMF-specific plotting beyond enrichment, parts may belong in a different subfolder — split accordingly
  - Extract shared enrichment patterns to utils.R stubs (note for Task 11)
  - Handle setwd() — `enrich_plot.R` reads from an external path (`EAC_Ref_all/00_merged/developmental/`)

  **Must NOT do**:
  - Do not lose the WNT-specific analysis if it’s unique
  - Do not merge scripts that use fundamentally different enrichment approaches (clusterProfiler vs custom overlap)

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: 5 scripts to read and intelligently merge, need to distinguish enrichment computation from plotting
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 4, 6, 7, 8)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 11, 12, 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/example_anno.R` — 315L. MP annotation via enrichment against Hallmark, GO, 3CA, Pan-Cancer, etc. Reads geneNMF_metaprograms_nMP_19.rds. Produces 9 enrichment heatmap PNGs.
  - `analysis/enrich_plot.R` — 270L. Reads enrich_dev.rds from external path. Produces Enrichment_Analysis.pdf. Uses clusterProfiler + pheatmap.
  - `analysis/terms_overlap.R` — 116L. Computes gene overlap enrichment → cluster_enrich.rds. Libraries: dplyr, tidyr, purrr, pheatmap.
  - `analysis/wnt_enrich.R` — 139L. WNT pathway enrichment. Libraries: pheatmap, dplyr.
  - `analysis/nmf_plot.R` — 304L. NMF result plotting with enrichment. Libraries: clusterProfiler, msigdbr, enrichplot, pheatmap.

  **Acceptance Criteria**:
  - [ ] All 5 original scripts fully read
  - [ ] `analysis/enrichment/` directory created with reorganized scripts
  - [ ] All unique enrichment analyses preserved
  - [ ] `Rscript -e "parse(file='<each_new_script>')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: All enrichment scripts parse correctly
    Tool: Bash
    Preconditions: Task 5 completed
    Steps:
      1. Run `find analysis/enrichment/ -name '*.R' -exec Rscript -e "parse(file='{}')" \;`
      2. Verify all exit 0
    Expected Result: All scripts parse without error
    Evidence: .sisyphus/evidence/task-5-enrichment-parse.txt

  Scenario: Key enrichment patterns preserved
    Tool: Bash (grep)
    Steps:
      1. Grep analysis/enrichment/ for 'enrichGO|enrichKEGG|enricher|clusterProfiler' — at least 1 match
      2. Grep for 'pheatmap|Heatmap\(' — at least 1 match
      3. Grep for 'Hallmark|hallmark' — at least 1 match
    Expected Result: Core enrichment functions present
    Evidence: .sisyphus/evidence/task-5-enrichment-functionality.txt
  ```

  **Commit**: YES (group with Wave 2)
  - Message: `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders`
  - Files: `analysis/enrichment/*`

---

- [ ] 6. Clinical/Survival — Read, Clean, Restructure (3 scripts, ~1.6K lines)

  **What to do**:
  - Read EVERY LINE of: `cibersort_result.R` (963L), `evaluate_clinical_MPs.R` (400L), `TCGA_data.R` (208L)
  - Create `analysis/clinical/` directory
  - Design structure:
    - `tcga_data_prep.R` — TCGA data download/processing (from TCGA_data.R). Produces tcga_esca_meta.rds + tcga_esca_tpm_matrix.rds
    - `survival_analysis.R` — KM survival analysis combining CIBERSORTx, GSVA, and MP-based approaches (merge cibersort_result.R + evaluate_clinical_MPs.R)
    - OR keep them separate if they do fundamentally different things — decide after reading
  - **setwd() warning**: cibersort_result.R uses paths to external project dirs. Must use full paths in merged version.
  - Both survival scripts read TCGA metadata — potential for shared data loading via utils.R

  **Must NOT do**:
  - Do not merge if survival approaches are fundamentally different (one may be CIBERSORTx deconvolution, other may be GSVA scoring)
  - Do not change TCGA data paths — these point to external data

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Complex survival analysis code, TCGA integration, must understand medical statistics context
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 4, 5, 7, 8)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 11, 12, 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/cibersort_result.R` — 963L. KM survival across multiple split methods. Reads tcga_esca_meta.rds, MP_outs_default.rds, EAC_Ref_epi.rds, state_temp.rds, states_degs.rds. Uses GSVA + survival + survminer.
  - `analysis/evaluate_clinical_MPs.R` — 400L. Clinical evaluation of MPs. Reads TCGA meta from external spatialtranscriptomics path. Uses survival + survminer + GSVA + UCell.
  - `analysis/TCGA_data.R` — 208L. TCGA data processing. Writes tcga_esca_meta.rds + tcga_esca_tpm_matrix.rds. Dependencies: readr, dplyr, stringr, purrr, org.Hs.eg.db.
  - AGENTS.md: "Always filter states_degs.rds with `group_by(cluster) %>% slice_max(n=100)`"

  **Acceptance Criteria**:
  - [ ] All 3 original scripts fully read
  - [ ] `analysis/clinical/` directory created
  - [ ] TCGA data prep script separated from survival analysis
  - [ ] External path references preserved correctly (no broken paths)
  - [ ] `Rscript -e "parse(file='<each_new_script>')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: Clinical scripts parse correctly
    Tool: Bash
    Preconditions: Task 6 completed
    Steps:
      1. Run `find analysis/clinical/ -name '*.R' -exec Rscript -e "parse(file='{}')" \;`
      2. Verify all exit 0
    Expected Result: All scripts parse without error
    Evidence: .sisyphus/evidence/task-6-clinical-parse.txt

  Scenario: TCGA paths preserved
    Tool: Bash (grep)
    Steps:
      1. Grep analysis/clinical/ for 'tcga_esca_meta' — must find reference
      2. Grep for 'survfit|Surv\(' — must find survival analysis code
    Expected Result: TCGA data and survival analysis preserved
    Evidence: .sisyphus/evidence/task-6-clinical-functionality.txt
  ```

  **Commit**: YES (group with Wave 2)
  - Message: `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders`
  - Files: `analysis/clinical/*`

---

- [ ] 7. CNV + Cell States — Read, Merge QC Pair, Restructure (6 scripts, ~1.5K lines)

  **What to do**:
  - Read EVERY LINE of: `cnv_profile_sc.R` (68L), `cnv_subset.R` (328L), `plot_CNV.R` (191L), `states_umap.R` (166L), `states_qccheck.R` (357L), `qc_unresolved_states.R` (354L)
  - Create `analysis/cnv/` and `analysis/cell_states/` directories
  - **CNV scripts** → `analysis/cnv/`:
    - Assess whether the 3 CNV scripts can merge into 1-2 files or stay separate (they may handle different stages: profiling → subsetting → plotting)
    - Rename descriptively: e.g., `cnv_profiling.R`, `cnv_subsetting.R`, `cnv_plotting.R` (or merged as appropriate)
  - **Cell state scripts** → `analysis/cell_states/`:
    - `states_qccheck.R` and `qc_unresolved_states.R` produce IDENTICAL output filenames (`states_status_quality_comparison.pdf`, `quality_comparison_by_original_state.pdf`) — likely near-duplicate. Merge into one `states_qc.R`
    - `states_umap.R` → `states_umap.R` (rename descriptively if needed)
  - All CNV scripts share identical library sets (data.table, dplyr, ComplexHeatmap, circlize, RColorBrewer, Seurat, infercna) — note for utils.R

  **Must NOT do**:
  - Do not merge CNV scripts blindly — they may represent a sequential workflow (profile → subset → plot)
  - Do not assume the QC pair are identical — read both fully, one may handle unresolved states specifically

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Must carefully compare QC pair and understand CNV workflow sequence
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 4, 5, 6, 8)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 11, 12, 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/cnv_profile_sc.R` — 68L. Reads EAC_Ref_filtered.rds (external). Writes all_outs.rds. Produces clustering_comparison_epithelial.png.
  - `analysis/cnv_subset.R` — 328L. Same library set as cnv_profile. No visible readRDS/saveRDS in grep — may use pre-loaded objects.
  - `analysis/plot_CNV.R` — 191L. Produces PDOs_cnv.pdf. Same library set.
  - `analysis/states_qccheck.R` — 357L. Produces states_status_quality_comparison.pdf + quality_comparison_by_original_state.pdf.
  - `analysis/qc_unresolved_states.R` — 354L. Produces SAME output filenames as states_qccheck.R. Likely a refined version.
  - `analysis/states_umap.R` — 166L. Produces columnar_states_UMAP.pdf.

  **Acceptance Criteria**:
  - [ ] All 6 original scripts fully read
  - [ ] `analysis/cnv/` and `analysis/cell_states/` directories created
  - [ ] QC pair merged into single script (or kept separate with justification)
  - [ ] All CNV stages preserved
  - [ ] `Rscript -e "parse(file='<each_new_script>')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: CNV and cell state scripts parse correctly
    Tool: Bash
    Preconditions: Task 7 completed
    Steps:
      1. Run `find analysis/cnv/ analysis/cell_states/ -name '*.R' -exec Rscript -e "parse(file='{}')" \;`
      2. Verify all exit 0
    Expected Result: All scripts parse without error
    Evidence: .sisyphus/evidence/task-7-cnv-states-parse.txt

  Scenario: QC merge preserves both quality check types
    Tool: Bash (grep)
    Steps:
      1. Grep analysis/cell_states/ for 'unresolved' — must be present (from qc_unresolved_states.R)
      2. Grep for 'quality_comparison' — must be present
    Expected Result: Both QC analysis types present in merged script
    Evidence: .sisyphus/evidence/task-7-qc-merge-check.txt
  ```

  **Commit**: YES (group with Wave 2)
  - Message: `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders`
  - Files: `analysis/cnv/*`, `analysis/cell_states/*`

---

- [ ] 8. Plotting + Standalone Scripts — Read, Rename, Move (~3.6K lines)

  **What to do**:
  - Read EVERY LINE of the remaining scripts not handled by other tasks. Based on Task 1 triage + Task 2 manifest, this likely includes:
    - `beaut_umap.R` (428L) — publication-quality UMAPs
    - `gene_expr_compare.R` (344L) — gene expression comparison heatmaps
    - `residual.R` (821L) — unclear from name, must read to determine
    - `plot_clinical.R` (342L) — clinical variable plotting
    - `celltyping.R` (198L) — cell type assignment
    - `annotation.R` (508L) — cell annotation
    - `summary.R` (696L) — cross-sample summary statistics
    - `sum_cancer.R` (183L) — cancer cell summary
    - Plus any files from Task 1 classified as KEEP
  - Assign each to the appropriate subfolder based on actual content:
    - General plotting scripts → `analysis/plotting/`
    - Cell annotation/typing → `analysis/cell_states/` (if Task 7 created it) or own subfolder
    - Summary statistics → `analysis/` root or a `analysis/summary/` subfolder
    - Expression analysis → decide based on content (may be metaprograms-adjacent or standalone)
  - Rename each file descriptively if current name is unclear
  - `residual.R` (821 lines!) is a major unknown — read carefully, may need its own category or belong in metaprograms/
  - For files classified as DELETE by Task 1: remove the original after confirming nothing unique is lost
  - For files classified as MERGE by Task 1: merge into the appropriate target script

  **Must NOT do**:
  - Do not create a "misc" or "other" subfolder — every file gets a purposeful home
  - Do not delete files without Task 1 confirming they're redundant
  - Do not move files that have unclear purpose without reading them first

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Many diverse files to process, but each is simpler than the MP merge batch
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 4, 5, 6, 7)
  - **Parallel Group**: Wave 2
  - **Blocks**: Tasks 12, 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `.sisyphus/evidence/task-1-triage-manifest.md` — Classification of unclear files from Task 1
  - `.sisyphus/evidence/task-2-file-manifest.md` — Master manifest with all file assignments
  - `analysis/beaut_umap.R` — 428L. Libraries: ggplot2, scattermore, pals, scales, cowplot.
  - `analysis/gene_expr_compare.R` — 344L. Produces Heatmap_Matched_Colors_Dediff.pdf. Libraries: ComplexHeatmap, Seurat.
  - `analysis/residual.R` — 821L. Produces MP_heatmap_states_subset_z_residual.pdf + 2 more heatmaps. Libraries: ComplexHeatmap, proxy, Seurat. Likely MP-adjacent.
  - `analysis/plot_clinical.R` — 342L. Produces CellState_PerVariable_ByStudy.pdf. Libraries: Seurat, readxl, patchwork.
  - `analysis/annotation.R` — 508L. Libraries: Seurat, dplyr, purrr, tidyr, ggplot2.
  - `analysis/celltyping.R` — 198L. Libraries: readxl, dplyr, purrr, stringr, Seurat.
  - `analysis/summary.R` — 696L. Cross-sample summary. Reads EAC_Ref_filtered.rds.
  - `analysis/sum_cancer.R` — 183L. Cancer cell summary.

  **Acceptance Criteria**:
  - [ ] All files in this batch fully read
  - [ ] Each file moved to appropriate themed subfolder with descriptive name
  - [ ] Files from Task 1 classified as DELETE removed (with justification logged)
  - [ ] Files from Task 1 classified as MERGE incorporated into target scripts
  - [ ] `Rscript -e "parse(file='<each_new_script>')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: All moved scripts parse correctly
    Tool: Bash
    Preconditions: Task 8 completed
    Steps:
      1. For each new/moved script in analysis/ subfolders, run parse check
      2. Verify all exit 0
    Expected Result: All scripts parse without error
    Evidence: .sisyphus/evidence/task-8-standalone-parse.txt

  Scenario: No orphaned files in analysis/ root
    Tool: Bash
    Steps:
      1. List all .R files remaining in analysis/ root (not in subfolders)
      2. Only `utils.R` should remain at root level (plus any justified exceptions)
    Expected Result: Only utils.R at analysis/ root
    Evidence: .sisyphus/evidence/task-8-root-cleanup.txt
  ```

  **Commit**: YES (group with Wave 2)
  - Message: `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders`
  - Files: All remaining moved/renamed scripts

---

- [ ] 9. Non-Malignant NMF — Parameterize Shells, Rename Folder

  **What to do**:
  - Rename `analysis/non_mali_nmf/` to `analysis/non_malignant_nmf/`
  - Read the 2 R scripts: `Auto_geneNMF_celltype.R` (196L) and `Auto_anno_celltype.R` (383L)
  - Read 2-3 representative shell scripts to confirm they are identical except for cell type and PBS resources
  - Read both master scripts: `Auto_geneNMF_master.sh` and `Auto_anno_master.sh`
  - **Create parameterized shell scripts**:
    - `run_geneNMF.sh` — single PBS script accepting cell type + resource args (replaces 7 per-cell-type scripts)
    - `run_annotation.sh` — single PBS script accepting cell type + resource args (replaces 7 per-cell-type scripts)
    - Preserve the varying PBS resource allocations by accepting them as arguments or using a case/switch for known cell types
  - **Update master scripts**:
    - `Auto_geneNMF_master.sh` → `geneNMF_master.sh` — loop submits `run_geneNMF.sh` with cell type + resources as vars
    - `Auto_anno_master.sh` → `annotation_master.sh` — same pattern
  - Delete the 14 per-cell-type .sh files after parameterized versions are created and verified
  - Rename R scripts descriptively if needed (current names are already clear)

  **Must NOT do**:
  - Do not change PBS resource allocations — different cell types need different memory/walltime
  - Do not modify the R scripts' logic — only touch shell wrappers
  - Do not break the `cd $WD` + relative Rscript path pattern — PBS jobs depend on this

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Shell scripting with PBS-specific patterns, must handle variable resource allocation
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 10, 11)
  - **Parallel Group**: Wave 3
  - **Blocks**: Task 13
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/non_mali_nmf/Auto_geneNMF_celltype.R` — 196L. Takes celltype arg via commandArgs(). Handles 7 cell types.
  - `analysis/non_mali_nmf/Auto_anno_celltype.R` — 383L. Takes celltype arg. Enrichment annotation.
  - `analysis/non_mali_nmf/Auto_geneNMF_master.sh` — Orchestrates PBS submissions for each cell type.
  - `analysis/non_mali_nmf/Auto_anno_master.sh` — Same pattern, checks for prerequisite MP files.
  - Per-cell-type .sh files: differ only in PBS `#PBS -l select=1:ncpus=N:mem=Mgb` and `#PBS -l walltime=HH:MM:SS`.
  - PBS template from AGENTS.md — standard structure to follow for parameterized script.

  **Acceptance Criteria**:
  - [ ] `analysis/non_malignant_nmf/` directory exists (renamed from non_mali_nmf)
  - [ ] 2 parameterized .sh scripts created (run_geneNMF.sh, run_annotation.sh)
  - [ ] 2 updated master scripts created
  - [ ] 14 per-cell-type .sh files removed
  - [ ] `bash -n` passes for all new .sh scripts
  - [ ] Master scripts correctly reference parameterized scripts
  - [ ] PBS resource allocations preserved per cell type

  **QA Scenarios:**

  ```
  Scenario: Parameterized shell scripts have valid syntax
    Tool: Bash
    Preconditions: Task 9 completed
    Steps:
      1. Run `bash -n analysis/non_malignant_nmf/run_geneNMF.sh`
      2. Run `bash -n analysis/non_malignant_nmf/run_annotation.sh`
      3. Run `bash -n analysis/non_malignant_nmf/geneNMF_master.sh`
      4. Run `bash -n analysis/non_malignant_nmf/annotation_master.sh`
    Expected Result: All 4 scripts pass syntax check
    Evidence: .sisyphus/evidence/task-9-shell-syntax.txt

  Scenario: Per-cell-type scripts removed
    Tool: Bash
    Steps:
      1. Count .sh files in analysis/non_malignant_nmf/
      2. Should be exactly 4 (2 runners + 2 masters), not 18
    Expected Result: 4 shell scripts total
    Evidence: .sisyphus/evidence/task-9-shell-count.txt

  Scenario: Resource allocations preserved
    Tool: Bash (grep)
    Steps:
      1. Grep parameterized scripts for cell type resource handling (case/if or variable substitution)
      2. Verify macrophage, fibroblast, endothelial, nk, plasma, cd4, cd8 all handled
    Expected Result: All 7 cell types have resource specifications
    Evidence: .sisyphus/evidence/task-9-resources-preserved.txt
  ```

  **Commit**: YES (group with Wave 3)
  - Message: `refactor(analysis): parameterize non-mali-nmf shells, populate utils.R`
  - Files: `analysis/non_malignant_nmf/*`

---

- [ ] 10. Developmental — Assess and Restructure

  **What to do**:
  - Read EVERY LINE of `analysis/developmental/developmental.R` (403L)
  - The script loads early embryogenesis markers from an Excel sheet and creates term2gene mappings
  - Rename `developmental.R` to something more descriptive if warranted (e.g., `embryogenesis_markers.R` or keep as `developmental.R` if name is clear enough)
  - Keep in `analysis/developmental/` subfolder (already correctly organized)
  - Check for any overlap with enrichment scripts (terms_overlap.R, example_anno.R) — this script may produce inputs used by enrichment

  **Must NOT do**:
  - Do not merge into enrichment unless the logic is truly duplicative

  **Recommended Agent Profile**:
  - **Category**: `quick`
    - Reason: Single file to read and potentially rename. Simple task.
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Tasks 9, 11)
  - **Parallel Group**: Wave 3
  - **Blocks**: Task 14
  - **Blocked By**: Tasks 1, 2

  **References**:

  **Pattern References**:
  - `analysis/developmental/developmental.R` — 403L. Loads early embryogenesis markers from Excel. Creates term2gene mappings for lineages (ICM, Epiblast, etc.).
  - `analysis/enrichment/` — enrichment scripts may consume the term2gene mappings this produces

  **Acceptance Criteria**:
  - [ ] Script fully read
  - [ ] Renamed descriptively if needed
  - [ ] Relationship to enrichment scripts documented
  - [ ] `Rscript -e "parse(file='<script>')"` exits 0

  **QA Scenarios:**

  ```
  Scenario: Developmental script parses correctly
    Tool: Bash
    Steps:
      1. Run `Rscript -e "parse(file='analysis/developmental/<final_name>.R')"` 
      2. Verify exit 0
    Expected Result: Script parses without error
    Evidence: .sisyphus/evidence/task-10-developmental-parse.txt
  ```

  **Commit**: YES (group with Wave 3)
  - Message: `refactor(analysis): parameterize non-mali-nmf shells, populate utils.R`
  - Files: `analysis/developmental/*`

---

- [ ] 11. Populate utils.R with Extracted Shared Functions

  **What to do**:
  - Now that all theme batches (Tasks 4-8) have been completed and all scripts have been read, identify the ACTUAL shared patterns
  - Populate the utils.R skeleton (from Task 3) with real implementations extracted from the merged/reorganized scripts:
    - `filter_silhouette(geneNMF_obj, threshold = 0)` — extract from whichever script has the cleanest implementation
    - `plot_correlation_heatmap(cor_matrix, ...)` — extract the ComplexHeatmap wrapper pattern used across MP scripts
    - `score_ucell(seurat_obj, gene_lists, ...)` — UCell scoring wrapper if the pattern is consistent
    - `load_epi_data(path = "ref_outs/EAC_Ref_epi.rds")` — standardized data loading
    - Any additional patterns that emerged during Tasks 4-8 (e.g., Fisher Z transform, Jaccard similarity)
  - **Only extract patterns used in 3+ scripts** — avoid premature DRY
  - Each function gets a roxygen-style comment block explaining parameters and return value
  - Document which scripts will source() this file

  **Must NOT do**:
  - Do not extract patterns used by fewer than 3 scripts
  - Do not create overly abstract wrappers that obscure what's happening
  - Do not change function signatures from how they're used in the source scripts

  **Recommended Agent Profile**:
  - **Category**: `deep`
    - Reason: Must identify common patterns across all reorganized scripts and extract clean implementations
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO (needs all Wave 2 tasks completed)
  - **Parallel Group**: Wave 3 (but depends on Wave 2 completion)
  - **Blocks**: Task 12
  - **Blocked By**: Tasks 3, 4, 5, 6, 7 (needs to see all scripts)

  **References**:

  **Pattern References**:
  - `analysis/utils.R` — Skeleton from Task 3 with stub functions
  - All reorganized scripts from Tasks 4-8 — source of actual implementations to extract
  - Research: 7 scripts use silhouette filtering, 11 use UCell, 8 read EAC_Ref_epi.rds

  **Acceptance Criteria**:
  - [ ] utils.R contains real implementations (not stubs) for 4+ shared functions
  - [ ] Each function has roxygen-style documentation
  - [ ] Only patterns used in 3+ scripts are extracted
  - [ ] `Rscript -e "parse(file='analysis/utils.R')"` exits 0
  - [ ] Functions match how they're called in the source scripts

  **QA Scenarios:**

  ```
  Scenario: utils.R parses and contains real functions
    Tool: Bash
    Preconditions: Task 11 completed
    Steps:
      1. Run `Rscript -e "parse(file='analysis/utils.R')"` — exit 0
      2. Grep for 'function\(' — at least 4 function definitions
      3. Grep for '#\' .*@param' — verify roxygen documentation exists
      4. Verify no remaining 'TODO: Extract' comments
    Expected Result: 4+ documented functions, no stubs remaining
    Evidence: .sisyphus/evidence/task-11-utils-populated.txt

  Scenario: Functions match usage in scripts
    Tool: Bash (grep)
    Steps:
      1. For each function name in utils.R, grep all analysis/ scripts for calls to that function
      2. Verify each function is called from at least 3 scripts
    Expected Result: All utils.R functions used by 3+ scripts
    Evidence: .sisyphus/evidence/task-11-utils-usage.txt
  ```

  **Commit**: YES (group with Wave 3)
  - Message: `refactor(analysis): parameterize non-mali-nmf shells, populate utils.R`
  - Files: `analysis/utils.R`

---

- [ ] 12. Wire utils.R into All Scripts That Need It

  **What to do**:
  - For every reorganized script that uses a function now in utils.R, add `source("analysis/utils.R")` at the top (after library() calls)
  - Replace the inline implementation with a call to the utils.R function
  - **CRITICAL setwd() handling**: Scripts that `setwd()` to ref_outs/ will break `source()` relative paths. Options:
    - Use absolute path: `source("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/analysis/utils.R")`
    - Or source before setwd(): place `source()` call before `setwd()` in each script
    - Choose one approach consistently across all scripts
  - Test each modified script with `Rscript -e "parse(file='<script>')"` after adding source() call
  - Update the manifest to reflect which scripts now depend on utils.R

  **Must NOT do**:
  - Do not force utils.R usage where the inline code is clearer or slightly different from the extracted version
  - Do not change function behavior when replacing inline with utils.R call
  - Do not break setwd() patterns

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Systematic find-and-replace across many files, must handle setwd() edge cases
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO (depends on Task 11)
  - **Parallel Group**: Wave 4 (sequential)
  - **Blocks**: Tasks 13, 14
  - **Blocked By**: Task 11

  **References**:

  **Pattern References**:
  - `analysis/utils.R` — Populated from Task 11 with real functions
  - All reorganized scripts in analysis/ subfolders
  - `.sisyphus/evidence/task-11-utils-usage.txt` — Which scripts use which utils functions

  **Acceptance Criteria**:
  - [ ] Every script that uses a utils.R function has `source()` call
  - [ ] No duplicate inline implementations remain for extracted functions
  - [ ] `Rscript -e "parse(file='<each_script>')"` exits 0 for all modified scripts
  - [ ] Consistent source() path approach across all scripts

  **QA Scenarios:**

  ```
  Scenario: All scripts with source() calls parse correctly
    Tool: Bash
    Preconditions: Task 12 completed
    Steps:
      1. Grep all analysis/ .R files for 'source.*utils'
      2. For each file found, run parse check
      3. Verify all exit 0
    Expected Result: All scripts parse, source() paths valid
    Evidence: .sisyphus/evidence/task-12-utils-wiring.txt

  Scenario: No remaining inline implementations of extracted functions
    Tool: Bash (grep)
    Steps:
      1. For each function name in utils.R (e.g., filter_silhouette), grep all analysis/ scripts
      2. Each function name should appear ONLY as a function CALL (not a function definition) in scripts that source utils.R
      3. Function DEFINITIONS should only exist in utils.R itself
    Expected Result: Zero duplicate function definitions outside utils.R
    Evidence: .sisyphus/evidence/task-12-no-duplicates.txt
  ```

  **Commit**: YES (group with Wave 4)
  - Message: `refactor(analysis): wire utils.R into scripts, fix all shell cross-references`
  - Files: All modified .R scripts

---

- [ ] 13. Audit All Shell Script Paths + Fix Cross-References

  **What to do**:
  - Grep EVERY `.sh` file in analysis/ (recursively) for `Rscript` calls
  - For each `Rscript` reference, verify the target `.R` file exists at the specified path
  - Fix any broken paths caused by file moves/renames
  - Key known cross-references to check:
    - `analysis/metaprograms/Auto_MP_correlation.sh` → should reference the renamed/moved R script
    - `analysis/non_malignant_nmf/run_geneNMF.sh` → should reference `Auto_geneNMF_celltype.R`
    - `analysis/non_malignant_nmf/run_annotation.sh` → should reference `Auto_anno_celltype.R`
    - Master scripts → should reference the parameterized runner scripts
  - Also check for any hardcoded paths in R scripts that reference old locations (e.g., `source("analysis/old_path.R")`)
  - Verify `cd $WD` patterns in shell scripts still work with new directory structure

  **Must NOT do**:
  - Do not change PBS resource allocations
  - Do not change working directory patterns in shell scripts

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
    - Reason: Systematic path auditing across many files
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES (can run alongside Task 12 if careful, but safer sequential)
  - **Parallel Group**: Wave 4
  - **Blocks**: Task 14
  - **Blocked By**: Tasks 4-9, 12

  **References**:

  **Pattern References**:
  - All `.sh` files in analysis/ — check every Rscript reference
  - `.sisyphus/evidence/task-2-file-manifest.md` — old → new path mapping
  - AGENTS.md PBS template — expected shell script structure

  **Acceptance Criteria**:
  - [ ] Every `Rscript` call in every `.sh` file resolves to an existing `.R` file
  - [ ] `bash -n` passes for all `.sh` files
  - [ ] No hardcoded references to old paths remain in any file

  **QA Scenarios:**

  ```
  Scenario: All shell-to-R references valid
    Tool: Bash
    Preconditions: Task 13 completed
    Steps:
      1. `grep -rn 'Rscript' analysis/ --include='*.sh'`
      2. For each match, extract the .R file path
      3. Verify each .R file exists
    Expected Result: 0 broken references
    Evidence: .sisyphus/evidence/task-13-shell-path-audit.txt

  Scenario: No old paths remain
    Tool: Bash (grep)
    Steps:
      1. Grep all analysis/ files for old paths: 'analysis/MP_analysis_sc', 'analysis/compare_pdos', 'analysis/Auto_MP_correlation_v2.R' (the old flat paths)
      2. Verify 0 matches (all references updated to new subfolder paths)
    Expected Result: 0 references to old flat paths
    Evidence: .sisyphus/evidence/task-13-no-old-paths.txt
  ```

  **Commit**: YES (group with Wave 4)
  - Message: `refactor(analysis): wire utils.R into scripts, fix all shell cross-references`
  - Files: All modified `.sh` scripts

---

- [ ] 14. Update AGENTS.md with Complete New Structure

  **What to do**:
  - Read current AGENTS.md sections: "Repository Structure", "Downstream Analysis", "Confirmed Analysis Scripts"
  - Rewrite the "Repository Structure" section to reflect the new analysis/ subfolder layout
  - Rewrite "Confirmed Analysis Scripts" section entirely:
    - For each new/merged script, document:
      - Purpose (1-2 sentences)
      - Inputs (which .rds files it reads)
      - Outputs (which .rds / .pdf / .png it produces)
      - Execution method (interactive vs PBS, which conda env)
      - If merged: list which original scripts were combined
  - Update "Key Shared Data Objects" table if any paths changed
  - Add new section documenting `analysis/utils.R` and its functions
  - Update the `non_mali_nmf` references to `non_malignant_nmf`
  - Preserve ALL existing documentation outside the analysis sections
  - Follow the "AGENTS.md Living Document Rule" — append/update, don't delete non-analysis content

  **Must NOT do**:
  - Do not rewrite sections unrelated to analysis/ reorganization
  - Do not remove existing documentation about pipeline steps 1-8
  - Do not change the structure/format of AGENTS.md beyond the analysis sections

  **Recommended Agent Profile**:
  - **Category**: `writing`
    - Reason: Documentation task requiring clear, structured technical writing
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO (final documentation, needs everything else done)
  - **Parallel Group**: Wave 5 (sequential, last implementation task)
  - **Blocks**: F1-F4
  - **Blocked By**: Tasks 12, 13

  **References**:

  **Pattern References**:
  - `AGENTS.md` — Current file, sections to update: Repository Structure, Downstream Analysis, Confirmed Analysis Scripts
  - `.sisyphus/evidence/task-2-file-manifest.md` — Complete old → new path mapping
  - All reorganized scripts in analysis/ — for documenting inputs/outputs/purpose
  - `analysis/utils.R` — Document shared functions

  **Acceptance Criteria**:
  - [ ] AGENTS.md "Repository Structure" reflects new analysis/ layout
  - [ ] Every script in reorganized analysis/ is documented in "Confirmed Analysis Scripts"
  - [ ] Each documented script has: purpose, inputs, outputs, execution method
  - [ ] `analysis/utils.R` documented with function list
  - [ ] `non_malignant_nmf` references updated from `non_mali_nmf`
  - [ ] No non-analysis documentation removed or damaged

  **QA Scenarios:**

  ```
  Scenario: AGENTS.md references match actual file structure
    Tool: Bash
    Preconditions: Task 14 completed
    Steps:
      1. Extract all file paths mentioned in AGENTS.md analysis sections
      2. For each path, verify the file exists on disk
      3. List all .R files in analysis/ recursively
      4. Verify every .R file is mentioned in AGENTS.md
    Expected Result: 100% path validity, 100% file coverage
    Evidence: .sisyphus/evidence/task-14-agents-md-audit.txt

  Scenario: No references to old paths in AGENTS.md
    Tool: Bash (grep)
    Steps:
      1. Grep AGENTS.md for old flat paths: 'analysis/MP_analysis_sc', 'analysis/example_anno', 'analysis/non_mali_nmf'
      2. Verify 0 matches for old paths
    Expected Result: All references use new subfolder paths
    Evidence: .sisyphus/evidence/task-14-no-old-refs.txt
  ```

  **Commit**: YES
  - Message: `docs(agents): update AGENTS.md with reorganized analysis/ structure`
  - Files: `AGENTS.md`

---

## Final Verification Wave (MANDATORY — after ALL implementation tasks)

> 4 review agents run in PARALLEL. ALL must APPROVE. Rejection → fix → re-run.

- [ ] F1. **Plan Compliance Audit** — `oracle`
  Read the plan end-to-end. For each "Must Have": verify implementation exists (check file locations, read merged scripts). For each "Must NOT Have": search codebase for violations — reject with file:line if found. Check evidence files exist in .sisyphus/evidence/. Compare deliverables against plan. Verify the file manifest (from Task 2) matches actual `find analysis/ -type f` output.
  Output: `Must Have [N/N] | Must NOT Have [N/N] | Tasks [N/N] | VERDICT: APPROVE/REJECT`

- [ ] F2. **Syntax Validation of All Reorganized Scripts** — `unspecified-high`
  For every `.R` file in `analysis/` (recursively): run `Rscript -e "parse(file='<path>')"` to verify syntax. For every `.sh` file: run `bash -n <path>` for syntax check. Report any failures with file path and error message. Check that `analysis/utils.R` can be sourced without error.
  Output: `R Scripts [N/N parse OK] | Shell Scripts [N/N syntax OK] | utils.R [PASS/FAIL] | VERDICT`

- [ ] F3. **File Accounting — Every Original File Mapped** — `unspecified-high`
  Compare the original file inventory (38 R + 17 sh) against the file manifest from Task 2. Every original file must be accounted for: moved (new path), merged (into which script + line range), or deleted (with justification). Flag any unaccounted files. Check no unexpected new files exist beyond what the plan specifies.
  Output: `Original Files [55] | Moved [N] | Merged [N] | Deleted [N] | Unaccounted [N] | VERDICT`

- [ ] F4. **Path Integrity Check** — `deep`
  Grep all `.sh` files in `analysis/` for `Rscript` calls. For each, verify the referenced `.R` file exists at the specified path. Grep all `.R` files for `source(` calls. For each, verify the sourced file exists. Check AGENTS.md for any remaining references to old `analysis/` paths. Check for broken `setwd()` calls in merged scripts.
  Output: `Shell→R refs [N/N valid] | source() refs [N/N valid] | AGENTS.md refs [N/N valid] | setwd() [N/N safe] | VERDICT`

---

## Commit Strategy

Commits happen per-wave to create clean rollback points:

| Wave | Commit Message | Key Files |
|------|---------------|-----------|
| 1 | `chore(analysis): triage scripts and design final structure` | `.sisyphus/` manifest |
| 2 | `refactor(analysis): reorganize metaprograms, enrichment, clinical, cnv, plotting into themed subfolders` | `analysis/metaprograms/`, `analysis/enrichment/`, `analysis/clinical/`, `analysis/cnv/`, `analysis/cell_states/`, `analysis/plotting/` |
| 3 | `refactor(analysis): parameterize non-mali-nmf shells, populate utils.R` | `analysis/non_malignant_nmf/`, `analysis/utils.R` |
| 4 | `refactor(analysis): wire utils.R into scripts, fix all shell cross-references` | All `.R` and `.sh` files |
| 5 | `docs(agents): update AGENTS.md with reorganized analysis/ structure` | `AGENTS.md` |

---

## Success Criteria

### Verification Commands
```bash
# All R scripts parse without error
find analysis/ -name "*.R" -exec sh -c 'Rscript -e "parse(file=\"{}\")" 2>&1 || echo "FAIL: {}"' \;

# All shell scripts have valid syntax
find analysis/ -name "*.sh" -exec bash -n {} \; 2>&1

# No orphaned Rscript references in shell scripts
grep -rn "Rscript" analysis/*.sh analysis/**/*.sh | while read line; do
  script=$(echo "$line" | grep -oP 'Rscript\s+\K\S+')
  [ ! -f "$script" ] && echo "BROKEN: $line"
done

# Directory structure matches plan
find analysis/ -type d | sort

# File count check (should be fewer files than original 55 due to merges)
find analysis/ -type f | wc -l
```

### Final Checklist
- [ ] All "Must Have" present
- [ ] All "Must NOT Have" absent
- [ ] Every original file accounted for in manifest
- [ ] All .sh → .R references resolve
- [ ] utils.R sourced where needed
- [ ] AGENTS.md fully updated
