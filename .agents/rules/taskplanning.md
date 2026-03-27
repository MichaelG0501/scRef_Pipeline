---
trigger: model_decision
description: Task Clarification Template
---

If my request for you to perform a task is unclear, you can ask to clarify any of the below points if highly relavent.

1. Target Scripts and Data Scope

Target Script(s): Are we modifying an existing script or creating a new one? (Please provide the exact /rds/general/... paths).

Data Inputs: Which exact .rds objects or metadata files are we reading from?

Methodology: Are we strictly sticking to Approach B and noreg (no cell cycle regression), or are there exceptions for this task?

2. Biological Logic and Definitions

MP/State Groupings: Are there any new collapsing rules for MPs into States? (e.g., Classic_Proliferative = c("MP5")).

Exclusions: Should we exclude, subset, or hide specific cell types, hybrid states, unresolved states, or cell cycle MPs (MP1, MP3, MP6, MP7) before running this analysis?

Thresholds: Are there specific numeric thresholds required for assignments (e.g., top 2 MPs, GSVA score minimums, gap < 0.3 for hybrids)?

3. Visualization and Plotting Strict Rules

Ordering: What is the exact sorting order? (Reminder: Default is strictly mp_tree_order or biological state order, never alphabetical unless explicitly requested).

Plot Layout: Do these need to be side-by-side, 2x2 grids, or separate pages in a single PDF?

Colors & Aesthetics: Are there specific color schemes required (e.g., matched patients needing distinct colors, with treated/untreated as light/dark shades)?

Legends & Labels: Do legends need to be off-plot? Do we need specific naming conventions (e.g., "MP number + description") or cell counts included in the legend?

Proportions: For stacked bar charts, should they be scaled strictly to 100% (no white space) based only on the subsets being plotted?

4. Outputs and Formatting

Directory: What is the exact output subfolder under ref_outs/ or analysis/?

Naming Convention: Is there a specific prefix (e.g., task number) or file name required for the output PDFs/PNGs?