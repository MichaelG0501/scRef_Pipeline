# Progress Update Workflow — Repeatable Instructions

**Purpose:** Generate supervisor-facing slides and progress markdown most Mondays and Thursdays when needed.  
**Location:** `updates/`  
**Discovery method:** `git status` — untracked files = new work since last commit.

---

## Quick Reference (Copy-Paste for New Session)

> Create my progress update slides and markdown. Follow the workflow in `updates/weekly_update_workflow.md`. Today is [DATE].

---

## Step 1: Discover This Cycle's Files

Untracked files in git represent all new work since the last commit:

```bash
WD=/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
cd "$WD"
git status --short
```

- `??` entries = new files (scripts, wrappers) → this cycle's work
- `M` entries = modified files → check `git diff <file>` for what changed

**For each untracked `.R` / `.sh` script**, read the header comment block (lines 1–20) to find:
- `# Output:` lines listing associated figures/data files in `ref_outs/`
- Brief description of what the script does

These output files are gitignored (`.png`, `.pdf`, `.rds` all in `ref_outs/`) but are the **associated results** of the untracked scripts.

```bash
# Example: extract output lines from a script header
head -20 analysis/cell_states/Auto_states_cluster.R | grep -i "output\|#.*ref_outs"
```

Also check git log for context on what was committed in the previous cycle:

```bash
git log --oneline -3
```

---

## Step 2: Select Key Figures

Pick **informative and critical figures** for a concise slide deck but do not miss information. Prioritise:

1. **The main result** (e.g., ComplexHeatmap at the chosen nMP resolution)
2. **Cell state / classification result** (e.g., stacked bar chart, state composition)
3. **Validation figure** (e.g., 3CA heatmap by state, external reference comparison)
4. **1–2 enrichment heatmaps** (Hallmark is usually the most informative)
5. **Correlation / scoring heatmap** (e.g., UCell, Spearman)
6. **1–2 TME figures** (if TME work was done)
7. **Survival / clinical** (if applicable)

**Known figure locations:**

| Type | Path pattern |
|------|-------------|
| NMF heatmaps | `ref_outs/Metaprogrammes_Results/*_complexheatmap*.png` |
| Enrichment | `ref_outs/enrich_*.png` |
| Cell states | `ref_outs/Auto_*_states_*.pdf` |
| Correlation | `ref_outs/Auto_MP_correlation_*.pdf` |
| UCell heatmap | `ref_outs/*UCellHeatmap*.pdf` |
| TME NMF | `ref_outs/nmf_summaries/metaprograms_heatmap_*.png` |
| TME enrichment | `ref_outs/nmf_{celltype}/enrich_*.png` |
| Annotation PDFs | `ref_outs/New_nMP*_anno.pdf` |
| Survival | `ref_outs/*survival*.pdf` or `ref_outs/*KM*.pdf` |
| 3CA validation | `ref_outs/3CA_*.pdf` |
| Presentation | `ref_outs/Presentation/*.pdf` |
| **PDO NMF / annotation** | `/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/` |

**Also check the PDO pipeline directory for new outputs each cycle:**

```
/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/
  PDOs_outs/          — NMF results, enrichment PDFs, correlation heatmaps
  analysis/           — PDO analysis R scripts
```

Key PDO figure types to include in updates:
- `PDOs_outs/optimal_nMP_metrics.png` — nMP selection plot
- `PDOs_outs/metaprograms_heatmap.png` — NMF factor heatmap
- `PDOs_outs/PDOs_enrichment_annotation.pdf` — MP annotation (all databases)
- `PDOs_outs/mp_jaccard_self_pdo.pdf` — Jaccard self-similarity
- `PDOs_outs/mp_correlation_pdo_meta_complexheatmap.pdf` — per-sample Spearman correlation
- `PDOs_outs/unresolved_states/*` — PDO unresolved cell pan-cancer relabeling results
- `PDOs_outs/Auto_PDO_mp_correlation_*.pdf` — PDO versus scRef MPs correlation and cross-validation plots
- `PDOs_outs/Auto_PDO_matched_*.pdf` — PDO Matched treated/untreated pairs temporal state analysis
Copy selected figures into this cycle's subfolder:

```bash
mkdir -p updates/DDmon/figures
cp <chosen_figures> updates/DDmon/figures/
```

Use short names (e.g., `cluster_heatmap.pdf`, `comparison_umap.pdf`).

### Machine-Readable Summaries

All future `Auto_` scripts **must** also save a small, machine-readable summary alongside their plot outputs. This allows AI agents to read results directly on the login node without loading heavy `.rds` files.

- Format: `.csv` or `.txt` (< 100 KB)
- Content: key metrics, counts, proportions — whatever summarises the visual output
- Location: saved directly into `updates/DDmon/summaries/` at script runtime
- Example: `Auto_states_comparison_summary.csv` contains ARI, NMI, bootstrap stats
- Scripts should accept an optional output directory argument, or default to the current cycle's update folder

---

## Step 3: Write the .tex Slides

Slides are **results-focused** — straight to figures, with an optional brief overview slide when comparing methods. Supervisor does not need methods detail.

### Rules
- **No title slide**. No `\titlepage`. Jump straight into result figures or a brief methods overview if needed.
- **Titles**: Should plainly state what the graph is (e.g. "TCGA Survival Analysis Comparison - States - Continuous Cox"). Never include things like "PartA", "page 1", "example", etc.
- **Footnotes & Spacing**: 
  - If self-explanatory (e.g. "State proportions"): Omit the footnote completely. Use figure height `0.92\textheight` or `0.94\textheight` and no `\vspace`.
  - If requires context (e.g. noting a biological finding or methodological detail): Provide a short, direct footnote. No fancy AI words. Use figure height `0.88\textheight` and add `\vspace{-3mm}` before the footnote.
  - For side-by-side or mutli-panel comparisons where a footnote is used, use `\vspace{1mm}` or `\vspace{2mm}` to ensure enough separation from the columns.
- **Cropping**: When placing multiple figures on one slide, crop white space from the source PDFs/PNGs so plots fill the available space for presentation.
- **Layouts**: 
  - **Single figure**: full width/height.
  - **Side-by-side**: Use `\begin{columns}` with `0.5\textwidth` each for direct comparisons (e.g., Correlation vs Jaccard, PDO vs scAtlas, reg vs noreg). Use `0.8\textheight` for figure height inside columns.
  - **Many panels (e.g., 4 plots)**: Use a 2x2 grid via `\begin{minipage}{0.48\textwidth}` macros with `0.48\textheight` or `0.5\textheight` per plot.
- **Text slides** Use only when a brief methods summary is essential (e.g. describing 'Five approaches to pairwise distance' or 'reg vs noreg'). Keep to bullet points, small font.
- **No next-steps slide**. Those go in the `.md` document only.
- **Never convert to .pptx until user confirms the content.**

### LaTeX Template

```latex
\documentclass[aspectratio=169,10pt]{beamer}

\usepackage{graphicx}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Times New Roman}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{array}

% Black theme
\usetheme{default}
\usecolortheme{default}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\setbeamertemplate{frametitle}[default][center]
\setbeamerfont{frametitle}{series=\bfseries,size=\large}
\setbeamerfont{title}{series=\bfseries,size=\large}
\setbeamersize{text margin left=10mm,text margin right=10mm}
\setbeamercolor{frametitle}{fg=black}
\setbeamercolor{title}{fg=black}
\setbeamercolor{structure}{fg=black}
\setbeamercolor{normal text}{fg=black}
\setbeamercolor{footline}{fg=black!60}
\setbeamertemplate{itemize items}{\textbullet}
\setbeamertemplate{itemize item}{\textbullet}
\setbeamertemplate{itemize subitem}{\textbullet}

\begin{document}

% Full-width figure slide (default pattern):
\begin{frame}{Descriptive Title}
\begin{center}
    \includegraphics[width=\textwidth,height=0.88\textheight,keepaspectratio]{figures/FILENAME}
\end{center}
\vspace{-3mm}
{\footnotesize One-line key finding or context.}
\end{frame}

% Side-by-side comparison (when needed):
\begin{frame}{Comparison Title}
\begin{columns}
\begin{column}{0.5\textwidth}
    \centering
    {\small \textbf{Left label}}\\[2mm]
    \includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{figures/LEFT}
\end{column}
\begin{column}{0.5\textwidth}
    \centering
    {\small \textbf{Right label}}\\[2mm]
    \includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{figures/RIGHT}
\end{column}
\end{columns}
\end{frame}

% Metrics table (when needed):
\begin{frame}{Metrics Title}
\vspace{5mm}
\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{>{\raggedright}p{5.5cm} c c}
\toprule
\textbf{Metric} & \textbf{Col A} & \textbf{Col B} \\
\midrule
Row 1 & val & val \\
Row 2 & val & val \\
\bottomrule
\end{tabular}
\end{center}
\end{frame}

% Brief methods bullet slide (only when essential):
\begin{frame}{Methods Note}
\begin{itemize}
    \setlength{\itemsep}{8pt}
    \item \textbf{Point 1}
    \begin{itemize}
        \item Sub-detail
    \end{itemize}
\end{itemize}
\end{frame}

\end{document}
```

### Known LaTeX Gotchas

- **DO NOT** use magic comments like `% !TEX program = lualatex` at the top of the file. This bypasses our custom wrapper and causes `ENOENT` errors.
- **DO NOT** use `\usepackage{enumitem}` — conflicts with beamer.
- **DO NOT** use `\setcounter{enumi}` — use `\item[N.]` for manual numbering.
- Figures in `figures/` subfolder are automatically copied to `/tmp` by the lualatex wrapper.
- PDF and PNG figures both work. PNG for heatmaps, PDF for R plot outputs.
- Use `10pt` document class option (not `12pt`) for cleaner density.
- Escape underscores in frame titles: `\_` or avoid them.

---

## Step 4: Compile LaTeX

Compile within VS Code by saving (Ctrl+S). If you need to compile manually in the terminal:

```bash
cd updates/DDmon
~/bin/lualatex --synctex=1 --interaction=nonstopmode update_DDmon.tex
```

The `~/bin/lualatex` wrapper compiles via `/tmp` for speed (~7s vs 82s on NFS).

---

## Step 5: Convert to PPTX

Convert the compiled PDF to high-resolution PNG slides, then assemble into PPTX:

```bash
module load tools/prod
module load texlive/20230313-GCC-13.2.0
module load poppler/24.04.0-GCC-13.2.0
cd updates/DDmon
mkdir -p pptx_slides_png
pdftoppm -png -r 600 update_DDmon.pdf pptx_slides_png/slide
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/py39
python ../png2pptx.py DDmon
```

The `png2pptx.py` script lives in the `updates/` root and accepts the update folder name as argument. It reads PNGs from `updates/DDmon/pptx_slides_png/` and writes `updates/DDmon/update_DDmon.pptx`.

---

## Step 6: Upload to Google Drive

Only upload the final PPTX:

```bash
module load rclone
rclone copyto updates/DDmon/update_DDmon.pptx gdrive:IMPERIAL/update_DDmon.pptx --progress
```

---

## Step 7: Write the .md Progress Document

Same git-based discovery as Step 1. Date range = "since last commit".

### Structure

```markdown
# Progress Update — DD Month YYYY

**Project:** ...

---

## Context
One paragraph: what the cycle's focus was and any key decisions.

## Summary of Work
Bullet points — results-oriented, not process-oriented.

## Key Outputs
Table of new/modified files with brief descriptions.
For modified scripts: describe ONLY what was added/changed, not the whole file.

## Decisions Made
Numbered list of key choices and rationale.

## Next Steps
Checkbox list.
```

---

## What NOT To Do

| Anti-pattern | Why it's slow | Do this instead |
|---|---|---|
| Multiple `find` passes | Redundant I/O; git is the source of truth | `git status --short` |
| Writing too much text on slides | Supervisor needs results, not methods | One sentence per key point, large figures |
| Inspect every figure with `look_at` | Slow, fails on PDFs | Refer to known figure locations table; you know what each output type shows |
| Full file reads of scripts | Wastes context reading unchanged code | Read only header block (lines 1–20) of untracked scripts |
| Write text-heavy slides then rewrite | Double work | Start from results-focused template |
| Query NotebookLM every cycle | ~60s browser automation | Optional — only query when project scope changes or when field context is needed |

---

## File Naming Convention

Each update cycle gets its own subfolder under `updates/`:

```
updates/
├── DDmon/                          # e.g., 02mar/
│   ├── figures/                    # copied figure files
│   ├── summaries/                  # machine-readable CSVs/TXTs saved directly by scripts
│   ├── pptx_slides_png/            # intermediate PNGs from pdftoppm (600 dpi)
│   ├── update_DDmon.tex            # LaTeX slides source
│   ├── update_DDmon.pdf            # compiled PDF (intermediate)
│   ├── update_DDmon.pptx           # final PPTX (uploaded to GDrive)
│   ├── update_DDmon.md             # progress markdown
│   └── update_DDmon_plan.md        # plan/prompt file (pre-execution)
├── png2pptx.py                     # shared PNG→PPTX converter
└── weekly_update_workflow.md       # this file (do not rename or move)
```

---

## Update after each update cycle, in case style or focus change

*Last updated: 12 March 2026*
