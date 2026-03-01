# Weekly Update Workflow — Repeatable Instructions

**Purpose:** Generate lab meeting slides and progress markdown each week.  
**Time target:** ~5 minutes total (vs 20+ min in the initial run).  
**Location:** `weekly_updates/`

---

## Quick Reference (Copy-Paste for New Session)

> Create my weekly update slides and markdown. Follow the workflow in `weekly_updates/weekly_update_workflow.md`. The date range is [START_DATE] to [END_DATE].

---

## Step 1: Discover This Week's Files (~30s)

Single command — finds everything modified since the start date:

```bash
WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline
find "$WD/ref_outs" "$WD/analysis" -newermt "YYYY-MM-DD" -type f \
  \( -name "*.png" -o -name "*.pdf" -o -name "*.rds" -o -name "*.R" -o -name "*.sh" \) \
  | sort | tee /tmp/weekly_files.txt
```

Also check git for what scripts changed and how:

```bash
cd "$WD"
git log --since="YYYY-MM-DD" --oneline --name-status
```

**Important:** For `.R` / `.sh` files marked as "Modified" in git, only the **added chunks** are this week's work — not the entire file. Use `git diff` to see what was actually added.

---

## Step 2: Select Key Figures (~1 min)

Pick **5–10 figures max** for a concise slide deck. Prioritise:

1. **The main result** (e.g., ComplexHeatmap at the chosen nMP resolution)
2. **Cell state / classification result** (e.g., stacked bar chart, state composition)
3. **Validation figure** (e.g., 3CA heatmap by state, external reference comparison)
4. **1–2 enrichment heatmaps** (Hallmark is usually the most informative)
5. **Correlation / scoring heatmap** (e.g., UCell, Spearman)
6. **1–2 TME figures** (if TME work was done that week)
7. **Survival / clinical** (if applicable)

**Known figure locations:**

| Type | Path pattern |
|------|-------------|
| NMF heatmaps | `ref_outs/Metaprogrammes_Results/*_complexheatmap*.png` |
| Enrichment | `ref_outs/enrich_*.png` |
| Cell states | `ref_outs/MP_state_*.pdf` |
| Correlation | `ref_outs/Auto_MP_correlation_*.pdf` |
| UCell heatmap | `ref_outs/*UCellHeatmap*.pdf` |
| TME NMF | `ref_outs/nmf_summaries/metaprograms_heatmap_*.png` |
| TME enrichment | `ref_outs/nmf_{celltype}/enrich_*.png` |
| Annotation PDFs | `ref_outs/New_nMP*_anno.pdf` |
| Survival | `ref_outs/*survival*.pdf` or `ref_outs/*KM*.pdf` |
| 3CA validation | `ref_outs/3CA_*.pdf` |
| Presentation | `ref_outs/Presentation/*.pdf` |

Copy selected figures:

```bash
mkdir -p weekly_updates/figures
cp <chosen_figures> weekly_updates/figures/
```

Use short names (e.g., `nmp19_heatmap.png`, `cell_states.pdf`, `enrich_hallmark.png`).

---

## Step 3: Write the .tex Slides (~2 min)

Use the template below. The slide deck should be **results-focused** — mostly full-page figures with minimal text.

### Template Structure

```
1. Title slide
2. Overview slide (5–6 bullet points of what was done)
3–N. One figure per slide (full width, large, clearly visible)
N+1. Next steps slide (4–5 bullet points)
```

### LaTeX Template

```latex
% !TEX program = lualatex
\documentclass[aspectratio=169,12pt]{beamer}

\usepackage{graphicx}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Times New Roman}
\usepackage{xcolor}

% Black theme
\usetheme{default}
\usecolortheme{default}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\setbeamertemplate{frametitle}[default][center]
\setbeamerfont{frametitle}{series=\bfseries,size=\Large}
\setbeamerfont{title}{series=\bfseries,size=\LARGE}
\setbeamersize{text margin left=10mm,text margin right=10mm}
\setbeamercolor{frametitle}{fg=black}
\setbeamercolor{title}{fg=black}
\setbeamercolor{structure}{fg=black}
\setbeamercolor{normal text}{fg=black}
\setbeamercolor{footline}{fg=black!60}
\setbeamertemplate{itemize items}{\textbullet}
\setbeamertemplate{itemize item}{\textbullet}
\setbeamertemplate{itemize subitem}{\textbullet}

\title{Weekly Update}
\subtitle{OAC scATLAS --- SUBTITLE\\Week of DD--DD Month YYYY}
\author{Shengbo Gong}
\date{DD Month YYYY}

\begin{document}

\begin{frame}
\titlepage
\end{frame}

\begin{frame}{This Week: Overview}
\large
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item Point 1
    \item Point 2
\end{itemize}
\end{frame}

% One figure per slide — use this pattern:
\begin{frame}{Slide Title}
\begin{center}
    \includegraphics[width=0.92\textwidth,height=0.82\textheight,keepaspectratio]{figures/FILENAME}
\end{center}
\end{frame}

% For figure + side text:
\begin{frame}{Slide Title}
\begin{columns}[c]
\begin{column}{0.55\textwidth}
    \begin{center}
        \includegraphics[width=\textwidth,height=0.78\textheight,keepaspectratio]{figures/FILENAME}
    \end{center}
\end{column}
\begin{column}{0.42\textwidth}
    \large
    \begin{itemize}
        \setlength{\itemsep}{8pt}
        \item Key point
    \end{itemize}
\end{column}
\end{columns}
\end{frame}

\begin{frame}{Next Steps}
\large
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item Next step 1
    \item Next step 2
\end{itemize}
\end{frame}

\end{document}
```

### Known LaTeX Gotchas

- **DO NOT** use `\usepackage{enumitem}` — conflicts with beamer.
- **DO NOT** use `\setcounter{enumi}` — use `\item[N.]` for manual numbering.
- Figures in `figures/` subfolder are automatically copied to `/tmp` by the lualatex wrapper.
- PDF and PNG figures both work. PNG for heatmaps, PDF for R plot outputs.

---

## Step 4: Compile (~10s)

```bash
module purge && module load tools/dev
~/bin/lualatex -interaction=nonstopmode weekly_updates/weekly_update_DDMMM.tex
# Second pass for page numbers:
~/bin/lualatex -interaction=nonstopmode weekly_updates/weekly_update_DDMMM.tex
```

The `~/bin/lualatex` wrapper compiles via `/tmp` for speed (~7s vs 82s on NFS).

---

## Step 5: Write the .md Progress Document (~2 min)

### Structure

```markdown
# Weekly Progress Update — DD–DD Month YYYY

**Project:** ...
**Author:** Shengbo Gong

---

## Context
One paragraph: what the week's focus was and any key decisions.

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

## What NOT To Do (Lessons from Week 1)

| Anti-pattern | Why it's slow | Do this instead |
|---|---|---|
| Query NotebookLM every week | ~60s browser automation | Only query when project scope changes; refer to previous `.md` for context |
| Launch explore/librarian agents | Unnecessary for known codebase | Use direct `find` + `git log` |
| Multiple `find` passes (PNG, PDF, RDS separately) | Redundant I/O | Single `find` with `-o` combining all extensions |
| Inspect every figure with `look_at` | Slow, fails on PDFs | Refer to known figure locations table above; you already know what each output type shows |
| Write text-heavy slides then rewrite | Double work | Start from results-focused template |
| Full file reads of scripts to understand changes | Wastes context reading unchanged code | Use `git diff --since` or check 20-hash comment blocks only |

---

## File Naming Convention

- Slides: `weekly_update_DDmon.tex` (e.g., `weekly_update_22feb.tex`)
- Markdown: `weekly_update_DDmon.md`
- Figures: `figures/` subfolder with short descriptive names
- This workflow doc: `weekly_update_workflow.md` (do not rename)

---

*Last updated: 1 March 2026*
