# Auto_instruction_latex.md

LuaLaTeX setup guide for compiling `.tex` files via VS Code LaTeX Workshop on the Imperial College HPC (SSH remote). Written for future agents/sessions to reproduce or debug the setup from scratch.

---

## Problem: `spawn lualatex ENOENT`

VS Code LaTeX Workshop tries to spawn `lualatex` as a subprocess. On HPC, `lualatex` is only available after loading HPC environment modules (`tools/prod` + `texlive/...`). The VS Code server process has a minimal PATH that does not include any module-loaded binaries, so the `spawn` call fails immediately with `ENOENT` (executable not found).

**Why `bash -c "module load ..."` in settings.json does not work:**
LaTeX Workshop substitutes its `%DOC%` / `%DOCFILE%` placeholders into the args array *before* spawning the process. When using `bash -c "... command"`, the placeholders end up inside the quoted string and are not expanded by LaTeX Workshop — they arrive as literal text. The extension effectively ignores the workspace `settings.json` tool definitions if they use `bash -c` and falls back to its built-in default tool (bare `lualatex`), which still fails with ENOENT.

---

## Solution Overview

Three files solve the problem completely:

| File | Location | Purpose |
|------|----------|---------|
| `lualatex` wrapper | `~/bin/lualatex` | Executable that sets up TexLive env and compiles via /tmp |
| `bibtex` wrapper | `~/bin/bibtex` | Same for bibtex |
| VS Code settings | `.vscode/settings.json` | Points LaTeX Workshop to the wrapper by absolute path |

---

## Step 1 — Create the lualatex wrapper at `~/bin/lualatex`

**Why this location:** `~/bin/` is already on VS Code server's PATH (`/rds/general/user/sg3723/home/.local/bin` is also on PATH as a fallback). The wrapper is kept at `~/bin/lualatex` so `settings.json` can reference it via a stable absolute path.

**Why skip module loading:** The `source /etc/profile.d/modules.sh && module load ...` approach is slow (~3-5s overhead per compile) and fragile inside subprocesses. Instead, the wrapper sets the required environment variables directly using absolute paths.

**Key finding:** `texlive/20230313-GCC-13.2.0` requires `tools/prod` as a prerequisite (NOT `tools/dev`). Attempting `tools/dev` silently fails. However the wrapper bypasses modules entirely using direct paths, making this moot.

**Why use /tmp for compilation:** LuaLaTeX reads hundreds of small `.sty`, `.fd`, and font files during compilation. On NFS (the `rds` filesystem), each file read is a network round-trip. Benchmarks showed 82 seconds on NFS vs 7 seconds on local `/tmp` — a 12x speedup. The wrapper copies the source directory to `/tmp`, compiles there, then copies the PDF and aux files back.

Create the file at `~/bin/lualatex`:

```bash
#!/bin/bash
# Fast lualatex wrapper for HPC
# Strategy: copy source to /tmp, compile locally (fast I/O), copy PDF back
TEXLIVE=/sw-eb/software/texlive/20230313-GCC-13.2.0
export LD_LIBRARY_PATH="/sw-eb/software/GCCcore/13.2.0/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$TEXLIVE/bin/x86_64-linux:$PATH"
export FONTCONFIG_PATH=/etc/fonts

# Redirect TeX caches to local /tmp
LOCALCACHE="/tmp/texcache_${USER}"
mkdir -p "$LOCALCACHE/texmf-var" "$LOCALCACHE/texmf-config" 2>/dev/null
export TEXMFVAR="$LOCALCACHE/texmf-var"
export TEXMFCONFIG="$LOCALCACHE/texmf-config"

# Parse args to find the .tex file (last argument)
ARGS=("$@")
TEXFILE="${ARGS[-1]}"

# If the tex file exists and is on NFS, compile in /tmp
if [[ -f "$TEXFILE" || -f "${TEXFILE}.tex" ]]; then
    # Resolve full path
    if [[ ! -f "$TEXFILE" && -f "${TEXFILE}.tex" ]]; then
        TEXFILE="${TEXFILE}.tex"
    fi
    SRCDIR="$(cd "$(dirname "$TEXFILE")" && pwd)"
    BASENAME="$(basename "$TEXFILE" .tex)"

    # Create local build dir
    BUILDDIR="/tmp/texbuild_${USER}_$$"
    mkdir -p "$BUILDDIR"

    # Copy entire source directory (tex, images, bib, etc.)
    cp -a "$SRCDIR"/* "$BUILDDIR"/ 2>/dev/null

    # Build remaining args (all except last)
    OTHERARGS=("${ARGS[@]:0:${#ARGS[@]}-1}")

    # Compile in /tmp
    cd "$BUILDDIR"
    "$TEXLIVE/bin/x86_64-linux/lualatex" "${OTHERARGS[@]}" "${BASENAME}.tex"
    EXITCODE=$?

    # Copy outputs back
    for ext in pdf synctex.gz aux log nav snm toc out; do
        [[ -f "$BUILDDIR/${BASENAME}.$ext" ]] && cp -f "$BUILDDIR/${BASENAME}.$ext" "$SRCDIR/"
    done

    # Cleanup
    rm -rf "$BUILDDIR"
    exit $EXITCODE
else
    # Fallback: run directly
    exec "$TEXLIVE/bin/x86_64-linux/lualatex" "$@"
fi
```

Run once to install:
```bash
mkdir -p ~/bin
# paste the above into ~/bin/lualatex, then:
chmod +x ~/bin/lualatex
```

Also copy to `~/.local/bin/` as a fallback (also on VS Code PATH):
```bash
cp ~/bin/lualatex ~/.local/bin/lualatex
chmod +x ~/.local/bin/lualatex
```

---

## Step 2 — Create the bibtex wrapper at `~/bin/bibtex`

```bash
#!/bin/bash
source /etc/profile.d/modules.sh
module purge 2>/dev/null
module load tools/prod 2>/dev/null
module load texlive/20230313-GCC-13.2.0 2>/dev/null
exec /sw-eb/software/texlive/20230313-GCC-13.2.0/bin/x86_64-linux/bibtex "$@"
```

```bash
chmod +x ~/bin/bibtex
cp ~/bin/bibtex ~/.local/bin/bibtex
chmod +x ~/.local/bin/bibtex
```

---

## Step 3 — Create `.vscode/settings.json`

Location: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/.vscode/settings.json`

**Key principle:** Set `command` to the **absolute path** of the wrapper. LaTeX Workshop spawns `command` directly as an executable — it does NOT go through a shell, so it must be a real file on disk. The `args` array uses standard LaTeX Workshop placeholders (`%DOC%`, `%DOCFILE%`) which the extension substitutes before spawning.

```jsonc
{
    "latex-workshop.latex.tools": [
        {
            "name": "lualatex",
            "command": "/rds/general/user/sg3723/home/bin/lualatex",
            "args": [
                "--synctex=1",
                "--interaction=nonstopmode",
                "%DOC%.tex"
            ],
            "env": {}
        },
        {
            "name": "bibtex",
            "command": "/rds/general/user/sg3723/home/bin/bibtex",
            "args": [
                "%DOCFILE%"
            ],
            "env": {}
        }
    ],

    "latex-workshop.latex.recipes": [
        {
            "name": "lualatex x2",
            "tools": ["lualatex", "lualatex"]
        },
        {
            "name": "lualatex",
            "tools": ["lualatex"]
        }
    ],

    "latex-workshop.latex.recipe.default": "first",
    "latex-workshop.latex.autoBuild.run": "onSave",
    "latex-workshop.view.pdf.viewer": "tab"
}
```

After editing settings.json: **Ctrl+Shift+P → "Developer: Reload Window"** to pick up changes.

---

## How to verify everything works

```bash
# 1. Check wrapper is executable and correct version
~/bin/lualatex --version    # should print: LuaHBTeX, Version 1.17.0 (TeX Live 2023)

# 2. Manual timed compile (should be ~7-10s)
cd /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Presentation
time ~/bin/lualatex --synctex=1 --interaction=nonstopmode methods_slides.tex

# 3. LaTeX Workshop: save the .tex file → auto-build fires, or:
#    Ctrl+Shift+P → "LaTeX Workshop: Build with recipe" → "lualatex x2"
```

---

## Performance notes

| Compile location | Time (22-page beamer .tex) | Notes |
|-----------------|---------------------------|-------|
| NFS directly | ~82 seconds | Hundreds of `.sty`/font reads over network |
| /tmp (wrapper) | ~7 seconds | All I/O local to login node |

**First compile in a new session:** may be 15-20s because the TeX font cache in `/tmp/texcache_sg3723/` was cleared when the login node rebooted. Subsequent compiles in the same session: ~7s.

**New `.tex` files:** work at full speed immediately — the wrapper copies any source directory to `/tmp` on the fly; no per-file setup needed.

**`/tmp` persistence:** `/tmp` is node-local and cleared on reboot. The build directory (`/tmp/texbuild_sg3723_*/`) is cleaned up by the wrapper after each compile. The font cache (`/tmp/texcache_sg3723/`) persists between compiles in the same session but is lost on reboot — this only causes a one-time slowdown.

---

## Troubleshooting

**`spawn lualatex ENOENT` still appears:**
- Reload VS Code window (Ctrl+Shift+P → "Developer: Reload Window") — settings.json changes require this
- Confirm the wrapper file exists: `ls -la ~/bin/lualatex`
- Confirm it is executable: `~/bin/lualatex --version`

**Compilation is slow again (~80s):**
- The `/tmp` font cache was cleared (node rebooted). First compile after reboot will be slow; subsequent ones will be fast automatically.
- If consistently slow, check the wrapper is the fast version: `head -3 ~/bin/lualatex` should show the `/tmp` strategy comment, not `module load`.

**`module load` errors in log:**
- The current wrapper does NOT use `module load` — it uses direct absolute paths. If you see module errors, the old wrapper is being used. Re-copy: `cp ~/.local/bin/lualatex ~/bin/lualatex`

**Missing image warnings (`File 'Figure_1.png' not found`):**
- These are non-fatal draft-mode placeholders — the PDF still generates. Place the actual image files in the same directory as the `.tex` file.

---

## HPC module reference (for manual terminal use only)

If you need to use lualatex manually in a terminal (not via VS Code), the correct module sequence is:

```bash
source /etc/profile.d/modules.sh
module purge
module load tools/prod                        # NOT tools/dev
module load texlive/20230313-GCC-13.2.0
export FONTCONFIG_PATH=/etc/fonts
lualatex --interaction=nonstopmode myfile.tex
```

The `texlive` module requires `tools/prod` as a prerequisite. Loading `tools/dev` first will fail with a silent module dependency error.
