#!/usr/bin/env python3
####################
# scatlas_velocity_scvelo_visualise.py
# Status: active
# Script: analysis/trajectory/scatlas_velocity_scvelo_visualise.py
# Methodology: analysis/methodology/trajectory/scatlas_velocity_methodology.md
# Inputs: live velocity cell metadata/sample manifest plus per-sample live or
#   ephemeral velocyto loom and versioned scVelo H5AD cache.
# Outputs: ref_outs/Auto_velocity_scATLAS/h5ad/, tables/ state/MP node and edge
#   CSVs, and figures/ per-sample/state-direction PDFs; persistent copies are live.
# Cache/replot: VELOCITY_CACHE_VERSION guards H5AD reuse; cached models may be
#   staged ephemeral but all replot tables/H5AD/final figures are live.
# Run: qsub analysis/trajectory/scatlas_velocity_run_scvelo_visualisation.sh
# Conda env: velocity
#
# Per-sample scVelo analysis with LOCAL per-sample UMAP recomputation
# and scATLAS state-transition summaries.
#
# Key changes from v1:
#   - UMAP is recomputed per sample (not injected from global coordinates)
#   - Uses centred refined noreg states only (no 3CA relabeling)
#   - MP-level velocity direction for Basal and SMG states
#   - Cache versioning (VELOCITY_CACHE_VERSION = 2)
####################

from pathlib import Path
from typing import List, Optional, Set, Tuple

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PathCollection
from matplotlib.patches import FancyArrowPatch
from matplotlib import patheffects as pe
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv


VELOCITY_CACHE_VERSION = 2

WD_EPH = Path("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline")
OUT_EPH = WD_EPH / "ref_outs" / "Auto_velocity_scATLAS"
WD_LIVE = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
OUT_LIVE = WD_LIVE / "ref_outs" / "Auto_velocity_scATLAS"

MAJOR_STATES = [
    "Classic proliferation",
    "Basal to intestinal metaplasia",
    "SMG to intestinal metaplasia",
    "Stress adaptive",
    "Cancer-cell immune mimicry",
]
STATE_COLORS = {
    "Classic proliferation": "#E41A1C",
    "Basal to intestinal metaplasia": "#4DAF4A",
    "SMG to intestinal metaplasia": "#FF7F00",
    "Stress adaptive": "#984EA3",
    "Cancer-cell immune mimicry": "#377EB8",
    "Unresolved": "grey",
    "Hybrid": "black",
    "Unassigned": "#A6A6A6",
}

BASAL_MPS = ["MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"]
SMG_MPS = ["MP8+", "MP8b", "MP16", "MP18b", "MP17"]

BASAL_MP_COLORS = {
    "MP14": "#4DAF4A", "MP3+": "#8DA0CB", "MP6+": "#66C2A5",
    "MP11+": "#FC8D62", "MP9+": "#A6D854", "MP10+": "#E78AC3",
    "Hybrid": "black", "Unassigned": "#A6A6A6",
}
SMG_MP_COLORS = {
    "MP8+": "#FF7F00", "MP8b": "#A65628", "MP16": "#F781BF",
    "MP18b": "#999999", "MP17": "#66C2A5",
    "Hybrid": "black", "Unassigned": "#A6A6A6",
}

MP_DESC_MAP = {
    "MP1": "G2/M cell cycle",
    "MP5": "G1/S cell cycle",
    "MP13+": "replication-stress-associated cell cycling",
    "MP2+": "MYC driven biosynthesis",
    "MP14": "Squamoid/basal transition",
    "MP3+": "Basal-columnar invasive epithelium",
    "MP6+": "Stress-reactive columnar epithelium",
    "MP11+": "Epithelial antiviral interferon response",
    "MP9+": "Metabolic columnar epithelium",
    "MP10+": "Intestinal metaplasia",
    "MP8+": "Glandular intestinal metaplasia",
    "MP8b": "Metabolic intestinal metaplasia",
    "MP16": "Mucous-secretory glandular epithelium",
    "MP18b": "Mucous-secretory differentiation",
    "MP17": "Immune-interactive glandular progenitor",
    "MP12": "Hypoxic inflammatory adaptive plasticity",
    "MP15": "T/NK-like cancer-cell immune mimicry",
}

LAYOUT = {
    "Classic proliferation": np.array([-1.18, 0.76]),
    "Basal to intestinal metaplasia": np.array([1.18, 0.76]),
    "SMG to intestinal metaplasia": np.array([1.18, -0.76]),
    "Stress adaptive": np.array([-1.18, -0.76]),
    "Cancer-cell immune mimicry": np.array([0.0, -1.34]),
}
EDGE_THRESHOLD = 0.35
CORE_PDF = "Auto_scatlas_velocity_per_sample_visualisations.pdf"
EXTENDED_PDF = "Auto_scatlas_velocity_per_sample_visualisations_extended.pdf"
DIRECTION_PDF = "Auto_scatlas_velocity_state_directions.pdf"
NODE_LABEL_FONTSIZE = 10
ALIGNMENT_LABEL_FONTSIZE = 10
ARROW_NODE_SHRINK = 54
ARROW_RAD = 0.24


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def comma_int(value: int | float) -> str:
    return f"{int(value):,}"


def total_cells_from_nodes(nodes: pd.DataFrame, sample_set: Optional[Set[str]] = None) -> int:
    if not len(nodes):
        return 0
    use_nodes = nodes
    if sample_set is not None:
        use_nodes = nodes[nodes["sample"].isin(sample_set)]
    if not len(use_nodes):
        return 0
    if "total_cells" in use_nodes.columns:
        return int(use_nodes[["sample", "total_cells"]].drop_duplicates()["total_cells"].sum())
    return int(use_nodes.groupby("sample")["cells"].sum().sum())


def title_with_n(title: str, nodes: pd.DataFrame, sample_set: Optional[Set[str]] = None) -> str:
    return f"{title} (n = {comma_int(total_cells_from_nodes(nodes, sample_set))} cells)"


def format_state_label(state: str) -> str:
    labels = {
        "Classic proliferation": "Classic\nproliferation",
        "Basal to intestinal metaplasia": "Basal to\nintestinal meta",
        "SMG to intestinal metaplasia": "SMG to\nintestinal meta",
        "Stress adaptive": "Stress\nadaptive",
        "Cancer-cell immune mimicry": "Cancer-cell\nimmune mimicry",
    }
    return labels.get(state, state)


def format_mp_label(mp: str) -> str:
    desc = MP_DESC_MAP.get(mp, mp)
    short = desc[:20] + "…" if len(desc) > 20 else desc
    return f"{mp}\n{short}"


def edge_curve_side(source: str, target: str) -> float:
    return 1.0 if MAJOR_STATES.index(source) < MAJOR_STATES.index(target) else -1.0


def extract_barcode(cell_id: str) -> str:
    bc = str(cell_id).split(":")[-1]
    if bc.endswith("x"):
        bc = bc[:-1]
    return bc


# ---------------------------------------------------------------------------
# Sample loading & velocity computation
# ---------------------------------------------------------------------------

def sample_has_velocity_input(sample: str) -> bool:
    h5ad_live = OUT_LIVE / "h5ad" / f"Auto_scvelo_{sample}.h5ad"
    h5ad_eph = OUT_EPH / "h5ad" / f"Auto_scvelo_{sample}.h5ad"
    loom_dir_live = OUT_LIVE / "looms" / sample
    loom_dir_eph = OUT_EPH / "looms" / sample
    return (
        h5ad_live.exists()
        or h5ad_eph.exists()
        or len(list(loom_dir_live.glob("*.loom"))) == 1
        or len(list(loom_dir_eph.glob("*.loom"))) == 1
    )


def load_loom_and_match(sample: str, meta: pd.DataFrame) -> ad.AnnData:
    """Load loom file, match barcodes to metadata, return raw AnnData."""
    loom_dir_live = OUT_LIVE / "looms" / sample
    loom_dir_eph = OUT_EPH / "looms" / sample
    looms = sorted(loom_dir_live.glob("*.loom"))
    if len(looms) != 1:
        looms = sorted(loom_dir_eph.glob("*.loom"))
    if len(looms) != 1:
        raise FileNotFoundError(f"Expected exactly one loom in {loom_dir_live} or {loom_dir_eph}, found {len(looms)}")

    sample_meta = meta[meta["sample"] == sample].copy()
    raw_to_cell = dict(zip(sample_meta["raw_barcode"], sample_meta["cell_id"]))
    raw_set = set(raw_to_cell)

    adata = sc.read_loom(str(looms[0]), sparse=True)
    adata.var_names_make_unique()

    raw_barcodes = []
    cell_ids = []
    keep = []
    for obs_name in adata.obs_names:
        bc = extract_barcode(obs_name)
        if bc not in raw_set and not bc.endswith("-1") and f"{bc}-1" in raw_set:
            bc = f"{bc}-1"
        raw_barcodes.append(bc)
        keep.append(bc in raw_set)
        cell_ids.append(raw_to_cell.get(bc, ""))

    adata.obs["raw_barcode"] = raw_barcodes
    adata = adata[np.array(keep), :].copy()
    adata.obs_names = [x for x, k in zip(cell_ids, keep) if k]
    adata.obs_names_make_unique()

    if adata.n_obs < 30:
        raise ValueError(f"{sample} has only {adata.n_obs} matched cells after loom/barcode join.")

    sample_meta = sample_meta.set_index("cell_id").loc[adata.obs_names]
    for col in sample_meta.columns:
        adata.obs[col] = sample_meta[col].values

    # State assignments from centred refined noreg
    adata.obs["state_noreg"] = adata.obs["state_noreg"].fillna("Unassigned").astype(str)
    # Direction states: only 5 primary states used for direction scoring
    adata.obs["state_direction"] = adata.obs["state_noreg"].apply(
        lambda x: x if x in MAJOR_STATES else "Unassigned"
    )
    # state_final = state_noreg (no 3CA relabeling)
    adata.obs["state_final"] = adata.obs["state_noreg"].astype(str)

    if "dominant_basal_mp" in adata.obs.columns:
        adata.obs["dominant_basal_mp"] = adata.obs["dominant_basal_mp"].fillna("Unassigned").astype(str)
    if "dominant_smg_mp" in adata.obs.columns:
        adata.obs["dominant_smg_mp"] = adata.obs["dominant_smg_mp"].fillna("Unassigned").astype(str)
    adata.obs["dataset"] = adata.obs.get("study", "unknown").astype(str)

    return adata


def sanitise_velocity_fields(adata: ad.AnnData) -> None:
    for key in ["X_umap", "velocity_umap"]:
        if key in adata.obsm:
            adata.obsm[key] = np.nan_to_num(adata.obsm[key], nan=0.0, posinf=0.0, neginf=0.0)
    for col in ["velocity_length", "velocity_confidence"]:
        if col in adata.obs:
            vals = pd.to_numeric(adata.obs[col], errors="coerce").to_numpy(dtype=float)
            adata.obs[col] = np.nan_to_num(vals, nan=0.0, posinf=0.0, neginf=0.0)


def run_velocity_local_umap(adata: ad.AnnData) -> ad.AnnData:
    """Run the full scVelo pipeline with per-sample UMAP recomputation."""
    scv.settings.verbosity = 2
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    if adata.n_vars > 3000:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True, flavor="seurat")
    n_pcs = max(2, min(30, adata.n_obs - 1, adata.n_vars - 1))
    n_neighbors = max(5, min(30, adata.n_obs - 1))
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Recompute UMAP locally for this sample
    sc.tl.umap(adata)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="umap")
    scv.tl.velocity_confidence(adata)
    sanitise_velocity_fields(adata)
    return adata


def _cache_is_v2(adata: ad.AnnData) -> bool:
    """Check if cached h5ad was produced with v2 (local UMAP) pipeline."""
    if "velocity_cache_version" not in adata.uns:
        return False
    return int(adata.uns["velocity_cache_version"]) >= VELOCITY_CACHE_VERSION


def load_or_run_velocity(sample: str, meta: pd.DataFrame, force_recompute: bool = False) -> ad.AnnData:
    h5ad_live = OUT_LIVE / "h5ad" / f"Auto_scvelo_{sample}.h5ad"
    h5ad_eph = OUT_EPH / "h5ad" / f"Auto_scvelo_{sample}.h5ad"

    if not force_recompute:
        for h5ad_path in [h5ad_live, h5ad_eph]:
            if h5ad_path.exists():
                adata = ad.read_h5ad(h5ad_path)
                required_obsm = {"X_umap", "velocity_umap"}
                required_obs = {"state_final", "state_direction", "sample", "dataset"}
                if (
                    required_obsm.issubset(adata.obsm.keys())
                    and required_obs.issubset(adata.obs.columns)
                    and _cache_is_v2(adata)
                ):
                    # Refresh MP assignments from metadata
                    meta_indexed = meta.set_index("cell_id")
                    for col in ["dominant_basal_mp", "dominant_smg_mp"]:
                        if col in meta_indexed.columns:
                            adata.obs[col] = meta_indexed[col].reindex(adata.obs_names).fillna("Unassigned").astype(str)
                    sanitise_velocity_fields(adata)
                    if not h5ad_live.exists():
                        adata.write(h5ad_live, compression="gzip")
                    return adata

    # Load from loom and run full pipeline
    adata = load_loom_and_match(sample, meta)
    adata = run_velocity_local_umap(adata)
    adata.uns["velocity_cache_version"] = VELOCITY_CACHE_VERSION

    # Save to both ephemeral (cache) and live (persistent)
    (OUT_EPH / "h5ad").mkdir(parents=True, exist_ok=True)
    (OUT_LIVE / "h5ad").mkdir(parents=True, exist_ok=True)
    adata.write(h5ad_eph, compression="gzip")
    adata.write(h5ad_live, compression="gzip")
    return adata


# ---------------------------------------------------------------------------
# State-level direction tables
# ---------------------------------------------------------------------------

def state_direction_tables(adata: ad.AnnData) -> Tuple[List[dict], List[dict]]:
    xy = adata.obsm["X_umap"]
    vv = adata.obsm["velocity_umap"]
    states = adata.obs["state_direction"].astype(str).to_numpy()
    sample = str(adata.obs["sample"].iloc[0])
    dataset = str(adata.obs["dataset"].iloc[0])

    total_direction_cells = int(len(states))
    total_major = int(np.isin(states, MAJOR_STATES).sum())
    node_rows = []
    centers = {}
    mean_velocities = {}
    for state in MAJOR_STATES:
        mask = states == state
        cells = int(mask.sum())
        pct_total = 100 * cells / total_direction_cells if total_direction_cells else 0.0
        node_rows.append({
            "sample": sample,
            "dataset": dataset,
            "source_group": dataset,
            "state": state,
            "cells": cells,
            "total_cells": total_direction_cells,
            "major_cells": total_major,
            "pct_major": pct_total,
            "pct_total_direction_states": pct_total,
            "pct_of_major_states": 100 * cells / total_major if total_major else 0.0,
        })
        if cells >= 5:
            centers[state] = xy[mask].mean(axis=0)
            mean_velocities[state] = vv[mask].mean(axis=0)

    edge_rows = []
    for source in MAJOR_STATES:
        if source not in centers:
            continue
        velocity = mean_velocities[source]
        velocity_norm = float(np.linalg.norm(velocity))
        if velocity_norm == 0:
            continue
        for target in MAJOR_STATES:
            if source == target or target not in centers:
                continue
            target_vec = centers[target] - centers[source]
            target_norm = float(np.linalg.norm(target_vec))
            if target_norm == 0:
                continue
            alignment = float(np.dot(velocity, target_vec) / (velocity_norm * target_norm))
            edge_rows.append({
                "sample": sample,
                "dataset": dataset,
                "source_group": dataset,
                "source": source,
                "target": target,
                "velocity_alignment": alignment,
                "source_velocity_norm": velocity_norm,
                "source_to_target_distance": target_norm,
                "source_cells": int((states == source).sum()),
                "target_cells": int((states == target).sum()),
            })
    return node_rows, edge_rows


# ---------------------------------------------------------------------------
# MP-level direction tables (Basal / SMG)
# ---------------------------------------------------------------------------

def mp_direction_tables(adata: ad.AnnData, obs_col: str, mp_list: List[str]) -> Tuple[List[dict], List[dict]]:
    """Compute velocity direction between dominant MPs within a state subset."""
    if obs_col not in adata.obs:
        return [], []

    states = adata.obs[obs_col].astype(str).to_numpy()
    if not isinstance(adata.obsm.get("velocity_umap"), np.ndarray):
        return [], []

    V = adata.obsm["velocity_umap"]
    U = adata.obsm["X_umap"]

    valid_states = [s for s in np.unique(states) if s in mp_list]

    sample = str(adata.obs["sample"].iloc[0])
    dataset = str(adata.obs["dataset"].iloc[0])
    total_direction_cells = int(np.isin(states, valid_states).sum())
    total_major = total_direction_cells

    node_rows = []
    centers = {}
    mean_velocities = {}
    for state in valid_states:
        mask = states == state
        cells = int(mask.sum())
        pct_total = 100 * cells / total_direction_cells if total_direction_cells else 0.0
        node_rows.append({
            "sample": sample,
            "dataset": dataset,
            "source_group": dataset,
            "state": state,
            "cells": cells,
            "total_cells": total_direction_cells,
            "major_cells": total_major,
            "pct_major": pct_total,
            "pct_total_direction_states": pct_total,
            "pct_of_major_states": pct_total,
        })
        if cells >= 5:
            centers[state] = U[mask].mean(axis=0)
            mean_velocities[state] = V[mask].mean(axis=0)

    edge_rows = []
    for source in valid_states:
        if source not in centers:
            continue
        velocity = mean_velocities[source]
        velocity_norm = float(np.linalg.norm(velocity))
        if velocity_norm == 0:
            continue
        for target in valid_states:
            if target == source or target not in centers:
                continue
            target_vec = centers[target] - centers[source]
            target_norm = float(np.linalg.norm(target_vec))
            if target_norm == 0:
                continue
            alignment = float(np.dot(velocity, target_vec) / (velocity_norm * target_norm))
            edge_rows.append({
                "sample": sample,
                "dataset": dataset,
                "source_group": dataset,
                "source": source,
                "target": target,
                "velocity_alignment": alignment,
                "source_velocity_norm": velocity_norm,
                "source_to_target_distance": target_norm,
                "source_cells": int((states == source).sum()),
                "target_cells": int((states == target).sum()),
            })
    return node_rows, edge_rows


# ---------------------------------------------------------------------------
# Categorical color helpers
# ---------------------------------------------------------------------------

def set_categorical_colors(adata: ad.AnnData, column: str, values: List[str], color_map: dict) -> None:
    observed = [str(x) for x in adata.obs[column].astype(str).unique()]
    categories = values + [x for x in observed if x not in values]
    adata.obs[column] = pd.Categorical(adata.obs[column].astype(str), categories=categories, ordered=True)
    adata.uns[f"{column}_colors"] = [color_map.get(value, "#808080") for value in categories]


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def drop_nonfinite_scatter_offsets(fig: plt.Figure) -> None:
    for ax in fig.axes:
        for artist in list(ax.collections):
            if not isinstance(artist, PathCollection):
                continue
            offsets = artist.get_offsets()
            if offsets is None or len(offsets) == 0:
                continue
            offsets_arr = np.asarray(offsets)
            if offsets_arr.ndim != 2 or offsets_arr.shape[1] < 2:
                continue
            keep = np.isfinite(offsets_arr[:, 0]) & np.isfinite(offsets_arr[:, 1])
            if np.all(keep):
                continue
            artist.set_offsets(offsets_arr[keep, :])
            values = artist.get_array()
            if values is not None and len(values) == len(keep):
                artist.set_array(np.asarray(values)[keep])


def point_size(adata: ad.AnnData) -> float:
    return float(max(10, min(72, 52000 / max(adata.n_obs, 1))))


def remove_axis_legend(ax: plt.Axes) -> None:
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()


def remove_axis_text_labels(ax: plt.Axes) -> None:
    for txt in list(ax.texts):
        txt.remove()


def save_raster_page(pdf: PdfPages, fig: plt.Figure, sample: str, suffix: str, figsize: Tuple[float, float]) -> None:
    drop_nonfinite_scatter_offsets(fig)
    tmp_dir = OUT_LIVE / "tmp_pages"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_png = tmp_dir / f"Auto_scatlas_velocity_{suffix}_{sample}.png"
    fig.savefig(tmp_png, dpi=190, bbox_inches="tight")
    plt.close(fig)

    page_img = plt.imread(tmp_png)
    raster_fig, raster_ax = plt.subplots(figsize=figsize)
    raster_ax.imshow(page_img)
    raster_ax.axis("off")
    raster_fig.tight_layout(pad=0)
    pdf.savefig(raster_fig, bbox_inches="tight")
    plt.close(raster_fig)


# ---------------------------------------------------------------------------
# Per-sample plotting (core + extended pages)
# ---------------------------------------------------------------------------

def plot_sample_core_page(pdf: PdfPages, adata: ad.AnnData) -> None:
    sample = str(adata.obs["sample"].iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8))
    axes = axes.ravel()

    final_values = [x for x in STATE_COLORS if x in set(adata.obs["state_final"])]
    direction_values = [x for x in MAJOR_STATES + ["Unassigned"] if x in set(adata.obs["state_direction"])]
    set_categorical_colors(adata, "state_final", final_values, STATE_COLORS)
    set_categorical_colors(adata, "state_direction", direction_values, STATE_COLORS)
    size = point_size(adata)

    sc.pl.embedding(
        adata,
        basis="umap",
        color="state_final",
        frameon=False,
        title="Cell states",
        ax=axes[0],
        show=False,
        size=size,
        legend_loc="right margin",
    )
    sc.pl.embedding(
        adata,
        basis="umap",
        color="state_direction",
        frameon=False,
        title="Direction states (5 primary)",
        ax=axes[1],
        show=False,
        size=size,
        legend_loc=None,
    )
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="state_final",
        legend_loc="none",
        title="Velocity stream",
        frameon=False,
        ax=axes[2],
        show=False,
        size=size,
        alpha=1.0,
    )
    remove_axis_legend(axes[1])
    remove_axis_legend(axes[2])
    remove_axis_text_labels(axes[2])

    fig.suptitle(sample, fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_raster_page(pdf, fig, sample, "core_page", (17, 5.8))


def plot_sample_extended_page(pdf: PdfPages, adata: ad.AnnData) -> None:
    sample = str(adata.obs["sample"].iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8))
    axes = axes.ravel()

    final_values = [x for x in STATE_COLORS if x in set(adata.obs["state_final"])]
    set_categorical_colors(adata, "state_final", final_values, STATE_COLORS)
    size = point_size(adata)

    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="state_final",
        arrow_color="black",
        arrow_size=3.8,
        arrow_length=6.5,
        density=0.65,
        scale=0.35,
        autoscale=False,
        alpha=1.0,
        size=size,
        legend_loc="right margin",
        title="Velocity grid",
        frameon=False,
        ax=axes[0],
        show=False,
    )
    sc.pl.embedding(
        adata,
        basis="umap",
        color="velocity_length",
        frameon=False,
        title="Velocity length",
        ax=axes[1],
        show=False,
        size=size,
        color_map="viridis",
    )
    sc.pl.embedding(
        adata,
        basis="umap",
        color="velocity_confidence",
        frameon=False,
        title="Velocity confidence",
        ax=axes[2],
        show=False,
        size=size,
        color_map="viridis",
    )

    fig.suptitle(sample, fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_raster_page(pdf, fig, sample, "extended_page", (17, 5.8))


# ---------------------------------------------------------------------------
# State-level direction network plots
# ---------------------------------------------------------------------------

def draw_direction_network(ax: plt.Axes, node_df: pd.DataFrame, edge_df: pd.DataFrame, title: str) -> None:
    ax.set_title(title, fontsize=14, fontweight="bold", pad=8)
    max_node = max(float(node_df["pct_major"].max()), 1.0) if len(node_df) else 1.0
    edge_plot = edge_df[edge_df["velocity_alignment"] >= EDGE_THRESHOLD].copy()
    for _, edge in edge_plot.iterrows():
        start = LAYOUT[edge["source"]]
        end = LAYOUT[edge["target"]]
        score = float(edge["velocity_alignment"])
        rad = ARROW_RAD * edge_curve_side(edge["source"], edge["target"])
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=24,
            linewidth=1.2 + 5.4 * score,
            color="black",
            alpha=0.95,
            shrinkA=ARROW_NODE_SHRINK,
            shrinkB=ARROW_NODE_SHRINK,
            connectionstyle=f"arc3,rad={rad}",
            zorder=8,
        )
        ax.add_patch(arrow)
        mid = (start + end) / 2
        edge_vec = end - start
        edge_len = float(np.linalg.norm(edge_vec))
        if edge_len > 0:
            perp = np.array([-edge_vec[1], edge_vec[0]]) / edge_len
            label_xy = mid + perp * rad * 0.82
        else:
            label_xy = mid
        ax.text(label_xy[0], label_xy[1], f"{score:.2f}", ha="center", va="center",
                fontsize=ALIGNMENT_LABEL_FONTSIZE, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="black", linewidth=0.35, alpha=0.92, pad=2.0),
                zorder=12)

    if not len(edge_plot):
        ax.text(0, 0, f"No state-pair alignment >= {EDGE_THRESHOLD:.2f}", ha="center", va="center",
                fontsize=10.5, color="#555555", zorder=2)

    for state in MAJOR_STATES:
        pct = 0.0
        hit = node_df[node_df["state"] == state]
        if len(hit):
            pct = float(hit["pct_major"].iloc[0])
        xy = LAYOUT[state]
        ax.scatter(
            xy[0], xy[1],
            s=300 + 1420 * pct / max_node,
            color=STATE_COLORS[state],
            edgecolor="black",
            linewidth=1.4,
            zorder=10,
        )
        label = format_state_label(state)
        text = ax.text(xy[0], xy[1], f"{label}\n{pct:.1f}%", ha="center", va="center",
                       fontsize=NODE_LABEL_FONTSIZE,
                       fontweight="bold", color="black", zorder=11)
        text.set_path_effects([pe.withStroke(linewidth=3.0, foreground="white")])

    ax.set_xlim(-1.85, 1.85)
    ax.set_ylim(-1.70, 1.22)
    ax.set_aspect("equal")
    ax.axis("off")


def aggregate_direction(nodes: pd.DataFrame, edges: pd.DataFrame, label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if len(nodes):
        denominator = total_cells_from_nodes(nodes)
        node_group = nodes.groupby("state", as_index=False).agg(cells=("cells", "sum"))
        node_group["total_cells"] = denominator
        node_group["pct_major"] = np.where(denominator > 0, 100 * node_group["cells"] / denominator, 0.0)
        node_group["pct_total_direction_states"] = node_group["pct_major"]
    else:
        node_group = pd.DataFrame(columns=["state", "cells", "total_cells", "pct_major", "pct_total_direction_states"])
    node_group = (
        pd.DataFrame({"state": MAJOR_STATES})
        .merge(node_group, on="state", how="left")
        .fillna({"cells": 0, "total_cells": 0, "pct_major": 0.0, "pct_total_direction_states": 0.0})
    )
    if len(edges):
        edge_group = (
            edges.groupby(["source", "target"], as_index=False)
            .agg(
                velocity_alignment=("velocity_alignment", "mean"),
                median_alignment=("velocity_alignment", "median"),
                positive_fraction=("velocity_alignment", lambda x: float((x > 0).mean())),
                sample_n=("sample", "nunique"),
            )
        )
    else:
        edge_group = pd.DataFrame(columns=["source", "target", "velocity_alignment", "median_alignment", "positive_fraction", "sample_n"])
    node_group["group"] = label
    edge_group["group"] = label
    return node_group, edge_group


def draw_direction_page(pdf: PdfPages, panels: List[Tuple[pd.DataFrame, pd.DataFrame, str]]) -> None:
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(7.4 * n, 6.7), squeeze=False)
    for ax, (node_df, edge_df, title) in zip(axes.ravel(), panels):
        draw_direction_network(ax, node_df, edge_df, title)
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_direction_pdf(nodes: pd.DataFrame, edges: pd.DataFrame) -> None:
    pdf_path = OUT_LIVE / "figures" / DIRECTION_PDF
    group_rows = []
    datasets = sorted(nodes["dataset"].dropna().astype(str).unique())

    def panel_for_samples(sample_set: List[str], label: str) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
        sample_set_obj = set(sample_set)
        sub_nodes = nodes[nodes["sample"].isin(sample_set_obj)]
        sub_edges = edges[edges["sample"].isin(sample_set_obj)]
        node_group, edge_group = aggregate_direction(sub_nodes, sub_edges, label)
        group_rows.append(edge_group)
        return node_group, edge_group, title_with_n(label, nodes, sample_set_obj)

    def panel_for_sample(sample: str) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
        sub_nodes = nodes[nodes["sample"] == sample]
        sub_edges = edges[edges["sample"] == sample]
        return sub_nodes, sub_edges, title_with_n(sample, nodes, {sample})

    with PdfPages(pdf_path) as pdf:
        all_samples = sorted(nodes["sample"].unique())
        node_all, edge_all = aggregate_direction(nodes, edges, "All raw-BAM scATLAS")
        group_rows.append(edge_all)
        draw_direction_page(pdf, [(node_all, edge_all, title_with_n("All raw-BAM scATLAS", nodes))])

        dataset_panels = []
        for dataset in datasets:
            sample_set = sorted(nodes.loc[nodes["dataset"].astype(str).eq(dataset), "sample"].unique())
            dataset_panels.append(panel_for_samples(sample_set, f"{dataset} {len(sample_set)} samples"))
        if dataset_panels:
            draw_direction_page(pdf, dataset_panels)

        for sample in all_samples:
            draw_direction_page(pdf, [panel_for_sample(sample)])

    if group_rows:
        pd.concat(group_rows, axis=0).to_csv(
            OUT_LIVE / "tables" / "Auto_scatlas_velocity_group_state_direction_edges.csv",
            index=False,
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    import os
    force_recompute = os.environ.get("SCATLAS_FORCE_VELOCITY_RECOMPUTE", "FALSE").upper() == "TRUE"

    for subdir in ["figures", "h5ad", "tables", "logs", "tmp_pages"]:
        (OUT_LIVE / subdir).mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_cell_metadata.csv")
    manifest = pd.read_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_sample_manifest.csv")
    samples = manifest.sort_values(["dataset", "sample"])["sample"].tolist()
    missing_inputs = [sample for sample in samples if not sample_has_velocity_input(sample)]
    if missing_inputs:
        missing_text = ", ".join(missing_inputs[:10])
        if len(missing_inputs) > 10:
            missing_text += f", ... ({len(missing_inputs)} total)"
        raise FileNotFoundError(
            "Missing both cached h5ad and loom input for scVelo visualisation samples: "
            f"{missing_text}. Run the velocity/velocyto generation chain before plotting."
        )

    all_obs = []
    all_nodes = []
    all_edges = []
    all_basal_nodes = []
    all_basal_edges = []
    all_smg_nodes = []
    all_smg_edges = []
    failures = []

    core_pdf_path = OUT_LIVE / "figures" / CORE_PDF
    extended_pdf_path = OUT_LIVE / "figures" / EXTENDED_PDF
    with PdfPages(core_pdf_path) as core_pdf, PdfPages(extended_pdf_path) as extended_pdf:
        for sample in samples:
            print(f"Loading/scVelo for {sample}", flush=True)
            try:
                adata = load_or_run_velocity(sample, meta, force_recompute=force_recompute)
                plot_sample_core_page(core_pdf, adata)
                plot_sample_extended_page(extended_pdf, adata)
                node_rows, edge_rows = state_direction_tables(adata)
                all_nodes.extend(node_rows)
                all_edges.extend(edge_rows)

                basal_node_rows, basal_edge_rows = mp_direction_tables(adata, "dominant_basal_mp", BASAL_MPS)
                all_basal_nodes.extend(basal_node_rows)
                all_basal_edges.extend(basal_edge_rows)

                smg_node_rows, smg_edge_rows = mp_direction_tables(adata, "dominant_smg_mp", SMG_MPS)
                all_smg_nodes.extend(smg_node_rows)
                all_smg_edges.extend(smg_edge_rows)

                obs = adata.obs.copy()
                obs["velocity_umap_1"] = adata.obsm["velocity_umap"][:, 0]
                obs["velocity_umap_2"] = adata.obsm["velocity_umap"][:, 1]
                all_obs.append(obs)
            except Exception as exc:
                failures.append({"sample": sample, "error": str(exc)})
                print(f"FAILED {sample}: {exc}", flush=True)

    if failures:
        pd.DataFrame(failures).to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_scvelo_failures.csv", index=False)
    if not all_nodes:
        raise RuntimeError("No samples completed scVelo; see Auto_scatlas_velocity_scvelo_failures.csv")

    obs_all = pd.concat(all_obs, axis=0)
    nodes = pd.DataFrame(all_nodes)
    edges = pd.DataFrame(all_edges)
    obs_all.to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_scvelo_cell_metadata.csv")
    nodes.to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_state_nodes.csv", index=False)
    edges.to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_state_direction_edges.csv", index=False)

    if all_basal_nodes:
        pd.DataFrame(all_basal_nodes).to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_basal_mp_nodes.csv", index=False)
        pd.DataFrame(all_basal_edges).to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_basal_mp_direction_edges.csv", index=False)

    if all_smg_nodes:
        pd.DataFrame(all_smg_nodes).to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_smg_mp_nodes.csv", index=False)
        pd.DataFrame(all_smg_edges).to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_smg_mp_direction_edges.csv", index=False)

    if len(edges):
        top_edges = (
            edges.sort_values(["sample", "source", "velocity_alignment"], ascending=[True, True, False])
            .groupby(["sample", "source"], as_index=False)
            .head(1)
        )
        top_edges.to_csv(OUT_LIVE / "tables" / "Auto_scatlas_velocity_top_state_direction_per_source.csv", index=False)
        positive_top = top_edges[top_edges["velocity_alignment"] > 0].copy()
        audit_rows = []
        for label, sample_set in {"All raw-BAM scATLAS": set(nodes["sample"].unique())}.items():
            sub_top = positive_top[positive_top["sample"].isin(sample_set)]
            sub_edges = edges[edges["sample"].isin(sample_set)]
            for target in MAJOR_STATES:
                audit_rows.append({
                    "group": label,
                    "target_state": target,
                    "top_positive_count": int((sub_top["target"] == target).sum()),
                    "mean_alignment_to_target": float(sub_edges.loc[sub_edges["target"] == target, "velocity_alignment"].mean()),
                    "median_alignment_to_target": float(sub_edges.loc[sub_edges["target"] == target, "velocity_alignment"].median()),
                })
        for dataset in sorted(nodes["dataset"].astype(str).unique()):
            sample_set = set(nodes.loc[nodes["dataset"].astype(str).eq(dataset), "sample"])
            sub_top = positive_top[positive_top["sample"].isin(sample_set)]
            sub_edges = edges[edges["sample"].isin(sample_set)]
            for target in MAJOR_STATES:
                audit_rows.append({
                    "group": dataset,
                    "target_state": target,
                    "top_positive_count": int((sub_top["target"] == target).sum()),
                    "mean_alignment_to_target": float(sub_edges.loc[sub_edges["target"] == target, "velocity_alignment"].mean()),
                    "median_alignment_to_target": float(sub_edges.loc[sub_edges["target"] == target, "velocity_alignment"].median()),
                })
        pd.DataFrame(audit_rows).to_csv(
            OUT_LIVE / "tables" / "Auto_scatlas_velocity_direction_audit_summary.csv",
            index=False,
        )

    plot_direction_pdf(nodes, edges)
    print(f"Wrote {core_pdf_path}", flush=True)
    print(f"Wrote {extended_pdf_path}", flush=True)


if __name__ == "__main__":
    main()
