#!/usr/bin/env python
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/map_scatlas_states_xenium.py
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs/outputs: documented in this header below and in the analysis map.
####################


import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns


MP_DESCRIPTIONS = {
    "MP1": "G2M Cell Cycle",
    "MP2": "MYC-related Proliferation",
    "MP5": "Epithelial IFN Resp.",
    "MP7": "DNA Damage Repair",
    "MP8": "Intestinal Diff.",
    "MP9": "G1S Cell Cycle",
    "MP10": "Columnar Diff.",
    "MP12": "Neuro-responsive Epi",
    "MP13": "Hypoxic Inflam. Epi.",
    "MP14": "Hypoxia Adapted Epi.",
    "MP15": "Immune Infiltration",
    "MP16": "Secretory Diff. (Gastric)",
    "MP17": "Basal-like Transition",
    "MP18": "Secretory Diff. (Intest.)",
}

CC_MPS = ["MP1", "MP7", "MP9"]

STATE_GROUPS = {
    "Classic Proliferative": ["MP2"],
    "Basal to Intestinal Metaplasia": ["MP17", "MP14", "MP5", "MP10", "MP8"],
    "SMG-like Metaplasia": ["MP18", "MP16"],
    "Stress-adaptive": ["MP13", "MP12"],
    "Immune Infiltrating": ["MP15"],
}

STATE_COLORS = {
    "Classic Proliferative": "#E41A1C",
    "Basal to Intestinal Metaplasia": "#4DAF4A",
    "SMG-like Metaplasia": "#FF7F00",
    "Stress-adaptive": "#984EA3",
    "Immune Infiltrating": "#377EB8",
    "Unresolved": "#8F8F8F",
    "Hybrid": "#111111",
    "Non-carcinoma": "#D9D9D9",
}

MP_COLORS = {
    "MP1": "#B0B0B0",
    "MP2": "#E41A1C",
    "MP5": "#66C2A5",
    "MP7": "#999999",
    "MP8": "#FC8D62",
    "MP9": "#C0C0C0",
    "MP10": "#A6D854",
    "MP12": "#E78AC3",
    "MP13": "#984EA3",
    "MP14": "#8DA0CB",
    "MP15": "#377EB8",
    "MP16": "#FFD92F",
    "MP17": "#4DAF4A",
    "MP18": "#FF7F00",
    "Unresolved": "#8F8F8F",
    "Non-carcinoma": "#D9D9D9",
}

STATE_ORDER = [
    "Classic Proliferative",
    "Basal to Intestinal Metaplasia",
    "SMG-like Metaplasia",
    "Stress-adaptive",
    "Immune Infiltrating",
    "Unresolved",
    "Hybrid",
    "Non-carcinoma",
]

STATE_SCORE_ORDER = [
    "Classic Proliferative",
    "Basal to Intestinal Metaplasia",
    "SMG-like Metaplasia",
    "Stress-adaptive",
    "Immune Infiltrating",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map scATLAS metaprograms and Approach B states onto the Yates et al. Xenium cohort."
    )
    parser.add_argument(
        "--dataset-root",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Zenodo_upload",
    )
    parser.add_argument(
        "--output-dir",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_xenium_states",
    )
    parser.add_argument(
        "--signature-dir",
        default=None,
        help="Directory containing Auto_scATLAS signature CSVs. Defaults to output-dir.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Use the top N ranked genes per retained MP for cell scoring.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=5,
        help="Drop an MP score if fewer than this many top-N genes overlap the Xenium panel.",
    )
    parser.add_argument(
        "--assignment-celltypes",
        default="Carcinoma",
        help="Comma-separated Cell_type values to assign scATLAS malignant states to.",
    )
    parser.add_argument(
        "--baseline-celltypes",
        default="Carcinoma,Epithelial,Epithelial (Pre-cancerous)",
        help="Comma-separated Cell_type values used for per-patient gene centering/scaling.",
    )
    parser.add_argument(
        "--threshold-control-celltypes",
        default="Epithelial",
        help="Comma-separated non-malignant Cell_type values used to calibrate the automatic threshold.",
    )
    parser.add_argument(
        "--threshold",
        default="auto",
        help="Use 'auto' or a numeric Approach B score threshold.",
    )
    parser.add_argument(
        "--threshold-quantile",
        type=float,
        default=0.95,
        help="For --threshold auto, use this quantile of control-cell best group scores.",
    )
    parser.add_argument(
        "--min-threshold",
        type=float,
        default=0.5,
        help="For --threshold auto, never use a threshold below this value.",
    )
    parser.add_argument(
        "--hybrid-gap",
        type=float,
        default=0.3,
        help="Approach B top-two group score gap below which a carcinoma cell is labelled Hybrid.",
    )
    return parser.parse_args()


def split_values(value):
    return [x.strip() for x in value.split(",") if x.strip()]


def decode_scalar(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def decode_array(values):
    if hasattr(values, "dtype") and values.dtype.kind in {"S", "O"}:
        return np.array([decode_scalar(x) for x in values], dtype=object)
    return np.array(values)


def read_series(node):
    if isinstance(node, h5py.Dataset):
        values = node[...]
        return pd.Series(decode_array(values))

    categories = decode_array(node["categories"][...])
    codes = node["codes"][...]
    values = np.empty(codes.shape[0], dtype=object)
    valid = codes >= 0
    values[valid] = categories[codes[valid]]
    values[~valid] = np.nan
    return pd.Series(values)


def read_dataframe(group):
    index_name = decode_scalar(group.attrs["_index"])
    column_order = [decode_scalar(x) for x in group.attrs["column-order"]]
    index = pd.Index(read_series(group[index_name]).astype(str), name=index_name)
    df = pd.DataFrame(index=index)
    for column in column_order:
        series = read_series(group[column])
        series.index = index
        df[column] = series.values
    return df


def read_sparse_matrix(group):
    encoding = decode_scalar(group.attrs["encoding-type"])
    shape = tuple(group.attrs["shape"])
    data = group["data"][...]
    indices = group["indices"][...]
    indptr = group["indptr"][...]
    if encoding == "csr_matrix":
        return sp.csr_matrix((data, indices, indptr), shape=shape)
    if encoding == "csc_matrix":
        return sp.csc_matrix((data, indices, indptr), shape=shape)
    raise ValueError(f"Unsupported sparse encoding: {encoding}")


def derive_stage(title):
    title_low = str(title).lower()
    if "metastasis" in title_low:
        return "metastasis"
    if "primary" in title_low:
        return "primary"
    return "other"


def label_mp(mp_name):
    if mp_name in MP_DESCRIPTIONS:
        return f"{mp_name}: {MP_DESCRIPTIONS[mp_name]}"
    return mp_name


def slug_state(state_name):
    return "Auto_group_" + state_name.lower().replace(" ", "_").replace("-", "_")


def load_signatures(signature_dir, top_n):
    ranked_path = signature_dir / "Auto_scATLAS_mp_gene_ranked.csv"
    state_groups_path = signature_dir / "Auto_scATLAS_state_groups.csv"
    summary_path = signature_dir / "Auto_scATLAS_mp_signature_summary.csv"
    order_path = signature_dir / "Auto_scATLAS_mp_order.csv"

    missing = [p for p in [ranked_path, state_groups_path, summary_path] if not p.exists()]
    if missing:
        missing_str = ", ".join(str(p) for p in missing)
        raise FileNotFoundError(f"Missing signature files: {missing_str}")

    ranked = pd.read_csv(ranked_path)
    ranked["rank"] = ranked["rank"].astype(int)
    ranked = ranked[ranked["rank"] <= top_n].copy()

    state_groups_df = pd.read_csv(state_groups_path)
    signature_summary = pd.read_csv(summary_path)
    if order_path.exists():
        mp_order = pd.read_csv(order_path).sort_values("plot_order")["mp"].astype(str).tolist()
    else:
        mp_order = list(MP_DESCRIPTIONS)
    return ranked, state_groups_df, signature_summary, mp_order


def build_signature_map(ranked, gene_index, min_genes, output_dir):
    signature_map = {}
    overlap_rows = []

    for mp_name, mp_df in ranked.groupby("mp", sort=False):
        genes = mp_df.sort_values("rank")["gene"].astype(str).tolist()
        idx = gene_index.get_indexer(genes)
        valid = idx >= 0
        used_idx = idx[valid].tolist()
        used_genes = [gene for gene, keep in zip(genes, valid) if keep]
        dropped = len(used_idx) < min_genes
        if not dropped:
            signature_map[mp_name] = used_idx
        overlap_rows.append(
            {
                "mp": mp_name,
                "description": MP_DESCRIPTIONS.get(mp_name, mp_name),
                "n_genes_requested_top_n": len(genes),
                "n_genes_overlap": int(valid.sum()),
                "min_genes_required": int(min_genes),
                "used_for_scoring": not dropped,
                "overlap_genes": ";".join(used_genes),
            }
        )

    overlap_df = pd.DataFrame(overlap_rows).sort_values("mp")
    overlap_df.to_csv(output_dir / "Auto_scATLAS_xenium_mp_overlap_summary.csv", index=False)
    if len(signature_map) == 0:
        raise ValueError("No MP signatures passed the Xenium overlap threshold.")
    return signature_map, overlap_df


def score_mp_signatures(norm_matrix, obs, signature_map, score_universe):
    patients = obs["patient"].astype(str).to_numpy()
    score_universe = np.asarray(score_universe, dtype=bool)
    mp_names = list(signature_map)
    union_idx = sorted(set(gene_idx for values in signature_map.values() for gene_idx in values))
    union_pos = {gene_idx: i for i, gene_idx in enumerate(union_idx)}
    signature_pos = {
        mp_name: [union_pos[gene_idx] for gene_idx in gene_indices]
        for mp_name, gene_indices in signature_map.items()
    }

    score_frames = []
    for patient in pd.unique(patients):
        row_idx = np.where((patients == patient) & score_universe)[0]
        if len(row_idx) == 0:
            continue

        dense = norm_matrix[row_idx][:, union_idx].toarray().astype(np.float32, copy=False)
        gene_mean = dense.mean(axis=0)
        gene_sd = dense.std(axis=0)
        gene_sd[~np.isfinite(gene_sd) | (gene_sd == 0)] = 1.0
        gene_z = (dense - gene_mean) / gene_sd

        sample_scores = np.empty((len(row_idx), len(mp_names)), dtype=np.float32)
        for col_idx, mp_name in enumerate(mp_names):
            positions = signature_pos[mp_name]
            if len(positions) == 0:
                sample_scores[:, col_idx] = np.nan
            else:
                sample_scores[:, col_idx] = gene_z[:, positions].mean(axis=1)

        sample_df = pd.DataFrame(sample_scores, index=obs.index[row_idx], columns=mp_names)
        score_frames.append(sample_df)

    if len(score_frames) == 0:
        raise ValueError("No cells were available for scATLAS Xenium scoring.")

    mp_raw_scored = pd.concat(score_frames, axis=0)
    scored_index = obs.index[score_universe]
    mp_raw_scored = mp_raw_scored.reindex(scored_index)
    mp_raw = pd.DataFrame(np.nan, index=obs.index, columns=mp_names, dtype=np.float32)
    mp_raw.loc[scored_index, :] = mp_raw_scored

    mp_adj = mp_raw.copy()
    for patient in pd.unique(patients):
        patient_mask = (patients == patient) & score_universe
        if patient_mask.sum() == 0:
            continue
        patient_index = obs.index[patient_mask]
        mp_adj.loc[patient_index, :] = mp_adj.loc[patient_index, :] - mp_adj.loc[patient_index, :].mean(axis=0)

    global_sd = mp_adj.loc[scored_index, :].std(axis=0, ddof=0)
    global_sd[~np.isfinite(global_sd) | (global_sd == 0)] = 1.0
    mp_adj = mp_adj.divide(global_sd, axis=1)
    return mp_raw, mp_adj


def compute_group_scores(mp_adj_noncc):
    group_scores = {}
    for state_name, mps in STATE_GROUPS.items():
        available = [mp for mp in mps if mp in mp_adj_noncc.columns]
        if len(available) == 0:
            group_scores[state_name] = np.nan
        elif len(available) == 1:
            group_scores[state_name] = mp_adj_noncc[available[0]]
        else:
            group_scores[state_name] = mp_adj_noncc[available].max(axis=1)
    return pd.DataFrame(group_scores, index=mp_adj_noncc.index)


def assign_states(group_scores, threshold, hybrid_gap):
    best_state = group_scores.idxmax(axis=1)
    best_value = group_scores.max(axis=1)
    ordered = np.sort(group_scores.to_numpy(dtype=float), axis=1)
    top1 = ordered[:, -1]
    top2 = ordered[:, -2] if ordered.shape[1] > 1 else np.zeros_like(top1)
    gap = top1 - top2

    state = best_state.astype(object)
    state[best_value < threshold] = "Unresolved"
    state[(gap < hybrid_gap) & (state != "Unresolved")] = "Hybrid"
    return pd.Series(state, index=group_scores.index, name="Auto_state_B"), pd.Series(gap, index=group_scores.index, name="Auto_state_gap")


def top_mp_labels(mp_adj_noncc, threshold, mp_order):
    available_order = [mp for mp in mp_order if mp in mp_adj_noncc.columns]
    extra_mps = [mp for mp in mp_adj_noncc.columns if mp not in available_order]
    mp_adj_noncc = mp_adj_noncc[available_order + extra_mps]
    top_mp = mp_adj_noncc.idxmax(axis=1)
    top_val = mp_adj_noncc.max(axis=1)
    top_mp = top_mp.astype(object)
    top_mp[top_val < threshold] = "Unresolved"
    return pd.Series(top_mp, index=mp_adj_noncc.index, name="Auto_top_mp")


def select_threshold(group_scores, assignment_mask, control_mask, args, output_dir):
    best_scores = group_scores.max(axis=1)
    assignment_best = best_scores[assignment_mask]
    control_best = best_scores[control_mask]

    if args.threshold.lower() == "auto":
        control_best_clean = control_best.dropna()
        if len(control_best_clean) > 0:
            control_quantile = float(control_best_clean.quantile(args.threshold_quantile))
            threshold = max(float(args.min_threshold), control_quantile)
            method = "max(min_threshold, threshold_control_quantile)"
        else:
            control_quantile = np.nan
            threshold = float(args.min_threshold)
            method = "min_threshold_no_control_cells"
    else:
        threshold = float(args.threshold)
        control_quantile = np.nan
        method = "manual"

    quantiles = [0.50, 0.75, 0.90, 0.95, 0.99]
    rows = [
        {
            "metric": "selected_threshold",
            "value": threshold,
            "method": method,
            "threshold_quantile": args.threshold_quantile,
            "min_threshold": args.min_threshold,
            "n_assignment_cells": int(assignment_mask.sum()),
            "n_control_cells": int(control_mask.sum()),
        },
        {
            "metric": "control_quantile",
            "value": control_quantile,
            "method": method,
            "threshold_quantile": args.threshold_quantile,
            "min_threshold": args.min_threshold,
            "n_assignment_cells": int(assignment_mask.sum()),
            "n_control_cells": int(control_mask.sum()),
        },
    ]
    for q in quantiles:
        rows.append(
            {
                "metric": f"assignment_best_q{int(q * 100)}",
                "value": float(assignment_best.quantile(q)),
                "method": method,
                "threshold_quantile": args.threshold_quantile,
                "min_threshold": args.min_threshold,
                "n_assignment_cells": int(assignment_mask.sum()),
                "n_control_cells": int(control_mask.sum()),
            }
        )
        rows.append(
            {
                "metric": f"control_best_q{int(q * 100)}",
                "value": float(control_best.quantile(q)) if len(control_best.dropna()) > 0 else np.nan,
                "method": method,
                "threshold_quantile": args.threshold_quantile,
                "min_threshold": args.min_threshold,
                "n_assignment_cells": int(assignment_mask.sum()),
                "n_control_cells": int(control_mask.sum()),
            }
        )

    diagnostics = pd.DataFrame(rows)
    diagnostics.to_csv(output_dir / "Auto_xenium_threshold_diagnostics.csv", index=False)
    return threshold, diagnostics


def make_threshold_sensitivity(group_scores, results, assignment_mask, selected_threshold, hybrid_gap, output_dir):
    candidate_thresholds = sorted(
        set([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, round(float(selected_threshold), 4)])
    )
    rows = []
    assignment_index = results.index[assignment_mask]
    patient_vec = results.loc[assignment_index, "patient"].astype(str)
    for threshold in candidate_thresholds:
        states, _ = assign_states(group_scores.loc[assignment_index], threshold, hybrid_gap)
        tmp = pd.DataFrame({"patient": patient_vec, "Auto_state_B": states})
        summary = tmp.groupby(["patient", "Auto_state_B"], observed=False).size().rename("cells").reset_index()
        totals = summary.groupby("patient")["cells"].transform("sum")
        summary["pct"] = 100.0 * summary["cells"] / totals
        summary["threshold"] = threshold
        rows.append(summary)
    sensitivity = pd.concat(rows, axis=0)
    sensitivity.to_csv(output_dir / "Auto_xenium_threshold_sensitivity_by_patient.csv", index=False)
    return sensitivity


def plot_categorical(ax, df, color_col, palette, title):
    values = df[color_col].astype(str)
    present = [cat for cat in palette if cat in set(values)]
    extras = [cat for cat in pd.unique(values) if cat not in present]
    order = [cat for cat in ["Non-carcinoma", "Unresolved", "Hybrid"] if cat in present]
    order += [cat for cat in present if cat not in order] + extras

    for cat in order:
        mask = values == cat
        if cat == "Non-carcinoma":
            point_size = 0.035 if len(df) > 100000 else 0.08
            alpha = 0.38
            zorder = 1
        elif cat in {"Unresolved", "Hybrid"}:
            point_size = 0.08 if len(df) > 100000 else 0.16
            alpha = 0.62
            zorder = 2
        else:
            point_size = 0.34 if len(df) > 100000 else 0.62
            alpha = 0.96
            zorder = 3
        ax.scatter(
            df.loc[mask, "spatial_x"],
            df.loc[mask, "spatial_y"],
            s=point_size,
            c=palette.get(cat, "#808080"),
            linewidths=0,
            alpha=alpha,
            label=cat,
            rasterized=True,
            zorder=zorder,
        )

    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    if order:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            frameon=False,
            fontsize=7,
            markerscale=10,
        )


def plot_top_mp(ax, df, title, mp_order):
    values = df["Auto_top_mp_label"].astype(str)
    mp_labels = [label_mp(mp) for mp in mp_order if label_mp(mp) in set(values)]
    order = ["Non-carcinoma"] + mp_labels + [x for x in ["Unresolved"] if x in set(values)]
    extras = [cat for cat in pd.unique(values) if cat not in order]
    order = order + extras
    palette = {"Non-carcinoma": MP_COLORS["Non-carcinoma"], "Unresolved": MP_COLORS["Unresolved"]}
    palette.update({label_mp(mp): MP_COLORS.get(mp, "#808080") for mp in mp_order})
    plot_categorical(ax, df.assign(Auto_top_mp_label=values), "Auto_top_mp_label", palette, title)


def plot_continuous(ax, df, color_col, title, vmax):
    non_assignment = ~df["Auto_is_assignment_cell"].astype(bool)
    if non_assignment.any():
        ax.scatter(
            df.loc[non_assignment, "spatial_x"],
            df.loc[non_assignment, "spatial_y"],
            s=0.035 if len(df) > 100000 else 0.08,
            c=STATE_COLORS["Non-carcinoma"],
            linewidths=0,
            alpha=0.35,
            label="Non-carcinoma",
            rasterized=True,
        )
    tumor_df = df.loc[~non_assignment].copy()
    values = tumor_df[color_col].to_numpy(dtype=float)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    sc = ax.scatter(
        tumor_df["spatial_x"],
        tumor_df["spatial_y"],
        s=0.22 if len(df) > 100000 else 0.46,
        c=values,
        cmap="RdBu_r",
        norm=norm,
        linewidths=0,
        alpha=0.95,
        rasterized=True,
    )
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    if non_assignment.any():
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=7, markerscale=10)
    return sc


def make_state_maps(results, output_dir, mp_order):
    state_pdf = output_dir / "Auto_xenium_spatial_states.pdf"
    score_pdf = output_dir / "Auto_xenium_spatial_group_scores.pdf"

    with PdfPages(state_pdf) as pdf:
        for patient in pd.unique(results["patient"]):
            sub = results.loc[results["patient"] == patient].copy()
            stage = sub["Auto_stage"].dropna().astype(str)
            stage_label = stage.iloc[0] if len(stage) > 0 else "unknown"
            sample_title = f"{patient} ({stage_label})"
            fig, axes = plt.subplots(1, 2, figsize=(15, 7))
            plot_categorical(axes[0], sub, "Auto_state_B", STATE_COLORS, f"{sample_title}: scATLAS state")
            plot_top_mp(axes[1], sub, f"{sample_title}: top non-CC MP", mp_order)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight", dpi=900)
            fig.savefig(output_dir / f"Auto_{patient}_xenium_state_map.png", dpi=900, bbox_inches="tight")
            plt.close(fig)

    with PdfPages(score_pdf) as pdf:
        for patient in pd.unique(results["patient"]):
            sub = results.loc[results["patient"] == patient].copy()
            tumor_values = sub.loc[sub["Auto_is_assignment_cell"], [slug_state(x) for x in STATE_SCORE_ORDER]].to_numpy(dtype=float)
            vmax = np.nanquantile(np.abs(tumor_values), 0.98)
            if not np.isfinite(vmax) or vmax == 0:
                vmax = 1.0

            fig, axes = plt.subplots(2, 3, figsize=(16, 10))
            axes_flat = axes.flatten()
            stage = sub["Auto_stage"].dropna().astype(str)
            stage_label = stage.iloc[0] if len(stage) > 0 else "unknown"
            sample_title = f"{patient} ({stage_label})"
            for i, state_name in enumerate(STATE_SCORE_ORDER):
                sc = plot_continuous(axes_flat[i], sub, slug_state(state_name), f"{sample_title}: {state_name}", vmax)
                fig.colorbar(sc, ax=axes_flat[i], fraction=0.046, pad=0.04)
            axes_flat[-1].axis("off")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight", dpi=900)
            fig.savefig(output_dir / f"Auto_{patient}_xenium_group_scores.png", dpi=900, bbox_inches="tight")
            plt.close(fig)


def summarise_counts(df, group_col, state_col):
    summary = df.groupby([group_col, state_col], observed=False).size().rename("cells").reset_index()
    totals = summary.groupby(group_col)["cells"].transform("sum")
    summary["pct"] = 100.0 * summary["cells"] / totals
    return summary


def make_assignment_score_matrix(sample_df):
    rows = [state for state in STATE_SCORE_ORDER + ["Unresolved", "Hybrid"] if (sample_df["Auto_state_B"] == state).any()]
    score_cols = [slug_state(state) for state in STATE_SCORE_ORDER]
    out = []
    for state in rows:
        vals = sample_df.loc[sample_df["Auto_state_B"] == state, score_cols].mean(axis=0)
        vals.index = STATE_SCORE_ORDER
        vals.name = state
        out.append(vals)
    total = sample_df[score_cols].mean(axis=0)
    total.index = STATE_SCORE_ORDER
    total.name = "Total carcinoma cells"
    out.append(total)
    return pd.DataFrame(out)


def make_summary_plots(results, output_dir, mp_order):
    tumor_results = results.loc[results["Auto_is_assignment_cell"]].copy()
    by_patient = summarise_counts(tumor_results, "patient", "Auto_state_B")
    mp_label_order = [label_mp(mp) for mp in mp_order if mp not in CC_MPS] + ["Unresolved"]
    by_top_mp = summarise_counts(tumor_results, "patient", "Auto_top_mp_label")
    by_top_mp["Auto_top_mp_label"] = pd.Categorical(by_top_mp["Auto_top_mp_label"], categories=mp_label_order, ordered=True)
    by_top_mp = by_top_mp.sort_values(["patient", "Auto_top_mp_label"])

    by_patient.to_csv(output_dir / "Auto_xenium_state_summary_by_patient.csv", index=False)
    by_top_mp.to_csv(output_dir / "Auto_xenium_top_mp_summary_by_patient.csv", index=False)

    score_cols = [slug_state(state) for state in STATE_SCORE_ORDER]
    state_means = tumor_results.groupby("patient", observed=False)[score_cols].mean().reset_index()
    state_means.columns = ["patient"] + STATE_SCORE_ORDER
    state_means.to_csv(output_dir / "Auto_xenium_group_score_means_by_patient.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    sample_pivot = (
        by_patient.pivot(index="patient", columns="Auto_state_B", values="pct")
        .reindex(columns=STATE_ORDER)
        .drop(columns=["Non-carcinoma"], errors="ignore")
        .fillna(0.0)
    )
    sample_pivot.plot(
        kind="bar",
        stacked=True,
        color=[STATE_COLORS[state] for state in sample_pivot.columns],
        ax=axes[0],
        width=0.8,
    )
    axes[0].set_title("State proportions by Xenium sample")
    axes[0].set_ylabel("% carcinoma cells")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].legend(frameon=False, title="State", bbox_to_anchor=(1.02, 1), loc="upper left")

    top_mp_pivot = (
        by_top_mp.pivot(index="patient", columns="Auto_top_mp_label", values="pct")
        .reindex(columns=mp_label_order)
        .fillna(0.0)
    )
    mp_palette = {label_mp(mp): MP_COLORS.get(mp, "#808080") for mp in mp_order if mp not in CC_MPS}
    mp_palette["Unresolved"] = MP_COLORS["Unresolved"]
    top_mp_pivot.plot(
        kind="bar",
        stacked=True,
        color=[mp_palette[col] for col in top_mp_pivot.columns],
        ax=axes[1],
        width=0.8,
    )
    axes[1].set_title("Top MP abundance by Xenium sample")
    axes[1].set_ylabel("% carcinoma cells")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(frameon=False, title="Top MP", bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(output_dir / "Auto_xenium_state_distribution_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "Auto_xenium_state_distribution_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    heatmap_pdf = output_dir / "Auto_xenium_state_score_heatmaps_by_patient.pdf"
    with PdfPages(heatmap_pdf) as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes_flat = axes.flatten()
        all_vals = []
        matrices = {}
        for patient in pd.unique(results["patient"]):
            sample_df = tumor_results.loc[tumor_results["patient"] == patient].copy()
            mat = make_assignment_score_matrix(sample_df)
            matrices[patient] = mat
            all_vals.append(mat.to_numpy().ravel())
        lim = np.nanquantile(np.abs(np.concatenate(all_vals)), 0.98)
        if not np.isfinite(lim) or lim == 0:
            lim = 1.0
        for i, patient in enumerate(pd.unique(results["patient"])):
            ax = axes_flat[i]
            sns.heatmap(
                matrices[patient],
                cmap="RdBu_r",
                center=0,
                vmin=-lim,
                vmax=lim,
                linewidths=0.4,
                linecolor="white",
                cbar=i == 2,
                ax=ax,
            )
            ax.set_title(patient)
            ax.set_xlabel("State marker score")
            ax.set_ylabel("Assigned state")
            ax.tick_params(axis="x", rotation=45)
            ax.tick_params(axis="y", rotation=0)
        for ax in axes_flat[len(matrices):]:
            ax.axis("off")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(output_dir / "Auto_xenium_state_score_heatmaps_by_patient.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main():
    args = parse_args()
    dataset_root = Path(args.dataset_root)
    output_dir = Path(args.output_dir)
    signature_dir = Path(args.signature_dir) if args.signature_dir else output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    xenium_path = dataset_root / "processed" / "Xenium.h5ad"
    sample_info_path = dataset_root / "sample_info.csv"

    assignment_celltypes = split_values(args.assignment_celltypes)
    baseline_celltypes = split_values(args.baseline_celltypes)
    threshold_control_celltypes = split_values(args.threshold_control_celltypes)

    ranked, state_groups_df, signature_summary, mp_order = load_signatures(signature_dir, args.top_n)

    with h5py.File(xenium_path, "r") as handle:
        obs = read_dataframe(handle["obs"])
        var = read_dataframe(handle["var"])
        spatial = handle["obsm"]["spatial"][...]
        norm_matrix = read_sparse_matrix(handle["X"]).tocsr()

    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)
    for column in ["cell_id", "Cell_type", "patient", "batch", "cluster_cellcharter"]:
        obs[column] = obs[column].astype(str)
    for column in ["area", "total_counts", "cNMF_1", "cNMF_2", "cNMF_3", "cNMF_4", "cNMF_5"]:
        if column in obs.columns:
            obs[column] = pd.to_numeric(obs[column], errors="coerce")

    results = obs.copy()
    results["spatial_x"] = spatial[:, 0]
    results["spatial_y"] = spatial[:, 1]

    sample_info = pd.read_csv(sample_info_path, encoding="utf-8-sig")
    sample_info = sample_info[sample_info["library name"].str.endswith("_Xenium")].copy()
    sample_info["patient"] = sample_info["library name"].str.replace("_Xenium", "", regex=False)
    sample_info["Auto_stage"] = sample_info["title"].apply(derive_stage)
    sample_info = sample_info[["patient", "title", "tissue", "description", "Auto_stage"]]
    results = results.merge(sample_info, on="patient", how="left")
    results.index = obs.index

    gene_index = pd.Index(var.index)
    signature_map, overlap_df = build_signature_map(ranked, gene_index, args.min_genes, output_dir)
    state_groups_df.to_csv(output_dir / "Auto_scATLAS_xenium_state_groups_resolved.csv", index=False)
    signature_summary.to_csv(output_dir / "Auto_scATLAS_xenium_mp_signature_summary_used.csv", index=False)

    assignment_mask = results["Cell_type"].isin(assignment_celltypes).to_numpy()
    score_universe = results["Cell_type"].isin(baseline_celltypes).to_numpy()
    control_mask = results["Cell_type"].isin(threshold_control_celltypes).to_numpy() & score_universe
    if int(assignment_mask.sum()) == 0:
        raise ValueError(f"No assignment cells found for Cell_type values: {assignment_celltypes}")
    if int(score_universe.sum()) == 0:
        raise ValueError(f"No baseline cells found for Cell_type values: {baseline_celltypes}")

    mp_raw, mp_adj = score_mp_signatures(norm_matrix, results, signature_map, score_universe)

    retained_noncc = [mp for mp in mp_order if mp in mp_adj.columns and mp not in CC_MPS]
    retained_noncc += [mp for mp in mp_adj.columns if mp not in retained_noncc and mp not in CC_MPS]
    mp_adj_noncc = mp_adj[retained_noncc].copy()
    group_scores = compute_group_scores(mp_adj_noncc)
    threshold, diagnostics = select_threshold(group_scores, assignment_mask, control_mask, args, output_dir)

    assignment_index = results.index[assignment_mask]
    state_assignment_tumor, state_gap_tumor = assign_states(
        group_scores.loc[assignment_index],
        threshold,
        args.hybrid_gap,
    )
    top_mp_tumor = top_mp_labels(mp_adj_noncc.loc[assignment_index], threshold, mp_order)

    state_assignment = pd.Series("Non-carcinoma", index=results.index, name="Auto_state_B", dtype=object)
    state_assignment.loc[assignment_index] = state_assignment_tumor
    state_gap = pd.Series(np.nan, index=results.index, name="Auto_state_gap", dtype=float)
    state_gap.loc[assignment_index] = state_gap_tumor
    top_mp = pd.Series("Non-carcinoma", index=results.index, name="Auto_top_mp", dtype=object)
    top_mp.loc[assignment_index] = top_mp_tumor
    best_group_score = group_scores.max(axis=1).rename("Auto_best_group_score")

    group_score_export = group_scores.rename(columns={state: slug_state(state) for state in group_scores.columns})
    results = pd.concat(
        [
            results,
            mp_raw.add_prefix("Auto_raw_"),
            mp_adj.add_prefix("Auto_adj_"),
            group_score_export,
            state_assignment,
            state_gap,
            top_mp,
            best_group_score,
        ],
        axis=1,
    )
    results["Auto_top_mp_label"] = results["Auto_top_mp"].map(label_mp)
    results.loc[results["Auto_top_mp"].eq("Non-carcinoma"), "Auto_top_mp_label"] = "Non-carcinoma"
    results["Auto_is_assignment_cell"] = assignment_mask
    results["Auto_score_universe"] = score_universe
    results["Auto_threshold_control"] = control_mask

    make_threshold_sensitivity(group_scores, results, assignment_mask, threshold, args.hybrid_gap, output_dir)
    results.to_csv(output_dir / "Auto_xenium_cell_annotations.csv.gz", index=True, compression="gzip")

    make_state_maps(results, output_dir, mp_order)
    make_summary_plots(results, output_dir, mp_order)

    with open(output_dir / "Auto_xenium_state_mapping_parameters.txt", "w", encoding="utf-8") as handle:
        handle.write(f"dataset_root={dataset_root}\n")
        handle.write(f"xenium_path={xenium_path}\n")
        handle.write(f"signature_dir={signature_dir}\n")
        handle.write(f"top_n={args.top_n}\n")
        handle.write(f"min_genes={args.min_genes}\n")
        handle.write(f"assignment_celltypes={','.join(assignment_celltypes)}\n")
        handle.write(f"baseline_celltypes={','.join(baseline_celltypes)}\n")
        handle.write(f"threshold_control_celltypes={','.join(threshold_control_celltypes)}\n")
        handle.write(f"threshold={threshold}\n")
        handle.write(f"threshold_request={args.threshold}\n")
        handle.write(f"threshold_quantile={args.threshold_quantile}\n")
        handle.write(f"min_threshold={args.min_threshold}\n")
        handle.write(f"hybrid_gap={args.hybrid_gap}\n")
        handle.write("expression_matrix=.X log1p(CP10K) from processed/Xenium.h5ad\n")
        handle.write("score_subset=baseline Cell_type set; state assignment subset=assignment Cell_type set\n")
        handle.write("platform=10x Xenium single-cell spatial transcriptomics\n")
        handle.write("python_env=/rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp\n")
        handle.write(f"n_cells_total={len(results)}\n")
        handle.write(f"n_assignment_cells={int(assignment_mask.sum())}\n")
        handle.write(f"n_score_universe_cells={int(score_universe.sum())}\n")
        handle.write(f"n_threshold_control_cells={int(control_mask.sum())}\n")


if __name__ == "__main__":
    main()
