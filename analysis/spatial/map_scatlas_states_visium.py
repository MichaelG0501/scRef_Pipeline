#!/usr/bin/env python
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/map_scatlas_states_visium.py
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
from matplotlib.patches import Patch
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
    "Hybrid": "black",
    "Normal/Mixed": "#E6E6E6",
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
    "Normal/Mixed": "#E6E6E6",
}

CNV_COLORS = {"Normal": "#F7F7F7", "Mixed": "#BDBDBD", "Tumor": "#111111"}
SMOOTHED_COLORS = {"Normal": "#D9D9D9", "Mixed": "#969696", "Tumor": "#08519C"}
STAGE_COLORS = {"primary": "#4C78A8", "metastasis": "#F58518", "other": "#54A24B"}

STATE_ORDER = [
    "Classic Proliferative",
    "Basal to Intestinal Metaplasia",
    "SMG-like Metaplasia",
    "Stress-adaptive",
    "Immune Infiltrating",
    "Unresolved",
    "Hybrid",
    "Normal/Mixed",
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
        description="Map scATLAS metaprograms and Approach B states onto the Yates et al. Visium cohort."
    )
    parser.add_argument(
        "--dataset-root",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Zenodo_upload",
    )
    parser.add_argument(
        "--output-dir",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_visium_initial",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Use the top N ranked genes per retained MP for spot scoring.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Approach B threshold for assigning a real state.",
    )
    parser.add_argument(
        "--hybrid-gap",
        type=float,
        default=0.3,
        help="Approach B top-two group score gap below which a spot is labelled Hybrid.",
    )
    return parser.parse_args()


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


def read_spatial_table(group):
    index = pd.Index(read_series(group["_index"]).astype(str), name="spot")
    df = pd.DataFrame(index=index)
    for column in group.keys():
        if column == "_index":
            continue
        values = group[column][...]
        df[column] = values
    return df


def normalise_log1p(matrix):
    row_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    scale = np.zeros_like(row_sums)
    valid = row_sums > 0
    scale[valid] = 1e4 / row_sums[valid]
    norm = matrix.multiply(scale[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    return norm


def derive_stage(title):
    title_low = str(title).lower()
    if "metastasis" in title_low:
        return "metastasis"
    if "primary" in title_low:
        return "primary"
    return "other"


def load_signatures(signature_dir, top_n):
    ranked = pd.read_csv(signature_dir / "Auto_scATLAS_mp_gene_ranked.csv")
    ranked["rank"] = ranked["rank"].astype(int)
    ranked = ranked[ranked["rank"] <= top_n].copy()

    state_groups_df = pd.read_csv(signature_dir / "Auto_scATLAS_state_groups.csv")
    signature_summary = pd.read_csv(signature_dir / "Auto_scATLAS_mp_signature_summary.csv")
    mp_order_path = signature_dir / "Auto_scATLAS_mp_order.csv"
    if mp_order_path.exists():
        mp_order = pd.read_csv(mp_order_path).sort_values("plot_order")["mp"].astype(str).tolist()
    else:
        mp_order = list(MP_DESCRIPTIONS)
    return ranked, state_groups_df, signature_summary, mp_order


def zscore_by_sample(norm_matrix, obs, signature_map, score_mask):
    samples = obs["sample"].astype(str).to_numpy()
    score_mask = np.asarray(score_mask, dtype=bool)
    score_frames = []

    for sample in pd.unique(samples):
        idx = np.where((samples == sample) & score_mask)[0]
        if len(idx) == 0:
            continue
        dense = norm_matrix[idx].toarray().astype(np.float32)
        gene_mean = dense.mean(axis=0)
        gene_sd = dense.std(axis=0)
        gene_sd[gene_sd == 0] = 1.0
        gene_z = (dense - gene_mean) / gene_sd

        sample_scores = {}
        for mp_name, gene_idx in signature_map.items():
            if len(gene_idx) == 0:
                sample_scores[mp_name] = np.full(len(idx), np.nan, dtype=np.float32)
            else:
                sample_scores[mp_name] = gene_z[:, gene_idx].mean(axis=1).astype(np.float32)

        sample_df = pd.DataFrame(sample_scores, index=obs.index[idx])
        score_frames.append(sample_df)

    if len(score_frames) == 0:
        raise ValueError("No tumor spots were available for scATLAS scoring.")

    tumor_index = obs.index[score_mask]
    mp_raw_tumor = pd.concat(score_frames, axis=0).loc[tumor_index]
    mp_raw = pd.DataFrame(np.nan, index=obs.index, columns=mp_raw_tumor.columns)
    mp_raw.loc[tumor_index, :] = mp_raw_tumor
    mp_adj = mp_raw.copy()

    for sample in pd.unique(samples):
        sample_mask = (obs["sample"].astype(str) == sample) & score_mask
        if sample_mask.sum() == 0:
            continue
        mp_adj.loc[sample_mask] = mp_adj.loc[sample_mask] - mp_adj.loc[sample_mask].mean(axis=0)

    global_sd = mp_adj.loc[score_mask].std(axis=0, ddof=0)
    global_sd[global_sd == 0] = 1.0
    mp_adj = mp_adj.divide(global_sd, axis=1)
    return mp_raw, mp_adj


def assign_states(mp_adj_noncc, threshold, hybrid_gap):
    group_scores = {}
    for state_name, mps in STATE_GROUPS.items():
        available = [mp for mp in mps if mp in mp_adj_noncc.columns]
        group_scores[state_name] = mp_adj_noncc[available].max(axis=1) if available else 0.0

    group_scores = pd.DataFrame(group_scores, index=mp_adj_noncc.index)
    best_state = group_scores.idxmax(axis=1)
    best_value = group_scores.max(axis=1)

    ordered = np.sort(group_scores.to_numpy(), axis=1)
    top1 = ordered[:, -1]
    top2 = ordered[:, -2] if ordered.shape[1] > 1 else np.zeros_like(top1)
    gap = top1 - top2

    state = best_state.astype(object)
    state[best_value < threshold] = "Unresolved"
    state[(gap < hybrid_gap) & (state != "Unresolved")] = "Hybrid"

    return group_scores, pd.Series(state, index=mp_adj_noncc.index, name="Auto_state_B"), pd.Series(gap, index=mp_adj_noncc.index, name="Auto_state_gap")


def top_mp_labels(mp_adj_noncc, threshold, mp_order):
    available_order = [mp for mp in mp_order if mp in mp_adj_noncc.columns]
    mp_adj_noncc = mp_adj_noncc[available_order]
    top_mp = mp_adj_noncc.idxmax(axis=1)
    top_val = mp_adj_noncc.max(axis=1)
    top_mp = top_mp.astype(object)
    top_mp[top_val < threshold] = "Unresolved"
    return pd.Series(top_mp, index=mp_adj_noncc.index, name="Auto_top_mp")


def label_mp(mp_name):
    if mp_name in MP_DESCRIPTIONS:
        return f"{mp_name}: {MP_DESCRIPTIONS[mp_name]}"
    return mp_name


def plot_categorical(ax, df, color_col, palette, title):
    values = df[color_col].astype(str)
    present = [cat for cat in palette if cat in set(values)]
    extras = [cat for cat in pd.unique(values) if cat not in present]
    order = present + extras

    for cat in order:
        mask = values == cat
        ax.scatter(
            df.loc[mask, "pxl_col_in_fullres"],
            df.loc[mask, "pxl_row_in_fullres"],
            s=24,
            c=palette.get(cat, "#808080"),
            linewidths=0,
            label=cat,
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
            fontsize=8,
            markerscale=1.3,
        )


def plot_top_mp(ax, df, title, mp_order):
    values = df["Auto_top_mp_label"].astype(str)
    mp_labels = [label_mp(mp) for mp in mp_order if label_mp(mp) in set(values)]
    order = ["Normal/Mixed"] + mp_labels + [x for x in ["Unresolved"] if x in set(values)]
    extras = [cat for cat in pd.unique(values) if cat not in order]
    order = order + extras
    palette = {"Normal/Mixed": MP_COLORS["Normal/Mixed"], "Unresolved": MP_COLORS["Unresolved"]}
    palette.update({label_mp(mp): MP_COLORS.get(mp, "#808080") for mp in mp_order})
    for cat in order:
        mask = values == cat
        ax.scatter(
            df.loc[mask, "pxl_col_in_fullres"],
            df.loc[mask, "pxl_row_in_fullres"],
            s=24,
            c=palette.get(cat, "#808080"),
            linewidths=0,
            label=cat,
        )
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=8, markerscale=1.3)


def plot_continuous(ax, df, color_col, title, vmax):
    non_tumor = ~df["Auto_is_tumor_bin"].astype(bool)
    if non_tumor.any():
        ax.scatter(
            df.loc[non_tumor, "pxl_col_in_fullres"],
            df.loc[non_tumor, "pxl_row_in_fullres"],
            s=24,
            c=STATE_COLORS["Normal/Mixed"],
            linewidths=0,
            label="Normal/Mixed",
        )
    tumor_df = df.loc[~non_tumor].copy()
    values = tumor_df[color_col].to_numpy(dtype=float)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    sc = ax.scatter(
        tumor_df["pxl_col_in_fullres"],
        tumor_df["pxl_row_in_fullres"],
        s=24,
        c=values,
        cmap="RdBu_r",
        norm=norm,
        linewidths=0,
    )
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    if non_tumor.any():
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=8)
    return sc


def make_state_maps(results, output_dir, mp_order):
    state_pdf = output_dir / "Auto_visium_spatial_states.pdf"
    score_pdf = output_dir / "Auto_visium_spatial_group_scores.pdf"

    with PdfPages(state_pdf) as pdf:
        for sample in pd.unique(results["sample"]):
            sub = results.loc[results["sample"] == sample].copy()
            sample_title = f"{sample} ({sub['Auto_stage'].iloc[0]})"
            fig, axes = plt.subplots(1, 2, figsize=(15, 7))
            plot_categorical(axes[0], sub, "Auto_state_B", STATE_COLORS, f"{sample_title}: scATLAS state")
            plot_top_mp(axes[1], sub, f"{sample_title}: top non-CC MP", mp_order)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            fig.savefig(output_dir / f"Auto_{sample}_state_map.png", dpi=300, bbox_inches="tight")
            plt.close(fig)

    with PdfPages(score_pdf) as pdf:
        for sample in pd.unique(results["sample"]):
            sub = results.loc[results["sample"] == sample].copy()
            sample_title = f"{sample} ({sub['Auto_stage'].iloc[0]})"
            vmax = np.nanquantile(np.abs(sub.loc[sub["Auto_is_tumor_bin"], STATE_SCORE_ORDER].to_numpy()), 0.98)
            if not np.isfinite(vmax) or vmax == 0:
                vmax = 1.0

            fig, axes = plt.subplots(2, 3, figsize=(16, 10))
            axes_flat = axes.flatten()
            for i, state_name in enumerate(STATE_SCORE_ORDER):
                sc = plot_continuous(axes_flat[i], sub, state_name, f"{sample_title}: {state_name}", vmax)
                fig.colorbar(sc, ax=axes_flat[i], fraction=0.046, pad=0.04)
            axes_flat[-1].axis("off")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            fig.savefig(output_dir / f"Auto_{sample}_group_scores.png", dpi=300, bbox_inches="tight")
            plt.close(fig)


def summarise_counts(df, group_col):
    summary = (
        df.groupby([group_col, "Auto_state_B"], observed=False)
        .size()
        .rename("spots")
        .reset_index()
    )
    totals = summary.groupby(group_col)["spots"].transform("sum")
    summary["pct"] = 100.0 * summary["spots"] / totals
    return summary


def summarise_top_mp(df, group_col, mp_label_order):
    summary = (
        df.groupby([group_col, "Auto_top_mp_label"], observed=False)
        .size()
        .rename("spots")
        .reset_index()
    )
    totals = summary.groupby(group_col)["spots"].transform("sum")
    summary["pct"] = 100.0 * summary["spots"] / totals
    summary["Auto_top_mp_label"] = pd.Categorical(summary["Auto_top_mp_label"], categories=mp_label_order, ordered=True)
    return summary.sort_values([group_col, "Auto_top_mp_label"])


def make_assignment_score_matrix(sample_df):
    rows = [state for state in STATE_SCORE_ORDER + ["Unresolved", "Hybrid"] if (sample_df["Auto_state_B"] == state).any()]
    out = []
    for state in rows:
        vals = sample_df.loc[sample_df["Auto_state_B"] == state, STATE_SCORE_ORDER].mean(axis=0)
        vals.name = state
        out.append(vals)
    total = sample_df[STATE_SCORE_ORDER].mean(axis=0)
    total.name = "Total tumor bins"
    out.append(total)
    return pd.DataFrame(out)


def make_summary_plots(results, output_dir, mp_order):
    tumor_results = results.loc[results["Auto_is_tumor_bin"]].copy()
    sample_labels = (
        results[["sample", "Auto_stage"]]
        .drop_duplicates()
        .assign(sample_label=lambda x: x["sample"] + np.where(x["Auto_stage"].eq("metastasis"), " (metastasis)", ""))
    )
    sample_label_map = dict(zip(sample_labels["sample"], sample_labels["sample_label"]))

    by_sample = summarise_counts(tumor_results, "sample")
    by_sample["sample_label"] = by_sample["sample"].map(sample_label_map)
    mp_label_order = [label_mp(mp) for mp in mp_order if mp not in CC_MPS] + ["Unresolved"]
    by_top_mp = summarise_top_mp(tumor_results, "sample", mp_label_order)
    by_top_mp["sample_label"] = by_top_mp["sample"].map(sample_label_map)

    by_sample.to_csv(output_dir / "Auto_visium_state_summary_by_sample.csv", index=False)
    by_top_mp.to_csv(output_dir / "Auto_visium_top_mp_summary_by_sample.csv", index=False)

    state_means = (
        tumor_results.groupby("sample", observed=False)[STATE_SCORE_ORDER]
        .mean()
        .reset_index()
    )
    state_means.to_csv(output_dir / "Auto_visium_group_score_means_by_sample.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    sample_pivot = (
        by_sample.pivot(index="sample_label", columns="Auto_state_B", values="pct")
        .reindex(columns=STATE_ORDER)
        .drop(columns=["Normal/Mixed"], errors="ignore")
        .fillna(0.0)
    )
    sample_pivot.plot(
        kind="bar",
        stacked=True,
        color=[STATE_COLORS[state] for state in sample_pivot.columns],
        ax=axes[0],
        width=0.8,
    )
    axes[0].set_title("State proportions by sample")
    axes[0].set_ylabel("% tumor spots")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].legend(frameon=False, title="State", bbox_to_anchor=(1.02, 1), loc="upper left")

    top_mp_pivot = (
        by_top_mp.pivot(index="sample_label", columns="Auto_top_mp_label", values="pct")
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
    axes[1].set_title("Top MP abundance by sample")
    axes[1].set_ylabel("% tumor spots")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(frameon=False, title="Top MP", bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(output_dir / "Auto_visium_state_distribution_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "Auto_visium_state_distribution_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    heatmap_pdf = output_dir / "Auto_visium_state_score_heatmaps_by_sample.pdf"
    with PdfPages(heatmap_pdf) as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes_flat = axes.flatten()
        all_vals = []
        matrices = {}
        for sample in pd.unique(results["sample"]):
            sample_df = tumor_results.loc[tumor_results["sample"] == sample].copy()
            mat = make_assignment_score_matrix(sample_df)
            matrices[sample] = mat
            all_vals.append(mat.to_numpy().ravel())
        lim = np.nanquantile(np.abs(np.concatenate(all_vals)), 0.98)
        if not np.isfinite(lim) or lim == 0:
            lim = 1.0
        for i, sample in enumerate(pd.unique(results["sample"])):
            ax = axes_flat[i]
            sample_label = sample_label_map.get(sample, sample)
            sns.heatmap(
                matrices[sample],
                cmap="RdBu_r",
                center=0,
                vmin=-lim,
                vmax=lim,
                linewidths=0.4,
                linecolor="white",
                cbar=i == 2,
                ax=ax,
            )
            ax.set_title(sample_label)
            ax.set_xlabel("State marker score")
            ax.set_ylabel("Assigned state")
            ax.tick_params(axis="x", rotation=45)
            ax.tick_params(axis="y", rotation=0)
        for ax in axes_flat[len(matrices):]:
            ax.axis("off")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(output_dir / "Auto_visium_state_score_heatmaps_by_sample.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main():
    args = parse_args()
    dataset_root = Path(args.dataset_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    visium_path = dataset_root / "processed" / "Visium.h5ad"
    sample_info_path = dataset_root / "sample_info.csv"

    ranked, state_groups_df, signature_summary, mp_order = load_signatures(output_dir, args.top_n)
    signature_summary = signature_summary.copy()

    with h5py.File(visium_path, "r") as handle:
        obs = read_dataframe(handle["obs"])
        var = read_dataframe(handle["var"])
        spatial = read_spatial_table(handle["obsm"]["spatial"])
        carcinoma = read_sparse_matrix(handle["layers"]["Carcinoma"]).tocsr()

    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)
    spatial.index = spatial.index.astype(str)

    for column in ["sample", "batch", "CNV_assignment", "Smoothed label", "spot_id"]:
        obs[column] = obs[column].astype(str)

    for column in ["array_row", "array_col", "Dist. to edge"]:
        obs[column] = pd.to_numeric(obs[column], errors="coerce")

    spatial = spatial[["pxl_col_in_fullres", "pxl_row_in_fullres"]].copy()
    results = obs.join(spatial, how="left")
    results["Carcinoma"] = pd.to_numeric(results["Carcinoma"], errors="coerce")

    sample_info = pd.read_csv(sample_info_path, encoding="utf-8-sig")
    sample_info = sample_info[sample_info["library name"].str.endswith("_Visium")].copy()
    sample_info["sample"] = sample_info["library name"].str.replace("_Visium", "", regex=False)
    sample_info["Auto_stage"] = sample_info["title"].apply(derive_stage)
    sample_info = sample_info[["sample", "title", "tissue", "description", "Auto_stage"]]
    results = results.merge(sample_info, on="sample", how="left")
    results["Auto_is_tumor_bin"] = results["Smoothed label"].astype(str).eq("Tumor")

    gene_index = pd.Index(var.index)
    signature_map = {}
    overlap_rows = []

    for mp_name, mp_df in ranked.groupby("mp", sort=False):
        genes = mp_df["gene"].astype(str).tolist()
        idx = gene_index.get_indexer(genes)
        valid = idx >= 0
        signature_map[mp_name] = idx[valid].tolist()
        overlap_rows.append(
            {
                "mp": mp_name,
                "description": MP_DESCRIPTIONS.get(mp_name, mp_name),
                "n_genes_requested": len(genes),
                "n_genes_overlap": int(valid.sum()),
            }
        )

    overlap_df = pd.DataFrame(overlap_rows).sort_values("mp")
    overlap_df.to_csv(output_dir / "Auto_scATLAS_mp_overlap_summary.csv", index=False)
    state_groups_df.to_csv(output_dir / "Auto_scATLAS_state_groups_resolved.csv", index=False)
    signature_summary.to_csv(output_dir / "Auto_scATLAS_mp_signature_summary_used.csv", index=False)

    norm_matrix = normalise_log1p(carcinoma)
    mp_raw, mp_adj = zscore_by_sample(norm_matrix, results, signature_map, results["Auto_is_tumor_bin"].to_numpy())

    retained_noncc = [mp for mp in mp_adj.columns if mp not in CC_MPS]
    mp_adj_noncc = mp_adj[retained_noncc].copy()
    tumor_index = results.index[results["Auto_is_tumor_bin"]]
    group_scores_tumor, state_assignment_tumor, state_gap_tumor = assign_states(
        mp_adj_noncc.loc[tumor_index],
        args.threshold,
        args.hybrid_gap,
    )
    top_mp_tumor = top_mp_labels(mp_adj_noncc.loc[tumor_index], args.threshold, mp_order)
    group_scores = pd.DataFrame(np.nan, index=results.index, columns=group_scores_tumor.columns)
    group_scores.loc[tumor_index, :] = group_scores_tumor
    state_assignment = pd.Series("Normal/Mixed", index=results.index, name="Auto_state_B", dtype=object)
    state_assignment.loc[tumor_index] = state_assignment_tumor
    state_gap = pd.Series(np.nan, index=results.index, name="Auto_state_gap", dtype=float)
    state_gap.loc[tumor_index] = state_gap_tumor
    top_mp = pd.Series("Normal/Mixed", index=results.index, name="Auto_top_mp", dtype=object)
    top_mp.loc[tumor_index] = top_mp_tumor

    results = pd.concat([results, mp_raw.add_prefix("Auto_raw_"), mp_adj.add_prefix("Auto_adj_"), group_scores, state_assignment, state_gap, top_mp], axis=1)

    results["Auto_top_mp_label"] = results["Auto_top_mp"].map(label_mp)

    results.to_csv(output_dir / "Auto_visium_spot_annotations.csv.gz", index=True, compression="gzip")

    make_state_maps(results, output_dir, mp_order)
    make_summary_plots(results, output_dir, mp_order)

    with open(output_dir / "Auto_visium_state_mapping_parameters.txt", "w", encoding="utf-8") as handle:
        handle.write(f"dataset_root={dataset_root}\n")
        handle.write(f"visium_path={visium_path}\n")
        handle.write(f"top_n={args.top_n}\n")
        handle.write(f"threshold={args.threshold}\n")
        handle.write(f"hybrid_gap={args.hybrid_gap}\n")
        handle.write("score_subset=Smoothed label == Tumor\n")
        handle.write("non_tumor_plot_label=Normal/Mixed\n")
        handle.write("platform=standard 10x Visium spot data, not Visium HD 8um/16um bins\n")
        handle.write("python_env=/rds/general/user/sg3723/home/miniforge3/envs/bidcell_temp\n")


if __name__ == "__main__":
    main()
