#!/usr/bin/env python

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


STATE_GROUPS = {
    "Classic Proliferative": ["MP2"],
    "Basal to Intestinal Metaplasia": ["MP17", "MP14", "MP5", "MP10", "MP8"],
    "SMG-like Metaplasia": ["MP18", "MP16"],
    "Stress-adaptive": ["MP13", "MP12"],
    "Immune Infiltrating": ["MP15"],
}

CC_MPS = ["MP1", "MP7", "MP9"]

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Temporary comparison of current Xenium scATLAS assignment versus direct UCell-style assignment."
    )
    parser.add_argument(
        "--xenium-h5ad",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Zenodo_upload/processed/Xenium.h5ad",
    )
    parser.add_argument(
        "--signature-dir",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_xenium_states",
    )
    parser.add_argument(
        "--current-annotations",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_xenium_states/Auto_xenium_cell_annotations.csv.gz",
    )
    parser.add_argument(
        "--output-dir",
        default="/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_xenium_state_method_comparison",
    )
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--max-rank", type=int, default=1500)
    parser.add_argument("--chunk-size", type=int, default=2000)
    parser.add_argument("--assignment-celltypes", default="Carcinoma")
    parser.add_argument("--control-celltypes", default="Epithelial")
    parser.add_argument(
        "--direct-gap-quantile",
        type=float,
        default=0.15,
        help="For the visual-only direct UCell filtered label, mark the lowest gap quantile as Hybrid.",
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
        return pd.Series(decode_array(node[...]))
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


def load_signatures(signature_dir, top_n, min_genes, gene_index):
    ranked = pd.read_csv(signature_dir / "Auto_scATLAS_mp_gene_ranked.csv")
    ranked["rank"] = ranked["rank"].astype(int)
    ranked = ranked[ranked["rank"] <= top_n].copy()

    signature_map = {}
    overlap_rows = []
    for mp_name, mp_df in ranked.groupby("mp", sort=False):
        genes = mp_df.sort_values("rank")["gene"].astype(str).tolist()
        idx = gene_index.get_indexer(genes)
        valid = idx >= 0
        used_idx = idx[valid].tolist()
        used_genes = [gene for gene, keep in zip(genes, valid) if keep]
        if len(used_idx) >= min_genes:
            signature_map[mp_name] = used_idx
        overlap_rows.append(
            {
                "mp": mp_name,
                "n_genes_requested_top_n": len(genes),
                "n_genes_overlap": int(valid.sum()),
                "used_for_direct_ucell": len(used_idx) >= min_genes,
                "overlap_genes": ";".join(used_genes),
            }
        )
    return signature_map, pd.DataFrame(overlap_rows)


def compute_ucell_scores(matrix, row_mask, signature_map, max_rank=1500, chunk_size=2000):
    row_mask = np.asarray(row_mask, dtype=bool)
    rows = np.where(row_mask)[0]
    mp_names = list(signature_map)
    scores = np.full((len(rows), len(mp_names)), np.nan, dtype=np.float32)
    n_genes = matrix.shape[1]
    max_rank = min(max_rank, n_genes)

    for start in range(0, len(rows), chunk_size):
        stop = min(start + chunk_size, len(rows))
        chunk_rows = rows[start:stop]
        dense = matrix[chunk_rows, :].toarray().astype(np.float32, copy=False)
        ranks = rankdata(-dense, axis=1, method="average").astype(np.float32, copy=False)

        for col_idx, mp_name in enumerate(mp_names):
            gene_idx = signature_map[mp_name]
            n_sig = len(gene_idx)
            rank_sub = ranks[:, gene_idx]
            rank_sub = np.minimum(rank_sub, max_rank)
            rank_sum = rank_sub.sum(axis=1)
            rank_sum_min = n_sig * (n_sig + 1) / 2.0
            denom = n_sig * max_rank - rank_sum_min
            scores[start:stop, col_idx] = 1.0 - (rank_sum - rank_sum_min) / denom

    return pd.DataFrame(scores, index=rows, columns=mp_names)


def state_group_scores(mp_scores):
    group_scores = {}
    for state_name, mps in STATE_GROUPS.items():
        available = [mp for mp in mps if mp in mp_scores.columns]
        if len(available) == 1:
            group_scores[state_name] = mp_scores[available[0]]
        else:
            group_scores[state_name] = mp_scores[available].max(axis=1)
    return pd.DataFrame(group_scores, index=mp_scores.index)


def top_state_from_group_scores(group_scores):
    top_state = group_scores.idxmax(axis=1).astype(object)
    top_score = group_scores.max(axis=1)
    ordered = np.sort(group_scores.to_numpy(dtype=float), axis=1)
    gap = ordered[:, -1] - ordered[:, -2]
    return top_state, top_score, pd.Series(gap, index=group_scores.index)


def summarise_states(df, state_col, out_col):
    out = df.groupby(["patient", state_col], observed=False).size().rename("cells").reset_index()
    totals = out.groupby("patient")["cells"].transform("sum")
    out["pct"] = 100.0 * out["cells"] / totals
    out = out.rename(columns={state_col: "state"})
    out["method"] = out_col
    return out


def plot_state_map(ax, df, state_col, title, point_size, rasterized):
    values = df[state_col].astype(str)
    present = [state for state in STATE_ORDER if state in set(values)]
    extras = [state for state in pd.unique(values) if state not in present]
    for state in present + extras:
        mask = values == state
        ax.scatter(
            df.loc[mask, "spatial_x"],
            df.loc[mask, "spatial_y"],
            s=point_size,
            c=STATE_COLORS.get(state, "#808080"),
            linewidths=0,
            alpha=0.92,
            rasterized=rasterized,
            label=state,
        )
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")


def add_legend(fig):
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=state, markerfacecolor=color, markersize=6)
        for state, color in STATE_COLORS.items()
        if state != "Non-carcinoma"
    ]
    fig.legend(handles=handles, loc="center right", frameon=False, fontsize=8)


def make_comparison_pdf(results, diagnostics, output_pdf):
    patients = list(pd.unique(results["patient"]))
    with PdfPages(output_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        lines = [
            "Temporary Xenium State Assignment Method Comparison",
            "",
            "Current: mean top-N signature z-score, centered per patient over carcinoma/epithelial/pre-cancerous epithelial cells,",
            "then Approach B state groups with threshold 0.5 and hybrid gap 0.3.",
            "",
            "Direct UCell top: UCell-style rank score on Xenium .X log1p(CP10K), no per-patient centering/scaling,",
            "state = max grouped non-CC MP score, no unresolved threshold.",
            "",
            "Direct UCell filtered: same direct UCell top state, but cells in the lowest direct top1-top2 gap quantile are Hybrid.",
            "",
            diagnostics.to_string(index=False),
        ]
        ax.text(0.03, 0.97, "\n".join(lines), va="top", ha="left", fontsize=9, family="monospace")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        carcinoma = results.loc[results["Auto_is_assignment_cell"]].copy()
        axes[0].hist(carcinoma["Auto_best_group_score"], bins=80, color="#4C78A8")
        axes[0].axvline(0.5, color="black", linestyle="--", linewidth=1)
        axes[0].set_title("Current best group score")
        axes[1].hist(carcinoma["Auto_direct_ucell_top_score"], bins=80, color="#F58518")
        axes[1].set_title("Direct UCell top score")
        axes[2].hist(carcinoma["Auto_direct_ucell_gap"], bins=80, color="#54A24B")
        axes[2].axvline(diagnostics.loc[diagnostics["metric"].eq("direct_gap_threshold"), "value"].iloc[0], color="black", linestyle="--", linewidth=1)
        axes[2].set_title("Direct UCell top1-top2 gap")
        for ax in axes:
            ax.set_ylabel("Cells")
            ax.set_xlabel("Score")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for patient in patients:
            sub = results.loc[results["patient"] == patient].copy()
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            plot_state_map(axes[0], sub, "Auto_state_B", f"{patient}: current normalized", 0.20, True)
            plot_state_map(axes[1], sub, "Auto_direct_ucell_state", f"{patient}: direct UCell top", 0.20, True)
            plot_state_map(axes[2], sub, "Auto_direct_ucell_state_filtered", f"{patient}: direct UCell gap-filtered", 0.20, True)
            add_legend(fig)
            fig.tight_layout(rect=(0, 0, 0.88, 1))
            pdf.savefig(fig, bbox_inches="tight", dpi=600)
            plt.close(fig)

            tumor = sub.loc[sub["Auto_is_assignment_cell"]].copy()
            if len(tumor) == 0:
                continue
            x_mid = tumor["spatial_x"].median()
            y_mid = tumor["spatial_y"].median()
            quadrants = [
                ("upper-left", tumor[(tumor["spatial_x"] <= x_mid) & (tumor["spatial_y"] <= y_mid)]),
                ("upper-right", tumor[(tumor["spatial_x"] > x_mid) & (tumor["spatial_y"] <= y_mid)]),
                ("lower-left", tumor[(tumor["spatial_x"] <= x_mid) & (tumor["spatial_y"] > y_mid)]),
                ("lower-right", tumor[(tumor["spatial_x"] > x_mid) & (tumor["spatial_y"] > y_mid)]),
            ]
            fig, axes = plt.subplots(2, 4, figsize=(18, 8))
            for j, (quad_name, quad_df) in enumerate(quadrants):
                plot_state_map(axes[0, j], quad_df, "Auto_state_B", f"{patient} {quad_name}: current", 3.0, False)
                plot_state_map(axes[1, j], quad_df, "Auto_direct_ucell_state", f"{patient} {quad_name}: direct UCell", 3.0, False)
            add_legend(fig)
            fig.tight_layout(rect=(0, 0, 0.88, 1))
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    current = pd.read_csv(args.current_annotations, compression="gzip", low_memory=False)
    current = current.rename(columns={current.columns[0]: "_index"})
    current["_index"] = current["_index"].astype(str)
    current = current.set_index("_index", drop=False)

    with h5py.File(args.xenium_h5ad, "r") as handle:
        obs = read_dataframe(handle["obs"])
        var = read_dataframe(handle["var"])
        matrix = read_sparse_matrix(handle["X"]).tocsr()

    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)
    current = current.loc[obs.index].copy()

    assignment_celltypes = split_values(args.assignment_celltypes)
    control_celltypes = split_values(args.control_celltypes)
    assignment_mask = current["Cell_type"].isin(assignment_celltypes).to_numpy()
    control_mask = current["Cell_type"].isin(control_celltypes).to_numpy()

    signature_map, overlap = load_signatures(Path(args.signature_dir), args.top_n, args.min_genes, pd.Index(var.index))
    overlap.to_csv(output_dir / "Auto_direct_ucell_signature_overlap.csv", index=False)

    score_mask = assignment_mask | control_mask
    ucell_scores_pos = compute_ucell_scores(
        matrix,
        score_mask,
        signature_map,
        max_rank=args.max_rank,
        chunk_size=args.chunk_size,
    )
    score_rows = np.where(score_mask)[0]
    ucell_scores_pos.index = obs.index[score_rows]

    non_cc_cols = [mp for mp in ucell_scores_pos.columns if mp not in CC_MPS]
    group_scores = state_group_scores(ucell_scores_pos[non_cc_cols])
    direct_state, direct_top_score, direct_gap = top_state_from_group_scores(group_scores)

    current["Auto_direct_ucell_state"] = "Non-carcinoma"
    current.loc[direct_state.index, "Auto_direct_ucell_state"] = direct_state
    current.loc[~assignment_mask, "Auto_direct_ucell_state"] = "Non-carcinoma"
    current["Auto_direct_ucell_state_filtered"] = current["Auto_direct_ucell_state"].copy()

    assignment_direct_gap = direct_gap.loc[current.index[assignment_mask]]
    direct_gap_threshold = float(assignment_direct_gap.quantile(args.direct_gap_quantile))
    low_gap_idx = assignment_direct_gap.index[assignment_direct_gap < direct_gap_threshold]
    current.loc[low_gap_idx, "Auto_direct_ucell_state_filtered"] = "Hybrid"

    current["Auto_direct_ucell_top_score"] = np.nan
    current.loc[direct_top_score.index, "Auto_direct_ucell_top_score"] = direct_top_score
    current["Auto_direct_ucell_gap"] = np.nan
    current.loc[direct_gap.index, "Auto_direct_ucell_gap"] = direct_gap

    summary = pd.concat(
        [
            summarise_states(current.loc[assignment_mask], "Auto_state_B", "current_normalized"),
            summarise_states(current.loc[assignment_mask], "Auto_direct_ucell_state", "direct_ucell_top"),
            summarise_states(current.loc[assignment_mask], "Auto_direct_ucell_state_filtered", "direct_ucell_gap_filtered"),
        ],
        axis=0,
    )
    summary.to_csv(output_dir / "Auto_xenium_state_method_comparison_summary.csv", index=False)

    pair = current.loc[assignment_mask, ["Auto_state_B", "Auto_direct_ucell_state"]].dropna()
    diagnostics = pd.DataFrame(
        [
            {"metric": "n_total_cells", "value": len(current)},
            {"metric": "n_assignment_cells", "value": int(assignment_mask.sum())},
            {"metric": "n_control_cells", "value": int(control_mask.sum())},
            {"metric": "direct_gap_quantile", "value": args.direct_gap_quantile},
            {"metric": "direct_gap_threshold", "value": direct_gap_threshold},
            {"metric": "direct_top_score_assignment_q10", "value": float(current.loc[assignment_mask, "Auto_direct_ucell_top_score"].quantile(0.10))},
            {"metric": "direct_top_score_assignment_q50", "value": float(current.loc[assignment_mask, "Auto_direct_ucell_top_score"].quantile(0.50))},
            {"metric": "direct_top_score_assignment_q90", "value": float(current.loc[assignment_mask, "Auto_direct_ucell_top_score"].quantile(0.90))},
            {"metric": "direct_top_score_control_q50", "value": float(current.loc[control_mask, "Auto_direct_ucell_top_score"].quantile(0.50))},
            {"metric": "direct_top_score_control_q90", "value": float(current.loc[control_mask, "Auto_direct_ucell_top_score"].quantile(0.90))},
            {"metric": "ari_current_vs_direct_top", "value": adjusted_rand_score(pair["Auto_state_B"], pair["Auto_direct_ucell_state"])},
            {"metric": "nmi_current_vs_direct_top", "value": normalized_mutual_info_score(pair["Auto_state_B"], pair["Auto_direct_ucell_state"])},
        ]
    )
    diagnostics.to_csv(output_dir / "Auto_xenium_state_method_comparison_diagnostics.csv", index=False)

    current[
        [
            "_index",
            "cell_id",
            "Cell_type",
            "patient",
            "spatial_x",
            "spatial_y",
            "Auto_state_B",
            "Auto_best_group_score",
            "Auto_direct_ucell_state",
            "Auto_direct_ucell_state_filtered",
            "Auto_direct_ucell_top_score",
            "Auto_direct_ucell_gap",
            "Auto_is_assignment_cell",
        ]
    ].to_csv(output_dir / "Auto_xenium_state_method_comparison_annotations.csv.gz", index=False, compression="gzip")

    make_comparison_pdf(
        current,
        diagnostics,
        output_dir / "Auto_xenium_state_method_comparison.pdf",
    )

    with open(output_dir / "Auto_xenium_state_method_comparison_parameters.txt", "w", encoding="utf-8") as handle:
        for key, value in vars(args).items():
            handle.write(f"{key}={value}\n")
        handle.write("direct_ucell_formula=UCell rank score with descending per-cell ranks, maxRank clipping, no per-patient centering/scaling\n")
        handle.write("direct_ucell_top=no unresolved threshold; direct_ucell_gap_filtered marks lowest gap quantile as Hybrid for visual comparison only\n")


if __name__ == "__main__":
    main()
