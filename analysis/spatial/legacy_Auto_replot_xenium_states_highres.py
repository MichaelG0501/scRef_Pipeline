#!/usr/bin/env python
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_Auto_replot_xenium_states_highres.py
#   Description: Replots cached historical uncentred Xenium state assignments.
#     It must not be used as a current centred-state analytical input.
#   Methodology: analysis/methodology/spatial/spatial_mapping_methodology.md
####################

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


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

BIO_STATE_ORDER = [
    "Classic Proliferative",
    "Basal to Intestinal Metaplasia",
    "SMG-like Metaplasia",
    "Stress-adaptive",
    "Immune Infiltrating",
]

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
        description="Replot cached Xenium scATLAS state assignments at higher resolution and run neighbourhood analysis."
    )
    parser.add_argument(
        "--annotations",
        default="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_Yates2025_EAC_spatial/Auto_scATLAS_xenium_states/Auto_xenium_cell_annotations.csv.gz",
    )
    parser.add_argument(
        "--output-dir",
        default="/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_xenium_state_replot_highres",
    )
    parser.add_argument("--knn-k", type=int, default=50)
    parser.add_argument("--dense-knn-k", type=int, default=100)
    parser.add_argument("--dense-target-cells", type=int, default=9000)
    return parser.parse_args()


def parse_bool(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "t", "1", "yes"])


def ordered_labels(values):
    seen = set(values.astype(str))
    ordered = [state for state in STATE_ORDER if state in seen]
    ordered += [state for state in pd.unique(values.astype(str)) if state not in ordered]
    return ordered


def point_style(label, n_points, mode):
    if mode == "context":
        return {
            "Non-carcinoma": (0.15, 0.18, 1),
            "Unresolved": (0.15, 0.48, 2),
            "Hybrid": (0.15, 0.54, 2),
        }.get(label, (0.6 if n_points > 100000 else 1.05, 0.96, 3))
    if mode == "full":
        return {
            "Non-carcinoma": (0.225, 0.24, 1),
            "Unresolved": (0.225, 0.54, 2),
            "Hybrid": (0.225, 0.60, 2),
        }.get(label, (0.85 if n_points > 100000 else 1.55, 0.98, 3))
    return {
        "Non-carcinoma": (5.0, 0.28, 1),
        "Unresolved": (5.0, 0.60, 2),
        "Hybrid": (5.0, 0.68, 2),
    }.get(label, (12.5, 0.98, 3))


def plot_state_map(ax, df, title, mode="full", rasterized=True, draw_legend=True):
    values = df["Auto_state_B"].astype(str)
    for label in ordered_labels(values):
        mask = values == label
        size, alpha, zorder = point_style(label, len(df), mode)
        ax.scatter(
            df.loc[mask, "spatial_x"],
            df.loc[mask, "spatial_y"],
            s=size,
            c=STATE_COLORS.get(label, "#808080"),
            alpha=alpha,
            linewidths=0,
            rasterized=rasterized,
            zorder=zorder,
            label=label,
        )
    ax.set_title(title, fontsize=11, pad=6)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    if draw_legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            frameon=False,
            fontsize=11,
            labelspacing=1.5,
            markerscale=10 if mode != "zoom" else 3,
        )


def load_annotations(path):
    keep = [
        "_index",
        "cell_id",
        "Cell_type",
        "patient",
        "spatial_x",
        "spatial_y",
        "Auto_state_B",
        "Auto_is_assignment_cell",
        "Auto_best_group_score",
        "Auto_state_gap",
    ]
    df = pd.read_csv(path, compression="gzip", low_memory=False, usecols=lambda col: col in keep)
    df["Auto_is_assignment_cell"] = parse_bool(df["Auto_is_assignment_cell"])
    df["Auto_state_B"] = df["Auto_state_B"].astype(str)
    df["patient"] = df["patient"].astype(str)
    df["spatial_x"] = pd.to_numeric(df["spatial_x"], errors="coerce")
    df["spatial_y"] = pd.to_numeric(df["spatial_y"], errors="coerce")
    df = df.dropna(subset=["spatial_x", "spatial_y", "patient", "Auto_state_B"])
    return df


def make_full_maps(df, output_dir):
    pdf_path = output_dir / "Auto_xenium_spatial_states_highres.pdf"
    with PdfPages(pdf_path) as pdf:
        for patient in pd.unique(df["patient"]):
            sub = df.loc[df["patient"] == patient].copy()
            fig, ax = plt.subplots(figsize=(10, 9))
            plot_state_map(ax, sub, f"{patient}: scATLAS state", mode="full", rasterized=True)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight", dpi=900)
            fig.savefig(output_dir / f"Auto_{patient}_xenium_state_map_highres.png", dpi=900, bbox_inches="tight")
            plt.close(fig)
    return pdf_path


def select_dense_region(sub, dense_knn_k, target_cells):
    all_xy = sub[["spatial_x", "spatial_y"]].to_numpy(dtype=float)
    state_mask = sub["Auto_state_B"].isin(BIO_STATE_ORDER)
    candidates = sub.loc[state_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
    if len(candidates) == 0:
        candidates = sub.loc[sub["Auto_is_assignment_cell"], ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
    if len(candidates) == 0:
        candidates = all_xy
    k = min(dense_knn_k + 1, len(all_xy))
    if k <= 2:
        center = np.median(all_xy, axis=0)
    else:
        n_state_neighbors = min(dense_knn_k + 1, len(candidates))
        nn_state = NearestNeighbors(n_neighbors=n_state_neighbors, algorithm="kd_tree").fit(candidates)
        state_radius = nn_state.kneighbors(candidates, return_distance=True)[0][:, -1]
        nn_all = NearestNeighbors(n_neighbors=k, algorithm="kd_tree").fit(all_xy)
        all_radius = nn_all.kneighbors(candidates, return_distance=True)[0][:, -1]
        score = state_radius + 0.25 * all_radius
        center = candidates[int(np.argmin(score))]

    dist_to_center = np.sqrt(((all_xy - center) ** 2).sum(axis=1))
    n_target = min(max(2500, target_cells), len(all_xy))
    selected = np.argsort(dist_to_center)[:n_target]
    selected_xy = all_xy[selected, :]
    xmin, ymin = selected_xy.min(axis=0)
    xmax, ymax = selected_xy.max(axis=0)
    width = xmax - xmin
    height = ymax - ymin
    pad = max(width, height) * 0.06
    if width < height:
        extra = (height - width) / 2
        xmin -= extra
        xmax += extra
    elif height < width:
        extra = (width - height) / 2
        ymin -= extra
        ymax += extra
    return xmin - pad, xmax + pad, ymin - pad, ymax + pad, center[0], center[1]


def make_dense_region_maps(df, output_dir, dense_knn_k, target_cells):
    rows = []
    pdf_path = output_dir / "Auto_xenium_dense_subregions_highres.pdf"
    with PdfPages(pdf_path) as pdf:
        for patient in pd.unique(df["patient"]):
            sub = df.loc[df["patient"] == patient].copy()
            xmin, xmax, ymin, ymax, cx, cy = select_dense_region(sub, dense_knn_k, target_cells)
            zoom = sub.loc[
                (sub["spatial_x"] >= xmin)
                & (sub["spatial_x"] <= xmax)
                & (sub["spatial_y"] >= ymin)
                & (sub["spatial_y"] <= ymax)
            ].copy()
            rows.append(
                {
                    "patient": patient,
                    "xmin": xmin,
                    "xmax": xmax,
                    "ymin": ymin,
                    "ymax": ymax,
                    "center_x": cx,
                    "center_y": cy,
                    "cells_in_region": len(zoom),
                    "carcinoma_cells_in_region": int(zoom["Auto_is_assignment_cell"].sum()),
                    "biological_state_cells_in_region": int(zoom["Auto_state_B"].isin(BIO_STATE_ORDER).sum()),
                }
            )

            fig, axes = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.0, 3.1]})
            plot_state_map(axes[0], sub, f"{patient}: selected dense region", mode="context", rasterized=True, draw_legend=False)
            axes[0].add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, fill=False, edgecolor="#111111", linewidth=1.2))
            axes[0].scatter([cx], [cy], s=10, c="#111111", linewidths=0, zorder=5)
            plot_state_map(axes[1], zoom, f"{patient}: dense subregion ({len(zoom):,} cells)", mode="zoom", rasterized=False)
            axes[1].set_xlim(xmin, xmax)
            axes[1].set_ylim(ymax, ymin)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight", dpi=900)
            fig.savefig(output_dir / f"Auto_{patient}_xenium_dense_subregion_highres.png", dpi=900, bbox_inches="tight")
            plt.close(fig)
    region_df = pd.DataFrame(rows)
    region_df.to_csv(output_dir / "Auto_xenium_dense_subregion_bounds.csv", index=False)
    return pdf_path


def neighbourhood_scores(df, k):
    rows = []
    tumor = df.loc[df["Auto_state_B"].isin(BIO_STATE_ORDER)].copy()
    for patient, sub in tumor.groupby("patient", sort=False):
        if len(sub) <= k:
            continue
        coords = sub[["spatial_x", "spatial_y"]].to_numpy(dtype=float)
        labels = sub["Auto_state_B"].astype(str).to_numpy()
        nn = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree").fit(coords)
        ind = nn.kneighbors(coords, return_distance=False)[:, 1:]
        same = (labels[ind] == labels[:, None]).mean(axis=1)
        rows.append(
            pd.DataFrame(
                {
                    "patient": patient,
                    "cell_id": sub.get("cell_id", sub.index).astype(str).to_numpy(),
                    "state": labels,
                    "same_neighbor_score": same,
                    "k": k,
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["patient", "cell_id", "state", "same_neighbor_score", "k"])
    return pd.concat(rows, axis=0, ignore_index=True)


def grouped_boxplot(ax, plot_df, states, title):
    data = [plot_df.loc[plot_df["state"] == state, "same_neighbor_score"].to_numpy() for state in states]
    positions = np.arange(len(states)) + 1
    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True, showfliers=False)
    for patch, state in zip(bp["boxes"], states):
        patch.set_facecolor(STATE_COLORS[state])
        patch.set_alpha(0.82)
        patch.set_edgecolor("#111111")
        patch.set_linewidth(0.8)
    for element in ["whiskers", "caps", "medians"]:
        for artist in bp[element]:
            artist.set_color("#111111")
            artist.set_linewidth(0.8)
    rng = np.random.default_rng(42)
    for pos, state in zip(positions, states):
        vals = plot_df.loc[plot_df["state"] == state, "same_neighbor_score"].to_numpy()
        if len(vals) > 2500:
            vals = rng.choice(vals, size=2500, replace=False)
        x = rng.normal(pos, 0.08, size=len(vals))
        ax.scatter(x, vals, s=4, c=STATE_COLORS[state], alpha=0.28, linewidths=0)
    ax.set_xticks(positions)
    ax.set_xticklabels(states, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Same-state neighbours")
    ax.set_title(title, fontsize=12)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_neighbourhood_plot(scores, output_dir):
    scores.to_csv(output_dir / "Auto_xenium_neighbourhood_same_state_scores.csv", index=False)
    summary = (
        scores.groupby(["patient", "state"], observed=False)["same_neighbor_score"]
        .agg(["count", "mean", "median"])
        .reset_index()
    )
    summary.to_csv(output_dir / "Auto_xenium_neighbourhood_same_state_summary.csv", index=False)

    pdf_path = output_dir / "Auto_xenium_neighbourhood_colocalisation.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(7.0, 5.5))
        states = [state for state in BIO_STATE_ORDER if state in set(scores["state"])]
        grouped_boxplot(ax, scores, states, "Xenium scATLAS state colocalisation")
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight", dpi=600)
        fig.savefig(output_dir / "Auto_xenium_neighbourhood_colocalisation.png", dpi=900, bbox_inches="tight")
        plt.close(fig)

        heat = summary.pivot(index="patient", columns="state", values="mean").reindex(columns=states)
        fig, ax = plt.subplots(figsize=(8.5, 4.2))
        im = ax.imshow(heat.to_numpy(dtype=float), vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index)
        ax.set_title("Mean same-state neighbour score by patient")
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight", dpi=600)
        plt.close(fig)
    return pdf_path


def make_proportion_plot(df, output_dir):
    prop_df = df[["patient", "Auto_state_B"]].copy()
    prop_df.columns = ["patient", "state"]
    prop_df = prop_df[prop_df["state"] != "Non-carcinoma"]
    
    overall = prop_df["state"].value_counts().reset_index()
    overall.columns = ["state", "count"]
    overall["patient"] = "Total"
    overall["pct"] = 100 * overall["count"] / overall["count"].sum()
    
    per_patient = prop_df.groupby(["patient", "state"]).size().reset_index(name="count")
    per_patient["pct"] = 100 * per_patient["count"] / per_patient.groupby("patient")["count"].transform("sum")
    
    plot_df = pd.concat([per_patient, overall], ignore_index=True)
    
    state_levels = [s for s in STATE_ORDER if s in plot_df["state"].unique()] + [s for s in plot_df["state"].unique() if s not in STATE_ORDER]
    
    patient_levels = sorted([p for p in plot_df["patient"].unique() if p != "Total"]) + ["Total"]
    plot_df["is_total"] = plot_df["patient"].apply(lambda x: "Total" if x == "Total" else "Patients")
    
    fig = plt.figure(figsize=(18, 8))
    gs = fig.add_gridspec(1, 3, width_ratios=[5, 1, 2])
    
    ax1 = fig.add_subplot(gs[0])
    patients_df = plot_df[plot_df["is_total"] == "Patients"]
    pivot_patients = patients_df.pivot(index="patient", columns="state", values="pct").reindex(columns=state_levels, fill_value=0)
    
    bottom = np.zeros(len(pivot_patients))
    for state in state_levels:
        values = pivot_patients[state].values
        color = STATE_COLORS.get(state, "#808080")
        ax1.bar(pivot_patients.index, values, bottom=bottom, color=color, edgecolor="black", linewidth=0.2, label=state)
        for i, val in enumerate(values):
            if val > 3:
                ax1.text(i, bottom[i] + val/2, f"{val:.1f}%", ha='center', va='center', fontsize=12, color="black")
        bottom += values
        
    ax1.set_ylabel("% of cells", fontsize=16, weight="bold")
    ax1.set_title("State proportions (Patients)", fontsize=18)
    ax1.set_xticks(range(len(pivot_patients)))
    ax1.set_xticklabels(pivot_patients.index, rotation=45, ha="right", fontsize=14)
    ax1.set_xlim(-0.6, len(pivot_patients)-0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    ax2 = fig.add_subplot(gs[1])
    total_df = plot_df[plot_df["is_total"] == "Total"]
    pivot_total = total_df.pivot(index="patient", columns="state", values="pct").reindex(columns=state_levels, fill_value=0)
    
    bottom = np.zeros(len(pivot_total))
    for state in state_levels:
        values = pivot_total[state].values
        color = STATE_COLORS.get(state, "#808080")
        ax2.bar(pivot_total.index, values, bottom=bottom, color=color, edgecolor="black", linewidth=0.2)
        for i, val in enumerate(values):
            if val > 3:
                ax2.text(i, bottom[i] + val/2, f"{val:.1f}%", ha='center', va='center', fontsize=12, color="black")
        bottom += values
        
    ax2.set_title("Total", fontsize=18)
    ax2.set_xticks([0])
    ax2.set_xticklabels(["Total"], rotation=45, ha="right", fontsize=14)
    ax2.set_xlim(-0.6, 0.6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.set_yticks([])
    
    ax3 = fig.add_subplot(gs[2])
    pie_data = pivot_total.loc["Total"]
    pie_data = pie_data[pie_data > 0]
    colors = [STATE_COLORS.get(s, "#808080") for s in pie_data.index]
    
    wedges, texts = ax3.pie(
        pie_data, 
        colors=colors,
        wedgeprops=dict(width=1, edgecolor="white")
    )
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        pct = pie_data.iloc[i]
        if pct > 3:
            ax3.annotate(f"{pie_data.index[i]}\n{pct:.1f}%", xy=(x*0.5, y*0.5), ha="center", va="center", color="black", weight="bold", fontsize=12)
            
    ax3.set_title("Overall pie", fontsize=18, weight="bold")
    
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc='center right', title="State", bbox_to_anchor=(0.98, 0.5), fontsize=14, title_fontsize=16)
    
    fig.tight_layout()
    fig.subplots_adjust(right=0.80, wspace=0.1)
    
    pdf_path = output_dir / "Auto_xenium_state_proportions.pdf"
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight", dpi=300)
    fig.savefig(output_dir / "Auto_xenium_state_proportions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_annotations(args.annotations)
    make_proportion_plot(df, output_dir)
    make_full_maps(df, output_dir)
    make_dense_region_maps(df, output_dir, args.dense_knn_k, args.dense_target_cells)
    scores = neighbourhood_scores(df, args.knn_k)
    make_neighbourhood_plot(scores, output_dir)

    with open(output_dir / "Auto_xenium_replot_parameters.txt", "w", encoding="utf-8") as handle:
        for key, value in vars(args).items():
            handle.write(f"{key}={value}\n")
        handle.write("colors=original Xenium state palette; non-biological states are plotted smaller than biological states\n")
        handle.write("neighbourhood_strategy=Visium-compatible biological states only; Hybrid, Unresolved and Non-carcinoma excluded\n")
        handle.write("full_map_pdf_dpi=900\n")
        handle.write("dense_zoom_points=rasterized_false\n")


if __name__ == "__main__":
    main()
