#!/usr/bin/env python
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_spatial_spillover_annotation.py
#   Methodology: not required (legacy spatial utility)
#   Map: analysis/ANALYSIS_MAP.md
####################
####################
# Analysis registry:
#   Status: active
#   Description: Annotate segmented Visium HD cells after calibrated spatial
#     marker-transcript spillover subtraction and global Leiden clustering.
#   Methodology: analysis/methodology/spatial/legacy_visiumhd_spatial_spillover_annotation_methodology.md
#   Inputs: Space Ranger segmented_outputs directories supplied by --inputs.
#   Outputs: ref_outs/visium_hd_outs/tables/Auto_<sample>_spatial_*.csv[.gz],
#     including raw-versus-corrected annotation transitions;
#     ref_outs/visium_hd_outs/intermediate/Auto_<sample>_spatial_*.npz;
#     ref_outs/visium_hd_outs/figures/Auto_<sample>_spatial_spillover_diagnostics.pdf;
#     updates/new_updates/summaries/visiumhd_spatial_annotation_summary.csv.
#   Cache/replot: final annotation and calibration tables are persistent; the
#     sparse selected kernel is saved for audit and downstream replotting.
#   Run: qsub -v SCREF_RUN_MODE=spatial_annotation analysis/spatial/run_visium_hd_states.sh
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

from process_visium_hd import MANUAL_MARKERS, read_10x_counts, read_spatial_coordinates


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--sample-names", nargs="+", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-counts", type=int, default=200)
    parser.add_argument("--max-mt", type=float, default=15.0)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--sigma-multipliers", nargs="+", type=float, default=[0.75, 1.25, 2.0])
    parser.add_argument("--contamination-fractions", nargs="+", type=float, default=[0.05, 0.10, 0.20, 0.30])
    parser.add_argument("--radius-sigma", type=float, default=3.0)
    parser.add_argument("--bootstrap-repeats", type=int, default=50)
    args = parser.parse_args()
    if len(args.inputs) != len(args.sample_names):
        parser.error("--inputs and --sample-names must have the same length")
    return args


def marker_matrix(adata):
    genes = sorted({g for panel in MANUAL_MARKERS.values() for g in panel if g in adata.var_names})
    if not genes:
        raise ValueError("No annotation marker genes were found")
    values = adata[:, genes].X
    values = values.toarray() if sp.issparse(values) else np.asarray(values)
    return genes, values.astype(np.float32, copy=False)


def median_nearest_distance(xy):
    model = NearestNeighbors(n_neighbors=2, algorithm="kd_tree").fit(xy)
    distances = model.kneighbors(return_distance=True)[0][:, 1]
    distances = distances[np.isfinite(distances) & (distances > 0)]
    if not len(distances):
        raise ValueError("Could not estimate a positive centroid nearest-neighbour distance")
    return float(np.median(distances))


def gaussian_donor_kernel(xy, sigma, radius_sigma):
    radius = sigma * radius_sigma
    graph = NearestNeighbors(radius=radius, algorithm="kd_tree").fit(xy).radius_neighbors_graph(
        xy, mode="distance"
    ).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    graph.data = np.exp(-0.5 * (graph.data / sigma) ** 2)
    donor_mass = np.asarray(graph.sum(axis=1)).ravel()
    valid = donor_mass > 0
    scale = np.zeros_like(donor_mass)
    scale[valid] = 1.0 / donor_mass[valid]
    return sp.diags(scale).dot(graph).tocsr(), radius, int((~valid).sum())


def log_cp10k_marker_expression(counts, totals):
    scale = np.zeros_like(totals, dtype=np.float64)
    valid = totals > 0
    scale[valid] = 1e4 / totals[valid]
    return np.log1p(counts * scale[:, None]).astype(np.float32)


def type_scores(gene_expression, genes):
    gene_index = {gene: i for i, gene in enumerate(genes)}
    values = {}
    for cell_type, panel in MANUAL_MARKERS.items():
        idx = [gene_index[g] for g in panel if g in gene_index]
        values[cell_type] = gene_expression[:, idx].mean(axis=1) if idx else np.zeros(gene_expression.shape[0])
    return pd.DataFrame(values)


def safe_spearman(x, y):
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 30 or np.nanstd(x[valid]) == 0 or np.nanstd(y[valid]) == 0:
        return 0.0
    value = spearmanr(x[valid], y[valid]).statistic
    return float(value) if np.isfinite(value) else 0.0


def calibration_metrics(raw_counts, expected_counts, corrected_counts, totals, genes, clusters):
    raw_scores = type_scores(log_cp10k_marker_expression(raw_counts, totals), genes)
    background_scores = type_scores(log_cp10k_marker_expression(expected_counts, totals), genes)
    corrected_scores = type_scores(log_cp10k_marker_expression(corrected_counts, totals), genes)
    leakage = []
    retention = []
    zero_increase = []
    for cell_type in corrected_scores.columns:
        raw = raw_scores[cell_type].to_numpy()
        bg = background_scores[cell_type].to_numpy()
        corrected = corrected_scores[cell_type].to_numpy()
        positive = raw > 0
        low_cut = np.quantile(raw[positive], 0.75) if positive.any() else 0
        low = positive & (raw <= low_cut) & (bg > 0)
        if low.sum() < 30:
            low = positive & (bg > 0)
        leakage.append(abs(safe_spearman(corrected[low], bg[low])))
        high_cut = np.quantile(raw[positive], 0.90) if positive.any() else 0
        high = positive & (raw >= high_cut)
        if high.sum() and raw[high].mean() > 0:
            retention.append(float(corrected[high].mean() / raw[high].mean()))
        zero_increase.append(float(np.mean(corrected <= 0) - np.mean(raw <= 0)))
    z = (corrected_scores - corrected_scores.mean(axis=0)) / corrected_scores.std(axis=0, ddof=0).replace(0, 1)
    z["cluster"] = clusters
    cluster_z = z.groupby("cluster", observed=True).mean()
    ordered = np.sort(cluster_z.to_numpy(), axis=1)
    margins = ordered[:, -1] - ordered[:, -2]
    median_retention = float(np.nanmedian(retention)) if retention else 1.0
    residual_leakage = float(np.nanmedian(leakage)) if leakage else 0.0
    median_zero_increase = float(np.nanmedian(zero_increase)) if zero_increase else 0.0
    cluster_margin = float(np.nanmedian(margins))
    objective = (
        residual_leakage
        + 2.0 * max(0.0, 0.75 - median_retention)
        + 1.5 * max(0.0, median_zero_increase - 0.20)
        - 0.10 * cluster_margin
    )
    return {
        "residual_background_abs_spearman": residual_leakage,
        "high_source_signal_retention": median_retention,
        "zero_fraction_increase": median_zero_increase,
        "cluster_top_second_margin": cluster_margin,
        "calibration_objective": objective,
    }


def cluster_evidence(scores, corrected_expression, genes, clusters):
    gene_index = {gene: i for i, gene in enumerate(genes)}
    score_sd = scores.std(axis=0, ddof=0).replace(0, 1)
    zscores = (scores - scores.mean(axis=0)) / score_sd
    zscores["cluster"] = clusters
    cluster_z = zscores.groupby("cluster", observed=True).mean()
    rows = []
    clusters_array = np.asarray(clusters).astype(str)
    for cluster in cluster_z.index.astype(str):
        inside = clusters_array == cluster
        outside = ~inside
        for cell_type, panel in MANUAL_MARKERS.items():
            present = [g for g in panel if g in gene_index]
            supported = []
            for gene in present:
                value = corrected_expression[:, gene_index[gene]]
                pct_in = np.mean(value[inside] > 0)
                pct_out = np.mean(value[outside] > 0) if outside.any() else 0
                mean_in = np.mean(value[inside])
                mean_out = np.mean(value[outside]) if outside.any() else 0
                if pct_in >= 0.05 and pct_in >= pct_out + 0.02 and mean_in > mean_out:
                    supported.append(gene)
            required = min(2, len(present))
            rows.append({
                "Auto_spatial_cluster": cluster,
                "cell_type": cell_type,
                "cluster_standardised_score": float(cluster_z.loc[cluster, cell_type]),
                "cluster_corrected_marker_score": float(scores.loc[inside, cell_type].mean()),
                "supported_markers": ";".join(supported),
                "n_supported_markers": len(supported),
                "required_markers": required,
                "passes_marker_support": len(supported) >= required,
            })
    return pd.DataFrame(rows)


def label_clusters(evidence, minimum_z, ambiguity_gap):
    labels = {}
    details = []
    for cluster, group in evidence.groupby("Auto_spatial_cluster", sort=False):
        eligible = group[(group["cluster_corrected_marker_score"] > 0) &
                         (group["cluster_standardised_score"] >= minimum_z)].copy()
        if eligible.empty:
            label = "unresolved"
            selected = eligible
        else:
            top_index = eligible["cluster_standardised_score"].idxmax()
            top = eligible.loc[top_index, "cluster_standardised_score"]
            selected = eligible[
                (eligible["cluster_standardised_score"] >= top - ambiguity_gap)
                & (eligible["passes_marker_support"] | (eligible.index == top_index))
            ].sort_values(
                "cluster_standardised_score", ascending=False
            )
            label = "|".join(selected["cell_type"].tolist())
        labels[str(cluster)] = label
        details.append({
            "Auto_spatial_cluster": str(cluster),
            "Auto_annotation_celltype": label,
            "Auto_annotation_minimum_z": minimum_z,
            "Auto_annotation_ambiguity_gap": ambiguity_gap,
            "Auto_annotation_selected_scores": ";".join(
                f"{r.cell_type}={r.cluster_standardised_score:.3f}" for r in selected.itertuples()
            ),
        })
    return labels, pd.DataFrame(details)


def calibrate_annotation_thresholds(scores, corrected_expression, genes, clusters, repeats):
    evidence = cluster_evidence(scores, corrected_expression, genes, clusters)
    rng = np.random.default_rng(20260722)
    cluster_values = np.asarray(clusters).astype(str)
    grid_rows = []
    for minimum_z in [0.0, 0.25, 0.5]:
        for gap in [0.15, 0.30, 0.50]:
            labels, _ = label_clusters(evidence, minimum_z, gap)
            stability_values = []
            for cluster in sorted(set(cluster_values)):
                idx = np.where(cluster_values == cluster)[0]
                original = set(labels[cluster].split("|"))
                if original == {"unresolved"}:
                    continue
                for _ in range(repeats):
                    sampled = rng.choice(idx, size=len(idx), replace=True)
                    bootstrap = scores.iloc[sampled].mean(axis=0)
                    supported = evidence[evidence["Auto_spatial_cluster"] == cluster]
                    z = (bootstrap - scores.mean(axis=0)) / scores.std(axis=0, ddof=0).replace(0, 1)
                    eligible = supported[supported["cell_type"].map(z) >= minimum_z].copy()
                    if eligible.empty:
                        selected = {"unresolved"}
                    else:
                        eligible["bootstrap_z"] = eligible["cell_type"].map(z)
                        top_index = eligible["bootstrap_z"].idxmax()
                        top = eligible.loc[top_index, "bootstrap_z"]
                        selected = set(eligible.loc[
                            (eligible["bootstrap_z"] >= top - gap)
                            & (eligible["passes_marker_support"] | (eligible.index == top_index)),
                            "cell_type",
                        ])
                    stability_values.append(len(original & selected) / len(original | selected))
            resolved = np.array([labels[c] != "unresolved" for c in cluster_values])
            multiplicity = np.mean([len(labels[c].split("|")) for c in cluster_values[resolved]]) if resolved.any() else 0
            stability = float(np.mean(stability_values)) if stability_values else 0
            coverage = float(np.mean(resolved))
            objective = stability + 0.75 * coverage - 0.05 * max(0, multiplicity - 1)
            grid_rows.append({
                "minimum_z": minimum_z,
                "ambiguity_gap": gap,
                "bootstrap_label_stability": stability,
                "resolved_cell_fraction": coverage,
                "mean_labels_per_resolved_cell": multiplicity,
                "threshold_objective": objective,
            })
    grid = pd.DataFrame(grid_rows).sort_values(
        ["threshold_objective", "minimum_z", "ambiguity_gap"], ascending=[False, False, True]
    )
    selected = grid.iloc[0]
    labels, label_details = label_clusters(evidence, selected.minimum_z, selected.ambiguity_gap)
    return labels, label_details, evidence, grid


def prepare_clusters(adata, resolution):
    work = adata.copy()
    sc.pp.filter_genes(work, min_cells=10)
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=min(3000, work.n_vars), span=0.3)
    cluster = work[:, work.var["highly_variable"]].copy()
    sc.pp.regress_out(cluster, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(cluster, max_value=10)
    n_comps = min(50, cluster.n_vars - 1, cluster.n_obs - 1)
    sc.tl.pca(cluster, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(cluster, n_neighbors=min(15, cluster.n_obs - 1), n_pcs=min(40, n_comps), metric="cosine")
    sc.tl.leiden(cluster, resolution=resolution, key_added="Auto_spatial_cluster", flavor="igraph", directed=False)
    sc.tl.umap(cluster, random_state=0)
    return cluster.obs["Auto_spatial_cluster"].astype(str).to_numpy(), cluster.obsm["X_umap"]


def annotation_comparison_page(pdf, annotation, x, y, page_title, invert_y=False):
    label_columns = [
        "Auto_annotation_celltype_uncorrected",
        "Auto_annotation_celltype",
    ]
    labels = sorted(set().union(*[
        set(annotation[column].dropna().astype(str)) for column in label_columns
    ]))
    colours = {
        label: plt.cm.tab20(index % 20)
        for index, label in enumerate(labels)
    }
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for axis, column, title in zip(axes, label_columns, ["Uncorrected", "Corrected"]):
        point_colours = annotation[column].astype(str).map(colours)
        axis.scatter(
            annotation[x],
            annotation[y],
            c=list(point_colours),
            s=2,
            rasterized=True,
        )
        if invert_y:
            axis.invert_yaxis()
        axis.set_title(title, fontsize=16, fontweight="bold")
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", markersize=8, color=colours[label], label=label)
        for label in labels
    ]
    fig.suptitle(page_title, fontsize=18, fontweight="bold")
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(5, max(1, len(handles))),
        frameon=False,
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.10, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def diagnostics_pdf(path, annotation, calibration, selected, cell_types):
    with PdfPages(path) as pdf:
        annotation_comparison_page(
            pdf,
            annotation,
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
            "Spatial annotation",
            invert_y=True,
        )
        annotation_comparison_page(
            pdf,
            annotation,
            "Auto_umap_1",
            "Auto_umap_2",
            "Whole-transcriptome UMAP annotation",
        )

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        pivot = calibration.pivot_table(index="sigma_multiplier", columns="contamination_fraction",
                                        values="residual_background_abs_spearman")
        im = axes[0].imshow(pivot, aspect="auto", cmap="viridis_r")
        axes[0].set_title("Residual neighbour dependence")
        axes[0].set_xticks(range(len(pivot.columns)), pivot.columns)
        axes[0].set_yticks(range(len(pivot.index)), pivot.index)
        axes[0].set_xlabel("Contamination fraction")
        axes[0].set_ylabel("Sigma / median NN distance")
        fig.colorbar(im, ax=axes[0], fraction=0.046)
        axes[1].scatter(calibration["high_source_signal_retention"], calibration["residual_background_abs_spearman"],
                        s=55, c=calibration["calibration_objective"], cmap="viridis_r")
        axes[1].axvline(0.75, color="black", linestyle="--")
        axes[1].set_xlabel("High-source signal retained")
        axes[1].set_ylabel("Residual neighbour dependence")
        axes[2].barh(cell_types.index.astype(str), cell_types.to_numpy())
        axes[2].set_xlabel("Cells")
        axes[2].set_title("Final annotations")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def annotate_sample(input_path, sample_name, args):
    output = Path(args.output_dir)
    for subdir in ["tables", "intermediate", "figures"]:
        (output / subdir).mkdir(parents=True, exist_ok=True)
    adata = read_10x_counts(Path(input_path))
    coordinates = read_spatial_coordinates(Path(input_path), "spatial", adata.obs_names)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    keep = (adata.obs["total_counts"] >= args.min_counts) & (adata.obs["pct_counts_mt"] <= args.max_mt)
    adata = adata[keep].copy()
    coordinates = coordinates.loc[adata.obs_names].copy()
    clusters, umap = prepare_clusters(adata, args.leiden_resolution)
    genes, raw_counts = marker_matrix(adata)
    totals = adata.obs["total_counts"].to_numpy(dtype=float)
    xy = coordinates[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    nn_distance = median_nearest_distance(xy)
    calibration_rows = []
    candidates = {}
    for sigma_multiplier in args.sigma_multipliers:
        sigma = nn_distance * sigma_multiplier
        kernel, radius, isolated = gaussian_donor_kernel(xy, sigma, args.radius_sigma)
        expected = np.asarray(kernel.T.dot(raw_counts), dtype=np.float32)
        for rho in args.contamination_fractions:
            corrected = np.maximum(0, raw_counts - rho * expected)
            metrics = calibration_metrics(raw_counts, rho * expected, corrected, totals, genes, clusters)
            row = {
                "sample": sample_name, "sigma_multiplier": sigma_multiplier,
                "sigma_coordinate_units": sigma, "radius_coordinate_units": radius,
                "contamination_fraction": rho, "median_nearest_neighbour_distance": nn_distance,
                "n_isolated_cells": isolated, **metrics,
            }
            calibration_rows.append(row)
            candidates[(sigma_multiplier, rho)] = (kernel, rho * expected, corrected)
    calibration = pd.DataFrame(calibration_rows).sort_values(
        ["calibration_objective", "sigma_multiplier", "contamination_fraction"]
    )
    best = calibration.iloc[0]
    kernel, expected, corrected = candidates[(best.sigma_multiplier, best.contamination_fraction)]
    corrected_expression = log_cp10k_marker_expression(corrected, totals)
    scores = type_scores(corrected_expression, genes)
    labels, label_details, evidence, threshold_grid = calibrate_annotation_thresholds(
        scores, corrected_expression, genes, clusters, args.bootstrap_repeats
    )
    ####################
    # Apply the corrected pipeline's selected z threshold and ambiguity gap to
    # uncorrected marker scores. Clusters, QC, marker-support rules, and every
    # annotation parameter are held fixed, isolating the effect of subtraction.
    raw_expression = log_cp10k_marker_expression(raw_counts, totals)
    raw_scores = type_scores(raw_expression, genes)
    raw_evidence = cluster_evidence(raw_scores, raw_expression, genes, clusters)
    selected_threshold = threshold_grid.iloc[0]
    raw_labels, raw_label_details = label_clusters(
        raw_evidence,
        selected_threshold.minimum_z,
        selected_threshold.ambiguity_gap,
    )
    ####################
    annotation = adata.obs.copy()
    annotation.insert(0, "barcode", annotation.index.astype(str))
    annotation["sample"] = sample_name
    annotation["Auto_spatial_cluster"] = clusters
    annotation["Auto_umap_1"] = umap[:, 0]
    annotation["Auto_umap_2"] = umap[:, 1]
    annotation["Auto_annotation_celltype_uncorrected"] = pd.Series(
        clusters, index=annotation.index
    ).map(raw_labels)
    annotation["Auto_annotation_celltype"] = pd.Series(clusters, index=annotation.index).map(labels)
    annotation["Auto_annotation_method"] = "spatial_spillover_corrected_global_cluster"
    annotation["Auto_annotation_pass_doublet_filter"] = True
    annotation["Auto_annotation_doublet_status"] = "not_applicable_segmented_cells"
    annotation["Auto_annotation_keep_epithelial"] = annotation["Auto_annotation_celltype"].eq("epithelial")
    annotation["Auto_spatial_sigma_multiplier"] = best.sigma_multiplier
    annotation["Auto_spatial_sigma_coordinate_units"] = best.sigma_coordinate_units
    annotation["Auto_spatial_contamination_fraction"] = best.contamination_fraction
    annotation = pd.concat([annotation, coordinates], axis=1)
    background_scores = type_scores(log_cp10k_marker_expression(expected, totals), genes)
    for cell_type in scores.columns:
        annotation[f"{cell_type}_score"] = scores[cell_type].to_numpy()
        annotation[f"{cell_type}_score2"] = raw_scores[cell_type].to_numpy()
        annotation[f"{cell_type}_expected_background"] = background_scores[cell_type].to_numpy()
    prefix = f"Auto_{sample_name}_spatial"
    annotation.to_csv(output / "tables" / f"{prefix}_cell_annotations.csv.gz", index=False, compression="gzip")
    calibration.to_csv(output / "tables" / f"{prefix}_spillover_calibration.csv", index=False)
    threshold_grid.insert(0, "sample", sample_name)
    threshold_grid.to_csv(output / "tables" / f"{prefix}_annotation_threshold_calibration.csv", index=False)
    evidence = evidence.merge(label_details, on="Auto_spatial_cluster", how="left")
    evidence.insert(0, "sample", sample_name)
    evidence.to_csv(output / "tables" / f"{prefix}_cluster_marker_evidence.csv", index=False)
    raw_evidence = raw_evidence.merge(
        raw_label_details, on="Auto_spatial_cluster", how="left"
    )
    raw_evidence.insert(0, "sample", sample_name)
    raw_evidence.to_csv(
        output / "tables" / f"{prefix}_uncorrected_cluster_marker_evidence.csv",
        index=False,
    )
    transitions = (
        annotation.groupby(
            ["Auto_annotation_celltype_uncorrected", "Auto_annotation_celltype"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    transitions.insert(0, "sample", sample_name)
    transitions.to_csv(
        output / "tables" / f"{prefix}_uncorrected_corrected_transitions.csv",
        index=False,
    )
    sp.save_npz(output / "intermediate" / f"{prefix}_selected_donor_kernel.npz", kernel)
    counts = annotation["Auto_annotation_celltype"].value_counts()
    diagnostics_pdf(output / "figures" / f"{prefix}_spillover_diagnostics.pdf", annotation, calibration, best, counts)
    changed = annotation["Auto_annotation_celltype_uncorrected"].ne(
        annotation["Auto_annotation_celltype"]
    )
    return {
        "sample": sample_name,
        "mode": "spatial",
        "n_annotated": len(annotation),
        "n_epithelial": int(annotation["Auto_annotation_keep_epithelial"].sum()),
        "n_unresolved": int(annotation["Auto_annotation_celltype"].eq("unresolved").sum()),
        "n_annotation_changed_by_correction": int(changed.sum()),
        "pct_annotation_changed_by_correction": float(100 * changed.mean()),
        "sigma_multiplier": best.sigma_multiplier,
        "contamination_fraction": best.contamination_fraction,
        "high_source_signal_retention": best.high_source_signal_retention,
        "residual_background_abs_spearman": best.residual_background_abs_spearman,
        "celltype_counts": ";".join(f"{k}={v}" for k, v in counts.items()),
    }


def main():
    args = parse_args()
    summaries = [annotate_sample(path, sample, args) for path, sample in zip(args.inputs, args.sample_names)]
    summary = pd.DataFrame(summaries)
    output = Path(args.output_dir)
    summary.to_csv(output / "tables" / "Auto_visiumhd_spatial_annotation_summary.csv", index=False)
    update_dir = Path("updates/new_updates/summaries")
    update_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(update_dir / "visiumhd_spatial_annotation_summary.csv", index=False)


if __name__ == "__main__":
    main()
