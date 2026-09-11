#!/usr/bin/env python
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_visiumhd_manual_annotation_threshold_calibration.py
#   Methodology: not required (legacy spatial utility)
#   Map: analysis/ANALYSIS_MAP.md
####################
####################
# Analysis registry:
#   Status: active
#   Description: Select one shared marker-enrichment threshold set for manual
#     segmented and custom 16 um singlet annotation.
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: ref_outs/visium_hd_outs/tables/Auto_<sample>_{custom,segmented}_cell_annotations.csv.gz;
#     ref_outs/visium_hd_outs/tables/Auto_<sample>_{custom,segmented}_cluster_gene_evidence.csv
#   Outputs: ref_outs/visium_hd_outs/tables/Auto_visiumhd_manual_threshold_*.csv;
#     updates/new_updates/summaries/visiumhd_manual_threshold_calibration_summary.csv
#   Cache/replot: Uses compact annotation/evidence tables; no count matrices are loaded.
#   Run: qsub analysis/spatial/visiumhd_manual_annotation_threshold_calibration.sh
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


REFERENCE_TYPES = {"endothelial", "macrophage", "fibroblast"}
RESIDUAL_TYPES = {"epithelial", "fibroblast"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-names", nargs="+", required=True)
    parser.add_argument("--spatial-max-distance", type=float, default=25.0)
    return parser.parse_args()


def as_bool(values):
    return values.astype(str).str.lower().isin({"true", "t", "1"})


####################
def canonical_cluster_id(values):
    # Custom tables contain NA for RCTD non-singlets, which makes pandas parse
    # otherwise integer Leiden IDs as floats (for example 6 -> 6.0).
    return values.astype(str).str.replace(r"\.0$", "", regex=True)
####################


def broad_label(values):
    values = values.fillna("unresolved").astype(str)
    return values.replace({"t.cell": "t_nk.cell", "nk.cell": "t_nk.cell"})


def cluster_labels(gene_evidence, module_scores, fallback_scores, min_log2fc, min_pct, min_pct_delta, fdr, min_cluster_cells=20):
    evidence = gene_evidence.copy()
    evidence["supported_grid"] = (
        (evidence["cluster_n"] >= min_cluster_cells)
        & (evidence["log2fc"] >= min_log2fc)
        & (evidence["pct_cluster"] >= min_pct)
        & (evidence["pct_delta"] >= min_pct_delta)
        & (evidence["p_adjusted"] <= fdr)
    )
    type_rows = []
    for (cluster, cell_type), group in evidence.groupby(["Auto_manual_cluster", "cell_type"], sort=False):
        supported = group.loc[group["supported_grid"]]
        n_available = len(group)
        ####################
        # Match the production rule: sparse one/two-gene panels require one
        # supported marker, while larger panels require two independent genes.
        required = 1 if n_available <= 2 else 2
        ####################
        type_rows.append({
            "Auto_manual_cluster": str(cluster),
            "cell_type": cell_type,
            "n_supported": len(supported),
            "fraction_supported": len(supported) / n_available,
            "median_log2fc": float(supported["log2fc"].median()) if len(supported) else np.nan,
            "passes": len(supported) >= required,
            "residual_module_score": module_scores.get((str(cluster), cell_type), np.nan),
            "residual_n_prevalent": int((group["pct_cluster"] >= 0.10).sum()),
            "residual_fraction_prevalent": float((group["pct_cluster"] >= 0.10).mean()),
            "fallback_raw_marker_score": fallback_scores.get((str(cluster), cell_type), np.nan),
        })
    type_table = pd.DataFrame(type_rows)
    labels = {}
    for cluster, group in type_table.groupby("Auto_manual_cluster", sort=False):
        passing = group.loc[group["passes"]]
        first_pass = passing.loc[~passing["cell_type"].isin(RESIDUAL_TYPES)]
        if not first_pass.empty:
            candidates = first_pass
            stage = "first"
        else:
            fibroblast_specific = passing.loc[passing["cell_type"].eq("fibroblast")]
            if not fibroblast_specific.empty:
                candidates = fibroblast_specific
                stage = "fibroblast_specific"
            else:
                candidates = group.loc[
                    group["cell_type"].eq("epithelial")
                    & (group["residual_module_score"] > 0)
                    & (group["residual_n_prevalent"] >= 2)
                ]
                stage = "epithelial_absolute"
        if candidates.empty:
            protected = group.loc[
                ~group["cell_type"].isin(RESIDUAL_TYPES)
                & (group["fallback_raw_marker_score"] > 0.10)
            ]
            if not protected.empty:
                candidates = protected
                stage = "protected_fallback"
            else:
                candidates = group.loc[group["cell_type"].isin(RESIDUAL_TYPES)]
                stage = "residual_forced_fallback"
        if stage in {"first", "fibroblast_specific"}:
            selected = candidates.sort_values(
                ["n_supported", "fraction_supported", "median_log2fc", "cell_type"],
                ascending=[False, False, False, True],
            ).iloc[0]
        elif stage == "epithelial_absolute":
            selected = candidates.sort_values(
                ["residual_fraction_prevalent", "residual_module_score", "cell_type"],
                ascending=[False, False, True],
            ).iloc[0]
        else:
            selected = candidates.sort_values(
                ["fallback_raw_marker_score", "cell_type"], ascending=[False, True]
            ).iloc[0]
        labels[str(cluster)] = selected["cell_type"]
    return pd.Series(labels, dtype=object)


def agreement_metrics(truth, prediction, min_class_n=20):
    truth = broad_label(pd.Series(truth).reset_index(drop=True))
    prediction = broad_label(pd.Series(prediction).reset_index(drop=True))
    valid_truth = truth.ne("unresolved") & truth.ne("nan") & truth.ne("")
    truth = truth.loc[valid_truth].reset_index(drop=True)
    prediction = prediction.loc[valid_truth].reset_index(drop=True)
    if truth.empty:
        return {"n_validation": 0, "accuracy": np.nan, "macro_recall": np.nan}
    classes = truth.value_counts()
    classes = classes.index[classes >= min_class_n]
    recalls = [float((prediction.loc[truth.eq(label)] == label).mean()) for label in classes]
    return {
        "n_validation": len(truth),
        "accuracy": float((truth == prediction).mean()),
        "macro_recall": float(np.mean(recalls)) if recalls else np.nan,
    }


def load_inputs(annotation_dir, sample_names):
    annotation_dir = Path(annotation_dir)
    annotations = {}
    evidence = {}
    for sample in sample_names:
        for mode in ["custom", "segmented"]:
            annotation_path = annotation_dir / f"Auto_{sample}_{mode}_cell_annotations.csv.gz"
            evidence_path = annotation_dir / f"Auto_{sample}_{mode}_cluster_gene_evidence.csv"
            if not annotation_path.exists() or not evidence_path.exists():
                raise FileNotFoundError(f"Missing annotation/evidence input for {sample} {mode}")
            annotations[(sample, mode)] = pd.read_csv(annotation_path, compression="gzip", low_memory=False)
            evidence[(sample, mode)] = pd.read_csv(evidence_path, low_memory=False)
    return annotations, evidence


def spatial_pair_indices(custom, segmented, max_distance):
    custom_xy = custom[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    segmented_xy = segmented[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    valid_custom = np.isfinite(custom_xy).all(axis=1)
    valid_segmented = np.isfinite(segmented_xy).all(axis=1)
    tree = cKDTree(segmented_xy[valid_segmented])
    distance, nearest = tree.query(custom_xy[valid_custom], k=1)
    retained = distance <= max_distance
    custom_index = np.flatnonzero(valid_custom)[retained]
    segmented_index = np.flatnonzero(valid_segmented)[nearest[retained]]
    return custom_index, segmented_index


def main():
    args = parse_args()
    annotation_dir = Path(args.annotation_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    annotations, evidence = load_inputs(annotation_dir, args.sample_names)
    pair_indices = {}
    module_scores = {}
    fallback_scores = {}
    for sample in args.sample_names:
        custom = annotations[(sample, "custom")]
        segmented = annotations[(sample, "segmented")]
        custom_manual = custom.loc[as_bool(custom["Auto_manual_qc_pass"])].reset_index(drop=True)
        annotations[(sample, "custom_manual")] = custom_manual
        pair_indices[sample] = spatial_pair_indices(custom_manual, segmented, args.spatial_max_distance)
        for mode, annotation in [("custom", custom_manual), ("segmented", segmented)]:
            cluster_ids = canonical_cluster_id(annotation["Auto_manual_cluster"])
            for cell_type in RESIDUAL_TYPES:
                score_col = f"{cell_type}_score"
                means = annotation.assign(Auto_cluster_id_calibration=cluster_ids).groupby(
                    "Auto_cluster_id_calibration", observed=True
                )[score_col].mean()
                for cluster, value in means.items():
                    module_scores[(sample, mode, str(cluster), cell_type)] = float(value)
            for score_col in [column for column in annotation.columns if column.endswith("_score2")]:
                cell_type = score_col[:-7]
                means = annotation.assign(Auto_cluster_id_calibration=cluster_ids).groupby(
                    "Auto_cluster_id_calibration", observed=True
                )[score_col].mean()
                for cluster, value in means.items():
                    fallback_scores[(sample, mode, str(cluster), cell_type)] = float(value)

    grid = list(itertools.product([0.25, 0.5, 0.75, 1.0], [0.02, 0.05, 0.10], [0.01, 0.02, 0.05]))
    metric_rows = []
    reference_rows = []
    label_cache = {}
    for min_log2fc, min_pct, min_pct_delta in grid:
        grid_id = f"lfc{min_log2fc:g}_pct{min_pct:g}_delta{min_pct_delta:g}"
        for sample in args.sample_names:
            for mode in ["custom", "segmented"]:
                mode_scores = {
                    (cluster, cell_type): value
                    for (score_sample, score_mode, cluster, cell_type), value in module_scores.items()
                    if score_sample == sample and score_mode == mode
                }
                mode_fallback_scores = {
                    (cluster, cell_type): value
                    for (score_sample, score_mode, cluster, cell_type), value in fallback_scores.items()
                    if score_sample == sample and score_mode == mode
                }
                labels = cluster_labels(
                    evidence[(sample, mode)], mode_scores, mode_fallback_scores,
                    min_log2fc, min_pct, min_pct_delta, 0.05
                )
                label_cache[(grid_id, sample, mode)] = labels
                annotation_key = (sample, "custom_manual") if mode == "custom" else (sample, "segmented")
                annotation = annotations[annotation_key]
                assigned = canonical_cluster_id(annotation["Auto_manual_cluster"]).map(labels).fillna("unresolved")
                counts = assigned.value_counts()
                reference_counts = {cell_type: int(counts.get(cell_type, 0)) for cell_type in REFERENCE_TYPES}
                reference_rows.append({
                    "grid_id": grid_id,
                    "marker_min_log2fc": min_log2fc,
                    "marker_min_pct": min_pct,
                    "marker_min_pct_delta": min_pct_delta,
                    "sample": sample,
                    "mode": mode,
                    "n_resolved": int(assigned.ne("unresolved").sum()),
                    "n_total": len(assigned),
                    "resolved_fraction": float(assigned.ne("unresolved").mean()),
                    "n_reference_types_ge20": sum(value >= 20 for value in reference_counts.values()),
                    **{f"n_{key}": value for key, value in reference_counts.items()},
                })

            custom = annotations[(sample, "custom_manual")]
            custom_labels = canonical_cluster_id(custom["Auto_manual_cluster"]).map(
                label_cache[(grid_id, sample, "custom")]
            ).fillna("unresolved")
            rctd_metrics = agreement_metrics(custom["Auto_rctd_first_type"], custom_labels)
            custom_index, segmented_index = pair_indices[sample]
            segmented = annotations[(sample, "segmented")]
            segmented_labels = canonical_cluster_id(segmented["Auto_manual_cluster"]).map(
                label_cache[(grid_id, sample, "segmented")]
            ).fillna("unresolved")
            spatial_metrics = agreement_metrics(
                segmented_labels.iloc[segmented_index], custom_labels.iloc[custom_index]
            )
            metric_rows.append({
                "grid_id": grid_id,
                "marker_min_log2fc": min_log2fc,
                "marker_min_pct": min_pct,
                "marker_min_pct_delta": min_pct_delta,
                "sample": sample,
                "rctd_n": rctd_metrics["n_validation"],
                "rctd_accuracy": rctd_metrics["accuracy"],
                "rctd_macro_recall": rctd_metrics["macro_recall"],
                "spatial_n": spatial_metrics["n_validation"],
                "spatial_accuracy": spatial_metrics["accuracy"],
                "spatial_macro_recall": spatial_metrics["macro_recall"],
            })

    metrics = pd.DataFrame(metric_rows)
    references = pd.DataFrame(reference_rows)
    pooled = metrics.groupby(["grid_id", "marker_min_log2fc", "marker_min_pct", "marker_min_pct_delta"], as_index=False).agg(
        mean_rctd_accuracy=("rctd_accuracy", "mean"),
        mean_rctd_macro_recall=("rctd_macro_recall", "mean"),
        mean_spatial_accuracy=("spatial_accuracy", "mean"),
        mean_spatial_macro_recall=("spatial_macro_recall", "mean"),
        min_spatial_macro_recall=("spatial_macro_recall", "min"),
    )
    reference_guard = references.groupby("grid_id", as_index=False).agg(
        min_reference_types_ge20=("n_reference_types_ge20", "min"),
        mean_resolved_fraction=("resolved_fraction", "mean"),
    )
    pooled = pooled.merge(reference_guard, on="grid_id", how="left")
    pooled["passes_reference_guard"] = pooled["min_reference_types_ge20"] >= 2
    pooled["validation_score"] = (
        0.6 * pooled["mean_spatial_macro_recall"] + 0.4 * pooled["mean_rctd_macro_recall"]
    )
    eligible = pooled.loc[pooled["passes_reference_guard"]].copy()
    if eligible.empty:
        raise RuntimeError("No shared threshold combination retains two reference types in every sample and mode")
    eligible = eligible.sort_values(
        ["validation_score", "min_spatial_macro_recall", "marker_min_log2fc", "marker_min_pct", "marker_min_pct_delta"],
        ascending=[False, False, False, False, False],
    )
    selected = eligible.iloc[[0]].copy()
    selected["marker_fdr"] = 0.05
    selected["marker_min_cluster_cells"] = 20
    selected["selection_rule"] = "max_0.6_spatial_macro_plus_0.4_rctd_macro_subject_to_two_reference_types_per_sample_mode"

    metrics.to_csv(output_dir / "Auto_visiumhd_manual_threshold_validation_by_sample.csv", index=False)
    references.to_csv(output_dir / "Auto_visiumhd_manual_threshold_reference_counts.csv", index=False)
    pooled.sort_values("validation_score", ascending=False).to_csv(
        output_dir / "Auto_visiumhd_manual_threshold_grid_summary.csv", index=False
    )
    selected.to_csv(output_dir / "Auto_visiumhd_manual_threshold_selected.csv", index=False)
    summary_dir = Path("updates/new_updates/summaries")
    summary_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(summary_dir / "visiumhd_manual_threshold_calibration_summary.csv", index=False)
    print(selected.to_string(index=False))


if __name__ == "__main__":
    main()
####################
