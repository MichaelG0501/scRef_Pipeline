#!/usr/bin/env python
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_binned_state_mapping.py
#   Description: Map current centred-refined scATLAS MPs and Approach-B states
#     onto level-1/level-2 malignant epithelial Visium HD 16 um bins only.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_state_mapping_methodology.md
#   Inputs:
#     analysis/spatial/visium_hd_samples.tsv
#     ref_outs/visium_hd_outs/malignancy/tables/
#       Auto_<sample>_binned_malignancy.csv.gz
#     ref_outs/visium_hd_outs/state_mapping/intermediate/signatures/
#     Space Ranger 16 um count matrices listed in the manifest
#   Outputs:
#     tables/: malignant-bin MP/state annotations, state abundance, parameters,
#       and signature coverage
#     figures/: per-sample state and top-MP spatial maps
#     logs/: run summary
#     updates/new_updates/summaries/: compact state abundance
#   Cache/replot: --replot-only rebuilds figures/summaries from the live mapped
#     annotation table without reading count matrices.
#   Run: python analysis/spatial/visium_hd_binned_state_mapping.py
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

####################
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
SAMPLES = ("SUR1231", "FFPEA1", "FFPED1")
MALIGNANT_LEVELS = ("malignant_level_1", "malignant_level_2")
CC_MPS = ("MP1", "MP5", "MP13+")
EXCLUDED_MPS: tuple[str, ...] = ()
STATE_GROUPS = {
    "Classic proliferation": ("MP2+",),
    "Basal to intestinal metaplasia": (
        "MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"
    ),
    "SMG to intestinal metaplasia": (
        "MP8+", "MP8b", "MP16", "MP18b", "MP17"
    ),
    "Stress adaptive": ("MP12",),
    "Cancer-cell immune mimicry": ("MP15",),
}
STATE_ORDER = tuple(STATE_GROUPS) + ("Unresolved", "Hybrid")
STATE_COLOURS = {
    "Classic proliferation": "#E41A1C",
    "Basal to intestinal metaplasia": "#4DAF4A",
    "SMG to intestinal metaplasia": "#FF7F00",
    "Stress adaptive": "#984EA3",
    "Cancer-cell immune mimicry": "#377EB8",
    "Unresolved": "#9E9E9E",
    "Hybrid": "#111111",
}
MP_COLOURS = {
    "MP2+": "#E41A1C", "MP14": "#1B9E77", "MP3+": "#4DAF4A",
    "MP6+": "#66C2A5", "MP11+": "#A6D854", "MP9+": "#8DA0CB",
    "MP10+": "#B3DE69", "MP8+": "#FF7F00", "MP8b": "#FC8D62",
    "MP16": "#FFD92F", "MP18b": "#E6AB02", "MP17": "#A6761D",
    "MP12": "#984EA3", "MP15": "#377EB8",
    "Unresolved": "#9E9E9E",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="",
    )
    parser.add_argument(
        "--malignancy-dir",
        default=str(WD / "ref_outs" / "visium_hd_outs" / "malignancy"),
    )
    parser.add_argument(
        "--signature-dir",
        default=str(
            WD / "ref_outs" / "visium_hd_outs" / "state_mapping"
            / "intermediate" / "signatures"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(WD / "ref_outs" / "visium_hd_outs" / "state_mapping"),
    )
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--hybrid-gap", type=float, default=0.3)
    parser.add_argument("--replot-only", action="store_true")
    return parser.parse_args()


def bool_values(values):
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).to_numpy(dtype=bool)
    return values.astype(str).str.lower().isin({"true", "t", "1"}).to_numpy()


def read_counts(input_dir):
    input_dir = Path(input_dir)
    for name in ("filtered_feature_bc_matrix.h5", "filtered_feature_cell_matrix.h5"):
        path = input_dir / name
        if path.exists():
            adata = sc.read_10x_h5(str(path), gex_only=True)
            adata.var_names_make_unique()
            return adata
    for name in ("filtered_feature_bc_matrix", "filtered_feature_cell_matrix"):
        path = input_dir / name
        if path.exists():
            adata = sc.read_10x_mtx(str(path), gex_only=True)
            adata.var_names_make_unique()
            return adata
    raise FileNotFoundError(f"No filtered expression matrix under {input_dir}")


def normalise_log1p(matrix):
    row_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    scale = np.zeros_like(row_sums)
    valid = row_sums > 0
    scale[valid] = 1e4 / row_sums[valid]
    normalised = matrix.multiply(scale[:, None]).tocsr()
    normalised.data = np.log1p(normalised.data)
    return normalised


def load_signatures(signature_dir, top_n):
    signature_dir = Path(signature_dir)
    ranked = pd.read_csv(
        signature_dir / "Auto_scATLAS_centred_refined_mp_gene_ranked.csv"
    )
    ranked = ranked.loc[ranked["rank"].astype(int) <= top_n].copy()
    order_table = pd.read_csv(
        signature_dir / "Auto_scATLAS_centred_refined_mp_order.csv"
    ).sort_values("plot_order")
    state_table = pd.read_csv(
        signature_dir / "Auto_scATLAS_centred_refined_state_groups.csv"
    )
    exported_groups = {
        state: tuple(group["mp"].astype(str))
        for state, group in state_table.groupby("state", sort=False)
    }
    if exported_groups != STATE_GROUPS:
        raise ValueError("Exported state groups do not match the production mapping constants")
    descriptions = dict(zip(order_table["mp"], order_table["description"]))
    mp_order = order_table["mp"].astype(str).tolist()
    return ranked, mp_order, descriptions


def load_malignant_inputs(manifest_path, malignancy_dir):
    manifest = pd.DataFrame({
        'sample': ['SUR1122', 'SUR1231', 'FFPEA1', 'FFPED1'],
        'binned_input': [
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/binned_outputs/square_016um',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/binned_outputs/square_016um',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/binned_outputs/square_016um',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/binned_outputs/square_016um'
        ],
        'segmented_input': [
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1122/outs/segmented_outputs',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/Frozen_batch/spaceranger_count/OCT-SUR1231/outs/segmented_outputs',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/A1/outs/segmented_outputs',
            '/rds/general/project/spatialtranscriptomics/live/Visium_HD/FFPE_batch/spaceranger_count/D1/outs/segmented_outputs'
        ]
    }).set_index("sample")
    loaded = []
    for sample in SAMPLES:
        if sample not in manifest.index:
            raise ValueError(f"Manifest lacks {sample}")
        table_path = (
            Path(malignancy_dir) / "tables"
            / f"Auto_{sample}_binned_malignancy.csv.gz"
        )
        columns = [
            "barcode", "sample", "is_epithelial_target", "Auto_malignancy",
            "Auto_malignancy_evidence", "Auto_cna_class",
            "Auto_cancer_signature_score", "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        table = pd.read_csv(
            table_path, compression="gzip", usecols=columns, low_memory=False
        )
        table["barcode"] = table["barcode"].astype(str)
        keep = bool_values(table["is_epithelial_target"]) & table[
            "Auto_malignancy"
        ].astype(str).isin(MALIGNANT_LEVELS).to_numpy()
        table = table.loc[keep].copy().set_index("barcode", drop=False)
        if table.empty:
            raise ValueError(f"No level-1/level-2 malignant epithelial bins for {sample}")
        if table.index.has_duplicates:
            raise ValueError(f"Duplicate malignant barcodes for {sample}")

        adata = read_counts(manifest.loc[sample, "binned_input"])
        missing = table.index.difference(adata.obs_names)
        if len(missing):
            raise ValueError(f"{sample}: {len(missing)} malignant bins absent from count matrix")
        adata = adata[table.index].copy()
        loaded.append((sample, adata, table.loc[adata.obs_names].copy()))
    return loaded


def score_mps(normalised, obs, signature_map):
    samples = obs["sample"].astype(str).to_numpy()
    frames = []
    for sample in SAMPLES:
        indices = np.where(samples == sample)[0]
        sample_matrix = normalised[indices]
        gene_mean = np.asarray(sample_matrix.mean(axis=0)).ravel()
        squared = sample_matrix.copy()
        squared.data **= 2
        gene_variance = np.asarray(squared.mean(axis=0)).ravel() - gene_mean**2
        gene_sd = np.sqrt(np.maximum(gene_variance, 0))
        gene_sd[gene_sd == 0] = 1
        sample_scores = {}
        for mp_name, gene_indices in signature_map.items():
            if not gene_indices:
                raise ValueError(f"No signature genes available for {mp_name}")
            inv_sd = 1.0 / gene_sd[gene_indices]
            scaled = sample_matrix[:, gene_indices].dot(sp.diags(inv_sd))
            centred_mean = np.mean(gene_mean[gene_indices] * inv_sd)
            sample_scores[mp_name] = (
                np.asarray(scaled.mean(axis=1)).ravel() - centred_mean
            ).astype(np.float32)
        frames.append(pd.DataFrame(sample_scores, index=obs.index[indices]))

    raw_scores = pd.concat(frames).loc[obs.index]
    adjusted = raw_scores.copy()
    for sample in SAMPLES:
        mask = obs["sample"].astype(str).eq(sample)
        adjusted.loc[mask] = adjusted.loc[mask] - adjusted.loc[mask].mean(axis=0)
    pooled_sd = adjusted.std(axis=0, ddof=0).replace(0, 1)
    return raw_scores, adjusted.divide(pooled_sd, axis=1)


def assign_states(adjusted, threshold, hybrid_gap):
    group_scores = pd.DataFrame(index=adjusted.index)
    for state, mps in STATE_GROUPS.items():
        group_scores[state] = adjusted[list(mps)].max(axis=1)
    state = group_scores.idxmax(axis=1).astype(object)
    top_value = group_scores.max(axis=1)
    ordered = np.sort(group_scores.to_numpy(), axis=1)
    gap = ordered[:, -1] - ordered[:, -2]
    state.loc[top_value < threshold] = "Unresolved"
    state.loc[(gap < hybrid_gap) & state.ne("Unresolved")] = "Hybrid"
    state_defining_mps = [mp for mps in STATE_GROUPS.values() for mp in mps]
    top_mp = adjusted[state_defining_mps].idxmax(axis=1).astype(object)
    top_mp.loc[adjusted[state_defining_mps].max(axis=1) < threshold] = "Unresolved"
    return group_scores, state.rename("Auto_state_B"), pd.Series(
        gap, index=adjusted.index, name="Auto_state_gap"
    ), top_mp.rename("Auto_top_mp")


def state_summary(results):
    counts = (
        results.groupby(["sample", "Auto_state_B"], observed=True)
        .size().rename("n_bins").reset_index()
    )
    counts["pct_malignant_bins"] = 100 * counts["n_bins"] / counts.groupby(
        "sample"
    )["n_bins"].transform("sum")
    counts["state_order"] = counts["Auto_state_B"].map(
        {state: index + 1 for index, state in enumerate(STATE_ORDER)}
    )
    return counts.sort_values(["sample", "state_order"])


def plot_category_map(
    data, category, palette, order, title, output_stem, background_data=None
):
    data = data.loc[
        np.isfinite(data["pxl_col_in_fullres"])
        & np.isfinite(data["pxl_row_in_fullres"])
    ].copy()
    present = [label for label in order if label in set(data[category].astype(str))]
    counts = data[category].astype(str).value_counts()
    technical = [label for label in ("Unresolved", "Hybrid") if label in present]
    draw_order = technical + [label for label in present if label not in technical]
    point_size = max(1.8, min(4.8, 50000 / max(len(data), 1)))
    figure = plt.figure(figsize=(16, 9), constrained_layout=True)
    grid = figure.add_gridspec(1, 2, width_ratios=(3.6, 2.4))
    axis = figure.add_subplot(grid[0, 0])
    legend_axis = figure.add_subplot(grid[0, 1])
    if background_data is not None and not background_data.empty:
        background_data = background_data.loc[
            np.isfinite(background_data["pxl_col_in_fullres"])
            & np.isfinite(background_data["pxl_row_in_fullres"])
        ]
        axis.scatter(
            background_data["pxl_col_in_fullres"],
            background_data["pxl_row_in_fullres"],
            s=max(0.8, point_size * 0.55), alpha=0.58, linewidths=0,
            color="#C7C7C7", rasterized=True,
        )
    for label in draw_order:
        selected = data[category].astype(str).eq(label)
        axis.scatter(
            data.loc[selected, "pxl_col_in_fullres"],
            data.loc[selected, "pxl_row_in_fullres"],
            s=point_size,
            alpha=0.78,
            linewidths=0,
            color=palette.get(label, "#666666"),
            rasterized=True,
        )
    axis.invert_yaxis()
    axis.set_aspect("equal", adjustable="datalim")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    handles = []
    if background_data is not None and not background_data.empty:
        handles.append(Line2D(
            [0], [0], marker="o", linestyle="", markersize=11,
            markerfacecolor="#C7C7C7", markeredgecolor="none",
            label=f"Other retained bins (n={len(background_data):,})",
        ))
    handles.extend([
        Line2D(
            [0], [0], marker="o", linestyle="", markersize=11,
            markerfacecolor=palette.get(label, "#666666"),
            markeredgecolor="none", label=f"{label} (n={counts[label]:,})",
        )
        for label in present
    ])
    legend_axis.axis("off")
    legend_axis.legend(
        handles=handles, loc="center left",
        frameon=False, fontsize=11, title=category.replace("Auto_", "").replace("_", " "),
        title_fontsize=12,
    )
    figure.suptitle(title, fontsize=16)
    figure.savefig(output_stem.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def write_outputs_and_plots(
    results, output_dir, descriptions, parameters, malignancy_dir
):
    output_dir = Path(output_dir)
    table_dir = output_dir / "tables"
    figure_dir = output_dir / "figures"
    log_dir = output_dir / "logs"
    for directory in (table_dir, figure_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    result_path = table_dir / "Auto_visiumhd_binned_malignant_state_annotations.csv.gz"
    results.to_csv(result_path, index=False, compression="gzip")
    summary = state_summary(results)
    summary.to_csv(
        table_dir / "Auto_visiumhd_binned_malignant_state_abundance.csv",
        index=False,
    )
    summary_dir = WD / "updates" / "new_updates" / "summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_dir / "visium_hd_binned_malignant_state_abundance.csv", index=False)
    pd.DataFrame(parameters.items(), columns=["parameter", "value"]).to_csv(
        table_dir / "Auto_visiumhd_binned_state_mapping_parameters.csv",
        index=False,
    )
    for sample in SAMPLES:
        sample_data = results.loc[results["sample"].eq(sample)].copy()
        background_path = (
            Path(malignancy_dir) / "tables"
            / f"Auto_{sample}_binned_malignancy.csv.gz"
        )
        background = pd.read_csv(
            background_path,
            compression="gzip",
            usecols=[
                "barcode", "Auto_postfilter_keep",
                "pxl_col_in_fullres", "pxl_row_in_fullres",
            ],
            low_memory=False,
        )
        retained = bool_values(background["Auto_postfilter_keep"])
        background = background.loc[
            retained & ~background["barcode"].astype(str).isin(
                sample_data["barcode"].astype(str)
            )
        ].copy()
        plot_category_map(
            sample_data, "Auto_state_B", STATE_COLOURS, STATE_ORDER,
            f"{sample} | malignant epithelial scATLAS state",
            figure_dir / f"Auto_{sample}_binned_malignant_scatlas_state_map",
            background_data=background,
        )
        mp_order = [mp for mps in STATE_GROUPS.values() for mp in mps] + ["Unresolved"]
        mp_labels = {mp: f"{mp}: {descriptions.get(mp, mp)}" for mp in mp_order if mp != "Unresolved"}
        mp_labels["Unresolved"] = "Unresolved"
        sample_data["Auto_top_mp_label"] = sample_data["Auto_top_mp"].map(mp_labels)
        label_palette = {mp_labels[key]: value for key, value in MP_COLOURS.items() if key in mp_labels}
        plot_category_map(
            sample_data, "Auto_top_mp_label", label_palette,
            [mp_labels[mp] for mp in mp_order],
            f"{sample} | malignant epithelial top scATLAS MP",
            figure_dir / f"Auto_{sample}_binned_malignant_scatlas_top_mp_map",
            background_data=background,
        )
    return summary


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    table_dir = output_dir / "tables"
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    ranked, mp_order, descriptions = load_signatures(args.signature_dir, args.top_n)
    parameters = {
        "samples": ";".join(SAMPLES),
        "representation": "binned_16um",
        "malignancy_levels": ";".join(MALIGNANT_LEVELS),
        "epithelial_gate": "is_epithelial_target_TRUE",
        "top_n_genes_per_mp": args.top_n,
        "normalisation": "log1p_CP10K_then_gene_z_within_sample",
        "mp_scale": "sample_center_then_pooled_spatial_sd",
        "state_threshold": args.threshold,
        "hybrid_gap": args.hybrid_gap,
        "cc_mps_excluded_from_states": ";".join(CC_MPS),
        "excluded_mps_not_state_defining": ";".join(EXCLUDED_MPS),
    }

    mapped_path = table_dir / "Auto_visiumhd_binned_malignant_state_annotations.csv.gz"
    if args.replot_only:
        if not mapped_path.exists():
            raise FileNotFoundError(f"Missing mapped state table: {mapped_path}")
        results = pd.read_csv(mapped_path, compression="gzip", low_memory=False)
        summary = write_outputs_and_plots(
            results, output_dir, descriptions, parameters, args.malignancy_dir
        )
        cache_status = "replot_only"
    else:
        loaded = load_malignant_inputs(args.manifest, args.malignancy_dir)
        common_genes = loaded[0][1].var_names
        for _, adata, _ in loaded[1:]:
            common_genes = common_genes.intersection(adata.var_names, sort=False)
        if len(common_genes) < 1000:
            raise ValueError(f"Only {len(common_genes)} common genes remain")

        matrices, metadata = [], []
        for sample, adata, table in loaded:
            matrices.append(adata[:, common_genes].X)
            table = table.copy()
            table.index = sample + "_" + table.index.astype(str)
            metadata.append(table)
        obs = pd.concat(metadata)
        matrix = sp.vstack(matrices).tocsr()
        normalised = normalise_log1p(matrix)
        gene_index = pd.Index(common_genes)
        signature_map = {}
        coverage_rows = []
        for mp, group in ranked.groupby("mp", sort=False):
            requested = group["gene"].astype(str).tolist()
            indices = gene_index.get_indexer(requested)
            available = indices[indices >= 0].tolist()
            signature_map[str(mp)] = available
            coverage_rows.append({
                "mp": mp, "n_requested_top_genes": len(requested),
                "n_available_common_genes": len(available),
                "available_fraction": len(available) / max(len(requested), 1),
            })
        required_state_mps = [mp for mps in STATE_GROUPS.values() for mp in mps]
        missing_signature = [mp for mp in required_state_mps if not signature_map.get(mp)]
        if missing_signature:
            raise ValueError(f"No available genes for state MP(s): {missing_signature}")

        raw_scores, adjusted = score_mps(normalised, obs, signature_map)
        group_scores, state, gap, top_mp = assign_states(
            adjusted, args.threshold, args.hybrid_gap
        )
        results = pd.concat(
            [
                obs,
                raw_scores.add_prefix("Auto_raw_"),
                adjusted.add_prefix("Auto_adj_"),
                group_scores.add_prefix("Auto_group_"),
                state,
                gap,
                top_mp,
            ],
            axis=1,
        ).reset_index(drop=True)
        table_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(coverage_rows).to_csv(
            table_dir / "Auto_visiumhd_binned_state_signature_coverage.csv",
            index=False,
        )
        summary = write_outputs_and_plots(
            results, output_dir, descriptions, parameters, args.malignancy_dir
        )
        cache_status = "rebuilt"

    lines = [
        f"status=complete", f"cache_status={cache_status}",
        f"samples={';'.join(SAMPLES)}", f"n_malignant_bins={len(results)}",
        f"state_counts={';'.join(f'{row.Auto_state_B}:{row.n_bins}' for row in summary.itertuples())}",
    ]
    (log_dir / "Auto_visiumhd_binned_state_mapping_run_summary.txt").write_text(
        "\n".join(lines) + "\n"
    )


if __name__ == "__main__":
    main()
####################
