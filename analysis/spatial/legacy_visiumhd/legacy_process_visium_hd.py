#!/usr/bin/env python
####################
# Analysis registry:
#   Status: legacy
#   Script: analysis/spatial/legacy_visiumhd/legacy_process_visium_hd.py
#   Methodology: not required (legacy spatial utility)
#   Map: analysis/ANALYSIS_MAP.md
####################
####################
# Analysis registry:
#   Status: active
#   Description: Annotate RCTD-binned, custom 16 um, or segmented Visium HD
#     observations and map scATLAS MPs/states after the malignancy gate.
#   Methodology: analysis/methodology/spatial/legacy_visium_hd_annotation_cnv_methodology.md
#   Inputs: Space Ranger square_016um or segmented output directories passed by
#     --inputs; ref_outs/visium_hd_outs/ exported scATLAS signature CSVs; for
#     binned/custom, ref_outs/visium_hd_outs/rctd/tables/Auto_<sample>_binned_rctd_annotations.csv.gz;
#     for mapping, ref_outs/visium_hd_outs/malignancy/tables/Auto_<sample>_<mode>_infercna_malignancy.csv.gz
#   Outputs: ref_outs/visium_hd_outs/tables/ annotation, marker-evidence,
#     threshold, state, and summary CSVs; ref_outs/visium_hd_outs/figures/ maps;
#     updates/new_updates/summaries/ compact annotation/state summaries.
#   Cache/replot: --stage annotate and --stage map are separable; map reuses the
#     persisted annotation and malignancy tables without rerunning annotation.
#   Run: qsub -v SCREF_RUN_MODE=<custom_annotation|segmented_annotation|custom|segmented_downstream>
#     analysis/spatial/run_visium_hd_states.sh
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

try:
    import geopandas as gpd
except ImportError:
    gpd = None

# New state configuration (noreg)
CC_MPS = ["MP1", "MP5", "MP13+"]

STATE_GROUPS = {
    "Classic proliferation": ["MP2+"],
    "Basal to intestinal metaplasia": ["MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"],
    "SMG to intestinal metaplasia": ["MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP2x"],
    "Stress adaptive": ["MP12"],
    "Cancer-cell immune mimicry": ["MP15"]
}

STATE_COLORS = {
    "Classic proliferation": "#E41A1C",
    "Basal to intestinal metaplasia": "#4DAF4A",
    "SMG to intestinal metaplasia": "#FF7F00",
    "Stress adaptive": "#984EA3",
    "Cancer-cell immune mimicry": "#377EB8",
    "Unresolved": "grey",
    "Hybrid": "black",
    "Normal/Mixed": "#E6E6E6",
}

MP_COLORS = {
    "MP1": "#B0B0B0",
    "MP5": "#C0C0C0",
    "MP13+": "#999999",
    "MP2+": "#E41A1C",
    "MP14": "#4DAF4A",
    "MP3+": "#66C2A5",
    "MP6+": "#A6D854",
    "MP11+": "#E78AC3",
    "MP9+": "#8DA0CB",
    "MP10+": "#FFD92F",
    "MP8+": "#FF7F00",
    "MP8b": "#FC8D62",
    "MP16": "#FFD92F",
    "MP18b": "#FF7F00",
    "MP17": "#4DAF4A",
    "MP2x": "#FF7F00",
    "MP12": "#984EA3",
    "MP15": "#377EB8",
    "Unresolved": "grey",
    "Normal/Mixed": "#E6E6E6",
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, help="Paths to square_016um/ or segmented_outputs/")
    parser.add_argument("--sample-names", nargs="+", required=True)
    parser.add_argument("--mode", choices=["binned", "segmented"], default="binned")
    parser.add_argument("--signature-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--hybrid-gap", type=float, default=0.3)
    return parser.parse_args()

def load_signatures(signature_dir, top_n):
    sig_dir = Path(signature_dir)
    ranked = pd.read_csv(sig_dir / "Auto_scATLAS_mp_gene_ranked.csv")
    ranked["rank"] = ranked["rank"].astype(int)
    ranked = ranked[ranked["rank"] <= top_n].copy()
    
    mp_order_path = sig_dir / "Auto_scATLAS_mp_order.csv"
    if mp_order_path.exists():
        mp_order_df = pd.read_csv(mp_order_path)
        mp_order = mp_order_df.sort_values("plot_order")["mp"].astype(str).tolist()
        mp_desc_map = dict(zip(mp_order_df["mp"], mp_order_df["description"]))
    else:
        mp_order = list(ranked["mp"].unique())
        mp_desc_map = {m: m for m in mp_order}
        
    return ranked, mp_order, mp_desc_map

def normalise_log1p(matrix):
    row_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    scale = np.zeros_like(row_sums)
    valid = row_sums > 0
    scale[valid] = 1e4 / row_sums[valid]
    norm = matrix.multiply(scale[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    return norm

def zscore_by_sample(norm_matrix, obs, signature_map):
    samples = obs["sample"].astype(str).to_numpy()
    score_frames = []

    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        if len(idx) == 0: continue
        
        sub_mat = norm_matrix[idx]
        gene_mean = np.asarray(sub_mat.mean(axis=0)).ravel()
        sub_mat_sq = sub_mat.copy()
        sub_mat_sq.data **= 2
        gene_var = np.asarray(sub_mat_sq.mean(axis=0)).ravel() - (gene_mean ** 2)
        gene_var[gene_var < 0] = 0
        gene_sd = np.sqrt(gene_var)
        gene_sd[gene_sd == 0] = 1.0

        sample_scores = {}
        for mp_name, gene_idx in signature_map.items():
            if len(gene_idx) == 0:
                sample_scores[mp_name] = np.full(len(idx), np.nan, dtype=np.float32)
            else:
                sub_genes = sub_mat[:, gene_idx]
                inv_sd = 1.0 / gene_sd[gene_idx]
                sub_scaled = sub_genes.dot(sp.diags(inv_sd))
                mean_scaled = gene_mean[gene_idx] / gene_sd[gene_idx]
                
                mp_val = np.asarray(sub_scaled.mean(axis=1)).ravel() - np.mean(mean_scaled)
                sample_scores[mp_name] = mp_val.astype(np.float32)

        score_frames.append(pd.DataFrame(sample_scores, index=obs.index[idx]))

    mp_raw = pd.concat(score_frames, axis=0).loc[obs.index]
    mp_adj = mp_raw.copy()

    for sample in pd.unique(samples):
        sample_mask = (obs["sample"].astype(str) == sample)
        mp_adj.loc[sample_mask] = mp_adj.loc[sample_mask] - mp_adj.loc[sample_mask].mean(axis=0)

    global_sd = mp_adj.std(axis=0, ddof=0)
    global_sd[global_sd == 0] = 1.0
    return mp_raw, mp_adj.divide(global_sd, axis=1)

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

def label_mp(mp_name, mp_desc_map):
    if mp_name in mp_desc_map:
        return f"{mp_name}: {mp_desc_map[mp_name]}"
    return mp_name

def plot_categorical(ax, df, color_col, palette, title):
    values = df[color_col].astype(str)
    
    # Bottom cats plotted first so they are under
    bottom_cats = ["Unresolved", "Hybrid", "Normal/Mixed"]
    
    present = [cat for cat in palette if cat in set(values)]
    extras = [cat for cat in pd.unique(values) if cat not in present]
    all_cats = present + extras
    
    order_bottom = [c for c in bottom_cats if c in all_cats]
    order_top = [c for c in all_cats if c not in bottom_cats]
    order = order_bottom + order_top

    for cat in order:
        mask = values == cat
        ax.scatter(
            df.loc[mask, "pxl_col_in_fullres"],
            df.loc[mask, "pxl_row_in_fullres"],
            s=1.0,
            c=palette.get(cat, "#808080"),
            linewidths=0,
            label=cat,
        )

    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    if order:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=8, markerscale=3)

def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ranked, mp_order, mp_desc_map = load_signatures(args.signature_dir, args.top_n)
    
    all_obs = []
    all_spatial = []
    matrices = []
    var_names = None
    
    for path_str, sample_name in zip(args.inputs, args.sample_names):
        path = Path(path_str)
        
        # Load count matrix
        h5_path_1 = path / "filtered_feature_bc_matrix.h5"
        h5_path_2 = path / "filtered_feature_cell_matrix.h5"
        
        if h5_path_1.exists():
            adata = sc.read_10x_h5(str(h5_path_1))
        elif h5_path_2.exists():
            adata = sc.read_10x_h5(str(h5_path_2))
        else:
            mtx_dir = path / "filtered_feature_cell_matrix"
            if not mtx_dir.exists():
                mtx_dir = path / "filtered_feature_bc_matrix"
            adata = sc.read_10x_mtx(mtx_dir)
            
        adata.var_names_make_unique()
        if var_names is None: var_names = adata.var_names
        
        umi_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["total_counts"] = umi_counts
        adata.obs["sample"] = sample_name
        adata = adata[umi_counts >= 200].copy()
        
        # Load spatial coords
        if args.mode == "binned":
            parquet_path = path / "spatial" / "tissue_positions.parquet"
            csv_path = path / "spatial" / "tissue_positions.csv"
            if parquet_path.exists():
                spatial_df = pd.read_parquet(parquet_path)
            elif csv_path.exists():
                spatial_df = pd.read_csv(csv_path)
            else:
                raise FileNotFoundError(f"No tissue_positions file found in {path / 'spatial'}")
            
            if "barcode" in spatial_df.columns:
                spatial_df = spatial_df.set_index("barcode")
            spatial_df = spatial_df.loc[adata.obs_names]
            
        else: # segmented
            geojson_path = path / "cell_segmentations.geojson"
            if not geojson_path.exists():
                raise FileNotFoundError(f"Missing {geojson_path}")
            if gpd is None:
                raise ImportError("geopandas is required for segmented mode")
            
            seg = gpd.read_file(geojson_path)
            seg["cell_id"] = seg["cell_id"].astype(str)
            seg = seg.set_index("cell_id", drop=False)
            
            # Map adata obs_names to cell_id integers
            import re
            extracted_ids = adata.obs_names.to_series().str.extract(r'(\d+)', expand=False).astype(int).astype(str)
            
            centroids = seg.loc[extracted_ids].geometry.centroid
            spatial_df = pd.DataFrame({
                "pxl_col_in_fullres": centroids.x.to_numpy(),
                "pxl_row_in_fullres": centroids.y.to_numpy()
            }, index=adata.obs_names)
            
        all_obs.append(adata.obs)
        all_spatial.append(spatial_df)
        matrices.append(adata.X)
    
    obs = pd.concat(all_obs)
    spatial = pd.concat(all_spatial)
    X = sp.vstack(matrices).tocsr()
    
    obs.index = obs["sample"] + "_" + obs.index.astype(str)
    spatial.index = obs.index
    
    norm_matrix = normalise_log1p(X)
    
    gene_index = pd.Index(var_names)
    signature_map = {mp: gene_index.get_indexer(mp_df["gene"])[gene_index.get_indexer(mp_df["gene"]) >= 0].tolist() 
                     for mp, mp_df in ranked.groupby("mp", sort=False)}
        
    mp_raw, mp_adj = zscore_by_sample(norm_matrix, obs, signature_map)
    mp_adj_noncc = mp_adj[[m for m in mp_adj.columns if m not in CC_MPS]].copy()
    
    group_scores, state_assignment, state_gap = assign_states(mp_adj_noncc, args.threshold, args.hybrid_gap)
    top_mp = top_mp_labels(mp_adj_noncc, args.threshold, mp_order)
    
    results = pd.concat([
        obs, spatial[["pxl_row_in_fullres", "pxl_col_in_fullres"]],
        mp_raw.add_prefix("Auto_raw_"), mp_adj.add_prefix("Auto_adj_"),
        group_scores, state_assignment, state_gap, top_mp
    ], axis=1)
    results["Auto_top_mp_label"] = results["Auto_top_mp"].apply(lambda x: label_mp(x, mp_desc_map))
    
    out_csv = out_dir / f"Auto_visiumhd_{args.mode}_spot_annotations.csv.gz"
    results.to_csv(out_csv, index=True, compression="gzip")
    
    for sample in args.sample_names:
        sub = results.loc[results["sample"] == sample].copy()
        fig, axes = plt.subplots(1, 2, figsize=(15, 7))
        plot_categorical(axes[0], sub, "Auto_state_B", STATE_COLORS, f"{sample} ({args.mode}): scATLAS state (noreg)")
        plot_categorical(axes[1], sub, "Auto_top_mp_label", 
                         {label_mp(m, mp_desc_map): MP_COLORS.get(m, "grey") for m in mp_order}, 
                         f"{sample} ({args.mode}): top non-CC MP")
        fig.tight_layout()
        fig.savefig(out_dir / f"Auto_{sample}_{args.mode}_state_map.png", dpi=400, bbox_inches="tight")
        with PdfPages(out_dir / f"Auto_{sample}_{args.mode}_state_map.pdf") as pdf:
            pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

####################
# Annotation and malignancy-gated Visium HD workflow. This supersedes the
# all-cell mapper above while retaining its scoring helpers and plot styling.
####################
MANUAL_MARKERS = {
    "fibroblast": ["COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"],
    "macrophage": ["CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"],
    "mast": ["MS4A2", "CPA3", "TPSB2", "TPSAB1"],
    "epithelial": ["KRT7", "MUC1", "KRT19", "EPCAM"],
    "t.cell": ["CD3E", "CD3D", "CD2", "CD3G"],
    "b.cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "nk.cell": ["GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"],
    "plasma": ["MZB1", "JCHAIN", "DERL3"],
    "dendritic": ["CLEC10A", "CCR7", "CD86"],
    "endothelial": ["ENG", "CLEC14A", "CLDN5", "VWF", "CDH5"],
    "lymph": ["CCL21"],
    "erythrocyte": ["HBA1", "HBA2", "HBB"],
    "keratinocyte": ["FLG", "IVL"],
    "neutrophil": ["CTSG", "ELANE", "MPO", "AZU1"],
}

####################
# Marker-enrichment annotation uses a one-sided detection-enrichment test and
# effect-size filters on CP10K expression. This avoids raw-UMI score thresholds,
# whose scale changes between segmented cells and 16 um bins.
from scipy.stats import hypergeom
####################


def gated_parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate Visium HD data and map scATLAS states only in malignant epithelial observations."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Paths to square_016um/ or segmented_outputs/")
    parser.add_argument("--sample-names", nargs="+", required=True)
    ####################
    # Spatial annotations are generated separately and mapped from segmented counts.
    parser.add_argument("--mode", choices=["binned", "segmented", "custom", "spatial"], required=True)
    ####################
    parser.add_argument("--stage", choices=["annotate", "map"], required=True)
    parser.add_argument("--signature-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rctd-dir", help="RCTD output directory; required for binned annotation")
    parser.add_argument("--malignancy-dir", help="InferCNA output directory; required for state mapping")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--hybrid-gap", type=float, default=0.3)
    parser.add_argument("--min-counts", type=int, default=200)
    parser.add_argument("--max-mt", type=float, default=15.0)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    ####################
    # Shared manual-annotation evidence thresholds. These apply identically to
    # segmented observations and custom RCTD-singlet bins.
    parser.add_argument("--marker-min-log2fc", type=float, default=1.0)
    parser.add_argument("--marker-min-pct", type=float, default=0.10)
    parser.add_argument("--marker-min-pct-delta", type=float, default=0.05)
    parser.add_argument("--marker-fdr", type=float, default=0.05)
    parser.add_argument("--marker-min-cluster-cells", type=int, default=20)
    parser.add_argument("--residual-min-module-score", type=float, default=0.0)
    parser.add_argument("--residual-min-pct", type=float, default=0.10)
    parser.add_argument("--residual-min-markers", type=int, default=2)
    parser.add_argument("--fallback-protected-min-score", type=float, default=0.10)
    ####################
    args = parser.parse_args()
    if len(args.inputs) != len(args.sample_names):
        parser.error("--inputs and --sample-names must have the same length")
    return args


def read_10x_counts(path):
    h5_path_1 = path / "filtered_feature_bc_matrix.h5"
    h5_path_2 = path / "filtered_feature_cell_matrix.h5"
    if h5_path_1.exists():
        adata = sc.read_10x_h5(str(h5_path_1), gex_only=True)
    elif h5_path_2.exists():
        adata = sc.read_10x_h5(str(h5_path_2), gex_only=True)
    else:
        mtx_dir = path / "filtered_feature_cell_matrix"
        if not mtx_dir.exists():
            mtx_dir = path / "filtered_feature_bc_matrix"
        adata = sc.read_10x_mtx(str(mtx_dir), gex_only=True)
    adata.var_names_make_unique()
    return adata


def read_spatial_coordinates(path, mode, obs_names):
    if mode in {"binned", "custom"}:
        parquet_path = path / "spatial" / "tissue_positions.parquet"
        csv_path = path / "spatial" / "tissue_positions.csv"
        if parquet_path.exists():
            spatial = pd.read_parquet(parquet_path)
        elif csv_path.exists():
            spatial = pd.read_csv(csv_path)
        else:
            raise FileNotFoundError(f"No tissue positions file in {path / 'spatial'}")
        barcode_col = next((c for c in ["barcode", "Barcode"] if c in spatial.columns), None)
        if barcode_col is None:
            raise ValueError(f"No barcode column in tissue positions for {path}")
        spatial = spatial.set_index(barcode_col)
        x_col = next((c for c in ["pxl_col_in_fullres", "array_col"] if c in spatial.columns), None)
        y_col = next((c for c in ["pxl_row_in_fullres", "array_row"] if c in spatial.columns), None)
        if x_col is None or y_col is None:
            raise ValueError(f"No supported coordinate columns in tissue positions for {path}")
        spatial = spatial.reindex(obs_names)
        return pd.DataFrame(
            {"pxl_col_in_fullres": spatial[x_col].to_numpy(), "pxl_row_in_fullres": spatial[y_col].to_numpy()},
            index=obs_names,
        )

    geojson_path = path / "cell_segmentations.geojson"
    if not geojson_path.exists():
        raise FileNotFoundError(f"Missing {geojson_path}")
    if gpd is None:
        raise ImportError("geopandas is required for segmented mode")
    seg = gpd.read_file(geojson_path)
    seg["cell_id"] = seg["cell_id"].astype(str)
    seg = seg.set_index("cell_id", drop=False)
    cell_ids = pd.Index(obs_names).to_series().str.extract(r"(\d+)", expand=False).astype(int).astype(str)
    if cell_ids.isna().any():
        raise ValueError(f"Could not extract numeric cell IDs for {cell_ids.isna().sum()} segmented observations")
    missing = ~cell_ids.isin(seg.index)
    if missing.any():
        raise ValueError(f"{missing.sum()} expression cells are missing from the segmentation GeoJSON")
    ####################
    # Space Ranger stores slide-coordinate polygons with geographic CRS
    # metadata. Calculate each Shapely centroid directly so GeoPandas does not
    # imply that a geographic reprojection is required; numeric coordinates
    # remain in the native slide system used by the count output.
    centroids = [geometry.centroid for geometry in seg.loc[cell_ids.to_numpy()].geometry]
    ####################
    return pd.DataFrame(
        {
            "pxl_col_in_fullres": np.asarray([centroid.x for centroid in centroids]),
            "pxl_row_in_fullres": np.asarray([centroid.y for centroid in centroids]),
        },
        index=obs_names,
    )


def raw_marker_mean(adata, genes):
    valid = [gene for gene in genes if gene in adata.var_names]
    if not valid:
        return np.full(adata.n_obs, np.nan, dtype=np.float32)
    values = adata[:, valid].X
    return np.asarray(values.mean(axis=1)).ravel().astype(np.float32)


####################
def benjamini_hochberg(p_values):
    p_values = np.asarray(p_values, dtype=float)
    if p_values.size == 0:
        return p_values
    order = np.argsort(p_values)
    ranked = p_values[order] * p_values.size / np.arange(1, p_values.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0.0, 1.0)
    return adjusted


def cluster_marker_evidence(adata, cluster_key, args):
    marker_genes = sorted({gene for genes in MANUAL_MARKERS.values() for gene in genes if gene in adata.var_names})
    if not marker_genes:
        raise ValueError("No manual marker genes are present after QC")
    marker_indices = adata.var_names.get_indexer(marker_genes)
    counts = adata.layers["counts"][:, marker_indices]
    if sp.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float64)
    totals = np.maximum(adata.obs["total_counts"].to_numpy(dtype=np.float64), 1.0)
    cp10k = counts * (1e4 / totals[:, None])
    detected = counts > 0
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    n_cells = adata.n_obs
    total_cp10k = cp10k.sum(axis=0)
    total_detected = detected.sum(axis=0)
    gene_to_index = {gene: idx for idx, gene in enumerate(marker_genes)}
    gene_rows = []
    for cluster in sorted(pd.unique(clusters)):
        in_cluster = clusters == cluster
        cluster_n = int(in_cluster.sum())
        rest_n = n_cells - cluster_n
        if rest_n <= 0:
            raise ValueError("Manual annotation requires more than one expression cluster")
        cluster_cp10k = cp10k[in_cluster].sum(axis=0)
        cluster_detected = detected[in_cluster].sum(axis=0)
        mean_cluster = cluster_cp10k / cluster_n
        mean_rest = (total_cp10k - cluster_cp10k) / rest_n
        pct_cluster = cluster_detected / cluster_n
        pct_rest = (total_detected - cluster_detected) / rest_n
        log2fc = np.log2((mean_cluster + 1.0) / (mean_rest + 1.0))
        p_values = hypergeom.sf(cluster_detected - 1, n_cells, total_detected, cluster_n)
        p_adjusted = benjamini_hochberg(p_values)
        for cell_type, genes in MANUAL_MARKERS.items():
            for gene in genes:
                if gene not in gene_to_index:
                    continue
                idx = gene_to_index[gene]
                supported = (
                    cluster_n >= args.marker_min_cluster_cells
                    and log2fc[idx] >= args.marker_min_log2fc
                    and pct_cluster[idx] >= args.marker_min_pct
                    and (pct_cluster[idx] - pct_rest[idx]) >= args.marker_min_pct_delta
                    and p_adjusted[idx] <= args.marker_fdr
                )
                gene_rows.append({
                    "Auto_manual_cluster": cluster,
                    "cell_type": cell_type,
                    "gene": gene,
                    "cluster_n": cluster_n,
                    "mean_cp10k_cluster": mean_cluster[idx],
                    "mean_cp10k_rest": mean_rest[idx],
                    "log2fc": log2fc[idx],
                    "pct_cluster": pct_cluster[idx],
                    "pct_rest": pct_rest[idx],
                    "pct_delta": pct_cluster[idx] - pct_rest[idx],
                    "p_value": p_values[idx],
                    "p_adjusted": p_adjusted[idx],
                    "supported": bool(supported),
                })
    gene_evidence = pd.DataFrame(gene_rows)
    type_rows = []
    for (cluster, cell_type), group in gene_evidence.groupby(["Auto_manual_cluster", "cell_type"], sort=False):
        supported = group.loc[group["supported"]]
        n_available = len(group)
        ####################
        # One marker is sufficient only for one- or two-gene panels (CCL21;
        # FLG/IVL), where requiring two would make dropout-driven false
        # negatives unavoidable. Panels with three or more genes require two.
        required = 1 if n_available <= 2 else 2
        ####################
        type_rows.append({
            "Auto_manual_cluster": cluster,
            "cell_type": cell_type,
            "cluster_n": int(group["cluster_n"].iloc[0]),
            "n_markers_available": n_available,
            "n_markers_supported": len(supported),
            "required_markers": required,
            "marker_fraction_supported": len(supported) / n_available,
            "marker_median_log2fc": float(supported["log2fc"].median()) if len(supported) else np.nan,
            "marker_max_padjusted": float(supported["p_adjusted"].max()) if len(supported) else np.nan,
            "supported_markers": ";".join(supported["gene"].tolist()),
            "passes_marker_evidence": len(supported) >= required,
        })
    type_evidence = pd.DataFrame(type_rows)
    ####################
    # Epithelial transcripts may be common across tumour clusters, so requiring
    # cluster-vs-rest enrichment for the final epithelial residual call creates
    # false negatives. Fibroblast remains a differential ECM-specificity call;
    # only epithelial uses module evidence and marker co-detection after the
    # protected first-stage and fibroblast-specific tests have failed.
    ####################
    cluster_score_means = adata.obs.groupby(cluster_key, observed=True)[
        ["epithelial_score", "fibroblast_score"]
    ].mean()
    fallback_score_means = adata.obs.groupby(cluster_key, observed=True)[
        [f"{cell_type}_score2" for cell_type in MANUAL_MARKERS]
    ].mean()
    residual_metrics = []
    for (cluster, cell_type), group in gene_evidence.groupby(["Auto_manual_cluster", "cell_type"], sort=False):
        if cell_type not in {"epithelial", "fibroblast"}:
            continue
        n_prevalent = int((group["pct_cluster"] >= args.residual_min_pct).sum())
        prevalent_markers = group.loc[group["pct_cluster"] >= args.residual_min_pct, "gene"].tolist()
        module_score = float(cluster_score_means.loc[str(cluster), f"{cell_type}_score"])
        residual_metrics.append({
            "Auto_manual_cluster": str(cluster),
            "cell_type": cell_type,
            "residual_module_score": module_score,
            "residual_n_markers_prevalent": n_prevalent,
            "residual_marker_fraction_prevalent": n_prevalent / len(group),
            "residual_prevalent_markers": ";".join(prevalent_markers),
            "residual_required_markers": args.residual_min_markers,
            "passes_residual_evidence": (
                cell_type == "epithelial"
                and
                int(group["cluster_n"].iloc[0]) >= args.marker_min_cluster_cells
                and module_score > args.residual_min_module_score
                and n_prevalent >= args.residual_min_markers
            ),
        })
    residual_metrics = pd.DataFrame(residual_metrics)
    type_evidence = type_evidence.merge(
        residual_metrics, on=["Auto_manual_cluster", "cell_type"], how="left"
    )
    type_evidence["fallback_raw_marker_score"] = [
        float(fallback_score_means.loc[str(cluster), f"{cell_type}_score2"])
        for cluster, cell_type in zip(type_evidence["Auto_manual_cluster"], type_evidence["cell_type"])
    ]
    type_evidence["passes_residual_evidence"] = (
        type_evidence["passes_residual_evidence"].astype("boolean").fillna(False).astype(bool)
    )
    ####################
    return gene_evidence, type_evidence


def assign_clusters_from_marker_evidence(type_evidence, fallback_protected_min_score=0.10):
    assignments = {}
    selected_rows = []
    for cluster, cluster_evidence in type_evidence.groupby("Auto_manual_cluster", sort=False):
        passing = cluster_evidence.loc[cluster_evidence["passes_marker_evidence"]].copy()
        first_pass = passing.loc[~passing["cell_type"].isin(["epithelial", "fibroblast"])]
        if not first_pass.empty:
            candidates = first_pass
            stage = "non_epithelial_first"
        else:
            ####################
            # Fibroblast remains a specificity call because stromal transcripts
            # spill broadly. Only when fibroblast lacks differential ECM-marker
            # support do residual clusters receive an absolute epithelial call.
            ####################
            fibroblast_specific = passing.loc[passing["cell_type"].eq("fibroblast")]
            if not fibroblast_specific.empty:
                candidates = fibroblast_specific
                stage = "fibroblast_specific_residual"
            else:
                candidates = cluster_evidence.loc[
                    cluster_evidence["cell_type"].eq("epithelial")
                    & cluster_evidence["passes_residual_evidence"]
                ]
                stage = "epithelial_absolute_residual"
            ####################
        if candidates.empty:
            ####################
            # The original notebook resolves low-evidence clusters by returning
            # to absolute marker scores. Preserve its protected-lineage-first
            # rule, then force the remaining residual cluster to the larger
            # epithelial/fibroblast raw marker mean so every cluster is labelled.
            ####################
            protected_fallback = cluster_evidence.loc[
                ~cluster_evidence["cell_type"].isin(["epithelial", "fibroblast"])
                & (cluster_evidence["fallback_raw_marker_score"] > fallback_protected_min_score)
            ]
            if not protected_fallback.empty:
                candidates = protected_fallback
                stage = "protected_absolute_fallback"
            else:
                candidates = cluster_evidence.loc[
                    cluster_evidence["cell_type"].isin(["epithelial", "fibroblast"])
                ]
                stage = "epithelial_fibroblast_forced_fallback"
            ####################
        if stage in {"non_epithelial_first", "fibroblast_specific_residual"}:
            selected = candidates.sort_values(
                ["n_markers_supported", "marker_fraction_supported", "marker_median_log2fc", "cell_type"],
                ascending=[False, False, False, True],
            ).iloc[0]
        elif stage == "epithelial_absolute_residual":
            selected = candidates.sort_values(
                ["residual_marker_fraction_prevalent", "residual_module_score", "cell_type"],
                ascending=[False, False, True],
            ).iloc[0]
        else:
            selected = candidates.sort_values(
                ["fallback_raw_marker_score", "cell_type"], ascending=[False, True]
            ).iloc[0]
        assignments[str(cluster)] = selected["cell_type"]
        if stage in {"non_epithelial_first", "fibroblast_specific_residual"}:
            selected_markers = selected["supported_markers"]
            selected_n_markers = int(selected["n_markers_supported"])
            selected_required_markers = int(selected["required_markers"])
            selected_effect = selected["marker_median_log2fc"]
            selected_padjusted = selected["marker_max_padjusted"]
        elif stage == "epithelial_absolute_residual":
            selected_markers = selected["residual_prevalent_markers"]
            selected_n_markers = int(selected["residual_n_markers_prevalent"])
            selected_required_markers = int(selected["residual_required_markers"])
            selected_effect = selected["residual_module_score"]
            selected_padjusted = np.nan
        else:
            selected_markers = ""
            selected_n_markers = 0
            selected_required_markers = 0
            selected_effect = selected["fallback_raw_marker_score"]
            selected_padjusted = np.nan
        selected_rows.append({
            "Auto_manual_cluster": str(cluster),
            "Auto_annotation_celltype": selected["cell_type"],
            "Auto_annotation_evidence_stage": stage,
            "Auto_annotation_supported_markers": selected_markers,
            "Auto_annotation_n_markers_supported": selected_n_markers,
            "Auto_annotation_required_markers": selected_required_markers,
            "Auto_annotation_marker_median_log2fc": selected_effect,
            "Auto_annotation_marker_max_padjusted": selected_padjusted,
        })
    return pd.Series(assignments, dtype=object), pd.DataFrame(selected_rows)
####################


def annotate_segmented(adata, coordinates, args):
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    keep = (adata.obs["total_counts"] >= args.min_counts) & (adata.obs["pct_counts_mt"] <= args.max_mt)
    adata = adata[keep].copy()
    coordinates = coordinates.loc[adata.obs_names].copy()
    if adata.n_obs < 20:
        raise ValueError(f"Only {adata.n_obs} segmented cells remain after QC")
    ####################
    # Raw marker means are retained as descriptive audit fields. Module scores
    # and all annotation evidence are computed after CP10K/log1p normalization,
    # so their interpretation does not depend directly on UMI depth.
    ####################
    score_cols = []
    for cell_type, genes in MANUAL_MARKERS.items():
        valid = [gene for gene in genes if gene in adata.var_names]
        score_col = f"{cell_type}_score"
        mean_col = f"{cell_type}_score2"
        adata.obs[mean_col] = raw_marker_mean(adata, valid)
        score_cols.append(score_col)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for cell_type, genes in MANUAL_MARKERS.items():
        valid = [gene for gene in genes if gene in adata.var_names]
        score_col = f"{cell_type}_score"
        if valid:
            sc.tl.score_genes(adata, gene_list=valid, score_name=score_col, ctrl_size=50, use_raw=False, random_state=0)
        else:
            adata.obs[score_col] = 0.0
    adata.obs["Auto_annotation_score_scale"] = "log1p_cp10k"
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=min(3000, adata.n_vars), span=0.3)
    if adata.var["highly_variable"].sum() < 20:
        raise ValueError("Too few highly variable genes for manual segmented annotation")
    cluster_adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.regress_out(cluster_adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(cluster_adata, max_value=10)
    n_comps = min(50, cluster_adata.n_vars - 1, cluster_adata.n_obs - 1)
    if n_comps < 2:
        raise ValueError("Too few cells/features after QC for manual segmented clustering")
    sc.tl.pca(cluster_adata, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(cluster_adata, n_neighbors=min(15, cluster_adata.n_obs - 1), n_pcs=min(40, n_comps), metric="cosine")
    sc.tl.leiden(cluster_adata, resolution=args.leiden_resolution, key_added="Auto_manual_cluster", flavor="igraph", directed=False)
    adata.obs["Auto_manual_cluster"] = cluster_adata.obs["Auto_manual_cluster"].astype(str).to_numpy()
    gene_evidence, type_evidence = cluster_marker_evidence(adata, "Auto_manual_cluster", args)
    cluster_labels, selected_evidence = assign_clusters_from_marker_evidence(
        type_evidence, args.fallback_protected_min_score
    )
    adata.obs["Auto_annotation_celltype"] = adata.obs["Auto_manual_cluster"].map(cluster_labels)
    selected_evidence = selected_evidence.set_index("Auto_manual_cluster")
    for column in selected_evidence.columns:
        if column == "Auto_annotation_celltype":
            continue
        adata.obs[column] = adata.obs["Auto_manual_cluster"].map(selected_evidence[column])
    adata.obs["Auto_annotation_method"] = "manual_cluster_marker_enrichment"
    adata.obs["Auto_annotation_pass_doublet_filter"] = True
    adata.obs["Auto_annotation_doublet_status"] = "not_applicable_manual_segmented"
    adata.obs["Auto_annotation_keep_epithelial"] = adata.obs["Auto_annotation_celltype"].eq("epithelial")
    annotation = pd.concat([adata.obs.copy(), coordinates], axis=1)
    annotation.insert(0, "barcode", annotation.index.astype(str))
    type_evidence = type_evidence.merge(selected_evidence.reset_index(), on="Auto_manual_cluster", how="left")
    gene_evidence = gene_evidence.merge(
        selected_evidence.reset_index()[["Auto_manual_cluster", "Auto_annotation_celltype", "Auto_annotation_evidence_stage"]],
        on="Auto_manual_cluster",
        how="left",
    )
    return annotation, type_evidence, gene_evidence


def rctd_annotation_table(sample_name, rctd_dir):
    path = Path(rctd_dir) / "tables" / f"Auto_{sample_name}_binned_rctd_annotations.csv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing RCTD annotation table: {path}")
    annotation = pd.read_csv(path, compression="gzip")
    required = {"barcode", "Auto_annotation_celltype", "Auto_annotation_pass_doublet_filter", "Auto_annotation_keep_epithelial"}
    missing = required.difference(annotation.columns)
    if missing:
        raise ValueError(f"RCTD table {path} lacks required columns: {sorted(missing)}")
    annotation["barcode"] = annotation["barcode"].astype(str)
    return annotation


####################
# Custom binned annotation uses RCTD only to identify and annotate non-singlet
# bins. RCTD singlets are passed through the same CP10K marker scoring, Leiden
# clustering, and hierarchical cluster-marker evidence rule as segmentation.
####################
def annotate_custom_binned(adata, coordinates, rctd_annotation, args):
    rctd = rctd_annotation.copy().set_index("barcode", drop=False)
    missing_counts = ~rctd.index.isin(adata.obs_names)
    if missing_counts.any():
        raise ValueError(f"{missing_counts.sum()} RCTD observations are missing from the binned count matrix")
    if "spot_class" not in rctd.columns:
        raise ValueError("RCTD annotation lacks spot_class required for custom singlet/doublet handling")
    rctd_singlet = rctd["spot_class"].astype(str).eq("singlet")
    custom = rctd.copy()
    custom["Auto_rctd_first_type"] = custom.get("first_type", pd.Series(index=custom.index, dtype=object)).astype(str)
    custom["Auto_rctd_second_type"] = custom.get("second_type", pd.Series(index=custom.index, dtype=object)).astype(str)
    custom["Auto_rctd_spot_class"] = custom["spot_class"].astype(str)
    custom["Auto_manual_qc_pass"] = False
    custom["Auto_annotation_celltype"] = "unresolved"
    custom["Auto_annotation_method"] = "manual_cluster_marker_enrichment_rctd_singlet"
    custom["Auto_annotation_pass_doublet_filter"] = rctd_singlet.to_numpy()
    custom["Auto_annotation_keep_epithelial"] = False
    custom["Auto_annotation_doublet_status"] = np.where(rctd_singlet, "RCTD_singlet", "RCTD_" + custom["spot_class"].astype(str))

    non_singlet = ~rctd_singlet
    doublet_label = custom.loc[non_singlet, "Auto_rctd_first_type"].copy()
    second_type = custom.loc[non_singlet, "Auto_rctd_second_type"]
    has_second = second_type.notna() & second_type.ne("") & second_type.ne("nan")
    doublet_label.loc[has_second] = doublet_label.loc[has_second] + "|" + second_type.loc[has_second]
    custom.loc[non_singlet, "Auto_annotation_celltype"] = doublet_label
    custom.loc[non_singlet, "Auto_annotation_method"] = "RCTD_doublet"

    singlet_barcodes = rctd.index[rctd_singlet]
    singlet_adata = adata[singlet_barcodes].copy()
    singlet_coordinates = coordinates.loc[singlet_barcodes].copy()
    ####################
    # The same marker-enrichment decision rule is used for custom singlet bins
    # and segmented cells. Only the Leiden resolution differs, because the bin
    # matrix contains many more observations and rare lineages require finer
    # clusters to remain detectable.
    ####################
    manual, type_evidence, gene_evidence = annotate_segmented(
        singlet_adata,
        singlet_coordinates,
        args,
    )
    manual = manual.set_index("barcode", drop=False)
    manual_barcodes = manual.index.intersection(custom.index)
    custom.loc[manual_barcodes, "Auto_manual_qc_pass"] = True
    custom.loc[manual_barcodes, "Auto_annotation_celltype"] = manual.loc[manual_barcodes, "Auto_annotation_celltype"]
    custom.loc[manual_barcodes, "Auto_annotation_method"] = "manual_cluster_marker_enrichment_rctd_singlet"
    custom.loc[manual_barcodes, "Auto_annotation_keep_epithelial"] = manual.loc[manual_barcodes, "Auto_annotation_keep_epithelial"].astype(bool)
    for column in [
        "Auto_manual_cluster", "Auto_annotation_score_scale", "Auto_annotation_evidence_stage",
        "Auto_annotation_supported_markers", "Auto_annotation_n_markers_supported",
        "Auto_annotation_required_markers", "Auto_annotation_marker_median_log2fc",
        "Auto_annotation_marker_max_padjusted",
    ]:
        if column in manual.columns:
            custom.loc[manual_barcodes, column] = manual.loc[manual_barcodes, column]
    for column in [col for col in manual.columns if col.endswith("_score") or col.endswith("_score2")]:
        custom.loc[manual_barcodes, column] = manual.loc[manual_barcodes, column]
    custom["sample"] = custom.get("sample", pd.Series(index=custom.index, dtype=object)).fillna("")
    return custom.reset_index(drop=True), type_evidence, gene_evidence
####################


def annotation_output_path(output_dir, sample_name, mode):
    return Path(output_dir) / "tables" / f"Auto_{sample_name}_{mode}_cell_annotations.csv.gz"


def run_annotation_stage(args):
    out_dir = Path(args.output_dir)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for path_str, sample_name in zip(args.inputs, args.sample_names):
        path = Path(path_str)
        type_evidence = None
        gene_evidence = None
        if args.mode == "binned":
            if not args.rctd_dir:
                raise ValueError("--rctd-dir is required for binned annotation")
            annotation = rctd_annotation_table(sample_name, args.rctd_dir)
        elif args.mode == "custom":
            if not args.rctd_dir:
                raise ValueError("--rctd-dir is required for custom annotation")
            adata = read_10x_counts(path)
            coordinates = read_spatial_coordinates(path, args.mode, adata.obs_names)
            annotation, type_evidence, gene_evidence = annotate_custom_binned(
                adata, coordinates, rctd_annotation_table(sample_name, args.rctd_dir), args
            )
            annotation["sample"] = sample_name
        else:
            adata = read_10x_counts(path)
            coordinates = read_spatial_coordinates(path, args.mode, adata.obs_names)
            annotation, type_evidence, gene_evidence = annotate_segmented(adata, coordinates, args)
            annotation["sample"] = sample_name
        annotation.to_csv(annotation_output_path(out_dir, sample_name, args.mode), index=False, compression="gzip")
        ####################
        # Persist cluster/type and cluster/gene evidence so threshold choices
        # can be audited or re-evaluated without rerunning count processing.
        ####################
        if type_evidence is not None:
            type_evidence.insert(0, "sample", sample_name)
            type_evidence.insert(1, "mode", args.mode)
            type_evidence.to_csv(
                out_dir / "tables" / f"Auto_{sample_name}_{args.mode}_cluster_marker_evidence.csv", index=False
            )
        if gene_evidence is not None:
            gene_evidence.insert(0, "sample", sample_name)
            gene_evidence.insert(1, "mode", args.mode)
            gene_evidence.to_csv(
                out_dir / "tables" / f"Auto_{sample_name}_{args.mode}_cluster_gene_evidence.csv", index=False
            )
        ####################
        labels = annotation["Auto_annotation_celltype"].fillna("NA").value_counts()
        summary_rows.append({
            "sample": sample_name,
            "mode": args.mode,
            "n_annotated": len(annotation),
            "n_pass_doublet_filter": int(annotation["Auto_annotation_pass_doublet_filter"].astype(bool).sum()),
            "n_epithelial": int(annotation["Auto_annotation_keep_epithelial"].astype(bool).sum()),
            "n_unresolved": int(annotation["Auto_annotation_celltype"].astype(str).eq("unresolved").sum()),
            "celltype_counts": ";".join(f"{label}={count}" for label, count in labels.items()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "tables" / f"Auto_visiumhd_{args.mode}_annotation_summary.csv", index=False)
    summary_dir = Path("updates/new_updates/summaries")
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_dir / f"visiumhd_{args.mode}_annotation_summary.csv", index=False)
    ####################
    threshold_summary = pd.DataFrame([{
        "mode": args.mode,
        "marker_min_log2fc": args.marker_min_log2fc,
        "marker_min_pct": args.marker_min_pct,
        "marker_min_pct_delta": args.marker_min_pct_delta,
        "marker_fdr": args.marker_fdr,
        "marker_min_cluster_cells": args.marker_min_cluster_cells,
        "residual_min_module_score": args.residual_min_module_score,
        "residual_min_pct": args.residual_min_pct,
        "residual_min_markers": args.residual_min_markers,
        "fallback_protected_min_score": args.fallback_protected_min_score,
        "minimum_supported_markers": "2_for_panels_ge3_genes_else_1",
        "decision_order": "non_epi_specific_then_fib_specific_then_epi_absolute_then_notebook_forced_fallback",
    }])
    threshold_summary.to_csv(out_dir / "tables" / f"Auto_visiumhd_{args.mode}_marker_thresholds.csv", index=False)
    ####################


def read_malignancy_table(malignancy_dir, sample_name, mode):
    path = Path(malignancy_dir) / "tables" / f"Auto_{sample_name}_{mode}_infercna_cells.csv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing InferCNA malignancy table: {path}")
    table = pd.read_csv(path, compression="gzip", low_memory=False)
    required = {"barcode", "Auto_malignant"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"InferCNA table {path} lacks required columns: {sorted(missing)}")
    table["barcode"] = table["barcode"].astype(str)
    table["Auto_malignant"] = table["Auto_malignant"].astype(str).str.lower().eq("true")
    return table.set_index("barcode")


def load_epithelial_data_for_mapping(args):
    loaded = []
    if not args.malignancy_dir:
        raise ValueError("--malignancy-dir is required for state mapping")
    for path_str, sample_name in zip(args.inputs, args.sample_names):
        path = Path(path_str)
        adata = read_10x_counts(path)
        coordinates = read_spatial_coordinates(path, args.mode, adata.obs_names)
        annotation = pd.read_csv(annotation_output_path(args.output_dir, sample_name, args.mode), compression="gzip")
        annotation["barcode"] = annotation["barcode"].astype(str)
        annotation = annotation.set_index("barcode")
        epithelial = annotation.index[annotation["Auto_annotation_keep_epithelial"].astype(bool)]
        epithelial = adata.obs_names.intersection(epithelial)
        if len(epithelial) == 0:
            raise ValueError(f"No annotated epithelial observations for {sample_name} ({args.mode})")
        malignancy = read_malignancy_table(args.malignancy_dir, sample_name, args.mode)
        missing_malignancy = epithelial.difference(malignancy.index)
        if len(missing_malignancy) > 0:
            raise ValueError(f"{sample_name} has {len(missing_malignancy)} epithelial observations without InferCNA calls")
        adata = adata[epithelial].copy()
        obs = annotation.loc[epithelial].copy()
        obs = obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], errors="ignore")
        obs["sample"] = sample_name
        obs["Auto_malignant"] = malignancy.loc[epithelial, "Auto_malignant"].to_numpy()
        for column in [
            "Auto_malignancy",
            "Auto_malignancy_evidence",
            "Auto_malignancy_before_keratinocyte_exclusion",
            "Auto_cna_class",
            "Auto_cancer_signature_score",
            "Auto_cancer_signature_status",
            "Auto_segmented_keratinocyte_n_cells",
            "Auto_segmented_keratinocyte_fraction",
            "Auto_segmented_keratinocyte_nearest_distance",
            ####################
            # Preserve full-CNA profile provenance in the state-mapping output.
            ####################
            "Auto_keratinocyte_cna_profile_class",
            "Auto_keratinocyte_tumour_profile_correlation",
            "Auto_keratinocyte_normal_profile_correlation",
            "Auto_keratinocyte_tumour_normal_margin",
            "Auto_keratinocyte_tumour_profile_threshold",
            "Auto_keratinocyte_margin_threshold",
            "Auto_malignancy_before_keratinocyte_profile",
            ####################
            "Auto_binned_malignant",
            "Auto_nearest_binned_barcode",
            "Auto_nearest_binned_distance",
        ]:
            if column in malignancy.columns:
                obs[column] = malignancy.loc[epithelial, column].to_numpy()
        loaded.append((adata, obs, coordinates.loc[epithelial].copy()))
    return loaded


def run_mapping_stage(args):
    out_dir = Path(args.output_dir)
    table_dir = out_dir / "tables"
    intermediate_dir = out_dir / "intermediate"
    figure_dir = out_dir / "figures"
    for directory in [table_dir, intermediate_dir, figure_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    ranked, mp_order, mp_desc_map = load_signatures(args.signature_dir, args.top_n)
    loaded = load_epithelial_data_for_mapping(args)
    common_genes = loaded[0][0].var_names
    for adata, _, _ in loaded[1:]:
        common_genes = common_genes.intersection(adata.var_names, sort=False)
    if len(common_genes) < 1000:
        raise ValueError(f"Only {len(common_genes)} common genes remain across epithelial inputs")
    matrices, obs_frames, spatial_frames = [], [], []
    for adata, obs, spatial in loaded:
        matrices.append(adata[:, common_genes].X)
        obs_frames.append(obs)
        spatial_frames.append(spatial)
    obs = pd.concat(obs_frames, axis=0)
    spatial = pd.concat(spatial_frames, axis=0)
    X = sp.vstack(matrices).tocsr()
    obs.index = obs["sample"].astype(str) + "_" + obs.index.astype(str)
    spatial.index = obs.index
    norm_matrix = normalise_log1p(X)
    gene_index = pd.Index(common_genes)
    signature_map = {}
    for mp_name, mp_df in ranked.groupby("mp", sort=False):
        indexer = gene_index.get_indexer(mp_df["gene"])
        signature_map[mp_name] = indexer[indexer >= 0].tolist()
    mp_raw, mp_adj = zscore_by_sample(norm_matrix, obs, signature_map)
    mp_adj_noncc = mp_adj[[mp for mp in mp_adj.columns if mp not in CC_MPS]].copy()
    group_scores, state_assignment, state_gap = assign_states(mp_adj_noncc, args.threshold, args.hybrid_gap)
    top_mp = top_mp_labels(mp_adj_noncc, args.threshold, mp_order)
    epithelial_scores = pd.concat([
        obs, spatial[["pxl_row_in_fullres", "pxl_col_in_fullres"]],
        mp_raw.add_prefix("Auto_raw_"), mp_adj.add_prefix("Auto_adj_"),
        group_scores, state_assignment, state_gap, top_mp,
    ], axis=1)
    epithelial_scores["Auto_top_mp_label"] = epithelial_scores["Auto_top_mp"].apply(lambda value: label_mp(value, mp_desc_map))
    epithelial_scores.to_csv(intermediate_dir / f"Auto_visiumhd_{args.mode}_epithelial_prestate_scores.csv.gz", compression="gzip")
    results = epithelial_scores.loc[epithelial_scores["Auto_malignant"].astype(bool)].copy()
    if results.empty:
        raise ValueError(f"No malignant epithelial observations passed InferCNA for {args.mode}")
    results.to_csv(table_dir / f"Auto_visiumhd_{args.mode}_malignant_epithelial_state_annotations.csv.gz", compression="gzip")
    summary = results.groupby("sample", observed=True).agg(
        n_malignant_epithelial=("Auto_malignant", "size"),
        n_states=("Auto_state_B", "nunique"),
    ).reset_index()
    summary.to_csv(table_dir / f"Auto_visiumhd_{args.mode}_malignant_epithelial_state_summary.csv", index=False)
    summary_dir = Path("updates/new_updates/summaries")
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_dir / f"visiumhd_{args.mode}_malignant_epithelial_state_summary.csv", index=False)
    for sample_name in args.sample_names:
        sub = results.loc[results["sample"] == sample_name].copy()
        if sub.empty:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(15, 7))
        plot_categorical(axes[0], sub, "Auto_state_B", STATE_COLORS, f"{sample_name} ({args.mode}): malignant epithelial scATLAS state")
        plot_categorical(
            axes[1], sub, "Auto_top_mp_label",
            {label_mp(mp, mp_desc_map): MP_COLORS.get(mp, "grey") for mp in mp_order},
            f"{sample_name} ({args.mode}): malignant epithelial top non-CC MP",
        )
        fig.tight_layout()
        fig.savefig(figure_dir / f"Auto_{sample_name}_{args.mode}_malignant_epithelial_state_map.png", dpi=400, bbox_inches="tight")
        with PdfPages(figure_dir / f"Auto_{sample_name}_{args.mode}_malignant_epithelial_state_map.pdf") as pdf:
            pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def gated_main():
    args = gated_parse_args()
    if args.stage == "annotate":
        run_annotation_stage(args)
    else:
        run_mapping_stage(args)
    log_dir = Path(args.output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / f"Auto_visiumhd_{args.mode}_{args.stage}_run_summary.txt").write_text(
        "\n".join([
            f"stage={args.stage}",
            f"mode={args.mode}",
            f"samples={','.join(args.sample_names)}",
            f"inputs={','.join(args.inputs)}",
            f"output_dir={args.output_dir}",
            f"rctd_dir={args.rctd_dir or ''}",
            f"malignancy_dir={args.malignancy_dir or ''}",
            f"threshold={args.threshold}",
            f"hybrid_gap={args.hybrid_gap}",
            ####################
            f"marker_min_log2fc={args.marker_min_log2fc}",
            f"marker_min_pct={args.marker_min_pct}",
            f"marker_min_pct_delta={args.marker_min_pct_delta}",
            f"marker_fdr={args.marker_fdr}",
            f"marker_min_cluster_cells={args.marker_min_cluster_cells}",
            f"residual_min_module_score={args.residual_min_module_score}",
            f"residual_min_pct={args.residual_min_pct}",
            f"residual_min_markers={args.residual_min_markers}",
            ####################
        ]) + "\n"
    )


if __name__ == "__main__":
    gated_main()
####################
