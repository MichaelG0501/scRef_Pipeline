#!/usr/bin/env python
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_celltype_annotation.py
#   Description: Apply one shared high-resolution marker-score annotation to all
#     >=100-UMI 16 um bins and Space Ranger segmented cells, with targeted
#     one-pass expression-graph refinement of raw-rejected ranked candidates.
#   Methodology: analysis/methodology/spatial/visium_hd_final_annotation_methodology.md
#   Inputs: analysis/spatial/visium_hd_samples.tsv;
#     analysis/shared/visium_hd_celltype_colours.tsv; Space Ranger count and
#     coordinate outputs; binned RCTD tables.
#   Outputs:
#     intermediate/: expression UMAP coordinates.
#     tables/: per-observation annotations, cluster/type evidence, marker-level
#       support, cluster score calls, and compact summaries.
#     figures/: matched annotation, cluster, and marker-score diagnostics.
#     logs/: lightweight run summary.
#   Cache/replot: annotation tables are reusable downstream; heavy annotation
#     is skipped only when --reuse-complete is supplied.
#   Run: python analysis/spatial/visium_hd_celltype_annotation.py
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

####################
import argparse
from pathlib import Path
import warnings

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
import scanpy as sc
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=PerformanceWarning)

WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")

MARKERS = {
    "erythrocyte": ["HBA1", "HBA2", "HBB"],
    "keratinocyte": ["FLG", "IVL"],
    "lymph": ["CCL21"],
    "neutrophil": ["CTSG", "ELANE", "MPO", "AZU1"],
    "endothelial": ["ENG", "CLEC14A", "CLDN5", "VWF", "CDH5"],
    "epithelial": ["KRT7", "MUC1", "KRT19", "EPCAM"],
    "fibroblast": ["COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"],
    "b.cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "plasma": ["MZB1", "JCHAIN", "DERL3"],
    "dendritic": ["CLEC10A", "CCR7", "CD86"],
    "macrophage": ["CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"],
    "mast": ["MS4A2", "CPA3", "TPSB2", "TPSAB1"],
    "nk.cell": ["GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"],
    "t.cell": ["CD3E", "CD3D", "CD2", "CD3G"],
}

_celltype_colour_table = pd.DataFrame({
    'celltype': ['epithelial', 'fibroblast', 'endothelial', 'macrophage', 'mast', 't.cell', 'b.cell', 'nk.cell', 'plasma', 'dendritic', 'lymph', 'erythrocyte', 'keratinocyte', 'neutrophil', 'unresolved', 'combined'],
    'colour': ['#D73027', '#8C564B', '#1F78B4', '#FF7F00', '#A65628', '#33A02C', '#377EB8', '#984EA3', '#E377C2', '#17BECF', '#6BAED6', '#7F7F7F', '#E6AB02', '#1B9E77', '#BDBDBD', '#555555']
})
CELLTYPE_COLOURS = dict(
    zip(_celltype_colour_table["celltype"], _celltype_colour_table["colour"])
)
COMBINED_CELLTYPE_COLOUR = CELLTYPE_COLOURS.pop("combined")


def celltype_colour(label):
    label = str(label)
    if "|" in label:
        return COMBINED_CELLTYPE_COLOUR
    return CELLTYPE_COLOURS.get(label, COMBINED_CELLTYPE_COLOUR)


# Replicated detection avoids assigning a high-resolution cluster from one
# isolated transcript without imposing a representation-specific percentage.
MARKER_MIN_CLUSTER_CELLS = 5
MARKER_MIN_POSITIVE_CELLS = 2
RAW_NORMALISED_SCORE_THRESHOLD = 0.25
STRUCTURAL_DOUBLET_STANDARDISED_GAP = 0.25
STRUCTURAL_TYPES = {"epithelial", "fibroblast"}
CLUSTERS_PER_HEATMAP_PAGE = 35
# Deliberate annotation overclustering limits rare-lineage dilution. Very
# small datasets retain the lower resolution to avoid median cluster sizes
# below five observations.
LEIDEN_RESOLUTION_LARGE = 10.0
LEIDEN_RESOLUTION_SMALL = 6.0
LEIDEN_LARGE_MIN_OBSERVATIONS = 1000
MIN_COUNTS = 100
MAX_MT_PERCENT = 15.0
GRAPH_MIN_GENE_CELLS = 10
GRAPH_SCALE_MAX = 10.0
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0

####################
# Targeted refinement records only the ranked cell types rejected before the
# first raw-score-passing assignment. A parent is reclustered once when any
# such candidate has enough individually high-scoring observations.
TARGETED_REFINEMENT_CELL_SCORE = 1.0
TARGETED_REFINEMENT_MIN_HIGH_CELLS = 20
TARGETED_REFINEMENT_LOCAL_RESOLUTION = 1.0
####################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Final two-method Visium HD cell-type annotation workflow."
    )
    parser.add_argument(
        "--manifest",
        default="",
    )
    parser.add_argument(
        "--output-dir",
        default=str(WD / "ref_outs" / "visium_hd_outs"),
    )
    parser.add_argument(
        "--rctd-dir",
        default=str(WD / "ref_outs" / "visium_hd_outs" / "rctd"),
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=["binned", "segmented"],
        default=["binned", "segmented"],
    )
    parser.add_argument("--reuse-complete", action="store_true")
    parser.add_argument("--force-samples", nargs="*", default=[])
    parser.add_argument("--samples", nargs="*", default=[])
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=None,
        help="Override the size-aware default resolution for audit runs.",
    )
    parser.add_argument(
        "--raw-score-threshold",
        type=float,
        default=RAW_NORMALISED_SCORE_THRESHOLD,
    )
    parser.add_argument(
        "--structural-doublet-gap",
        type=float,
        default=STRUCTURAL_DOUBLET_STANDARDISED_GAP,
    )
    return parser.parse_args()


def read_counts(input_dir):
    candidates = [
        input_dir / "filtered_feature_bc_matrix.h5",
        input_dir / "filtered_feature_cell_matrix.h5",
    ]
    for path in candidates:
        if path.exists():
            adata = sc.read_10x_h5(str(path), gex_only=True)
            adata.var_names_make_unique()
            return adata
    for name in ["filtered_feature_bc_matrix", "filtered_feature_cell_matrix"]:
        path = input_dir / name
        if path.exists():
            adata = sc.read_10x_mtx(str(path), gex_only=True)
            adata.var_names_make_unique()
            return adata
    raise FileNotFoundError(f"No filtered expression matrix in {input_dir}")


def read_binned_coordinates(input_dir, barcodes):
    parquet_path = input_dir / "spatial" / "tissue_positions.parquet"
    csv_path = input_dir / "spatial" / "tissue_positions.csv"
    if parquet_path.exists():
        positions = pd.read_parquet(parquet_path)
    elif csv_path.exists():
        positions = pd.read_csv(csv_path)
    else:
        raise FileNotFoundError(f"No tissue positions in {input_dir / 'spatial'}")
    barcode_col = next(
        (column for column in ["barcode", "Barcode"] if column in positions.columns),
        None,
    )
    if barcode_col is None:
        raise ValueError(f"No barcode column in tissue positions for {input_dir}")
    x_col = next(
        (
            column
            for column in ["pxl_col_in_fullres", "array_col"]
            if column in positions.columns
        ),
        None,
    )
    y_col = next(
        (
            column
            for column in ["pxl_row_in_fullres", "array_row"]
            if column in positions.columns
        ),
        None,
    )
    if x_col is None or y_col is None:
        raise ValueError(f"No supported coordinates for {input_dir}")
    positions[barcode_col] = positions[barcode_col].astype(str)
    positions = positions.set_index(barcode_col).reindex(barcodes)
    return pd.DataFrame(
        {
            "pxl_col_in_fullres": positions[x_col].to_numpy(dtype=float),
            "pxl_row_in_fullres": positions[y_col].to_numpy(dtype=float),
        },
        index=barcodes,
    )


def read_segmented_coordinates(input_dir, barcodes):
    geojson_path = input_dir / "cell_segmentations.geojson"
    if not geojson_path.exists():
        raise FileNotFoundError(f"Missing {geojson_path}")
    segmentations = gpd.read_file(geojson_path)
    segmentations["cell_id"] = segmentations["cell_id"].astype(str)
    segmentations = segmentations.set_index("cell_id", drop=False)
    cell_ids = (
        pd.Index(barcodes)
        .to_series(index=barcodes)
        .str.extract(r"(\d+)", expand=False)
    )
    if cell_ids.isna().any():
        raise ValueError("Could not parse one or more segmented cell IDs")
    cell_ids = cell_ids.astype(int).astype(str)
    missing = ~cell_ids.isin(segmentations.index)
    if missing.any():
        raise ValueError(
            f"{int(missing.sum())} expression cells are absent from segmentation geometry"
        )
    centroids = [
        geometry.centroid
        for geometry in segmentations.loc[cell_ids.to_numpy()].geometry
    ]
    return pd.DataFrame(
        {
            "pxl_col_in_fullres": [centroid.x for centroid in centroids],
            "pxl_row_in_fullres": [centroid.y for centroid in centroids],
        },
        index=barcodes,
    )


def marker_required_count(n_available):
    return 1 if n_available <= 2 else 2


def normalised_marker_scores(adata):
    scores = {}
    for cell_type, genes in MARKERS.items():
        valid = [gene for gene in genes if gene in adata.var_names]
        if not valid:
            scores[cell_type] = np.zeros(adata.n_obs, dtype=np.float32)
            continue
        values = adata[:, valid].X
        scores[cell_type] = np.asarray(values.mean(axis=1)).ravel().astype(
            np.float32
        )
    return pd.DataFrame(scores, index=adata.obs_names)


def cluster_marker_support(
    adata,
    cluster_key,
    score_frame,
    standardised_score_frame,
):
    marker_genes = sorted(
        {
            gene
            for genes in MARKERS.values()
            for gene in genes
            if gene in adata.var_names
        }
    )
    marker_indices = adata.var_names.get_indexer(marker_genes)
    counts = adata.layers["counts"][:, marker_indices]
    if sp.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float64)
    expression = adata[:, marker_genes].X
    if sp.issparse(expression):
        expression = expression.toarray()
    expression = np.asarray(expression, dtype=np.float64)
    detected = counts > 0
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    gene_index = {gene: index for index, gene in enumerate(marker_genes)}
    gene_rows = []

    for cluster in sorted(pd.unique(clusters)):
        selected = clusters == cluster
        cluster_n = int(selected.sum())
        cluster_detected = detected[selected].sum(axis=0)
        mean_expression = expression[selected].mean(axis=0)
        for cell_type, genes in MARKERS.items():
            for gene in genes:
                if gene not in gene_index:
                    continue
                index = gene_index[gene]
                supported = (
                    cluster_n >= MARKER_MIN_CLUSTER_CELLS
                    and cluster_detected[index] >= MARKER_MIN_POSITIVE_CELLS
                )
                gene_rows.append(
                    {
                        "Auto_manual_cluster": cluster,
                        "cell_type": cell_type,
                        "gene": gene,
                        "cluster_n": cluster_n,
                        "positive_cells": int(cluster_detected[index]),
                        "detection_fraction": (
                            float(cluster_detected[index]) / cluster_n
                        ),
                        "mean_log1p_cp10k": float(mean_expression[index]),
                        "supported": bool(supported),
                    }
                )

    gene_evidence = pd.DataFrame(gene_rows)
    cluster_scores = score_frame.assign(
        Auto_manual_cluster=adata.obs[cluster_key].astype(str).to_numpy()
    ).groupby("Auto_manual_cluster", observed=True).mean()
    cluster_standardised_scores = standardised_score_frame.assign(
        Auto_manual_cluster=adata.obs[cluster_key].astype(str).to_numpy()
    ).groupby("Auto_manual_cluster", observed=True).mean()
    type_rows = []
    for (cluster, cell_type), group in gene_evidence.groupby(
        ["Auto_manual_cluster", "cell_type"], sort=False
    ):
        supported = group.loc[group["supported"]]
        required = marker_required_count(len(group))
        type_rows.append(
            {
                "Auto_manual_cluster": str(cluster),
                "cell_type": cell_type,
                "cluster_n": int(group["cluster_n"].iloc[0]),
                "n_markers_available": len(group),
                "n_markers_supported": len(supported),
                "required_markers": required,
                "marker_fraction_supported": len(supported) / len(group),
                "marker_minimum_positive_cells": (
                    int(supported["positive_cells"].min())
                    if len(supported)
                    else 0
                ),
                "marker_mean_detection_fraction": (
                    float(supported["detection_fraction"].mean())
                    if len(supported)
                    else 0.0
                ),
                "supported_markers": ";".join(supported["gene"]),
                "passes_marker_evidence": len(supported) >= required,
                "cluster_normalised_marker_score": float(
                    cluster_scores.loc[str(cluster), cell_type]
                ),
                "cluster_standardised_marker_score": float(
                    cluster_standardised_scores.loc[str(cluster), cell_type]
                ),
            }
        )
    return gene_evidence, pd.DataFrame(type_rows)


def choose_cluster_labels(type_evidence):
    assignments = {}
    selected_rows = []
    for cluster, evidence in type_evidence.groupby(
        "Auto_manual_cluster", sort=False
    ):
        ranked = evidence.sort_values(
            [
                "cluster_standardised_marker_score",
                "cluster_normalised_marker_score",
                "cell_type",
            ],
            ascending=[False, False, True],
        ).reset_index(drop=True)
        ranked["standardised_rank"] = np.arange(1, len(ranked) + 1)
        ranked["passes_raw_threshold"] = (
            ranked["cluster_normalised_marker_score"]
            > RAW_NORMALISED_SCORE_THRESHOLD
        )
        passing = ranked.loc[ranked["passes_raw_threshold"]].copy()

        initial = ranked.iloc[0]
        if passing.empty:
            final_label = "unresolved"
            selected = initial.copy()
            selected["cell_type"] = "unresolved"
            selected_rank = np.nan
            stage = "all_raw_scores_below_threshold"
            second = None
            score_gap = np.nan
            doublet_partner = ""
            doublet_gap = np.nan
        else:
            selected = passing.iloc[0]
            selected_rank = int(selected["standardised_rank"])
            final_label = str(selected["cell_type"])
            stage = (
                "top_standardised_score_passed_raw_threshold"
                if selected_rank == 1
                else f"rank_{selected_rank}_passed_raw_threshold"
            )
            second = passing.iloc[1] if len(passing) > 1 else None
            score_gap = (
                float(
                    selected["cluster_standardised_marker_score"]
                    - second["cluster_standardised_marker_score"]
                )
                if second is not None
                else np.nan
            )
            doublet_partner = ""
            doublet_gap = np.nan

            if final_label in STRUCTURAL_TYPES:
                non_structural = passing.loc[
                    ~passing["cell_type"].isin(STRUCTURAL_TYPES)
                ]
                if not non_structural.empty:
                    partner = non_structural.iloc[0]
                    doublet_gap = float(
                        selected["cluster_standardised_marker_score"]
                        - partner["cluster_standardised_marker_score"]
                    )
                    if (
                        doublet_gap
                        <= STRUCTURAL_DOUBLET_STANDARDISED_GAP
                    ):
                        doublet_partner = str(partner["cell_type"])
                        final_label = (
                            f"{doublet_partner}|{selected['cell_type']}"
                        )
                        stage = "structural_nonstructural_score_doublet"

        excluded = ranked.loc[
            ~ranked["passes_raw_threshold"],
            "cell_type",
        ].tolist()
        rejected_before_selection = ranked.iloc[
            : (len(ranked) if passing.empty else int(selected_rank) - 1)
        ].copy()
        rejected_candidates = rejected_before_selection["cell_type"].tolist()
        rejected_candidate_details = ";".join(
            (
                f"{row.cell_type}:z="
                f"{row.cluster_standardised_marker_score:.4f},raw="
                f"{row.cluster_normalised_marker_score:.4f}"
            )
            for row in rejected_before_selection.itertuples()
        )
        ranking = ";".join(
            (
                f"{row.cell_type}:z="
                f"{row.cluster_standardised_marker_score:.4f},raw="
                f"{row.cluster_normalised_marker_score:.4f},"
                f"pass={bool(row.passes_raw_threshold)}"
            )
            for row in ranked.itertuples()
        )
        assignments[str(cluster)] = final_label
        selected_rows.append(
            {
                "Auto_manual_cluster": str(cluster),
                "Auto_annotation_cluster_n": int(initial["cluster_n"]),
                "Auto_annotation_celltype": final_label,
                "Auto_annotation_primary_celltype": (
                    str(selected["cell_type"])
                    if final_label != "unresolved"
                    else "unresolved"
                ),
                "Auto_annotation_evidence_stage": stage,
                "Auto_annotation_initial_celltype": str(initial["cell_type"]),
                "Auto_annotation_initial_standardised_score": float(
                    initial["cluster_standardised_marker_score"]
                ),
                "Auto_annotation_initial_raw_score": float(
                    initial["cluster_normalised_marker_score"]
                ),
                "Auto_annotation_selected_rank": selected_rank,
                "Auto_annotation_cluster_normalised_score": selected[
                    "cluster_normalised_marker_score"
                ],
                "Auto_annotation_cluster_standardised_score": selected[
                    "cluster_standardised_marker_score"
                ],
                "Auto_annotation_second_celltype": (
                    second["cell_type"] if second is not None else ""
                ),
                "Auto_annotation_score_gap": score_gap,
                "Auto_annotation_doublet_partner": doublet_partner,
                "Auto_annotation_doublet_gap": doublet_gap,
                "Auto_annotation_excluded_below_raw_threshold": ";".join(
                    excluded
                ),
                "Auto_annotation_rejected_candidate_celltypes": ";".join(
                    rejected_candidates
                ),
                "Auto_annotation_rejected_candidate_details": (
                    rejected_candidate_details
                ),
                "Auto_annotation_score_ranking": ranking,
                "Auto_annotation_raw_threshold": (
                    RAW_NORMALISED_SCORE_THRESHOLD
                ),
            }
        )
    return pd.Series(assignments, dtype=object), pd.DataFrame(selected_rows)


####################
def targeted_refine_weak_top_clusters(
    adata,
    marker_scores,
    initial_calls,
):
    parent_labels = adata.obs["Auto_manual_cluster"].astype(str).to_numpy()
    refined_labels = parent_labels.copy().astype(object)
    refinement_rows = []
    for call in initial_calls.itertuples(index=False):
        parent = str(call.Auto_manual_cluster)
        parent_selected = parent_labels == parent
        parent_n = int(parent_selected.sum())
        parent_marker_scores = marker_scores.loc[parent_selected]
        rejected_candidates = [
            value
            for value in str(
                call.Auto_annotation_rejected_candidate_celltypes
            ).split(";")
            if value
        ]
        high_counts = {
            cell_type: int(
                (
                    parent_marker_scores[cell_type]
                    > TARGETED_REFINEMENT_CELL_SCORE
                ).sum()
            )
            for cell_type in rejected_candidates
        }
        eligible_candidates = [
            cell_type
            for cell_type in rejected_candidates
            if high_counts[cell_type] >= TARGETED_REFINEMENT_MIN_HIGH_CELLS
        ]
        eligible = bool(eligible_candidates)
        base_row = {
            "parent_cluster": parent,
            "parent_n": parent_n,
            "rejected_candidate_celltypes": ";".join(rejected_candidates),
            "rejected_candidate_details": str(
                call.Auto_annotation_rejected_candidate_details
            ),
            "candidate_high_score_counts": ";".join(
                f"{cell_type}={high_counts[cell_type]}"
                for cell_type in rejected_candidates
            ),
            "eligible_candidate_celltypes": ";".join(eligible_candidates),
            "minimum_high_score_observations": (
                TARGETED_REFINEMENT_MIN_HIGH_CELLS
            ),
            "eligible": eligible,
            "selected_local_resolution": (
                TARGETED_REFINEMENT_LOCAL_RESOLUTION if eligible else np.nan
            ),
            "n_child_clusters": 1,
            "child_clusters": parent,
            "status": "not_eligible_no_rejected_candidate_with_20_high_scores",
        }
        if not eligible:
            refinement_rows.append(base_row)
            continue

        parent_indices = np.flatnonzero(parent_selected)
        ####################
        # Recompute an expression graph within the eligible parent. Reusing
        # the induced global graph failed to expose the C120 T-cell axis
        # because global HVG/PCA selection is dominated by between-lineage
        # variation. This local graph still contains no marker-score features.
        parent_graph = adata[parent_indices].copy()
        sc.pp.filter_genes(
            parent_graph,
            min_cells=min(GRAPH_MIN_GENE_CELLS, parent_graph.n_obs),
        )
        sc.pp.highly_variable_genes(
            parent_graph,
            flavor="seurat",
            n_top_genes=min(3000, parent_graph.n_vars),
            span=0.3,
        )
        if int(parent_graph.var["highly_variable"].sum()) < 20:
            base_row["status"] = "eligible_too_few_local_hvgs"
            refinement_rows.append(base_row)
            continue
        parent_graph = parent_graph[
            :, parent_graph.var["highly_variable"]
        ].copy()
        sc.pp.regress_out(
            parent_graph,
            ["total_counts", "pct_counts_mt"],
        )
        sc.pp.scale(parent_graph, max_value=GRAPH_SCALE_MAX)
        local_n_comps = min(
            50,
            parent_graph.n_vars - 1,
            parent_graph.n_obs - 1,
        )
        sc.tl.pca(
            parent_graph,
            n_comps=local_n_comps,
            svd_solver="arpack",
        )
        sc.pp.neighbors(
            parent_graph,
            n_neighbors=min(15, parent_graph.n_obs - 1),
            n_pcs=min(40, local_n_comps),
            metric="cosine",
        )
        ####################
        key = "Auto_local_cluster"
        sc.tl.leiden(
            parent_graph,
            resolution=TARGETED_REFINEMENT_LOCAL_RESOLUTION,
            key_added=key,
            flavor="igraph",
            directed=False,
            random_state=0,
        )
        local_labels = parent_graph.obs[key].astype(str).to_numpy()
        local_communities = sorted(pd.unique(local_labels))
        if len(local_communities) < 2:
            base_row["status"] = "eligible_local_resolution_returned_one_cluster"
            refinement_rows.append(base_row)
            continue

        ####################
        # Retain every resolution-1 community. The usual standardized ranking
        # and raw-score validation are rerun once after all eligible parents
        # have been split; no child can trigger another refinement pass.
        suffixes = {
            community: (
                chr(ord("a") + index)
                if index < 26
                else f"s{index + 1}"
            )
            for index, community in enumerate(local_communities)
        }
        child_ids = np.asarray(
            [f"{parent}_{suffixes[value]}" for value in local_labels],
            dtype=object,
        )
        refined_labels[parent_indices] = child_ids
        ####################
        base_row.update(
            {
                "n_child_clusters": len(local_communities),
                "child_clusters": ";".join(sorted(pd.unique(child_ids))),
                "status": "refined",
            }
        )
        refinement_rows.append(base_row)

    adata.obs["Auto_manual_cluster_parent"] = parent_labels
    adata.obs["Auto_manual_cluster"] = refined_labels
    adata.obs["Auto_manual_cluster_refined"] = refined_labels != parent_labels
    return pd.DataFrame(refinement_rows)
####################


def cluster_and_annotate(adata, leiden_resolution):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    marker_scores = normalised_marker_scores(adata)
    ####################
    # Fit the standardized marker-score scale once across all observations.
    # Targeted child clusters only receive new means on this fixed scale; they
    # never alter its centre or standard deviation.
    marker_score_sd = marker_scores.std(axis=0, ddof=0).replace(0, 1)
    standardised_marker_scores = (
        marker_scores - marker_scores.mean(axis=0)
    ) / marker_score_sd
    ####################
    for column in marker_scores:
        adata.obs[f"{column}_score"] = marker_scores[column].to_numpy()

    graph_adata = adata.copy()
    sc.pp.filter_genes(graph_adata, min_cells=GRAPH_MIN_GENE_CELLS)
    sc.pp.highly_variable_genes(
        graph_adata,
        flavor="seurat",
        n_top_genes=min(3000, graph_adata.n_vars),
        span=0.3,
    )
    if int(graph_adata.var["highly_variable"].sum()) < 20:
        raise ValueError("Too few highly variable genes for clustering")
    graph_adata = graph_adata[:, graph_adata.var["highly_variable"]].copy()
    sc.pp.regress_out(graph_adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(graph_adata, max_value=GRAPH_SCALE_MAX)
    n_comps = min(50, graph_adata.n_vars - 1, graph_adata.n_obs - 1)
    sc.tl.pca(graph_adata, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(
        graph_adata,
        n_neighbors=min(15, graph_adata.n_obs - 1),
        n_pcs=min(40, n_comps),
        metric="cosine",
    )
    sc.tl.leiden(
        graph_adata,
        resolution=leiden_resolution,
        key_added="Auto_manual_cluster",
        flavor="igraph",
        directed=False,
    )
    sc.tl.umap(
        graph_adata,
        random_state=0,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
    )
    adata.obs["Auto_manual_cluster"] = (
        graph_adata.obs["Auto_manual_cluster"].astype(str).to_numpy()
    )
    adata.obs["Auto_manual_cluster_parent"] = adata.obs[
        "Auto_manual_cluster"
    ].astype(str)
    adata.obs["Auto_manual_cluster_refined"] = False
    adata.obsm["X_umap"] = graph_adata.obsm["X_umap"].copy()
    adata.obs["Auto_umap_display_pass_qc"] = True

    ####################
    # The initial calls record ranked candidates rejected before the first raw
    # score passes. Eligible parents are locally reclustered once at resolution
    # 1.0, then final evidence and labels are recomputed from scratch.
    _, initial_type_evidence = cluster_marker_support(
        adata,
        "Auto_manual_cluster",
        marker_scores,
        standardised_marker_scores,
    )
    _, initial_calls = choose_cluster_labels(initial_type_evidence)
    refinement_audit = targeted_refine_weak_top_clusters(
        adata,
        marker_scores,
        initial_calls,
    )
    adata.uns["Auto_targeted_refinement_audit"] = refinement_audit
    ####################

    gene_evidence, type_evidence = cluster_marker_support(
        adata,
        "Auto_manual_cluster",
        marker_scores,
        standardised_marker_scores,
    )
    labels, selected = choose_cluster_labels(type_evidence)
    selected["Auto_annotation_leiden_resolution"] = leiden_resolution
    parent_lookup = (
        adata.obs[["Auto_manual_cluster", "Auto_manual_cluster_parent"]]
        .drop_duplicates()
        .set_index("Auto_manual_cluster")["Auto_manual_cluster_parent"]
    )
    selected["Auto_annotation_parent_cluster"] = selected[
        "Auto_manual_cluster"
    ].map(parent_lookup)
    selected["Auto_annotation_targeted_refinement"] = selected[
        "Auto_manual_cluster"
    ].ne(selected["Auto_annotation_parent_cluster"])
    adata.obs["Auto_annotation_celltype"] = (
        adata.obs["Auto_manual_cluster"].map(labels).to_numpy()
    )
    selected = selected.set_index("Auto_manual_cluster")
    for column in selected.columns:
        if column == "Auto_annotation_celltype":
            continue
        adata.obs[column] = adata.obs["Auto_manual_cluster"].map(
            selected[column]
        ).to_numpy()

    adata.obs["Auto_annotation_celltype_pre_filter"] = adata.obs[
        "Auto_annotation_celltype"
    ]
    adata.obs["Auto_annotation_coexpression"] = (
        "not_applied_cluster_score_doublets_retained"
    )
    adata.obs["Auto_annotation_active_marker_types"] = ""

    type_evidence = type_evidence.merge(
        selected.reset_index(),
        on="Auto_manual_cluster",
        how="left",
    )
    gene_evidence = gene_evidence.merge(
        selected.reset_index()[
            [
                "Auto_manual_cluster",
                "Auto_annotation_celltype",
                "Auto_annotation_evidence_stage",
            ]
        ],
        on="Auto_manual_cluster",
        how="left",
    )
    return adata, type_evidence, gene_evidence, selected.reset_index()


def rctd_non_singlet_label(frame):
    first = frame["first_type"].fillna("unresolved").astype(str)
    second = frame.get(
        "second_type", pd.Series("", index=frame.index, dtype=object)
    ).fillna("").astype(str)
    return np.where(second.ne("") & second.ne("nan"), first + "|" + second, first)


def prepare_sample(row, method, rctd_dir):
    input_dir = Path(row[f"{method}_input"] if method == "segmented" else row["binned_input"])
    adata = read_counts(input_dir)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    if method == "binned":
        rctd_path = (
            rctd_dir
            / "tables"
            / f"Auto_{row['sample']}_binned_rctd_annotations.csv.gz"
        )
        if not rctd_path.exists():
            raise FileNotFoundError(f"Missing RCTD table: {rctd_path}")
        rctd = pd.read_csv(rctd_path, compression="gzip")
        rctd["barcode"] = rctd["barcode"].astype(str)
        rctd = rctd.set_index("barcode", drop=False)
        missing = ~rctd.index.isin(adata.obs_names)
        if missing.any():
            raise ValueError(
                f"{int(missing.sum())} RCTD bins are absent from {input_dir}"
            )
        ####################
        # RCTD remains an auditable singlet/doublet call but is no longer an
        # exclusion gate. Every >=100-UMI RCTD bin enters the same custom
        # score and graph annotation used for segmented cells.
        analysed_barcodes = rctd.index
        ####################
        coordinates = read_binned_coordinates(input_dir, analysed_barcodes)
        analysed = adata[analysed_barcodes].copy()
        ####################
        rctd_columns = {
            "spot_class": "Auto_rctd_spot_class",
            "first_type": "Auto_rctd_first_type",
            "second_type": "Auto_rctd_second_type",
            "Auto_rctd_is_singlet": "Auto_rctd_is_singlet",
        }
        for source, destination in rctd_columns.items():
            if source in rctd.columns:
                analysed.obs[destination] = rctd.loc[
                    analysed.obs_names, source
                ].to_numpy()
        return analysed, coordinates, None, input_dir
        ####################

    keep = (adata.obs["total_counts"] >= MIN_COUNTS) & (
        adata.obs["pct_counts_mt"] < MAX_MT_PERCENT
    )
    analysed = adata[keep].copy()
    if analysed.n_obs < MARKER_MIN_CLUSTER_CELLS:
        raise ValueError(
            f"Only {analysed.n_obs} segmented cells passed QC for {row['sample']}"
        )
    coordinates = read_segmented_coordinates(input_dir, analysed.obs_names)
    return analysed, coordinates, None, input_dir


def build_annotation(adata, coordinates, sample, method):
    obs = adata.obs.copy()
    obs.insert(0, "barcode", obs.index.astype(str))
    obs.insert(1, "sample", sample)
    obs.insert(2, "method", method)
    obs["Auto_annotation_method"] = (
        "unified_umap_ranked_standardised_raw_threshold"
    )
    obs["Auto_annotation_pass_doublet_filter"] = True
    obs["Auto_annotation_doublet_status"] = np.where(
        obs["Auto_annotation_celltype"].astype(str).str.contains(
            "|",
            regex=False,
        ),
        "cluster_score_doublet_retained",
        "singlet",
    )
    obs["Auto_annotation_keep_epithelial"] = (
        obs["Auto_annotation_pass_doublet_filter"]
        & obs["Auto_annotation_celltype"].eq("epithelial")
    )
    obs["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    obs["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    coordinates = coordinates.copy()
    coordinates["barcode"] = coordinates.index.astype(str)
    return obs.reset_index(drop=True).merge(
        coordinates.reset_index(drop=True),
        on="barcode",
        how="left",
        validate="one_to_one",
    )


def plot_annotations(
    annotation,
    sample,
    method,
    figure_dir,
    output_name=None,
    title=None,
    point_size_reference_n=None,
):
    plot_data = annotation.loc[
        annotation["Auto_annotation_pass_doublet_filter"].astype(bool)
    ].copy()
    labels = sorted(plot_data["Auto_annotation_celltype"].dropna().unique())
    label_counts = plot_data["Auto_annotation_celltype"].value_counts()
    colours = {
        label: celltype_colour(label)
        for label in labels
    }
    n_points = max(
        point_size_reference_n
        if point_size_reference_n is not None
        else len(plot_data),
        1,
    )
    point_size = max(1.8, min(7.0, 80000 / n_points))
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(19, 9),
        gridspec_kw={"width_ratios": [1, 1, 0.28]},
        constrained_layout=True,
    )
    panels = [
        ("pxl_col_in_fullres", "pxl_row_in_fullres", "Spatial", True),
        ("UMAP_1", "UMAP_2", "UMAP", False),
    ]
    for axis, (x_col, y_col, panel_title, reverse_y) in zip(axes[:2], panels):
        panel_data = plot_data
        if panel_title == "UMAP":
            panel_data = plot_data.loc[
                np.isfinite(plot_data[x_col]) & np.isfinite(plot_data[y_col])
            ]
        for label in labels:
            selected = panel_data["Auto_annotation_celltype"].eq(label)
            axis.scatter(
                panel_data.loc[selected, x_col],
                panel_data.loc[selected, y_col],
                s=point_size,
                alpha=0.78,
                linewidths=0,
                color=colours[label],
                rasterized=True,
            )
        axis.set_title(panel_title, fontsize=15)
        axis.set_aspect("equal", adjustable="datalim")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
        if reverse_y:
            axis.invert_yaxis()
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=15,
            markerfacecolor=colours[label],
            markeredgecolor="none",
            label=f"{label} (n={int(label_counts[label]):,})",
        )
        for label in labels
    ]
    axes[2].axis("off")
    axes[2].legend(
        handles=handles,
        loc="center left",
        frameon=False,
        title="Cell type",
        fontsize=13,
        title_fontsize=14,
    )
    fig.suptitle(
        title if title is not None else f"{sample} | {method}",
        fontsize=17,
    )
    base = figure_dir / (
        output_name
        if output_name is not None
        else f"Auto_{sample}_{method}_annotation_diagnostics"
    )
    fig.savefig(base.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_assignment_umap(annotation, sample, method, figure_dir):
    plot_data = annotation.loc[
        annotation["Auto_annotation_pass_doublet_filter"].astype(bool)
        & np.isfinite(annotation["UMAP_1"])
        & np.isfinite(annotation["UMAP_2"])
    ].copy()
    labels = sorted(plot_data["Auto_annotation_celltype"].dropna().unique())
    output_path = (
        figure_dir
        / f"Auto_{sample}_{method}_cluster_assignment_umap.pdf"
    )
    with PdfPages(output_path) as pdf:
        for label in labels:
            selected = plot_data.loc[
                plot_data["Auto_annotation_celltype"].eq(label)
            ]
            fig, axis = plt.subplots(
                figsize=(11, 9),
                constrained_layout=True,
            )
            axis.scatter(
                plot_data["UMAP_1"],
                plot_data["UMAP_2"],
                s=1.2,
                alpha=0.10,
                linewidths=0,
                color="#BDBDBD",
                rasterized=True,
            )
            axis.scatter(
                selected["UMAP_1"],
                selected["UMAP_2"],
                s=4.0,
                alpha=0.82,
                linewidths=0,
                color=celltype_colour(label),
                rasterized=True,
            )
            for cluster, cluster_data in selected.groupby(
                "Auto_manual_cluster", observed=True
            ):
                axis.text(
                    cluster_data["UMAP_1"].median(),
                    cluster_data["UMAP_2"].median(),
                    f"C{cluster}\nn={len(cluster_data):,}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                    color="black",
                    bbox={
                        "boxstyle": "round,pad=0.18",
                        "facecolor": "white",
                        "edgecolor": "#636363",
                        "linewidth": 0.4,
                        "alpha": 0.86,
                    },
                )
            axis.set_title(
                f"{sample} | {method} | {label}"
                f" | n={len(selected):,}",
                fontsize=16,
            )
            axis.set_aspect("equal", adjustable="datalim")
            axis.set_xticks([])
            axis.set_yticks([])
            for spine in axis.spines.values():
                spine.set_visible(False)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


def plot_cluster_score_heatmaps(
    type_evidence,
    cluster_calls,
    sample,
    method,
    figure_dir,
):
    call_table = cluster_calls.copy()
    call_table["Auto_manual_cluster"] = call_table[
        "Auto_manual_cluster"
    ].astype(str)
    call_table = call_table.sort_values(
        [
            "Auto_annotation_celltype",
            "Auto_annotation_cluster_n",
            "Auto_manual_cluster",
        ],
        ascending=[True, False, True],
    )
    cluster_order = call_table["Auto_manual_cluster"].tolist()
    celltype_order = list(MARKERS)
    display_labels = {
        row.Auto_manual_cluster: (
            f"C{row.Auto_manual_cluster} | {row.Auto_annotation_celltype}"
            f" | n={int(row.Auto_annotation_cluster_n)}"
        )
        for row in call_table.itertuples()
    }
    evidence = type_evidence.copy()
    evidence["Auto_manual_cluster"] = evidence[
        "Auto_manual_cluster"
    ].astype(str)
    raw = evidence.pivot(
        index="cell_type",
        columns="Auto_manual_cluster",
        values="cluster_normalised_marker_score",
    ).reindex(index=celltype_order, columns=cluster_order)
    standardised = evidence.pivot(
        index="cell_type",
        columns="Auto_manual_cluster",
        values="cluster_standardised_marker_score",
    ).reindex(index=celltype_order, columns=cluster_order)
    raw_threshold_pass = raw > RAW_NORMALISED_SCORE_THRESHOLD

    height = max(9.0, 0.48 * len(celltype_order) + 3.0)
    pdf_path = (
        figure_dir
        / f"Auto_{sample}_{method}_cluster_marker_score_heatmaps.pdf"
    )
    with PdfPages(pdf_path) as pdf:
        for matrix, title, colour_map, symmetric, legend_title in [
            (
                raw,
                "Raw normalized marker score",
                "Blues",
                False,
                "Mean log1p CP10K",
            ),
            (
                standardised,
                "Standardized marker score",
                "RdBu_r",
                True,
                "Mean within-sample z-score",
            ),
        ]:
            full_values = matrix.to_numpy(dtype=float)
            finite = full_values[np.isfinite(full_values)]
            if symmetric:
                limit = max(float(np.quantile(np.abs(finite), 0.99)), 0.1)
                vmin, vmax = -limit, limit
            else:
                vmin = 0.0
                vmax = max(float(np.quantile(finite, 0.99)), 0.01)
            cluster_pages = [
                cluster_order[index:index + CLUSTERS_PER_HEATMAP_PAGE]
                for index in range(0, len(cluster_order), CLUSTERS_PER_HEATMAP_PAGE)
            ]
            for page_index, page_clusters in enumerate(cluster_pages, start=1):
                page_matrix = matrix.loc[:, page_clusters]
                page_threshold_pass = raw_threshold_pass.loc[:, page_clusters]
                values = page_matrix.to_numpy(dtype=float)
                width = max(15.0, 0.48 * len(page_clusters) + 5.0)
                fig, axis = plt.subplots(
                    figsize=(width, height),
                    constrained_layout=True,
                )
                image = axis.imshow(
                    values,
                    aspect="auto",
                    interpolation="nearest",
                    cmap=colour_map,
                    vmin=vmin,
                    vmax=vmax,
                )
                support_y, support_x = np.where(
                    page_threshold_pass.to_numpy(dtype=bool)
                )
                axis.scatter(
                    support_x,
                    support_y,
                    s=18,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=0.65,
                )
                axis.set_xticks(np.arange(len(page_clusters)))
                axis.set_xticklabels(
                    [display_labels[cluster] for cluster in page_clusters],
                    rotation=90,
                    ha="center",
                    fontsize=8,
                )
                axis.set_yticks(np.arange(len(celltype_order)))
                axis.set_yticklabels(celltype_order, fontsize=10)
                axis.set_xlabel(
                    "Cluster | final assignment | cell count",
                    fontsize=12,
                )
                axis.set_ylabel("")
                axis.set_title(
                    f"{sample} | {method} | {title}"
                    f" | page {page_index}/{len(cluster_pages)}",
                    fontsize=15,
                )
                colour_bar = fig.colorbar(
                    image,
                    ax=axis,
                    fraction=0.022,
                    pad=0.015,
                )
                colour_bar.set_label(legend_title, fontsize=10)
                axis.text(
                    1.0,
                    -0.19,
                    (
                        "Black outline: raw normalized score > "
                        f"{RAW_NORMALISED_SCORE_THRESHOLD}"
                    ),
                    transform=axis.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                )
                pdf.savefig(fig, dpi=300, bbox_inches="tight")
                if symmetric:
                    png_name = (
                        f"Auto_{sample}_{method}_cluster_marker_score_heatmap"
                        + ("" if page_index == 1 else f"_page{page_index}")
                        + ".png"
                    )
                    fig.savefig(
                        figure_dir / png_name,
                        dpi=300,
                        bbox_inches="tight",
                    )
                plt.close(fig)


def run_method(
    manifest,
    method,
    output_dir,
    rctd_dir,
    reuse_complete,
    force_samples,
    leiden_resolution_override,
):
    table_dir = output_dir / "tables"
    intermediate_dir = output_dir / "intermediate"
    figure_dir = output_dir / "figures" / "annotation_diagnostics"
    log_dir = output_dir / "logs"
    for directory in [table_dir, intermediate_dir, figure_dir, log_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    log_lines = [f"method={method}"]

    for _, row in manifest.iterrows():
        sample = row["sample"]
        annotation_path = table_dir / f"Auto_{sample}_{method}_cell_annotations.csv.gz"
        evidence_path = (
            table_dir
            / f"Auto_{sample}_{method}_cluster_celltype_evidence.csv"
        )
        marker_support_path = (
            table_dir
            / f"Auto_{sample}_{method}_cluster_marker_support.csv"
        )
        cluster_calls_path = (
            table_dir
            / f"Auto_{sample}_{method}_cluster_score_calls.csv"
        )
        ####################
        refinement_path = (
            table_dir
            / f"Auto_{sample}_{method}_targeted_refinement.csv"
        )
        ####################
        if (
            reuse_complete
            and annotation_path.exists()
            and evidence_path.exists()
            and cluster_calls_path.exists()
            and refinement_path.exists()
            and sample not in force_samples
        ):
            annotation = pd.read_csv(annotation_path, compression="gzip")
            type_evidence = pd.read_csv(evidence_path)
            cluster_calls = pd.read_csv(cluster_calls_path)
            refinement_audit = (
                pd.read_csv(refinement_path)
                if refinement_path.exists()
                else pd.DataFrame()
            )
            log_lines.append(f"sample={sample}; status=reused")
        else:
            analysed, coordinates, excluded, _ = prepare_sample(
                row, method, rctd_dir
            )
            leiden_resolution = (
                leiden_resolution_override
                if leiden_resolution_override is not None
                else (
                    LEIDEN_RESOLUTION_LARGE
                    if analysed.n_obs >= LEIDEN_LARGE_MIN_OBSERVATIONS
                    else LEIDEN_RESOLUTION_SMALL
                )
            )
            analysed, type_evidence, gene_evidence, cluster_calls = (
                cluster_and_annotate(
                    analysed,
                    leiden_resolution,
                )
            )
            ####################
            refinement_audit = analysed.uns[
                "Auto_targeted_refinement_audit"
            ].copy()
            refinement_audit.insert(0, "sample", sample)
            refinement_audit.insert(1, "method", method)
            refinement_audit.to_csv(refinement_path, index=False)
            ####################
            annotation = build_annotation(analysed, coordinates, sample, method)
            if excluded is not None:
                excluded["sample"] = sample
                excluded["method"] = method
                for column in annotation.columns:
                    if column not in excluded.columns:
                        excluded[column] = np.nan
                annotation = pd.concat(
                    [annotation, excluded[annotation.columns]],
                    ignore_index=True,
                )
            annotation.to_csv(annotation_path, index=False, compression="gzip")
            type_evidence.insert(0, "sample", sample)
            type_evidence.insert(1, "method", method)
            gene_evidence.insert(0, "sample", sample)
            gene_evidence.insert(1, "method", method)
            cluster_calls.insert(0, "sample", sample)
            cluster_calls.insert(1, "method", method)
            type_evidence.to_csv(
                evidence_path,
                index=False,
            )
            gene_evidence.to_csv(
                marker_support_path,
                index=False,
            )
            cluster_calls.to_csv(
                cluster_calls_path,
                index=False,
            )
            annotation.loc[
                annotation["Auto_annotation_pass_doublet_filter"].astype(bool),
                ["barcode", "UMAP_1", "UMAP_2"],
            ].to_csv(
                intermediate_dir / f"Auto_{sample}_{method}_annotation_umap.csv.gz",
                index=False,
                compression="gzip",
            )
            log_lines.append(
                f"sample={sample}; status=calculated; n_analysed={analysed.n_obs}"
            )

        plot_annotations(annotation, sample, method, figure_dir)
        plot_cluster_assignment_umap(
            annotation,
            sample,
            method,
            figure_dir,
        )
        plot_cluster_score_heatmaps(
            type_evidence,
            cluster_calls,
            sample,
            method,
            figure_dir,
        )
        passing = annotation["Auto_annotation_pass_doublet_filter"].astype(bool)
        label_counts = (
            annotation.loc[passing, "Auto_annotation_celltype"]
            .astype("string")
            .fillna("unresolved")
            .value_counts()
        )
        if "Auto_annotation_leiden_resolution" in annotation.columns:
            resolution_values = (
                pd.to_numeric(
                    annotation.loc[
                        passing,
                        "Auto_annotation_leiden_resolution",
                    ],
                    errors="coerce",
                )
                .dropna()
                .unique()
            )
        else:
            resolution_values = np.array([], dtype=float)
        summary_rows.append(
            {
                "sample": sample,
                "method": method,
                "n_output": len(annotation),
                "n_pass_doublet_filter": int(passing.sum()),
                "n_rctd_non_singlet_excluded": int((~passing).sum()),
                ####################
                "n_rctd_non_singlet_annotated": int(
                    annotation.get(
                        "Auto_rctd_spot_class",
                        pd.Series("singlet", index=annotation.index),
                    )
                    .astype(str)
                    .ne("singlet")
                    .sum()
                ),
                "n_targeted_refinement_parents": int(
                    refinement_audit.get(
                        "status", pd.Series(dtype=object)
                    ).eq("refined").sum()
                ),
                ####################
                "n_cluster_score_doublet_retained": int(
                    annotation.loc[
                        passing,
                        "Auto_annotation_celltype",
                    ]
                    .astype(str)
                    .str.contains("|", regex=False)
                    .sum()
                ),
                "leiden_resolution": ";".join(
                    f"{value:g}" for value in sorted(resolution_values)
                ),
                "n_epithelial": int(
                    annotation["Auto_annotation_keep_epithelial"]
                    .fillna(False)
                    .astype(bool)
                    .sum()
                ),
                "celltype_counts": ";".join(
                    f"{label}={count}" for label, count in label_counts.items()
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        table_dir / f"Auto_visiumhd_{method}_annotation_summary.csv",
        index=False,
    )
    (log_dir / f"Auto_visiumhd_{method}_annotation_run_summary.txt").write_text(
        "\n".join(log_lines) + "\n"
    )
    return summary
def plot_summary_proportions(summary_df, output_path):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    records = []
    for _, row in summary_df.iterrows():
        sample = row["sample"]
        method = row["method"]
        counts_str = row["celltype_counts"]
        if pd.isna(counts_str) or not counts_str:
            continue
        
        total = 0
        parsed = {}
        for item in str(counts_str).split(";"):
            if "=" in item:
                ct, count = item.split("=")
                parsed[ct] = int(count)
                total += int(count)
        
        for ct, count in parsed.items():
            records.append({
                "sample": sample,
                "method": method,
                "celltype": ct,
                "count": count,
                "proportion": count / total if total > 0 else 0
            })
            
    long_df = pd.DataFrame(records)
    if long_df.empty:
        print("No celltype counts found for plotting summary.")
        return
        
    samples = sorted(long_df["sample"].unique())
    
    with PdfPages(output_path) as pdf:
        for sample in samples:
            sample_df = long_df[long_df["sample"] == sample].copy()
            pivot_df = sample_df.pivot(index="method", columns="celltype", values="proportion").fillna(0)
            
            available_methods = [m for m in ["binned", "segmented"] if m in pivot_df.index]
            pivot_df = pivot_df.loc[available_methods]
            
            cols = sorted(pivot_df.columns)
            if "unresolved" in cols:
                cols.remove("unresolved")
                cols.append("unresolved")
            pivot_df = pivot_df[cols]
            
            fig, ax = plt.subplots(figsize=(8, 6))
            bottoms = np.zeros(len(pivot_df))
            for ct in pivot_df.columns:
                values = pivot_df[ct].values
                color = celltype_colour(ct)
                
                ax.bar(pivot_df.index, values, bottom=bottoms, label=ct, color=color, edgecolor="white", width=0.6)
                bottoms += values
                
            ax.set_ylabel("Proportion of Cells", fontsize=12)
            ax.set_title(f"Cell Type Proportions - {sample}", fontsize=14)
            ax.set_ylim(0, 1.05)
            
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], title="Cell Type", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)


def main():
    global RAW_NORMALISED_SCORE_THRESHOLD
    global STRUCTURAL_DOUBLET_STANDARDISED_GAP

    args = parse_args()
    if args.leiden_resolution is not None and args.leiden_resolution <= 0:
        raise ValueError("--leiden-resolution must be positive")
    if args.raw_score_threshold < 0:
        raise ValueError("--raw-score-threshold cannot be negative")
    if args.structural_doublet_gap < 0:
        raise ValueError("--structural-doublet-gap cannot be negative")
    RAW_NORMALISED_SCORE_THRESHOLD = args.raw_score_threshold
    STRUCTURAL_DOUBLET_STANDARDISED_GAP = args.structural_doublet_gap

    manifest_path = Path(args.manifest)
    output_dir = Path(args.output_dir)
    rctd_dir = Path(args.rctd_dir)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    required = {"sample", "binned_input", "segmented_input"}
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest lacks columns: {sorted(missing)}")
    if manifest["sample"].duplicated().any():
        raise ValueError("Manifest sample names must be unique")
    if args.samples:
        unknown_samples = sorted(set(args.samples).difference(manifest["sample"]))
        if unknown_samples:
            raise ValueError(
                f"Requested samples are absent from manifest: {unknown_samples}"
            )
        manifest = manifest.loc[manifest["sample"].isin(args.samples)].copy()
    for column in ["binned_input", "segmented_input"]:
        absent = [path for path in manifest[column] if not Path(path).is_dir()]
        if absent:
            raise FileNotFoundError(f"Missing {column} paths: {absent}")

    summaries = [
        run_method(
            manifest,
            method,
            output_dir,
            rctd_dir,
            args.reuse_complete,
            set(args.force_samples),
            args.leiden_resolution,
        )
        for method in args.methods
    ]
    combined = pd.concat(summaries, ignore_index=True)
    combined.to_csv(
        output_dir / "tables" / "Auto_visiumhd_annotation_summary.csv",
        index=False,
    )
    parameters = pd.DataFrame(
        [
            {
                "minimum_counts": MIN_COUNTS,
                "maximum_mitochondrial_percent": MAX_MT_PERCENT,
                "graph_minimum_cells_per_gene": GRAPH_MIN_GENE_CELLS,
                "graph_regressed_covariates": "total_counts;pct_counts_mt",
                "graph_scale_maximum": GRAPH_SCALE_MAX,
                "umap_minimum_distance": UMAP_MIN_DIST,
                "umap_spread": UMAP_SPREAD,
                "normalisation": "log1p_CP10K",
                "highly_variable_genes": 3000,
                "neighbours": 15,
                "leiden_resolution_rule": (
                    "10_if_at_least_1000_analysed_observations_else_6"
                    if args.leiden_resolution is None
                    else "command_line_override"
                ),
                "leiden_resolution_large": LEIDEN_RESOLUTION_LARGE,
                "leiden_resolution_small": LEIDEN_RESOLUTION_SMALL,
                "leiden_large_minimum_observations": (
                    LEIDEN_LARGE_MIN_OBSERVATIONS
                ),
                "leiden_resolution_override": args.leiden_resolution,
                ####################
                "binned_rctd_use": (
                    "metadata_only_all_rctd_bins_enter_custom_annotation"
                ),
                "local_cluster_refinement": (
                    "one_pass_rejected_rank_candidate_expression_graph_refinement"
                ),
                "targeted_refinement_cell_score": (
                    TARGETED_REFINEMENT_CELL_SCORE
                ),
                "targeted_refinement_minimum_high_cells": (
                    TARGETED_REFINEMENT_MIN_HIGH_CELLS
                ),
                "targeted_refinement_local_resolution": (
                    TARGETED_REFINEMENT_LOCAL_RESOLUTION
                ),
                ####################
                "whole_transcriptome_dge": "none",
                "marker_min_cluster_cells": MARKER_MIN_CLUSTER_CELLS,
                "marker_minimum_positive_cells": MARKER_MIN_POSITIVE_CELLS,
                "marker_support_use": "diagnostic_only_not_used_for_assignment",
                "raw_normalised_score_threshold": (
                    RAW_NORMALISED_SCORE_THRESHOLD
                ),
                "assignment_rule": (
                    "rank_all_celltypes_by_cluster_mean_standardised_score_"
                    "then_select_first_with_raw_normalised_score_above_threshold"
                ),
                "structural_doublet_types": "epithelial;fibroblast",
                "structural_doublet_standardised_gap": (
                    STRUCTURAL_DOUBLET_STANDARDISED_GAP
                ),
                "structural_doublet_rule": (
                    "when_structural_is_selected_and_best_nonstructural_"
                    "passing_type_is_within_gap_label_nonstructural|structural"
                ),
                "clustering_purpose": (
                    "size_aware_high_resolution_clusters_for_score_aggregation"
                ),
                "clustering_and_umap": (
                    "same_regressed_scaled_pca_cosine_neighbor_graph"
                ),
                "targeted_reclustering": (
                    "one_pass_resolution_1_all_local_children_retained"
                ),
                "cell_level_coexpression_filter": "none",
            }
        ]
    )
    parameters.to_csv(
        output_dir / "tables" / "Auto_visiumhd_annotation_parameters.csv",
        index=False,
    )
    (output_dir / "logs" / "Auto_visiumhd_superseded_outputs.txt").write_text(
        "The current annotation does not use DGE or blanket local refinement.\n"
        "Targeted local expression-graph refinement is used once when a "
        "raw-rejected ranked candidate has at least 20 observations above "
        "individual score 1.0, and is audited in "
        "*_targeted_refinement.csv.\n"
        "Pre-existing files ending cluster_marker_dge.csv or "
        "cluster_top_dge.csv.gz are retained for file safety but are "
        "superseded.\n"
        "Use cluster_marker_support.csv, cluster_celltype_evidence.csv, and "
        "cluster_score_calls.csv.\n"
    )
    canonical_output_dir = WD / "ref_outs" / "visium_hd_outs"
    if output_dir.resolve() == canonical_output_dir.resolve():
        update_dir = WD / "updates" / "new_updates" / "summaries"
        update_dir.mkdir(parents=True, exist_ok=True)
        combined.to_csv(
            update_dir / "visiumhd_final_annotation_summary.csv",
            index=False,
        )
    
    # Generate the multi-page PDF summary for celltype proportions
    plot_summary_proportions(
        combined,
        output_dir / "figures" / "annotation_diagnostics" / "Auto_visiumhd_celltype_proportions.pdf"
    )

if __name__ == "__main__":
    main()
####################
