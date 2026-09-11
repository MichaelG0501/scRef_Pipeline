#!/usr/bin/env python
####################
# Analysis registry:
#   Status: active
#   Script: analysis/spatial/visium_hd_postfilter_diagnostics.py
#   Description: Replot final binned post-filter annotations with the exact
#     production annotation-diagnostic layout and point styling.
#   Methodology:
#     analysis/methodology/spatial/visium_hd_binned_filter_malignancy_methodology.md
#   Inputs:
#     analysis/shared/visium_hd_celltype_colours.tsv
#     ref_outs/visium_hd_outs/post_annotation_filter/tables/
#       Auto_<sample>_binned_filtered_annotations.csv.gz
#   Outputs:
#     figures/: matched post-filter PDF/PNG diagnostics under
#       ref_outs/visium_hd_outs/figures/annotation_diagnostics/after_filtering/
#   Cache/replot: plot-only; always rebuilt from live filtered tables.
#   Run: python analysis/spatial/visium_hd_postfilter_diagnostics.py
#   Environment: /rds/general/user/sg3723/home/miniforge3/envs/jupyter
####################

####################
from pathlib import Path

import pandas as pd

from visium_hd_celltype_annotation import plot_annotations


WD = Path("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline")
SAMPLES = ("SUR1231", "FFPEA1", "FFPED1")
FILTER_TABLE_DIR = (
    WD / "ref_outs" / "visium_hd_outs" / "post_annotation_filter" / "tables"
)
FIGURE_DIR = (
    WD
    / "ref_outs"
    / "visium_hd_outs"
    / "figures"
    / "annotation_diagnostics"
    / "after_filtering"
)
def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for sample in SAMPLES:
        input_path = (
            FILTER_TABLE_DIR
            / f"Auto_{sample}_binned_filtered_annotations.csv.gz"
        )
        if not input_path.exists():
            raise FileNotFoundError(
                f"Missing filtered annotation table: {input_path}"
            )
        annotation = pd.read_csv(input_path, compression="gzip")
        required = {
            "Auto_postfilter_keep",
            "Auto_postfilter_celltype",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
            "UMAP_1",
            "UMAP_2",
        }
        missing = required.difference(annotation.columns)
        if missing:
            raise ValueError(
                f"{input_path.name} is missing columns: {sorted(missing)}"
            )

        keep_values = annotation["Auto_postfilter_keep"]
        keep_mask = (
            keep_values
            if pd.api.types.is_bool_dtype(keep_values)
            else keep_values.astype(str).str.lower().isin({"true", "t", "1"})
        )
        retained = annotation.loc[keep_mask].copy()
        retained["Auto_annotation_pass_doublet_filter"] = True
        retained["Auto_annotation_celltype"] = retained[
            "Auto_postfilter_celltype"
        ].astype(str)
        plot_annotations(
            retained,
            sample,
            "binned",
            FIGURE_DIR,
            output_name=(
                f"Auto_{sample}_binned_annotation_diagnostics_after_filtering"
            ),
            title=f"{sample} | binned | after filtering",
            point_size_reference_n=len(annotation),
        )


if __name__ == "__main__":
    main()
####################
