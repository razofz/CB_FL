# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import fl_utils
import pandas as pd
import scarf
from scarf.plots import plot_heatmap

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"
zarr_dir = out_dir + "/zarr_files"

# +
in_files = {
    "zarr_merged": zarr_dir + "/RNA_yBM_and_FL_hpc_merged",
}

out_files = {
    "plot_hox_heatmap": figures_dir + "/hox_heatmap.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(out_adt_matrices_dir)
# -

# ---

ds = scarf.DataStore(in_files["zarr_merged"], default_assay="RNA")
# ds = scarf.DataStore("_zarr_files/RNA_yBM_and_FL_hpc_merged", default_assay="RNA")

ds.cells.head().columns

ds.cells.to_pandas_dataframe(columns=["names"]).value_counts()

ds.cells.insert(
    column_name="is_mapped",
    values=ds.cells.fetch("clusters") != "bm_NA",
    overwrite=True,
    fill_value=False,
)

ds.RNA.feats.insert(
    "is_hox",
    ds.RNA.feats.index_to_bool(
        ds.RNA.feats.get_index_by(ds.RNA.feats.grep("^HOX"), "names")
    ),
    fill_value=False,
    overwrite=True,
)
c = ds.RNA.feats.fetch("nCells", key="is_hox") > 10
ds.RNA.feats.insert(
    column_name="is_filtered_hox",
    values=c,
    key="is_hox",
    fill_value=False,
    overwrite=True,
)
ds.RNA.z.create_dataset(
    name="markers/is_mapped__clusters/bm_HSC/names",
    data=ds.RNA.feats.active_index(key="is_filtered_hox"),
    overwrite=True,
)

ds.plot_marker_heatmap(
    group_key="clusters",
    cell_key="is_mapped",
    topn=len(c[c]),
    savename=out_files["plot_hox_heatmap"],
    # savename="./images/hox_heatmap.svg",
)

ds.RNA.feats.fetch_all("ids")[ds.RNA.feats.fetch_all("is_filtered_hox")]
