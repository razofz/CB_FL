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

# +
# %config InlineBackend.figure_format = 'retina'

import glob

import fl_utils
import pandas as pd
import scarf

# + tags=[]
fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
figures_dir = out_dir + "/Figures"
external_dir = "data/external"
roy_dir = external_dir + "/Roy_et_al"
roy_analysis_dir = out_dir + "/Roy_et_al_analysis"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_roy": roy_analysis_dir + "/fetal_roy.zarr",
}

out_files = {}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(roy_dir)
# -

# ---

# +
ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")

# ds_yBM = scarf.DataStore('./_zarr_files/RNA_yBM_hpc_merged', default_assay='RNA')
# ds_FL = scarf.DataStore('./_zarr_files/RNA_FL_hpc_merged', default_assay='RNA')
# -

ds_roy = scarf.DataStore(in_files["zarr_roy"])
# ds_roy = scarf.DataStore('./_zarr_files/CB_fetal_roy.zarr')

# +
color_palette_BM = {}
x = pd.crosstab(
    ds_yBM.cells.fetch("RNA_cluster_colors"),
    ds_yBM.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_BM[i] = x[i][x[i] > 0].index[0]

color_palette_FL = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"),
    ds_FL.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_FL[i] = x[i][x[i] > 0].index[0]
# -

ds_FL.run_mapping(
    target_assay=ds_roy.RNA,
    target_name="roy_data",
    target_feat_key="I__hvgs",
    save_k=9,
)

ds_yBM.run_mapping(
    target_assay=ds_roy.RNA,
    target_name="roy_data",
    target_feat_key="I__hvgs",
    save_k=9,
)
