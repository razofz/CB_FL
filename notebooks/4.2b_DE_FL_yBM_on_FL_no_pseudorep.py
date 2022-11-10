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

# + tags=[]
# %config InlineBackend.figure_format = 'retina'

import scarf
import json
import pandas as pd
import numpy as np
import fl_utils
from distutils.dir_util import copy_tree


# + tags=[]
fl_utils.set_project_path()

# + tags=[]
raw_dir = "data/raw/paper_specific/from_seurat"
zarr_dir = "data/processed/notebooks/zarr_files"
preds_dir = "data/processed/notebooks/mapping_predictions"
out_dir = "data/processed/notebooks"
deseq2_dir = out_dir + "/DEseq2/FLcoreNoPseudorep"
deseq2_clw_dir = deseq2_dir + "/cluster_wise"

# + tags=[]
fl_utils.create_dir(deseq2_dir)
fl_utils.create_dir(deseq2_clw_dir)
fl_utils.create_dir(zarr_dir + "/FLcoreNoPseudorep")

# + tags=[]
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "FL_target_preds": preds_dir + "/RNA_FL_target_preds.json",
    # "metadata_FL": 'data/processed/notebooks/FL_seurat_metadata.csv',
    "cluster_markers_FL": raw_dir + "/FL_combined_cluster_markers.csv",
}

out_files = {
    "zarr_yBM_FLcoreNoPseudorep": zarr_dir + "/FLcoreNoPseudorep/RNA_yBM_hpc_merged",
    "zarr_FL_FLcoreNoPseudorep": zarr_dir + "/FLcoreNoPseudorep/RNA_FL_hpc_merged",
    "deseq2_dir": deseq2_dir,
    "deseq2_clw_dir": deseq2_dir + "/cluster_wise",
}

fl_utils.create_dir(out_files["zarr_yBM_FLcoreNoPseudorep"])
fl_utils.create_dir(out_files["zarr_FL_FLcoreNoPseudorep"])

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# + tags=[]
threads = 16
# -

# ---

# + tags=[]
preds = json.load(open(in_files["FL_target_preds"]))
preds.keys()

# + tags=[]
dfs = []
# -

_ = copy_tree(src=in_files["zarr_yBM"], dst=out_files["zarr_yBM_FLcoreNoPseudorep"],
        dry_run=0, verbose=0, update=1)
_ = copy_tree(src=in_files["zarr_FL"], dst=out_files["zarr_FL_FLcoreNoPseudorep"],
        dry_run=0, verbose=0, update=1)

ds = scarf.DataStore(out_files["zarr_FL_FLcoreNoPseudorep"], nthreads=threads)

ds.cells.insert(column_name="sample", values=[x.split("__")[0] for x in ds.cells.fetch_all("ids")], overwrite=True)

x = np.char.add(
    ds.cells.fetch_all("sample"),
    ["_" + x for x in ds.cells.fetch_all("RNA_cluster_names")],
)
ds.cells.insert("split_sample_clusts", x, fill_value="NA", overwrite=True)
df = ds.make_bulk(
    from_assay="RNA", group_key="split_sample_clusts", pseudo_reps=1
)
df = df.drop(
    columns=[
        "fl1_NA_Rep1",
        # "fl1_NA_Rep2",
        "fl2_NA_Rep1",
        # "fl2_NA_Rep2",
        "names",
    ]
)
df.columns = [f"FL_hpc_{x}" for x in df.columns]
dfs.append(df.T)

# -

ds = scarf.DataStore(out_files["zarr_yBM_FLcoreNoPseudorep"], nthreads=threads)

ds.cells.head()

ds.cells.insert(column_name="sample", values=[x.split("__")[0] for x in ds.cells.fetch_all("ids")], overwrite=True)

ds.cells.head()

# + tags=[]
x = (
    pd.Series(preds["yBM_HPC"])
    .reindex(ds.cells.fetch_all("ids"))
    .fillna("NA")
    .values
)
ds.cells.insert("predicted_celltype_FL", x, fill_value="NA", overwrite=True)
x = np.char.add(
    ds.cells.fetch_all("sample"),
    ["_" + x for x in ds.cells.fetch_all("predicted_celltype_FL")],
)
ds.cells.insert("split_sample_pred_clusts", x, fill_value="NA", overwrite=True)
df = ds.make_bulk(
    from_assay="RNA", group_key="split_sample_pred_clusts", pseudo_reps=1
)
df = df.drop(
    columns=[
        "ybm1_NA_Rep1",
        # "ybm1_NA_Rep2",
        "ybm2_NA_Rep1",
        # "ybm2_NA_Rep2",
        "names",
    ]
)
df.columns = [f"yBM_hpc_{x}" for x in df.columns]
dfs.append(df.T)

dfs = pd.concat(dfs).fillna(0).T
dfs.columns = [x.replace("-", ".") for x in dfs.columns]

# + tags=[]

metadata = pd.DataFrame(
    [
        (x, x.split("_", 1)[0], x.split("_")[3], x.rsplit("_", 1)[1])
        for x in dfs.columns
    ],
    columns=["sample_id", "sample", "cluster_name", "pseudo_rep"],
).set_index("sample_id")

for i in metadata.cluster_name.unique():
    print(i)
    idx = metadata.cluster_name == i
    x = dfs[metadata[idx].index]
    x.to_csv(f"{deseq2_clw_dir}/{i}_cts_all_hpc.csv")
    metadata[idx].to_csv(f"{deseq2_clw_dir}/{i}_coldata_all_hpc.csv")

# +
ds = scarf.DataStore(out_files["zarr_FL_FLcoreNoPseudorep"],
                     nthreads=threads)
df = ds.make_bulk(
    from_assay="RNA", group_key="split_sample_clusts", pseudo_reps=1
)
df = df.drop(
    columns=[
        "fl1_NA_Rep1",
        # "fl1_NA_Rep2",
        "fl2_NA_Rep1",
        # "fl2_NA_Rep2",
        "names",
    ]
)
df.columns = [f"FL_hpc_{x}" for x in df.columns]

df.columns = [x.replace("-", ".") for x in df.columns]
metadata = pd.DataFrame(
    [
        ((x, x.split("_", 1)[0], x.split("_")[3], x.rsplit("_", 1)[1]))
        for x in df.columns
    ],
    columns=["sample_id", "sample", "cluster_name", "pseudo_rep"],
).set_index("sample_id")


df.to_csv(f"{deseq2_dir}/FL_cts_hpc.csv")
metadata.to_csv(f"{deseq2_dir}/FL_coldata.csv")
