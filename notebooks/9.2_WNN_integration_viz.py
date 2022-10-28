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
import glob
import os

import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scarf as scf

# + tags=[]
fl_utils.set_project_path()

# +
threads = 16

out_dir = "data/processed/notebooks"
raw_dir = "data/raw/paper_specific/from_seurat"
figures_dir = "data/processed/notebooks/Figures"
zarr_dir = out_dir + "/zarr_files"

# + tags=[]
in_files = {
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "metadata": raw_dir + "/WNN_FL_combined_metadata.csv",
    "umap": raw_dir + "/WNN_FL_combined_umap.csv",
}

out_files = {
    "umap_wnn_clusters": figures_dir + "/WNN_UMAP_WNN_Clusters.svg",
    "umap_org_clusters": figures_dir + "/WNN_UMAP_org_Clusters.svg",
    "barplot": figures_dir + "/bar_plot_WNN_org_clust_comparison.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

# ---

# +
# %load_ext autotime
# %config InlineBackend.figure_format = 'retina'

plt.style.use("fivethirtyeight")
plt.rcParams["svg.fonttype"] = "none"


def clean_axis(ax, ts=11, ga=0.4):
    ax.xaxis.set_tick_params(labelsize=ts)
    ax.yaxis.set_tick_params(labelsize=ts)
    for i in ["top", "bottom", "left", "right"]:
        ax.spines[i].set_visible(False)
    ax.grid(which="major", linestyle="--", alpha=ga)
    ax.figure.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    return True


# -

metadata = pd.read_csv(in_files["metadata"])
umap = pd.read_csv(in_files["umap"])

ds_FL = scf.DataStore(in_files["zarr_FL"], default_assay="RNA")

metadata["WNN_u1"] = umap["wnnUMAP_1"]
metadata["WNN_u2"] = umap["wnnUMAP_2"] * (
    -1
)  # multipled by -1 to flip umap so HSCs are on top
metadata.index = [x.split(" ")[0] for x in metadata.index]
metadata.index = [
    x.split("_")[0].lower() + "__" + x.split("_")[2] for x in metadata.index
]

metadata = metadata.drop(
    columns=[
        "RNA_snn_res.0.8",
        "seurat_clusters",
        "RNA_snn_res.0.7",
        "orig_seurat_clusters",
        "clusters_annotated",
        "old.ident",
        "integrated_snn_res.0.7",
        "clust.names",
    ]
)


def import_seurat_column_FL(data, colname, save_key, fill_val):
    col = data[colname]
    col = col.reindex(ds_FL.cells.fetch("ids")).values
    ds_FL.cells.insert(save_key, col, overwrite=True, fill_value=fill_val)
    return None


import_seurat_column_FL(metadata, "WNN_u1", "WNN_UMAP1", -1)
import_seurat_column_FL(metadata, "WNN_u2", "WNN_UMAP2", -1)
import_seurat_column_FL(metadata, "wsnn_res.0.5", "WNN_clusts", -1)

ds_FL.plot_layout(
    layout_key="WNN_UMAP",
    color_by="WNN_clusts",
    height=7,
    width=7,
    legend_onside=True,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    cmap="Paired",
    save_dpi=300,
    savename=out_files["umap_wnn_clusters"],
)

color_palette_FL = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"), ds_FL.cells.fetch("RNA_cluster_names")
)
for i in x:
    color_palette_FL[i] = x[i][x[i] > 0].index[0]

ds_FL.plot_layout(
    layout_key="WNN_UMAP",
    color_by="RNA_cluster_names",
    color_key=color_palette_FL,
    height=7,
    width=7,
    legend_onside=True,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    save_dpi=300,
    savename=out_files["umap_org_clusters"],
)

df = ds_FL.cells.to_pandas_dataframe(columns=["RNA_cluster_names", "WNN_clusts"])
df = pd.crosstab(df["RNA_cluster_names"], df["WNN_clusts"])
df = df / df.sum() * 100
df.drop(columns=-1)
df.drop(index="NA")

# + tags=[]
clust_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.sort_values().cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_FL[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of projected cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot"], dpi=300)
plt.show()
