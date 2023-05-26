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
#     display_name: notebooks_FL
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
import seaborn as sns

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
    "zarr": zarr_dir + "/RNA_yBM_and_FL_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_BM": zarr_dir + "/RNA_yBM_hpc_merged",
    "metadata": raw_dir + "/FL_BM_combined_metadata.csv",
    "umap": raw_dir + "/FL_BM_combined_umap.csv",
}

out_files = {
    "umap_harmony_clusters": figures_dir + "/harmony_UMAP_harmony_clusters.svg",
    "umap_harmony_clusters_BM_only": figures_dir
    + "/harmony_UMAP_harmony_clusters_BM_only.svg",
    "umap_harmony_clusters_FL_only": figures_dir
    + "/harmony_UMAP_harmony_clusters_FL_only.svg",
    "umap_FL_harmony_clusters_FL_only": figures_dir
    + "/harmony_UMAP_FL_clusters_FL_only.svg",
    "harmony_UMAP_predicted_BM_clusters_BM_only": figures_dir
    + "/harmony_UMAP_predicted_BM_clusters_BM_only.svg",
    "harmony_UMAP_BM_cell_clusters": figures_dir + "/harmony_UMAP_BM_cell_clusters.svg",
    "umap_harmony_sample_distribution": figures_dir
    + "/harmony_UMAP_sample_distribution.svg",
    "barplot_org_BM": figures_dir + "/bar_plot_harmony_org_BM_clust_comparison.svg",
    "barplot_FL": figures_dir + "/bar_plot_harmony_FL_clust_comparison.svg",
    "barplot_predicted_BM": figures_dir
    + "/bar_plot_harmony_predicted_BM_clust_comparison.svg",
    "barplot_sample_contrib": figures_dir + "/bar_plot_harmony_sample_contribution.svg",
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

ds = scf.DataStore(in_files["zarr"], default_assay="RNA")

metadata["harmony_u1"] = umap["UMAPharmonyRNA_1"]
metadata["harmony_u2"] = umap["UMAPharmonyRNA_2"] * (
    -1
)  # multipled by -1 to flip umap so HSCs are on top
metadata.index = [x.split(" ")[0] for x in metadata.index]
metadata.index = [
    x.split("_")[1].lower() + "__" + x.split("_")[-1] for x in metadata.index
]


def import_seurat_column_FL(data, colname, save_key, fill_val):
    col = data[colname]
    col = col.reindex(ds.cells.fetch("ids")).values
    ds.cells.insert(save_key, col, overwrite=True, fill_value=fill_val)
    return None


import_seurat_column_FL(metadata, "harmony_u1", "harmony_UMAP1", -1)
import_seurat_column_FL(metadata, "harmony_u2", "harmony_UMAP2", -1)
import_seurat_column_FL(metadata, "clusters_RNA_harmony", "harmony_clusts", -1)

import_seurat_column_FL(metadata, "paper_clusters", "bm_clusts", -1)

ds.cells.insert(
    "BM_cells",
    [x[:-1] == "ybm" for x in ds.cells.fetch_all("sample")],
    overwrite=True,
    fill_value="NA",
)
ds.cells.insert(
    "FL_cells",
    [x[:-1] == "fl" for x in ds.cells.fetch_all("sample")],
    overwrite=True,
    fill_value="NA",
)

colors = {
    0: "#a6cee3",
    1: "#1f78b4",
    2: "#b2df8a",
    3: "#33a02c",
    4: "#fb9a99",
    5: "#e31a1c",
    6: "#fdbf6f",
    7: "#ff7f00",
    8: "#cab2d6",
    9: "#6a3d9a",
    10: "#ffff99",
    11: "#b15928",
    12: "#a6cee3",
    -1: "#000000",
}

# + tags=[]
ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="harmony_clusts",
    height=7,
    width=7,
    legend_onside=True,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    color_key=colors,
    save_dpi=300,
    savename=out_files["umap_harmony_clusters"],
)

# + tags=[]
ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="harmony_clusts",
    subselection_key="BM_cells",
    height=7,
    width=7,
    legend_onside=True,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    color_key=colors,
    save_dpi=300,
    savename=out_files["umap_harmony_clusters_BM_only"],
)
# -

ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="harmony_clusts",
    subselection_key="FL_cells",
    height=7,
    width=7,
    legend_onside=True,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    color_key=colors,
    save_dpi=300,
    savename=out_files["umap_harmony_clusters_FL_only"],
)

# +
ds_FL = scf.DataStore(in_files["zarr_FL"], default_assay="RNA", nthreads=threads)
color_palette_FL = {}
color_palette_pred_BM = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"), ds_FL.cells.fetch("RNA_cluster_names")
)
for i in x:
    color_palette_FL["fl_" + i] = x[i][x[i] > 0].index[0]
    color_palette_pred_BM["bm_" + i] = x[i][x[i] > 0].index[0]

color_palette_pred_BM["bm_NA"] = "#D3D3D3"
# -

ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="clusters",
    subselection_key="FL_cells",
    height=7,
    width=7,
    legend_onside=True,
    color_key=color_palette_FL,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    save_dpi=300,
    savename=out_files["umap_FL_harmony_clusters_FL_only"],
)

ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="clusters",
    subselection_key="BM_cells",
    height=7,
    width=7,
    legend_onside=True,
    color_key=color_palette_pred_BM,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    save_dpi=300,
    savename=out_files["harmony_UMAP_predicted_BM_clusters_BM_only"],
)

# +
ds_BM = scf.DataStore(in_files["zarr_BM"], nthreads=threads)

df_BM = ds_BM.cells.to_pandas_dataframe(["ids", "RNA_cluster_names"])
df_BM.index = df_BM["ids"]
df_BM["RNA_cluster_names"] = ["bm_" + x for x in df_BM["RNA_cluster_names"]]
df_BM.columns = ["ids", "clusters"]

cell_clusters = pd.DataFrame()
cell_clusters = df_BM

cell_clusters = cell_clusters["clusters"].reindex(ds.cells.fetch("ids")).values
ds.cells.insert("bm_clusters", cell_clusters, overwrite=True, fill_value="NA")

# +
# ds_FL = scf.DataStore(zarr_path + "RNA_FL_hpc_merged", nthreads=2)

df_FL = ds_FL.cells.to_pandas_dataframe(["ids", "RNA_cluster_names"])
df_FL.index = df_FL["ids"]
df_FL["RNA_cluster_names"] = ["fl_" + x for x in df_FL["RNA_cluster_names"]]
df_FL.columns = ["ids", "clusters"]

cell_clusters = pd.DataFrame()
cell_clusters = df_FL

cell_clusters = cell_clusters["clusters"].reindex(ds.cells.fetch("ids")).values
cell_clusters = np.array(cell_clusters, dtype="str")
ds.cells.insert("fl_clusters", cell_clusters, overwrite=True, fill_value="NA")
# -

color_palette_BM = {}
x = pd.crosstab(
    ds_BM.cells.fetch("RNA_cluster_colors"), ds_BM.cells.fetch("RNA_cluster_names")
)
for i in x:
    color_palette_BM["bm_" + i] = x[i][x[i] > 0].index[0]

ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="bm_clusters",
    subselection_key="BM_cells",
    height=7,
    width=7,
    legend_onside=True,
    color_key=color_palette_BM,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    save_dpi=300,
    savename=out_files["harmony_UMAP_BM_cell_clusters"],
)

df = ds.cells.to_pandas_dataframe(columns=["bm_clusters", "harmony_clusts"])
df = pd.crosstab(df["bm_clusters"], df["harmony_clusts"])
df = df.drop(columns=-1)
df = df.drop(index="NA")
df = df.drop(index="nan")
df = df / df.sum() * 100

color_palette_BM_slight_change = color_palette_BM
color_palette_BM_slight_change["MEP-II"] = "#cc4a4a"
color_palette_BM_slight_change["MPP-II"] = "#ffc4d4"

# +
clust_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.sort_values().cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_BM_slight_change[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of BM cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_org_BM"], dpi=300)
plt.show()
# -

df = ds.cells.to_pandas_dataframe(columns=["clusters", "harmony_clusts"])
df = pd.crosstab(df["clusters"], df["harmony_clusts"])
df = df.drop(columns=-1)
df = df.drop(index="NA")
df = df.drop(
    [
        "bm_Cyc",
        "bm_DC-I",
        "bm_DC-Mono",
        "bm_GMP",
        "bm_HSC",
        "bm_Ly-I",
        "bm_Ly-II",
        "bm_MEP",
        "bm_MPP-I",
        "bm_MPP-II",
        "bm_NA",
        "bm_Ly-III",
    ]
)
df = df / df.sum() * 100

# +
clust_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.sort_values().cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_FL[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of FL cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_FL"], dpi=300)
plt.show()
# -

df = ds.cells.to_pandas_dataframe(columns=["clusters", "harmony_clusts"])
df = pd.crosstab(df["clusters"], df["harmony_clusts"])
df = df.drop(columns=-1)
df = df.drop(index="NA")
df = df.drop(
    [
        "fl_Cyc",
        "fl_DC-I",
        "fl_DC-Mono",
        "fl_GMP",
        "fl_HSC",
        "fl_Ly-I",
        "fl_Ly-II",
        "fl_MEP",
        "fl_MPP-I",
        "fl_MPP-II",
        "fl_Ly-III",
    ]
)
df = df / df.sum() * 100

# +
clust_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.sort_values().cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_pred_BM[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of predicted BM cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_predicted_BM"], dpi=300)
plt.show()
# -

df = ds.cells.to_pandas_dataframe(columns=["sample", "harmony_clusts"])
df = pd.crosstab(df["sample"], df["harmony_clusts"])
df = df.drop(columns=-1)
df = df / df.sum() * 100

pal = sns.color_palette("tab20", 4)
print(pal.as_hex())

# +
clust_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

sample_colors = {
    "fl1": "#1f77b4",
    "fl2": "#aec7e8",
    "ybm1": "#ff7f0e",
    "ybm2": "#ffbb78",
}

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.sort_values().cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = sample_colors[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of sample in each cluster", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_sample_contrib"], dpi=300)
plt.show()
# -

ds.plot_layout(
    layout_key="harmony_UMAP",
    color_by="sample",
    height=7,
    width=7,
    legend_onside=True,
    color_key=sample_colors,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    shuffle_df=True,
    save_dpi=300,
    savename=out_files["umap_harmony_sample_distribution"],
)
