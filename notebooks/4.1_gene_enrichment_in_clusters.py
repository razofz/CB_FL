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

import scarf
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
import fl_utils

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


fl_utils.set_project_path()

raw_dir = "data/raw/paper_specific/from_seurat"
zarr_dir = "data/processed/notebooks/zarr_files"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "metadata_yBM": raw_dir + "/yBM_hpc_cite_seq_celldata_seurat.csv",
    "metadata_FL": "data/processed/notebooks/FL_seurat_metadata.csv",
    "cluster_markers_FL": raw_dir + "/FL_combined_cluster_markers.csv",
}

out_files = {
    "dot_plot": figures_dir + "/dot_plot_diff_genes_FL.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

ds_yBM = scarf.DataStore(in_files["zarr_yBM"])
ds_FL = scarf.DataStore(in_files["zarr_FL"])

for i in [
    "Ly-III",
    "DC-I",
    "Cyc",
    "MEP",
    "DC-Mono",
    "GMP",
    "Ly-II",
    "Ly-I",
    "MPP-II",
    "MPP-I",
    "HSC",
]:
    print(
        i,
        len(
            ds_FL.cells.fetch("RNA_cluster_names")[
                ds_FL.cells.fetch("RNA_cluster_names") == i
            ]
        ),
    )


metadata_FL = pd.read_csv(in_files["metadata_FL"], index_col=0)
metadata_BM = pd.read_csv(in_files["metadata_yBM"], index_col=0)
cluster_markers_FL = pd.read_csv(in_files["cluster_markers_FL"])

cluster_markers_FL.cluster.unique()

clust_names_to_num = {
    "Ly-III": 10,
    "DC-I": 9,
    "Cyc": 8,
    "MEP": 3,
    "DC-Mono": 7,
    "GMP": 4,
    "Ly-II": 6,
    "Ly-I": 2,
    "MPP-II": 1,
    "MPP-I": 0,
    "HSC": 5,
}

# +
fig, ax = plt.subplots(1, 1, figsize=(17, 4))
clust_order = [
    "Ly-III",
    "DC-I",
    "Cyc",
    "MEP",
    "DC-Mono",
    "GMP",
    "Ly-II",
    "Ly-I",
    "MPP-II",
    "MPP-I",
    "HSC",
]  # [::-1]
pos = 0
names = []
for i in clust_order[::-1]:
    genes = cluster_markers_FL[
        cluster_markers_FL.cluster == clust_names_to_num[i]
    ]["gene"][:10].values
    print(i, genes)

    feat_idx = ds_FL.RNA.feats.get_index_by(genes, column="names")
    df = pd.DataFrame(
        ds_FL.RNA.normed(feat_idx=feat_idx).compute(),
        columns=ds_FL.RNA.feats.fetch_all("names")[feat_idx],
    )
    df.index = [
        x.split("__")[1] + "-" + x[2] + "_FL_hpc"
        for x in ds_FL.cells.fetch("ids")
    ]
    df = df.reindex(metadata_FL.index)

    mdf = df.apply(lambda x: (x - x.mean()) / x.std())
    mdf["cluster"] = metadata_FL.cluster_names.copy()
    mdf = mdf.groupby("cluster").mean()
    mdf = mdf.reindex(clust_order)

    bdf = df.copy()
    bdf[bdf > 0] = 1
    bdf["cluster"] = metadata_FL.cluster_names.copy()
    bdf = bdf.groupby("cluster").mean()
    bdf = bdf.reindex(clust_order)

    n = 11
    for j in mdf:
        if j in names:
            continue
        if (bdf[j] > 0.9).sum() < 5:
            ax.scatter(
                np.repeat(pos, n),
                range(n),
                s=120 * bdf[j],
                c=mdf[j],
                vmin=-2,
                vmax=2,
                cmap="bwr",
                alpha=0.7,
                lw=0.3,
                edgecolors="k",
                rasterized=True,
            )
            names.append(j)
            pos += 1

clean_axis(ax, ga=0.1)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=90, fontsize=10)
ax.set_yticks(range(11))
ax.set_yticklabels(clust_order, fontsize=12)
ax.set_xlim((-1, len(names)))
ax.set_ylim((-1, 11))
plt.tight_layout()
plt.savefig(out_files["dot_plot"], dpi=300)
plt.show()


# +
fig, ax = plt.subplots(1, 1, figsize=(15, 4))
clust_order = [
    "Ly-III",
    "DC-I",
    "Cyc",
    "MEP",
    "DC-Mono",
    "GMP",
    "Ly-II",
    "Ly-I",
    "MPP-II",
    "MPP-I",
    "HSC",
]
pos = 0
names = []

genes = ds_FL.ADT.feats.fetch_all("names")

feat_idx = ds_FL.ADT.feats.get_index_by(genes, column="names")
df = pd.DataFrame(
    ds_FL.ADT.normed(feat_idx=feat_idx).compute(),
    columns=ds_FL.ADT.feats.fetch_all("names")[feat_idx],
)
df.index = [
    x.split("__")[1] + "-" + x[2] + "_FL_hpc" for x in ds_FL.cells.fetch("ids")
]
df = df.reindex(metadata_FL.index)

mdf = df.copy()
mdf["cluster"] = metadata_FL.cluster_names.copy()
mdf = mdf.groupby("cluster").mean()
mdf = mdf.reindex(clust_order)
mdf = mdf.T.sort_values(by=["HSC"], ascending=False).T

adt_gene_names = {
    "BAH1": "MPL",
    "CD10": "MME",
    "CD105": "ENG",
    "CD107A": "LAMP1",
    "CD117": "KIT",
    "CD11A": "ITGAL",
    "CD11C": "ITGAX",
    "CD123": "IL3RA",
    "CD135": "FLT3",
    "CD155": "PVR",
    "CD18": "ITGB2",
    "CD196": "CCR6",
    "CD201": "PROCR",
    "CD25": "IL2RA",
    "CD26": "DPP4",
    "CD32": "FCGR2A",
    "CD33": "CD33",
    "CD34": "CD34",
    "CD35": "CR1",
    "CD36": "CD36",
    "CD38": "CD38",
    "CD4": "CD4",
    "CD41": "ITGA2B",
    "CD42B": "GP1BA",
    "CD44": "CD44",
    "CD45RA": "PTPRC",
    "CD45RO": "PTPRC",
    "CD48": "CD48",
    "CD49D": "ITGA4",
    "CD49F": "ITGA6",
    "CD52": "CD52",
    "CD54": "ICAM1",
    "CD56": "NCAM1",
    "CD61": "ITGB3",
    "CD62L": "SELL",
    "CD7": "CD7",
    "CD70": "CD70",
    "CD71": "TFRC",
    "CD79B": "CD79B",
    "CD9": "CD9",
    "CD90": "THY1",
    "CD93": "CD93",
    "CD95": "FAS",
    "IL7RA": "IL7R",
}

bdf = df.copy()
bdf["cluster"] = metadata_FL.cluster_names.copy()
bdf = bdf.groupby("cluster").mean()
bdf = bdf.reindex(clust_order)
bdf = bdf[list(mdf.columns)]

print(mdf.shape, bdf.shape)

n = 11
for j1, j2 in zip(mdf, bdf):
    if len(bdf[j2].shape) == 2:
        c = bdf[j2].mean(axis=1).values
    else:
        c = bdf[j2].values
    c = (c - c.mean()) / c.std()
    ax.scatter(
        np.repeat(pos, n),
        range(n),
        s=100 * mdf[j1],
        c=c,
        vmin=-2,
        vmax=2,
        cmap="bwr",
        alpha=0.7,
        lw=0.3,
        edgecolors="k",
        rasterized=True,
    )
    names.append(j1)
    pos += 1
clean_axis(ax, ga=0.1)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(
    [f"{x[4:]} ({adt_gene_names[x[4:]]})" for x in names],
    rotation=45,
    fontsize=10,
    ha="right",
)
ax.set_yticks(range(11))
ax.set_yticklabels(clust_order, fontsize=12)
ax.set_xlim((-1, len(names)))
ax.set_ylim((-1, 11))
plt.tight_layout()
# plt.savefig('./images/dot_plot_ADTs_FL.svg', dpi=300)
plt.show()
