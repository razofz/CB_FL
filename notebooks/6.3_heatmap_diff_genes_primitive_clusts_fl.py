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
from collections import defaultdict

import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scarf
import scipy
import seaborn as sns

# + tags=[]
fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "cluster_markers": glob.glob(
        seurat_dir + "/Individual_cluster_markers/*.csv"
    ),
}

out_files = {
    "plot": figures_dir + "/heatmap_prim_FL_clusts.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(out_adt_matrices_dir)
# fl_utils.create_dir(out_adt_matrices_dir + "/FCS_files")
# -

# ---

# +
plt.style.use("fivethirtyeight")
plt.rcParams["svg.fonttype"] = "none"


def clean_axis(ax, ts=11, ga=0.6):
    ax.xaxis.set_tick_params(labelsize=ts)
    ax.yaxis.set_tick_params(labelsize=ts)
    for i in ["top", "bottom", "left", "right"]:
        ax.spines[i].set_visible(False)
    ax.grid(which="major", linestyle="--", alpha=ga)
    ax.figure.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    return True


# -

ds = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
# ds = scarf.DataStore("./_zarr_files/RNA_FL_hpc_merged", default_assay="RNA")

feat_idx = np.where(ds.RNA.feats.fetch_all("nCells") >= 10)[0]
ban_genes = ds.RNA.feats.get_index_by(
    ds.RNA.feats.grep("^RPS|^RPL|^MRPS|^MRPL|^MT-"), "names"
)
feat_idx = np.array(list(set(feat_idx).difference(ban_genes)))
data = ds.RNA.normed(feat_idx=feat_idx, renormalize_subset=True).compute()
data.shape

data = pd.DataFrame(data, columns=ds.RNA.feats.fetch_all("names")[feat_idx])
data.index = ds.cells.fetch("ids")
data.shape


# +
# clust_order = [
#     "T",
#     "DC-I",
#     "Cyc",
#     "MEP",
#     "DC-Mono",
#     "GMP",
#     "Ly-II",
#     "Ly-I",
#     "MPP-II",
#     "MPP-I",
#     "HSC",
# ][::-1]
# clust_colors = [
#     ds.cells.fetch("RNA_cluster_colors")[ds.cells.fetch("RNA_cluster_names") == x][0]
#     for x in clust_order
# ]
# -


def get_goi_values(genes):
    df = {}
    for i in genes:
        if i in data:
            df[i] = data[i]
    df = pd.DataFrame(df)
    return df


genes_up_dn = []
for fn in in_files["cluster_markers"]:
    # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
    name = fn.split("/")[-1].split(".csv")[0]
    if (
        (name == "HSC_cluster_markers")
        | (name == "MPP-I_cluster_markers")
        | (name == "MPP-II_cluster_markers")
    ):
        print(name)
        res2 = pd.read_csv(fn)
        res2["avg_log2FC"] = res2["avg_log2FC"].fillna(0)
        res2["p_val_adj"] = res2["p_val_adj"].fillna(1)
        genes = list(
            res2[
                (abs(res2["avg_log2FC"]) > 0.5) & (res2["p_val_adj"] < 1e-3)
            ].index
        )
        genes_up_dn += genes
        print(len(genes_up_dn), len(genes))
genes_up_dn = pd.unique(genes_up_dn)

print(len(pd.unique(genes_up_dn)))

df_goi = get_goi_values(genes_up_dn)

df_goi_prim = df_goi[
    (ds.cells.fetch("RNA_cluster_names") == "HSC")
    | (ds.cells.fetch("RNA_cluster_names") == "MPP-I")
    | (ds.cells.fetch("RNA_cluster_names") == "MPP-II")
]
colors = col_colors = ds.cells.fetch("RNA_cluster_colors")[
    (ds.cells.fetch("RNA_cluster_names") == "HSC")
    | (ds.cells.fetch("RNA_cluster_names") == "MPP-I")
    | (ds.cells.fetch("RNA_cluster_names") == "MPP-II")
]

len(df_goi_prim.columns[(df_goi_prim.sum() > 0)])

cmg = sns.clustermap(
    df_goi_prim.T,
    col_colors=colors,
    robust=True,
    method="ward",
    xticklabels=[],
    yticklabels=df_goi.columns,
    z_score=0,
    cmap="bwr",
    figsize=(20, 20),
    rasterized=True,
)
clean_axis(cmg.ax_heatmap)
plt.show()

den = scipy.cluster.hierarchy.dendrogram(
    cmg.dendrogram_row.linkage, labels=df_goi_prim.T.index, color_threshold=95
)


def get_cluster_classes(den, label="ivl"):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den["color_list"], den["icoord"]):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes


clusters = get_cluster_classes(den)

# + tags=[]
for x in clusters.keys():
    print(x)
    for y in clusters[x]:
        print(y)

# +
cluster = []
for i in df_goi_prim.T.index:
    included = False
    for j in clusters.keys():
        if i in clusters[j]:
            cluster.append(j)
            included = True
    if not included:
        cluster.append(None)

len(cluster)
# -

gene_cluster_cols = {
    "C1": "#e74c3c",
    "C2": "#9b59b6",
    "C3": "#3498db",
    "C4": "#1abc9c",
    "C5": "#27ae60",
    "C6": "#f1c40f",
    "C7": "#e67e22",
    "C8": "#839192",
    "C9": "#5d6d7e",
}

row_colors = [gene_cluster_cols[x] for x in cluster]

cmg = sns.clustermap(
    df_goi_prim.T,
    col_colors=colors,
    row_colors=row_colors,
    robust=True,
    method="ward",
    xticklabels=[],
    yticklabels=df_goi.columns,
    z_score=0,
    cmap="bwr",
    figsize=(20, 20),
    rasterized=True,
)
clean_axis(cmg.ax_heatmap)
plt.savefig(out_files["plot"], dpi=300)
# plt.savefig("./images/heatmap_prim_FL_clusts.svg", dpi=300)
plt.show()
