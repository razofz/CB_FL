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

# + tags=[]
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

# + tags=[]
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
cs22_dir = "data/processed/seurat/CS22"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
figures_dir = out_dir + "/Figures"

# + tags=[]
in_files = {
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "cluster_markers": glob.glob(seurat_dir + "/Individual_cluster_markers/*.csv"),
    "new_cluster_markers": {
        "hsc": cs22_dir + "/cc_mpp1_v_hsc.csv",
        "mpp2": cs22_dir + "/cc_mpp1_v_mpp2.csv",
    },
}

out_files = {
    "plot": figures_dir + "/heatmap_prim_FL_clusts.svg",
    "genes_p5": out_dir + "/fig3f_tables_point5.csv",
    "genes_p25": out_dir + "/fig3f_tables_point25.csv",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(out_adt_matrices_dir)
# fl_utils.create_dir(out_adt_matrices_dir + "/FCS_files")
# -

# ---

# + tags=[]
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


# + tags=[]
ds = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
# ds = scarf.DataStore("./_zarr_files/RNA_FL_hpc_merged", default_assay="RNA")

# + tags=[]
feat_idx = np.where(ds.RNA.feats.fetch_all("nCells") >= 10)[0]
ban_genes = ds.RNA.feats.get_index_by(
    ds.RNA.feats.grep("^RPS|^RPL|^MRPS|^MRPL|^MT-"), "names"
)
feat_idx = np.array(list(set(feat_idx).difference(ban_genes)))
data = ds.RNA.normed(feat_idx=feat_idx, renormalize_subset=True).compute()
data.shape

# + tags=[]
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


# + tags=[]
def get_goi_values(genes):
    df = {}
    for i in genes:
        if i in data:
            df[i] = data[i]
    df = pd.DataFrame(df)
    return df


# + tags=[]
# genes_up_dn = {"point" + str(threshold): [] for threshold in [25, 5]}
# for fn in in_files["new_cluster_markers"].items():
# # for fn in in_files["cluster_markers"]:
#     # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
#     print(fn)
#     fn = fn[1]
#     res = pd.read_csv(fn)
#     res["avg_log2FC"] = res["avg_log2FC"].fillna(0)
#     res["p_val_adj"] = res["p_val_adj"].fillna(1)
#     genes = list(
#         res[(abs(res["avg_log2FC"]) > 0.25) & (res["p_val_adj"] < 1e-3)].index
#         # res[(abs(res["avg_log2FC"]) > 0.5) & (res["p_val_adj"] < 1e-3)].index
#     )
#     genes_up_dn["point25"] += genes
#     print(len(genes_up_dn), len(genes))
# genes_up_dn["point25"] = pd.unique(genes_up_dn["point25"])

# + tags=[]
# for fn in in_files["new_cluster_markers"].items():
# # for fn in in_files["cluster_markers"]:
#     # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
#     print(fn)
#     fn = fn[1]
#     res = pd.read_csv(fn)
#     res["avg_log2FC"] = res["avg_log2FC"].fillna(0)
#     res["p_val_adj"] = res["p_val_adj"].fillna(1)
#     genes = list(
#         res[(abs(res["avg_log2FC"]) > 0.25) & (res["p_val_adj"] < 1e-3)].index
#         # res[(abs(res["avg_log2FC"]) > 0.5) & (res["p_val_adj"] < 1e-3)].index
#     )
#     genes_up_dn["point5"] += genes
#     print(len(genes_up_dn), len(genes))
# genes_up_dn["point5"] = pd.unique(genes_up_dn["point5"])
# -

pd.read_csv(in_files["new_cluster_markers"]["hsc"], index_col=0)

pd.read_csv(in_files["new_cluster_markers"]["mpp2"], index_col=0)

threshold = 0.25
# threshold = 0.5

# + tags=[]
genes_up_dn = []
for fn in in_files["new_cluster_markers"].items():
# for fn in in_files["cluster_markers"]:
    # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
    print(fn)
    fn = fn[1]
    res = pd.read_csv(fn, index_col=0)
    res["avg_log2FC"] = res["avg_log2FC"].fillna(0)
    res["p_val_adj"] = res["p_val_adj"].fillna(1)
    genes = list(
        res[(abs(res["avg_log2FC"]) > threshold) & (res["p_val_adj"] < 1e-3)].index
        # res[(abs(res["avg_log2FC"]) > 0.25) & (res["p_val_adj"] < 1e-3)].index
        # res[(abs(res["avg_log2FC"]) > 0.5) & (res["p_val_adj"] < 1e-3)].index
    )
    genes_up_dn += genes
    print(len(genes_up_dn), len(genes))
genes_up_dn = pd.unique(genes_up_dn)

# + tags=[]
genes_up_dn
# -

pd.DataFrame(genes_up_dn).to_csv(out_files["genes_p25"])

threshold = 0.5

genes_up_dn = []
for fn in in_files["new_cluster_markers"].items():
# for fn in in_files["cluster_markers"]:
    # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
    print(fn)
    fn = fn[1]
    res = pd.read_csv(fn, index_col=0)
    res["avg_log2FC"] = res["avg_log2FC"].fillna(0)
    res["p_val_adj"] = res["p_val_adj"].fillna(1)
    genes = list(
        res[(abs(res["avg_log2FC"]) > threshold) & (res["p_val_adj"] < 1e-3)].index
        # res[(abs(res["avg_log2FC"]) > 0.25) & (res["p_val_adj"] < 1e-3)].index
        # res[(abs(res["avg_log2FC"]) > 0.5) & (res["p_val_adj"] < 1e-3)].index
    )
    genes_up_dn += genes
    print(len(genes_up_dn), len(genes))
genes_up_dn = pd.unique(genes_up_dn)

pd.DataFrame(genes_up_dn).to_csv(out_files["genes_p5"])

# + tags=[]
# genes_up_dn = []
# for fn in in_files["new_cluster_markers"].items():
# # for fn in in_files["cluster_markers"]:
#     # for fn in glob.glob("./from_seurat/Individual_cluster_markers/*.csv"):
#     print(fn)
#     fn = fn[1]
#     name = fn.split("/")[-1].split(".csv")[0]
#     if (
#         (name == "HSC_cluster_markers")
#         | (name == "MPP-I_cluster_markers")
#         | (name == "MPP-II_cluster_markers")
#     ):
#         print(name)
#         res2 = pd.read_csv(fn)
#         res2["avg_log2FC"] = res2["avg_log2FC"].fillna(0)
#         res2["p_val_adj"] = res2["p_val_adj"].fillna(1)
#         genes = list(
#             res2[(abs(res2["avg_log2FC"]) > 0.5) & (res2["p_val_adj"] < 1e-3)].index
#         )
#         genes_up_dn += genes
#         print(len(genes_up_dn), len(genes))
# genes_up_dn = pd.unique(genes_up_dn)

# + tags=[]
# [len(pd.unique(genes)) for genes in list(genes_up_dn.values())]
# -

len(list(genes_up_dn))

# + tags=[]
df_goi = get_goi_values(genes_up_dn)

# + tags=[]
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

# + tags=[]
len(df_goi_prim.columns[(df_goi_prim.sum() > 0)])

# + tags=[]
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
# plt.savefig(out_files["plot"], dpi=300)
# plt.savefig(f"tmp/heatmap_fig3f_threshold_{threshold}.svg", dpi=300)
plt.savefig(out_files["plot"], dpi=300)
plt.show()
# -

# ---

# +
# cmg = sns.clustermap(
#     df_goi_prim.T,
#     col_colors=colors,
#     robust=True,
#     method="ward",
#     xticklabels=[],
#     yticklabels=df_goi.columns,
#     z_score=0,
#     cmap="bwr",
#     figsize=(20, 20),
#     rasterized=True,
# )
# clean_axis(cmg.ax_heatmap)
# plt.show()

# + tags=[]
# den = scipy.cluster.hierarchy.dendrogram(
#     cmg.dendrogram_row.linkage, labels=df_goi_prim.T.index, color_threshold=95
# )

# +
# den = scipy.cluster.hierarchy.dendrogram(
#     cmg.dendrogram_row.linkage, labels=df_goi_prim.T.index, color_threshold=95
# )

# +
# def get_cluster_classes(den, label="ivl"):
#     cluster_idxs = defaultdict(list)
#     for c, pi in zip(den["color_list"], den["icoord"]):
#         for leg in pi[1:3]:
#             i = (leg - 5.0) / 10.0
#             if abs(i - int(i)) < 1e-5:
#                 cluster_idxs[c].append(int(i))

#     cluster_classes = {}
#     for c, l in cluster_idxs.items():
#         i_l = [den[label][i] for i in l]
#         cluster_classes[c] = i_l

#     return cluster_classes

# +
# clusters = get_cluster_classes(den)

# + tags=[]
# for x in clusters.keys():
#     print(x)
#     for y in clusters[x]:
#         print(y)

# +
# cluster = []
# for i in df_goi_prim.T.index:
#     included = False
#     for j in clusters.keys():
#         if i in clusters[j]:
#             cluster.append(j)
#             included = True
#     if not included:
#         cluster.append(None)

# len(cluster)

# +
# gene_cluster_cols = {
#     "C1": "#e74c3c",
#     "C2": "#9b59b6",
#     "C3": "#3498db",
#     "C4": "#1abc9c",
#     "C5": "#27ae60",
#     "C6": "#f1c40f",
#     "C7": "#e67e22",
#     "C8": "#839192",
#     "C9": "#5d6d7e",
# }

# +
# row_colors = [gene_cluster_cols[x] for x in cluster]

# +
# cmg = sns.clustermap(
#     df_goi_prim.T,
#     col_colors=colors,
#     row_colors=row_colors,
#     robust=True,
#     method="ward",
#     xticklabels=[],
#     yticklabels=df_goi.columns,
#     z_score=0,
#     cmap="bwr",
#     figsize=(20, 20),
#     rasterized=True,
# )
# clean_axis(cmg.ax_heatmap)
# plt.savefig(out_files["plot"], dpi=300)
# # plt.savefig("./images/heatmap_prim_FL_clusts.svg", dpi=300)
# plt.show()
