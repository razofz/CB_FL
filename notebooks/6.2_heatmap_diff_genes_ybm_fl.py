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

import glob

import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scarf
import seaborn as sns


# + tags=[]
fl_utils.set_project_path()

# + tags=[]
threads = 16

raw_dir = "data/raw"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
sup_tables_raw_dir = raw_dir + "/paper_specific/Supplemental_tables"
figures_dir = out_dir + "/Figures"

# + tags=[]
in_files = {
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_merged": zarr_dir + "/RNA_yBM_and_FL_hpc_merged",
    "fetal_signature": sup_tables_raw_dir + "/Fetal_signature_p005.txt",
}

out_files = {
    "plot_FL": figures_dir + "/heatmap_yBM_FL_fetal_core_no_BM_Gene.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# + tags=[]
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
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
ds = scarf.DataStore(in_files["zarr_merged"], default_assay="RNA")

# + tags=[]
len(ds.cells.fetch("ids"))

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


# + tags=[]
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
#     ds_FL.cells.fetch("RNA_cluster_colors")[
#         ds_FL.cells.fetch("RNA_cluster_names") == x
#     ][0]
#     if len(ds_FL.cells.fetch("RNA_cluster_colors")[
#         ds_FL.cells.fetch("RNA_cluster_names") == x
#     ]) > 0
#     else []
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
venn_res = pd.read_csv(in_files["fetal_signature"], delimiter="\t", header=None)
genes = [x for x in venn_res[0].values]

# + tags=[]
df_goi = get_goi_values(genes)

# + tags=[]
cmg = sns.clustermap(
    df_goi.T,
    col_colors=ds.cells.fetch("sample_colors"),
    # row_colors=row_colors,
    robust=True,
    method="ward",
    row_cluster=True,
    xticklabels=[],
    yticklabels=df_goi.columns,
    z_score=0,
    cmap="bwr",
    figsize=(20, 20),
    rasterized=True,
)
clean_axis(cmg.ax_heatmap)
plt.savefig(out_files["plot_FL"], dpi=300)
# plt.savefig("./images/heatmap_yBM_FL_fetal_core_no_BM_Gene.svg", dpi=300)
plt.show()
