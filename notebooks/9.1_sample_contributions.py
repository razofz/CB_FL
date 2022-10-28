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
import json
import os
from collections import Counter

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
tables_dir = "data/processed/notebooks/Supplemental_tables"

# + tags=[]
in_files = {
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_BM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarrs": {
        x: f"{zarr_dir}/RNA_{x}_merged"
        for x in ["FL_CS16", "FL_hpc", "FL_W9", "CB_hpc", "yBM_hpc"]
    },
}

out_files = {
    "umap_sample_contrib": figures_dir + "/FL_UMAP_sample_contribution.svg",
    "bar_plot_sample_contrib": figures_dir
    + "/bar_plot_FL_sample_contribution.svg",
    "BM_preds": tables_dir + "/BM_preds_sample_contribution.csv",
    "FL_preds": tables_dir + "/FL_preds_sample_contribution.csv",
    "FL_gating": tables_dir + "/sample_contribution_FL_gating.csv",
    "BM_gating": tables_dir + "/sample_contribution_BM_gating.csv",
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

ds_BM = scf.DataStore(in_files["zarr_BM"], default_assay="RNA")
ds_FL = scf.DataStore(in_files["zarr_FL"], default_assay="RNA")

ds_FL.cells.insert(
    "sample",
    [x.split("_")[0] for x in ds_FL.cells.fetch_all("ids")],
    overwrite=True,
)
ds_BM.cells.insert(
    "sample",
    [x.split("_")[0] for x in ds_BM.cells.fetch_all("ids")],
    overwrite=True,
)

sample_colors = {
    "fl1": "#1f77b4",
    "fl2": "#aec7e8",
    "ybm1": "#ff7f0e",
    "ybm2": "#ffbb78",
}

# +
color_palette_BM = {}
x = pd.crosstab(
    ds_BM.cells.fetch("RNA_cluster_colors"),
    ds_BM.cells.fetch("RNA_cluster_names"),
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

ds_FL.plot_layout(
    layout_key="RNA_seurat_UMAP",
    color_by="sample",
    height=7,
    width=7,
    legend_onside=True,
    color_key=sample_colors,
    scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.9},
    point_size=15,
    save_dpi=300,
    savename=out_files["umap_sample_contrib"],
)

df = ds_FL.cells.to_pandas_dataframe(columns=["sample", "RNA_cluster_names"])
df = pd.crosstab(df["sample"], df["RNA_cluster_names"])
df = df.drop(columns="NA")
df = df / df.sum() * 100
df

# +
clust_order = [
    "HSC",
    "MPP-I",
    "MPP-II",
    "MEP",
    "GMP",
    "DC-Mono",
    "Ly-I",
    "Ly-II",
    "Cyc",
    "DC-I",
    "Ly-III",
]
fig, ax = plt.subplots(1, 1, figsize=(10, 4))

for x in clust_order:
    df_sub = df[x]
    df_sub = df_sub.cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = sample_colors[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order])
ax.set_ylabel("% of sample in each cluster", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["bar_plot_sample_contrib"], dpi=300)
plt.show()
# -

classified_cells_BM = {}
for sample in sorted(list(ds_BM.z.RNA.projections.keys())):
    classified_cells_BM[sample] = ds_BM.get_target_classes(
        target_name=sample,
        threshold_fraction=0.4,
        reference_class_group="RNA_cluster_names",
    )
    print(sample, len(classified_cells_BM[sample]))

sample_oi = ["FL_CS16", "FL_hpc", "FL_W9", "CB_hpc"]
# data_stores = {
#     "FL_CS16": "RNA_FL_CS16_merged",
#     "FL_HPC": "RNA_FL_hpc_merged",
#     "FL_W9": "RNA_FL_W9_merged",
#     "CB_HPC": "RNA_CB_hpc_merged",
# }
preds_BM = {}
for i in sample_oi:
    # print(i, data_stores[i])
    print(i, in_files["zarrs"][i])
    # data_stores[i]
    ds2 = scf.DataStore(in_files["zarrs"][i])
    if "hpc" in i:
        i = i.upper()
    preds_BM[i] = dict(
        zip(ds2.cells.fetch("ids"), classified_cells_BM[i].values)
    )

# + tags=[]
sample_dist = {}
for i in preds_BM.keys():
    print(i)
    df = pd.DataFrame.from_dict(
        preds_BM[i], orient="index", columns=["cluster"]
    )
    df["sample"] = [i + "_" + x.split("__")[0] for x in df.index]
    df = pd.crosstab(df["sample"], df["cluster"])
    df = df.T / df.T.sum() * 100
    for y in df.columns:
        sample_dist[y] = dict(df[y])
# -

sample_dist = pd.DataFrame.from_dict(sample_dist, orient="index")
sample_dist = sample_dist.fillna(0)
sample_dist = sample_dist[
    [
        "HSC-I",
        "HSC-II",
        "HSC-III",
        "MPP-I",
        "MPP-II",
        "MPP-III",
        "MEP-I",
        "MEP-II",
        "GMP",
        "MB",
        "DC-I",
        "DC-II",
        "Ly-I",
        "Ly-II",
        "Cyc-I",
        "Cyc-II",
        "NA",
    ]
]
sample_dist.to_csv(out_files["BM_preds"])

classified_cells_FL = {}
for sample in sorted(list(ds_FL.z.RNA.projections.keys())):
    classified_cells_FL[sample] = ds_FL.get_target_classes(
        target_name=sample,
        threshold_fraction=0.4,
        reference_class_group="RNA_cluster_names",
    )
    print(sample, len(classified_cells_FL[sample]))

sample_oi = ["FL_CS16", "FL_W9", "CB_hpc", "yBM_hpc"]
# data_stores = {
#     "FL_CS16": "RNA_FL_CS16_merged",
#     "yBM_HPC": "RNA_yBM_hpc_merged",
#     "FL_W9": "RNA_FL_W9_merged",
#     "CB_HPC": "RNA_CB_hpc_merged",
# }
preds_FL = {}
for i in sample_oi:
    print(i, in_files["zarrs"][i])
    # data_stores[i]
    ds2 = scf.DataStore(in_files["zarrs"][i])
    if "CB_hpc" in i:
        i = i.upper()
    elif "yBM_hpc" in i:
        i = i[:3] + i[3:].upper()
    preds_FL[i] = dict(
        zip(ds2.cells.fetch("ids"), classified_cells_FL[i].values)
    )

sample_dist = {}
for i in preds_FL.keys():
    print(i)
    df = pd.DataFrame.from_dict(
        preds_FL[i], orient="index", columns=["cluster"]
    )
    df["sample"] = [i + "_" + x.split("__")[0] for x in df.index]
    df = pd.crosstab(df["sample"], df["cluster"])
    df = df.T / df.T.sum() * 100
    for y in df.columns:
        sample_dist[y] = dict(df[y])

sample_dist = pd.DataFrame.from_dict(sample_dist, orient="index")
sample_dist = sample_dist.fillna(0)
sample_dist = sample_dist[
    [
        "HSC",
        "MPP-I",
        "MPP-II",
        "MEP",
        "GMP",
        "DC-Mono",
        "Ly-I",
        "Ly-II",
        "Cyc",
        "DC-I",
        "Ly-III",
        "NA",
    ]
]
sample_dist.to_csv(out_files["FL_preds"])
sample_dist

# +
gates_oi = [
    "gated_CD10_CLP",
    "gated_CD38neg_MPP",
    "gated_CD49Fpos_HSC",
    "gated_CD7_CLP",
    "gated_CD90pos_HSC",
    "gated_CMP_123",
    "gated_CMP_135",
    "gated_GMP_123",
    "gated_GMP_135",
    "gated_IL7RA_CLP",
    "gated_LMPP",
    "gated_MEP_123",
    "gated_MEP_135",
]

cells_in_gates = {}

for i in gates_oi:
    print(i)
    c = Counter(ds_FL.cells.fetch(i))
    del c["False"]
    c = {k[:3]: v for k, v in c.items()}
    cells_in_gates[i] = c

    # c = ds_FL.cells.fetch("ids")[ds_FL.cells.fetch(i)]
    # c = [x.split("__")[0] for x in c]
    # c = Counter(c)
    # cells_in_gates[i] = c
cells_in_gates = pd.DataFrame.from_dict(cells_in_gates, orient="index")
cells_in_gates = cells_in_gates.T
cells_in_gates = cells_in_gates / cells_in_gates.sum() * 100
cells_in_gates = cells_in_gates.T
cells_in_gates.to_csv(out_files["FL_gating"])


cells_in_gates = {}

for i in gates_oi:
    print(i)
    c = Counter(ds_BM.cells.fetch(i))
    del c["False"]
    c = {k[:4]: v for k, v in c.items()}
    cells_in_gates[i] = c
    # c = ds_BM.cells.fetch("ids")[ds_BM.cells.fetch(i)]
    # c = [x.split("__")[0] for x in c]
    # c = Counter(c)
    # cells_in_gates[i] = c
cells_in_gates = pd.DataFrame.from_dict(cells_in_gates, orient="index")
cells_in_gates = cells_in_gates.T
cells_in_gates = cells_in_gates / cells_in_gates.sum() * 100
cells_in_gates = cells_in_gates.T
cells_in_gates.to_csv(out_files["BM_gating"])
