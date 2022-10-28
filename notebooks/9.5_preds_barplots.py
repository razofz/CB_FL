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

# + tags=[]
threads = 16

figures_dir = "data/processed/notebooks/Figures"
tables_dir = "data/processed/notebooks/Supplemental_tables"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"

# + tags=[]
in_files = {
    "BM_preds": tables_dir + "/BM_preds_sample_contribution.csv",
    "FL_preds": tables_dir + "/FL_preds_sample_contribution.csv",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_BM": zarr_dir + "/RNA_yBM_hpc_merged",
}

out_files = {
    "BM_preds_plot": figures_dir + "/bar_plot_BM_preds_contribution.svg",
    "FL_preds_plot": figures_dir + "/bar_plot_FL_preds_contribution.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# + tags=[]
plt.style.use("fivethirtyeight")
plt.rcParams["svg.fonttype"] = "none"


# + tags=[]
def clean_axis(ax, ts=11, ga=0.4):
    ax.xaxis.set_tick_params(labelsize=ts)
    ax.yaxis.set_tick_params(labelsize=ts)
    for i in ["top", "bottom", "left", "right"]:
        ax.spines[i].set_visible(False)
    ax.grid(which="major", linestyle="--", alpha=ga)
    ax.figure.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    return True


# + tags=[]
ds_BM = scf.DataStore(in_files["zarr_BM"], default_assay="RNA")
ds_FL = scf.DataStore(in_files["zarr_FL"], default_assay="RNA")

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

# fig, ax = plt.subplots(1, 1, figsize=(10, 4))

bm = pd.read_csv(in_files["BM_preds"], index_col=0)
fl = pd.read_csv(in_files["FL_preds"], index_col=0)

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

fl = fl.T

fl

# labels = list(fl.index)
labels = fl.columns
fig, ax = plt.subplots()
clean_axis(ax)
fig.set_figheight(8)
fig.set_figwidth(8)
for x in fl.columns:
    df_sub = fl[x]
    df_sub = df_sub.cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_FL[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k") #, label=k)
        # ax.bar([x], [v], color=color, lw=1, edgecolor="k")
# for row in fl.index:
#     ax.bar(labels, fl.T[row]/sum(fl.T[row]), label=row)
ax.set_ylabel("Contribution")
ax.set_xticklabels(labels, rotation=90)
# ax.legend(fl.index)
plt.savefig(out_files["FL_preds_plot"])

bm = bm.T

bm

labels = bm.columns
fig, ax = plt.subplots()
clean_axis(ax)
fig.set_figheight(8)
fig.set_figwidth(8)
for x in bm.columns:
    df_sub = bm[x]
    df_sub = df_sub.cumsum()[::-1]
    for k, v in df_sub.to_dict().items():
        color = color_palette_BM[k] if k != "NA" else "#000000"
        ax.bar([x], [v], color=color, lw=1, edgecolor="k")  # , label=k)
        # ax.bar([x], [v], color=color, lw=1, edgecolor="k")
# for row in bm.index:
#     ax.bar(labels, bm.T[row]/sum(bm.T[row]), label=row)
ax.set_ylabel("Contribution")
ax.set_xticklabels(labels, rotation=90)
# ax.legend(bm.index)
plt.savefig(out_files["BM_preds_plot"])

# # labels = bm.columns
# fig, ax = plt.subplots()
# fig.set_figheight(8)
# fig.set_figwidth(8)
# for row in bm.index:
#     ax.bar(labels, bm.T[row], label=row)
# ax.set_ylabel("Contribution")
# ax.set_xticklabels(labels, rotation = 90)
# ax.legend()
# plt.savefig(out_files["BM_preds_plot"])
