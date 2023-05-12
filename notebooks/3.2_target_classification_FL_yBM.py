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
# %load_ext autotime
# %config InlineBackend.figure_format = 'retina'

import scarf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import json
import fl_utils

plt.rcParams["svg.fonttype"] = "none"


def clean_axis(ax, ts=11, ga=0.4):
    ax.xaxis.set_tick_params(labelsize=ts)
    ax.yaxis.set_tick_params(labelsize=ts)
    ax.grid(which="major", linestyle="--", alpha=ga)
    ax.figure.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    ax.set_axisbelow(True)
    return True


# + tags=[]
fl_utils.set_project_path()

# + tags=[]
raw_dir = "data/raw/paper_specific/from_seurat/"
zarr_dir = "data/processed/notebooks/zarr_files"
out_dir = "data/processed/notebooks"
sup_tables_dir = out_dir + "/Supplemental_tables"
figures_dir = out_dir + "/Figures"

# + tags=[]
fl_utils.create_dir(sup_tables_dir)
fl_utils.create_dir(out_dir + "/mapping_predictions")

# + tags=[]
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
}

out_files = {
    "preds_BM": out_dir + "/mapping_predictions/RNA_BM_target_preds.json",
    "preds_FL": out_dir + "/mapping_predictions/RNA_FL_target_preds.json",
    "best_nn_BM": out_dir + "/mapping_predictions/RNA_BM_best_nn.json",
    "best_nn_FL": out_dir + "/mapping_predictions/RNA_FL_best_nn.json",
    "table_FL_on_yBM": sup_tables_dir + "/predictions_FL_hpc_on_yBM_hpc.csv",
    "table_yBM_on_FL": sup_tables_dir + "/predictions_yBM_hpc_on_FL_hpc.csv",
    "classified_cells": out_dir + "/mapping_predictions/classified_cells.json",
    "plot_bar_ybm_hpc_on_FL_hpc": figures_dir
    + "/prediction_barplots_ybm_hpc_on_FL_hpc.svg",
    "plot_bar_FL_HPC_on_ybm_hpc": figures_dir
    + "/prediction_barplots_FL_HPC_on_ybm_hpc.svg",
}


fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

# ---

# + tags=[]
ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
# -

Counter(ds_FL.cells.fetch("RNA_cluster_names"))

# + tags=[]
color_palette_BM = {}
x = pd.crosstab(
    ds_yBM.cells.fetch("RNA_cluster_colors"),
    ds_yBM.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_BM[i] = x[i][x[i] > 0].index[0]

# + tags=[]
color_palette_FL = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"),
    ds_FL.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_FL[i] = x[i][x[i] > 0].index[0]

# + tags=[]
classified_cells_BM = {}
for sample in sorted(list(ds_yBM.z.RNA.projections.keys())):
    classified_cells_BM[sample] = ds_yBM.get_target_classes(
        target_name=sample,
        threshold_fraction=0.4,
        reference_class_group="RNA_cluster_names",
    )
    print(sample, len(classified_cells_BM[sample]))

# + tags=[]
classified_cells_FL = {}
for sample in sorted(list(ds_FL.z.RNA.projections.keys())):
    classified_cells_FL[sample] = ds_FL.get_target_classes(
        target_name=sample,
        threshold_fraction=0.4,
        reference_class_group="RNA_cluster_names",
    )
    print(sample, len(classified_cells_FL[sample]))

# + tags=[]
samples_pred_table_BM = {}
targets_reodered = ["FL_CS16", "FL_HPC", "FL_W9", "CB_HPC", "yBM_HPC"]
for n, sample in enumerate(targets_reodered):
    df = classified_cells_BM[sample].value_counts(ascending=True).drop("NA")
    samples_pred_table_BM[sample] = 100 * df / df.sum()
samples_pred_table_BM = pd.DataFrame(samples_pred_table_BM).fillna(0)
samples_pred_table_BM.to_csv(out_files["table_FL_on_yBM"])

clust_order = [
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
]
cmap = ["#FFFFFF", "#D3D3D3", "#A8A8A8", "#808080", "#000000"]

ax = samples_pred_table_BM[samples_pred_table_BM.columns[::-1]].reindex(clust_order).iloc[::-1].plot(
    kind="barh", figsize=(5, 10), color=cmap[::-1], width=0.7, lw=1, edgecolor="k"
)
ax.legend(frameon=False, bbox_to_anchor=(1, 0), loc="lower right")
ax.set_xlabel("% cells in cluster", fontsize=13)
clean_axis(ax)
plt.tight_layout()
# plt.savefig(f'./images/prediction_barplots_FL_HPC_on_ybm_hpc.svg', dpi=300)
plt.savefig(out_files["plot_bar_FL_HPC_on_ybm_hpc"], dpi=300)
plt.show()

# + tags=[]
samples_pred_table_FL = {}
targets_reodered = ["FL_CS16", "FL_HPC", "FL_W9", "CB_HPC", "yBM_HPC"]
for n, sample in enumerate(targets_reodered):
    df = classified_cells_FL[sample].value_counts(ascending=True).drop("NA")
    samples_pred_table_FL[sample] = 100 * df / df.sum()
samples_pred_table_FL = pd.DataFrame(samples_pred_table_FL).fillna(0)
# samples_pred_table_FL.to_csv(out_files["table_yBM_on_FL"])

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
cmap = ["#FFFFFF", "#D3D3D3", "#A8A8A8", "#808080", "#000000"]

ax = samples_pred_table_FL[samples_pred_table_FL.columns[::-1]].reindex(clust_order).iloc[::-1].plot(
    kind="barh", figsize=(5, 10), color=cmap[::-1], width=0.7, lw=1, edgecolor="k",
    orientation="horizontal"
)
ax.legend(frameon=False, bbox_to_anchor=(1, 0), loc="lower right")
ax.set_xlabel("% cells in cluster", fontsize=13)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["plot_bar_ybm_hpc_on_FL_hpc"], dpi=300)
# plt.savefig(f'./images/prediction_barplots_ybm_hpc_on_FL_hpc.svg', dpi=300)
plt.show()
# -

classified_cells_FL

cc_fl = {k: list(v) for k, v in classified_cells_FL.items()}

import json
with open(out_files["classified_cells"], "w") as f:
    f.write(json.dumps(cc_fl))

{k: Counter(v) for k, v in classified_cells_FL.items()}

ds_FL

Counter(ds_FL.cells.fetch("RNA_cluster_names"))

samples_pred_table_FL

classified_cells_BM

# + tags=[]
preds_BM = {}
for i in ["FL_HPC", "CB_HPC"]:
    print(i)
    ds2 = scarf.DataStore(
        zarr_dir
        + f"/RNA_{i.split('_')[0] + '_' + i.split('_')[1].lower()}_merged"
    )
    preds_BM[i] = dict(
        zip(ds2.cells.fetch("ids"), classified_cells_BM[i].values)
    )
# with open(zarr_dir + "/RNA_BM_target_preds.json", 'w') as OUT:
with open(out_files["preds_BM"], "w") as OUT:
    json.dump(preds_BM, OUT, indent=2)

best_nn_BM = {}
for target in ["FL_HPC", "CB_HPC"]:
    print(target)
    proj = ds_yBM.z.RNA.projections[target]
    ref_names = ds_yBM.cells.fetch("ids")[proj["indices"][:, 0]]
    target_names = ds2.cells.fetch("ids")
    best_nn_BM[target] = dict(zip(target_names, ref_names))
with open(out_files["best_nn_BM"], "w") as OUT:
    json.dump(best_nn_BM, OUT, indent=2)

# +
preds_FL = {}
for i in ["yBM_HPC", "CB_HPC"]:
    print(i)
    ds2 = scarf.DataStore(
        zarr_dir
        + f"/RNA_{i.split('_')[0] + '_' + i.split('_')[1].lower()}_merged"
    )
    preds_FL[i] = dict(
        zip(ds2.cells.fetch("ids"), classified_cells_FL[i].values)
    )
with open(out_files["preds_FL"], "w") as OUT:
    json.dump(preds_FL, OUT, indent=2)

best_nn_FL = {}
for target in ["yBM_HPC", "CB_HPC"]:
    print(target)
    proj = ds_FL.z.RNA.projections[target]
    ref_names = ds_FL.cells.fetch("ids")[proj["indices"][:, 0]]
    target_names = ds2.cells.fetch("ids")
    best_nn_FL[target] = dict(zip(target_names, ref_names))
with open(out_files["best_nn_FL"], "w") as OUT:
    json.dump(best_nn_FL, OUT, indent=2)
