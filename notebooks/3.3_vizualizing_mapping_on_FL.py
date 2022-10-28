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
# %load_ext autotime
# %config InlineBackend.figure_format = 'retina'

import scarf
import pandas as pd
import glob
import numpy as np
import alphashape
from descartes import PolygonPatch
import matplotlib.pyplot as plt
import fl_utils

# -

fl_utils.set_project_path()

raw_dir = "data/raw/paper_specific/from_seurat/"
zarr_dir = "data/processed/notebooks/zarr_files"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

fl_utils.create_dir(figures_dir)

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
}

out_files = {
    "plot_yBM_on_FL": figures_dir + "/bar_plot_yBM_preds_on_FL.svg",
    "proj_plots": [
        figures_dir + f"/Figure_3B_{sample}.svg"
        for sample in ["FL_CS16", "FL_W9", "CB_HPC", "yBM_HPC"]
    ],
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

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

ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")

# +
color_palette_BM = {}
x = pd.crosstab(
    ds_yBM.cells.fetch("RNA_cluster_colors"),
    ds_yBM.cells.fetch("RNA_cluster_names"),
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

targets_reodered = ["FL_CS16", "FL_W9", "CB_HPC", "yBM_HPC"]
for sample in targets_reodered:
    print(sample)
    for g, ms in ds_FL.get_mapping_score(
        target_name=sample, log_transform=False, weighted=False, multiplier=1e5
    ):
        q = np.percentile(ms[ms > 0], 95)
        ms[ms > q] = q
        fn = [x for x in out_files["proj_plots"] if sample in x][0]
        ds_FL.plot_layout(
            layout_key="RNA_seurat_UMAP",
            color_by="RNA_cluster_names",
            color_key=color_palette_FL,
            size_vals=ms,
            height=7,
            width=7,
            legend_onside=False,
            scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.8},
            savename=fn,
            save_dpi=300,
        )

yBM_classified_clusts = {}
yBM_clusts = ds_yBM.cells.fetch("RNA_cluster_names")
for i in sorted(set(yBM_clusts)):
    yBM_classified_clusts[i] = ds_FL.get_target_classes(
        target_name="yBM_HPC",
        reference_class_group="RNA_cluster_names",
        threshold_fraction=0.4,
        target_subset=list(np.where(yBM_clusts == i)[0]),
    )

# +
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

fig, ax = plt.subplots(1, 1, figsize=(10, 4))
yBM_clusts_pred_table = {}
for n, i in enumerate(clust_order):
    df = yBM_classified_clusts[i].value_counts(ascending=True)
    vals = 100 * df / df.sum()
    yBM_clusts_pred_table[i] = vals
    vals = vals.cumsum()[::-1]
    for k, v in vals.to_dict().items():
        color = color_palette_FL[k] if k != "NA" else "#000000"
        ax.bar([n], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order], rotation=90)
ax.set_ylabel("% of projected cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["plot_yBM_on_FL"], dpi=300)

plt.show()
yBM_clusts_pred_table = pd.DataFrame(yBM_clusts_pred_table).fillna(0)
yBM_clusts_pred_table
