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

fl_utils.set_project_path()

zarr_dir = "data/processed/notebooks/zarr_files"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
}

out_files = {
    "sup_fig_2e": figures_dir + "/Sup_Figure_2e.svg",
    "sup_fig_2c": figures_dir + "/Sup_Figure_2c.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

ds_yBM = scarf.DataStore(in_files["zarr_yBM"])
ds_FL = scarf.DataStore(in_files["zarr_FL"])

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

for sample in ["FL_CS16", "FL_HPC", "FL_W9", "CB_HPC"]:
    print(sample)
    for g, ms in ds_yBM.get_mapping_score(
        target_name=sample, log_transform=False, weighted=False, multiplier=1e5
    ):
        q = np.percentile(ms[ms > 0], 95)
        ms[ms > q] = q
        fn = (
            out_files["sup_fig_2c"].split(".")[0]
            + f"_{sample}."
            + out_files["sup_fig_2c"].split(".")[1]
        )
        print(fn)
        ds_yBM.plot_layout(
            layout_key="RNA_seurat_UMAP",
            color_by="RNA_cluster_names",
            color_key=color_palette_BM,
            size_vals=ms,
            height=7,
            width=7,
            legend_onside=False,
            scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.8},
            savename=fn,
            save_dpi=300,
        )
        # savename=f'./images/Sup_Figure_2c_{sample}.svg', save_dpi=300)

FL_classified_clusts = {}
FL_clusts = ds_FL.cells.fetch("RNA_cluster_names")
for i in sorted(set(FL_clusts)):
    FL_classified_clusts[i] = ds_yBM.get_target_classes(
        target_name="FL_HPC",
        reference_class_group="RNA_cluster_names",
        threshold_fraction=0.4,
        target_subset=list(np.where(FL_clusts == i)[0]),
    )

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

fig, ax = plt.subplots(1, 1, figsize=(6, 3))
FL_clusts_pred_table = {}
for n, i in enumerate(clust_order):
    df = FL_classified_clusts[i].value_counts(ascending=True)
    vals = 100 * df / df.sum()
    FL_clusts_pred_table[i] = vals
    vals = vals.cumsum()[::-1]
    for k, v in vals.to_dict().items():
        color = color_palette_BM[k] if k != "NA" else "#000000"
        ax.bar([n], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order], rotation=90)
ax.set_ylabel("% of projected cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["sup_fig_2e"], format="svg", dpi=300)
plt.show()
FL_clusts_pred_table = pd.DataFrame(FL_clusts_pred_table).fillna(0)
FL_clusts_pred_table
