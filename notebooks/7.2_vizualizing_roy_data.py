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

import alphashape
import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scarf
from descartes import PolygonPatch

# + tags=[]
fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
figures_dir = out_dir + "/Figures"
external_dir = "data/external"
roy_dir = external_dir + "/Roy_et_al"
roy_analysis_dir = out_dir + "/Roy_et_al_analysis"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_roy": roy_analysis_dir + "/fetal_roy.zarr",
}

out_files = {
    "plots_proj": {
        g: figures_dir + f"/umap_roy_proj_{g}_FLCS22_ref.svg"
        for g in ["ABM", "FBM", "FL", "PBM", "eFL"]
    },
    "plot_preds": figures_dir + f"/bar_plot_roy_preds_on_FLCS22.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(roy_dir)
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

ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
ds_roy = scarf.DataStore(in_files["zarr_roy"])

# +
# ds_yBM = scf.DataStore("./_zarr_files/RNA_yBM_hpc_merged")
# ds_FL = scf.DataStore("./_zarr_files/RNA_FL_hpc_merged")
# ds_roy = scf.DataStore("./_zarr_files/CB_fetal_roy.zarr")

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

# + tags=[]
score_gen = ds_FL.get_mapping_score(
    target_name="roy_data",
    target_groups=ds_roy.cells.fetch("roy_tissue_type"),
    log_transform=False,
    weighted=False,
    multiplier=1e5,
)
for g, ms in score_gen:
    print(g)
    ds_FL.plot_layout(
        layout_key="RNA_seurat_UMAP",
        color_by="RNA_cluster_names",
        color_key=color_palette_FL,
        size_vals=ms,
        height=7,
        width=7,
        legend_onside=False,
        scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.8},
        savename=out_files["plots_proj"][g],
        # savename=f"./images/SVG/umap_roy_proj_{g}_FLCS22_ref.svg",
        save_dpi=300,
    )
# -

roy_classified_type = {}
roy_type = ds_roy.cells.fetch("roy_tissue_type")
for i in sorted(set(roy_type)):
    roy_classified_type[i] = ds_FL.get_target_classes(
        target_name="roy_data",
        reference_class_group="RNA_cluster_names",
        threshold_fraction=0.4,
        target_subset=list(np.where(roy_type == i)[0]),
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
cmap = ["#FFFFFF", "#D3D3D3", "#A8A8A8", "#808080", "#000000"]
sample_order = ["eFL", "FL", "FBM", "PBM", "ABM"]

samples_pred_table_roy = {}
for n, sample in enumerate(sample_order):
    df = roy_classified_type[sample].value_counts(ascending=True).drop("NA")
    samples_pred_table_roy[sample] = 100 * df / df.sum()
samples_pred_table_roy = pd.DataFrame(samples_pred_table_roy).fillna(0)
ax = samples_pred_table_roy.reindex(clust_order).plot(
    kind="bar", figsize=(10, 5), color=cmap, width=0.7, lw=1, edgecolor="k"
)
ax.legend(frameon=False)
ax.set_ylabel("% cells in cluster", fontsize=13)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["plot_preds"], dpi=300)
# plt.savefig(f"./images/bar_plot_roy_preds_on_FLCS22.svg", dpi=300)
plt.show()
