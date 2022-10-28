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
import ast
import json
import pickle
from glob import glob

import fl_utils
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scarf
import seaborn as sns
from descartes import PolygonPatch
from numba import jit
from scipy.stats import gaussian_kde
from tqdm import tqdm

# -

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
gated_dir = raw_dir + "/paper_specific/ADT_gated_populations/conventional_gates"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
adt_matrices_dir = out_dir + "/adt_matrices"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "FL_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_FL_hpc.csv",
    "yBM_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_yBM_hpc.csv",
    "csvs": {
        sample: f"{adt_matrices_dir}/FCS_files/adt_matrices_{sample}.csv"
        for sample in ["yBM_hpc", "FL_hpc"]
    },
    "shapes": out_dir + "/shapes.pickle",
    "adt_gates": glob(f"{gated_dir}/*.csv"),
}

out_files = {
    "gated": out_dir + "/RNA_identified_gate_cell.json",
    "plot_contours": figures_dir + "/FL_gated_cells_contours_conventional.svg",
    "plot_pies": figures_dir + "/FL_gated_cells_pies_conventional.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(out_adt_matrices_dir)
# fl_utils.create_dir(out_adt_matrices_dir + "/FCS_files")
# -

# ---

# +
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


@jit(nopython=True)
def matrix_diff(X, Y):
    n_adts, n_cells = Y.shape
    min_poses = np.zeros(n_cells)
    min_vals = np.zeros(n_cells)
    for i in range(n_cells):
        y = Y[:, i]
        diff = np.abs(X - y).sum(axis=1)
        mp = int(np.argsort(diff)[0])
        min_poses[i] = mp
        min_vals[i] = diff[mp]
    return min_poses, min_vals


ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")

# +
# ds_yBM = scarf.DataStore("./_zarr_files/RNA_yBM_hpc_merged")
# ds_FL = scarf.DataStore("./_zarr_files/RNA_FL_hpc_merged")
# -

cellindices = {}
cellindices["FL_hpc"] = dict(
    zip(
        [
            x.split("__")[1] + "_" + x.split("__")[0].upper() + "_hpc"
            for x in ds_FL.cells.fetch_all("ids")
        ],
        ds_FL.cells.fetch_all("ids"),
    )
)
cellindices["yBM_hpc"] = dict(
    zip(
        [
            x.split("__")[1]
            + "_"
            + x.split("__")[0].replace("ybm", "yBM")
            + "_hpc"
            for x in ds_yBM.cells.fetch_all("ids")
        ],
        ds_yBM.cells.fetch_all("ids"),
    )
)

adt_order = pd.read_csv(
    in_files["FL_adt_matrix"],
    # "./seurat_normalized_adt_matrices/ADT_scaled_scarf_naming/seurat_adt_matrices_FL_hpc.csv",
    index_col=0,
).columns


# + tags=[]
def identify_cells(diff_tol=2):
    gated_cells = {}
    for fn in tqdm(sorted(in_files["adt_gates"])):
        sample, gate = fn.split("/")[-1][:-4].split("_null_")
        sample = sample + "_hpc"
        ref_adts = pd.read_csv(
            in_files["csvs"][sample],
            # f"./seurat_normalized_adt_matrices/FCS_files/adt_matrices_{sample}.csv",
            index_col=0,
        )
        cell_names = ref_adts.index
        ref_adts = ref_adts[adt_order].values
        df = pd.read_csv(fn).T
        df = df.reindex([x.replace("_", "") for x in adt_order]).values
        cells = []
        missed = 0
        for mp, mv in zip(*matrix_diff(ref_adts, df)):
            if mv < diff_tol:
                cells.append(cell_names[mp])
            else:
                missed += 1
                # 'print' here might give too many warnings
                pass
        if len(cells) != df.shape[1]:
            print(
                f"WARNING: {missed}/{df.shape[1]} cells not found in {sample} {gate}"
            )

        if gate not in gated_cells:
            gated_cells[gate] = {}
        if sample not in gated_cells[gate]:
            gated_cells[gate][sample] = []
        missed = 0
        none_vals = 0
        for i in cells:
            if i in cellindices[sample]:
                if cellindices[sample][i] is None:
                    none_vals += 1
                else:
                    gated_cells[gate][sample].append(cellindices[sample][i])
            else:
                missed += 1
        if missed > 0:
            print(
                f"INFO {sample} {gate}: {missed}/{len(cells)} were not found in cellindices"
            )
        if none_vals > 0:
            print(
                f"INFO {sample} {gate}: {none_vals}/{len(cells)} had None values in cellindices"
            )
    return gated_cells


gated_cells = identify_cells()
with open(out_files["gated"], "w") as OUT:
    # with open("./_zarr_files/RNA_identified_gate_cell.json", "w") as OUT:
    json.dump(gated_cells, OUT, indent=2)
gated_cells.keys()
# -

gated_cells = json.load(open(out_files["gated"]))
# gated_cells = json.load(open("./_zarr_files/RNA_identified_gate_cell.json"))
gated_cells.keys()

ashapes_FL = pickle.load(open(in_files["shapes"], "rb"))
# ashapes_FL = pickle.load(open("_zarr_files/shapes.pickle", "rb"))
# ashapes_BM = pickle.load(open("../data/_zarr_files/shapes.pickle", "rb" ))

# +
def density_estimation(x, y):
    deltaX = max(x) - min(x)
    deltaY = max(y) - min(y)
    xmin = min(x) - deltaX
    xmax = max(x) + deltaX
    ymin = min(y) - deltaY
    ymax = max(y) + deltaY
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    z = np.reshape(kernel(positions).T, xx.shape)
    return xx, yy, z


def contour_plot(ax, xvals, yvals, levels=12, ls=":", c="k", lw=2, zorder=10):
    x, y, z = density_estimation(xvals, yvals)
    cset = ax.contour(x, y, z, levels=levels, linewidths=0)
    for j in range(len(cset.allsegs)):
        color = c
        linew = lw
        for ii, seg in enumerate(cset.allsegs[j]):
            ax.plot(
                seg[:, 0], seg[:, 1], ls=ls, c=color, lw=linew, zorder=zorder
            )


def hull_plot(ax, ashape, disp_alpha, color, lw=1):
    ax.add_patch(
        PolygonPatch(ashape, alpha=disp_alpha, color=color, ec="grey", lw=lw)
    )


# -

metadata_FL = ds_FL.cells.to_pandas_dataframe(
    [
        "ids",
        "RNA_seurat_UMAP1",
        "RNA_seurat_UMAP2",
        "RNA_seurat_clusters",
        "RNA_cluster_names",
        "RNA_cluster_colors",
    ]
)
metadata_FL.index = metadata_FL.ids

# + tags=[]
gates = [
    "CD49Fpos_HSC",
    "CD90pos_HSC",
    "CD38neg_MPP",
    "LMPP",
    "CMP_123",
    "CMP_135",
    "GMP_123",
    "GMP_135",
    "MEP_123",
    "MEP_135",
    "CD10_CLP",
    "CD7_CLP",
    "IL7RA_CLP",
]

fig = plt.figure(figsize=(12, 12))
for n, gate in enumerate(gates):

    cells = gated_cells[gate]["FL_hpc"]
    if len(cells) < 1:
        print(f"WARNING: No cells found in gate {gate} for sample {sample}")
        continue
    x = metadata_FL.reindex(
        cells
    ).RNA_seurat_UMAP1  # Even if there are duplicate cell names they are included
    y = metadata_FL.reindex(cells).RNA_seurat_UMAP2
    ax = fig.add_subplot(4, 4, n + 1)
    title = f"{gate.replace('_', ' ')} (n={len(cells)})"
    ax.set_title(title, fontsize=12)

    for i in metadata_FL.RNA_seurat_clusters.unique():
        if i == -1:
            continue
        idx = metadata_FL.RNA_seurat_clusters == i
        hull_plot(
            ax,
            ashapes_FL[i],
            0.2,
            metadata_FL[metadata_FL.RNA_seurat_clusters == i][
                "RNA_cluster_colors"
            ].values[0],
            1.5,
        )
    if len(cells) > 5:
        contour_plot(ax, x, y, ls="--", lw=0.8, zorder=5, levels=4, c="k")
    gate_clean_name = gate.replace("n", "-").replace("p", "+").replace("_", "")

    padding = 2
    ax.set_xlim(
        (
            min(metadata_FL.RNA_seurat_UMAP1) - padding,
            max(metadata_FL.RNA_seurat_UMAP1) + padding,
        )
    )
    ax.set_ylim(
        (
            min(metadata_FL.RNA_seurat_UMAP2) - padding,
            max(metadata_FL.RNA_seurat_UMAP2) + padding,
        )
    )
    clean_axis(ax)
    ax.set_axis_off()

plt.tight_layout()
plt.savefig(out_files["plot_contours"], dpi=300)
# plt.savefig(f"./images/PNG/FL_gated_cells_contours_conventional.png", dpi=300)
plt.show()

# + tags=[]
gates = [
    "CD49Fpos_HSC",
    "CD90pos_HSC",
    "CD38neg_MPP",
    "LMPP",
    "CMP_123",
    "CMP_135",
    "GMP_123",
    "GMP_135",
    "MEP_123",
    "MEP_135",
    "CD10_CLP",
    "CD7_CLP",
    "IL7RA_CLP",
]

fig = plt.figure(figsize=(12, 12))
for n, gate in enumerate(gates):
    cells = gated_cells[gate]["FL_hpc"]
    if len(cells) < 1:
        print(f"WARNING: No cells found in gate {gate} for sample {sample}")
        continue
    ax = fig.add_subplot(4, 4, n + 1)
    title = f"{gate.replace('_', ' ')} (n={len(cells)})"
    ax.set_title(title, fontsize=12)
    v = metadata_FL.reindex(cells).RNA_cluster_names.value_counts()
    idx = v / v.sum() < 0.05
    v["Rest"] = v[idx].sum()
    v = v.drop(idx[idx].index)
    colors = [
        metadata_FL.RNA_cluster_colors[
            metadata_FL.RNA_cluster_names == x
        ].values[0]
        for x in v.index[:-1]
    ] + [[0.5, 0.5, 0.5]]
    explode = [0.05 if x == 0 else 0.02 for x in range(len(colors))]
    ax.pie(
        v,
        labels=v.index,
        autopct="%1.1f%%",
        startangle=90,
        colors=colors,
        explode=explode,
        wedgeprops={"alpha": 0.7},
        textprops={"fontsize": 10},
    )
    clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["plot_pies"], dpi=300)
# plt.savefig(f"./images/PNG/FL_gated_cells_pies_conventional.png", dpi=300)
plt.show()
