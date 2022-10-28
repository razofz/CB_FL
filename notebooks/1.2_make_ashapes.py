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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import alphashape
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pickle
from descartes import PolygonPatch
import fl_utils

fl_utils.set_project_path()

in_dir = "data/processed/notebooks"
out_dir = in_dir
figures_dir = out_dir + "/Figures"

fl_utils.create_dir(figures_dir)

# +
in_files = {
    "metadata": out_dir + "/FL_seurat_metadata.csv",
}

out_files = {
    "shapes": out_dir + "/shapes.pickle",
    "plot_num": figures_dir + "/umap_hulls_labelled_num.svg",
    "plot_cluster": figures_dir + "/umap_hulls_labelled_cluster.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")


# -


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


def make_shape(xvals, yvals, levels, ashape_alpha=0.15):
    x, y, z = density_estimation(xvals, yvals)
    cset = plt.contour(x, y, z, levels=levels, linewidths=0)
    seg = cset.allsegs[1]  # choose the outermost contourset
    seg = seg[
        np.argmax(np.array([len(x) for x in seg]))
    ]  # choose the biggest of outer closed loops
    poly = Polygon(np.vstack([seg[:, 0], seg[:, 1]]).T)
    newx, newy = [], []
    for i, j in zip(xvals, yvals):
        p = Point(i, j)
        if poly.contains(p):
            newx.append(i)
            newy.append(j)
    ashape = alphashape.alphashape(np.vstack([newx, newy]).T, ashape_alpha)
    return ashape


def make_hulls(metadata):
    levels = {
        0: 2,
        1: 5,
        2: 5,
        3: 4,
        4: 6,
        5: 20,
        6: 4,
        7: 8,
        8: 10,
        9: 0,
        10: 50,
    }
    ashapes = {}
    for i in sorted(metadata.cluster.unique()):
        idx = metadata.cluster == i
        x = metadata[idx]["u1"].values
        y = metadata[idx]["u2"].values
        ashapes[i] = make_shape(x, y, levels=levels[i])
    plt.close()
    return ashapes


metadata = pd.read_csv(in_files["metadata"], index_col=0)
ashapes = make_hulls(metadata)
pickle.dump(ashapes, open(out_files["shapes"], "wb"))
# ashapes = pickle.load(open(out_files["shapes"], "rb"))

# +
def initfig():
    return plt.subplots(1, 1, figsize=(6, 6))


def cleanfig(ax):
    padding = 5
    ax.set_xlim((metadata.u1.min() - padding, metadata.u1.max() + padding))
    ax.set_ylim((metadata.u2.min() - padding, metadata.u2.max() + padding))
    ax.set_axis_off()
    clean_axis(ax)


# -

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


# +
def make_hull_clust_colored_labeled():
    fig, ax = initfig()
    for i in sorted(metadata.cluster.unique()):
        idx = metadata.cluster == i
        x = metadata[idx]["u1"].values
        y = metadata[idx]["u2"].values
        color = metadata["clust_color"][idx].values[0]
        ax.add_patch(
            PolygonPatch(ashapes[i], alpha=0.7, color=color, ec="k", lw=0.6)
        )
        ax.text(
            x.mean(),
            y.mean(),
            metadata[idx].cluster.values[0],
            ha="center",
            va="center",
            fontsize=11,
        )
    cleanfig(ax)
    plt.savefig(out_files["plot_num"], dpi=300)
    plt.show()


make_hull_clust_colored_labeled()


# +
def make_hull_clust_colored_labeled():
    fig, ax = initfig()
    for i in sorted(metadata.cluster.unique()):
        idx = metadata.cluster == i
        x = metadata[idx]["u1"].values
        y = metadata[idx]["u2"].values
        color = metadata["clust_color"][idx].values[0]
        ax.add_patch(
            PolygonPatch(ashapes[i], alpha=0.7, color=color, ec="k", lw=0.6)
        )
        ax.text(
            x.mean(),
            y.mean(),
            metadata[idx].cluster_names.values[0],
            ha="center",
            va="center",
            fontsize=11,
        )
    cleanfig(ax)
    plt.savefig(out_files["plot_cluster"], dpi=300)
    plt.show()


make_hull_clust_colored_labeled()
