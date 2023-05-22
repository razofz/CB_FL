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
# %load_ext autotime
# %config InlineBackend.figure_format = 'retina'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fl_utils
import seaborn as sns


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
fl_utils.set_project_path()

# + tags=[]
in_dir = "data/processed/DESeq2/results/sub_setter/FLcorePseudotech/"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

fl_clusters = [
    "Cyc",
    "DC.Mono",
    "GMP",
    "HSC",
    "Ly.I",
    "Ly.II",
    "MEP",
    "MPP.I",
    "MPP.II",
]

in_files = {
    "csvs": {
        cl: {"pos": in_dir + cl + "_pos" + ".csv", "neg": in_dir + cl + "_neg" + ".csv"}
        for cl in fl_clusters
    }
}

out_files = {
    "plot_bar_n_diff_genes_all_clust_up_dn": figures_dir
    + "/bar_n_diff_genes_all_clust_up_dn.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

sizes = {cl: {"pos": 0, "neg": 0} for cl in fl_clusters}
for cl in fl_clusters:
    for reg in ["pos", "neg"]:
        with open(in_files["csvs"][cl][reg], "r") as f:
            sizes[cl][reg] = len(f.readlines()) - 1

n_diff = pd.DataFrame(sizes).T
n_diff["tot"] = n_diff.pos + n_diff.neg
n_diff.rename(columns={"neg": "dn", "pos": "up"}, inplace=True)
n_diff

# + tags=[]
fig, ax = plt.subplots(1, 1, figsize=(5, 3))
colors = {"tot": "r", "dn": "b"}
n_diff_ord = n_diff.sort_values(by=["tot"], ascending=False)
ax.bar(n_diff_ord.index, n_diff_ord["tot"], color=colors["tot"], edgecolor="k")
ax.bar(n_diff_ord.index, n_diff_ord["dn"], color=colors["dn"], edgecolor="k")

ax.set_xticks(np.array(range(len(n_diff_ord.index))) + 0.25)
ax.set_xticklabels(
    [x.replace(".", "-") for x in n_diff_ord.index],
    rotation=45,
    ha="right",
    fontsize=12,
)
ax.set_ylabel("n genes", fontsize=12)
labels = list(colors.keys())
handles = [plt.Rectangle((0, 0), 1, 1, color=colors[label]) for label in labels]
labels = ["up", "down"]
plt.legend(handles, labels, loc="upper right")
clean_axis(ax)
plt.savefig(out_files["plot_bar_n_diff_genes_all_clust_up_dn"], dpi=300)
plt.tight_layout()
plt.show()
