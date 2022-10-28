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
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

in_files = {}

out_files = {
    "plot_bar_n_diff_genes_all_clust_up_dn": figures_dir
    + "/bar_n_diff_genes_all_clust_up_dn.svg",
}

fl_utils.print_files(out_files, "out")

# + tags=[]
numbers = """547 Cyc_neg
436 Cyc_pos
758 DC.Mono_neg
395 DC.Mono_pos
1157 GMP_neg
581 GMP_pos
1306 HSC_neg
973 HSC_pos
502 Ly.II_neg
358 Ly.II_pos
1001 Ly.I_neg
616 Ly.I_pos
1190 MEP_neg
606 MEP_pos
1130 MPP.II_neg
739 MPP.II_pos
1205 MPP.I_neg
941 MPP.I_pos"""

# + tags=[]
tmp = [x.split(" ") for x in numbers.split("\n")]
values = {x[1]: int(x[0]) - 1 for x in tmp}

tot = []
up = []
dn = []
cluster_names = []
for i in range(0, len(values.keys()), 2):
    # print(i)
    cluster_name = list(values.keys())[i].split("_")[0]
    print(cluster_name)
    cluster_names.append(cluster_name)
    # print(list(values.items())[i:i+2])
    cluster_counts = list(values.values())[i: i + 2]
    print(f"{cluster_counts=}")
    dn.append(cluster_counts[0])
    up.append(cluster_counts[1])
    tot.append(sum(cluster_counts))

print(dn, up, tot)

# + tags=[]
n_diff = pd.DataFrame(
    data={"tot": tot, "up": up, "dn": dn}, index=cluster_names
)
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
# plt.savefig(out_dir + "raz_plot_bar_n_diff_genes_all_clust_up_dn.svg", dpi=300)
plt.tight_layout()
plt.show()
