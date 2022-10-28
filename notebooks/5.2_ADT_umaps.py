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

import fl_utils
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"
adt_matrices_dir = out_dir + "/adt_matrices"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "FL_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_FL_hpc.csv",
    "yBM_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_yBM_hpc.csv",
    "metadata_yBM": seurat_dir + "/yBM_hpc_cite_seq_celldata_seurat.csv",
    "metadata_FL": out_dir + "/FL_seurat_metadata.csv",
}

out_files = {
    "plot_markers": f"{figures_dir}/adt_umaps_conventional_markers.svg",
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


# +
def mmnorm(v: pd.Series, f: int) -> pd.Series:
    return (f * ((v - v.min()) / (v.max() - v.min())).round(1)).astype(int)


def quantile_trim(v: pd.Series) -> pd.Series:
    q1: float = v.quantile(q=0.05)
    q2: float = v.quantile(q=0.95)
    v[v < q1] = q1
    v[v > q2] = q2
    return v


def make_scatter(ax, data, cmap):
    v = mmnorm(quantile_trim(data), len(cmap) - 1).sort_values()
    # Reindex to bring values to the top
    ax.scatter(
        metadata.u1.reindex(v.index),
        metadata.u2.reindex(v.index),
        c=[cmap[x] for x in v.values],
        lw=0,
        s=10,
        rasterized=True,
    )
    clean_axis(ax)
    ax.set_axis_off()
    return True


def make_scatter_bm(ax, data, cmap):
    v = mmnorm(quantile_trim(data), len(cmap) - 1).sort_values()
    # Reindex to bring values to the top
    ax.scatter(
        metadata_bm.u1.reindex(v.index),
        metadata_bm.u2.reindex(v.index),
        c=[cmap[x] for x in v.values],
        lw=0,
        s=10,
        rasterized=True,
    )
    clean_axis(ax)
    ax.set_axis_off()
    return True


# -

adt = pd.read_csv(
    in_files["FL_adt_matrix"],
    # "./seurat_normalized_adt_matrices/ADT_scaled_scarf_naming/seurat_adt_matrices_FL_hpc.csv",
    index_col=0,
)
metadata = pd.read_csv(in_files["metadata_FL"], index_col=0)
# metadata = pd.read_csv("./from_seurat/FL_seurat_metadata.csv", index_col=0)
adt.index = [
    x.split("_")[0] + "-" + x.rsplit("_", 1)[0][-1] + "_FL_hpc"
    for x in adt.index
]

adt_bm = pd.read_csv(
    in_files["yBM_adt_matrix"],
    # "./seurat_normalized_adt_matrices/ADT_scaled_scarf_naming/seurat_adt_matrices_yBM_hpc.csv",
    index_col=0,
)
metadata_bm = pd.read_csv(in_files["metadata_yBM"], index_col=0)
# metadata_bm = pd.read_csv(
#     "../data/cell_metadata/yBM_hpc_cite_seq_celldata_seurat.csv", index_col=0
# )
adt_bm.index = [
    x.split("_")[0] + "-" + x.rsplit("_", 1)[0][-1] + "_yBM_hpc"
    for x in adt_bm.index
]

# +
cmap = sns.color_palette("coolwarm", n_colors=11)

fig = plt.figure(figsize=(14, 9))
adt_oi = [
    "ADT-CD38",
    "ADT-CD90",
    "ADT-CD49F",
    "ADT-CD45RA",
    "ADT-CD10",
    "ADT-IL7RA",
    "ADT-CD123",
    "ADT-CD135",
    "ADT-CD71",
    "ADT-CD201",
    "ADT-CD7",
]
for n, adt_name in enumerate(adt_oi):
    ax = fig.add_subplot(3, 4, (n + 1))
    make_scatter(ax, adt[adt_name], cmap)
    ax.set_title(adt_name.upper().replace("_", " "), fontsize=12)
plt.tight_layout()
plt.savefig(out_files["plot_markers"], dpi=300)
# plt.savefig(f"./images/SVG/adt_umaps_conventional_markers.svg", dpi=300)
# plt.savefig(f"./images/PNG/adt_umaps_conventional_markers.png", dpi=300)
plt.show()
