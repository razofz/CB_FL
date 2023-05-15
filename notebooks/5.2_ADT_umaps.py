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

import fl_utils
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scarf
import numpy as np

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"
zarr_dir = "data/processed/notebooks/zarr_files"
adt_matrices_dir = out_dir + "/adt_matrices"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "FL_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_FL_hpc.csv",
    "yBM_adt_matrix": f"{adt_matrices_dir}/seurat_adt_matrices_yBM_hpc.csv",
    "metadata_yBM": seurat_dir + "/yBM_hpc_cite_seq_celldata_seurat.csv",
    "metadata_FL": out_dir + "/FL_seurat_metadata.csv",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
}

out_files = {
    "plot_markers": f"{figures_dir}/adt_umaps_conventional_markers.svg",
    "dot_plot_ADT": f"{figures_dir}/adt_dot_plot.svg",
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

ds_FL = scarf.DataStore(in_files["zarr_FL"], nthreads=16)

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

# -

# ---

adt_gene_names = {
    'CD35': 'CR1',
    'CD90': 'THY1',
    'CD201': 'PROCR',
    'CD32': 'FCGR2A',
    'CD4': 'CD4',
    'CD11C': 'ITGAX',
    'CD93': 'CD93',
    'CD62L': 'SELL',
    'CD49F': 'ITGA6',
    'CD54': 'ICAM1',
    'CD105': 'ENG',
    'CD117': 'KIT',
    'CD155': 'PVR',
    'CD52': 'CD52',
    'CD95': 'FAS',
    'CD33': 'CD33',
    'CD9': 'CD9',
    'CD36': 'CD36',
    'CD107A': 'LAMP1',
    'CD196': 'CCR6',
    'CD70': 'CD70',
    'CD123': 'IL3RA',
    'CD45RO': 'PTPRC',
    'CD25': 'IL2RA',
    'CD26': 'DPP4',
    'CD44': 'CD44',
    'CD48': 'CD48',
    'CD135': 'FLT3',
    'CD42B': 'GP1BA',
    'CD79B': 'CD79B',
    'CD41': 'ITGA2B',
    'IL7RA': 'IL7R',    
    'CD56': 'NCAM1',
    'CD61': 'ITGB3',
    'CD10': 'MME',
    'CD38': 'CD38',
    'BAH1': 'MPL',
    'CD49D': 'ITGA4',
    'CD71': 'TFRC',
    'CD18': 'ITGB2',
    'CD45RA': 'PTPRC',
    'CD7': 'CD7',
    'CD11A': 'ITGAL',
    'CD34': 'CD34',
}


adt_gene_names = {k: v for k, v in adt_gene_names.items().__reversed__()}


# +
fig, ax = plt.subplots(1, 1, figsize=(6, 12))
clust_order = ["Ly-III", "DC-I", "Cyc", "MEP", "DC-Mono", "GMP", "Ly-II", "Ly-I", "MPP-II", "MPP-I", "HSC"]
clust_order.reverse()
pos = 0
names = []

genes = ds_FL.ADT.feats.fetch_all('names')

feat_idx = ds_FL.ADT.feats.get_index_by(genes, column='names')
df = pd.DataFrame(ds_FL.ADT.normed(feat_idx=feat_idx).compute(), columns=genes)
df.index = [x.split('__')[1] + '-' + x[2] + '_FL_hpc' for x in ds_FL.cells.fetch('ids')]
df = df.reindex(metadata.index)

mdf = df.copy()
mdf['cluster'] = metadata.cluster_names.copy()
mdf = mdf.groupby('cluster').mean()
mdf = mdf.reindex(clust_order)
mdf = mdf[["ADT_" + x for x in adt_gene_names if x != "CD34"]]

bdf = df.copy()
bdf['cluster'] = metadata.cluster_names.copy()
bdf = bdf.groupby('cluster').mean()
bdf = bdf.reindex(clust_order)
bdf = bdf[["ADT_" + x for x in adt_gene_names if x != "CD34"]]

print (mdf.shape, bdf.shape)

n = 11
for j1, j2 in zip(mdf, bdf):
    if len(bdf[j2].shape) == 2:
        c = bdf[j2].mean(axis=1).values
    else:
        c = bdf[j2].values
    c = (c - c.mean()) / c.std()
    ax.scatter(range(n), np.repeat(pos, n), s=100*mdf[j1], c=c, vmin=-2, vmax=2,
               cmap='bwr', alpha=0.7, lw=0.3, edgecolors='k', rasterized=True)
    names.append(j1)
    pos += 1
clean_axis(ax, ga=0.1)        
ax.set_yticks(range(len(names)))
ax.set_yticklabels([ f"{x[4:]} ({adt_gene_names[x[4:]]})" for x in names], fontsize=10, ha='right')
ax.set_xticks(range(11))
ax.set_xticklabels(clust_order, fontsize=12, rotation = 45)
ax.set_ylim((-1, len(names)))
ax.set_xlim((-1, 11))
plt.tight_layout()
plt.savefig(out_files["dot_plot_ADT"], dpi=300)
plt.show()
