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

import glob
import os

import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scarf as scf
from IPython.display import display

scf.__version__

# + tags=[]
fl_utils.set_project_path()

# +
threads = 16

out_dir = "data/processed/notebooks"
external_dir = "data/external/E-MTAB-9067"
raw_dir = "data/raw/paper_specific/from_seurat"
figures_dir = "data/processed/notebooks/Figures"
zarr_dir = out_dir + "/zarr_files"

# + tags=[]
in_files = {
    "h5": external_dir + "/MergedAllSamples_annotated.h5ad",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "civij_clusters": raw_dir + "/../civij_clusters.txt",
}

# + tags=[]
with open(in_files["civij_clusters"]) as f:
    civij_clusters = f.readlines()
    civij_clusters = [x.rstrip() for x in civij_clusters]

# + tags=[]
out_files = {
    "zarr_civic": zarr_dir + "/RNA_civicdata",
    "umap_FL_ref": figures_dir + "/umap_ms_FL_ref_FL_Civij.svg",
    "barplot_civij_preds": figures_dir + "/bar_plot_civij_preds_on_FLCS22.svg",
    "barplot_civij_preds_on_FL": figures_dir + "/bar_plot_civij_preds_on_FL.svg",
    "roy_projs": {
        g: figures_dir + f"/umap_roy_proj_{g}_FLCS22_ref.svg" for g in civij_clusters
    },
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
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

reader = scf.H5adReader(
    in_files["h5"],
    matrix_key="raw/X",
    cell_ids_key="CELL_NAME",
    feature_ids_key="Ensembl",
    feature_name_key="_index",
)

reader.groupCodes["X"] = reader.groupCodes["raw/X"]

writer = scf.H5adToZarr(reader, zarr_fn=out_files["zarr_civic"])
writer.dump()

ds = scf.DataStore(out_files["zarr_civic"], default_assay="RNA")

ds.filter_cells(
    attrs=["RNA_nCounts", "RNA_nFeatures", "RNA_percentMito"],
    highs=[40000, 8000, 20],
    lows=[1000, 500, 1],
)

ds.mark_hvgs(min_cells=20, top_n=2000, max_mean=2, min_mean=-3, max_var=10)

ds.make_graph(feat_key="hvgs", k=21, dims=15, n_centroids=100)

# + tags=[]
ds.run_umap(n_epochs=200, spread=2, min_dist=2, parallel=True)
# -

ds.plot_layout(layout_key="RNA_UMAP", color_by="gate")
ds.plot_layout(layout_key="RNA_UMAP", color_by="CD34")
ds.plot_layout(layout_key="RNA_UMAP", color_by="SPINK2")

ds.plot_layout(layout_key="RNA_UMAP", color_by="gate")

ds.run_leiden_clustering(resolution=1)

ds.plot_layout(layout_key="RNA_UMAP", color_by="RNA_leiden_cluster")

ds.run_marker_search(group_key="RNA_leiden_cluster", threshold=0.25)

# + tags=[]
ds.plot_marker_heatmap(group_key="RNA_leiden_cluster", topn=5)
# -

ds_FL = scf.DataStore(in_files["zarr_FL"], default_assay="RNA")

x = ds.RNA.feats.fetch_all("names")
len(x), len(set(x))

temp = {}
for i in ds.RNA.feats.fetch_all("names"):
    if i not in temp:
        temp[i] = 0
    temp[i] += 1
new_ids = []
for i in temp:
    new_ids.append(f"{i}")
    for j in range(1, temp[i]):
        new_ids.append(f"{i}_{j+1}")
len(temp), len(new_ids)

ds.RNA.feats.insert("new_ids", new_ids, overwrite=True)

# rm -rf {out_files["zarr_civic"] + "/RNA/featureData/ids"}

# !mv {out_files["zarr_civic"] + "/RNA/featureData/new_ids"} {out_files["zarr_civic"] + "/RNA/featureData/ids"}

# + tags=[]
ds.run_mapping(
    target_assay=ds_FL.RNA,
    target_name="FL_cs22",
    target_feat_key="I__hvgs",
    save_k=11,
    run_coral=True,
)

# + tags=[]
ds_FL.run_mapping(
    target_assay=ds.RNA,
    target_name="FL_Civij",
    target_feat_key="I__hvgs",
    save_k=11,
    run_coral=True,
)
# -

color_palette_FL = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"), ds_FL.cells.fetch("RNA_cluster_names")
)
for i in x:
    color_palette_FL[i] = x[i][x[i] > 0].index[0]

targets_reodered = ["FL_Civij"]
for sample in targets_reodered:
    print(sample)
    for g, ms in ds_FL.get_mapping_score(
        target_name=sample, log_transform=False, weighted=False, multiplier=1e5
    ):
        # fig, ax = plt.subplots(1, 1, figsize=(5,5.25))
        q = np.percentile(ms[ms > 0], 95)
        ms[ms > q] = q
        ds_FL.plot_layout(
            layout_key="RNA_seurat_UMAP",
            color_by="RNA_cluster_names",
            color_key=color_palette_FL,
            size_vals=ms,
            height=7,
            width=7,
            legend_onside=False,
            scatter_kwargs={"lw": 0.2, "zorder": 2, "alpha": 0.8},
        )
        plt.savefig(out_files["umap_FL_ref"], dpi=300)
        plt.show()

# + tags=[]
score_gen = ds_FL.get_mapping_score(
    target_name="FL_Civij",
    target_groups=ds.cells.fetch("gate"),
    log_transform=False,
    weighted=False,
    multiplier=1e5,
)
for i in score_gen:
    print(i[0])
score_gen = ds_FL.get_mapping_score(
    target_name="FL_Civij",
    target_groups=ds.cells.fetch("gate"),
    log_transform=False,
    weighted=False,
    multiplier=1e5,
)

# + tags=[]
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
        savename=out_files["roy_projs"][g],
        show_fig=True,
        save_dpi=300,
    )
# -

np.unique(ds.cells.fetch_all("manual_annotation"))

# +
Civij_classified_type = {}
Civij_type = ds.cells.fetch("gate")
for i in sorted(set(Civij_type)):
    Civij_classified_type[i] = ds_FL.get_target_classes(
        target_name="FL_Civij",
        reference_class_group="RNA_cluster_names",
        threshold_fraction=0.4,
        target_subset=list(np.where(Civij_type == i)[0]),
    )

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
sample_order = ["HSC", "MPP", "CMP", "MEP", "GMP", "CLP"]

samples_pred_table_civij = {}
for n, sample in enumerate(sample_order):
    df = Civij_classified_type[sample].value_counts(ascending=True).drop("NA")
    samples_pred_table_civij[sample] = 100 * df / df.sum()
samples_pred_table_civij = pd.DataFrame(samples_pred_table_civij).fillna(0)
ax = samples_pred_table_civij.reindex(clust_order).plot(
    kind="bar", figsize=(10, 5), width=0.7, lw=1, edgecolor="k"
)
ax.legend(frameon=False)
ax.set_ylabel("% cells in cluster", fontsize=13)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_civij_preds"], dpi=300)
plt.show()
samples_pred_table_civij

# +
clust_order = ["HSC", "MPP", "CMP", "MEP", "GMP", "CLP"]

fig, ax = plt.subplots(1, 1, figsize=(5, 6))
clust_pred_table_civij = {}
for n, i in enumerate(clust_order):
    df = Civij_classified_type[i].value_counts(ascending=True)
    vals = 100 * df / df.sum()
    clust_pred_table_civij[i] = vals
    vals = vals.cumsum()[::-1]
    for k, v in vals.to_dict().items():
        color = color_palette_FL[k] if k != "NA" else "#000000"
        ax.bar([n], [v], color=color, lw=1, edgecolor="k")
ax.set_xticks(range(len(clust_order)))
ax.set_xticklabels([x for x in clust_order], rotation=90)
ax.set_ylabel("% of projected cells", fontsize=12)
clean_axis(ax)
plt.tight_layout()
plt.savefig(out_files["barplot_civij_preds_on_FL"], dpi=300)
plt.show()
