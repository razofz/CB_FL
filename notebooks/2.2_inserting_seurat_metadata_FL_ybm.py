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

import scarf
import pandas as pd
import fl_utils

fl_utils.set_project_path()

raw_dir = "data/raw/paper_specific/from_seurat/"
in_dir = "data/processed/notebooks/zarr_files"
out_dir = "data/processed/notebooks"
figures_dir = out_dir + "/Figures"

# +
in_files = {
    "zarr_yBM": in_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": in_dir + "/RNA_FL_hpc_merged",
    "metadata_yBM": raw_dir + "/yBM_hpc_cite_seq_celldata_seurat.csv",
    "metadata_FL": "data/processed/notebooks/FL_seurat_metadata.csv",
}

out_files = {
    "umap_FL": figures_dir + "/UMAP_FL.svg",
    "umap_yBM": figures_dir + "/UMAP_yBM.svg",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")


def fix_index_BM(indexv):
    ret_val = []
    for i in indexv:
        x = i.rsplit("-", 1)
        ret_val.append(f"ybm{x[1][0]}__{x[0]}")
    return ret_val


def fix_index_FL(indexv):
    ret_val = []
    for i in indexv:
        x = i.rsplit("-", 1)
        ret_val.append(f"fl{x[1][0]}__{x[0]}")
    return ret_val


seurat_md_fn_BM = in_files["metadata_yBM"]
seurat_filtered_cells_BM = fix_index_BM(
    pd.read_csv(seurat_md_fn_BM, index_col=0).index
)
len(seurat_filtered_cells_BM)

seurat_md_fn_FL = in_files["metadata_FL"]
seurat_filtered_cells_FL = fix_index_FL(
    pd.read_csv(seurat_md_fn_FL, index_col=0).index
)
len(seurat_filtered_cells_FL)

ds_yBM.cells.update_key(
    ds_yBM.cells.index_to_bool(
        ds_yBM.cells.get_index_by(seurat_filtered_cells_BM, column="ids")
    ),
    key="I",
)
ds_yBM.cells.fetch_all("I").sum()

ds_FL.cells.update_key(
    ds_FL.cells.index_to_bool(
        ds_FL.cells.get_index_by(seurat_filtered_cells_FL, column="ids")
    ),
    key="I",
)
ds_FL.cells.fetch_all("I").sum()


def import_seurat_column_BM(fn, colname, save_key, fill_val):
    col = pd.read_csv(fn, index_col=0)[colname]
    col.index = fix_index_BM(col.index)
    col = col.reindex(ds_yBM.cells.fetch("ids")).values
    ds_yBM.cells.insert(save_key, col, overwrite=True, fill_value=fill_val)
    return None


def import_seurat_column_FL(fn, colname, save_key, fill_val):
    col = pd.read_csv(fn, index_col=0)[colname]
    col.index = fix_index_FL(col.index)
    col = col.reindex(ds_FL.cells.fetch("ids")).values
    ds_FL.cells.insert(save_key, col, overwrite=True, fill_value=fill_val)
    return None


import_seurat_column_BM(seurat_md_fn_BM, "clusts", "RNA_seurat_clusters", -1)
import_seurat_column_BM(
    seurat_md_fn_BM, "clust_names", "RNA_cluster_names", "NA"
)
import_seurat_column_BM(
    seurat_md_fn_BM, "clust_color", "RNA_cluster_colors", "NA"
)
import_seurat_column_BM(seurat_md_fn_BM, "u1", "RNA_seurat_UMAP1", -1)
import_seurat_column_BM(seurat_md_fn_BM, "u2", "RNA_seurat_UMAP2", -1)

import_seurat_column_FL(seurat_md_fn_FL, "cluster", "RNA_seurat_clusters", -1)
import_seurat_column_FL(
    seurat_md_fn_FL, "cluster_names", "RNA_cluster_names", "NA"
)
import_seurat_column_FL(
    seurat_md_fn_FL, "clust_color", "RNA_cluster_colors", "NA"
)
import_seurat_column_FL(seurat_md_fn_FL, "u1", "RNA_seurat_UMAP1", -1)
import_seurat_column_FL(seurat_md_fn_FL, "u2", "RNA_seurat_UMAP2", -1)

# plot UMAPs
color_palette_FL = {}
x = pd.crosstab(
    ds_FL.cells.fetch("RNA_cluster_colors"),
    ds_FL.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_FL[i] = x[i][x[i] > 0].index[0]

color_palette_yBM = {}
x = pd.crosstab(
    ds_yBM.cells.fetch("RNA_cluster_colors"),
    ds_yBM.cells.fetch("RNA_cluster_names"),
)
for i in x:
    color_palette_yBM[i] = x[i][x[i] > 0].index[0]

ds_FL.plot_layout(
    from_assay="RNA",
    layout_key="RNA_seurat_UMAP",
    color_by="RNA_cluster_names",
    color_key=color_palette_FL,
    savename=out_files["umap_FL"],
)
ds_yBM.plot_layout(
    from_assay="RNA",
    layout_key="RNA_seurat_UMAP",
    color_by="RNA_cluster_names",
    color_key=color_palette_yBM,
    savename=out_files["umap_yBM"],
)
