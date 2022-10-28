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

import anndata as ad
import fl_utils
import numpy as np
import scanpy as sc
import scarf

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
in_files = {}

out_files = {
    "roy_h5ad": roy_dir + "/Roy_et_al.h5ad",
    "roy_h5ad_analysis": roy_analysis_dir + "/fetal_roy.h5ad",
    "zarr_roy": roy_analysis_dir + "/fetal_roy.zarr",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

fl_utils.create_dir(roy_dir)
fl_utils.create_dir(roy_analysis_dir)

# ---

# Download their provided AnnData file:

# !wget "https://data.mendeley.com/public-files/datasets/phfgms85x2/files/f14ae771-4aeb-4871-a7ec-5dedc1ee6679/file_downloaded" --output-document {out_files["roy_h5ad"]}.gz

# !gunzip -c {out_files["roy_h5ad"]}.gz > {out_files["roy_h5ad"]}

# !ls {out_files["roy_h5ad"]}
# !file {out_files["roy_h5ad"]}

# ---

adata = sc.read_h5ad(out_files["roy_h5ad"])

adata.X = adata.X.toarray()

adata.write_h5ad(out_files["roy_h5ad_analysis"])  # requires a bit of memory

reader = scarf.H5adReader(
    out_files["roy_h5ad_analysis"], feature_name_key="name"
)
writer = scarf.H5adToZarr(
    reader,
    zarr_fn=out_files["zarr_roy"],
    chunk_size=(1000, 1000),
    assay_name="RNA",
)
writer.dump()
reader.nCells, reader.nFeatures

# _Note:_ Change the `nthreads` parameter to something suitable for your computer. Also, depending on how much RAM you have, might want/have to change the `chunk_size` parameter above (can just remove it completely and it'll use a default value).

ds = scarf.DataStore(
    out_files["zarr_roy"], nthreads=threads, min_features_per_cell=10
)

ds

adata

ds.cells.head()

ds.cells.columns

for col in [
    "UMI_count",
    "cell_type",
    "cluster_id",
    "detectedGenesPerCell",
    "donor",
    "sampleID",
    "tissue_type",
]:
    ds.cells.insert(
        f"roy_{col}",
        np.array(ds.cells.to_pandas_dataframe(columns=[col])[col]),
        overwrite=True,
    )
    ds.cells.drop(column=col)

# Here I renamed the metadata columns that came from their AnnData object as `roy_{col}`, to make it clear which metadata column were provided by them.

roy_umap = adata.obsm["X_umap"].T
roy_umap

ds.cells.insert("roy_UMAP1", roy_umap[0], overwrite=True)
ds.cells.insert("roy_UMAP2", roy_umap[1], overwrite=True)

# Also, add their UMAP coordinates to our cell metadata table. (There is also t-SNE and force-directed graph coordinates in the AnnData object, if you want them you can add them in the same way)

ds.cells.columns

ds.cells.head()

# Notice that the suffix after the barcodes (e. g. `AAACCTGAGAATGTTG-1_eFL`) signifies the tissue type each cell came from. Notice also that the suffix matches with the `roy_tissue_type` column value.

ds.RNA.feats.head()

ds.plot_cells_dists()

# _Note:_ Probably already preprocessed, e. g. doesn't look to be a lot of outliers. So I don't filter out anything here, feel free if you want to set any thresholds.

# _Preprocessing done._
#
# ---
#
# Look at their UMAP with their metadata:

ds.plot_layout(layout_key="roy_UMAP")

ds.plot_layout(layout_key="roy_UMAP", color_by="roy_cluster_id")

ds.plot_layout(layout_key="roy_UMAP", color_by="roy_tissue_type")

ds.plot_layout(layout_key="roy_UMAP", color_by="roy_cell_type")

ds.plot_layout(layout_key="roy_UMAP", color_by="roy_donor")

ds.plot_layout(layout_key="roy_UMAP", color_by="roy_donor", shuffle_df=True)

# ---
#
# And now, do our own hvg identification, umap etc.

ds.mark_hvgs(max_mean=1, max_var=6)

ds.make_graph(feat_key="hvgs")

ds.run_umap(parallel=True)

ds.plot_cells_dists(group_key="roy_donor")

ds.plot_cells_dists(group_key="roy_tissue_type")

ds.plot_layout(layout_key="RNA_UMAP")

# Using their provided metadata:

ds.plot_layout(layout_key="RNA_UMAP", color_by="roy_cluster_id")

ds.plot_layout(layout_key="RNA_UMAP", color_by="roy_tissue_type")

ds.plot_layout(layout_key="RNA_UMAP", color_by="roy_cell_type")

ds.plot_layout(layout_key="RNA_UMAP", color_by="roy_donor")

ds.plot_layout(layout_key="RNA_UMAP", color_by="roy_donor", shuffle_df=True)

# ### Do own clustering etc

ds

ds.run_leiden_clustering()

ds.plot_layout(layout_key="RNA_UMAP", color_by="RNA_leiden_cluster")

# Look at our clustering on their UMAP:

ds.plot_layout(layout_key="roy_UMAP", color_by="RNA_leiden_cluster")

ds.plot_layout(
    layout_key="roy_UMAP", color_by="RNA_leiden_cluster", shuffle_df=True
)

# Looks pretty smoshed, might be that their UMAP has more of a 3D structure.

# ---
#
# Okay, all done for me. Good luck!
#
# /Rasmus
