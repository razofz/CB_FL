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

import scarf
import pandas as pd
import glob
import fl_utils

# -

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
seurat_dir = raw_dir + "/paper_specific/from_seurat"
gated_dir = raw_dir + "/paper_specific/ADT_gated_populations/conventional_gates"
out_dir = "data/processed/notebooks"
zarr_dir = out_dir + "/zarr_files"
external_dir = "data/external/osf-vdf42/cell_metadata"

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_split_yBM": {i: zarr_dir + f"/RNA_yBM{i}_hpc" for i in [1, 2]},
    "zarr_split_FL": {i: zarr_dir + f"/RNA_FL{i}_hpc" for i in [1, 2]},
    "metadata_FL": out_dir + "/FL_seurat_metadata.csv",
    "metadata_yBM": external_dir + "/yBM_hpc_cite_seq_celldata_seurat.csv",
    "hvgs_FL": seurat_dir + "/FL_combined_hvgs.csv",
    "hvgs_yBM": external_dir + "/yBM_hpc_seurat_hvgs.csv",
}

out_files = {
    "zarr_merged": zarr_dir + "/RNA_yBM_and_FL_hpc_merged",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# +
# fl_utils.create_dir(out_adt_matrices_dir)
# fl_utils.create_dir(out_adt_matrices_dir + "/FCS_files")
# -

# ---

ds1 = scarf.DataStore(
    in_files["zarr_split_yBM"][1],
    default_assay="RNA",
    assay_types={"HTO": "ADT"},
)
ds2 = scarf.DataStore(
    in_files["zarr_split_yBM"][2],
    default_assay="RNA",
    assay_types={"HTO": "ADT"},
)
ds3 = scarf.DataStore(
    in_files["zarr_split_FL"][1],
    default_assay="RNA",
    assay_types={"HTO": "ADT"},
)
ds4 = scarf.DataStore(
    in_files["zarr_split_FL"][2],
    default_assay="RNA",
    assay_types={"HTO": "ADT"},
)
_ = scarf.ZarrMerge(
    out_files["zarr_merged"],
    [ds1.RNA, ds2.RNA, ds3.RNA, ds4.RNA],
    ["ybm1", "ybm2", "fl1", "fl2"],
    "RNA",
    overwrite=True,
).dump(nthreads=threads)

del ds1
del ds2
del ds3
del ds4

# +
# ds1 = scarf.DataStore('_zarr_files/RNA_yBM1_hpc', default_assay='RNA', assay_types={'HTO': 'ADT'})
# ds2 = scarf.DataStore('_zarr_files/RNA_yBM2_hpc', default_assay='RNA', assay_types={'HTO': 'ADT'})
# ds3 = scarf.DataStore('_zarr_files/RNA_FL1_hpc', default_assay='RNA', assay_types={'HTO': 'ADT'})
# ds4 = scarf.DataStore('_zarr_files/RNA_FL2_hpc/', default_assay='RNA', assay_types={'HTO': 'ADT'})
# _ = scarf.ZarrMerge('_zarr_files/RNA_yBM_and_FL_hpc_merged', [ds1.RNA, ds2.RNA, ds3.RNA, ds4.RNA], ['ybm1', 'ybm2', 'fl1', 'fl2'], 'RNA').write()
# -

ds = scarf.DataStore(
    out_files["zarr_merged"], default_assay="RNA", nthreads=threads
)


# +
# ds = scarf.DataStore('_zarr_files/RNA_yBM_and_FL_hpc_merged', default_assay='RNA')

# +
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


# -

seurat_md_fn_BM = in_files["metadata_yBM"]
# seurat_md_fn_BM = '../data/cell_metadata/yBM_hpc_cite_seq_celldata_seurat.csv'
seurat_filtered_cells_BM = fix_index_BM(
    pd.read_csv(seurat_md_fn_BM, index_col=0).index
)
len(seurat_filtered_cells_BM)

seurat_md_fn_FL = in_files["metadata_FL"]
# seurat_md_fn_FL = './from_seurat/FL_seurat_metadata.csv'
seurat_filtered_cells_FL = fix_index_FL(
    pd.read_csv(seurat_md_fn_FL, index_col=0).index
)
len(seurat_filtered_cells_FL)

# + tags=[]
seurat_filtered_cells = []
seurat_filtered_cells = seurat_filtered_cells_BM + seurat_filtered_cells_FL
# -

ds.cells.update_key(
    ds.cells.index_to_bool(
        ds.cells.get_index_by(seurat_filtered_cells, column="ids")
    ),
    key="I",
)
ds.cells.fetch_all("I").sum()

ds.cells.insert(
    "sample",
    [x.split("__")[0] for x in ds.cells.fetch_all("ids")],
    overwrite=True,
)
for sample in ["ybm1", "ybm2", "fl1", "fl2"]:
    v = ds.cells.fetch_all("I") & (ds.cells.fetch_all("sample") == sample)
    ds.cells.insert(sample, v, fill_value=False, overwrite=True)

pd.unique([x.split("__")[0] for x in ds.cells.fetch_all("ids")])

hvg_fl = pd.read_csv(in_files["hvgs_FL"])["x"].values
hvg_bm = pd.read_csv(in_files["hvgs_yBM"])["x"].values
# hvg_fl = pd.read_csv('./from_seurat/FL_combined_hvgs.csv')['x'].values
# hvg_bm = pd.read_csv('../data/cell_metadata/yBM_hpc_seurat_hvgs.csv')['x'].values
print(len(hvg_fl), len(hvg_bm))

# + tags=[]
hvg_FL_bm = set(hvg_fl).intersection(hvg_bm)
len(set(hvg_fl).intersection(hvg_bm))
# -

ds_FL = scarf.DataStore(in_files["zarr_FL"], nthreads=threads)
ds_BM = scarf.DataStore(in_files["zarr_yBM"], nthreads=threads)
# ds_FL = scf.DataStore('./_zarr_files/RNA_FL_hpc_merged', nthreads=2)
# ds_BM = scf.DataStore('./_zarr_files/RNA_yBM_hpc_merged', nthreads=2)

# +
df_BM = ds_BM.cells.to_pandas_dataframe(["ids", "predicted_celltype_FL"])
df_BM.index = df_BM["ids"]
df_BM["predicted_celltype_FL"] = [
    "bm_" + x for x in df_BM["predicted_celltype_FL"]
]
df_BM.columns = ["ids", "clusters"]

df_FL = ds_FL.cells.to_pandas_dataframe(["ids", "RNA_cluster_names"])
df_FL.index = df_FL["ids"]
df_FL["RNA_cluster_names"] = ["fl_" + x for x in df_FL["RNA_cluster_names"]]
df_FL.columns = ["ids", "clusters"]

cell_clusters = pd.DataFrame()
cell_clusters = df_BM
cell_clusters = cell_clusters.append(df_FL)

cell_clusters = cell_clusters["clusters"].reindex(ds.cells.fetch("ids")).values
ds.cells.insert("clusters", cell_clusters, overwrite=True, fill_value="NA")
# -

ds.RNA.feats.insert(
    "I__hvgs",
    ds.RNA.feats.index_to_bool(
        ds.RNA.feats.get_index_by(hvg_FL_bm, column="ids")
    ),
    fill_value=False,
    overwrite=True,
)

# + tags=[]
ds.make_graph(
    feat_key="hvgs", k=21, dims=30, n_centroids=100, renormalize_subset=True
)

# + tags=[]
ds.run_umap(n_epochs=200, spread=5, min_dist=2, parallel=True)
# -

ds.plot_layout(layout_key="RNA_UMAP", color_by="sample")

sample_colors = {
    "ybm1": "#ff7f0e",
    "ybm2": "#ff7f0e",
    "fl1": "#1f77b4",
    "fl2": "#1f77b4",
}
ds.cells.insert(
    "sample_colors",
    [sample_colors[x.split("__")[0]] for x in ds.cells.fetch("ids")],
    overwrite=True,
)
