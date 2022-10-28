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
external_dir = "data/external/"
zarr_dir = "data/processed/notebooks/zarr_files"
# out_dir = 'data/processed/notebooks/zarr_files'

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "zarr_CS16": zarr_dir + "/RNA_FL_CS16_merged",
    "zarr_W9": zarr_dir + "/RNA_FL_W9_merged",
    "zarr_CB": zarr_dir + "/RNA_CB_hpc_merged",
    "hvgs_yBM": external_dir
    + "/osf-vdf42/cell_metadata/yBM_hpc_seurat_hvgs.csv",
    "hvgs_FL": raw_dir + "/FL_combined_hvgs.csv",
}

out_files = {}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")
ds_CB = scarf.DataStore(in_files["zarr_CB"], default_assay="RNA")
ds_FL9 = scarf.DataStore(in_files["zarr_W9"], default_assay="RNA")
ds_FL16 = scarf.DataStore(in_files["zarr_CS16"], default_assay="RNA")

ds_yBM.RNA.feats.insert(
    "I__hvgs",
    ds_yBM.RNA.feats.index_to_bool(
        ds_yBM.RNA.feats.get_index_by(
            pd.read_csv(in_files["hvgs_yBM"])["x"].values,
            column="ids",
        )
    ),
    fill_value=False,
    overwrite=True,
)

len(
    ds_FL.cells.fetch("RNA_cluster_names")[
        ds_FL.cells.fetch("RNA_cluster_names") == "MEP"
    ]
)

ds_FL.RNA.feats.insert(
    "I__hvgs",
    ds_FL.RNA.feats.index_to_bool(
        ds_FL.RNA.feats.get_index_by(
            pd.read_csv(in_files["hvgs_FL"])["x"].values, column="ids"
        )
    ),
    fill_value=False,
    overwrite=True,
)

print(
    ds_FL.RNA.feats.fetch("I__hvgs").sum(),
    ds_yBM.RNA.feats.fetch("I__hvgs").sum(),
)

ds_yBM.make_graph(feat_key="hvgs", k=9, dims=15, n_centroids=100)
ds_FL.make_graph(feat_key="hvgs", k=9, dims=15, n_centroids=100)

ds_yBM.run_mapping(
    target_assay=ds_FL.RNA,
    target_name="FL_HPC",
    target_feat_key="hvgs_ybmhpc",
    save_k=9,
)
ds_yBM.run_mapping(
    target_assay=ds_CB.RNA,
    target_name="CB_HPC",
    target_feat_key="hvgs_ybmhpc",
    save_k=9,
)
ds_yBM.run_mapping(
    target_assay=ds_yBM.RNA,
    target_name="yBM_HPC",
    target_feat_key="hvgs_ybmhpc",
    save_k=9,
)
ds_yBM.run_mapping(
    target_assay=ds_FL9.RNA,
    target_name="FL_W9",
    target_feat_key="hvgs_ybmhpc",
    save_k=9,
)
ds_yBM.run_mapping(
    target_assay=ds_FL16.RNA,
    target_name="FL_CS16",
    target_feat_key="hvgs_ybmhpc",
    save_k=9,
)

ds_FL.run_mapping(
    target_assay=ds_yBM.RNA,
    target_name="yBM_HPC",
    target_feat_key="hvgs_flhpc",
    save_k=9,
)
ds_FL.run_mapping(
    target_assay=ds_CB.RNA,
    target_name="CB_HPC",
    target_feat_key="hvgs_flhpc",
    save_k=9,
)
ds_FL.run_mapping(
    target_assay=ds_FL.RNA,
    target_name="FL_HPC",
    target_feat_key="hvgs_flhpc",
    save_k=9,
)
ds_FL.run_mapping(
    target_assay=ds_FL9.RNA,
    target_name="FL_W9",
    target_feat_key="hvgs_flhpc",
    save_k=9,
)
ds_FL.run_mapping(
    target_assay=ds_FL16.RNA,
    target_name="FL_CS16",
    target_feat_key="hvgs_flhpc",
    save_k=9,
)
