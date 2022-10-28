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
import numpy as np
import os
import shutil
import fl_utils

fl_utils.set_project_path()

raw_dir = "data/raw"
zarr_dir = "data/processed/notebooks/zarr_files"

# +
in_files = {
    "matrices": raw_dir + "/cite_cell_matrices/FL_CS16_and_W9",
    "metadata": raw_dir + "/paper_specific/from_seurat/FL_W9_metadata.csv",
}

out_files = {
    "zarr_CS16_and_W9": zarr_dir + "/RNA_FL_CS16_and_W9",
    "zarr_W9": zarr_dir + "/RNA_FL_W9_merged",
    "zarr_CS16": zarr_dir + "/RNA_FL_CS16_merged",
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

assay_renamer = {
    "Gene Expression": "RNA",
    "Antibody capture": "ADT",
    "Antibody capture 2": "HTO",
}
fn = in_files["matrices"]
reader = scarf.CrDirReader(fn, "rna")
assay_names = {
    k: assay_renamer[v] for k, v in reader.assayFeats.T["type"].items()
}
reader.rename_assays(assay_names)
_ = scarf.CrToZarr(
    reader, zarr_fn=out_files["zarr_CS16_and_W9"], chunk_size=(1000, 1000)
).dump()

ds = scarf.DataStore(out_files["zarr_CS16_and_W9"], default_assay="RNA")


def fix_index_FL(indexv):
    ret_val = []
    for i in indexv:
        x = i + "-1"
        ret_val.append(f"{x}")
    return ret_val


seurat_md_fn_FL = in_files["metadata"]
seurat_filtered_cells_FL = fix_index_FL(
    pd.read_csv(seurat_md_fn_FL, index_col=0).index
)
len(seurat_filtered_cells_FL)

ds.cells.update_key(
    ds.cells.index_to_bool(
        ds.cells.get_index_by(seurat_filtered_cells_FL, column="ids")
    ),
    key="I",
)
ds.cells.fetch_all("I").sum()

names = ds.RNA.feats.fetch_all("names")
name_counts = pd.Series(names).value_counts()

new_ids = []
for i in names:
    if name_counts[i] > 0:
        if name_counts[i] == 1:
            new_ids.append(i)
        else:
            new_ids.append(f"{i}_{name_counts[i]-1}")
        name_counts[i] = name_counts[i] - 1

len(names), len(name_counts), len(new_ids)

ds.RNA.feats.insert("new_ids", new_ids, overwrite=True)

shutil.rmtree(out_files["zarr_CS16_and_W9"] + "/RNA/featureData/ids")
os.rename(
    src=out_files["zarr_CS16_and_W9"] + "/RNA/featureData/new_ids",
    dst=out_files["zarr_CS16_and_W9"] + "/RNA/featureData/ids",
)


def import_seurat_column_FL(fn, colname, save_key, fill_val):
    col = pd.read_csv(fn, index_col=0)[colname]
    col.index = fix_index_FL(col.index)
    col = col.reindex(ds.cells.fetch("ids")).values
    ds.cells.insert(save_key, col, overwrite=True, fill_value=fill_val)
    return None


import_seurat_column_FL(seurat_md_fn_FL, "hash.ID", "hash_ID", "NA")


y = np.array([x for x in ds.cells.fetch("hash_ID")])

np.unique(y)

for sample in ["HTO-FL1", "HTO-FL2", "HTO-FL3", "HTO-FL4"]:
    v = ds.cells.fetch_all("I") & (ds.cells.fetch_all("hash_ID") == sample)
    ds.cells.insert(sample, v, fill_value=False, overwrite=True)

# zarr_dir = os.path.dirname(out_files["zarr_CS16_and_W9"])
print(zarr_dir)

for subsample in list(range(1, 4 + 1)):
    scarf.SubsetZarr(
        in_zarr=out_files["zarr_CS16_and_W9"],
        out_zarr=f"{zarr_dir}/RNA_FL{subsample}_"
        + ("CS16" if subsample < 3 else "W9"),
        cell_key=f"HTO-FL{subsample}",
    ).dump()

ds1 = scarf.DataStore(f"{zarr_dir}/RNA_FL3_W9/", default_assay="RNA")
ds2 = scarf.DataStore(f"{zarr_dir}/RNA_FL4_W9/", default_assay="RNA")
_ = scarf.ZarrMerge(
    out_files["zarr_W9"],
    [ds1.RNA, ds2.RNA],
    ["FL3", "FL4"],
    "RNA",
    overwrite=True,
).dump()

ds1 = scarf.DataStore(f"{zarr_dir}/RNA_FL1_CS16/", default_assay="RNA")
ds2 = scarf.DataStore(f"{zarr_dir}/RNA_FL2_CS16/", default_assay="RNA")
_ = scarf.ZarrMerge(
    out_files["zarr_CS16"],
    [ds1.RNA, ds2.RNA],
    ["FL1", "FL2"],
    "RNA",
    overwrite=True,
).dump()
