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

import scarf
import glob
from IPython.display import display
import fl_utils

# -

fl_utils.set_project_path()

in_dir = "data/raw/cite_cell_matrices"
out_dir = "data/processed/notebooks/zarr_files"

# +
in_files = {
    "matrices": glob.glob(in_dir + "/*"),
}

out_files = {
    "zarr_files": [
        f"{out_dir}/RNA_{fn.split('/')[-1]}" for fn in in_files["matrices"][1:]
    ],
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

assay_renamer = {
    "Gene Expression": "RNA",
    "Antibody capture": "ADT",
    "Antibody capture 2": "HTO",
}
for fn in in_files["matrices"]:
    if fn.split("/")[-1] == "FL_CS16_and_W9":
        continue
    print(fn)
    reader = scarf.CrDirReader(fn, "rna", mtx_separator="\t")
    assay_names = {
        k: assay_renamer[v] for k, v in reader.assayFeats.T["type"].items()
    }
    reader.rename_assays(assay_names)
    display(reader.assayFeats)
    _ = scarf.CrToZarr(
        reader,
        zarr_fn=f"{out_dir}/RNA_{fn.split('/')[-1]}",
        chunk_size=(1000, 1000),
    ).dump()
    print()

samples_list = {
    "yBM": 2,
    "FL": 2,
    "CB": 3,
}

for sample in samples_list:
    print(sample)
    datastores = [
        scarf.DataStore(
            out_dir + f"/RNA_{sample}{i}_hpc",
            default_assay="RNA",
            assay_types={"HTO": "ADT"},
        )
        for i in range(1, 1 + samples_list[sample])
    ]
    _ = scarf.ZarrMerge(
        out_dir + f"/RNA_{sample}_hpc_merged",
        [ds.RNA for ds in datastores],
        [sample.lower() + str(i) for i in range(1, 1 + samples_list[sample])],
        "RNA",
        overwrite=True,
    ).dump()
    _ = scarf.ZarrMerge(
        out_dir + f"/RNA_{sample}_hpc_merged",
        [ds.ADT for ds in datastores],
        [sample.lower() + str(i) for i in range(1, 1 + samples_list[sample])],
        "ADT",
        overwrite=True,
    ).dump()
