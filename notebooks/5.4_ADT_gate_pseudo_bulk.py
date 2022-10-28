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
import json

import numpy as np
import pandas as pd
import scarf
import fl_utils

# -

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
out_dir = "data/processed/notebooks"
deseq2_dir = out_dir + "/DEseq2"
deseq2_gatewise_dir = deseq2_dir + "/gate_wise"
zarr_dir = out_dir + "/zarr_files"
# -

gates = [
    "CD49Fpos_HSC",
    "CD90pos_HSC",
    "CD38neg_MPP",
    "LMPP",
    "CMP_123",
    "CMP_135",
    "GMP_123",
    "GMP_135",
    "MEP_123",
    "MEP_135",
    "CD10_CLP",
    "CD7_CLP",
    "IL7RA_CLP",
]

# +
in_files = {
    "zarr_yBM": zarr_dir + "/RNA_yBM_hpc_merged",
    "zarr_FL": zarr_dir + "/RNA_FL_hpc_merged",
    "gated": out_dir + "/RNA_identified_gate_cell.json",
}

out_files = {
    "gatewise_cts": {
        gate: deseq2_gatewise_dir + f"/{gate}_cts_all_hpc.csv" for gate in gates
    },
    "gatewise_coldata": {
        gate: deseq2_gatewise_dir + f"/{gate}_coldata_all_hpc.csv"
        for gate in gates
    },
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

fl_utils.create_dir(deseq2_gatewise_dir)

# ---

ds_yBM = scarf.DataStore(in_files["zarr_yBM"], default_assay="RNA")
ds_FL = scarf.DataStore(in_files["zarr_FL"], default_assay="RNA")

gated_cells = json.load(open(in_files["gated"]))
gated_cells.keys()

# +
gates = [
    "CD49Fpos_HSC",
    "CD90pos_HSC",
    "CD38neg_MPP",
    "LMPP",
    "CMP_123",
    "CMP_135",
    "GMP_123",
    "GMP_135",
    "MEP_123",
    "MEP_135",
    "CD10_CLP",
    "CD7_CLP",
    "IL7RA_CLP",
]

for n, gate in enumerate(gates):
    cells = gated_cells[gate]["FL_hpc"]
    ds_FL.cells.insert(
        f"gated_{gate}",
        [i.split("__")[0] + gate if i in cells else False for i in
         ds_FL.cells.fetch_all("ids")],
        overwrite=True,
    )

x = np.char.add(
    ds_FL.cells.fetch_all("sample"),
    ["_" + x for x in ds_FL.cells.fetch_all("RNA_cluster_names")],
)

# + tags=[]

for n, gate in enumerate(gates):
    cells = gated_cells[gate]["yBM_hpc"]
    ds_yBM.cells.insert(
        f"gated_{gate}",
        [i.split("__")[0] + gate if i in cells else False for i in
         ds_yBM.cells.fetch_all("ids")],
        overwrite=True,
    )
# -

dfs = []

# + tags=[]
for gate in gates:
    df = ds_FL.make_bulk(
        from_assay="RNA", group_key=f"gated_{gate}", pseudo_reps=2
    )
    df = df.drop(columns=["False_Rep1", "False_Rep2", "names"])
    df.columns = [f"FL_{x[:3]}_{gate}_" + x.split("_")[-1] for x in df.columns]
    dfs.append(df.T)

# + tags=[]
for gate in gates:
    df = ds_yBM.make_bulk(
        from_assay="RNA", group_key=f"gated_{gate}", pseudo_reps=2
    )
    df = df.drop(columns=["False_Rep1", "False_Rep2", "names"])
    df.columns = [f"yBM_{x[:4]}_{gate}_" + x.split("_")[-1] for x in df.columns]
    dfs.append(df.T)
# -

dfs = pd.concat(dfs).fillna(0).T
dfs.columns = [x.replace("-", ".") for x in dfs.columns]

metadata = pd.DataFrame(
    [
        (
            x,
            x.split("_", 1)[0],
            x.rsplit("_", 1)[0].split("_", 2)[2],
            x.rsplit("_", 1)[1],
        )
        for x in dfs.columns
    ],
    columns=["sample_id", "sample", "gate_name", "pseudo_rep"],
).set_index("sample_id")

for i in metadata.gate_name.unique():
    print(i)
    idx = metadata.gate_name == i
    x = dfs[metadata[idx].index]
    x.to_csv(out_files["gatewise_cts"][i])
    metadata[idx].to_csv(out_files["gatewise_coldata"][i])
