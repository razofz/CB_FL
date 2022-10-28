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
import glob
import json

import fcswrite
import fl_utils
import numpy as np
import pandas as pd

# -

fl_utils.set_project_path()

# +
threads = 16

raw_dir = "data/raw"
adt_matrices_dir = raw_dir + "/paper_specific/from_seurat/adt_matrices"
out_dir = "data/processed/notebooks"
out_adt_matrices_dir = out_dir + "/adt_matrices"

# +
in_files = {
    "adt_matrices": glob.glob(adt_matrices_dir + "/*"),
}

out_files = {
    "csvs": [
        f"{out_adt_matrices_dir}/FCS_files/adt_matrices_{sample}.csv"
        for sample in ["yBM_hpc", "FL_hpc"]
    ],
    "fcss": [
        f"{out_adt_matrices_dir}/FCS_files/adt_matrices_{sample}.fcs"
        for sample in ["yBM_hpc", "FL_hpc"]
    ],
    "seurat_matrices": [
        f"{out_adt_matrices_dir}/seurat_adt_matrices_{sample}.csv"
        for sample in ["yBM_hpc", "FL_hpc"]
    ],
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")
# -

fl_utils.create_dir(out_adt_matrices_dir)
fl_utils.create_dir(out_adt_matrices_dir + "/FCS_files")


# ---


def read_csv_mod_idx(fn, suf):
    df = pd.read_csv(fn).T
    if df.index[0].find("-") == -1:
        suf = "-1" + suf
    df.index = [x + suf for x in df.index]
    return df


def load_rep_adt(name):
    sample, pop = name.split("_")
    fn = f"{adt_matrices_dir}/{sample}%s_{pop}_adt_scaled_values.csv"
    in_files["adt_matrices"]
    # fn = f'./seurat_normalized_adt_matrices/{sample}%s_{pop}_adt_scaled_values.csv'
    return read_csv_mod_idx(fn % "1", f"_{sample}1_{pop}").append(
        read_csv_mod_idx(fn % "2", f"_{sample}2_{pop}"), sort=False
    )


# + tags=[]
seurat_adts = {
    "yBM_hpc": load_rep_adt("yBM_hpc"),
    "FL_hpc": load_rep_adt("FL_hpc"),
}
seurat_adts

# + tags=[]
names = []
for i in list(seurat_adts["FL_hpc"].index):
    x = i.split("_", 2)[2]
    names.append(x)
seurat_adts["FL_hpc"].index = names
# -

for i in seurat_adts:
    df = seurat_adts[i]
    df.to_csv(f"{out_adt_matrices_dir}/seurat_adt_matrices_{i}.csv")
    # df.to_csv(f"./seurat_normalized_adt_matrices/ADT_scaled_scarf_naming/seurat_adt_matrices_{i}.csv")

for i in seurat_adts:
    df = np.e ** seurat_adts[i] * 1e3
    df.to_csv(f"{out_adt_matrices_dir}/FCS_files/adt_matrices_{i}.csv")
    fcswrite.write_fcs(
        f"{out_adt_matrices_dir}/FCS_files/adt_matrices_{i}.fcs",
        list(df.columns),
        df.values,
    )
    # df.to_csv(f"./seurat_normalized_adt_matrices/FCS_files/adt_matrices_{i}.csv")
    # fcswrite.write_fcs(f"./seurat_normalized_adt_matrices/FCS_files/adt_matrices_{i}.fcs", list(df.columns), df.values)
