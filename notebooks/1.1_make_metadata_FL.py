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

import pandas as pd
import fl_utils

fl_utils.set_project_path()

in_dir = "data/raw/paper_specific/from_seurat"
out_dir = "data/processed/notebooks"

fl_utils.create_dir(out_dir)

# N. B. we use metadata and umap coordinates here to reproduce the results in
# the paper. If you would like to use the output from the seurat part instead,
# use files from the `data/processed/seurat/CS22` dir. And keep in mind that the
# results will probably be slightly different from ours due to randomness in the
# algorithms (e. g. cluster numbers might be different, resulting in wrong
# automatic naming of clusters).
metadata = pd.read_csv(in_dir + "/FL_combined_metadata.csv")
umap = pd.read_csv(in_dir + "/FL_combined_umap.csv")

metadata["u1"] = umap["UMAP_1"]
metadata["u2"] = umap["UMAP_2"] * (
    -1
)  # multipled by -1 to flip umap so HSCs are on top
metadata.index = [x.split(" ")[0] for x in metadata.index]
metadata.index = [
    x.split("_")[2] + "-" + x[2] + "_FL_hpc" for x in metadata.index
]

cluster_names = {
    1: "MPP-II",
    2: "Ly-I",
    3: "MEP",
    4: "GMP",
    5: "HSC",
    6: "Ly-II",
    7: "DC-Mono",
    8: "Cyc",
    9: "DC-I",
    10: "Ly-III",
    0: "MPP-I",
}

metadata["cluster"] = metadata["integrated_snn_res.0.7"]
metadata["cluster_names"] = [cluster_names[x] for x in metadata.cluster]

# + tags=[]
metadata = metadata.drop(
    columns=[
        "RNA_snn_res.0.8",
        "seurat_clusters",
        "RNA_snn_res.0.7",
        "orig_seurat_clusters",
        "clusters_annotated",
        "old.ident",
        "integrated_snn_res.0.7",
        "clust.names",
    ]
)
# -

colors = {
    0: "#ffb4b4",
    1: "#addfff",
    2: "#d18ce2",
    3: "#ff5a5a",
    4: "#6469ff",
    5: "#fff148",
    6: "#dd00ff",
    7: "#5dcbd6",
    8: "#efefd9",
    9: "#6af8c5",
    10: "#a96bff",
}
metadata["clust_color"] = [colors[x] for x in metadata.cluster]

metadata.to_csv(out_dir + "/FL_seurat_metadata.csv")
