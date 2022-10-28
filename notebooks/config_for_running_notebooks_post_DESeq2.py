# execute all notebooks with this command:
# `jupyter nbconvert --inplace --execute --config config_for_running_notebooks_post_DESeq2.py`
# note: if only the .py files exist, run this command first to generate the .ipynb files:
# `jupytext --sync *.py`

c = get_config()
c.NbConvertApp.notebooks = [
    "6.1_combined_FL_yBM.ipynb",
    "6.2_heatmap_diff_genes_ybm_fl.ipynb",
    "6.3_heatmap_diff_genes_primitive_clusts_fl.ipynb",
    "7.0_Roy_data_load.ipynb",
    "7.1_projection_roy_data.ipynb",
    "7.2_vizualizing_roy_data.ipynb",
    "8.1_hox_genes_FL_BM.ipynb",
    "9.1_sample_contributions.ipynb",
    "9.2_WNN_integration_viz.ipynb",
    "9.3_harmony_integration_viz.ipynb",
    "9.4_z_imported_Civij_data-OSF.ipynb",
    "9.5_preds_barplots.ipynb"
]
c.NbConvertApp.export_format = "notebook"
