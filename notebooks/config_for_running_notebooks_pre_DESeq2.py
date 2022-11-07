# execute all notebooks with this command:
# `jupyter nbconvert --inplace --execute --config config_for_running_notebooks_pre_DESeq2.py`
# note: if only the .py files exist, run this command first to generate the .ipynb files:
# `jupytext --sync *.py`

c = get_config()
c.NbConvertApp.notebooks = [
    "1.1_make_metadata_FL.ipynb",
    "1.2_make_ashapes.ipynb",
    "2.1_load_CITEseq.ipynb",
    "2.2_inserting_seurat_metadata_FL_ybm.ipynb",
    "2.3_load_CITEseq_FL_W9.ipynb",
    "3.1_mapping_FL_yBM.ipynb",
    "3.2_target_classification_FL_yBM.ipynb",
    "3.3_vizualizing_mapping_on_FL.ipynb",
    "3.4_vizualizing_mapping_on_BM.ipynb",
    "4.1_gene_enrichment_in_clusters.ipynb",
    "4.2_DE_FL_yBM_on_FL.ipynb",
    "4.3_DE_gene_analysis_cluster.ipynb",
    "4.4_bar_plot.ipynb",
    "5.1_make_adt_fcs.ipynb",
    "5.2_ADT_umaps.ipynb",
    "5.3_ADT_gate_analysis.ipynb",
    "5.4_ADT_gate_pseudo_bulk.ipynb",
]
c.NbConvertApp.export_format = "notebook"
