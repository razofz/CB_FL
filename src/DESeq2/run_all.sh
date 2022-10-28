which Rscript
export file=1.1_All_clustered_pops.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=1.2_All_gated_pops.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=1.3_produce_FL_core.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=2.0_runner.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=3.1_investigate_FL_BM_clusters_in_PCA.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=3.2_investigate_FL_BM_gates_in_PCA.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=4.0_export_excel_sheet_to_tsv.py &&
echo "> Executing $file:" &&
python $file &&
export file=4.1_fetal_signature_in_leukemia.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=4.2_fetal_signature_in_leukemia_ind_PCA_fetalcore.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=4.3_fetal_signature_in_leukemia_randomgenes.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=5.1_dev_cell_heatmap.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=6.1_rev_stat_test_distributions_PC.R &&
echo "> Executing $file:" &&
Rscript $file &&
export file=6.2_rev_fetal_core_boxplots_final.R &&
echo "> Executing $file:" &&
Rscript $file 
