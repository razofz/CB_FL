configfile: "config.yaml"


DESEQ2_DIR = config["processed_dir"] + "DESeq2/results/"
DESEQ2_PLOT_DIR = config["processed_dir"] + "DESeq2/images/"


rule all:
    input:
        expand(
            DESEQ2_DIR
            + "Normed_counts/{fl_core_version}/{cluster}_DESeq_normed.txt",
            cluster=config["fl_clusters_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.csv",
            cluster=config["fl_clusters_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR
            + "Normed_counts/{fl_core_version}/gated_{gate}_DESeq_normed.txt",
            gate=config["gates_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/gated_{gate}_FL_yBM.csv",
            gate=config["gates_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/Fetal_signature_p005.txt",
            fl_core_version=config["fl_cores"],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/Adult_signature_p005.txt",
            fl_core_version=config["fl_cores"],
        ),
        DESEQ2_PLOT_DIR + "PC12_singel_cell_clusters_top500.pdf",
        DESEQ2_PLOT_DIR + "PC12_FL_BM_singel_cell_gates_top500.pdf",
        DESEQ2_PLOT_DIR + "PC1Age_All_leuk_FL_core.pdf",
        DESEQ2_PLOT_DIR + "IndividualPC1_PC2_fetalcore_MLLAF4.pdf",
        DESEQ2_PLOT_DIR + "Main_Figure_5a.pdf",


rule deseq_all_clustered_pops:
    input:
        cts_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_cts_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        coldata_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_coldata_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
    output:
        normed_counts=expand(
            DESEQ2_DIR
            + "Normed_counts/{fl_core_version}/{cluster}_DESeq_normed.txt",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        deseq_results=expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        fl_only_normed_counts=expand(
            DESEQ2_DIR
            + "Normed_counts/{fl_core_version}/only_FL_DESeq_normed.txt",
            allow_missing=True,
        ),
        gsea_fc_rnk=expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.rnk",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        mpp1_hsc_rnk=expand(
            DESEQ2_DIR + "{fl_core_version}/FL_MPPI_HSC.rnk",
            allow_missing=True,
        ),
        mpp1_mpp2_rnk=expand(
            DESEQ2_DIR + "{fl_core_version}/FL_MPPI_MPPII.rnk",
            allow_missing=True,
        ),
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.1_All_clustered_pops.R"


rule deseq_all_gated_pops:
    input:
        cts_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_cts_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
        coldata_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_coldata_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
    output:
        normed_counts=expand(
            DESEQ2_DIR
            + "Normed_counts/{fl_core_version}/gated_{gate}_DESeq_normed.txt",
            gate=config["gates_to_use"],
            allow_missing=True,
        ),
        deseq_results=expand(
            DESEQ2_DIR + "{fl_core_version}/gated_{gate}_FL_yBM.csv",
            gate=config["gates_to_use"],
            allow_missing=True,
        ),
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.2_All_gated_pops.R"


# rule gather_cores:
#     input:
#         fetal_signature=expand(
#             DESEQ2_DIR + "{fl_core_version}/Fetal_signature_p005.txt",
#             fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
#         ),
#         adult_signature=expand(
#             DESEQ2_DIR + "{fl_core_version}/Adult_signature_p005.txt",
#             fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
#         ),


rule deseq_produce_fl_core:
    input:
        deseq_results=expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
    output:
        sub_setter_pos=expand(
            DESEQ2_DIR + "sub_setter/{fl_core_version}/{cluster}_pos.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        sub_setter_neg=expand(
            DESEQ2_DIR + "sub_setter/{fl_core_version}/{cluster}_neg.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        de_genes=DESEQ2_DIR + "{fl_core_version}/DE_genes.csv",
        adult_signature=DESEQ2_DIR
        + "{fl_core_version}/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR
        + "{fl_core_version}/Fetal_signature_p005.txt",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.3a_produce_FL_core.R"


# not properly utilising smk features, due to time constraints
rule single_cell_deg:
    input:
        cite_cell_matrices_dir=config["raw_dir"] + "cite_cell_matrices/",
        combined_metadata=config["raw_dir"]
        + "paper_specific/from_seurat/FL_BM_combined_metadata.csv",
        mapping_predictions=config["processed_dir"]
        + "notebooks/mapping_predictions/RNA_FL_target_preds.json",
        fl_metadata=config["processed_dir"]
        + "seurat/CS22/FL_combined_metadata.csv",
    output:
        cells=expand(
            config["interim_dir"] + "single_cell_deg/{sample}_cells.csv",
            sample=config["samples"],
        ),
        deg_results_bm=expand(
            DESEQ2_DIR + "FLcoreSC/{cluster}_BM_specific_markers.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        deg_results_fl=expand(
            DESEQ2_DIR + "FLcoreSC/{cluster}_FL_specific_markers.csv",
            cluster=config["fl_clusters_to_use"],
        ),
    conda:
        "envs/seurat.yaml"
    threads: 4
    script:
        "src/DESeq2/0.1_sc_DEG.R"


rule produce_sc_fl_core:
    input:
        deg_results_bm=rules.single_cell_deg.output.deg_results_bm,
        deg_results_fl=rules.single_cell_deg.output.deg_results_fl,
    output:
        sub_setter_pos=expand(
            DESEQ2_DIR + "sub_setter/FLcoreSC/{cluster}_pos.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        sub_setter_neg=expand(
            DESEQ2_DIR + "sub_setter/FLcoreSC/{cluster}_neg.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        de_genes=DESEQ2_DIR + "FLcoreSC/DE_genes.csv",
        adult_signature=DESEQ2_DIR + "FLcoreSC/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcoreSC/Fetal_signature_p005.txt",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.3b_produce_sc_FL_core.R"


rule plot_pca:
    input:
        cts_files_cluster=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_cts_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        coldata_files_cluster=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_coldata_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        cts_files_gate=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_cts_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
        coldata_files_gate=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_coldata_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
    output:
        plot_pc12_cluster=DESEQ2_PLOT_DIR
        + "PC12_singel_cell_clusters_top500.pdf",
        plot_pc23_cluster=DESEQ2_PLOT_DIR
        + "PC23_singel_cell_clusters_top500.pdf",
        plot_pc12_gate=DESEQ2_PLOT_DIR + "PC12_singel_cell_gates_top500.pdf",
        plot_pc23_gate=DESEQ2_PLOT_DIR + "PC23_singel_cell_gates_top500.pdf",
        plot_pc12_both=DESEQ2_PLOT_DIR
        + "PC12_singel_cell_clusters_gates_top500.pdf",
        plot_pc23_both=DESEQ2_PLOT_DIR
        + "PC23_singel_cell_clusters_gates_top500.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/2.1_plot_pca.R"


rule plot_fl_bm_clusters_pca:
    input:
        cts_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_cts_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        coldata_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/cluster_wise/"
            + "{cluster}_coldata_all_hpc.csv",
            cluster=config["fl_clusters_to_use"],
        ),
    output:
        plot_pc12=DESEQ2_PLOT_DIR + "PC12_FL_BM_singel_cell_clusters_top500.pdf",
        plot_pc23=DESEQ2_PLOT_DIR + "PC23_FL_BM_singel_cell_clusters_top500.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/3.1_investigate_FL_BM_clusters_in_PCA.R"


rule plot_fl_bm_gates_pca:
    input:
        cts_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_cts_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
        coldata_files=expand(
            config["raw_dir"]
            + "paper_specific/DESeq2/gate_wise/"
            + "{gate}_coldata_all_hpc.csv",
            gate=config["gates_to_use"],
        ),
    output:
        plot_pc12=DESEQ2_PLOT_DIR + "PC12_FL_BM_singel_cell_gates_top500.pdf",
        plot_pc23=DESEQ2_PLOT_DIR + "PC23_FL_BM_singel_cell_gates_top500.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/3.2_investigate_FL_BM_gates_in_PCA.R"


rule export_ball_excel_to_tsv:
    input:
        ball_xlsx=config["external_dir"]
        + "BALL-1988S-HTSeq/B-ALL-subtyping.xlsx",
    output:
        ball_tsv=config["external_dir"] + "BALL-1988S-HTSeq/subtypes.tsv",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/4.0_export_excel_sheet_to_tsv.py"


rule fetal_signature_in_leukemia:
    input:
        adult_signature=DESEQ2_DIR + "FLcorePseudotech/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcorePseudotech/Fetal_signature_p005.txt",
        ball_tsv=rules.export_ball_excel_to_tsv.output.ball_tsv,
        ball_dir=config["external_dir"] + "BALL-1988S-HTSeq/",
    output:
        plot_pc12_prebatch=DESEQ2_PLOT_DIR + "PC12_top500_preBatch_ALL_leuk.pdf",
        plot_pc12_postbatch=DESEQ2_PLOT_DIR
        + "PC12_top500_postBatch_ALL_leuk.pdf",
        plot_pc1_age=DESEQ2_PLOT_DIR + "PC1Age_All_leuk_FL_core.pdf",
        plot_pc2_age=DESEQ2_PLOT_DIR + "PC2Age_All_leuk_FL_core.pdf",
        down_inf_leukemia=DESEQ2_DIR + "down_inf_luek.csv",
        up_inf_leukemia=DESEQ2_DIR + "up_inf_luek.csv",
        plot_heatmap_fl_sig=DESEQ2_PLOT_DIR + "heatmap_A_fetal_sig_ALL_AF4.pdf",
        plot_pc12_af4_fl_core=DESEQ2_PLOT_DIR + "PC12_AF4_FL_core.pdf",
        plot_pc1age_af4_fl_core=DESEQ2_PLOT_DIR + "PC1Age_AF4_FL_core.pdf",
        down_inf_leukemia_af4=DESEQ2_DIR + "down_inf_luek_AF4.csv",
        up_inf_leukemia_af4=DESEQ2_DIR + "up_inf_luek_AF4.csv",
        plot_heatmap_fl_sig_af4=DESEQ2_PLOT_DIR
        + "heatmap_B_fetal_sig_ALL_AF4.pdf",
        plot_heatmap_hox_genes_af4=DESEQ2_PLOT_DIR + "heatmap_Hox_genes_AF4.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/4.1_fetal_signature_in_leukemia.R"


rule fetal_signature_in_leukemia_pca_fl_core:
    input:
        adult_signature=DESEQ2_DIR + "FLcorePseudotech/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcorePseudotech/Fetal_signature_p005.txt",
        ball_tsv=rules.export_ball_excel_to_tsv.output.ball_tsv,
        ball_dir=config["external_dir"] + "BALL-1988S-HTSeq/",
    output:
        plots_pc1=expand(
            DESEQ2_PLOT_DIR + "IndividualPC1_fetalcore_{leuk_type}.pdf",
            leuk_type=config["leukemia_types"],
        ),
        plots_pc1_gender=expand(
            DESEQ2_PLOT_DIR + "IndividualPC1_fetalcore_{leuk_type}_gender.pdf",
            leuk_type=config["leukemia_types"],
        ),
        plot_fl_core_same_pc1=DESEQ2_PLOT_DIR + "fetalcore_samePC1_all.pdf",
        plot_fl_core_same_pc1_gender=DESEQ2_PLOT_DIR
        + "fetalcore_samePC1_all_gender.pdf",
        plots_type=expand(
            DESEQ2_PLOT_DIR + "{leuk_type}.pdf",
            leuk_type=config["leukemia_types"],
        ),
        plots_type_gender=expand(
            DESEQ2_PLOT_DIR + "{leuk_type}_gender.pdf",
            leuk_type=config["leukemia_types"],
        ),
        plot_fl_core_mllaf4=DESEQ2_PLOT_DIR
        + "IndividualPC1_PC2_fetalcore_MLLAF4.pdf",
        plot_fl_core_mllaf4_gender=DESEQ2_PLOT_DIR
        + "IndividualPC1_PC2_fetalcore_MLLAF4_gender.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/4.2_fetal_signature_in_leukemia_ind_PCA_fetalcore.R"


rule fetal_signature_in_leukemia_random_genes:
    input:
        adult_signature=DESEQ2_DIR + "FLcorePseudotech/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcorePseudotech/Fetal_signature_p005.txt",
        ball_tsv=rules.export_ball_excel_to_tsv.output.ball_tsv,
        ball_dir=config["external_dir"] + "BALL-1988S-HTSeq/",
    output:
        random_genes=expand(
            DESEQ2_PLOT_DIR + "random_{i}.csv", i=list(range(1, 6))
        ),
        plots_random_pc1_all=expand(
            DESEQ2_PLOT_DIR + "random_samePC1_all_{i}.pdf",
            i=list(range(1, 6))
        ),
        plots_random_pc1_type=expand(
            DESEQ2_PLOT_DIR + "random_PC1_{leuk_type}_{i}.pdf",
            i=list(range(1, 6)),
            leuk_type=config["leukemia_types"],
        ),
        plots_random_pc1_type_type=expand(
            DESEQ2_PLOT_DIR + "{leuk_type}_random_samePC1_{leuk_type}_{i}.pdf",
            i=list(range(1, 6)),
            leuk_type=config["leukemia_types"],
        ),
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/4.3_fetal_signature_in_leukemia_randomgenes.R"


rule plot_main_fig5:
    input:
        deg=DESEQ2_DIR + "FLcorePseudotech/DE_genes.csv",
        samples=config["external_dir"] + "iPS_ETVRUNX/dev_cell_samples.txt",
        fpkm=config["external_dir"] + "iPS_ETVRUNX/fpkm.tsv",
    output:
        plot_fig5=DESEQ2_PLOT_DIR + "Main_Figure_5a.pdf",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/5.1_dev_cell_heatmap.R"
