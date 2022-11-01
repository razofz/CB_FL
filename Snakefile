configfile: "config.yaml"


DESEQ2_DIR = config["processed_dir"] + "DESeq2/results/"


rule gather_pops:
    input:
        expand(
            DESEQ2_DIR + "Normed_counts/{fl_core_version}/{cluster}_DESeq_normed.txt",
            cluster=config["fl_clusters_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.csv",
            cluster=config["fl_clusters_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "Normed_counts/{fl_core_version}/gated_{gate}_DESeq_normed.txt",
            gate=config["gates_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        expand(
            DESEQ2_DIR + "{fl_core_version}/gated_{gate}_FL_yBM.csv",
            gate=config["gates_to_use"],
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),


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
            DESEQ2_DIR + "Normed_counts/{fl_core_version}/{cluster}_DESeq_normed.txt",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        deseq_results=expand(
            DESEQ2_DIR + "{fl_core_version}/{cluster}_FL_yBM.csv",
            cluster=config["fl_clusters_to_use"],
            allow_missing=True,
        ),
        fl_only_normed_counts=expand(
            DESEQ2_DIR + "Normed_counts/{fl_core_version}/only_FL_DESeq_normed.txt",
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
            DESEQ2_DIR + "Normed_counts/{fl_core_version}/gated_{gate}_DESeq_normed.txt",
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


rule gather_cores:
    input:
        fetal_signature=expand(
            DESEQ2_DIR + "{fl_core_version}/Fetal_signature_p005.txt",
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),
        adult_signature=expand(
            DESEQ2_DIR + "{fl_core_version}/Adult_signature_p005.txt",
            fl_core_version=[v for v in config["fl_cores"] if "eudo" in v],
        ),


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
        adult_signature=DESEQ2_DIR + "{fl_core_version}/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "{fl_core_version}/Fetal_signature_p005.txt",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.3a_produce_FL_core.R"


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
        adult_signature=DESEQ2_DIR + "{fl_core_version}/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "{fl_core_version}/Fetal_signature_p005.txt",
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/1.3b_produce_sc_FL_core.R"
