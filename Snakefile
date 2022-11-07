configfile: "config.yaml"


import json
import shutil


DESEQ2_DIR = config["processed_dir"] + "DESeq2/results/"
DESEQ2_PLOT_DIR = config["processed_dir"] + "DESeq2/images/"
NOTEBOOKS_PLOT_DIR = config["processed_dir"] + "/notebooks/Figures/"
GLOBAL_FIGURES_DIR = "figures/figures/"

NOTEBOOKS_PLOT_DIR = "data/processed/notebooks/Figures/"
DESEQ2_PLOT_DIR = "data/processed/DESeq2/images/"
GLOBAL_FIGURES_DIR = "figures/figures/"


figures = json.load(open("figures/figures.json", "r"))

def get_paper_figures_indices(figures_dict):
    pdfs = []
    svgs = []
    for k, v in figures_dict.items():
        if v["origin"] == "py":
            svgs.append(k)
        elif v["origin"] == "DEseq":
            pdfs.append(k)

    return()

in_out_pairs = []
ins = []
outs = []
ins_dict = {}
outs_dict = {}
indices = []
for figure in list(figures.keys()):
    # print(figure)
    if figures[figure]["origin"] not in ["py", "DEseq"]:
        ...
        # print(f"Skipping: {figure} {figures[figure]=}")
    else:
        if figures[figure]["origin"] == "py":
            path = NOTEBOOKS_PLOT_DIR
            ft = ".svg"
        elif figures[figure]["origin"] == "DEseq":
            path = DESEQ2_PLOT_DIR
            ft = ".pdf"
        i = 1
        for fn in figures[figure]["filenames"]:
            idx = ""
            infile = path + fn + ft
            idx = figure + ( str(i) if len(figures[figure]["filenames"]) > 1 else "")
            outfile = GLOBAL_FIGURES_DIR + idx + "_" + fn + ft
            # print(infile, "\t" + outfile)
            in_out_pairs.append((infile, outfile))
            ins.append(infile)
            outs.append(outfile)
            indices.append(idx)
            ins_dict[idx] = infile
            outs_dict[idx] = outfile
            i += 1

# print(ins)
# print(outs)

# key_use = ""
# for k, v in outs_dict.items():
#     # print(k, v)
#     if indices[2] in v:
#         key_use = k
#         print(key_use)
# print(ins_dict[key_use])
# print(outs_dict[key_use])
# assert ins_dict[key_use].split("/")[-1] in outs_dict[key_use]

# paper_figures =

# rule gather_figures:
#     input:
#         outs,


def get_file(wildcards):
    return ins_dict[wildcards.figure]


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
        DESEQ2_PLOT_DIR + "PC12_FL_BM_singel_cell_gates_top500.pdf",
        DESEQ2_PLOT_DIR + "PC1Age_All_leuk_FL_core.pdf",
        DESEQ2_PLOT_DIR + "IndividualPC1_PC2_fetalcore_MLLAF4.pdf",
        DESEQ2_PLOT_DIR + "Main_Figure_5a.pdf",
        expand(
            DESEQ2_DIR + "FLcoreSC/{cluster}_FL_specific_markers.csv",
            cluster=config["fl_clusters_to_use"],
        ),
        expand(
            DESEQ2_PLOT_DIR + "random_{i}.csv", i=list(range(1, 6))
        ),
        expand(
            DESEQ2_PLOT_DIR
            + "PC1_variance_test_all_leuk_in_PCA_its_{iters}.pdf",
            iters=config["random_params"]["iterations"],
        ),
        DESEQ2_DIR + "Pvals_anova_FL_up2.csv",
        outs,


rule gather_figures:
    input:
        # source_file=get_file,
        # source=ins_dict["{figure,Figure[0-9]+[a-z]+[0-9]+\_}"],
        ins=ins,
    output:
        outs=outs
    run:
        for out_file in output.outs:
            key_use = ""
            for k, v in outs_dict.items():
                # print(k, v)
                if out_file in v:
                    key_use = k
                    print(key_use)
            source = ins_dict[key_use]
            destination = outs_dict[key_use]
            assert source.split("/")[-1] in destination
            shutil.copy(source, destination)

        # destination=expand(outs_dict[figure],
        #                    figure=indices)
        # destination_file="figures/figures/{figure,Figure[0-9]+[a-z]+[0-9]+}_{etc}",
    # params:
    #     figure = indices
    # shell:
    #     "cp {input.source_file} {output.destination_file}"


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


rule single_cell_deg:
    input:
        cite_cell_matrices=expand(
            config["raw_dir"]
            + "cite_cell_matrices/{donor}_hpc/",
            donor=config["donors"],
        ),
        combined_metadata=config["raw_dir"]
        + "paper_specific/from_seurat/FL_BM_combined_metadata.csv",
        fl_metadata=config["raw_dir"]
        + "paper_specific/from_seurat/FL_combined_metadata.csv",
        mapping_predictions=config["processed_dir"]
        + "notebooks/mapping_predictions/RNA_FL_target_preds.json",
    output:
        cells=expand(
            config["interim_dir"] + "single_cell_deg/{donor}_cells.csv",
            donor=config["donors"],
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
            DESEQ2_PLOT_DIR + "random_samePC1_all_{i}.pdf", i=list(range(1, 6))
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


rule rev_stat_test_distributions_pc:
    input:
        adult_signature=DESEQ2_DIR + "FLcorePseudotech/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcorePseudotech/Fetal_signature_p005.txt",
        ball_tsv=rules.export_ball_excel_to_tsv.output.ball_tsv,
        ball_dir=config["external_dir"] + "BALL-1988S-HTSeq/",
    output:
        plots_random_pc1_variance=expand(
            DESEQ2_PLOT_DIR + "PC1_variance_test_{leuk_type}_its_{iters}.pdf",
            leuk_type=config["leukemia_types"],
            iters=config["random_params"]["iterations"],
        ),
        plot_same_pc1=expand(
            DESEQ2_PLOT_DIR
            + "PC1_variance_test_all_leuk_in_PCA_its_{iters}.pdf",
            iters=config["random_params"]["iterations"],
        ),
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/6.1_rev_stat_test_distributions_PC.R"


types_oi = ["KMT2A-AFF1", "BCR-ABL1", "High hyperdiploid"]

rule rev_fetal_core_box_plots:
    input:
        adult_signature=DESEQ2_DIR + "FLcorePseudotech/Adult_signature_p005.txt",
        fetal_signature=DESEQ2_DIR + "FLcorePseudotech/Fetal_signature_p005.txt",
        ball_tsv=rules.export_ball_excel_to_tsv.output.ball_tsv,
        ball_dir=config["external_dir"] + "BALL-1988S-HTSeq/",
    output:
        plots_up=expand(
            DESEQ2_PLOT_DIR + "Fetal_core_up_box_plot_{leuk_type}.pdf",
            leuk_type=types_oi,
        ),
        plots_down=expand(
            DESEQ2_PLOT_DIR + "Fetal_core_down_box_plot_{leuk_type}.pdf",
            leuk_type=types_oi,
        ),
        pvals_up=DESEQ2_DIR + "Pvals_anova_FL_up2.csv",
        pvals_down=DESEQ2_DIR + "Pvals_anova_FL_down2.csv",
    params:
        types_oi=types_oi
    conda:
        "envs/DESeq2.yaml"
    script:
        "src/DESeq2/6.2_rev_fetal_core_boxplots_final.R"
