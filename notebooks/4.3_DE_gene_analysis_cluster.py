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

# + tags=[]
# %load_ext autotime
# %config InlineBackend.figure_format = 'retina'

import glob
import fl_utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use("fivethirtyeight")
plt.rcParams["svg.fonttype"] = "none"


def clean_axis(ax, ts=11, ga=0.6):
    ax.xaxis.set_tick_params(labelsize=ts)
    ax.yaxis.set_tick_params(labelsize=ts)
    for i in ["top", "bottom", "left", "right"]:
        ax.spines[i].set_visible(False)
    ax.grid(which="major", linestyle="--", alpha=ga)
    ax.figure.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    return True


# + tags=[]
fl_utils.set_project_path()

# + tags=[]
threads = 16

raw_dir = "data/raw"
out_dir = "data/processed/notebooks"
deseq2_dir = out_dir + "/DEseq2"
deseq2_clw_dir = deseq2_dir + "/cluster_wise"
deseq2_scripts_dir = deseq2_dir + "/scripts"
deseq2_results_dir = deseq2_dir + "/results"
deseq2_comparison_dir = raw_dir + "/paper_specific/DESeq2/FL_yBM_comparison/clusters"
# deseq2_comparison_dir = "data/processed/DEseq2/results"
sup_tables_dir = raw_dir + "/paper_specific/Supplemental_tables"
figures_dir = out_dir + "/Figures"

# + tags=[]
in_files = {
    "deseq2_comparison_dir": [y
        for y in [x for x in glob.glob(deseq2_comparison_dir + "/*.csv")
        if "Ly.III_" not in x] if "DC.I" not in y]
    # "deseq2_comparison_dir": glob.glob(deseq2_comparison_dir + "/clusters/*"),
}

out_files = {
    "plots": [
        f'{figures_dir}/DEseq_{fn.split("/")[-1].split(".csv")[0]}.svg'
        for fn in in_files["deseq2_comparison_dir"]
    ],
}

fl_utils.print_files(in_files, "in")
fl_utils.print_files(out_files, "out")

# + tags=[]
fl_utils.create_dir(deseq2_results_dir)


# -

# ---

# + tags=[]
def read_diff_data(fn):
    i = fn.split("/")[-1][:-4]
    df = pd.read_csv(fn)
    # print(f"read_diff_data: {df=}")
    return i, df


def filter_genes(names):
    # print(f"filter_genes: {names=}")
    return [
        x
        for x in names
        if not x.startswith("RPS")
        and not x.startswith("RPL")
        and not x.startswith("MT-")
    ]


def get_up_genes(df, p_thresh=0.05, fc_thresh=1):
    xdf = df[df.padj < p_thresh].sort_values("padj").log2FoldChange
    # print(f"get_up_genes: {df=}")
    return filter_genes(list(xdf[xdf > fc_thresh].index))


def get_down_genes(df, p_thresh=0.05, fc_thresh=-1):
    # print(f"get_down_genes: {df=}")
    xdf = df[df.padj < p_thresh].sort_values("padj").log2FoldChange
    return filter_genes(list(xdf[xdf < -1].index))


def get_gene_sets(sample_pair):
    up_set = []
    down_set = []
    for fn in in_files["deseq2_comparison_dir"]:
        # print(f"get_gene_sets: {fn=}")
        sample, df = read_diff_data(fn)
        if sample.endswith(sample_pair) is False:
            continue
        up_genes = get_up_genes(df)
        up_set.extend(up_genes)
        down_genes = get_down_genes(df)
        down_set.extend(down_genes)
    up_set = sorted(set(up_set))
    down_set = sorted(set(down_set))
    return up_set, down_set


# def merge_sets(*args):
#     tot = []
#     for i in args:
#         tot.extend(i)
#     return sorted(set(i))


def build_set_matrix(gene_set, sample_pair):
    data = {}
    # for fn in glob.glob('./DEseq2/results/FL_yBM_comparrison/clusters/*.csv'):
    # for fn in glob.glob(f"{deseq2_comparison_dir}/clusters/*.csv"):
    for fn in in_files["deseq2_comparison_dir"]:
        sample, df = read_diff_data(fn)
        if sample.endswith(sample_pair) is False:
            continue
        df = df[df.padj < 0.05]
        data[sample.rsplit("_", 1)[0]] = (
            df["log2FoldChange"].reindex(gene_set).fillna(0)
        )
    data = pd.DataFrame(data)
    # print(data)
    return data


# + tags=[]
fl_ybm_up, fl_ybm_down = get_gene_sets("FL_yBM")
print("FL_yBM", len(fl_ybm_up), len(fl_ybm_down))
len(fl_ybm_up + fl_ybm_down)

# + tags=[]
gene_sets = [fl_ybm_up, fl_ybm_down]
sample_pairs = ["FL_yBM", "FL_yBM"]
names = ["FL_yBM_up", "FL_yBM_down"]
for i, j, k in zip(gene_sets, sample_pairs, names):
    print(k)
    # print(f"{i=} {j=} {k=}")
    mat = build_set_matrix(i, j)
    if k.endswith("down"):
        sorted_genes = pd.DataFrame(
            [mat.idxmin(axis=1).sort_values(), (mat < 0).sum(axis=1)]
        ).T.sort_values(by=[1, 0])
        mat = mat.reindex(sorted_genes.index)
        mat[mat > 0] = 0
        cgx = sns.clustermap(
            mat,
            row_cluster=False,
            robust=True,
            figsize=(12, 5),
            xticklabels=mat.columns,
            cmap="OrRd_r",
            rasterized=True,
        )
    elif k.endswith("up"):
        sorted_genes = pd.DataFrame(
            [mat.idxmax(axis=1).sort_values(), (mat > 0).sum(axis=1)]
        ).T.sort_values(by=[1, 0])
        mat = mat.reindex(sorted_genes.index)
        mat[mat < 0] = 0
        cgx = sns.clustermap(
            mat,
            row_cluster=False,
            robust=True,
            figsize=(12, 5),
            xticklabels=mat.columns,
            cmap="BuPu",
            rasterized=True,
        )
    else:
        print("oh no!")
    # mat.to_csv(f"./deseq2_analysis/downstream/{k}.csv")
    # plt.savefig(f"./images/deseq_diff_heatmaps/{k}.svg", dpi=300)
    plt.show()


# + tags=[]
def plot_vulc(name, res, ciel=20, minp=3, mino=2, display_names=None):
    x, y = res["log2FoldChange"].values, -np.log10(res["padj"].values)
    y[y > ciel] = ciel
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    ax.scatter(
        x[y > -np.log10(0.05)],
        y[y > -np.log10(0.05)],
        rasterized=True,
        s=5,
        c="crimson",
    )
    ax.scatter(
        x[y <= -np.log10(0.05)],
        y[y <= -np.log10(0.05)],
        rasterized=True,
        s=3,
        c="k",
    )
    for i, j, k in zip(x, y, res.index):
        if display_names is not None:
            if k in display_names:
                ax.text(i, j, k, fontsize=10)
        elif j > minp or i > mino:
            ax.text(i, j, k, fontsize=10)
    ax.axhline(-np.log10(0.05), lw=1, ls="--", c="k")
    ax.axvline(1, lw=1, ls="--", c="k")
    ax.axvline(-1, lw=1, ls="--", c="k")
    clean_axis(ax)
    ax.set_xlabel("Log2 FC")
    ax.set_ylabel("-Log10 p-value")
    # plt.savefig(f"{figures_dir}/SVG/DEseq_{name}.svg", dpi=300)
    plt.savefig(f"{figures_dir}/DEseq_{name}.svg", dpi=300)
    plt.show()


# + tags=[]
# for fn in glob.glob('./DEseq2/results/FL_yBM_comparrison/clusters/*.csv'):
# for fn in glob.glob(f"{deseq2_results_dir}/FL_yBM_comparison/clusters/*.csv"):
for fn in in_files["deseq2_comparison_dir"]:
    name = fn.split("/")[-1].split(".csv")[0]
    print(name)
    res2 = pd.read_csv(fn)
    res2["log2FoldChange"] = res2["log2FoldChange"].fillna(0)
    res2["padj"] = res2["padj"].fillna(1)
    res2["padj"][res2["padj"] == 0] = res2.padj[res2.padj > 0].min()
    plot_vulc(
        name,
        res=res2,
        ciel=100,
        display_names={
            x: None
            for x in res2[
                (abs(res2["log2FoldChange"]) > 1) & (res2["padj"] < 1e-10)
            ].index
        },
    )

# + tags=[]
gene_sets = [fl_ybm_up, fl_ybm_down]
sample_pairs = ["FL_yBM", "FL_yBM"]
mat = build_set_matrix(fl_ybm_up + fl_ybm_down, "FL_yBM")

# + tags=[]
n_diff_tot = []
n_diff_up = []
n_diff_dn = []
for i, n in enumerate(mat.columns):
    print(i, n)
    n_diff_tot.append(len(mat[n][abs(mat[n]) >= 1]))
    n_diff_up.append(len(mat[n][mat[n] >= 1]))
    n_diff_dn.append(len(mat[n][mat[n] <= -1]))
n_diff = pd.DataFrame()
n_diff["tot"] = n_diff_tot
n_diff["up"] = n_diff_up
n_diff["dn"] = n_diff_dn
n_diff.index = [x.split("_")[0] for x in mat.columns]

# + tags=[]
for i in n_diff.columns:
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    print(i)
    ax.bar(
        n_diff.sort_values(by=[i], ascending=False).index,
        n_diff.sort_values(by=[i], ascending=False)[i],
        color="grey",
        edgecolor="k",
    )
    ax.set_xticks(np.array(range(len(n_diff.index))) + 0.25)
    ax.set_xticklabels(
        [
            x.replace(".", "-")
            for x in n_diff.sort_values(by=[i], ascending=False).index
        ],
        rotation=45,
        ha="right",
        fontsize=12,
    )
    ax.set_ylabel("n genes", fontsize=12)
    clean_axis(ax)
    plt.tight_layout()
    # plt.savefig(f'./images/SVG/bar_n_diff_genes_{i}.svg', dpi=300)
    plt.show()

# +
n_diff_tot = {}
n_diff_up = {}
n_diff_dn = {}
for n in mat.columns:
    print(n)
    n_diff_tot[n.split("_")[0]] = list(mat[n][abs(mat[n]) >= 1].index)
    n_diff_up[n.split("_")[0]] = list(mat[n][mat[n] >= 1].index)
    n_diff_dn[n.split("_")[0]] = list(mat[n][mat[n] <= -1].index)

ly_genes = []
my_genes = []
prim_genes = []
for x in [n_diff_tot, n_diff_up, n_diff_dn]:
    ly_genes.append(len(pd.unique(x["Ly.I"] + x["Ly.II"])))
    my_genes.append(len(pd.unique(x["Ly.I"] + x["DC.Mono"] + x["GMP"])))
    prim_genes.append(
        len(pd.unique(x["HSC"] + x["MPP.I"] + x["MPP.II"] + x["Cyc"]))
    )
n_diff = pd.DataFrame()
n_diff["ly_genes"] = ly_genes
n_diff["my_genes"] = my_genes
n_diff["prim_genes"] = prim_genes
n_diff.index = ["tot", "up", "dn"]
n_diff = n_diff.T
n_diff

# + tags=[]
for i in n_diff_up.keys():
    print(i)
    for x in n_diff_up[i]:
        print(x)

# +
all_diff_gene_up = pd.DataFrame.from_dict(n_diff_up, orient="index")
all_diff_gene_up = all_diff_gene_up.T

# all_diff_gene_up.to_csv(out_files["diff_genes_up"], index=False)
# all_diff_gene_up.to_csv(f"{sup_tables_dir}/all_diff_genes_up.csv", index=False)
# all_diff_gene_up.to_csv('./Supplemental_tables/all_diff_genes_up.csv',index=False)

# +
all_diff_gene_dn = pd.DataFrame.from_dict(n_diff_dn, orient="index")
all_diff_gene_dn = all_diff_gene_dn.T

# all_diff_gene_dn.to_csv(out_files["diff_genes_dn"], index=False)
# all_diff_gene_dn.to_csv(f"{sup_tables_dir}/all_diff_genes_dn.csv", index=False)
# all_diff_gene_dn.to_csv('./Supplemental_tables/all_diff_genes_dn.csv',index=False)

# + tags=[]
for i in n_diff.columns:
    fig, ax = plt.subplots(1, 1, figsize=(2.5, 3))
    print(i)
    ax.bar(
        n_diff.sort_values(by=[i], ascending=False).index,
        n_diff.sort_values(by=[i], ascending=False)[i],
        color="grey",
        edgecolor="k",
    )
    ax.set_xticks(np.array(range(len(n_diff.index))) + 0.25)
    ax.set_xticklabels(
        [
            x.replace(".", "-")
            for x in n_diff.sort_values(by=[i], ascending=False).index
        ],
        rotation=45,
        ha="right",
        fontsize=12,
    )
    ax.set_ylabel("n genes", fontsize=12)
    clean_axis(ax)
    plt.tight_layout()
    # plt.savefig(f'./images/SVG/bar_n_diff_genes_{i}.svg', dpi=300)
    plt.show()

# + tags=[]
fig, ax = plt.subplots(1, 1, figsize=(2.5, 3))
colors = {"tot": "r", "dn": "b"}
n_diff_ord = n_diff.sort_values(by=["tot"], ascending=False)
ax.bar(n_diff_ord.index, n_diff_ord["tot"], color=colors["tot"], edgecolor="k")
ax.bar(n_diff_ord.index, n_diff_ord["dn"], color=colors["dn"], edgecolor="k")

ax.set_xticks(np.array(range(len(n_diff_ord.index))) + 0.25)
ax.set_xticklabels(
    [x.replace(".", "-") for x in n_diff_ord.index],
    rotation=45,
    ha="right",
    fontsize=12,
)
ax.set_ylabel("n genes", fontsize=12)
labels = list(colors.keys())
handles = [plt.Rectangle((0, 0), 1, 1, color=colors[label]) for label in labels]
labels = ["up", "down"]
plt.legend(handles, labels, loc="upper right")
clean_axis(ax)
# plt.savefig(f'./images/SVG/bar_n_diff_genes_prim_my_ly_up_dn.svg', dpi=300)
plt.tight_layout()
# plt.savefig(f'./images/PNG/bar_n_diff_genes_prim_my_ly_up_dn.png', dpi=300)
plt.show()

# +
# mat = mat.drop("DC.I_FL", axis=1)
# mat = mat.drop("T_FL", axis=1)
# -


