invisible(lapply(list(
  "stringr",
  "DESeq2",
  "biomaRt",
  "limma",
  "ggplot2",
  "pheatmap"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

clusters <- snakemake@config[["fl_clusters_to_use"]]
gates <- snakemake@config[["gates_to_use"]]

cts_files_cluster <- snakemake@input[["cts_files_cluster"]]
coldata_files_cluster <- snakemake@input[["coldata_files_cluster"]]
cts_files_gate <- snakemake@input[["cts_files_gate"]]
coldata_files_gate <- snakemake@input[["coldata_files_gate"]]

make_named_version <- function(ids, file_list, pattern = "_") {
  named_file_list <- unlist(
    lapply(ids,
      FUN = function(id) {
        file_list[str_detect(string = file_list, pattern = str_c(id, pattern))]
      }
    )
  )
  names(named_file_list) <- ids
  return(named_file_list)
}

named_cts_files_cluster <- make_named_version(
  ids = clusters, file_list = cts_files_cluster
)
named_cts_files_gate <- make_named_version(
  ids = gates, file_list = cts_files_gate
)
named_coldata_files_cluster <- make_named_version(
  ids = clusters, file_list = coldata_files_cluster
)
named_coldata_files_gate <- make_named_version(
  ids = gates, file_list = coldata_files_gate
)

###################
#  colour scales  #
###################

colors_scale_cluster <- c(
  "#efefd9", # Cyc
  "#5dcbd6", # DC.mono
  "#6469ff", # GMP
  "#fff148", # HSC
  "#d18ce2", # LyI
  "#dd00ff", # LyII
  "#ff5a5a", # MEP
  "#ffb4b4", # MPPI
  "#addfff" # MPPII
)

colors_scale_gate <- c(
  "#dd00ff", # CD10
  "#addfff", # MPP
  "#fff148", # 49F
  "#F3A3E3", # CD7
  "#F9CF3F", # 90
  "#5dcbd6", # CMP123
  "#76F5F2", # CMP135
  "#6469ff", # GMP123
  "#AEB1FC", # GMP135
  "#EDC5F7", # IL7
  "#d18ce2", # LMPP
  "#ffb4b4", # MEP123
  "#ff5a5a"
) # MEP135

colors_scale_both <- c(
  "#dd00ff", # CD10
  "#addfff", # MPP
  "#fff148", # 49F
  "#F3A3E3", # CD7
  "#F9CF3F", # 90
  "#5dcbd6", # CMP123
  "#76F5F2", # CMP135
  "#efefd9", # Cyc
  "#5dcbd6", # DC.mono
  "#6469ff", # GMP
  "#6469ff", # GMP123
  "#AEB1FC", # GMP135
  "#fff148", # HSC
  "#EDC5F7", # IL7
  "#d18ce2", # LMPP
  "#d18ce2", # LyI
  "#dd00ff", # LyII
  "#ff5a5a", # MEP
  "#ffb4b4", # MEP123
  "#ff5a5a", # MEP135
  "#ffb4b4", # MPPI
  "#addfff" # MPPII
)

################################################################################
#                                  Load files                                  #
################################################################################

##############
#  clusters  #
##############

cts_combined_cluster <- data.frame()
coldata_combined_cluster <- data.frame()
for (cluster_idx in seq_len(length(clusters))) {
  cts <- read.csv(
    named_cts_files_cluster[clusters[cluster_idx]],
    row.names = 1
  )
  coldata <- read.csv(
    named_coldata_files_cluster[clusters[cluster_idx]],
    row.names = 1
  )
  if (cluster_idx == 1) {
    cts_combined_cluster <- cts
    coldata_combined_cluster <- coldata
  } else {
    cts_combined_cluster <- cbind(cts_combined_cluster, cts)
    coldata_combined_cluster <- rbind(coldata_combined_cluster, coldata)
  }
}

cts_combined_cluster <- cts_combined_cluster[!grepl(
  pattern = "DC.I",
  x = colnames(cts_combined_cluster)
)]
cts_combined_cluster <- cts_combined_cluster[!grepl(
  pattern = "_Ly.III_",
  x = colnames(cts_combined_cluster)
)]
coldata_combined_cluster <- coldata_combined_cluster[!grepl(
  pattern = "DC.I",
  x = rownames(coldata_combined_cluster)
), ]
coldata_combined_cluster <- coldata_combined_cluster[!grepl(
  pattern = "_Ly.III_",
  x = rownames(coldata_combined_cluster)
), ]
cts_combined_cluster <- cts_combined_cluster[!grepl(
  pattern = "yBM_",
  x = colnames(cts_combined_cluster)
)]
coldata_combined_cluster <- coldata_combined_cluster[!grepl(
  pattern = "yBM_",
  x = rownames(coldata_combined_cluster)
), ]

###########
#  gates  #
###########

cts_combined_gate <- data.frame()
coldata_combined_gate <- data.frame()
for (gate_idx in seq_len(length(gates))) {
  cts <- read.csv(
    named_cts_files_gate[gates[gate_idx]],
    row.names = 1
  )
  coldata <- read.csv(
    named_coldata_files_gate[gates[gate_idx]],
    row.names = 1
  )
  if (gate_idx == 1) {
    cts_combined_gate <- cts
    coldata_combined_gate <- coldata
  } else {
    cts_combined_gate <- cbind(cts_combined_gate, cts)
    coldata_combined_gate <- rbind(coldata_combined_gate, coldata)
  }
}

cts_combined_gate <- cts_combined_gate[!grepl(
  pattern = "yBM_",
  x = colnames(cts_combined_gate)
)]
coldata_combined_gate <- coldata_combined_gate[!grepl(
  pattern = "yBM_",
  x = rownames(coldata_combined_gate)
), ]

##########
#  both  #
##########

coldata_combined_both_cluster <- coldata_combined_cluster
coldata_combined_both_gate <- coldata_combined_gate

colnames(coldata_combined_both_cluster) <- c(
  "sample", "population", "pseudo_rep"
)
coldata_combined_both_cluster$group_type <- rep(
  x = "Cluster",
  dim(coldata_combined_both_cluster)[1]
)
colnames(coldata_combined_both_gate) <- c("sample", "population", "pseudo_rep")
coldata_combined_both_gate$group_type <- rep(
  x = "Gated",
  dim(coldata_combined_both_gate)[1]
)
coldata_combined_gate_cluster <- rbind.data.frame(
  coldata_combined_both_cluster,
  coldata_combined_both_gate
)

cts_combined_gate_cluster <- cbind.data.frame(
  cts_combined_cluster,
  cts_combined_gate
)

################################################################################
#                                 DESeq2 part                                  #
################################################################################

########################
#  deseq for clusters  #
########################

dds_cluster <- DESeqDataSetFromMatrix(
  countData = cts_combined_cluster,
  colData = coldata_combined_cluster,
  design = ~1
)
dds_cluster <- DESeq(dds_cluster)
vsd_FL_cluster <- vst(dds_cluster, blind = TRUE)

rv_FL_cluster <- rowVars(assay(vsd_FL_cluster)) # calculate variance
select_FL_cluster <- order(rv_FL_cluster, decreasing = T)[
  seq_len(min(500, length(rv_FL_cluster)))
] # find the top 500 genes by variance
pca_FL_cluster <- prcomp(t(assay(vsd_FL_cluster)[select_FL_cluster, ])) # calculate PCA

# calc PCA contribution
percentVar_FL_cluster <- pca_FL_cluster$sdev^2 / sum(pca_FL_cluster$sdev^2)

intgroup_FL_cluster <- c("sample", "cluster_name")
intgroup_FL.df_cluster <- as.data.frame(
  colData(vsd_FL_cluster)[, intgroup_FL_cluster, drop = FALSE]
) # add grouping if wanted

group_cluster <- if (length(intgroup_FL_cluster) > 1) {
  factor(apply(intgroup_FL.df_cluster, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL_cluster)[[intgroup_FL_cluster]]
}
d_FL_cluster <- data.frame(
  PC1 = pca_FL_cluster$x[, 1],
  PC2 = pca_FL_cluster$x[, 2],
  PC3 = pca_FL_cluster$x[, 3],
  group = group_cluster, intgroup_FL.df_cluster, name = colnames(vsd_FL_cluster)
)
loading_FL_cluster <- data.frame(
  LO1 = pca_FL_cluster$rotation[, 1],
  LO2 = pca_FL_cluster$rotation[, 2],
  LO3 = pca_FL_cluster$rotation[, 3],
  name = rownames(pca_FL_cluster$rotation)
)

#####################
#  deseq for gates  #
#####################

dds_gate <- DESeqDataSetFromMatrix(
  countData = cts_combined_gate,
  colData = coldata_combined_gate,
  design = ~1
)
dds_gate <- DESeq(dds_gate)
vsd_FL_gate <- vst(dds_gate, blind = TRUE)

rv_FL_gate <- rowVars(assay(vsd_FL_gate)) # calculate variance
select_FL_gate <- order(rv_FL_gate, decreasing = T)[
  seq_len(min(500, length(rv_FL_gate)))
] # find the top 500 genes by variance
pca_FL_gate <- prcomp(t(assay(vsd_FL_gate)[select_FL_gate, ]))

percentVar_FL_gate <- pca_FL_gate$sdev^2 / sum(pca_FL_gate$sdev^2) # calc PCA contribution

intgroup_FL_gate <- c("sample", "gate_name")
# add grouping if wanted
intgroup_FL.df_gate <- as.data.frame(
  colData(vsd_FL_gate)[, intgroup_FL_gate, drop = FALSE]
)

group_gate <- if (length(intgroup_FL_gate) > 1) {
  factor(apply(intgroup_FL.df_gate, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL_gate)[[intgroup_FL_gate]]
}

d_FL_gate <- data.frame(
  PC1 = pca_FL_gate$x[, 1],
  PC2 = pca_FL_gate$x[, 2],
  PC3 = pca_FL_gate$x[, 3],
  group = group_gate, intgroup_FL.df_gate, name = colnames(vsd_FL_gate)
)
loading_FL_gate <- data.frame(
  LO1 = pca_FL_gate$rotation[, 1],
  LO2 = pca_FL_gate$rotation[, 2],
  LO3 = pca_FL_gate$rotation[, 3],
  name = rownames(pca_FL_gate$rotation)
)

####################
#  deseq for both  #
####################

dds_both <- DESeqDataSetFromMatrix(
  countData = cts_combined_gate_cluster,
  colData = coldata_combined_gate_cluster,
  design = ~1
)
dds_both <- DESeq(dds_both)
vsd_FL_both <- vst(dds_both, blind = TRUE)

rv_FL_both <- rowVars(assay(vsd_FL_both)) # calculate variance
# find the top 500 genes by variance
select_FL_both <- order(rv_FL_both, decreasing = T)[
  seq_len(min(500, length(rv_FL_both)))
]
pca_FL_both <- prcomp(t(assay(vsd_FL_both)[select_FL_both, ])) # calculate PCA
percentVar_FL_both <- pca_FL_both$sdev^2 / sum(pca_FL_both$sdev^2) # calc PCA contribution
intgroup_FL_both <- c("sample", "population", "group_type")
# add grouping if wanted
intgroup_FL.df_both <- as.data.frame(
  colData(vsd_FL_both)[, intgroup_FL_both,
    drop = FALSE
  ]
)
group_both <- if (length(intgroup_FL_both) > 1) {
  factor(apply(intgroup_FL.df_both, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL_both)[[intgroup_FL_both]]
}
d_FL_both <- data.frame(
  PC1 = pca_FL_both$x[, 1],
  PC2 = pca_FL_both$x[, 2],
  PC3 = pca_FL_both$x[, 3],
  group = group_both,
  intgroup_FL.df_both,
  name = colnames(vsd_FL_both)
)
loading_FL_both <- data.frame(
  LO1 = pca_FL_both$rotation[, 1],
  LO2 = pca_FL_both$rotation[, 2],
  LO3 = pca_FL_both$rotation[, 3],
  name = rownames(pca_FL_both$rotation)
)

################################################################################
#                                    Plots                                     #
################################################################################

#################
#  for cluster  #
#################

for (pcs in list(1:2, 2:3)) {
  p <- ggplot(
    data = d_FL_cluster,
    aes_string(
      x = str_c("PC", pcs[[1]]),
      y = str_c("PC", pcs[[2]]),
      shape = "sample"
    )
  ) +
    geom_point(
      aes_string(fill = "cluster_name"),
      color = "black", shape = 21, size = 5
    ) +
    xlab(paste0(
      "PC", pcs[[1]], ": ", round(percentVar_FL_cluster[pcs[[1]]] * 100),
      "% variance"
    )) +
    ylab(paste0(
      "PC", pcs[[2]], ": ", round(percentVar_FL_cluster[pcs[[2]]] * 100),
      "% variance"
    )) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed() +
    scale_fill_manual(values = colors_scale_cluster) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed()
  ggsave(
    snakemake@output[[str_c("plot_pc", str_flatten(pcs), "_cluster")]],
    plot = p
  )
}

###############
#  for gates  #
###############

for (pcs in list(1:2, 2:3)) {
  p <- ggplot(
    data = d_FL_gate,
    aes_string(
      x = str_c("PC", pcs[[1]]),
      y = str_c("PC", pcs[[2]]),
      shape = "sample"
    )
  ) +
    geom_point(
      aes_string(fill = "gate_name"),
      color = "black", shape = 21, size = 5
    ) +
    xlab(paste0(
      "PC", pcs[[1]], ": ", round(percentVar_FL_gate[pcs[[1]]] * 100),
      "% variance"
    )) +
    ylab(paste0(
      "PC", pcs[[2]], ": ", round(percentVar_FL_gate[pcs[[2]]] * 100),
      "% variance"
    )) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed() +
    scale_fill_manual(values = colors_scale_gate) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed()
  ggsave(
    snakemake@output[[str_c("plot_pc", str_flatten(pcs), "_gate")]],
    plot = p
  )
}

##############
#  for both  #
##############

for (pcs in list(1:2, 2:3)) {
  p <- ggplot(
    data = d_FL_both,
    aes_string(
      x = str_c("PC", pcs[[1]]),
      y = str_c("PC", pcs[[2]]),
      shape = "sample"
    )
  ) +
    geom_point(
      aes_string(fill = "population", shape = "group_type"),
      size = 5
    ) +
    xlab(paste0(
      "PC", pcs[[1]], ": ", round(percentVar_FL_both[pcs[[1]]] * 100),
      "% variance"
    )) +
    ylab(paste0(
      "PC", pcs[[2]],": ", round(percentVar_FL_both[pcs[[2]]] * 100),
      "% variance"
    )) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed() +
    scale_colour_manual(
      values = colors_scale_both,
      aesthetics = c("color", "fill")
    ) +
    scale_shape_manual(values = c(21, 24)) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21)),
      shape = guide_legend(override.aes = list(fill = "black"))
    ) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed()
  ggsave(
    snakemake@output[[str_c("plot_pc", str_flatten(pcs), "_both")]],
    plot = p,
    width = 7,
    height = 7,
    units = "in"
  )
}
