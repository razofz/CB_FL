# test of combined gating and clusters pca on deseq2 normalized values
# must run Cluster and gate analysis first

# N. B.: scripts 2.1, 2.2 and 2.3 must all be run in the same namespace.
# E. g. start an R console and `source()` the scripts one after another. - raz

source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("limma"))

out_dir <- "data/processed/DEseq2"
images_dir <- paste0(out_dir, "/images")
if (dir.exists(images_dir) == FALSE) {
  dir.create(images_dir, showWarnings = TRUE, recursive = TRUE)
}

colnames(coldata_combined_clust) <- c("sample", "population", "pseudo_rep")
coldata_combined_clust$group_type <- rep(
  x = "Cluster",
  dim(coldata_combined_clust)[1]
)
colnames(coldata_combined_gate) <- c("sample", "population", "pseudo_rep")
coldata_combined_gate$group_type <- rep(
  x = "Gated",
  dim(coldata_combined_gate)[1]
)
coldata_combined_gate_clust <- rbind.data.frame(
  coldata_combined_clust,
  coldata_combined_gate
)

cts_combined_gate_clust <- cbind.data.frame(
  cts_combined_clust,
  cts_combined_gate
)

dds <- DESeqDataSetFromMatrix(
  countData = cts_combined_gate_clust, colData = coldata_combined_gate_clust,
  design = ~1
)
dds <- DESeq(dds)
vsd_FL <- vst(dds, blind = TRUE)

rv_FL <- rowVars(assay(vsd_FL)) # calculate variance
# find the top 500 genes by variance
select_FL <- order(rv_FL, decreasing = T)[seq_len(min(500, length(rv_FL)))]
pca_FL <- prcomp(t(assay(vsd_FL)[select_FL, ])) # calculate PCA
percentVar_FL <- pca_FL$sdev^2 / sum(pca_FL$sdev^2) # calc PCA contribution
intgroup_FL <- c("sample", "population", "group_type")
# add grouping if wanted
intgroup_FL.df <- as.data.frame(colData(vsd_FL)[, intgroup_FL, drop = FALSE])
group <- if (length(intgroup_FL) > 1) {
  factor(apply(intgroup_FL.df, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL)[[intgroup_FL]]
}
d_FL <- data.frame(PC1 = pca_FL$x[, 1], PC2 = pca_FL$x[, 2], PC3 = pca_FL$x[
  ,
  3
], group = group, intgroup_FL.df, name = colnames(vsd_FL))
loading_FL <- data.frame(
  LO1 = pca_FL$rotation[, 1], LO2 = pca_FL$rotation[
    ,
    2
  ], LO3 = pca_FL$rotation[, 3], name =
    rownames(pca_FL$rotation)
)

colors_scale3 <- c(
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

ggplot(data = d_FL, aes_string(x = "PC1", y = "PC2")) +
  geom_point(aes_string(fill = "population", shape = "group_type"), size = 5) +
  xlab(paste0("PC1: ", round(percentVar_FL[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # geom_text(aes(label=population), size=3, nudge_x = 1, nudge_y = 1)+
  coord_fixed() +
  scale_colour_manual(values = colors_scale3, aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(21, 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "black"))
  )
  # scale_colour_manual(values = colors_scale3, aesthetics = c("color", "fill"), guide = "none") +
  # scale_shape_manual(values = c(21, 24), guide = "none") #+
ggsave(paste0(images_dir, "/PC12_singel_cell_gates_clusters_top500.pdf"),
  plot = last_plot(),
  width = 7,
  height = 7,
  units = "in"
)

ggplot(data = d_FL, aes_string(x = "PC2", y = "PC3")) +
  geom_point(aes_string(fill = "population", shape = "group_type"), size = 5) +
  xlab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar_FL[3] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # geom_text(aes(label=population), size=3, nudge_x = 1, nudge_y = 1)+
  coord_fixed() +
  scale_colour_manual(values = colors_scale3, aesthetics = c("color", "fill")) +
  scale_shape_manual(values = c(21, 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "black"))
  )
  # scale_colour_manual(values = colors_scale3, aesthetics = c("color", "fill"), guide = "none") +
  # scale_shape_manual(values = c(21, 24), guide = "none") #+
ggsave(paste0(images_dir, "/PC23_singel_cell_gates_clusters_top500.pdf"),
  plot = last_plot(),
  width = 7,
  height = 7,
  units = "in"
)
