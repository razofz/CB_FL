# N. B.: scripts 2.1 (this), 2.2 and 2.3 must all be run in the same namespace.
# E. g. start an R console and `source()` the scripts one after another. - raz

source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))

out_dir <- "data/processed/DEseq2"
# path_to_files <- "data/processed/notebooks/DEseq2/cluster_wise/"
path_to_files <- "data/raw/paper_specific/DEseq2/cluster_wise"
images_dir <- paste0(out_dir, "/images")
if (dir.exists(images_dir) == FALSE) {
  dir.create(images_dir, showWarnings = TRUE, recursive = TRUE)
}
raw_dir <- "data/raw/paper_specific/Supplemental_tables"

files <- list.files(
  path = path_to_files, pattern = "*_coldata_all_hpc.csv",
  full.names = TRUE, recursive = FALSE
)
cts_combined <- data.frame()
coldata_combined <- data.frame()
for (x in 1:length(files)) {
  spl <- str_split(files[x], "/", simplify = T)
  names <- str_split(spl[length(spl)], "_coldata_all_hpc.csv", simplify = T)[1]
  print(paste0("> ", names))
  cts <- read.csv(paste0(path_to_files, "/", names, "_cts_all_hpc.csv"),
    row.names = 1
  )
  coldata <- read.csv(files[x], row.names = 1)
  if (x == 1) {
    cts_combined <- cts
    coldata_combined <- coldata
  } else {
    cts_combined <- cbind(cts_combined, cts)
    coldata_combined <- rbind(coldata_combined, coldata)
  }
}

cts_combined <- cts_combined[!grepl(
  pattern = "DC.I", x =
    colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "_Ly.III_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "DC.I", x =
    rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "_Ly.III_", x =
    rownames(coldata_combined)
), ]
cts_combined <- cts_combined[!grepl(
  pattern = "yBM_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "yBM_", x =
    rownames(coldata_combined)
), ]
cts_combined_clust <- cts_combined
coldata_combined_clust <- coldata_combined


dds <- DESeqDataSetFromMatrix(
  countData = cts_combined, colData = coldata_combined,
  design = ~1
)
dds <- DESeq(dds)
vsd_FL <- vst(dds, blind = TRUE)

rv_FL <- rowVars(assay(vsd_FL)) # calculate variance
select_FL <- order(rv_FL, decreasing = T)[seq_len(min(500, length(rv_FL)))] # find the top 500 genes by variance
pca_FL <- prcomp(t(assay(vsd_FL)[select_FL, ])) # calculate PCA
percentVar_FL <- pca_FL$sdev^2 / sum(pca_FL$sdev^2) # calc PCA contribution
intgroup_FL <- c("sample", "cluster_name")
intgroup_FL.df <- as.data.frame(colData(vsd_FL)[, intgroup_FL, drop = FALSE]) # add grouping if wanted
group <- if (length(intgroup_FL) > 1) {
  factor(apply(intgroup_FL.df, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL)[[intgroup_FL]]
}
d_FL <- data.frame(
  PC1 = pca_FL$x[, 1], PC2 = pca_FL$x[, 2], PC3 = pca_FL$x[ , 3
  ], group = group, intgroup_FL.df, name = colnames(vsd_FL)
)
loading_FL <- data.frame(
  LO1 = pca_FL$rotation[, 1], LO2 = pca_FL$rotation[, 2], LO3 = pca_FL$rotation[, 3],
  name = rownames(pca_FL$rotation)
)

colors_scale <- c(
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

ggplot(data = d_FL, aes_string(x = "PC1", y = "PC2", shape = "sample")) +
  geom_point(aes_string(fill = "cluster_name"), color = "black", shape = 21, size = 5) +
  xlab(paste0("PC1: ", round(percentVar_FL[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # geom_text(aes(label=cluster_name), size=3, nudge_x = 1, nudge_y = 1)+
  coord_fixed() +
  scale_fill_manual(values = colors_scale)
ggsave(paste0(images_dir, "/PC12_singel_cell_clusters_top500.pdf"), plot = last_plot())

# ggplot(data = loading_FL, aes_string(x = "LO1", y = "LO2")) +
#   geom_point(size = 3) +
#   geom_text(aes(label = name), size = 3, nudge_x = 0.005, nudge_y = 0.005) +
#   theme_bw()

ggplot(data = d_FL, aes_string(x = "PC2", y = "PC3", shape = "sample")) +
  geom_point(aes_string(fill = "cluster_name"), color = "black", shape = 21, size = 5) +
  xlab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar_FL[3] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # geom_text(aes(label=cluster_name), size=3, nudge_x = 1, nudge_y = 1)+
  coord_fixed() +
  scale_fill_manual(values = colors_scale)
ggsave(paste0(images_dir, "/PC23_singel_cell_clusters_top500.pdf"), plot = last_plot())
