invisible(lapply(list(
  "stringr",
  "biomaRt",
  "ggplot2",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

clusters <- snakemake@config[["fl_clusters_to_use"]]
cts_files <- snakemake@input[["cts_files"]]
coldata_files <- snakemake@input[["coldata_files"]]

print(clusters)

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

named_cts_files <- make_named_version(ids = clusters, file_list = cts_files)
named_coldata_files <- make_named_version(
  ids = clusters, file_list = coldata_files
)

cts_combined <- data.frame()
coldata_combined <- data.frame()
for (cluster_idx in seq_len(length(clusters))) {
  cts <- read.csv(named_cts_files[clusters[cluster_idx]],
    row.names = 1
  )
  coldata <- read.csv(named_coldata_files[clusters[cluster_idx]],
    row.names = 1
  )
  if (cluster_idx == 1) {
    cts_combined <- cts
    coldata_combined <- coldata
  } else {
    cts_combined <- cbind(cts_combined, cts)
    coldata_combined <- rbind(coldata_combined, coldata)
  }
}

cts_combined <- cts_combined[!grepl(
  pattern = "Ly.III", x = colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "DC.I", x = colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "_T_", x = colnames(cts_combined)
)]

coldata_combined <- coldata_combined[!grepl(
  pattern = "DC.I", x = rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "Ly.III", x = rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "_T_", x = rownames(coldata_combined)
), ]


dds <- DESeqDataSetFromMatrix(
  countData = cts_combined,
  colData = coldata_combined,
  design = ~1
)
dds <- DESeq(dds)
vsd_fl <- vst(dds, blind = TRUE)

rv_fl <- rowVars(assay(vsd_fl)) # calculate variance
select_fl <- order(rv_fl, decreasing = T)[seq_len(min(500, length(rv_fl)))]
# find the top 500 genes by variance
pca_fl <- prcomp(t(assay(vsd_fl)[select_fl, ])) # calculate PCA
percentVar_fl <- pca_fl$sdev^2 / sum(pca_fl$sdev^2) # calc PCA contribution
intgroup_fl <- c("sample", "cluster_name")
intgroup_fl.df <- as.data.frame(colData(vsd_fl)[, intgroup_fl, drop = FALSE])
# add grouping if wanted
group <- if (length(intgroup_fl) > 1) {
  factor(apply(intgroup_fl.df, 1, paste, collapse = ":"))
} else {
  colData(vsd_fl)[[intgroup_fl]]
}
d_fl <- data.frame(
  PC1 = pca_fl$x[, 1],
  PC2 = pca_fl$x[, 2],
  PC3 = pca_fl$x[, 3],
  group = group, intgroup_fl.df, name = colnames(vsd_fl)
)
loading_fl <- data.frame(
  LO1 = pca_fl$rotation[, 1],
  LO2 = pca_fl$rotation[, 2],
  LO3 = pca_fl$rotation[, 3],
  name = rownames(pca_fl$rotation)
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
  "#addfff"
) # MPPII

for (pcs in list(1:2, 2:3)) {
  p <- ggplot(
    data = d_fl,
    aes_string(
      x = str_c("PC", pcs[[1]]),
      y = str_c("PC", pcs[[2]])
    )
  ) +
    geom_point(
      aes_string(fill = "cluster_name", shape = "sample"),
      size = 5
    ) +
    xlab(paste0(
      "PC", pcs[[1]], ": ", round(percentVar_fl[pcs[[1]]] * 100),
      "% variance"
    )) +
    ylab(paste0(
      "PC", pcs[[2]], ": ", round(percentVar_fl[pcs[[2]]] * 100),
      "% variance"
    )) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    coord_fixed() +
    scale_colour_manual(
      values = colors_scale,
      aesthetics = c("color", "fill")
    ) +
    scale_shape_manual(values = c(21, 24)) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21)),
      shape = guide_legend(override.aes = list(fill = "black"))
    )
  ggsave(
    snakemake@output[[str_c("plot_pc", str_flatten(pcs))]],
    plot = p
  )
}
