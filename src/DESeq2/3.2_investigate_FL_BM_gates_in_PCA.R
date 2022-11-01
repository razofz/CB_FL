invisible(lapply(list(
  "stringr",
  "biomaRt",
  "ggplot2",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

gates <- snakemake@config[["gates_to_use"]]
cts_files <- snakemake@input[["cts_files"]]
coldata_files <- snakemake@input[["coldata_files"]]

print(gates)

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

named_cts_files <- make_named_version(ids = gates, file_list = cts_files)
named_coldata_files <- make_named_version(
  ids = gates, file_list = coldata_files
)

cts_combined <- data.frame()
coldata_combined <- data.frame()
for (gate_idx in seq_len(length(gates))) {
  cts <- read.csv(named_cts_files[gates[gate_idx]],
    row.names = 1
  )
  coldata <- read.csv(named_coldata_files[gates[gate_idx]],
    row.names = 1
  )
  if (gate_idx == 1) {
    cts_combined <- cts
    coldata_combined <- coldata
  } else {
    cts_combined <- cbind(cts_combined, cts)
    coldata_combined <- rbind(coldata_combined, coldata)
  }
}

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
intgroup_fl <- c("sample", "gate_name")
intgroup_fl_df <- as.data.frame(colData(vsd_fl)[, intgroup_fl, drop = FALSE])
# add grouping if wanted
group <- if (length(intgroup_fl) > 1) {
  factor(apply(intgroup_fl_df, 1, paste, collapse = ":"))
} else {
  colData(vsd_fl)[[intgroup_fl]]
}
d_fl <- data.frame(
  PC1 = pca_fl$x[, 1],
  PC2 = pca_fl$x[, 2],
  PC3 = pca_fl$x[, 3],
  group = group, intgroup_fl_df, name = colnames(vsd_fl)
)
loading_fl <- data.frame(
  LO1 = pca_fl$rotation[, 1],
  LO2 = pca_fl$rotation[, 2],
  LO3 = pca_fl$rotation[, 3],
  name = rownames(pca_fl$rotation)
)

colors_scale <- c(
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
  "#ff5a5a" # MEP135
)

for (pcs in list(1:2, 2:3)) {
  p <- ggplot(
    data = d_fl,
    aes_string(
      x = str_c("PC", pcs[[1]]),
      y = str_c("PC", pcs[[2]])
    )
  ) +
    geom_point(
      aes_string(fill = "gate_name", shape = "sample"),
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
    plot = p,
    width = 7,
    height = 7,
    units = "in"
  )
}
