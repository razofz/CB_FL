# N. B.: scripts 2.1, 2.2 and 2.3 must all be run in the same namespace.
# E. g. start an R console and `source()` the scripts one after another. - raz

source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))

out_dir <- "data/processed/DEseq2"
path_to_files <- "data/raw/paper_specific/DEseq2/gate_wise"
# path_to_files <- "data/processed/notebooks/DEseq2/gate_wise/"
images_dir <- paste0(out_dir, "/images")
if (dir.exists(images_dir) == FALSE) {
  dir.create(images_dir, showWarnings = TRUE, recursive = TRUE)
}

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
  pattern = "yBM_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "yBM_", x =
    rownames(coldata_combined)
), ]
cts_combined_gate <- cts_combined
coldata_combined_gate <- coldata_combined

dds <- DESeqDataSetFromMatrix(
  countData = cts_combined, colData = coldata_combined,
  design = ~1
)
dds <- DESeq(dds)
vsd_FL_gate <- vst(dds, blind = TRUE)


rv_FL <- rowVars(assay(vsd_FL_gate)) # calculate variance
# find the top 500 genes by variance
select_FL <- order(rv_FL, decreasing = T)[seq_len(min(500, length(rv_FL)))]
pca_FL <- prcomp(t(assay(vsd_FL_gate)[select_FL, ])) # calculate PCA
percentVar_FL <- pca_FL$sdev^2 / sum(pca_FL$sdev^2) # calc PCA contribution
intgroup_FL <- c("sample", "gate_name")
# add grouping if wanted
intgroup_FL.df <- as.data.frame(colData(vsd_FL_gate)[, intgroup_FL, drop = FALSE])
group <- if (length(intgroup_FL) > 1) {
  factor(apply(intgroup_FL.df, 1, paste, collapse = ":"))
} else {
  colData(vsd_FL)[[intgroup_FL]]
}
d_FL <- data.frame(
  PC1 = pca_FL$x[, 1], PC2 = pca_FL$x[, 2], PC3 = pca_FL$x[
    ,
    3
  ], group = group, intgroup_FL.df, name =
    colnames(vsd_FL_gate)
)
loading_FL <- data.frame(
  LO1 = pca_FL$rotation[, 1], LO2 = pca_FL$rotation[
    ,
    2
  ], LO3 = pca_FL$rotation[, 3], name =
    rownames(pca_FL$rotation)
)

colors_scale2 <- c(
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

ggplot(data = d_FL, aes_string(x = "PC1", y = "PC2", shape = "sample")) +
  geom_point(aes_string(fill = "gate_name"),
    color = "black", shape = 21, size = 5
  ) +
  xlab(paste0("PC1: ", round(percentVar_FL[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  theme_bw() +
  # geom_text(aes(label=gate_name), size=3, nudge_x = 1, nudge_y = 1)+
  theme(aspect.ratio = 1) +
  coord_fixed() +
  scale_fill_manual(values = colors_scale2) #, guide = "none")
ggsave(paste0(images_dir, "/PC12_singel_cell_gates_top500.pdf"),
  plot = last_plot(),
  width = 7,
  height = 7,
  units = "in"
)

ggplot(data = d_FL, aes_string(x = "PC2", y = "PC3", shape = "sample")) +
  geom_point(aes_string(fill = "gate_name"),
    color = "black", shape = 21, size =
      5
  ) +
  xlab(paste0("PC2: ", round(percentVar_FL[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar_FL[3] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  # geom_text(aes(label=gate_name), size=3, nudge_x = 0.5, nudge_y = 0.5)+
  coord_fixed() +
  scale_fill_manual(values = colors_scale2) #, guide = "none")
ggsave(paste0(images_dir, "/PC23_singel_cell_gates_top500.pdf"),
  plot = last_plot(),
  width = 7,
  height = 7,
  units = "in"
)
