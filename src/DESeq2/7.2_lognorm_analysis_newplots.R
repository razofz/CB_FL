library(Seurat)
library(ggplot2)
library(patchwork)
# library(SeuratWrappers)

set.seed(snakemake@config[["seed"]])

# roy <- LoadH5Seurat(snakemake@input[["h5seurat"]])
roy <- readRDS(snakemake@input[["roy_rds"]])
# roy <- readRDS('roy.rds')
meta <- read.csv(snakemake@input[["roy_metadata"]])

rownames(meta) <- meta$ids
meta$X <- NULL
meta$ids <- NULL

roy@meta.data <- meta

roy$roy_tissue_type <- factor(
  roy$roy_tissue_type,
  levels = c("eFL", "FL", "FBM", "PBM", "ABM")
)

Idents(roy) <- roy$roy_tissue_type

# DimPlot(roy, reduction = "umap", group.by = "roy_tissue_type") + coord_fixed()
#ggsave(paste0(fig_path, 'umap_sampletype.pdf'))

# DimPlot(roy, reduction = 'umap', group.by = 'roy_cell_type') + coord_fixed()
#ggsave(paste0(fig_path, 'umap_celltype.pdf'), width = 11)

# load genes of interest (GOI)
GOI_FL <- read.table(
  file = snakemake@input[["fetal_signature"]], sep = "\t")

GOI_AD <- read.table(
  file = snakemake@input[["adult_signature"]], sep = "\t")

GOI_AD <- GOI_AD$V1
GOI_FL <- GOI_FL$V1

roy <- NormalizeData(
  roy,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
roy <- ScaleData(roy, features = rownames(roy))

roy <- AddModuleScore(roy, features = list(GOI_FL), name = "FL")
roy <- AddModuleScore(roy, features = list(GOI_AD), name = "AD")

plot_cols <- c("#90CAF9", "#1E88E5", "#0D47A1", "#EF5350", "#7A0E0E")
legend_labs <- c(
  "eFL" = "Early FL",
  "FL" = "FL",
  "FBM" = "Fetal BM",
  "PBM" = "Pediatric BM",
  "ABM" = "Adult BM"
)

VlnPlot(roy,
  features = "FL1",
  split.by = "roy_tissue_type",
  group.by = "roy_tissue_type"
) +
  labs(title = "Fetal signature") +
  ylab("Module score") +
  xlab("") +
  scale_x_discrete(labels = legend_labs) +
  scale_fill_manual(values = plot_cols, labels = legend_labs)
ggsave(snakemake@output[["plot_fl"]], width = 8, height = 8)

VlnPlot(roy,
  features = "AD1",
  split.by = "roy_tissue_type",
  group.by = "roy_tissue_type",
  cols = plot_cols
) +
  labs(title = "Adult signature") +
  ylab("Module score") +
  xlab("") +
  scale_x_discrete(labels = legend_labs) +
  scale_fill_manual(values = plot_cols, labels = legend_labs)
ggsave(snakemake@output[["plot_ad"]], width = 8, height = 8)

score_table <- roy@meta.data[
  ,
  c(
    "roy_sampleID",
    "roy_tissue_type",
    "FL1",
    "AD1"
  )
]

write.csv(score_table, snakemake@output[["scores"]])
