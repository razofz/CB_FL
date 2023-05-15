library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

source("env.R")

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/CS22")
output.dir <- paste0(input.dir, "/")
plot.dir <- paste0(input.dir, "/figures/")
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

seed <- 1234
set.seed(seed)

# Seurat speedup variables
block.size <- 1e6

# Plotting variables
pt.size <- 1.5
label.size <- 5
w_gg <- 8.3
h_gg <- 6

options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed = TRUE)


# Load data ---------------------------------------------------------------

sobj <- readRDS(paste0(input.dir, "/FL.combined.rds"))

DefaultAssay(sobj) <- "RNA"

colours_features <- c("moccasin", "darkslategray")
(DimPlot(sobj) + coord_fixed()) %>% ggsave(.,
  filename = str_c(plot.dir, "cc_umap_clusters.svg"),
  device = "svg", units = "in", width = 7,
  height = 7
)
(FeaturePlot(sobj,
  features = "G2M.Score", coord.fixed = T, cols = colours_features
)) %>% ggsave(.,
  filename = str_c(plot.dir, "cc_umap_G2M.svg"),
  device = "svg", units = "in", width = 7, height = 7
)
(FeaturePlot(sobj,
  features = "S.Score", coord.fixed = T, cols = colours_features
)) %>% ggsave(.,
  filename = str_c(plot.dir, "cc_umap_S.svg"),
  device = "svg", units = "in", width = 7, height = 7
)
(VlnPlot(sobj, features = "S.Score")) %>% ggsave(.,
  filename = str_c(plot.dir, "cc_violin_S.svg"), device = "svg",
  units = "in", width = 7, height = 7
)
(VlnPlot(sobj, features = "G2M.Score")) %>% ggsave(.,
  filename = str_c(plot.dir, "cc_violin_G2M.svg"),
  device = "svg", unit = "in", width = 7, height = 7
)

markers <- FindMarkers(sobj,
  ident.1 = "MPP-I", ident.2 = "HSC",
  only.pos = F
) %>% arrange(desc(avg_log2FC))
markers %>% write.csv(., file = str_c(output.dir, "cc_mpp1_v_hsc.csv"))

markers <- FindMarkers(sobj,
  ident.1 = "MPP-I", ident.2 = "MPP-II",
  only.pos = F
) %>% arrange(desc(avg_log2FC))
markers %>% write.csv(., file = str_c(output.dir, "cc_mpp1_v_mpp2.csv"))

clusters <- unique(sobj[["clust.names"]])[["clust.names"]] %>% levels()

counts <- list()
counts <- table(sobj$clust.names, sobj$Phase, sobj$orig.ident)

for (donor in c("FL1_hpc", "FL2_hpc")) {
  write.csv(counts[, , donor], file = str_c(
    output.dir, "cc_phase_counts_donor_", donor, ".csv"
  ))
}
