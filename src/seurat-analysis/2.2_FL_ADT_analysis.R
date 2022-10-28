library(Seurat)
library(ggplot2)
library(future)

source("env.R")

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/CS22")

output.dir <- input.dir
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

figures.dir <- paste0(output.dir, "/figures")
if (!dir.exists(figures.dir)) {
  dir.create(figures.dir, recursive = TRUE)
}

raw.dir <- paste0("data/raw/from_seurat")

# plan("multicore", workers = 16)
plan("sequential")
plan()

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

FL.combined <- readRDS(paste0(input.dir, "/FL.combined.rds"))

adt_markers <- FindAllMarkers(FL.combined, assay = "ADT")

# markers here (rownames) for this to work
feats <- c(
  "adt_ADT-CD90", "adt_ADT-CD38", "adt_ADT-CD45RA", "adt_ADT-CD10", "adt_ADT-CD71",
  "adt_ADT-CD123", "adt_ADT-CD135", "adt_ADT-IL7RA", "adt_ADT-CD201"
)
FeaturePlot(FL.combined, features = feats, pt.size = 1, cols = c("deepskyblue1", "orangered2"), order = TRUE)
FeaturePlot(FL.combined,
  features = unique(adt_markers$gene)[1:16], pt.size = 1,
  cols = c("deepskyblue1", "orangered2"), order = TRUE
)
FeaturePlot(FL.combined,
  features = unique(adt_markers$gene)[17:32], pt.size = 1,
  cols = c("deepskyblue1", "orangered2"), order = TRUE
)
DimPlot(FL.combined, reduction = "umap", pt.size = 2, label = TRUE, label.size = 12)
