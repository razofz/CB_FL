library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(future)

source("env.R")

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/CS22")
output.dir <- paste0(input.dir, "/cluster_markers")
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

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

for (i in unique(Idents(FL.combined))) {
  markers <- FindMarkers(FL.combined, ident.1 = i)
  hist(markers$avg_log2FC[markers$p_val < 0.05],
    breaks = 20, main =
      paste0(str_replace(i, "/", "-"))
  )
  write.table(markers, paste0(
    output.dir, "/", str_replace(i, "/", "-"),
    "_cluster_markers.csv"
  ), sep = ",")
}
