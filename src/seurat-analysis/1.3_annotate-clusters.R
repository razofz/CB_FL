library(Seurat)
library(ggplot2)
library(future)

# N.B. use 1.3_runner.R to run this file, i. e. do not run this file directly

source("citeseq-functions.R") # needs to be before env.R
source("env.R")

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/CS22/", sample)

output.dir <- paste0(input.dir)
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

figures.dir <- paste0(output.dir, "/figures")
if (!dir.exists(figures.dir)) {
  dir.create(figures.dir, recursive = TRUE)
}
raw.dir <- paste0("data/raw/paper_specific/from_seurat")

plan("multicore", workers = 16)
plan()

# Seurat speedup variables
block.size <- 1e6

# Plotting variables
pt.size <- 1.5
label.size <- 5
w_gg <- 8.3
h_gg <- 6

seed <- 1234
set.seed(seed)

options(future.globals.maxSize = 5000 * 1024 ^ 2)
options(future.seed = TRUE)

# Helper functions --------------------------------------------------------

save_plot <-
  function(name,
           plot = last_plot(),
           filetype = 'svg',
           width = w_gg,
           height = h_gg) {
    ggsave(
      filename = paste0(figures.dir, '/', name, '.', filetype),
      plot = plot,
      device = filetype,
      units = 'in',
      width = width,
      height = height
    )
  }

# Load data ---------------------------------------------------------------

obj <- readRDS(paste0(input.dir, "/", sample, "_post_step_1.2.rds"))

# Change default umap -----------------------------------------------------

obj@reductions$umap.pre.cc <- obj@reductions$umap
obj@reductions$umap <- obj@reductions$umap.post.cc
obj@reductions$umap@key <- 'UMAP_'
colnames(obj@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")

# Annotate clusters -------------------------------------------------------

obj[['orig_seurat_clusters']] <- obj[['seurat_clusters']]

if (sample == "FL2_hpc") {
  print(table(obj$seurat_clusters))
  # note: R indexing from 1, seurat cluster levels from 0
  levels(obj$seurat_clusters)[7] <- 1
  print(table(obj$seurat_clusters))

  # merge 6 and 1
  # 1 hsc
  # 4 mpp1
  # 5 mpp2
  # 7 cyc 8
  # 2 ery
  # 3 my
  # 0 mpp3
  # 6 ly 7

  obj$clusters_annotated <- obj$seurat_clusters
  print(levels(obj$clusters_annotated))
  levels(obj$clusters_annotated) <- c(
    "MPP2",
    "HSC",
    "Ery",
    "My",
    "MPP1",
    "MPP3",
    "Ly",
    "Cyc"
  )
  print(levels(obj$clusters_annotated))
  print(table(obj$clusters_annotated))
} else if (sample == "FL1_hpc") {
  print(table(obj$seurat_clusters))
  levels(obj$seurat_clusters)[9] <- 0
  levels(obj$seurat_clusters)[8] <- 5
  print(table(obj$seurat_clusters))

  # merge 7 and 5 to 5
  # merge 8 and 0 to 0
  # 3 hsc
  # 4 mpp1
  # 1 mpp2
  # 2 ery
  # 5 my
  # 0 mpp3
  # 6 ly

  obj$clusters_annotated <- obj$seurat_clusters
  print(levels(obj$clusters_annotated))
  levels(obj$clusters_annotated) <- c(
    "MPP3",
    "MPP2",
    "Ery",
    "HSC",
    "MPP1",
    "My",
    "Ly"
    # "Cyc"
  )
  levels(obj$clusters_annotated)
  table(obj$clusters_annotated)
}

# Visualise annotations ---------------------------------------------------

FL2_cols <- read.csv(paste0(raw.dir, '/', 'cluster_colours_', samples[2], '.csv'), row.names = 'X')

FL2_cols <- FL2_cols[order(match(FL2_cols[[1]], levels(obj$clusters_annotated))), ]
print(FL2_cols)

cl_cols <- FL2_cols$colour

DimPlot(
  obj,
  reduction = "umap",
  group.by = "orig_seurat_clusters",
  label = TRUE,
  label.box = T,
  repel = T,
  pt.size = pt.size,
  label.size = label.size
) + labs(subtitle = 'Seurat clusters',
         caption = paste0('[', obj@project.name, ']'))
save_plot("RNA_UMAP_Seurat-clusters")

DimPlot(
  obj,
  reduction = "umap",
  group.by = "clusters_annotated",
  label = TRUE,
  label.box = T,
  repel = T,
  cols = cl_cols,
  pt.size = pt.size,
  label.size = label.size
) + labs(subtitle = 'Annotated clusters',
         caption = paste0('[', obj@project.name, ']'))
save_plot("RNA_UMAP_annotated-clusters")

# Save data ---------------------------------------------------------

meta.data <-
  cbind(obj@meta.data,
        obj@reductions$umap@cell.embeddings)
write.csv(meta.data, paste0(output.dir, '/metadata_step_1.3.csv'))

saveRDS(obj, file = paste0(output.dir, '/', sample, '_post_step_1.3.rds'))
