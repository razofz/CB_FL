library(Seurat)
library(Matrix)
library(dplyr)
library(dict)
library(ggplot2)
library(tictoc)

source("env.R")

# Global variables --------------------------------------------------------

input.dir <- "data/raw/cite_cell_matrices/FL_CS16_and_W9"

output.dir <- paste0("data/processed/seurat/W9")
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

# figures.dir <- paste0(output.dir, "/figures")
# if (!dir.exists(figures.dir)) {
#   dir.create(figures.dir, recursive = TRUE)
# }

seed <- 1234
set.seed(seed)

# Utility functions -------------------------------------------------------

# w_gg <- 8.3
# h_gg <- 6

# save.plot <- function(name, plot = last_plot(), filetype = "svg", width = w_gg,
#                       height = h_gg) {
#   ggsave(
#     filename = paste0(figures_dir, "/", name, ".", filetype),
#     plot = plot,
#     device = filetype,
#     units = "in",
#     width = width,
#     height = height
#   )
# }

# Load data ---------------------------------------------------------------

min.cells <- 3
min.features <- 200

load.data <- function(project.name) {
  citerun.data <- Read10X(data.dir = input.dir)
  rownames(x = citerun.data[["Antibody capture 2"]]) <- gsub(
    pattern = "_[control_]*TotalSeqB", replacement = "",
    x = rownames(x = citerun.data[["Antibody capture 2"]])
  )

  loaded.object <- CreateSeuratObject(
    counts = citerun.data[["Gene Expression"]],
    min.cells = min.cells,
    min.features = min.features,
    project = project.name
  )
  loaded.object[["ADT"]] <- CreateAssayObject(citerun.data[["Antibody capture"]][, colnames(x = loaded.object)])
  loaded.object[["HTO"]] <- CreateAssayObject(citerun.data[["Antibody capture 2"]][, colnames(x = loaded.object)])

  return(loaded.object)
}

citerun <- load.data("W9")

# Filtering data ----------------------------------------------------------

citerun[["percent.mt"]] <- PercentageFeatureSet(citerun, pattern = "^MT-")

filtering_params <- dict()
filtering_params[["nFeature_RNA_min"]] <- 1500
filtering_params[["nFeature_RNA_max"]] <- 6500
filtering_params[["nCount_RNA_min"]] <- 0
filtering_params[["nCount_RNA_max"]] <- 6e4
filtering_params[["percent.mt_min"]] <- 1
filtering_params[["percent.mt_max"]] <- 4

citerun <- subset(citerun,
  subset = nFeature_RNA > filtering_params[["nFeature_RNA_min"]] &
    nFeature_RNA < filtering_params[["nFeature_RNA_max"]] &
    percent.mt > filtering_params[["percent.mt_min"]] &
    percent.mt < filtering_params[["percent.mt_max"]] &
    nCount_RNA > filtering_params[["nCount_RNA_min"]] &
    nCount_RNA < filtering_params[["nCount_RNA_max"]]
)

# ADT and HTO part --------------------------------------------------------

citerun <- NormalizeData(citerun, assay = "ADT", normalization.method = "CLR")
citerun <- ScaleData(citerun, assay = "ADT")
citerun <- NormalizeData(citerun, assay = "HTO", normalization.method = "CLR")
citerun <- ScaleData(citerun, assay = "HTO")

citerun <- HTODemux(citerun, assay = "HTO")
saveRDS(citerun, file = paste0(output.dir, "/FL_CS16_and_W9_subset.rds"))
