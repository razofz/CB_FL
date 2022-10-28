library(Seurat)
library(Matrix)
library(dplyr)
library(dict)
library(ggplot2)
library(future)

source("citeseq-functions.R") # needs to be before env.R
source("env.R")

# Global variables --------------------------------------------------------

samples <- c("FL1_hpc", "FL2_hpc")

input_dir <- paste0("data/raw/cite_cell_matrices")
output_dir <- paste0("data/processed/seurat/CS22")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
figures_dir <- paste0(output_dir, "/figures")
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

plan("multicore", workers = 16)
plan()

seed <- 1234
set.seed(seed)

# Helper functions --------------------------------------------------------

w_gg <- 8.3
h_gg <- 6

save_plot <- function(name, plot = last_plot(), filetype = "svg", width = w_gg,
                      height = h_gg) {
  ggsave(
    filename = paste0(figures_dir, "/", name, ".", filetype),
    plot = plot,
    device = filetype,
    units = "in",
    width = width,
    height = height
  )
}

# Load data ---------------------------------------------------------------

min_cells <- 3
min_features <- 200

citerun <-
  get.Seurat.obj.from.data(
    project.name = paste0("CITE-CS22"),
    input.dir = input_dir,
    samples = samples,
    min.cells = min_cells,
    min.features = min_features
  )

# Filter data ----------------------------------------------------------

VlnPlot(citerun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
save_plot("pre_filtering")

plot1 <- FeatureScatter(citerun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(citerun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
save_plot("feature-scatter-pre-filtering")

filtering_params <- dict()
filtering_params[["nFeature_RNA_min"]] <- 1000
filtering_params[["nFeature_RNA_max"]] <- 5000
filtering_params[["nCount_RNA_min"]] <- 0
filtering_params[["nCount_RNA_max"]] <- 2.5e4
filtering_params[["percent.mt_min"]] <- 1
filtering_params[["percent.mt_max"]] <- 3.5

citerun <- subset(citerun,
  subset = nFeature_RNA > filtering_params[["nFeature_RNA_min"]] &
    nFeature_RNA < filtering_params[["nFeature_RNA_max"]] &
    percent.mt > filtering_params[["percent.mt_min"]] &
    percent.mt < filtering_params[["percent.mt_max"]] &
    nCount_RNA > filtering_params[["nCount_RNA_min"]] &
    nCount_RNA < filtering_params[["nCount_RNA_max"]]
)
plot1 <- FeatureScatter(citerun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(citerun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
save_plot("feature-scatter-post-filtering")

VlnPlot(citerun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
save_plot("post_filtering")

# ADT and HTO part --------------------------------------------------------

# Normalizes and scales ADT and HTO data, as well as demultiplexes
citerun <- demultiplex(citerun)

cells_per_HTO <- table(citerun@meta.data[["HTO_maxID"]])
cells_per_HTO
HTOs <- names(cells_per_HTO)

colour.palette <- "Accent"

cells.per.HTO.plot(citerun, colour.palette = colour.palette)
save_plot("cells-per-HTO-only")

# Global classification results
nbr_cells <- table(citerun$HTO_classification.global)
nbr_cells

Idents(citerun) <- "HTO_maxID"

RidgePlot(citerun, assay = "HTO", features = rownames(citerun[["HTO"]])[1:4], ncol = 2) +
  labs(caption = paste0("[", citerun@project.name, "]"))
save_plot("enrichment-HTOs")

# Compare number of UMIs for singlets, doublets and negative cells:
starlet.distribution.plot(citerun)
save_plot("distribution-singlet-doublet-negative")

# Barplot for doublets, singlets (per HTOs) and negatives
cells.per.HTO.and.starlet.plot(citerun, HTOs, colour.palette = colour.palette)
save_plot("cells-per-HTO-doublet-negative", width = 10, height = 8)

# Subsetting, removing negatives
Idents(citerun) <- "HTO_classification.global"
citerun.subset <- subset(citerun, idents = "Negative", invert = TRUE)
Idents(citerun.subset) <- "HTO_maxID"

# Calculate a distance matrix using HTOs
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = citerun.subset, assay = "HTO"))))
HTOHeatmap(citerun.subset, assay = "HTO")
save_plot("HTOHeatmap")

# Extract the singlets, only use these onwards:
citerun.singlet <- take.only.singlets(citerun)

HTOs

citerun.FL1 <- subset.HTO(citerun.singlet, HTOs, HTOs[2])
citerun.FL2 <- subset.HTO(citerun.singlet, HTOs, HTOs[3])

# subset on the two samples and save RDS's, then handle them in separate files
# after that

# Save data ---------------------------------------------------------------

saveRDS(citerun.FL1, paste0(
  output_dir, "/", samples[[1]],
  "_joint_preprocessing.rds"
))
saveRDS(citerun.FL2, paste0(
  output_dir, "/", samples[[2]],
  "_joint_preprocessing.rds"
))
