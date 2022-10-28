library(Seurat)
library(Matrix)
library(dplyr)
library(dict)
library(ggplot2)
library(future)

source("citeseq-functions.R") # needs to be before env.R
source("env.R")

# N.B. use 1.2_runner.R to run this file, i. e. do not run this file directly

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/CS22")

output.dir <- paste0(input.dir, "/", sample)
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

figures.dir <- paste0(output.dir, "/figures")
if (!dir.exists(figures.dir)) {
  dir.create(figures.dir, recursive = TRUE)
}

plan("multicore", workers = availableCores())
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
           filetype = "svg",
           width = w_gg,
           height = h_gg) {
    ggsave(
      filename = paste0(figures.dir, "/", name, ".", filetype),
      plot = plot,
      device = filetype,
      units = "in",
      width = width,
      height = height
    )
  }

# Load data ---------------------------------------------------------------

obj <- readRDS(paste0(input.dir, "/", sample, "_joint_preprocessing.rds"))
obj@project.name <- paste0(obj@project.name, "-", sample)

# RNA-seq ----------------------------------------------------------

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, assay = "RNA", block.size = block.size)

val <- 1500
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = val,
  # nfeatures = 2500,
  assay = "RNA"
)
write.csv(
  VariableFeatures(obj),
  file = paste0(output.dir, "/", val, "-nfeatures.csv"),
  row.names = FALSE,
  col.names = NULL,
  quote = TRUE
)

obj <- ScaleData(obj,
  features = VariableFeatures(obj),
  block.size = block.size,
  assay = "RNA"
)

variable.features.plot(obj)
save_plot("variable-features", width = 12)

obj <-
  CellCycleScoring(obj, cc.genes$s.genes, cc.genes$g2m.genes)
head(obj[[]])

npcs <- 50
if (length(obj$orig.ident) < npcs) {
  npcs <- length(obj$orig.ident) - 1
}

obj <- RunPCA(
  obj,
  assay = "RNA",
  npcs = npcs,
  features = VariableFeatures(obj),
  verbose = F
)
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

ElbowPlot(obj) + labs(
  title = "PCA Elbowplot",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("PCA-elbowplot")
pc <- 9

obj <- clustering.and.UMAP(obj,
  pcs.to.use = pc
)
table(obj[["seurat_clusters"]])

DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.box = T,
  pt.size = pt.size,
  label.size = label.size
  # cols = wes_palette("Darjeeling1",
  #                    n = as.numeric(
  #                      levels(obj@meta.data$seurat_clusters)[max(as.numeric(obj@meta.data$seurat_clusters))]
  #                    ) + 1,
  #                    type = "continuous")
) + labs(
  subtitle = "Seurat clusters",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_UMAP_clusters_pre-cc-removal")

pt.size <- 1
DimPlot(
  obj,
  group.by = "HTO_maxID",
  reduction = "umap",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "HTO",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_UMAP_HTOs_pre-cc-removal")

DimPlot(
  obj,
  group.by = "Phase",
  reduction = "umap",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Cell cycle phase",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_UMAP_cc-phase_pre-cc-removal")

DimPlot(
  obj,
  group.by = "HTO_maxID",
  reduction = "pca",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "HTO",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_PCA_HTOs_pre-cc-removal")

DimPlot(
  obj,
  group.by = "Phase",
  reduction = "pca",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Cell cycle phase",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_PCA_cc-phase_pre-cc-removal")

# Removing cell cycle effects ---------------------------------------------

# saveRDS(obj,
#   file = paste0(output.dir, "/obj_pre-cc-effect-removal.rds")
# )
# obj <-
#   readRDS(paste0(output.dir, "/obj_pre-cc-effect-removal.rds"))

obj <- remove.cell.cycle.effects(obj,
  pcs.to.use = pc,
  features = rownames(obj),
  cluster.resolution = .7,
  block.size = 1e6
  # block.size = block.size
)

saveRDS(obj,
  file = paste0(output.dir, "/", sample, "_post_step_1.2.rds")
)
# obj <-
#   readRDS(paste0(output.dir, "/cite-singlet-post-cc-effect-removal.rds"))


plot <- DimPlot(
  obj,
  reduction = "umap.post.cc",
  group.by = "seurat_clusters",
  label = TRUE,
  label.box = T,
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Seurat clusters",
  caption = paste0("[", obj@project.name, "]")
)
plot
save_plot("RNA_UMAP_clusters_post-cc-removal")
g <- ggplot_build(plot)
cols <- g$data[[1]]
cl_cols <- list()
for (i in levels(obj@meta.data$seurat_clusters)) {
  cl_cols[i] <- unique(cols[cols$group == as.character(as.integer(i) + 1), ]$colour)
}
write.csv(t(t(cl_cols)), file = paste0(output.dir, "/", "cluster_colours.csv"))
# can use this to extract the cluster colours used here, to be able to use them in the bar plots in Nabo.
# Will be helpful for recognizing clusters in bar plot.

DimPlot(
  obj,
  group.by = "Phase",
  reduction = "umap.post.cc",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Cell cycle phase",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_UMAP_cc-phase_post-cc-removal")

DimPlot(
  obj,
  group.by = "Phase",
  reduction = "pca.post.cc",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Cell cycle phase",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_PCA_cc-phase_post-cc-removal")

DimPlot(obj,
  group.by = "Phase",
  reduction = "pca.post.cc.on.cc.genes",
  pt.size = pt.size,
  label.size = label.size
) + labs(
  subtitle = "Cell cycle phase",
  caption = paste0("[", obj@project.name, "]")
)
save_plot("RNA_PCA-on-cc-genes_cc-phase_post-cc-removal")

# Export metadata ---------------------------------------------------------

meta.data <-
  cbind(
    obj@meta.data,
    obj@reductions$umap.post.cc@cell.embeddings
  )
colnames(meta.data)[which(names(meta.data) == "UMAP_1")] <-
  "UMAP.POST.CC_1"
colnames(meta.data)[which(names(meta.data) == "UMAP_2")] <-
  "UMAP.POST.CC_2"
write.csv(meta.data, paste0(output.dir, "/metadata_step_1.2.csv"))
