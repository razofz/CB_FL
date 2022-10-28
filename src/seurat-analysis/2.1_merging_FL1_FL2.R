library(Seurat)
library(dplyr)
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

FL1 <- readRDS(paste0(input.dir, "/FL1_hpc/FL1_hpc_post_step_1.3.rds"))
FL2 <- readRDS(paste0(input.dir, "/FL2_hpc/FL2_hpc_post_step_1.3.rds"))

# integrate RNA
selfeats <- SelectIntegrationFeatures(c(FL1, FL2), nfeatures = 2000, assay = c("RNA", "RNA"))
FL.anchors <- FindIntegrationAnchors(c(FL1, FL2), anchor.features = selfeats)

FL.combined <- IntegrateData(FL.anchors)
DefaultAssay(FL.combined) <- "integrated"

FL.combined <- ScaleData(FL.combined, verbose = FALSE)
FL.combined <- ScaleData(FL.combined, verbose = FALSE, assay = "RNA")
FL.combined <- RunPCA(FL.combined, npcs = 30, verbose = FALSE)
ElbowPlot(FL.combined)
FL.combined <- RunUMAP(FL.combined, reduction = "pca", dims = 1:14)
DimPlot(FL.combined, reduction = "umap", pt.size = 2, label = TRUE, label.size = 12)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

FL.combined <- CellCycleScoring(FL.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
FL.combined <- ScaleData(FL.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(FL.combined))
FL.combined <- RunPCA(FL.combined, npcs = 30, verbose = FALSE)

ElbowPlot(FL.combined)
FL.combined <- RunUMAP(FL.combined, reduction = "pca", dims = 1:14)
FL.combined@reductions$umap@cell.embeddings[,'UMAP_2'] <- -Embeddings(FL.combined, reduction='umap')[,'UMAP_2'] # flip UMAP to be similar to sample-specific umaps
DimPlot(FL.combined, reduction = "umap", pt.size = 2, label = TRUE, label.size = 12)

FL.combined <- FindNeighbors(FL.combined, reduction = "pca", dims = 1:30)
FL.combined <- FindClusters(FL.combined, resolution = 0.7)

DimPlot(FL.combined, reduction = "umap", pt.size = 2, label = TRUE, label.size = 12)

Idents(FL.combined) <- "orig.ident"
DimPlot(FL.combined, reduction = "umap", pt.size = 2, label = TRUE, label.size = 8)

Idents(FL.combined) <- "integrated_snn_res.0.7"
FL.markers <- FindAllMarkers(FL.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
genes.var <- FL.combined@assays$integrated@var.features

top10 <- FL.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(FL.combined, features = top10$gene, assay = "integrated")
DotPlot(FL.combined, features = unique(top10$gene), assay = "integrated") + RotatedAxis()
DoHeatmap(FL.combined, features = top10$gene, assay = "RNA")
DotPlot(FL.combined, features = unique(top10$gene), assay = "RNA") + RotatedAxis()

# 0 MPP-II
# 1 MPP-I
# 2 Ly-I
# 3 MEP
# 4 GMP
# 5 Ly-II
# 6 HSC
# 7 DC-Mono
# 8 Cyc
# 9 DC-I
# 10 Ly-III (T)


new.cluster.ids <- c("MPP-II", "MPP-I", "Ly-I", "MEP", "GMP", "Ly-II", "HSC",
                     "DC-Mono", "Cyc", "DC-I", "Ly-III")
names(new.cluster.ids) <- levels(FL.combined)
FL.combined <- RenameIdents(FL.combined, new.cluster.ids)
# levels(FL.combined) <- c("Ly-III", "DC-I", "Cyc", "MEP", "DC-Mono", "GMP", "Ly-II", "Ly-I", "MPP-II", "MPP-I", "HSC")
FL.combined[["clust.names"]] <- Idents(object = FL.combined)

DefaultAssay(FL.combined) <- "ADT"
FL.combined <- ScaleData(FL.combined, assay = "ADT")
FL.combined <- NormalizeData(FL.combined, normalization.method = "CLR", margin = 2)

# Save data ---------------------------------------------------------------

saveRDS(FL.combined, file = paste0(output.dir, "/FL.combined.rds"))

write.table(FL1@assays$ADT@scale.data, paste0(output.dir, "/FL1_adt_scaled_values.csv"), sep = ",")
write.table(FL2@assays$ADT@scale.data, paste0(output.dir, "/FL2_adt_scaled_values.csv"), sep = ",")

write.table(genes.var, paste0(output.dir, "/FL_combined_hvgs.csv"), sep = ",")
write.table(FL.combined@meta.data, paste0(output.dir, "/FL_combined_metadata.csv"), sep = ",")
write.table(FL.combined@reductions$umap@cell.embeddings, paste0(output.dir, "/FL_combined_umap.csv"), sep = ",")
write.table(FL.markers, paste0(output.dir, "/FL_combined_cluster_markers.csv"), sep = ",")
