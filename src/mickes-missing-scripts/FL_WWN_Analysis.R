# WNN test FL
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


setwd('C:/Users/mikae/Documents/Featal liver project/Bioinformatics/From_Rasmus/CS22-CITE-FL_Seurat-processing_210208-selected')
FL.combined <- readRDS('./FL_combined/FL.combined.rds')

#normalize ADTs and perform PCA
DefaultAssay(FL.combined) <- 'ADT'
VariableFeatures(FL.combined) <- rownames(FL.combined[["ADT"]])
FL.combined <- ScaleData(FL.combined, assay = 'ADT')
FL.combined <- NormalizeData(FL.combined, normalization.method = "CLR", margin = 2)
FL.combined <- RunPCA(FL.combined, reduction.name = 'apca',approx=FALSE)

#Perform WNN
FL.combined <- FindMultiModalNeighbors(
  FL.combined, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

FL.combined[['pca']]@assay.used
FL.combined[['apca']]@assay.used

# Maker WNN UMAP
FL.combined <- RunUMAP(FL.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
FL.combined <- FindClusters(FL.combined, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
# Save WNN analysis
saveRDS(FL.combined, file = "../CS22-CITE-FL_Seurat-processing_210208-selected/FL_combined/FL.combined.WNN.rds")

#reload WNN analysis. Start from this step later
FL.combined <- readRDS('./FL_combined/FL.combined.WNN.rds')


#vizualization of new UMAP, new clusters and old clusters
DimPlot(
  FL.combined,
  reduction = 'umap',
  group.by = 'clust.names',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)
DimPlot(
  FL.combined,
  reduction = 'wnn.umap',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)
DimPlot(
  FL.combined,
  reduction = 'wnn.umap',
  group.by = 'clust.names',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)
DimPlot(
  FL.combined,
  reduction = 'wnn.umap',
  group.by = 'orig.ident',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2,
  shuffle = TRUE
)
VlnPlot(
  FL.combined,
  features = "integrated.weight",
  group.by = 'clust.names',
  sort = TRUE,
  pt.size = 0.1
) +
  NoLegend()
VlnPlot(
  FL.combined,
  features = "ADT.weight",
  group.by = 'clust.names',
  sort = TRUE,
  pt.size = 0.1
) +
  NoLegend()


adt_markers <- FindAllMarkers(FL.combined, assay = "ADT", only.pos = TRUE)
DoHeatmap(FL.combined, features=adt_markers$gene)

RidgePlot(FL.combined, features = c("ADT-CD90", "ADT-CD45RA", "ADT-CD117", "ADT-CD32"), ncol = 2)
VlnPlot(FL.combined, features = unique(adt_markers$gene), 
        group.by = 'wsnn_res.0.5', sort = TRUE, pt.size = 0.1) +
  NoLegend()

FindMarkers(FL.combined, ident.1 = 5, ident.2 = 2, min.pct = 0.25)
FeaturePlot(FL.combined, features = unique(adt_markers$gene))

Orgindet_dis <- as.data.frame.matrix(table(FL.combined@meta.data$orig.ident, FL.combined@meta.data$clust.names))
Orgindet_prop <- NULL
for (i in 1:dim(Orgindet_dis)[1]) {
  a <- Orgindet_dis[i,]/rowSums(Orgindet_dis)[i]*100
  Orgindet_prop <- rbind(Orgindet_prop, a)
}
View(Orgindet_prop)

Comp_WNN_orgClusts <- as.data.frame.matrix(table(FL.combined@meta.data$wsnn_res.0.5, FL.combined@meta.data$clust.names))

OrgClust_prop <- NULL
for (i in 1:dim(Comp_WNN_orgClusts)[1]) {
  a <- Comp_WNN_orgClusts[i,]/colSums(Comp_WNN_orgClusts)*100
  OrgClust_prop <- rbind(OrgClust_prop, a)
}
View(OrgClust_prop)

Comp_orgClusts_WNN <- as.data.frame.matrix(table(FL.combined@meta.data$clust.names, FL.combined@meta.data$wsnn_res.0.5))
WNNClust_prop <- NULL
for (i in 1:dim(Comp_orgClusts_WNN)[1]) {
  a <- Comp_orgClusts_WNN[i,]/colSums(Comp_orgClusts_WNN)*100
  WNNClust_prop <- rbind(WNNClust_prop, a)
}
View(WNNClust_prop)

#exporting UMAP and Metadata
write.table(FL.combined@meta.data, './output/WNN_FL_combined_metadata.csv', sep=',')
write.table(FL.combined@reductions$wnn.umap@cell.embeddings, './output/WNN_FL_combined_umap.csv', sep=',')



