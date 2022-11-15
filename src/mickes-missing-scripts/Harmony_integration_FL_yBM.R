#harmony analysis of FL yBM 

setwd('C:/Users/mikae/Documents/Featal liver project/Bioinformatics/From_Rasmus/CS22-CITE-FL_Seurat-processing_210208-selected')
library(Seurat)
library(harmony)
library(stringr)

FL <- readRDS('./FL_combined/FL.combined.WNN.rds')
yBM <- readRDS('./FL_combined/yBM.rds')

FL_BM.combined <- merge(FL, yBM, add.cell.ids = c('FL', 'BM'), project = 'FL_BM_merge')

DefaultAssay(FL_BM.combined) <- "RNA"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

FL_BM.combined <- CellCycleScoring(FL_BM.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#FL_BM.combined <- ScaleData(FL_BM.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(FL_BM.combined))


FL_BM.combined <- NormalizeData(FL_BM.combined, verbose = F)
FL_BM.combined <- FindVariableFeatures(FL_BM.combined, verbose = F)
FL_BM.combined <- ScaleData(FL_BM.combined, vars.to.regress = c("S.Score", "G2M.Score"), verbose = F) #run if CC regression
#FL_BM.combined <- ScaleData(FL_BM.combined, verbose = F) #run if no CC regression
FL_BM.combined <- RunPCA(FL_BM.combined, verbose = F)

table(FL_BM.combined@meta.data$orig.ident)

#Harmony analysis

FL_BM.combined <- RunHarmony(
  FL_BM.combined,
  group.by.vars = "orig.ident",
)

FL_BM.combined@reductions$harmony_RNA <- FL_BM.combined@reductions$harmony

FL_BM.combined <- FindNeighbors(
  FL_BM.combined,
  reduction = "harmony_RNA",
  dims = 1:15,
  assay = "RNA",
  graph.name = c("NNharmonyRNA", "SNNharmonyRNA")
)
FL_BM.combined <- FindClusters(
  FL_BM.combined,
  graph.name = "SNNharmonyRNA",resolution = 0.8
)
FL_BM.combined[["clusters_RNA_harmony"]] <- Idents(FL_BM.combined)
# FL_BM.combined <- RunUMAP(  #Used for analysis without CC regression
#   FL_BM.combined,
#   reduction = "harmony_RNA",
#   dims = 1:15,
#   reduction.name = "UMAPharmonyRNA",
#   reduction.key = "UMAPharmonyRNA_",
#   nn.name = "NNharmonyRNA",
#   assay = "RNA", n.neighbors = 10, min.dist = 0.005, spread = 6
# )

FL_BM.combined <- RunUMAP( #Used with CC regression
  FL_BM.combined,
  reduction = "harmony_RNA",
  dims = 1:15,
  reduction.name = "UMAPharmonyRNA",
  reduction.key = "UMAPharmonyRNA_",
  nn.name = "NNharmonyRNA",
  assay = "RNA", n.neighbors = 15, min.dist = 0.005, spread = 6
)


DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', group.by = 'orig.ident', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2, shuffle = TRUE)
DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', group.by = 'clusters_RNA_harmony', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2)

#without harmony

FL_BM.combined <- FindNeighbors(
  FL_BM.combined,
  reduction = "pca",
  dims = 1:30,
  assay = "RNA",
  graph.name = c("NNunintegratedRNA", "SNNunintegratedRNA")
)
FL_BM.combined <- FindClusters(
  FL_BM.combined,
  graph.name = "SNNunintegratedRNA"
)
FL_BM.combined[["clusters_RNA_unintegrated"]] <- Idents(FL_BM.combined)
FL_BM.combined <- RunUMAP(
  FL_BM.combined,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "UMAPunintegratedRNA",
  reduction.key = "UMAPunintegratedRNA_",
  nn.name = "NNunintegratedRNA",
  assay = "RNA",
)

#saveRDS(FL_BM.combined, file = "../CS22-CITE-FL_Seurat-processing_210208-selected/FL_combined/FL_BM_combined_CC.rds")
#saveRDS(FL_BM.combined, file = "../CS22-CITE-FL_Seurat-processing_210208-selected/FL_combined/FL_BM_combined.rds")

#if you want only to read again
#FL_BM.combined <- readRDS('./FL_combined/FL_BM_combined.rds') #without CC regression
FL_BM.combined <- readRDS('./FL_combined/FL_BM_combined_CC.rds') #with CC regression

FL_BM.combined@meta.data$combined_clusts <- c(
  paste('FL', FL_BM.combined@meta.data$clust.names[!is.na(FL_BM.combined@meta.data$clust.names)], sep = '_'),
  paste('BM', FL_BM.combined@meta.data$paper_clusters[!is.na(FL_BM.combined@meta.data$paper_clusters)], sep = '_')
)

#plotting harmony
DimPlot(
  FL_BM.combined,
  reduction = 'UMAPharmonyRNA',
  group.by = 'clusters_RNA_harmony',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)
DimPlot(
  FL_BM.combined,
  reduction = 'UMAPharmonyRNA',
  group.by = 'orig.ident',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2,
  shuffle = TRUE
)

DimPlot(
  FL_BM.combined,
  reduction = 'UMAPharmonyRNA',
  group.by = 'combined_clusts',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)

DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', split.by =  'orig.ident' ,group.by = 'clusters_RNA_harmony', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2, ncol = 2)

DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', split.by =  'orig.ident' ,group.by = 'combined_clusts', 
        label = FALSE, repel = TRUE, label.size = 5,pt.size = 2, ncol = 2)

DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', split.by =  'orig.ident' ,group.by = 'clusters_RNA_harmony', 
        label = FALSE, repel = TRUE, label.size = 5,pt.size = 2, ncol = 2)

#plotting unintegrated
DimPlot(
  FL_BM.combined,
  reduction = 'UMAPunintegratedRNA',
  group.by = 'clusters_RNA_unintegrated',
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 2
)
DimPlot(FL_BM.combined, reduction = 'UMAPunintegratedRNA', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 5,pt.size = 2)
DimPlot(FL_BM.combined, reduction = 'UMAPunintegratedRNA', split.by =  'orig.ident' ,group.by = 'clusters_RNA_unintegrated', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2)

#comparing old clusters 
head(FL_BM.combined@meta.data$clust.names)
tail(FL_BM.combined@meta.data$clust.names)
head(FL_BM.combined@meta.data$paper_clusters)
tail(FL_BM.combined@meta.data$paper_clusters)

#compare sample distribution in harmony clusts
Orgindet_dis <- as.data.frame.matrix(table(FL_BM.combined@meta.data$orig.ident, 
                                           FL_BM.combined@meta.data$clusters_RNA_harmony))
Orgindet_prop <- NULL
for (i in 1:dim(Orgindet_dis)[1]) {
  a <- Orgindet_dis[i,]/rowSums(Orgindet_dis)[i]*100
  Orgindet_prop <- rbind(Orgindet_prop, a)
}
View(Orgindet_prop)

#compare FL cluster distribution in harmony clusts
Comp_Har_FLClusts <- as.data.frame.matrix(table(FL_BM.combined@meta.data$clusters_RNA_harmony, 
                                                 FL_BM.combined@meta.data$clust.names))

Comp_Har_FLClusts_prop <- NULL
for (i in 1:dim(Comp_Har_FLClusts)[1]) {
  a <- Comp_Har_FLClusts[i,]/colSums(Comp_Har_FLClusts)*100
  Comp_Har_FLClusts_prop <- rbind(Comp_Har_FLClusts_prop, a)
}
View(Comp_Har_FLClusts_prop)

#compare yBM cluster distribution in harmony clusts
Comp_Har_BMClusts <- as.data.frame.matrix(table(FL_BM.combined@meta.data$clusters_RNA_harmony, 
                                                FL_BM.combined@meta.data$paper_clusters))

Comp_Har_BMClusts_prop <- NULL
for (i in 1:dim(Comp_Har_BMClusts)[1]) {
  a <- Comp_Har_BMClusts[i,]/colSums(Comp_Har_BMClusts)*100
  Comp_Har_BMClusts_prop <- rbind(Comp_Har_BMClusts_prop, a)
}
View(Comp_Har_BMClusts_prop)

#compare FL and yBM cluster distribution in harmony clusts
Comp_Har_FL_BMClusts <- as.data.frame.matrix(table(FL_BM.combined@meta.data$combined_clusts, 
                                                FL_BM.combined@meta.data$clusters_RNA_harmony))

Comp_Har_FL_BMClusts_prop <- NULL
for (i in 1:dim(Comp_Har_FL_BMClusts)[1]) {
  a <- Comp_Har_FL_BMClusts[i,]/colSums(Comp_Har_FL_BMClusts)*100
  Comp_Har_FL_BMClusts_prop <- rbind(Comp_Har_FL_BMClusts_prop, a)
}
View(Comp_Har_FL_BMClusts_prop)


#compare harmony clusts distribution in FL clusters
Comp_FLClusts_Har <- as.data.frame.matrix(table(FL_BM.combined@meta.data$clust.names, 
                                                FL_BM.combined@meta.data$clusters_RNA_harmony))

Comp_FLClusts_Har_prop <- NULL
for (i in 1:dim(Comp_FLClusts_Har)[1]) {
  a <- Comp_FLClusts_Har[i,]/colSums(Comp_FLClusts_Har)*100
  Comp_FLClusts_Har_prop <- rbind(Comp_FLClusts_Har_prop, a)
}
View(Comp_FLClusts_Har_prop)

#compare harmony clusts distribution in yBM cluster
Comp_BMClusts_Har <- as.data.frame.matrix(table(FL_BM.combined@meta.data$paper_clusters, 
                                                FL_BM.combined@meta.data$clusters_RNA_harmony))

Comp_BMClusts_Har_prop <- NULL
for (i in 1:dim(Comp_BMClusts_Har)[1]) {
  a <- Comp_BMClusts_Har[i,]/colSums(Comp_BMClusts_Har)*100
  Comp_BMClusts_Har_prop <- rbind(Comp_BMClusts_Har_prop, a)
}
View(Comp_BMClusts_Har_prop)

#FL_BM.combined <- SetIdent(FL_BM.combined, value = 'clusters_RNA_harmony')
#prop.table(table(Idents(FL_BM.combined), FL_BM.combined$orig.ident),margin = 2)*100

#FL_BM.combined <- SetIdent(FL_BM.combined, value = 'clust.names')
#prop.table(table(Idents(FL_BM.combined), FL_BM.combined$clusters_RNA_harmony),margin = 2)*100


DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', group.by = 'clust.names', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2, 
        cells =Cells(FL_BM.combined)[!is.na(FL_BM.combined@meta.data$clust.names)])

DimPlot(FL_BM.combined, reduction = 'UMAPharmonyRNA', group.by = 'paper_clusters', 
        label = TRUE, repel = TRUE, label.size = 5,pt.size = 2, 
        cells =Cells(FL_BM.combined)[!is.na(FL_BM.combined@meta.data$paper_clusters)])


# Export of data for python analysis
write.table(FL_BM.combined@assays$RNA@var.features, './output/FL_BM_combined_hvgs.csv', sep=',')
write.table(FL_BM.combined@meta.data, './output/FL_BM_combined_metadata.csv', sep=',')
write.table(FL_BM.combined@reductions$UMAPharmonyRNA@cell.embeddings, './output/FL_BM_combined_umap.csv', sep=',')







