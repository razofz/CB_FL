# author: "[Rasmus Olofzon](https://github.com/razofz)"

library(Seurat)
library(Matrix)
library(dplyr)
library(dict)
library(ggplot2)
library(tictoc)

# Utility functions -------------------------------------------------------

save.plot <- function(plot, plotname, analysis_dir, width = 1000, height = 600, filetype = "png") {
  if (filetype == "png") {
    png(
      filename = paste0(analysis_dir, "/", plotname, ".png"),
      width = width, height = height
    )
  } else if (filetype == "svg") {
    svg(
      filename = paste0(analysis_dir, "/", plotname, ".svg"),
      width = width, height = height
    )
  }
  print(plot)
  dev.off()
}

write_report <- function(entry, analysis_dir, append = TRUE) {
  cat(entry, file = paste0(analysis_dir, "/analysis-report_", what, ".md"), sep = "\n\n", append = append)
}

save_output_start <- function(analysis_dir) {
  log.file <- file(paste0(analysis_dir, "/analysis-report_", what, ".md"), open = "a+")
  sink(log.file, append = TRUE)
  sink(log.file, append = TRUE, type = "message")
}

save_output_stop <- function() {
  sink()
  sink(type = "message")
}

print_dict <- function(d) {
  for (val in d$items()) {
    cat(paste0(val[2], " \t ", val[1], "\n"))
  }
}

# Load data ---------------------------------------------------------------

get.Seurat.obj.from.data.samplewise <- function(sample,
                                                input.dir,
                                                exact.input.dir = FALSE,
                                                report.data = NULL,
                                                min.cells = 3,
                                                min.features = 200) {
  print(paste0("Loading ", sample))
  if (exact.input.dir) {
    fn <- paste0(input.dir)
  } else {
    fn <- paste0(input.dir, "/", sample)
  }

  citerun.data <- Read10X(data.dir = fn)
  # rownames(x = citerun.data[["Antibody Capture"]]) <- gsub(
  #   pattern = "_[control_]*TotalSeqB", replacement = "",
  #   x = rownames(x = citerun.data[["Antibody Capture"]])
  # )
  # page(citerun.data)

  if (!is.null(report.data)) {
    report.data[[paste0("raw-nbr-cells-", sample)]] <-
      citerun.data$`Gene Expression`@Dim[2]
    report.data[[paste0("raw-nbr-features-", sample)]] <-
      citerun.data$`Gene Expression`@Dim[1]
  }

  loaded.object <- CreateSeuratObject(
    counts = citerun.data[["Gene Expression"]],
    min.cells = min.cells,
    min.features = min.features,
    project = sample
  )
  loaded.object[["ADT"]] <- CreateAssayObject(citerun.data[["Antibody capture"]][, colnames(x = loaded.object)])
  loaded.object[["HTO"]] <- CreateAssayObject(citerun.data[["Antibody capture 2"]][, colnames(x = loaded.object)])

  if (!is.null(report.data)) {
    report.data[[paste0("seurat-nbr-cells-", sample)]] <-
      loaded.object@assays$RNA@counts@Dim[2]
    report.data[[paste0("seurat-nbr-features-", sample)]] <-
      loaded.object@assays$RNA@counts@Dim[1]
  }

  loaded.object[["percent.mt"]] <- PercentageFeatureSet(loaded.object, pattern = "^MT-")

  return(loaded.object)
}

get.Seurat.obj.from.data <- function(project.name,
                                     input.dir,
                                     exact.input.dir = FALSE,
                                     report.data = NULL,
                                     samples = NULL,
                                     min.cells = 3,
                                     min.features = 200) {
  if (exists("samples") && length(samples) > 1) {
    print(paste0("Doing ", project.name, " for samples ", samples))
    loaded.object <- get.Seurat.obj.from.data.samplewise(
      sample = samples[1],
      input.dir = input.dir,
      exact.input.dir = exact.input.dir,
      report.data = report.data,
      min.cells = min.cells,
      min.features = min.features
    )

    for (i in 2:length(samples)) {
      print(c(samples[1:i]))
      loaded.object <- merge(loaded.object,
        y = get.Seurat.obj.from.data.samplewise(
          sample = samples[i],
          input.dir = input.dir,
          report.data = report.data,
          min.cells = min.cells,
          min.features = min.features
        ),
        add.cell.ids = c(samples[1:i]),
        project = project.name
      )
    }
  } else {
    if (exact.input.dir) {
      fn <- paste0(input.dir)
    } else {
      fn <- paste0(input.dir, "/", sample)
    }
    citerun.data <- Read10X(data.dir = fn)
    # rownames(x = citerun.data[["Antibody Capture"]]) <- gsub(
    #   pattern = "_[control_]*TotalSeqB", replacement = "",
    #   x = rownames(x = citerun.data[["Antibody Capture"]])
    # )

    if (!is.null(report.data)) {
      report.data[[paste0("raw-nbr-cells")]] <-
        citerun.data$`Gene Expression`@Dim[2]
      report.data[[paste0("raw-nbr-features")]] <-
        citerun.data$`Gene Expression`@Dim[1]
    }

    loaded.object <- CreateSeuratObject(
      counts = citerun.data[["Gene Expression"]],
      min.cells = min.cells,
      min.features = min.features,
      project = project.name
    )
    loaded.object[["ADT"]] <- CreateAssayObject(citerun.data[["Antibody Capture"]][, colnames(x = loaded.object)])
    loaded.object[["HTO"]] <- CreateAssayObject(citerun.data[["Antibody Capture 2"]][, colnames(x = loaded.object)])

    if (!is.null(report.data)) {
      report.data[["seurat-nbr-cells"]] <-
        loaded.object@assays$RNA@counts@Dim[2]
      report.data[["seurat-nbr-features"]] <-
        loaded.object@assays$RNA@counts@Dim[1]
    }

    loaded.object[["percent.mt"]] <- PercentageFeatureSet(loaded.object, pattern = "^MT-")
  }

  return(loaded.object)
}

# Filtering data ----------------------------------------------------------

filter.data <- function(object,
                        nFeature_RNA_min = 1500,
                        nFeature_RNA_max = 6500,
                        nCount_RNA_min = 0,
                        nCount_RNA_max = 6e4,
                        percent.mt_min = 1,
                        percent.mt_max = 4,
                        pt.size = 0.01) {
  print(nFeature_RNA_min)
  print(nFeature_RNA_max)
  print(nCount_RNA_min)
  print(nCount_RNA_max)
  print(percent.mt_min)
  print(percent.mt_max)




  object <- subset(object,
    subset = nFeature_RNA > nFeature_RNA_min &
      nFeature_RNA < nFeature_RNA_max &
      percent.mt > percent.mt_min &
      percent.mt < percent.mt_max &
      nCount_RNA > nCount_RNA_min &
      nCount_RNA < nCount_RNA_max
  )

  # return(VlnPlot(citerun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = pt.size))
  return(object)
}

# ADT and HTO part --------------------------------------------------------

demultiplex <- function(object) {
  object <- NormalizeData(object, assay = "ADT", normalization.method = "CLR")
  object <- ScaleData(object, assay = "ADT")
  object <- NormalizeData(object, assay = "HTO", normalization.method = "CLR")
  object <- ScaleData(object, assay = "HTO")
  return(HTODemux(object, assay = "HTO"))
}

cells.per.HTO.plot <- function(object, colour.palette = "Paired") {
  cells <- table(object@meta.data[["HTO_maxID"]])

  df <- data.frame(cells)
  p <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), vjust = 1.6, color = "white", size = 6) +
    labs(title = "Number of cells per HTO:", subtitle = paste0("[", object@project.name, "]"))
  p <- p + scale_fill_brewer(palette = colour.palette) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  return(p)
}


# Global classification results
# write_report("")
# nbr.cells <- table(citerun$HTO_classification.global)
# nbr.cells
# for (i in 1:length(nbr.cells)) {
#   write_report(paste0( sapply(dimnames(table(citerun$HTO_classification.global)), "[", i),
#                        's: \t',
#                        table(citerun$HTO_classification.global)[i]))
# }

viz.enrichment.for.HTOs.plot <- function(object, HTOs) {
  Idents(object) <- "HTO_maxID"
  return(RidgePlot(object, assay = "HTO", features = rownames(object[["HTO"]])[1:length(HTOs)], ncol = 2))
}

# Compare number of UMIs for singlets, doublets and negative cells:
starlet.distribution.plot <- function(object, pt.size = 0.5) {
  Idents(object) <- "HTO_classification.global"
  return(
    VlnPlot(object, features = "nCount_RNA", pt.size = pt.size, log = TRUE) +
      labs(caption = paste0("[", citerun@project.name, "]"))
  )
}

# Barplot for doublets, singlets (per HTOs) and negatives
cells.per.HTO.and.starlet.plot <- function(object,
                                           HTOlist,
                                           colour.palette = "Paired") {
  df <- data.frame(summary(object$hash.ID))
  colnames(df) <- c("cell.count")
  df$type <- rownames(df)
  df$type <- factor(df$type, levels = c(HTOlist, "Doublet", "Negative"))
  # df

  p <- ggplot(df, aes(x = type, y = cell.count, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = cell.count), vjust = -.5, color = "black", size = 6)
  p <- p + scale_fill_brewer(palette = colour.palette) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

  p <- p + labs(
    title = "Number of cells per HTO, doublets and negatives:",
    caption = paste0("[", citerun@project.name, "]")
  )
  return(p)
}

# RNA-seq -----------------------------------------------------------------

# Cluster and visualize cells using the usual scRNA-seq workflow,
# and examine for the potential presence of batch effects.
# Extract the singlets, only use these onwards:
take.only.singlets <- function(object) {
  Idents(object) <- "HTO_classification.global"
  return(subset(object, idents = "Singlet"))
}

subset.HTO <- function(object, HTO.list, HTOs.to.use) {
  Idents(object) <- "HTO_maxID"
  for (i in HTO.list[!HTO.list %in% HTOs.to.use]) {
    object <- subset(object, idents = i, invert = TRUE)
  }
  return(object)
}

norm.find.hvgs.scale.data <- function(object, nfeatures = 2000, assay = "RNA",
                                      verbose = TRUE) {
  object <- NormalizeData(object,
    assay = assay,
    verbose = verbose
  )
  object <- FindVariableFeatures(object,
    selection.method = "mean.var.plot",
    nfeatures = nfeatures,
    verbose = verbose
  )
  object <- ScaleData(object,
    features = VariableFeatures(object),
    verbose = verbose
  )
  return(object)
}

# plot variable features with and without labels
variable.features.plot <- function(object) {
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(object), 10)

  plot1 <- VariableFeaturePlot(object)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  return(plot1 + plot2 +
    labs(title = "Variable features", caption = paste0("[", object@project.name, "]")))
}

cells.per.cluster.plot <- function(object, colour.palette = "Paired") {
  cells <- table(object@meta.data[["seurat_clusters"]])

  df <- data.frame(cells)
  p <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), vjust = 1.6, color = "white", size = 6) +
    ggtitle("Number of cells per cluster:")
  # p <- p + scale_fill_manual(values=colour.palette) + theme(axis.title.x=element_blank())
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  return(p)
}

# RNA-seq, continued ------------------------------------------------------

clustering.and.UMAP <- function(object,
                                pcs.to.use,
                                n.neighbours = 30L,
                                reduction.name = "umap",
                                reduction.key = "UMAP_",
                                reduction.to.use = "pca",
                                cluster.resolution = .8,
                                assay = "RNA",
                                verbose = TRUE) {
  object <- FindNeighbors(object,
    reduction = reduction.to.use,
    dims = 1:pcs.to.use,
    features = VariableFeatures(object),
    verbose = verbose
  )
  object <- FindClusters(object,
    resolution = cluster.resolution,
    group.singletons = T,
    verbose = verbose
  )
  object <- RunUMAP(object,
    dims = 1:pcs.to.use,
    n.neighbors = n.neighbours,
    reduction = reduction.to.use,
    assay = assay,
    umap.method = "umap-learn",
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    verbose = verbose
  )
  return(object)
}

# Removing cell cycle effects ---------------------------------------------

remove.cell.cycle.effects <- function(
  object, pcs.to.use, block.size = 1e6, cluster.resolution = .8,
    features = rownames(object)) {
  tic("regressing out cell cycle effects")
  object <- ScaleData(object,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = features,
    block.size = block.size
  )
  toc()

  # Now, a PCA on the variable genes no longer returns components associated with cell cycle
  object <- RunPCA(object,
    features = VariableFeatures(object),
    # nfeatures.print = 10,
    npcs = npcs,
    verbose = F,
    reduction.name = "pca.post.cc"
  )

  # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
  object <- RunPCA(object,
    features = c(cc.genes$s.genes, cc.genes$g2m.genes),
    npcs = npcs,
    reduction.name = "pca.post.cc.on.cc.genes"
  )

  object <- clustering.and.UMAP(object,
    pcs.to.use,
    cluster.resolution = cluster.resolution,
    reduction.to.use = "pca.post.cc",
    reduction.name = "umap.post.cc",
    reduction.key = "UMAP.POST.CC"
  )

  return(object)
}

# Save data ---------------------------------------------------------------


export.metadata <- function(object, output.dir, reduction.to.use = "umap.post.cc") {
  meta.data <- cbind(object@meta.data, object@reductions[[reduction.to.use]]@cell.embeddings)
  write.csv(meta.data, paste0(output.dir, "/metadata.csv"))
}
