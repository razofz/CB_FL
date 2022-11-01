invisible(lapply(list(
  "stringr",
  "future.apply",
  "ggplot2",
  "jsonlite",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

plan(multisession, workers = 16)
options(future.globals.maxSize = 20e3 * 1024^2)
set.seed(12345)

out_dir <- str_c(
  Sys.getenv("PROJECT_PATH"),
  "/data/processed/DEseq2/"
)
plots_dir <- str_c(out_dir, "/images/")

donors <- c(str_c("yBM", 1:2), str_c("FL", 1:2))
samples <- c("FL", "BM")

metadata <- read.table(
  str_c(
    Sys.getenv("PROJECT_PATH"),
    "/data/raw/paper_specific/from_seurat/FL_BM_combined_metadata.csv"
  ),
  sep = ","
)
rownames(metadata) <- str_replace(
  rownames(metadata),
  pattern = "hpc_",
  replacement = ""
)
metadata$orig.ident <- str_replace(
  metadata$orig.ident,
  pattern = "_hpc",
  replacement = ""
)

for (sample in samples) {
  if (str_detect(pattern = "BM", string = sample)) {
    formatted_sample <- str_to_lower(sample)
  } else if (str_detect(pattern = "FL", string = sample)) {
    formatted_sample <- sample
  }
  write.table(
    rownames(metadata[metadata$orig.ident == formatted_sample, ]),
    file = str_c(out_dir, sample, "_cells.csv")
  )
}

object_list <- future_lapply(samples,
  FUN = function(sample) {
    mat <- Read10X(str_c(
      Sys.getenv("PROJECT_PATH"),
      "/data/raw/cite_cell_matrices/",
      sample,
      "_hpc/"
    ))$`Gene Expression`
    if (str_detect(pattern = "BM", string = sample)) {
      colnames_mat <- str_c(
        "BM_", str_to_lower(sample), "_",
        colnames(mat)
      )
    } else if (str_detect(pattern = "FL", string = sample)) {
      colnames_mat <- str_c(
        "FL_", sample, "_",
        colnames(mat)
      )
    }
    colnames(mat) <- colnames_mat
    colnames(mat) <- str_replace(
      colnames(mat),
      pattern = "hpc_",
      replacement = ""
    )
    obj <- CreateSeuratObject(
      counts = mat,
      project = sample,
      assay = "RNA"
    )
    cells_to_keep <- read.table(str_c(out_dir, sample, "_cells.csv"))$x
    obj <- subset(obj, cells = cells_to_keep)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(
      obj,
      selection.method = "vst",
      nfeatures = 2000
    )
    obj[["orig.sample"]] <- obj[["orig.ident"]]
    return(obj)
  },
  future.seed = TRUE
)

features <- SelectIntegrationFeatures(object.list = object_list)
anchors <- FindIntegrationAnchors(
  object.list = object_list, anchor.features = features
)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
# combined <- FindClusters(combined, resolution = 0.5)

combined[["donor"]] <- sapply(rownames(combined[[]]), FUN = function(x) {
  unlist(str_split(x, pattern = "_"))[2]
})

# combined[["sample"]] <- lapply(combined[["orig.sample"]], FUN = function(x) {
#   str_sub(x, end = -2)
# })

Idents(combined) <- "donor"
p <- DimPlot(combined, reduction = "umap") + coord_fixed()
p
ggsave(p,
  filename = str_c(plots_dir, "umap.svg"),
  device = "svg", width = 7, height = 7
)

p <- DimPlot(combined, reduction = "umap", split.by = "orig.sample") +
  coord_fixed()
p
ggsave(p,
  filename = str_c(plots_dir, "umap_split.svg"),
  device = "svg", width = 11, height = 5
)

# bm both boys
# fl one boy one girl

path <- str_c(
  Sys.getenv("PROJECT_PATH"),
  "/data/processed/notebooks/mapping_predictions/RNA_FL_target_preds.json"
)

fl_metadata <- read.table(
  str_c(
    Sys.getenv("PROJECT_PATH"),
    "/data/processed/seurat/CS22/FL_combined_metadata.csv"
  ),
  sep = ","
)

mapping_predictions <- unlist(fromJSON(txt = path)$yBM_HPC)
head(mapping_predictions)
new_names <- str_replace(names(mapping_predictions),
  pattern = "__", replacement = "_"
)
new_names <- str_c("BM_", new_names)
names(mapping_predictions) <- new_names
head(mapping_predictions)

new_names <- str_replace(rownames(fl_metadata),
  pattern = "hpc_", replacement = ""
)
new_names <- str_c("FL_", new_names)
rownames(fl_metadata) <- new_names

md <- combined[[]]
md["fl_cluster"] <- "NA"
md[names(mapping_predictions), "fl_cluster"] <- mapping_predictions

fl_metadata <- fl_metadata[rownames(md[md$orig.sample == "FL", ]), ]
md[rownames(fl_metadata), "fl_cluster"] <- fl_metadata$clust.names

combined@meta.data <- md

DimPlot(combined,
  reduction = "umap", group.by = "fl_cluster",
  split.by = "orig.sample"
) + coord_fixed()

DimPlot(combined,
  reduction = "umap", group.by = "fl_cluster",
  split.by = "donor"
) + coord_fixed()

DefaultAssay(combined) <- "RNA"
combined[["comparees"]] <- str_c(
  combined[[]]$"orig.sample", "_",
  combined[[]]$"fl_cluster"
)
Idents(combined) <- "comparees"

fl_clusters <- unique(combined[[]][, "fl_cluster"])
fl_clusters <- fl_clusters[fl_clusters != "NA"]
fl_clusters <- fl_clusters[fl_clusters != "Ly-III"]
fl_clusters <- fl_clusters[fl_clusters != "DC-I"]

# for (cluster in fl_clusters) {
#   markers <- FindConservedMarkers(combined,
#     ident.1 = cluster,
#     ident.2 = cluster,
#     grouping.var = "donor", verbose = FALSE
#   )
#   print(head(markers))
#   # write.table(markers, file = str_c(out_dir, cluster, "_markers.csv"))
# }

ordering <- c(
  "FL_HSC", "BM_HSC", "FL_MPP-I", "BM_MPP-I", "FL_MPP-II",
  "BM_MPP-II", "FL_Ly-I", "BM_Ly-I", "FL_Ly-II", "BM_Ly-II",
  "FL_Ly-III", "BM_Ly-III", "FL_DC-I", "BM_DC-I", "FL_DC-Mono",
  "BM_DC-Mono", "FL_GMP", "BM_GMP", "FL_MEP", "BM_MEP", "FL_Cyc",
  "BM_Cyc", "BM_NA"
)
tmp <- as.factor(combined[[]][, "comparees"])
levels(tmp) <- ordering
combined[["comparees"]] <- tmp

simple_dir <- str_c(out_dir, "markers_simple_way/")
dir.create(simple_dir, recursive = T)

DefaultAssay(combined) <- "RNA"
Idents(combined) <- "comparees"
for (cluster in fl_clusters) {
# for (cluster in c("HSC")) {
  print(cluster)
  fc_threshold <- .25
  markers <- FindMarkers(combined,
    ident.1 = str_c("BM_", cluster),
    ident.2 = str_c("FL_", cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  write.table(markers, file = str_c(
    simple_dir, cluster, "_BM_specific",
    "_markers.csv"
  ))
  markers <- FindMarkers(combined,
    ident.1 = str_c("FL_", cluster),
    ident.2 = str_c("BM_", cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  write.table(markers, file = str_c(
    simple_dir, cluster, "_FL_specific",
    "_markers.csv"
  ))
}

Idents(combined) <- "comparees"
p <- DotPlot(combined,
  features = feats, cols = c("blue", "red"),
  group.by = "comparees", split.by = "orig.sample"
) + RotatedAxis() + scale_y_discrete(limits = rev)
p
ggsave(p,
  filename = str_c(out_dir, "dotplot_sample.svg"),
  device = "svg", width = 10, height = 11
)

p <- DotPlot(combined,
  features = feats, cols = c(
    "magenta", "blue", "green", "brown"
  ), group.by = "comparees", split.by = "donor"
) + RotatedAxis() + scale_y_discrete(limits = rev)
p
ggsave(p,
  filename = str_c(out_dir, "dotplot_donors.svg"),
  device = "svg", width = 10, height = 11
)

################################################################################
#                          Finding conserved markers                           #
################################################################################

find_d_i_m_in_sample <- function(cluster,
                                 sample,
                                 seurat_object) {
  markers <- FindConservedMarkers(
    subset(seurat_object, subset = orig.sample == sample),
    ident.1 = cluster, only.pos = T, grouping.var = "donor", verbose = F
  )
  if (sample == "BM") sample <- "ybm"
  markers <- markers[order(
    markers[, str_c(sample, "1_avg_log2FC")],
    markers[, str_c(sample, "2_avg_log2FC")],
    decreasing = T
  ), ]
  return(markers)
}

# find donor-inspecific markers (d.i.m.)
find_d_i_m_in_cluster <- function(cluster, seurat_object) {
  print(cluster)
  samples <- c("FL", "BM")
  cluster_result <- lapply(samples, FUN = function(sample) {
    find_d_i_m_in_sample(
      cluster = cluster, sample = sample,
      seurat_object = seurat_object
    )
  })
  names(cluster_result) <- samples
  return(cluster_result)
}

find_donor_inspecific_markers <- function(seurat_object, cluster_list) {
  result <- lapply(cluster_list, FUN = function(cluster) {
    find_d_i_m_in_cluster(
      cluster = cluster,
      seurat_object = seurat_object
    )
  })
  names(result) <- cluster_list
  return(result)
}

conserved_markers_list <- find_donor_inspecific_markers(
  seurat_object = combined, cluster_list = fl_clusters
)

print(lapply(conserved_markers_list, function(x) lapply(x, function(y) dim(y))))

cm_dir <- str_c(out_dir, "conserved_markers/")
dir.create(cm_dir, recursive = T)
write(toJSON(conserved_markers_list), str_c(
  out_dir,
  "conserved_markers_list.json"
))

for (cluster in fl_clusters) {
  for (sample in samples) {
    mlist <- getElement(
      getElement(conserved_markers_list, name = cluster),
      name = sample
    )
    print(str_c(cluster, ": ", sample))
    print(dim(mlist))
    write.table(
      mlist,
      file = str_c(
        cm_dir, cluster, "_", sample, "_specific_markers.csv"
      )
    )
  }
}

conserved_markers_list <- fromJSON(str_c(
  out_dir, "conserved_markers_list.json"
))

################################################################################
#          Using conserved markers to find DE markers between samples          #
################################################################################

find_markers_in_cluster <- function(cluster,
                                    seurat_object) {
  print(cluster)
  result <- lapply(
    list(
      "combo_1" = c("BM", "FL"),
      "combo_2" = c("FL", "BM")
    ),
    function(sample_pair) {
      markers <- FindMarkers(
        object = seurat_object,
        ident.1 = str_c(sample_pair[1], "_", cluster),
        ident.2 = str_c(sample_pair[2], "_", cluster),
        features = rownames(getElement(
          getElement(conserved_markers_list, name = cluster),
          name = sample_pair[1]
        )),
        only.pos = T,
        # logfc.threshold = fc_threshold,
        verbose = F
      )
      markers <- markers[order(
        markers[, "avg_log2FC"],
        decreasing = T
      ), ]
      return(markers)
    }
  )
  names(result) <- c("BM", "FL")
  return(result)
}

find_markers <- function(seurat_object, cluster_list) {
  Idents(seurat_object) <- "comparees"
  result <- lapply(cluster_list, FUN = function(cluster) {
    find_markers_in_cluster(
      cluster = cluster,
      seurat_object = seurat_object
    )
  })
  names(result) <- cluster_list
  return(result)
}

markers_list <- find_markers(
  seurat_object = combined, cluster_list = fl_clusters
)

print(lapply(markers_list, function(x) lapply(x, function(y) dim(y))))

write(toJSON(conserved_markers_list), str_c(
  out_dir,
  "markers_list_from_conserved_lists_only.json"
))

mcl_dir <- str_c(out_dir, "markers_with_genes_from_conserved_lists_only/")
dir.create(mcl_dir, recursive = T)

for (cluster in fl_clusters) {
  for (sample in samples) {
    mlist <- getElement(
      getElement(markers_list, name = cluster),
      name = sample
    )
    print(str_c(cluster, ": ", sample))
    print(dim(mlist))
    write.table(
      mlist,
      file = str_c(
        mcl_dir, cluster, "_", sample, "_specific_markers.csv"
      )
    )
  }
}

feats <- unlist(lapply(markers_list, function(x) {
  lapply(x, function(y) {
    head(rownames(y),
      n = 2
    )
  })
}))
names(feats) <- NULL
feats <- unique(feats)
feats

Idents(combined) <- "comparees"

p <- DotPlot(combined,
  features = feats, cols = c("blue", "red"),
  group.by = "comparees", split.by = "orig.sample"
) + RotatedAxis() + scale_y_discrete(limits = rev)
p
ggsave(p,
  filename = str_c(mcl_dir, "dotplot_sample.svg"),
  device = "svg", width = 10, height = 11
)

p <- DotPlot(combined,
  features = feats, cols = c(
    "magenta", "blue", "green", "brown"
  ), group.by = "comparees", split.by = "donor"
) + RotatedAxis() + scale_y_discrete(limits = rev)
p
ggsave(p,
  filename = str_c(mcl_dir, "dotplot_donors.svg"),
  device = "svg", width = 10, height = 11
)

