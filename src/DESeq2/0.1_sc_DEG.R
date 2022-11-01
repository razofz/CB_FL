invisible(lapply(list(
  "stringr",
  "future.apply",
  "ggplot2",
  "jsonlite",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])
threads <- snakemake@threads
# print(threads)
plan(multisession, workers = threads)
options(future.globals.maxSize = 16e3 * 1024^2)

clusters <- snakemake@config[["fl_clusters_to_use"]]

# out_dir <- str_c(
#   Sys.getenv("PROJECT_PATH"),
#   "/data/processed/DEseq2/"
# )
# plots_dir <- str_c(out_dir, "/images/")

deg_results_bm_files <- snakemake@output[["deg_results_bm"]]
deg_results_fl_files <- snakemake@output[["deg_results_fl"]]
named_deg_results_bm <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      deg_results_bm_files[str_detect(deg_results_bm_files, str_c(cl, "_"))]
    }
  )
)
named_deg_results_fl <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      deg_results_fl_files[str_detect(deg_results_fl_files, str_c(cl, "_"))]
    }
  )
)
names(named_deg_results_bm) <- clusters
names(named_deg_results_fl) <- clusters

donors <- c(str_c("yBM", 1:2), str_c("FL", 1:2))
samples <- c("FL", "BM")
config_samples <- snakemake@config[["samples"]]

metadata <- read.table(
  snakemake@input[["combined_metadata"]],
  # str_c(
  #   Sys.getenv("PROJECT_PATH"),
  #   "/data/raw/paper_specific/from_seurat/FL_BM_combined_metadata.csv"
  # ),
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

for (sample in config_samples) {
  if (str_detect(pattern = "BM", string = sample)) {
    formatted_sample <- str_to_lower(sample)
  } else if (str_detect(pattern = "FL", string = sample)) {
    formatted_sample <- sample
  }
  path <- str_c(
    dirname(
      snakemake@output[["cells"]][[1]]
    ),
    "/", sample, "_cells.csv"
  )
  # print(path)
  write.table(
    rownames(metadata[str_detect(metadata$orig.ident, formatted_sample), ]),
    file = path
    # file = str_c(out_dir, sample, "_cells.csv")
  )
}


# print(snakemake@input[["cite_cell_matrices_dir"]])
# print(dir(snakemake@input[["cite_cell_matrices_dir"]]))

object_list <- future_lapply(donors,
  FUN = function(donor) {
    path <- str_c(
      snakemake@input[["cite_cell_matrices_dir"]],
      "/",
      donor,
      "_hpc/"
    )
    print(str_c("> Loading ", donor, " data.."))
    mat <- Read10X(path)$`Gene Expression`
    if (str_detect(pattern = "BM", string = donor)) {
      colnames_mat <- str_c(
        "BM_", str_to_lower(donor), "_",
        colnames(mat)
      )
    } else if (str_detect(pattern = "FL", string = donor)) {
      colnames_mat <- str_c(
        "FL_", donor, "_",
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
      project = donor,
      assay = "RNA"
    )
    sample <- str_sub(donor, end = -2)
    # if (str_detect(pattern = "BM", string = sample)) {
    #   formatted_sample <- str_sub(sample, start = 2)
    # } else if (str_detect(pattern = "FL", string = sample)) {
    #   formatted_sample <- sample
    # }
    cells_to_keep <- read.table(
      str_c(
        dirname(
          snakemake@output[["cells"]][[1]]
        ), "/", sample, "_cells.csv"
      )
    )$x
    # cells_to_keep <- read.table(str_c(out_dir, donor, "_cells.csv"))$x
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

# Idents(combined) <- "donor"
# p <- DimPlot(combined, reduction = "umap") + coord_fixed()
# p
# ggsave(p,
#   filename = str_c(plots_dir, "umap.svg"),
#   device = "svg", width = 7, height = 7
# )

# p <- DimPlot(combined, reduction = "umap", split.by = "orig.sample") +
#   coord_fixed()
# p
# ggsave(p,
#   filename = str_c(plots_dir, "umap_split.svg"),
#   device = "svg", width = 11, height = 5
# )

# bm both boys
# fl one boy one girl

# path <- str_c(
#   Sys.getenv("PROJECT_PATH"),
#   "/data/processed/notebooks/mapping_predictions/RNA_FL_target_preds.json"
# )
path <- snakemake@input[["mapping_predictions"]]

fl_metadata <- read.table(
  snakemake@input[["fl_metadata"]],
  # str_c(
  #   Sys.getenv("PROJECT_PATH"),
  #   "/data/processed/seurat/CS22/FL_combined_metadata.csv"
  # ),
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

# DimPlot(combined,
#   reduction = "umap", group.by = "fl_cluster",
#   split.by = "orig.sample"
# ) + coord_fixed()

# DimPlot(combined,
#   reduction = "umap", group.by = "fl_cluster",
#   split.by = "donor"
# ) + coord_fixed()

DefaultAssay(combined) <- "RNA"
combined[["comparees"]] <- str_c(
  combined[[]]$"orig.sample", "_",
  combined[[]]$"fl_cluster"
)
Idents(combined) <- "comparees"

fl_clusters <- unique(combined[[]][, "fl_cluster"])
fl_clusters <- fl_clusters[fl_clusters != "NA"]
fl_clusters <- fl_clusters[!is.na(fl_clusters)]
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
  "BM_Cyc", "BM_NA", "FL_NA"
)
tmp <- as.factor(combined[[]][, "comparees"])
levels(tmp) <- ordering
combined[["comparees"]] <- tmp

# simple_dir <- str_c(out_dir, "markers_simple_way/")
# dir.create(simple_dir, recursive = T)

DefaultAssay(combined) <- "RNA"
Idents(combined) <- "comparees"
# for (cluster in fl_clusters) {
for (cluster in clusters) {
  # for (cluster in c("HSC")) {
  print(cluster)
  fc_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["log2foldchange"]]
  formatted_cluster <- str_replace(
    pattern = "\\.",
    replacement = "-",
    string = cluster
  )
  markers <- FindMarkers(combined,
    ident.1 = str_c("BM_", formatted_cluster),
    ident.2 = str_c("FL_", formatted_cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  # write.table(markers, file = named_deg_results_bm[[str_replace(
  #   pattern = "-",
  #   replacement = "\\.",
  #   string = cluster
  # )]])
  print(named_deg_results_bm[[cluster]])
  write.table(markers, file = named_deg_results_bm[[cluster]])
  markers <- FindMarkers(combined,
    ident.1 = str_c("FL_", formatted_cluster),
    ident.2 = str_c("BM_", formatted_cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  print(named_deg_results_fl[[cluster]])
  write.table(markers, file = named_deg_results_fl[[cluster]])
}

# Idents(combined) <- "comparees"
# p <- DotPlot(combined,
#   features = feats, cols = c("blue", "red"),
#   group.by = "comparees", split.by = "orig.sample"
# ) + RotatedAxis() + scale_y_discrete(limits = rev)
# p
# ggsave(p,
#   filename = str_c(out_dir, "dotplot_sample.svg"),
#   device = "svg", width = 10, height = 11
# )

# p <- DotPlot(combined,
#   features = feats, cols = c(
#     "magenta", "blue", "green", "brown"
#   ), group.by = "comparees", split.by = "donor"
# ) + RotatedAxis() + scale_y_discrete(limits = rev)
# p
# ggsave(p,
#   filename = str_c(out_dir, "dotplot_donors.svg"),
#   device = "svg", width = 10, height = 11
# )
