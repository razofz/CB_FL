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

make_named_version <- function(ids, file_list, pattern = "_") {
  named_file_list <- unlist(
    lapply(ids,
      FUN = function(id) {
        file_list[str_detect(string = file_list, pattern = str_c(id, pattern))]
      }
    )
  )
  names(named_file_list) <- ids
  return(named_file_list)
}

deg_results_bm_files <- snakemake@output[["deg_results_bm"]]
deg_results_fl_files <- snakemake@output[["deg_results_fl"]]
named_deg_results_bm <- make_named_version(clusters, deg_results_bm_files)
named_deg_results_fl <- make_named_version(clusters, deg_results_fl_files)
names(named_deg_results_bm) <- clusters
names(named_deg_results_fl) <- clusters

# donors <- c(str_c("yBM", 1:2), str_c("FL", 1:2))
donors <- snakemake@config[["donors"]]
# samples <- c("FL", "BM")
samples <- snakemake@config[["samples"]]

cells_files <- snakemake@output[["cells"]]
named_cells_files <- make_named_version(donors, cells_files)
names(named_cells_files) <- donors

matrices_files <- snakemake@input[["cite_cell_matrices"]]
named_matrices_files <- make_named_version(donors, matrices_files)
names(named_matrices_files) <- donors


metadata <- read.table(snakemake@input[["combined_metadata"]], sep = ",")
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

for (donor in donors) {
  if (str_detect(pattern = "BM", string = donor)) {
    formatted_donor <- str_to_lower(donor)
  } else if (str_detect(pattern = "FL", string = donor)) {
    formatted_donor <- donor
  }
  write.table(
    rownames(metadata[str_detect(metadata$orig.ident, formatted_donor), ]),
    file = named_cells_files[[donor]]
  )
}


object_list <- future_lapply(donors,
  FUN = function(donor) {
    path <- named_matrices_files[[donor]]
    sample <- str_sub(donor, end = -2)
    if (sample == "yBM") sample <- "BM"
    print(str_c("> Loading ", donor, " data.."))
    print(str_c("> from ", named_matrices_files[[donor]]))
    mat <- Read10X(path)$`Gene Expression`
    # if (str_detect(pattern = "BM", string = donor)) {
    #   colnames_mat <- str_c(
    #     "BM_", str_to_lower(donor), "_",
    #     colnames(mat)
    #   )
    # } else if (str_detect(pattern = "FL", string = donor)) {
    #   colnames_mat <- str_c(
    #     "FL_", donor, "_",
    #     colnames(mat)
    #   )
    # }
    if (str_detect(pattern = "BM", string = donor)) {
      formatted_donor <- str_to_lower(donor)
    } else if (str_detect(pattern = "FL", string = donor)) {
      formatted_donor <- donor
    }
    colnames(mat) <- str_c(
      sample, "_", formatted_donor, "_",
      colnames(mat)
    )
    # colnames(mat) <- colnames_mat
    # colnames(mat) <- str_replace(
    #   colnames(mat),
    #   pattern = "hpc_",
    #   replacement = ""
    # )
    obj <- CreateSeuratObject(
      counts = mat,
      project = donor,
      assay = "RNA"
    )
    cells_to_keep <- read.table(named_cells_files[[donor]])$x
    obj <- subset(obj, cells = cells_to_keep)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(
      obj,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = F
    )
    obj[["orig.sample"]] <- sample
    obj[["orig.donor"]] <- donor
    return(obj)
  },
  future.seed = TRUE
)
names(object_list) <- donors

stopifnot(all(sapply(donors, FUN = function(donor) {
  donor == unique(object_list[[donor]][[]][["orig.donor"]])
})))

features <- SelectIntegrationFeatures(object.list = object_list)
anchors <- FindIntegrationAnchors(
  object.list = object_list, anchor.features = features
)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

# combined[["donor"]] <- sapply(rownames(combined[[]])[1:10], FUN = function(x) {
#   unlist(str_split(x, pattern = "_"))[2]
# })
combined[["donor"]] <- combined[["orig.donor"]]

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

fl_metadata <- read.table(
  snakemake@input[["fl_metadata"]],
  sep = ","
)

# foo <- read.table(
#   "data/processed/seurat/CS22/FL_combined_metadata.csv",
#   sep = ","
# )


mapping_predictions <- unlist(
  fromJSON(txt = snakemake@input[["mapping_predictions"]])$yBM_HPC
)
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
fl_metadata[fl_metadata$clust.names == "ERY", "clust.names"] <- "MEP"
md[rownames(fl_metadata), "fl_cluster"] <- fl_metadata$clust.names
table(md$fl_cluster)
md$fl_cluster <- str_replace(
  pattern = "\\/", replacement = "\\-", string =
    md$fl_cluster
)

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
fl_clusters <- fl_clusters[fl_clusters != "T"]
fl_clusters <- fl_clusters[fl_clusters != "ERY"]

print(table(combined[["comparees"]]))


DefaultAssay(combined) <- "RNA"
Idents(combined) <- "comparees"
for (cluster in fl_clusters) {
  # for (cluster in c("HSC")) {
  print(cluster)
  fc_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["log2foldchange"]]
  padj_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["adjustedpvalue"]]
  formatted_cluster <- str_replace(
    pattern = "-",
    replacement = ".",
    string = cluster
  )
  markers <- FindMarkers(combined,
    ident.1 = str_c("BM_", cluster),
    ident.2 = str_c("FL_", cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  markers <- markers[markers$p_val_adj < padj_threshold, ]
  print(named_deg_results_bm[[formatted_cluster]])
  write.table(markers, file = named_deg_results_bm[[formatted_cluster]])
  markers <- FindMarkers(combined,
    ident.1 = str_c("FL_", cluster),
    ident.2 = str_c("BM_", cluster),
    only.pos = T,
    logfc.threshold = fc_threshold,
    verbose = FALSE
  )
  markers <- markers[order(markers$avg_log2FC, decreasing = T), ]
  markers <- markers[markers$p_val_adj < padj_threshold, ]
  print(named_deg_results_fl[[formatted_cluster]])
  write.table(markers, file = named_deg_results_fl[[formatted_cluster]])
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
