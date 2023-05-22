invisible(lapply(list(
  "stringr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

clusters <- snakemake@config[["fl_clusters_to_use"]]
cts_files <- snakemake@input[["cts_files"]]
coldata_files <- snakemake@input[["coldata_files"]]
normed_counts_files <- snakemake@output[["normed_counts"]]
deseq_results_files <- snakemake@output[["deseq_results"]]
# gsea_fc_rnk_files <- snakemake@output[["gsea_fc_rnk"]]

stopifnot(all(
  as.logical(
    lapply(
      list(
        clusters,
        cts_files,
        coldata_files,
        normed_counts_files,
        deseq_results_files
        # gsea_fc_rnk_files
      ),
      FUN = function(x) !is.null(x)
    )
  )
))

is_pseudotech <- F
if (!is.null(snakemake@wildcards[["fl_core_version"]])) {
  if (snakemake@wildcards[["fl_core_version"]] == "FLcorePseudotech") {
    deseq_design <- ~ sample
    is_pseudotech <- T
  }
} else {
  deseq_design <- ~sample
}

print(clusters)

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

named_cts_files <- make_named_version(clusters, cts_files)
named_coldata_files <- make_named_version(clusters, coldata_files)
named_normed_counts_files <- make_named_version(clusters, normed_counts_files)
named_deseq_results_files <- make_named_version(clusters, deseq_results_files)
# named_gsea_fc_rnk_files <- make_named_version(clusters, gsea_fc_rnk_files)

stopifnot(all(
  as.logical(
    lapply(
      list(
        named_cts_files,
        named_coldata_files,
        named_normed_counts_files,
        named_deseq_results_files
        # named_gsea_fc_rnk_files
      ),
      FUN = function(x) length(x) == length(clusters)
    )
  )
))

for (cluster in clusters) {
  cts <- read.csv(named_cts_files[cluster],
    row.names = 1
  )
  coldata <- read.csv(named_coldata_files[cluster], row.names = 1)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = deseq_design
  )
  dds <- DESeq(dds, test = "Wald")
  res <- results(dds, contrast = c("sample", "FL", "yBM"))
  write.table(res,
    named_deseq_results_files[cluster],
    sep = ","
  )
  # create dataframe containing log2 fold change and no NAs for GSEA
  res_table <- data.frame(res)
  res_table <- subset.data.frame(res_table, select = "log2FoldChange")
  res_table <- na.omit(res_table)
  # write.table(res_table,
  #   file = named_gsea_fc_rnk_files[[cluster]],
  #   sep = "\t"
  # )
  norm_counts <- counts(dds, normalized = T)
  write.table(norm_counts,
    file = named_normed_counts_files[cluster],
    append = F, sep = "\t",
    row.names = TRUE, col.names = TRUE, quote = FALSE
  )
}

cts_combined <- data.frame()
coldata_combined <- data.frame()
for (cluster_idx in seq_len(length(clusters))) {
  cts <- read.csv(named_cts_files[clusters[cluster_idx]],
    row.names = 1
  )
  coldata <- read.csv(named_coldata_files[clusters[cluster_idx]],
    row.names = 1
  )
  if (cluster_idx == 1) {
    cts_combined <- cts
    coldata_combined <- coldata
  } else {
    cts_combined <- cbind(cts_combined, cts)
    coldata_combined <- rbind(coldata_combined, coldata)
  }
}

cts_combined <- cts_combined[!grepl(
  pattern = "Ly.III", x = colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "DC.I", x = colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "_T_", x = colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "DC.I", x = rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "Ly.III", x = rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "_T_", x = rownames(coldata_combined)
), ]
cts_combined <- cts_combined[!grepl(
  pattern = "yBM_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "yBM_", x =
    rownames(coldata_combined)
), ]

dds <- DESeqDataSetFromMatrix(
  countData = cts_combined, colData = coldata_combined,
  design = ~1
  # design = ~pseudo_rep + cluster_name
  # design = deseq_design
)

# dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = T)

# norm_counts <- counts(dds, normalized = F)
# norm_counts <- counts(dds, normalized = T)
write.table(norm_counts,
  file = snakemake@output[["fl_only_normed_counts"]],
  append = F, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
)

# GSEA stuff

# dds@colData$cluster_name <- as.factor(dds@colData$cluster_name)

# res <- results(dds, contrast = c("cluster_name", "MPP.I", "HSC"))

# # create dataframe containing log2 fold change and no NAs for GSEA
# res_table <- data.frame(res)
# res_table <- subset.data.frame(res_table, select = "log2FoldChange")
# res_table <- na.omit(res_table)
# write.table(res_table,
#   # paste0(out_dir, "/results/", "FL_MPPI_HSC.rnk"),
#   snakemake@output[["mpp1_hsc_rnk"]],
#   sep = "\t"
# )

# res <- results(dds, contrast = c("cluster_name", "MPP.I", "MPP.II"))

# res_table <- data.frame(res)
# res_table <- subset.data.frame(res_table, select = "log2FoldChange")
# res_table <- na.omit(res_table)

# write.table(res_table,
#   # paste0(out_dir, "/results/", "FL_MPPI_MPPII.rnk"),
#   snakemake@output[["mpp1_mpp2_rnk"]],
#   sep = "\t"
# )
