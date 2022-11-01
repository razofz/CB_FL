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
gsea_fc_rnk_files <- snakemake@output[["gsea_fc_rnk"]]

print(clusters)

named_cts_files <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      cts_files[str_detect(cts_files, str_c(cl, "_"))]
    }
  )
)
stopifnot(length(named_cts_files) == length(clusters))
named_coldata_files <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      coldata_files[str_detect(coldata_files, str_c(cl, "_"))]
    }
  )
)
names(named_cts_files) <- clusters
names(named_coldata_files) <- clusters

named_normed_counts_files <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      normed_counts_files[str_detect(normed_counts_files, str_c(cl, "_"))]
    }
  )
)
named_deseq_results_files <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      deseq_results_files[str_detect(deseq_results_files, str_c(cl, "_"))]
    }
  )
)
named_gsea_fc_rnk_files <- unlist(
  lapply(clusters,
    FUN = function(cl) {
      gsea_fc_rnk_files[str_detect(gsea_fc_rnk_files, str_c(cl, "_"))]
    }
  )
)
names(named_deseq_results_files) <- clusters
names(named_normed_counts_files) <- clusters
names(named_gsea_fc_rnk_files) <- clusters

if (snakemake@wildcards[["fl_core_version"]] == "FLcorePseudotech") {
  deseq_design <- ~ pseudo_rep + sample
} else if (snakemake@wildcards[["fl_core_version"]] == "FLcoreNoPseudotech") {
  deseq_design <- ~sample
}

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
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  res <- results(dds, contrast = c("sample", "FL", "yBM"))
  write.table(res,
    named_deseq_results_files[cluster],
    sep = ","
  )

  # create dataframe containing log2 fold change and no NAs for GSEA
  res_table <- data.frame(res)
  res_table <- subset.data.frame(res_table, select = "log2FoldChange")
  res_table <- na.omit(res_table)
  write.table(res_table,
    file = named_gsea_fc_rnk_files[[cluster]],
    sep = "\t"
  )

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

dds <- DESeqDataSetFromMatrix(
  countData = cts_combined, colData = coldata_combined,
  design = deseq_design
)

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = T)
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
