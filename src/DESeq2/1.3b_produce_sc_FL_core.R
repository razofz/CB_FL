invisible(lapply(list(
  "stringr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

# comparison DEseq2 res Pseudo bulking approaches

# in_dir <- "data/processed/DEseq2/results/"
# out_dir <- "data/processed/DEseq2/results/"
# sup_table_dir <- "data/raw/paper_specific/Supplemental_tables/"

set.seed(snakemake@config[["seed"]])
clusters <- snakemake@config[["fl_clusters_to_use"]]
deseq_results_files <- snakemake@input[["deseq_results"]]
sub_setter_pos_files <- snakemake@output[["sub_setter_pos"]]
sub_setter_neg_files <- snakemake@output[["sub_setter_neg"]]

fc_threshold <- snakemake@config$deseq_cutoffs$log2foldchange
padj_threshold <- snakemake@config$deseq_cutoffs$adjustedpvalue

named_deseq_results_files <- unlist(
  lapply(clusters,
    FUN = function(cluster) {
      deseq_results_files[str_detect(deseq_results_files, str_c(cluster, "_"))]
    }
  )
)
stopifnot(length(named_deseq_results_files) == length(clusters))
names(named_deseq_results_files) <- clusters

named_sub_setter_pos_files <- unlist(
  lapply(clusters,
    FUN = function(cluster) {
      sub_setter_pos_files[
        str_detect(sub_setter_pos_files, str_c(cluster, "_"))
      ]
    }
  )
)
stopifnot(length(named_sub_setter_pos_files) == length(clusters))
names(named_sub_setter_pos_files) <- clusters
named_sub_setter_neg_files <- unlist(
  lapply(clusters,
    FUN = function(cluster) {
      sub_setter_neg_files[
        str_detect(sub_setter_neg_files, str_c(cluster, "_"))
      ]
    }
  )
)
stopifnot(length(named_sub_setter_neg_files) == length(clusters))
names(named_sub_setter_neg_files) <- clusters

deseq_results <- lapply(clusters, FUN = function(cluster) {
  read.csv(
    file = named_deseq_results_files[cluster],
    header = TRUE, sep = ","
  )
})
names(deseq_results) <- clusters

sub_setter_pos <- function(df) {
  df_sub <- df[!is.na(df$padj), ]
  df_sub <- df_sub[
    df_sub$log2FoldChange > fc_threshold & df_sub$padj < padj_threshold,
  ]
}

sub_setter_neg <- function(df) {
  df_sub <- df[!is.na(df$padj), ]
  df_sub <- df_sub[
    df_sub$log2FoldChange < -fc_threshold & df_sub$padj < padj_threshold,
  ]
}

write_sub_setter <- function(cluster) {
  write.table(sub_setter_pos(deseq_results[[cluster]]),
    # file = paste0(sub_setter_dir, cluster, "_pos.csv"),
    file = named_sub_setter_pos_files[[cluster]],
    sep = ","
  )
  write.table(sub_setter_neg(deseq_results[[cluster]]),
    # file = paste0(sub_setter_dir, cluster, "_neg.csv"),
    file = named_sub_setter_neg_files[[cluster]],
    sep = ","
  )
}

invisible(sapply(clusters, write_sub_setter))

prim_pos <- intersect(
  intersect(
    intersect(
      rownames(sub_setter_pos(deseq_results[["HSC"]])),
      rownames(sub_setter_pos(deseq_results[["MPP.I"]]))
    ),
    rownames(sub_setter_pos(deseq_results[["MPP.II"]]))
  ), rownames(sub_setter_pos(deseq_results[["Cyc"]]))
)

lymph_pos <- intersect(
  rownames(sub_setter_pos(deseq_results[["Ly.I"]])),
  rownames(sub_setter_pos(deseq_results[["Ly.II"]]))
)

mye_pos <- intersect(
  intersect(
    rownames(sub_setter_pos(deseq_results[["DC.Mono"]])),
    rownames(sub_setter_pos(deseq_results[["GMP"]]))
  ),
  rownames(sub_setter_pos(deseq_results[["MEP"]]))
)

all_core_pos <- unique(c(prim_pos, lymph_pos, mye_pos))
fl_core_pos <- intersect(intersect(prim_pos, lymph_pos), mye_pos)

prim_neg <- intersect(
  intersect(
    intersect(
      rownames(sub_setter_neg(deseq_results[["HSC"]])),
      rownames(sub_setter_neg(deseq_results[["MPP.I"]]))
    ),
    rownames(sub_setter_neg(deseq_results[["MPP.II"]]))
  ),
  rownames(sub_setter_neg(deseq_results[["Cyc"]]))
)

lymph_neg <- intersect(
  rownames(sub_setter_neg(deseq_results[["Ly.I"]])),
  rownames(sub_setter_neg(deseq_results[["Ly.II"]]))
)
mye_neg <- intersect(
  intersect(
    rownames(sub_setter_neg(deseq_results[["DC.Mono"]])),
    rownames(sub_setter_neg(deseq_results[["GMP"]]))
  ),
  rownames(sub_setter_neg(deseq_results[["MEP"]]))
)

all_core_neg <- unique(c(prim_neg, lymph_neg, mye_neg))
fl_core_neg <- intersect(intersect(prim_neg, lymph_neg), mye_neg)

core <- list(
  prim_pos, lymph_pos, mye_pos, fl_core_pos, all_core_pos,
  prim_neg, lymph_neg, mye_neg, fl_core_neg, all_core_neg
)

max_length <- max(sapply(core, length))
core_filled <- sapply(core, function(x) {
  c(x, rep(NA, max_length - length(x)))
})

core_filled <- data.frame(core_filled)
colnames(core_filled) <- c(
  "prim_pos", "lymph_pos", "mye_pos", "FL_core_pos",
  "all_core_pos", "prim_neg", "lymph_neg", "mye_neg",
  "FL_core_neg", "all_core_neg"
)

write.table(
  x = core_filled,
  file = snakemake@output[["de_genes"]],
  # file = paste0(out_dir, "pseudorep_DE_genes.csv"),
  quote = FALSE, sep = ","
)

write.table(all_core_neg,
  file = snakemake@output[["adult_signature"]],
  # file = paste0(
  #   # sup_table_dir,
  #   out_dir,
  #   "Adult_signature_p005.txt"),
  row.names = F, col.names = F, quote = F
)

write.table(all_core_pos,
  file = snakemake@output[["fetal_signature"]],
  # file = paste0(
  #   # sup_table_dir,
  #   out_dir,
  #   "Fetal_signature_p005.txt"
  # ),
  row.names = F, col.names = F, quote = F
)
