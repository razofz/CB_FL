source("env.R")
raw_dir <- "data/raw/paper_specific/DEseq2/gate_wise"
# raw_dir <- "data/processed/notebooks/DEseq2/gate_wise"
out_dir <- "data/processed/DEseq2"
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))

files <- list.files(
  path = raw_dir, pattern = "*_coldata_all_hpc.csv",
  full.names = TRUE, recursive = FALSE
)

for (x in files) {
  names <- str_split(str_split(x, "/", simplify = T)[6], "_coldata_all_hpc.csv",
    simplify = T
  )[1]
  print(paste0("> names: ", names))

  cts <- read.csv(paste0(raw_dir, "/", names, "_cts_all_hpc.csv"),
    row.names = 1
  )
  coldata <- read.csv(x, row.names = 1)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~sample
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  res <- results(dds, contrast = c("sample", "FL", "yBM"))
  write.table(res, paste0(out_dir, "/results/gated_", names, "_FL_yBM.csv"),
    sep = ","
  )
  norm_counts <- counts(dds, normalized = T)
  write.table(norm_counts,
    file = paste0(
      out_dir,
      "/results/Normed_counts/gated_", names,
      "_DEseq_normed.txt"
    ), append = TRUE,
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
  )
}
