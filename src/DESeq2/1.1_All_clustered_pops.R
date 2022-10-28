source("env.R")
raw_dir <- "data/raw/paper_specific/DEseq2/cluster_wise"
# raw_dir <- "data/processed/notebooks/DEseq2/cluster_wise"
out_dir <- "data/processed/DEseq2"
if (dir.exists(out_dir) == FALSE) {
  dir.create(out_dir, showWarnings = TRUE, recursive = TRUE)
}
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))

files <- list.files(
  path = raw_dir, pattern = "*_coldata_all_hpc.csv",
  full.names = TRUE, recursive = FALSE
)
files <- files[-grep("DC.I", files)]
files <- files[-grep("Ly.III", files)]
print(files)
if (dir.exists(paste0(out_dir, "/results/Normed_counts/")) == FALSE) {
  dir.create(paste0(out_dir, "/results/Normed_counts/"),
    showWarnings = TRUE,
    recursive = TRUE
  )
}

for (x in files) {
  names <- str_split(str_split(x, "/", simplify = T)[6], "_coldata_all_hpc.csv",
    simplify = T
  )[1]
  print(paste0("> names: ", names))

  cts <- read.csv(paste0(raw_dir, "/", names, "_cts_all_hpc.csv"),
    row.names =
      1
  )
  coldata <- read.csv(x, row.names = 1)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~pseudo_rep + sample
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  res <- results(dds, contrast = c("sample", "FL", "yBM"))
  write.table(res, paste0(out_dir, "/results/", names, "_FL_yBM.csv"),
    sep = ","
  )
  norm_counts <- counts(dds, normalized = T)
  write.table(norm_counts,
    file = paste0(
      out_dir, "/results/Normed_counts/",
      names, "_DEseq_normed.txt"
    ), append = F, sep = "\t",
    row.names = TRUE, col.names = TRUE, quote = FALSE
  )
}


cts_combined <- data.frame()
coldata_combined <- data.frame()
names <- NULL

for (x in 1:length(files)) {
  if (x == 1) {
  names <- str_split(str_split(files[x], "/", simplify = T)[6],
    "_coldata_all_hpc.csv",
    simplify = T
  )[1]
    print(names)
  cts <- read.csv(paste0(raw_dir, "/", names, "_cts_all_hpc.csv"),
    row.names = 1
  )
    cts_combined <- cts
    coldata <- read.csv(files[x], row.names = 1)
    coldata_combined <- coldata
  } else {
  names <- str_split(str_split(files[x], "/", simplify = T)[6],
    "_coldata_all_hpc.csv",
    simplify = T
  )[1]
    print(names)
  cts <- read.csv(paste0(raw_dir, "/", names, "_cts_all_hpc.csv"),
    row.names
    = 1
  )
    cts_combined <- cbind(cts_combined, cts)
    coldata <- read.csv(files[x], row.names = 1)
    coldata_combined <- rbind(coldata_combined, coldata)
  }
}

cts_combined <- cts_combined[!grepl(
  pattern = "Ly.III", x =
    colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "DC.I", x =
    colnames(cts_combined)
)]
cts_combined <- cts_combined[!grepl(
  pattern = "_T_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "DC.I", x =
    rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "Ly.III", x =
    rownames(coldata_combined)
), ]
coldata_combined <- coldata_combined[!grepl(
  pattern = "_T_", x =
    rownames(coldata_combined)
), ]
cts_combined <- cts_combined[!grepl(
  pattern = "yBM_", x =
    colnames(cts_combined)
)]
coldata_combined <- coldata_combined[!grepl(
  pattern = "yBM_", x =
    rownames(coldata_combined)
), ]
cts_combined_clust <- cts_combined
coldata_combined_clust <- coldata_combined

dds <- DESeqDataSetFromMatrix(
  countData = cts_combined, colData = coldata_combined,
  design = ~pseudo_rep + sample
)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = T)
write.table(norm_counts,
 file = paste0(out_dir, "/results/Normed_counts/", "only_FL_DEseq_normed.txt"),
  append = F, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
)
