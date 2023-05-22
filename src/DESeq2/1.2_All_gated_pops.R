invisible(lapply(list(
  "stringr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

gates <- snakemake@config[["gates_to_use"]]
cts_files <- snakemake@input[["cts_files"]]
coldata_files <- snakemake@input[["coldata_files"]]
normed_counts_files <- snakemake@output[["normed_counts"]]
deseq_results_files <- snakemake@output[["deseq_results"]]

named_cts_files <- unlist(
  lapply(gates,
    FUN = function(gate) {
      cts_files[str_detect(cts_files, str_c(gate, "_"))]
    }
  )
)
stopifnot(length(named_cts_files) == length(gates))
names(named_cts_files) <- gates
named_coldata_files <- unlist(
  lapply(gates,
    FUN = function(gate) {
      coldata_files[str_detect(coldata_files, str_c(gate, "_"))]
    }
  )
)
names(named_coldata_files) <- gates

named_normed_counts_files <- unlist(
  lapply(gates,
    FUN = function(gate) {
      normed_counts_files[str_detect(normed_counts_files, str_c(gate, "_"))]
    }
  )
)
named_deseq_results_files <- unlist(
  lapply(gates,
    FUN = function(gate) {
      deseq_results_files[str_detect(deseq_results_files, str_c(gate, "_"))]
    }
  )
)
names(named_normed_counts_files) <- gates
names(named_deseq_results_files) <- gates

  deseq_design <- ~ sample
# if (snakemake@wildcards[["fl_core_version"]] == "FLcorePseudotech") {
#   deseq_design <- ~ sample
# } else if (snakemake@wildcards[["fl_core_version"]] == "FLcoreNoPseudotech") {
#   deseq_design <- ~sample
# }

for (gate in gates) {
  cts <- read.csv(named_cts_files[gate],
    row.names = 1
  )
  coldata <- read.csv(named_coldata_files[gate], row.names = 1)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = deseq_design
  )
  dds <- DESeq(dds, test = "Wald")
  res <- results(dds, contrast = c("sample", "FL", "yBM"))
  write.table(res,
    file = named_deseq_results_files[gate],
    sep = ","
  )
  norm_counts <- counts(dds, normalized = T)
  write.table(norm_counts,
    file = named_normed_counts_files[gate],
    append = TRUE, sep = "\t", row.names = TRUE,
    col.names = TRUE, quote = FALSE
  )
}
