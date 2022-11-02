invisible(lapply(list(
  "stringr",
  "DESeq2",
  "biomaRt",
  "limma",
  "ggplot2",
  "pheatmap"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

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

leukemia_types <- snakemake@config[["leukemia_types"]]
print(leukemia_types)

named_plots_filenames <- make_named_version(
  ids = leukemia_types,
  file_list = snakemake@output[["plots_random_pc1_variance"]]
)

sample_table <- read.table(snakemake@input[["ball_tsv"]],
  header = TRUE, dec = ",", sep = "\t"
)


colums_oi <- c(
  "patient", "fusion", "age", "gender", "primary.subtype",
  "RNA.seq.library"
)
sample_table <- sample_table[colums_oi]

fusions_oi <- c(
  "KMT2A-AFF1", "KMT2A-MLLT1", "KMT2A-MLLT3", "ETV6-RUNX1",
  "NoFusion", "BCR-ABL1"
)
sample_table <- sample_table[sample_table$fusion %in% fusions_oi, ]

subtype_oi <- c(
  "KMT2A", "ETV6-RUNX1", "Ph", "High hyperdiploid",
  "Low hyperdiploid"
)
sample_table <- sample_table[sample_table$primary.subtype %in% subtype_oi, ]

sample_table <- subset.data.frame(sample_table, !age == "NA")
sample_table <- subset.data.frame(sample_table, !gender == "NA")

for (i in 1:dim(sample_table)[1]) {
  if (sample_table$fusion[i] == "NoFusion") {
    sample_table$fusion[i] <- sample_table$primary.subtype[i]
  }
}

sample_table$age_group <- cut(sample_table$age, c(
  0, 2, 16, 40,
  max(sample_table$age)
))
levels(sample_table$age_group) <- c("0-2", "2-16", "16-40", ">40")

# get file names
file_names <- c()
for (i in 1:dim(sample_table)[1]) {
  file_names[i] <- list.files(
    path = snakemake@input[["ball_dir"]],
    # path = str_c(external_dir, "/BALL-1988S-HTSeq/"),
    pattern = as.character(sample_table$patient[i]),
    ignore.case = TRUE, full.names = TRUE
  )[1]
  # hard-coded to only take the first sample,
  # works for this one but might not for others
}
sample_table$file_name <- file_names

#Load samples into DEseq2

de_sample_table <- data.frame(
  sampleNames = sample_table$patient,
  fileName = as.character(sample_table$file_name),
  age_group = sample_table$age_group,
  type = sample_table$fusion,
  age = sample_table$age, gender = sample_table$gender,
  RNA_lib = sample_table$RNA.seq.library
)

counts <- DESeqDataSetFromHTSeqCount(
  sampleTable = de_sample_table, design = ~RNA_lib
)
counts <- DESeq(counts)

vsd <- vst(counts, blind = FALSE)
# plotPCA(vsd, "RNA_lib") + theme_bw() + coord_fixed(ratio = 8 / 6.3)

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)
# plotPCA(vsd, "RNA_lib") + theme_bw() + coord_fixed(ratio = 5.6 / 7.2)

# load genes of interest (GOI)
goi_fl <- read.table(
  file = snakemake@input[["fetal_signature"]],
 sep = "\t"
)
goi_ad <- read.table(
  file = snakemake@input[["adult_signature"]],
  sep = "\t"
)
goi_fl <- goi_fl$V1
goi_ad <- goi_ad$V1
goi <- data.frame(
  gene_name = append(goi_fl, goi_ad),
  fl_or_bm = append(rep("FL", length(goi_fl)), rep("BM", length(goi_ad)))
)

ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
  path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl"
)

ids <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = goi$gene_name,
  mart = ensembl
)

# n_genes <- sum(rownames(vsd) %in% ids$ensembl_gene_id)

select_fl_core <- rownames(vsd) %in% ids$ensembl_gene_id
n_genes <- as.numeric(table(select_fl_core)["TRUE"])

iter <- snakemake@config[["random_params"]][["iterations"]]

#loop over types PC1
# sup6d (keep names)
for (leuk in unique(vsd$type)) {
  print(leuk)
  vsd_sub <- vsd[, vsd$type %in% leuk]

  pca <- prcomp(t(assay(vsd_sub)[select_fl_core, ])) # calculate PCA
  percent_var_fl_core <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
  comp1 <- NULL
  for (rounds in 1:iter) {
    pca <- prcomp(t(assay(vsd_sub)[sample(nrow(assay(vsd_sub)), n_genes), ]))
    # calculate PCA
    percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution

    comp1 <- c(comp1, percent_var[1])
  }
  pval <- length(which(comp1 >= percent_var_fl_core[1])) / iter
  high_x <- max(c(comp1, percent_var_fl_core[1])) + 0.02
  if (leuk %in% c("High hyperdiploid", "ETV6-RUNX1")) {
    high_x <- .25
  }
  low_x <- min(c(comp1, percent_var_fl_core[1])) - 0.02
  pdf(
    named_plots_filenames[[leuk]],
    # paste0(images_dir, "/PC1_variance_test_", leuk, "_its_", iter, ".pdf"),
    width = 10, height = 5
  )
  hist(comp1, 100,
    main = paste0("Histogram of variance in ", leuk),
    xlab = "PC1 variance",
    xlim = c(low_x, high_x)
  )
  abline(v = percent_var_fl_core[1], col = "red", lwd = 3, lty = 2)
  text(
    x = percent_var_fl_core[1] + 0.008, y = 300,
    labels = paste0("p-val: ", pval)
  )
  dev.off()
}

# compare with all leuk subtype in same PCA

pca <- prcomp(t(assay(vsd)[select_fl_core, ])) # calculate PCA
percent_var_fl_core <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
comp1 <- NULL
for (rounds in 1:iter) {
  print(rounds)
  # genes_use <- norm_counts_use[sample(
  #   nrow(norm_counts_use), n_genes,
  #   replace = FALSE
  # ), ]
  pca <- prcomp(t(assay(vsd)[sample(nrow(assay(vsd)), n_genes), ]))
  # calculate PCA
  percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution

  comp1 <- c(comp1, percent_var[1])
}
pval <- length(which(comp1 >= percent_var_fl_core[1])) / iter
high_x <- max(c(comp1, percent_var_fl_core[1])) + 0.02
low_x <- min(c(comp1, percent_var_fl_core[1])) - 0.02
# fig5c
pdf(
  file = snakemake@output[["plot_same_pc1"]],
  # file = paste0(
  #   images_dir, "/PC1_variance_test_all_leuk_in_PCA_its_",
  #   iter, ".pdf"
  # ),
  width = 10, height = 5
)
hist(comp1, 100,
  main = paste0("Histogram of variance in all"),
  xlab = "PC1 variance",
  xlim = c(low_x, high_x)
)
abline(v = percent_var_fl_core[1], col = "red", lwd = 3, lty = 2)
text(x = percent_var_fl_core[1] + 0.008, y = 250, labels = paste0(
  "p-val: ",
  pval
))
dev.off()
