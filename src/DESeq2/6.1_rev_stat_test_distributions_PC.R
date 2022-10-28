source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("pheatmap"))

raw_dir <- "data/raw/paper_specific/Supplemental_tables"
external_dir <- "data/external"
out_dir <- "data/processed/DEseq2"
images_dir <- paste0(out_dir, "/images")
if (dir.exists(images_dir) == FALSE) {
  dir.create(images_dir, showWarnings = TRUE, recursive = TRUE)
}
set.seed(12345)

sample_table3 <- read.table(str_c(
  external_dir,
  "/BALL-1988S-HTSeq/B-ALL-subtyping.txt"
  ),
  header = TRUE, dec = ",", sep = "\t"
)
colums_oi <- c(
  "patient", "fusion", "age", "gender", "primary.subtype",
  "RNA.seq.library"
)
sample_table3 <- sample_table3[colums_oi]

fusions_oi3 <- c(
  "KMT2A-AFF1", "KMT2A-MLLT1", "KMT2A-MLLT3", "ETV6-RUNX1",
  "NoFusion", "BCR-ABL1"
)
sample_table3 <- sample_table3[sample_table3$fusion %in% fusions_oi3, ]

subtype_oi <- c(
  "KMT2A", "ETV6-RUNX1", "Ph", "High hyperdiploid",
  "Low hyperdiploid"
)
sample_table3 <- sample_table3[sample_table3$primary.subtype %in% subtype_oi, ]

sample_table3 <- subset.data.frame(sample_table3, !age == "NA")
sample_table3 <- subset.data.frame(sample_table3, !gender == "NA")

for (i in 1:dim(sample_table3)[1]) {
  if (sample_table3$fusion[i] == "NoFusion") {
    sample_table3$fusion[i] <- sample_table3$primary.subtype[i]
  }
}

sample_table3$age_group <- cut(sample_table3$age, c(
  0, 2, 16, 40,
  max(sample_table3$age)
))
levels(sample_table3$age_group) <- c("0-2", "2-16", "16-40", ">40")

# get file names
file_names3 <- c()
for (i in 1:dim(sample_table3)[1]) {
  file_names3[i] <- list.files(
    path = str_c(external_dir, "/BALL-1988S-HTSeq/"),
    pattern = as.character(sample_table3$patient[i]),
    ignore.case = TRUE, full.names = TRUE
  )[1] 
  # hard-coded to only take the first sample,
  # works for this one but might not for others
}
sample_table3$file_name <- file_names3

#Load samples into DEseq2

DEsampleTable3 <- data.frame(
  sampleNames = sample_table3$patient,
  fileName = as.character(sample_table3$file_name),
  age_group = sample_table3$age_group,
  type = sample_table3$fusion,
  age = sample_table3$age, gender = sample_table3$gender,
  RNA_lib = sample_table3$RNA.seq.library
)

counts3 <- DESeqDataSetFromHTSeqCount(
  sampleTable = DEsampleTable3, design =
    ~RNA_lib
)
counts3 <- DESeq(counts3)

vsd3 <- vst(counts3, blind = FALSE)
assay(vsd3) <- limma::removeBatchEffect(assay(vsd3), vsd3$RNA_lib)

# load genes of interest (GOI)
GOI_FL <- read.table(
  file = str_c(raw_dir, "/Fetal_signature_p005.txt"),
  sep = "\t"
)
GOI_AD <- read.table(
  file = str_c(raw_dir, "/Adult_signature_p005.txt"),
  sep = "\t"
)
GOI_AD <- GOI_AD$V1
GOI_FL <- GOI_FL$V1

GOI <- data.frame(
  gene_name = append(GOI_FL, GOI_AD),
  FLorBM = append(rep("FL", length(GOI_FL)), rep("BM", length(GOI_AD)))
)

ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
  path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl"
)
ids <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = GOI$gene_name,
  mart = ensembl
)
select_FL_Core <- rownames(vsd3) %in% ids$ensembl_gene_id
# n_genes <- length(ids$ensembl_gene_id) #when wrong length
n_genes <- as.numeric(table(select_FL_Core)["TRUE"])
#colors and categories

iter <- 10000

#loop over types PC1
# sup6d (keep names)
for (Luek in unique(vsd3$type)) {
  vsd3.sub <- vsd3[, vsd3$type %in% Luek]

  pca <- prcomp(t(assay(vsd3.sub)[select_FL_Core, ])) # calculate PCA
  percentVar_FL_core <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
  comp1 <- NULL
  for (rounds in 1:iter) {
    pca <- prcomp(t(assay(vsd3.sub)[sample(nrow(assay(vsd3.sub)), n_genes), ]))
    # calculate PCA
    percentVar <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution

    comp1 <- c(comp1, percentVar[1])
  }
  pval <- length(which(comp1 >= percentVar_FL_core[1])) / iter
  high_x <- max(c(comp1, percentVar_FL_core[1])) + 0.02
  if (Luek %in% c("High hyperdiploid", "ETV6-RUNX1")) {
    high_x <- .25
  }
  low_x <- min(c(comp1, percentVar_FL_core[1])) - 0.02
  pdf(paste0(images_dir, "/PC1_variance_test_", Luek, "_its_", iter, ".pdf"),
    width = 10, height = 5
  )
  hist(comp1, 100,
    main = paste0("Histogram of variance in ", Luek),
    xlab = "PC1 variance",
    xlim = c(low_x, high_x)
  )
  abline(v = percentVar_FL_core[1], col = "red", lwd = 3, lty = 2)
  text(x = percentVar_FL_core[1] + 0.008, y = 300, labels = paste0(
    "p-val: ",
    pval
  ))
  dev.off()
}

# compare with all leuk subtype in same PCA

pca <- prcomp(t(assay(vsd3)[select_FL_Core, ])) # calculate PCA
percentVar_FL_core <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
comp1 <- NULL
for (rounds in 1:iter) {
  print(rounds)
  pca <- prcomp(t(assay(vsd3)[sample(nrow(assay(vsd3)), n_genes), ]))
  # calculate PCA
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution

  comp1 <- c(comp1, percentVar[1])
}
pval <- length(which(comp1 >= percentVar_FL_core[1])) / iter
high_x <- max(c(comp1, percentVar_FL_core[1])) + 0.02
low_x <- min(c(comp1, percentVar_FL_core[1])) - 0.02
# fig5c
pdf(file = paste0(
  images_dir, "/PC1_variance_test_all_leuk_in_PCA_its_",
  iter, ".pdf"
), width = 10, height = 5)
hist(comp1, 100,
  main = paste0("Histogram of variance in all"),
  xlab = "PC1 variance",
  xlim = c(low_x, high_x)
)
abline(v = percentVar_FL_core[1], col = "red", lwd = 3, lty = 2)
text(x = percentVar_FL_core[1] + 0.008, y = 250, labels = paste0(
  "p-val: ",
  pval
))
dev.off()

