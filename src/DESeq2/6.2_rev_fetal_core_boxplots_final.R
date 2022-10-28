source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("fdrtool"))

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
subtype_oi <- c(
  "KMT2A", "ETV6-RUNX1", "Ph", "High hyperdiploid",
  "Low hyperdiploid"
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
fusions_oi3 <- c(
  "KMT2A-AFF1", "KMT2A-MLLT1", "KMT2A-MLLT3", "ETV6-RUNX1",
  "NoFusion", "BCR-ABL1"
)
sample_table3 <- sample_table3[sample_table3$fusion %in% fusions_oi3, ]

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
GOI_FL <- GOI_FL$V1
GOI_AD <- GOI_AD$V1

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
ids_FL <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = GOI_FL,
  mart = ensembl
)
ids_BM <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = GOI_AD,
  mart = ensembl
)
select_FL_Core <- rownames(vsd3) %in% ids$ensembl_gene_id
select_FL_Core_FL <- rownames(vsd3) %in% ids_FL$ensembl_gene_id
select_FL_Core_BM <- rownames(vsd3) %in% ids_BM$ensembl_gene_id

age_group_col <- c(
  "0-2" = "#F8766D", "2-16" = "#7CAE00",
  "16-40" = "#00BFC4", ">40" = "#C77CFF"
)
types_oi <- c("KMT2A-AFF1", "BCR-ABL1", "High hyperdiploid")
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib")

sig_levels_FL <- NULL
sig_levels_BM <- NULL
vsd3.sub_FL <- NULL
vsd3.sub_BM <- NULL

# Fetal_core_up_box_plot_KMT2A-AFF1,Fetal_core_up_box_plot_BCR-ABL1,Fetal_core_up_box_plot_High hyperdiploid
# Fetal_core_down_box_plot_KMT2A-AFF1,Fetal_core_down_box_plot_BCR-ABL1,Fetal_core_down_box_plot_High hyperdiploid

# fig5e (all of them)
for (leuk in types_oi) {
  vsd3.sub <- vsd3[, vsd3$type %in% leuk]
  # add grouping if wanted
  intgroup.df <- as.data.frame(colData(vsd3.sub)[, intgroup, drop = FALSE])

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(vsd3.sub)[[intgroup]]
  }

  vsd3.sub_FL <- as.data.frame(t(assay(vsd3.sub)[select_FL_Core_FL, ]))
  vsd3.sub_FL <- sapply(vsd3.sub_FL, function(vsd3.sub_FL) {
    (vsd3.sub_FL -
      mean(vsd3.sub_FL)) /
      sd(vsd3.sub_FL)
  })
  MEANS_FL <- rowMeans(vsd3.sub_FL)
  vsd3.sub_FL <- cbind.data.frame(intgroup.df, MEANS_FL)
  p1 <- ggplot(data = vsd3.sub_FL, aes_string(
    y = "MEANS_FL", x = "age_group",
    fill = "age_group"
  )) +
    geom_boxplot() +
    scale_fill_manual(values = age_group_col) +
    coord_fixed() +
    theme_bw() +
    theme(aspect.ratio = 2) +
    ylab("Z-scored mean expression") +
    xlab("") +
    ggtitle(paste0(leuk, " Fetal core up"))
  print(p1)
  ggsave(paste0(images_dir, "/Fetal_core_up_box_plot_", leuk, ".pdf"),
    plot =
      last_plot()
  )
  vsd3.sub_BM <- as.data.frame(t(assay(vsd3.sub)[select_FL_Core_BM, ]))
  vsd3.sub_BM <- sapply(vsd3.sub_BM, function(vsd3.sub_BM) {
    (vsd3.sub_BM -
      mean(vsd3.sub_BM)) /
      sd(vsd3.sub_BM)
  })
  MEANS_BM <- rowMeans(vsd3.sub_BM)
  vsd3.sub_BM <- cbind.data.frame(intgroup.df, MEANS_BM)
  p2 <- ggplot(data = vsd3.sub_BM, aes_string(
    y = "MEANS_BM", x = "age_group",
    fill = "age_group"
  )) +
    geom_boxplot() +
    scale_fill_manual(values = age_group_col) +
    coord_fixed() +
    theme_bw() +
    theme(aspect.ratio = 2) +
    ylab("Z-scored mean expression") +
    xlab("") +
    ggtitle(paste0(leuk, " Fetal core down"))

  print(p2)
  ggsave(paste0(images_dir, "/Fetal_core_down_box_plot_", leuk, ".pdf"),
    plot =
      last_plot()
  )

  res.aov_FL <- aov(MEANS_FL ~ age_group, data = vsd3.sub_FL)
  sig_levels_FL <- rbind.data.frame(
    sig_levels_FL,
    TukeyHSD(res.aov_FL)[["age_group"]][, 4]
  )

  res.aov_BM <- aov(MEANS_BM ~ age_group, data = vsd3.sub_BM)
  sig_levels_BM <- rbind.data.frame(
    sig_levels_BM,
    TukeyHSD(res.aov_BM)[["age_group"]][, 4]
  )
}

colnames(sig_levels_FL) <- rownames(TukeyHSD(res.aov_FL)[["age_group"]])
rownames(sig_levels_FL) <- types_oi

colnames(sig_levels_BM) <- rownames(TukeyHSD(res.aov_BM)[["age_group"]])
rownames(sig_levels_BM) <- types_oi

write.csv(sig_levels_FL,
  file = paste0(out_dir, "/Pvals_anova_FL_up2.csv"),
  quote = F
)
write.csv(sig_levels_BM,
  file = paste0(out_dir, "/Pvals_anova_FL_down2.csv"),
  quote = F
)
