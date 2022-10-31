#source("env.R")

setwd("~/OneDrive - Lund University/ALL/R scripts/CMLS manuscript (recent things)/RUNX1 PCA 20221028")

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))

set.seed(12345)

sample_table <- read.table("dev_cell_samples.txt", header = TRUE, sep = "\t")

states_oi <- c(
  "Fetal Liver CS21-22",
  "Adult BM",
  "Fetal Liver CS17",
  "hPSC D31"
)

sample_table <- sample_table[sample_table$State %in% states_oi, ]


for (i in 1:dim(sample_table)[1]) {
  if (str_detect(string = sample_table$State[i], pattern = "Fetal Liver")) {
    sample_table$State[i] <- "Fetal Liver"
  }
  if (str_detect(string = sample_table$Cell_type[i], pattern = "IL7R")) {
    sample_table$Cell_type[i] <- "IL7R+"
  } else if (str_detect(string = sample_table$Cell_type[i], pattern = "LIN-")) {
    sample_table$Cell_type[i] <- "HSC-like"
  }
}

colnames(sample_table) <- c(
  "sample_name", "sample_type", "timing", "cell_type",
  "proposed.sample.names"
)

sample_colors <- c(
  "Adult BM" = "#7A0E0E",
  "Fetal Liver" = "#08298A",
  "hPSC D31" = "#21610B"
)

#read fpkm table and subset for samples of interest
geneexp <- read.table("fpkm.tsv", header = TRUE, sep = "\t")

rownames(geneexp) <- geneexp$X
geneexp$X <- NULL

geneexp <- geneexp[colnames(geneexp) %in% sample_table$proposed.sample.names]



# load genes of interest (GOI)
GOI_FL <- read.table(
  file = "Fetal_signature_p005.txt", sep = "\t") 

GOI_AD <- read.table(
  file = "Adult_signature_p005.txt",
  sep = "\t"
)
GOI_AD <- GOI_AD$V1
GOI_FL <- GOI_FL$V1

GOI <- data.frame(
  gene_name = append(GOI_FL, GOI_AD),
  FLorBM = append(rep("FL", length(GOI_FL)), rep("BM", length(GOI_AD))))



# remove genes with lower mean expression than the lowest GOI to avoid having
# non-expressed genes affect the data later on
geneexp <- geneexp[
  rowMeans(geneexp) >= min(rowMeans(geneexp[rownames(geneexp) %in%
    GOI$gene_name, ])),
]

geneexp <- log2(geneexp + 1)

rv <- rowVars(as.matrix(geneexp))

# select top 500 genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(geneexp[select, ]))
percent_var <- pca$sdev^2 / sum(pca$sdev^2)
genes_in_PCA <- rownames(geneexp[select, ])

pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)


ggplot(
  data = pcatable,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  coord_fixed() +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  xlim(-42, 55) +
  ylim(-43, 43)

ggsave('20221028_PCA_top500.pdf')


# PCA with fetal core genes removed from the top 500 variable genes

genes_in_top500 <- geneexp[select, ]

genes_use <- genes_in_top500[!(rownames(genes_in_top500) %in% GOI$gene_name), ]

pca <- prcomp(t(genes_use))

percent_var <- pca$sdev^2 / sum(pca$sdev^2)


pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)


ggplot(
  data = pcatable,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  coord_fixed() +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  xlim(-42, 55) +
  ylim(-43, 43)

ggsave('20221028_PCA_pseudorep_core_rem_after_top500.pdf')


# Check which genes overlap between top 500 variable and fetal core genes,
# should be 0
sum(rownames(genes_use) %in% GOI$gene_name)


# PCA with fetal core

select_core <- rownames(geneexp) %in% GOI$gene_name # find our gene sets

pca <- prcomp(t(geneexp[select_core, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2)

pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)


ggplot(
  data = pcatable,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  coord_fixed() +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  xlim(-42, 55) +
  ylim(-43, 43)

ggsave('20221028_PCA_pseudoreps_core.pdf')

overlapping_genes <- genes_in_PCA[genes_in_PCA %in% rownames(pca$rotation)]
write.csv(overlapping_genes, "20221028_overlapping_genes_pseudoreps_core.csv")


#random genes PCA

n_genes <- length(rownames(pca$rotation))

for (i in 1:5) {
  print(i)
  pca <- prcomp(t(geneexp[sample(nrow(geneexp), n_genes), ]))
  percent_var <- pca$sdev^2 / sum(pca$sdev^2)
  write.csv(rownames(pca$rotation),
    file =
      paste0("20221028_random_PCA_genes_pseudoreps_core_", i, ".csv")
  )
  pcatable <-
    data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      PC3 = pca$x[, 3],
      sample_table
    )
  p <-
    ggplot(
      data = pcatable,
      aes_string(
        x = "PC1",
        y = "PC2",
        color = "sample_type",
        shape = "cell_type"
      )
    ) +
    geom_point(size = 4) +
    theme_bw() +
    coord_fixed() +
    scale_color_manual(values = sample_colors) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
    xlim(-42, 55) +
    ylim(-43, 43)
  print(p)
  ggsave(paste0("20221028_PCA_random_pseudoreps_core_", i, ".pdf"))
}
