invisible(lapply(list(
  "stringr",
  "ggplot2",
  "stringr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

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

sample_table <- read.table(
  snakemake@input[["samples"]],
  header = TRUE,
  sep = "\t"
)

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

geneexp <- read.table(snakemake@input[["fpkm"]], header = TRUE, sep = "\t")

rownames(geneexp) <- geneexp$X
geneexp$X <- NULL
geneexp <- geneexp[colnames(geneexp) %in% sample_table$proposed.sample.names]

# remove genes with lower mean expression than the lowest GOI to avoid having
# non-expressed genes affect the data later on
geneexp <- geneexp[
  rowMeans(geneexp) >= min(rowMeans(geneexp[rownames(geneexp) %in%
    goi$gene_name, ])),
]

geneexp <- log2(geneexp + 1)

rv <- rowVars(as.matrix(geneexp))

# select top 500 genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(geneexp[select, ]))
percent_var <- pca$sdev^2 / sum(pca$sdev^2)
genes_in_pca <- rownames(geneexp[select, ])

pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)

p <- ggplot(
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
ggsave(snakemake@output[["plot_top500"]], plot = p)
# ggsave("20221028_PCA_top500.pdf", plot = p)

# PCA with fetal core genes removed from the top 500 variable genes
genes_in_top500 <- geneexp[select, ]
genes_use <- genes_in_top500[!(rownames(genes_in_top500) %in% goi$gene_name), ]
pca <- prcomp(t(genes_use))
percent_var <- pca$sdev^2 / sum(pca$sdev^2)

pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)

p <- ggplot(
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
ggsave(snakemake@output[["plot_subset500"]], plot = p)
# ggsave("20221028_PCA_pseudorep_core_rem_after_top500.pdf", plot = p)

# Check which genes overlap between top 500 variable and fetal core genes,
# should be 0
stopifnot(sum(rownames(genes_use) %in% goi$gene_name) == 0)

# PCA with fetal core
select_core <- rownames(geneexp) %in% goi$gene_name # find our gene sets

pca <- prcomp(t(geneexp[select_core, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2)

pcatable <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_table
)

p <- ggplot(
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

ggsave(snakemake@output[["plot_core"]], plot = p)
# ggsave("20221028_PCA_pseudoreps_core.pdf", plot = p)

overlapping_genes <- genes_in_pca[genes_in_pca %in% rownames(pca$rotation)]
write.csv(overlapping_genes, snakemake@output[["overlapping_genes"]])
# write.csv(overlapping_genes, "20221028_overlapping_genes_pseudoreps_core.csv")

#random genes PCA
n_genes <- length(rownames(pca$rotation))

for (i in 1:5) {
  print(i)
  pca <- prcomp(t(geneexp[sample(nrow(geneexp), n_genes), ]))
  percent_var <- pca$sdev^2 / sum(pca$sdev^2)
  write.csv(rownames(pca$rotation),
    file = snakemake@output[["csvs_random"]][[i]]
    # file = paste0("20221028_random_PCA_genes_pseudoreps_core_", i, ".csv")
  )
  pcatable <-
    data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      PC3 = pca$x[, 3],
      sample_table
    )
  p <- ggplot(
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
  ggsave(snakemake@output[["plots_random"]][[i]], plot = p)
  # ggsave(paste0("PCA_random_pseudoreps_core_", i, ".pdf"), plot = p)
}
