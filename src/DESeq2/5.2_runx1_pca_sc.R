invisible(lapply(list(
  "stringr",
  "ggplot2",
  "stringr",
  "dplyr",
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
pca_1 <- prcomp(t(geneexp[select, ]))
percent_var_1 <- pca_1$sdev^2 / sum(pca_1$sdev^2)
genes_in_pca_1 <- rownames(geneexp[select, ])

pca_1table <- data.frame(
  PC1 = pca_1$x[, 1], PC2 = pca_1$x[, 2], PC3 = pca_1$x[, 3],
  sample_table
)

# PCA with fetal core genes removed from the top 500 variable genes
genes_in_top500 <- geneexp[select, ]
genes_use <- genes_in_top500[!(rownames(genes_in_top500) %in% goi$gene_name), ]
pca_2 <- prcomp(t(genes_use))
percent_var_2 <- pca_2$sdev^2 / sum(pca_2$sdev^2)

pca_2table <- data.frame(
  PC1 = pca_2$x[, 1], PC2 = pca_2$x[, 2], PC3 = pca_2$x[, 3],
  sample_table
)

# Check which genes overlap between top 500 variable and fetal core genes,
# should be 0
stopifnot(sum(rownames(genes_use) %in% goi$gene_name) == 0)

# PCA with fetal core
select_core <- rownames(geneexp) %in% goi$gene_name # find our gene sets

pca_3 <- prcomp(t(geneexp[select_core, ])) # calculate PCA
percent_var_3 <- pca_3$sdev^2 / sum(pca_3$sdev^2)

pca_3table <- data.frame(
  PC1 = pca_3$x[, 1], PC2 = pca_3$x[, 2], PC3 = pca_3$x[, 3],
  sample_table
)

overlapping_genes <- genes_in_pca_1[genes_in_pca_1 %in% rownames(pca_3$rotation)]
write.csv(overlapping_genes, snakemake@output[["overlapping_genes"]])

min_x <- min(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
  FUN = function(x) min(x["PC1"])
)))
min_y <- min(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
  FUN = function(x) min(x["PC2"])
)))
max_x <- max(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
  FUN = function(x) max(x["PC1"])
)))
max_y <- max(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
  FUN = function(x) max(x["PC2"])
)))
width_x <- max_x - min_x
width_y <- max_y - min_y
margin_x <- width_x * .05
margin_y <- width_y * .05

p1 <- ggplot(
  data = pca_1table,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1 / 1) +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var_1[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var_1[2] * 100), "% variance")) +
  xlim(min_x - margin_x, max_x + margin_x) +
  ylim(min_y - margin_y, max_y + margin_y)
  # xlim(-42, 55) +
  # ylim(-43, 43)
p2 <- ggplot(
  data = pca_2table,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1 / 1) +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var_2[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var_2[2] * 100), "% variance")) +
  xlim(min_x - margin_x, max_x + margin_x) +
  ylim(min_y - margin_y, max_y + margin_y)
p3 <- ggplot(
  data = pca_3table,
  aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
) +
  geom_point(size = 4) +
  theme_bw() +
  theme(aspect.ratio = 1 / 1) +
  coord_fixed() +
  scale_color_manual(values = sample_colors) +
  xlab(paste0("PC1: ", round(percent_var_3[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var_3[2] * 100), "% variance")) +
  xlim(min_x - margin_x, max_x + margin_x) +
  ylim(min_y - margin_y, max_y + margin_y)

ggsave(snakemake@output[["plot_top500"]], plot = p1)
ggsave(snakemake@output[["plot_core"]], plot = p3)
ggsave(snakemake@output[["plot_subset500"]], plot = p2)

#random genes PCA
n_genes <- length(rownames(pca_3$rotation))

for (i in 1:5) {
  print(str_c("> i:", i))
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
    theme(aspect.ratio = 1 / 1) +
    scale_color_manual(values = sample_colors) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
    xlim(min_x - margin_x, max_x + margin_x) +
    ylim(min_y - margin_y, max_y + margin_y)
    # xlim(-42, 55) +
    # ylim(-43, 43)
  ggsave(snakemake@output[["plots_random"]][[i]], plot = p)
  # ggsave(paste0("PCA_random_pseudoreps_core_", i, ".pdf"), plot = p)
}


######################################################
# Do this for the extra cutoffs as well
######################################################

get_cutoff_files <- function(struct, cutoff) {
  struct[
    str_detect(
      pattern = cutoff %>% as.character(),
      struct
    )
  ]
}

for (cutoff in
     snakemake@config[["seurat_deg_cutoffs"]][["log2foldchange"]][["extra"]]) {
  print(str_c("> cutoff: ", cutoff))
  # load genes of interest (GOI)
  goi_fl <- get_cutoff_files(
    struct = snakemake@input[["fetal_signature_extra"]],
    cutoff = cutoff
  )
  goi_ad <- get_cutoff_files(
    struct = snakemake@input[["adult_signature_extra"]],
    cutoff = cutoff
  )

  stopifnot(all(c(
    !is.null(goi_fl),
    !is.null(goi_ad)
  )))

  goi_fl <- read.table(file = goi_fl, sep = "\t")
  goi_ad <- read.table(file = goi_ad, sep = "\t")
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
    } else if (
      str_detect(string = sample_table$Cell_type[i], pattern = "LIN-")) {
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
  pca_1 <- prcomp(t(geneexp[select, ]))
  percent_var_1 <- pca_1$sdev^2 / sum(pca_1$sdev^2)
  genes_in_pca_1 <- rownames(geneexp[select, ])

  pca_1table <- data.frame(
    PC1 = pca_1$x[, 1], PC2 = pca_1$x[, 2], PC3 = pca_1$x[, 3],
    sample_table
  )

  # PCA with fetal core genes removed from the top 500 variable genes
  genes_in_top500 <- geneexp[select, ]
  genes_use <- genes_in_top500[
    !(rownames(genes_in_top500) %in% goi$gene_name),
  ]
  pca_2 <- prcomp(t(genes_use))
  percent_var_2 <- pca_2$sdev^2 / sum(pca_2$sdev^2)

  pca_2table <- data.frame(
    PC1 = pca_2$x[, 1], PC2 = pca_2$x[, 2], PC3 = pca_2$x[, 3],
    sample_table
  )

  # Check which genes overlap between top 500 variable and fetal core genes,
  # should be 0
  stopifnot(sum(rownames(genes_use) %in% goi$gene_name) == 0)

  # PCA with fetal core
  select_core <- rownames(geneexp) %in% goi$gene_name # find our gene sets

  pca_3 <- prcomp(t(geneexp[select_core, ])) # calculate PCA
  percent_var_3 <- pca_3$sdev^2 / sum(pca_3$sdev^2)

  pca_3table <- data.frame(
    PC1 = pca_3$x[, 1], PC2 = pca_3$x[, 2], PC3 = pca_3$x[, 3],
    sample_table
  )

  overlapping_genes <- genes_in_pca_1[
    genes_in_pca_1 %in% rownames(pca_3$rotation)
  ]
  filename <- get_cutoff_files(
    struct = snakemake@output[["overlapping_genes_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  write.csv(overlapping_genes, filename)

  min_x <- min(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
    FUN = function(x) min(x["PC1"])
  )))
  min_y <- min(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
    FUN = function(x) min(x["PC2"])
  )))
  max_x <- max(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
    FUN = function(x) max(x["PC1"])
  )))
  max_y <- max(unlist(lapply(list(pca_1table, pca_2table, pca_3table),
    FUN = function(x) max(x["PC2"])
  )))
  width_x <- max_x - min_x
  width_y <- max_y - min_y
  margin_x <- width_x * .05
  margin_y <- width_y * .05

  p1 <- ggplot(
    data = pca_1table,
    aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
  ) +
    geom_point(size = 4) +
    theme_bw() +
    coord_fixed() +
    theme(aspect.ratio = 1 / 1) +
    scale_color_manual(values = sample_colors) +
    xlab(paste0("PC1: ", round(percent_var_1[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var_1[2] * 100), "% variance")) +
    xlim(min_x - margin_x, max_x + margin_x) +
    ylim(min_y - margin_y, max_y + margin_y)
    # xlim(-42, 55) +
    # ylim(-43, 43)
  p2 <- ggplot(
    data = pca_2table,
    aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
  ) +
    geom_point(size = 4) +
    theme_bw() +
    coord_fixed() +
    theme(aspect.ratio = 1 / 1) +
    scale_color_manual(values = sample_colors) +
    xlab(paste0("PC1: ", round(percent_var_2[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var_2[2] * 100), "% variance")) +
    xlim(min_x - margin_x, max_x + margin_x) +
    ylim(min_y - margin_y, max_y + margin_y)
  p3 <- ggplot(
    data = pca_3table,
    aes_string(x = "PC1", y = "PC2", color = "sample_type", shape = "cell_type")
  ) +
    geom_point(size = 4) +
    theme_bw() +
    theme(aspect.ratio = 1 / 1) +
    coord_fixed() +
    scale_color_manual(values = sample_colors) +
    xlab(paste0("PC1: ", round(percent_var_3[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percent_var_3[2] * 100), "% variance")) +
    xlim(min_x - margin_x, max_x + margin_x) +
    ylim(min_y - margin_y, max_y + margin_y)

  filename <- get_cutoff_files(
    struct = snakemake@output[["plot_top500_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  ggsave(filename, plot = p1)

  filename <- get_cutoff_files(
    struct = snakemake@output[["plot_core_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  ggsave(filename, plot = p3)

  filename <- get_cutoff_files(
    struct = snakemake@output[["plot_subset500_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  ggsave(filename, plot = p2)


  #random genes PCA
  n_genes <- length(rownames(pca_3$rotation))
  

  # super annoying to have to do this, a result of R's implicit conversion of
  # integer doubles to actual integers (e.g. 1.0 becomes 1).
  formatted_cutoff <- ""
  if (cutoff %>% as.character() %>% str_length() < 3) {
    formatted_cutoff <- str_c(
      cutoff %>% as.character(),
      ".0"
    )
    print(str_c("> formatted_cutoff: ", formatted_cutoff))
  }
  
  plots <- NULL
  plots <- get_cutoff_files(
    struct = snakemake@output[["plots_random_extra"]],
    cutoff = formatted_cutoff
  )
  stopifnot(!is.null(plots))
  # print(plots)

  csvs <- NULL
  csvs <- get_cutoff_files(
    struct = snakemake@output[["csvs_random_extra"]],
    cutoff = formatted_cutoff
  )
  stopifnot(!is.null(csvs))
  # print(csvs)

  for (i in 1:5) {
    print(str_c("> i: ", i))
    pca <- prcomp(t(geneexp[sample(nrow(geneexp), n_genes), ]))
    percent_var <- pca$sdev^2 / sum(pca$sdev^2)
    # print(csvs[[i]])
    write.csv(rownames(pca$rotation),
      file = csvs[[i]]
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
      theme(aspect.ratio = 1 / 1) +
      scale_color_manual(values = sample_colors) +
      xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
      xlim(min_x - margin_x, max_x + margin_x) +
      ylim(min_y - margin_y, max_y + margin_y)
    # xlim(-42, 55) +
    # ylim(-43, 43)
    # print(plots[[i]])
    ggsave(filename = plots[[i]], plot = p)
    # ggsave(paste0("PCA_random_pseudoreps_core_", i, ".pdf"), plot = p)
  }
}
