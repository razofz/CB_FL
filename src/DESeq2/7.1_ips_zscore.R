invisible(lapply(list(
  "stringr",
  "ggplot2",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sample_table <- read.table(
  snakemake@input[["samples"]],
  header = TRUE,
  sep = "\t"
)
geneexp <- read.table(snakemake@input[["fpkm"]], header = TRUE, sep = "\t")

states_oi <- c(
  "Fetal Liver CS21-22",
  "Adult BM",
  "Fetal Liver CS17",
  "hPSC D31",
  "ETV6RUNX1 hIPS"
)

sample_table <- sample_table[sample_table$State %in% states_oi, ]

# rename some sample and cell types
for (i in 1:dim(sample_table)[1]) {
  if (str_detect(string = sample_table$State[i], pattern = "Fetal Liver")) {
    sample_table$State[i] <- "Fetal Liver"
  } else if (str_detect(string = sample_table$State[i], pattern = "hPSC D31")) {
    sample_table$State[i] <- "IPSCs"
  } else if (
    str_detect(string = sample_table$State[i], pattern = "ETV6RUNX1 hIPS")
  ) {
    sample_table$State[i] <- "ETV6-RUNX1 IPSCs"
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

sample_table$sample_type <- factor(sample_table$sample_type, levels = c(
  "Adult BM",
  "Fetal Liver",
  "IPSCs",
  "ETV6-RUNX1 IPSCs"
))

annotation_colors <- list(
  "sample_type" = c(
    "Adult BM" = "#7A0E0E",
    "Fetal Liver" = "#08298A",
    "IPSCs" = "#21610B",
    "ETV6-RUNX1 IPSCs" = "#4C0B5F"
  ), "cell_type" = c(
    "HSC-like" = "#D8D8D8",
    "IL7R+" = "#585858",
    "ProB" = "#151515"
  )
)

rownames(geneexp) <- geneexp$X
geneexp$X <- NULL
geneexp <- geneexp[colnames(geneexp) %in% sample_table$proposed.sample.names]

# load genes of interest (GOI)

goi_fl <- read.table(
  file = snakemake@input[["deg"]],
  sep = ","
)
goi_fl <- goi_fl$FL_core_pos
goi_fl <- goi_fl[!is.na(goi_fl)]

# subset for GOI, log-transform and scale the fpkm table
geneexp <- geneexp[rownames(geneexp) %in% goi_fl, ]

geneexp <- log2(geneexp + 1)
geneexp <- data.frame(t(geneexp))
geneexp <- scale(geneexp, center = TRUE, scale = TRUE)
geneexp <- data.frame(geneexp)

# calculate mean z score of all the genes per sample
geneexp$mean_expr <- rowMeans(geneexp)

sample_table <- cbind(sample_table, geneexp$mean_expr)

colnames(sample_table)[6] <- "mean_FL_score"

ggplot(data = sample_table, aes_string(
  y = "mean_FL_score", x = "sample_type",
  fill = "cell_type"
)) +
  geom_boxplot(width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = annotation_colors[["cell_type"]]) +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab("Z-scored mean expression") +
  xlab("") +
  ggtitle("Fetal signature - universally upregulated")

ggsave(snakemake@output[["plot"]])
