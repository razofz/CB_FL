invisible(lapply(list(
  "stringr",
  "DESeq2",
  "biomaRt",
  "limma",
  "ggplot2",
  "fdrtool",
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

types_oi <- snakemake@params[["types_oi"]]
print(types_oi)
# leukemia_types <- snakemake@config[["leukemia_types"]]
# print(leukemia_types)

named_plots_up <- make_named_version(
  ids = types_oi,
  file_list = snakemake@output[["plots_up"]],
  pattern = "\\."
)
named_plots_down <- make_named_version(
  ids = types_oi,
  file_list = snakemake@output[["plots_down"]],
  pattern = "\\."
)

sample_table <- read.table(snakemake@input[["ball_tsv"]],
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
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)

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
ids_fl <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = goi_fl,
  mart = ensembl
)
ids_bm <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = goi_ad,
  mart = ensembl
)
select_fl_core <- rownames(vsd) %in% ids$ensembl_gene_id
select_fl_core_fl <- rownames(vsd) %in% ids_fl$ensembl_gene_id
select_fl_core_bm <- rownames(vsd) %in% ids_bm$ensembl_gene_id

age_group_col <- c(
  "0-2" = "#F8766D", "2-16" = "#7CAE00",
  "16-40" = "#00BFC4", ">40" = "#C77CFF"
)
# types_oi <- c("KMT2A-AFF1", "BCR-ABL1", "High hyperdiploid")
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib")

sig_levels_fl <- NULL
sig_levels_bm <- NULL
vsd_sub_fl <- NULL
vsd_sub_m <- NULL

# Fetal_core_up_box_plot_KMT2A-AFF1
# Fetal_core_up_box_plot_BCR-ABL1
# Fetal_core_up_box_plot_High hyperdiploid
# Fetal_core_down_box_plot_KMT2A-AFF1
# Fetal_core_down_box_plot_BCR-ABL1
# Fetal_core_down_box_plot_High hyperdiploid

# fig5e (all of them)
for (leuk in types_oi) {
  vsd_sub <- vsd[, vsd$type %in% leuk]
  # add grouping if wanted
  intgroup_df <- as.data.frame(colData(vsd_sub)[, intgroup, drop = FALSE])

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = ":"))
  } else {
    colData(vsd_sub)[[intgroup]]
  }

  vsd_sub_fl <- as.data.frame(t(assay(vsd_sub)[select_fl_core_fl, ]))
  vsd_sub_fl <- sapply(vsd_sub_fl, function(vsd_sub_fl) {
    (vsd_sub_fl -
      mean(vsd_sub_fl)) /
      sd(vsd_sub_fl)
  })
  means_fl <- rowMeans(vsd_sub_fl)
  vsd_sub_fl <- cbind.data.frame(intgroup_df, means_fl)
  p1 <- ggplot(data = vsd_sub_fl, aes_string(
    y = "means_fl", x = "age_group",
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
  # print(p1)
  ggsave(
    named_plots_up[[leuk]],
    # paste0(images_dir, "/Fetal_core_up_box_plot_", leuk, ".pdf"),
    plot = p1
  )
  vsd_sub_bm <- as.data.frame(t(assay(vsd_sub)[select_fl_core_bm, ]))
  vsd_sub_bm <- sapply(vsd_sub_bm, function(vsd_sub_bm) {
    (vsd_sub_bm -
      mean(vsd_sub_bm)) /
      sd(vsd_sub_bm)
  })
  means_bm <- rowMeans(vsd_sub_bm)
  vsd_sub_bm <- cbind.data.frame(intgroup_df, means_bm)
  p2 <- ggplot(data = vsd_sub_bm, aes_string(
    y = "means_bm", x = "age_group",
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

  # print(p2)
  ggsave(
    named_plots_down[[leuk]],
    # paste0(images_dir, "/Fetal_core_down_box_plot_", leuk, ".pdf"),
    plot = p2
  )

  res_aov_fl <- aov(means_fl ~ age_group, data = vsd_sub_fl)
  sig_levels_fl <- rbind.data.frame(
    sig_levels_fl,
    TukeyHSD(res_aov_fl)[["age_group"]][, 4]
  )

  res_aov_bm <- aov(means_bm ~ age_group, data = vsd_sub_bm)
  sig_levels_bm <- rbind.data.frame(
    sig_levels_bm,
    TukeyHSD(res_aov_bm)[["age_group"]][, 4]
  )
}

colnames(sig_levels_fl) <- rownames(TukeyHSD(res_aov_fl)[["age_group"]])
rownames(sig_levels_fl) <- types_oi

colnames(sig_levels_bm) <- rownames(TukeyHSD(res_aov_bm)[["age_group"]])
rownames(sig_levels_bm) <- types_oi

write.csv(
  sig_levels_fl,
  file = snakemake@output[["pvals_up"]],
  # file = paste0(out_dir, "/Pvals_anova_FL_up2.csv"),
  quote = F
)
write.csv(
  sig_levels_bm,
  file = snakemake@output[["pvals_down"]],
  # file = paste0(out_dir, "/Pvals_anova_FL_down2.csv"),
  quote = F
)
