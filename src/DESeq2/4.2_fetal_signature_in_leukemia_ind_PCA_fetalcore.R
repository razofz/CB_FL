invisible(lapply(list(
  "stringr",
  "biomaRt",
  "limma",
  "ggplot2",
  "pheatmap",
  "stringr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

# images_dir <- dirname(snakemake@output[["plot_fl_core_mllaf4"]])

plots_pc1 <- snakemake@output[["plots_pc1"]]
plots_pc1_gender <- snakemake@output[["plots_pc1_gender"]]
plots_type <- snakemake@output[["plots_type"]]
plots_type_gender <- snakemake@output[["plots_type_gender"]]
leukemia_types <- snakemake@config[["leukemia_types"]]
print(leukemia_types)
print(plots_pc1)

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

named_plots_pc1 <- make_named_version(
  ids = leukemia_types,
  file_list = plots_pc1,
  pattern = "."
)
print(named_plots_pc1)
named_plots_pc1_gender <- make_named_version(
  ids = leukemia_types,
  file_list = plots_pc1_gender,
  pattern = "_"
)
print(named_plots_pc1_gender)
named_plots_type <- make_named_version(
  ids = leukemia_types,
  file_list = plots_type,
  pattern = "."
)
print(named_plots_type)
named_plots_type_gender <- make_named_version(
  ids = leukemia_types,
  file_list = plots_type_gender,
  pattern = "_"
)
print(named_plots_type_gender)

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

subtype_oi <- c("KMT2A", "ETV6-RUNX1", "Ph", "High hyperdiploid",
                "Low hyperdiploid")
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

# Load samples into DEseq2

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
  values = goi,
  mart = ensembl
)

type_colours <-
  c(
    "BCR-ABL1" = "#F8766D",
    "ETV6-RUNX1" = "#C49A00",
    "High hyperdiploid" = "#53B400",
    "KMT2A-AFF1" = "#00C094",
    "KMT2A-MLLT1" = "#00B6EB",
    "KMT2A-MLLT3" = "#A58AFF",
    "Low hyperdiploid" = "#FB61D7"
  )

age_colours <- c(
  "0-2" = "#F8766D",
  "2-16" = "#7CAE00",
  "16-40" = "#00BFC4",
  ">40" = "#C77CFF"
)

# IndividualPC1_fetalcore_KMT2A-MLLT1  IndividualPC1_fetalcore_BCR-ABL1
# IndividualPC1_fetalcore_KMT2A-AFF1  IndividualPC1_fetalcore_High hyperdiploid
# IndividualPC1_fetalcore_Low hyperdiploid  IndividualPC1_fetalcore_ETV6-RUNX1
# IndividualPC1_fetalcore_KMT2A-MLLT3

stopifnot(all(unique(vsd$type) %in% leukemia_types))

# subset vsd object and perform PCA individually using the GOI
for (x in 1:length(unique(vsd$type))) {
  print(x)
  print(unique(vsd$type)[x])
  vsd_sub <- vsd[, vsd$type %in% unique(vsd$type)[x]]
  plotPCA(vsd_sub, "type")

  # breaking the PCA analysis into its parts
  intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
  select <- rownames(vsd_sub) %in% ids$ensembl_gene_id # find the our gene sets
  pca <- prcomp(t(assay(vsd_sub)[select, ])) # calculate PCA
  percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
  intgroup_df <- as.data.frame(colData(vsd_sub)[, intgroup, drop = FALSE])
  # add grouping if wanted
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = ":"))
  } else {
    colData(vsd_sub)[[intgroup]]
  }
  princpial3 <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
    group = group, intgroup_df, name =
      colnames(vsd_sub)
  )

  ensm_to_sym <- data.frame(
    ensblID = rownames(pca$rotation),
    gene_sym = getBM(
      attributes = "hgnc_symbol", filters = "ensembl_gene_id",
      values = rownames(pca$rotation), mart = ensembl
    )
  )

  p <- ggplot(data = princpial3, aes_string(
    x = "PC1", y = "age", color =
      "age_group", shape = "gender"
  )) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    scale_color_manual(values = age_colours) +
    coord_fixed()
  print(p)
  ggsave(named_plots_pc1_gender[[unique(vsd$type)[x]]],
    # ggsave(paste0(
    #   images_dir,
    #   "/IndividualPC1_fetalcore_", unique(vsd$type)[x], "_gender.pdf"
    # ),
    plot = last_plot()
  )

  # KMT2A-AFF1 is Fig 5d, left, the others are fig S6c, upper row
  p <- ggplot(data = princpial3, aes_string(
    x = "PC1", y = "age", color =
      "age_group"
  )) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    scale_color_manual(values = age_colours) +
    coord_fixed()
  print(p)
  ggsave(named_plots_pc1[[unique(vsd$type)[x]]],
    # ggsave(paste0(
    #   images_dir, "/",
    #   "/IndividualPC1_fetalcore_", unique(vsd$type)[x], ".pdf"
    # ),
    plot = last_plot()
  )
}

# PCA on all samples combined
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
select <- rownames(vsd) %in% ids$ensembl_gene_id # find the our gene sets
pca <- prcomp(t(assay(vsd)[select, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
# add grouping if wanted
intgroup_df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = ":"))
} else {
  colData(vsd)[[intgroup]]
}
principal <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  group = group, intgroup_df, name = colnames(vsd)
)

ensm_to_sym <- data.frame(
  ensblID = rownames(pca$rotation),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(pca$rotation), mart = ensembl
  )
)

fl_or_bm <- NULL
for (i in 1:length(ensm_to_sym$hgnc_symbol)) {
  x <- goi$fl_or_bm[goi$gene_name %in%
    as.character(lapply(ensm_to_sym$hgnc_symbol, toupper))[i]]
  fl_or_bm[i] <- x
}

loading <- data.frame(
  LO1 = pca$rotation[, 1], LO2 = pca$rotation[, 2], LO3 = pca$rotation[, 3],
  name = ensm_to_sym$ensblID, symb = ensm_to_sym$hgnc_symbol,
  origin = fl_or_bm
)


f <- ggplot(data = principal, aes_string(
  x = "PC1", y = "age", color = "type",
  shape = "gender"
)) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("Age")) +
  theme_bw() +
  theme(aspect.ratio = 1.5) +
  scale_color_manual(values = type_colours) +
  xlim(-38, 38) +
  ylim(-1, 82) +
  coord_fixed()
ggsave(snakemake@output[["plot_fl_core_same_pc1_gender"]],
  # ggsave(str_c(images_dir, "/fetalcore_samePC1_all_gender.pdf"),
  plot = f
)

# fig 5b, left
f <- ggplot(data = principal, aes_string(
  x = "PC1", y = "age",
  color = "type"
)) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("Age")) +
  theme_bw() +
  theme(aspect.ratio = 1.5) +
  scale_color_manual(values = type_colours) +
  xlim(-38, 38) +
  ylim(-1, 82) +
  coord_fixed()
ggsave(snakemake@output[["plot_fl_core_same_pc1"]], plot = f)
# ggsave(str_c(images_dir, "/fetalcore_samePC1_all.pdf"), plot = last_plot())

# loop over the subtypes and plot them individually in the same PCA space
for (x in 1:length(unique(principal$type))) {
  print(x)
  f <- ggplot(
    data = subset.data.frame(principal, type == unique(principal$type)[x]),
    aes_string(x = "PC1", y = "age", color = "type", shape = "gender")
  ) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    scale_color_manual(values = type_colours) +
    coord_fixed()
  print(f)
  ggsave(named_plots_type_gender[[unique(principal$type)[x]]],
    plot = last_plot()
  )
  # ggsave(paste0(
  #   images_dir, "/",
  #   unique(principal$type)[x], "_gender.pdf"
  # ), plot = last_plot())

  f <- ggplot(
    data = subset.data.frame(principal, type == unique(principal$type)[x]),
    aes_string(x = "PC1", y = "age", color = "type")
  ) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    scale_color_manual(values = type_colours) +
    coord_fixed()
  print(f)
  ggsave(named_plots_type[[unique(principal$type)[x]]],
    plot = last_plot()
  )
  # ggsave(paste0(
  #   images_dir, "/",
  #   unique(principal$type)[x], ".pdf"
  # ), plot = last_plot())
}

# subset vsd object and perform PCA individually on MLL-AF4 using the GOI
vsd_sub <- vsd[, vsd$type %in% "KMT2A-AFF1"]
# plotPCA(vsd_sub, "type")

# breaking the PCA analysis into its parts
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
select <- rownames(vsd_sub) %in% ids$ensembl_gene_id # find the our gene sets
pca <- prcomp(t(assay(vsd_sub)[select, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
 # add grouping if wanted
intgroup_df <- as.data.frame(colData(vsd_sub)[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = ":"))
} else {
  colData(vsd_sub)[[intgroup]]
}
principal <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  group = group, intgroup_df, name = colnames(vsd_sub)
)


p2 <- ggplot(data = principal, aes_string(
  x = "PC1", y = "PC2", color =
    "age_group", shape = "gender"
)) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1.5) +
  xlim(-25, 27) +
  ylim(-24, 22) +
  coord_fixed()
ggsave(snakemake@output[["plot_fl_core_mllaf4_gender"]],
  # ggsave(paste0(images_dir, "/IndividualPC1_PC2_fetalcore_MLLAF4_gender.pdf"),
  plot = p2
)

# Fig S6b
p2 <- ggplot(data = principal, aes_string(
  x = "PC1", y = "PC2", color = "age_group"
)) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  theme_bw() +
  theme(aspect.ratio = 1.5) +
  xlim(-25, 27) +
  ylim(-24, 22) +
  coord_fixed()
ggsave(snakemake@output[["plot_fl_core_mllaf4"]],
  # ggsave(paste0(images_dir, "/IndividualPC1_PC2_fetalcore_MLLAF4.pdf"),
  plot = p2
)
