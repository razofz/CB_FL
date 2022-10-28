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

sample_table <- read.table(str_c(
  external_dir,
  "/BALL-1988S-HTSeq/B-ALL-subtyping.txt"
),
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
    path = str_c(external_dir, "/BALL-1988S-HTSeq/"),
    pattern = as.character(sample_table$patient[i]),
    ignore.case = TRUE, full.names = TRUE
  )[1] 
  # hard-coded to only take the first sample,
  # works for this one but might not for others
}
sample_table$file_name <- file_names

# Load samples into DEseq2

DEsampleTable <- data.frame(
  sampleNames = sample_table$patient,
  fileName = as.character(sample_table$file_name),
  age_group = sample_table$age_group,
  type = sample_table$fusion,
  age = sample_table$age, gender = sample_table$gender,
  RNA_lib = sample_table$RNA.seq.library
)

counts <- DESeqDataSetFromHTSeqCount(
  sampleTable = DEsampleTable, design =
    ~RNA_lib
)
counts <- DESeq(counts)
vsd <- vst(counts, blind = FALSE)
# plotPCA(vsd3, "RNA_lib") + theme_bw() + coord_fixed(ratio = 8 / 6.3)

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)
# plotPCA(vsd3, "RNA_lib") + theme_bw() + coord_fixed(ratio = 5.6 / 7.2)

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
  values = GOI,
  # values = GOI$gene_name,
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

# IndividualPC1_fetalcore_KMT2A-MLLT1  IndividualPC1_fetalcore_BCR-ABL1  IndividualPC1_fetalcore_KMT2A-AFF1  IndividualPC1_fetalcore_High hyperdiploid  IndividualPC1_fetalcore_Low hyperdiploid  IndividualPC1_fetalcore_ETV6-RUNX1  IndividualPC1_fetalcore_KMT2A-MLLT3

# subset vsd3 object and perform PCA individually using the GOI
for (x in 1:length(unique(vsd$type))) {
  print(x)
  print(unique(vsd$type)[x])
  vsd3_sub <- vsd[, vsd$type %in% unique(vsd$type)[x]]
  plotPCA(vsd3_sub, "type")

  # breaking the PCA analysis into its parts
  intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
  select <- rownames(vsd3_sub) %in% ids$ensembl_gene_id # find the our gene sets
  pca <- prcomp(t(assay(vsd3_sub)[select, ])) # calculate PCA
  percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
  intgroup_df <- as.data.frame(colData(vsd3_sub)[, intgroup, drop = FALSE])
  # add grouping if wanted
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = ":"))
  } else {
    colData(vsd3_sub)[[intgroup]]
  }
  princpial3 <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
    group = group, intgroup_df, name =
      colnames(vsd3_sub)
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
  ggsave(paste0(
    images_dir,
    "/IndividualPC1_fetalcore_", unique(vsd$type)[x], "_gender.pdf"
  ),
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
  ggsave(paste0(
    images_dir, "/",
    "/IndividualPC1_fetalcore_", unique(vsd$type)[x], ".pdf"
  ),
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
princpial3 <- data.frame(
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

FLor_BM <- NULL
for (i in 1:length(ensm_to_sym$hgnc_symbol)) {
  x <- GOI$FLorBM[GOI$gene_name %in%
    as.character(lapply(ensm_to_sym$hgnc_symbol, toupper))[i]]
  FLor_BM[i] <- x
}
# for (i in 1:length(ensm_to_sym$hgnc_symbol)) {
#   x <- GOI$FLorBM[GOI$gene_name %in%
#     as.character(lapply(ensm_to_sym$hgnc_symbol, toupper))[i]]
#   FLor_BM[i] <- x
# }

loading <- data.frame(
  LO1 = pca$rotation[, 1], LO2 = pca$rotation[, 2], LO3 = pca$rotation[, 3],
  name = ensm_to_sym$ensblID, symb = ensm_to_sym$hgnc_symbol,
  origin = FLor_BM
)


f <- ggplot(data = princpial3, aes_string(
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
print(f)
ggsave(str_c(images_dir, "/fetalcore_samePC1_all_gender.pdf"),
  plot =
    last_plot()
)

# fig 5b, left
f <- ggplot(data = princpial3, aes_string(
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
print(f)
ggsave(str_c(images_dir, "/fetalcore_samePC1_all.pdf"), plot = last_plot())


# loop over the subtypes and plot them individually in the same PCA space
for (x in 1:length(unique(princpial3$type))) {
  print(x)
  f <- ggplot(
    data = subset.data.frame(princpial3, type == unique(princpial3$type)[x]),
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
  ggsave(paste0(
    images_dir, "/",
    unique(princpial3$type)[x], "_gender.pdf"
  ), plot = last_plot())

  f <- ggplot(
    data = subset.data.frame(princpial3, type == unique(princpial3$type)[x]),
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
  ggsave(paste0(
    images_dir, "/",
    unique(princpial3$type)[x], ".pdf"
  ), plot = last_plot())
}


# subset vsd3 object and perform PCA individually on MLL-AF4 using the GOI
vsd3_sub <- vsd[, vsd$type %in% "KMT2A-AFF1"]
plotPCA(vsd3_sub, "type")

# breaking the PCA analysis into its parts
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
select <- rownames(vsd3_sub) %in% ids$ensembl_gene_id # find the our gene sets
pca <- prcomp(t(assay(vsd3_sub)[select, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
 # add grouping if wanted
intgroup_df <- as.data.frame(colData(vsd3_sub)[, intgroup, drop = FALSE])
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = ":"))
} else {
  colData(vsd3_sub)[[intgroup]]
}
princpial3 <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  group = group, intgroup_df, name = colnames(vsd3_sub)
)


p2 <- ggplot(data = princpial3, aes_string(
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
print(p2)
ggsave(paste0(images_dir, "/IndividualPC1_PC2_fetalcore_MLLAF4_gender.pdf"),
  plot = last_plot()
)

# Fig S6b
p2 <- ggplot(data = princpial3, aes_string(
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
print(p2)
ggsave(paste0(images_dir, "/IndividualPC1_PC2_fetalcore_MLLAF4.pdf"),
  plot = last_plot()
)
