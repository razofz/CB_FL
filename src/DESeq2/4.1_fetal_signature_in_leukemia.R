# loading KMT2A-MLLT1, KMT2A-MLLT3, ETV6-RUNX1, BCR-ABL1
source("env.R")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("limma"))

external_dir <- "data/external"
raw_dir <- "data/raw/paper_specific/Supplemental_tables"
out_dir <- "data/processed/DEseq2"
images_dir <- paste0(out_dir, "/images")
if (dir.exists(images_dir) == FALSE) {
  dir.create(images_dir, showWarnings = TRUE, recursive = TRUE)
}

if (dir.exists(paste0(external_dir, "/BALL-1988S-HTSeq")) == FALSE) {
  stop("BALL dir not found")
}

# load genes of interest (GOI)
GOI_FL <- read.table(
  file = paste0(raw_dir, "/Fetal_signature_p005.txt"),
  sep = "\t"
)
GOI_AD <- read.table(
  file = paste0(raw_dir, "/Adult_signature_p005.txt"),
  sep = "\t"
)
GOI_FL <- GOI_FL$V1
GOI_AD <- GOI_AD$V1
GOI <- data.frame(
  gene_name = append(GOI_FL, GOI_AD),
  FLorBM = append(rep("FL", length(GOI_FL)), rep("BM", length(GOI_AD)))
)

sample_table3 <- read.table(
  paste0(external_dir, "/BALL-1988S-HTSeq/subtypes.tsv"),
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
) # , 'High hyperdiploid', 'Low hyperdiploid'
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
    path = paste0(external_dir, "/BALL-1988S-HTSeq"),
    pattern = as.character(sample_table3$patient[i]),
    ignore.case = TRUE, full.names = TRUE
  )[1] # hard-coded to only take the first sample, works for this one but might not for others
}
sample_table3$file_name <- file_names3

# Load samples into DEseq2
DEsampleTable3 <- data.frame(
  sampleNames = sample_table3$patient,
  fileName = as.character(sample_table3$file_name),
  age_group = sample_table3$age_group,
  type = sample_table3$fusion,
  age = sample_table3$age, gender = sample_table3$gender,
  RNA_lib = sample_table3$RNA.seq.library
)

counts3 <- DESeqDataSetFromHTSeqCount(
  sampleTable = DEsampleTable3, design = ~RNA_lib
)
counts3 <- DESeq(counts3)

vsd3 <- vst(counts3, blind = FALSE)
plotPCA(vsd3, "RNA_lib") + theme_bw()
ggsave(paste0(images_dir, "/PC12_top500_preBatch_ALL_leuk.pdf"),
  plot = last_plot()
)
assay(vsd3) <- limma::removeBatchEffect(assay(vsd3), vsd3$RNA_lib)
plotPCA(vsd3, "RNA_lib") + theme_bw()
ggsave(paste0(images_dir, "/PC12_top500_postBatch_ALL_leuk.pdf"),
  plot = last_plot()
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


# PCA with our FL/BM genes
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
select <- rownames(vsd3) %in% ids$ensembl_gene_id # find the our gene sets
length(select)
table(select)
pca <- prcomp(t(assay(vsd3)[select, ])) # calculate PCA
percentVar <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
intgroup.df <- as.data.frame(colData(vsd3)[, intgroup, drop = FALSE]) 
# add grouping if wanted
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = ":"))
} else {
  colData(vsd3)[[intgroup]]
}
princpial3 <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  group = group, intgroup.df, name = colnames(vsd3)
)

Ensm_to_sym <- data.frame(
  ensblID = rownames(pca$rotation),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(pca$rotation), mart = ensembl
  )
)

FLor_BM <- NULL
for (i in 1:length(Ensm_to_sym$hgnc_symbol)) {
  x <- GOI$FLorBM[GOI$gene_name %in% as.character(lapply(Ensm_to_sym$hgnc_symbol, toupper))[i]]
  FLor_BM[i] <- x
}

loading <- data.frame(
  LO1 = pca$rotation[, 1], LO2 = pca$rotation[, 2], LO3 = pca$rotation[, 3],
  name = Ensm_to_sym$ensblID, symb = Ensm_to_sym$hgnc_symbol,
  origin = FLor_BM
)

ggplot(data = princpial3, aes_string(x = "PC1", y = "age", color = "type")) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  xlim(-35, 35) +
  ylim(-1, 82) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(paste0(images_dir, "/PC1Age_All_leuk_FL_core.pdf"), plot = last_plot())

ggplot(data = princpial3, aes_string(x = "PC2", y = "age", color = "type")) +
  geom_point(size = 3) +
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(paste0(images_dir, "/PC2Age_All_leuk_FL_core.pdf"), plot = last_plot())


for (x in 1:length(unique(princpial3$type))) {
  print(x)
  f <- ggplot(
    data = subset.data.frame(princpial3, type == unique(princpial3$type)[x]),
    aes_string(x = "PC1", y = "age", color = "type")
  ) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-35, 35) +
    ylim(-1, 82) +
    coord_fixed()
  print(f)
  ggsave(
    paste0(paste0(
      images_dir, "/PC1Age_", unique(princpial3$type)[x],
      "_FL_core.pdf"
    )),
    plot = last_plot()
  )
}

ngenes <- 30
# the genes contributing most to fetal signature
FL_up <- head(loading[order(loading$LO1, decreasing = TRUE), ], ngenes)
sum(str_count(FL_up$origin, pattern = "FL")) / ngenes * 100
# the genes contributing most to adult signature
FL_dn <- head(loading[order(loading$LO1, decreasing = FALSE), ], ngenes)
sum(str_count(FL_dn$origin, pattern = "BM")) / ngenes * 100
write.csv(x = FL_dn, file = paste0(out_dir, "/results/down_inf_luek.csv"))
write.csv(x = FL_up, file = paste0(out_dir, "/results/up_inf_luek.csv"))

Categories_to_show <- as.data.frame(colData(vsd3)[, c(
  "age_group", "gender", "type"
)])
Categories_color <- list(
  gender = c("Male" = "#F8766D", "Female" = "#7CAE00"),
  age_group = c(
    "0-2" = "#F8766D", "2-16" = "#7CAE00", "16-40" = "#00BFC4",
    ">40" = "#C77CFF"
  )
)

select2 <- rownames(vsd3) %in% c(FL_dn$name, FL_up$name)

pheatmap(
  assay(vsd3)[select2, ],
  cluster_rows = TRUE,
  labels_row = Ensm_to_sym$hgnc_symbol[
    Ensm_to_sym$ensblID %in%
      rownames(assay(vsd3)[select2, ])
  ],
  cluster_cols = TRUE,
  annotation_col = Categories_to_show,
  annotation_colors = Categories_color,
  show_colnames = FALSE,
  cellwidth = 2,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  cutree_cols = 3,
  filename = paste0(images_dir, "/heatmap_A_fetal_sig_ALL_AF4.pdf")
)

# pheatmap(assay(vsd3)[select2, ],
#   cluster_rows = TRUE,
#   show_rownames = FALSE,
#   cluster_cols = TRUE,
#   annotation_col = Catagories_to_show,
#   show_colnames = TRUE,
#   filename = paste0(images_dir, "/heatmap_A_from_script_4.1.pdf"),
# )


# PCA with our FL/BM genes
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib")
select <- rownames(vsd3) %in% ids$ensembl_gene_id # find the our gene sets
cdata <- colData(vsd3)
select4 <- cdata$type == "KMT2A-AFF1"
pca2 <- prcomp(t(assay(vsd3)[select, select4])) # calculate PCA
percentVar <- pca2$sdev^2 / sum(pca2$sdev^2) # calc PCA contribution
intgroup.df <- as.data.frame(colData(vsd3)[select4, intgroup, drop = FALSE]) # add grouping if wanted
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = ":"))
} else {
  colData(vsd3)[[intgroup]]
}
princpial3 <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2], PC3 = pca2$x[, 3], group = group, intgroup.df, name = colnames(vsd3)[select4])

Ensm_to_sym <- data.frame(
  ensblID = rownames(pca2$rotation),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(pca2$rotation), mart = ensembl
  )
)

FLor_BM <- NULL
for (i in 1:length(Ensm_to_sym$hgnc_symbol)) {
  x <- GOI$FLorBM[GOI$gene_name %in% as.character(lapply(Ensm_to_sym$hgnc_symbol, toupper))[i]]
  FLor_BM[i] <- x
}

loading <- data.frame(
  LO1 = pca2$rotation[, 1], LO2 = pca2$rotation[, 2], LO3 = pca2$rotation[, 3],
  name = Ensm_to_sym$ensblID, symb = Ensm_to_sym$hgnc_symbol,
  origin = FLor_BM
)

ggplot(data = princpial3, aes_string(x = "PC1", y = "PC2", color = "age_group")) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(paste0(images_dir, "/PC12_AF4_FL_core.pdf"), plot = last_plot())


ggplot(data = princpial3, aes_string(x = "PC1", y = "age", color = "age_group")) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("Age")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(paste0(images_dir, "/PC1Age_AF4_FL_core.pdf"), plot = last_plot())

ngenes <- 25
FL_up <- head(loading[order(loading$LO1, decreasing = TRUE), ], ngenes) # the genes contributing most to fetal signature
sum(str_count(FL_up$origin, pattern = "BM")) / ngenes * 100
FL_dn <- head(loading[order(loading$LO1, decreasing = FALSE), ], ngenes) # the genes contributing most to adult signature
sum(str_count(FL_dn$origin, pattern = "FL")) / ngenes * 100
write.csv(x = FL_dn, file = paste0(out_dir, "/results/down_inf_luek_AF4.csv"))
write.csv(x = FL_up, file = paste0(out_dir, "/results/up_inf_luek_AF4.csv"))

select2 <- rownames(vsd3) %in% c(FL_dn$name, FL_up$name)

Categories_to_show_AF4 <- as.data.frame(colData(vsd3)[, c("age_group", "gender")])
Categories_color_AF4 <- list(
  gender = c("Male" = "#F8766D", "Female" = "#7CAE00"),
  age_group = c(
    "0-2" = "#F8766D", "2-16" = "#7CAE00",
    "16-40" = "#00BFC4", ">40" = "#C77CFF"
  )
)


pheatmap(assay(vsd3)[select2, select4],
  cluster_rows = TRUE,
  labels_row = Ensm_to_sym$hgnc_symbol[Ensm_to_sym$ensblID %in% rownames(assay(vsd3)[select2, ])],
  cluster_cols = TRUE,
  annotation_col = Categories_to_show_AF4,
  annotation_colors = Categories_color_AF4,
  show_colnames = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  cutree_cols = 5,
  filename = paste0(images_dir, "/heatmap_B_fetal_sig_ALL_AF4.pdf")
)

# HOX analysis
hox_OI <- c(
  "HOXA-AS3", "HOXA1", "HOXA10", "HOXA10-AS", "HOXA11-AS", "HOXA2",
  "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB-AS1",
  "HOXB-AS3", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7",
  "HOXB8", "HOXB9"
)

HOX_ids <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = hox_OI,
  mart = ensembl
)

select_hox <- rownames(vsd3) %in% HOX_ids$ensembl_gene_id

Ensm_to_sym_hox <- data.frame(
  ensblID = rownames(assay(vsd3)[select_hox, ]),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(assay(vsd3)[select_hox, ]),
    mart = ensembl
  )
)


pheatmap(assay(vsd3)[select_hox, select4],
  cluster_rows = TRUE,
  labels_row = Ensm_to_sym_hox$hgnc_symbol,
  cluster_cols = TRUE,
  annotation_col = Categories_to_show_AF4,
  annotation_colors = Categories_color_AF4,
  show_colnames = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  filename = paste0(images_dir, "/heatmap_Hox_genes_AF4.pdf")
)
