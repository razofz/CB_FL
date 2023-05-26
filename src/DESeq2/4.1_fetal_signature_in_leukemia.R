# loading KMT2A-MLLT1, KMT2A-MLLT3, ETV6-RUNX1, BCR-ABL1
invisible(lapply(list(
  "stringr",
  "biomaRt",
  "ggplot2",
  "pheatmap",
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
  snakemake@input[["ball_tsv"]],
  header = TRUE, dec = ",", sep = "\t", quote = '""'
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
) # , 'High hyperdiploid', 'Low hyperdiploid'
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
plotPCA(vsd, "RNA_lib") + theme_bw()
ggsave(
  snakemake@output[["plot_pc12_prebatch"]],
  # paste0(images_dir, "/PC12_top500_preBatch_ALL_leuk.pdf"),
  plot = last_plot()
)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)
plotPCA(vsd, "RNA_lib") + theme_bw()
ggsave(
  snakemake@output[["plot_pc12_postbatch"]],
  # paste0(images_dir, "/PC12_top500_postBatch_ALL_leuk.pdf"),
  plot = last_plot()
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

# PCA with our FL/BM genes
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib") #
select <- rownames(vsd) %in% ids$ensembl_gene_id # find the our gene sets
length(select)
table(select)
pca <- prcomp(t(assay(vsd)[select, ])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
intgroup_df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
# add grouping if wanted
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
  x <- goi$fl_or_bm[
    goi$gene_name %in% as.character(lapply(ensm_to_sym$hgnc_symbol, toupper))[i]
  ]
  fl_or_bm[i] <- x
}

loading <- data.frame(
  LO1 = pca$rotation[, 1], LO2 = pca$rotation[, 2], LO3 = pca$rotation[, 3],
  name = ensm_to_sym$ensblID, symb = ensm_to_sym$hgnc_symbol,
  origin = fl_or_bm
)

ggplot(data = principal, aes_string(x = "PC1", y = "age", color = "type")) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  xlim(-35, 35) +
  ylim(-1, 82) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(snakemake@output[["plot_pc1_age"]], plot = last_plot())
# ggsave(paste0(images_dir, "/PC1Age_All_leuk_FL_core.pdf"), plot = last_plot())

ggplot(data = principal, aes_string(x = "PC2", y = "age", color = "type")) +
  geom_point(size = 3) +
  xlab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(snakemake@output[["plot_pc2_age"]], plot = last_plot())
# ggsave(paste0(images_dir, "/PC2Age_All_leuk_FL_core.pdf"), plot = last_plot())

images_dir <- dirname(snakemake@output[["plot_pc1_age"]])

for (x in 1:length(unique(principal$type))) {
  print(x)
  f <- ggplot(
    data = subset.data.frame(principal, type == unique(principal$type)[x]),
    aes_string(x = "PC1", y = "age", color = "type")
  ) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    xlim(-35, 35) +
    ylim(-1, 82) +
    coord_fixed()
  print(f)
  ggsave(
    paste0(paste0(
      images_dir, "/PC1Age_", unique(principal$type)[x],
      "_FL_core.pdf"
    )),
    plot = last_plot()
  )
}

ngenes <- 30
# the genes contributing most to fetal signature
fl_up <- head(loading[order(loading$LO1, decreasing = TRUE), ], ngenes)
sum(str_count(fl_up$origin, pattern = "FL")) / ngenes * 100
# the genes contributing most to adult signature
fl_dn <- head(loading[order(loading$LO1, decreasing = FALSE), ], ngenes)
sum(str_count(fl_dn$origin, pattern = "BM")) / ngenes * 100

write.csv(x = fl_dn, file = snakemake@output[["down_inf_leukemia"]])
write.csv(x = fl_dn, file = snakemake@output[["up_inf_leukemia"]])
# write.csv(x = fl_up, file = paste0(out_dir, "/results/up_inf_luek.csv"))
# write.csv(x = fl_up, file = paste0(out_dir, "/results/up_inf_luek.csv"))

categories_to_show <- as.data.frame(colData(vsd)[, c(
  "age_group", "gender", "type"
)])
categories_color <- list(
  gender = c("Male" = "#2E3192", "Female" = "#BE1E2D"),
  # gender = c("Male" = "#F8766D", "Female" = "#7CAE00"),
  age_group = c(
    "0-2" = "#F8766D", "2-16" = "#7CAE00", "16-40" = "#00BFC4",
    ">40" = "#C77CFF"
  )
)

select <- rownames(vsd) %in% c(fl_dn$name, fl_up$name)

pheatmap(
  assay(vsd)[select, ],
  cluster_rows = TRUE,
  labels_row = ensm_to_sym$hgnc_symbol[
    ensm_to_sym$ensblID %in%
      rownames(assay(vsd)[select, ])
  ],
  cluster_cols = TRUE,
  annotation_col = categories_to_show,
  annotation_colors = categories_color,
  show_colnames = FALSE,
  cellwidth = 2,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  cutree_cols = 3,
  filename = snakemake@output[["plot_heatmap_fl_sig"]]
  # filename = paste0(images_dir, "/heatmap_A_fetal_sig_ALL_AF4.pdf")
)

# PCA with our FL/BM genes
intgroup <- c("age_group", "type", "age", "gender", "RNA_lib")
select <- rownames(vsd) %in% ids$ensembl_gene_id # find the our gene sets
cdata <- colData(vsd)
select_aff4 <- cdata$type == "KMT2A-AFF1"
pca <- prcomp(t(assay(vsd)[select, select_aff4])) # calculate PCA
percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
intgroup_df <- as.data.frame(colData(vsd)[select_aff4, intgroup, drop = FALSE])
# add grouping if wanted
group <- if (length(intgroup) > 1) {
  factor(apply(intgroup_df, 1, paste, collapse = ":"))
} else {
  colData(vsd)[[intgroup]]
}
principal <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  group = group, intgroup_df, name =
    colnames(vsd)[select_aff4]
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

ggplot(
  data = principal,
  aes_string(x = "PC1", y = "PC2", color = "age_group")
) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(snakemake@output[["plot_pc12_af4_fl_core"]], plot = last_plot())
# ggsave(paste0(images_dir, "/PC12_AF4_FL_core.pdf"), plot = last_plot())


ggplot(
  data = principal,
  aes_string(x = "PC1", y = "age", color = "age_group")
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("Age")) +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  theme_bw()
ggsave(snakemake@output[["plot_pc1age_af4_fl_core"]], plot = last_plot())
# ggsave(paste0(images_dir, "/PC1Age_AF4_FL_core.pdf"), plot = last_plot())

ngenes <- 25
fl_up <- head(loading[order(loading$LO1, decreasing = TRUE), ], ngenes)
# the genes contributing most to fetal signature
# sum(str_count(fl_up$origin, pattern = "BM")) / ngenes * 100
fl_dn <- head(loading[order(loading$LO1, decreasing = FALSE), ], ngenes)
# the genes contributing most to adult signature
# sum(str_count(fl_dn$origin, pattern = "FL")) / ngenes * 100
write.csv(x = fl_dn, file = snakemake@output[["down_inf_leukemia_af4"]])
write.csv(x = fl_up, file = snakemake@output[["up_inf_leukemia_af4"]])
# write.csv(x = fl_dn, file = paste0(out_dir, "/results/down_inf_luek_AF4.csv"))
# write.csv(x = fl_up, file = paste0(out_dir, "/results/up_inf_luek_AF4.csv"))

select <- rownames(vsd) %in% c(fl_dn$name, fl_up$name)

categories_to_show_af4 <- as.data.frame(colData(vsd)[, c(
  "age_group",
  "gender"
)])
categories_color_af4 <- list(
  # gender = c("Male" = "#F8766D", "Female" = "#7CAE00"),
  gender = c("Male" = "#2E3192", "Female" = "#BE1E2D"),
  age_group = c(
    "0-2" = "#F8766D", "2-16" = "#7CAE00",
    "16-40" = "#00BFC4", ">40" = "#C77CFF"
  )
)


pheatmap(assay(vsd)[select, select_aff4],
  cluster_rows = TRUE,
  labels_row = ensm_to_sym$hgnc_symbol[ensm_to_sym$ensblID %in% rownames(assay(vsd)[select, ])],
  cluster_cols = TRUE,
  annotation_col = categories_to_show_af4,
  annotation_colors = categories_color_af4,
  show_colnames = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  cutree_cols = 5,
  filename = snakemake@output[["plot_heatmap_fl_sig_af4"]]
  # filename = paste0(images_dir, "/heatmap_B_fetal_sig_ALL_AF4.pdf")
)

# HOX analysis
hox_oi <- c(
  "HOXA-AS3", "HOXA1", "HOXA10", "HOXA10-AS", "HOXA11-AS", "HOXA2",
  "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB-AS1",
  "HOXB-AS3", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7",
  "HOXB8", "HOXB9"
)

hox_ids <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = hox_oi,
  mart = ensembl
)

select_hox <- rownames(vsd) %in% hox_ids$ensembl_gene_id

ensm_to_sym_hox <- data.frame(
  ensblID = rownames(assay(vsd)[select_hox, ]),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(assay(vsd)[select_hox, ]),
    mart = ensembl
  )
)

pheatmap(assay(vsd)[select_hox, select_aff4],
  cluster_rows = TRUE,
  labels_row = ensm_to_sym_hox$hgnc_symbol,
  cluster_cols = TRUE,
  annotation_col = categories_to_show_af4,
  annotation_colors = categories_color_af4,
  show_colnames = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize_row = 8,
  border_color = "black",
  filename = snakemake@output[["plot_heatmap_hox_genes_af4"]]
  # filename = paste0(images_dir, "/heatmap_Hox_genes_AF4.pdf")
)
