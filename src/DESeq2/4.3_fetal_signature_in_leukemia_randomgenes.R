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

sample_table3 <- read.table(str_c(
  external_dir,
  "/BALL-1988S-HTSeq/B-ALL-subtyping.txt"
  ),
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
)
sample_table3 <- sample_table3[sample_table3$primary.subtype %in% subtype_oi, ]

sample_table3 <- subset.data.frame(sample_table3, !age == "NA")
sample_table3 <- subset.data.frame(sample_table3, !gender == "NA")

for (i in 1:dim(sample_table3)[1]) {
  if (sample_table3$fusion[i] == 'NoFusion') {
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

counts3 <- DESeqDataSetFromHTSeqCount(sampleTable = DEsampleTable3, design = ~RNA_lib)
counts3 <- DESeq(counts3)

vsd3 <- vst(counts3, blind = FALSE)
plotPCA(vsd3, "RNA_lib") + theme_bw() + coord_fixed(ratio = 8 / 6.3)

assay(vsd3) <- limma::removeBatchEffect(assay(vsd3), vsd3$RNA_lib)
plotPCA(vsd3, "RNA_lib") + theme_bw() + coord_fixed(ratio = 5.6 / 7.2)

ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
  path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl"
)

intgroup <- c("age_group", "type", "age", "gender", "RNA_lib")

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

# ran 5 randoms for paper, only one as representative now
for (i in c(1:5)) {
# for (i in c(1)) {
  print(i)

  # the size is the number of genes matching the fetal core genes
  # present in the DESeq object in script 4.1_fetal_signature_in_leukemia.R
  random_num <- sample.int(length(rownames(assay(vsd3))),
    size = 691, replace = FALSE
  )
  select <- logical(length = length(rownames(assay(vsd3))))
  select[random_num] <- TRUE
  
  # make a combined PCA for all samples
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
  
  gene_ids <-getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = rownames(pca$rotation), 
                   mart = ensembl)
  
  write.csv(gene_ids, file = paste0(images_dir, '/random_', i, '.csv'))
  
  
  # Fig 5b, right (representative figure shown in manuscript)
  f <- ggplot(data = princpial3,
              aes_string(x = "PC1", y = "age", color = "type")) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    scale_color_manual(values = type_colours) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    coord_fixed()
  print(f)
  ggsave(paste0(images_dir, "/random_samePC1_all_", i, ".pdf"), plot = last_plot())
  
  # plot each subtype individually in the same PCA space
  for (y in 1:length(unique(princpial3$type))) {
    print(y)
    f <- ggplot(
      data = subset.data.frame(princpial3, type == unique(princpial3$type)[y]),
      aes_string(x = "PC1", y = "age", color = "type")
    ) +
      geom_point(size = 3.5) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("Age")) +
      theme_bw() +
      theme(aspect.ratio = 1.5) +
      xlim(-38, 38) +
      ylim(-1, 82) +
      scale_color_manual(values = type_colours) +
      coord_fixed()
    print(f)
    ggsave(paste0(
      images_dir, "/", unique(princpial3$type)[y],
      "_random_samePC1_", unique(princpial3$type)[y], "_", i,
      ".pdf"
    ), plot = last_plot())
  }
  
  #subset vsd3 object and perform PCA individually using the random genes
  for (z in 1:length(unique(vsd3$type))) {
    print(z)
    print(unique(vsd3$type)[z])
    vsd3.sub <- vsd3[, vsd3$type %in% unique(vsd3$type)[z]]
    plotPCA(vsd3.sub, "type")
    
    # breaking the PCA analysis into its parts
    pca <- prcomp(t(assay(vsd3.sub)[select, ])) # calculate PCA
    percentVar <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
    intgroup.df <- as.data.frame(colData(vsd3.sub)[, intgroup, drop = FALSE]) # add grouping if wanted
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
      colData(vsd3.sub)[[intgroup]]
    }
    princpial3 <- data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[
        ,
        3
      ], group = group, intgroup.df, name =
        colnames(vsd3.sub)
    )
    
    # KMT2A-AFF1 in Fig 5d, right, rest in fig S6c, lower row (representative
    # figures shown in manuscript)
    p <- ggplot(data = princpial3, aes_string(x = "PC1", y = "age", color = "age_group")) +
      geom_point(size = 3.5) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("Age")) +
      theme_bw() +
      theme(aspect.ratio = 1.5) +
      scale_color_manual(values = age_colours) +
      xlim(-38, 38) +
      ylim(-1, 82) +
      coord_fixed()
    print(p)
    ggsave(paste0(
      images_dir, "/",
      "random_PC1_", unique(vsd3$type)[z], "_", i, ".pdf"
    ),
    plot = last_plot()
    )
  }
}

# random_PC1_KMT2A-MLLT1_1,random_PC1_BCR-ABL1_1,random_PC1_KMT2A-AFF1_1,random_PC1_High hyperdiploid_1,random_PC1_Low hyperdiploid_1,random_PC1_ETV6-RUNX1_1,random_PC1_KMT2A-MLLT3_1





