invisible(lapply(list(
  "stringr",
  "biomaRt",
  "ggplot2",
  "limma",
  "pheatmap",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

leukemia_types <- snakemake@config[["leukemia_types"]]
print(leukemia_types)

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

make_named_version_nested <- function(ids_level1,
                                      ids_level2,
                                      file_list,
                                      pattern_level1 = "\\.",
                                      pattern_level2 = "_") {
  named_file_list <- unlist(
    lapply(ids_level1,
      FUN = function(id) {
        file_list[str_detect(
          string = file_list,
          pattern = str_c(id, pattern_level1)
        )]
      }
    )
  )
  # print(named_file_list)
  new_list <- split(
    named_file_list,
    ceiling(seq_along(named_file_list) / length(ids_level2))
  )
  # print(new_list)
  nested_list <- lapply(new_list, FUN = function(sub_list) {
    make_named_version(
      ids = ids_level2, file_list = sub_list, pattern = pattern_level2
    )
  })
  # print(nested_list)
  # names(named_file_list) <- ids_level1
  return(nested_list)
}

plots_random_pc1_all <- snakemake@output[["plots_random_pc1_all"]]
named_plots_random_pc1_all <- make_named_version(
  ids = as.character(1:5),
  file_list = plots_random_pc1_all,
  pattern = "\\."
)
plots_random_pc1_type <- snakemake@output[["plots_random_pc1_type"]]
named_plots_random_pc1_type <- make_named_version_nested(
  ids_level1 = as.character(1:5),
  ids_level2 = leukemia_types,
  file_list = plots_random_pc1_type,
  pattern_level1 = "\\.",
  pattern_level2 = "_"
)
plots_random_pc1_type_type <- snakemake@output[["plots_random_pc1_type_type"]]
named_plots_random_pc1_type_type <- make_named_version_nested(
  ids_level1 = as.character(1:5),
  ids_level2 = leukemia_types,
  file_list = plots_random_pc1_type_type,
  pattern_level1 = "\\.",
  pattern_level2 = "_"
)
# print("> plots_random_pc1_all:")
# print(named_plots_random_pc1_all)
# print("> plots_random_pc1_type:")
# print(named_plots_random_pc1_type)
# print("> plots_random_pc1_type_type:")
# print(named_plots_random_pc1_type_type)

random_csvs <- snakemake@output[["random_genes"]]
named_random_csvs <- make_named_version(
  ids = as.character(1:5),
  file_list = random_csvs,
  pattern = "\\."
)
# print("> random_csvs:")
# print(named_random_csvs)

sample_table <- read.table(snakemake@input[["ball_tsv"]],
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
# plotPCA(vsd, "RNA_lib") + theme_bw() + coord_fixed(ratio = 8 / 6.3)

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)
# plotPCA(vsd, "RNA_lib") + theme_bw() + coord_fixed(ratio = 5.6 / 7.2)

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

n_genes <- sum(rownames(vsd) %in% ids$ensembl_gene_id)
# select number of genes matching the number of fetal cores genes present in the
# object

# subset assay object to not contain genes with lower expression than the
# lowest-expressed gene in the fetal core
norm_counts_use <- assay(vsd)[rowMeans(assay(vsd)) >=
  min(rowMeans(assay(vsd)[rownames(assay(vsd))
  %in% ids$ensembl_gene_id, ])), ]

# ran 5 randoms for paper, only one as representative now
for (i in c(1:5)) {
# for (i in c(1)) {
  print(i)

  genes_use <- norm_counts_use[sample(nrow(norm_counts_use), n_genes,
    replace = FALSE
  ), ]

  # make a combined PCA for all samples
  pca <- prcomp(t(assay(vsd)[
    rownames(assay(vsd)) %in% rownames(genes_use),
  ]))
  # calculate PCA

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

  gene_ids <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(pca$rotation),
    mart = ensembl
  )

  write.csv(gene_ids, file = snakemake@output[["random_genes"]][[i]])
  # write.csv(gene_ids, file = paste0(images_dir, "/random_", i, ".csv"))

  # Fig 5b, right (representative figure shown in manuscript)
  f <- ggplot(data = principal,
              aes_string(x = "PC1", y = "age", color = "type")) +
    geom_point(size = 3.5) +
    xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
    ylab(paste0("Age")) +
    theme_bw() +
    theme(aspect.ratio = 1.5) +
    scale_color_manual(values = type_colours) +
    xlim(-38, 38) +
    ylim(-1, 82) +
    coord_fixed()
  # print(f)
  ggsave(named_plots_random_pc1_all[[i]], plot = f)
  # ggsave(paste0(images_dir, "/random_samePC1_all_", i, ".pdf"), plot = last_plot())

  # plot each subtype individually in the same PCA space
  for (y in 1:length(unique(principal$type))) {
    print(y)
    f <- ggplot(
      data = subset.data.frame(principal, type == unique(principal$type)[y]),
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
    # print(f)
    ggsave(named_plots_random_pc1_type_type[[i]][[unique(vsd$type)[y]]], plot = f)
    # ggsave(paste0(
    #   images_dir, "/", unique(princpial3$type)[y],
    #   "_random_samePC1_", unique(princpial3$type)[y], "_", i,
    #   ".pdf"
    # ), plot = last_plot())
  }

  #subset vsd object and perform PCA individually using the random genes
  for (z in 1:length(unique(vsd$type))) {
    print(z)
    print(unique(vsd$type)[z])
    vsd_sub <- vsd[, vsd$type %in% unique(vsd$type)[z]]
    # plotPCA(vsd_sub, "type")

    # breaking the PCA analysis into its parts:
    # calculate PCA
    pca <- prcomp(t(assay(vsd_sub)[
      rownames(assay(vsd_sub)) %in% rownames(genes_use),
    ]))
    percent_var <- pca$sdev^2 / sum(pca$sdev^2) # calc PCA contribution
    intgroup_df <- as.data.frame(colData(vsd_sub)[, intgroup, drop = FALSE])
    # add grouping if wanted
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup_df, 1, paste, collapse = ":"))
    } else {
      colData(vsd_sub)[[intgroup]]
    }
    principal <- data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
      group = group, intgroup_df, name =
        colnames(vsd_sub)
    )

    # KMT2A-AFF1 in Fig 5d, right, rest in fig S6c, lower row (representative
    # figures shown in manuscript)
    p <- ggplot(
      data = principal,
      aes_string(
        x = "PC1", y = "age", color = "age_group"
      )
    ) +
      geom_point(size = 3.5) +
      xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
      ylab(paste0("Age")) +
      theme_bw() +
      theme(aspect.ratio = 1.5) +
      scale_color_manual(values = age_colours) +
      xlim(-38, 38) +
      ylim(-1, 82) +
      coord_fixed()
    # print(p)
    # ggsave(paste0(
    #   images_dir, "/",
    #   "random_PC1_", unique(vsd$type)[z], "_", i, ".pdf"
    # ),
    ggsave(named_plots_random_pc1_type[[i]][[unique(vsd$type)[z]]],
      plot = p
    )
  }
}

# random_PC1_KMT2A-MLLT1_1,random_PC1_BCR-ABL1_1,random_PC1_KMT2A-AFF1_1,
# random_PC1_High
# hyperdiploid_1,random_PC1_Low
# hyperdiploid_1,random_PC1_ETV6-RUNX1_1,random_PC1_KMT2A-MLLT3_1
