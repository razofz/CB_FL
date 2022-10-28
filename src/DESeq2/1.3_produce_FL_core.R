# comparison DEseq2 res Pseudo bulking approaches
source("env.R")
in_dir <- "data/processed/DEseq2/results/"
out_dir <- "data/processed/DEseq2/results/"
sup_table_dir <- "data/raw/paper_specific/Supplemental_tables/"
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("stringr"))

HSC <- read.csv(
  file = paste0(in_dir, "HSC_FL_yBM.csv"), header = TRUE,
  sep = ","
)
Cyc <- read.csv(
  file = paste0(in_dir, "Cyc_FL_yBM.csv"), header = TRUE,
  sep = ","
)
DCmono <- read.csv(
  file = paste0(in_dir, "DC.Mono_FL_yBM.csv"), header =
    TRUE, sep = ","
)
GMP <- read.csv(
  file = paste0(in_dir, "GMP_FL_yBM.csv"), header = TRUE,
  sep = ","
)
LyI <- read.csv(
  file = paste0(in_dir, "Ly.I_FL_yBM.csv"), header = TRUE,
  sep = ","
)
LyII <- read.csv(
  file = paste0(in_dir, "Ly.II_FL_yBM.csv"), header = TRUE,
  sep = ","
)
MEP <- read.csv(
  file = paste0(in_dir, "MEP_FL_yBM.csv"), header = TRUE,
  sep = ","
)
MPPI <- read.csv(
  file = paste0(in_dir, "MPP.I_FL_yBM.csv"), header = TRUE,
  sep = ","
)
MPPII <- read.csv(
  file = paste0(in_dir, "MPP.II_FL_yBM.csv"), header =
    TRUE, sep = ","
)

Sub_setter_pos <- function(DF) {
  DF_sub <- DF[!is.na(DF$padj),]
  DF_sub <- DF_sub[DF_sub$log2FoldChange > 1 & DF_sub$padj < 0.05, ]
}

Sub_setter_neg <- function(DF) {
  DF_sub <- DF[!is.na(DF$padj),]
  DF_sub <- DF_sub[DF_sub$log2FoldChange < -1 & DF_sub$padj < 0.05, ]
}

clusters <- c(
  "HSC",
  "Cyc",
  "DC.Mono",
  "GMP",
  "Ly.I",
  "Ly.II",
  "MEP",
  "MPP.I",
  "MPP.II"
)

sub_setter_dir <- paste0(out_dir, "sub_setter/")
if (dir.exists(sub_setter_dir) == FALSE) {
  dir.create(sub_setter_dir,
    showWarnings = TRUE,
    recursive = TRUE
  )
}

write_sub_setter <- function(cluster, cluster_df) {
  write.table(Sub_setter_pos(cluster_df),
    file = paste0(sub_setter_dir, cluster, "_pos.csv"),
    sep = ","
  )
  write.table(Sub_setter_neg(cluster_df),
    file = paste0(sub_setter_dir, cluster, "_neg.csv"),
    sep = ","
  )
}

write_sub_setter("HSC", HSC)
write_sub_setter("Cyc", Cyc)
write_sub_setter("DC.Mono", DCmono)
write_sub_setter("GMP", GMP)
write_sub_setter("Ly.I", LyI)
write_sub_setter("Ly.II", LyII)
write_sub_setter("MEP", MEP)
write_sub_setter("MPP.I", MPPI)
write_sub_setter("MPP.II", MPPII)

prim_pos <- intersect(intersect(intersect(
    rownames(Sub_setter_pos(HSC)), rownames(Sub_setter_pos(MPPI))),
    rownames(Sub_setter_pos(MPPII))), rownames(Sub_setter_pos(Cyc)))

lymph_pos <- intersect(
  rownames(Sub_setter_pos(LyI)),
  rownames(Sub_setter_pos(LyII))
)

mye_pos <- intersect(
  intersect(
    rownames(Sub_setter_pos(DCmono)),
    rownames(Sub_setter_pos(GMP))
  ),
  rownames(Sub_setter_pos(MEP))
)

all_core_pos <- unique(c(prim_pos, lymph_pos, mye_pos))
FL_core_pos <- intersect(intersect(prim_pos, lymph_pos), mye_pos)

prim_neg <- intersect(
  intersect(
    intersect(
      rownames(Sub_setter_neg(HSC)),
      rownames(Sub_setter_neg(MPPI))
    ),
    rownames(Sub_setter_neg(MPPII))
  ),
  rownames(Sub_setter_neg(Cyc))
)

lymph_neg <- intersect(
  rownames(Sub_setter_neg(LyI)),
  rownames(Sub_setter_neg(LyII))
)
mye_neg <- intersect(
  intersect(
    rownames(Sub_setter_neg(DCmono)),
    rownames(Sub_setter_neg(GMP))
  ),
  rownames(Sub_setter_neg(MEP))
)

all_core_neg <- unique(c(prim_neg, lymph_neg, mye_neg))
FL_core_neg <- intersect(intersect(prim_neg, lymph_neg), mye_neg)

core <- list(
  prim_pos, lymph_pos, mye_pos, FL_core_pos, all_core_pos, prim_neg,
  lymph_neg, mye_neg, FL_core_neg, all_core_neg
)

max_length <- max(sapply(core, length))
core_filled <- sapply(core, function(x) {
  c(x, rep(NA, max_length - length(x)))
})

core_filled <- data.frame(core_filled)
colnames(core_filled) <- c(
  "prim_pos", "lymph_pos", "mye_pos", "FL_core_pos",
  "all_core_pos", "prim_neg", "lymph_neg", "mye_neg",
  "FL_core_neg", "all_core_neg"
)

write.table(
  x = core_filled,
  file = paste0(out_dir, "pseudorep_DE_genes.csv"),
  quote = FALSE, sep = ","
)

write.table(all_core_neg,
  file = paste0(
    # sup_table_dir,
    out_dir,
    "Adult_signature_p005.txt"
  ), row.names = F,
  col.names = F, quote = F
)

write.table(all_core_pos,
  file = paste0(
    # sup_table_dir,
    out_dir,
    "Fetal_signature_p005.txt"
  ), row.names = F,
  col.names = F, quote = F
)
