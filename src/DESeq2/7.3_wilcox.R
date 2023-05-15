library(ggplot2)

set.seed(snakemake@config[["seed"]])

mod_scores <- read.csv(snakemake@input[["scores"]], row.names = 1)

mod_scores$roy_tissue_type <- factor(mod_scores$roy_tissue_type,
  levels = c("eFL", "FL", "FBM", "PBM", "ABM")
)

w_FL <- pairwise.wilcox.test(
  x = mod_scores$FL1, g = mod_scores$roy_tissue_type,
  p.adjust.method = "bonferroni", alternative = "two.sided"
)
write.csv(w_FL$p.value, file = snakemake@output[["pvals_fl"]])

w_AD <- pairwise.wilcox.test(
  x = mod_scores$AD1, g = mod_scores$roy_tissue_type,
  p.adjust.method = "bonferroni", alternative = "two.sided"
)
write.csv(w_AD$p.value, file = snakemake@output[["pvals_ad"]])
