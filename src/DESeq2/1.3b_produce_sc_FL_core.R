invisible(lapply(list(
  "stringr",
  "dplyr",
  "DESeq2"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])
clusters <- snakemake@config[["fl_clusters_to_use"]]
deg_results_bm_files <- snakemake@input[["deg_results_bm"]]
deg_results_fl_files <- snakemake@input[["deg_results_fl"]]
sub_setter_pos_files <- snakemake@output[["sub_setter_pos"]]
sub_setter_neg_files <- snakemake@output[["sub_setter_neg"]]

fc_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["log2foldchange"]][["main"]]
padj_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["adjustedpvalue"]]

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


named_deg_results_bm_files <- make_named_version(clusters, deg_results_bm_files)
named_deg_results_fl_files <- make_named_version(clusters, deg_results_fl_files)
stopifnot(length(named_deg_results_bm_files) == length(clusters))
stopifnot(length(named_deg_results_fl_files) == length(clusters))

named_sub_setter_pos_files <- make_named_version(clusters, sub_setter_pos_files)
named_sub_setter_neg_files <- make_named_version(clusters, sub_setter_neg_files)
stopifnot(length(named_sub_setter_pos_files) == length(clusters))
stopifnot(length(named_sub_setter_neg_files) == length(clusters))

deg_results_bm <- lapply(clusters, FUN = function(cluster) {
  read.table(
    file = named_deg_results_bm_files[[cluster]],
  )
})
names(deg_results_bm) <- clusters
deg_results_fl <- lapply(clusters, FUN = function(cluster) {
  read.table(
    file = named_deg_results_fl_files[[cluster]],
  )
})
names(deg_results_fl) <- clusters

sub_setter <- function(df) {
  df_sub <- df[!is.na(df$p_val_adj), ]
  df_sub <- df_sub[
    df_sub$avg_log2FC > fc_threshold & df_sub$p_val_adj < padj_threshold,
  ]
  return(df_sub)
}

write_sub_setter <- function(cluster) {
  write.table(sub_setter(deg_results_fl[[cluster]]),
    file = named_sub_setter_pos_files[[cluster]],
    sep = ","
  )
  write.table(sub_setter(deg_results_bm[[cluster]]),
    file = named_sub_setter_neg_files[[cluster]],
    sep = ","
  )
}

for (cluster in clusters) {
  print(cluster)
  write_sub_setter(cluster)
}

prim_pos <- intersect(
  intersect(
    intersect(
      rownames(sub_setter(deg_results_fl[["HSC"]])),
      rownames(sub_setter(deg_results_fl[["MPP.I"]]))
    ),
    rownames(sub_setter(deg_results_fl[["MPP.II"]]))
  ), rownames(sub_setter(deg_results_fl[["Cyc"]]))
)

lymph_pos <- intersect(
  rownames(sub_setter(deg_results_fl[["Ly.I"]])),
  rownames(sub_setter(deg_results_fl[["Ly.II"]]))
)

mye_pos <- intersect(
  intersect(
    rownames(sub_setter(deg_results_fl[["DC.Mono"]])),
    rownames(sub_setter(deg_results_fl[["GMP"]]))
  ),
  rownames(sub_setter(deg_results_fl[["MEP"]]))
)

all_core_pos <- unique(c(prim_pos, lymph_pos, mye_pos))
fl_core_pos <- intersect(intersect(prim_pos, lymph_pos), mye_pos)
print(lapply(list(prim_pos, lymph_pos, mye_pos), length))

prim_neg <- intersect(
  intersect(
    intersect(
      rownames(sub_setter(deg_results_bm[["HSC"]])),
      rownames(sub_setter(deg_results_bm[["MPP.I"]]))
    ),
    rownames(sub_setter(deg_results_bm[["MPP.II"]]))
  ),
  rownames(sub_setter(deg_results_bm[["Cyc"]]))
)

lymph_neg <- intersect(
  rownames(sub_setter(deg_results_bm[["Ly.I"]])),
  rownames(sub_setter(deg_results_bm[["Ly.II"]]))
)
mye_neg <- intersect(
  intersect(
    rownames(sub_setter(deg_results_bm[["DC.Mono"]])),
    rownames(sub_setter(deg_results_bm[["GMP"]]))
  ),
  rownames(sub_setter(deg_results_bm[["MEP"]]))
)

all_core_neg <- unique(c(prim_neg, lymph_neg, mye_neg))
fl_core_neg <- intersect(intersect(prim_neg, lymph_neg), mye_neg)
print(lapply(list(prim_neg, lymph_neg, mye_neg), length))

core <- list(
  prim_pos, lymph_pos, mye_pos, fl_core_pos, all_core_pos,
  prim_neg, lymph_neg, mye_neg, fl_core_neg, all_core_neg
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
  file = snakemake@output[["de_genes"]],
  quote = FALSE, sep = ","
)

write.table(all_core_neg,
  file = snakemake@output[["adult_signature"]],
  row.names = F, col.names = F, quote = F
)

write.table(all_core_pos,
  file = snakemake@output[["fetal_signature"]],
  row.names = F, col.names = F, quote = F
)


######################################################
# Do this for the extra cutoffs as well
######################################################

get_cutoff_files <- function(struct, cutoff) {
  struct[
    str_detect(
      pattern = cutoff %>% as.character(),
      struct
    )
  ]
}

for (cutoff in
  snakemake@config[["seurat_deg_cutoffs"]][["log2foldchange"]][["extra"]]) {
  print(cutoff)

  deg_results_bm_files <- get_cutoff_files(
    struct = snakemake@input[["deg_results_bm_extra"]],
    cutoff = cutoff
  )
  deg_results_fl_files <- get_cutoff_files(
    struct = snakemake@input[["deg_results_fl_extra"]],
    cutoff = cutoff
  )
  sub_setter_pos_files <- get_cutoff_files(
    struct = snakemake@output[["sub_setter_pos_extra"]],
    cutoff = cutoff
  )
  sub_setter_neg_files <- get_cutoff_files(
    struct = snakemake@output[["sub_setter_neg_extra"]],
    cutoff = cutoff
  )

  # print(sub_setter_pos_files)
  # print(sub_setter_neg_files)

  # print(is.null(deg_results_bm_files))
  # print(is.null(deg_results_fl_files))
  # print(is.null(sub_setter_pos_files))
  # print(is.null(sub_setter_neg_files))

  stopifnot(
    all(
      c(
        !is.null(deg_results_bm_files),
        !is.null(deg_results_fl_files),
        !is.null(sub_setter_pos_files),
        !is.null(sub_setter_neg_files)
      )
    )
  )

  fc_threshold <- cutoff
  padj_threshold <- snakemake@config[["seurat_deg_cutoffs"]][["adjustedpvalue"]]

  make_named_version <- function(ids, file_list, pattern = "_") {
    named_file_list <- unlist(
      lapply(ids,
        FUN = function(id) {
          file_list[str_detect(string = file_list, pattern = str_c(
            id,
            pattern
          ))]
        }
      )
    )
    names(named_file_list) <- ids
    return(named_file_list)
  }

  named_deg_results_bm_files <- make_named_version(
    clusters,
    deg_results_bm_files
  )
  named_deg_results_fl_files <- make_named_version(
    clusters,
    deg_results_fl_files
  )
  stopifnot(length(named_deg_results_bm_files) == length(clusters))
  stopifnot(length(named_deg_results_fl_files) == length(clusters))

  named_sub_setter_pos_files <- make_named_version(
    clusters,
    sub_setter_pos_files
  )
  named_sub_setter_neg_files <- make_named_version(
    clusters,
    sub_setter_neg_files
  )
  stopifnot(length(named_sub_setter_pos_files) == length(clusters))
  stopifnot(length(named_sub_setter_neg_files) == length(clusters))

  deg_results_bm <- lapply(clusters, FUN = function(cluster) {
    read.table(
      file = named_deg_results_bm_files[[cluster]],
    )
  })
  names(deg_results_bm) <- clusters
  deg_results_fl <- lapply(clusters, FUN = function(cluster) {
    read.table(
      file = named_deg_results_fl_files[[cluster]],
    )
  })
  names(deg_results_fl) <- clusters

  sub_setter <- function(df) {
    df_sub <- df[!is.na(df$p_val_adj), ]
    df_sub <- df_sub[
      df_sub$avg_log2FC > fc_threshold & df_sub$p_val_adj < padj_threshold,
    ]
    return(df_sub)
  }

  write_sub_setter <- function(cluster) {
    write.table(sub_setter(deg_results_fl[[cluster]]),
      file = named_sub_setter_pos_files[[cluster]],
      sep = ","
    )
    write.table(sub_setter(deg_results_bm[[cluster]]),
      file = named_sub_setter_neg_files[[cluster]],
      sep = ","
    )
  }

  for (cluster in clusters) {
    print(cluster)
    write_sub_setter(cluster)
  }

  prim_pos <- intersect(
    intersect(
      intersect(
        rownames(sub_setter(deg_results_fl[["HSC"]])),
        rownames(sub_setter(deg_results_fl[["MPP.I"]]))
      ),
      rownames(sub_setter(deg_results_fl[["MPP.II"]]))
    ), rownames(sub_setter(deg_results_fl[["Cyc"]]))
  )

  lymph_pos <- intersect(
    rownames(sub_setter(deg_results_fl[["Ly.I"]])),
    rownames(sub_setter(deg_results_fl[["Ly.II"]]))
  )

  mye_pos <- intersect(
    intersect(
      rownames(sub_setter(deg_results_fl[["DC.Mono"]])),
      rownames(sub_setter(deg_results_fl[["GMP"]]))
    ),
    rownames(sub_setter(deg_results_fl[["MEP"]]))
  )

  all_core_pos <- unique(c(prim_pos, lymph_pos, mye_pos))
  fl_core_pos <- intersect(intersect(prim_pos, lymph_pos), mye_pos)
  print(lapply(list(prim_pos, lymph_pos, mye_pos), length))

  prim_neg <- intersect(
    intersect(
      intersect(
        rownames(sub_setter(deg_results_bm[["HSC"]])),
        rownames(sub_setter(deg_results_bm[["MPP.I"]]))
      ),
      rownames(sub_setter(deg_results_bm[["MPP.II"]]))
    ),
    rownames(sub_setter(deg_results_bm[["Cyc"]]))
  )

  lymph_neg <- intersect(
    rownames(sub_setter(deg_results_bm[["Ly.I"]])),
    rownames(sub_setter(deg_results_bm[["Ly.II"]]))
  )
  mye_neg <- intersect(
    intersect(
      rownames(sub_setter(deg_results_bm[["DC.Mono"]])),
      rownames(sub_setter(deg_results_bm[["GMP"]]))
    ),
    rownames(sub_setter(deg_results_bm[["MEP"]]))
  )

  all_core_neg <- unique(c(prim_neg, lymph_neg, mye_neg))
  fl_core_neg <- intersect(intersect(prim_neg, lymph_neg), mye_neg)
  print(lapply(list(prim_neg, lymph_neg, mye_neg), length))

  core <- list(
    prim_pos, lymph_pos, mye_pos, fl_core_pos, all_core_pos,
    prim_neg, lymph_neg, mye_neg, fl_core_neg, all_core_neg
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

  filename <- get_cutoff_files(
    struct = snakemake@output[["de_genes_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  write.table(
    x = core_filled,
    file = filename,
    quote = FALSE, sep = ","
  )

  filename <- get_cutoff_files(
    struct = snakemake@output[["adult_signature_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  write.table(all_core_neg,
    file = filename,
    row.names = F, col.names = F, quote = F
  )

  filename <- get_cutoff_files(
    struct = snakemake@output[["fetal_signature_extra"]],
    cutoff = cutoff
  )
  stopifnot(!is.null(filename))
  write.table(all_core_pos,
    file = filename,
    row.names = F, col.names = F, quote = F
  )
}
