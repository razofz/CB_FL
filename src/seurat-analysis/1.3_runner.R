curdir <- getwd()
path <- Sys.getenv("PROJECT_PATH")

if (path == "") {
  source("env.R")
  path <- getwd()
  setwd("src/seurat-analysis")
} else {
  setwd(path)
  setwd("src/seurat-analysis")
}

samples <- c("FL1_hpc", "FL2_hpc")

sample <- samples[[1]]
source("1.3_annotate-clusters.R")

setwd(path)
setwd("src/seurat-analysis")

sample <- samples[[2]]
source("1.3_annotate-clusters.R")

setwd(path)
setwd("src/seurat-analysis")
