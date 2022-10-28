curdir <- getwd()
path <- Sys.getenv("PROJECT_PATH")

if (path == "") {
  source("env.R")
  path <- getwd()
  setwd("src/DESeq2")
} else {
  setwd(path)
  setwd("src/DESeq2")
}

script <- "2.1_investigate_FL_clusters_in_PCA.R"
print("")
print(paste0("> Running script ", script))
source(script)
setwd(path)
setwd("src/DESeq2")
script <- "2.2_investigate_FL_gated_in_PCA.R"
print("")
print(paste0("> Running script ", script))
source(script)
setwd(path)
setwd("src/DESeq2")
script <- "2.3_investigate_FL_gates_n_clusters_in_PCA.R"
print("")
print(paste0("> Running script ", script))
source(script)
