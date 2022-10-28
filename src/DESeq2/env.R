library(readr)

path <- Sys.getenv("PROJECT_PATH")
if (nchar(path) == 0) {
  path <- sub("\n", "", read_file("../../project_path"))
}
# print(getwd())
setwd(path)
print(paste0("> changed working dir to ", getwd()))
