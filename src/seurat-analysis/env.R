library(readr)

path <- Sys.getenv("PROJECT_PATH")
if (nchar(path) == 0) {
	path <- sub('\n', '', read_file("../../project_path"))
}
# print(getwd())
setwd(path)
print(paste0("> setting working directory to ", getwd()))

