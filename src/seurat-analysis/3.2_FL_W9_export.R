library(Seurat)
library(dplyr)
library(ggplot2)

source("env.R")

# Global variables --------------------------------------------------------

input.dir <- paste0("data/processed/seurat/W9")

output.dir <- input.dir
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

seed <- 1234
set.seed(seed)

# Load data ---------------------------------------------------------------

FL_W9 <- readRDS(paste0(input.dir, "/FL_CS16_and_W9_subset.rds"))

# Save data ---------------------------------------------------------------

write.table(
            FL_W9@meta.data,
            paste0(input.dir, "/FL_W9_metadata.csv"), sep = ","
)
write.table(
            FL_W9@assays$ADT@scale.data,
            paste0(input.dir, "/FL_W9_adt_scaled_values.csv"), sep = ","
)
