
.libPaths(c(
  "/home/yang-z/R/bioinfo/4.4",
  "/usr/lib64/R/library",
  "/usr/share/R/library"
))

cat(crlf, "LOADING REQUIRED PACKAGES ...", crlf)

suppressPackageStartupMessages({
  require(crayon)
  require(argparse)
  require(dplyr)
  require(stringr)
  require(tibble)

  require(SingleCellExperiment)
  require(Seurat)
})

genes <- readRDS("~/Documents/bioinfo/refseq/9606/genes-ext.rds")
