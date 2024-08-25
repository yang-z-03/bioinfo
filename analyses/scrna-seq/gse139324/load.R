
library(Seurat)
library(SingleCellExperiment)

patients <- c(
  "geo/GSE139324/GSM4138111_HNSCC_1/TIL/",
  "geo/GSE139324/GSM4138113_HNSCC_2/TIL/",
  "geo/GSE139324/GSM4138115_HNSCC_3/TIL/",
  "geo/GSE139324/GSM4138117_HNSCC_4/TIL/",
  "geo/GSE139324/GSM4138119_HNSCC_5/TIL/",
  "geo/GSE139324/GSM4138121_HNSCC_6/TIL/",
  "geo/GSE139324/GSM4138123_HNSCC_7/TIL/",
  "geo/GSE139324/GSM4138125_HNSCC_8/TIL/",
  "geo/GSE139324/GSM4138127_HNSCC_9/TIL/",
  "geo/GSE139324/GSM4138129_HNSCC_10/TIL/",
  "geo/GSE139324/GSM4138131_HNSCC_11/TIL/",
  "geo/GSE139324/GSM4138133_HNSCC_12/TIL/",
  "geo/GSE139324/GSM4138135_HNSCC_13/TIL/",
  "geo/GSE139324/GSM4138137_HNSCC_14/TIL/",
  "geo/GSE139324/GSM4138139_HNSCC_15/TIL/",
  "geo/GSE139324/GSM4138141_HNSCC_16/TIL/",
  "geo/GSE139324/GSM4138143_HNSCC_17/TIL/",
  "geo/GSE139324/GSM4138145_HNSCC_18/TIL/",
  "geo/GSE139324/GSM4138147_HNSCC_19/TIL/",
  "geo/GSE139324/GSM4138149_HNSCC_20/TIL/",
  "geo/GSE139324/GSM4138151_HNSCC_21/TIL/",
  "geo/GSE139324/GSM4138153_HNSCC_22/TIL/",
  "geo/GSE139324/GSM4138155_HNSCC_23/TIL/",
  "geo/GSE139324/GSM4138157_HNSCC_24/TIL/",
  "geo/GSE139324/GSM4138159_HNSCC_25/TIL/",
  "geo/GSE139324/GSM4138161_HNSCC_26/TIL/"
)

read_files_10x <- function(patients) {
  seurat_list <- list()

  id <- 1
  for (sample in patients) {

    matr <- Seurat::Read10X(

      # Specify which column of genes.tsv or features.tsv to use
      # for gene names; default is 2 (gene symbol), or 1 (ensembl id)
      gene.column = 2,

      # Specify which column of barcodes.tsv to use for cell names;
      cell.column = 1,
      data.dir = sample
    )

    seurat_list[[id]] <- Seurat::CreateSeuratObject(
      counts = matr, project = paste0("rt-", id),
      min.cells = 3, min.features = 200
    )

    cat("reading sample", id, "\n")
    id <- id + 1
  }

  return(seurat_list)
}

plist <- read_files_10x(patients)
