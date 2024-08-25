
library(Seurat)
library(scater)

f <- function(x) {
  paste("geo", "GSE202374", x, sep = "/")
}

srat <- readRDS(f("merged-dataset.rds"))

# Seurat integration creates a unified object that contains both original data
# (‘RNA’ assay) as well as integrated data (‘integrated’ assay).
# if set the assay to RNA, we may directly access to raw data, and if set to
# integrated, we will access transformed data.

trim_seurat <- function(srat, assay_name, layer_name) {
  idents_df <- srat @ meta.data |> select(seurat_clusters, orig.ident,
                                          percent_mt, percent_rb,
                                          nCount_RNA, nFeature_RNA)

  export <- Seurat::CreateSeuratObject(
    counts = SeuratObject::LayerData(
      srat, assay = assay_name,
      layer = layer_name
    ),
    project = "integrated",
    meta.data = idents_df
  )

  export @ reductions $ pca <- srat @ reductions $ pca
  export @ reductions $ umap <- srat @ reductions $ umap
  return(export)
}

DefaultAssay(srat) <- "integrated"
srat <- ScaleData(srat)
srat <- RunPCA(srat, npcs = 30, verbose = F)
srat <- RunUMAP(srat, reduction = "pca", dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, k.param = 10, verbose = TRUE)
srat <- FindClusters(srat, verbose = TRUE)

# convert seurat back to single cell experiment
# note that this data is already log2-based. and that data scaling is only for
# better pca performance. when we come to differential expression, we should
# use unscaled data.
sce <- as.SingleCellExperiment(srat |> trim_seurat("integrated", "data"))
# see if the scale is successful
plotRLE(sce, exprs_values = "counts")

scedf <- counts(sce) |> as.data.frame()
rowmeds <- apply(scedf, 1, median)

global_de <- map(0:30, function(id) {
  sel <- (sce |> colData()) $ seurat_clusters == id
  subsets <- scedf[, sel]
  delta <- apply(subsets, 1, median) - rowmeds
  return(delta)
}) |> data.frame()

colnames(global_de) = paste("cluster", 0:30, sep = "")
saveRDS(global_de, f("merged-de.rds"))
