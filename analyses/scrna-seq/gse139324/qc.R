
library(scater)
source("scrublet/scrublet.R")

qc_violin <- function(current_data) {
  Seurat::VlnPlot(
    current_data,
    features = c(
      "reads",
      "percent_mt",
      "percent_rb"
    ),
    ncol = 3, pt.size = 0.1
  )
}

qc <- function(srat) {
  srat <- Seurat::NormalizeData(srat)

  assaydata <- GetAssayData(srat)
  genes <- srat |> rownames()
  mtgenes <- genes[grep("^MT-", genes, ignore.case = TRUE)]
  srat $ percent_mt <- colSums(assaydata[mtgenes,]) / colSums(assaydata)

  rbgenes <- genes[c(
    grep("^RPL", genes, ignore.case = TRUE),
    grep("^RPS", genes, ignore.case = TRUE)
  )]
  srat $ percent_rb <- colSums(assaydata[rbgenes,]) / colSums(assaydata)
  srat $ reads <- colSums(assaydata)

  srat <- scrublet_matrix(seurat_obj = srat,
                          matrix = srat @ assays $ RNA $ data)

  # list(qc_violin(srat), srat)
  return(srat)
}

for (i in seq_along(plist)) {
  cat("calculating qc metrics for", i, "\n")
  plist[[i]] <- qc(plist[[i]])
}

genes <- rownames(plist[[1]])
genes[grep("cd3", genes, ignore.case = TRUE)]

for (i in seq_along(plist)) {
  cat("preparing and scaling for", i, "\n")
  len <- rownames(plist[[i]]) |> length()
  plist[[i]] $ patient <- rep(i, len)
  plist[[i]] <- Seurat::ScaleData(plist[[i]])
  plist[[i]] <- FindVariableFeatures(plist[[i]], nfeatures = 10000,
                                     selection.method = "vst")
  plist[[i]] <- Seurat::RunPCA(plist[[i]])
}

merge1 <- list()
for (i in 1:13) {
  merge1[[i]] <- plist[[i]]
}

features <- Seurat::SelectIntegrationFeatures(plist, nfeatures = 3000,
                                              fvf.nfeatures = 3000)
anchors <- Seurat::FindIntegrationAnchors(
  plist, k.anchor = 5, k.filter = 25, dims = 1:10,
  anchor.features = c(
    features
    #, "CD3D", "CD3E", "CD3G", "IL17A", "IL17F", "IL23A",
    # genes[grep("hla", genes, ignore.case = TRUE)],
    # genes[grep("ifn", genes, ignore.case = TRUE)]
  )
)

saveRDS(anchors, "anchors.rds")
anchors <- readRDS("anchors.rds")
genes <- readRDS("features.rds")
features <- c(genes, "CD3D", "CD3E", "CD3G", "IL17A", "IL17F", "IL23A")

integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:10,
                                    features = features,
                                    features.to.integrate = features)

saveRDS(integrated, "integrated.rds")
saveRDS(genes, "features.rds")
