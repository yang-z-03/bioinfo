
library(Seurat)

# 创建六个样本的 Seurat

get_sample_dir <- function(x) {
  paste("geo", "GSE202374", x, sep = "/")
}

f <- function(x) {
  paste("geo", "GSE202374", x, sep = "/")
}

read_sample <- function(x) {
  df <- Seurat::Read10X(
    # gene.column: Specify which column of genes.tsv or features.tsv to use
    # for gene names; default is 2 (gene symbol), or 1 (ensembl id)
    gene.column = 2,
    # cell.column: Specify which column of barcodes.tsv to use for cell names;
    # default is 1
    cell.column = 1,
    data.dir = get_sample_dir(x)
  )

  srat <- Seurat::CreateSeuratObject(
    counts = df, project = x,
    min.cells = 3, min.features = 200
  )

  return(srat)
}

c1 <- read_sample("gsm6112169_c1")
c23 <- read_sample("gsm6112170_c23")

p6 <- read_sample("gsm6112171_p6")
p10 <- read_sample("gsm6112172_p10")
p14 <- read_sample("gsm6112173_p14")
p15 <- read_sample("gsm6112174_p15")

calc_mt_rb <- function(srat, verbose = FALSE) {
  srat[["percent_mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
  srat[["percent_rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

  if (verbose) {
    VlnPlot(srat, features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent_mt",
      "percent_rb"
    ), ncol = 4)
  }

  return(srat)
}

normalize_srat <- function(srat, lower_feature, upper_feature, upper_mt) {
  srat <- subset(srat,
    subset = nFeature_RNA > lower_feature &
      nFeature_RNA < upper_feature &
      percent_mt < upper_mt
  )

  srat <- Seurat::NormalizeData(srat, verbose = TRUE)
  srat <- Seurat::FindVariableFeatures(
    srat, selection.method = "vst", nfeatures = 2000, verbose = TRUE
  )

  return(srat)
}

c1  <- calc_mt_rb(c1)
c23 <- calc_mt_rb(c23)
p6  <- calc_mt_rb(p6)
p10 <- calc_mt_rb(p10)
p14 <- calc_mt_rb(p14)
p15 <- calc_mt_rb(p15)

c1  <- normalize_srat(c1, 700, 5000, 15)
c23 <- normalize_srat(c23, 500, 5000, 15)
p6  <- normalize_srat(p6, 700, 5000, 15)
p10 <- normalize_srat(p10, 500, 5000, 15)
p14 <- normalize_srat(p14, 500, 5000, 12)
p15 <- normalize_srat(p15, 500, 5000, 12)

sample_list <- list(c1, c23, p6, p10, p14, p15)

anchors <- Seurat::FindIntegrationAnchors(
  object.list = sample_list, dims = 1:30
)

srat <- Seurat::IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(srat, f("merged-dataset.rds")) # ======================================
srat <- readRDS(f("merged-dataset.rds"))

# clean up
rm(c1, c23, p10, p14, p15, p6, anchors)
gc()

DefaultAssay(srat) <- "RNA"
srat <- Seurat::NormalizeData(srat, verbose = TRUE)
srat <- Seurat::FindVariableFeatures(
  srat, selection.method = "vst",
  nfeatures = 2000, verbose = TRUE
)
srat <- Seurat::ScaleData(srat, verbose = TRUE)
srat <- Seurat::RunPCA(srat, npcs = 30, verbose = TRUE)

srat <- Seurat::RunUMAP(
  srat, reduction = "pca",
  dims = 1:30,
  verbose = TRUE,
  n.neighbors = 10,
  min.dist = 0.5
)

# Now let's change the assay to integrated and do the same do the same thing in
# the integrated assay (it's already normalized and HVGs are selected):

DefaultAssay(srat) <- "integrated"
srat <- Seurat::ScaleData(srat, verbose = TRUE)
srat <- Seurat::RunPCA(srat, npcs = 30, verbose = TRUE)

srat <- Seurat::RunUMAP(
  srat, reduction = "pca",
  dims = 1:30,
  verbose = TRUE,
  n.neighbors = 10,
  min.dist = 0.5
)

merged_dim <- Seurat::DimPlot(
  srat, reduction = "umap"
) +
  patchwork::plot_annotation(title = "merged dataset (c2, p4)") +
  unify_theme_font()

ggsave(
  f("merged-umap.png"), dpi = 300, plot = merged_dim,
  width = 9, height = 9, units = "in"
)

# 使用 Seurat 自动聚类

srat <- FindNeighbors(srat, dims = 1:30, k.param = 10)
srat <- FindClusters(srat)
table(srat @ meta.data $ seurat_clusters)

dim_empirical <- DimPlot(srat, label = TRUE) +
  unify_theme_font()

ggsave(
  f("umap-cluster.png"), dpi = 300, plot = dim_empirical,
  width = 9, height = 9, units = "in"
)

# 获取 DE 列表，以备后用

all_markers <- Seurat::FindAllMarkers(
  srat, only.pos = TRUE,
  min.pct = 0.5, logfc.threshold = 0.5
)

dim(all_markers)
table(all_markers $ cluster)
top3_markers <- as.data.frame(
  all_markers |> group_by(cluster) |>
    top_n(n = 3, wt = avg_log2FC)
)

saveRDS(all_markers, f("merged-differential-expressions.rds"))

# =============================================================================

DefaultAssay(srat) <- "integrated"
monaco_ref <- celldex::MonacoImmuneData()

# 我想只保留 integrated 中的 scale.data 层，它本可以使用这句话
# srat |> DietSeurat(layers = "scale.data") 实现
# 然而，Seurat 5.1.0 这个有 bug，这样并不能去掉层，我只能

seurat_scale_data_only <- function(srat) {
  # the meta.data items i want to keep
  idents_df <- data.frame(
    orig.ident = srat @ meta.data $ orig.ident
  )

  export <- Seurat::CreateSeuratObject(
    counts = SeuratObject::LayerData(
      srat, assay = "integrated",
      layer = "scale.data"
    ),
    project = "integrated",
    meta.data = idents_df
  )

  export @ reductions $ pca <- srat @ reductions $ pca
  export @ reductions $ umap <- srat @ reductions $ umap

  return(export)
}

sce <- as.SingleCellExperiment(srat |> seurat_scale_data_only())

monaco_main <- SingleR::SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco_ref,
  labels = monaco_ref $ label.main
)

monaco_fine <- SingleR::SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco_ref,
  labels = monaco_ref $ label.fine
)

monaco_main $ pruned.labels |> length()
monaco_main $ pruned.labels |> table()

srat @ meta.data $ monaco_main <- monaco_main $ pruned.labels
srat @ meta.data $ monaco_fine <- monaco_fine $ pruned.labels

srat <- Seurat::SetIdent(srat, value = "monaco_fine")

annotation <- Seurat::DimPlot(
  srat,
  label = TRUE, repel = TRUE, label.size = 4
) + unify_theme_font()

ggsave(
  f("annotation.png"), dpi = 300, plot = annotation,
  width = 15, height = 9, units = "in"
)

srat <- Seurat::SetIdent(srat, value = "monaco_main")

annotation <- Seurat::DimPlot(
  srat,
  label = TRUE, repel = TRUE, label.size = 4
) + unify_theme_font()

ggsave(
  f("annotation-coarse.png"), dpi = 300, plot = annotation,
  width = 9, height = 9, units = "in"
)

# the sample groupings
srat @ meta.data $ orig.ident |> table()

grouping <- srat @ meta.data $ orig.ident
is_control <- grep("_c[0-9]*$", grouping)
is_patient <- grep("_p[0-9]*$", grouping)
grouping[is_control] <- "control"
grouping[is_patient] <- "patient"
grouping |> table()

# the control/patient grouping tag
srat @ meta.data $ group.ident <- grouping

srat <- Seurat::SetIdent(srat, value = "group.ident")

annotation <- Seurat::DimPlot(
  srat,
  label = TRUE, repel = TRUE, label.size = 4
) + unify_theme_font()

ggsave(
  f("annotation-group.png"), dpi = 300, plot = annotation,
  width = 9, height = 9, units = "in"
)
