
integrated <- Seurat::ScaleData(integrated)
integrated <- Seurat::RunPCA(integrated)
gc()

integrated <- Seurat::RunTSNE(integrated, perplexity = 25)
integrated <- Seurat::RunUMAP(integrated, n.neighbors = 20, min.dist = 0.3)
TSNEPlot(integrated)

integrated <- Seurat::FindNeighbors(integrated, dims = 1:20, reduction = "pca")
integrated <- Seurat::FindClusters(integrated, resolution = 0.1)
integrated @ meta.data |> str()

library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
source("singler/init.R")

all_markers <- Seurat::FindAllMarkers(
  integrated, only.pos = TRUE,
  min.pct = 0.5, logfc.threshold = 0.5
)

dim(all_markers)
table(all_markers $ cluster)
top3_markers <- as.data.frame(
  all_markers |> group_by(cluster) |> top_n(n = 3, wt = avg_log2FC)
)
top1_markers <- as.data.frame(
  all_markers |> group_by(cluster) |> top_n(n = 1, wt = avg_log2FC)
)

Seurat::VlnPlot(object = integrated, features = top1_markers $ gene)

immgen <- celldex::ImmGenData()
upper_names <- str_to_upper(immgen @ NAMES)
gene_names <- rownames(integrated)
gene_names %in% upper_names |> table()
immgen @ NAMES <- upper_names

pred_annot <- singler(
  integrated @ assays $ integrated $ data,
  ref = immgen, labels = immgen $ label.fine,
  fine.tune = TRUE, clusters = integrated @ active.ident,
  assay.type.test = "data"
)

monaco <- celldex::MonacoImmuneData()
gene_names %in% monaco @ NAMES |> table()
pred_annot <- singler(
  integrated @ assays $ integrated $ data,
  ref = monaco, labels = monaco $ label.main,
  fine.tune = TRUE, clusters = integrated @ active.ident,
  assay.type.test = "data"
)

pred_annot

levels(integrated)
integrated <- RenameIdents(
  integrated,
  "0" = "T cells",
  "1" = "T cells",
  "2" = "Monocyte or macrophages",
  "3" = "B cells",
  "4" = "T cells",
  "5" = "Undetermined lymphocytes",
  "6" = "Dendritic cells",
  "7" = "Dendritic cells",
  "8" = "Basophils"
)

TSNEPlot(integrated)

catdata <- table(Idents(integrated), integrated $ patient)
catdata <- as.data.frame(catdata)

draw_stacked_bars <- function(data) {
  ggplot(data, aes(x = .data $ Var2, y = .data $ Freq, fill = .data $ Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    scale_fill_manual(values = brewer.pal(12, "Paired")) +
    theme(legend.title = element_blank())
}

draw_stacked_bars(catdata)

# the T cell subclones ========================================================

tsub <- subset(integrated, seurat_clusters == "0" |
                 seurat_clusters == "1" | seurat_clusters == "4")
TSNEPlot(tsub)
Idents(tsub) |> table()

tsub <- Seurat::ScaleData(tsub)
tsub <- Seurat::RunPCA(tsub)
tsub <- Seurat::RunTSNE(tsub, perplexity = 30)

tsub <- Seurat::FindNeighbors(tsub, dims = 1:30, reduction = "pca")
tsub <- Seurat::FindClusters(tsub, resolution = 0.5)

# here, if you create the seurat object with seurat v5, this method won't work
# for you. i choose to manually export the data here
tsce <- as.SingleCellExperiment(tsub)

seurat_scale_data_only <- function(srat) {
  # the meta.data items i want to keep
  idents_df <- data.frame(
    orig.ident = srat @ meta.data $ orig.ident,
    patient = srat @ meta.data $ patient,
    cluster = srat @ meta.data $ seurat_clusters
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
  return(export)
}

# here, we generate dataset only for scaled data. thus only tended to
# run clustering methods here. not for de.
tsce <- seurat_scale_data_only(tsub)
sce <- as.SingleCellExperiment(tsce)
sce <- scater::runUMAP(sce, n_neighbors = 10, min_dist = 0.1)
scater::plotUMAP(sce)
tsce <- as.Seurat(sce)

rm(tsce, sce)
gc()

tsub @ reductions $ umap <- tsce @ reductions $ UMAP
tsub @ reductions $ umap @ assay.used <- "integrated"
tsub[['umap']] |> str()
Seurat::UMAPPlot(tsub)
Seurat::TSNEPlot(tsub)
Seurat::DimPlot(tsub)

fine_annot_imm <- singler(
  tsub @ assays $ integrated $ data,
  ref = immgen, labels = immgen $ label.fine,
  fine.tune = TRUE, clusters = tsub @ active.ident,
  assay.type.test = "data"
)

fine_annot_mon <- singler(
  tsub @ assays $ integrated $ data,
  ref = monaco, labels = monaco $ label.fine,
  fine.tune = TRUE, clusters = tsub @ active.ident,
  assay.type.test = "data"
)

tsub $ seurat_clusters |> str()
fine_annot <- data.frame(
  immgen = fine_annot_imm $ labels,
  monaco = fine_annot_mon $ labels
)

#                                    immgen                      monaco
# 1  T cells (T.8MEMKLRG1-CD127+.D8.LISOVA)              Vd2 gd T cells
# 2        T cells (T.8EFF.OT1.12HR.LISOVA)          T regulatory cells
# 3  T cells (T.8MEMKLRG1-CD127+.D8.LISOVA) Effector memory CD8 T cells
# 4        T cells (T.8EFF.OT1.12HR.LISOVA)          T regulatory cells
# 5                       T cells (T.Tregs)          T regulatory cells
# 6                       T cells (T.Tregs)          T regulatory cells
# 7                       ILC (LPL.NCR+CNK)        Natural killer cells
# 8          T cells (T.8EFF.OT1.D5.VSVOVA) Effector memory CD8 T cells
# 9        T cells (T.8EFF.OT1.12HR.LISOVA) Effector memory CD8 T cells
# 10        T cells (T.8MEM.OT1.D45.LISOVA) Effector memory CD8 T cells
# 11       T cells (T.8EFF.TBET-.OT1LISOVA)          T regulatory cells
# 12       T cells (T.8EFF.OT1.12HR.LISOVA)   Follicular helper T cells
# 13       T cells (T.8EFF.OT1.48HR.LISOVA)                Plasmablasts
# 14                      T cells (T.Tregs)          T regulatory cells

levels(tsub)
tsub <- RenameIdents(
  tsub,
  "0"  = "Vd2 gd T cells",
  "1"  = "T regulatory cells",
  "2"  = "Effector memory CD8 T cells",
  "3"  = "T regulatory cells",
  "4"  = "T regulatory cells",
  "5"  = "T regulatory cells",
  "6"  = "Natural killer cells",
  "7"  = "Effector memory CD8 T cells",
  "8"  = "Effector memory CD8 T cells",
  "9"  = "Effector memory CD8 T cells",
  "10" = "T regulatory cells",
  "11" = "Follicular helper T cells",
  "12" = "Undetermined lymphocytes",
  "13" = "T regulatory cells"
)

t_all_markers <- Seurat::FindAllMarkers(
  tsub, only.pos = TRUE,
  min.pct = 0.5, logfc.threshold = 0.5
)

dim(t_all_markers)
table(t_all_markers $ cluster)
t_top3_markers <- as.data.frame(
  t_all_markers |> group_by(cluster) |> top_n(n = 3, wt = avg_log2FC)
)
t_top1_markers <- as.data.frame(
  t_all_markers |> group_by(cluster) |> top_n(n = 1, wt = avg_log2FC)
)

t_catdata <- table(Idents(tsub), tsub $ patient)
t_catdata <- as.data.frame(t_catdata)

ggplot(t_catdata,
       aes(x = .data $ Var2, y = .data $ Freq, fill = .data $ Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

Seurat::FeaturePlot(object = tsub, features = c("IL17A"))
Seurat::VlnPlot(object = tsub, features = c("IL17A"))

Seurat::FeaturePlot(object = tsub, features = t_top1_markers $ gene)
Seurat::VlnPlot(object = tsub, features = t_top1_markers $ gene)

saveRDS(integrated, "geo/GSE139324/integrated.rds")
saveRDS(tsub, "geo/GSE139324/t-cells.rds")
