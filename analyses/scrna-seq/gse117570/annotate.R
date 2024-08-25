
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
source("singler/init.R")
source("utils/themes.R")

f <- function(t) {
  paste("geo", "GSE117570", t, sep = "/")
}

srat <- readRDS(f("integrated.rds"))

# re-cluster programs.
srat <- Seurat::FindNeighbors(srat, dims = 1:2, reduction = "tsne")
srat <- Seurat::FindClusters(srat, resolution = 0.1)
srat @ meta.data |> str()

# find markers according to cluster
all_markers <- Seurat::FindAllMarkers(
  srat, only.pos = TRUE,
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

immgen <- celldex::ImmGenData()
upper_names <- str_to_upper(immgen @ NAMES)
gene_names <- rownames(srat)
gene_names %in% upper_names |> table()

immgen @ NAMES <- upper_names

pred_annot <- singler(
  srat @ assays $ integrated $ data,
  ref = immgen, labels = immgen $ label.fine,
  fine.tune = TRUE, clusters = srat @ active.ident,
  assay.type.test = "data"
)

monaco <- celldex::MonacoImmuneData()
gene_names %in% monaco @ NAMES |> table()
pred_annot <- singler(
  srat @ assays $ integrated $ data,
  ref = monaco, labels = monaco $ label.main,
  fine.tune = TRUE, clusters = srat @ active.ident,
  assay.type.test = "data"
)

hca <- celldex::HumanPrimaryCellAtlasData()
gene_names %in% hca @ NAMES |> table()
pred_annot <- singler(
  srat @ assays $ integrated $ data,
  ref = hca, labels = hca $ label.main,
  fine.tune = TRUE, clusters = srat @ active.ident,
  assay.type.test = "data"
)

TSNEPlot(srat)
Idents(srat) |> table()

levels(srat)
srat <- RenameIdents(
  srat,
  "0" = "T cells",
  "1" = "Monocytes",
  "2" = "Macrophages",
  "3" = "Epithelial cells",
  "4" = "Epithelial cells",
  "5" = "Epithelial cells",
  "6" = "Monocytes",
  "7" = "Epithelial cells",
  "8" = "Stem cell-likes",
  "9" = "B cells"
)

catdata <- table(Idents(srat), srat $ patient)
catdata <- as.data.frame(catdata)

bar_all <- ggplot(catdata,
                  aes(x = .data $ Var2, y = .data $ Freq, fill = .data $ Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

# we may remove the epithelial cell fraction. and watch only the immune cell
# partition and its distributions.

immune_catdata <- catdata[catdata $ Var1 != "Epithelial cells", ]

bar_immune <- ggplot(immune_catdata,
                     aes(x = .data $ Var2, y = .data $ Freq, fill = .data $ Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

tsub <- subset(srat, seurat_clusters == "0")
TSNEPlot(tsub)
Idents(tsub) |> table()

# recluster of the T cell portions.

tsub <- Seurat::ScaleData(tsub)
tsub <- Seurat::RunPCA(tsub)
tsub <- Seurat::RunTSNE(tsub, perplexity = 20)

tsub <- Seurat::FindNeighbors(tsub, dims = 1:30, reduction = "pca")
tsub <- Seurat::FindClusters(tsub, resolution = 1)

TSNEPlot(tsub)

saveRDS(tsub, f('t-cells.rds'))
tsub <- readRDS(f('t-cells.rds'))

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

t_catdata <- table(Idents(tsub), tsub $ patient)
t_catdata <- as.data.frame(t_catdata)

bar_t <- ggplot(t_catdata,
                aes(x = .data $ Var2, y = .data $ Freq, fill = .data $ Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

Seurat::FeaturePlot(object = tsub, features = c("IL17RA"))

tsub <- RenameIdents(
  tsub,
  "0" = "T Th1/Th17",
  "1" = "T regulatory",
  "2" = "T Switched memory",
  "3" = "T CD8+",
  "4" = "ILC (LPL.NCR+CNK)",
  "5" = "ILC (ILC2)"
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

Seurat::FeaturePlot(object = tsub, features = t_top1_markers $ gene)
Seurat::VlnPlot(object = tsub, features = t_top1_markers $ gene)

violins <- list()
expr_umaps <- list()
# genes_to_draw <- t_top1_markers $ gene
genes_to_draw <- c("ANXA1", "DUSP4", "ARHGAP24", "SLC2A3", "ITM2C", "PTGS2")
n <- 1
for (g in genes_to_draw) {
  violins[[n]] <- Seurat::VlnPlot(object = tsub, features = g, pt.size = 0.01,
                                  alpha = 0.3, combine = FALSE) [[1]] +
    unify_theme_font(9, "Arial") +
    Seurat::NoLegend() +
    ggplot2::scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
    ggplot2::xlab("Cluster")

  expr_umaps[[n]] <- Seurat::FeaturePlot(object = tsub, features = g,
                                         combine = FALSE) [[1]] +
    unify_theme_font(9, "Arial") +
    # Seurat::NoLegend() +
    ggplot2::xlab("tSNE 1") +
    ggplot2::ylab("tSNE 2") +
    scale_color_gradient(low = "white", "high" = "red")

  n <- n + 1
}

violin_leg <- Seurat::VlnPlot(object = tsub, features = g, pt.size = 0.01,
                              alpha = 0.3, combine = FALSE) [[1]] +
  unify_theme_font(9, "Arial") +
  ggplot2::scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
  ggplot2::xlab("Cluster")

pw_violin <-
  (violins[[1]] | violins[[2]] | violins[[3]]) /
  (violins[[4]] | violins[[5]] | violins[[6]])

pw_umaps <-
  (expr_umaps[[1]] | expr_umaps[[2]] | expr_umaps[[3]]) /
  (expr_umaps[[4]] | expr_umaps[[5]] | expr_umaps[[6]])


ggsave(f("violins.pdf"), pw_violin, width = 7.5, height = 4, units = "in")
ggsave(f("violins.png"), pw_violin, width = 7.5, height = 4, units = "in", dpi = 600)
ggsave(f("violin-legend.pdf"), violin_leg, width = 5, height = 2, units = "in")
ggsave(f("tsnes.pdf"), pw_umaps, width = 10, height = 6, units = "in")
ggsave(f("tsnes.png"), pw_umaps, width = 10, height = 6, units = "in", dpi = 600)

tsne_int <- DimPlot(srat)
umap_tsub <- DimPlot(tsub)

ggsave(f("tsne-all.pdf"), tsne_int, width = 6, height = 4, units = "in")
ggsave(f("tsne-all.png"), tsne_int, width = 6, height = 4, units = "in", dpi = 600)

ggsave(f("tsne-t.pdf"), umap_tsub, width = 6, height = 4, units = "in")
ggsave(f("tsne-t.png"), umap_tsub, width = 6, height = 4, units = "in", dpi = 600)

ggsave(f("stack-t.pdf"), bar_t, width = 5, height = 4, units = "in")
ggsave(f("stack-immune.pdf"), bar_immune, width = 5, height = 4, units = "in")
ggsave(f("stack-all.pdf"), bar_all, width = 5, height = 4, units = "in")
