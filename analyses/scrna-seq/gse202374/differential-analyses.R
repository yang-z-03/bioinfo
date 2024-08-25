
library(Seurat)
library(RColorBrewer)
library(ggplot2)

devtools::install_github("immunogenomics/presto")

gse <- "GSE202374"
gsm <- "GSM6112169_C1"
base_path <- paste("geo", gse, gsm, sep = "/")

f <- function(a) {
  paste(base_path, a, sep = "/")
}

current_data <- readRDS(f("normalized-data.rds"))

current_data <- Seurat::RunUMAP(
  current_data,
  dims = 1:30,
  verbose = TRUE,
  n.neighbors = 10,
  min.dist = 0.5
)

graph <- Seurat::DimPlot(
  current_data,
  label.size = 4,
  repel = TRUE, label = TRUE
) + unify_theme_font()

ggsave(f("umap.png"), graph, width = 8, height = 5, dpi = 100)

draw_single_gene_expr <- function(gene, description) {
  graph <- Seurat::FeaturePlot(current_data, gene, label = TRUE) +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = "Spectral"))
    ) +
    ggtitle(paste0(gene, ": ", description)) +
    unify_theme_font()

  ggsave(f(paste0("expr-", gene |> tolower(), ".png")), graph,
         width = 5, height = 5, dpi = 100)
}

draw_single_gene_expr("PPBP", "Platlets") # no such gene
draw_single_gene_expr("LILRA4", "Dendritic cells")
draw_single_gene_expr("MS4A1", "B cells")
draw_single_gene_expr("LYZ", "Monocytes")
draw_single_gene_expr("NKG7", "NK cells")
draw_single_gene_expr("CD8B", "CD8+ T cells")
draw_single_gene_expr("IL7R", "CD4+ T cells")

# 使用 SingleR 进行细胞类型标注

# 鉴于我们已经定义的标志物，我们可以挖掘文献并识别每种观察到的细胞类型（这可能是 PBMC
# 最容易的）。 然而，我们可以尝试自动注释， SingleR 是工作流不敏感的（它可以与 Seurat，
# SCE 等一起使用）。详细的 singleR 手册与先进的用法可以在这里找到。

# https://bioconductor.org/books/release/SingleRBook/introduction.html

install_packages(c("celldex", "SingleR"))

# 下面是人类细胞的几个常用参考标记集

monaco_ref <- celldex::MonacoImmuneData()
# hpca_ref <- celldex::HumanPrimaryCellAtlasData() # nolint
# dice_ref <- celldex::DatabaseImmuneCellExpressionData() # nolint

# 为了方便起见，让我们将 Seurat 对象转换为单细胞实验（SCE）

sce <- as.SingleCellExperiment(DietSeurat(current_data))

monaco_main <- SingleR::SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco_ref,
  labels = monaco_ref $ label.main
)

monaco_ref $ label.fine

#   [1] "Naive CD8 T cells"             "Central memory CD8 T cells"
#   [3] "Effector memory CD8 T cells"   "Terminal effector CD8 T cells"
#   [5] "MAIT cells"                    "Vd2 gd T cells"
#   [7] "Non-Vd2 gd T cells"            "Follicular helper T cells"
#   [9] "T regulatory cells"            "Th1 cells"
#  [11] "Th1/Th17 cells"                "Th17 cells"
#  [13] "Th2 cells"                     "Naive CD4 T cells"
#  [15] "Progenitor cells"              "Naive B cells"
#  [17] "Non-switched memory B cells"   "Exhausted B cells"
#  [19] "Switched memory B cells"       "Plasmablasts"
#  [21] "Classical monocytes"           "Intermediate monocytes"
#  [23] "Non classical monocytes"       "Natural killer cells"
#  [25] "Plasmacytoid dendritic cells"  "Myeloid dendritic cells"
#  [27] "Low-density neutrophils"       "Low-density basophils"
#  [29] "Naive CD8 T cells"             "Central memory CD8 T cells"
#  [31] "Effector memory CD8 T cells"   "Terminal effector CD8 T cells"
#  [33] "MAIT cells"                    "Vd2 gd T cells"
#  [35] "Non-Vd2 gd T cells"            "Follicular helper T cells"
#  [37] "T regulatory cells"            "Th1 cells"
#  [39] "Th1/Th17 cells"                "Th17 cells"
#  [41] "Th2 cells"                     "Naive CD4 T cells"
#  [43] "Terminal effector CD4 T cells" "Progenitor cells"
#  [45] "Naive B cells"                 "Non-switched memory B cells"
#  [47] "Exhausted B cells"             "Switched memory B cells"
#  [49] "Plasmablasts"                  "Classical monocytes"
#  [51] "Intermediate monocytes"        "Non classical monocytes"
#  [53] "Natural killer cells"          "Plasmacytoid dendritic cells"
#  [55] "Myeloid dendritic cells"       "Low-density neutrophils"
#  [57] "Low-density basophils"         "Naive CD8 T cells"
#  [59] "Central memory CD8 T cells"    "Effector memory CD8 T cells"
#  [61] "Terminal effector CD8 T cells" "MAIT cells"
#  [63] "Vd2 gd T cells"                "Non-Vd2 gd T cells"
#  [65] "Follicular helper T cells"     "T regulatory cells"
#  [67] "Th1 cells"                     "Th1/Th17 cells"
#  [69] "Th17 cells"                    "Th2 cells"
#  [71] "Naive CD4 T cells"             "Terminal effector CD4 T cells"
#  [73] "Progenitor cells"              "Naive B cells"
#  [75] "Non-switched memory B cells"   "Exhausted B cells"
#  [77] "Switched memory B cells"       "Plasmablasts"
#  [79] "Classical monocytes"           "Intermediate monocytes"
#  [81] "Non classical monocytes"       "Natural killer cells"
#  [83] "Plasmacytoid dendritic cells"  "Myeloid dendritic cells"
#  [85] "Low-density neutrophils"       "Low-density basophils"
#  [87] "Naive CD8 T cells"             "Central memory CD8 T cells"
#  [89] "Effector memory CD8 T cells"   "Terminal effector CD8 T cells"
#  [91] "MAIT cells"                    "Vd2 gd T cells"
#  [93] "Non-Vd2 gd T cells"            "Follicular helper T cells"
#  [95] "T regulatory cells"            "Th1 cells"
#  [97] "Th1/Th17 cells"                "Th17 cells"
#  [99] "Th2 cells"                     "Naive CD4 T cells"
# [101] "Progenitor cells"              "Naive B cells"
# [103] "Non-switched memory B cells"   "Exhausted B cells"
# [105] "Switched memory B cells"       "Plasmablasts"
# [107] "Classical monocytes"           "Intermediate monocytes"
# [109] "Non classical monocytes"       "Natural killer cells"
# [111] "Plasmacytoid dendritic cells"  "Myeloid dendritic cells"
# [113] "Low-density neutrophils"       "Low-density basophils"

monaco_fine <- SingleR::SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco_ref,
  labels = monaco_ref $ label.fine
)

table(monaco_main $ pruned.labels)
table(monaco_fine $ pruned.labels)

# 用 Seurat 画 UMAP 图，将分组标签导入 Seurat
current_data @ meta.data $ monaco_main <- monaco_main $ pruned.labels
current_data @ meta.data $ monaco_fine <- monaco_fine $ pruned.labels

current_data <- Seurat::SetIdent(current_data, value = "monaco_fine")

annotation <- Seurat::DimPlot(
  current_data,
  label = TRUE, repel = TRUE, label.size = 4
) + unify_theme_font()

ggsave(
  f("annotation.png"), dpi = 300, plot = annotation,
  width = 15, height = 9, units = "in"
)

current_data <- Seurat::SetIdent(current_data, value = "monaco_main")

annotation <- Seurat::DimPlot(
  current_data,
  label = TRUE, repel = TRUE, label.size = 4
) + unify_theme_font()

ggsave(
  f("annotation-coarse.png"), dpi = 300, plot = annotation,
  width = 15, height = 9, units = "in"
)
