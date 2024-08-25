
install_packages(c("Seurat", "glmGamPoi"))
devtools::install_github("Moonerss/scrubletR")
devtools::install_github("constantAmateur/SoupX", ref = "devel")

library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)
library(scales)

library(SingleCellExperiment)
library(scater)
library(scran)

library(Seurat)
library(patchwork)
library(fs)
library(devtools)
library(reticulate)
library(RColorBrewer)

reticulate::use_condaenv(
  condaenv = "bioinfo",
  conda = "D:/applications/yang-z/conda/condabin/conda.bat"
)

reticulate::py_config()

series_id <- "GSE202374"

# download_table <- geo_get_supplemental_links(series_id) # nolint
# geo_get_supplemental_files(series_id)                   # nolint

expr_set <- geo_get_series_annot(series_id)
platform <- geo_get_platform_annot(expr_set)

# 读取文件，创建 Seurat 对象

root_dir <- paste("geo", series_id, sep = "/")

arrange_files_to_10x <- function(root_dir) {
  files <- list.files(root_dir)
  matrix <- strsplit(files, split = "_")

  for (x in matrix) {
    if (length(x) != 3) next

    attempt_dir <- paste(root_dir, paste(x[1], x[2], sep = "_"), sep = "/")
    if (!dir.exists(attempt_dir)) {
      dir.create(attempt_dir)
    }

    file_move(
      paste(root_dir, paste(x[1], x[2], x[3], sep = "_"), sep = "/"),
      paste(root_dir, paste(x[1], x[2], sep = "_"), x[3], sep = "/")
    )
  }
}

arrange_files_to_10x(root_dir)

read_files_10x <- function(root_dir) {
  rt_seurat_list <- list()
  sample_list <- list.dirs(root_dir, recursive = FALSE)

  id <- 1
  for (sample in sample_list) {
    rt <- Seurat::Read10X(
      # gene.column: Specify which column of genes.tsv or features.tsv to use
      # for gene names; default is 2 (gene symbol), or 1 (ensembl id)
      gene.column = 2,
      # cell.column: Specify which column of barcodes.tsv to use for cell names;
      # default is 1
      cell.column = 1,
      data.dir = sample
    )

    soup_channel <- SoupX::SoupChannel(tod = rt, toc = rt,
                                       calcSoupProfile = FALSE)

    # SoupChannel 对象的构成：

    #  $ toc        :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #   ..@ i       : int [1:11293999] 39 53 ...
    #   ..@ p       : int [1:8219] 0 1488 ...
    #   ..@ Dim     : int [1:2] 33538 8218
    #   ..@ Dimnames:List of 2
    #   .. ..$ : chr [1:33538] "MIR1302-2HG" "FAM138A" ...
    #   .. ..$ : chr [1:8218] "AAACCCAAGCGATGGT-1" ...
    #   ..@ x       : num [1:11293999] 1 1 2 1 3 9 1 1 2 1 ...
    #   ..@ factors : list() ·
    #  $ metaData   :'data.frame':    8218 obs. of  1 variable:
    #   ..$ nUMIs: num [1:8218] 4254 2238 4009 6294 1584 ...
    #  $ nDropUMIs  : Named num [1:8218] 4254 2238 4009 6294 1584 ...
    #   ..- attr(*, "names")= chr [1:8218] "AAACCCAAGCGATGGT-1" ...
    #  $ soupProfile:'data.frame':    33538 obs. of  2 variables:
    #   ..$ est   : num [1:33538] NaN NaN 。..
    #   ..$ counts: num [1:33538] 0 0 0 0 0 0 0 0 0 0 ...
    #  - attr(*, "class")= chr [1:2] "list" "SoupChannel"

    soup_profile <- data.frame(
      row.names = rownames(soup_channel $ toc),
      est = rowSums(soup_channel $ toc) / sum(soup_channel $ toc),
      counts = rowSums(soup_channel $ toc)
    )
    soup_channel <- setSoupProfile(soup_channel, soup_profile)

    # soupx 要求我们必须有一个分群，但是我连 QC 都没做我怎么分群？
    # 我先随便分一个，喂给它试试

    qclust <- scran::quickCluster(rt, min.size = 30)
    soup_channel <- SoupX::setClusters(soup_channel, setNames(
      qclust, rt |> colnames()
    ))

    # estimate and clean the data using SoupX
    soup_channel <- SoupX::autoEstCont(soup_channel)
    soup_channel <- SoupX::adjustCounts(soup_channel)

    # List of 6
    #  $ tod        :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #   .. ..@ i        : int [1:11293999]
    #   .. ..@ p        : int [1:8219]
    #   .. ..@ Dim      : int [1:2]
    #   .. ..@ Dimnames :List of 2
    #   .. .. ..$ : chr [1:33538]
    #   .. .. ..$ : chr [1:8218]
    #   .. ..@ x       : num [1:11293999]
    #   .. ..@ factors : list()
    #  $ toc        :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #   .. ..@ i       : int [1:11293999]
    #   .. ..@ p       : int [1:8219]
    #   .. ..@ Dim     : int [1:2]
    #   .. ..@ Dimnames:List of 2
    #   .. .. ..$ : chr [1:33538]
    #   .. .. ..$ : chr [1:8218]
    #   .. ..@ x       : num [1:11293999]
    #   .. ..@ factors : list()
    #  $ metaData   :'data.frame':    8218 obs. of  3 variables:
    #   ..$ nUMIs   : num [1:8218] ·
    #   ..$ clusters: chr [1:8218] ·
    #   ..$ rho     : num [1:8218] ·
    #  $ nDropUMIs  : Named num [1:8218]
    #   ..- attr(*, "names")= chr [1:8218]
    #  $ soupProfile:'data.frame':    33538 obs. of  2 variables:
    #   ..$ est   : num [1:33538] ·
    #   ..$ counts: num [1:33538] ·
    #  $ fit        :List of 7
    #   ..$ dd            :'data.frame':      3000 obs. of  14 variables:
    #   .. ..$ gene          : chr [1:3000] "CPA3" "NDUFA4L2" ...
    #   .. ..$ passNonExp    : logi [1:3000] TRUE TRUE TRUE TRUE TRUE TRUE ...
    #   .. ..$ rhoEst        : num [1:3000] 0 0.03278 0.00555 ...
    #   .. ..$ rhoIdx        : int [1:3000] 1 17 12 11 15 1 24 15 13 1 ...
    #   .. ..$ obsCnt        : num [1:3000] 0 3 1 1 1 0 3 1 3 0 ...
    #   .. ..$ expCnt        : num [1:3000] 76.8 91.5 180.2 197.4 46.1 ...
    #   .. ..$ isExpressedFDR: num [1:3000] 1.16e-26 5.95e-27 3.90e-60 ...
    #   .. ..$ geneIdx       : int [1:3000] 5 7 13 17 19 21 22 30 36 39 ...
    #   .. ..$ tfidf         : num [1:3000] 3.5 3.49 3.43 3.39 3.38 ...
    #   .. ..$ soupIdx       : int [1:3000] 1365 1110 487 419 2380 866 ...
    #   .. ..$ soupExp       : num [1:3000] 9.29e-05 1.11e-04 2.18e-04 ...
    #   .. ..$ useEst        : logi [1:3000] TRUE TRUE TRUE TRUE TRUE TRUE ...
    #   .. ..$ rhoHigh       : num [1:3000] 0.0481 0.0958 0.0309 0.0282 ...
    #   .. ..$ rhoLow        : num [1:3000] 0 0.00676 0.00014 0.000128 ...
    #   ..$ priorRho      : num 0.05
    #   ..$ priorRhoStdDev: num 0.1
    #   ..$ posterior     : num [1:1001] 0 8.57 9.7 10.47 11.14 ...
    #   ..$ rhoEst        : num 0.01
    #   ..$ rhoFWHM       : num [1:2] 0.01 0.037
    #   ..$ markersUsed   :'data.frame':      3450 obs. of  10 variables:
    #   .. ..$ gene                       : chr [1:3450] "S100A1" "CMBL" ...
    #   .. ..$ cluster                    : chr [1:3450] "29" "29" "29" "2" ...
    #   .. ..$ geneFrequency              : num [1:3450] 0.804 0.891 0.826 ...
    #   .. ..$ geneFrequencyOutsideCluster: num [1:3450] 0.00649 0.01346 ...
    #   .. ..$ geneFrequencySecondBest    : num [1:3450] 0.559 0.603 0.603 ...
    #   .. ..$ geneFrequencyGlobal        : num [1:3450] 0.011 0.0184 0.0138 ...
    #   .. ..$ secondBestClusterName      : chr [1:3450] "12" "12" "12" "10" ...
    #   .. ..$ tfidf                      : num [1:3450] 3.63 3.56 3.54 ...
    #   .. ..$ idf                        : num [1:3450] 4.51 4 4.29 3.65 ...
    #   .. ..$ qval                       : num [1:3450] 6.25e-64 4.03e-64 ...
    #  - attr(*, "class")= chr [1:2] "list" "SoupChannel"

    segs <- strsplit(sample, split = "/")[[1]]
    rt_seurat_list[[id]] <- Seurat::CreateSeuratObject(
      counts = soup_channel, project = paste0("rt-", segs[3]),
      min.cells = 3, min.features = 200
    )

    id <- id + 1
  }

  return(rt_seurat_list)
}

data <- read_files_10x(root_dir)

# 注意到，读入 soup_channel 中的矩阵大小和 Seurat 中的矩阵大小不一致？
# 可能在 Seurat::CreateSeuratObject 中忽略掉了一些行和列！

# [[1]]
# An object of class Seurat
# 20140 features across 8126 samples within 1 assay
# Active assay: RNA (20140 features, 0 variable features)
#  1 layer present: counts
#
# [[2]]
# An object of class Seurat
# 20723 features across 8528 samples within 1 assay
# Active assay: RNA (20723 features, 0 variable features)
#  1 layer present: counts
#
# [[3]]
# An object of class Seurat
# 20459 features across 10637 samples within 1 assay
# Active assay: RNA (20459 features, 0 variable features)
#  1 layer present: counts
#
# [[4]]
# An object of class Seurat
# 19742 features across 4088 samples within 1 assay
# Active assay: RNA (19742 features, 0 variable features)
#  1 layer present: counts
#
# [[5]]
# An object of class Seurat
# 19478 features across 5025 samples within 1 assay
# Active assay: RNA (19478 features, 0 variable features)
#  1 layer present: counts
#
# [[6]]
# An object of class Seurat
# 19672 features across 4972 samples within 1 assay
# Active assay: RNA (19672 features, 0 variable features)
#  1 layer present: counts

# 每一个对象的结构如下

str(data[[1]])
saveRDS(data, paste(root_dir, "data-soupx-corrected.rds", sep = "/"))

# Formal class 'Seurat' [package "SeuratObject"] with 13 slots
# ..@ assays      :List of 1
# .. ..$ RNA:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
# .. .. .. ..@ layers    :List of 1
# .. .. .. .. ..$ counts:Formal class 'dgCMatrix' [package "Matrix"] with 6 slot
# .. .. .. .. .. .. ..@ i       : int [1:11278290] 22 36 ...
# .. .. .. .. .. .. ..@ p       : int [1:8127] 0 1488 ...
# .. .. .. .. .. .. ..@ Dim     : int [1:2] 20140 8126
# .. .. .. .. .. .. ..@ Dimnames:List of 2
# .. .. .. .. .. .. .. ..$ : NULL
# .. .. .. .. .. .. .. ..$ : NULL
# .. .. .. .. .. .. ..@ x       : num [1:11278290] 1 1 2 1 3 9 1 1 2 1 ...
# .. .. .. .. .. .. ..@ factors : list()
# .. .. .. ..@ cells     :Formal class 'LogMap' ["SeuratObject"] with 1 slot
# .. .. .. .. .. ..@ .Data: logi [1:8126, 1] TRUE TRUE TRUE TRUE TRUE TRUE ...
# .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. .. .. .. .. ..$ : chr [1:8126] "AAACCCAAGCGATGGT-1" ...
# .. .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. .. .. ..$ dim     : int [1:2] 8126 1
# .. .. .. .. .. ..$ dimnames:List of 2
# .. .. .. .. .. .. ..$ : chr [1:8126] "AAACCCAAGCGATGGT-1" ...
# .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. ..@ features  :Formal class 'LogMap' ["SeuratObject"] with 1 slot
# .. .. .. .. .. ..@ .Data: logi [1:20140, 1] TRUE TRUE TRUE TRUE TRUE TRUE ...
# .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. .. .. .. .. ..$ : chr [1:20140] "AL627309.1" ...
# .. .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. .. .. ..$ dim     : int [1:2] 20140 1
# .. .. .. .. .. ..$ dimnames:List of 2
# .. .. .. .. .. .. ..$ : chr [1:20140] "AL627309.1" ...
# .. .. .. .. .. .. ..$ : chr "counts"
# .. .. .. ..@ default   : int 1
# .. .. .. ..@ assay.orig: chr(0)
# .. .. .. ..@ meta.data :'data.frame': 20140 obs. of  0 variables
# .. .. .. ..@ misc      :List of 1
# .. .. .. .. ..$ calcN: logi TRUE
# .. .. .. ..@ key       : chr "rna_"
# ..@ meta.data   :'data.frame':        8126 obs. of  3 variables:
# .. ..$ orig.ident  : Factor w/ 1 level "rt": 1 1 1 1 1 1 1 1 1 1 ...
# .. ..$ nCount_RNA  : num [1:8126] 4254 2237 4009 6294 1584 ...
# .. ..$ nFeature_RNA: int [1:8126] 1488 1114 ...
# ..@ active.assay: chr "RNA"
# ..@ active.ident: Factor w/ 1 level "rt": 1 1 1 1 1 1 1 1 1 1 ...
# .. ..- attr(*, "names")= chr [1:8126] "AAACCCAAGCGATGGT-1" ...
# ..@ graphs      : list() ·
# ..@ neighbors   : list() ·
# ..@ reductions  : list() ·
# ..@ images      : list() ·
# ..@ project.name: chr "rt"
# ..@ misc        : list() ·
# ..@ version     :Classes 'package_version', 'numeric_version' hidden list of 1
# .. ..$ : int [1:3] 5 0 2
# ..@ commands    : list() ·
# ..@ tools       : list() ·

# meta.data 是下一步最重要的字段。它可以使用 @ 和 [[ 操作符访问。
# 目前，每个细胞有 3 个字段：数据集 ID、每个细胞检测到的 UMI 读数数（nCount_RNA）和
# 每个相同细胞表达（检测到）的基因数（nFeature_RNA）。

data <- readRDS(paste(root_dir, "data-soupx-corrected.rds", sep = "/"))
x <- 1

for (x in seq_along(data)) {
  current_data <- data[[x]]

  meta <- current_data @ meta.data
  head(meta)

  summary(meta $ nCount_RNA)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  （未使用 SoupX）
  #   500    1998    3276    3799    4522   44178
  #
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  （使用 SoupX）
  # 497.4  1979.1  3244.5  3760.8  4478.1 43686.9

  summary(meta $ nFeature_RNA)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 200.0   901.2  1343.0  1387.9  1710.0  5522.0

  # 计算线粒体基因（用名称作为索引）和核糖体基因的占比
  current_data[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(current_data, pattern = "^MT-")
  current_data[["percent.rb"]] <-
    Seurat::PercentageFeatureSet(current_data, pattern = "^RP[SL]")

  # 无论是 python 版本还是 R 版本的 scrublet 程序在我的电脑上都无法运行；
  # 半路导致不可逆的崩溃，在 R 中甚至能把 REPL 都崩掉 ...

  # scrublet_meta <- scrublet(seurat_obj = current_data)
  current_data[["is_doublet"]] <- rep(FALSE, ncol(current_data))

  current_data @ active.ident

  # 绘制各项 QC 指标的分布状况；一个点代表一个细胞
  qc_violin <- function(current_data) {
    Seurat::VlnPlot(
      current_data,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.rb"
      ),
      ncol = 4, pt.size = 0.1
    ) +
      ggplot2::theme(plot.title = element_text(size = 10))
  }

  qc_violin(current_data)

  # 观察各个变量之间存在什么关系？
  # 测序深度 v.s. 线粒体基因（核糖体基因）比例：
  # - 深度越大，线粒体基因越少，核糖体基因不甚变化或略减少
  Seurat::FeatureScatter(current_data, feature1 = "nCount_RNA",
                         feature2 = "percent.mt")
  Seurat::FeatureScatter(current_data, feature1 = "nCount_RNA",
                         feature2 = "percent.rb")
  # 测序深度 v.s. 检出基因数：正相关
  Seurat::FeatureScatter(current_data, feature1 = "nCount_RNA",
                         feature2 = "nFeature_RNA")
  # 线粒体基因和核糖体基因检出之间：一般线粒体基因检出率高的，核糖体检出率低（质量差）
  Seurat::FeatureScatter(current_data, feature1 = "percent.mt",
                         feature2 = "percent.rb")

  # 上面的图清楚地表明，高 MT 百分比与低 UMI 计数强烈相关，并且通常被解释为死细胞。
  # 然而，高核糖体蛋白含量与 MT 强烈反相关，并且似乎包含生物信号。 双联体分数和表达基因的
  # 数量之间也有很强的相关性。 让我们在元数据中设置 QC 列，并以信息性的方式定义它。

  # QC

  # 如果 QC 合格，最终显示为 pass.

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ is_doublet == TRUE,
    "doublet", "pass"
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 500 &
      current_data @ meta.data $ qc == "pass",
    "low_n_feature", current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 500 &
      current_data @ meta.data $ qc != "pass" &
      current_data @ meta.data $ qc != "low_n_feature",
    paste("low_n_feature", current_data @ meta.data $ qc, sep = ":"),
    current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ percent.mt > 15 &
      current_data@ meta.data $ qc == "pass", # nolint
    "high_mt", current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 500 &
      current_data @ meta.data $ qc != "pass" &
      current_data @ meta.data $ qc != "high_mt",
    paste("high_mt", current_data @ meta.data $ qc, sep = ":"),
    current_data @ meta.data $ qc
  )

  table(current_data[["qc"]])
  qc_violin(current_data |> subset(qc == "pass"))

  # 归一化

  # 为了做进一步的分析，我们需要将数据标准化以考虑测序深度。传统的方法是将其缩放到 10，000
  # （就好像所有单元格总共有 10 k 个 UMI，CPM），并对获得的值进行 log2 变换。
  # 标准化数据储存在 srat[['RNA']] @ data 中

  current_data <- Seurat::NormalizeData(current_data)

  # 下一步发现全部细胞间最可变的特征（基因）- 这些通常是下游分析最感兴趣的。
  current_data <- Seurat::FindVariableFeatures(current_data,
                                               selection.method = "vst",
                                               nfeatures = 2000)

  # 标在图上
  # geo/gse202374/var-genes.bmp # nolint

  top10 <- head(Seurat::VariableFeatures(current_data), 10)
  top10_plot <- Seurat::VariableFeaturePlot(current_data)
  top10_plot <- Seurat::LabelPoints(
    plot = top10_plot,
    points = top10,
    repel = TRUE,
    xnudge = 0, ynudge = 0,
    size = 4.0
  ) + unify_theme_font()

  # ScaleData 将标准化的基因表达转换为 z 分数（以 0 为中心的值，方差为 1）。
  # 存储在 srat[['RNA']] @ scale.data 中，用于后续 PCA；默认值是仅对可变基因运行缩放。

  all_genes <- rownames(current_data)
  current_data <- Seurat::ScaleData(current_data, features = all_genes)

  # 我们只取前两千个高变异基因做 PCA

  current_data <- RunPCA(
    current_data,
    features = Seurat::VariableFeatures(object = current_data)
  )

  # Prinicpal component “loadings” should match markers of distinct populations
  # for well behaved datasets.

  pc_loadings <- Seurat::VizDimLoadings(
    current_data, dims = 1:9, reduction = "pca"
  ) + unify_theme_font()

  ggsave(
    "pc-loadings.png", pc_loadings, dpi = 100,
    width = 12, height = 18, units = "in"
  )

  pc_heatmap <- Seurat::DimHeatmap(
    current_data, dims = 1:9, reduction = "pca",
    balanced = TRUE, cells = 50,
    fast = FALSE, combine = TRUE
  ) + unify_theme_font() + theme(axis.text.x = element_blank()) +
    scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"))

  ggsave(
    "pc-heatmap.png", dpi = 100,
    width = 12, height = 18, units = "in"
  )

  # 画出主成分分析的图，以及各 PC 的占比

  dim_plot_pca <- Seurat::DimPlot(current_data, reduction = "pca") +
    unify_theme_font()
  elbow_plot_pca <- Seurat::ElbowPlot(current_data, reduction = "pca") +
    unify_theme_font()

  # 我们现在可以进行聚类。更高的分辨率会导致更多的聚类（默认值为 0.8）。这将是非常重要的，
  # 以找到正确的集群分辨率，因为细胞类型的标志物取决于集群的定义。
  #
  # 我们可能需要更换 resolution 做许多实验

  current_data <- Seurat::FindNeighbors(current_data, dims = 1:10)
  current_data <- Seurat::FindClusters(current_data, resolution = 0.5)
  table(current_data @ meta.data $ seurat_clusters)

  # 出于可视化的目的，我们还需要生成 UMAP 降维表示：

  current_data <- Seurat::RunUMAP(current_data, dims = 1:10, verbose = FALSE)

  # DimPlot 默认使用 UMAP，Seurat 集群作为标识：

  dim_plot_umap <- Seurat::DimPlot(
    current_data,
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  # 我们可以在 UMAP 图上表示出某个基因表达的亚群位置

  gene_names <- current_data @ assays $ RNA @ features |> rownames()
  ils <- gene_names[grep("^IL[0-9]", gene_names)]
  expr_maps_ils <- Seurat::FeaturePlot(
    current_data, features = ils,
    combine = FALSE
  )

  cds <- gene_names[grep("^CD[0-9]", gene_names)]
  expr_maps_cds <- Seurat::FeaturePlot(
    current_data, features = cds,
    combine = FALSE
  )

  patch <- expr_maps_cds[[1]] + unify_theme_font()
  for (x in seq_along(cds)) {
    if (x == 1) next
    graph <- expr_maps_cds[[x]] + unify_theme_font()
    patch <- patch + graph
  }

  ggsave(
    "expression-umap-cd.png", dpi = 300, plot = patch,
    width = 40, height = 18, units = "in"
  )

  # 让我们删除未通过 QC 的细胞并比较图。我们现在可以看到更多定义的集群。

  dim_plot_umap_qc <- Seurat::DimPlot(
    current_data |> subset(qc == "pass"),
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  ggsave(
    "umap.png", dpi = 300, plot = dim_plot_umap,
    width = 5, height = 5, units = "in"
  )

  ggsave(
    "umap-qc.png", dpi = 300, plot = dim_plot_umap_qc,
    width = 5, height = 5, units = "in"
  )

  # 细胞周期修正

  cc_genes <- Seurat::cc.genes
  cc_genes_2019 <- Seurat::cc.genes.updated.2019

  current_data <- Seurat::CellCycleScoring(
    current_data,
    s.features = cc_genes_2019 $ s.genes,
    g2m.features = cc_genes_2019 $ g2m.genes
  )

  table(current_data[[]] $ Phase)

  # 最后，我们看一下 QC 指标在 UMAP 图上的分布

  mt_feature <- Seurat::FeaturePlot(
    current_data, features = "percent.mt",
    label.size = 4, repel = TRUE, label = TRUE
  ) + unify_theme_font()

  rb_feature <- Seurat::FeaturePlot(
    current_data, features = "percent.rb",
    label.size = 4, repel = TRUE, label = TRUE
  ) + unify_theme_font()

  ggsave(
    "umap-mt.png", dpi = 300, plot = mt_feature,
    width = 5, height = 5, units = "in"
  )

  ggsave(
    "umap-rb.png", dpi = 300, plot = rb_feature,
    width = 5, height = 5, units = "in"
  )

  # 观察一个指标在不同细胞分群中的分布状况

  mt_violin <- Seurat::VlnPlot(
    current_data, features = "percent.mt"
  ) + unify_theme_font()

  rb_violin <- Seurat::VlnPlot(
    current_data, features = "percent.rb"
  ) + unify_theme_font()

  depth_violin <- Seurat::VlnPlot(
    current_data, features = "nCount_RNA"
  ) + unify_theme_font()

  feature_violin <- Seurat::VlnPlot(
    current_data, features = "nFeature_RNA"
  ) + unify_theme_font()

  ggsave(
    "violin-mt.png", dpi = 300, plot = mt_violin,
    width = 6, height = 4, units = "in"
  )

  ggsave(
    "violin-rb.png", dpi = 300, plot = rb_violin,
    width = 6, height = 4, units = "in"
  )

  ggsave(
    "violin-depth.png", dpi = 300, plot = depth_violin,
    width = 6, height = 4, units = "in"
  )

  ggsave(
    "violin-feature.png", dpi = 300, plot = feature_violin,
    width = 6, height = 4, units = "in"
  )

  # 由于我们已经进行了大量的 QC，去除了双联体和空细胞，我们现在可以应用 SCTransform
  # 归一化，这被证明有利于通过提高信噪比来发现稀有细胞群。单个 SCTransform 命令替换
  # NormalizeData、ScaleData 和 FindVariableFeatures。
  #
  # 我们还将使用 vars.to.regress 变量校正 %mt 基因和细胞周期分数；我们先前的探索
  # 已经表明，细胞周期分数和 %mt 在簇之间都没有非常显著的变化，因此我们不会去除生物信号，
  # 而只是一些不需要的变化。

  current_data <- Seurat::SCTransform(
    current_data,
    method = "glmGamPoi",
    ncells = sum(current_data @ meta.data $ qc == "pass"),
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    verbose = FALSE
  )

  current_data <- Seurat::RunPCA(current_data, verbose = FALSE)
  current_data <- Seurat::RunUMAP(current_data, dims = 1:30, verbose = FALSE)
  current_data <- Seurat::FindNeighbors(current_data, dims = 1:30,
                                        verbose = FALSE)
  current_data <- Seurat::FindClusters(current_data, verbose = FALSE)
  table(current_data[[]] $ seurat_clusters)

  normalized_umap <- Seurat::DimPlot(
    current_data,
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  normalized_umap_qc <- Seurat::DimPlot(
    current_data |> subset(qc == "pass"),
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  ggsave(
    "umap-norm.png", dpi = 300, plot = normalized_umap,
    width = 5, height = 5, units = "in"
  )

  ggsave(
    "umap-qc-norm.png", dpi = 300, plot = normalized_umap_qc,
    width = 5, height = 5, units = "in"
  )

  # 差异表达和标记选择

  # 差异表达使我们能够定义每个簇的特异性基因标记。根据定义，它受到聚类定义方式的影响，
  # 因此在定义标记之前找到聚类的正确分辨率非常重要。如果某些聚类缺少任何显著的标记，
  # 请调整聚类。建议在 RNA 检测试剂盒上进行差异表达，而不是 SCTransform 后的数据集。
  # 差异表达可以在两个特定簇之间以及簇和所有其他细胞之间进行。

  Seurat::DefaultAssay(current_data) <- "RNA"
  current_data <- Seurat::NormalizeData(current_data)
  current_data <- Seurat::FindVariableFeatures(
    current_data, selection.method = "vst",
    nfeatures = 2000
  )

  all_genes <- rownames(current_data)
  current_data <- Seurat::ScaleData(current_data, features = all_genes)

  # 下面的函数允许通过将每个聚类与所有剩余的细胞进行比较来找到每个聚类的标记，
  # 同时仅报告阳性细胞。有许多测试可以用来定义标记，包括非常快速和直观的 tf-idf。

  # 默认情况下，使用 Wilcoxon 秩和检验。这需要一段时间！为了提高速度，我们增加了默认的
  # 最小百分比和 log 2FC 截止值；这些应该进行调整以适应数据集！

  all_markers <- Seurat::FindAllMarkers(
    current_data, only.pos = TRUE,
    min.pct = 0.5, logfc.threshold = 0.5
  )

  dim(all_markers)
  table(all_markers $ cluster)
  top3_markers <- as.data.frame(
    all_markers |> group_by(cluster) |>
      top_n(n = 3, wt = avg_log2FC)
  )

  #            p_val avg_log2FC pct.1 pct.2     p_val_adj cluster    gene
  # 1   0.000000e+00   3.261608 0.884 0.220  0.000000e+00       0    GZMK (1)
  # 2   0.000000e+00   2.777491 0.729 0.130  0.000000e+00       0    CD8A
  # 3   0.000000e+00   2.920525 0.624 0.098  0.000000e+00       0    CD8B
  # 4  2.172191e-226   1.787039 0.580 0.161 4.374793e-222       1  CD40LG (2)
  # 5  1.905731e-176   1.687896 0.885 0.546 3.838142e-172       1     LTB
  # 6  2.470413e-165   1.543393 0.656 0.275 4.975411e-161       1    CAPG
  # 7   0.000000e+00   1.126051 1.000 0.859  0.000000e+00       2   RPS12 (3)
  # 8  1.501876e-137   1.156144 0.862 0.495 3.024779e-133       2    LDHB
  # 9   3.684642e-76   1.237884 0.546 0.267  7.420870e-72       2   NOSIP
  # 10  0.000000e+00   1.363160 1.000 0.998  0.000000e+00       3  MALAT1 (4)
  # 11  8.052117e-99   1.592294 0.753 0.757  1.621696e-94       3    XIST
  # 12  7.283491e-50   1.495279 0.606 0.602  1.466895e-45       3    NKTR
  # 13 1.150099e-255   1.569820 0.979 0.556 2.316299e-251       4    IL7R (5)
  # 14 2.097075e-161   1.142827 0.980 0.864 4.223508e-157       4   TXNIP
  # 15  4.461819e-70   1.129163 0.545 0.263  8.986103e-66       4  GPR171
  # 16  0.000000e+00   4.118931 0.753 0.092  0.000000e+00       5    GNLY (6)
  # 17  0.000000e+00   4.321742 0.686 0.045  0.000000e+00       5    TRDC
  # 18  0.000000e+00   7.347585 0.512 0.005  0.000000e+00       5   KRT81
  # 19 3.791806e-268   3.483729 0.912 0.188 7.636697e-264       6   KLRD1 (7)
  # 20  2.708023e-65   2.613728 0.558 0.224  5.453958e-61       6    ZEB2
  # 21  1.203032e-58   2.881647 0.542 0.249  2.422907e-54       6   DIP2A
  # 22  0.000000e+00   5.227008 0.878 0.054  0.000000e+00       7    GZMB (8)
  # 23  0.000000e+00   7.653717 0.779 0.008  0.000000e+00       7  FGFBP2
  # 24  0.000000e+00   4.563473 0.775 0.039  0.000000e+00       7   KLRF1
  # 25 5.456180e-132   2.490954 1.000 0.970 1.098875e-127       8  MT-CO3 (9)
  # 26  1.317216e-82   3.180156 0.662 0.225  2.652872e-78       8    ZEB2
  # 27  1.578227e-28   2.762834 0.585 0.390  3.178549e-24       8    SAT1
  # 28  0.000000e+00  11.694427 1.000 0.003  0.000000e+00       9   TPSB2 (10)
  # 29  0.000000e+00  12.129890 1.000 0.003  0.000000e+00       9  TPSAB1
  # 30  0.000000e+00  11.816134 0.991 0.002  0.000000e+00       9    CPA3
  # 31  0.000000e+00   8.075968 0.604 0.005  0.000000e+00      10  LGALS2 (11)
  # 32  0.000000e+00   7.855680 0.595 0.008  0.000000e+00      10    CD14
  # 33  0.000000e+00   8.339411 0.515 0.002  0.000000e+00      10 CLEC10A
  # 34  0.000000e+00   6.301742 0.606 0.010  0.000000e+00      11   FOXP3 (12)
  # 35 1.577843e-304   4.826458 0.571 0.031 3.177776e-300      11   CTLA4
  # 36 1.718181e-107   4.065246 0.581 0.105 3.460417e-103      11 TNFRSF4
  # 37  0.000000e+00   6.981980 0.995 0.017  0.000000e+00      12  SMIM25 (13)
  # 38  0.000000e+00   6.874433 0.913 0.026  0.000000e+00      12  CDKN1C
  # 39  0.000000e+00   9.120403 0.579 0.004  0.000000e+00      12   LYPD2
  # 40 9.119791e-177   6.413925 0.538 0.047 1.836726e-172      13  CCL4L2 (14)
  # 41 1.986613e-174   5.723539 1.000 0.268 4.001038e-170      13    CCL4
  # 42  2.070554e-31   2.814468 0.532 0.206  4.170096e-27      13    XCL2
  # 43  0.000000e+00   8.524827 0.754 0.007  0.000000e+00      14   CD79A (15)
  # 44  0.000000e+00   9.418402 0.626 0.004  0.000000e+00      14    IGHM
  # 45  0.000000e+00  13.110377 0.544 0.004  0.000000e+00      14    IGKC
  # 46  0.000000e+00  11.844500 0.984 0.006  0.000000e+00      15    RGS5 (16)
  # 47  0.000000e+00  12.059344 0.762 0.001  0.000000e+00      15   KCNE4
  # 48  0.000000e+00  10.506028 0.648 0.000  0.000000e+00      15   EDNRA
  # 49 4.634498e-215   5.973344 0.550 0.025 9.333879e-211      16   TIMP3 (17)
  # 50 8.198994e-123   5.037129 0.550 0.049 1.651277e-118      16  IGFBP4
  # 51 1.895203e-111   6.164919 0.595 0.066 3.816939e-107      16     GSN
  # 52  0.000000e+00  10.274761 0.871 0.001  0.000000e+00      17  S100A1 (18)
  # 53  0.000000e+00  10.541747 0.613 0.001  0.000000e+00      17   NPTX2
  # 54  0.000000e+00  10.179475 0.505 0.000  0.000000e+00      17  ACSM2B
  # 55 9.369878e-176   8.010029 0.900 0.067 1.887093e-171      18    MT1E (19)
  # 56  7.601643e-95   7.084123 0.971 0.175  1.530971e-90      18    MT1X
  # 57  8.241566e-59   6.255310 1.000 0.404  1.659851e-54      18    MT2A

  # 如果我们知道一些先验知识，可以先大致看出某些细胞分群

  draw_single_gene_expr <- function(gene, description) {
    Seurat::FeaturePlot(current_data, gene, label = TRUE) +
      scale_colour_gradientn(
        colours = rev(brewer.pal(n = 11, name = "Spectral"))
      ) +
      ggtitle(paste0(gene, ": ", description)) +
      unify_theme_font()
  }

  draw_single_gene_expr("PPBP", "Platlets") # no such gene
  draw_single_gene_expr("LILRA4", "Dendritic cells")
  draw_single_gene_expr("MS4A1", "B cells")
  draw_single_gene_expr("LYZ", "Monocytes")
  draw_single_gene_expr("NKG7", "NK cells")
  draw_single_gene_expr("CD8B", "CD8+ T cells")
  draw_single_gene_expr("IL7R", "CD4+ T cells")

  saveRDS(current_data, "normalized-data.rds")
}