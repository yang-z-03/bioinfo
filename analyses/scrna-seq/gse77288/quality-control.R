
library(GEOquery)

library(tidyverse)
library(magrittr)
library(dplyr)
library(data.table)
library(scales)

install_packages(c("SingleCellExperiment"))
library(SingleCellExperiment)

install_packages(c("scater"))
library(scater)

install_packages(c("scran"))
library(scran)

## Download GEO database ------------------------------------------------------

# 从原始数据（DNA reads）分析单细胞测序结果
series_id <- "GSE77288"
download_table <- geo_get_supplemental_links(series_id)
geo_get_supplemental_files(series_id)

expr_set <- geo_get_series_annot(series_id)
platform <- geo_get_platform_annot(expr_set)

file_name <- "GSE77288_molecules-raw-single-per-sample.txt.gz"
root_dir <- paste("geo", series_id, sep = "/")

gunzip(paste("geo", series_id, file_name, sep = "/"),
       destname = paste("geo", series_id, "txt", sep = "/"))

## Read database content ------------------------------------------------------

molecules_raw <- read.table(paste("geo", series_id, "txt", sep = "/"),
                            sep = "\t", header = TRUE) |> t() # transpose it.

# a data table with 864 rows and 20422 columns.
# 每行表示一个细胞，每列表示一个检测基因，
# 但是我们更希望每行表示一个基因，因为我们往往需要根据基因筛选行

# # A tibble: 864 × 20,422
#    individual replicate well  ENSG00000186092 ENSG00000237683 ENSG00000235249
#    <chr>      <chr>     <chr>           <int>           <int>           <int>
#  1 NA19098    r1        A01                 0               0               0
#  2 NA19098    r1        A02                 0               0               0
#  3 NA19098    r1        A03                 0               0               0
#  4 NA19098    r1        A04                 0               1               0
#  5 NA19098    r1        A05                 0               0               0
#  6 NA19098    r1        A06                 0               0               0
#  7 NA19098    r1        A07                 0               0               0
#  8 NA19098    r1        A08                 0               0               0
#  9 NA19098    r1        A09                 0               0               0
# 10 NA19098    r1        A10                 0               0               0

# 现在我们已经获得了数据，构造行名称
# 注意， df[1,] 选择第一行， df[,1] 选择第一列； 而 df[1] 选择第一列返回一个 data.frame

dataset <- molecules_raw[4:nrow(molecules_raw), ] |> as_tibble()
sample_names <- paste(molecules_raw[1, ],
                      molecules_raw[2, ], molecules_raw[3, ], sep = "-")

# note that the data passed to the assay slot has to be a matrix!
tung <- SingleCellExperiment(

  # 这个 SCE 中只有一个实验，实验名称叫 counts
  assays = list(counts = as.matrix(dataset) |> apply(2, as.numeric))
)

colnames(tung) <- sample_names
colData(tung) $ individual <- molecules_raw[1, ]
colData(tung) $ replicate <- molecules_raw[2, ]
colData(tung) $ well <- molecules_raw[3, ]
colData(tung) $ batch <-
  paste(molecules_raw[1, ], molecules_raw[2, ], sep = "-")
colData(tung) $ id <- sample_names
rownames(tung) <- rownames(molecules_raw[4:nrow(molecules_raw), ])

rm(dataset, sample_names, molecules_raw)

## Statistical analysis -------------------------------------------------------

# 基本的描述统计学处理

# 计算每个基因在不同细胞中的表达均数和标准差
colData(tung) $ mean_counts <- colMeans(counts(tung))
colData(tung) $ median_counts <- colMedians(counts(tung))
colData(tung) $ sd_counts <- colSds(counts(tung))
colData(tung) $ total_counts <- colSums(counts(tung))

# 计算 CPM 矩阵， sweep 函数 1 为按行， 2 为按列
assay(tung, "cpm") <- sweep(counts(tung), 2, tung $ total_counts / 1e6,
                            function(x, y) {
                              return(x / y)
                            })

# 我们可以使用上述统计量筛选基因行
by_row <- tung[rowMeans(counts(tung)) > 0.01, ]
by_col <- tung[, colData(tung) $ mean_counts > 0.01]

## Descriptive visualizations -------------------------------------------------

cell_info <- as.data.frame(colData(tung))

# 每一批次样本检测到的总基因数（批次效应）
ggplot(data = cell_info, aes(x = batch, y = total_counts)) +
  geom_violin(fill = "brown") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Quality control -------------------------------------------------------------

# 此前步骤中获得的表达矩阵和注释信息
# tung

# 从主对象中删除 ERCC
altExp(tung, "ercc") <-
  tung[grep("^ercc.", rownames(tung), ignore.case = TRUE), ]

tung <- tung[grep("^ercc.", rownames(tung),
                  invert = TRUE, ignore.case = TRUE), ]
rowData(tung) $ src_id <- rownames(tung)
ensembl_genes <- readRDS("ensembl/9606/genes.rds")

appending <- ensembl_genes[, c("seqname", "feature", "start", "end", "score",
                               "strand", "frame", "name", "biotype",
                               "description", "version", "gene_id")] |>
  `colnames<-`(c("ens_seqname", "ens_feature", "ens_start", "ens_end",
                 "ens_score", "ens_strand", "ens_frame", "ens_name",
                 "ens_biotype", "ens_desc", "ens_v", "gene_id")) |>
  data.table()

# If all.x is true, all the non matching cases of x are appended to the result
# as well, with NA filled in the corresponding columns of y.
rowData(tung) <- merge(x = rowData(tung), sort = FALSE,
                       y = appending,
                       all.x = TRUE, by.x = "src_id", by.y = "gene_id")

# 删除没有检索到基因名的行
tung <- tung[!is.na(rowData(tung) $ ens_name), ]

# 线粒体质控：下面的 scater 函数允许我们添加对数据集评估有用的每细胞和每基因指标。 单细胞
# 实验最常见的度量是计数总数（UMI）、检测到的基因总数、线粒体计数总数、线粒体计数百分比等。
mitochondrial_genes <- appending[ens_seqname == "MT", ]

is_mito <- rowData(tung) $ ens_name %in% mitochondrial_genes $ ens_name
is_ribosomal <- rowData(tung) $ ens_name %in%
  grep("^RP[LS]", rowData(tung) $ ens_name, value = TRUE)

# 使用正则表达式可以简单的查看线粒体和核糖体基因
grep("^RP[LS]", rowData(tung) $ ens_name, value = TRUE)
grep("^MT-", rowData(tung) $ ens_name, value = TRUE)

#  [1] "sum"                        "detected"
#  [3] "subsets_mito_sum"           "subsets_mito_detected"
#  [5] "subsets_mito_percent"       "subsets_ribosomal_sum"
#  [7] "subsets_ribosomal_detected" "subsets_ribosomal_percent"
#  [9] "altexps_ercc_sum"           "altexps_ercc_detected"
# [11] "altexps_ercc_percent"       "total"

subset_metrics <- perCellQCMetrics(tung,
                                   subsets = list(mito = is_mito,
                                                  ribosomal = is_ribosomal))

# [1] mean      [2] detected
feature_metrics <- perFeatureQCMetrics(tung)

# 这样可以直接将这些 QC 加到元数据中
tung <- addPerCellQC(tung, subsets = list(mito = is_mito,
                                          ribosomal = is_ribosomal))
tung <- addPerFeatureQC(tung)

# 分布图：一个细胞有多少个分子被检测到了
draw_total_expression_per_cell <- function(data) {
  hist(
    data $ total,
    breaks = 100
  )
  abline(v = 25000, col = "red")
}

# 分布图：一个细胞有多少个不同的基因被检测到了
draw_total_detected_per_cell <- function(data) {
  hist(
    data $ detected,
    breaks = 100
  )
  abline(v = 7000, col = "red")
}

draw_total_expression_per_cell(tung)
draw_total_detected_per_cell(tung)

# 有时很难给出一个明显的过滤截止值。  在这种情况下，自适应阈值可以帮助我们识别与我们用于 QC
# 的任何变量中的中位数相差超过3个中位数绝对偏差（MAD）的点。  注意指定偏差的正确方向：
# 检测到的基因数量少，但 MT 基因百分比高，是低质量细胞的标志

qc_low_sum <-
  isOutlier(subset_metrics $ sum, log = TRUE, type = "lower")
qc_low_detect <-
  isOutlier(subset_metrics $ detected, log = TRUE, type = "lower")

# ERCC 评估
#
# ERCC（外部RNA控制联盟）开发了一套 RNA 标准品，用于微阵列、qPCR和测序应用中的质量控制。
# 这些RNA标准品是已知浓度和组成（即序列长度和GC含量）的加标 RNA。它们可用于评估 RNA-seq
# 数据的灵敏度和准确性。
#
# the altexps_ercc_percent column contains the percentage of reads mapped
# to ERCC transcripts.
#
# ERCC 相当于一种用于定量的标准曲线。在每个样本中加入固定的剂量。  当其百分比很高时，我们应该
# 怀疑样本的质量不好。

qc_ercc <- isOutlier(subset_metrics $ altexps_ercc_percent, type = "higher")
qc_mito <- isOutlier(subset_metrics $ subsets_mito_percent, type = "higher")

discard_crit <- qc_low_detect | qc_low_sum | qc_ercc | qc_mito

# 手动显示计算的离群值数量：
DataFrame(low_lib_size = sum(qc_low_sum),
          low_n_features = sum(qc_low_detect),
          high_ercc_percent = sum(qc_ercc),
          high_mito_percent = sum(qc_mito),
          discard = sum(discard_crit))

reasons <- quickPerCellQC(subset_metrics,
                          sub.fields = c("subsets_mito_percent",
                                         "altexps_ercc_percent"))
# 自动计算指定的两列离群值舍掉的样本个数
colSums(as.matrix(reasons))

# 自动的或手动的将筛选条件塞进元数据
# tung $ discard <- discard_crit # nolint
tung $ discard <- reasons $ discard

# 可视化两群细胞的筛选量
# 参数：sum, detected, mito-percent, ercc-percent
plotColData(tung, x = "sum", y = "detected", colour_by = "discard")
plotColData(tung, x = "sum", y = "subsets_mito_percent", colour_by = "discard")
plotColData(tung, x = "altexps_ercc_percent",
            y = "subsets_mito_percent", colour_by = "discard")

# 我们可以分批次的观察这些筛选量，来发现潜在的批次效应
plotColData(tung, x = "sum", y = "detected", colour_by = "discard",
            other_fields = "replicate") +
  facet_wrap(~replicate) +
  scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))

# 至此，我们完成了样本的筛选

# Expression visualizations ---------------------------------------------------

# 整个数据集中表达最多的基因
# 我们看到的大多数基因都是线粒体或核糖体蛋白，这对于大多数 scRNA-seq 数据集来说是非常典型的

expr_rank <- plotHighestExprs(tung, exprs_values = "counts",
                              feature_names_to_plot = "ens_name",
                              colour_cells_by = "detected")

ggsave(paste0(root_dir, "/", "expression-rank.png"), expr_rank,
       width = 4, height = 6, units = "in", dpi = 600)

keep_feature <- nexprs(tung, byrow = TRUE, detection_limit = 1) >= 2
expr_crit <- !keep_feature
table(expr_crit)

# 我们完成基因行的筛选
rowData(tung) $ discard_row <- expr_crit

assay(tung, "logcounts_raw") <- log2(counts(tung) + 1)
saveRDS(tung, file = paste0(root_dir, "/", "dataset.rds"))

# Getting perspects on what QC has done ---------------------------------------

tung <- readRDS(paste0(root_dir, "/", "dataset.rds"))
tung_qc <- tung[!rowData(tung) $ discard_row, !colData(tung) $ discard]

# 在清洗前后的数据上运行 PCA
tung <- runPCA(tung, exprs_values = "counts")
plotPCA(tung, colour_by = "batch",
        size_by = "detected", shape_by = "individual")

tung <- runPCA(tung, exprs_values = "logcounts_raw")
plotPCA(tung, colour_by = "batch",
        size_by = "detected", shape_by = "individual")

tung_qc <- runPCA(tung_qc, exprs_values = "logcounts_raw")
plotPCA(tung_qc, colour_by = "batch",
        size_by = "detected", shape_by = "individual")

# 默认时，我们只选了最大的 500 个基因做主成分分析的数据，以获得一个速度的妥协；我们可以
# 指定使用更少或者更多的基因。

tung_qc <- runPCA(tung_qc, exprs_values = "logcounts_raw", ntop = 50)
plotPCA(tung_qc, colour_by = "batch",
        size_by = "detected", shape_by = "individual")

tung_qc <- runPCA(tung_qc, exprs_values = "logcounts_raw", ntop = nrow(tung_qc))
plotPCA(tung_qc, colour_by = "batch",
        size_by = "detected", shape_by = "individual")

# t-SNE

set.seed(42)
tung <- runTSNE(tung, exprs_values = "logcounts_raw", perplexity = 130)
plotTSNE(tung, colour_by = "batch", size_by = "detected",
         shape_by = "individual")

set.seed(42)
tung_qc <- runTSNE(tung_qc, exprs_values = "logcounts_raw", perplexity = 130)
plotTSNE(tung_qc, colour_by = "batch", size_by = "detected",
         shape_by = "individual")

# scRNA-seq数据中存在大量潜在的混杂因素、伪影和偏倚。分析 scRNA-seq 数据的主要挑战之一
# 是很难进行真正的技术复制（为什么？）来区分生物学和技术上的差异。在前面的章节中，我们考虑了
# 批量效应，在本章中，我们将继续探索如何识别和删除实验性伪影

# scater 允许识别与感兴趣的实验和QC变量相关的主成分（它根据线性模型回归 PC 值对感兴趣的
# 变量的主成分进行 R^2 排名）

plotExplanatoryVariables(tung_qc, exprs_values = "logcounts_raw",
                         variables = c("detected", "sum", "batch",
                                       "individual", "altexps_ercc_percent",
                                       "subsets_mito_percent"))

# 我们可以发现，sum 变量（测序深度）、detected 变量（检测到的基因数）解释力很强，这是单细胞
# 测序最常见的问题。
# 此外，batch R2 > individual，这是典型的批次效应表现。

# 除了对批次进行校正之外，还有其他因素可能需要补偿。  与批量校正一样，这些调整需要外部信息。
# 一种流行的方法是 scLVM，它允许您识别并减去细胞周期或凋亡等过程的影响。
# 此外，不同的方案可能在对每个转录物的覆盖范围、它们基于 A/T 核苷酸的平均含量产生的偏倚或
# 它们捕获短转录物的能力方面有所不同。理想情况下，我们希望补偿所有这些差异和偏倚。

# Normalization ---------------------------------------------------------------

## 原始 logrank ----------------------------------------------------------------

tung_qc <- runPCA(tung_qc, exprs_values = "logcounts_raw")
plotPCA(tung_qc, colour_by = "batch", size_by = "detected",
        shape_by = "individual")

## 使用 CPM 矫正 ---------------------------------------------------------------

logcounts(tung_qc) <- log2(calculateCPM(tung_qc) + 1)
tung_qc <- runPCA(tung_qc)
plotPCA(tung_qc, colour_by = "batch", size_by = "detected",
        shape_by = "individual")

# a relative log expression (RLE) plots can be very useful assessing whether
# normalization procedure was successful.

plotRLE(tung_qc, exprs_values = "logcounts_raw", colour_by = "batch") +
  ggtitle("RLE plot for logcounts_raw")

plotRLE(tung_qc, exprs_values = "logcounts", colour_by = "batch") +
  ggtitle("RLE plot for log2(CPM) counts")

## size factor 矫正 ------------------------------------------------------------

# 基于 CPM 和其他类似的文库缩放方法假设所有细胞含有相似量的 RNA，因此应该产生相似的
# UMI 计数。  这并不总是正确的。  下面的方法，在 scran 和其他几个 bioconductoR 包中可用，
# 使用聚类来进行规范化，这有时被称为通过去卷积的归一化。  这些集群看起来很像我们的批次

qclust <- quickCluster(tung_qc, min.size = 30)

# 函数向 colData 添加一个名为 sizeFactor 的列。这些值然后由 logNormCounts 使用

tung_qc <- computeSumFactors(tung_qc, clusters = qclust)

# 有时 scran 产生负的或零的尺寸因子。这将完全扭曲规范化的表达式矩阵。
# 我们可以这样检查scran计算出的尺寸因子：如果您发现 scran 计算的大小因子为负，请尝试
# 增加群集和池大小，直到它们都为正。

summary(sizeFactors(tung_qc))

#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.3938  0.7614  0.9551  1.0000  1.1566  3.2914

tung_qc <- logNormCounts(tung_qc)

# 一样的方法检测效果
tung_qc <- runPCA(tung_qc)
plotPCA(tung_qc, colour_by = "batch", size_by = "detected",
        shape_by = "individual")
plotRLE(tung_qc, exprs_values = "logcounts", colour_by = "batch")

## downsample 矫正 -------------------------------------------------------------

logcounts(tung_qc) <- log2(downsampleMatrix(counts(tung_qc), prop = 0.5) + 1)

# 一样的方法检测效果
tung_qc <- runPCA(tung_qc)
plotPCA(tung_qc, colour_by = "batch", size_by = "detected",
        shape_by = "individual")
plotRLE(tung_qc, exprs_values = "logcounts", colour_by = "batch")

# Adjusting correlations for other parameters (Batch effect) ------------------

# 我们对库大小进行了标准化，有效地消除了它作为混淆因素的影响。现在我们将考虑从我们的数据中
# 删除其他定义不太明确的混杂因素。技术混杂因素（又称批次效应）可能源于试剂、分离方法、
# 进行实验的实验室、实验员，甚至是进行实验的日期或时间的差异。
#
# 解释技术混杂因素，特别是批量效应，是一个涉及实验设计原则的大课题。  在这里，我们解决的方法，
# 可以考虑到混杂因素时，必须假定实验设计是适当的
#
# 从根本上讲，解释技术混杂因素涉及识别并理想地去除表达数据中与感兴趣的生物信号不相关（即混淆）
# 的变异来源。有多种方法，其中一些使用加标、或管家基因，并且其中一些使用内源基因

# 使用 spike-in 作为对照基因是有吸引力的，因为在我们的实验中，将相同量的 ERCC（或其他）
# spike-in 添加到每个细胞中。原则上，我们观察到的这些基因的所有变异性都是由于技术噪音造成的;
# 而内源性基因则受到技术噪音和生物变异性的影响。技术噪声可以通过将模型拟合到加标物并从
# 内源基因中“减去”来去除。
#
# 有几种方法可以基于这个前提（例如 BASiCS、scLVM、RUVg）; 每一个都使用不同的噪声模型和
# 不同的拟合过程。  或者，人们可以识别出表现出超出技术噪音的显著变异的基因
# （例如，与中位数的距离，高度可变的基因）
#
# However, there are issues with the use of spike-ins for normalisation
# (particularly ERCCs, derived from bacterial sequences), including that their
# variability can, for various reasons, actually be higher than that of
# endogenous genes.

# scRNA-seq 数据集的集成中有两种情况。
#
# 在第一种情况下，细胞组成预计是相同的，这时为批量 RNA-seq 开发的方法（例如 ComBat）
# 表现出良好的性能。对于同一实验的生物重复，这通常是正确的；对于 tung 数据集中的批次也是如此。
#
# 在第二种情况下，数据集之间的重叠是部分的 —— 例如，如果数据集表示健康和患病组织，
# 其在细胞类型组成上基本上不同。在这种情况下，基于互最近邻（MNN）的方法往往表现得更好。

# 前述规范化：消除掉测序深度的混杂

umi    <- readRDS(paste0(root_dir, "/", "dataset.rds"))
umi_qc <- umi[!rowData(umi) $ discard_row, !colData(umi)$discard]
qclust <- quickCluster(umi_qc, min.size = 30)
umi_qc <- computeSumFactors(umi_qc, clusters = qclust)
umi_qc <- logNormCounts(umi_qc)

install_packages(c("batchelor", "sva"))
library(batchelor)
library(sva)
library(kBET)
set.seed(42)

# 以 replicate 分组矫正批次效应
assay(umi_qc, "combat") <- sva::ComBat(logcounts(umi_qc),
                                       batch = umi_qc $ replicate)

# 以 replicate 分组矫正批次效应
assay(umi_qc, "combat_detects") <- sva::ComBat(logcounts(umi_qc),
                                               batch = umi_qc $ detected)

# 使用 mnn 矫正
mnn_out <- batchelor::fastMNN(umi_qc, batch = umi_qc $ replicate)
assay(umi_qc, "mnn") <- assay(mnn_out, "reconstructed")

# 最后，我们显示以上所有处理的效果
for (n in assayNames(umi_qc)) {
  tmp <- runPCA(umi_qc, exprs_values = n, ncomponents = 20)

  print(
    plotPCA(
      tmp,
      colour_by = "batch",
      size_by = "detected",
      shape_by = "individual"
    ) +
      ggtitle(n)
  )

  print(
    plotRLE(umi_qc, exprs_values = n, colour_by = "batch") +
      ggtitle(n)
  )
}

# 我们还可以使用全批次的相对对数表达式（RLE）来检查校正的有效性
res <- list()
rle_all <- function(mat) {
  return(mat - median(mat))
}
for (n in assayNames(umi_qc)) {
  res[[n]] <- suppressWarnings(rle_all(assay(umi_qc, n)))
}
par(mar = c(6, 4, 1, 1))
boxplot(res, las = 2)

# 使用 kBET 包确认批次之间分布的差异性水平
compare_kbet_results <- function(sce) {
  indiv <- unique(as.character(sce $ individual))
  norms <- assayNames(sce) # Get all normalizations
  results <- list()

  for (i in indiv){ 
    for (j in norms){
      tmp <- kBET(df = assay(sce[, sce $ individual == i], j) |> t(),
                  batch = sce$batch[sce $ individual == i],
                  heuristic = TRUE,
                  verbose = FALSE,
                  addTest = FALSE,
                  plot = FALSE)
      results[[i]][[j]] <- tmp $ summary $ kBET.observed[1]
    }
  }

  return(do.call(rbind.data.frame, results))
}

eff_debatching <- compare_kbet_results(umi_qc)
eff_debatching

saveRDS(umi_qc, file = paste0(root_dir, "/", "dataset-qc.rds"))
