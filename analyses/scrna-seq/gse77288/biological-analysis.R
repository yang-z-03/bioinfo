
install_packages(c("pcaMethods", "SC3", "scater", "SingleCellExperiment",
                   "pheatmap", "mclust"))
library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
library(ggplot2)

set.seed(42)

series_id <- "GSE77288"
root_dir <- paste("geo", series_id, sep = "/")
data <- readRDS(paste0(root_dir, "/", "dataset-qc.rds"))
data <- data[!rowData(data)["ens_name"] |> duplicated()]

data |> rowData() |> colnames()
data |> colData() |> colnames()
data |> metadata() |> names()

# 最简单的分群，单纯 PCA
data <- runPCA(data, exprs_values = "combat")
plotPCA(data, color_by = "individual")

# SC3 -------------------------------------------------------------------------

data <- sc3_estimate_k(data)
k_est <- metadata(data) $ sc3 $ k_estimation

rowData(data) $ feature_symbol <- rowData(data) $ ens_name
data <- sc3(data, ks = k_est, biology = TRUE, n_cores = 23)

# for more details of sc3, see
# ://www.bioconductor.org/packages/devel/bioc/vignettes/SC3/inst/doc/SC3.html

save_to_png <- function(x, name, width, height) {
  ggsave(paste0(root_dir, "/", name, ".png"), x,
         width = width, height = height, units = "in", dpi = 600)
}

# 共识矩阵是一个 N × N 矩阵，其中 N 是输入数据集中的细胞数。它表示基于来自聚类参数的所有
# 组合的聚类结果的平均值的细胞之间的相似性。相似度 0（蓝色）意味着两个细胞总是被分配到不同
# 的簇。相比之下，相似性 1（红色）意味着两个细胞总是被分配到同一个簇。共识矩阵采用层次聚类
# 的方法进行聚类，具有对角块结构。直观地说，当所有对角块都是完全红色，所有非对角元素都是
# 完全蓝色时，就实现了完美的聚类。

x <- sc3_plot_consensus(data, k = k_est)
save_to_png(x, "sc3-consensus", 4, 4)

# 轮廓是一致性矩阵对角性的定量度量。平均轮廓宽度（在轮廓图的左下角示出）从 0 到 1 变化，
# 其中 1 表示完美的块对角一致矩阵，并且 0 表示不存在块对角结构的情况。当平均轮廓宽度
# 接近 1 时，实现最佳聚类。
x <- sc3_plot_silhouette(data, k = k_est)
save_to_png(x, "sc3-silhouette", 4, 4)

# 表达面板表示经过细胞和基因过滤器后的原始输入表达矩阵（列中的细胞和行中的基因）。
# 基因通过 k 均值聚类，k = 100（左侧的树状图），热图表示 log2 缩放后基因簇中心的表达水平。

x <- sc3_plot_expression(data, k = k_est)
save_to_png(x, "sc3-expression", 8, 16)

# 稳定性指数显示每个群集在选定的 k 范围内的稳定程度。稳定性指数在 0 和 1 之间变化，
# 其中 1 意味着相同的簇出现在不同 k 的每个解中。

x <- sc3_plot_cluster_stability(data, k = k_est)
save_to_png(x, "sc3-cluster-stability", 4, 4)

# 使用非参数 Kruskal-Wallis 检验计算差异表达。显著的 p 值表明至少一个簇中的基因表达随机
# 支配另一个簇。 SC3 提供了调整后 p 值 < 0.01 的所有差异表达基因的列表，并绘制了具有
# 最低 p 值的 50 个基因的基因表达谱。请注意，聚类后差异表达的计算可能会在 p 值的分布中
# 引入偏差，因此我们建议仅使用 p 值对基因进行排名。

x <- sc3_plot_de_genes(data, k = k_est)
save_to_png(x, "sc3-differential-expression", 8, 16)

# 为了找到标记基因，对于每个基因，基于平均聚类表达值构建二元分类器。然后使用基因表达等级
# 计算分类器预测。受试者工作特征（ROC）曲线下面积用于量化预测的准确性。 通过使用 Wilcoxon
# 符号秩检验将 p 值分配给每个基因。 默认情况下，选择 ROC 曲线下面积（AUROC）> 0.85
# 且 p 值 < 0.01 的基因，并且每个聚类的前 10 个标记基因在该热图中可视化。

x <- sc3_plot_markers(data, k = k_est)
save_to_png(x, "sc3-markers", 8, 16)

# Differential expression: Statistical methods --------------------------------

assay(data, "combat") |> class() # a matrix
individual_group <- colData(data)["individual"] |> as.data.frame()
individual_group <- individual_group $ individual

## Kolmogorov-Smirnov Test (non-parametric) -----------------------------------

ks_pvals <- apply(
  assay(data, "combat"), 1, function(x) { # for each row
    ks.test(x[individual_group == "NA19101"],
            x[individual_group == "NA19239"]) $ p.value
  }
)

# multiple testing correction
ks_pvals <- p.adjust(ks_pvals, method = "fdr")

## Wilcox / Mann-Whitney U Test -----------------------------------------------

wilcox_pvals <- apply(
  assay(data, "combat"), 1, function(x) { # for each row
    wilcox.test(x[individual_group == "NA19101"],
                x[individual_group == "NA19239"]) $ p.value
  }
)

# multiple testing correction
wilcox_pvals <- p.adjust(wilcox_pvals, method = "fdr")

## EdgeR ----------------------------------------------------------------------

# We’ve already used edgeR for differential expression. edgeR is based on a
# negative binomial model of gene expression and uses a generalized linear
# model (GLM) framework, the enables us to include other factors such as
# batch to the model.

library(edgeR)

dge <- DGEList(
  counts = assay(data, "counts"),
  norm.factors = rep(1, length(assay(data, "counts")[1, ])),
  group = individual_group
)

edge_group <- factor(individual_group)
design <- model.matrix(~ edge_group)

# error: non-comfortable arrays. ?
dge <- estimateDisp(dge, design = design, trend.method = "none")

fit <- glmFit(dge, design)
res <- glmLRT(fit)
pvals <- res $ table[, 4]
names(pvals) <- rownames(res $ table)

pvals <- p.adjust(pvals, method = "fdr")

## MAST -----------------------------------------------------------------------

install_packages("MAST")
library(MAST)

# MAST is based on a zero-inflated negative binomial model. It tests for
# differential expression using a hurdle model to combine tests of discrete
# (0 vs not zero) and continuous (non-zero values) aspects of gene expression.
# Again this uses a linear modelling framework to enable complex models
# to be considered.

log_counts <- logcounts(data)
df <- data.frame(names = rownames(log_counts))
rownames(df) <- rownames(log_counts)

cond_data <- data.frame(cond = individual_group)
rownames(cond_data) <- colnames(log_counts)

obj <- FromMatrix(as.matrix(log_counts), cond_data, df)
colData(obj) $ cngeneson <- scale(colSums(assay(obj) > 0))
cond <- factor(colData(obj) $ cond)

# Model expression as function of condition & number of detected genes
zlm_cond <- zlm(~ cond + cngeneson, obj)

summary_cond <- summary(zlm_cond, doLRT = "condNA19239")
summary_dt <- summary_cond $ datatable

summary_dt <- as.data.frame(summary_dt)
pvals <- unlist(summary_dt[summary_dt $ component == "H", 4]) # H = hurdle model
names(pvals) <- unlist(summary_dt[summary_dt $ component == "H", 1])

pvals <- p.adjust(pvals, method = "fdr")
