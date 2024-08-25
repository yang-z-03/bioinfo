
library(reticulate)
library(fs)

library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)
library(patchwork)
library(SoupX)
library(RColorBrewer)

reticulate::py_config()

series_id <- "GSE210132"

# download_table <- geo_get_supplemental_links(series_id) # nolint
# geo_get_supplemental_files(series_id)                   # nolint

expr_set <- geo_get_series_annot(series_id)
platform <- geo_get_platform_annot(expr_set)

arrange_files_to_10x <- function(root_dir) {
  files <- list.files(root_dir)
  matrix <- strsplit(files, split = "_")

  i <- 1
  for (x in matrix) {
    attempt_dir <- paste(root_dir, x[1], sep = "/")
    if (!dir.exists(attempt_dir)) {
      dir.create(attempt_dir)
    }

    print(i)
    print(files[i])
    print(x[1])
    print(last(x))

    file_move(
      paste(root_dir, files[i], sep = "/"),
      paste(root_dir, x[1], last(x), sep = "/")
    )

    i <- i + 1
  }
}

root_dir <- paste("geo", series_id, sep = "/")
arrange_files_to_10x(root_dir)

sample_id_s <- function(id) {
  paste("geo/GSE210132/GSM6422", id, sep = "")
}

sample_list <- list.dirs(root_dir, recursive = FALSE)
# choose one sample first. this one, #133, is ntCZC with NALCN-wt.

load_sample <- function(id) {
  sample <- sample_id_s(id)

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

  soup_profile <- data.frame(
    row.names = rownames(soup_channel $ toc),
    est = rowSums(soup_channel $ toc) / sum(soup_channel $ toc),
    counts = rowSums(soup_channel $ toc)
  )
  soup_channel <- setSoupProfile(soup_channel, soup_profile)

  # provide a temp cluster
  qclust <- scran::quickCluster(rt, min.size = 30)
  soup_channel <- SoupX::setClusters(soup_channel, setNames(
    qclust, rt |> colnames()
  ))

  # estimate and clean the data using SoupX
  soup_channel <- SoupX::autoEstCont(soup_channel)
  soup_channel <- SoupX::adjustCounts(soup_channel)

  segs <- strsplit(sample, split = "/")[[1]]
  result <- Seurat::CreateSeuratObject(
    counts = soup_channel, project = paste0("sample-", segs[3]),
    min.cells = 3, min.features = 200
  )

  return(result)
}

ntczc_wt <- load_sample(133)
current_data <- ntczc_wt

analyse_data <- function(current_data) {

  # take a look at the metadata
  meta <- current_data @ meta.data
  head(meta)
  summary(meta $ nCount_RNA)
  summary(meta $ nFeature_RNA)

  # gene sets vary in nomenclature. this is rather annoying. we can see the
  # names are camel/pascal cases here.
  # there are only 13 mitochondrial genes in human. you may just memorize that.

  # x_id <- grep('mt-', current_data |> rownames(), ignore.case = T)
  # (current_data |> rownames())[x_id]

  #  [1] "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6" "mt-Co3"
  #  [8] "mt-Nd3"  "mt-Nd4l" "mt-Nd4"  "mt-Nd5"  "mt-Nd6"  "mt-Cytb"

  current_data[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(current_data, pattern = "^mt-")
  current_data[["percent.rb"]] <-
    Seurat::PercentageFeatureSet(current_data, pattern = "^Rp[sl]")

  current_data <- scrublet(seurat_obj = current_data)
  current_data[["is_doublet"]] <- current_data[["predicted_doublets"]]

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

  # QC

  # 如果 QC 合格，最终显示为 pass.

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ is_doublet == TRUE,
    "doublet", "pass"
  )

  # 注意，这里的值要根据小提琴图来看
  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 1000 &
      current_data @ meta.data $ qc == "pass",
    "low_n_feature", current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 1000 &
      current_data @ meta.data $ qc != "pass" &
      current_data @ meta.data $ qc != "low_n_feature",
    paste("low_n_feature", current_data @ meta.data $ qc, sep = ":"),
    current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ percent.mt > 6 &
      current_data@ meta.data $ qc == "pass", # nolint
    "high_mt", current_data @ meta.data $ qc
  )

  current_data[["qc"]] <- ifelse(
    current_data @ meta.data $ nFeature_RNA < 1000 &
      current_data @ meta.data $ qc != "pass" &
      current_data @ meta.data $ qc != "high_mt",
    paste("high_mt", current_data @ meta.data $ qc, sep = ":"),
    current_data @ meta.data $ qc
  )

  table(current_data[["qc"]])
  qc_violin(current_data |> subset(qc == "pass"))

  current_data <- Seurat::NormalizeData(current_data)
  current_data <- Seurat::FindVariableFeatures(current_data,
                                               selection.method = "vst",
                                               nfeatures = 2000)

  # the top10 variable genes.
  top10 <- head(Seurat::VariableFeatures(current_data), 10)
  top10_plot <- Seurat::VariableFeaturePlot(current_data)
  top10_plot <- Seurat::LabelPoints(
    plot = top10_plot,
    points = top10,
    repel = TRUE,
    xnudge = 0, ynudge = 0,
    size = 4.0
  ) + unify_theme_font()

  # variance pca plots
  all_genes <- rownames(current_data)
  current_data <- Seurat::ScaleData(current_data, features = all_genes)
  current_data <- RunPCA(
    current_data,
    features = Seurat::VariableFeatures(object = current_data)
  )

  # the top and bottom pca in graphs
  pc_loadings <- Seurat::VizDimLoadings(
    current_data, dims = 1:9, reduction = "pca"
  )

  for(x in seq_along(pc_loadings)) {
    pc_loadings[[x]] <- pc_loadings[[x]] + unify_theme_font()
  }

  # pca heatmap
  pc_heatmap <- Seurat::DimHeatmap(
    current_data, dims = 1:9, reduction = "pca",
    balanced = TRUE, cells = 50,
    fast = FALSE, combine = TRUE
  )

  for(x in seq_along(pc_heatmap)) {
    pc_heatmap[[x]] <- pc_heatmap[[x]] +
      unify_theme_font() +
      theme(axis.text.x = element_blank()) +
      scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"))
  }

  dim_plot_pca <- Seurat::DimPlot(current_data, reduction = "pca") +
    unify_theme_font()
  elbow_plot_pca <- Seurat::ElbowPlot(current_data, reduction = "pca") +
    unify_theme_font()

  # by default now, we did not remove the items that doesn't pass UC.
  # clustering, but using primitive automatic methods.
  current_data <- Seurat::FindNeighbors(current_data, dims = 1:10)
  current_data <- Seurat::FindClusters(current_data, resolution = 0.5)
  table(current_data @ meta.data $ seurat_clusters)

  # run umap. it performs better than pca.
  current_data <- Seurat::RunUMAP(current_data, dims = 1:10, verbose = FALSE)
  dim_plot_umap <- Seurat::DimPlot(
    current_data,
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  # remove those failed with QC.
  dim_plot_umap_qc <- Seurat::DimPlot(
    current_data |> subset(qc == "pass"),
    label.size = 4,
    repel = TRUE, label = TRUE
  ) + unify_theme_font()

  # see the distribution of qc-failed cells.
  mt_feature <- Seurat::FeaturePlot(
    current_data, features = "percent.mt",
    label.size = 4, repel = TRUE, label = TRUE
  ) + unify_theme_font()

  rb_feature <- Seurat::FeaturePlot(
    current_data, features = "percent.rb",
    label.size = 4, repel = TRUE, label = TRUE
  ) + unify_theme_font()

  # cell cycle
  cc_genes <- Seurat::cc.genes
  cc_genes_2019 <- Seurat::cc.genes.updated.2019

  # adapt the gene table (hopefully) to the wierd nomenclature
  current_nomen <- current_data |> rownames()
  cc_genes_2019[[1]] <- map(cc_genes_2019[[1]], function(x) {
    current_nomen[grep(paste("^", x, "$", sep = ""),
                       current_nomen, ignore.case = TRUE)]
  }) |> unlist()
  cc_genes_2019[[2]] <- map(cc_genes_2019[[2]], function(x) {
    current_nomen[grep(paste("^", x, "$", sep = ""),
                       current_nomen, ignore.case = TRUE)]
  }) |> unlist()

  current_data <- Seurat::CellCycleScoring(
    current_data,
    s.features = cc_genes_2019 $ s.genes,
    g2m.features = cc_genes_2019 $ g2m.genes
  )

  table(current_data[[]] $ Phase)

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

  # re-run the pca, umap and cluster
  current_data <- Seurat::RunPCA(current_data, verbose = FALSE)
  current_data <- Seurat::RunUMAP(current_data, dims = 1:30, verbose = FALSE)
  current_data <- Seurat::FindNeighbors(current_data, dims = 1:30,
                                        verbose = FALSE)
  current_data <- Seurat::FindClusters(current_data, verbose = FALSE)
  table(current_data[[]] $ seurat_clusters)

  # normalized umap can be compared between datasets.
  # and it looks less wierd for border conditions.
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
  # 同时仅报告阳性细胞。有许多测试可以用来定义标记，包括非常快速和直观的 tf-idf

  # 默认情况下，使用 Wilcoxon 秩和检验。这需要一段时间！为了提高速度，我们增加了默认的
  # 最小百分比和 log 2FC 截止值；这些应该进行调整以适应数据集！

  all_markers <- Seurat::FindAllMarkers(
    current_data, only.pos = TRUE,
    min.pct = 0.5, logfc.threshold = 0.5
  )

  dim(all_markers)
  table(all_markers $ cluster)
  top10_markers <- as.data.frame(
    all_markers |> group_by(cluster) |>
      top_n(n = 10, wt = avg_log2FC)
  )

  draw_single_gene_expr <- function(gene, description) {
    Seurat::FeaturePlot(current_data, gene, label = TRUE) +
      scale_colour_gradientn(
        colours = rev(brewer.pal(n = 11, name = "Spectral"))
      ) +
      ggtitle(paste0(gene, ": ", description)) +
      unify_theme_font()
  }

  # see valid gene names ...
  all_genes[grep('mmp', all_genes, ignore.case = TRUE)]
  draw_single_gene_expr("Mmp9", "MMP9")

  return(current_data)
}