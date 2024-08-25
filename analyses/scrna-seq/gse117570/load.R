
# processed data is contained in the gz files in tab-delimited format.

f <- function(t) {
  paste("geo", "GSE117570", t, sep = "/")
}

p1t <- read.delim(f("gsm3304007-p1-tumor.txt.gz"), header = TRUE, sep = "\t")
p1n <- read.delim(f("gsm3304008-p1-normal.txt.gz"), header = TRUE, sep = "\t") #
p2t <- read.delim(f("gsm3304009-p2-tumor.txt.gz"), header = TRUE, sep = "\t")
p2n <- read.delim(f("gsm3304010-p2-normal.txt.gz"), header = TRUE, sep = "\t") #
p3t <- read.delim(f("gsm3304011-p3-tumor.txt.gz"), header = TRUE, sep = "\t")
p3n <- read.delim(f("gsm3304012-p3-normal.txt.gz"), header = TRUE, sep = "\t") #
p4t <- read.delim(f("gsm3304013-p4-tumor.txt.gz"), header = TRUE, sep = "\t")
p4n <- read.delim(f("gsm3304014-p4-normal.txt.gz"), header = TRUE, sep = "\t") #

library(SingleCellExperiment)
library(scater)
library(Seurat)
library(reticulate)
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

qc <- function(table) {
  sce1t <- SingleCellExperiment(assays = list(counts = table |> as.matrix()))
  sce1t <- scater::logNormCounts(sce1t)

  genes <- sce1t |> rownames()
  mtgenes <- genes[grep("^MT-", genes, ignore.case = TRUE)]
  colData(sce1t) $ percent_mt <- colSums(table[mtgenes,]) / colSums(table)
  rbgenes <- genes[c(
    grep("^RPL", genes, ignore.case = TRUE),
    grep("^RPS", genes, ignore.case = TRUE)
  )]
  colData(sce1t) $ percent_rb <- colSums(table[rbgenes,]) / colSums(table)
  colData(sce1t) $ reads <- colSums(table)

  srat1t <- as.Seurat(sce1t, counts = "counts", data = "logcounts")
  srat1t <- RenameAssays(srat1t, assay.name = "originalexp",
                         new.assay.name = "RNA")

  srat1t <- scrublet_matrix(seurat_obj = srat1t,
                            matrix = srat1t @ assays $ RNA $ counts)
  srat1t[["is_doublet"]] <- srat1t[["predicted_doublets"]]

  list(qc_violin(srat1t), srat1t)
}

prepare_merge <- function(srat1t, patient) {
  len <- rownames(srat1t) |> length()
  srat1t[["patient"]] <- rep(patient, len)
  srat1t <- Seurat::ScaleData(srat1t)
  srat1t <- FindVariableFeatures(srat1t, nfeatures = 10000,
                                 selection.method = "vst")
  srat1t <- Seurat::RunPCA(srat1t)
  srat1t
}

qc1 <- qc(p1t)
srat1t <- qc1[[2]]
srat1t <- prepare_merge(srat1t, "p1")

qc2 <- qc(p2t)
srat2t <- qc2[[2]]
srat2t <- prepare_merge(srat2t, "p2")

qc3 <- qc(p3t)
srat3t <- qc3[[2]]
srat3t <- prepare_merge(srat3t, "p3")

qc4 <- qc(p4t)
srat4t <- qc4[[2]]
srat4t <- prepare_merge(srat4t, "p4")

samples <- list(srat1t, srat2t, srat4t)
features <- Seurat::SelectIntegrationFeatures(samples, nfeatures = 10000,
                                              fvf.nfeatures = 10000)
anchors <- Seurat::FindIntegrationAnchors(
  samples, k.anchor = 5, k.filter = 25, dims = 1:20,
  anchor.features = features
)
integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20)

var1 <- VariableFeatures(srat1t)
var2 <- VariableFeatures(srat2t)
var3 <- VariableFeatures(srat3t)
var4 <- VariableFeatures(srat4t)
vars <- c(var1, var2, var4)
vars <- vars[!duplicated(vars)]

integrated <- Seurat::ScaleData(integrated, features = vars)
integrated <- Seurat::RunPCA(integrated, features = vars)

# Only one parameter among 'dims', 'nn.name', 'graph', or 'features' should be
# used at a time to run UMAP
integrated <- Seurat::RunUMAP(
  integrated,
  n.neighbors = 20,
  min.dist = 0.3,
  dims = NULL
)

integrated <- Seurat::RunTSNE(integrated, perplexity = 50)
TSNEPlot(integrated)

integrated <- Seurat::FindNeighbors(integrated, dims = 1:2, reduction = "tsne")
integrated <- Seurat::FindClusters(integrated,
                                   resolution = 0.1)
integrated @ meta.data |> str()
integrated[["patient"]] |> table()

plot_tsne <- function(srat) {
  data <- srat[["tsne"]] @ cell.embeddings |> as.data.frame()
  data $ clusters <- srat[["seurat_clusters"]][[1]] |> factor()
  data $ patient <- srat[["patient"]][[1]] |> factor()
  ggplot(data, aes(x = .data $ tSNE_1,
                   y = .data $ tSNE_2)) +
    geom_point(aes(color = factor(.data $ patient)))
}

plot_expression_tsne <- function(srat, gene) {
  gene_row <- GetAssayData(integrated, layer = "data")[gene, ]
  names(gene_row) <- NULL
  data <- srat[["tsne"]] @ cell.embeddings |> as.data.frame()
  data $ clusters <- srat[["seurat_clusters"]][[1]]
  data $ expr <- gene_row
  ggplot(data, aes(x = .data $ tSNE_1,
                   y = .data $ tSNE_2)) +
    geom_point(aes(color = .data $ expr)) +
    scale_color_gradient(low = "white", "high" = "red")
}

plot_tsne(integrated)
plot_expression_tsne(integrated, "IL17RA")
genes <- rownames(integrated)
genes[grep("IL17", genes)]

# annotation ===================================================================

integrated <- readRDS(f("integrated.rds"))

BiocManager::install("scAnnotatR", Ncpus = 44)
BiocManager::install("SingleR", Ncpus = 44)

# library(scAnnotatR)
source("init.R")
install_packages(c("celldex", "SingleR"), j = 44)

library(dplyr)
library(Matrix)
library(SingleR)
library(celldex)

# well, our gene features are too small ...

# default_models <- load_models("default")
# integrated <- classify_cells(classify_obj = integrated,
#                              assay = "integrated", slot = "data",
#                              cell_types = "all",
#                              path_to_models = "default")

all_markers <- Seurat::FindAllMarkers(
  integrated, only.pos = TRUE,
  min.pct = 0.5, logfc.threshold = 0.5
)

dim(all_markers)
table(all_markers $ cluster)
top3_markers <- as.data.frame(
  all_markers |> group_by(cluster) |>
    top_n(n = 3, wt = avg_log2FC)
)
top1_markers <- as.data.frame(
  all_markers |> group_by(cluster) |>
    top_n(n = 1, wt = avg_log2FC)
)

Seurat::DoHeatmap(object = integrated, features = top3_markers $ gene)
Seurat::VlnPlot(object = integrated, features = top1_markers $ gene)
Seurat::FeaturePlot(object = integrated, features = top1_markers $ gene)
Seurat::DotPlot(object = integrated, features = top1_markers $ gene)

immgen <- celldex::ImmGenData()
sce_int <- as.SingleCellExperiment(DietSeurat(integrated))
pred_annot <- SingleR::SingleR(
  # SeuratObject::GetAssayData(integrated, layer = "data"),
  # integrated @ assays $ integrated $ data,
  test = sce_int,
  ref = immgen, labels = immgen $ label.main
)
