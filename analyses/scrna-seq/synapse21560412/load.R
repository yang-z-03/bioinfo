
library(Seurat)
library(ggplot2)

p1 <- readRDS("synapse/21560412/droplet-lung-blood-p1.rds")

names(p1)
dim(p1)

p1 @ meta.data |> str()

draw_single_gene_expr <- function(gene, description, data) {
  Seurat::FeaturePlot(data, gene, label = TRUE) +
    ggplot2::scale_colour_gradientn(
      colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
    ) +
    ggplot2::ggtitle(paste0(gene, ": ", description))
}

human_genes <- rownames(p1)
human_genes[grep("irf7", human_genes, ignore.case = TRUE)]

tsne <- TSNEPlot(p1)
pirg1 <- draw_single_gene_expr("IRG1", "ACOD1", p1)

ggsave(
  "synapse/21560412/p1.tsne.png", tsne,
  width = 10, height = 8, units = "in", dpi = 600
)

ggsave(
  "synapse/21560412/p1.expression.irg1.png", pirg1,
  width = 10, height = 8, units = "in", dpi = 600
)

################################################################################

p3 <- readRDS("synapse/21560412/droplet-lung-blood-p3.rds")

tsne3 <- TSNEPlot(p3)
pirg3 <- draw_single_gene_expr("IRG1", "ACOD1", p3)

ggsave(
  "synapse/21560412/p3.tsne.png", tsne3,
  width = 15, height = 8, units = "in", dpi = 600
)

ggsave(
  "synapse/21560412/p3.expression.irg1.png", pirg3,
  width = 10, height = 8, units = "in", dpi = 600
)

ggsave(
  "synapse/21560412/p3.expression.sftpc.png",
  draw_single_gene_expr("SFTPC", "SFTPC, Marker of AT2", p3),
  width = 10, height = 8, units = "in", dpi = 600
)

################################################################################

p2 <- readRDS("synapse/21560412/droplet-lung-p2.rds")

load("synapse/21560412/facs-lung-blood-p1.robj")
load("synapse/21560412/facs-lung-blood-p2.robj")
load("synapse/21560412/facs-lung-blood-p3.robj")

f1 <- UpdateSeuratObject(ntiss.P1.anno.gencode)
f2 <- UpdateSeuratObject(ntiss.P2.anno.gencode)
f3 <- UpdateSeuratObject(ntiss.P3.anno.gencode)

f1 <- FindVariableFeatures(f1, nfeatures = 10000, selection.method = "vst")
f1 <- RunPCA(f1)
f1 <- RunTSNE(f1, perplexity = 30)

human_genes <- rownames(f2)
human_genes[grep("acod", human_genes, ignore.case = TRUE)]

reveal_acod1 <- function(name, data) {
  ggsave(
    paste("synapse/21560412/", name, ".tsne.png", sep = ""),
    TSNEPlot(data),
    width = 12, height = 8, units = "in", dpi = 600
  )
  ggsave(
    paste("synapse/21560412/", name, ".expression.acod1.png", sep = ""),
    draw_single_gene_expr("ACOD1", "ACOD1", data),
    width = 10, height = 8, units = "in", dpi = 600
  )
}

reveal_acod1("facs1", f1)
reveal_acod1("facs2", f2)
reveal_acod1("facs3", f3)

# for now, all f1, f2, and f3 datasets contains normalized data (data).
# and have picked variable features. now let's merge them into a whole.

samples <- list()
samples[[1]] <- f1
samples[[2]] <- f2
samples[[3]] <- f3

Seurat::SelectIntegrationFeatures(samples)
anchors <- Seurat::FindIntegrationAnchors(
  samples, k.anchor = 5, k.filter = 25, dims = 1:20
)
integrated <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20)

var1 <- VariableFeatures(f1)
var2 <- VariableFeatures(f2)
var3 <- VariableFeatures(f3)
vars <- c(var1, var2, var3)
vars <- vars[!duplicated(vars)]

integrated <- Seurat::ScaleData(integrated, features = vars)
integrated <- Seurat::RunPCA(integrated, features = vars)
integrated <- Seurat::RunTSNE(integrated, perplexity = 30)

reveal_acod1("facs", integrated)
