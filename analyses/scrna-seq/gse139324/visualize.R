
require(Seurat)
require(ggplot2)
require(patchwork)
require(RColorBrewer)

f <- function(x) {
  paste("geo", "GSE139324", x, sep = "/")
}

integrated <- readRDS("geo/GSE139324/integrated.rds")
tsub <- readRDS("geo/GSE139324/t-cells.rds")

t_all_markers <- Seurat::FindAllMarkers(
  tsub, only.pos = TRUE,
  min.pct = 0.5, logfc.threshold = 0.5
)

dim(t_all_markers)
table(t_all_markers $ cluster)
t_top5_markers <- as.data.frame(
  t_all_markers |> group_by(cluster) |> top_n(n = 5, wt = avg_log2FC)
)

# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[1:5])   # KLRG1, GZMK, C1orf21
# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[6:10])  # ATP6V0A1, GNG8, IL1R1
# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[11:15]) # TYMS, TK1, CDT1, PRR11
# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[16:20]) # IGFBP7, KLRF1
# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[21:25]) # ICOS, CCR7
# Seurat::VlnPlot(object = tsub, features = t_top5_markers $ gene[26:30]) # CD180
# Seurat::VlnPlot(object = tsub, features = c("CXCR3"))
#
# violins <- Seurat::VlnPlot(object = tsub, features = c(
#   "GZMK", "KLRG1", "IL1R1", "GNG8", "TYMS", "KLRF1", "CCR7", "CD180"
# ))
#
# ggsave("violin.png", violins, width = 8, height = 10, units = "in", dpi = 600)
#
# Seurat::FeaturePlot(object = tsub, features = c(
#   "GZMK", "KLRG1", "IL1R1", "GNG8", "TYMS", "KLRF1", "CCR7", "CD180"
# ))
#
# Seurat::VlnPlot(object = tsub, features = c(
#   "IL17A", "IL17F", "IL23A"
# ))
#
# Seurat::FeaturePlot(object = tsub, features = c(
#   "IL17A", "IL17F", "IL23A"
# ))

genes_to_draw <- c("GZMK", "KLRG1", "IL1R1", "GNG8", "TYMS", "IFNG",
                   "KLRF1", "CCR7", "CD180", "IL17A", "IL17F", "IL23A")

violins <- list()
expr_umaps <- list()
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
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    scale_color_gradient(low = "white", "high" = "red")

  n <- n + 1
}

violin_leg <- Seurat::VlnPlot(object = tsub, features = g, pt.size = 0.01,
                              alpha = 0.3, combine = FALSE) [[1]] +
  unify_theme_font(9, "Arial") +
  ggplot2::scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
  ggplot2::xlab("Cluster")

violins[[1]]
expr_umaps[[1]]

pw_violin <-
  (violins[[1]] | violins[[2]] | violins[[3]] | violins[[4]]) /
  (violins[[5]] | violins[[6]] | violins[[7]] | violins[[8]]) /
  (violins[[9]] | violins[[10]] | violins[[11]] | violins[[12]])

pw_umaps <-
  (expr_umaps[[1]] | expr_umaps[[2]] | expr_umaps[[3]] | expr_umaps[[4]]) /
  (expr_umaps[[5]] | expr_umaps[[6]] | expr_umaps[[7]] | expr_umaps[[8]]) /
  (expr_umaps[[9]] | expr_umaps[[10]] | expr_umaps[[11]] | expr_umaps[[12]])

ggsave(f("violins.pdf"), pw_violin, width = 10, height = 6, units = "in")
ggsave(f("violins.png"), pw_violin, width = 10, height = 6, units = "in", dpi = 600)
ggsave(f("violin-legend.pdf"), violin_leg, width = 5, height = 2, units = "in")
ggsave(f("umaps.pdf"), pw_umaps, width = 12, height = 8, units = "in")
ggsave(f("umaps.png"), pw_umaps, width = 12, height = 8, units = "in", dpi = 600)

classif <- rep(0, length(tsub $ patient))
classif[tsub $ patient == 4 | tsub $ patient == 7 | tsub $ patient == 10 |
          tsub $ patient == 11 | tsub $ patient == 12 | tsub $ patient == 13 |
          tsub $ patient == 14 | tsub $ patient == 15 | tsub $ patient == 16 |
          tsub $ patient == 17 | tsub $ patient == 18 | tsub $ patient == 19] <- 1

integrated $ class <- classif
tsub $ class <- classif

t_catdata <- table(Idents(tsub), tsub $ class)
t_catdata <- as.data.frame(t_catdata)
t_stack <- draw_stacked_bars(t_catdata) +
  ggplot2::scale_x_discrete(labels = c("CT", "RR")) +
  ggplot2::xlab("Group")

ggsave(f("stack-t.pdf"), t_stack, width = 5, height = 4, units = "in")
ggsave(f("stack-t.png"), t_stack, width = 5, height = 4, units = "in", dpi = 600)

catdata <- table(Idents(integrated), integrated $ class)
catdata <- as.data.frame(catdata)
all_stack <- draw_stacked_bars(catdata) +
  ggplot2::scale_x_discrete(labels = c("CT", "RR")) +
  ggplot2::xlab("Group")

ggsave(f("stack-all.pdf"), all_stack, width = 5, height = 4, units = "in")
ggsave(f("stack-all.png"), all_stack, width = 5, height = 4, units = "in", dpi = 600)

tsne_int <- DimPlot(integrated)
umap_tsub <- DimPlot(tsub)

ggsave(f("tsne-all.pdf"), tsne_int, width = 7, height = 4, units = "in")
ggsave(f("tsne-all.png"), tsne_int, width = 7, height = 4, units = "in", dpi = 600)

ggsave(f("umap-t.pdf"), umap_tsub, width = 7, height = 4, units = "in")
ggsave(f("umap-t.png"), umap_tsub, width = 7, height = 4, units = "in", dpi = 600)
