
library(SingleCellExperiment)
library(scater)
library(Seurat)

mouse <- readRDS("geo/GSE127465/sce.mouse.rds")

mouse <- scater::runPCA(mouse, exprs_values = "normlog")
mconf <- plotExplanatoryVariables(mouse, exprs_values = "normlog",
                                  variables = c("rep", "tissue", "batch",
                                                "reads", "percent_mt", "annot_maj"))

set.seed(42)
mouse <- scater::runUMAP(
  mouse, n_threads = 40,
  n_neighbors = 8, pca = 50,
  exprs_values = "normlog",
  BPPARAM = BiocParallel::MulticoreParam())

mumap <- scater::plotUMAP(
  mouse,
  colour_by = "annot_maj"
)

m_gene_names <- mouse |> rownames()
acod1 <- grep('irg', m_gene_names, ignore.case = TRUE)
m_gene_names[acod1]

# in mice, Acod1 is named as Irg1?
# plot expression of Irg1
mirg <- plotExpression(mouse, "Irg1", exprs_values = "normlog", log2_values = TRUE)
mirg <- scater::plotUMAP(
  mouse, dimred = "UMAP",
  colour_by = "Irg1", by_exprs_values = "normlog"
)

ggsave(
  "geo/GSE127465/mouse.confounders.png",
  plot = mconf + unify_theme_font(),
  width = 6, height = 4, units = "in", dpi = 600
)

ggsave(
  "geo/GSE127465/mouse.umap.png",
  plot = mumap + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

ggsave(
  "geo/GSE127465/mouse.expression.irg1.png",
  plot = mirg + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

################################################################################

human <- readRDS("geo/GSE127465/sce.human.rds")

# calculate differential expression.

human $ tissue |> table()
hb <- human[, human $ tissue == "blood"]
ht <- human[, human $ tissue == "tumor"]
ht <- ht[, ht $ is_tumor_immune == "True"]
ht <- ht[, ht $ annot_maj != "tRBC"]

h_gene_names <- human |> rownames()
acod1 <- grep('acod1', h_gene_names, ignore.case = TRUE)
h_gene_names[acod1]

# run for human blood sample
hb <- scater::runPCA(hb, exprs_values = "normlog")
ht <- scater::runPCA(ht, exprs_values = "normlog")

hbconf <- plotExplanatoryVariables(hb, exprs_values = "normlog",
                                   variables = c("sample", "reads",
                                                 "percent_mt", "annot_maj"))

htconf <- plotExplanatoryVariables(ht, exprs_values = "normlog",
                                   variables = c("sample", "reads",
                                                 "percent_mt", "annot_maj"))

ggsave(
  "geo/GSE127465/hb.confounders.png",
  plot = hbconf + unify_theme_font(),
  width = 6, height = 4, units = "in", dpi = 600
)

ggsave(
  "geo/GSE127465/ht.confounders.png",
  plot = htconf + unify_theme_font(),
  width = 6, height = 4, units = "in", dpi = 600
)

set.seed(42)
hb <- scater::runUMAP(
  hb, n_threads = 40,
  n_neighbors = 5, pca = 10,
  exprs_values = "normlog",
  verbose = TRUE,
  min_dist = 0.05,
  BPPARAM = BiocParallel::MulticoreParam()
)

set.seed(42)
ht <- scater::runUMAP(
  ht, n_threads = 44,
  n_neighbors = 20, pca = 30,
  exprs_values = "normlog",
  verbose = TRUE,
  min_dist = 0.05,
  BPPARAM = BiocParallel::MulticoreParam()
)

# the umap plot.
hbumap <- scater::plotUMAP(
  hb,
  colour_by = "annot_maj"
)

ggsave(
  "geo/GSE127465/hb.umap.png",
  plot = hbumap + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

htumap <- scater::plotUMAP(
  ht,
  colour_by = "annot_maj"
)

ggsave(
  "geo/GSE127465/ht.umap.png",
  plot = htumap + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

# plot expression of Acod1
hbirg <- scater::plotUMAP(
  hb, dimred = "UMAP",
  colour_by = "ACOD1", by_exprs_values = "normlog"
)

ggsave(
  "geo/GSE127465/hb.expression.acod1.png",
  plot = hbirg + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

htirg <- scater::plotUMAP(
  ht, dimred = "UMAP",
  colour_by = "ACOD1", by_exprs_values = "normlog"
)

ggsave(
  "geo/GSE127465/ht.expression.acod1.png",
  plot = htirg + unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

library(RColorBrewer)

# plot other marker genes in human.
ggsave(
  "geo/GSE127465/ht.expression.cd19.png",
  plot = scater::plotUMAP(
    ht, dimred = "UMAP",
    colour_by = "CD19", by_exprs_values = "normlog"
  ) +
    scale_fill_continuous(low = "black", high = "red", limits = c(0, 1)) +
    scale_fill_brewer(palette="Set2") +
    unify_theme_font(),
  width = 10, height = 8, units = "in", dpi = 600
)

rescale_umap <- function(gene, dataname, sce) {
  subset <- sce[gene, ]
  subset_data <- assay(subset, "normlog")[gene, ]
  names(subset_data) <- NULL
  nonzero <- subset_data > 0
  nzpart <- subset_data[nonzero]
  ord <- order(nzpart) / length(nzpart)
  assay(subset, "rescale") <- assay(subset, "normlog")
  assay(subset, "rescale")[gene, nonzero] <- ord

  p <- scater::plotUMAP(
    subset, dimred = "UMAP",
    colour_by = gene, by_exprs_values = "rescale"
  ) + unify_theme_font()

  xrange <- (layer_scales(p) $ x $ range $ range)
  yrange <- (layer_scales(p) $ y $ range $ range)

  p2 <- scater::plotUMAP(
    subset[, nonzero], dimred = "UMAP",
    colour_by = gene, by_exprs_values = "rescale"
  ) + unify_theme_font() +
    scale_x_continuous(limits = xrange) +
    scale_y_continuous(limits = yrange)

  ggsave(
    paste("geo/GSE127465/", dataname, ".expression.",
          gene |> tolower(), ".png", sep = ""),
    plot = p + p2,
    width = 16, height = 8, units = "in", dpi = 600
  )

  p
}

# human tumor marker genes
rescale_umap("CD19", "ht", ht)
rescale_umap("CD3E", "ht", ht)
rescale_umap("NCR1", "ht", ht)
rescale_umap("CSF1R", "ht", ht)
rescale_umap("CSF3R", "ht", ht)
rescale_umap("LILRA4", "ht", ht)
rescale_umap("ACOD1", "ht", ht)

# human blood
rescale_umap("CD19", "hb", hb)
rescale_umap("CD3E", "hb", hb)
rescale_umap("NCR1", "hb", hb)
rescale_umap("CSF1R", "hb", hb)
rescale_umap("CSF3R", "hb", hb)
rescale_umap("LILRA4", "hb", hb)
rescale_umap("ACOD1", "hb", hb)

# mouse
rescale_umap("Irg1", "mouse", mouse)
rescale_umap("Cd19", "mouse", mouse)
rescale_umap("Cd3e", "mouse", mouse)
rescale_umap("Ncr1", "mouse", mouse)
rescale_umap("Csf1r", "mouse", mouse)
rescale_umap("Csf3r", "mouse", mouse)
rescale_umap("Siglech", "mouse", mouse)
