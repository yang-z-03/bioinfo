
parser <- argparse::ArgumentParser(
  prog = "cluster",
  description = "clustering according to dimension reduction"
)

parser $ add_argument(
  "-m", dest = "method", type = "character", default = "",
  help = paste("method of dimension reduction, accepts 'sc3' or 'seurat'")
)

parser $ add_argument(
  "-k", dest = "k", type = "integer", default = 20,
  help = paste("number of nearest neighbors to find")
)

parser $ add_argument(
  "--expected", dest = "c", type = "integer", nargs = "+", default = 6:10,
  help = paste("the expected cluster number")
)

parser $ add_argument(
  "-a", dest = "alg", type = "integer", default = 2,
  help = paste("algorithm for modularity optimization (1 = original louvain",
               "algorithm; 2 = louvain algorithm with multilevel refinement;",
               "3 = slm algorithm; 4 = leiden algorithm). leiden requires",
               "the leidenalg python package to be installed")
)

parser $ add_argument(
  "-r", dest = "resolution", type = "double", default = 0.8,
  help = paste("resolution, the smaller, the finer. commonly between 0.5-1")
)

if (length(vargs) == 0 ||
      (length(vargs) == 1 && (vargs[1] == "-h" || vargs[1] == "--help"))) {
  parser $ print_help()
  stop()
} else {
  pargs <- parser $ parse_args(vargs)
}

switch(
  pargs $ method,
  sc3 = {
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(GetAssayData(shared[["seurat"]], layer = "counts")),
        logcounts = as.matrix(GetAssayData(shared[["seurat"]], layer = "data"))
      )
    )

    rowData(sce) $ feature_symbol <- rownames(sce)
    sce <- sc3(sce, ks = pargs $ c, biology = TRUE)
  },
  seurat = {

    # requires that you run PCA in dimreduc first.
    shared[["seurat"]] <- Seurat::FindNeighbors(
      shared[["seurat"]],
      k.param = pargs $ k
    )

    shared[["seurat"]] <- Seurat::FindClusters(
      shared[["seurat"]],
      resolution = pargs $ resolution,
      algorithm = pargs $ alg
    )

    # stores in the metadata column: seurat_clusters
    seurat @ active.ident <- seurat $ seurat_clusters
  }
)