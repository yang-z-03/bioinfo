
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
  "--expected", dest = "c", type = "integer", nargs = "+", default = 8,
  help = paste("the expected cluster number")
)

parser $ add_argument(
  "--ncore", dest = "core", type = "integer", default = 100,
  help = paste("cpu cores to run sc3")
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
        counts = GetAssayData(shared[["seurat"]], layer = "counts") |> as("sparseMatrix"),
        logcounts = GetAssayData(shared[["seurat"]], layer = "data") |> as("sparseMatrix")
      )
    )

    gene_meta <- shared[["meta_gene"]]
    sample_meta <- shared[["meta_sample"]]
    rowData(sce) $ feature_symbol <- gene_meta $ name

    # one opened parallel worker consumes one connection. and R is hard-coded
    # with a maximum of 128 connections.

    sce <- sc3_estimate_k(sce)
    k <- (metadata(sce) $ sc3) $ k_estimation

    cat(blue("the estimation of k:"), k, crlf)
    cat(blue("running sc3"), crlf)

    # run sc3.
    suppressMessages({
      sce <- sc3(sce, ks = pargs $ c, biology = TRUE, n_cores = pargs $ core)
    })

    # plotting consensus matrix. the consensus matrix is a diagonal matrix with
    # ncells * ncells in size. when the cell number is extremely high (e.g.
    # greater than 5000), you should not attempt to plot this matrix. since the
    # process will be extremely slow.

    #. sc3_plot_consensus(sce, k)

    plotk <- paste("sc3", pargs $ c[1], "clusters", sep = "_")
    if (!dir.exists("sc3"))
      dir.create("sc3")

    # this object is also extremely large :(

    #. cat(blue("plotting cluster stability ..."), crlf)
    #. sc3stab <- sc3_plot_cluster_stability(sce, k = pargs $ c[1])
    #. saveRDS(sc3stab, "sc3/cluster-stability.rds")

    cat(blue("plotting cluster markers ..."), crlf)
    sc3mark <- sc3_plot_markers(sce, k = pargs $ c[1], show_pdata = plotk)
    saveRDS(sc3mark, "sc3/markers.rds")

    cat(blue("plotting differential expressions ..."), crlf)
    sc3de <- sc3_plot_de_genes(sce, k = pargs $ c[1], show_pdata = plotk)
    saveRDS(sc3de, "sc3/de.rds")

    # extract the biological information back to seurat.
    # the sample-specific metadata can be directly stored in seurat

    for (cold in sce |> colData() |> colnames()) {
      shared[["seurat"]][[cold]] <- colData(sce)[[cold]]
    }

    for (cold in sce |> colData() |> colnames()) {
      sample_meta[[cold]] <- colData(sce)[[cold]]
    }

    # the gene specific (e.g. marker genes, differential expressions etc.)
    # will be stored in append to the genes-meta.rds

    for (rowd in sce |> rowData() |> colnames()) {
      if (rowd == "feature_symbol") {
        next
      } else {
        gene_meta[[rowd]] <- rowData(sce)[[rowd]]
      }
    }

    #. saveRDS(gene_meta, "norm/genes-meta.rds")
    #. saveRDS(sample_meta, "norm/samples-meta.rds")
    #. saveRDS(shared[["seurat"]], "norm/seurat.rds")
    
    # we should avoid too much file saving operations, since when the dimensions
    # of the experiments grow large, these disk operation take too much time
    # at an unnecessary frequency.
    
    shared[["meta_gene"]] <- gene_meta
    shared[["meta_sample"]] <- sample_meta
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
    
    # we should avoid too much file saving operations, since when the dimensions
    # of the experiments grow large, these disk operation take too much time
    # at an unnecessary frequency.
    
    Idents(shared[["seurat"]]) <- "seurat_clusters"
    shared[["meta_sample"]] $ seurat_clusters <-
      shared[["seurat"]] $ seurat_clusters
    
    #. saveRDS(shared[["seurat"]], "norm/seurat.rds")
    #. saveRDS(shared[["meta_sample"]], "norm/samples-meta.rds")
  }
)