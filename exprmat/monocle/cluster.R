
# 2.  SUPERVISED CLUSTERING USING MONOCLE ======================================

# the functions are organized into a small data structure called a
# celltypehierarchy, that monocle uses to classify the cells. you first
# initialize a new celltypehierarchy object, then register your gating
# functions within it. once the data structure is set up, you can use it to
# classify all the cells in the experiment:

myf5 <- row.names(subset(fData(expr), gene_short_name == "MYF5"))
anpep <- row.names(subset(fData(expr), gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()

cth <- addCellType(
  cth, "Myoblast",
  classify_func = function(x) {
    x[myf5, ] >= 1
  }
)

cth <- addCellType(
  cth, "Fibroblast",
  classify_func = function(x) {
    x[myf5, ] < 1 & x[anpep, ] > 1
  }
)

# first, we'll select a different set of genes to use for clustering the cells.
# before we just picked genes that were highly expressed and highly variable.
# now, we'll pick genes that co-vary with our markers. in a sense, we'll be
# building a large list of genes to use as markers, so that even if a cell
# doesn't have myf5, it might be recognizable as a myoblast based on other
# genes.

marker_diff <- markerDiffTable(
  expr[expressed_genes, ],
  cth, residualModelFormulaStr = "~ Media + num_genes_expressed",
  cores = 1
)

# the function markerdifftable takes a celldataset and a celltypehierarchy and
# classifies all the cells into types according to your provided functions.
# it then removes all the "unknown" and "ambiguous" functions before
# identifying genes that are differentially expressed between the types. again,
# you can provide a residual model of effects to exclude from this test. the
# function then returns a data frame of test results, and you can use this to
# pick the genes you want to use for clustering. often it's best to pick the
# top 10 or 20 genes that are most specific for each cell type. this ensures
# that the clustering genes aren't dominated by markers for one cell type.
# you generally want a balanced panel of markers for each type if possible.
# monocle provides a handy function for ranking genes by how restricted their
# expression is for each type.

candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))

marker_spec <-
  calculateMarkerSpecificity(expr[candidate_clustering_genes, ], cth)

head(selectTopMarkers(marker_spec, 3))

# to cluster the cells, we'll choose the top 500 markers for each of these
# cell types:

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500) $ gene_id)
expr <- setOrderingFilter(expr, semisup_clustering_genes)

plot_ordering_genes(expr)

# here, we pass the selected filter to cell clustering.

expr <- clusterCells(expr, num_clusters = 2)


## 2.2  imputing cell types ----------------------------------------------------
##      (tagging whole clusters, or cluster annotation)

# Note that we've reduced the number of "contaminating" fibroblasts in the
# myoblast cluster, and vice versa. But what about the "Unknown" cells? If you
# provide clusterCells with a the CellTypeHierarcy, Monocle will use it
# classify whole clusters, rather than just individual cells. Essentially,
# clusterCells works exactly as before, except after the clusters are built,
#it counts the frequency of each cell type in each cluster. When a cluster is
# composed of more than a certain percentage (in this case, 10%) of a certain
# type, all the cells in the cluster are set to that type. If a cluster is
# composed of more than one cell type, the whole thing is marked "Ambiguous".
# If there's no cell type thats above the threshold, the cluster is marked
# "Unknown". Thus, Monocle helps you impute the type of each cell even in the
# presence of missing marker data.

expr <- clusterCells(
  expr, num_clusters = 2,
  frequency_thresh = 0.1, # > 10%
  cell_type_hierarchy = cth
)

plot_cell_clusters(
  expr, 1, 2, color = "CellType",
  markers = c("MYF5", "ANPEP")
)

# finally, we should pick out a cluster of a certain cell type, and run our
# trajectory and pseudo-time analysis then.
