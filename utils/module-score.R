
# features are lists of vectors.
# this function accepts log-normed data (but not scaled) as input (slot 'data').

mod.score <- function(
    assay.data, features,
    pool = NULL,
    nbin = 24, ctrl = 100,
    name = 'cluster', seed = 1
) {
  
  # set the random seed.
  if (!is.null(x = seed)) set.seed(seed = seed)
  features.old <- features
  
  if (is.null(x = features))
    stop("missing input feature list")
  cluster.length <- length(x = features)
  
  pool <- pool %||% rownames(assay.data)
  data.avg <- Matrix::rowMeans(assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  
  # add a random small value, to ensure the average data should not be
  # completely identical.
  
  data.cut <- ggplot2::cut_number(
    x = data.avg + rnorm(n = length(data.avg)) / 1e30,
    n = nbin, labels = FALSE, right = FALSE
  )
  
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  
  # determine the control features: we will first split the genes by expression
  # levels, and we will find for each module feature, a set of control genes
  # with the same level of expression (sized ctrl = 100). we may happen to include
  # the gene itself as its reference.
  
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(sample(
          data.cut[which(data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(ctrl.use), # number of scores
    ncol = ncol(assay.data) # number of cells
  )
  
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    # calculate the mean log-normed data (of control genes) for each cell
    ctrl.scores[i, ] <- Matrix::colMeans(assay.data[features.use, ])
  }
  
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(assay.data)
  )
  
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    # features.use may only be 1.
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- names(features) %||% paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(assay.data)
  
  return (features.scores.use)
}
