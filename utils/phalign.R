
# the inputs should follow normal distribution, since we will run pca analysis
# directly, we will expect 0-centered scaled data matrix.

phalign <- function (
    ref, query, hvgs,
    ref.tag, query.cluster,
    ndim.pca = 50,
    n.neighbors = 50, valid.nn.r = 0.50, valid.nn.p = 0.05,
    disting.gap = 0.01, top.tag = 3,
    verbose = TRUE,
    permutation.test = FALSE, shuffle = 100
) {
  
  ref.tag <- as.character(ref.tag)
  query.cluster <- as.character(query.cluster)
  
  # suppose the rows are genes and columns are cells, and that the gene and cell
  # matrix are well ordered.
  
  q.gnames <- rownames(query)
  r.gnames <- rownames(ref)
  intersects <- r.gnames %in% q.gnames
  if (verbose) {
    cat('found', sum(intersects), 'genes (out of', length(r.gnames), 
        ') matching between query and reference.', '\n')
  }
  
  # make sure that the two matrix is in order.
  
  filt.genes <- r.gnames[intersects]
  query <- query[filt.genes, ]
  ref <- ref[filt.genes, ]
  
  # make sure that the two matrix do not have sharing cells.
  
  q.cnames <- colnames(query)
  r.cnames <- colnames(ref)
  matches <- match(q.cnames, r.cnames)
  dup <- (!(matches |> is.na())) |> sum()
  if (dup > 0) {
    cat('there are', dup, 'cells identical from query',
        'set to the reference set. this should not happen, you should check',
        'your data and make sure the cells are unique.')
    return(NULL)
  }
  
  merged.scaled <- cbind(query, ref)
  
  # here, we should only use highly-variable genes to perform the pca analysis
  # since this will save a great many time. notice that the prcomp function regards
  # rows as observations (cells), but not genes.
  
  hvg.scaled <- merged.scaled[hvgs, ]
  pca <- prcomp(hvg.scaled |> t(), scale. = FALSE) # we have already scaled it.
  
  # and we only utilize the first several components of the pca.
  
  pcs <- pca $ x[, 1:ndim.pca]
  pc.ref <- pcs[r.cnames, ]
  pc.query <- pcs[q.cnames, ]
  
  # query the k nearest neighbors within the pcs of reference.
  
  knn <- BiocNeighbors::queryKNN(pc.ref, pc.query, k = n.neighbors)
  
  # for each cell's neighbors (here n = 50 by default), we will calculate the
  # heterogeneity within them using gene-space linear correlation.
  
  corrs <- matrix(nrow = nrow(pc.query), ncol = n.neighbors)
  corrs.p <- matrix(nrow = nrow(pc.query), ncol = n.neighbors)
  knn.i <- matrix(nrow = nrow(pc.query), ncol = n.neighbors)
  knn.c <- matrix(nrow = nrow(pc.query), ncol = n.neighbors)
  for (i in 1:nrow(pc.query)) {
    nns <- knn $ index[i, ]
    neigh.corr <- rep(0, length(nns))
    neigh.corr.p <- rep(0, length(nns))
    for (j in seq_along(nns)) {
      neigh.corr[j] <- cor(pc.query[i, ], pc.ref[nns[j], ], method = 'pearson')
      neigh.corr.p[j] <- cor.test(pc.query[i, ], pc.ref[nns[j], ], method = 'pearson') $ p.value
    }
    
    # plot(pc.query[1, ], pc.ref[nns[1], ])
    
    # for those with principles components not linear correlates well. it suggest
    # that the two points (although picked out as neighbors), actually varies
    # greatly in the transcriptome landscape. these are not actual neighbors,
    # just picking the best from all the bads.
    
    corrs[i, ] <- neigh.corr
    corrs.p[i, ] <- neigh.corr.p
    
    # what is adjusted p here? why adjust?
    filt <- (neigh.corr > valid.nn.r) & (neigh.corr.p < valid.nn.p)
    nns[!filt] <- NA
    knn.i[i, ] <- nns
    knn.c[i, ] <- rep(i, n.neighbors)
  }
  
  # now, we will query the tag distribution in reference of selected nearest
  # neighbors, invalid neighbors are set to NA in index.
  
  # convert the matrix to long form data.
  # some cells have left without a single valid neighbors, and these cells will
  # be completely thrown out in the later analyses. they will not present in the
  # long form data, as well as all others.
  
  na.nums <- knn.i |> is.na() |> rowSums()
  
  # hist(na.nums)
  
  longf <- data.frame(
    cell.index = as.vector(knn.c |> t()),
    nn.index = as.vector(knn.i |> t()),
    corr = as.vector(corrs |> t()),
    p.val = as.vector(corrs.p |> t())
  )
  
  real <- phalign.stats(
    query, longf, ref.tag, query.cluster, top.tag, na.nums, n.neighbors, disting.gap,
    verbose = TRUE
  )
  
  order.row <- table(query.cluster) |> names()
  order.col <- table(ref.tag) |> names()
  result.matrix = matrix(nrow = length(order.row), ncol = length(order.col))
  rownames(result.matrix) <- order.row
  colnames(result.matrix) <- paste('label', order.col, sep = '.')
  result.matrix[
    real $ cluster |> rownames(),
    real $ cluster |> colnames()
  ] <- real $ cluster
  
  if (! permutation.test) {
    return (list(
      cell = real $ cell,
      cluster = result.matrix
    ))
  }
  
  # permute the ref.tag for a certain time, to generate an array of null distributions
  
  perm.matrix = array(dim = c(length(order.row), length(order.col), shuffle))
  
  if (verbose) pbar <- progress::progress_bar $ new(total = shuffle)
  
  for (z in 1:shuffle) {
    perm <- sample(ref.tag)
    cls <- phalign.stats(
      query, longf, perm, query.cluster, top.tag, na.nums, n.neighbors, disting.gap,
      verbose = FALSE
    ) $ cluster
    
    pm = matrix(nrow = length(order.row), ncol = length(order.col))
    rownames(pm) <- order.row
    colnames(pm) <- paste('label', order.col, sep = '.')
    pm[
      cls |> rownames(),
      cls |> colnames()
    ] <- cls
    
    perm.matrix[, , z] <- pm
    if (verbose) pbar $ tick()
  }
  
  res <- list(
    cell = real $ cell,
    cluster = result.matrix,
    permutation = perm.matrix
  )
  
  res <- phalign.significance(res)
  return (res)
}


# internal function, reused during permutation test.

phalign.stats <- function(
    query, longf,
    ref.tag, query.cluster, 
    top.tag, na.nums, n.neighbors = 50, disting.gap = 0.01,
    verbose = TRUE
) {
  
  longf $ tag <- ref.tag[longf $ nn.index]
  
  # discrimination test: tell if the cell match uniformly to all the tags,
  # if the correlation values for each type of tags is not discriminating, this
  # cell identity is thus undetermined. they just ignore these cells in the
  # following procedures.
  
  # they indicates that we should use the single maximum correlation when removing
  # cells, but we now choose the top-n cell's median.
  
  longf <- longf[! is.na(longf $ nn.index), ] |> tibble::tibble()
  
  # extract top-n neighbors aligned to each label:
  
  top.n.mean <- longf |>
    dplyr::group_by(cell.index, tag) |>
    dplyr::top_n(corr, n = top.tag) |>
    dplyr::summarise(mean.top = mean(corr), .groups = 'keep')
  
  top.1 <- longf |>
    dplyr::group_by(cell.index, tag) |>
    dplyr::top_n(corr, n = 1)
  
  top.1 $ uid <- 1:nrow(top.1)
  
  # boxplot(corr ~ tag, data = subset(longf, subset = cell.index == 1))
  # boxplot(corr ~ tag, data = subset(top.n, subset = cell.index == 1))
  
  top.case <- top.1 |>
    dplyr::select(uid, cell.index, tag, corr) |>
    dplyr::group_by(cell.index) |>
    dplyr::top_n(corr, n = 1)
  
  n.identical.top <- duplicated(top.case $ cell.index) |> sum()
  dup.cell.id <- unique(top.case $ cell.index[duplicated(top.case $ cell.index)])
  if (n.identical.top > 0 && verbose) {
    cat('there are', n.identical.top, 'cells with identical maximum correlation',
        'scores between reference tags. these cells are prone to be ambiguous.', '\n')
  }
  
  # remove identical cells
  
  top.2 <- top.1[!(top.1 $ uid %in% top.case $ uid), ]
  top.case <- top.case[!(top.case $ cell.index %in% dup.cell.id), ]
  
  # ensure that cells have two possible tags, if the cell only have one. it is
  # aligned exactly to that tag.
  
  second.case <- top.2 |>
    dplyr::select(uid, cell.index, tag, corr) |>
    dplyr::group_by(cell.index) |>
    dplyr::top_n(corr, n = 1)
  
  # have a secondary cell index
  unique.map <- !(top.case $ cell.index %in% second.case $ cell.index)
  if (unique.map |> sum() > 0 && verbose) {
    cat('there are', unique.map |> sum(), 'cells with exact match', '\n')
  }
  
  unique.id <- top.case[unique.map, ]
  
  # non-unique mappings will then calculate a gap metric.
  
  nu.1 <- top.case[!unique.map, ]
  nu.2 <- second.case
  colnames(nu.1) <- c('uid.1', 'cell.index', 'tag.1', 'corr.1')
  colnames(nu.2) <- c('uid.2', 'cell.index', 'tag.1', 'corr.2')
  
  # assertion here:
  # identical(nu.1 $ cell.index, nu.2 $ cell.index)
  
  nu <- merge(nu.1, nu.2, by = 'cell.index')
  nu $ delta <- nu $ corr.1 - nu $ corr.2
  
  # finally, we will construct the cell table:
  
  wide <- tidyr::pivot_wider(
    top.n.mean, names_from = 'tag', names_prefix = 'tag.', values_from = 'mean.top'
  ) |> as.data.frame()
  
  rownames(wide) <- wide $ cell.index
  cell.table <- data.frame(cell.index = 1:ncol(query))
  cell.table <- merge(cell.table, wide, by = 'cell.index', all.x = TRUE)
  cell.table $ .cluster <- query.cluster
  cell.table $ .status <- rep('.', ncol(query))
  cell.table $ .label <- rep('.', ncol(query))
  
  # no valid neighbors status.
  
  cell.table $ .na.nums <- na.nums
  cell.table $ .status[na.nums == n.neighbors] <- 'no.valid.nn'
  
  # unique matching status
  
  cell.table $ .status[unique.id $ cell.index |> unique()] <- 'unique'
  unique.id <- unique.id[!duplicated(unique.id $ cell.index), ]
  cell.table $ .label[unique.id $ cell.index] <- unique.id $ tag
  
  # not significant cells
  
  notsig <- nu $ cell.index[nu $ delta < disting.gap]
  cell.table $ .status[notsig] <- 'not.sig'
  cell.table $ .label[nu $ cell.index] <- nu $ tag.1.x
  cell.table $ .label[notsig] <- '.'
  
  # by now, we have constructed the cell table successfully. the only thing left
  # is to summarize this cell table according to the cluster, and calculate
  # the proportion of cluster identity.
  
  cluster.ident <- cell.table[cell.table $ .label != '.', c('.cluster', '.label')]
  stats <- cluster.ident |> tibble() |>
    dplyr::group_by(.cluster, .label) |>
    dplyr::summarise(n = dplyr::n(), .groups = 'keep') |>
    tidyr::pivot_wider(
      names_from = '.label', names_prefix = 'label.', values_from = 'n'
    ) |> as.data.frame()
  
  rownames(stats) <- stats $ .cluster
  stats $ .cluster <- NULL
  stats.sum <- rowSums(stats, na.rm = TRUE)
  stats <- as.matrix(stats) / stats.sum
  
  # finalize
  
  cell.table $ cell.name <- colnames(query)
  rownames(cell.table) <- cell.table $ cell.name
  
  return (list(
    cell = cell.table,
    cluster = stats
  ))
}

phalign.significance <- function(phalign.res) {
  dims <- dim(phalign.res $ permutation)
  pvals <- matrix(nrow = dims[1], ncol = dims[2])
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      null.dist <- phalign.res $ permutation[i, j, ]
      x <- phalign.res $ cluster[i, j]
      pvals[i, j] <- 1 - pnorm(
        x, mean = mean(null.dist, na.rm = TRUE), 
        sd = stats::sd(null.dist, na.rm = TRUE)
      )
    }
  }
  
  phalign.res $ p.values <- pvals
  return (phalign.res)
}

phalign.plot.single <- function(res, cluster, tag) {
  cluster.id <- match(cluster, rownames(res $ cluster))
  label.id <- match(paste('label', tag, sep = '.'), colnames(res $ cluster))
  
  plot.new()
  pval <- res $ p.values[cluster.id, label.id]
  if (pval > 0.5) pval <- 1 - pval
  hist(
    res $ permutation[cluster.id, label.id, ],
    main = paste('p =', pval),
    xlab = paste('Cluster:', cluster, 'Tag:', tag)
  )
  abline(v = res $ cluster[cluster.id, label.id], col = 'red', lwd = 3, lty = 2)
}

phalign.plot.cluster <- function(res, cluster) {
  cluster.id <- match(cluster, rownames(res $ cluster))
  lab.names <- colnames(res $ cluster)
  nlab <- lab.names |> length()
  
  .label <- c()
  .perm <- c()
  .res.label <- c()
  .res.value <- c()
  .res.p <- c()
  for (i in 1:nlab) {
    values <- res $ permutation[cluster.id, i, ]
    label.name <- lab.names[i] |> stringr::str_replace('label.', '')
    .label <- c(.label, rep(label.name, length(values)))
    .perm <- c(.perm, values)
    .res.label <- c(.res.label, label.name)
    .res.value <- c(.res.value, res $ cluster[cluster.id, i])
    .res.p <- c(.res.p, res $ p.values[cluster.id, i])
  }
  
  df <- data.frame(label = .label, perm = .perm)
  df.points <- data.frame(label = .res.label, value = .res.value, p = .res.p)
  
  df <- df[!is.na(df $ perm), ]
  df.points <- df.points[!is.na(df.points $ value), ]
  
  df.points[df.points $ p > 0.5, ] $ p <- 1 - df.points[df.points $ p > 0.5, ] $ p
  df.points $ p.text <- rep('', nrow(df.points))
  if (df.points[df.points $ p < 0.10, ] |> nrow() > 0)
    df.points[df.points $ p < 0.10, ] $ p.text <- ''
  if (df.points[df.points $ p < 0.05, ] |> nrow() > 0)
    df.points[df.points $ p < 0.05, ] $ p.text <- '*'
  if (df.points[df.points $ p < 0.01, ] |> nrow() > 0)
    df.points[df.points $ p < 0.01, ] $ p.text <- '**'
  if (df.points[df.points $ p < 0.001, ] |> nrow() > 0)
    df.points[df.points $ p < 0.001, ] $ p.text <- '***'
  
  ggplot2::ggplot(ggplot2::aes(y = perm, x = label), data = df) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(0.1)) +
    ggplot2::geom_point(ggplot2::aes(y = value, x = label), data = df.points, color = 'red') +
    ggplot2::geom_text(
      ggplot2::aes(y = value + 0.02, x = label, label = p.text),
      data = df.points, color = 'red', size = 10
    ) +
    ggplot2::labs(
      y = 'Similarity Score (PhenoAligner)',
      x = 'Label', title = cluster
    )
}

# dec 16, 2024
# yang-z. a rewrite of phenoaligner algorithm (zhang, 2022)
