
aggregate_reference <- function(
    ref, labels, ncenters = NULL, power = 0.5, ntop = 1000,
    assay.type = "logcounts", rank = 20, subset.row = NULL,
    check.missing = TRUE,
    BPPARAM = BiocParallel::SerialParam(), BSPARAM = BiocSingular::bsparam()
  ) {

  by.label <- split(seq_along(labels), labels)
  if (is(ref, "SummarizedExperiment")) {
    ref <- assay(ref, i = assay.type)
  }

  fragmenter <- function(x, i, j) {
    x <- x[, j, drop = FALSE] # this is usually smaller, so do this first.
    if (!is.null(i)) { x <- x[i, , drop = FALSE] }
    x
  }

  all.seeds <- .define_seeds(length(by.label))

  if (is.null(BPPARAM) ||
      is(BPPARAM, "BiocParallel::SerialParam") || 
      is(BPPARAM, "MulticoreParam")) {

    output.vals <- bpmapply(
      chosen = by.label, .seed = all.seeds,
      FUN = function(chosen, x, ...) {
        .aggregate_internal(fragmenter(x, i = subset.row, j = chosen), ...)
      },
      MoreArgs = list(
        x = ref, ncenters = ncenters, power = power, rank = rank,
        check.missing = check.missing, ntop = ntop, BSPARAM = BSPARAM
      ),
      BPPARAM = BPPARAM, SIMPLIFY = FALSE, USE.NAMES = FALSE
    )

  } else {
    by.mat <- lapply(by.label, fragmenter, x = ref, i = subset.row)
    output.vals <- bpmapply(by.mat,
      .seed = all.seeds, FUN = .aggregate_internal,
      MoreArgs = list(
        ncenters = ncenters, power = power, rank = rank,
        check.missing = check.missing, ntop = ntop, BSPARAM = BSPARAM
      ),
      BPPARAM = BPPARAM, SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
  }

  if (length(output.vals) == 0L) {
    output.vals[[1]] <- matrix(0, nrow(ref), 0, 
                               dimnames = list(rownames(ref), NULL))
  }

  first <- labels[vapply(by.label, function(i) i[1], 0L)]
  num <- vapply(output.vals, ncol, 0L)
  output.labels <- rep(first, num)

  output <- SummarizedExperiment(list(
    logcounts = do.call(cbind, output.vals)),
    colData = DataFrame(label = output.labels)
  )
  colnames(output) <- sprintf("%s.%s", output.labels, sequence(num))
  output
}

.aggregate_internal <- function(
    current, ncenters, power, rank, ntop, check.missing, BSPARAM, .seed
  ) {

  oldseed <- .get_seed()
  on.exit(.set_seed(oldseed))

  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1]), add = TRUE, after = FALSE)
  assign(".Random.seed", .seed, envir = .GlobalEnv)

  old.bp <- DelayedArray::getAutoBPPARAM()
  DelayedArray::setAutoBPPARAM(BiocParallel::SerialParam())
  on.exit(DelayedArray::setAutoBPPARAM(old.bp), add = TRUE, after = FALSE)

  cur.ncenters <- ncenters
  if (is.null(cur.ncenters)) {
    cur.ncenters <- floor(ncol(current)^power)
  }

  if (cur.ncenters <= 1) {
    val <- matrix(rowMeans(current), dimnames = list(rownames(current), NULL))
  } else if (cur.ncenters >= ncol(current)) {
    val <- current
  } else {
    to.use <- current

    if (ntop <= nrow(current)) {
      o <- order(DelayedMatrixStats::rowVars(to.use), decreasing = TRUE)
      to.use <- to.use[head(o, ntop), , drop = FALSE]
    }

    to.use <- t(to.use)
    to.use <- beachmat::realizeFileBackedMatrix(to.use)

    # Identifying the top PCs to avoid realizing the entire matrix.
    if (rank <= min(dim(to.use)) - 1L) {
      pcs <- BiocSingular::runPCA(to.use, rank = rank,
                    get.rotation = FALSE, BSPARAM = BSPARAM) $ x
    } else {
      pcs <- as.matrix(to.use)
    }

    clustered <- kmeans(pcs, centers = cur.ncenters)
    val <- DelayedArray::colsum(DelayedArray::DelayedArray(current),
                                clustered $ cluster)
    tab <- table(clustered $ cluster)[colnames(val)]
    val <- sweep(val, 2, tab, "/")
  }

  val
}

.define_seeds <- function(n) {
  if (!n) {
    return(list())
  }

  # bumping the RNG so that repeated calls to this function
  # generate different results.
  runif(1)

  oldseed <- .get_seed()
  on.exit(.set_seed(oldseed))

  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1]), add = TRUE, after = FALSE)

  seeds <- vector("list", n)
  seeds[[1L]] <- .Random.seed
  for (i in seq_len(n - 1L)) {
    seeds[[i + 1L]] <- parallel::nextRNGStream(seeds[[i]])
  }

  seeds
}

.get_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

.set_seed <- function(oldseed) {
  if (!is.null(oldseed)) {
    assign(".Random.seed", oldseed, envir = .GlobalEnv)
  } else {
    rm(.Random.seed, envir = .GlobalEnv)
  }
}
