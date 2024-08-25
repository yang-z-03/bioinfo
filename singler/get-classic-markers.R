
get_classic_markers <- function(
    ref, labels, assay.type = "logcounts", check.missing = TRUE, de.n = NULL,
    num.threads = BiocParallel::bpnworkers(BPPARAM),
    BPPARAM = BiocParallel::SerialParam()
  ) {

  if (!BiocParallel::bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
    BiocParallel::bpstart(BPPARAM)
    on.exit(BiocParallel::bpstop(BPPARAM))
  }

  if (!.is_list(ref)) {
    ref <- list(ref)
    labels <- list(labels)
  }

  # Setting up references.
  common <- Reduce(intersect, lapply(ref, rownames))
  if (length(common) == 0L && any(vapply(ref, nrow, 0L) > 0L)) {
    stop("no common row names across 'ref'")
  }
  common <- as.character(common) # avoid problems with NULL rownames for zero-row inputs.

  for (i in seq_along(ref)) {
    current <- ref[[i]][common, , drop = FALSE]
    current <- .to_clean_matrix(current, assay.type, check.missing,
                                msg = "ref", BPPARAM = BPPARAM)
    curptr <- beachmat::initializeCpp(current)

    flabels <- factor(labels[[i]])
    gm <- grouped_medians(curptr, as.integer(flabels) - 1L, 
                          nlevels(flabels), nthreads = num.threads)
    gm <- t(gm)
    colnames(gm) <- levels(flabels)
    ref[[i]] <- gm
  }

  ulabels <- .get_levels(unlist(lapply(ref, colnames)))
  labels <- list()
  for (i in seq_along(ref)) {
    m <- match(colnames(ref[[i]]), ulabels)
    labels[[i]] <- m - 1L
  }

  # identify top hits based on the average (or sum, it doesn't matter)
  # of the log-fold changes between labels across references.
  if (is.null(de.n)) {
    de.n <- -1L
  } else {
    stopifnot(de.n > 0)
  }

  out <- find_classic_markers(
    nlabels = length(ulabels),
    ngenes = length(common),
    labels = labels,
    ref = lapply(ref, beachmat::initializeCpp),
    de_n = de.n,
    nthreads = num.threads
  )

  names(out) <- ulabels
  for (i in seq_along(out)) {
    names(out[[i]]) <- ulabels
  }

  relist(common[unlist(out)], out)
}
