
singler <- function(
    test, ref, labels, method = NULL,
    clusters = NULL,
    genes = "de", sd.thresh = 1,
    de.method = "classic", de.n = NULL, de.args = list(),
    aggr.ref = FALSE, aggr.args = list(),
    recompute = TRUE,
    restrict = NULL,
    quantile = 0.8,
    fine.tune = TRUE, tune.thresh = 0.05, prune = TRUE,
    assay.type.test = "logcounts", assay.type.ref = "logcounts",
    check.missing = TRUE,
    num.threads = BiocParallel::bpnworkers(BPPARAM), BNPARAM = NULL,
    BPPARAM = BiocParallel::SerialParam()) {

  if (!BiocParallel::bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
    BiocParallel::bpstart(BPPARAM)
    on.exit(BiocParallel::bpstop(BPPARAM))
  }

  test <- .to_clean_matrix(test, assay.type.test, check.missing, msg = "test",
                           BPPARAM = BPPARAM)

  # Converting to a common list format for ease of data munging.
  if (single.ref <- !.is_list(ref)) {
    ref <- list(ref)
  }

  ref <- lapply(ref,
    FUN = .to_clean_matrix, assay.type = assay.type.ref,
    check.missing = check.missing, msg = "ref", BPPARAM = BPPARAM
  )
  refnames <- Reduce(intersect, lapply(ref, rownames))

  keep <- intersect(rownames(test), refnames)
  if (length(keep) == 0) {
    stop("no common genes between 'test' and 'ref'")
  }
  if (!identical(keep, rownames(test))) {
    test <- test[keep, ]
  }

  for (i in seq_along(ref)) {
    if (!identical(keep, rownames(ref[[i]]))) {
      ref[[i]] <- ref[[i]][keep, , drop = FALSE]
    }
  }

  # Converting back.
  if (single.ref) {
    ref <- ref[[1]]
  }

  trained <- train_singler(
    ref, labels, genes = genes,
    sd.thresh = sd.thresh,
    de.method = de.method,
    de.n = de.n,
    de.args = de.args,
    aggr.ref = aggr.ref,
    aggr.args = aggr.args,
    recompute = recompute,
    restrict = restrict,
    check.missing = FALSE,
    BNPARAM = BNPARAM,
    num.threads = num.threads,
    BPPARAM = BPPARAM
  )

  if (!is.null(method)) {
    .Deprecated(msg = paste(
      "method = 'cluster' is no longer necessary when argument 'cluster'",
      "is set to non-null value.", sep = " "))
  }

  if (!is.null(clusters)) {
    oldp <- DelayedArray::getAutoBPPARAM()
    DelayedArray::setAutoBPPARAM(BPPARAM)
    on.exit(DelayedArray::setAutoBPPARAM(oldp), add = TRUE)

    if (test |> dim() |> is.null()) {
      # here, test should be a matrix. however, it may shrink into an named
      # vector when their is no genes in the reference marker. we should
      # post an error message -- yang-z, aug. 17
      stop("there is no overlapping genes between your data and the reference!")
    }
    test <- DelayedArray::colsum(DelayedArray::DelayedArray(test), clusters)
  }

  classify_singler(
    test, trained,
    quantile = quantile,
    fine.tune = fine.tune,
    tune.thresh = tune.thresh,
    prune = prune,
    check.missing = FALSE,
    num.threads = num.threads,
    BPPARAM = BPPARAM
  )
}
