
classify_singler <- function(
    test,
    trained,
    quantile = 0.8,
    fine.tune = TRUE,
    tune.thresh = 0.05,
    sd.thresh = NULL,
    prune = TRUE,
    assay.type = "logcounts",
    check.missing = TRUE,
    num.threads = BiocParallel::bpnworkers(BPPARAM),
    BPPARAM = BiocParallel::SerialParam()
  ) {

  test <- .to_clean_matrix(test, assay.type, check.missing,
                           msg = "test", BPPARAM = BPPARAM)

  solo <- .is_solo(trained)
  if (solo) {
    trained <- list(trained)
  }

  results <- lapply(trained,
    FUN = .classify_internals,
    test = test,
    quantile = quantile,
    fine.tune = fine.tune,
    tune.thresh = tune.thresh,
    prune = prune,
    num.threads = num.threads
  )

  if (solo) {
    results[[1]]
  } else {
    combine_recomputed(
      results,
      test = test,
      trained = trained,
      check.missing = FALSE,
      quantile = quantile
    )
  }
}

.classify_internals <- function(test, trained, quantile, fine.tune,
                                tune.thresh = 0.05, prune = TRUE,
                                num.threads = 1) {

  m <- match(trained $ markers $ unique, rownames(test))
  if (anyNA(m)) {
    stop("'rownames(test)' does not contain all genes used in 'trained'")
  }

  trained <- rebuild_index(trained, num.threads = num.threads)

  parsed <- beachmat::initializeCpp(test)
  out <- run(parsed, m - 1L, trained $ built,
    quantile = quantile,
    use_fine_tune = fine.tune,
    fine_tune_threshold = tune.thresh,
    nthreads = num.threads
  )

  colnames(out $ scores) <- trained $ labels $ unique
  output <- DataFrame(
    scores = I(out $ scores),
    labels = trained $ labels $ unique[out $ best + 1L],
    delta.next = out $ delta,
    check.names = FALSE
  )

  if (prune) {
    output $ pruned.labels <- output $ labels
    output $ pruned.labels[prune_scores(output)] <- NA_character_
  }

  rownames(output) <- colnames(test)
  metadata(output) $ common.genes <- trained $ markers $ unique
  metadata(output) $ de.genes <- trained $ markers $ full

  output
}
