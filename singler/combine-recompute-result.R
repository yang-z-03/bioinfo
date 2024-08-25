
combine_recomputed <- function(
    results,
    test,
    trained,
    quantile = 0.8,
    assay.type.test = "logcounts",
    check.missing = TRUE,
    allow.lost = FALSE,
    warn.lost = TRUE,
    num.threads = BiocParallel::bpnworkers(BPPARAM),
    BPPARAM = BiocParallel::SerialParam()
  ) {
  
  all.names <- c(list(colnames(test)), lapply(results, rownames))
  if (length(unique(all.names)) != 1) {
    stop("cell/cluster names in 'results' are not identical")
  }
  all.nrow <- c(ncol(test), vapply(results, nrow, 0L))
  if (length(unique(all.nrow)) != 1) {
    stop("numbers of cells/clusters in 'results' are not identical")
  }

  # Checking the marker consistency.
  all.refnames <- lapply(trained, function(x) rownames(x $ ref))
  intersected <- Reduce(intersect, all.refnames)
  for (i in seq_along(trained)) {
    if (!all(trained[[i]] $ markers $ unique %in% rownames(test))) {
      stop("all markers stored in 'results' should be present in 'test'")
    } else if (warn.lost && !all(trained[[i]] $ markers $ unique %in% intersected)) {
      warning("entries of 'trained' differ in the universe of available markers")
    }
  }

  # Applying the integration.
  universe <- Reduce(union, c(list(rownames(test)), all.refnames))
  ibuilt <- integrate_build(
    match(rownames(test), universe) - 1L,
    lapply(trained, function(x) beachmat::initializeCpp(x $ ref)),
    lapply(trained, function(x) match(rownames(x $ ref), universe) - 1L),
    lapply(trained, function(x) match(x $ labels $ full,
                                      x $ labels $ unique) - 1L),
    lapply(trained, function(x) x $ built),
    nthreads = num.threads
  )

  test <- .to_clean_matrix(test, assay.type = assay.type.test,
                           check.missing = check.missing, msg = "test",
                           BPPARAM = BPPARAM)

  collated <- vector("list", length(trained))
  for (i in seq_along(collated)) {
    collated[[i]] <- match(results[[i]] $ labels,
                           trained[[i]] $ labels $ unique) - 1L
  }

  parsed <- beachmat::initializeCpp(test)
  irun <- integrate_run(parsed, collated, ibuilt,
                        quantile = quantile, nthreads = num.threads)
  scores <- irun $ scores

  # Organizing the outputs.
  base.scores <- vector("list", length(results))
  for (r in seq_along(base.scores)) {
    mat <- results[[r]] $ scores
    mat[] <- NA_real_
    idx <- cbind(seq_len(nrow(mat)), collated[[r]] + 1L)
    mat[idx] <- scores[, r]
    base.scores[[r]] <- mat
  }

  all.scores <- do.call(cbind, base.scores)
  output <- DataFrame(scores = I(all.scores), row.names = rownames(results[[1]]))
  S4Vectors::metadata(output) $ label.origin <- .create_label_origin(base.scores)

  chosen <- irun $ best + 1L
  cbind(output, .combine_result_frames(chosen, results))
}

#' @importFrom S4Vectors DataFrame
.combine_result_frames <- function(chosen, results) {
  has.pruned <- !is.null(results[[1]] $ pruned.labels)

  # Organizing the statistics based on the chosen results.
  chosen.label <- chosen.first <- chosen.pruned <- rep(NA_character_, nrow(results[[1]]))

  for (u in unique(chosen)) {
    current <- chosen == u
    res <- results[[u]]
    chosen.label[current] <- res $ labels[current]

    if (has.pruned) { # same for pruned.
      chosen.pruned[current] <- res $ pruned.labels[current]
    }
  }

  output <- DataFrame(labels = chosen.label, row.names = rownames(results[[1]]))

  if (has.pruned) {
    output $ pruned.labels <- chosen.pruned
  }

  output $ reference <- chosen

  if (is.null(names(results))) {
    names(results) <- sprintf("ref%i", seq_along(results))
  }

  output $ orig.results <- do.call(DataFrame, lapply(results, I))

  output
}

#' @importFrom S4Vectors DataFrame
.create_label_origin <- function(collected.scores) {
  DataFrame(
    label = unlist(lapply(collected.scores, colnames)),
    reference = rep(
      seq_along(collected.scores), 
      vapply(collected.scores, ncol, 0L)
    )
  )
}
