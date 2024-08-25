
combine_common <- function(results) {
  .Deprecated(new = "combine_recomputed")

  if (length(unique(lapply(results, rownames))) != 1) {
    stop("cell/cluster names in 'results' are not identical")
  }
  if (length(unique(vapply(results, nrow, 0L))) != 1) {
    stop("numbers of cells/clusters in 'results' are not identical")
  }

  all.common <- lapply(results, function(x) sort(metadata(x) $ common.genes))
  if (length(unique(all.common)) != 1) {
    # This should be changed to 'stop' before release/after merge with PR #60.
    warning("common genes are not identical")
  }

  ncells <- nrow(results[[1]])
  collected.scores <- collected.best <- vector("list", length(results))
  for (i in seq_along(results)) {
    scores <- results[[i]] $ scores
    collected.best[[i]] <- scores[cbind(seq_len(ncells), max.col(scores))]
    collected.scores[[i]] <- scores
  }

  all.scores <- do.call(cbind, collected.scores)
  output <- DataFrame(scores = I(all.scores),
                      row.names = rownames(results[[1]]))

  metadata(output) $ common.genes <- all.common[[1]]
  metadata(output) $ label.origin <- .create_label_origin(collected.scores)

  chosen <- max.col(do.call(cbind, collected.best))
  cbind(output, .combine_result_frames(chosen, results))
}
