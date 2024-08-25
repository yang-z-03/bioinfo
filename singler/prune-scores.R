
prune_scores <- function(
    results, nmads = 3,
    min.diff.med = -Inf, min.diff.next = 0,
    get.thresholds = FALSE
  ) {

  delta <- delta_from_median(results)
  keep <- delta >= min.diff.med

  dn <- results$delta.next
  if (!is.null(dn)) {
    keep <- keep & dn >= min.diff.next
  }

  # ignoring the fine-tuning when allocating cells to labels.
  labels <- results$labels
  by.label <- split(which(keep), labels[keep])

  thresholds <- list()
  for (l in names(by.label)) {
    idx <- by.label[[l]]
    current <- delta[idx]

    med <- median(current)
    MAD <- mad(current, center = med)
    thresholds[[l]] <- curthresh <- med - nmads * MAD

    keep[idx] <- keep[idx] & current >= curthresh
  }

  if (get.thresholds) {
    unlist(thresholds)
  } else {
    !keep
  }
}
