
delta_from_median <- function(results) {
  scores <- results$scores
  labels <- results$labels
  assigned <- scores[cbind(seq_along(labels), match(labels, colnames(scores)))]
  assigned - DelayedMatrixStats::rowMedians(DelayedArray::DelayedArray(scores))
}
