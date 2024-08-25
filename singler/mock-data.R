
.mock_ref_data <- function(ngroups = 5, nreps = 4, ngenes = 1000, prop = 0.5) {
  nmarkers <- ngenes * prop
  nmarkers.per.group <- ceiling(nmarkers / ngroups)

  means <- matrix(0,
    ncol = ngroups, nrow = ngenes,
    dimnames = list(
      sprintf("GENE_%i", seq_len(ngenes)),
      LETTERS[seq_len(ngroups)]
    )
  )

  counter <- 0L
  for (i in seq_len(ngroups)) {
    means[counter + seq_len(nmarkers.per.group), i] <- rnorm(nmarkers.per.group)
    counter <- counter + nmarkers.per.group
  }

  g <- rep(colnames(means), each = nreps)
  mat <- matrix(rpois(ngenes * length(g),
                      lambda = 10 * 2 ^ means[, g]),
                ncol = length(g))

  rownames(mat) <- rownames(means)

  SummarizedExperiment(
    list(counts = mat),
    colData = DataFrame(label = g),
    metadata = list(means = means)
  )
}
